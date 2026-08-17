"""Generate word clouds for TCR/AIRR CDR3 datasets.

This module provides a CLI-backed `wordcloud(args)` entrypoint and several
helpers used to create a wordcloud grouped by chain + repertoire.

The implementation is kept intentionally small and exposes helpers to make
unit testing feasible without requiring matplotlib rendering.
"""

import json
from argparse import Namespace
from pathlib import Path
from typing import Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from natsort import natsorted
from wordcloud import WordCloud

import tcrcloud.colours
import tcrcloud.format
from tcrcloud.errors import TCRcloudError

# No fixed species default here; per-chain color lookup is dynamic via tcrcloud.colours helpers.


def _vcall_color(
    vcall: str, species: str = "homo_sapiens", default: str = "grey"
) -> str:
    """Resolve a V-gene call to a color string using the requested species palette."""
    return tcrcloud.colours.get_color_from_vcall(vcall, species, default)


# This colours the wordclouds
class SimpleGroupedColorFunc:
    """Create a color function object which assigns EXACT colors
    to certain words based on the color to words mapping

    Parameters
    ----------
    color_to_words : dict(str -> list(str))
        A dictionary that maps a color to the list of words.

    default_color : str
        Color that will be assigned to a word that's not a member
        of any value from color_to_words.
    """

    def __init__(self, color_to_words: dict[str, list[str]], default_color: str):
        self.word_to_color = {
            word: color for (color, words) in color_to_words.items() for word in words
        }

        self.default_color = default_color

    def __call__(self, word: str, **kwargs: Any) -> str:
        return self.word_to_color.get(word, self.default_color)


def handle_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure CDR3 sequences are unique for wordcloud generation.

    WordCloud requires unique word keys. If the input dataframe contains
    repeated `junction_aa` sequences, it would collapse them into a single
    word with combined frequencies.

    To avoid that, we prepend increasing numbers of spaces to duplicates.
    For example, if `CASSIRSSYEQYF` occurs three times, it becomes:
    - "CASSIRSSYEQYF"
    - " CASSIRSSYEQYF"
    - "  CASSIRSSYEQYF"
    """

    duplicated = df["junction_aa"].duplicated(keep=False)
    if not duplicated.any():
        return df

    # Prepend a growing number of spaces to each duplicate to ensure uniqueness.
    # The first occurrence of a repeated sequence keeps its original spelling
    # (0 spaces); each subsequent occurrence gets one more leading space.
    counts = df.loc[duplicated].groupby("junction_aa").cumcount()

    # `counts` is a pandas Series, so we use a Python-level mapping to avoid
    # pandas/numpy string ufunc multiplication issues.
    prefixes = counts.apply(lambda c: " " * int(c))
    df.loc[duplicated, "junction_aa"] = prefixes + df.loc[duplicated, "junction_aa"]
    return df


def _ensure_required_columns(df: pd.DataFrame) -> None:
    required_columns = {"junction_aa", "v_call", "counts", "chain", "repertoire_id"}
    missing = required_columns - set(df.columns)
    if missing:
        raise TCRcloudError(f"missing required columns: {sorted(missing)}")


def _extract_family_and_text(df: pd.DataFrame) -> tuple[dict, dict]:
    """Extract the mapping used to build the wordcloud.

    WordCloud wants a mapping of word->weight. Here, `junction_aa` is treated
    as the word and `counts` is treated as the weight.

    Additionally, we keep a `family` mapping that remembers which V-gene call is
    associated with each `junction_aa`, so that we can colour words consistently.

    `iloc[:, 0]` (not the deprecated `DataFrame.squeeze()`) selects the single
    remaining column, and works uniformly for groups of one or many rows.
    """

    family = df[["junction_aa", "v_call"]].set_index("junction_aa").iloc[:, 0].to_dict()
    text = df[["junction_aa", "counts"]].set_index("junction_aa").iloc[:, 0].to_dict()
    return family, text


def _load_colour_mapping(colours_path: str) -> dict[str, list[str]]:
    try:
        with open(colours_path) as json_file:
            return json.load(json_file)
    except FileNotFoundError as exc:
        raise TCRcloudError(f"{colours_path} doesn't seem to exist") from exc
    except json.decoder.JSONDecodeError as exc:
        raise TCRcloudError(
            f"{colours_path} doesn't seem properly formatted. Check https://github.com/eriicdesousa/TCRcloud for more information"
        ) from exc


def _build_color_to_words(
    family: dict, colours_path: str | None, species: str = "homo_sapiens"
) -> dict:
    """Build mapping from colours -> list of words for WordCloud.

    If `colours_path` is set, the JSON file is treated as the source of truth.
    Otherwise, the mapping is derived from the V-gene calls.

    This function exists so unit tests can validate both behaviour modes.
    """

    if colours_path is not None:
        return _load_colour_mapping(colours_path)

    color_to_words: dict[str, list[str]] = {}
    for aa, vcall in family.items():
        colour = _vcall_color(vcall, species)
        color_to_words.setdefault(colour, []).append(aa)
    return color_to_words


def _add_legend(colour_map: dict[str, str]) -> None:
    # Natural-sort so gene numbers order numerically (TRBV2 before TRBV10).
    sorted_legend = natsorted(colour_map)

    patch_list = [
        mpatches.Patch(color=colour_map[key], label=key) for key in sorted_legend
    ]

    plt.legend(
        handles=patch_list,
        bbox_to_anchor=(0.5, -0.01),
        loc="upper center",
        ncol=4,
        prop={"size": 6},
    )


def wordcloud(args: Namespace) -> None:
    """Main entrypoint for the `TCRcloud cloud` command.

    This function is intentionally small; it delegates most of the work to
    helpers so that tests can cover behaviour without rendering plots.
    """

    # Convert legacy string boolean values into a real boolean.
    legend = args.legend
    if isinstance(legend, str):
        if legend.lower() in ("true", "t", "1", "yes", "y"):
            legend = True
        elif legend.lower() in ("false", "f", "0", "no", "n"):
            legend = False
        else:
            raise TCRcloudError("please indicate legend True or False")

    legend = bool(legend)

    # Determine the output image format (defaults to "png"). Validated here
    # too since `wordcloud()` may be called directly (e.g. in tests) rather
    # than only via the argparse CLI, whose `choices=["svg", "png"]` wouldn't
    # otherwise catch a bad value.
    output_format = (getattr(args, "format", None) or "png").strip().lower()
    if output_format not in ("svg", "png"):
        raise TCRcloudError(
            f"unsupported output format '{output_format}'. "
            "Please choose 'svg' or 'png'"
        )

    # Format and validate the input AIRR CDR3 data.
    samples_df = tcrcloud.format.format_data(args)
    formatted_samples = tcrcloud.format.format_cloud(samples_df)
    _ensure_required_columns(formatted_samples)

    if not formatted_samples["junction_aa"].is_unique:
        handle_duplicates(formatted_samples)

    # Use the base name of the input file to generate output filenames.
    input_stem = Path(args.rearrangements).stem

    for (chain, repertoire_id), df in formatted_samples.groupby(
        ["chain", "repertoire_id"]
    ):
        family, text = _extract_family_and_text(df)

        # Build the wordcloud object using the token frequency map.
        wordcloud_obj = WordCloud(
            width=1000,
            height=args.size,
            background_color="white",
            relative_scaling=0.7,
            prefer_horizontal=1.0,
            scale=2.5,
            max_font_size=3000,
            max_words=len(df),
        ).generate_from_frequencies(text)

        # Determine the colors used for each token.
        species = getattr(args, "species", "homo_sapiens") or "homo_sapiens"
        color_to_words = _build_color_to_words(family, args.colours, species)
        try:
            grouped_color_func = SimpleGroupedColorFunc(color_to_words, "grey")
        except TypeError as exc:
            raise TCRcloudError(
                f"{args.colours} doesn't seem properly formatted. Check https://github.com/eriicdesousa/TCRcloud for more information"
            ) from exc
        wordcloud_obj.recolor(color_func=grouped_color_func)

        # Plot the wordcloud in a fixed square area (5x5) and put legend below.
        fig = plt.figure(figsize=(7, 7), dpi=300)
        # wordcloud axis as large as possible while keeping legend just underneath
        ax_wordcloud = fig.add_axes([0.05, 0.28, 0.9, 0.68])
        ax_wordcloud.imshow(wordcloud_obj, interpolation="bilinear")
        ax_wordcloud.set_xticks([])
        ax_wordcloud.set_yticks([])

        if legend and args.colours is None:
            # The legend maps V-gene calls to colors, which only makes sense
            # when colors are derived from V-gene calls (the built-in
            # species palettes). When a custom `--colours` JSON file is used,
            # colors are assigned directly to arbitrary CDR3 words instead of
            # V-genes, so there is no V-gene/color mapping to show and the
            # legend is intentionally skipped even if `--legend` is True.
            colour_map = {v: _vcall_color(v, species) for v in set(family.values())}
            ax_legend = fig.add_axes([0.05, 0.25, 0.9, 0.22])
            ax_legend.axis("off")
            with plt.rc_context({"axes.axisbelow": True}):
                plt.sca(ax_legend)
                _add_legend(colour_map)

        outputname = f"{input_stem}_{repertoire_id}_{chain}.{output_format}"
        fig.savefig(outputname, dpi=300, bbox_inches="tight", pad_inches=0.2)
        plt.close(fig)
        print("Word cloud saved as " + outputname)
