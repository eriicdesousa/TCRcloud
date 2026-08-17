"""Radar plot generation for TCRcloud.

This module computes per-repertoire diversity metrics and renders them as a
radar/spider plot. Some metrics span many orders of magnitude, so a mix of
scaling strategies is used to map each metric onto the shared [0, 1] axis:
linear (most metrics), log (Distinct CDR3, Chao1), and a narrow-range
("tail") linear scale for the Gini-Simpson index, whose axis range is
restricted to [0.7, 1.0] so that values close to 1 are easier to tell
apart. Note that "tail-log" is a linear scale over that narrow range, not
an actual logarithmic transform - see _METRIC_SCALES below.

"""

import json
import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skbio

import tcrcloud.format
from tcrcloud.errors import TCRcloudError

logger = logging.getLogger(__name__)

# For large-range metrics (Distinct CDR3, Chao1), the axis uses log scaling.
# For the Gini-Simpson metric, the axis range is restricted to its upper
# "tail" ([0.7, 1.0]) rather than the full [0, 1] range, so that values
# close to 1 (where most real repertoires fall) are easier to distinguish
# visually.
_METRIC_RANGES = np.array(
    [
        [0.0, 50.0],  # D50 Index (linear)
        [0.0, 1.0],  # 1 - Gini Coefficient (linear)
        [0.0, 15],  # Shannon Index (linear)
        [0.7, 1.0],  # Gini-Simpson Index (tail-log)
        [1.0, 250000.0],  # Distinct CDR3 (log)
        [1.0, 250000.0],  # Chao1 Index (log)
        [0.0, 1.0],  # Pielou Evenness (linear)
    ],
    dtype=float,
)

# Scale type per metric (used when mapping raw values into [0,1]).
# NOTE: "tail-log" is scaled identically to "linear" (see the `is_log` mask
# in calculate_metrics and the `scale == "log"` check in _scale_value_to_01
# below); the visual "zoom" for Gini-Simpson comes entirely from its narrow
# _METRIC_RANGES bounds ([0.7, 1.0]), not from a log transform. The distinct
# label is kept mainly for documentation purposes.
_METRIC_SCALES = [
    "linear",  # D50
    "linear",  # Gini
    "linear",  # Shannon
    "tail-log",  # Gini-Simpson (narrow linear range, emphasizes values close to 1)
    "log",  # Distinct CDR3
    "log",  # Chao1
    "linear",  # Pielou Evenness
]


def calculate_dfifty(df: pd.DataFrame, length: int) -> float:
    """Compute D50 (% clones needed to reach 50% of total counts).

    D50 is the percentage of unique clones required to accumulate 50% of the
    total counts. This is computed on the top 10,000 clones to keep runtime
    bounded for large repertoires.

    `df["counts"]` must already be sorted in descending order (as returned by
    `tcrcloud.format.format_metrics`), since the top-10,000 truncation below
    and the cumulative-sum search both assume the largest clones come first.
    """

    counts = np.asarray(df["counts"], dtype=float)
    if length >= 10000:
        counts = counts[:10000]

    total = counts.sum()
    if total <= 0:
        return 0.0

    # Find how many top clones are needed to exceed 50% cumulative count.
    idx = int(np.searchsorted(np.cumsum(counts), total / 2, side="right"))
    return (idx + 1) * 100.0 / (10000 if length >= 10000 else length)


def _load_legend_mapping(legend_file: str) -> dict[str, str]:
    """Load the repertoire_id -> legend label mapping from a JSON file.

    Raises `tcrcloud.errors.TCRcloudError` (instead of exiting the process
    directly) so that callers - and TCRcloud.py's centralized error
    handling in particular - can surface a clean "TCRcloud error: ..."
    message, consistent with tcrcloud.cloud._load_colour_mapping.
    """

    try:
        with open(legend_file) as json_file:
            return json.load(json_file)
    except FileNotFoundError as exc:
        raise TCRcloudError(f"{legend_file} doesn't seem to exist") from exc
    except json.decoder.JSONDecodeError as exc:
        raise TCRcloudError(
            f"{legend_file} doesn't seem properly formatted. "
            "Check https://github.com/eriicdesousa/TCRcloud for more information"
        ) from exc


def calculate_metrics(
    keys: Sequence[tuple[str, str]],
    samples: pd.core.groupby.DataFrameGroupBy,
    legend_file: str | None,
    export: bool,
    filename: str,
) -> tuple[list[list[Any]], np.ndarray, np.ndarray, list[str]]:
    """Compute all radar metrics for each repertoire.

    This function computes a fixed set of repertoire diversity metrics and
    then scales them into the [0, 1] range according to the fixed ranges
    defined in _METRIC_RANGES/_METRIC_SCALES.

    Returns:
        (datasets, min_vals, max_vals, scales) where `datasets` is a list of
        `[label, *scaled_metric_values]` (one entry per repertoire/chain
        key), and `min_vals`/`max_vals`/`scales` mirror the rows of
        _METRIC_RANGES/_METRIC_SCALES (used by `radar()` to draw axis ticks).
    """

    legend_dict = {}
    if legend_file is not None:
        legend_dict = _load_legend_mapping(legend_file)

    chain_names = {
        "A": "α chain",
        "B": "β chain",
        "G": "γ chain",
        "D": "δ chain",
        "H": "Heavy chain",
        "K": "Kappa chain",
        "L": "Lambda chain",
    }

    def _format_label(key: tuple[str, str]) -> str:
        prefix = legend_dict.get(key[1], key[1])
        suffix = chain_names.get(key[0], "")
        return f"{prefix} {suffix}".strip()

    raw_metrics = []
    labels = []

    for key in keys:
        df = samples.get_group(key)

        df = tcrcloud.format.format_metrics(df)
        length = len(df)
        counts = df["counts"].to_numpy(dtype=float)

        raw_metrics.append(
            [
                calculate_dfifty(df, length),
                # Gini index measures inequality (higher = less diverse), so
                # it's inverted here (1 - Gini) to keep every plotted axis
                # consistent: higher value always means more diversity.
                1 - skbio.diversity.alpha_diversity("gini_index", counts)[0],
                skbio.diversity.alpha_diversity("shannon", counts)[0],
                skbio.diversity.alpha_diversity("simpson", counts)[0],
                length,
                skbio.diversity.alpha_diversity("chao1", counts)[0],
                skbio.diversity.alpha_diversity("pielou_e", counts)[0],
            ]
        )
        labels.append(_format_label(key))

    if not raw_metrics:
        # No repertoires to plot; the min/max placeholders below are never
        # used since `radar()` bails out early when `datasets` is empty.
        return [], np.zeros(7), np.zeros(7), _METRIC_SCALES

    raw_arr = np.asarray(raw_metrics, dtype=float)

    # Scale each raw metric into [0, 1] using the fixed per-metric ranges
    # defined in _METRIC_RANGES above.
    min_vals = _METRIC_RANGES[:, 0]
    max_vals = _METRIC_RANGES[:, 1]

    # Vectorize scaling for speed.
    scaled_arr = np.empty_like(raw_arr, dtype=float)

    # Apply a log transform for "log"-scaled metrics; everything else
    # ("linear" and "tail-log") is scaled linearly within its own range.
    is_log = np.array([s == "log" for s in _METRIC_SCALES], dtype=bool)
    if np.any(is_log):
        # Clip to avoid log of zero/negative values.
        raw_log = np.clip(raw_arr[:, is_log], min_vals[is_log], None)
        raw_log = np.log10(raw_log)
        min_log = np.log10(min_vals[is_log])
        max_log = np.log10(max_vals[is_log])
        scaled_arr[:, is_log] = (raw_log - min_log) / (max_log - min_log)

    if np.any(~is_log):
        raw_lin = raw_arr[:, ~is_log]
        min_lin = min_vals[~is_log]
        max_lin = max_vals[~is_log]
        scaled_arr[:, ~is_log] = (raw_lin - min_lin) / (max_lin - min_lin)

    scaled_arr = np.clip(scaled_arr, 0.0, 1.0)

    datasets = [[label, *row.tolist()] for label, row in zip(labels, scaled_arr, strict=True)]

    if export:
        metrics_filename = (
            str(Path(filename).with_suffix("")) + "_repertoire_metrics.txt"
        )
        with open(metrics_filename, "w") as fileout:
            for label, row in zip(labels, raw_arr, strict=True):
                print("Repertoire:", label, file=fileout)
                print(
                    "D50 Index:",
                    float(np.format_float_positional(row[0], precision=3)),
                    file=fileout,
                )
                print(
                    "1 - Gini Coefficient:",
                    float(np.format_float_positional(row[1], precision=3)),
                    file=fileout,
                )
                print(
                    "Shannon Index:",
                    float(np.format_float_positional(row[2], precision=3)),
                    file=fileout,
                )
                print(
                    "Gini-Simpson Index:",
                    float(np.format_float_positional(row[3], precision=3)),
                    file=fileout,
                )
                print(
                    "Distinct CDR3:",
                    float(np.format_float_positional(row[4], precision=3)),
                    file=fileout,
                )
                print(
                    "Chao1 Index:",
                    float(np.format_float_positional(row[5], precision=3)),
                    file=fileout,
                )
                print(
                    "Pielou Evenness:",
                    float(np.format_float_positional(row[6], precision=3)),
                    file=fileout,
                )
                print(file=fileout)
        logger.info("Repertoire metrics saved as " + metrics_filename)

    return datasets, min_vals, max_vals, _METRIC_SCALES


def _scale_value_to_01(
    value: float, min_val: float, max_val: float, scale: str
) -> float:
    """Map a raw metric value onto the shared [0, 1] radar axis.

    "tail-log" metrics fall through to the linear branch below (see the
    NOTE on _METRIC_SCALES above).
    """
    if scale == "log":
        min_t = np.log10(min_val)
        max_t = np.log10(max_val)
        value = np.clip(value, min_val, max_val)
        return (np.log10(value) - min_t) / (max_t - min_t)
    v = np.clip(value, min_val, max_val)
    return (v - min_val) / (max_val - min_val)


def _format_tick(value: float, is_int: bool) -> Any:
    """Format a single axis tick value as a display label."""
    if is_int:
        return int(round(value))
    if value >= 1000:
        # Avoid compact formats like 10k/1M; show full integer with commas.
        return f"{int(round(value)):,}"
    # Avoid scientific notation for very small values.
    if 0 < abs(value) < 1e-3:
        return f"{value:.6f}".rstrip("0").rstrip(".")
    return round(value, 5)


def radar(
    rearrangements: str,
    custom_legend: str | None = None,
    legend: bool = True,
    export: bool = False,
    output_format: str = "png",
) -> None:
    """Entry point for the `TCRcloud radar` command."""

    # Determine the output image format (defaults to "png"). Validated here
    # too since `radar()` may be called directly as a library function rather
    # than only via the argparse CLI, whose `choices=["svg", "png"]` wouldn't
    # otherwise catch a bad value.
    output_format = output_format.strip().lower()
    if output_format not in ("svg", "png"):
        raise TCRcloudError(
            f"unsupported output format '{output_format}'. "
            "Please choose 'svg' or 'png'"
        )

    # The radar plot categories (one per metric) + repeat first category for closure.
    categories = [
        "D50\nIndex",
        "1 - Gini\nCoefficient",
        "Shannon\nIndex",
        "Gini-Simpson\nIndex",
        "Distinct\nCDR3",
        "Chao1\nIndex",
        "Pielou\nEvenness",
    ]
    categories = [*categories, categories[0]]

    # Load and filter the input repertoire TSV file, then group by chain and repertoire.
    samples_df = tcrcloud.format.format_data(rearrangements)
    samples = samples_df.groupby(["chain", "repertoire_id"])
    keys = list(samples.groups.keys())

    datasets, min_vals, max_vals, scales = calculate_metrics(
        keys,
        samples,
        custom_legend,
        export,
        rearrangements,
    )

    if not datasets:
        raise TCRcloudError("no repertoires found for plotting")

    label_loc = np.linspace(start=0, stop=2 * np.pi, num=len(datasets[0]))

    plt.figure(figsize=(10, 14))
    plt.subplot(polar=True)

    # Predefined palette for the radar lines; will cycle if there are more
    # repertoires than colors.
    radar_colours = [
        "#f0e442",
        "#0072b2",
        "#cc79a7",
        "#009e73",
        "#e69f00",
        "#56b4e9",
        "#8a6bbf",
    ]

    for i in datasets:
        dataset = i[1:]
        dataset = [*dataset, dataset[0]]

        try:
            thecolour = radar_colours.pop(0)
        except IndexError:
            thecolour = "#BBBBBB"

        plt.plot(
            label_loc,
            dataset,
            # No label here: plt.fill (below) draws the same data and is
            # given the legend label instead, so each repertoire only gets
            # a single legend entry rather than one per plot()+fill() pair.
            linewidth=4.0,
            alpha=0.4,
            color=thecolour,
        )

        plt.fill(
            label_loc, dataset, label=i[0], linewidth=4.0, alpha=0.6, color=thecolour
        )

    # Whether each metric's tick labels should be rendered as whole numbers;
    # positions match the metric order in categories/_METRIC_RANGES above
    # (only Distinct CDR3 and Chao1 are counts, so only they round to int).
    integer_metrics = [False, False, False, False, True, True, False]

    # Draw tick labels for each axis based on the fixed default ranges.
    for idx, (min_val, max_val, scale, is_int) in enumerate(
        zip(min_vals, max_vals, scales, integer_metrics, strict=True)
    ):
        if scale == "log":
            lo = int(np.floor(np.log10(min_val)))
            hi = int(np.ceil(np.log10(max_val)))
            tick_values = [10**e for e in range(lo, hi + 1)]
        else:
            tick_values = [
                min_val + (max_val - min_val) * f for f in (0.1, 0.3, 0.5, 0.7, 0.9)
            ]

        # Optionally suppress the first (minimum) tick for cleaner axis labeling.
        if len(tick_values) > 1 and abs(tick_values[0] - min_val) < 1e-12:
            tick_values = tick_values[1:]

        # For Distinct CDR3 and Chao1 (idx 4, 5), the log-scale tick range is
        # rounded outward to the nearest power of 10 (`hi = ceil(log10(max_val))`),
        # so the last generated tick can overshoot max_val (e.g. 1,000,000
        # instead of the actual 250,000 max). That overshot tick would be
        # clipped to the same axis position as max_label below, drawing a
        # misleading duplicate/overlapping value there, so it's dropped.
        if idx in (4, 5) and len(tick_values) > 1:
            tick_values = tick_values[:-1]

        for tick in tick_values:
            pos = _scale_value_to_01(tick, min_val, max_val, scale)
            label = _format_tick(tick, is_int)
            plt.text(
                label_loc[idx],
                pos,
                label,
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=12,
                fontweight="bold",
            )

        max_label = _format_tick(max_val, is_int)
        # Strip trailing .0 for axis max values to keep labels clean
        if isinstance(max_label, float) and max_label.is_integer():
            max_label = int(max_label)
        plt.text(
            label_loc[idx],
            1.0,
            max_label,
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
            fontweight="bold",
        )

    plt.ylim(0, 1.01)
    plt.yticks([0.1, 0.3, 0.5, 0.7, 0.9], [])
    plt.tick_params(pad=32, labelsize=16)
    lines, labels = plt.thetagrids(np.degrees(label_loc), labels=categories)
    outputname = (
        str(Path(rearrangements).with_suffix("")) + "_radar" + "." + output_format
    )
    if legend:
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), fontsize=16)
    plt.savefig(outputname, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Radar saved as " + outputname)
