"""Unit tests for tcrcloud.cloud."""

import json

import pandas as pd
import pytest

import tcrcloud.cloud as cloud
import tcrcloud.format as tformat
from tcrcloud.errors import TCRcloudError

# ---------------------------------------------------------------------------
# legend natural ordering (via natsort)
# ---------------------------------------------------------------------------


def test_add_legend_orders_vgene_labels_naturally(monkeypatch):
    # cloud._add_legend should natural-sort V-gene labels so that gene
    # numbers order numerically (TRBV2 before TRBV10), not lexicographically.
    captured = {}

    def fake_legend(*args, **kwargs):
        captured["handles"] = kwargs.get("handles") or args[0]

    monkeypatch.setattr(cloud.plt, "legend", fake_legend)

    cloud._add_legend({"TRBV10": "#111111", "TRBV2": "#222222", "TRBV1": "#333333"})

    labels = [h.get_label() for h in captured["handles"]]
    assert labels == ["TRBV1", "TRBV2", "TRBV10"]


# ---------------------------------------------------------------------------
# SimpleGroupedColorFunc
# ---------------------------------------------------------------------------


def test_simple_grouped_color_func_returns_mapped_color():
    color_func = cloud.SimpleGroupedColorFunc(
        {"#ff0000": ["CASSA"], "#00ff00": ["CASSB"]}, default_color="grey"
    )
    assert color_func("CASSA") == "#ff0000"
    assert color_func("CASSB") == "#00ff00"


def test_simple_grouped_color_func_returns_default_for_unknown_word():
    color_func = cloud.SimpleGroupedColorFunc(
        {"#ff0000": ["CASSA"]}, default_color="grey"
    )
    assert color_func("UNKNOWN") == "grey"


# ---------------------------------------------------------------------------
# handle_duplicates
# ---------------------------------------------------------------------------


def test_handle_duplicates_is_noop_when_junction_aa_is_unique():
    df = pd.DataFrame({"junction_aa": ["CASSA", "CASSB"]})
    result = cloud.handle_duplicates(df)
    assert result["junction_aa"].tolist() == ["CASSA", "CASSB"]


def test_handle_duplicates_keeps_first_occurrence_unprefixed():
    # The first occurrence of a repeated CDR3 must keep its original spelling
    # (0 leading spaces); later occurrences get 1, 2, ... leading spaces so
    # that WordCloud treats them as distinct tokens.
    df = pd.DataFrame({"junction_aa": ["CASSA", "CASSA", "CASSA", "CASSB"]})
    result = cloud.handle_duplicates(df.copy())
    assert result["junction_aa"].tolist() == ["CASSA", " CASSA", "  CASSA", "CASSB"]


# ---------------------------------------------------------------------------
# _ensure_required_columns
# ---------------------------------------------------------------------------


def test_ensure_required_columns_passes_when_all_present():
    df = pd.DataFrame(
        columns=["junction_aa", "v_call", "counts", "chain", "repertoire_id"]
    )
    cloud._ensure_required_columns(df)  # should not raise


def test_ensure_required_columns_raises_when_missing():
    df = pd.DataFrame(columns=["junction_aa", "v_call"])
    with pytest.raises(TCRcloudError, match="missing required columns"):
        cloud._ensure_required_columns(df)


# ---------------------------------------------------------------------------
# _extract_family_and_text
# ---------------------------------------------------------------------------


def test_extract_family_and_text_with_multiple_rows():
    df = pd.DataFrame(
        {
            "junction_aa": ["CASSA", "CASSB"],
            "v_call": ["TRBV1", "TRBV2"],
            "counts": [3, 5],
        }
    )
    family, text = cloud._extract_family_and_text(df)
    assert family == {"CASSA": "TRBV1", "CASSB": "TRBV2"}
    assert text == {"CASSA": 3, "CASSB": 5}


def test_extract_family_and_text_with_single_row():
    df = pd.DataFrame({"junction_aa": ["CASSA"], "v_call": ["TRBV1"], "counts": [7]})
    family, text = cloud._extract_family_and_text(df)
    assert family == {"CASSA": "TRBV1"}
    assert text == {"CASSA": 7}


# ---------------------------------------------------------------------------
# _load_colour_mapping
# ---------------------------------------------------------------------------


def test_load_colour_mapping_reads_valid_json(tmp_path):
    colours_file = tmp_path / "colours.json"
    colours_file.write_text(json.dumps({"#000000": ["CASSA"]}))
    assert cloud._load_colour_mapping(str(colours_file)) == {"#000000": ["CASSA"]}


def test_load_colour_mapping_missing_file_raises_with_its_own_name(tmp_path):
    missing_path = tmp_path / "missing_colours.json"
    with pytest.raises(TCRcloudError, match=r"doesn't seem to exist"):
        cloud._load_colour_mapping(str(missing_path))
    # The message must reference the colours file itself, not some other
    # filename, so that main() (see TCRcloud.py) can surface it verbatim.
    try:
        cloud._load_colour_mapping(str(missing_path))
    except TCRcloudError as exc:
        assert str(missing_path) in str(exc)


def test_load_colour_mapping_invalid_json_raises_tcrcloud_error(tmp_path):
    bad_file = tmp_path / "bad_colours.json"
    bad_file.write_text("{not valid json")
    with pytest.raises(TCRcloudError, match="doesn't seem properly formatted"):
        cloud._load_colour_mapping(str(bad_file))


# ---------------------------------------------------------------------------
# _build_color_to_words
# ---------------------------------------------------------------------------


def test_build_color_to_words_uses_colours_path_when_provided(tmp_path):
    colours_file = tmp_path / "colours.json"
    colours_file.write_text(json.dumps({"#ff0000": ["CASSA"]}))
    result = cloud._build_color_to_words({}, str(colours_file))
    assert result == {"#ff0000": ["CASSA"]}


def test_build_color_to_words_derives_colors_from_vcall_when_no_path(monkeypatch):
    monkeypatch.setattr(
        cloud,
        "_vcall_color",
        lambda vcall, species="homo_sapiens", default="grey": f"color-{vcall}",
    )
    family = {"CASSA": "TRBV1", "CASSB": "TRBV1", "CASSC": "TRBV2"}
    result = cloud._build_color_to_words(family, colours_path=None)
    assert result == {
        "color-TRBV1": ["CASSA", "CASSB"],
        "color-TRBV2": ["CASSC"],
    }


# ---------------------------------------------------------------------------
# wordcloud() output format (-f/--format svg or png)
# ---------------------------------------------------------------------------
#
# `tcrcloud.format.format_data`/`format_cloud` are mocked out so these tests
# exercise only the output-format handling, not the AIRR loading pipeline.


def _make_kwargs(tmp_path, **overrides):
    rearrangements = tmp_path / "repertoire.tsv"
    rearrangements.write_text(
        "junction_aa\tv_call\tj_call\tjunction\trepertoire_id\tproductive\n"
    )
    defaults = dict(
        rearrangements=str(rearrangements),
        colours=None,
        species="homo_sapiens",
        legend=True,
        size=1000,
        output_format="png",
    )
    defaults.update(overrides)
    return defaults


def _fake_formatted_samples():
    return pd.DataFrame(
        {
            "junction_aa": ["CASSA", "CASSB"],
            "v_call": ["TRBV1", "TRBV2"],
            "counts": [3, 5],
            "chain": ["B", "B"],
            "repertoire_id": ["rep1", "rep1"],
        }
    )


def _patch_format_pipeline(monkeypatch):
    monkeypatch.setattr(tformat, "format_data", lambda args: pd.DataFrame())
    monkeypatch.setattr(tformat, "format_cloud", lambda df: _fake_formatted_samples())


def test_wordcloud_writes_png_by_default(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    cloud.wordcloud(**_make_kwargs(tmp_path))

    outputs = list(tmp_path.glob("*.png"))
    assert [p.name for p in outputs] == ["repertoire_rep1_B.png"]
    assert not list(tmp_path.glob("*.svg"))


def test_wordcloud_writes_svg_when_requested(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    cloud.wordcloud(**_make_kwargs(tmp_path, output_format="svg"))

    outputs = list(tmp_path.glob("*.svg"))
    assert [p.name for p in outputs] == ["repertoire_rep1_B.svg"]
    assert not list(tmp_path.glob("*.png"))


def test_wordcloud_format_is_case_insensitive(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    cloud.wordcloud(**_make_kwargs(tmp_path, output_format="SVG"))

    assert [p.name for p in tmp_path.glob("*.svg")] == ["repertoire_rep1_B.svg"]


def test_wordcloud_rejects_unsupported_format(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)

    with pytest.raises(TCRcloudError, match="unsupported output format"):
        cloud.wordcloud(**_make_kwargs(tmp_path, output_format="pdf"))


def test_wordcloud_defaults_to_png_when_format_not_given(tmp_path, monkeypatch):
    # `output_format` is an optional parameter that defaults to "png".
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)
    rearrangements = tmp_path / "repertoire.tsv"
    rearrangements.write_text(
        "junction_aa\tv_call\tj_call\tjunction\trepertoire_id\tproductive\n"
    )

    cloud.wordcloud(str(rearrangements))

    assert [p.name for p in tmp_path.glob("*.png")] == ["repertoire_rep1_B.png"]
