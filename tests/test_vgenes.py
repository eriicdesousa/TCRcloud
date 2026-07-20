"""Unit tests for tcrcloud.vgenes."""

import argparse

import numpy as np
import pandas as pd
import pytest

import tcrcloud.colours
import tcrcloud.format as tformat
import tcrcloud.vgenes as vgenes

# ---------------------------------------------------------------------------
# _palette_for_chain
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "chain_letter, expected_vgene_type",
    [
        ("A", "TRAV"),
        ("B", "TRBV"),
        ("G", "TRGV"),
        ("D", "TRDV"),
        ("H", "IGHV"),
        ("K", "IGKV"),
        ("L", "IGLV"),
    ],
)
def test_palette_for_chain_maps_to_correct_vgene_type(
    monkeypatch, chain_letter, expected_vgene_type
):
    calls = []
    monkeypatch.setattr(
        tcrcloud.colours,
        "get_vgene_palette",
        lambda vgene_type, species: calls.append((vgene_type, species)) or {"X": "c"},
    )

    result = vgenes._palette_for_chain(chain_letter, species="mus_musculus")

    assert calls == [(expected_vgene_type, "mus_musculus")]
    assert result == {"X": "c"}


def test_palette_for_chain_unknown_letter_returns_empty_dict():
    assert vgenes._palette_for_chain("Z") == {}


# ---------------------------------------------------------------------------
# get_table
# ---------------------------------------------------------------------------
#
# `samples` in these tests mimics the output of
# `tcrcloud.format.format_vgene(...).groupby(["chain", "repertoire_id"])`, i.e.
# one row per distinct CDR3 sequence with its assigned v_call/CDR3_length.


def _patch_palette(monkeypatch, palettes):
    def fake_get_vgene_palette(vgene_type, species=None):
        return dict(palettes.get(vgene_type, {}))

    monkeypatch.setattr(tcrcloud.colours, "get_vgene_palette", fake_get_vgene_palette)


def _make_vgene_df(chain="B", repertoire_id="rep1", rows=None):
    if rows is None:
        rows = [
            {"junction_aa": "CASSAF", "v_call": "TRBV1", "CDR3_length": 12},
            {"junction_aa": "CASSBF", "v_call": "TRBV1", "CDR3_length": 12},
            {"junction_aa": "CASSCF", "v_call": "TRBV2", "CDR3_length": 13},
            {"junction_aa": "CASSDF", "v_call": "TRBV2", "CDR3_length": 13},
        ]
    df = pd.DataFrame(rows)
    df["repertoire_id"] = repertoire_id
    df["chain"] = chain
    df["counts"] = 1
    return df


def _make_vgene_samples(chain="B", repertoire_id="rep1", rows=None):
    df = _make_vgene_df(chain=chain, repertoire_id=repertoire_id, rows=rows)
    samples = df.groupby(["chain", "repertoire_id"])
    return samples, list(samples.groups.keys())


def test_get_table_builds_one_dataset_per_key(monkeypatch):
    _patch_palette(monkeypatch, {"TRBV": {"TRBV1": "#111", "TRBV2": "#222"}})
    samples, keys = _make_vgene_samples()
    args = argparse.Namespace(
        rearrangements="repertoire.tsv", export=False, species="homo_sapiens"
    )

    datasets = vgenes.get_table(keys, samples, args)

    assert len(datasets) == 1
    dataset = datasets[0]
    assert dataset["chain"] == "B"
    assert dataset["repertoire_id"] == "rep1"
    assert dataset["x_axis_names"] == ["TRBV1", "TRBV2"]
    assert dataset["x_axis_ticks"] == [0, 1]
    assert dataset["plot_aspect"] == vgenes._CHAIN_PLOT_SETTINGS["B"]["plot_aspect"]
    assert dataset["x_size"] == vgenes._CHAIN_PLOT_SETTINGS["B"]["x_size"]

    # `frequency` is a percentage of total reads, so it must sum to ~100
    # across the whole V-gene x CDR3-length grid.
    total_frequency = sum(np.sum(row) for row in dataset["z"])
    assert total_frequency == pytest.approx(100.0)


def test_get_table_skips_unsupported_chain_letters():
    args = argparse.Namespace(
        rearrangements="repertoire.tsv", export=False, species="homo_sapiens"
    )

    datasets = vgenes.get_table([("Z", "rep1")], None, args)

    assert datasets == []


def test_get_table_returns_multiple_datasets_for_multiple_keys(monkeypatch):
    _patch_palette(
        monkeypatch,
        {"TRBV": {"TRBV1": "#111", "TRBV2": "#222"}},
    )
    df1 = _make_vgene_df(repertoire_id="rep1")
    df2 = _make_vgene_df(repertoire_id="rep2")
    df = pd.concat([df1, df2], ignore_index=True)
    samples = df.groupby(["chain", "repertoire_id"])
    keys = list(samples.groups.keys())
    args = argparse.Namespace(
        rearrangements="repertoire.tsv", export=False, species="homo_sapiens"
    )

    datasets = vgenes.get_table(keys, samples, args)

    assert {d["repertoire_id"] for d in datasets} == {"rep1", "rep2"}
    assert all(d["chain"] == "B" for d in datasets)


def test_get_table_exports_csv_when_requested(tmp_path, monkeypatch):
    _patch_palette(monkeypatch, {"TRBV": {"TRBV1": "#111", "TRBV2": "#222"}})
    samples, keys = _make_vgene_samples()
    rearrangements = tmp_path / "repertoire.tsv"
    args = argparse.Namespace(
        rearrangements=str(rearrangements), export=True, species="homo_sapiens"
    )

    vgenes.get_table(keys, samples, args)

    expected = tmp_path / "repertoire_vgenes_tablerep1_B.csv"
    assert expected.exists()
    assert "v_call" in expected.read_text()


def test_get_table_does_not_export_csv_by_default(tmp_path, monkeypatch):
    _patch_palette(monkeypatch, {"TRBV": {"TRBV1": "#111", "TRBV2": "#222"}})
    samples, keys = _make_vgene_samples()
    rearrangements = tmp_path / "repertoire.tsv"
    args = argparse.Namespace(
        rearrangements=str(rearrangements), export=False, species="homo_sapiens"
    )

    vgenes.get_table(keys, samples, args)

    assert not list(tmp_path.glob("*.csv"))


def test_get_table_renames_ambiguous_alpha_delta_vgene_for_delta_chain(monkeypatch):
    _patch_palette(monkeypatch, {"TRDV": {"TRDV1": "#111", "TRDV4": "#222"}})
    rows = [
        {"junction_aa": "CASSAF", "v_call": "TRAV14/DV4", "CDR3_length": 10},
        {"junction_aa": "CASSBF", "v_call": "TRDV1", "CDR3_length": 11},
    ]
    samples, keys = _make_vgene_samples(chain="D", rows=rows)
    args = argparse.Namespace(
        rearrangements="repertoire.tsv", export=False, species="homo_sapiens"
    )

    datasets = vgenes.get_table(keys, samples, args)

    assert len(datasets) == 1
    dataset = datasets[0]
    # "TRAV14/DV4" must have been folded into the "TRDV4" bucket rather than
    # appearing as its own (unrecognized) V gene.
    assert dataset["x_axis_names"] == ["TRDV1", "TRDV4"]
    total_frequency = sum(np.sum(row) for row in dataset["z"])
    assert total_frequency == pytest.approx(100.0)


def test_get_table_pops_ambiguous_vgenes_from_delta_palette(monkeypatch):
    _patch_palette(
        monkeypatch,
        {
            "TRDV": {
                "TRDV1": "#111",
                # A raw palette should never expose these ambiguous
                # alpha/delta calls directly, but the pop() is defensive.
                "TRAV14/DV4": "#999",
            }
        },
    )
    rows = [{"junction_aa": "CASSAF", "v_call": "TRDV1", "CDR3_length": 10}]
    samples, keys = _make_vgene_samples(chain="D", rows=rows)
    args = argparse.Namespace(
        rearrangements="repertoire.tsv", export=False, species="homo_sapiens"
    )

    datasets = vgenes.get_table(keys, samples, args)

    assert datasets[0]["x_axis_names"] == ["TRDV1"]


# ---------------------------------------------------------------------------
# barplot() output format (-f/--format svg or png)
# ---------------------------------------------------------------------------
#
# `tcrcloud.format.format_data`/`format_vgene` are mocked out so these tests
# exercise only barplot()'s output-format/export handling, not the AIRR
# loading pipeline.


def _fake_formatted_vgene_samples():
    return pd.DataFrame(
        {
            "junction_aa": ["CASSAF", "CASSBF", "CASSCF", "CASSDF"],
            "v_call": ["TRBV1", "TRBV1", "TRBV2", "TRBV2"],
            "repertoire_id": ["rep1"] * 4,
            "chain": ["B"] * 4,
            "counts": [1, 1, 1, 1],
            "CDR3_length": [12, 12, 13, 13],
        }
    )


def _patch_format_pipeline(monkeypatch):
    monkeypatch.setattr(tformat, "format_data", lambda args: pd.DataFrame())
    monkeypatch.setattr(
        tformat, "format_vgene", lambda df: _fake_formatted_vgene_samples()
    )


def _patch_palette_for_barplot(monkeypatch):
    monkeypatch.setattr(
        tcrcloud.colours,
        "get_vgene_palette",
        lambda vgene_type, species=None: (
            {"TRBV1": "#111111", "TRBV2": "#222222"} if vgene_type == "TRBV" else {}
        ),
    )


def _make_args(tmp_path, **overrides):
    rearrangements = tmp_path / "repertoire.tsv"
    rearrangements.write_text(
        "junction_aa\tv_call\tj_call\tjunction\trepertoire_id\tproductive\n"
    )
    defaults = dict(
        rearrangements=str(rearrangements),
        export=False,
        species="homo_sapiens",
        format="png",
    )
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


def test_barplot_writes_png_by_default(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    _patch_palette_for_barplot(monkeypatch)
    monkeypatch.chdir(tmp_path)

    vgenes.barplot(_make_args(tmp_path))

    outputs = list(tmp_path.glob("*.png"))
    assert [p.name for p in outputs] == ["repertoire_vgenes_rep1_B.png"]
    assert not list(tmp_path.glob("*.svg"))
    assert [p.name for p in tmp_path.glob("*.html")] == [
        "repertoire_vgenes_rep1_B.html"
    ]


def test_barplot_writes_svg_when_requested(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    _patch_palette_for_barplot(monkeypatch)
    monkeypatch.chdir(tmp_path)

    vgenes.barplot(_make_args(tmp_path, format="svg"))

    outputs = list(tmp_path.glob("*.svg"))
    assert [p.name for p in outputs] == ["repertoire_vgenes_rep1_B.svg"]
    assert not list(tmp_path.glob("*.png"))


def test_barplot_html_export_is_independent_of_static_image_format(
    tmp_path, monkeypatch
):
    # The interactive HTML export name/content must not depend on the
    # chosen static image format.
    _patch_format_pipeline(monkeypatch)
    _patch_palette_for_barplot(monkeypatch)
    monkeypatch.chdir(tmp_path)

    vgenes.barplot(_make_args(tmp_path, format="svg"))

    assert [p.name for p in tmp_path.glob("*.html")] == [
        "repertoire_vgenes_rep1_B.html"
    ]


def test_barplot_format_is_case_insensitive(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    _patch_palette_for_barplot(monkeypatch)
    monkeypatch.chdir(tmp_path)

    vgenes.barplot(_make_args(tmp_path, format="PNG"))

    assert [p.name for p in tmp_path.glob("*.png")] == ["repertoire_vgenes_rep1_B.png"]


def test_barplot_rejects_unsupported_format(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    _patch_palette_for_barplot(monkeypatch)

    with pytest.raises(ValueError, match="unsupported output format"):
        vgenes.barplot(_make_args(tmp_path, format="pdf"))


def test_barplot_defaults_to_png_when_format_missing(tmp_path, monkeypatch):
    # `barplot()` may be called with an args object that predates the
    # --format option (e.g. older scripts); it should default to png.
    _patch_format_pipeline(monkeypatch)
    _patch_palette_for_barplot(monkeypatch)
    monkeypatch.chdir(tmp_path)
    args = _make_args(tmp_path)
    del args.format

    vgenes.barplot(args)

    assert [p.name for p in tmp_path.glob("*.png")] == ["repertoire_vgenes_rep1_B.png"]


def test_barplot_accepts_string_export_flag(tmp_path, monkeypatch):
    # `barplot()` may be called directly (e.g. from older scripts/tests)
    # with a plain "true"/"false" string rather than a real bool.
    _patch_format_pipeline(monkeypatch)
    _patch_palette_for_barplot(monkeypatch)
    monkeypatch.chdir(tmp_path)

    vgenes.barplot(_make_args(tmp_path, export="true"))

    assert list(tmp_path.glob("*_vgenes_table*.csv"))


def test_barplot_raises_when_no_repertoires_found(tmp_path, monkeypatch):
    monkeypatch.setattr(tformat, "format_data", lambda args: pd.DataFrame())
    monkeypatch.setattr(
        tformat,
        "format_vgene",
        lambda df: pd.DataFrame(
            columns=[
                "junction_aa",
                "v_call",
                "repertoire_id",
                "chain",
                "counts",
                "CDR3_length",
            ]
        ),
    )

    with pytest.raises(ValueError, match="no repertoires found"):
        vgenes.barplot(_make_args(tmp_path))
