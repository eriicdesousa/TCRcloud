"""Unit tests for tcrcloud.aminoacids."""

import argparse

import pandas as pd
import pytest

import tcrcloud.aminoacids as aminoacids
import tcrcloud.format as tformat
from tcrcloud.errors import TCRcloudError

# ---------------------------------------------------------------------------
# generate_mesh
# ---------------------------------------------------------------------------


def test_generate_mesh_inflates_upper_bounds_by_half():
    # x_max/y_max are inflated by 0.5 inside generate_mesh so adjacent bars
    # touch without gaps; the returned mesh's x/y coordinates must reflect
    # that inflated value (1.5), not the raw input (1).
    mesh = aminoacids.generate_mesh(
        x_min=0, x_max=1, y_min=0, y_max=1, z_min=0, z_max=5, color_value="#ff0000"
    )
    assert max(mesh.x) == pytest.approx(1.5)
    assert max(mesh.y) == pytest.approx(1.5)
    assert min(mesh.x) == 0
    assert min(mesh.y) == 0
    assert max(mesh.z) == 5
    assert mesh.color == "#ff0000"


def test_generate_mesh_returns_mesh3d_with_expected_vertex_count():
    mesh = aminoacids.generate_mesh(0, 1, 0, 1, 0, 1, "blue")
    assert len(mesh.x) == 8
    assert len(mesh.y) == 8
    assert len(mesh.z) == 8


# ---------------------------------------------------------------------------
# desired_order / colours consistency
# ---------------------------------------------------------------------------
#
# `normalized.reindex(desired_order[::-1])` and the padding loop in
# aminoacids() assume every amino acid in `colours` also appears (exactly
# once) in `desired_order`, and vice versa.


def test_desired_order_matches_colours_keys_exactly():
    assert sorted(aminoacids.desired_order) == sorted(aminoacids.colours.keys())


def test_desired_order_has_no_duplicates():
    assert len(aminoacids.desired_order) == len(set(aminoacids.desired_order)) == 20


# ---------------------------------------------------------------------------
# aminoacids()
# ---------------------------------------------------------------------------
#
# `tcrcloud.format.format_data`/`format_aminoacids` are mocked out so these
# tests exercise only aminoacids()'s plotting/export/format handling, not the
# AIRR loading pipeline.


def _fake_formatted_aminoacids_samples():
    return pd.DataFrame(
        {
            "junction_aa": ["CASSAF", "CASSGF", "CASSLF"],
            "repertoire_id": ["rep1"] * 3,
            "chain": ["B"] * 3,
            "counts": [3, 2, 1],
        }
    )


def _patch_format_pipeline(monkeypatch):
    monkeypatch.setattr(tformat, "format_data", lambda args: pd.DataFrame())
    monkeypatch.setattr(
        tformat, "format_aminoacids", lambda df: _fake_formatted_aminoacids_samples()
    )


def _make_args(tmp_path, **overrides):
    rearrangements = tmp_path / "repertoire.tsv"
    rearrangements.write_text(
        "junction_aa\tv_call\tj_call\tjunction\trepertoire_id\tproductive\n"
    )
    defaults = dict(
        rearrangements=str(rearrangements),
        export=False,
        threeD=False,
        format="png",
    )
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


def test_aminoacids_writes_2d_png_by_default(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    aminoacids.aminoacids(_make_args(tmp_path))

    outputs = list(tmp_path.glob("*.png"))
    assert [p.name for p in outputs] == ["repertoire_aminoacids_rep1_B.png"]
    assert not list(tmp_path.glob("*.svg"))
    assert not list(tmp_path.glob("*.html"))


def test_aminoacids_writes_2d_svg_when_requested(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    aminoacids.aminoacids(_make_args(tmp_path, format="svg"))

    outputs = list(tmp_path.glob("*.svg"))
    assert [p.name for p in outputs] == ["repertoire_aminoacids_rep1_B.svg"]
    assert not list(tmp_path.glob("*.png"))


def test_aminoacids_format_is_case_insensitive(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    aminoacids.aminoacids(_make_args(tmp_path, format="PNG"))

    assert [p.name for p in tmp_path.glob("*.png")] == [
        "repertoire_aminoacids_rep1_B.png"
    ]


def test_aminoacids_rejects_unsupported_format(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)

    with pytest.raises(TCRcloudError, match="unsupported output format"):
        aminoacids.aminoacids(_make_args(tmp_path, format="pdf"))


def test_aminoacids_defaults_to_png_when_format_missing(tmp_path, monkeypatch):
    # `aminoacids()` may be called with an args object that predates the
    # --format option (e.g. older scripts); it should default to png.
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)
    args = _make_args(tmp_path)
    del args.format

    aminoacids.aminoacids(args)

    assert [p.name for p in tmp_path.glob("*.png")] == [
        "repertoire_aminoacids_rep1_B.png"
    ]


def test_aminoacids_writes_3d_png_and_html_when_threed_requested(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    aminoacids.aminoacids(_make_args(tmp_path, threeD=True))

    png_outputs = list(tmp_path.glob("*.png"))
    html_outputs = list(tmp_path.glob("*.html"))
    assert [p.name for p in png_outputs] == ["repertoire_aminoacids3D_rep1_B.png"]
    assert [p.name for p in html_outputs] == ["repertoire_aminoacids3D_rep1_B.html"]


def test_aminoacids_3d_svg_still_writes_html_extension(tmp_path, monkeypatch):
    # The interactive HTML export name must not depend on the chosen static
    # image format (previously derived via outputname.replace(".png", ".html"),
    # which silently produced no ".html" file for non-png formats).
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    aminoacids.aminoacids(_make_args(tmp_path, threeD=True, format="svg"))

    svg_outputs = list(tmp_path.glob("*.svg"))
    html_outputs = list(tmp_path.glob("*.html"))
    assert [p.name for p in svg_outputs] == ["repertoire_aminoacids3D_rep1_B.svg"]
    assert [p.name for p in html_outputs] == ["repertoire_aminoacids3D_rep1_B.html"]


def test_aminoacids_accepts_string_threed_and_export_flags(tmp_path, monkeypatch):
    # aminoacids() may be called directly (e.g. from older scripts/tests)
    # with plain "true"/"false" strings rather than real bools.
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    aminoacids.aminoacids(_make_args(tmp_path, threeD="true", export="true"))

    assert list(tmp_path.glob("*_aminoacids3D_*.png"))
    assert list(tmp_path.glob("*_aminoacids_table*.csv"))


def test_aminoacids_exports_csv_table_when_requested(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    aminoacids.aminoacids(_make_args(tmp_path, export=True))

    outputs = list(tmp_path.glob("*_aminoacids_table*.csv"))
    assert [p.name for p in outputs] == ["repertoire_aminoacids_tablerep1_B.csv"]


def test_aminoacids_does_not_export_csv_by_default(tmp_path, monkeypatch):
    _patch_format_pipeline(monkeypatch)
    monkeypatch.chdir(tmp_path)

    aminoacids.aminoacids(_make_args(tmp_path))

    assert not list(tmp_path.glob("*.csv"))


def test_aminoacids_writes_one_plot_per_chain_repertoire(tmp_path, monkeypatch):
    def fake_format_aminoacids(df):
        return pd.DataFrame(
            {
                "junction_aa": ["CASSAF", "CASSGF", "CAVSAF", "CAVSGF"],
                "repertoire_id": ["rep1", "rep1", "rep2", "rep2"],
                "chain": ["B", "B", "A", "A"],
                "counts": [3, 2, 4, 1],
            }
        )

    monkeypatch.setattr(tformat, "format_data", lambda args: pd.DataFrame())
    monkeypatch.setattr(tformat, "format_aminoacids", fake_format_aminoacids)
    monkeypatch.chdir(tmp_path)

    aminoacids.aminoacids(_make_args(tmp_path))

    outputs = {p.name for p in tmp_path.glob("*.png")}
    assert outputs == {
        "repertoire_aminoacids_rep1_B.png",
        "repertoire_aminoacids_rep2_A.png",
    }
