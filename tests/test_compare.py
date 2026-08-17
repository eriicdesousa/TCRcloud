"""Unit tests for tcrcloud.compare."""

import pandas as pd
import pytest

import tcrcloud.colours
import tcrcloud.compare as compare
import tcrcloud.format as tformat
from tcrcloud.errors import TCRcloudError

# ---------------------------------------------------------------------------
# _load_combined
# ---------------------------------------------------------------------------


def _make_kwargs(rearrangements="repertoire.tsv", rearrangements2=None, **overrides):
    defaults = dict(
        rearrangements=rearrangements,
        rearrangements2=rearrangements2,
        species="homo_sapiens",
        export=False,
        output_format="png",
    )
    defaults.update(overrides)
    return defaults


def test_load_combined_single_file_keeps_repertoire_ids_unchanged(monkeypatch):
    df = pd.DataFrame({"repertoire_id": ["rep1", "rep2"], "chain": ["B", "B"]})
    monkeypatch.setattr(tformat, "format_data", lambda args: df)

    combined, prefix = compare._load_combined("repertoire.tsv")

    assert combined["repertoire_id"].tolist() == ["rep1", "rep2"]
    assert prefix == "repertoire"


def test_load_combined_two_files_prefixes_repertoire_ids(monkeypatch):
    def fake_format_data(rearrangements):
        if rearrangements == "fileA.tsv":
            return pd.DataFrame({"repertoire_id": ["rep1"], "chain": ["B"]})
        return pd.DataFrame({"repertoire_id": ["rep1"], "chain": ["B"]})

    monkeypatch.setattr(tformat, "format_data", fake_format_data)

    combined, prefix = compare._load_combined("fileA.tsv", rearrangements2="fileB.tsv")

    assert sorted(combined["repertoire_id"].tolist()) == ["fileA__rep1", "fileB__rep1"]
    assert prefix == "fileA_vs_fileB"


# ---------------------------------------------------------------------------
# _vgene_frequency_tables
# ---------------------------------------------------------------------------


def _patch_palette(monkeypatch, palette):
    monkeypatch.setattr(
        tcrcloud.colours,
        "get_vgene_palette",
        lambda vgene_type, species=None: dict(palette) if vgene_type == "TRBV" else {},
    )


def _make_vgene_df():
    rows = [
        {"chain": "B", "repertoire_id": "rep1", "v_call": "TRBV1", "CDR3_length": 12},
        {"chain": "B", "repertoire_id": "rep1", "v_call": "TRBV1", "CDR3_length": 12},
        {"chain": "B", "repertoire_id": "rep1", "v_call": "TRBV2", "CDR3_length": 13},
        {"chain": "B", "repertoire_id": "rep2", "v_call": "TRBV1", "CDR3_length": 13},
        {"chain": "B", "repertoire_id": "rep2", "v_call": "TRBV2", "CDR3_length": 13},
    ]
    return pd.DataFrame(rows)


def test_vgene_frequency_tables_align_columns_and_sum_to_100(monkeypatch):
    _patch_palette(monkeypatch, {"TRBV1": "#111", "TRBV2": "#222"})
    df = _make_vgene_df()

    tableA, tableB, v_genes, lengths = compare._vgene_frequency_tables(
        df, "B", "rep1", "rep2", "homo_sapiens"
    )

    assert v_genes == ["TRBV1", "TRBV2"]
    # Both repertoires' CDR3 lengths (12 for rep1, 13 for rep2) must appear
    # in the aligned column set for both tables.
    assert set(lengths) == {12, 13}
    assert tableA.shape == tableB.shape == (2, 2)
    assert tableA.values.sum() == pytest.approx(100.0)
    assert tableB.values.sum() == pytest.approx(100.0)


def test_vgene_frequency_tables_returns_none_for_unsupported_chain(monkeypatch):
    _patch_palette(monkeypatch, {"TRBV1": "#111"})
    df = _make_vgene_df()

    result = compare._vgene_frequency_tables(df, "Z", "rep1", "rep2", "homo_sapiens")

    assert result is None


# ---------------------------------------------------------------------------
# _aminoacid_position_table
# ---------------------------------------------------------------------------


def test_aminoacid_position_table_columns_sum_to_100():
    df = pd.DataFrame(
        {
            "junction_aa": ["CASSAF", "CASSGF", "CASSLF"],
            "repertoire_id": ["rep1"] * 3,
            "chain": ["B"] * 3,
            "counts": [1, 1, 1],
        }
    )

    table = compare._aminoacid_position_table(df)

    assert list(table.index) == compare.aminoacids.desired_order
    # Every CDR3-position column should sum to 100% across all amino acids.
    for col in table.columns:
        assert table[col].sum() == pytest.approx(100.0)


# ---------------------------------------------------------------------------
# compare() end-to-end (format_data mocked out)
# ---------------------------------------------------------------------------


def _fake_raw_rows():
    rows = []
    for rep, vcall in [("rep1", "TRBV1"), ("rep2", "TRBV2")]:
        for i in range(6):
            rows.append(
                {
                    "junction_aa": "CASSA" + ("A" * (i % 3)) + "F",
                    "v_call": vcall,
                    "repertoire_id": rep,
                    "chain": "B",
                }
            )
    return pd.DataFrame(rows)


def test_compare_writes_output_files_for_two_repertoires(tmp_path, monkeypatch):
    _patch_palette(monkeypatch, {"TRBV1": "#111111", "TRBV2": "#222222"})
    monkeypatch.setattr(tformat, "format_data", lambda args: _fake_raw_rows())
    monkeypatch.chdir(tmp_path)

    compare.compare(**_make_kwargs(str(tmp_path / "repertoire.tsv"), export=True))

    outputs = {p.name for p in tmp_path.glob("*")}
    assert any(name.endswith("_vgenes_diff_rep1_vs_rep2_B.png") for name in outputs)
    assert any(name.endswith("_vgenes_diff_rep1_vs_rep2_B.html") for name in outputs)
    assert any(
        name.endswith("_aminoacids3D_diff_rep1_vs_rep2_B.png") for name in outputs
    )
    assert any(
        name.endswith("_aminoacids3D_diff_rep1_vs_rep2_B.html") for name in outputs
    )
    assert any(
        name.endswith("_aminoacids2D_diff_rep1_vs_rep2_B.png") for name in outputs
    )
    assert any(
        name.endswith("_aminoacids2D_diff_squashedAA_rep1_vs_rep2_B.png")
        for name in outputs
    )
    # Export=True should also produce the underlying comparison tables.
    assert any(
        name.endswith("_vgenes_diff_table_rep1_vs_rep2_B.csv") for name in outputs
    )


def test_compare_raises_when_no_chain_has_multiple_repertoires(tmp_path, monkeypatch):
    single_rep = _fake_raw_rows()
    single_rep = single_rep[single_rep["repertoire_id"] == "rep1"]
    monkeypatch.setattr(tformat, "format_data", lambda args: single_rep)

    with pytest.raises(TCRcloudError, match="no chain had 2 or more repertoires"):
        compare.compare(**_make_kwargs(str(tmp_path / "repertoire.tsv")))


def test_compare_rejects_unsupported_format(tmp_path, monkeypatch):
    monkeypatch.setattr(tformat, "format_data", lambda args: _fake_raw_rows())

    with pytest.raises(TCRcloudError, match="unsupported output format"):
        compare.compare(**_make_kwargs(str(tmp_path / "repertoire.tsv"), output_format="pdf"))
