"""Unit tests for tcrcloud.format."""

from types import SimpleNamespace

import pandas as pd
import pytest

import tcrcloud.format as format

# ---------------------------------------------------------------------------
# _is_valid_cdr3
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cdr3, expected",
    [
        ("CASSLGTDTQYF", True),
        ("CASSLGTDTQYW", True),
        ("", False),
        ("ASSLGTDTQYF", False),  # doesn't start with C
        ("CASSLGTDTQYA", False),  # doesn't end with F/W
        ("CASSLGXDTQYF", False),  # contains invalid character X
        ("CASSLG*DTQYF", False),  # contains invalid character *
    ],
)
def test_is_valid_cdr3(cdr3, expected):
    assert format._is_valid_cdr3(cdr3) is expected


# ---------------------------------------------------------------------------
# _clean_v_calls
# ---------------------------------------------------------------------------


def test_clean_v_calls_extracts_chain_letter():
    # For a typical 8+ character V call, both the primary (index 2) and
    # secondary (index -6) positions land on the locus letter.
    assert format._clean_v_calls("TRBV20-1") == ("B", "B")


def test_clean_v_calls_extracts_secondary_locus_for_dv_genes():
    # "TRAV1-2/DV8*01" is an ambiguous alpha/delta V-gene; the secondary
    # character should land on the "D" of the "/DV" suffix.
    v_call, v_call2 = format._clean_v_calls("TRAV1-2/DV8*01")
    assert v_call == "A"
    assert v_call2 == "D"


def test_clean_v_calls_handles_empty_string():
    assert format._clean_v_calls("") == ("", "")


def test_clean_v_calls_handles_short_string():
    assert format._clean_v_calls("TR") == ("", "")


# ---------------------------------------------------------------------------
# format_data
# ---------------------------------------------------------------------------
#
# `airr.validate_rearrangement`/`airr.read_rearrangement` are mocked out so
# these tests exercise only tcrcloud's own filtering/cleaning logic, not the
# `airr` package's schema validation.


def _write_header(tmp_path, columns):
    path = tmp_path / "repertoire.airr.rearrangements.tsv"
    path.write_text("\t".join(columns) + "\n")
    return path


def _patch_airr(monkeypatch, rows):
    monkeypatch.setattr(format.airr, "validate_rearrangement", lambda *a, **k: None)
    monkeypatch.setattr(format.airr, "read_rearrangement", lambda *a, **k: iter(rows))


def test_format_data_keeps_only_productive_valid_matching_rows(tmp_path, monkeypatch):
    path = _write_header(
        tmp_path,
        ["junction_aa", "v_call", "j_call", "junction", "repertoire_id", "productive"],
    )
    rows = [
        # kept: productive, valid CDR3, matching V/J calls (beta chain)
        {
            "junction_aa": "CASSLGTDTQYF",
            "v_call": "TRBV20-1*01",
            "j_call": "TRBJ2-3*01",
            "junction": "TGTGCCAGCAGCTTAGGGACAGATACGCAGTATTTT",
            "repertoire_id": "rep1",
            "productive": "T",
        },
        # dropped: not productive
        {
            "junction_aa": "CASSLGTDTQYF",
            "v_call": "TRBV20-1*01",
            "j_call": "TRBJ2-3*01",
            "junction": "TGTGCCAGCAGCTTAGGGACAGATACGCAGTATTTT",
            "repertoire_id": "rep1",
            "productive": "F",
        },
        # dropped: invalid CDR3 (doesn't start with C)
        {
            "junction_aa": "XASSLGTDTQYF",
            "v_call": "TRBV20-1*01",
            "j_call": "TRBJ2-3*01",
            "junction": "TGTGCCAGCAGCTTAGGGACAGATACGCAGTATTTT",
            "repertoire_id": "rep1",
            "productive": "T",
        },
        # dropped: V/J calls disagree (alpha V-gene, beta J-gene)
        {
            "junction_aa": "CASSLGTDTQYF",
            "v_call": "TRAV1-2*01",
            "j_call": "TRBJ2-3*01",
            "junction": "TGTGCCAGCAGCTTAGGGACAGATACGCAGTATTTT",
            "repertoire_id": "rep1",
            "productive": "T",
        },
    ]
    _patch_airr(monkeypatch, rows)

    df = format.format_data(SimpleNamespace(rearrangements=str(path)))

    assert len(df) == 1
    assert df.iloc[0]["junction_aa"] == "CASSLGTDTQYF"
    assert df.iloc[0]["v_call"] == "TRBV20-1"  # allele suffix stripped
    assert df.iloc[0]["chain"] == "B"


def test_format_data_assigns_delta_chain_for_plain_trdv_trdj_pair(
    tmp_path, monkeypatch
):
    # Regression test: a genuine (non-ambiguous) delta V-gene like "TRDV2-1"
    # contains "DV" as a substring (from "TRDV"), which used to be
    # misdetected as the ambiguous "TRAV.../DV8" suffix pattern and made the
    # chain column read a garbage character (e.g. "2") instead of "D".
    path = _write_header(
        tmp_path,
        ["junction_aa", "v_call", "j_call", "junction", "repertoire_id", "productive"],
    )
    rows = [
        {
            "junction_aa": "CASSLGTDTQYF",
            "v_call": "TRDV2-1*01",
            "j_call": "TRDJ1*01",
            "junction": "TGTGCCAGCAGCTTAGGGACAGATACGCAGTATTTT",
            "repertoire_id": "rep1",
            "productive": "T",
        },
    ]
    _patch_airr(monkeypatch, rows)

    df = format.format_data(SimpleNamespace(rearrangements=str(path)))

    assert len(df) == 1
    assert df.iloc[0]["chain"] == "D"


def test_format_data_assigns_delta_chain_for_ambiguous_dv_suffix(tmp_path, monkeypatch):
    # The genuinely ambiguous alpha/delta notation ("TRAV.../DV8") paired
    # with a TRDJ gene must still resolve to the delta chain.
    path = _write_header(
        tmp_path,
        ["junction_aa", "v_call", "j_call", "junction", "repertoire_id", "productive"],
    )
    rows = [
        {
            "junction_aa": "CASSLGTDTQYF",
            "v_call": "TRAV1-2/DV8*01",
            "j_call": "TRDJ2*01",
            "junction": "TGTGCCAGCAGCTTAGGGACAGATACGCAGTATTTT",
            "repertoire_id": "rep1",
            "productive": "T",
        },
    ]
    _patch_airr(monkeypatch, rows)

    df = format.format_data(SimpleNamespace(rearrangements=str(path)))

    assert len(df) == 1
    assert df.iloc[0]["chain"] == "D"


def test_format_data_raises_clear_error_when_no_rows_pass_filters(
    tmp_path, monkeypatch
):
    path = _write_header(
        tmp_path,
        ["junction_aa", "v_call", "j_call", "junction", "repertoire_id", "productive"],
    )
    _patch_airr(monkeypatch, [])  # no rows at all

    with pytest.raises(ValueError, match="no productive rearrangements"):
        format.format_data(SimpleNamespace(rearrangements=str(path)))


def test_format_data_detects_duplicate_count_column_by_exact_header_match(
    tmp_path, monkeypatch
):
    path = _write_header(
        tmp_path,
        [
            "junction_aa",
            "v_call",
            "j_call",
            "junction",
            "repertoire_id",
            "duplicate_count",
            "productive",
        ],
    )
    row = {
        "junction_aa": "CASSLGTDTQYF",
        "v_call": "TRBV20-1*01",
        "j_call": "TRBJ2-3*01",
        "junction": "TGTGCCAGCAGCTTAGGGACAGATACGCAGTATTTT",
        "repertoire_id": "rep1",
        "duplicate_count": "5",
        "productive": "T",
    }
    _patch_airr(monkeypatch, [row])

    df = format.format_data(SimpleNamespace(rearrangements=str(path)))

    assert "duplicate_count" in df.columns
    assert df.iloc[0]["duplicate_count"] == "5"


def test_format_data_does_not_falsely_detect_duplicate_count_substring(
    tmp_path, monkeypatch
):
    # Regression test: detection used to be a substring check (`"duplicate_count"
    # in first_line`), which could spuriously match unrelated columns whose
    # name merely contains that substring.
    path = _write_header(
        tmp_path,
        [
            "junction_aa",
            "v_call",
            "j_call",
            "junction",
            "repertoire_id",
            "notes_duplicate_count_ref",
            "productive",
        ],
    )
    row = {
        "junction_aa": "CASSLGTDTQYF",
        "v_call": "TRBV20-1*01",
        "j_call": "TRBJ2-3*01",
        "junction": "TGTGCCAGCAGCTTAGGGACAGATACGCAGTATTTT",
        "repertoire_id": "rep1",
        "notes_duplicate_count_ref": "foo",
        "productive": "T",
    }
    _patch_airr(monkeypatch, [row])

    df = format.format_data(SimpleNamespace(rearrangements=str(path)))

    assert "duplicate_count" not in df.columns


def test_format_data_defaults_repertoire_id_when_column_missing(tmp_path, monkeypatch):
    # Regression test: when the input file has no `repertoire_id` column at
    # all, `row[k]` used to raise a `KeyError` for every row. Instead, the
    # whole file should be treated as a single repertoire named after the
    # input file.
    path = _write_header(
        tmp_path,
        ["junction_aa", "v_call", "j_call", "junction", "productive"],
    )
    row = {
        "junction_aa": "CASSLGTDTQYF",
        "v_call": "TRBV20-1*01",
        "j_call": "TRBJ2-3*01",
        "junction": "TGTGCCAGCAGCTTAGGGACAGATACGCAGTATTTT",
        "productive": "T",
    }
    _patch_airr(monkeypatch, [row])

    df = format.format_data(SimpleNamespace(rearrangements=str(path)))

    assert "repertoire_id" in df.columns
    assert df.iloc[0]["repertoire_id"] == path.stem


# ---------------------------------------------------------------------------
# _aggregate_counts
# ---------------------------------------------------------------------------


def test_aggregate_counts_sums_duplicate_count_when_present():
    df = pd.DataFrame(
        {
            "junction_aa": ["A", "A", "B"],
            "repertoire_id": ["r1", "r1", "r1"],
            "chain": ["B", "B", "B"],
            "duplicate_count": [2, 3, 1],
        }
    )
    result = format._aggregate_counts(df, ["junction_aa", "repertoire_id", "chain"])
    row_a = result[result["junction_aa"] == "A"].iloc[0]
    assert row_a["counts"] == 5


def test_aggregate_counts_counts_rows_when_no_duplicate_count():
    df = pd.DataFrame(
        {
            "junction_aa": ["A", "A", "B"],
            "repertoire_id": ["r1", "r1", "r1"],
            "chain": ["B", "B", "B"],
        }
    )
    result = format._aggregate_counts(df, ["junction_aa", "repertoire_id", "chain"])
    row_a = result[result["junction_aa"] == "A"].iloc[0]
    assert row_a["counts"] == 2


def test_aggregate_counts_sorted_descending():
    df = pd.DataFrame(
        {
            "junction_aa": ["A", "B", "B"],
            "repertoire_id": ["r1"] * 3,
            "chain": ["B"] * 3,
        }
    )
    result = format._aggregate_counts(df, ["junction_aa", "repertoire_id", "chain"])
    assert result["counts"].is_monotonic_decreasing


# ---------------------------------------------------------------------------
# format_metrics / format_cloud / format_vgene / format_aminoacids
# ---------------------------------------------------------------------------


def test_format_vgene_adds_cdr3_length():
    df = pd.DataFrame(
        {
            "junction_aa": ["CASSLGTDTQYF"],
            "v_call": ["TRBV20-1"],
            "repertoire_id": ["r1"],
            "chain": ["B"],
        }
    )
    result = format.format_vgene(df)
    # CDR3_length excludes the conserved flanking C and F/W residues.
    assert result.iloc[0]["CDR3_length"] == len("CASSLGTDTQYF") - 2
