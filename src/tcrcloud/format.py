"""Formatting utilities for TCRcloud.

This module provides helpers to load and clean AIRR rearrangement data and to
produce summary tables used by the downstream TCRcloud analysis pipeline.

Expected input is an AIRR rearrangement file (TSV) that may include a
`duplicate_count` column. The downstream pipeline assumes the output contains
cleaned `junction_aa`, `v_call`, `j_call`, and an inferred `chain` value.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import pandas as pd

import airr

INVALID_CDR3_CHARS = {"X", "*", "B", "Z", "J", "_"}


def _clean_v_calls(v_call: str) -> tuple[str, str]:
    """Extract the locus letter(s) from a V call string for V/J matching.

    Returns a `(primary, secondary)` pair of single characters used to check
    that a V call agrees with a J call:

    - `primary` is `v_call[2]`, the locus letter in a standard call, e.g.
      "TRBV20-1" -> "B", matching `j_call[2]` of "TRBJ2-3" -> "B".
    - `secondary` is `v_call[-6]`, which lands on the "D" of the "/DV"
      suffix for ambiguous alpha/delta calls like "TRAV1-2/DV8*01", so they
      can also match a "TRDJ..." j_call.

    NOTE: both offsets assume a fixed-width allele suffix ("*01") and a
    single-digit DV gene number; unusual formatting (e.g. two-digit DV gene
    numbers, or no allele suffix) can shift these positions.
    """

    if not v_call:
        return "", ""

    return (
        v_call[2] if len(v_call) > 2 else "",
        v_call[-6] if len(v_call) > 5 else "",
    )


def _is_valid_cdr3(cdr3: str) -> bool:
    """Return True if CDR3 sequence meets expected quality filters."""

    if not cdr3:
        return False
    if not cdr3.startswith("C"):
        return False
    if cdr3[-1] not in {"F", "W"}:
        return False
    if INVALID_CDR3_CHARS.intersection(cdr3):
        return False
    return True


def format_data(args):
    """Load, validate and filter an AIRR rearrangements TSV file.

    Rows are kept only if they are marked `productive`, have a CDR3
    (`junction_aa`) that passes `_is_valid_cdr3`, and have V/J calls that
    agree with each other (see `_clean_v_calls`).

    Returns a DataFrame with cleaned `junction_aa`/`v_call` columns and an
    inferred `chain` column, ready for the `format_*` aggregation helpers.
    """

    # Determine which columns are present in the input rearrangement file.
    # Some AIRR exports include a `duplicate_count` column; we need to keep it if
    # present so that downstream aggregates can sum properly.
    # `repertoire_id` is optional too: some single-repertoire exports omit it
    # entirely, in which case `row[k]` would raise a `KeyError` for every row.
    # When that happens we simply don't request the column from `row` here and
    # fill it in afterwards with a default derived from the input filename.
    with open(args.rearrangements) as f:
        first_line = f.readline()
        header_columns = first_line.rstrip("\n").split("\t")
        has_repertoire_id = "repertoire_id" in header_columns

        keys = ["junction_aa", "v_call", "j_call", "junction"]
        if has_repertoire_id:
            keys.append("repertoire_id")
        if "duplicate_count" in header_columns:
            keys.append("duplicate_count")
        keys.append("productive")

    # Validate the file against the AIRR rearrangement schema; the second
    # positional argument (`debug=True`) makes `airr` raise on the first
    # validation error instead of only logging warnings.
    airr.validate_rearrangement(args.rearrangements, True)
    reader = airr.read_rearrangement(args.rearrangements)

    # Collect filtered rows in a list so we can build a DataFrame at the end.
    valid_rows: list[dict] = []
    # Iterate over each rearrangement record and keep only the rows that
    # satisfy a set of biological quality filters.
    for row in reader:
        # Keep only productive rearrangements (per AIRR definition).
        productive = row.get("productive")
        if productive is not True and str(productive).lower() not in (
            "true",
            "t",
            "1",
        ):
            continue

        # Filter invalid or low-quality CDR3 sequences.
        cdr3 = row.get("junction_aa") or ""
        if not _is_valid_cdr3(cdr3):
            continue

        # Extract the key matching characters from the V and J calls.
        v_call_raw = row.get("v_call") or ""
        j_call_raw = row.get("j_call") or ""

        v_call, v_call2 = _clean_v_calls(v_call_raw)
        j_call = j_call_raw[2] if len(j_call_raw) > 2 else ""
        if not (v_call and j_call):
            continue

        # Keep only rows where the V and J calls agree (either exact or via
        # the secondary position in the V call string).
        if v_call != j_call and v_call2 != j_call:
            continue

        valid_rows.append({k: row[k] for k in keys})

    if not valid_rows:
        # `pd.DataFrame([])` has no columns, which would turn the next few
        # lines into a confusing `KeyError: 'v_call'` instead of explaining
        # that the file simply had no rows passing the quality filters
        # (productive, valid CDR3, matching V/J calls).
        raise ValueError(
            "TCRcloud error: no productive rearrangements with a valid CDR3 "
            f"and matching V/J calls were found in {args.rearrangements}"
        )

    # Build a DataFrame from the filtered records.
    df = pd.DataFrame(valid_rows)

    # When the input file has no `repertoire_id` column at all, treat the
    # whole file as a single repertoire, named after the input file itself.
    if not has_repertoire_id:
        df["repertoire_id"] = Path(args.rearrangements).stem

    # If multiple V gene assignments are present (comma-separated), keep only the
    # first one.
    df["v_call"] = df["v_call"].str.split(",", n=1, expand=True)[0]

    # Drop allele suffixes (e.g. '*01' etc.) from v_call for consistent grouping.
    if df["v_call"].astype(str).str.contains(r"\*", regex=True).any():
        df["v_call"] = df["v_call"].str.replace(r"\*.*$", "", regex=True)

    # Infer chain information (TCR alpha/beta/gamma/delta etc.) from V/J calls.
    # Genes like "TRAV1-2/DV8" are ambiguously alpha or delta; when paired
    # with a "TRDJ..." j_call (checked above via `_clean_v_calls`), treat
    # them as delta by reading the "D" from the "/DV" suffix (3rd-from-last
    # character, e.g. "TRAV1-2/DV8" -> "D"). Otherwise the chain is simply
    # the locus letter at index 2 of the (allele-stripped) v_call.
    #
    # NOTE: this must check for the literal "/DV" suffix (with slash), not
    # just "DV" - ordinary delta V genes like "TRDV2-1" also contain "DV" as
    # a substring (the "D" of "TRDV" followed by "V"), so a plain "DV" check
    # would misfire on genuine TRDV genes too, reading a garbage character
    # (e.g. "2" from "TRDV2-1") from `v_call[-3]` instead of "D".
    is_dv_dj = df["v_call"].str.contains("/DV", na=False, regex=False) & df[
        "j_call"
    ].str.contains("DJ", na=False)

    df["chain"] = pd.NA
    df.loc[is_dv_dj, "chain"] = df.loc[is_dv_dj, "v_call"].str[-3]
    df.loc[~is_dv_dj, "chain"] = df.loc[~is_dv_dj, "v_call"].str[2]

    return df


def _aggregate_counts(df: pd.DataFrame, group_by: Sequence[str]) -> pd.DataFrame:
    """Aggregate counts by group, using duplicate_count when available."""

    # If the input carries explicit counts per row, sum them; otherwise count
    # each row as a single unit.
    if "duplicate_count" in df.columns:
        dup = df.loc[:, list(group_by) + ["duplicate_count"]].copy()
        dup["duplicate_count"] = pd.to_numeric(
            dup["duplicate_count"], errors="coerce"
        ).fillna(1)
        agg = (
            dup.rename(columns={"duplicate_count": "counts"})
            .groupby(list(group_by), dropna=False)
            .sum()
            .reset_index()
        )
    else:
        agg = (
            df.pivot_table(index=list(group_by), aggfunc="size")
            .reset_index()
            .rename(columns={0: "counts"})
        )

    return agg.sort_values(by="counts", ascending=False)


def format_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate sequence counts for metrics."""

    return _aggregate_counts(df, ["junction_aa", "repertoire_id", "chain"])


def format_cloud(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate sequence counts for cloud visualization."""

    return _aggregate_counts(df, ["junction_aa", "v_call", "repertoire_id", "chain"])


def format_vgene(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate v-gene counts and annotate CDR3 length."""

    agg = _aggregate_counts(df, ["junction_aa", "v_call", "repertoire_id", "chain"])
    agg["CDR3_length"] = agg["junction_aa"].str[1:-1].str.len()
    return agg


def format_aminoacids(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate amino-acid composition counts."""

    return _aggregate_counts(df, ["junction_aa", "repertoire_id", "chain"])


# print("""
#                   _____ ____ ____      _                 _
#                  |_   _/ ___|  _ \ ___| | ___  _   _  __| |
#                    | || |   | |_) / __| |/ _ \| | | |/ _` |
#                    | || |___|  _ < (__| | (_) | |_| | (_| |
#                    |_| \____|_| \_\___|_|\___/ \__,_|\__,_|
#    """)
