"""Formatting utilities for TCRcloud.

This module provides helpers to load and clean AIRR rearrangement data and to
produce summary tables used by the downstream TCRcloud analysis pipeline.

Expected input is an AIRR rearrangement file (TSV) that may include a
`duplicate_count` column. The downstream pipeline assumes the output contains
cleaned `junction_aa`, `v_call`, `j_call`, and an inferred `chain` value.
"""

from __future__ import annotations

from typing import Sequence

import pandas as pd

import airr


INVALID_CDR3_CHARS = {"X", "*", "B", "Z", "J", "_"}


def _clean_v_calls(v_call: str) -> tuple[str, str]:
    """Extract two key characters from a V call string for matching."""

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
    # Determine which columns are present in the input rearrangement file.
    # Some AIRR exports include a `duplicate_count` column; we need to keep it if
    # present so that downstream aggregates can sum properly.
    with open(args.rearrangements) as f:
        first_line = f.readline()
        if "duplicate_count" in first_line:
            keys = [
                "junction_aa",
                "v_call",
                "j_call",
                "junction",
                "repertoire_id",
                "duplicate_count",
                "productive",
            ]
        else:
            keys = [
                "junction_aa",
                "v_call",
                "j_call",
                "junction",
                "repertoire_id",
                "productive",
            ]

    # Validate the file in-place with AIRR schema rules and open a streaming reader.
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

    # Build a DataFrame from the filtered records.
    df = pd.DataFrame(valid_rows)

    # If multiple V gene assignments are present (comma-separated), keep only the
    # first one.
    df["v_call"] = df["v_call"].str.split(",", n=1, expand=True)[0]

    # Drop allele suffixes (e.g. '*01' etc.) from v_call for consistent grouping.
    if df["v_call"].astype(str).str.contains(r"\*", regex=True).any():
        df["v_call"] = df["v_call"].str.replace(r"\*.*$", "", regex=True)

    # Infer chain information (TCR alpha/beta etc.) from V/J calls.
    is_dv_dj = df["v_call"].str.contains("DV", na=False) & df["j_call"].str.contains(
        "DJ", na=False
    )

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
        dup["duplicate_count"] = pd.to_numeric(dup["duplicate_count"], errors="coerce").fillna(1)
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


def format_convergence(df: pd.DataFrame) -> pd.DataFrame:
    """Compute convergence (# unique junctions per CDR3/repertoire/chain)."""

    # Count how many unique junctions exist for each (CDR3, repertoire, chain).
    per_junction = (
        df.groupby(
            ["junction_aa", "junction", "repertoire_id", "chain"]
        )  # noqa: WPS221
        .size()
        .reset_index(name="count")
    )

    # Then count how many distinct junctions each CDR3 has across repertoires.
    return (
        per_junction.groupby(["junction_aa", "repertoire_id", "chain"])  # noqa: WPS221
        .size()
        .reset_index(name="counts")
        .sort_values(by="counts", ascending=False)
    )


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
