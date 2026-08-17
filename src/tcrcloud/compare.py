"""Repertoire comparison plots for TCRcloud.

This module compares pairs of repertoires that share the same chain/locus,
producing difference plots analogous to the single-repertoire plots from:

- `tcrcloud.vgenes`: V gene x CDR3 length surface plot (here, the
  percentage-point difference between the two repertoires).
- `tcrcloud.aminoacids`: amino acid composition by CDR3 position (again, the
  percentage-point difference), shown as a 3D surface as well as a 2D
  diverging stacked bar chart (per CDR3 position) and a squashed
  per-amino-acid summary bar chart.

If a single rearrangements file contains more than one repertoire for a
given chain, every pair of repertoires for that chain is compared. If a
second rearrangements file is given, repertoires from both files are pooled
per chain before pairing, so repertoires can also be compared across files
(always restricted to repertoires sharing the same chain/locus).
"""

import copy
import itertools
import logging
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from natsort import natsort_keygen

import tcrcloud.aminoacids as aminoacids
import tcrcloud.format
import tcrcloud.vgenes as vgenes
from tcrcloud.errors import TCRcloudError

logger = logging.getLogger(__name__)

_HTML_CONFIG = {"responsive": True, "displayModeBar": True, "displaylogo": False}


def _load_combined(
    rearrangements: str, rearrangements2: str | None = None
) -> tuple[pd.DataFrame, str]:
    """Load one or two rearrangements files into a single tagged DataFrame.

    Returns `(combined_df, output_prefix)`. When a second file is given,
    each file's `repertoire_id` values are prefixed with that file's stem
    (e.g. "fileA__rep1") so repertoires from different files never collide,
    even if they happen to share a raw repertoire_id.
    """

    df1 = tcrcloud.format.format_data(rearrangements).copy()
    stem1 = Path(rearrangements).stem

    if not rearrangements2:
        return df1, str(Path(rearrangements).with_suffix(""))

    df2 = tcrcloud.format.format_data(rearrangements2).copy()
    stem2 = Path(rearrangements2).stem

    df1["repertoire_id"] = stem1 + "__" + df1["repertoire_id"].astype(str)
    df2["repertoire_id"] = stem2 + "__" + df2["repertoire_id"].astype(str)

    combined = pd.concat([df1, df2], ignore_index=True)
    return combined, f"{stem1}_vs_{stem2}"


# ---------------------------------------------------------------------------
# V gene usage comparison (like tcrcloud.vgenes)
# ---------------------------------------------------------------------------


def _vgene_frequency_tables(
    formatted_vgene: pd.DataFrame,
    chain_letter: str,
    repA: str,
    repB: str,
    species: str,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[Any]] | None:
    """Build aligned V-gene x CDR3-length frequency (%) tables for repA/repB."""

    if chain_letter == "D":
        x_axis = dict(vgenes._palette_for_chain("D", species))
        for gene in vgenes._ALPHA_DELTA_VGENES:
            x_axis.pop(gene, None)
    else:
        x_axis = vgenes._palette_for_chain(chain_letter, species)

    if not x_axis:
        return None

    v_genes = sorted(x_axis.keys(), key=natsort_keygen())

    tables = {}
    for repertoire_id in (repA, repB):
        df = formatted_vgene[
            (formatted_vgene["chain"] == chain_letter)
            & (formatted_vgene["repertoire_id"] == repertoire_id)
        ].copy()

        if chain_letter == "D":
            df.replace(vgenes._ALPHA_DELTA_VGENES, vgenes._DELTA_VGENES, inplace=True)

        counts = (
            df.pivot_table(index=["v_call", "CDR3_length"], aggfunc="size")
            .reset_index()
            .rename(columns={0: "counts"})
        )
        total = counts["counts"].sum()
        pivot = (
            counts.pivot(index="v_call", columns="CDR3_length", values="counts")
            .reindex(index=v_genes, fill_value=0)
            .fillna(0)
        )
        tables[repertoire_id] = (pivot * 100 / total) if total else pivot * 0.0

    all_lengths = sorted(set().union(*(t.columns for t in tables.values())))
    aligned = {
        rid: t.reindex(columns=all_lengths, fill_value=0.0) for rid, t in tables.items()
    }
    return aligned[repA], aligned[repB], v_genes, all_lengths


def _vgene_comparison(
    formatted_vgene: pd.DataFrame,
    chain_letter: str,
    repA: str,
    repB: str,
    species: str,
    prefix: str,
    output_format: str,
    export: bool,
) -> None:
    tables = _vgene_frequency_tables(formatted_vgene, chain_letter, repA, repB, species)
    if tables is None:
        logger.warning(
            f"Skipping V genes comparison for {repA} vs {repB} ({chain_letter}): "
            "no built-in V-gene palette for this chain/species"
        )
        return
    tableA, tableB, v_genes, lengths = tables

    diff = tableA - tableB

    if export:
        diff_filename = (
            f"{prefix}_vgenes_diff_table_{repA}_vs_{repB}_{chain_letter}.csv"
        )
        diff.to_csv(diff_filename)
        logger.info("V genes comparison table saved as " + diff_filename)

    settings = vgenes._CHAIN_PLOT_SETTINGS.get(
        chain_letter, {"plot_aspect": (3, 1, 1), "x_size": 6}
    )
    x = list(range(len(v_genes)))
    z = [diff[col].values for col in diff.columns]
    max_abs = max(abs(diff.values.min()), abs(diff.values.max())) or 1.0

    fig = go.Figure(
        go.Surface(
            x=x,
            y=lengths,
            z=z,
            colorscale="Portland_r",
            cmin=-max_abs,
            cmax=max_abs,
            cmid=0,
        )
    )
    camera = dict(eye=dict(x=2.5, y=-3.5, z=2.5))
    sc = dict(
        aspectratio=dict(
            x=settings["plot_aspect"][0],
            y=settings["plot_aspect"][1],
            z=settings["plot_aspect"][2],
        ),
        xaxis_title=v_genes[0][:4],
        yaxis_title="CDR3 Length",
        zaxis_title=f"{repA} minus {repB} (percentage points)",
        xaxis=dict(
            tickmode="array",
            ticktext=v_genes,
            tickvals=x,
            tickfont=dict(size=settings["x_size"]),
            title=dict(font=dict(size=10)),
        ),
        yaxis=dict(tickfont=dict(size=8), title=dict(font=dict(size=10))),
        zaxis=dict(tickfont=dict(size=8), title=dict(font=dict(size=10))),
    )
    fig.update_layout(
        width=700,
        margin=dict(r=10, l=10, b=10, t=10),
        scene_camera=camera,
        scene=sc,
        template="plotly_white",
    )
    outputname = f"{prefix}_vgenes_diff_{repA}_vs_{repB}_{chain_letter}.{output_format}"
    fig.write_image(outputname, scale=6)
    logger.info("V genes comparison plot saved as " + outputname)

    fig_html = copy.deepcopy(fig)
    sc_html = dict(sc)
    sc_html["xaxis"] = dict(
        tickmode="array",
        ticktext=v_genes,
        tickvals=x,
        tickfont=dict(size=13),
        title=dict(font=dict(size=18)),
    )
    sc_html["yaxis"] = dict(tickfont=dict(size=14), title=dict(font=dict(size=18)))
    sc_html["zaxis"] = dict(tickfont=dict(size=14), title=dict(font=dict(size=18)))
    fig_html.update_layout(
        width=1280,
        height=800,
        scene=sc_html,
        template="plotly_white",
    )
    html_outputname = f"{prefix}_vgenes_diff_{repA}_vs_{repB}_{chain_letter}.html"
    fig_html.write_html(html_outputname, config=_HTML_CONFIG)
    logger.info("Interactive HTML plot saved as " + html_outputname)


# ---------------------------------------------------------------------------
# Amino acid composition comparison (like tcrcloud.aminoacids)
# ---------------------------------------------------------------------------


def _aminoacid_position_table(df: pd.DataFrame) -> pd.DataFrame:
    """Build a (20 amino acids) x (CDR3 position) frequency (%) table."""

    df = df.reset_index(drop=True)
    splitted = df["junction_aa"].str.split("", expand=True)
    new_df = pd.concat([df, splitted], axis=1)

    names = list(new_df.columns[5:])
    long_df = new_df.melt(
        id_vars=["counts"], value_vars=names, var_name="position", value_name="aa"
    )
    final = (
        long_df.groupby(["aa", "position"])["counts"]
        .sum()
        .unstack("position")
        .reindex(columns=names)
    )

    for aa in aminoacids.colours:
        if aa not in final.index:
            final.loc[aa] = np.nan
    final = final.sort_index()
    final.drop(final.head(1).index, inplace=True)  # drop blank row from split
    final.columns = names
    final.drop(columns=final.columns[-1:], inplace=True)  # drop trailing blank column
    final = final.replace(np.nan, 0)

    totals = final.sum(axis=0)
    normalized = (final * 100).div(totals.replace(0, np.nan), axis=1).fillna(0)
    return normalized.reindex(aminoacids.desired_order)


def _aminoacid_comparison(
    formatted_aminoacids: pd.DataFrame,
    chain_letter: str,
    repA: str,
    repB: str,
    prefix: str,
    output_format: str,
    export: bool,
) -> None:
    tables = {}
    for repertoire_id in (repA, repB):
        df = formatted_aminoacids[
            (formatted_aminoacids["chain"] == chain_letter)
            & (formatted_aminoacids["repertoire_id"] == repertoire_id)
        ]
        if df.empty:
            logger.warning(
                f"Skipping amino acids comparison for {repA} vs {repB} "
                f"({chain_letter}): no data for {repertoire_id}"
            )
            return
        tables[repertoire_id] = _aminoacid_position_table(df)

    all_positions = sorted(
        set().union(*(t.columns for t in tables.values())), key=lambda v: int(v)
    )
    tableA = tables[repA].reindex(columns=all_positions, fill_value=0.0)
    tableB = tables[repB].reindex(columns=all_positions, fill_value=0.0)
    diff = tableA - tableB

    if export:
        diff_filename = (
            f"{prefix}_aminoacids_diff_table_{repA}_vs_{repB}_{chain_letter}.csv"
        )
        diff.to_csv(diff_filename)
        logger.info("Amino acids comparison table saved as " + diff_filename)

    _aminoacid_3d_diff_plot(diff, repA, repB, chain_letter, prefix, output_format)
    _aminoacid_2d_diff_plots(diff, repA, repB, chain_letter, prefix, output_format)


def _aminoacid_3d_diff_plot(
    diff: pd.DataFrame,
    repA: str,
    repB: str,
    chain_letter: str,
    prefix: str,
    output_format: str,
) -> None:
    x_ordered = [aa for aa in aminoacids.desired_order if aa in diff.index]
    y_ordered = sorted(diff.columns, key=lambda v: int(v))
    diff_reordered = diff.loc[x_ordered, y_ordered]
    z_values = diff_reordered.values

    step = 1
    mesh_list = []
    for i, aa in enumerate(x_ordered):
        for j, _pos in enumerate(y_ordered):
            color_value = aminoacids.colours.get(aa, "gray")
            x_min = i * 2 * step
            x_max = x_min + step
            y_min = j * 2 * step
            y_max = y_min + step
            z_max = z_values[i, j]
            mesh_list.append(
                aminoacids.generate_mesh(
                    x_min, x_max, y_min, y_max, 0, z_max, color_value
                )
            )

    x_tickvals = [i * 2 * step + step / 2 for i in range(len(x_ordered))]
    y_tickvals = [j * 2 * step + step / 2 for j in range(len(y_ordered))]

    sc = dict(
        aspectratio=dict(x=1, y=1, z=1),
        xaxis_title="Amino acids",
        yaxis_title="CDR3 Length",
        zaxis_title=f"{repA} minus {repB} (%)",
        xaxis=dict(
            tickmode="array",
            ticktext=x_ordered,
            tickvals=x_tickvals,
            tickfont=dict(size=7),
            title=dict(font=dict(size=8)),
        ),
        yaxis=dict(
            tickmode="array",
            ticktext=[str(pos) for pos in y_ordered],
            tickvals=y_tickvals,
            tickfont=dict(size=7),
            title=dict(font=dict(size=8)),
        ),
        zaxis=dict(tickfont=dict(size=8), title=dict(font=dict(size=8))),
    )
    camera = dict(eye=dict(x=2.0, y=2.0, z=0.5))
    fig = go.Figure(mesh_list)
    fig.update_layout(
        width=700,
        margin=dict(r=10, l=10, b=10, t=10),
        scene_camera=camera,
        scene=sc,
        template="plotly_white",
    )
    outputname = (
        f"{prefix}_aminoacids3D_diff_{repA}_vs_{repB}_{chain_letter}.{output_format}"
    )
    fig.write_image(outputname, scale=6)
    logger.info("Tridimensional amino acids comparison plot saved as " + outputname)

    fig_html = copy.deepcopy(fig)
    sc_html = dict(
        aspectratio=dict(x=1, y=1, z=1),
        xaxis_title="Amino acids",
        yaxis_title="CDR3 Length",
        zaxis_title=f"{repA} minus {repB} (%)",
        xaxis=dict(
            tickmode="array",
            ticktext=x_ordered,
            tickvals=x_tickvals,
            tickfont=dict(size=13),
            title=dict(font=dict(size=18)),
        ),
        yaxis=dict(
            tickmode="array",
            ticktext=[str(pos) for pos in y_ordered],
            tickvals=y_tickvals,
            tickfont=dict(size=14),
            title=dict(font=dict(size=18)),
        ),
        zaxis=dict(tickfont=dict(size=14), title=dict(font=dict(size=18))),
    )
    fig_html.update_layout(
        width=1280, height=800, scene=sc_html, template="plotly_white"
    )
    html_outputname = f"{prefix}_aminoacids3D_diff_{repA}_vs_{repB}_{chain_letter}.html"
    fig_html.write_html(html_outputname, config=_HTML_CONFIG)
    logger.info("Interactive HTML plot saved as " + html_outputname)


def _aminoacid_2d_diff_plots(
    diff: pd.DataFrame,
    repA: str,
    repB: str,
    chain_letter: str,
    prefix: str,
    output_format: str,
) -> None:
    # (A) Diverging stacked bar chart across CDR3 positions: for each
    # position, positive percentage-point differences stack upward and
    # negative ones stack downward, colored per amino acid.
    comp_T = diff.transpose()  # rows = positions, columns = amino acids
    fig, ax = plt.subplots(figsize=(10, 6))
    positions = comp_T.index.tolist()
    cum_pos = np.zeros(len(comp_T))
    cum_neg = np.zeros(len(comp_T))

    for aa in aminoacids.desired_order:
        if aa in comp_T.columns:
            values = comp_T[aa].values
            pos_values = np.where(values > 0, values, 0)
            neg_values = np.where(values < 0, values, 0)
            ax.bar(
                positions,
                pos_values,
                bottom=cum_pos,
                color=aminoacids.colours[aa],
                label=aa,
            )
            cum_pos += pos_values
            ax.bar(positions, neg_values, bottom=cum_neg, color=aminoacids.colours[aa])
            cum_neg += neg_values

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("CDR3 Position")
    ax.set_ylabel("Difference in % (" + repA + " minus " + repB + ")")
    ax.set_title(f"Amino acid composition: {repA} vs {repB} ({chain_letter})")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    ax.set_xticks(positions)
    ax.set_xticklabels(positions)
    outname = (
        f"{prefix}_aminoacids2D_diff_{repA}_vs_{repB}_{chain_letter}.{output_format}"
    )
    plt.tight_layout()
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Amino acids 2D comparison plot saved as " + outname)

    # (B) Squashed per-amino-acid bar chart: sum the difference across all
    # CDR3 positions to give one net difference value per amino acid.
    sum_by_aa = diff.sum(axis=1).reindex(aminoacids.desired_order).fillna(0)
    fig2, ax2 = plt.subplots(figsize=(8, 4))
    colors_aa = [aminoacids.colours[a] for a in sum_by_aa.index]
    ax2.bar(sum_by_aa.index, sum_by_aa.values, color=colors_aa)
    ax2.axhline(0, color="black", linewidth=0.8)
    ax2.set_xlabel("Amino acids")
    ax2.set_ylabel("Difference in % (" + repA + " minus " + repB + ")")
    ax2.set_title(f"Sum across positions: {repA} vs {repB} ({chain_letter})")
    outname2 = (
        f"{prefix}_aminoacids2D_diff_squashedAA_{repA}_vs_{repB}_{chain_letter}"
        f".{output_format}"
    )
    plt.tight_layout()
    plt.savefig(outname2, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Amino acids squashed comparison plot saved as " + outname2)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def compare(
    rearrangements: str,
    rearrangements2: str | None = None,
    species: str = "homo_sapiens",
    export: bool = False,
    output_format: str = "png",
) -> None:
    """Entry point for the `TCRcloud compare` command."""

    # Determine the output image format (defaults to "png"). Validated here
    # too since `compare()` may be called directly as a library function
    # rather than only via the argparse CLI, whose `choices=["svg", "png"]`
    # wouldn't otherwise catch a bad value.
    output_format = output_format.strip().lower()
    if output_format not in ("svg", "png"):
        raise TCRcloudError(
            f"unsupported output format '{output_format}'. "
            "Please choose 'svg' or 'png'"
        )

    combined_df, prefix = _load_combined(rearrangements, rearrangements2)

    formatted_vgene = tcrcloud.format.format_vgene(combined_df)
    formatted_aminoacids = tcrcloud.format.format_aminoacids(combined_df)

    chain_groups = combined_df.groupby("chain")["repertoire_id"].unique()

    compared_any = False
    for chain_letter, reps in chain_groups.items():
        reps = sorted(reps)
        if len(reps) < 2:
            continue
        for repA, repB in itertools.combinations(reps, 2):
            compared_any = True
            logger.info(f"Comparing {repA} vs {repB} ({chain_letter} chain)")
            _vgene_comparison(
                formatted_vgene,
                str(chain_letter),
                repA,
                repB,
                species,
                prefix,
                output_format,
                export,
            )
            _aminoacid_comparison(
                formatted_aminoacids,
                str(chain_letter),
                repA,
                repB,
                prefix,
                output_format,
                export,
            )

    if not compared_any:
        raise TCRcloudError(
            "no chain had 2 or more repertoires to compare. "
            "Provide a rearrangements file with multiple repertoires of the "
            "same chain, or a second rearrangements file with a matching "
            "chain to compare across files"
        )
