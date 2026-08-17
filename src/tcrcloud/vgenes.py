"""V gene surface plot generation for TCRcloud.

This module computes, per chain/repertoire, a table of V-gene usage broken
down by CDR3 length and renders it as a 3D surface plot (V gene x CDR3
length x percentage of reads).
"""

import copy
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from natsort import natsort_keygen

import tcrcloud.colours
import tcrcloud.format
from tcrcloud.errors import TCRcloudError

# V genes ambiguously shared between the alpha and delta loci; when a
# repertoire's delta chain is being plotted, these are renamed to their
# delta-locus equivalent so they group with the rest of the TRDV genes.
_ALPHA_DELTA_VGENES = [
    "TRAV14/DV4",
    "TRAV29/DV5",
    "TRAV23/DV6",
    "TRAV36/DV7",
    "TRAV38-2/DV8",
]
_DELTA_VGENES = ["TRDV4", "TRDV5", "TRDV6", "TRDV7", "TRDV8"]

# Per-chain plot styling: 3D aspect ratio and x-axis tick font size, tuned
# so that each chain's V-gene axis stays readable regardless of how many
# V genes it has.
_CHAIN_PLOT_SETTINGS = {
    "A": {"plot_aspect": (3.5, 1, 1), "x_size": 6},
    "B": {"plot_aspect": (3, 1, 1), "x_size": 4},
    "G": {"plot_aspect": (2, 1, 1), "x_size": 10},
    "D": {"plot_aspect": (1.5, 1, 1), "x_size": 8},
    "H": {"plot_aspect": (3.5, 1, 1), "x_size": 4},
    "K": {"plot_aspect": (3.5, 1, 1), "x_size": 4},
    "L": {"plot_aspect": (3.5, 1, 1), "x_size": 4},
}


def _palette_for_chain(chain_letter: str, species: str = "homo_sapiens"):
    mapping = {
        "A": "TRAV",
        "B": "TRBV",
        "G": "TRGV",
        "D": "TRDV",
        "H": "IGHV",
        "K": "IGKV",
        "L": "IGLV",
    }
    vgene_type = mapping.get(chain_letter.upper())
    if vgene_type is None:
        return {}
    return tcrcloud.colours.get_vgene_palette(vgene_type, species)


def get_table(keys, samples, args):
    """Build the per-chain/repertoire V-gene usage tables used for plotting."""

    species = getattr(args, "species", "homo_sapiens") or "homo_sapiens"

    datasets = []
    for chain_letter, repertoire_id in keys:
        settings = _CHAIN_PLOT_SETTINGS.get(chain_letter)
        if settings is None:
            continue

        if chain_letter == "D":
            x_axis = dict(_palette_for_chain("D", species))
            for gene in _ALPHA_DELTA_VGENES:
                x_axis.pop(gene, None)
        else:
            x_axis = _palette_for_chain(chain_letter, species)

        x_axis_ticks = list(range(len(x_axis)))

        df = samples.get_group((chain_letter, repertoire_id))
        jic = df.copy()

        if "D" in jic["chain"].values:
            jic.replace(_ALPHA_DELTA_VGENES, _DELTA_VGENES, inplace=True)

        new_df = jic.pivot_table(
            index=["v_call", "CDR3_length"], aggfunc="size"
        ).reset_index()
        new_df.rename(columns={0: "counts"}, inplace=True)

        # Create an empty df (0 counts for every V gene / CDR3 length
        # combination in range) to serve as the base grid for the surface
        # plot, so every V gene has an entry at every CDR3 length even if no
        # reads were observed there.
        empty_df = pd.DataFrame(columns=["v_call", "CDR3_length", "counts"])
        x_axis_names = []
        limitymin = new_df.loc[new_df["CDR3_length"].idxmin()].iloc[1] - 1
        limitymax = new_df.loc[new_df["CDR3_length"].idxmax()].iloc[1] + 1
        for v_genes in x_axis:
            x_axis_names.append(v_genes)
            for c in range(limitymin, limitymax):
                df_new_row = pd.DataFrame(
                    {"v_call": [v_genes], "CDR3_length": [c], "counts": [0]}
                )
                empty_df = pd.concat([empty_df, df_new_row])

        df_merged = pd.concat([new_df, empty_df], ignore_index=True, sort=True)

        final_df = df_merged.groupby(["v_call", "CDR3_length"]).sum().reset_index()
        final_df["frequency"] = 100 * (final_df["counts"] / final_df["counts"].sum())
        final_df = final_df.sort_values(
            by=["CDR3_length", "v_call"], key=natsort_keygen()
        )
        df_grouped = final_df.groupby(["v_call", "CDR3_length"]).sum()
        df_reset = df_grouped.reset_index()
        df_reformat = (
            df_reset.pivot(index="v_call", columns="CDR3_length", values="frequency")
            .reset_index()
            .rename_axis(index=None, columns=None)
        )
        df_sorted = df_reformat.sort_values(by=["v_call"], key=natsort_keygen())

        if args.export:
            df_filename = (
                str(Path(args.rearrangements).with_suffix(""))
                + "_vgenes_table"
                + repertoire_id
                + "_"
                + chain_letter
                + ".csv"
            )
            df_sorted.to_csv(df_filename, index=False)

        x = df_sorted["v_call"].factorize()[0]
        y = np.array(df_sorted.columns.values.tolist()[1:])
        df_transpose = df_sorted.transpose()
        z = np.array(df_transpose.values.tolist()[1])
        for i in range(2, len(df_sorted.columns)):
            z = np.append(z, np.array(df_transpose.values.tolist()[i]))
        z = np.array_split(z, len(df_sorted.columns) - 1)

        datasets.append(
            {
                "x": x,
                "y": y,
                "z": z,
                "plot_aspect": settings["plot_aspect"],
                "x_size": settings["x_size"],
                "x_axis_ticks": x_axis_ticks,
                "x_axis_names": x_axis_names,
                "chain": chain_letter,
                "repertoire_id": repertoire_id,
            }
        )

    return datasets


def barplot(args):
    # Normalize boolean-style CLI flags (allow strings like "true"/"false").
    # argparse already handles this via str2bool when barplot() is invoked
    # through the CLI, but barplot() may also be called directly (e.g. in
    # tests) with plain strings, so we re-check defensively here.
    if isinstance(args.export, str):
        args.export = args.export.lower() in ("yes", "true", "t", "y", "1")

    # Determine the output image format (defaults to "png"). Validated here
    # too since `barplot()` may be called directly (e.g. in tests) rather
    # than only via the argparse CLI, whose `choices=["svg", "png"]` wouldn't
    # otherwise catch a bad value.
    output_format = (getattr(args, "format", None) or "png").strip().lower()
    if output_format not in ("svg", "png"):
        raise TCRcloudError(
            f"unsupported output format '{output_format}'. "
            "Please choose 'svg' or 'png'"
        )

    samples_df = tcrcloud.format.format_data(args)
    formatted_samples = tcrcloud.format.format_vgene(samples_df)

    samples = formatted_samples.groupby(["chain", "repertoire_id"])
    keys = list(samples.groups.keys())
    datasets = get_table(keys, samples, args)

    if not datasets:
        raise TCRcloudError("no repertoires found for plotting")

    for d in datasets:
        fig = go.Figure(go.Surface(x=d["x"], y=d["y"], z=d["z"], colorscale="Turbo"))
        camera = dict(eye=dict(x=2.5, y=-3.5, z=2.5))

        sc = dict(
            aspectratio=dict(
                x=d["plot_aspect"][0], y=d["plot_aspect"][1], z=d["plot_aspect"][2]
            ),
            xaxis_title=d["x_axis_names"][0][:4],
            yaxis_title="CDR3 Length",
            zaxis_title="Percentage of reads",
            xaxis=dict(
                tickmode="array",
                ticktext=d["x_axis_names"],
                tickvals=d["x_axis_ticks"],
                tickfont=dict(size=d["x_size"]),
                title=dict(font=dict(size=10)),
            ),
            yaxis=dict(
                tickfont=dict(size=8),
                title=dict(font=dict(size=10)),
            ),
            zaxis=dict(
                tickfont=dict(size=8),
                title=dict(font=dict(size=10)),
            ),
        )

        fig.update_layout(
            width=700,
            margin=dict(r=10, l=10, b=10, t=10),
            scene_camera=camera,
            scene=sc,
            template="plotly_white",
        )
        outputname = (
            str(Path(args.rearrangements).with_suffix(""))
            + "_vgenes_"
            + d["repertoire_id"]
            + "_"
            + d["chain"]
            + "."
            + output_format
        )
        fig.write_image(outputname, scale=6)
        print("V genes plot saved as " + outputname)

        # Also make an interactive HTML version with bigger fonts, independent
        # of the chosen static image format.
        fig_html = copy.deepcopy(fig)
        sc_html = dict(
            aspectratio=dict(
                x=d["plot_aspect"][0], y=d["plot_aspect"][1], z=d["plot_aspect"][2]
            ),
            xaxis_title=d["x_axis_names"][0][:4],
            yaxis_title="CDR3 Length",
            zaxis_title="Percentage of reads",
            xaxis=dict(
                tickmode="array",
                ticktext=d["x_axis_names"],
                tickvals=d["x_axis_ticks"],
                tickfont=dict(size=13),
                title=dict(font=dict(size=18)),
            ),
            yaxis=dict(
                tickfont=dict(size=14),
                title=dict(font=dict(size=18)),
            ),
            zaxis=dict(
                tickfont=dict(size=14),
                title=dict(font=dict(size=18)),
            ),
        )
        fig_html.update_layout(
            width=1280,
            height=800,
            scene_camera=camera,
            scene=sc_html,
            template="plotly_white",
        )
        html_outputname = (
            str(Path(args.rearrangements).with_suffix(""))
            + "_vgenes_"
            + d["repertoire_id"]
            + "_"
            + d["chain"]
            + ".html"
        )
        fig_html.write_html(
            html_outputname,
            config={"responsive": True, "displayModeBar": True, "displaylogo": False},
        )
        print("Interactive HTML plot saved as " + html_outputname)
