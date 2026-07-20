import copy
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import plotly.graph_objects as go

import tcrcloud.format
import tcrcloud.colours

colours = tcrcloud.colours.aminoacids

desired_order = [
    "H",
    "K",
    "R",  # Polar Positively Charged (Basic)
    "D",
    "E",  # Polar Negatively Charged (Acidic)
    "C",
    "N",
    "Q",
    "S",
    "T",  # Polar Uncharged
    "Y",
    "W",
    "F",  # Aromatic Hydrophobic
    "G",
    "P",
    "A",
    "M",
    "V",
    "I",
    "L",  # Aliphatic Hydrophobic
]


def generate_mesh(x_min, x_max, y_min, y_max, z_min, z_max, color_value):
    y_max = y_max + 0.5
    x_max = x_max + 0.5
    mesh = go.Mesh3d(
        x=[x_min, x_min, x_max, x_max, x_min, x_min, x_max, x_max],
        y=[y_min, y_max, y_max, y_min, y_min, y_max, y_max, y_min],
        z=[z_min, z_min, z_min, z_min, z_max, z_max, z_max, z_max],
        color=color_value,
        i=[7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2],
        j=[3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
        k=[0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6],
        opacity=1.0,
        flatshading=True,
    )
    return mesh


def aminoacids(args):
    all_aa = pd.DataFrame(
        np.zeros((20, 1)),
        columns=["just_empty"],
        index=[
            "H",
            "K",
            "R",  # Polar Positively Charged (Basic)
            "D",
            "E",  # Polar Negatively Charged (Acidic)
            "C",
            "N",
            "Q",
            "S",
            "T",  # Polar Uncharged
            "Y",
            "W",
            "F",  # Aromatic Hydrophobic
            "G",
            "P",
            "A",
            "M",
            "V",
            "I",
            "L",  # Aliphatic Hydrophobic
        ],
    )

    # Get full dataset from the formatting module
    samples_df = tcrcloud.format.format_data(args)
    formatted_samples = tcrcloud.format.format_aminoacids(samples_df)
    samples = formatted_samples.groupby(["chain", "repertoire_id"])
    keys = [key for key, _ in samples]

    # Normalize boolean-style CLI flags (allow strings like "true"/"false").
    # argparse already handles this via str2bool, but aminoacids() may also
    # be called directly (e.g. in tests) with plain strings, so we
    # re-check defensively here.
    export = args.export
    if isinstance(export, str):
        export = export.lower() in ("yes", "true", "t", "y", "1")

    three_d = args.threeD
    if isinstance(three_d, str):
        three_d = three_d.lower() in ("yes", "true", "t", "y", "1")

    # Determine the output image format (defaults to "png"). Validated
    # here too since `aminoacids()` may be called directly (e.g. in
    # tests) rather than only via the argparse CLI, whose
    # `choices=["svg", "png"]` wouldn't otherwise catch a bad value.
    output_format = (getattr(args, "format", None) or "png").strip().lower()
    if output_format not in ("svg", "png"):
        raise ValueError(
            f"TCRcloud error: unsupported output format '{output_format}'. "
            "Please choose 'svg' or 'png'"
        )

    for j in keys:
        df = samples.get_group(j)
        splitted = df["junction_aa"].str.split("", expand=True)
        new_df = pd.concat([df, splitted], axis=1)

        positions = []
        names = []
        for i in new_df.columns[5:]:
            positions.append(new_df.groupby(i)[["counts"]].sum())
            names.append(i)
        final = pd.concat(positions, axis=1)

        # Ensure every amino acid is present, even if zero
        for i in colours.keys():
            if i not in list(final.index):
                final.loc[i] = np.nan
        final = final.sort_index()
        final.drop(final.head(1).index, inplace=True)  # drop blank from split
        final.columns = names
        final.drop(columns=final.columns[-1:], inplace=True)

        # Replace missing with zero
        final = final.replace(np.nan, 0)

        # Normalize to a percentage of reads per CDR3 position; the
        # CDR3-length axis always adapts to the data.
        normalized = final.apply(lambda x: x * 100 / sum(x), axis=0).fillna(0)
        normalized = normalized.reindex(
            all_aa.index, fill_value=0
        )  # fill missing a.a value to 0 if absent

        # If user wants a 2D stacked bar plot
        if not three_d:
            transposed_df = (
                normalized.transpose()
            )  # Convert for plotting, for 2D bar plot
            transposed_df.plot(
                kind="bar", stacked=True, color=colours, figsize=(10, 14)
            )
            outputname = (
                args.rearrangements[:-4]
                + "_aminoacids_"
                + j[1]
                + "_"
                + j[0]
                + "."
                + output_format
            )
            plt.legend(bbox_to_anchor=(1.01, 1), reverse=False, loc="upper left")
            plt.savefig(outputname, dpi=300, bbox_inches="tight")
            plt.close()
            print("Amino acids plot saved as " + outputname)

        # If user wants a 3D bar plot
        else:
            # If fewer than 20 rows, pad
            if len(normalized) < 20:
                pad_aa = pd.DataFrame(
                    np.zeros((20, 1)),
                    columns=["just_empty"],
                    index=[
                        "L",
                        "I",
                        "V",
                        "M",
                        "A",
                        "P",
                        "G",  # Aliphatic Hydrophobic
                        "F",
                        "W",
                        "Y",  # Aromatic Hydrophobic
                        "T",
                        "S",
                        "Q",
                        "N",
                        "C",  # Polar Uncharged
                        "E",
                        "D",  # Polar Negatively Charged (Acidic)
                        "R",
                        "K",
                        "H",  # Polar Positively Charged (Basic)
                    ],
                )
                result = pd.concat([normalized, pad_aa], axis=1)
                normalized = result.drop("just_empty", axis=1).fillna(0)
            # Flip so top row is last in the final plot
            normalized = normalized.reindex(desired_order[::-1])

            # build 3D bars
            x = normalized.index.tolist()  # Each amino acid
            y = np.array(normalized.columns)  # CDR3 positions
            z_vals = np.array(normalized.iloc[0])
            for i in range(1, len(normalized)):
                z_vals = np.append(z_vals, normalized.iloc[i].values)

            x_df = pd.Series(y)
            y_df = pd.Series(x)
            z_df = pd.Series(z_vals)

            x_min = 0
            y_min = 0
            z_min = 0
            step = 1
            mesh_list = []

            x_unique = x_df.unique()
            y_unique = y_df.unique()
            len_x_df_uniq = len(x_unique)

            for idx, x_data in enumerate(x_unique):
                for idx2, y_data in enumerate(y_unique):
                    color_value = colours[y_data]
                    x_max = x_min + step
                    y_max = y_min + step
                    z_max = z_df[idx + idx2 * len_x_df_uniq]
                    mesh_list.append(
                        generate_mesh(
                            x_min,
                            x_max,
                            y_min,
                            y_max,
                            z_min,
                            z_max,
                            color_value,
                        )
                    )
                    x_min += 2 * step
                y_min += 2 * step
                x_min = 0

            fig = go.Figure(mesh_list)
            camera = dict(eye=dict(x=2.0, y=2.0, z=2.0))
            sc = dict(
                aspectratio=dict(x=1, y=1, z=1),
                xaxis_title="Amino acids",
                yaxis_title="CDR3 Length",
                zaxis_title="Percentage of reads",
                xaxis=dict(
                    tickmode="array",
                    ticktext=desired_order,  # use the global order
                    tickvals=[v for v in range(1, 41) if v % 2 != 0],
                    tickfont=dict(size=7),
                    title=dict(font=dict(size=8)),
                ),
                yaxis=dict(
                    tickmode="array",
                    ticktext=[str(v) for v in range(1, y[-1] + 1)],
                    tickvals=[v for v in range(1, y[-1] * 2 + 1) if v % 2 != 0],
                    tickfont=dict(size=7),
                    title=dict(font=dict(size=8)),
                ),
                zaxis=dict(
                    tickfont=dict(size=8),
                    title=dict(font=dict(size=8)),
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
                args.rearrangements[:-4]
                + "_aminoacids3D_"
                + j[1]
                + "_"
                + j[0]
                + "."
                + output_format
            )
            fig.write_image(outputname, scale=6)

            # Also make interactive HTML version with bigger fonts
            fig_html = copy.deepcopy(fig)

            x_ordered = (
                [aa for aa in desired_order if aa in x]
                if "x" in locals()
                else desired_order
            )
            y_ordered = sorted(y, key=lambda v: int(v)) if "y" in locals() else []

            num_x = len(x_ordered)
            num_y = len(y_ordered)

            # Compute tick centers for the reordered data (if not already computed)
            x_tickvals = [i * 2 * step + step / 2 for i in range(num_x)]
            y_tickvals = [j * 2 * step + step / 2 for j in range(num_y)]

            # Define a scene for the HTML export with larger fonts:
            sc_html = dict(
                aspectratio=dict(x=1, y=1, z=1),
                xaxis_title="Amino acids",
                yaxis_title="CDR3 Length",
                zaxis_title="Percentage of reads",
                xaxis=dict(
                    tickmode="array",
                    ticktext=desired_order,
                    tickvals=[v for v in range(1, 41) if v % 2 != 0],
                    tickfont=dict(size=13),
                    title=dict(font=dict(size=18)),
                ),
                yaxis=dict(
                    tickmode="array",
                    ticktext=[str(v) for v in range(1, y[-1] + 1)],
                    tickvals=[v for v in range(1, y[-1] * 2 + 1) if v % 2 != 0],
                    tickfont=dict(size=14),
                    title=dict(font=dict(size=18)),
                ),
                zaxis=dict(
                    tickfont=dict(size=14),
                    title=dict(font=dict(size=18)),
                ),
            )
            camera_html = dict(eye=dict(x=2.0, y=2.0, z=2.1))
            fig_html.update_layout(
                width=1280, height=800, scene=sc_html, template="plotly_white"
            )

            # Export interactive HTML version (independent of the
            # chosen static image format).
            html_outputname = (
                args.rearrangements[:-4]
                + "_aminoacids3D_"
                + j[1]
                + "_"
                + j[0]
                + ".html"
            )
            fig_html.write_html(
                html_outputname,
                config={
                    "responsive": True,
                    "displayModeBar": True,
                    "displaylogo": False,
                },
            )

            print("Tridimensional Amino acids plot saved as " + outputname)
            print("Interactive HTML plot saved as", html_outputname)

        # Export the processed data to a CSV file
        if export:
            df_filename = (
                args.rearrangements[:-4]
                + "_aminoacids_table"
                + j[1]
                + "_"
                + j[0]
                + ".csv"
            )
            normalized.to_csv(df_filename, index=True)
