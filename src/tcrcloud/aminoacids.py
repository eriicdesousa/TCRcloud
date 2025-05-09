import sys
import copy
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio

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

    # If by_length flag is set, process each CDR3 length group separately.
    if args.by_length:
        # Compute CDR3 length from the junction_aa column
        samples_df["cdr3_length"] = samples_df["junction_aa"].apply(len)
        # Use a "counts" column if available; otherwise use the number of rows
        total_reads = (
            samples_df["counts"].sum()
            if "counts" in samples_df.columns
            else len(samples_df)
        )
        # Process each length group (filtering groups with <2.5% of total reads)
        for length, group in samples_df.groupby("cdr3_length"):
            group_reads = (
                group["counts"].sum() if "counts" in group.columns else len(group)
            )
            if group_reads / total_reads < 0.025:
                sys.stderr.write(
                    f"Skipping CDR3 length {length} due to low read count.\n"
                )
                continue
            print(f"Processing CDR3 length {length}")

            # Format the subgroup as for grouped
            formatted_samples = tcrcloud.format.format_aminoacids(group)
            samples = formatted_samples.groupby(["chain", "repertoire_id"])
            keys = [key for key, _ in samples]

            # Prepare dictionary for comparing 2 repertoires (if -c True)
            for_comparison = {
                "A": [],
                "B": [],
                "G": [],
                "D": [],
                "H": [],
                "K": [],
                "L": [],
            }

            for j in keys:
                df = samples.get_group(j)
                splitted = df["junction_aa"].str.split("", expand=True)
                new_df = pd.concat([df, splitted], axis=1)

                positions = []
                names = []
                for i in new_df.columns[5:]:
                    positions.append(new_df.groupby(i).sum("counts"))
                    names.append(i)
                final = pd.concat(positions, axis=1)

                # Ensure every amino acid is present, even if zero
                for i in colours.keys():
                    if i not in list(final.index):
                        final.loc[i] = np.nan
                final = final.sort_index()
                final.drop(final.head(1).index, inplace=True)  # drop blank from split
                final.columns = names
                final.drop(columns=final.columns[-1:], axis=1, inplace=True)

                # If user specified a max CDR3 length, trim/pad
                if args.length is not None:
                    if len(final.columns) > args.length:
                        diff = len(final.columns) - args.length
                        final.drop(columns=final.columns[-diff:], axis=1, inplace=True)
                    elif len(final.columns) < args.length:
                        diff = args.length - len(final.columns)
                        for number in range(1, diff + 1):
                            final[len(final.columns) + 1] = np.nan

                # Replace missing with zero
                final = final.replace(np.nan, 0)

                # --------------------------
                # Single-sample code (no comparison):
                # --------------------------
                if args.compare.lower() == "false":
                    # Make single-sample normalized data
                    normalized = final.apply(lambda x: x * 100 / sum(x), axis=0).fillna(
                        0
                    )
                    normalized = normalized.reindex(
                        all_aa.index, fill_value=0
                    )  # fill missing a.a value to 0 if absent
                    transposed_df = (
                        normalized.transpose()
                    )  # Convert for plotting, for 2D bar plot

                    # If user wants single-sample 2D stacked bar
                    if args.threeD.lower() == "false":
                        transposed_df.plot(
                            kind="bar", stacked=True, color=colours, figsize=(10, 14)
                        )
                        outputname = (
                            args.rearrangements[:-4]
                            + "_aminoacids_"
                            + j[1]
                            + "_"
                            + j[0]
                            + f"_L{length}.png"
                        )
                        plt.legend(
                            bbox_to_anchor=(1.01, 1), reverse=False, loc="upper left"
                        )
                        plt.savefig(outputname, dpi=300, bbox_inches="tight")
                        plt.close()
                        print("Amino acids plot saved as " + outputname)

                    # If user wants single sample 3D plot
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
                                tickvals=[
                                    v for v in range(1, y[-1] * 2 + 1) if v % 2 != 0
                                ],
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

                        # Save PNG
                        outputname = (
                            args.rearrangements[:-4]
                            + "_aminoacids3D_"
                            + j[1]
                            + "_"
                            + j[0]
                            + f"_L{length}.png"
                        )
                        fig.write_image(outputname, scale=6)

                        # Also make interactive HTML version with bigger fonts
                        fig_html = copy.deepcopy(fig)

                        x_ordered = (
                            [aa for aa in desired_order if aa in x]
                            if "x" in locals()
                            else desired_order
                        )
                        y_ordered = (
                            sorted(y, key=lambda v: int(v)) if "y" in locals() else []
                        )

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
                                tickvals=[
                                    v for v in range(1, y[-1] * 2 + 1) if v % 2 != 0
                                ],
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
                            width=1920,
                            height=1080,
                            scene=sc_html,
                            template="plotly_white",
                        )

                        # Export interactive HTML version
                        html_outputname = outputname.replace(".png", ".html")
                        fig_html.write_html(html_outputname)

                        print("Tridimensional Amino acids plot saved as " + outputname)
                        print("Interactive HTML plot saved as", html_outputname)

                    # --------------------------
                    # If compare == True => store normalized data for difference
                    # (and skip the single-sample plots)
                    # --------------------------
                else:  # args.compare.lower() == "true"
                    normalized = final.apply(lambda x: x * 100 / sum(x), axis=0).fillna(
                        0
                    )
                    normalized = normalized.reindex(all_aa.index, fill_value=0)
                    # Just store it for later difference
                    for_comparison[j[0]].append([normalized, j[1] + "_" + j[0]])

                # exports the processed dara to a CSV file
                if args.export.lower() == "true":
                    df_filename = (
                        args.rearrangements[:-4]
                        + "_aminoacids_table"
                        + j[1]
                        + "_"
                        + j[0]
                        + f"_L{length}.csv"
                    )
                    normalized.to_csv(df_filename, index=True)

            # End of loop over keys

            # --------------------------
            # 2-repertoires sample comparison if requested
            # --------------------------
            if args.compare.lower() == "true":
                comparisons = []  # initialize comparisons list
                for chain_id, data_list in for_comparison.items():
                    if len(data_list) < 2:
                        sys.stderr.write(
                            f"Less than 2 repertoires from the {chain_id} chain were detected in the rearrangements file\n"
                        )
                    elif len(data_list) > 2:
                        sys.stderr.write(
                            f"More than 2 repertoires from the {chain_id} chain were detected in the rearrangements file\n"
                        )
                    elif len(data_list) == 2:
                        # data_list[0] = [normalized_df, "SampleA"], data_list[1] = [normalized_df, "SampleB"]
                        dfA, labelA = data_list[0]
                        dfB, labelB = data_list[1]
                        comparison1 = dfA - dfB
                        comparison2 = dfB - dfA
                        comparisons.append([comparison1, labelA])  # A minus B
                        comparisons.append([comparison2, labelB])  # B minus A

                # If user wants a 3D difference plot:
                if args.threeD.lower() == "true":
                    for diff_df, label in comparisons:
                        # Order the axes explicitly from diff_df:
                        x = diff_df.index.tolist()  # amino acid labels
                        y = diff_df.columns.tolist()  # CDR3 positions

                        # Order x according to desired_order (preserving only those present)
                        x_ordered = [aa for aa in desired_order if aa in x]
                        # Sort y numerically if possible (or alphabetically otherwise)
                        try:
                            y_ordered = sorted(y, key=lambda v: int(v))
                        except:
                            y_ordered = sorted(y)

                        # Reorder the diff_df accordingly:
                        diff_df_reordered = diff_df.loc[x_ordered, y_ordered]
                        Z = (
                            diff_df_reordered.values
                        )  # shape: (num_amino, num_positions)
                        step = 1
                        num_x = len(x_ordered)
                        num_y = len(y_ordered)
                        mesh_list = []

                        # Build a mesh for each cell of the matrix:
                        for i in range(num_x):
                            for j in range(num_y):
                                z_val = Z[i, j]
                                # Color based on the amino acid (x_ordered[i])
                                cval = colours.get(x_ordered[i], "gray")
                                x_min = i * 2 * step
                                x_max = x_min + step
                                y_min = j * 2 * step
                                y_max = y_min + step
                                z_min = 0
                                z_max = z_val
                                mesh_list.append(
                                    generate_mesh(
                                        x_min, x_max, y_min, y_max, z_min, z_max, cval
                                    )
                                )

                        # Compute tick values for each axis (center of each cell)
                        x_tickvals = [i * 2 * step + step / 2 for i in range(num_x)]
                        y_tickvals = [j * 2 * step + step / 2 for j in range(num_y)]

                        sc_png = dict(
                            aspectratio=dict(x=1, y=1, z=1),
                            xaxis_title="Amino acids",
                            yaxis_title="CDR3 Length",
                            zaxis_title="Difference in %",
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
                            zaxis=dict(
                                tickfont=dict(size=8),
                                title=dict(font=dict(size=8)),
                            ),
                        )
                        camera_png = dict(eye=dict(x=2.0, y=2.0, z=0.5))
                        fig = go.Figure(mesh_list)
                        fig.update_layout(
                            width=700,
                            margin=dict(r=10, l=10, b=10, t=10),
                            scene_camera=camera_png,
                            scene=sc_png,
                            template="plotly_white",
                        )
                        outputname = (
                            args.rearrangements[:-4]
                            + "_aminoacids3D_"
                            + label
                            + "_comparison"
                            + f"_L{length}.png"
                        )
                        fig.write_image(outputname, scale=6)
                        print("Tridimensional Amino acids plot saved as", outputname)

                        # --- HTML export --- larger fonts
                        sc_html = dict(
                            aspectratio=dict(x=1, y=1, z=1),
                            xaxis_title="Amino acids",
                            yaxis_title="CDR3 Length",
                            zaxis_title="Difference in %",
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
                            zaxis=dict(
                                tickfont=dict(size=14),
                                title=dict(font=dict(size=18)),
                            ),
                        )
                        # For HTML camera is slightly shifted (if desired, adjust here)
                        camera_html = dict(eye=dict(x=2.0, y=2.0, z=2.1))
                        fig_html = copy.deepcopy(fig)
                        fig_html.update_layout(
                            width=1920,
                            height=1080,
                            scene_camera=camera_html,
                            scene=sc_html,
                            template="plotly_white",
                        )
                        html_outputname = outputname.replace(".png", ".html")
                        fig_html.write_html(html_outputname)
                        print("Interactive HTML plot saved as", html_outputname)

                # If user wants a 2D difference
                else:
                    for diff_df, label in comparisons:
                        # diff_df is the difference in % (rows=aa, cols=positions)

                        # (A) "Diverging stacked bar" across positions
                        # Exactly like 3D but flattened onto the "position" axis
                        # one bar per position, stacked by amino acids
                        comp_T = diff_df.transpose()  # now rows=positions, cols=aa
                        comp_T = comp_T.loc[
                            :, comp_T.columns.str.strip() != ""
                        ]  # Drop empty column if present

                        # Export the underlying table to CSV only if export flag is True
                        if args.export.lower() == "true":
                            csv_filename = (
                                args.rearrangements[:-4]
                                + "_aminoacids2D_"
                                + label
                                + "_comparison_table"
                                + f"_L{length}.csv"
                            )
                            comp_T.to_csv(csv_filename, index_label="Position")
                            print("CSV table saved as", csv_filename)

                        fig, ax = plt.subplots(figsize=(10, 6))
                        positions = comp_T.index.tolist()
                        cum_pos = np.zeros(len(comp_T))
                        cum_neg = np.zeros(len(comp_T))

                        for aa in desired_order:
                            if aa in comp_T.columns:
                                values = comp_T[aa].values
                                pos_values = np.where(values > 0, values, 0)
                                neg_values = np.where(values < 0, values, 0)

                                # Positive bars
                                ax.bar(
                                    positions,
                                    pos_values,
                                    bottom=cum_pos,
                                    color=colours[aa],
                                    label=(
                                        aa
                                        if aa not in ax.get_legend_handles_labels()[1]
                                        else ""
                                    ),
                                )
                                cum_pos += pos_values

                                # Negative bars
                                ax.bar(
                                    positions,
                                    neg_values,
                                    bottom=cum_neg,
                                    color=colours[aa],
                                )
                                cum_neg += neg_values

                        ax.axhline(0, color="black", linewidth=0.8)
                        ax.set_xlabel("CDR3 Position")
                        ax.set_ylabel("Difference in %")
                        ax.set_title("2D Comparison (stacked): " + label)
                        ax.legend(
                            bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8
                        )
                        ax.set_xticks(positions)
                        ax.set_xticklabels(positions)
                        outname = (
                            args.rearrangements[:-4]
                            + "_aminoacids2D_"
                            + label
                            + "_comparison_stacked"
                            + f"_L{length}.png"
                        )
                        plt.tight_layout()
                        plt.savefig(outname, dpi=300, bbox_inches="tight")
                        plt.close()
                        print("2D stacked difference plot saved as", outname)

                        # (B) Squash across positions => one bar per amino acid
                        # sum of all positions for each aa
                        sum_by_aa = diff_df.sum(axis=1)  # row-wise sum
                        sum_by_aa = sum_by_aa.reindex(desired_order).fillna(0)

                        # Export the squashed differences to CSV only if export flag is True
                        if args.export.lower() == "true":
                            df_sum_by_aa = sum_by_aa.to_frame(name="Difference")
                            df_sum_by_aa.index.name = "Amino Acid"
                            csv_filename2 = (
                                args.rearrangements[:-4]
                                + "_aminoacids2D_"
                                + label
                                + "_squashedAA_table"
                                + f"_L{length}.csv"
                            )
                            df_sum_by_aa.to_csv(csv_filename2)
                            print("CSV table saved as", csv_filename2)

                        fig2, ax2 = plt.subplots(figsize=(8, 4))
                        colors_aa = [colours[a] for a in sum_by_aa.index]
                        ax2.bar(sum_by_aa.index, sum_by_aa.values, color=colors_aa)
                        ax2.axhline(0, color="black", linewidth=0.8)
                        ax2.set_xlabel("Amino acids")
                        ax2.set_ylabel("Difference in %")
                        ax2.set_title("Sum across positions: " + label)
                        outname2 = (
                            args.rearrangements[:-4]
                            + "_aminoacids2D_"
                            + label
                            + "_squashedAA"
                            + f"_L{length}.png"
                        )
                        plt.tight_layout()
                        plt.savefig(outname2, dpi=300, bbox_inches="tight")
                        plt.close()
                        print("2D squashed-by-amino-acid plot saved as", outname2)

    # End processing for --by_length
    else:
        # Original processing when --by_length flag is not set
        samples_df = tcrcloud.format.format_data(args)
        formatted_samples = tcrcloud.format.format_aminoacids(samples_df)
        samples = formatted_samples.groupby(["chain", "repertoire_id"])
        keys = [key for key, _ in samples]

        # Prepare dictionary for comparing 2 repertoires (if -c True)
        for_comparison = {"A": [], "B": [], "G": [], "D": [], "H": [], "K": [], "L": []}

        for j in keys:
            df = samples.get_group(j)
            splitted = df["junction_aa"].str.split("", expand=True)
            new_df = pd.concat([df, splitted], axis=1)

            positions = []
            names = []
            for i in new_df.columns[5:]:
                positions.append(new_df.groupby(i).sum("counts"))
                names.append(i)
            final = pd.concat(positions, axis=1)

            # Ensure every amino acid is present, even if zero
            for i in colours.keys():
                if i not in list(final.index):
                    final.loc[i] = np.nan
            final = final.sort_index()
            final.drop(final.head(1).index, inplace=True)  # drop blank from split
            final.columns = names
            final.drop(columns=final.columns[-1:], axis=1, inplace=True)

            # If user specified a max CDR3 length, trim/pad
            if args.length is not None:
                if len(final.columns) > args.length:
                    diff = len(final.columns) - args.length
                    final.drop(columns=final.columns[-diff:], axis=1, inplace=True)
                elif len(final.columns) < args.length:
                    diff = args.length - len(final.columns)
                    for number in range(1, diff + 1):
                        final[len(final.columns) + 1] = np.nan

            # Replace missing with zero
            final = final.replace(np.nan, 0)

            # --------------------------
            # Single-sample code (no compare):
            # --------------------------
            if args.compare.lower() == "false":
                # Make single-sample normalized data
                normalized = final.apply(lambda x: x * 100 / sum(x), axis=0).fillna(0)
                normalized = normalized.reindex(
                    all_aa.index, fill_value=0
                )  # fill missing a.a value to 0 if absent
                transposed_df = (
                    normalized.transpose()
                )  # Convert for plotting, for 2D bar plot

                # If user wants single-sample 2D stacked bar
                if args.threeD.lower() == "false":
                    transposed_df.plot(
                        kind="bar", stacked=True, color=colours, figsize=(10, 14)
                    )
                    outputname = (
                        args.rearrangements[:-4]
                        + "_aminoacids_"
                        + j[1]
                        + "_"
                        + j[0]
                        + ".png"
                    )
                    plt.legend(
                        bbox_to_anchor=(1.01, 1), reverse=False, loc="upper left"
                    )
                    plt.savefig(outputname, dpi=300, bbox_inches="tight")
                    plt.close()
                    print("Amino acids plot saved as " + outputname)

                # If user wants single sample 3D plot
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

                    # Save PNG
                    outputname = (
                        args.rearrangements[:-4]
                        + "_aminoacids3D_"
                        + j[1]
                        + "_"
                        + j[0]
                        + ".png"
                    )
                    fig.write_image(outputname, scale=6)

                    # Also make interactive HTML version with bigger fonts
                    fig_html = copy.deepcopy(fig)

                    x_ordered = (
                        [aa for aa in desired_order if aa in x]
                        if "x" in locals()
                        else desired_order
                    )
                    y_ordered = (
                        sorted(y, key=lambda v: int(v)) if "y" in locals() else []
                    )

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
                        width=1920, height=1080, scene=sc_html, template="plotly_white"
                    )

                    # Export interactive HTML version
                    html_outputname = outputname.replace(".png", ".html")
                    fig_html.write_html(html_outputname)

                    print("Tridimensional Amino acids plot saved as " + outputname)
                    print("Interactive HTML plot saved as", html_outputname)

                # --------------------------
                # If compare == True => store normalized data for difference
                # (and skip the single-sample plots)
                # --------------------------
            else:  # args.compare.lower() == "true"
                normalized = final.apply(lambda x: x * 100 / sum(x), axis=0).fillna(0)
                normalized = normalized.reindex(all_aa.index, fill_value=0)
                # Just store it for later difference
                for_comparison[j[0]].append([normalized, j[1] + "_" + j[0]])

            # exports the processed dara to a CSV file
            if args.export.lower() == "true":
                df_filename = (
                    args.rearrangements[:-4]
                    + "_aminoacids_table"
                    + j[1]
                    + "_"
                    + j[0]
                    + ".csv"
                )
                normalized.to_csv(df_filename, index=True)

        # End of loop over keys

        # --------------------------
        # 2-repertoires sample comparison if requested
        # --------------------------
        if args.compare.lower() == "true":
            comparisons = []  # initialize comparisons list
            for chain_id, data_list in for_comparison.items():
                if len(data_list) < 2:
                    sys.stderr.write(
                        f"Less than 2 repertoires from the {chain_id} chain were detected in the rearrangements file\n"
                    )
                elif len(data_list) > 2:
                    sys.stderr.write(
                        f"More than 2 repertoires from the {chain_id} chain were detected in the rearrangements file\n"
                    )
                elif len(data_list) == 2:
                    # data_list[0] = [normalized_df, "SampleA"], data_list[1] = [normalized_df, "SampleB"]
                    dfA, labelA = data_list[0]
                    dfB, labelB = data_list[1]
                    comparison1 = dfA - dfB
                    comparison2 = dfB - dfA
                    comparisons.append([comparison1, labelA])  # A minus B
                    comparisons.append([comparison2, labelB])  # B minus A

            # If user wants a 3D difference plot:
            if args.threeD.lower() == "true":
                for diff_df, label in comparisons:
                    # Order the axes explicitly from diff_df:
                    x = diff_df.index.tolist()  # amino acid labels
                    y = diff_df.columns.tolist()  # CDR3 positions

                    # Order x according to desired_order (preserving only those present)
                    x_ordered = [aa for aa in desired_order if aa in x]
                    # Sort y numerically if possible (or alphabetically otherwise)
                    try:
                        y_ordered = sorted(y, key=lambda v: int(v))
                    except:
                        y_ordered = sorted(y)

                    # Reorder the diff_df accordingly:
                    diff_df_reordered = diff_df.loc[x_ordered, y_ordered]
                    Z = diff_df_reordered.values  # shape: (num_amino, num_positions)
                    step = 1
                    num_x = len(x_ordered)
                    num_y = len(y_ordered)
                    mesh_list = []

                    # Build a mesh for each cell of the matrix:
                    for i in range(num_x):
                        for j in range(num_y):
                            z_val = Z[i, j]
                            # Color based on the amino acid (x_ordered[i])
                            cval = colours.get(x_ordered[i], "gray")
                            x_min = i * 2 * step
                            x_max = x_min + step
                            y_min = j * 2 * step
                            y_max = y_min + step
                            z_min = 0
                            z_max = z_val
                            mesh_list.append(
                                generate_mesh(
                                    x_min, x_max, y_min, y_max, z_min, z_max, cval
                                )
                            )

                    # Compute tick values for each axis (center of each cell)
                    x_tickvals = [i * 2 * step + step / 2 for i in range(num_x)]
                    y_tickvals = [j * 2 * step + step / 2 for j in range(num_y)]

                    sc_png = dict(
                        aspectratio=dict(x=1, y=1, z=1),
                        xaxis_title="Amino acids",
                        yaxis_title="CDR3 Length",
                        zaxis_title="Difference in %",
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
                        zaxis=dict(
                            tickfont=dict(size=8),
                            title=dict(font=dict(size=8)),
                        ),
                    )
                    camera_png = dict(eye=dict(x=2.0, y=2.0, z=0.5))
                    fig = go.Figure(mesh_list)
                    fig.update_layout(
                        width=700,
                        margin=dict(r=10, l=10, b=10, t=10),
                        scene_camera=camera_png,
                        scene=sc_png,
                        template="plotly_white",
                    )
                    outputname = (
                        args.rearrangements[:-4]
                        + "_aminoacids3D_"
                        + label
                        + "_comparison.png"
                    )
                    fig.write_image(outputname, scale=6)
                    print("Tridimensional Amino acids plot saved as", outputname)

                    # --- HTML export --- larger fonts
                    sc_html = dict(
                        aspectratio=dict(x=1, y=1, z=1),
                        xaxis_title="Amino acids",
                        yaxis_title="CDR3 Length",
                        zaxis_title="Difference in %",
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
                        zaxis=dict(
                            tickfont=dict(size=14),
                            title=dict(font=dict(size=18)),
                        ),
                    )
                    # For HTML camera is slightly shifted (if desired, adjust here)
                    camera_html = dict(eye=dict(x=2.0, y=2.0, z=2.1))
                    fig_html = copy.deepcopy(fig)
                    fig_html.update_layout(
                        width=1920,
                        height=1080,
                        scene_camera=camera_html,
                        scene=sc_html,
                        template="plotly_white",
                    )
                    html_outputname = outputname.replace(".png", ".html")
                    fig_html.write_html(html_outputname)
                    print("Interactive HTML plot saved as", html_outputname)

            # If user wants a 2D difference
            else:
                for diff_df, label in comparisons:
                    # diff_df is the difference in % (rows=aa, cols=positions)

                    # (A) "Diverging stacked bar" across positions
                    # Exactly like 3D but flattened onto the "position" axis
                    # one bar per position, stacked by amino acids
                    comp_T = diff_df.transpose()  # now rows=positions, cols=aa
                    comp_T = comp_T.loc[
                        :, comp_T.columns.str.strip() != ""
                    ]  # Drop empty column if present

                    # Export the underlying table to CSV only if export flag is True
                    if args.export.lower() == "true":
                        csv_filename = (
                            args.rearrangements[:-4]
                            + "_aminoacids2D_"
                            + label
                            + "_comparison_table.csv"
                        )
                        comp_T.to_csv(csv_filename, index_label="Position")
                        print("CSV table saved as", csv_filename)

                    fig, ax = plt.subplots(figsize=(10, 6))
                    positions = comp_T.index.tolist()
                    cum_pos = np.zeros(len(comp_T))
                    cum_neg = np.zeros(len(comp_T))

                    for aa in desired_order:
                        if aa in comp_T.columns:
                            values = comp_T[aa].values
                            pos_values = np.where(values > 0, values, 0)
                            neg_values = np.where(values < 0, values, 0)

                            # Positive bars
                            ax.bar(
                                positions,
                                pos_values,
                                bottom=cum_pos,
                                color=colours[aa],
                                label=(
                                    aa
                                    if aa not in ax.get_legend_handles_labels()[1]
                                    else ""
                                ),
                            )
                            cum_pos += pos_values

                            # Negative bars
                            ax.bar(
                                positions, neg_values, bottom=cum_neg, color=colours[aa]
                            )
                            cum_neg += neg_values

                    ax.axhline(0, color="black", linewidth=0.8)
                    ax.set_xlabel("CDR3 Position")
                    ax.set_ylabel("Difference in %")
                    ax.set_title("2D Comparison (stacked): " + label)
                    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
                    ax.set_xticks(positions)
                    ax.set_xticklabels(positions)
                    outname = (
                        args.rearrangements[:-4]
                        + "_aminoacids2D_"
                        + label
                        + "_comparison_stacked.png"
                    )
                    plt.tight_layout()
                    plt.savefig(outname, dpi=300, bbox_inches="tight")
                    plt.close()
                    print("2D stacked difference plot saved as", outname)

                    # (B) Squash across positions => one bar per amino acid
                    # sum of all positions for each aa
                    sum_by_aa = diff_df.sum(axis=1)  # row-wise sum
                    sum_by_aa = sum_by_aa.reindex(desired_order).fillna(0)

                    # Export the squashed differences to CSV only if export flag is True
                    if args.export.lower() == "true":
                        df_sum_by_aa = sum_by_aa.to_frame(name="Difference")
                        df_sum_by_aa.index.name = "Amino Acid"
                        csv_filename2 = (
                            args.rearrangements[:-4]
                            + "_aminoacids2D_"
                            + label
                            + "_squashedAA_table.csv"
                        )
                        df_sum_by_aa.to_csv(csv_filename2)
                        print("CSV table saved as", csv_filename2)

                    fig2, ax2 = plt.subplots(figsize=(8, 4))
                    colors_aa = [colours[a] for a in sum_by_aa.index]
                    ax2.bar(sum_by_aa.index, sum_by_aa.values, color=colors_aa)
                    ax2.axhline(0, color="black", linewidth=0.8)
                    ax2.set_xlabel("Amino acids")
                    ax2.set_ylabel("Difference in %")
                    ax2.set_title("Sum across positions: " + label)
                    outname2 = (
                        args.rearrangements[:-4]
                        + "_aminoacids2D_"
                        + label
                        + "_squashedAA.png"
                    )
                    plt.tight_layout()
                    plt.savefig(outname2, dpi=300, bbox_inches="tight")
                    plt.close()
                    print("2D squashed-by-amino-acid plot saved as", outname2)
