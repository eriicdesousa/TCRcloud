import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import plotly.graph_objects as go

import tcrcloud.format
import tcrcloud.colours

colours = tcrcloud.colours.aminoacids


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
    samples_df = tcrcloud.format.format_data(args)
    formatted_samples = tcrcloud.format.format_aminoacids(samples_df)
    samples = formatted_samples.groupby(["chain", "repertoire_id"])
    keys = [key for key, _ in samples]
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
        final = final.sort_index()
        final.drop(final.head(1).index, inplace=True)
        final.columns = names
        final.drop(columns=final.columns[-1:], axis=1, inplace=True)
        if args.length is not None:
            if len(final.columns) > args.length:
                diff = len(final.columns) - args.length
                final.drop(columns=final.columns[-diff:], axis=1, inplace=True)
            if len(final.columns) < args.length:
                diff = args.length - len(final.columns)
                for number in range(1, diff + 1):
                    final[len(final.columns) + 1] = np.nan
        final = final.replace(np.nan, 0)
        normalized = final.apply(lambda x: x * 100 / sum(x), axis=0)
        normalized = normalized.replace(np.nan, 0)
        aa_rank = [
            16,
            9,
            2,
            1,
            20,
            17,
            3,
            12,
            5,
            14,
            13,
            7,
            8,
            6,
            4,
            11,
            10,
            15,
            18,
            19,
        ]
        normalized["rank"] = aa_rank
        normalized.sort_values(by=["rank"], inplace=True)
        normalized.drop("rank", inplace=True, axis=1)
        transposed_df = normalized.transpose()

        if args.export.lower() == "true":
            df_filename = args.rearrangements[:-4] + "_aminoacids_table.csv"
            normalized.to_csv(df_filename, index=True)

        if args.threeD.lower() == "false":
            transposed_df.plot(
                kind="bar", stacked=True, color=colours, figsize=(10, 14)
            )
            outputname = (
                args.rearrangements[:-4] + "_aminoacids_" + j[1] + "_" + j[0] + ".png"
            )
            plt.legend(bbox_to_anchor=(1.01, 1), reverse=True, loc="upper left")
            plt.savefig(outputname, dpi=300, bbox_inches="tight")
            print("Amino acids plot saved as " + outputname)

        elif args.threeD.lower() == "true":
            if len(normalized) < 20:
                all_aa = pd.DataFrame(
                    np.zeros((20, 1)),
                    columns=["just_empty"],
                    index=[
                        "F",
                        "Y",
                        "W",
                        "G",
                        "A",
                        "V",
                        "L",
                        "M",
                        "I",
                        "S",
                        "T",
                        "C",
                        "P",
                        "N",
                        "Q",
                        "K",
                        "R",
                        "H",
                        "D",
                        "E",
                    ],
                )
                result = pd.concat([normalized, all_aa], axis=1)
                normalized = result.drop("just_empty", axis=1)
                normalized = normalized.fillna(0)
            normalized = normalized.iloc[::-1]
            x = normalized.index.factorize()[0]
            # for i in range(len(normalized.columns) - 1):
            #     x = np.append(x, normalized.index.factorize()[0])
            # x = np.sort(x)
            y = np.array(normalized.columns.values.tolist())
            # for i in range(len(normalized) - 1):
            #     y = np.append(y, np.array(normalized.columns.values.tolist()))
            # y = np.sort(y)
            # df_transpose = normalized.transpose()
            z = np.array(normalized.values.tolist()[0])
            for i in range(1, len(normalized)):
                z = np.append(z, np.array(normalized.values.tolist()[i]))
            # fig = plt.figure(figsize=(10, 14))
            # ax = fig.add_subplot(projection="3d")
            # counter = 0
            # for i in colours:
            #     ax.bar3d(
            #         x[int(len(z) / 20) * counter : int(len(z) / 20) * (counter + 1)],
            #         y[int(len(z) / 20) * counter : int(len(z) / 20) * (counter + 1)],
            #         z[int(len(z) / 20) * counter : int(len(z) / 20) * (counter + 1)]
            #         * 0,
            #         0.7,
            #         0.7,
            #         z[int(len(z) / 20) * counter : int(len(z) / 20) * (counter + 1)],
            #         shade=True,
            #         color=colours[i],
            #         alpha=0.5,
            #     )
            #     counter += 1

            # ax.set_xticks(
            #     [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            # )
            # ax.set_xticklabels([key for key in colours], ha="right")
            # ax.set_ylim(1.5, len(z) / 20 + 2)
            # ax.yaxis.get_major_locator().set_params(integer=True)
            # ax.elev = 35
            # ax.azim = 40

            x_df = pd.Series(y)
            y_df = pd.Series(x)
            z_df = pd.Series(z)

            x_min = 0
            y_min = 0
            z_min = 0.8 * min(z_df)
            step = 1

            mesh_list = []
            colors = list(colours.values())
            color_value = 0

            x_df_uniq = x_df.unique()
            y_df_uniq = y_df.unique()
            len_x_df_uniq = len(x_df_uniq)

            for idx, x_data in enumerate(x_df_uniq):
                for idx2, y_data in enumerate(y_df_uniq):
                    color_value = colors[idx2 % 22]
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
                        ),
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
                    ticktext=[key for key in colours],
                    tickvals=[v for v in range(1, 40 + 1) if v % 2 != 0],
                    tickfont=dict(size=7),
                    titlefont=dict(size=8),
                ),
                yaxis=dict(
                    tickmode="array",
                    ticktext=[str(v) for v in range(1, y[-1] + 1)],
                    tickvals=[v for v in range(1, y[-1] * 2 + 1) if v % 2 != 0],
                    tickfont=dict(size=7),
                    titlefont=dict(size=8),
                ),
                zaxis=dict(
                    tickfont=dict(size=8),
                    titlefont=dict(size=8),
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
                args.rearrangements[:-4] + "_aminoacids3D_" + j[1] + "_" + j[0] + ".png"
            )
            # plt.savefig(outputname, dpi=300, bbox_inches="tight")
            fig.write_image(outputname, scale=6)
            print("Tridimensional Amino acids plot saved as " + outputname)
