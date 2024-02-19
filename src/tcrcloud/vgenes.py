import sys
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import plotly.graph_objects as go
from natsort import natsort_keygen

import tcrcloud.format
import tcrcloud.colours

# Import default dicts that contain V gene information
TRAV = tcrcloud.colours.TRAV
TRBV = tcrcloud.colours.TRBV
TRGV = tcrcloud.colours.TRGV
TRDV = tcrcloud.colours.TRDV
IGHV = tcrcloud.colours.IGHV
IGKV = tcrcloud.colours.IGKV
IGLV = tcrcloud.colours.IGLV


def get_table(keys, samples, args):
    if args.compare.lower() != "true":
        if args.compare.lower() != "false":
            sys.stderr.write(
                "TCRcloud error: please indicate \
True or False\n"
            )
            exit()

    datasets = []
    for_comparison = {}
    for_comparison["A"] = []
    for_comparison["B"] = []
    for_comparison["G"] = []
    for_comparison["D"] = []
    for_comparison["H"] = []
    for_comparison["K"] = []
    for_comparison["L"] = []
    for j in keys:
        if j[0] == "A":
            x_axis = TRAV
            plot_aspect = (3.5, 1, 1)
            x_size = 6
            ymax = args.yhighalpha
            ymin = args.ylowalpha
            zmax = args.zhighalpha
            zmin = args.zlowalpha
        if j[0] == "B":
            x_axis = TRBV
            plot_aspect = (3, 1, 1)
            x_size = 4
            ymax = args.yhighbeta
            ymin = args.ylowbeta
            zmax = args.zhighbeta
            zmin = args.zlowbeta
        if j[0] == "G":
            x_axis = TRGV
            plot_aspect = (2, 1, 1)
            x_size = 10
            ymax = args.yhighgamma
            ymin = args.ylowgamma
            zmax = args.zhighgamma
            zmin = args.zlowgamma
        if j[0] == "D":
            x_axis = TRDV
            plot_aspect = (1.5, 1, 1)
            x_size = 10
            ymax = args.yhighdelta
            ymin = args.ylowdelta
            zmax = args.zhighdelta
            zmin = args.zlowdelta
        if j[0] == "H":
            x_axis = IGHV
            plot_aspect = (5, 1, 1)
            x_size = 2
            ymax = args.yhighheavy
            ymin = args.ylowheavy
            zmax = args.zhighheavy
            zmin = args.zlowheavy
        if j[0] == "K":
            x_axis = IGKV
            plot_aspect = (3.5, 1, 1)
            x_size = 2
            ymax = args.yhighkappa
            ymin = args.ylowhkappa
            zmax = args.zhighkappa
            zmin = args.zlowkappa
        if j[0] == "L":
            x_axis = IGLV
            plot_aspect = (3.5, 1, 1)
            x_size = 4
            ymax = args.yhighlambda
            ymin = args.ylowhlambda
            zmax = args.zhighlambda
            zmin = args.zlowlambda

        x_axis_ticks = []
        for i in range(0, len(x_axis)):
            x_axis_ticks.append(i)

        df = samples.get_group(j)

        new_df = df.pivot_table(
            index=["v_call", "CDR3_length"], aggfunc="size"
        ).reset_index()
        new_df.rename(columns={0: "counts"}, inplace=True)

        # create an empty df to serve as base
        empty_df = pd.DataFrame(columns=["v_call", "CDR3_length", "counts"])
        x_axis_names = []
        for v_genes, colour in x_axis.items():
            x_axis_names.append(v_genes)
            if ymin is None:
                limitymin = new_df.loc[new_df["CDR3_length"].idxmin()].iloc[1] - 1
            else:
                limitymin = ymin + 1
            if ymax is None:
                limitymax = new_df.loc[new_df["CDR3_length"].idxmax()].iloc[1] + 1
            else:
                limitymax = ymax - 1
            for c in range(limitymin, limitymax):
                df_new_row = pd.DataFrame(
                    {"v_call": [v_genes], "CDR3_length": [c], "counts": [0]}
                )
                empty_df = pd.concat([empty_df, df_new_row])

        df_merged = pd.concat([new_df, empty_df], ignore_index=True, sort=True)

        final_df = df_merged.groupby(["v_call", "CDR3_length"]).sum().reset_index()
        final_df["frequency"] = 100 * (
            final_df["counts"] / final_df["counts"].sum()
        )  # .round(3)
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

        if args.export.lower() == "true":
            df_filename = (
                args.rearrangements[:-4] + "_vgenes_table" + j[1] + "_" + j[0] + ".csv"
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
            [
                x,
                y,
                z,
                plot_aspect,
                x_size,
                ymin,
                ymax,
                zmin,
                zmax,
                x_axis_ticks,
                x_axis_names,
                j[0],
                j[1],
                False,
            ]
        )
        for_comparison[j[0]].append(
            [
                df_sorted,
                plot_aspect,
                x_size,
                ymin,
                ymax,
                zmin,
                zmax,
                x_axis_ticks,
                x_axis_names,
                j[0] + "_comparison",
                j[1],
            ]
        )

    if args.compare.lower() == "false":
        return datasets
    elif args.compare.lower() == "true":
        comparisons = []
        for m in for_comparison:
            if len(for_comparison[m]) < 2:
                sys.stderr.write(
                    "Less than 2 repertoires from the "
                    + m
                    + " chain were detected in the rearragements file\n"
                )
            if len(for_comparison[m]) > 2:
                sys.stderr.write(
                    "More than 2 repertoires from the "
                    + m
                    + " chain were detected in the rearragements file\n"
                )
            if len(for_comparison[m]) == 2:
                comb = for_comparison[m]
                num1 = comb[0][0].copy()
                num1 = num1.drop("v_call", axis=1)
                num2 = comb[1][0].copy()
                num2 = num2.drop("v_call", axis=1)
                comparison1 = num1 - num2
                comparison1 = comparison1.fillna(0)
                comparison2 = num2 - num1
                comparison2 = comparison2.fillna(0)
                comparison1.insert(0, "v_call", comb[0][0]["v_call"])
                comparison2.insert(0, "v_call", comb[0][0]["v_call"])

                if comb[0][3] is not None:
                    ymin = min(comb[0][3], comb[1][3])
                else:
                    ymin = None
                if comb[0][4] is not None:
                    ymax = max(comb[0][4], comb[1][4])
                else:
                    ymax = None
                if comb[0][5] is not None:
                    zmin = max(comb[0][5], comb[1][5])
                else:
                    zmin = None
                if comb[0][6] is not None:
                    zmax = max(comb[0][6], comb[1][6])
                else:
                    zmax = None

                x = comparison1["v_call"].factorize()[0]
                y = np.array(comparison1.columns.values.tolist()[1:])
                df_transpose = comparison1.transpose()
                z = np.array(df_transpose.values.tolist()[1])
                for i in range(2, len(comparison1.columns)):
                    z = np.append(z, np.array(df_transpose.values.tolist()[i]))
                z = np.array_split(z, len(comparison1.columns) - 1)

                comparisons.append(
                    [
                        x,
                        y,
                        z,
                        comb[0][1],
                        comb[0][2],
                        ymin,
                        ymax,
                        zmin,
                        zmax,
                        comb[0][7],
                        comb[0][8],
                        comb[0][9],
                        comb[0][10],
                        True,
                    ]
                )

                x = comparison2["v_call"].factorize()[0]
                y = np.array(comparison2.columns.values.tolist()[1:])
                df_transpose = comparison2.transpose()
                z = np.array(df_transpose.values.tolist()[1])
                for i in range(2, len(comparison2.columns)):
                    z = np.append(z, np.array(df_transpose.values.tolist()[i]))
                z = np.array_split(z, len(comparison2.columns) - 1)

                comparisons.append(
                    [
                        x,
                        y,
                        z,
                        comb[1][1],
                        comb[1][2],
                        ymin,
                        ymax,
                        zmin,
                        zmax,
                        comb[1][7],
                        comb[1][8],
                        comb[1][9],
                        comb[1][10],
                        True,
                    ]
                )
    return comparisons


def barplot(args):
    samples_df = tcrcloud.format.format_data(args)

    formatted_samples = tcrcloud.format.format_vgene(samples_df)

    samples = formatted_samples.groupby(["chain", "repertoire_id"])
    keys = [key for key, _ in samples]
    datasets = get_table(keys, samples, args)

    for i in datasets:
        fig = plt.figure(figsize=(10, 14))
        dataset = i[1:]
        dataset = [*dataset, dataset[0]]

        if i[13] is False:
            fig = go.Figure(
                go.Surface(
                    x=i[0], y=i[1], z=i[2], colorscale="Turbo", cmin=i[7], cmax=i[8]
                )
            )
            camera = dict(eye=dict(x=2.5, y=-3.5, z=2.5))

            sc = dict(
                aspectratio=dict(x=i[3][0], y=i[3][1], z=i[3][2]),
                xaxis_title=i[10][0][:4],
                yaxis_title="CDR3 Length",
                zaxis_title="Percentage of reads",
                xaxis=dict(
                    tickmode="array",
                    ticktext=i[10],
                    tickvals=i[9],
                    tickfont=dict(size=i[4]),
                    titlefont=dict(size=10),
                ),
                yaxis=dict(
                    tickfont=dict(size=8), titlefont=dict(size=10), range=[i[5], i[6]]
                ),
                zaxis=dict(
                    tickfont=dict(size=8), titlefont=dict(size=10), range=[i[7], i[8]]
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
                args.rearrangements[:-4] + "_vgenes_" + i[12] + "_" + i[11] + ".png"
            )
            fig.write_image(outputname, scale=6)
            print("V genes plot saved as " + outputname)

        if i[13] is True:
            i = datasets[0]
            fig = go.Figure(
                go.Surface(
                    x=i[0], y=i[1], z=i[2], colorscale="Portland", cmin=i[7], cmax=i[8]
                )
            )
            camera = dict(eye=dict(x=2.5, y=-5, z=0.5))

            sc = dict(
                aspectratio=dict(x=i[3][0], y=i[3][1], z=i[3][2]),
                xaxis_title=i[10][0][:4],
                yaxis_title="CDR3 Length",
                zaxis_title="Percentage of reads",
                xaxis=dict(
                    tickmode="array",
                    ticktext=i[10],
                    tickvals=i[9],
                    tickfont=dict(size=i[4]),
                    tickangle=-45,
                    titlefont=dict(size=10),
                ),
                yaxis=dict(
                    tickfont=dict(size=6), titlefont=dict(size=10), range=[i[5], i[6]]
                ),
                zaxis=dict(
                    tickfont=dict(size=8), titlefont=dict(size=10), range=[i[7], i[8]]
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
                args.rearrangements[:-4] + "_vgenes_" + i[12] + "_" + i[11] + ".png"
            )
            fig.write_image(outputname, scale=6)
            print("V genes plot saved as " + outputname)

            i = datasets[1]
            fig = go.Figure(
                go.Surface(
                    x=i[0],
                    y=i[1],
                    z=i[2],
                    colorscale="Portland_r",
                    cmin=i[7],
                    cmax=i[8],
                )
            )
            camera = dict(eye=dict(x=2.5, y=-5, z=0.5))

            sc = dict(
                aspectratio=dict(x=i[3][0], y=i[3][1], z=i[3][2]),
                xaxis_title=i[10][0][:4],
                yaxis_title="CDR3 Length",
                zaxis_title="Percentage of reads",
                xaxis=dict(
                    tickmode="array",
                    ticktext=i[10],
                    tickvals=i[9],
                    tickfont=dict(size=i[4]),
                    tickangle=-45,
                    titlefont=dict(size=10),
                ),
                yaxis=dict(
                    tickfont=dict(size=6), titlefont=dict(size=10), range=[i[5], i[6]]
                ),
                zaxis=dict(
                    tickfont=dict(size=8), titlefont=dict(size=10), range=[i[7], i[8]]
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
                args.rearrangements[:-4] + "_vgenes_" + i[12] + "_" + i[11] + ".png"
            )
            fig.write_image(outputname, scale=6)
            print("V genes plot saved as " + outputname)
            break
