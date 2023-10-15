import sys
import pandas as pd
import numpy as np
import itertools

import matplotlib.pyplot as plt
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
            sys.stderr.write("TCRcloud error: please indicate \
True or False\n")
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

        new_df = df.pivot_table(index=["v_call", "CDR3_length"],
                                aggfunc="size").reset_index()
        new_df.rename(columns={0: "counts"}, inplace=True)

        # create an empty df to serve as base
        empty_df = pd.DataFrame(columns=["v_call",
                                         "CDR3_length",
                                         "counts"])
        x_axis_names = []
        for v_genes, colour in x_axis.items():
            x_axis_names.append(v_genes)
            if ymin is None:
                limitymin = new_df.loc[new_df["CDR3_length"].idxmin()][1] - 1
            else:
                limitymin = ymin + 1
            if ymax is None:
                limitymax = new_df.loc[new_df["CDR3_length"].idxmax()][1] + 1
            else:
                limitymax = ymax - 1
            for c in range(limitymin, limitymax):
                df_new_row = pd.DataFrame({"v_call": [v_genes],
                                           "CDR3_length": [c],
                                           "counts": [0]})
                empty_df = pd.concat([empty_df, df_new_row])

        df_merged = pd.concat([new_df, empty_df],
                              ignore_index=True,
                              sort=True)

        final_df = df_merged.groupby(["v_call",
                                      "CDR3_length"]).sum().reset_index()
        final_df["frequency"] = 100 * (final_df["counts"]
                                       / final_df["counts"].sum()).round(3)
        final_df = final_df.sort_values(by=["CDR3_length", "v_call"],
                                        key=natsort_keygen())
        df_grouped = final_df.groupby(["v_call", "CDR3_length"]).sum()
        df_reset = df_grouped.reset_index()
        df_reformat = df_reset.pivot("v_call",
                                     "CDR3_length",
                                     "frequency").reset_index().rename_axis(index=None, columns=None)
        df_sorted = df_reformat.sort_values(by=["v_call"],
                                            key=natsort_keygen())
        breakpoint()
        x = df_sorted["v_call"].factorize()[0]
        for i in range(len(df_sorted.columns) - 2):
            x = np.append(x, df_sorted["v_call"].factorize()[0])
        y = np.array(df_sorted.columns.values.tolist()[1:])
        for i in range(len(df_sorted) - 1):
            y = np.append(y, np.array(df_sorted.columns.values.tolist()[1:]))
        y = np.sort(y)
        df_transpose = df_sorted.transpose()
        z = np.array(df_transpose.values.tolist()[1])
        for i in range(2, len(df_sorted.columns)):
            z = np.append(z, np.array(df_transpose.values.tolist()[i]))

        datasets.append([x, y, z, plot_aspect, x_size, ymin, ymax,
                         zmin, zmax, x_axis_ticks, x_axis_names, j[0],
                         j[1], False])
        for_comparison[j[0]].append([df_sorted, plot_aspect, x_size, ymin,
                                     ymax, zmin, zmax, x_axis_ticks,
                                     x_axis_names, j[0] + "_comparison",
                                     j[1]])

    if args.compare.lower() == "false":
        return datasets
    elif args.compare.lower() == "true":
        comparisons = []
        for m in for_comparison:
            for comb in itertools.combinations(for_comparison[m], 2):
                num1 = comb[0][0].copy()
                num1 = num1._get_numeric_data()
                num2 = comb[1][0].copy()
                num2 = num2._get_numeric_data()
                comparison1 = num1 - num2
                comparison1 = comparison1.fillna(0)
                comparison2 = num2 - num1
                comparison2 = comparison2.fillna(0)

                x = comb[0][0]["v_call"].factorize()[0]
                for i in range(len(comparison1.columns) - 1):
                    x = np.append(x, comb[0][0]["v_call"].factorize()[0])
                y = np.array(comparison1.columns.values.tolist())
                for i in range(len(comparison1) - 1):
                    y = np.append(y, np.array(
                        comparison1.columns.values.tolist()))
                y = np.sort(y)
                df_transpose = comparison1.transpose()
                z = np.array(df_transpose.values.tolist()[0])
                for i in range(1, len(comparison1.columns)):
                    z = np.append(z, np.array(df_transpose.values.tolist()[i]))

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
                    zmin = z.min() - 0.1
                if comb[0][6] is not None:
                    zmax = max(comb[0][6], comb[1][6])
                else:
                    zmax = None
                comparisons.append([x, y, z, comb[0][1], comb[0][2], ymin,
                                    ymax, zmin, zmax, comb[0][7], comb[0][8],
                                    comb[0][9], comb[0][10], True])
                x = comb[0][0]["v_call"].factorize()[0]
                for i in range(len(comparison2.columns) - 1):
                    x = np.append(x, comb[0][0]["v_call"].factorize()[0])
                y = np.array(comparison2.columns.values.tolist())
                for i in range(len(comparison2) - 1):
                    y = np.append(y, np.array(
                        comparison2.columns.values.tolist()))
                y = np.sort(y)
                df_transpose = comparison2.transpose()
                z = np.array(df_transpose.values.tolist()[0])
                for i in range(1, len(comparison2.columns)):
                    z = np.append(z, np.array(df_transpose.values.tolist()[i]))

                comparisons.append([x, y, z, comb[1][1], comb[1][2], ymin,
                                    ymax, zmin, zmax, comb[1][7], comb[1][8],
                                    comb[1][9], comb[1][10], True])
    return comparisons


def barplot(args):
    samples_df = tcrcloud.format.format_data(args)

    formatted_samples = tcrcloud.format.format_vgene(samples_df)

    samples = formatted_samples.groupby(["chain", "repertoire_id"])
    keys = [key for key, _ in samples]
    datasets = get_table(keys, samples, args)
    figure_colours = ["#847AB7",
                      "#8FCCC1",
                      "#B77A99",
                      "#E0A3AD",
                      "#70AD84",
                      "#C1C184",
                      "#B7E0F4",
                      "#EAE0AD"]
    alternate = []
    for i in datasets:
        fig = plt.figure(figsize=(10, 14))
        dataset = i[1:]
        dataset = [*dataset, dataset[0]]

        if i[13] is False:
            try:
                thecolour = figure_colours.pop(0)
            except IndexError:
                thecolour = "#BBBBBB"
            ax = fig.add_subplot(projection="3d")
            x = i[0]
            y = i[1]
            z = i[2]
            plot_aspect = i[3]
            x_size = i[4]
            ymin = i[5]
            ymax = i[6]
            zmin = i[7]
            zmax = i[8]
            x_axis_ticks = i[9]
            x_axis_names = i[10]
            ax.set_box_aspect(aspect=plot_aspect)
            ax.bar3d(x, y, z * 0, 1, 1, z, shade=True, color=thecolour)
            ax.set_ylim(ymin, ymax)
            ax.set_zlim(zmin, zmax)
            ax.set_xticks(x_axis_ticks)
            ax.set_xticklabels(x_axis_names,
                               rotation=45,
                               ha="right",
                               fontsize=x_size)
            ax.yaxis.get_major_locator().set_params(integer=True)
            ax.zaxis.get_major_locator().set_params(prune="lower")

        if i[13] is True:
            try:
                thecolour = figure_colours.pop(0)
                alternate.append(thecolour)
            except IndexError:
                thecolour = "#BBBBBB"
            ax = fig.add_subplot(projection="3d")
            x = i[0]
            y = i[1]
            z = i[2]
            plot_aspect = (i[3][0], 2, i[3][2])
            x_size = i[4]
            ymin = i[5]
            ymax = i[6]
            zmin = i[7]
            zmax = i[8]
            x_axis_ticks = i[9]
            x_axis_names = i[10]
            ax.set_box_aspect(aspect=plot_aspect)
            z1 = z.copy()
            z1[z1 > 0] = 0
            if len(alternate) == 1:
                ax.bar3d(x, y, z1 * 0, 1, 1, z1, shade=True,
                         color=figure_colours[0])
            else:
                ax.bar3d(x, y, z1 * 0, 1, 1, z1,
                         shade=True, color=alternate[-2])
            z[z < 0] = 0
            ax.bar3d(x, y, z * 0, 1, 1, z, shade=True, color=thecolour)
            ax.set_ylim(ymin, ymax)
            ax.set_zlim(zmin, zmax)
            ax.set_xticks(x_axis_ticks)
            ax.set_xticklabels(x_axis_names,
                               rotation=60,
                               ha="right",
                               fontsize=x_size)
            ax.yaxis.get_major_locator().set_params(integer=True)
            ax.zaxis.get_major_locator().set_params(prune="lower")
            ax.elev = 5

        outputname = args.rearrangements[:-4] + \
            "_vgenes_" + i[12] + "_" + i[11] + ".png"
        plt.savefig(outputname, dpi=300, bbox_inches="tight")
        print("V genes plot saved as " + outputname)
