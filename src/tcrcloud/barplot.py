import pandas as pd

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
    datasets = []
    for j in keys:
        if j[0] == "A":
            x_axis = TRAV
            plot_aspect = (3.5, 1, 1)
            x_size = 6
            ymax = args.yhighalpha
            ymin = args.ylowalpha
            zmax = args.zhighalpha
        if j[0] == "B":
            x_axis = TRBV
            plot_aspect = (3, 1, 1)
            x_size = 4
            ymax = args.yhighbeta
            ymin = args.ylowbeta
            zmax = args.zhighbeta
        if j[0] == "G":
            x_axis = TRGV
            plot_aspect = (2, 1, 1)
            x_size = 10
            ymax = args.yhighgamma
            ymin = args.ylowgamma
            zmax = args.zhighgamma
        if j[0] == "D":
            x_axis = TRDV
            plot_aspect = (1.5, 1, 1)
            x_size = 10
            ymax = args.yhighdelta
            ymin = args.ylowdelta
            zmax = args.zhighdelta
        if j[0] == "H":
            x_axis = IGHV
            plot_aspect = (5, 1, 1)
            x_size = 2
            ymax = args.yhighheavy
            ymin = args.ylowheavy
            zmax = args.zhighheavy
        if j[0] == "K":
            x_axis = IGKV
            plot_aspect = (3.5, 1, 1)
            x_size = 2
            ymax = args.yhighkappa
            ymin = args.ylowhkappa
            zmax = args.zhighkappa
        if j[0] == "L":
            x_axis = IGLV
            plot_aspect = (3.5, 1, 1)
            x_size = 4
            ymax = args.yhighlambda
            ymin = args.ylowhlambda
            zmax = args.zhighlambda

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

        x = final_df["v_call"].factorize()[0]
        y = final_df["CDR3_length"].to_numpy()
        z = final_df["frequency"].to_numpy()
        datasets.append([x, y, z, plot_aspect, x_size, ymin, ymax,
                         zmax, x_axis_ticks, x_axis_names, j[0], j[1]])
    return datasets


def barplot(args):

    samples_df = tcrcloud.format.format_data(args)

    formatted_samples = tcrcloud.format.format_barplot(samples_df)

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

    for i in datasets:
        fig = plt.figure(figsize=(10, 14))
        dataset = i[1:]
        dataset = [*dataset, dataset[0]]

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
        zmax = i[7]
        x_axis_ticks = i[8]
        x_axis_names = i[9]
        ax.set_box_aspect(aspect=plot_aspect)
        ax.bar3d(x, y, z * 0, 1, 1, z, shade=True, color=thecolour)
        ax.set_ylim(ymin, ymax)
        ax.set_zlim(0, zmax)
        ax.set_xticks(x_axis_ticks)
        ax.set_xticklabels(x_axis_names,
                           rotation=45,
                           ha="right",
                           fontsize=x_size)
        ax.yaxis.get_major_locator().set_params(integer=True)
        ax.zaxis.get_major_locator().set_params(prune="lower")

        outputname = args.rearrangements[:-4] + "_barplot_" + i[11] + "_" + i[10] + ".png"
        plt.savefig(outputname, dpi=300, bbox_inches="tight")
        print("Barplot saved as " + outputname)
