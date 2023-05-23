import pandas as pd

import matplotlib.pyplot as plt

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


def surface(args):

    samples_df = tcrcloud.format.format_data(args)

    formatted_samples = tcrcloud.format.format_surface(samples_df)

    samples = formatted_samples.groupby(["chain", "repertoire_id"])
    keys = [key for key, _ in samples]

    for j in keys:
        if j[0] == "A":
            x_axis = TRAV
            plot_aspect = (3.5, 1, 1)
            x_size = 6
        if j[0] == "B":
            x_axis = TRBV
            plot_aspect = (3, 1, 1)
            x_size = 4
        if j[0] == "G":
            x_axis = TRGV
            plot_aspect = (2, 1, 1)
            x_size = 10
        if j[0] == "D":
            x_axis = TRDV
            plot_aspect = (1.5, 1, 1)
            x_size = 10
        if j[0] == "H":
            x_axis = IGHV
            plot_aspect = (5, 1, 1)
            x_size = 2
        if j[0] == "K":
            x_axis = IGKV
            plot_aspect = (3.5, 1, 1)
            x_size = 2
        if j[0] == "L":
            x_axis = IGLV
            plot_aspect = (3.5, 1, 1)
            x_size = 4

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
            if args.ymin is None:
                limitymin = new_df.loc[new_df["CDR3_length"].idxmin()][1] - 1
            else:
                limitymin = args.ymin + 1
            if args.ymax is None:
                limitymax = new_df.loc[new_df["CDR3_length"].idxmax()][1] + 1
            else:
                limitymax = args.ymax - 1
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
        final_df = final_df.sort_values(by=["CDR3_length", "v_call"])

        x = final_df["v_call"].factorize()[0]
        y = final_df["CDR3_length"].to_numpy()
        z = final_df["frequency"].to_numpy()
        fig = plt.figure(figsize=(10, 14))
        ax = fig.add_subplot(projection="3d")
        ax.set_box_aspect(aspect=plot_aspect)
        ax.plot_trisurf(x, y, z, cmap="viridis")
        ax.set_ylim(args.ymin, args.ymax)
        ax.set_zlim(0, args.zmax)
        ax.set_xticks(x_axis_ticks)
        ax.set_xticklabels(x_axis_names,
                           rotation=45,
                           ha="right",
                           fontsize=x_size)
        ax.yaxis.get_major_locator().set_params(integer=True)
        ax.zaxis.get_major_locator().set_params(prune="lower")
        outputname = args.rearrangements[:-4] + "_" + j[1] + "_surface_" + j[0] + ".png"
        plt.savefig(outputname, dpi=300, bbox_inches="tight")
        print("Surface plot saved as " + outputname)
