import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import tcrcloud.format
import tcrcloud.colours

colours = tcrcloud.colours.aminoacids


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
        for i in list(new_df.columns)[5:]:
            positions.append(new_df.groupby(i).sum())
            names.append(i)
        final = pd.concat(positions, axis=1)
        final = final.sort_index()
        final.drop(final.head(1).index, inplace=True)
        final.set_axis(names, axis="columns", inplace=True)
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
        aa_rank = [4, 1, 15, 16, 9, 2, 13, 6, 14,
                   7, 8, 19, 3, 20, 12, 17, 18, 5, 11, 10]
        normalized["rank"] = aa_rank
        normalized.sort_values(by=["rank"], inplace=True)
        normalized.drop("rank", inplace=True, axis=1)
        transposed_df = normalized.transpose()

        if args.threeD.lower() == "false":
            transposed_df.plot(kind="bar", stacked=True,
                               color=colours, figsize=(10, 14))
            outputname = args.rearrangements[:-4] + \
                "_aminoacids_" + j[1] + "_" + j[0] + ".png"
            plt.legend(bbox_to_anchor=(1.01, 1),
                       reverse=True, loc="upper left")
            plt.savefig(outputname, dpi=300, bbox_inches="tight")
            print("Amino acids plot saved as " + outputname)

        elif args.threeD.lower() == "true":
            if len(normalized) < 20:
                all_aa = pd.DataFrame(np.zeros((20, 1)),
                                      columns=["just_empty"],
                                      index=["C", "G", "P", "A", "V", "I", "L",
                                             "M", "F", "Y", "W", "R", "H", "K",
                                             "D", "E", "S", "T", "N", "Q"])
                result = pd.concat([normalized, all_aa], axis=1)
                normalized = result.drop("just_empty", axis=1)
                normalized = normalized.fillna(0)

            x = normalized.index.factorize()[0]
            for i in range(len(normalized.columns) - 1):
                x = np.append(x, normalized.index.factorize()[0])
            x = np.sort(x)
            y = np.array(normalized.columns.values.tolist())
            for i in range(len(normalized) - 1):
                y = np.append(y, np.array(
                    normalized.columns.values.tolist()))
            # y = np.sort(y)
            # df_transpose = normalized.transpose()
            z = np.array(normalized.values.tolist()[0])
            for i in range(1, len(normalized)):
                z = np.append(z, np.array(normalized.values.tolist()[i]))
            fig = plt.figure(figsize=(10, 14))
            ax = fig.add_subplot(projection="3d")
            counter = 0
            for i in colours:
                ax.bar3d(x[int(len(z)/20)*counter:int(len(z)/20)*(counter+1)], y[int(len(z)/20)*counter:int(len(z)/20)*(counter+1)], z[int(
                    len(z)/20)*counter:int(len(z)/20)*(counter+1)] * 0, 0.7, 0.7, z[int(len(z)/20)*counter:int(len(z)/20)*(counter+1)], shade=True, color=colours[i], alpha=0.5)
                counter += 1

            ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                          11, 12, 13, 14, 15, 16, 17, 18, 19])
            ax.set_xticklabels([key for key in colours], ha="right")
            ax.set_ylim(1.5, len(z)/20 + 2)
            ax.yaxis.get_major_locator().set_params(integer=True)
            ax.elev = 35
            ax.azim = 40
            outputname = args.rearrangements[:-4] + \
                "_aminoacids3D_" + j[1] + "_" + j[0] + ".png"
            plt.savefig(outputname, dpi=300, bbox_inches="tight")
            print("Tridimensional Amino acids plot saved as " + outputname)
