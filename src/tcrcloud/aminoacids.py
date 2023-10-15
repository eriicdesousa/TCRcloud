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
        splitted = df['junction_aa'].str.split('', expand=True)
        new_df = pd.concat([df, splitted], axis=1)
        positions = []
        names = []
        for i in list(new_df.columns)[5:]:
            positions.append(new_df.groupby(i).sum())
            names.append(i)
        final = pd.concat(positions, axis=1)
        final = final.sort_index()
        final.drop(final.head(1).index, inplace=True)
        final.set_axis(names, axis='columns', inplace=True)
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
        aa_rank = [4, 1, 15, 16, 9, 2, 13, 6, 14,
                   7, 8, 19, 3, 20, 12, 17, 18, 5, 11, 10]
        normalized["rank"] = aa_rank
        normalized.sort_values(by=['rank'], inplace=True)
        normalized.drop("rank", inplace=True, axis=1)
        transposed_df = normalized.transpose()
        # breakpoint()
        # x = df_sorted["v_call"].factorize()[0]
        # for i in range(len(df_sorted.columns) - 2):
        #     x = np.append(x, df_sorted["v_call"].factorize()[0])
        # y = np.array(df_sorted.columns.values.tolist()[1:])
        # for i in range(len(df_sorted) - 1):
        #     y = np.append(y, np.array(df_sorted.columns.values.tolist()[1:]))
        # y = np.sort(y)
        # df_transpose = df_sorted.transpose()
        # z = np.array(df_transpose.values.tolist()[1])
        # for i in range(2, len(df_sorted.columns)):
        #     z = np.append(z, np.array(df_transpose.values.tolist()[i]))

        transposed_df.plot(kind="bar", stacked=True,
                           color=colours, figsize=(10, 14))
        outputname = args.rearrangements[:-4] + \
            "_aminoacids_" + j[1] + "_" + j[0] + ".png"
        plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
        plt.savefig(outputname, dpi=300, bbox_inches="tight")
        print("Amino acids plot saved as " + outputname)
