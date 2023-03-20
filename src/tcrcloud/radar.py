from sklearn import preprocessing
import json
import sys

import numpy as np
import matplotlib.pyplot as plt
import skbio

import tcrcloud.format


def calculate_convergence(df):
    counts = df["counts"]
    total = counts.sum()
    conv_value = 0
    for i in df["counts"]:
        if i > 1:
            conv_value += i
    return conv_value / total


def calculate_dfifty(df, length):
    counts = df["counts"]
    counter = 0
    if length < 10000:
        total = counts.sum()
        for i in counts.cumsum():
            counter += 1
            if i > total / 2:
                break
        return (counter * 100) / length
    elif length >= 10000:
        counts = counts.head(10000)
        total = counts.sum()
        for i in counts.cumsum():
            counter += 1
            if i > total / 2:
                break
        return (counter * 100) / 10000


def calculate_metrics(keys, samples, legend_file, export, filename):
    minmax_scale = preprocessing.MinMaxScaler(feature_range=(0, 1))

    patients = []
    list_for_printing = []

    for j in keys:
        df = samples.get_group(j)
        convergence = np.array([0, calculate_convergence(tcrcloud.format.
                                                         format_convergence(df)
                                                         ), 1]).reshape(-1, 1)
        df = tcrcloud.format.format_metrics(df)
        length = len(df)
        distinct = np.array([0, length, 55000]).reshape(-1, 1)
        counts = df["counts"].tolist()
        dfifty = np.array([0, calculate_dfifty(df, length), 50]).reshape(-1, 1)
        del df
        shannon = np.array([0, skbio.diversity.
                            alpha_diversity("shannon", counts)[0], 15]
                           ).reshape(-1, 1)
        simpson = np.array([0.75, skbio.diversity.
                            alpha_diversity("simpson", counts)[0], 1]
                           ).reshape(-1, 1)
        chao = np.array([0, skbio.diversity.
                         alpha_diversity("chao1", counts)[0], 90000]
                        ).reshape(-1, 1)
        gini = np.array([0, skbio.diversity.
                         alpha_diversity("gini_index", counts)[0], 1]
                        ).reshape(-1, 1)
        metrics = []

        if legend_file is not None:
            try:
                with open(legend_file) as json_file:
                    legend_dict = json.load(json_file)
            except FileNotFoundError:
                sys.stderr.write("TCRcloud error: " + legend_file
                                 + " doesn't seem to exist\n")
                exit()
            except json.decoder.JSONDecodeError:
                sys.stderr.write("TCRcloud error: " + legend_file
                                 + " doesn't seem properly formatted. Check \
https://github.com/oldguyeric/TCRcloud for more information\n")
                exit()
        else:
            legend_dict = {}

        if j[0] == "A":
            metrics.append(legend_dict.get(j[1], j[1]) + " α chain")
        elif j[0] == "B":
            metrics.append(legend_dict.get(j[1], j[1]) + " β chain")
        elif j[0] == "G":
            metrics.append(legend_dict.get(j[1], j[1]) + " γ chain")
        elif j[0] == "D":
            metrics.append(legend_dict.get(j[1], j[1]) + " δ chain")
        elif j[0] == "H":
            metrics.append(legend_dict.get(j[1], j[1]) + " Heavy chain")
        elif j[0] == "K":
            metrics.append(legend_dict.get(j[1], j[1]) + " Kappa chain")
        elif j[0] == "L":
            metrics.append(legend_dict.get(j[1], j[1]) + " Light chain")

        if export.lower() == "true":
            metrics_for_printing = metrics.copy()
            metrics_for_printing.append(float(np.format_float_positional(
                                        dfifty[1], precision=3)))
            metrics_for_printing.append(float(np.format_float_positional(
                                        distinct[1], precision=3)))
            metrics_for_printing.append(float(np.format_float_positional(
                                        chao[1], precision=3)))
            metrics_for_printing.append(float(np.format_float_positional(
                                        simpson[1], precision=3)))
            metrics_for_printing.append(float(np.format_float_positional(
                                        shannon[1], precision=3)))
            metrics_for_printing.append(float(np.format_float_positional(
                                        gini[1], precision=3)))
            metrics_for_printing.append(float(np.format_float_positional(
                                        convergence[1], precision=3)))
            list_for_printing.append(metrics_for_printing)

        metrics.append(minmax_scale.fit_transform(dfifty)[1].
                       astype(float))
        metrics.append(minmax_scale.fit_transform(distinct)[1].
                       astype(float))
        metrics.append(minmax_scale.fit_transform(chao)[1].
                       astype(float))
        metrics.append(minmax_scale.fit_transform(simpson)[1].
                       astype(float))
        metrics.append(minmax_scale.fit_transform(shannon)[1].
                       astype(float))
        metrics.append(minmax_scale.fit_transform(gini)[1].
                       astype(float))
        metrics.append(minmax_scale.fit_transform(convergence)[1].
                       astype(float))
        patients.append(metrics)
    # patients = []
    # patients.append(["test1", np.array(0.1), np.array(0.9), np.array(0.1),
    #                  np.array(0.9), np.array(0.1), np.array(0.9),
    #                  np.array(0.1)])
    # patients.append(["test2", np.array(0.2), np.array(0.8), np.array(0.2),
    #                  np.array(0.8), np.array(0.2), np.array(0.8),
    #                  np.array(0.2)])
    # patients.append(["test3", np.array(0.3), np.array(0.7), np.array(0.3),
    #                  np.array(0.7), np.array(0.3), np.array(0.7),
    #                  np.array(0.3)])
    # patients.append(["test4", np.array(0.4), np.array(0.6), np.array(0.4),
    #                  np.array(0.6), np.array(0.4), np.array(0.6),
    #                  np.array(0.4)])
    # patients.append(["test5", np.array(0.5), np.array(0.5), np.array(0.5),
    #                  np.array(0.5), np.array(0.5), np.array(0.5),
    #                  np.array(0.5)])
    # patients.append(["test6", np.array(0.6), np.array(0.4), np.array(0.6),
    #                  np.array(0.4), np.array(0.6), np.array(0.4),
    #                  np.array(0.6)])
    # patients.append(["test7", np.array(0.7), np.array(0.3), np.array(0.7),
    #                  np.array(0.3), np.array(0.7), np.array(0.3),
    #                  np.array(0.7)])
    # patients.append(["test8", np.array(0.8), np.array(0.2), np.array(0.8),
    #                  np.array(0.2), np.array(0.8), np.array(0.2),
    #                  np.array(0.8)])
    # patients.append(["test9", np.array(0.9), np.array(0.1), np.array(0.9),
    #                  np.array(0.1), np.array(0.9), np.array(0.1),
    #                  np.array(0.9)])
    if export.lower() == "true":
        metrics_filename = filename[:-4] + "_repertoire_metrics.txt"
        with open(metrics_filename, "w") as fileout:
            for i in list_for_printing:
                print("Repertoire:", i[0], file=fileout)
                print("D50 Index:", i[1], file=fileout)
                print("Distinct CDR3:", int(i[2]), file=fileout)
                print("Chao1 Index:", i[3], file=fileout)
                print("Simpson Index:", i[4], file=fileout)
                print("Shannon Index:", i[5], file=fileout)
                print("Gini Index:", i[6], file=fileout)
                print("Convergence:", i[7], file=fileout)
                print(file=fileout)
        print("Repertoire metrics saved as " + metrics_filename)
    return patients


def radar(args):
    if args.legend.lower() != "true":
        if args.legend.lower() != "false":
            sys.stderr.write("TCRcloud error: please indicate \
True or False\n")
            exit()
    if args.export.lower() != "true":
        if args.export.lower() != "false":
            sys.stderr.write("TCRcloud error: please indicate \
True or False\n")
            exit()

    categories = ["D50\nIndex",
                  "Distinct\nCDR3",
                  "Chao1\nIndex",
                  "Simpson\nIndex",
                  "Shannon\nIndex",
                  "Gini\nindex",
                  "Convergence"
                  ]

    categories = [*categories, categories[0]]
    samples_df = tcrcloud.format.format_data(args)
    samples = samples_df.groupby(["chain", "repertoire_id"])
    keys = [key for key, _ in samples]
    patients = calculate_metrics(keys, samples, args.custom_legend,
                                 args.export, args.rearrangements)

    label_loc = np.linspace(start=0, stop=2 * np.pi, num=len(patients[0]))

    plt.figure(figsize=(10, 14))
    plt.subplot(polar=True)

    radar_colours = ["#332288",
                     "#44AA99",
                     "#882255",
                     "#CC6677",
                     "#117733",
                     "#999933",
                     # "#AA4499",
                     "#88CCEE",
                     "#DDCC77"
                     ]

    for i in patients:
        patient = i[1:]
        patient = [*patient, patient[0]]

        try:
            thecolour = radar_colours.pop(0)
        except IndexError:
            thecolour = "#BBBBBB"

        plt.plot(label_loc,
                 patient,
                 # label=i[0],
                 linewidth=4.0,
                 alpha=0.4,
                 color=thecolour)

        plt.fill(label_loc,
                 patient,
                 label=i[0],
                 linewidth=4.0,
                 alpha=0.6,
                 color=thecolour)

    plt.text(label_loc[0],
             0.11,
             "5",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[0],
             0.26,
             "15",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[0],
             0.46,
             "25",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[0],
             0.66,
             "35",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[0],
             0.86,
             "45",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[0],
             1.00,
             "50",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")

    plt.text(label_loc[1],
             0.1,
             "5500",
             horizontalalignment="left",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[1],
             0.26,
             "16500",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[1],
             0.46,
             "27500",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[1],
             0.66,
             "38500",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[1],
             0.86,
             "49500",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[1],
             1.00,
             "55000",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")

    plt.text(label_loc[2],
             0.11,
             "9000",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[2],
             0.26,
             "27000",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[2],
             0.46,
             "45000",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[2],
             0.66,
             "63000",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[2],
             0.86,
             "81000",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[2],
             1.00,
             "90000",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")

    plt.text(label_loc[3],
             0.11,
             "0.78",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[3],
             0.26,
             "0.83",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[3],
             0.46,
             "0.88",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[3],
             0.66,
             "0.92",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[3],
             0.86,
             "0.97",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[3],
             1.01,
             "1",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")

    plt.text(label_loc[4],
             0.11,
             "1.5",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[4],
             0.26,
             "4.5",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[4],
             0.46,
             "7.5",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[4],
             0.66,
             "10.5",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[4],
             0.86,
             "13.5",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[4],
             1.00,
             "15",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")

    plt.text(label_loc[5],
             0.11,
             "0.1",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[5],
             0.26,
             "0.3",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[5],
             0.46,
             "0.5",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[5],
             0.66,
             "0.7",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[5],
             0.86,
             "0.9",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[5],
             1.01,
             "1",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")

    plt.text(label_loc[6],
             0.11,
             "0.1",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[6],
             0.26,
             "0.3",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[6],
             0.46,
             "0.5",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[6],
             0.66,
             "0.7",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[6],
             0.86,
             "0.9",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")
    plt.text(label_loc[6],
             1.01,
             "1",
             horizontalalignment="center",
             verticalalignment="center",
             fontsize=12,
             fontweight="bold")

    plt.ylim(0, 1.01)
    plt.yticks([0.1, 0.3, 0.5, 0.7, 0.9], [])
    plt.tick_params(pad=32, labelsize=16)
    lines, labels = plt.thetagrids(np.degrees(label_loc), labels=categories)
    outputname = args.rearrangements[:-4] + "_radar" + ".png"
    if args.legend.lower() == "true":
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1),
                   fontsize=16)
    plt.savefig(outputname, dpi=300, bbox_inches="tight")
    print("Radar saved as " + outputname)
