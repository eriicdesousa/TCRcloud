"""Radar plot generation for TCRcloud.

This module computes per-repertoire diversity metrics and renders them as a
radar/spider plot. Some metrics span many orders of magnitude, so a mix of
linear, log, and tail-log scaling is used for display.

"""

import json
import sys

import numpy as np
import matplotlib.pyplot as plt
import skbio

import tcrcloud.format

# For large-range metrics (Distinct CDR3, Chao1), the axis uses log scaling.
# For the Gini-Simpson metric, we use a "tail-log" scale to better separate
# values close to 1.
_METRIC_RANGES = np.array(
    [
        [0.0, 50.0],  # D50 Index (linear)
        [0.0, 1.0],  # Gini Coefficient (linear)
        [0.0, 15],  # Shannon Index (linear)
        [0.7, 1.0],  # Gini-Simpson Index (tail-log)
        [1.0, 250000.0],  # Distinct CDR3 (log)
        [1.0, 250000.0],  # Chao1 Index (log)
        [0.000001, 1.0],  # Convergence (log)
    ],
    dtype=float,
)

# Scale type per metric (used when mapping raw values into [0,1]).
_METRIC_SCALES = [
    "linear",  # D50
    "linear",  # Gini
    "linear",  # Shannon
    "tail-log",  # Gini-Simpson (emphasize values close to 1)
    "log",  # Distinct CDR3
    "log",  # Chao1
    "log",  # Convergence
]


def calculate_convergence(df):
    """Compute the convergence metric.

    Convergence is defined as the proportion of total counts that come from
    CDR3 sequences that have more than one distinct junction (i.e., convergent
    recombination events).
    """

    counts = np.asarray(df["counts"], dtype=float)
    total = counts.sum()
    if total <= 0:
        return 0.0

    return float(counts[counts > 1].sum()) / float(total)


def calculate_dfifty(df, length):
    """Compute D50 (% clones needed to reach 50% of total counts).

    D50 is the percentage of unique clones required to accumulate 50% of the
    total counts. This is computed on the top 10,000 clones to keep runtime
    bounded for large repertoires.
    """

    counts = np.asarray(df["counts"], dtype=float)
    if length >= 10000:
        counts = counts[:10000]

    total = counts.sum()
    if total <= 0:
        return 0.0

    # Find how many top clones are needed to exceed 50% cumulative count.
    idx = int(np.searchsorted(np.cumsum(counts), total / 2, side="right"))
    return (idx + 1) * 100.0 / (10000 if length >= 10000 else length)


def calculate_metrics(keys, samples, legend_file, export, filename):
    """Compute all radar metrics for each repertoire.

    This function computes a fixed set of repertoire diversity metrics and
    then scales them into the [0, 1] range according to the fixed ranges
    defined in _METRIC_RANGES/_METRIC_SCALES.

    Returns:
        (datasets, min_vals, max_vals, scales, raw_arr, scaled_arr)
    """

    legend_dict = {}
    if legend_file is not None:
        try:
            with open(legend_file) as json_file:
                legend_dict = json.load(json_file)
        except FileNotFoundError:
            sys.stderr.write(f"TCRcloud error: {legend_file} doesn't seem to exist\n")
            exit(1)
        except json.decoder.JSONDecodeError:
            sys.stderr.write(
                "TCRcloud error: "
                + legend_file
                + " doesn't seem properly formatted. Check https://github.com/oldguyeric/TCRcloud for more information\n"
            )
            exit(1)

    chain_names = {
        "A": "α chain",
        "B": "β chain",
        "G": "γ chain",
        "D": "δ chain",
        "H": "Heavy chain",
        "K": "Kappa chain",
        "L": "Lambda chain",
    }

    def _format_label(key):
        prefix = legend_dict.get(key[1], key[1])
        suffix = chain_names.get(key[0], "")
        return f"{prefix} {suffix}".strip()

    raw_metrics = []
    labels = []

    for key in keys:
        df = samples.get_group(key)

        convergence = calculate_convergence(tcrcloud.format.format_convergence(df))

        df = tcrcloud.format.format_metrics(df)
        length = len(df)
        counts = df["counts"].to_numpy(dtype=float)

        raw_metrics.append(
            [
                calculate_dfifty(df, length),
                skbio.diversity.alpha_diversity("gini_index", counts)[0],
                skbio.diversity.alpha_diversity("shannon", counts)[0],
                skbio.diversity.alpha_diversity("simpson", counts)[0],
                length,
                skbio.diversity.alpha_diversity("chao1", counts)[0],
                convergence,
            ]
        )
        labels.append(_format_label(key))

    if not raw_metrics:
        return [], np.zeros(7), np.zeros(7), _METRIC_SCALES

    raw_arr = np.asarray(raw_metrics, dtype=float)

    # Scale using fixed ranges (pre-refactor behavior).
    min_vals = _METRIC_RANGES[:, 0]
    max_vals = _METRIC_RANGES[:, 1]

    # Vectorize scaling for speed.
    scaled_arr = np.empty_like(raw_arr, dtype=float)

    # Apply log transform for log-scaled metrics; keep all others linear.
    is_log = np.array([s == "log" for s in _METRIC_SCALES], dtype=bool)
    if np.any(is_log):
        # Clip to avoid log of zero/negative values.
        raw_log = np.clip(raw_arr[:, is_log], min_vals[is_log], None)
        raw_log = np.log10(raw_log)
        min_log = np.log10(min_vals[is_log])
        max_log = np.log10(max_vals[is_log])
        scaled_arr[:, is_log] = (raw_log - min_log) / (max_log - min_log)

    if np.any(~is_log):
        raw_lin = raw_arr[:, ~is_log]
        min_lin = min_vals[~is_log]
        max_lin = max_vals[~is_log]
        scaled_arr[:, ~is_log] = (raw_lin - min_lin) / (max_lin - min_lin)

    scaled_arr = np.clip(scaled_arr, 0.0, 1.0)

    datasets = [[label, *row.tolist()] for label, row in zip(labels, scaled_arr)]

    if export:
        metrics_filename = filename[:-4] + "_repertoire_metrics.txt"
        with open(metrics_filename, "w") as fileout:
            for label, row in zip(labels, raw_arr):
                print("Repertoire:", label, file=fileout)
                print(
                    "D50 Index:",
                    float(np.format_float_positional(row[0], precision=3)),
                    file=fileout,
                )
                print(
                    "Gini Coefficient:",
                    float(np.format_float_positional(row[1], precision=3)),
                    file=fileout,
                )
                print(
                    "Shannon Index:",
                    float(np.format_float_positional(row[2], precision=3)),
                    file=fileout,
                )
                print(
                    "Gini-Simpson Index:",
                    float(np.format_float_positional(row[3], precision=3)),
                    file=fileout,
                )
                print(
                    "Distinct CDR3:",
                    float(np.format_float_positional(row[4], precision=3)),
                    file=fileout,
                )
                print(
                    "Chao1 Index:",
                    float(np.format_float_positional(row[5], precision=3)),
                    file=fileout,
                )
                print(
                    "Convergence:",
                    float(np.format_float_positional(row[6], precision=3)),
                    file=fileout,
                )
                print(file=fileout)
        print("Repertoire metrics saved as " + metrics_filename)

    return datasets, min_vals, max_vals, _METRIC_SCALES


def radar(args):
    # argparse handles boolean conversion when using str2bool.
    # Normalize boolean-style CLI flags (allow strings like "true"/"false").
    if isinstance(args.legend, str):
        args.legend = args.legend.lower() in ("yes", "true", "t", "y", "1")
    if isinstance(args.export, str):
        args.export = args.export.lower() in ("yes", "true", "t", "y", "1")

    # The radar plot categories (one per metric) + repeat first category for closure.
    categories = [
        "D50\nIndex",
        "Gini\nCoefficient",
        "Shannon\nIndex",
        "Gini-Simpson\nIndex",
        "Distinct\nCDR3",
        "Chao1\nIndex",
        "Convergence",
    ]
    categories = [*categories, categories[0]]

    # Load and filter the input repertoire TSV file, then group by chain and repertoire.
    samples_df = tcrcloud.format.format_data(args)
    samples = samples_df.groupby(["chain", "repertoire_id"])
    keys = list(samples.groups.keys())

    datasets, min_vals, max_vals, scales = calculate_metrics(
        keys,
        samples,
        args.custom_legend,
        args.export,
        args.rearrangements,
    )

    if not datasets:
        sys.stderr.write("TCRcloud error: no repertoires found for plotting\n")
        return

    label_loc = np.linspace(start=0, stop=2 * np.pi, num=len(datasets[0]))

    plt.figure(figsize=(10, 14))
    plt.subplot(polar=True)

    # Predefined palette for the radar lines; will cycle if there are more
    # repertoires than colors.
    radar_colours = [
        "#f0e442",
        "#0072b2",
        "#cc79a7",
        "#009e73",
        "#e69f00",
        "#56b4e9",
    ]

    for i in datasets:
        dataset = i[1:]
        dataset = [*dataset, dataset[0]]

        try:
            thecolour = radar_colours.pop(0)
        except IndexError:
            thecolour = "#BBBBBB"

        plt.plot(
            label_loc,
            dataset,
            # label=i[0],
            linewidth=4.0,
            alpha=0.4,
            color=thecolour,
        )

        plt.fill(
            label_loc, dataset, label=i[0], linewidth=4.0, alpha=0.6, color=thecolour
        )

    integer_metrics = [False, False, False, False, True, True, False]

    # Draw tick labels for each axis based on the fixed default ranges.
    def _scale_value_to_01(value, min_val, max_val, scale):
        if scale == "log":
            min_t = np.log10(min_val)
            max_t = np.log10(max_val)
            value = np.clip(value, min_val, max_val)
            return (np.log10(value) - min_t) / (max_t - min_t)
        v = np.clip(value, min_val, max_val)
        return (v - min_val) / (max_val - min_val)

    def _format_tick(value, is_int):
        if is_int:
            return int(round(value))
        if value >= 1000:
            # Avoid compact formats like 10k/1M; show full integer with commas.
            return f"{int(round(value)):,}"
        # Avoid scientific notation for very small values (e.g. convergence on log scale)
        if 0 < abs(value) < 1e-3:
            return f"{value:.6f}".rstrip("0").rstrip(".")
        return round(value, 5)

    for idx, (min_val, max_val, scale, is_int) in enumerate(
        zip(min_vals, max_vals, scales, integer_metrics)
    ):
        if scale == "log":
            lo = int(np.floor(np.log10(min_val)))
            hi = int(np.ceil(np.log10(max_val)))
            tick_values = [10**e for e in range(lo, hi + 1)]
        else:
            tick_values = [
                min_val + (max_val - min_val) * f for f in (0.1, 0.3, 0.5, 0.7, 0.9)
            ]

        # Optionally suppress the first (minimum) tick for cleaner axis labeling.
        if len(tick_values) > 1 and abs(tick_values[0] - min_val) < 1e-12:
            tick_values = tick_values[1:]

        # For Distinct CDR3 and Chao1, skip rendering the last tick value
        if idx in (4, 5) and len(tick_values) > 1:
            tick_values = tick_values[:-1]

        for tick in tick_values:
            pos = _scale_value_to_01(tick, min_val, max_val, scale)
            label = _format_tick(tick, is_int)
            plt.text(
                label_loc[idx],
                pos,
                label,
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=12,
                fontweight="bold",
            )

        max_label = _format_tick(max_val, is_int)
        # Strip trailing .0 for axis max values to keep labels clean
        if isinstance(max_label, float) and max_label.is_integer():
            max_label = int(max_label)
        plt.text(
            label_loc[idx],
            1.0,
            max_label,
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
            fontweight="bold",
        )

    plt.ylim(0, 1.01)
    plt.yticks([0.1, 0.3, 0.5, 0.7, 0.9], [])
    plt.tick_params(pad=32, labelsize=16)
    lines, labels = plt.thetagrids(np.degrees(label_loc), labels=categories)
    outputname = args.rearrangements[:-4] + "_radar" + ".svg"
    if args.legend:
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), fontsize=16)
    plt.savefig(outputname, dpi=300, bbox_inches="tight")
    print("Radar saved as " + outputname)
