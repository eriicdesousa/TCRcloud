"""Unit tests for tcrcloud.radar."""

import argparse
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import skbio

import tcrcloud.format as tformat
import tcrcloud.radar as radar

# ---------------------------------------------------------------------------
# calculate_dfifty
# ---------------------------------------------------------------------------


def test_calculate_dfifty_basic_case():
    # Sorted descending counts [4, 2, 1]; cumulative sum first exceeds
    # total/2 (3.5) at the first clone -> 1 of 3 clones -> 33.33%.
    df = pd.DataFrame({"counts": [4, 2, 1]})
    assert radar.calculate_dfifty(df, len(df)) == pytest.approx(100 / 3)


def test_calculate_dfifty_returns_zero_when_total_is_zero():
    df = pd.DataFrame({"counts": []})
    assert radar.calculate_dfifty(df, 0) == 0.0


def test_calculate_dfifty_truncates_to_top_10000_when_length_is_larger():
    # Two large clones followed by many singletons; only the first 10000
    # rows are considered so the extra singleton beyond that must not
    # affect the result.
    counts = np.concatenate([[6000, 6000], np.ones(9999)])
    df = pd.DataFrame({"counts": counts})
    result = radar.calculate_dfifty(df, len(counts))
    assert result == pytest.approx(2 * 100.0 / 10000)


# ---------------------------------------------------------------------------
# _load_legend_mapping
# ---------------------------------------------------------------------------


def test_load_legend_mapping_reads_valid_json(tmp_path):
    legend_file = tmp_path / "legend.json"
    legend_file.write_text(json.dumps({"rep1": "Patient 1"}))
    assert radar._load_legend_mapping(str(legend_file)) == {"rep1": "Patient 1"}


def test_load_legend_mapping_missing_file_raises_filenotfound(tmp_path):
    missing_path = tmp_path / "missing_legend.json"
    with pytest.raises(
        FileNotFoundError, match=r"TCRcloud error:.*doesn't seem to exist"
    ):
        radar._load_legend_mapping(str(missing_path))


def test_load_legend_mapping_invalid_json_raises_value_error(tmp_path):
    bad_file = tmp_path / "bad_legend.json"
    bad_file.write_text("{not valid json")
    with pytest.raises(ValueError, match="doesn't seem properly formatted"):
        radar._load_legend_mapping(str(bad_file))


# ---------------------------------------------------------------------------
# calculate_metrics
# ---------------------------------------------------------------------------
#
# Synthetic "raw" rearrangement rows (one row per rearrangement, pre
# aggregation) for a single (chain, repertoire_id) group, aggregating (by
# junction_aa) to counts of 4 (CASSA), 2 (CASSB), 1 (CASSC).


def _make_group_df(chain="B"):
    return pd.DataFrame(
        {
            "junction_aa": [
                "CASSA",
                "CASSA",
                "CASSA",
                "CASSA",
                "CASSB",
                "CASSB",
                "CASSC",
            ],
            "junction": [
                "TGTA",
                "TGTA",
                "TGTA",
                "TGTB",
                "TGTC",
                "TGTC",
                "TGTD",
            ],
            "repertoire_id": ["rep1"] * 7,
            "chain": [chain] * 7,
        }
    )


def _make_samples(chain="B"):
    df = _make_group_df(chain=chain)
    samples = df.groupby(["chain", "repertoire_id"])
    return samples, list(samples.groups.keys())


def test_calculate_metrics_returns_one_dataset_per_key():
    samples, keys = _make_samples()
    datasets, min_vals, max_vals, scales = radar.calculate_metrics(
        keys, samples, None, False, "sample.tsv"
    )

    assert len(datasets) == 1
    label, *scaled_values = datasets[0]
    assert label == "rep1 β chain"
    assert len(scaled_values) == 7
    # All scaled metric values must fall within the plotted [0, 1] range.
    assert all(0.0 <= v <= 1.0 for v in scaled_values)


def test_calculate_metrics_applies_custom_legend(tmp_path):
    legend_file = tmp_path / "legend.json"
    legend_file.write_text(json.dumps({"rep1": "Patient 1"}))

    samples, keys = _make_samples()
    datasets, *_ = radar.calculate_metrics(
        keys, samples, str(legend_file), False, "sample.tsv"
    )

    assert datasets[0][0] == "Patient 1 β chain"


def test_calculate_metrics_missing_legend_file_raises(tmp_path):
    samples, keys = _make_samples()
    missing_path = tmp_path / "missing.json"

    with pytest.raises(FileNotFoundError, match="doesn't seem to exist"):
        radar.calculate_metrics(keys, samples, str(missing_path), False, "sample.tsv")


def test_calculate_metrics_exports_repertoire_metrics_file(tmp_path):
    samples, keys = _make_samples()
    filename = str(tmp_path / "sample.tsv")

    radar.calculate_metrics(keys, samples, None, True, filename)

    metrics_path = tmp_path / "sample_repertoire_metrics.txt"
    assert metrics_path.exists()
    content = metrics_path.read_text()
    assert "Repertoire: rep1 β chain" in content
    assert "D50 Index:" in content

    expected_pielou = float(
        np.format_float_positional(
            skbio.diversity.alpha_diversity("pielou_e", [4.0, 2.0, 1.0])[0],
            precision=3,
        )
    )
    assert f"Pielou Evenness: {expected_pielou}" in content


def test_calculate_metrics_returns_empty_when_no_keys():
    samples, _ = _make_samples()
    datasets, min_vals, max_vals, scales = radar.calculate_metrics(
        [], samples, None, False, "sample.tsv"
    )
    assert datasets == []


def test_calculate_metrics_scaled_values_match_expected_scaling():
    # Verifies the vectorized scaling in calculate_metrics agrees with the
    # single-value _scale_value_to_01 helper used to draw axis ticks, for
    # every metric (not just that values fall within [0, 1]).
    samples, keys = _make_samples()
    datasets, min_vals, max_vals, scales = radar.calculate_metrics(
        keys, samples, None, False, "sample.tsv"
    )
    _, *scaled_values = datasets[0]

    counts = [4.0, 2.0, 1.0]
    raw_values = [
        radar.calculate_dfifty(pd.DataFrame({"counts": counts}), len(counts)),
        1 - skbio.diversity.alpha_diversity("gini_index", counts)[0],
        skbio.diversity.alpha_diversity("shannon", counts)[0],
        skbio.diversity.alpha_diversity("simpson", counts)[0],
        len(counts),
        skbio.diversity.alpha_diversity("chao1", counts)[0],
        skbio.diversity.alpha_diversity("pielou_e", counts)[0],
    ]

    for value, min_val, max_val, scale, actual_scaled in zip(
        raw_values, min_vals, max_vals, scales, scaled_values, strict=True
    ):
        expected_scaled = np.clip(
            radar._scale_value_to_01(value, min_val, max_val, scale), 0.0, 1.0
        )
        assert actual_scaled == pytest.approx(expected_scaled)


@pytest.mark.parametrize(
    "chain,expected_suffix",
    [
        ("A", "α chain"),
        ("B", "β chain"),
        ("G", "γ chain"),
        ("D", "δ chain"),
        ("H", "Heavy chain"),
        ("K", "Kappa chain"),
        ("L", "Lambda chain"),
        ("Z", ""),  # unrecognized chain code falls back to no suffix
    ],
)
def test_calculate_metrics_labels_all_chain_suffixes(chain, expected_suffix):
    samples, keys = _make_samples(chain=chain)
    datasets, *_ = radar.calculate_metrics(keys, samples, None, False, "sample.tsv")

    expected_label = f"rep1 {expected_suffix}".strip()
    assert datasets[0][0] == expected_label


# ---------------------------------------------------------------------------
# _scale_value_to_01 / _format_tick (axis tick helpers used by radar())
# ---------------------------------------------------------------------------


def test_scale_value_to_01_linear_midpoint():
    assert radar._scale_value_to_01(25, 0.0, 50.0, "linear") == pytest.approx(0.5)


def test_scale_value_to_01_linear_clips_out_of_range_values():
    assert radar._scale_value_to_01(-10, 0.0, 50.0, "linear") == pytest.approx(0.0)
    assert radar._scale_value_to_01(1000, 0.0, 50.0, "linear") == pytest.approx(1.0)


def test_scale_value_to_01_log_endpoints():
    assert radar._scale_value_to_01(1.0, 1.0, 250000.0, "log") == pytest.approx(0.0)
    assert radar._scale_value_to_01(250000.0, 1.0, 250000.0, "log") == pytest.approx(
        1.0
    )


def test_scale_value_to_01_log_midpoint_matches_log10_formula():
    value, min_val, max_val = 500.0, 1.0, 250000.0
    expected = (np.log10(value) - np.log10(min_val)) / (
        np.log10(max_val) - np.log10(min_val)
    )
    assert radar._scale_value_to_01(value, min_val, max_val, "log") == pytest.approx(
        expected
    )


def test_format_tick_integer_metric_rounds_to_int():
    result = radar._format_tick(4.6, True)
    assert result == 5
    assert isinstance(result, int)


def test_format_tick_large_value_uses_comma_grouping():
    assert radar._format_tick(250000, False) == "250,000"


def test_format_tick_small_value_avoids_scientific_notation():
    assert radar._format_tick(0.00012345, False) == "0.000123"


def test_format_tick_default_rounds_to_five_decimals():
    assert radar._format_tick(3.14159265, False) == pytest.approx(3.14159)


# ---------------------------------------------------------------------------
# radar() entry point
# ---------------------------------------------------------------------------


def _make_args(tmp_path, **overrides):
    rearrangements = tmp_path / "repertoire.tsv"
    rearrangements.write_text(
        "junction_aa\tjunction\trepertoire_id\tchain\tproductive\n"
    )
    defaults = dict(
        rearrangements=str(rearrangements),
        custom_legend=None,
        legend=True,
        export=False,
        format="png",
    )
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


def _patch_format_data(monkeypatch):
    monkeypatch.setattr(tformat, "format_data", lambda args: _make_group_df())


def test_radar_writes_png_by_default(tmp_path, monkeypatch):
    _patch_format_data(monkeypatch)
    monkeypatch.chdir(tmp_path)

    radar.radar(_make_args(tmp_path))

    outputs = list(tmp_path.glob("*_radar.png"))
    assert [p.name for p in outputs] == ["repertoire_radar.png"]
    assert not list(tmp_path.glob("*_radar.svg"))


def test_radar_writes_svg_when_requested(tmp_path, monkeypatch):
    _patch_format_data(monkeypatch)
    monkeypatch.chdir(tmp_path)

    radar.radar(_make_args(tmp_path, format="svg"))

    outputs = list(tmp_path.glob("*_radar.svg"))
    assert [p.name for p in outputs] == ["repertoire_radar.svg"]
    assert not list(tmp_path.glob("*_radar.png"))


def test_radar_format_is_case_insensitive(tmp_path, monkeypatch):
    _patch_format_data(monkeypatch)
    monkeypatch.chdir(tmp_path)

    radar.radar(_make_args(tmp_path, format="SVG"))

    assert [p.name for p in tmp_path.glob("*_radar.svg")] == ["repertoire_radar.svg"]


def test_radar_rejects_unsupported_format(tmp_path, monkeypatch):
    _patch_format_data(monkeypatch)

    with pytest.raises(ValueError, match="unsupported output format"):
        radar.radar(_make_args(tmp_path, format="pdf"))


def test_radar_defaults_to_png_when_format_missing(tmp_path, monkeypatch):
    # `radar()` may be called with an args object that predates the
    # --format option (e.g. older scripts); it should default to png.
    _patch_format_data(monkeypatch)
    monkeypatch.chdir(tmp_path)
    args = _make_args(tmp_path)
    del args.format

    radar.radar(args)

    assert (tmp_path / "repertoire_radar.png").exists()


def test_radar_accepts_string_booleans_for_legend_and_export(tmp_path, monkeypatch):
    # Older callers may pass string booleans instead of real bools.
    _patch_format_data(monkeypatch)
    monkeypatch.chdir(tmp_path)

    radar.radar(_make_args(tmp_path, legend="false", export="true"))

    assert (tmp_path / "repertoire_radar.png").exists()
    assert (tmp_path / "repertoire_repertoire_metrics.txt").exists()


def test_radar_raises_when_no_repertoires_found(tmp_path, monkeypatch):
    empty_df = pd.DataFrame(
        columns=["junction_aa", "junction", "repertoire_id", "chain"]
    )
    monkeypatch.setattr(tformat, "format_data", lambda args: empty_df)
    with pytest.raises(ValueError, match="no repertoires found"):
        radar.radar(_make_args(tmp_path))


def test_radar_closes_figure_after_saving(tmp_path, monkeypatch):
    # Regression test: radar() must close the figure it creates, otherwise
    # repeated calls (e.g. batch processing) leak matplotlib figures.
    _patch_format_data(monkeypatch)
    monkeypatch.chdir(tmp_path)
    plt.close("all")

    radar.radar(_make_args(tmp_path))

    assert plt.get_fignums() == []


def test_radar_calls_legend_only_when_enabled(tmp_path, monkeypatch):
    _patch_format_data(monkeypatch)
    monkeypatch.chdir(tmp_path)
    legend_calls = []
    monkeypatch.setattr(plt, "legend", lambda *a, **k: legend_calls.append((a, k)))

    radar.radar(_make_args(tmp_path, legend=True))
    assert len(legend_calls) == 1

    radar.radar(_make_args(tmp_path, legend=False))
    assert len(legend_calls) == 1  # unchanged: plt.legend() not called again
