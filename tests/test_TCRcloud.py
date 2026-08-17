"""Unit tests for tcrcloud.TCRcloud (the CLI entry point)."""

import sys

import pytest

import tcrcloud.cloud as cloud
import tcrcloud.TCRcloud as TCRcloud
from tcrcloud.errors import TCRcloudError

# ---------------------------------------------------------------------------
# main() error handling
# ---------------------------------------------------------------------------
#
# Regression tests: main() used to unconditionally rebuild the error message
# from `args.rearrangements`, discarding any custom, already well-formatted
# message raised by the subcommand (e.g. cloud.py reporting that a --colours
# file doesn't exist). That made a missing colours.json report the
# rearrangements filename instead. Subcommands now raise TCRcloudError with
# the complete message (no "TCRcloud error: " prefix - main() adds it).


def test_main_preserves_custom_error_message(monkeypatch, capsys):
    def fake_wordcloud(*args, **kwargs):
        # Library modules raise TCRcloudError with a complete message (no
        # prefix); main() must surface it verbatim instead of rebuilding a
        # message from the rearrangements argument.
        raise TCRcloudError("colours.json doesn't seem to exist")

    monkeypatch.setattr(cloud, "wordcloud", fake_wordcloud)
    monkeypatch.setattr(sys, "argv", ["TCRcloud", "cloud", "rearrangements.tsv"])

    with pytest.raises(SystemExit) as excinfo:
        TCRcloud.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "TCRcloud error: colours.json doesn't seem to exist" in captured.err
    assert "rearrangements.tsv" not in captured.err


def test_main_falls_back_to_argument_filename_for_generic_filenotfound(
    monkeypatch, capsys
):
    def fake_wordcloud(*args, **kwargs):
        # A plain FileNotFoundError with no custom "TCRcloud error:" message,
        # e.g. raised by a bare `open()` call somewhere downstream.
        raise FileNotFoundError()

    monkeypatch.setattr(cloud, "wordcloud", fake_wordcloud)
    monkeypatch.setattr(sys, "argv", ["TCRcloud", "cloud", "missing.tsv"])

    with pytest.raises(SystemExit) as excinfo:
        TCRcloud.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "missing.tsv doesn't seem to exist" in captured.err


def test_main_exits_nonzero_on_tcrcloud_error(monkeypatch, capsys):
    """TCRcloudError raised by a subcommand must surface as a clean message
    plus a non-zero exit code (not a swallowed print)."""

    def fake_wordcloud(*args, **kwargs):
        raise TCRcloudError("something went wrong")

    monkeypatch.setattr(cloud, "wordcloud", fake_wordcloud)
    monkeypatch.setattr(sys, "argv", ["TCRcloud", "cloud", "rearrangements.tsv"])

    with pytest.raises(SystemExit) as excinfo:
        TCRcloud.main()

    assert excinfo.value.code == 1
    captured = capsys.readouterr()
    assert "TCRcloud error: something went wrong" in captured.err


# ---------------------------------------------------------------------------
# main() must NOT mask programming errors
# ---------------------------------------------------------------------------
#
# Regression tests: main() used to catch ValueError, KeyError and TypeError
# broadly, which conflated genuine programming bugs with bad user input and
# reported them as "you did not indicate a properly formatted file". Those
# handlers are gone; user-input errors are validated and converted to
# TCRcloudError at the module boundary instead, so internal errors must now
# propagate uncaught (with a traceback) rather than being masked.


def test_main_propagates_keyerror_from_subcommand(monkeypatch, capsys):
    def fake_wordcloud(*args, **kwargs):
        raise KeyError("internal bug")

    monkeypatch.setattr(cloud, "wordcloud", fake_wordcloud)
    monkeypatch.setattr(sys, "argv", ["TCRcloud", "cloud", "rearrangements.tsv"])

    with pytest.raises(KeyError):
        TCRcloud.main()

    # The bug must NOT be masked as a malformed-input message.
    assert "properly formatted" not in capsys.readouterr().err


def test_main_propagates_valueerror_from_subcommand(monkeypatch, capsys):
    def fake_wordcloud(*args, **kwargs):
        raise ValueError("internal bug")

    monkeypatch.setattr(cloud, "wordcloud", fake_wordcloud)
    monkeypatch.setattr(sys, "argv", ["TCRcloud", "cloud", "rearrangements.tsv"])

    with pytest.raises(ValueError):
        TCRcloud.main()

    assert "properly formatted" not in capsys.readouterr().err


# ---------------------------------------------------------------------------
# main() argument validation
# ---------------------------------------------------------------------------
#
# Invalid argument values should abort with an argparse error (exit code 2)
# and a message on stderr instead of silently running with wrong assumptions.


def _assert_invalid_args_exit_with_error(monkeypatch, capsys, argv, message_part):
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit) as excinfo:
        TCRcloud.main()

    assert excinfo.value.code == 2
    captured = capsys.readouterr()
    assert message_part in captured.err
    assert "usage:" in captured.err


def test_main_rejects_invalid_species_for_cloud(monkeypatch, capsys):
    _assert_invalid_args_exit_with_error(
        monkeypatch,
        capsys,
        ["TCRcloud", "cloud", "rearrangements.tsv", "-p", "mus_musculs"],
        "invalid choice: 'mus_musculs'",
    )


def test_main_rejects_invalid_species_for_vgenes(monkeypatch, capsys):
    _assert_invalid_args_exit_with_error(
        monkeypatch,
        capsys,
        ["TCRcloud", "vgenes", "rearrangements.tsv", "-p", "dog"],
        "invalid choice: 'dog'",
    )


def test_main_rejects_invalid_species_for_compare(monkeypatch, capsys):
    _assert_invalid_args_exit_with_error(
        monkeypatch,
        capsys,
        ["TCRcloud", "compare", "rearrangements.tsv", "-p", "human"],
        "invalid choice: 'human'",
    )


def test_main_rejects_invalid_boolean_for_aminoacids_export(monkeypatch, capsys):
    _assert_invalid_args_exit_with_error(
        monkeypatch,
        capsys,
        ["TCRcloud", "aminoacids", "rearrangements.tsv", "-e", "banana"],
        "Expected a boolean value",
    )


def test_main_rejects_invalid_boolean_for_aminoacids_threed(monkeypatch, capsys):
    _assert_invalid_args_exit_with_error(
        monkeypatch,
        capsys,
        ["TCRcloud", "aminoacids", "rearrangements.tsv", "-t", "yes please"],
        "Expected a boolean value",
    )


def test_main_false_legend_value_turns_off_legend(monkeypatch, capsys):
    captured = {}

    def fake_wordcloud(*args, **kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(cloud, "wordcloud", fake_wordcloud)
    monkeypatch.setattr(
        sys, "argv", ["TCRcloud", "cloud", "rearrangements.tsv", "-l", "False"]
    )

    TCRcloud.main()

    assert captured["legend"] is False


def test_main_true_export_value_turns_on_export(monkeypatch):
    captured = {}

    def fake_radar(*args, **kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(TCRcloud.tcrcloud.radar, "radar", fake_radar)
    monkeypatch.setattr(
        sys, "argv", ["TCRcloud", "radar", "rearrangements.tsv", "-e", "True"]
    )

    TCRcloud.main()

    assert captured["export"] is True
