"""Unit tests for tcrcloud.TCRcloud (the CLI entry point)."""

import sys

import pytest

import tcrcloud.cloud as cloud
import tcrcloud.TCRcloud as TCRcloud

# ---------------------------------------------------------------------------
# main() FileNotFoundError handling
# ---------------------------------------------------------------------------
#
# Regression tests: main() used to unconditionally rebuild the error message
# from `args.rearrangements`, discarding any custom, already well-formatted
# "TCRcloud error: ..." message raised by the subcommand (e.g. cloud.py
# reporting that a --colours file doesn't exist). That made a missing
# colours.json report the rearrangements filename instead.


def test_main_preserves_custom_filenotfound_message(monkeypatch, capsys):
    def fake_wordcloud(args):
        raise FileNotFoundError("TCRcloud error: colours.json doesn't seem to exist")

    monkeypatch.setattr(cloud, "wordcloud", fake_wordcloud)
    monkeypatch.setattr(sys, "argv", ["TCRcloud", "cloud", "rearrangements.tsv"])

    TCRcloud.main()

    captured = capsys.readouterr()
    assert "colours.json doesn't seem to exist" in captured.err
    assert "rearrangements.tsv" not in captured.err


def test_main_falls_back_to_argument_filename_for_generic_filenotfound(
    monkeypatch, capsys
):
    def fake_wordcloud(args):
        # A plain FileNotFoundError with no custom "TCRcloud error:" message,
        # e.g. raised by a bare `open()` call somewhere downstream.
        raise FileNotFoundError()

    monkeypatch.setattr(cloud, "wordcloud", fake_wordcloud)
    monkeypatch.setattr(sys, "argv", ["TCRcloud", "cloud", "missing.tsv"])

    TCRcloud.main()

    captured = capsys.readouterr()
    assert "missing.tsv doesn't seem to exist" in captured.err


def test_main_preserves_custom_message_and_does_not_raise(monkeypatch):
    """main() should swallow the FileNotFoundError, not propagate it."""

    def fake_wordcloud(args):
        raise FileNotFoundError("TCRcloud error: colours.json doesn't seem to exist")

    monkeypatch.setattr(cloud, "wordcloud", fake_wordcloud)
    monkeypatch.setattr(sys, "argv", ["TCRcloud", "cloud", "rearrangements.tsv"])

    TCRcloud.main()  # should not raise


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
