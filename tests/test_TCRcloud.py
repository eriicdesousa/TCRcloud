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
