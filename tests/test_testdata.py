"""Unit tests for tcrcloud.testdata."""

import json
from unittest.mock import Mock

import pytest
import requests

import tcrcloud.testdata as testdata
from tcrcloud.errors import TCRcloudError

# ---------------------------------------------------------------------------
# _make_query
# ---------------------------------------------------------------------------


def test_make_query_without_schema():
    query = testdata._make_query("su008", "TRA")

    filters = query["filters"]["content"]
    assert query["filters"]["op"] == "and"
    assert len(filters) == 2
    assert filters[0] == {
        "op": "contains",
        "content": {"field": "subject.subject_id", "value": "su008"},
    }
    assert filters[1] == {
        "op": "in",
        "content": {"field": "sample.pcr_target.pcr_target_locus", "value": ["TRA"]},
    }


def test_make_query_with_schema():
    query = testdata._make_query("su008", "TRB", require_schema=True)

    filters = query["filters"]["content"]
    assert len(filters) == 3
    assert filters[2] == {
        "op": "contains",
        "content": {
            "field": "study.keywords_study",
            "value": "contains_schema_rearrangement",
        },
    }


# ---------------------------------------------------------------------------
# _download_repertoire
# ---------------------------------------------------------------------------


def _make_response(json_payload=None, json_error=None, raise_for_status_error=None):
    """Build a fake `requests.Response`-like object for mocking `session.post`."""

    resp = Mock()
    if raise_for_status_error is not None:
        resp.raise_for_status.side_effect = raise_for_status_error
    else:
        resp.raise_for_status.return_value = None

    if json_error is not None:
        resp.json.side_effect = json_error
    else:
        resp.json.return_value = json_payload
    return resp


def test_download_repertoire_success(tmp_path, monkeypatch, capsys):
    output_path = tmp_path / "alpharepertoire.airr.json"
    payload = {"Repertoire": [{"repertoire_id": "r1"}], "Info": {"title": "test"}}

    session = Mock()
    session.post.return_value = _make_response(json_payload=payload)

    write_repertoire = Mock()
    monkeypatch.setattr(testdata.airr, "write_repertoire", write_repertoire)

    testdata._download_repertoire(session, {"filters": {}}, str(output_path))

    session.post.assert_called_once_with(
        f"{testdata.HOST_URL}/repertoire", json={"filters": {}}, timeout=30
    )
    write_repertoire.assert_called_once_with(
        str(output_path), payload["Repertoire"], info=payload["Info"]
    )

    captured = capsys.readouterr()
    assert "Received 1 repertoires" in captured.out
    assert str(output_path) in captured.out


def test_download_repertoire_network_error_exits(monkeypatch, capsys):
    session = Mock()
    session.post.side_effect = requests.exceptions.ConnectionError("boom")

    with pytest.raises(TCRcloudError) as exc_info:
        testdata._download_repertoire(session, {}, "out.json")

    assert "could not reach" in str(exc_info.value)
    assert testdata.HOST_URL in str(exc_info.value)


def test_download_repertoire_http_error_exits(monkeypatch, capsys):
    session = Mock()
    session.post.return_value = _make_response(
        raise_for_status_error=requests.exceptions.HTTPError("500 error")
    )

    with pytest.raises(TCRcloudError) as exc_info:
        testdata._download_repertoire(session, {}, "out.json")

    assert "could not reach" in str(exc_info.value)


def test_download_repertoire_invalid_json_exits(monkeypatch, capsys):
    session = Mock()
    session.post.return_value = _make_response(json_error=ValueError("bad json"))

    with pytest.raises(TCRcloudError) as exc_info:
        testdata._download_repertoire(session, {}, "out.json")

    assert "invalid JSON" in str(exc_info.value)


def test_download_repertoire_empty_repertoire_exits(monkeypatch, capsys):
    session = Mock()
    session.post.return_value = _make_response(
        json_payload={"Repertoire": [], "Info": {}}
    )
    write_repertoire = Mock()
    monkeypatch.setattr(testdata.airr, "write_repertoire", write_repertoire)

    with pytest.raises(TCRcloudError) as exc_info:
        testdata._download_repertoire(session, {}, "out.json")

    assert "no repertoires were returned" in str(exc_info.value)
    write_repertoire.assert_not_called()


# ---------------------------------------------------------------------------
# download (CLI entry point)
# ---------------------------------------------------------------------------


def test_download_calls_helper_for_both_chains_and_writes_legend(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    fake_session = Mock()
    monkeypatch.setattr(testdata, "get_session", Mock(return_value=fake_session))

    calls = []

    def fake_download_repertoire(session, query, output_path):
        calls.append((session, query, output_path))

    monkeypatch.setattr(testdata, "_download_repertoire", fake_download_repertoire)

    testdata.download(args=None)

    assert len(calls) == 2

    alpha_session, alpha_query, alpha_output = calls[0]
    assert alpha_session is fake_session
    assert alpha_output == "alpharepertoire.airr.json"
    assert alpha_query == testdata._make_query("su008", "TRA")

    beta_session, beta_query, beta_output = calls[1]
    assert beta_session is fake_session
    assert beta_output == "betarepertoire.airr.json"
    assert beta_query == testdata._make_query("su008", "TRB", require_schema=True)

    legend_path = tmp_path / "legend.json"
    assert legend_path.exists()
    legend = json.loads(legend_path.read_text())
    assert legend == {
        "PRJNA509910-su008_pre-TRA": "Subject 8 pre-treatment",
        "PRJNA509910-su008_post-TRA": "Subject 8 post-treatment",
        "PRJNA509910-su008_pre-TRB": "Subject 8 pre-treatment",
        "PRJNA509910-su008_post-TRB": "Subject 8 post-treatment",
    }
