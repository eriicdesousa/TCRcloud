"""Unit tests for tcrcloud.download."""

import logging
from unittest.mock import Mock

import airr
import pytest
import requests

import tcrcloud.download as download
from tcrcloud.errors import TCRcloudError


def _response(json_payload=None, json_error=None, raise_for_status_error=None):
    """Build a fake `requests.Response`-like object."""

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


@pytest.fixture(autouse=True)
def _clear_discover_repositories_cache():
    """`discover_repositories` is `lru_cache`d; reset it around every test."""

    download.discover_repositories.cache_clear()
    yield
    download.discover_repositories.cache_clear()


# ---------------------------------------------------------------------------
# get_session
# ---------------------------------------------------------------------------


def test_get_session_configures_retries_for_http_and_https():
    session = download.get_session()

    assert isinstance(session, requests.Session)
    for prefix in ("http://", "https://"):
        adapter = session.get_adapter(prefix + "example.com")
        assert adapter.max_retries.total == 5
        assert 429 in adapter.max_retries.status_forcelist


# ---------------------------------------------------------------------------
# discover_repositories
# ---------------------------------------------------------------------------


def test_discover_repositories_collects_urls_from_list_payload(monkeypatch):
    monkeypatch.setenv(
        "TCRCLOUD_REPOSITORY_SEEDS",
        "https://seed-a.example/airr/v1,https://seed-b.example/airr/v1",
    )

    def fake_get(url, timeout=10):
        if "seed-a" in url:
            return _response(
                json_payload=[
                    {"url": "https://repo1.example/airr/v1"},
                    {"baseUrl": "https://repo2.example/airr/v1/"},
                ]
            )
        raise requests.exceptions.ConnectionError("unreachable")

    fake_session = Mock()
    fake_session.get.side_effect = fake_get
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    result = download.discover_repositories()

    assert result == [
        "https://repo1.example/airr/v1",
        "https://repo2.example/airr/v1",
    ]


def test_discover_repositories_handles_dict_payload_with_repositories_key(
    monkeypatch,
):
    monkeypatch.setenv("TCRCLOUD_REPOSITORY_SEEDS", "https://seed-a.example/airr/v1")

    fake_session = Mock()
    fake_session.get.return_value = _response(
        json_payload={"repositories": [{"base_url": "https://repo3.example/airr/v1"}]}
    )
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    result = download.discover_repositories()

    assert result == ["https://repo3.example/airr/v1"]


def test_discover_repositories_returns_empty_list_when_all_seeds_fail(monkeypatch):
    monkeypatch.setenv(
        "TCRCLOUD_REPOSITORY_SEEDS",
        "https://seed-a.example/airr/v1,https://seed-b.example/airr/v1",
    )

    fake_session = Mock()
    fake_session.get.side_effect = requests.exceptions.ConnectionError("unreachable")
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    assert download.discover_repositories() == []


def test_discover_repositories_ignores_invalid_json(monkeypatch):
    monkeypatch.setenv("TCRCLOUD_REPOSITORY_SEEDS", "https://seed-a.example/airr/v1")

    fake_session = Mock()
    fake_session.get.return_value = _response(json_error=ValueError("bad json"))
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    assert download.discover_repositories() == []


# ---------------------------------------------------------------------------
# testserver
# ---------------------------------------------------------------------------


def _repertoire_data(repertoire_id="rep1", study_id="study1"):
    return {
        "Repertoire": [
            {"repertoire_id": repertoire_id, "study": {"study_id": study_id}}
        ]
    }


def test_testserver_returns_first_matching_repository(monkeypatch, caplog):
    monkeypatch.setattr(
        download,
        "discover_repositories",
        Mock(
            return_value=[
                "https://repo1.example/airr/v1",
                "https://repo2.example/airr/v1",
            ]
        ),
    )

    fake_session = Mock()
    fake_session.post.return_value = _response(
        json_payload={"Repertoire": [{"repertoire_id": "rep1"}]}
    )
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    with caplog.at_level(logging.INFO, logger="tcrcloud.download"):
        host_url = download.testserver(_repertoire_data())

    assert host_url == "https://repo1.example/airr/v1"
    fake_session.post.assert_called_once()
    assert "Your repertoire was found" in caplog.text


def test_testserver_falls_back_to_next_repository_on_request_exception(
    monkeypatch, caplog
):
    monkeypatch.setattr(
        download,
        "discover_repositories",
        Mock(
            return_value=[
                "https://repo1.example/airr/v1",
                "https://repo2.example/airr/v1",
            ]
        ),
    )

    responses = [
        requests.exceptions.ConnectionError("unreachable"),
        _response(json_payload={"Repertoire": [{"repertoire_id": "rep1"}]}),
    ]

    def fake_post(url, json=None, timeout=10):
        result = responses.pop(0)
        if isinstance(result, Exception):
            raise result
        return result

    fake_session = Mock()
    fake_session.post.side_effect = fake_post
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    host_url = download.testserver(_repertoire_data())

    assert host_url == "https://repo2.example/airr/v1"
    assert "could not reach" in caplog.text


def test_testserver_skips_repository_with_invalid_json(monkeypatch, caplog):
    monkeypatch.setattr(
        download,
        "discover_repositories",
        Mock(
            return_value=[
                "https://repo1.example/airr/v1",
                "https://repo2.example/airr/v1",
            ]
        ),
    )

    fake_session = Mock()
    fake_session.post.side_effect = [
        _response(json_error=ValueError("bad json")),
        _response(json_payload={"Repertoire": [{"repertoire_id": "rep1"}]}),
    ]
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    host_url = download.testserver(_repertoire_data())

    assert host_url == "https://repo2.example/airr/v1"
    assert "invalid JSON" in caplog.text


def test_testserver_exits_when_repertoire_not_found_anywhere(monkeypatch):
    monkeypatch.setattr(
        download,
        "discover_repositories",
        Mock(return_value=["https://repo1.example/airr/v1"]),
    )

    fake_session = Mock()
    fake_session.post.return_value = _response(json_payload={"Repertoire": []})
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    with pytest.raises(TCRcloudError, match="could not reach any"):
        download.testserver(_repertoire_data())


# ---------------------------------------------------------------------------
# airrdownload
# ---------------------------------------------------------------------------


def _airr_metadata():
    return {
        "Repertoire": [{"repertoire_id": "rep1"}],
        "Info": {"title": "Test Study", "version": "1.0", "description": "desc"},
    }


def _patch_common_airrdownload_deps(monkeypatch, read_airr_return=None):
    """Patch the parts of airrdownload unrelated to the pagination behaviour
    being tested: file validation, metadata loading and repository lookup."""

    monkeypatch.setattr(download.airr, "validate_repertoire", Mock())
    monkeypatch.setattr(
        download.airr,
        "read_airr",
        Mock(return_value=read_airr_return or _airr_metadata()),
    )
    monkeypatch.setattr(
        download, "testserver", Mock(return_value="https://fakehost.example/airr/v1")
    )


def test_airrdownload_writes_rearrangements_and_closes_file(monkeypatch, caplog):
    _patch_common_airrdownload_deps(monkeypatch)

    rows = [{"junction_aa": "CASSX", "v_call": "TRBV1"}]
    fake_session = Mock()
    fake_session.post.return_value = _response(json_payload={"Rearrangement": rows})
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    fake_writer = Mock()
    create_rearrangement = Mock(return_value=fake_writer)
    monkeypatch.setattr(download.airr, "create_rearrangement", create_rearrangement)

    with caplog.at_level(logging.INFO, logger="tcrcloud.download"):
        download.airrdownload("sample.airr.json")

    create_rearrangement.assert_called_once_with(
        "sample.airr.rearrangements.tsv", fields=rows[0].keys()
    )
    fake_writer.write.assert_called_once_with(rows[0])
    fake_writer.close.assert_called_once()
    assert "Saved as sample.airr.rearrangements.tsv" in caplog.text


def test_airrdownload_handles_empty_rearrangements_without_crashing(
    monkeypatch, caplog
):
    """A repertoire with zero productive rearrangements must not raise
    IndexError, and no file should be created (regression test)."""

    _patch_common_airrdownload_deps(monkeypatch)

    fake_session = Mock()
    fake_session.post.return_value = _response(json_payload={"Rearrangement": []})
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    create_rearrangement = Mock()
    monkeypatch.setattr(download.airr, "create_rearrangement", create_rearrangement)

    download.airrdownload("sample.airr.json")  # should not raise

    create_rearrangement.assert_not_called()
    assert "no rearrangements were found" in caplog.text


def test_airrdownload_paginates_until_short_page(monkeypatch):
    _patch_common_airrdownload_deps(monkeypatch)

    first_page = [{"junction_aa": f"CASS{i}"} for i in range(1000)]
    second_page = [{"junction_aa": "LAST"}]

    fake_session = Mock()
    fake_session.post.side_effect = [
        _response(json_payload={"Rearrangement": first_page}),
        _response(json_payload={"Rearrangement": second_page}),
    ]
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    fake_writer = Mock()
    monkeypatch.setattr(
        download.airr, "create_rearrangement", Mock(return_value=fake_writer)
    )

    download.airrdownload("sample.airr.json")

    assert fake_session.post.call_count == 2
    # The second request should page using `from` = number fetched so far.
    second_call_kwargs = fake_session.post.call_args_list[1].kwargs
    assert second_call_kwargs["json"]["from"] == 1000
    assert fake_writer.write.call_count == len(first_page) + len(second_page)
    fake_writer.close.assert_called_once()


def test_airrdownload_invalid_repertoire_file_exits(monkeypatch):
    monkeypatch.setattr(download.airr, "validate_repertoire", Mock())
    monkeypatch.setattr(download.airr, "read_airr", Mock(side_effect=TypeError()))

    with pytest.raises(TCRcloudError, match="properly formatted AIRR repertoire file"):
        download.airrdownload("sample.airr.json")


def test_airrdownload_validation_error_is_converted_to_tcrcloud_error(monkeypatch):
    """A file failing airr.validate_repertoire must surface as a
    TCRcloudError naming the file, not as an internal ValidationError that
    the CLI would otherwise misreport."""

    monkeypatch.setattr(
        download.airr,
        "validate_repertoire",
        Mock(side_effect=airr.ValidationError("bad schema")),
    )

    with pytest.raises(
        TCRcloudError, match="not a valid AIRR repertoire metadata file"
    ):
        download.airrdownload("sample.airr.json")


def test_airrdownload_missing_repertoire_key_raises_tcrcloud_error(monkeypatch):
    """A metadata dict without the top-level "Repertoire" key must be
    reported as a malformed input file, not surface as a bare KeyError."""

    _patch_common_airrdownload_deps(
        monkeypatch,
        read_airr_return={"Info": {"title": "t", "version": 1, "description": "d"}},
    )

    with pytest.raises(TCRcloudError, match="properly formatted AIRR repertoire file"):
        download.airrdownload("sample.airr.json")


def test_airrdownload_raises_on_request_exception(monkeypatch):
    _patch_common_airrdownload_deps(monkeypatch)

    fake_session = Mock()
    fake_session.post.side_effect = requests.exceptions.ConnectionError("boom")
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    create_rearrangement = Mock()
    monkeypatch.setattr(download.airr, "create_rearrangement", create_rearrangement)

    with pytest.raises(TCRcloudError, match="failed to download rearrangements"):
        download.airrdownload("sample.airr.json")

    create_rearrangement.assert_not_called()


def test_airrdownload_raises_on_invalid_json(monkeypatch):
    _patch_common_airrdownload_deps(monkeypatch)

    fake_session = Mock()
    fake_session.post.return_value = _response(json_error=ValueError("bad json"))
    monkeypatch.setattr(download, "get_session", Mock(return_value=fake_session))

    create_rearrangement = Mock()
    monkeypatch.setattr(download.airr, "create_rearrangement", create_rearrangement)

    with pytest.raises(TCRcloudError, match="invalid JSON response"):
        download.airrdownload("sample.airr.json")

    create_rearrangement.assert_not_called()
