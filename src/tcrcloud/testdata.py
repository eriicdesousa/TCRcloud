"""Helpers to download and prepare example AIRR test data.

This module is used to download a pair of repertoire files (alpha/beta) from the
iReceptor AIRR API and to generate a small legend file used by the TCRcloud
example workflows.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import airr
import requests

from tcrcloud.download import get_session
from tcrcloud.errors import TCRcloudError

# Base URL for the iReceptor AIRR API.
# This specific node is hardcoded (rather than discovered, unlike download.py)
# because the TCRcloud example dataset (subject "su008" from Yost et al. 2019,
# study PRJNA509910) is known to live on this repository only.
HOST_URL = "https://ipa6.ireceptor.org/airr/v1"


def _make_query(
    subject_id: str, locus: str, require_schema: bool = False
) -> Mapping[str, Any]:
    """Build the AIRR search query payload.

    The query selects repertoires for a given subject and PCR locus (TRA or TRB).
    When ``require_schema`` is True, it additionally filters for repertoires that
    contain the AIRR rearrangement schema.

    The returned structure matches the iReceptor API filter syntax.
    """

    # Base filters required for every query.
    filters: list[Mapping[str, Any]] = [
        {
            "op": "contains",
            "content": {"field": "subject.subject_id", "value": subject_id},
        },
        {
            "op": "in",
            "content": {
                "field": "sample.pcr_target.pcr_target_locus",
                "value": [locus],
            },
        },
    ]

    # Optionally require the presence of the rearrangement schema.
    if require_schema:
        filters.append(
            {
                "op": "contains",
                "content": {
                    "field": "study.keywords_study",
                    "value": "contains_schema_rearrangement",
                },
            }
        )

    # Combine the filters using a top-level AND operation.
    return {"filters": {"op": "and", "content": filters}}


def _download_repertoire(
    session: requests.Session, query: Mapping[str, Any], output_path: str
) -> None:
    """Download a single repertoire and write it to disk using the AIRR library.

    ``session`` should be a retry-enabled `requests.Session` (see
    `tcrcloud.download.get_session`) so transient failures (rate limits,
    5xx errors) are retried automatically instead of failing the whole
    `testdata` command on the first hiccup.

    Any failure (network error, bad/non-JSON response, or an empty
    ``Repertoire`` list) is raised as a `tcrcloud.errors.TCRcloudError`,
    which the CLI (`tcrcloud.TCRcloud.main`) surfaces as a clean
    "TCRcloud error: ..." message with a non-zero exit code.
    """

    # Query iReceptor for the requested repertoire. `timeout` prevents the
    # CLI from hanging indefinitely if the server stops responding.
    try:
        resp = session.post(f"{HOST_URL}/repertoire", json=query, timeout=30)
        resp.raise_for_status()
    except requests.exceptions.RequestException as exc:
        raise TCRcloudError(f"could not reach {HOST_URL} ({exc})") from exc

    # Parse JSON response and write using the airr library.
    try:
        data = resp.json()
    except ValueError as exc:
        raise TCRcloudError(f"invalid JSON response from {HOST_URL} ({exc})") from exc

    # Guard against a "successful" response that has no actual data, which
    # would otherwise silently produce an empty repertoire file.
    repertoires = data.get("Repertoire") or []
    if not repertoires:
        raise TCRcloudError(
            f"no repertoires were returned by {HOST_URL} for {output_path}"
        )

    airr.write_repertoire(output_path, repertoires, info=data["Info"])

    print(f"Received {len(repertoires)} repertoires. Saved as {output_path}")


def download(args):
    """Download example test-data repertoires and generate a legend file.

    This is the entry point for the `TCRcloud testdata` CLI command (see
    the ``testdata`` subparser in `tcrcloud.TCRcloud`). ``args`` is accepted
    for consistency with the other subcommand entry points but is currently
    unused, since this command always fetches the same fixed example dataset.
    """

    # Reuse a single retry-enabled session for both requests so connections
    # are pooled and transient failures are retried consistently.
    session = get_session()

    # Download alpha and beta repertoires using pre-defined query parameters.
    # The TRA query is a basic locus filter; the TRB query also requires the
    # presence of the rearrangement schema so we get full AIRR rearrangement data.
    _download_repertoire(
        session, _make_query("su008", "TRA"), "alpharepertoire.airr.json"
    )
    _download_repertoire(
        session,
        _make_query("su008", "TRB", require_schema=True),
        "betarepertoire.airr.json",
    )

    # A small legend mapping the output identifiers used by the example pipeline
    # to human-friendly labels (used for plot legends etc.).
    legend = {
        "PRJNA509910-su008_pre-TRA": "Subject 8 pre-treatment",
        "PRJNA509910-su008_post-TRA": "Subject 8 post-treatment",
        "PRJNA509910-su008_pre-TRB": "Subject 8 pre-treatment",
        "PRJNA509910-su008_post-TRB": "Subject 8 post-treatment",
    }

    Path("legend.json").write_text(json.dumps(legend, indent=4))
    print("json file for legend saved as legend.json")
