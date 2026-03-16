"""Helpers to download and prepare example AIRR test data.

This module is used to download a pair of repertoire files (alpha/beta) from the
iReceptor AIRR API and to generate a small legend file used by the TCRcloud
example workflows.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import airr
import requests


# Base URL for the iReceptor AIRR API.
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


def _download_repertoire(query: Mapping[str, Any], output_path: str) -> None:
    """Download a repertoire and write it to disk using the AIRR library."""

    # Query iReceptor for the requested repertoire.
    resp = requests.post(f"{HOST_URL}/repertoire", json=query)
    resp.raise_for_status()

    # Parse JSON response and write using the airr library.
    data = resp.json()
    airr.write_repertoire(output_path, data["Repertoire"], info=data["Info"])

    print(f"Received {len(data['Repertoire'])} repertoires. Saved as {output_path}")


def download(args):
    """Download example test-data repertoires and generate a legend file."""

    # Download alpha and beta repertoires using pre-defined query parameters.
    # The TRA query is a basic locus filter; the TRB query also requires the
    # presence of the rearrangement schema so we get full AIRR rearrangement data.
    _download_repertoire(_make_query("su008", "TRA"), "alpharepertoire.airr.json")
    _download_repertoire(
        _make_query("su008", "TRB", require_schema=True), "betarepertoire.airr.json"
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
