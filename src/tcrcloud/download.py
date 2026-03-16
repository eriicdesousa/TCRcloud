"""Download AIRR rearrangements for a local repertoire metadata file.

This module is invoked by the `TCRcloud download` CLI command.

Workflow:
1. Read the local AIRR repertoire metadata file.
2. Discover available AIRR Data Commons repositories via `/repositories`.
3. Find the repository that contains the requested repertoire via `/repertoire` search.
4. Page through `/rearrangement` results and write them to a TSV file.

If repository discovery fails, a built-in fallback list of known AIRR nodes is used.
"""

# Standard library imports
import os
import sys
import time
from functools import lru_cache

# Third-party imports
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

import requests
import airr


def get_session() -> requests.Session:
    """Return a requests.Session configured with sensible default retries.

    This session is shared across the whole download process to:
    - reuse HTTP connections for better performance
    - automatically retry transient failures (rate limits, 5xx, etc.)
    - avoid flooding the server with repeated immediate retries
    """

    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "POST"),
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


# Repository discovery
# --------------------
# The AIRR Data Commons defines a standard `/repositories` endpoint that
# lists all known AIRR repository instances. We query a small set of well-known
# seed servers and gather their listings, so we don't need a hard-coded repo list.


@lru_cache(maxsize=1)
def discover_repositories() -> list[str]:
    """Discover AIRR Data Commons repository URLs via /repositories endpoints.

    Queries a small set of known AIRR hosts for their `/repositories`
    registry listing and returns the base URLs found.

    If discovery fails (e.g., network/DNS errors), this returns an empty list.
    """

    session = get_session()

    # Seeds come from a small set of reliably-maintained AIRR commons.
    # Can be overridden with the TCRCLOUD_REPOSITORY_SEEDS env var.
    # Example: export TCRCLOUD_REPOSITORY_SEEDS="https://vdjserver.org/airr/v1,https://ipa1.ireceptor.org/airr/v1"
    seeds = os.environ.get(
        "TCRCLOUD_REPOSITORY_SEEDS",
        "https://vdjserver.org/airr/v1,https://ipa1.ireceptor.org/airr/v1,https://airr-seq.vdjbase.org/airr/v1",
    ).split(",")

    found: set[str] = set()
    # Query each seed server for their `/repositories` response.
    # The first successful response contributes URLs to the candidate set.
    for seed in seeds:
        seed = seed.strip()
        if not seed:
            continue
        try:
            resp = session.get(seed.rstrip("/") + "/repositories", timeout=10)
            resp.raise_for_status()
        except requests.exceptions.RequestException:
            # Skip seeds that cannot be reached or return errors.
            continue

        try:
            payload = resp.json()
        except ValueError:
            continue

        # Depending on the server, `/repositories` may return a list or a dict.
        candidates = (
            payload
            if isinstance(payload, list)
            else payload.get("repositories") or payload.get("Repositories")
        )
        if not candidates:
            continue

        for entry in candidates:
            if not isinstance(entry, dict):
                continue
            url = (
                entry.get("url")
                or entry.get("baseUrl")
                or entry.get("base_url")
                or entry.get("baseURL")
            )
            if url:
                found.add(url.rstrip("/"))

    # Return a stable (sorted) list for deterministic behavior.
    return sorted(found)


# Fallback repository list
# ------------------------
# If discovery fails (e.g., if `/repositories` is unsupported), we fall back to a hard-coded list of known nodes.
FALLBACK_REPOSITORIES = [
    "https://vdjserver.org/airr/v1",
    "https://ipa1.ireceptor.org/airr/v1",
    "https://ipa2.ireceptor.org/airr/v1",
    "https://ipa3.ireceptor.org/airr/v1",
    "https://ipa4.ireceptor.org/airr/v1",
    "https://ipa5.ireceptor.org/airr/v1",
    "https://ipa6.ireceptor.org/airr/v1",
    "https://covid19-1.ireceptor.org/airr/v1",
    "https://covid19-2.ireceptor.org/airr/v1",
    "https://covid19-3.ireceptor.org/airr/v1",
    "https://covid19-4.ireceptor.org/airr/v1",
    "https://scireptor.dkfz.de/airr/v1",
    "https://airr-seq.vdjbase.org/airr/v1",
    "https://agschwab.uni-muenster.de/airr/v1",
    "https://roche-airr.ireceptor.org/airr/v1",
    "https://t1d-1.ireceptor.org/airr/v1",
    "https://t1d-2.ireceptor.org/airr/v1",
    "https://t1d-3.ireceptor.org/airr/v1",
    "https://hpap.ireceptor.org/airr/v1/",
    "https://greifflab-1.ireceptor.org/airr/v1",
]


def testserver(data):
    # Determine which AIRR repository holds the requested repertoire.
    # We do this by issuing a simple `POST /repertoire` search for the
    # repertoire_id and study_id from the metadata file.
    query = {
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {
                        "field": "repertoire_id",
                        "value": data["Repertoire"][0]["repertoire_id"],
                    },
                },
                {
                    "op": "=",
                    "content": {
                        "field": "study.study_id",
                        "value": data["Repertoire"][0]["study"]["study_id"],
                    },
                },
            ],
        }
    }

    session = get_session()
    print("Looking for your repertoire in the AIRR Data Commons repositories...")
    repositories = discover_repositories() or FALLBACK_REPOSITORIES

    host_url = None
    for i in repositories:
        test_url = i
        try:
            resp = session.post(
                test_url.rstrip("/") + "/repertoire", json=query, timeout=10
            )
            resp.raise_for_status()
        except requests.exceptions.RequestException:
            sys.stderr.write(
                f"TCRcloud warning: could not reach {test_url}. Trying the next repository.\n"
            )
            continue

        try:
            repertoires = resp.json().get("Repertoire", [])
        except ValueError as e:
            sys.stderr.write(
                f"TCRcloud warning: invalid JSON from {test_url} ({e}). Trying the next repository.\n"
            )
            continue

        if len(repertoires) > 0:
            host_url = i
            print("Your repertoire was found at " + host_url)
            break

    if host_url is None:
        sys.stderr.write(
            "TCRcloud error: could not reach any of the configured AIRR repositories. "
            "Please check your network connection and try again.\n"
        )
        sys.exit(1)

    return host_url


def airrdownload(args):
    # Validate the input file path and load the AIRR repertoire metadata.
    airr.validate_repertoire(args.repertoire)
    repertoire_file = args.repertoire
    rearrangements_file = repertoire_file[:-4] + "rearrangements.tsv"

    try:
        data = airr.read_airr(args.repertoire)
    except TypeError:
        sys.stderr.write(
            "TCRcloud error: It seems you did not indicate a \
properly formatted AIRR repertoire file\n"
        )
        exit()

    # The metadata file may contain multiple repertoires; we will iterate them.
    repertoires = data["Repertoire"]

    # Find the correct AIRR repository that holds this repertoire.
    host_url = testserver(data)

    # Use a shared HTTP session for all subsequent requests.
    session = get_session()

    print(f"Using AIRR repository: {host_url}")

    # Retrieve repository-level metadata from the AIRR file.
    # Some versions nest the metadata under `Info.Info`.
    try:
        info = data["Info"]["Info"]
    except KeyError:
        info = data["Info"]

    # Print out basic metadata for user visibility.
    print("       Info: " + info["title"])
    print("    version: " + str(info["version"]))
    print("description: " + info["description"])
    print(
        "Found "
        + str(len(data["Repertoire"]))
        + " repertoires in \
repertoire metadata file."
    )

    # Build a reusable query template for the rearrangements endpoint.
    # We will replace `repertoire_id` inside the loop for each repertoire.
    # We also filter for only productive rearrangements.
    query = {
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {"field": "repertoire_id", "value": "random_value"},
                },
                {"op": "=", "content": {"field": "productive", "value": True}},
            ],
        },
        "size": 1000,
        "from": 0,
    }

    # Loop through each repertoire and download rearrangements in pages.
    # Many servers limit how many records can be returned per request, so
    # we page through results using `from` and `size`.

    first = True
    for idx, r in enumerate(repertoires, start=1):
        print(
            f"Retrieving rearrangements for repertoire {idx}/{len(repertoires)}: {r['repertoire_id']}"
        )
        print(
            "This process may take some time depending on the number of "
            "rearrangements you are downloading"
        )
        query["filters"]["content"][0]["content"]["value"] = r["repertoire_id"]
        query["size"] = 1000
        query["from"] = 0

        cnt = 0
        start_time = time.perf_counter()
        while True:
            # send the request
            try:
                resp = session.post(
                    host_url.rstrip("/") + "/rearrangement", json=query, timeout=30
                )
                resp.raise_for_status()
            except requests.exceptions.RequestException as e:
                sys.stderr.write(
                    f"TCRcloud error: failed to download rearrangements from {host_url} ({e})\n"
                )
                return

            try:
                data = resp.json()
                rearrangements = data.get("Rearrangement", [])
            except ValueError as e:
                sys.stderr.write(
                    f"TCRcloud error: invalid JSON response from {host_url} ({e})\n"
                )
                return

            # progress report
            cnt += len(rearrangements)
            elapsed = time.perf_counter() - start_time
            print(
                f"  fetched {cnt} rows (last chunk {len(rearrangements)}) "
                f"in {elapsed:.1f}s",
                end="\r",
                flush=True,
            )

            # Open a file for writing the rearrangements. We do this here
            # because we need to know the full set of fields being
            # returned from the data repository, otherwise by default only
            # the required fields will be written to the file.
            if first:
                out_file = airr.create_rearrangement(
                    rearrangements_file, fields=rearrangements[0].keys()
                )
                first = False

            # save the rearrangements to a file
            for row in rearrangements:
                out_file.write(row)

            # looping until zero rearrangements are returned from the query.
            if len(rearrangements) < 1000:
                break

            # Need to update the from parameter to get the next chunk
            query["from"] = cnt

        print()
        print(
            "Retrieved "
            + str(cnt)
            + " rearrangements for repertoire: "
            + r["repertoire_id"]
        )
    print("Saved as " + rearrangements_file)
