import requests
import airr


def testserver(data):
    # Find which server has the matching rearragements
    query = {
        "filters": {
            "op": "=",
            "content": {
                "field": "repertoire_id",
                "value": data["Repertoire"][0]["repertoire_id"]
            }
        }
    }
    repositories = ["https://vdjserver.org/airr/v1",
                    "https://ipa1.ireceptor.org/airr/v1",
                    "https://ipa2.ireceptor.org/airr/v1",
                    "https://ipa3.ireceptor.org/airr/v1",
                    "https://ipa4.ireceptor.org/airr/v1",
                    "http://ipa5.ireceptor.org/airr/v1",
                    "http://covid19-1.ireceptor.org/airr/v1",
                    "http://covid19-2.ireceptor.org/airr/v1",
                    "http://covid19-3.ireceptor.org/airr/v1",
                    "http://covid19-4.ireceptor.org/airr/v1",
                    "https://scireptor.dkfz.de/airr/v1",
                    "http://airr-seq.vdjbase.org/airr/v1"]

    host_url = "https://vdjserver.org/airr/v1"
    for i in repositories:
        test_url = i
        resp = requests.post(test_url + "/rearrangement", json=query)
        if len(resp.json()["Rearrangement"]) > 0:
            host_url = i
    return host_url


def airrdownload(args):
    try:
        airr.validate_repertoire(args.repertoire)
    except TypeError:
        print("TCRcloud error: unknown file type for repertoire. Supported \
file extensions are 'yaml', 'yml' or 'json'")
    repertoire_file = args.repertoire
    rearrangements_file = repertoire_file[:-4] + "rearrangements.tsv"
    data = airr.load_repertoire(args.repertoire)
    repertoires = data["Repertoire"]
    host_url = testserver(data)

    # Print out some Info
    print("       Info: " + data["Info"]["title"])
    print("    version: " + str(data["Info"]["version"]))
    print("description: " + data["Info"]["description"])
    print("Found " + str(len(data["Repertoire"])) + " repertoires in \
repertoire metadata file.")

    # Query the rearrangement endpoint
    # Define a generic query object, and we will replace the repertoire_id
    # within the loop. We also only request productive rearrangements as
    # an additional filter.

    query = {
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {
                        "field": "repertoire_id",
                        "value": "XXX"
                    }
                },
                {
                    "op": "=",
                    "content": {
                        "field": "productive",
                        "value": True
                    }
                }
            ]
        },
        "size": 1000,
        "from": 0
    }

    # Loop through each repertoire and query rearrangement data for
    # each. We download in chunks of 10000 because of the server
    # limitations using the from and size parameters.

    first = True
    for r in repertoires:
        print("Retrieving rearrangements for repertoire: "
              + r["repertoire_id"])
        query["filters"]["content"][0]["content"]["value"] = r["repertoire_id"]
        query["size"] = 1000
        query["from"] = 0

        cnt = 0
        while True:
            # send the request
            resp = requests.post(host_url + "/rearrangement", json=query)
            data = resp.json()
            rearrangements = data["Rearrangement"]

            # Open a file for writing the rearrangements. We do this here
            # because we need to know the full set of fields being
            # returned from the data repository, otherwise by default only
            # the required fields will be written to the file.
            if first:
                out_file = airr.create_rearrangement(
                    rearrangements_file,
                    fields=rearrangements[0].keys())
                first = False

            # save the rearrangements to a file
            for row in rearrangements:
                out_file.write(row)

            # looping until zero rearrangements are returned from the query.
            cnt += len(rearrangements)
            if len(rearrangements) < 1000:
                break

            # Need to update the from parameter to get the next chunk
            query["from"] = cnt

        print("Retrieved " + str(cnt) + " rearrangements for repertoire: "
                           + r["repertoire_id"])
        print("Saved as " + rearrangements_file)
