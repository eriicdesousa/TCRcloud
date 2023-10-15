import requests
import airr


def download(args):
    host_url = " https://ipa6.ireceptor.org/airr/v1"

    query1 = {
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "contains",
                    "content": {
                        "field": "subject.subject_id",
                        "value": "su008"
                    }
                },
                {
                    "op": "in",
                    "content": {
                        "field": "sample.pcr_target.pcr_target_locus",
                        "value": ["TRA"]
                    }
                }
            ]
        }
    }

    query2 = {
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "contains",
                    "content": {
                        "field": "subject.subject_id",
                        "value": "su008"
                    }
                },
                {
                    "op": "in",
                    "content": {
                        "field": "sample.pcr_target.pcr_target_locus",
                        "value": ["TRB"]
                    }
                },
                {
                    "op": "contains",
                    "content": {
                        "field": "study.keywords_study",
                        "value": "contains_schema_rearrangement"
                    }
                }
            ]
        }
    }
    # Send the query
    resp = requests.post(host_url + "/repertoire", json=query1)

    # The data is returned as JSON, use AIRR library to write out data
    data = resp.json()
    airr.write_repertoire("alpharepertoire.airr.json",
                          data["Repertoire"], info=data["Info"])

    # Print out some Info
    print("       Info: " + data["Info"]["Info"]["title"])
    print("    version: " + str(data["Info"]["Info"]["version"]))
    print("description: " + data["Info"]["Info"]["description"])

    # Save them out to a file
    print("Received " + str(len(data["Repertoire"]))
          + " alpha repertoires. Saved as alpharepertoire.airr.json")

    # Send the query
    resp = requests.post(host_url + "/repertoire", json=query2)

    # The data is returned as JSON, use AIRR library to write out data
    data = resp.json()
    airr.write_repertoire("betarepertoire.airr.json",
                          data["Repertoire"], info=data["Info"])

    # Save them out to a file
    print("Received " + str(len(data["Repertoire"]))
          + " beta repertoires. Saved as betarepertoire.airr.json")

    print()

    with open("legend.json", "w") as fileout:
        print("{", file=fileout)
        print('    "PRJNA509910-su008_pre-TRA":"Subject 8 pre-treatment",',
              file=fileout)
        print('    "PRJNA509910-su008_post-TRA":"Subject 8 post-treatment",',
              file=fileout)
        print('    "PRJNA509910-su008_pre-TRB":"Subject 8 pre-treatment",',
              file=fileout)
        print('    "PRJNA509910-su008_post-TRB":"Subject 8 post-treatment"',
              file=fileout)
        print("}", file=fileout)
    print("json file for legend saved as legend.json")
