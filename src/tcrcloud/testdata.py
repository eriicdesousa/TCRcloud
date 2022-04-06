import requests
import airr


def download(args):
    # This study is stored at VDJServer data repository
    host_url = "https://vdjserver.org/airr/v1"

    # POST data is sent with the query. Here we construct an object for
    # the query ((study_id == "PRJNA300878") AND (locus == "TRB"))

    query = {
        "filters": {
            "op": "and",
            "content": [
                {
                    "op": "=",
                    "content": {
                        "field": "study.study_id",
                        "value": "PRJNA300878"
                    }
                },
                {
                    "op": "=",
                    "content": {
                        "field": "sample.pcr_target.pcr_target_locus",
                        "value": "TRB"
                    }
                },
                {
                    "op": "or",
                    "content": [
                        {
                            "op": "=",
                            "content": {
                                "field": "sample.sample_id",
                                "value": "TW02B_T_memory_CD8"
                            }
                        },
                        {
                            "op": "=",
                            "content": {
                                "field": "sample.sample_id",
                                "value": "TW02A_T_memory_CD8"
                            }
                        }
                    ]
                }
            ]
        }
    }

    # Send the query
    resp = requests.post(host_url + "/repertoire", json=query)

    # The data is returned as JSON, use AIRR library to write out data
    data = resp.json()
    airr.write_repertoire("testdata.airr.json",
                          data["Repertoire"], info=data["Info"])

    # Print out some Info
    print("       Info: " + data["Info"]["title"])
    print("    version: " + str(data["Info"]["version"]))
    print("description: " + data["Info"]["description"])

    # Save them out to a file
    print("Received " + str(len(data["Repertoire"]))
          + " repertoires. Saved as testdata.airr.json")
    with open("legendfortestdata.json", "w") as fileout:
        print("{", file=fileout)
        print('    "2839362682105696746-242ac113-0001-012":"Twin 2A",',
              file=fileout)
        print('    "2939134772391776746-242ac113-0001-012":"Twin 2B"',
              file=fileout)
        print("}", file=fileout)
    print("json file for legend saved as testdata.json")
