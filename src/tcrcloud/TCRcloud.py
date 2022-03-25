#!/usr/bin/python3
import argparse
import yaml
import json

import matplotlib.pyplot as plt

import tcrcloud.cloud
import tcrcloud.radar
import tcrcloud.download
import tcrcloud.testdata

plt.rcParams["font.family"] = "serif"


def main():
    # create the top-level parser
    parser = argparse.ArgumentParser(
        description="Create a wordcloud of CDR3 \
    sequences from TCR AIRR-seq data or a radar plot with diversity metrics.",
        prog="TCRcloud")

    subparsers = parser.add_subparsers(
        title="Command Options",
        help="The program has 3 modes: cloud, radar or download",
        dest="cloud, radar, download, testdata or -h/--help",
        required=True)

    # create subparser for making the wordcloud
    parser_cloud = subparsers.add_parser("cloud", help="Create a wordcloud \
        from TCR CDR3 data")

    required_group = parser_cloud.add_argument_group("required arguments")

    # required_group.add_argument("-r","--repertoire",
    #     type= tcrcloud.base.jsonfile,
    #     help= "indicate the name of the AIRR Standards repertoire file",
    #     metavar= "repertoires.airr.json",
    #     required= True)

    required_group.add_argument("-a", "--rearrangements", type=str,
                                help="indicate the name of the AIRR Standards \
                                rearrangements file",
                                metavar="rearrangements.tsv", required=True)

    parser_cloud.set_defaults(func=tcrcloud.cloud.wordcloud)

    # create subparser for making the radar
    parser_radar = subparsers.add_parser("radar", help="Create a radar plot \
        with diversity metrics")
    required_group = parser_radar.add_argument_group("required arguments")
    required_group.add_argument("-a", "--rearrangements", type=str,
                                help="indicate the name of the AIRR \
                                Standards rearrangements file",
                                metavar="rearrangements.tsv", required=True)
    parser_radar.set_defaults(func=tcrcloud.radar.radar)

    # create subparser for downloading the rearregement data
    parser_download = subparsers.add_parser("download",
                                            help="Download TCR AIRR-seq \
                                            rearrangements data matching a \
                                            given repertoire metadata file")

    required_group = parser_download.add_argument_group("required arguments")

    required_group.add_argument("-r", "--repertoire", type=str,
                                help="indicate the name of the \
                                AIRR Standards repertoire file",
                                metavar="repertoires.airr.json",
                                required=True)

    parser_download.set_defaults(func=tcrcloud.download.airrdownload)

    # create subparser for downloading the test repertoire

    parser_testdata = subparsers.add_parser("testdata",
                                            help="Download TCR AIRR-seq \
                                            repertoire file to test TCRcloud")

    parser_testdata.set_defaults(func=tcrcloud.testdata.download)

    args = parser.parse_args()
    try:
        args.func(args)
    except FileNotFoundError:
        if dir(args)[-1] == "repertoire":
            print("TCRcloud error: " + args.repertoire
                  + " doesn't seem to exist")
        elif dir(args)[-1] == 'rearrangements':
            print("TCRcloud error: " + args.rearrangements
                  + " doesn't seem to exist")
    except (yaml.scanner.ScannerError, json.decoder.JSONDecodeError):
        print("TCRcloud error: It seems you did not indicate a properly \
formatted AIRR repertoire file")


if __name__ == "__main__":
    main()
