#!/usr/bin/env python3
import argparse
import yaml
import json
import sys

import matplotlib.pyplot as plt

import tcrcloud.cloud
import tcrcloud.radar
import tcrcloud.download
import tcrcloud.testdata

plt.rcParams["font.family"] = "serif"


def main():
    # create the top-level parser
    parser = argparse.ArgumentParser(
        description="Create a word cloud of CDR3 \
    sequences from TCR AIRR-seq data or a radar plot with diversity metrics.",
        prog="TCRcloud")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s 1.4.1")
    subparsers = parser.add_subparsers(
        title="command options",
        help="The program has 4 options: cloud, radar, download or testdata",
        dest="cloud, radar, download or testdata",
        required=True)

    # create subparser for making the wordcloud
    parser_cloud = subparsers.add_parser("cloud", help="Create a wordcloud \
        from TCR CDR3 data")

    # required_group = parser_cloud.add_argument_group("arguments")

    # required_group.add_argument("-r","--repertoire",
    #     type= tcrcloud.base.jsonfile,
    #     help= "indicate the name of the AIRR Standards repertoire file",
    #     metavar= "repertoires.airr.json",
    #     required= True)

    parser_cloud.add_argument("rearrangements", type=str,
                              help="indicate the name of the AIRR Standards \
                              rearrangements file",
                              metavar="rearrangements.tsv")
    parser_cloud.add_argument("-c", "--colours", type=str,
                              help="indicate the name of a json file to \
                              change the colours of the wordcloud",
                              metavar="colours.json", required=False)
    parser_cloud.add_argument("-l", "--legend", type=str,
                              help="indicate if legend should be included, \
                              default = True",
                              metavar="True or False", default="True",
                              required=False)
    parser_cloud.add_argument("-s", "--size", type=int,
                              help=argparse.SUPPRESS,
                              default=1000,
                              required=False)

    parser_cloud.set_defaults(func=tcrcloud.cloud.wordcloud)

    # create subparser for making the radar
    parser_radar = subparsers.add_parser("radar", help="Create a radar plot \
        with diversity metrics")

    parser_radar.add_argument("rearrangements", type=str,
                              help="indicate the name of the AIRR \
                              Standards rearrangements file",
                              metavar="rearrangements.tsv")
    parser_radar.add_argument("-c", "--custom_legend", type=str,
                              help="indicate the name of a json \
                              file to convert repertoire_id to what you \
                              want to appear in the legend",
                              metavar="legend.json", required=False)
    parser_radar.add_argument("-l", "--legend", type=str,
                              help="indicate if legend should be included, \
                              default = True",
                              metavar="True or False", default="True",
                              required=False)
    parser_radar.add_argument("-e", "--export", type=str,
                              help="indicate if the metrics from the radar \
                              should be exported to a text file, \
                              default = False",
                              metavar="True or False", default="False",
                              required=False)
    parser_radar.set_defaults(func=tcrcloud.radar.radar)

    # create subparser for downloading the rearregement data
    parser_download = subparsers.add_parser("download",
                                            help="Download TCR AIRR-seq \
                                            rearrangements data matching a \
                                            given repertoire metadata file")

    parser_download.add_argument("repertoire", type=str,
                                 help="indicate the name of the \
                                 AIRR Standards repertoire file",
                                 metavar="repertoires.airr.json")

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
        if dir(args)[-2] == "repertoire":
            sys.stderr.write("TCRcloud error: " + args.repertoire
                             + " doesn't seem to exist\n")
        elif dir(args)[-2] == "rearrangements":
            sys.stderr.write("TCRcloud error: " + args.rearrangements
                             + " doesn't seem to exist\n")
    except (yaml.scanner.ScannerError, json.decoder.JSONDecodeError):
        sys.stderr.write("TCRcloud error: It seems you did not indicate a \
properly formatted AIRR repertoire file\n")
    except (KeyError, TypeError):
        sys.stderr.write("TCRcloud error: It seems you did not indicate a \
properly formatted AIRR rearrangements file\n")


if __name__ == "__main__":
    main()
