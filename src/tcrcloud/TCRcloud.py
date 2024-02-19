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
import tcrcloud.vgenes
import tcrcloud.aminoacids

plt.rcParams["font.family"] = "serif"


def main():
    # create the top-level parser
    parser = argparse.ArgumentParser(
        description="Create visualizations from AIRR-seq data", prog="TCRcloud"
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.4.1")
    subparsers = parser.add_subparsers(
        title="command options",
        help="The program has 6 options: cloud, radar, vgenes, aminoacids, \
             download or testdata",
        dest="cloud, radar, vgenes, aminoacids, download or testdata",
        required=True,
    )

    # create subparser for making the word cloud
    parser_cloud = subparsers.add_parser(
        "cloud",
        help="Create a word cloud \
        from AIRR CDR3 data",
    )

    # required_group = parser_cloud.add_argument_group("arguments")

    # required_group.add_argument("-r","--repertoire",
    #     type= tcrcloud.base.jsonfile,
    #     help= "indicate the name of the AIRR Standards repertoire file",
    #     metavar= "repertoires.airr.json",
    #     required= True)

    parser_cloud.add_argument(
        "rearrangements",
        type=str,
        help="indicate the name of the AIRR \
                              rearrangements file",
        metavar="rearrangements.tsv",
    )
    parser_cloud.add_argument(
        "-c",
        "--colours",
        type=str,
        help="indicate the name of a json file to \
                              change the colours of the word cloud",
        metavar="colours.json",
        required=False,
    )
    parser_cloud.add_argument(
        "-l",
        "--legend",
        type=str,
        help="indicate if legend should be included, \
                              default = True",
        metavar="True or False",
        default="True",
        required=False,
    )
    parser_cloud.add_argument(
        "-s", "--size", type=int, help=argparse.SUPPRESS, default=1000, required=False
    )

    parser_cloud.set_defaults(func=tcrcloud.cloud.wordcloud)

    # create subparser for making the radar
    parser_radar = subparsers.add_parser(
        "radar",
        help="Create a radar plot \
        with diversity metrics from AIRR CDR3 data",
    )

    parser_radar.add_argument(
        "rearrangements",
        type=str,
        help="indicate the name of the AIRR \
                              rearrangements file",
        metavar="rearrangements.tsv",
    )
    parser_radar.add_argument(
        "-c",
        "--custom_legend",
        type=str,
        help="indicate the name of a json \
                              file to convert repertoire_id to what you \
                              want to appear in the legend",
        metavar="legend.json",
        required=False,
    )
    parser_radar.add_argument(
        "-l",
        "--legend",
        type=str,
        help="indicate if legend should be included, \
                              default = True",
        metavar="True or False",
        default="True",
        required=False,
    )
    parser_radar.add_argument(
        "-e",
        "--export",
        type=str,
        help="indicate if the metrics from the radar \
                              should be exported to a text file, \
                              default = False",
        metavar="True or False",
        default="False",
        required=False,
    )
    parser_radar.add_argument(
        "-fh",
        "--fmax",
        type=float,
        help="define the max value for the D50 \
                              Index axis, default = 50",
        metavar="float",
        default=50,
        required=False,
    )
    parser_radar.add_argument(
        "-fl",
        "--fmin",
        type=float,
        help="define the min value for the D50 \
                              Index axis, default = 0",
        metavar="float",
        default=0,
        required=False,
    )
    parser_radar.add_argument(
        "-ch",
        "--cmax",
        type=float,
        help="define the max value for the Convergence \
                              axis, default = 1",
        metavar="float",
        default=1,
        required=False,
    )
    parser_radar.add_argument(
        "-cl",
        "--cmin",
        type=float,
        help="define the min value for the Convergence \
                              axis, default = 0",
        metavar="float",
        default=0,
        required=False,
    )
    parser_radar.add_argument(
        "-chh",
        "--chmax",
        type=int,
        help="define the max value for the Chao1 \
                              Index axis, default = 17000",
        metavar="int",
        default=17000,
        required=False,
    )
    parser_radar.add_argument(
        "-chl",
        "--chmin",
        type=int,
        help="define the min value for the Chao1 \
                              Index axis, default = 0",
        metavar="int",
        default=0,
        required=False,
    )
    parser_radar.add_argument(
        "-uh",
        "--umax",
        type=int,
        help="define the max value for the Distinct \
                              CDR3 axis, default = 10000",
        metavar="int",
        default=10000,
        required=False,
    )
    parser_radar.add_argument(
        "-ul",
        "--umin",
        type=int,
        help="define the min value for the Distinct \
                              CDR3 axis, default = 0",
        metavar="int",
        default=0,
        required=False,
    )
    parser_radar.add_argument(
        "-gsh",
        "--gsmax",
        type=float,
        help="define the max value for the Gini-Simpson \
                              Index axis, default = 1",
        metavar="float",
        default=1,
        required=False,
    )
    parser_radar.add_argument(
        "-gsl",
        "--gsmin",
        type=float,
        help="define the min value for the Gini-Simpson \
                              Index axis, default = 0.90",
        metavar="float",
        default=0.90,
        required=False,
    )
    parser_radar.add_argument(
        "-sh",
        "--smax",
        type=float,
        help="define the max value for the Shannon \
                              Index axis, default = 15",
        metavar="float",
        default=15,
        required=False,
    )
    parser_radar.add_argument(
        "-sl",
        "--smin",
        type=float,
        help="define the min value for the Shannon \
                              Index axis, default = 0",
        metavar="float",
        default=0,
        required=False,
    )
    parser_radar.add_argument(
        "-gh",
        "--gmax",
        type=float,
        help="define the max value for the Gini \
                              Coefficient axis, default = 1",
        metavar="float",
        default=1,
        required=False,
    )
    parser_radar.add_argument(
        "-gl",
        "--gmin",
        type=float,
        help="define the min value for the Gini \
                              Coefficent axis, default = 0",
        metavar="float",
        default=0,
        required=False,
    )

    parser_radar.set_defaults(func=tcrcloud.radar.radar)

    # create subparser for making the V gene plot
    parser_vgenes = subparsers.add_parser(
        "vgenes",
        help="Create a V gene surface plot \
        from AIRR CDR3 data",
    )

    parser_vgenes.add_argument(
        "rearrangements",
        type=str,
        help="indicate the name of the AIRR \
                              rearrangements file",
        metavar="rearrangements.tsv",
    )
    parser_vgenes.add_argument(
        "-e",
        "--export",
        type=str,
        help="indicate if the metrics from the plot \
                              should be exported to a csv file, \
                              default = False",
        metavar="True or False",
        default="False",
        required=False,
    )
    parser_vgenes.add_argument(
        "-yha",
        "--yhighalpha",
        type=int,
        help="indicate the max value for the y axis \
                              for alpha chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-yla",
        "--ylowalpha",
        type=int,
        help="indicate the min value for the y axis \
                              for alpha chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zha",
        "--zhighalpha",
        type=float,
        help="indicate the max value for the z axis \
                              for alpha chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zla",
        "--zlowalpha",
        type=float,
        help="indicate the min value for the z axis \
                              for alpha chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-yhb",
        "--yhighbeta",
        type=int,
        help="indicate the max value for the y axis \
                              for beta chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-ylb",
        "--ylowbeta",
        type=int,
        help="indicate the min value for the y axis \
                              for beta chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zhb",
        "--zhighbeta",
        type=float,
        help="indicate the max value for the z axis \
                              for beta chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zlb",
        "--zlowbeta",
        type=float,
        help="indicate the min value for the z axis \
                              for beta chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-yhg",
        "--yhighgamma",
        type=int,
        help="indicate the max value for the y axis \
                              for gamma chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-ylg",
        "--ylowgamma",
        type=int,
        help="indicate the min value for the y axis \
                              for gamma chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zhg",
        "--zhighgamma",
        type=float,
        help="indicate the max value for the z axis \
                              for gamma chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zlg",
        "--zlowgamma",
        type=float,
        help="indicate the min value for the z axis \
                              for gamma chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-yhd",
        "--yhighdelta",
        type=int,
        help="indicate the max value for the y axis \
                              for delta chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-yld",
        "--ylowdelta",
        type=int,
        help="indicate the min value for the y axis \
                              for delta chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zhd",
        "--zhighdelta",
        type=float,
        help="indicate the max value for the z axis \
                              for delta chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zld",
        "--zlowdelta",
        type=float,
        help="indicate the min value for the z axis \
                              for delta chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-yhh",
        "--yhighheavy",
        type=int,
        help="indicate the max value for the y axis \
                              for heavy chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-ylh",
        "--ylowheavy",
        type=int,
        help="indicate the min value for the y axis \
                              for heavy chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zhh",
        "--zhighheavy",
        type=float,
        help="indicate the max value for the z axis \
                              for heavy chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zlh",
        "--zlowheavy",
        type=float,
        help="indicate the min value for the z axis \
                              for heavy chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-yhk",
        "--yhighkappa",
        type=int,
        help="indicate the max value for the y axis \
                              for kappa chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-ylk",
        "--ylowkappa",
        type=int,
        help="indicate the min value for the y axis \
                              for kappa chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zhk",
        "--zhighkappa",
        type=float,
        help="indicate the max value for the z axis \
                              for kappa chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zlk",
        "--zlowkappa",
        type=float,
        help="indicate the min value for the z axis \
                              for kappa chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-yhl",
        "--yhighlambda",
        type=int,
        help="indicate the max value for the y axis \
                              for lambda chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-yll",
        "--ylowlambda",
        type=int,
        help="indicate the min value for the y axis \
                              for lambda chain, default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zhl",
        "--zhighlambda",
        type=float,
        help="indicate the max value for the z axis \
                              for lambda chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-zll",
        "--zlowlambda",
        type=float,
        help="indicate the min value for the z axis \
                              for lambda chain, default = adapts to the data",
        metavar="float",
        required=False,
    )
    parser_vgenes.add_argument(
        "-c",
        "--compare",
        type=str,
        help="indicate if you want to compare  \
                              surface plots, default = False",
        metavar="True or False",
        default="False",
        required=False,
    )
    parser_vgenes.set_defaults(func=tcrcloud.vgenes.barplot)

    # create subparser for making the amino acids plot
    parser_aminoacids = subparsers.add_parser(
        "aminoacids",
        help="Create a \
        amino acids plot from AIRR CDR3 data",
    )

    parser_aminoacids.add_argument(
        "rearrangements",
        type=str,
        help="indicate the name of the AIRR \
                                   Standards rearrangements file",
        metavar="rearrangements.tsv",
    )
    parser_aminoacids.add_argument(
        "-l",
        "--length",
        type=int,
        help="indicate the value for the axis \
                                   representing the length of the CDR3, \
                                   default = adapts to the data",
        metavar="integer",
        required=False,
    )
    parser_aminoacids.add_argument(
        "-t",
        "--threeD",
        type=str,
        help="indicate if you want a tridimensional  \
                                   bar plot, default = False",
        metavar="True or False",
        default="False",
        required=False,
    )
    parser_aminoacids.add_argument(
        "-c",
        "--compare",
        type=str,
        help="indicate if you want to compare  \
                              3D barplot plots, default = False",
        metavar="True or False",
        default="False",
        required=False,
    )
    parser_aminoacids.add_argument(
        "-e",
        "--export",
        type=str,
        help="indicate if the metrics from the plot \
                              should be exported to a csv file, \
                              default = False",
        metavar="True or False",
        default="False",
        required=False,
    )
    parser_aminoacids.set_defaults(func=tcrcloud.aminoacids.aminoacids)

    # create subparser for downloading the rearrangement data
    parser_download = subparsers.add_parser(
        "download",
        help="Download TCR AIRR-seq \
                                            rearrangements data matching a \
                                            given repertoire metadata file",
    )

    parser_download.add_argument(
        "repertoire",
        type=str,
        help="indicate the name of the \
                                 AIRR repertoire file",
        metavar="repertoires.airr.json",
    )

    parser_download.set_defaults(func=tcrcloud.download.airrdownload)

    # create subparser for downloading the test repertoire
    parser_testdata = subparsers.add_parser(
        "testdata",
        help="Download example TCR \
                                            AIRR-seq repertoire data  \
                                            to test TCRcloud",
    )

    parser_testdata.set_defaults(func=tcrcloud.testdata.download)

    args = parser.parse_args()
    try:
        args.func(args)
    except FileNotFoundError:
        if dir(args)[-2] == "repertoire":
            sys.stderr.write(
                "TCRcloud error: " + args.repertoire + " doesn't seem to exist\n"
            )
        elif dir(args)[-2] == "rearrangements":
            sys.stderr.write(
                "TCRcloud error: " + args.rearrangements + " doesn't seem to exist\n"
            )
    except (yaml.scanner.ScannerError, json.decoder.JSONDecodeError):
        sys.stderr.write(
            "TCRcloud error: It seems you did not indicate a \
properly formatted AIRR repertoire file\n"
        )
    except (KeyError, TypeError):
        sys.stderr.write(
            "TCRcloud error: It seems you did not indicate a \
properly formatted AIRR rearrangements file\n"
        )


if __name__ == "__main__":
    main()
