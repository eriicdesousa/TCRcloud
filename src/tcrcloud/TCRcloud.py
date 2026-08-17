#!/usr/bin/env python3
import argparse
import yaml
import json
import sys
from importlib.metadata import PackageNotFoundError, version

import matplotlib.pyplot as plt

import tcrcloud.cloud
import tcrcloud.radar
import tcrcloud.download
import tcrcloud.testdata
import tcrcloud.vgenes
import tcrcloud.aminoacids
import tcrcloud.compare
from tcrcloud.errors import TCRcloudError

try:
    __version__ = version("TCRcloud")
except PackageNotFoundError:
    # Package isn't installed (e.g. running directly from a source checkout).
    __version__ = "0.0.0.dev"


def str2bool(v):
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        if v.lower() in ("yes", "true", "t", "y", "1"):
            return True
        if v.lower() in ("no", "false", "f", "n", "0"):
            return False
    raise argparse.ArgumentTypeError("Expected a boolean value (true/false)")


plt.rcParams["font.family"] = "serif"


def main():
    # create the top-level parser
    parser = argparse.ArgumentParser(
        description="Create visualizations from AIRR-seq data", prog="TCRcloud"
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    subparsers = parser.add_subparsers(
        title="command options",
        help="The program has 7 options: cloud, radar, vgenes, aminoacids, compare, download or testdata",
        dest="command",
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
        "-p",
        "--species",
        type=str,
        choices=[
            "homo_sapiens",
            "mus_musculus",
            "macaca_mulatta",
            "macaca_fascicularis",
        ],
        help="Species to use for built-in V-gene colours (homo_sapiens, mus_musculus, macaca_mulatta, macaca_fascicularis)",
        metavar="species",
        default="homo_sapiens",
        required=False,
    )
    parser_cloud.add_argument(
        "-l",
        "--legend",
        type=str2bool,
        help="Include a legend in the word cloud (True/False)",
        metavar="True or False",
        default=True,
        required=False,
    )
    parser_cloud.add_argument(
        "-s", "--size", type=int, help=argparse.SUPPRESS, default=1000, required=False
    )
    parser_cloud.add_argument(
        "-f",
        "--format",
        type=str,
        choices=["svg", "png"],
        help="Output image format for the word cloud (svg or png)",
        metavar="svg or png",
        default="png",
        required=False,
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
        type=str2bool,
        help="Include legend in the radar plot (True/False)",
        metavar="True or False",
        default=True,
        required=False,
    )
    parser_radar.add_argument(
        "-e",
        "--export",
        type=str2bool,
        help="Export computed radar metrics to a text file (True/False)",
        metavar="True or False",
        default=False,
        required=False,
    )
    parser_radar.add_argument(
        "-f",
        "--format",
        type=str,
        choices=["svg", "png"],
        help="Output image format for the radar plot (svg or png)",
        metavar="svg or png",
        default="png",
        required=False,
    )
    # Note: min/max scaling for radar metrics is fixed / computed automatically.
    # Removed CLI options for explicit min/max values to simplify usage.

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
        type=str2bool,
        help="Export the V gene usage table to a csv file (True/False)",
        metavar="True or False",
        default=False,
        required=False,
    )
    parser_vgenes.add_argument(
        "-p",
        "--species",
        type=str,
        choices=[
            "homo_sapiens",
            "mus_musculus",
            "macaca_mulatta",
            "macaca_fascicularis",
        ],
        help="Species to use for built-in V-gene (homo_sapiens, mus_musculus, macaca_mulatta, macaca_fascicularis)",
        metavar="species",
        default="homo_sapiens",
        required=False,
    )
    parser_vgenes.add_argument(
        "-f",
        "--format",
        type=str,
        choices=["svg", "png"],
        help="Output image format for the V gene plot (svg or png)",
        metavar="svg or png",
        default="png",
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
        "-t",
        "--threeD",
        type=str2bool,
        help="indicate if you want a tridimensional  \
                                   bar plot, default = False",
        metavar="True or False",
        default=False,
        required=False,
    )
    parser_aminoacids.add_argument(
        "-e",
        "--export",
        type=str2bool,
        help="indicate if the metrics from the plot \
                              should be exported to a csv file, \
                              default = False",
        metavar="True or False",
        default=False,
        required=False,
    )
    parser_aminoacids.add_argument(
        "-f",
        "--format",
        type=str,
        choices=["svg", "png"],
        help="Output image format for the amino acids plot (svg or png)",
        metavar="svg or png",
        default="png",
        required=False,
    )
    parser_aminoacids.set_defaults(func=tcrcloud.aminoacids.aminoacids)

    # create subparser for comparing 2 repertoires
    parser_compare = subparsers.add_parser(
        "compare",
        help="Compare repertoires that share the same chain \
        (V genes and amino acids usage)",
    )

    parser_compare.add_argument(
        "rearrangements",
        type=str,
        help="indicate the name of the AIRR \
                              rearrangements file",
        metavar="rearrangements.tsv",
    )
    parser_compare.add_argument(
        "rearrangements2",
        type=str,
        nargs="?",
        default=None,
        help="optionally indicate the name of a second AIRR rearrangements \
                              file, to compare repertoires across two files \
                              (still restricted to matching chains); if \
                              omitted, every pair of repertoires that share \
                              a chain within the single file is compared",
        metavar="rearrangements2.tsv",
    )
    parser_compare.add_argument(
        "-p",
        "--species",
        type=str,
        choices=[
            "homo_sapiens",
            "mus_musculus",
            "macaca_mulatta",
            "macaca_fascicularis",
        ],
        help="Species to use for built-in V-gene (homo_sapiens, mus_musculus, macaca_mulatta, macaca_fascicularis)",
        metavar="species",
        default="homo_sapiens",
        required=False,
    )
    parser_compare.add_argument(
        "-e",
        "--export",
        type=str2bool,
        help="Export the comparison tables to csv files (True/False)",
        metavar="True or False",
        default=False,
        required=False,
    )
    parser_compare.add_argument(
        "-f",
        "--format",
        type=str,
        choices=["svg", "png"],
        help="Output image format for the comparison plots (svg or png)",
        metavar="svg or png",
        default="png",
        required=False,
    )
    parser_compare.set_defaults(func=tcrcloud.compare.compare)

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
    except TCRcloudError as exc:
        # Expected, user-facing failures raised by library modules. Print a
        # clean message and exit non-zero (so the failure is visible to
        # scripts/pipelines), instead of letting helpers call sys.exit().
        sys.stderr.write(f"TCRcloud error: {exc}\n")
        sys.exit(1)
    except FileNotFoundError as exc:
        message = str(exc)
        if not message.startswith("TCRcloud error:"):
            filename = getattr(args, "repertoire", None) or getattr(
                args, "rearrangements", None
            )
            message = (
                f"TCRcloud error: {filename} doesn't seem to exist"
                if filename
                else message
            )
        sys.stderr.write(message + "\n")
        sys.exit(1)
    except ValueError as exc:
        sys.stderr.write(str(exc) + "\n")
        sys.exit(1)
    except (yaml.scanner.ScannerError, json.decoder.JSONDecodeError):
        sys.stderr.write("TCRcloud error: It seems you did not indicate a \
properly formatted AIRR repertoire file\n")
        sys.exit(1)
    except (KeyError, TypeError):
        sys.stderr.write("TCRcloud error: It seems you did not indicate a \
properly formatted AIRR rearrangements file\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
