#!/usr/bin/env python3
import argparse
import logging
import sys
from importlib.metadata import PackageNotFoundError, version

import tcrcloud.aminoacids
import tcrcloud.cloud
import tcrcloud.compare
import tcrcloud.download
import tcrcloud.radar
import tcrcloud.testdata
import tcrcloud.vgenes
from tcrcloud.errors import TCRcloudError

try:
    __version__ = version("TCRcloud")
except PackageNotFoundError:
    # Package isn't installed (e.g. running directly from a source checkout).
    __version__ = "0.0.0.dev"


# Species with a built-in V-gene colour palette (see tcrcloud.colours*).
_SPECIES = (
    "homo_sapiens",
    "mus_musculus",
    "macaca_mulatta",
    "macaca_fascicularis",
)


def str2bool(v: bool | str) -> bool:
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        if v.lower() in ("yes", "true", "t", "y", "1"):
            return True
        if v.lower() in ("no", "false", "f", "n", "0"):
            return False
    raise argparse.ArgumentTypeError("Expected a boolean value (true/false)")


def main() -> None:
    # Surface library status/progress messages (logged at INFO) on the CLI.
    # Library modules log via `logging.getLogger(__name__)` instead of
    # printing, so importers of the package stay free of stdout noise.
    logging.basicConfig(level=logging.INFO, format="%(message)s")

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
        choices=_SPECIES,
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
        choices=_SPECIES,
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
        choices=_SPECIES,
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

    # create subparser for downloading the test repertoire
    subparsers.add_parser(
        "testdata",
        help="Download example TCR \
                                            AIRR-seq repertoire data  \
                                            to test TCRcloud",
    )

    args = parser.parse_args()
    try:
        # Explicit per-command dispatch: the CLI layer translates the parsed
        # argparse namespace into typed function calls, so library modules
        # never see argparse objects and keep real, type-checkable
        # signatures.
        if args.command == "cloud":
            tcrcloud.cloud.wordcloud(
                args.rearrangements,
                colours=args.colours,
                species=args.species,
                legend=args.legend,
                size=args.size,
                output_format=args.format,
            )
        elif args.command == "radar":
            tcrcloud.radar.radar(
                args.rearrangements,
                custom_legend=args.custom_legend,
                legend=args.legend,
                export=args.export,
                output_format=args.format,
            )
        elif args.command == "vgenes":
            tcrcloud.vgenes.barplot(
                args.rearrangements,
                export=args.export,
                species=args.species,
                output_format=args.format,
            )
        elif args.command == "aminoacids":
            tcrcloud.aminoacids.aminoacids(
                args.rearrangements,
                three_d=args.threeD,
                export=args.export,
                output_format=args.format,
            )
        elif args.command == "compare":
            tcrcloud.compare.compare(
                args.rearrangements,
                rearrangements2=args.rearrangements2,
                species=args.species,
                export=args.export,
                output_format=args.format,
            )
        elif args.command == "download":
            tcrcloud.download.airrdownload(args.repertoire)
        elif args.command == "testdata":
            tcrcloud.testdata.download()
    except TCRcloudError as exc:
        # Expected, user-facing failures raised by library modules. Print a
        # clean message and exit non-zero (so the failure is visible to
        # scripts/pipelines), instead of letting helpers call sys.exit().
        sys.stderr.write(f"TCRcloud error: {exc}\n")
        sys.exit(1)
    except FileNotFoundError as exc:
        # A bare FileNotFoundError, e.g. from open() on a missing input
        # file. Library modules raise TCRcloudError for expected failures,
        # so this only handles unexpected ones: report the input filename
        # instead of the raw errno message.
        filename = getattr(args, "repertoire", None) or getattr(
            args, "rearrangements", None
        )
        message = f"{filename} doesn't seem to exist" if filename else str(exc)
        sys.stderr.write(f"TCRcloud error: {message}\n")
        sys.exit(1)

    # Deliberately no catch-all for ValueError/KeyError/TypeError/etc.
    # Expected user-input failures are validated and converted to
    # TCRcloudError at the module boundary (see format.py/download.py); a
    # broad catch here would conflate genuine programming bugs with bad
    # user input, misreporting them as "not properly formatted" and
    # hiding the traceback needed to fix the actual bug. Unexpected
    # exceptions should propagate with a traceback.


if __name__ == "__main__":
    main()
