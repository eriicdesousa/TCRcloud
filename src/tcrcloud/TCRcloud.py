#!/usr/bin/python3
import argparse

import matplotlib.pyplot as plt

import tcrcloud.cloud
import tcrcloud.radar
import tcrcloud.download

plt.rcParams["font.family"] = "serif"


def main():
    # create the top-level parser
    parser = argparse.ArgumentParser(
        description='Create a wordcloud of CDR3 \
    sequences from TCR AIRR-seq data or a radar plot with diversity metrics.',
        prog='TCRcloud')

    subparsers = parser.add_subparsers(
        title='Command Options',
        help='The program has 3 modes: cloud, radar or download',
        dest='cloud, radar, download or -h/--help',
        required=True)

    # create subparser for making the wordcloud
    parser_cloud = subparsers.add_parser('cloud', help='Create a wordcloud \
        from TCR CDR3 data')

    required_group = parser_cloud.add_argument_group('required arguments')

    # required_group.add_argument('-r','--repertoire', 
    #     type= tcrcloud.base.jsonfile, 
    #     help= 'indicate the name of the AIRR Standards repertoire file',
    #     metavar= 'repertoires.airr.json', 
    #     required= True)

    required_group.add_argument('-s', '--rearrangements', type=str,
                                help='indicate the name of the AIRR Standards \
                                rearrangements file',
                                metavar='rearrangements.tsv', required=True)

    parser_cloud.set_defaults(func=tcrcloud.cloud.wordcloud)

    # create subparser for making the radar
    parser_radar = subparsers.add_parser('radar', help='Only calculate \
    the metrics and print to file')
    required_group = parser_radar.add_argument_group('required arguments')
    required_group.add_argument('-s', '--rearrangements', type=str,
                                help='indicate the name of the AIRR \
                                Standards rearrangements file',
                                metavar='rearrangements.tsv', required=True)
    parser_radar.set_defaults(func=tcrcloud.radar.radar)

    # create subparser for downloading the data
    parser_download = subparsers.add_parser('download',
                                            help='Download TCR AIRR-seq \
                                            rearrangements data matching a \
                                            given repertoire metadata file')

    required_group = parser_download.add_argument_group('required arguments')

    required_group.add_argument('-r', '--repertoire', type=str,
                                help='indicate the name of the \
                                AIRR Standards repertoire file',
                                metavar='repertoires.airr.json',
                                required=True)

    parser_download.set_defaults(func=tcrcloud.download.airrdownload)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
