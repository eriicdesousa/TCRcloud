#!/usr/bin/python3
import argparse

import base
import cloud
import download

# create the top-level parser
parser= argparse.ArgumentParser(
    description= 'Create a wordcloud of CDR3 \
sequences from TCR AIRR-seq data and calculate some basic metrics.',
    prog= 'TCRcloud')

subparsers= parser.add_subparsers(
    title= 'Command Options', 
    help= 'The program has 3 modes: cloud, radar or download', 
    dest= 'cloud, radar, download or -h/--help', 
    required= True)

# create subparser for making the wordcloud
parser_cloud= subparsers.add_parser('cloud', 
    help='Create a wordcloud from TCR CDR3 data')

required_group= parser_cloud.add_argument_group('required arguments')

# required_group.add_argument('-r','--repertoire', 
#     type= base.jsonfile, 
#     help= 'indicate the name of the AIRR Standards repertoire file',
#     metavar= 'repertoires.airr.json', 
#     required= True)

required_group.add_argument('-s','--rearrangements', 
    type= str, 
    help= 'indicate the name of the AIRR Standards rearrangements file', 
    metavar= 'rearrangements.tsv', 
    required= True)

parser_cloud.set_defaults(func=cloud.wordcloud)

# # create subparser for printing the metrics
# parser_metrics = subparsers.add_parser('metrics', help='Only calculate the metrics and print to file')
# required_group = parser_metrics.add_argument_group('required arguments')
# required_group.add_argument('-j','--repertoire', type=airrdownload.jsonfile, help='indicate the name of the AIRR Standards repertoire file', metavar='repertoires.airr.json', required=True)
# required_group.add_argument('-r','--rearrangements', type=str, help='indicate the name of the AIRR Standards rearrangements file', metavar='rearrangements.tsv', required=True)
# parser_metrics.set_defaults(func=metrics.metrics)

# create subparser for downloading the data
parser_testdata= subparsers.add_parser('download', 
    help= 'Download TCR AIRR-seq rearrangements data matching a given \
repertoire metadata file')

required_group= parser_testdata.add_argument_group('required arguments')

parser_testdata.add_argument('-r','--repertoire', 
    type= base.jsonfile, 
    help= 'indicate the name of the AIRR Standards repertoire file',  
    metavar= 'repertoires.airr.json')

parser_testdata.set_defaults(func=download.airrdownload)

args = parser.parse_args()
args.func(args)                                 