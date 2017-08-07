#!/usr/bin/env python
# Create json files to split dataset into several parts.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
from sets import Set
from FWCore.PythonUtilities.LumiList import LumiList
from dbs.apis.dbsClient import DbsApi

parser = argparse.ArgumentParser(description='Create json files to split dataset into several parts.',
                  formatter_class = lambda prog: argparse.HelpFormatter(prog,width=90))
parser.add_argument('--dataset', required=True, dest='dataset', type=str, help="Dataset name")
parser.add_argument('--output-prefix', required=True, dest='output_prefix', type=str,
                    help="Prefix for output splitted json files")
parser.add_argument('--output-suffix', required=False, dest='output_suffix', type=str, default='sub',
                    help="Prefix for output splitted json files")
parser.add_argument('--n-splits', required=True, dest='n_splits', type=int, help="Number of splits")
args = parser.parse_args()

if args.n_splits < 1:
    raise RuntimeError('Number of splits should be >= 1.')

def FindMaxLumi(dbs, dataset):
    blocks = dbs.listBlocks(dataset=dataset)
    max_lumi = 0
    for block_entry in blocks:
        block_lumis = dbs.listFileLumis(block_name=block_entry['block_name'])
        for file_entry in block_lumis:
            file_lumis = file_entry['lumi_section_num']
            max_file_lumi = max(file_lumis)
            max_lumi = max(max_lumi, max_file_lumi)
    return max_lumi

def GetRunList(dbs, dataset):
    runs = dbs.listRuns(dataset=dataset)
    run_list = []
    for run in runs:
        run_list.extend(run['run_num'])
    run_set = Set(run_list)
    return list(run_set)

def SaveLumis(file_name, lumis):
    lumi_file = open(file_name, 'w')
    lumi_file.write(str(lumis))
    lumi_file.close()

dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')

print("Loading runs...")
runs = GetRunList(dbs, args.dataset)
if len(runs) != 1:
    raise RuntimeError('Only datasets with one run are currently supported.')

print("Loading lumis...")
max_lumi = FindMaxLumi(dbs, args.dataset)
splits = [ int(float(n + 1) / args.n_splits * max_lumi) for n in range(0, args.n_splits) ]

print("Max lumi: {}".format(max_lumi))
print("Lumi splits: {}".format(splits))

last_lumi = 0
for split_number in range(0, len(splits)):
    split = splits[split_number]
    lumis = {}
    lumis[runs[0]] = []
    lumis[runs[0]].append([last_lumi + 1, split])
    file_name = '{}_{}{}.json'.format(args.output_prefix, args.output_suffix, split_number + 1)
    SaveLumis(file_name, LumiList(compactList=lumis))
    last_lumi = split

print("Dataset lumis are split into {} parts.".format(args.n_splits))
