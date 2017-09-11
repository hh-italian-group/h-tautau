#!/usr/bin/env python
# Create json files to split dataset into several parts.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
import re
import subprocess
import json
import sys
from LumiList import LumiList
from multiprocessing import Process, Queue

parser = argparse.ArgumentParser(description='Create json files to split dataset into several parts.',
                  formatter_class = lambda prog: argparse.HelpFormatter(prog,width=90))
parser.add_argument('--dataset', required=True, dest='dataset', type=str, help="Dataset name")
parser.add_argument('--output-prefix', required=True, dest='output_prefix', type=str,
                    help="Prefix for output splitted json files")
parser.add_argument('--output-suffix', required=False, dest='output_suffix', type=str, default='sub',
                    help="Prefix for output splitted json files")
parser.add_argument('--n-splits', required=True, dest='n_splits', type=int, help="Number of splits")
parser.add_argument('--n-parallel', required=False, dest='n_parallel', type=int, default=20,
                    help="Number of parallel request to the DAS server.")
args = parser.parse_args()

if args.n_splits < 1:
    raise RuntimeError('Number of splits should be >= 1.')

def RunDasClient(query):
    cmd = 'das_client.py --limit=0 --format=json --query="{}"'.format(query)
    print '>>', cmd
    out = str(subprocess.check_output([cmd], shell=True))
    return json.loads(out)

def LoadFiles(dataset):
    files = []
    full_json = RunDasClient('file dataset={}'.format(dataset))
    for data_entry in full_json['data']:
        files.append(data_entry['file'][0]['name'])
    return files

def LoadFileLumis(queue, file_name, file_number):
    lumis = {}
    full_json = RunDasClient('lumi file={}'.format(file_name))
    for data_entry in full_json['data']:
        for entry in data_entry['lumi']:
            key = entry['run_number']
            if key not in lumis:
                lumis[key] = []
            lumis[key].extend(entry['number'])
    queue.put((LumiList(compactList=lumis), file_number))

def SaveLumis(file_name, lumis):
    lumi_file = open(file_name, 'w')
    lumi_file.write(str(lumis))
    lumi_file.close()

files = LoadFiles(args.dataset)
n_files = len(files)
splits = [ n_files / args.n_splits ] * args.n_splits
splits[-1] = n_files - splits[0] * (args.n_splits - 1)

print "Total number of files:", n_files
print "Number of files per split:", splits

split_number = 0
file_number = 0
all_lumis = LumiList()
for split_number in range(0, args.n_splits):
    print "Processing split {}...".format(split_number + 1)
    split_lumis = LumiList()
    queue = Queue()
    processes = {}
    n = 0
    n_processed = 0
    while n_processed < splits[split_number]:
        if n >= args.n_parallel:
            file_lumis,file_number = queue.get()
            processes[file_number].join()
            split_lumis += file_lumis
            n_processed += 1
        if n < splits[split_number]:
            p = Process(target=LoadFileLumis, args=(queue, files[file_number], file_number))
            p.start()
            processes[file_number] = p
            n += 1
            file_number += 1

    all_lumis += split_lumis
    file_name = '{}_{}{}.json'.format(args.output_prefix, args.output_suffix, split_number + 1)
    SaveLumis(file_name, split_lumis)

file_name = '{}.json'.format(args.output_prefix)
SaveLumis(file_name, all_lumis)

print "All splits have been processed."
