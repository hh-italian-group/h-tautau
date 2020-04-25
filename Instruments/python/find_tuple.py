#!/usr/bin/env python
# Find tuple file.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse

parser = argparse.ArgumentParser(description='Filter tuple events.')
parser.add_argument('--input', required=True, type=str, help="input path")
parser.add_argument('--tree', required=True, type=str, help="tree name")
parser.add_argument('--event-id', required=True, type=str, help="event id")
parser.add_argument('--n-threads', required=False, type=int, default=4, help="number of threads")
args = parser.parse_args()

import os
import ctypes
import ROOT
from tqdm import tqdm

if args.n_threads > 1:
    ROOT.ROOT.EnableImplicitMT(args.n_threads)
ROOT.gROOT.SetBatch(True)

event_id = [ int(s) for s in args.event_id.split(':') ]
if len(event_id) != 3:
    raise RuntimeError("Invalid event id = {}".format(args.event_id))

filter_line = 'run == {} && lumi == {} && evt == {}ULL'.format(event_id[0], event_id[1], event_id[2])
all_files = []
for root, subdirs, files in os.walk(args.input):
    for file in files:
        if file.endswith('.root'):
            file_path = os.path.join(root, file)
            all_files.append(file_path)

for file in tqdm(all_files):
    df = ROOT.RDataFrame(args.tree, file)
    df_filtered = df.Filter(filter_line)
    cnt = df_filtered.Count()
    if cnt.GetValue() > 0:
        print('\n{}'.format(file))
