#!/usr/bin/env python
# Filter tuple events.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse

parser = argparse.ArgumentParser(description='Filter tuple events.')
parser.add_argument('--input', required=True, type=str, help="input file")
parser.add_argument('--output', required=True, type=str, help="output file")
parser.add_argument('--trees', required=True, type=str, help="comma separated list of TTrees")
parser.add_argument('--filter-file', required=True, type=str, help=".h file with the filter definition")
parser.add_argument('--filter', required=True, type=str, help="filter string")
parser.add_argument('--other-trees', required=False, type=str, default=None,
                    help="comma separated list of TTrees that will be copied without applying the filter")
parser.add_argument('--n-threads', required=False, type=int, default=4, help="number of threads")
parser.add_argument('--comp-algo', required=False, type=str, default='kLZ4',
                    help="compression algorithm for the output file")
parser.add_argument('--comp-level', required=False, type=int, default=5,
                    help="compression level for the output file")
args = parser.parse_args()

import os
import ctypes
import ROOT

if args.n_threads > 1:
    ROOT.ROOT.EnableImplicitMT(args.n_threads)
ROOT.gROOT.SetBatch(True)

if not os.path.exists(args.filter_file):
    raise RuntimeError("File with filter definition does not exists.")
print("Loading filter file...")
error_code = ctypes.c_int(17)
ROOT.gInterpreter.ProcessLine('.L {}+'.format(args.filter_file), ctypes.addressof(error_code))
if error_code.value != ROOT.TInterpreter.kNoError:
    raise RuntimeError("Error while loading the filter file. Error code = {}".format(error_code.value) )

if os.path.exists(args.output):
    print('Removing existing output file "{}"...'.format(args.output))
    os.remove(args.output)

opt = ROOT.RDF.RSnapshotOptions()
opt.fCompressionAlgorithm = getattr(ROOT.ROOT, args.comp_algo)
opt.fCompressionLevel = args.comp_level
opt.fMode = 'UPDATE'

trees = args.trees.split(',')
for tree_name in trees:
    print("Processing {}...".format(tree_name))
    df = ROOT.RDataFrame(tree_name, args.input)
    df_filtered = df.Filter(args.filter)
    cnt_before = df.Count()
    cnt_after = df_filtered.Count()
    print('Number of events before the filter: {: >8}'.format(cnt_before.GetValue()))
    print('Number of events after the filter:  {: >8}'.format(cnt_after.GetValue()))
    df_filtered.Snapshot(tree_name, args.output, '.*', opt)

if args.other_trees is not None:
    other_trees = args.other_trees.split(',')
    for tree_name in other_trees:
        print("Copying {}...".format(tree_name))
        df = ROOT.RDataFrame(tree_name, args.input)
        df.Snapshot(tree_name, args.output, '.*', opt)
