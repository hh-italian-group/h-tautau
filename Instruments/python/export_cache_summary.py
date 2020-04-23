#!/usr/bin/env python
# Export cache summary into csv file.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
import os
import ROOT

parser = argparse.ArgumentParser(description='Export cache summary into csv file')
parser.add_argument('input', type=str, help="Input path")
args = parser.parse_args()

branches = [ 'exeTime', 'n_orig_events', 'n_stored_events', 'n_SVfit', 'n_KinFit', 'n_HHbtag' ]
header = [ 'channel', 'sample_name', 'n_entries', 'entry_id' ]
header.extend(branches)
print(','.join(header))
for root, subdirs, files in os.walk(args.input):
    for file in files:
        channel = os.path.basename(root)
        sample = os.path.splitext(file)[0]
        root_file = ROOT.TFile(os.path.join(root, file), 'READ')
        tree = root_file.Get('summary')
        if tree is None or type(tree) != ROOT.TTree: continue
        entry_id = 0
        n_entries = int(tree.GetEntries())
        for event in tree:
            entry = [ channel, sample, str(n_entries), str(entry_id) ]
            for branch in branches:
                entry.append(str(int(getattr(event, branch))))
            print(','.join(entry))
            entry_id += 1
