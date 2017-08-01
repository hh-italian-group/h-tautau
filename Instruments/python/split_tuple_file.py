#!/usr/bin/env python
# Split TTree into a separate file.
# This file is part of https://github.com/hh-italian-group/AnalysisTools.

import sys
import os
import subprocess
import re
import argparse
import shutil
import ROOT

parser = argparse.ArgumentParser(description='Split tuple ROOT file.',
                                 formatter_class = lambda prog: argparse.HelpFormatter(prog,width=90))
parser.add_argument('--input', required=True, dest='input', type=str, metavar='FILE', help="original ROOT file")
parser.add_argument('--output', required=True, dest='output', type=str, metavar='FILE', help="output ROOT file")
parser.add_argument('--max-size', required=False, dest='max_size', type=float, metavar='SIZE', default=8000,
                    help="maximal file size after the split in MiB")
parser.add_argument('--min-split-size', required=False, dest='min_split_size', type=float, metavar='SIZE', default=10,
                    help="minimal size of a split of a tuple in MiB")
args = parser.parse_args()

def IsComplete(tuple, entry_id):
    tuple.GetEntry(entry_id)
    storage_mode = tuple.GetListOfLeaves().FindObject('storageMode').GetValue()
    return storage_mode == 0

def FindSplitPoint(tuple, entry_id):
    if entry_id + 1 >= tuple.GetEntries(): return entry_id
    while entry_id >= 0:
        if IsComplete(tuple, entry_id + 1):
            return entry_id
        entry_id -= 1
    raise RuntimeError("Can't find complete entry")

MiB_factor = 2 ** 20
accuracy_margin = 0.95

class TupleIter:
    def __init__(self, tuple, pos = 0):
        self.tuple = tuple
        self.pos = pos
        self.n_entries = tuple.GetEntries()
        self.size = float(tuple.GetZipBytes()) / MiB_factor
        self.size_per_entry = self.size / self.n_entries

    def HasUnwrittenEvents(self):
        return self.pos < self.n_entries

    def Write(self, file, max_size):
        n_entries = min(self.n_entries - self.pos, int(max_size / self.size_per_entry))
        n_entries = FindSplitPoint(self.tuple, self.pos + n_entries - 1) - self.pos + 1
        file.cd()
        new_tuple = self.tuple.CopyTree('', '', n_entries, self.pos)
        new_tuple.Write()
        print '{}: {} entries has been written.'.format(self.tuple.GetName(), n_entries)
        self.pos += n_entries

class TupleList:
    def __init__(self, input_file, tuple_names, max_size, min_size):
        self.tuple_iters = [ TupleIter(input_file.Get(name)) for name in tuple_names ]
        self.current_iter = 0
        self.max_size = max_size
        self.min_size = min_size

    def HasUnwrittenEvents(self):
        return self.current_iter < len(self.tuple_iters)

    def Write(self, file):
        while self.HasUnwrittenEvents():
            max_write_size = self.max_size - float(file.GetSize()) / MiB_factor
            if max_write_size < self.min_size: break
            iter = self.tuple_iters[self.current_iter]
            iter.Write(file, max_write_size)
            if not iter.HasUnwrittenEvents():
                self.current_iter += 1


tuple_names = [ 'eTau', 'muMu', 'muTau', 'tauTau' ]
name_prefix = re.sub('\.[^\.]*$', '', args.output)
name_format = '{}_part{}.root'

print 'Creating a copy of the original file...'
name = name_format.format(name_prefix, 1)
shutil.copyfile(args.input, name)

input_file = ROOT.TFile(args.input, 'READ')
tuple_list = TupleList(input_file, tuple_names, args.max_size * accuracy_margin, args.min_split_size)

print 'Removing tuples from the copy...'
output_file = ROOT.TFile(name, 'UPDATE')
for tuple_name in tuple_names:
    output_file.Delete(tuple_name + ';*')
output_file.Close()
subprocess.check_output([ 'hadd -f9 {} {}'.format(args.output, name) ], shell=True)
os.remove(name)

print 'Writing part 1...'
output_file = ROOT.TFile(args.output, 'UPDATE')
tuple_list.Write(output_file)

n = 2
while tuple_list.HasUnwrittenEvents():
    print 'Writing part {}...'.format(n)
    name = name_format.format(name_prefix, n)
    output_file = ROOT.TFile(name, 'RECREATE')
    tuple_list.Write(output_file)
    output_file.Close()
    n += 1
