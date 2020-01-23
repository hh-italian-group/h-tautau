#!/usr/bin/env python
# Submit CacheTupleProducer jobs.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
import os
import re
import shutil
import subprocess

parser = argparse.ArgumentParser(description='Submit CacheTupleProducer jobs.')
parser.add_argument('--input', required=True, type=str, help="path with original root-tuples")
parser.add_argument('--output', required=True, type=str, help="path where to store job outputs")
parser.add_argument('--central-output', required=True, type=str, help="path where to store the final outputs")
parser.add_argument('--pattern', required=True, type=str, help="regex for file names")
parser.add_argument('--unc-sources', required=True, type=str, help="list of uncertainty sources")
parser.add_argument('--period', required=True, type=str, help="data taking period")
parser.add_argument('--queue', required=True, type=str, help="queue to submit the job")
parser.add_argument('--long-jobs', required=False, type=str, default=None, help="file with the list of long jobs")
parser.add_argument('--channels', required=False, type=str, default='eTau,muTau,tauTau', help="list of channels")
parser.add_argument('--selections', required=False, type=str, default='HH_legacy,HH', help="list of selections")
parser.add_argument('--jet-orderings', required=False, type=str, default='CSV,DeepCSV,DeepFlavour',
                    help="list of jet orderings")
parser.add_argument('--data-file-pattern', required=False, type=str, default='.*_Run201[678][A-H].*',
                    help="pattern to determine whatever the input file is data")
parser.add_argument('--cache-tuple-producer', required=False, type=str, default='../build/h-tautau/CacheTupleProducer',
                    help="path to the CacheTupleProducer executable")
parser.add_argument('--hasBjetPair', required=False, type=int, default=1,
                    help="require 2 b jet candidates to be present in the event")
parser.add_argument('--runSVFit', required=False, type=int, default=1, help="run SVfit")
parser.add_argument('--runKinFit', required=False, type=int, default=1, help="run KinFit")
parser.add_argument('--verbose', action="store_true", help="print verbose output")
args = parser.parse_args()

class Job:
    def __init__(self, file_name, channel, is_data, unc_sources, output, central_output, entry_range=None):
        self.input_file = file_name
        self.channel = channel
        self.is_data = is_data
        self.unc_sources = unc_sources
        self.input_basename = os.path.basename(file_name)
        self.input_name = os.path.splitext(self.input_basename)[0]

        if entry_range is None:
            self.range_index = None
            self.begin_entry_index = None
            self.end_entry_index = None
            self.name = '{}_{}'.format(self.input_name, channel)
            range_suffix = ''
        else:
            self.range_index = entry_range[0]
            self.begin_entry_index = entry_range[1]
            self.end_entry_index = entry_range[2]
            range_suffix = '_cache{}'.format(self.range_index)

        self.name = '{}_{}{}'.format(self.input_name, channel, range_suffix)
        self.output = os.path.abspath(os.path.join(output, ch, '{}{}'.format(self.input_name, range_suffix)))
        output_file_name = '{}{}.root'.format(self.input_name, range_suffix)
        self.output_file = os.path.join(self.output, output_file_name)
        self.central_output_dir = os.path.join(central_output, ch)
        self.central_output_file = os.path.join(self.central_output_dir, output_file_name)

    def __repr__(self):
        return self.name

def makedirs(path, mode):
    if not os.path.exists(path):
        original_umask = os.umask(0)
        try:
            os.makedirs(path, mode)
        finally:
            os.umask(original_umask)

def getNumberOfEntries(file_name, tree_name):
    import ROOT
    file = ROOT.TFile(file_name)
    tree = file.Get(tree_name)
    n_entries = tree.GetEntries()
    file.Close()
    return n_entries

def createSplits(n_entries, n_splits):
    entries_per_split = n_entries // n_splits
    if entries_per_split <= 0:
        raise RuntimeError("Invalid splitting: n_entries={}, n_splits={}".format(n_entries, n_splits))
    n_taken = 0
    splits = []
    for i in range(n_splits - 1):
        splits.append( (i, n_taken, n_taken + entries_per_split) )
        n_taken += entries_per_split
    splits.append( (n_splits - 1, n_taken, n_entries) )
    return splits

channels = args.channels.split(',')
unc_sources = args.unc_sources.split(',')

long_jobs = {}
if args.long_jobs is not None and os.path.exists(args.long_jobs):
    with open(args.long_jobs, 'r') as f:
        for line in f.readlines():
            if len(line) == 0 or line[0] == '#': continue
            split_line = [ s for s in line.split(' ') if len(s) > 0 ]
            if len(split_line) != 2:
                raise RuntimeError("Invalid format in long_jobs file.")
            long_jobs[split_line[0]] = int(split_line[1])

jobs = []
for f in sorted(os.listdir(args.input)):
    file_name = os.path.abspath(os.path.join(args.input, f))
    if os.path.isfile(file_name) and re.match(args.pattern, f):
        is_data = re.match(args.data_file_pattern, f) is not None
        if is_data:
            if 'None' in unc_sources:
                job_unc_sources = [ 'None' ]
            else:
                continue
        else:
            job_unc_sources = unc_sources
        for ch in channels:
            ref_job = Job(file_name, ch, is_data, job_unc_sources, args.output, args.central_output)
            if ref_job.name in long_jobs:
                n_entries = getNumberOfEntries(file_name, ch)
                n_splits = long_jobs[ref_job.name]
                splits = createSplits(n_entries, n_splits)
                for entry_range in splits:
                    job = Job(file_name, ch, is_data, job_unc_sources, args.output, args.central_output, entry_range)
                    jobs.append(job)
            else:
                jobs.append(ref_job)

n_central = 0
n_transfered = 0
n_failed = 0
n_submitted = 0
n_running = 0
for job in jobs:
    if os.path.exists(job.central_output_file):
        print('{}: already in the central storage'.format(job.name))
        n_central += 1
    elif os.path.isdir(job.output):
        job_successfully_ended = False
        job_failed = False
        log_file = os.path.join(job.output, 'run_job.log')
        if os.path.isfile(job.output_file) and os.path.isfile(log_file):
            with open(log_file, 'r') as f:
                last_line = f.readlines()[-1]
            if re.match('^Job successfully ended at.*', last_line):
                job_successfully_ended = True
            elif re.match('^Job failed at.*', last_line):
                job_failed = True
        if job_successfully_ended:
            print('{}: successfully ended'.format(job.name))
            os.chmod(job.output_file, 0660)
            makedirs(job.central_output_dir, 0770)
            shutil.move(job.output_file, job.central_output_file)
            print('{}: transfered into the central storage'.format(job.name))
            n_central += 1
            n_transfered += 1
        elif job_failed:
            print('{}: failed.'.format(job.name))
            n_failed += 1
        else:
            print('{}: submitted, but not ended yet'.format(job.name))
            n_running += 1
    else:
        print('Submitting "{}"...'.format(job.name))
        submit_cmd = './AnalysisTools/Run/submit_job.sh {} {} {} {}'.format(args.queue, job.name, job.output,
                                                                            os.path.abspath(args.cache_tuple_producer))
        submit_cmd += ' --input_file {} --output_file {}'.format(job.input_file, job.output_file)
        submit_cmd += ' --channels {} --period {} --selections {}'.format(job.channel, args.period, args.selections)
        submit_cmd += ' --unc_sources {} --jet_orderings {}'.format(','.join(job.unc_sources), args.jet_orderings)
        submit_cmd += ' --isData {} --hasBjetPair {}'.format(int(job.is_data), args.hasBjetPair)
        submit_cmd += ' --runSVFit {} --runKinFit {}'.format(args.runSVFit, args.runKinFit)
        submit_cmd += ' --working_path {}'.format(os.path.abspath('.'))
        if job.range_index is not None:
            submit_cmd += ' --begin_entry_index {}'.format(job.begin_entry_index)
            submit_cmd += ' --end_entry_index {}'.format(job.end_entry_index)
        result = subprocess.call([submit_cmd], shell=True)
        if result != 0:
            raise RuntimeError("Failed to submit job {}".format(job.name))
        n_submitted += 1

print('\nJOB SUMMARY\n\t{} in the central storage, out of which {} are newly transfered' \
      '\n\t{} newly submitted\n\t{} submitted, but not ended yet\n\t{} failed' \
      .format(n_central, n_transfered, n_submitted, n_running, n_failed))
