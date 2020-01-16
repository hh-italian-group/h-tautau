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
    def __init__(self, file_name, channel, is_data, unc_sources):
        self.input_file = file_name
        self.channel = channel
        self.is_data = is_data
        self.unc_sources = unc_sources
        self.input_basename = os.path.basename(file_name)
        self.input_name = os.path.splitext(self.input_basename)[0]
        self.name = '{}_{}'.format(self.input_name, channel)

    def __repr__(self):
        return self.name

channels = args.channels.split(',')
unc_sources = args.unc_sources.split(',')

def makedirs(path, mode):
    if not os.path.exists(path):
        try:
            original_umask = os.umask(0)
            os.makedirs(path, mode)
        finally:
            os.umask(original_umask)

jobs = []
for f in os.listdir(args.input):
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
            job = Job(file_name, ch, is_data, job_unc_sources)
            job.output = os.path.abspath(os.path.join(args.output, ch, job.input_name))
            job.output_file = os.path.join(job.output, job.input_basename)
            job.central_output_dir = os.path.join(args.central_output, ch)
            job.central_output_file = os.path.join(job.central_output_dir, job.input_basename)
            jobs.append(job)

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
        result = subprocess.call([submit_cmd], shell=True)
        if result != 0:
            raise RuntimeError("Failed to submit job {}".format(job.name))
        n_submitted += 1

print('\nJOB SUMMARY\n\t{} in the central storage, out of which {} are newly transfered' \
      '\n\t{} newly submitted\n\t{} submitted, but not ended yet\n\t{} failed' \
      .format(n_central, n_transfered, n_submitted, n_running, n_failed))
