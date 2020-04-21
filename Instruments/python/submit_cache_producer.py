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
parser.add_argument('--btaggers', required=False, type=str, default='DeepCSV,DeepFlavour,HHbtag',
                    help="list of b taggers")
parser.add_argument('--data-file-pattern', required=False, type=str, default='.*Run201[678][A-H].*',
                    help="pattern to determine whatever the input file is data")
parser.add_argument('--cache-tuple-producer', required=False, type=str, default='../build/h-tautau/CacheTupleProducer',
                    help="path to the CacheTupleProducer executable")
parser.add_argument('--cache-merger-cpp', required=False, type=str, default='../build/h-tautau/CacheMerger',
                    help="path to the CacheMerger executable")
parser.add_argument('--cache-merger-py', required=False, type=str,
                    default='./h-tautau/Instruments/python/cache_merger.py',
                    help="path to the cache_merger")
parser.add_argument('--hasBjetPair', required=False, type=int, default=1,
                    help="require 2 b jet candidates to be present in the event")
parser.add_argument('--runSVFit', required=False, type=int, default=1, help="run SVfit")
parser.add_argument('--runKinFit', required=False, type=int, default=1, help="run KinFit")
parser.add_argument('--verbose', action="store_true", help="print verbose output")
args = parser.parse_args()

def EnumValues(cl):
    return [k for k in dir(cl) if k[0] != '_' ]

def EnumToString(cl, value):
    keys = EnumValues(cl)
    for key in keys:
        if getattr(cl, key) == value:
            return key
    raise RuntimeError('Value "{}" is not part of the enum "{}".'.format(value, cl.__name__))

class JobStatus:
    Unknown = 0
    Submitted = 1
    Running = 2
    Finished = 4
    Failed = 8
    WaitingForSubJobs = 16
    NewlyTransfered = 32

class JobBase(object):
    def __init__(self, file_name, channel):
        self.status = JobStatus.Unknown
        self.channel = channel
        self.input_basename = os.path.basename(file_name)
        self.input_name = os.path.splitext(self.input_basename)[0]

    def __repr__(self):
        return self.name

    def UpdateStatus(self, status_bits):
        self.status = self.status | status_bits
    def HasStatus(self, status_bits):
        return self.status & status_bits == status_bits

class CalcJob(JobBase):
    def __init__(self, file_name, channel, is_data, unc_sources, output, central_output, entry_range=None):
        super(CalcJob, self).__init__(file_name, channel)
        self.input_file = file_name
        self.is_data = is_data
        self.unc_sources = unc_sources
        self.merge_job = None

        if entry_range is None:
            self.range_index = None
            self.begin_entry_index = None
            self.end_entry_index = None
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

class MergeJob(JobBase):
    def __init__(self, file_name, channel, output, central_output, calc_jobs):
        super(MergeJob, self).__init__(file_name, channel)
        self.calc_jobs = calc_jobs
        for job in calc_jobs:
            job.merge_job = self

        self.name = '{}_{}_merge'.format(self.input_name, channel)
        self.output = os.path.abspath(os.path.join(output, ch, '{}_merge'.format(self.input_name)))
        output_file_name = '{}.root'.format(self.input_name)
        self.output_file = os.path.join(self.output, output_file_name)
        self.central_output_dir = os.path.join(central_output, ch)
        self.central_output_file = os.path.join(self.central_output_dir, output_file_name)

    def AllCalcJobsFinished(self):
        for job in self.calc_jobs:
            if not job.HasStatus(JobStatus.Finished):
                return False
        return True

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
        for line_number, line in enumerate(f.readlines()):
            line = line.strip()
            if len(line) == 0 or line[0] == '#': continue
            split_line = [ s.strip() for s in line.split(' ') if len(s.strip()) > 0 ]
            if len(split_line) != 2:
                raise RuntimeError("Invalid format in long_jobs file in line {}.".format(line_number + 1))
            long_jobs[split_line[0]] = int(split_line[1])

calc_jobs = []
merge_jobs = []
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
            ref_job = CalcJob(file_name, ch, is_data, job_unc_sources, args.output, args.central_output)
            if ref_job.name in long_jobs:
                n_entries = getNumberOfEntries(file_name, ch)
                n_splits = long_jobs[ref_job.name]
                splits = createSplits(n_entries, n_splits)
                file_calc_jobs = []
                for entry_range in splits:
                    job = CalcJob(file_name, ch, is_data, job_unc_sources, args.output, args.central_output,
                                  entry_range)
                    file_calc_jobs.append(job)
                merge_jobs.append(MergeJob(file_name, ch, args.output, args.central_output, file_calc_jobs))
                calc_jobs.extend(file_calc_jobs)
            else:
                calc_jobs.append(ref_job)

def ProcessExistingJob(job):
    job_successfully_ended = False
    job_failed = False
    log_file = os.path.join(job.output, 'run_job.log')
    if os.path.isfile(log_file):
        with open(log_file, 'r') as f:
            log_file_lines = f.readlines()
            #last_line = log_file_lines[-1] if len(log_file_lines) > 0 else ""
        has_errors = False
        job_successfully_ended = False

        for line in log_file_lines:
            if re.match('.*(error|warning).*', line):
                has_errors = True
            if re.match('^Job successfully ended at.*', line) is not None:
                job_successfully_ended = True
            if re.match('^Job failed at.*', line) is not None:
                job_failed = True
        if job_failed:
            job_successfully_ended = False
        if (job_successfully_ended or job_failed) and (has_errors or not os.path.isfile(job.output_file)):
            job_successfully_ended = False
            job_failed = True
    else:
        job_failed = True
    if job_successfully_ended:
        print('{}: successfully ended'.format(job.name))
        os.chmod(job.output_file, 0660)
        makedirs(job.central_output_dir, 0770)
        shutil.move(job.output_file, job.central_output_file)
        print('{}: transfered into the central storage'.format(job.name))
        job.UpdateStatus(JobStatus.Finished | JobStatus.NewlyTransfered)
    elif job_failed:
        print('{}: failed.'.format(job.name))
        job.UpdateStatus(JobStatus.Failed)
    else:
        print('{}: submitted, but not ended yet'.format(job.name))
        job.UpdateStatus(JobStatus.Running)

for job in merge_jobs:
    if os.path.exists(job.central_output_file):
        print('{}: already in the central storage'.format(job.name))
        job.UpdateStatus(JobStatus.Finished)
        for calc_job in job.calc_jobs:
            calc_job.UpdateStatus(JobStatus.Finished)
            if os.path.exists(calc_job.central_output_file):
                os.remove(calc_job.central_output_file)


for job in calc_jobs:
    if job.HasStatus(JobStatus.Finished):
        print('{}: already merged'.format(job.name))
    elif os.path.exists(job.central_output_file):
        print('{}: already in the central storage'.format(job.name))
        job.UpdateStatus(JobStatus.Finished)
    elif os.path.isdir(job.output):
        ProcessExistingJob(job)
    else:
        print('Submitting "{}"...'.format(job.name))
        submit_cmd = './AnalysisTools/Run/submit_job.sh {} {} {} {}'.format(args.queue, job.name, job.output,
                                                                            os.path.abspath(args.cache_tuple_producer))
        submit_cmd += ' --input_file {} --output_file {}'.format(job.input_file, job.output_file)
        submit_cmd += ' --channels {} --period {} --selections {}'.format(job.channel, args.period, args.selections)
        submit_cmd += ' --unc_sources {} --btaggers {}'.format(','.join(job.unc_sources), args.btaggers)
        submit_cmd += ' --isData {} --hasBjetPair {}'.format(int(job.is_data), args.hasBjetPair)
        submit_cmd += ' --runSVFit {} --runKinFit {}'.format(args.runSVFit, args.runKinFit)
        submit_cmd += ' --working_path {}'.format(os.path.abspath('.'))
        if job.range_index is not None:
            submit_cmd += ' --begin_entry_index {}'.format(job.begin_entry_index)
            submit_cmd += ' --end_entry_index {}'.format(job.end_entry_index)
        result = subprocess.call([submit_cmd], shell=True)
        if result != 0:
            raise RuntimeError("Failed to submit job {}".format(job.name))
        job.UpdateStatus(JobStatus.Submitted)

for job in merge_jobs:
    if job.HasStatus(JobStatus.Finished): continue
    if job.AllCalcJobsFinished():
        if os.path.exists(job.central_output_file):
            print('{}: already in the central storage'.format(job.name))
            job.UpdateStatus(JobStatus.Finished)
        elif os.path.isdir(job.output):
            ProcessExistingJob(job)
            if job.HasStatus(JobStatus.Finished):
                for calc_job in job.calc_jobs:
                    if os.path.exists(calc_job.central_output_file):
                        os.remove(calc_job.central_output_file)
        else:
            print('Submitting "{}"...'.format(job.name))
            submit_cmd = './AnalysisTools/Run/submit_job.sh {} {} {}'.format(args.queue, job.name, job.output)
            submit_cmd += ' python -u {} '.format(os.path.abspath(args.cache_merger_py))
            submit_cmd += ' --channel {} --output {} '.format(job.channel, job.output_file)
            submit_cmd += ' --cache-merger {} '.format(os.path.abspath(args.cache_merger_cpp))
            for calc_job in job.calc_jobs:
                submit_cmd += ' {} '.format(os.path.abspath(calc_job.central_output_file))
            result = subprocess.call([submit_cmd], shell=True)
            if result != 0:
                raise RuntimeError("Failed to submit job {}".format(job.name))
            job.UpdateStatus(JobStatus.Submitted)
    else:
        print('{}: waiting for producer jobs to finish'.format(job.name))
        job.UpdateStatus(JobStatus.WaitingForSubJobs)

class JobCounter:
    def __init__(self, title, jobs=None, report_n_waiting=False):
        self.title = title
        self.report_n_waiting = report_n_waiting
        self.counts = { }
        for key in EnumValues(JobStatus):
            self.counts[getattr(JobStatus, key)] = 0
        if jobs is not None:
            self.AddJobs(jobs)

    def Add(self, job_status):
        for key in EnumValues(JobStatus):
            val = getattr(JobStatus, key)
            if job_status & val != 0:
                self.counts[val] += 1

    def AddJobs(self, jobs):
        for job in jobs:
            self.Add(job.status)

    def __str__(self):
        s = '\n{} JOB SUMMARY'.format(self.title)
        s += '\n\t{} are finished'.format(self.counts[JobStatus.Finished])
        s += ', out of which {} are newly transfered'.format(self.counts[JobStatus.NewlyTransfered])
        if self.report_n_waiting:
            s += '\n\t{} waiting for calc jobs to finish'.format(self.counts[JobStatus.WaitingForSubJobs])
        s += '\n\t{} newly submitted'.format(self.counts[JobStatus.Submitted])
        s += '\n\t{} submitted, but not ended yet'.format(self.counts[JobStatus.Running])
        s += '\n\t{} failed'.format(self.counts[JobStatus.Failed])
        return s

print(str(JobCounter('PLAIN', [job for job in calc_jobs if job.merge_job is None])))
print(JobCounter('SPLIT', [job for job in calc_jobs if job.merge_job is not None]))
print(JobCounter('MERGE', [job for job in merge_jobs], report_n_waiting=True))
