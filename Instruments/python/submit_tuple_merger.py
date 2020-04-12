#!/usr/bin/env python
# Submit TupleMerger jobs.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
import os
import re
import shutil
import subprocess

parser = argparse.ArgumentParser(description='Submit TupleMerger jobs.')
parser.add_argument('--crab-outputs', required=True, type=str, help="path to CRAB outputs")
parser.add_argument('--finished-tasks', required=True, type=str,
                    help="file with the list tasks that are finished on CRAB")
parser.add_argument('--output', required=True, type=str, help="path where to store job outputs")
parser.add_argument('--central-output', required=True, type=str, help="path where to store the final outputs")
parser.add_argument('--queue', required=True, type=str, help="queue to submit the job")
parser.add_argument('--tuple-merger', required=False, type=str, default='../build/h-tautau/TupleMerger',
                    help="path to the TupleMerger executable")
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
    NewlyTransfered = 16

class CrabTask:
    def __init__(self, full_id, crab_outputs):
        self.full_id = full_id
        id_parts = re.search('^([0-9]+_[0-9]+):(.*_)(crab_)([^ ]*)$', full_id)
        if id_parts is None:
            raise RuntimeError('Invalid crab task {}.'.format(full_id))
        self.number = id_parts.group(1)
        self.name = id_parts.group(4)
        self.root_directory = id_parts.group(3) + self.name
        self.path = None
        paths = []
        if os.path.exists(crab_outputs):
            for d in os.listdir(crab_outputs):
                path = os.path.abspath(os.path.join(crab_outputs, d, self.root_directory, self.number))
                if os.path.exists(path):
                    paths.append(path)
        if len(paths) > 1:
            raise RuntimeError("Found more than one path for crab task {}.".format(full_id))
        if len(paths) > 0:
            self.path = paths[0]
        name_parts = re.search('^(.*)_recovery[0-9]+', self.name)
        self.is_recovery = name_parts is not None
        self.dataset_name = name_parts.group(1) if self.is_recovery else self.name

    def __repr__(self):
        return self.full_id

class MergeJob:
    def __init__(self, output, central_output, crab_tasks):
        self.status = JobStatus.Unknown
        if len(crab_tasks) == 0:
            raise RuntimeError("Unable to create merge job without crab tasks.")
        self.crab_tasks = crab_tasks
        self.name = crab_tasks[0].dataset_name
        for task in crab_tasks:
            if task.dataset_name != self.name:
                raise RuntimeError('Inconsistent input crab tasks for the merge job "{}".'.format(self.name))
        self.output = os.path.abspath(os.path.join(output, self.name))
        output_file_name = '{}.root'.format(self.name)
        self.output_file = os.path.join(self.output, output_file_name)
        self.central_output_dir = central_output
        self.central_output_file = os.path.join(central_output, output_file_name)

    def UpdateStatus(self, status_bits):
        self.status = self.status | status_bits

    def HasStatus(self, status_bits):
        return self.status & status_bits == status_bits

def makedirs(path, mode):
    if not os.path.exists(path):
        original_umask = os.umask(0)
        try:
            os.makedirs(path, mode)
        finally:
            os.umask(original_umask)


crab_tasks = {}
with open(args.finished_tasks, 'r') as f:
    for line_number, line in enumerate(f.readlines()):
        line = line.strip()
        if len(line) == 0 or line[0] == '#': continue
        task = CrabTask(line, args.crab_outputs)
        if task.dataset_name not in crab_tasks:
            crab_tasks[task.dataset_name] = []
        crab_tasks[task.dataset_name].append(task)

def ProcessExistingJob(job):
    job_successfully_ended = False
    job_failed = False
    log_file = os.path.join(job.output, 'run_job.log')
    if os.path.isfile(log_file):
        with open(log_file, 'r') as f:
            log_file_lines = f.readlines()
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

merge_jobs = []

for dataset_name in sorted(crab_tasks):
    job = MergeJob(args.output, args.central_output, crab_tasks[dataset_name])
    merge_jobs.append(job)
    if os.path.exists(job.central_output_file):
        print('{}: already in the central storage'.format(job.name))
        job.UpdateStatus(JobStatus.Finished)
    elif os.path.isdir(job.output):
        ProcessExistingJob(job)
    else:
        tasks_str = ', '.join([task.full_id for task in job.crab_tasks])
        print('Submitting {} ({})...'.format(job.name, tasks_str))
        submit_cmd = './AnalysisTools/Run/submit_job.sh {} {} {} {}'.format(args.queue, job.name, job.output,
                                                                            os.path.abspath(args.tuple_merger))
        submit_cmd += ' --output {} '.format(job.output_file)
        for task in job.crab_tasks:
            if task.path is None:
                raise RuntimeError('Unable to find outputs for crab task {}.'.format(task.full_id))
            submit_cmd += ' --input-dir {} '.format(task.path)
        result = subprocess.call([submit_cmd], shell=True)
        if result != 0:
            raise RuntimeError("Failed to submit job {}".format(job.name))
        job.UpdateStatus(JobStatus.Submitted)

class JobCounter:
    def __init__(self, title, jobs=None):
        self.title = title
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
        s += '\n\t{} newly submitted'.format(self.counts[JobStatus.Submitted])
        s += '\n\t{} submitted, but not ended yet'.format(self.counts[JobStatus.Running])
        s += '\n\t{} failed'.format(self.counts[JobStatus.Failed])
        return s

print(JobCounter('MERGE', merge_jobs))
