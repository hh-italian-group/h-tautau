#!/usr/bin/env python
# Create an extended list of CRAB jobs.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
import re
import os
import shutil
import tarfile
import json
import sys
from FWCore.PythonUtilities.LumiList import LumiList

parser = argparse.ArgumentParser(description='Create an extended list of CRAB jobs.',
                  formatter_class = lambda prog: argparse.HelpFormatter(prog,width=90))
parser.add_argument('--job-list', required=True, dest='job_list', type=str, help="simple job list.")
parser.add_argument('--work-area', required=True, dest='work_area', type=str, help="mulit-CRAB work area.")
parser.add_argument('--output', required=False, dest='output', type=str, default='',
                    help="file where to store the extended job list.")
parser.add_argument('--prev-output', required=False, dest='prev_output', type=str, default='',
                    help="file with the extended job list from the previous check.")
parser.add_argument('--crab-results-out', required=True, dest='crab_results_out', type=str,
                    help="directory where to store the processed crab results.")
args = parser.parse_args()

if len(args.output) == 0:
    args.output = args.job_list

def ReadJobs(file_name):
    jobs = {}
    with open(file_name, 'r') as jobs_file:
        for line in jobs_file.readlines():
            line = line.strip()
            items = filter(lambda s: len(s) != 0, re.split(" |\t", line))
            n_items = len(items)
            if n_items < 1 or n_items > 2:
                raise RuntimeError('Invalid job list in file "{}"'.format(args.job_list))
            name = items[0]
            jobs[name] = []
            if n_items > 1:
                id_list = [ int(s) for s in re.split(",", items[1]) ]
                jobs[name] = id_list
    return jobs

jobs = ReadJobs(args.job_list)
if os.path.isfile(args.prev_output):
    prev_jobs = ReadJobs(args.prev_output)
    for job, id_list in prev_jobs.iteritems():
        if job in jobs and len(jobs[job]) == 0:
            jobs[job] = id_list

lumi_per_job_file = "run_and_lumis.tar.gz"
processed_lumis_file = "processedLumis.json"
crab_result_files = [ "inputDatasetLumis.json", "input_dataset_duplicate_lumis.json", "input_dataset_lumis.json",
                      "lumisToProcess.json", "notFinishedLumis.json", processed_lumis_file, lumi_per_job_file ]
jobs_ex = {}
for job,prev_failed_ids in sorted(jobs.iteritems()):
    job_name = re.sub(r'[0-9]+_[0-9]+:[^_]*_', '', job)
    job_path = os.path.join(args.work_area, job_name)
    results_tmp = os.path.join(args.crab_results_out, job_name)
    job_lumis_tmp = os.path.join(results_tmp, 'job_lumis')
    results_out = os.path.join(args.crab_results_out, job_name + '.tar.bz2')
    crab_results_path = os.path.join(job_path, 'results')
    processed_lumis_file_path = os.path.join(results_tmp, processed_lumis_file)
    jobs_ex[job] = prev_failed_ids
    if os.path.exists(results_out):
        print('Skipping job {}, because the file "{}" with the processed' \
              ' CRAB results already exists.'.format(job_name, results_out))
        continue
    print("Processing {}...".format(job_name))
    if os.path.exists(results_tmp):
        print('WARNING: temporary directory "{}" needed to process the job already exists.'\
              ' Removing it...'.format(results_tmp))
        shutil.rmtree(results_tmp)
    if not os.path.isdir(job_path):
        raise RuntimeError('CRAB working area not found for job "{}".'.format(job_name))
    if not os.path.isdir(crab_results_path):
        raise RuntimeError('CRAB results not found for job "{}".'.format(job_name))
    for file in crab_result_files:
        file_path = os.path.join(crab_results_path, file)
        if not os.path.isfile(file_path):
            raise RuntimeError('CRAB results are not complete for job "{}".' \
                               ' File "{}" not found.'.format(job_name, file))

    print("Copying files into temporary area...")
    os.makedirs(job_lumis_tmp)
    for file in crab_result_files:
        file_path = os.path.join(crab_results_path, file)
        if file == lumi_per_job_file:
            tar = tarfile.open(file_path)
            tar.extractall(job_lumis_tmp)
            tar.close()
        else:
            shutil.copy(file_path, results_tmp)
    print("Creating list of not processed jobs...")

    job_lumis_files = os.listdir(job_lumis_tmp)
    failed_job_ids = []
    processed_lumis = LumiList(filename=processed_lumis_file_path)
    for job_lumis_file in job_lumis_files:
        job_id = int(re.sub(r'job_lumis_([0-9]*)\.json', r'\1', job_lumis_file))
        job_lumis = LumiList(filename=os.path.join(job_lumis_tmp, job_lumis_file))
        if len(job_lumis - processed_lumis):
            failed_job_ids.append(job_id)
    failed_job_ids = sorted(failed_job_ids)
    jobs_ex[job] = failed_job_ids
    print("Not processed jobs: {}".format(jobs_ex[job]))

    print("Creating results archive...")
    tar = tarfile.open(results_out, 'w:bz2')
    tar.add(results_tmp, arcname=job_name)
    tar.close()
    print("Removing temporary files...")
    shutil.rmtree(results_tmp)

print('Writting extended list into "{}"...'.format(args.output))
with open(args.output, 'w') as f:
    for key, value in sorted(jobs_ex.iteritems()):
        line = key
        n_jobs = len(value)
        if len(value) > 0:
            line += ' {}'.format(value[0])
            for n in range(1, n_jobs):
                line += ',{}'.format(value[n])
        line += '\n'
        f.write(line)
