#!/usr/bin/env python
# Iteratively merge cache files using CacheMerger.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
import math
import os
import re
import shutil
import subprocess

parser = argparse.ArgumentParser(description='Iteratively merge cache files using CacheMerger.')
parser.add_argument('--channel', required=True, type=str, help="Analysis channel")
parser.add_argument('--output', required=True, type=str, help="output file name")
parser.add_argument('--max-files-per-step', required=False, type=int, default=2,
                    help="maximal number of inputs to merge per iteration")
parser.add_argument('--cache-merger', required=False, type=str, default='../build/h-tautau/CacheMerger',
                    help="path to the CacheMerger executable")
parser.add_argument('input', nargs='+', type=str, help="list of input cache files")
args = parser.parse_args()

class Job:
    def __init__(self, inputs, output):
        self.inputs = inputs
        self.output = output

class Iteration:
    def __init__(self, inputs, final_output, max_files_per_step, iteration_index):
        self.jobs = []
        if len(inputs) <= 1:
            raise RuntimeError("Inputs length should be > 1.")
        if max_files_per_step <= 1:
            raise RuntimeError("maximal number of inputs to merge per iteration should be > 1.")
        n_splits = int(math.ceil(len(inputs) / float(max_files_per_step)))
        output_path, final_output_basename = os.path.split(final_output)
        output_name_prefix = os.path.splitext(final_output_basename)[0]
        for n in range(n_splits):
            if n_splits > 1:
                job_output = '{}_iter{}_part{}.root'.format(output_name_prefix, iteration_index, n)
            else:
                job_output = final_output
            self.jobs.append(Job(inputs[n * max_files_per_step : (n + 1) * max_files_per_step], job_output))

def CreateJobs(inputs, final_output, max_files_per_step):
    iterations = []
    while len(iterations) == 0 or len(iterations[-1].jobs) > 1:
        iter_idx = len(iterations)
        if iter_idx == 0:
            current_inputs = inputs
        else:
            current_inputs = [ job.output for job in iterations[-1].jobs ]
        iterations.append(Iteration(current_inputs, final_output, max_files_per_step, iter_idx))
    return iterations

iterations = CreateJobs(sorted(args.input), args.output, args.max_files_per_step)

for iter_idx, iteration in enumerate(iterations):
    print("Starting merge iteration {} of {}".format(iter_idx + 1, len(iterations)))
    for job_idx, job in enumerate(iteration.jobs):
        print("Starting merging part {} of {} for iteration {}".format(job_idx + 1, len(iteration.jobs), iter_idx + 1))
        if len(job.inputs) > 1:
            job_cmd = '{} {} {}'.format(args.cache_merger, args.channel, job.output)
            for input in job.inputs:
                job_cmd += ' {} '.format(input)
            print("> {}".format(job_cmd))
            result = subprocess.call([job_cmd], shell=True)
            if result != 0:
                raise RuntimeError("Failed to merge part {} for iteration {}".format(job_idx, iter_idx))
        else:
            if iter_idx == 0:
                shutil.copyfile(job.inputs[0], job.output)
            else:
                shutil.move(job.inputs[0], job.output)
    if iter_idx != 0:
        print("Removing temporary files from iteration {}...".format(iter_idx))
        prev_iteration = iterations[iter_idx - 1]
        for job in prev_iteration.jobs:
            if os.path.exists(job.output):
                os.remove(job.output)

print("All files are successfully merged.")
