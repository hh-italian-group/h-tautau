#!/usr/bin/env python
# Create splitting for cache jobs.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
import os
import numpy as np
import pandas

parser = argparse.ArgumentParser(description='Create splitting for cache jobs')
parser.add_argument('--input', type=str, required=True, help="Input csv file with summary")
parser.add_argument('--max-exe-time', type=float, required=True, help="Maximal desired execution time in hours")
args = parser.parse_args()

def create_splits(df, max_exe_time):
    jobs = {}
    n_splits = np.ceil(df.exeTime_h.values / max_exe_time)
    for n in range(df.shape[0]):
        if n_splits[n] > 1:
            job_name = df.index.values[n]
            jobs[job_name] = int(n_splits[n])
    return jobs

df = pandas.read_csv(args.input)
df['exeTime_h'] = df.exeTime / (60. * 60.)
df['job_name'] = df.sample_name + '_' + df.channel
df_gr = df.groupby(['job_name'])

jobs = create_splits(df_gr.sum(), args.max_exe_time)
for job_name in sorted(jobs):
    print('{} {}'.format(job_name, jobs[job_name]))
