#!/usr/bin/env python
# Create json files to split dataset into several parts.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
import re
import subprocess
import json
import sys
from LumiList import LumiList
from multiprocessing import Process, Queue

parser = argparse.ArgumentParser(description='Create json files to split dataset into several parts.',
                  formatter_class = lambda prog: argparse.HelpFormatter(prog,width=90))
parser.add_argument('job_file', type=str, nargs='+', help="json files")
args = parser.parse_args()

lumi1 = LumiList(filename=args.job_file[0])
lumi2 = LumiList(filename=args.job_file[1])

lumi = lumi1 & lumi2

print lumi
