#!/usr/bin/env python
# Submit jobs on CRAB.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import argparse
import sys
import re

parser = argparse.ArgumentParser(description='Submit jobs on CRAB.',
                  formatter_class = lambda prog: argparse.HelpFormatter(prog,width=90))
parser.add_argument('--cfg', required=True, dest='cfg', type=str, help="CMSSW configuration file")
parser.add_argument('--site', required=True, dest='site', type=str, help="Site for stage out.")
parser.add_argument('--output', required=True, dest='output', type=str,
                    help="output path after /store/user/USERNAME")
parser.add_argument('job_file', type=str, nargs='+', help="text file with jobs descriptions")
args = parser.parse_args()

from CRABClient.UserUtilities import config, ClientException, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand
from httplib import HTTPException

config = config()

config.General.workArea = 'work_area'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = args.cfg

config.Data.inputDBS = 'global'
config.General.transferOutputs = True
config.General.transferLogs = True
config.Data.publication = False

config.Site.storageSite = args.site
config.Data.outLFNDirBase = "/store/user/{}/{}".format(getUsernameFromSiteDB(), args.output)

from crab_tools import JobCollection

for job_file in args.job_file:
    job_collection = JobCollection(job_file)
    print job_file
    print job_collection
    job_collection.submit(config)
