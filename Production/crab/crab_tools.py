# Definition of tools for CRAB job submission.
# This file is part of https://github.com/hh-italian-group/h-tautau.

import re
from sets import Set
from CRABClient.UserUtilities import ClientException
from CRABAPI.RawCommand import crabCommand
from httplib import HTTPException

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print str(hte)
        print "\n{}\nERROR: failed to submit task due to HTTPException.\n{}".format(hte, hte.headers)
    except ClientException as cle:
        print "ERROR: failed to submit task due to ClientException.\n{}".format(cle)

class Job:
    def __init__(self, line, jobNameSuffix = '', unitsPerJob = 1000):
        items = filter(lambda s: len(s) != 0, re.split(" |\t", line))
        n_items = len(items)
        if n_items < 2 or n_items > 3:
            raise RuntimeError("invalid job description = '{}'.".format(line))
        self.jobName = items[0]
        self.requestName = self.jobName + jobNameSuffix
        self.unitsPerJob = unitsPerJob
        self.inputDataset = items[1]
        if n_items > 2:
            self.lumiMask = items[2]
        else:
            self.lumiMask = None


    def __str__(self):
        str = "requestName = '{}', unitsPerJob = {}, inputDataset = '{}'".format(self.requestName, self.unitsPerJob,
                                                                                 self.inputDataset)
        if self.lumiMask is not None:
            str += ", lumiMask = '{}'".format(self.lumiMask)
        return str

    def submit(self, config):
        config.General.requestName = self.requestName
        config.Data.inputDataset = self.inputDataset
        config.Data.unitsPerJob = self.unitsPerJob
        if self.lumiMask is not None:
            config.Data.lumiMask = self.lumiMask
        submit(config)

class JobCollection:
    def __init__(self, file_name, job_names = '', lumi_mask = '', jobNameSuffix = '', unitsPerJob = 1000):
        self.jobs = []
        self.jobNames = job_names
        input_file = open(file_name, 'r')
        lines = [ s.strip() for s in input_file.readlines() ]
        lines = filter(lambda s: len(s) != 0 and s[0] != '#', lines)
        if len(lines) <= 2:
            raise RuntimeError("file '{}' is empty".format(file_name))
        header_items = filter(lambda s: len(s) != 0, re.split(" |\n", lines[0]))
        index_line = 0
        if header_items[0].startswith("lumiMask"):
            index_line = 1
            lumi = filter(lambda s: len(s) != 0, re.split("=", header_items[0]))
            self.lumiMask = lumi[1]
        else:
            self.lumiMask =  ''
        self.pyCfgParams = filter(lambda s: len(s) != 0, re.split(" |\t", lines[index_line]))
        #if len(header_items) == 0 or len(header_items) > 1:
        #    raise RuntimeError("invalid jobs file header '{}' in file '{}'".format(lines[0], file_name))
        self.splitting = 'Automatic'
        #known_splittings = Set(['FileBased', 'LumiBased', 'EventAwareLumiBased','Automatic'])
        #if not self.splitting in known_splittings:
        #    raise RuntimeError("unknown splitting = '{}' in file '{}'".format(self.splitting, file_name))

        if len(header_items) > 0:
            if header_items[0].lower() == "signal":
                if len(lines) < 4:
                    raise RuntimeError("invalid signal jobs definition in file '{}'".format(file_name))
                masses = filter(lambda s: len(s) != 0, re.split(" |\t", lines[2]))
                template = lines[3]
                for mass in masses:
                    line = template.format(M = mass)
                    self.jobs.append(Job(line))
                return
        if len(lumi_mask) != 0:
            self.lumiMask = lumi_mask

        for line in lines[2:]:
            self.jobs.append(Job(line, jobNameSuffix,unitsPerJob))
        input_file.close()

    def __str__(self):
        result = "Splitting = '{}', cfgParams = {}, lumiMask = '{}'".format(self.splitting, self.pyCfgParams,
                                                                              self.lumiMask)
        for job in self.jobs:
            if len(self.jobNames) == 0 or job.jobName in self.jobNames:
                result += "\n" + str(job)
        return result

    def submit(self, config):
        config.Data.splitting = self.splitting
        config.JobType.pyCfgParams = self.pyCfgParams
        for job in self.jobs:
            if len(self.jobNames) == 0 or job.jobName in self.jobNames:
                config.Data.lumiMask = self.lumiMask
                job.submit(config)
