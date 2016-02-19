/*! Definition of RunReport class that produces report about processed runs in JSON format.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <map>
#include <set>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "EventId.h"

namespace analysis {

class RunReport {
public:
    typedef std::pair<unsigned, unsigned> LumiBlockRange;
    typedef std::set<unsigned> LumiBlockSet;
    typedef std::map<unsigned, LumiBlockSet> RunMap;
    typedef std::set<LumiBlockRange> LumiBlockRangeSet;

    RunReport(const std::string& _outputFileName)
        : outputFileName(_outputFileName) {}

    void AddEvent(const EventId& eventId)
    {
        runMap[eventId.runId].insert(eventId.lumiBlock);
    }

    void Report()
    {
        using namespace boost::property_tree;
        ptree report;

        for(auto& runDescriptor : runMap) {
            std::ostringstream runName;
            runName << runDescriptor.first;

            ptree rangesPtree;
            const LumiBlockRangeSet ranges = CreateRanges(runDescriptor.second);
            for(const auto& range : ranges) {
                AddRange(rangesPtree, range.first, range.second);
            }
            report.add_child(runName.str(), rangesPtree);
        }
        std::ofstream output_stream(outputFileName.c_str());
        json_parser::write_json(output_stream, report, false);
    }

private:
    static LumiBlockRangeSet CreateRanges(const LumiBlockSet& lumiBlocks)
    {
        if(!lumiBlocks.size())
            throw std::runtime_error("no lumis found.");
        LumiBlockRangeSet result;
        LumiBlockSet::const_iterator iter = lumiBlocks.begin();
        unsigned startBlock = *iter, previousBlock = *iter;
        ++iter;
        for(; iter != lumiBlocks.end(); ++iter) {
            const unsigned lumiBlock = *iter;
            if(previousBlock != lumiBlock - 1) {
                result.insert(LumiBlockRange(startBlock, previousBlock));
                startBlock = lumiBlock;
            }
            previousBlock = lumiBlock;
        }
        result.insert(LumiBlockRange(startBlock, previousBlock));
        return result;
    }

    static void AddRange(boost::property_tree::ptree& ranges, unsigned minLumiBlock, unsigned maxLumiBlock)
    {
        boost::property_tree::ptree range, minLumiBlockBranch, maxLumiBlockBranch;
        minLumiBlockBranch.put("", minLumiBlock);
        maxLumiBlockBranch.put("", maxLumiBlock);
        range.push_back(std::make_pair("", minLumiBlockBranch));
        range.push_back(std::make_pair("", maxLumiBlockBranch));
        ranges.push_back(std::make_pair("", range));
    }

private:
    std::string outputFileName;
    RunMap runMap;
};

} // analysis
