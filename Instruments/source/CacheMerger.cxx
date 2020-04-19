/*! Merge multiple CacheTuples files into a single file.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <iostream>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Core/include/CacheTuple.h"
#include "h-tautau/Analysis/include/EventCacheProvider.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"

struct Arguments {
    REQ_ARG(std::string, channel);
    REQ_ARG(std::string, outputFile);
    REQ_ARG(std::vector<std::string>, inputs);
};

namespace analysis {

class CacheMerger {
public:
    CacheMerger(const Arguments& _args) :
            args(_args), cache_reader(args.inputs(), args.channel()),
            output(root_ext::CreateRootFile(args.outputFile())),
            output_cache(args.channel(), output.get(), false), output_summary("summary", output.get(), false),
            progressReporter(10, std::cout)
    {
        output_cache.SetAutoFlush(1000);
        output_cache.SetMaxVirtualSize(10000000);
        progressReporter.SetTotalNumberOfEvents(cache_reader.GetTotalNumberOfEntries());

        std::cout << "Inputs:\n";
        for(const auto& input : args.inputs())
            std::cout << "\t" << input << "\n";
        std::cout << "Output: " << args.outputFile() << std::endl;
    }

    void Run()
    {
        std::cout << "Merging cache summaries... ";
        for(const auto& summary : cache_reader.GetSummary()) {
            output_summary() = summary;
            output_summary.Fill();
        }
        output_summary.Write();

        std::cout << "done.\nMerging caches..." << std::endl;

        size_t n_steps = 0;
        while(true) {
            const auto entry_index = cache_reader.GetCurrentEntryIndex();
            if(!entry_index) break;
            const auto provider = cache_reader.Read(*entry_index);
            if(!provider.IsEmpty()) {
                provider.FillEvent(output_cache());
                output_cache.Fill();
            }
            if(++n_steps % 100 == 0) {
                const size_t n_processed = cache_reader.GetTotalNumberOfEntries()
                                           - cache_reader.GetRemainingNumberOfEntries();
                progressReporter.Report(n_processed, false);
            }
        }
        output_cache.Write();
        progressReporter.Report(cache_reader.GetTotalNumberOfEntries(), true);
    }
private:
    Arguments args;
    EventCacheReader cache_reader;
    std::shared_ptr<TFile> output;
    cache_tuple::CacheTuple output_cache;
    cache_tuple::CacheSummaryTuple output_summary;
    analysis::tools::ProgressReporter progressReporter;
};
}// namespace analysis

PROGRAM_MAIN(analysis::CacheMerger, Arguments)
