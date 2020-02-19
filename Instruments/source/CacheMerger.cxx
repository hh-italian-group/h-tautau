/*! Merge multiple CacheTuples files into a single file.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <iostream>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Core/include/CacheTuple.h"
#include "h-tautau/Analysis/include/EventCacheProvider.h"

struct Arguments {
    REQ_ARG(std::vector<std::string>, inputs);
    REQ_ARG(std::string, outputFile);
    REQ_ARG(std::string, channel);
};

namespace analysis {

class CacheMerger {
public:
    CacheMerger(const Arguments& _args): args(_args), output(root_ext::CreateRootFile(args.outputFile())),
                                         output_summary("summary", output.get(), false) {}

    void Run()
    {
        std::vector<std::shared_ptr<TFile>> all_cache;
        for (size_t n = 0; n < args.inputs().size(); ++n){
            std::cout << args.inputs().at(n) << "\n";
            all_cache.push_back(root_ext::OpenRootFile(args.inputs().at(n)));
        }

        std::vector<std::shared_ptr<cache_tuple::CacheTuple>> cacheTuples;
        for (size_t n = 0; n < all_cache.size(); ++n){
            auto cacheFile = all_cache.at(n);
            auto cacheTuple = std::make_shared<cache_tuple::CacheTuple>(args.channel(), cacheFile.get(), true);
            cacheTuples.push_back(cacheTuple);
        }
        cache_tuple::CacheTuple cache_out(args.channel(), output.get(), false);

        const Long64_t n_entries = cacheTuples.at(0)->GetEntries();

        for (size_t i = 0; i < cacheTuples.size() ; ++i){
            if(cacheTuples.at(0)->GetEntries() != cacheTuples.at(i)->GetEntries())
                throw exception ("The cache ntuple number '%1%' has an incorrect number of events '%2%' ('%3%').") % i
                                %cacheTuples.at(i)->GetEntries() %cacheTuples.at(0)->GetEntries();
        }

         for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
             EventCacheProvider eventCacheProvider;
             for (size_t i = 0; i < cacheTuples.size(); ++i){
                 cacheTuples.at(i)->GetEntry(current_entry);
                 const cache_tuple::CacheEvent& cache_event = cacheTuples.at(i)->data();
                 eventCacheProvider.AddEvent(cache_event);

            }
            eventCacheProvider.FillEvent(cache_out());
            cache_out.Fill();
         }
         cache_out.Write();

         auto cache_tuple_summary = std::make_shared<cache_ntuple::CacheSummaryTuple>("summary", all_cache.at(0).get(),
                                                                                     true);
         for(const auto& summary : *cache_tuple_summary) {
            output_summary() = summary;
            output_summary.Fill();
        }
         output_summary.Write();
         if(n_entries != cache_out.GetEntries())
            throw exception ("The cache merged output ntuple has '%1%' events, while the input as '%2%'.")
                             %cache_out.GetEntries() %n_entries;

    }
private:
    Arguments args;
    std::shared_ptr<TFile> output;
    cache_ntuple::CacheSummaryTuple output_summary;


};
}// namespace analysis

PROGRAM_MAIN(analysis::CacheMerger, Arguments)
