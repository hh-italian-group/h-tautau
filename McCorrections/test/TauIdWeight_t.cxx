/*! TauId weight test.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "h-tautau/McCorrections/include/TauIdWeight.h"
#include <iostream>


struct Arguments {
    run::Argument<std::string> input_file{"input_file", "input file"};
    run::Argument<std::string> json_file{"json_file", "json file"};
    run::Argument<std::string> iso_type{"iso_type", "iso type"};
};

namespace analysis {

class TauIdWeight_t {
public:
    using TauIdWeight = analysis::mc_corrections::TauIdWeight;

    TauIdWeight_t(const Arguments& _args) : args(_args),
        tauId_weight(args.json_file(), args.iso_type())
    {
        std::cout << "Ciao" << std::endl;
    }

public:

    void Run()
    {
        std::string filename = args.input_file();
        auto inputFile = root_ext::OpenRootFile(filename);
        ntuple::EventTuple eventTuple("tauTau", inputFile.get(), true, {}, {});
        const Long64_t n_entries = eventTuple.GetEntries();

        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
            eventTuple.GetEntry(current_entry);
            const ntuple::Event& event = eventTuple.data();

            std::cout << "TauID weight: " << tauId_weight.GetWeight(event) << std::endl;
        }

    }


private:
    Arguments args;
    TauIdWeight tauId_weight;

};

} //namespace analysis

PROGRAM_MAIN(analysis::TauIdWeight_t, Arguments)
