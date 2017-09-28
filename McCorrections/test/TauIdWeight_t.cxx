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
        tauId_weight(args.json_file(), Parse<DiscriminatorWP>(args.iso_type()))
    {
        std::cout << "Ciao" << std::endl;
    }

public:

    void Run()
    {
        auto inputFile = root_ext::OpenRootFile(args.input_file());
        auto eventTuple = ntuple::CreateEventTuple("tauTau",inputFile.get(),true,ntuple::TreeState::Full);

        for(const Event& event : *eventTuple) {
            std::cout << "TauID weight: " << tauId_weight.Get(event) << std::endl;
        }

    }


private:
    Arguments args;
    TauIdWeight tauId_weight;

};

} //namespace analysis

PROGRAM_MAIN(analysis::TauIdWeight_t, Arguments)
