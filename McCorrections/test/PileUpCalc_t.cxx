/*! TauId weight test.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "h-tautau/McCorrections/include/PileUpWeight.h"
#include <iostream>


struct Arguments {
    run::Argument<std::string> input_MC_file{"input_MC_file", "input MC file"};
    run::Argument<std::string> input_data_file{"input_data_file", "input data file"};
    run::Argument<std::string> input_weight_file{"input_weight_file", "input weight root file"};
    run::Argument<std::string> output_file{"output_file", "Output root file"};
};

namespace analysis {

class PileUpCalc_t {
public:
    using PileUpWeight = analysis::mc_corrections::PileUpWeight;

    PileUpCalc_t(const Arguments& _args) : args(_args),
        tauId_weight(args.json_file(), Parse<DiscriminatorWP>(args.iso_type()))
    {

    }

public:

    void Run()
    {
        auto inputFile = root_ext::OpenRootFile(args.input_file());
        auto eventTuple = ntuple::CreateEventTuple("tauTau",inputFile.get(),true,ntuple::TreeState::Full);

        for(const ntuple::Event& event : *eventTuple) {
            std::cout << "TauID weight: " << tauId_weight.Get(event) << std::endl;
        }

    }


private:
    Arguments args;
    TauIdWeight tauId_weight;

};

} //namespace analysis

PROGRAM_MAIN(analysis::PileUpCalc_t, Arguments)
