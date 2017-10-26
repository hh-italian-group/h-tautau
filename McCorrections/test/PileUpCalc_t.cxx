/*! TauId weight test.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "h-tautau/McCorrections/include/PileUpWeight.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include <iostream>


struct Arguments {
    run::Argument<std::string> input_weight_file{"input_weight_file", "input weight root file"};
    run::Argument<std::string> input_MC_file{"input_MC_file", "input MC file"};
    run::Argument<std::string> input_data_file{"input_data_file", "input data file"};
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> output_file{"output_file", "Output root file"};
};

namespace analysis {

class PileUpCalcData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(npv, 100, 0, 100)
};

class PileUpCalc_t {
public:
    using PileUpWeight = analysis::mc_corrections::PileUpWeight;

    PileUpCalc_t(const Arguments& _args) : args(_args),
        pu_weight(args.input_weight_file(), "pileup_weight", 60, 0)
    {

    }

public:

    void Run()
    {
        auto inputFile_mc = root_ext::OpenRootFile(args.input_MC_file());
        auto inputFile_data = root_ext::OpenRootFile(args.input_data_file());
        auto eventTuple_mc =
                ntuple::CreateEventTuple(args.tree_name(),inputFile_mc.get(),true,ntuple::TreeState::Full,GetEnabledBranches());
        auto eventTuple_data =
                ntuple::CreateEventTuple(args.tree_name(),inputFile_data.get(),true,ntuple::TreeState::Full,GetEnabledBranches());

        for(const ntuple::Event& event : *eventTuple_mc) {
            anaData.npv("mc").Fill(event.npv, pu_weight.Get(event));
        }

        for(const ntuple::Event& event : *eventTuple_data) {
            anaData.npv("data").Fill(event.npv, pu_weight.Get(event));
        }

        const double chi2_test = anaData.npv("mc").Chi2Test(&anaData.npv("data"), "WW");
        const double ks_test = anaData.npv("mc").KolmogorovTest(&anaData.npv("data"));

        std::cout << "Chi2 test: " << chi2_test << ", KS test: " << ks_test << std::endl;
    }


private:
    Arguments args;
    PileUpWeight pu_weight;
    std::shared_ptr<TFile> output;
    PileUpCalcData anaData;

    static const std::set<std::string>& GetEnabledBranches()
    {
        static const std::set<std::string> EnabledBranches_read = {
            "npv", "npu"
        };
        return EnabledBranches_read;
    }


};

} //namespace analysis

PROGRAM_MAIN(analysis::PileUpCalc_t, Arguments)
