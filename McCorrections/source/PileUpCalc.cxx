/*! Merge BSM files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"


struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> MC_input_file{"MC_input_file", "MC input file"};
    run::Argument<std::string> data_pileup_file{"data_pileup_file", "Pileup file for data"};
    run::Argument<std::string> output_weight_file{"output_weight_file", "Output weight root file"};
};

namespace analysis {

class PileUpCalcData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(n_pu_mc, 200, 0, 200)
    TH1D_ENTRY(pileup, 200, 0, 200)
};


class PileUpCalc {
public:
    PileUpCalc(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.output_weight_file())), anaData(output)
    {
        LoadInputs();
    }

    void Run() {}

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    PileUpCalcData anaData;

    void LoadInputs()
    {

        auto inputFile = root_ext::OpenRootFile(args.MC_input_file());
        auto data_pileup_file = root_ext::OpenRootFile(args.data_pileup_file());
        ntuple::ExpressTuple tuple(args.tree_name(), inputFile.get(), true);
        for(const auto& event : tuple)
            anaData.n_pu_mc().Fill(event.npu);

        auto n_pu_data = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*data_pileup_file, "pileup"));
        anaData.pileup().CopyContent(*n_pu_data);
        RenormalizeHistogram(anaData.pileup(), 1, true);
        RenormalizeHistogram(anaData.n_pu_mc(), 1, true);
        anaData.pileup().Divide(&anaData.n_pu_mc());
    }
};

} //namespace analysis

PROGRAM_MAIN(analysis::PileUpCalc, Arguments)
