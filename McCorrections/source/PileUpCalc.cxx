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
    TH1D_ENTRY(weight_pu, 200, 0, 200)
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
        std::shared_ptr<ntuple::ExpressTuple> AllEventTuple(new ntuple::ExpressTuple(args.tree_name(), inputFile.get(), true));
        const Long64_t n_entries = AllEventTuple->GetEntries();

        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries

            AllEventTuple->GetEntry(current_entry);
            std::shared_ptr<ntuple::ExpressEvent> event(new ntuple::ExpressEvent(AllEventTuple->data()));

            anaData.n_pu_mc().Fill(event->npu);
        } //end loop on entries

        TH1D* n_pu_data = (TH1D*)data_pileup_file->Get("pileup");
        anaData.weight_pu().CopyContent(*n_pu_data);
        anaData.weight_pu().Divide(&anaData.n_pu_mc());
    }
};

} //namespace analysis

PROGRAM_MAIN(analysis::PileUpCalc, Arguments)
