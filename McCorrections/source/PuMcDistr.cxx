/*! Pileup distribution for MC.
This file is part of https://github.com/hh-italian-group/h-tautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"


struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> output_weight_file{"output_weight_file", "Output weight root file"};
    run::Argument<std::vector<std::string>> MC_input_files{"MC_input_files", "MC input files"};
};

namespace analysis {

class PileUpCalcData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(n_pu_mc, 1000, 0, 100)
    TH1D_ENTRY(n_pu_mc_norm, 1000, 0, 100)
};


class PuMcDistr {
public:
    PuMcDistr(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.output_weight_file())), anaData(output)
    {

    }

    void Run()
    {
        for(const auto& file_name : args.MC_input_files()){
            auto inputFile = root_ext::OpenRootFile(file_name);
            std::cout << "File opened: " << file_name << std::endl;
            try{
                ntuple::ExpressTuple tuple(args.tree_name(), inputFile.get(), true);
                for(const auto& event : tuple)
                    anaData.n_pu_mc().Fill(event.npu);

                anaData.n_pu_mc_norm().CopyContent(anaData.n_pu_mc());
                RenormalizeHistogram(anaData.n_pu_mc_norm(), 1, true);
            }catch(std::exception&) {}

        }
    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    PileUpCalcData anaData;
};

} //namespace analysis

PROGRAM_MAIN(analysis::PuMcDistr, Arguments)
