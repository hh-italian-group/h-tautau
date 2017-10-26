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
    run::Argument<std::string> weight_hist_name{"weight_hist_name", "weight histogram name"};
    run::Argument<std::string> input_path{"input_path", "input path of files"};
    run::Argument<std::vector<std::string>> input_MC_files{"input_MC_files", "input MC files"};
    run::Argument<std::vector<std::string>> input_data_files{"input_data_files", "input data files"};
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
        pu_weight(args.input_weight_file(), args.weight_hist_name(), 60, 0)
    {

    }

public:

    void Run()
    {
        for(const auto& file_name : args.input_MC_files()){
            const std::string filename = args.input_path()  + "/" + file_name;
            auto inputFile_mc = root_ext::OpenRootFile(filename);
            ntuple::EventTuple eventTuple_mc(args.tree_name(), inputFile_mc.get(), true, {}, GetEnabledBranches());
            for(const ntuple::Event& event : eventTuple_mc) {
                if (event.eventEnergyScale != 0 ) continue;
                anaData.npv("mc").Fill(event.npv, pu_weight.Get(event));
            }
        }

        for(const auto& file_name : args.input_data_files()){
            const std::string filename = args.input_path()  + "/" + file_name;
            auto inputFile_data = root_ext::OpenRootFile(filename);
            ntuple::EventTuple eventTuple_data(args.tree_name(), inputFile_data.get(), true, {}, GetEnabledBranches());
            for(const ntuple::Event& event : eventTuple_data) {
                anaData.npv("data").Fill(event.npv);
            }
        }

        anaData.npv("mc_norm").CopyContent(anaData.npv("mc"));
        RenormalizeHistogram(anaData.npv("mc_norm"), 1, true);

        anaData.npv("data_norm").CopyContent(anaData.npv("data"));
        RenormalizeHistogram(anaData.npv("data_norm"), 1, true);

        const double chi2_test = anaData.npv("mc_norm").Chi2Test(&anaData.npv("data_norm"), "WW");
        const double ks_test = anaData.npv("mc_norm").KolmogorovTest(&anaData.npv("data_norm"), "X");

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
            "npv", "npu", "eventEnergyScale"
        };
        return EnabledBranches_read;
    }


};

} //namespace analysis

PROGRAM_MAIN(analysis::PileUpCalc_t, Arguments)
