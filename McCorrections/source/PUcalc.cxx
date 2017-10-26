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
    run::Argument<std::string> MC_input_file{"MC_input_file", "MC input file"};
    run::Argument<std::string> data_pileup_file{"data_pileup_file", "Pileup file for data"};
    run::Argument<int> max_n_pu{"max_n_pu","maximum value for n_pu"};
    run::Argument<std::string> output_weight_file{"output_weight_file", "Output weight root file"};

};

namespace analysis {

class PileUpCalcData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    ANA_DATA_ENTRY(TH1D, pileup)
};


class PileUpCalc {
public:
    PileUpCalc(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.output_weight_file())), anaData(output)
    {

    }

    void Run()
    {

        auto mc_pileup_file = root_ext::OpenRootFile(args.MC_input_file());
        auto data_pileup_file = root_ext::OpenRootFile(args.data_pileup_file());
        auto pu_data = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*data_pileup_file, "pileup"));
        anaData.pileup.SetMasterHist(pu_data->GetNbinsX(),pu_data->GetXaxis()->GetBinLowEdge(1),
                                     pu_data->GetXaxis()->GetBinUpEdge(pu_data->GetNbinsX()));

        auto pu_mc = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*mc_pileup_file, "n_pu_mc"));
        pu_mc->Rebin(10);//only pileupdata obs
        anaData.pileup("mc").CopyContent(*pu_mc);
        anaData.pileup("mc_norm").CopyContent(anaData.pileup("mc"));
        const int max_bin = anaData.pileup().FindBin(args.max_n_pu());
        anaData.pileup("data").CopyContent(*pu_data);
        anaData.pileup("data_norm").CopyContent(anaData.pileup("data"));

        for(unsigned i = 0; i < anaData.pileup().GetNbinsX(); ++i){
            if(i<=max_bin) continue;
            anaData.pileup("mc_norm").SetBinContent(i,0);
            anaData.pileup("data_norm").SetBinContent(i,0);
            anaData.pileup("mc_norm").SetBinError(i,0);
            anaData.pileup("data_norm").SetBinError(i,0);
        }

        const double integral_mc = anaData.pileup("mc_norm").Integral(1,max_bin);
        if (integral_mc == 0)
            throw analysis::exception("Integral mc is zero.");
        anaData.pileup("mc_norm").Scale(1 /integral_mc);

        anaData.pileup("mc_norm_full").CopyContent(anaData.pileup("mc"));
        RenormalizeHistogram(anaData.pileup("mc_norm_full"), 1, true);


        const double integral_data = anaData.pileup("data_norm").Integral(1,max_bin);
        if (integral_data == 0)
            throw analysis::exception("Integral data is zero.");
        anaData.pileup("data_norm").Scale(1 /integral_data);

        anaData.pileup("data_norm_full").CopyContent(anaData.pileup("data"));
        RenormalizeHistogram(anaData.pileup("data_norm_full"), 1, true);

        anaData.pileup("weight").CopyContent(anaData.pileup("data_norm"));

        for(unsigned i = 0; i < anaData.pileup("weight").GetNbinsX(); ++i){
            if(i<=max_bin) continue;
            anaData.pileup("weight").SetBinContent(i,0);
            anaData.pileup("weight").SetBinError(i,0);
        }

        anaData.pileup("weight").Divide(&anaData.pileup("mc_norm"));

    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    PileUpCalcData anaData;
};

} //namespace analysis

PROGRAM_MAIN(analysis::PileUpCalc, Arguments)
