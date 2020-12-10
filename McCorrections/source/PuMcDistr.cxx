/*! Pileup distribution for MC.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"

#include <TKey.h>

struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> output_weight_file{"output_weight_file", "Output weight root file"};
    run::Argument<bool> use_previous_file{"use_previous_file", "Boolean to select if use some input file"};
    run::Argument<std::string> input_weight_file{"input_weight_file", "Optional input file.", ""}; //In case use_previous_file is false, pass any string in the command line
    run::Argument<std::vector<std::string>> MC_input_files{"MC_input_files", "MC input files"};
};

namespace analysis {

class PileUpCalcData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(n_pu_mc, 100, 0, 100)
    TH1D_ENTRY(n_pu_mc_norm, 100, 0, 100)
};

class PuMcDistr {
public:
    PuMcDistr(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.output_weight_file())), anaData(output)
    {
        //Creates list of sample names contained in the input file
        if(args.use_previous_file()){
            prev_input = root_ext::OpenRootFile(args.input_weight_file());
            TIter nextkey(prev_input->GetListOfKeys());
            for(TKey* t_key; (t_key = dynamic_cast<TKey*>(nextkey()));){
                std::string dir_name = t_key->GetName();
                if(dir_name.substr(0,13) == "n_pu_mc_norm_") continue;
                if(dir_name.substr(0,8) == "n_pu_mc_"){
                    dir_name.erase(0,8);
                    prev_hists.push_back(dir_name);
                }
            }
        }
    }

    void Run()
    {
        for(const auto& file_name : args.MC_input_files()){
            if(file_name.find("Embedding") != std::string::npos) continue;
            auto inputFile = root_ext::OpenRootFile(file_name);
            std::cout << "File opened: " << file_name << std::endl;
            auto input_file = GetFileNameWithoutPath(file_name);
            auto name = RemoveFileExtension(input_file);
            std::cout << "name : " << name << '\n';

            //Check if exist the histogram in the input file and copy it to the new output file
            if(prev_hists.size() > 0 && std::count(prev_hists.begin(), prev_hists.end(), name)){
                std::cout << "histogram for sample " << name  << " already found. Copying histogram to new output file. "<< '\n';

                std::ostringstream n_pu_mc;
                n_pu_mc << "n_pu_mc_" << name;
                auto pu_mc = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*prev_input, n_pu_mc.str()));
                anaData.n_pu_mc(name).CopyContent(*pu_mc);

                std::ostringstream n_pu_mc_norm;
                n_pu_mc_norm << "n_pu_mc_norm_" << name;
                auto pu_mc_norm = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*prev_input, n_pu_mc_norm.str()));
                anaData.n_pu_mc_norm(name).CopyContent(*pu_mc_norm);
            }
            else{
                try{
                    ntuple::ExpressTuple tuple(args.tree_name(), inputFile.get(), true);
                    for(const auto& event : tuple)
                        anaData.n_pu_mc(name).Fill(event.npu);

                    anaData.n_pu_mc_norm(name).CopyContent(anaData.n_pu_mc(name));
                    RenormalizeHistogram(anaData.n_pu_mc_norm(name), 1, true);
                }catch(std::exception&) {}
            }
        }
    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    PileUpCalcData anaData;
    std::shared_ptr<TFile> prev_input;
    std::vector<std::string> prev_hists;
};

} //namespace analysis

PROGRAM_MAIN(analysis::PuMcDistr, Arguments)
