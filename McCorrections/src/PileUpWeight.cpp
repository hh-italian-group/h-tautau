/*! The pile up weight.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/McCorrections/include/PileUpWeight.h"

#include <fstream>
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {
namespace mc_corrections {

PileUpWeight::PileUpWeight(const std::string& pu_reweight_file_name, const std::string& hist_name,
                           double _max_available_pu, double _default_pu_weight) :
    max_available_pu(_max_available_pu), default_pu_weight(_default_pu_weight),
    pu_weights(LoadPUWeights(pu_reweight_file_name, hist_name))
{
}

PileUpWeight::HistPtr PileUpWeight::LoadPUWeights(const std::string& pu_reweight_file_name,
                                                  const std::string& hist_name)
{
    auto file = root_ext::OpenRootFile(pu_reweight_file_name);
    return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
}

double PileUpWeight::Get(EventInfo& eventInfo) const
{
    return Get(eventInfo->npu);
}
double PileUpWeight::Get(const ntuple::ExpressEvent& event) const { return Get(event.npu); }

double PileUpWeight::Get(double nPU) const
{
    const Int_t bin = pu_weights->FindBin(nPU);
    const Int_t maxBin = pu_weights->FindBin(max_available_pu);
    const bool goodBin = bin >= 1 && bin <= maxBin;
    return goodBin ? pu_weights->GetBinContent(bin) : default_pu_weight;
}

PileUpWeightEx::PileUpWeightEx(const std::string& pu_data_file_name, const std::string& pu_data_file_up_name,
                               const std::string& pu_data_file_down_name, const std::string& pu_mc_file_name,
                               const std::string& cfg_file_name, double _max_available_pu,
                               double _default_pu_weight) :
    max_available_pu(_max_available_pu), default_pu_weight(_default_pu_weight)
{
    data_files[UncertaintyScale::Central] = pu_data_file_name;
    data_files[UncertaintyScale::Up] = pu_data_file_up_name;
    data_files[UncertaintyScale::Down] = pu_data_file_down_name;

    LoadPUWeights(pu_mc_file_name, cfg_file_name);
}

void PileUpWeightEx::SetActiveDataset(const std::string& active_dataset)
{
    const std::string active_dataset_full_name = "n_pu_mc_" + active_dataset;

    if(!datasets.count(active_dataset_full_name))
        throw exception("PileUp weight for dataset '%1%' isn't found") % active_dataset;
    active_group = datasets.at(active_dataset_full_name);
}

void PileUpWeightEx::LoadPUWeights(const std::string& pu_mc_file_name, const std::string& cfg_file_name)
{
    auto mc_pileup_file = root_ext::OpenRootFile(pu_mc_file_name);

    std::ifstream cfg (cfg_file_name);
    if(cfg.fail())
        throw exception("Failed to open config file '%1%'.") % cfg_file_name;
    std::vector<std::string> lines;
    while (cfg.good()){
        std::string line;
        std::getline(cfg, line);
        if(!line.empty())
            lines.push_back(line);
    }

    for(const auto& [unc_scale, data_file] : data_files) {
        auto data_pileup_file = root_ext::OpenRootFile(data_file);
        auto pu_data = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*data_pileup_file, "pileup"));
        auto data_norm = std::shared_ptr<TH1D>(root_ext::CloneObject(*pu_data));
        const int max_bin = pu_data->FindBin(max_available_pu);

        for(int i = max_bin + 1; i <= pu_data->GetNbinsX(); ++i){
            data_norm->SetBinContent(i,0);
            data_norm->SetBinError(i,0);
        }

        RenormalizeHistogram(*data_norm, 1, true);

        std::vector<HistPtr> pu_weights;
        for(size_t id = 0; id < lines.size(); ++id) {
            const std::string& line = lines.at(id);
            auto vector_samples = SplitValueList(line, false, " ");
            if(vector_samples.size() == 0)
                throw exception("The line is empty");
            auto pu_mc = std::make_shared<TH1D>("", "", pu_data->GetNbinsX(),pu_data->GetXaxis()->GetBinLowEdge(1),
                                                  pu_data->GetXaxis()->GetBinUpEdge(pu_data->GetNbinsX()));
            for (size_t i = 0; i < vector_samples.size(); i++) {
                datasets[vector_samples.at(i)] = id;

                auto hist =  std::shared_ptr<TH1>(root_ext::ReadObject<TH1D>(*mc_pileup_file, vector_samples.at(i)));
                pu_mc->Add(&(*hist), 1);
            }
            auto mc_norm = std::shared_ptr<TH1D>(root_ext::CloneObject(*pu_mc));
            for(int i = max_bin + 1; i <= pu_data->GetNbinsX(); ++i){
                mc_norm->SetBinContent(i,0);
                mc_norm->SetBinError(i,0);
            }

            RenormalizeHistogram(*mc_norm, 1, true);
            auto weight = std::shared_ptr<TH1D>(root_ext::CloneObject(*data_norm));
            weight->Divide(&(*mc_norm));
            pu_weights.emplace_back(weight);
        }
        pu_weights_map[unc_scale] = pu_weights;
    }
}

double PileUpWeightEx::Get(EventInfo& eventInfo) const
{
    return Get(eventInfo->npu);
}
double PileUpWeightEx::Get(const ntuple::ExpressEvent& event) const { return Get(event.npu); }

double PileUpWeightEx::Get(const ntuple::ExpressEvent& event, UncertaintyScale unc_scale) const { return Get(event.npu, unc_scale); }

double PileUpWeightEx::Get(double nPU, UncertaintyScale unc_scale) const
{
    if(!active_group.is_initialized())
         throw exception("active group isn't initialized");
    auto w_hist = pu_weights_map.at(unc_scale).at(*active_group);
    const Int_t bin = w_hist->FindBin(nPU);
    const Int_t maxBin = w_hist->FindBin(max_available_pu);
    const bool goodBin = bin >= 1 && bin <= maxBin;
    return goodBin ? w_hist->GetBinContent(bin) : default_pu_weight;
}

double PileUpWeightEx::Get(EventInfo& eventInfo, UncertaintyScale unc_scale) const
{
    return Get(eventInfo->npu, unc_scale);
}

} // namespace mc_corrections
} // namespace analysis
