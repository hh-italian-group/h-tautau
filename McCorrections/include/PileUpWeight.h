/*! The pile up weight.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {

class PileUpWeight : public IWeightProvider {
public:
    using Event = ntuple::Event;
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    PileUpWeight(const std::string& pu_reweight_file_name, const std::string& hist_name, double _max_available_pu,
                 double _default_pu_weight) :
        max_available_pu(_max_available_pu), default_pu_weight(_default_pu_weight),
        pu_weights(LoadPUWeights(pu_reweight_file_name, hist_name)) { }

    virtual double Get(const Event& event) const override { return GetT(event); }
    virtual double Get(const ntuple::ExpressEvent& event) const override { return GetT(event); }

private:
	template<typename Event>
    double GetT(const Event& event) const
    {
        const double nPU = event.npu;
        const Int_t bin = pu_weights->FindBin(nPU);
        const Int_t maxBin = pu_weights->FindBin(max_available_pu);
        const bool goodBin = bin >= 1 && bin <= maxBin;
        return goodBin ? pu_weights->GetBinContent(bin) : default_pu_weight;
    }

    static HistPtr LoadPUWeights(const std::string& pu_reweight_file_name, const std::string& hist_name)
    {
        auto file = root_ext::OpenRootFile(pu_reweight_file_name);
        return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
    }

private:
    double max_available_pu, default_pu_weight;
    HistPtr pu_weights;
};

class PileUpWeightEx : public IWeightProvider {
public:
    using Event = ntuple::Event;
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    PileUpWeightEx(const std::string& pu_data_file_name, const std::string& pu_mc_file_name,
                 const std::string& cfg_file_name, double _max_available_pu,
                 double _default_pu_weight) :
        max_available_pu(_max_available_pu), default_pu_weight(_default_pu_weight)
    {
        LoadPUWeights(pu_data_file_name, pu_mc_file_name, cfg_file_name);
    }

    virtual double Get(const Event& event) const override { return GetT(event); }
    virtual double Get(const ntuple::ExpressEvent& event) const override { return GetT(event); }

    void SetActiveDataset(const std::string& active_dataset)
    {
        const std::string active_dataset_full_name = "n_pu_mc_" + active_dataset;
        
        if(!datasets.count(active_dataset_full_name))
            throw exception("PileUp weight for dataset '%1%' isn't found") % active_dataset;
        active_group = datasets.at(active_dataset_full_name);
    }

private:
	template<typename Event>
    double GetT(const Event& event) const
    {
        if(!active_group.is_initialized())
             throw exception("active group isn't initialized");
        const double nPU = event.npu;
        auto w_hist = pu_weights.at(*active_group);
        const Int_t bin = w_hist->FindBin(nPU);
        const Int_t maxBin = w_hist->FindBin(max_available_pu);
        const bool goodBin = bin >= 1 && bin <= maxBin;
        return goodBin ? w_hist->GetBinContent(bin) : default_pu_weight;
    }

    void LoadPUWeights(const std::string& pu_data_file_name, const std::string& pu_mc_file_name,
                                              const std::string& cfg_file_name)
    {
        auto mc_pileup_file = root_ext::OpenRootFile(pu_mc_file_name);
        auto data_pileup_file = root_ext::OpenRootFile(pu_data_file_name);
        auto pu_data = std::shared_ptr<TH1D>(root_ext::ReadObject<TH1D>(*data_pileup_file, "pileup"));

        auto data_norm = std::shared_ptr<TH1D>(root_ext::CloneObject(*pu_data));
        const int max_bin = pu_data->FindBin(max_available_pu);

        for(int i = max_bin + 1; i <= pu_data->GetNbinsX(); ++i){
            data_norm->SetBinContent(i,0);
            data_norm->SetBinError(i,0);
        }

        RenormalizeHistogram(*data_norm, 1, true);

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
    }

private:
    double max_available_pu, default_pu_weight;
    std::map<std::string, size_t> datasets;
    std::vector<HistPtr> pu_weights;
    boost::optional<size_t> active_group;
};

} // namespace mc_corrections
} // namespace analysis
