/*! The pile up weight.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {

class PileUpWeight : public IWeightProvider {
public:
    using Event = ntuple::Event;
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    PileUpWeight(const std::string& pu_reweight_file_name, const std::string& hist_name, double _max_available_pu,
                 double _default_pu_weight);

    virtual double Get(EventInfoBase& eventInfo) const override;
    virtual double Get(const ntuple::ExpressEvent& event) const override;

private:
    double Get(double nPU) const;
    static HistPtr LoadPUWeights(const std::string& pu_reweight_file_name, const std::string& hist_name);

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
                 double _default_pu_weight);

    virtual double Get(EventInfoBase& eventInfo) const override;
    virtual double Get(const ntuple::ExpressEvent& event) const override;

    void SetActiveDataset(const std::string& active_dataset);

private:
    double Get(double nPU) const;
    void LoadPUWeights(const std::string& pu_data_file_name, const std::string& pu_mc_file_name,
                       const std::string& cfg_file_name);

private:
    double max_available_pu, default_pu_weight;
    std::map<std::string, size_t> datasets;
    std::vector<HistPtr> pu_weights;
    boost::optional<size_t> active_group;
};

} // namespace mc_corrections
} // namespace analysis
