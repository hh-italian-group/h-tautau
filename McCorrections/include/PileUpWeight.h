/*! The pile up weight.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"


namespace analysis {
namespace mc_corrections {

class PileUpWeight {
public:
    using Event = ntuple::Event;
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    PileUpWeight(const std::string& pu_reweight_file_name, const std::string& hist_name, double _max_available_pu,
                 double _default_pu_weight) :
        max_available_pu(_max_available_pu), default_pu_weight(_default_pu_weight),
        pu_weights(LoadPUWeights(pu_reweight_file_name, hist_name)) { }

    double Get(const Event& event) const
    {
        const double nPU = event.npu;
        const Int_t bin = pu_weights->FindBin(nPU);
        const Int_t maxBin = pu_weights->FindBin(max_available_pu);
        const bool goodBin = bin >= 1 && bin <= maxBin;
        return goodBin ? pu_weights->GetBinContent(bin) : default_pu_weight;
    }

private:
    static HistPtr LoadPUWeights(const std::string& pu_reweight_file_name, const std::string& hist_name)
    {
        auto file = root_ext::OpenRootFile(pu_reweight_file_name);
        return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
    }

private:
    double max_available_pu, default_pu_weight;
    HistPtr pu_weights;
};

} // namespace mc_corrections
} // namespace analysis
