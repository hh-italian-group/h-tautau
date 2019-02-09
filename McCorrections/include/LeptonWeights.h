/* Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */


#pragma once

#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "WeightProvider.h"
#include "TauTriggerSFs2017.h"
#include "TauIdWeight.h"
#include "TauTriggerWeight.h"

namespace analysis {
namespace mc_corrections {
namespace detail {

class LeptonScaleFactors {
public:
    using SF = htt_utilities::ScaleFactor;

    LeptonScaleFactors(const std::string& idIsoInput, const std::string& triggerInputSingle, const std::string& triggerInputCross = "");

    template<typename LorentzVector>
    double GetIdIsoSF(const LorentzVector& p4) const
    {
        return idIso->get_ScaleFactor(p4.pt(), p4.eta());
    }

    template<typename LorentzVector>
    double GetTriggerEff(const LorentzVector& p4, bool isData) const
    {
        if(isData)
            return triggerSingle->get_EfficiencyData(p4.pt(), p4.eta());
        else
            return triggerSingle->get_EfficiencyMC(p4.pt(), p4.eta());
    }

    template<typename LorentzVector>
    double GetTriggerEffCross(const LorentzVector& p4, bool isData ) const
    {
        if(isData)
            return triggerCross->get_EfficiencyData(p4.pt(), p4.eta());
        else
            return triggerCross->get_EfficiencyMC(p4.pt(), p4.eta());
    }

    bool HasCrossTriggers() const;

private:
    std::shared_ptr<SF> idIso, triggerSingle, triggerCross;
};

class ElectronScaleFactorPOG {
public:
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    ElectronScaleFactorPOG(const std::string& idInput, const std::string& isoInput);

    double GetTriggerSF() const;

    template<typename LorentzVector>
    double GetIsoSF(const LorentzVector& p4) const
    {
        const Int_t bin_eta = iso_hist->GetXaxis()->FindBin(std::abs(p4.eta()));
        const Int_t bin_et = iso_hist->GetYaxis()->FindBin(p4.Et());
        double sf = iso_hist->GetBinContent(bin_eta,bin_et);

        return sf;
    }

    template<typename LorentzVector>
    double GetIdSF(const LorentzVector& p4) const
    {
        const Int_t bin_eta = id_hist->GetXaxis()->FindBin(std::abs(p4.eta()));
        const Int_t bin_et = id_hist->GetYaxis()->FindBin(p4.Et());
        double sf = id_hist->GetBinContent(bin_eta,bin_et);

        return sf;
    }

    template<typename LorentzVector>
    double GetIdIsoSF(const LorentzVector& p4) const { return GetIdSF(p4) * GetIsoSF(p4); }

    static HistPtr LoadWeight(const std::string& weight_file_name, const std::string& hist_name);

private:
    HistPtr id_hist, iso_hist;
 };

 class MuonScaleFactorPOG {
 public:
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    MuonScaleFactorPOG(const std::string& idInput, const std::string& isoInput, const std::string& triggerInput);

    template<typename LorentzVector>
    double GetTriggerSF(const LorentzVector& p4) const
    {
        const Int_t bin_pt = trigger_hist->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta = trigger_hist->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf = trigger_hist->GetBinContent(bin_pt,bin_eta);

        return sf;
    }

    template<typename LorentzVector>
    double GetIsoSF(const LorentzVector& p4) const
    {
        const Int_t bin_pt = iso_hist->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta = iso_hist->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf = iso_hist->GetBinContent(bin_pt,bin_eta);

        return sf;
    }

    template<typename LorentzVector>
    double GetIdSF(const LorentzVector& p4) const
    {
        const Int_t bin_pt = id_hist->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta = id_hist->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf = id_hist->GetBinContent(bin_pt,bin_eta);

        return sf;
    }

    template<typename LorentzVector>
    double GetIdIsoSF(const LorentzVector& p4) const { return GetIdSF(p4) * GetIsoSF(p4); }

    static HistPtr LoadWeight(const std::string& weight_file_name, const std::string& hist_name);

private:
    HistPtr id_hist, iso_hist, trigger_hist;
};

} // namespace detail

class LeptonWeights : public IWeightProvider {
public:
    using Event = ntuple::Event;

    LeptonWeights(const std::string& electron_idIsoInput, const std::string& electron_SingletriggerInput,
                  const std::string& electron_CrossTriggerInput, const std::string& muon_idIsoInput,
                  const std::string& muon_SingletriggerInput, const std::string& muon_CrossTriggerInput,
                  const std::string& tauTriggerInput, const std::string& tauTriggerInputOld, Period period, DiscriminatorWP _tau_iso_wp);

    double GetIdIsoWeight(const Event& event) const;
    double GetTriggerWeight(const Event& event) const;

    virtual double Get(const Event& event) const override;
    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override;

private:
    double GetTriggerEfficiency(const Event& event, bool isData) const;

private:
    detail::LeptonScaleFactors electronSF, muonSF;
    std::shared_ptr<TauIdWeight> tauIdWeight;
    std::shared_ptr<TauTriggerWeight> tauTriggerWeight;
    DiscriminatorWP tau_iso_wp;
};

} // namespace mc_corrections
} // namespace analysis
