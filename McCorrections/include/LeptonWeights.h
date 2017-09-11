/*! Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {

namespace detail {
class LeptonScaleFactors {
public:
    using SF = htt_utilities::ScaleFactor;

    LeptonScaleFactors(const std::string& idIsoInput, const std::string& triggerInput) :
        idIso(new SF()), trigger(new SF())
    {
        idIso->init_ScaleFactor(idIsoInput);
        trigger->init_ScaleFactor(triggerInput);
    }

    template<typename LorentzVector>
    double GetIdIsoSF(const LorentzVector& p4) const
    {
        return idIso->get_ScaleFactor(p4.pt(), p4.eta());
    }

    template<typename LorentzVector>
    double GetTriggerSF(const LorentzVector& p4) const
    {
        return trigger->get_ScaleFactor(p4.pt(), p4.eta());
    }

    template<typename LorentzVector>
    double GetTotalSF(const LorentzVector& p4) const { return GetIsoSF(p4) * GetTriggerSF(p4); }

private:
    std::shared_ptr<SF> idIso, trigger;
};

class MuonScaleFactorPOG {
public:
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    MuonScaleFactorPOG(const std::string& idInput_B_F, const std::string& isoInput_B_F,
                       const std::string& triggerInput_B_F, const std::string& idInput_G_H,
                       const std::string& isoInput_G_H, const std::string& triggerInput_G_H) :
        id_hist_B_F(LoadWeight(idInput_B_F,"MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")),
        id_hist_G_H(LoadWeight(idInput_G_H,"MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")),
        iso_hist_B_F(LoadWeight(isoInput_B_F,"LooseISO_TightID_pt_eta/pt_abseta_ratio")),
        iso_hist_G_H(LoadWeight(isoInput_G_H,"LooseISO_TightID_pt_eta/pt_abseta_ratio")),
        trigger_hist_B_F(LoadWeight(triggerInput_B_F,"IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio")),
        trigger_hist_G_H(LoadWeight(triggerInput_G_H,"IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio")),
        lumi_B_F(19.72), lumi_G_H(15.931)
    {
//        file_idInput_B_F = new TFile(idInput_B_F, "read");
    }

    template<typename LorentzVector>
    double GetIdSF(const LorentzVector& p4) const
    {
//        TH2F* hist_B_F = (TH2F*)file_idInput_B_F->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
        const Int_t bin_pt_B_F = id_hist_B_F->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta_B_F = id_hist_B_F->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf_B_F = id_hist_B_F->GetBinContent(bin_pt_B_F,bin_eta_B_F);

        const Int_t bin_pt_G_H = id_hist_G_H->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta_G_H = id_hist_G_H->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf_G_H = id_hist_G_H->GetBinContent(bin_pt_G_H,bin_eta_G_H);

        return ((sf_B_F * lumi_B_F) + (sf_G_H * lumi_G_H))/(lumi_B_F + lumi_G_H);
    }

    template<typename LorentzVector>
    double GetIsoSF(const LorentzVector& p4) const
    {
        const Int_t bin_pt_B_F = iso_hist_B_F->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta_B_F = iso_hist_B_F->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf_B_F = iso_hist_B_F->GetBinContent(bin_pt_B_F,bin_eta_B_F);

        const Int_t bin_pt_G_H = iso_hist_G_H->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta_G_H = iso_hist_G_H->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf_G_H = iso_hist_G_H->GetBinContent(bin_pt_G_H,bin_eta_G_H);

        return ((sf_B_F * lumi_B_F) + (sf_G_H * lumi_G_H))/(lumi_B_F + lumi_G_H);
    }

    template<typename LorentzVector>
    double GetTriggerSF(const LorentzVector& p4) const
    {
        const Int_t bin_pt_B_F = trigger_hist_B_F->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta_B_F = trigger_hist_B_F->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf_B_F = trigger_hist_B_F->GetBinContent(bin_pt_B_F,bin_eta_B_F);

        const Int_t bin_pt_G_H = trigger_hist_G_H->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta_G_H = trigger_hist_G_H->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf_G_H = trigger_hist_G_H->GetBinContent(bin_pt_G_H,bin_eta_G_H);

        return ((sf_B_F * lumi_B_F) + (sf_G_H * lumi_G_H))/(lumi_B_F + lumi_G_H);
    }

    template<typename LorentzVector>
    double GetTotalSF(const LorentzVector& p4) const { return GetIdSF(p4) * GetIsoSF(p4) * GetTriggerSF(p4); }

    static HistPtr LoadWeight(const std::string& weight_file_name, const std::string& hist_name)
    {
        auto file = root_ext::OpenRootFile(weight_file_name);
        return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
    }

private:
    HistPtr id_hist_B_F, id_hist_G_H, iso_hist_B_F, iso_hist_G_H, trigger_hist_B_F, trigger_hist_G_H;
    double lumi_B_F, lumi_G_H;

};

} // namespace detail

class LeptonWeights : public IWeightProvider {
public:
    using Event = ntuple::Event;

    LeptonWeights(const std::string& electron_idIsoInput, const std::string& electron_triggerInput,
                  const std::string& muon_idIsoInput, const std::string& muon_triggerInput) :
        electronSF(electron_idIsoInput, electron_triggerInput),
        muonSF(muon_idIsoInput, muon_triggerInput)
    {
    }

    double GetIdIsoWeight(const Event& event) const
    {
        const Channel channel = static_cast<Channel>(event.channelId);
        if(channel == Channel::ETau) return electronSF.GetIdIsoSF(event.p4_1);
        if(channel == Channel::MuTau) return muonSF.GetIdIsoSF(event.p4_1);
        return 1.0;
    }

    double GetTriggerWeight(const Event& event) const
    {
        const Channel channel = static_cast<Channel>(event.channelId);
        if(channel == Channel::ETau) return electronSF.GetTriggerSF(event.p4_1);
        if(channel == Channel::MuTau) return muonSF.GetTriggerSF(event.p4_1);
        return 1.0;
    }

    virtual double Get(const Event& event) const override { return GetIdIsoWeight(event) * GetTriggerWeight(event); }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
    }

private:
    detail::LeptonScaleFactors electronSF, muonSF;
//    detail::MuonScaleFactorPOG muonSF;
};

} // namespace mc_corrections
} // namespace analysis
