/*! Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */


#pragma once

#include "AnalysisTools/Core/include/EventIdentifier.h"
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

class ElectronScaleFactorPOG {
public:
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    ElectronScaleFactorPOG(const std::string& idInput, const std::string& isoInput) :
        id_hist(LoadWeight(idInput,"EGamma_SF2D")),
        iso_hist(LoadWeight(isoInput,"EGamma_SF2D"))
        {}

    double GetTriggerSF() const { return 0.991; }

    template<typename LorentzVector>
    double GetIsoSF(const LorentzVector& p4) const
    {
        const Int_t bin_eta = iso_hist->GetXaxis()->FindBin(std::abs(p4.eta()));
        const Int_t bin_et = iso_hist->GetYaxis()->FindBin(p4.Et());
        double sf = iso_hist->GetBinContent(bin_eta,bin_et);

         std::cout << "Iso (EleSF function) : " << sf << '\n';

        return sf;
    }

    template<typename LorentzVector>
    double GetIdSF(const LorentzVector& p4) const
    {
        const Int_t bin_eta = id_hist->GetXaxis()->FindBin(std::abs(p4.eta()));
        const Int_t bin_et = id_hist->GetYaxis()->FindBin(p4.Et());
        double sf = id_hist->GetBinContent(bin_eta,bin_et);
         std::cout << "Id (EleSF function) : " << sf << '\n';
        return sf;
    }

    template<typename LorentzVector>
    double GetIdIsoSF(const LorentzVector& p4) const {
//         std::cout << "GetIdSF(p4) * GetIsoSF(p4) : " << GetIdSF(p4) * GetIsoSF(p4) << '\n';
        return GetIdSF(p4) * GetIsoSF(p4); }

    static HistPtr LoadWeight(const std::string& weight_file_name, const std::string& hist_name)
    {
        auto file = root_ext::OpenRootFile(weight_file_name);
        return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
    }

private:
    HistPtr id_hist, iso_hist, trigger_hist;
 };

 class MuonScaleFactorPOG {
 public:
     using Hist = TH1;
     using HistPtr = std::shared_ptr<Hist>;

     MuonScaleFactorPOG(const std::string& idInput, const std::string& isoInput, const std::string& triggerInput) :
         id_hist(LoadWeight(idInput,"NUM_MediumID_DEN_genTracks_pt_abseta")),
         iso_hist(LoadWeight(isoInput,"NUM_TightRelIso_DEN_MediumID_pt_abseta")), //waht point
         trigger_hist(LoadWeight( triggerInput,"IsoMu27_PtEtaBins/pt_abseta_ratio"))
         {}

     template<typename LorentzVector>
     double GetTriggerSF(const LorentzVector& p4) const
     {

         const Int_t bin_pt = trigger_hist->GetXaxis()->FindBin(p4.pt());
         const Int_t bin_eta = trigger_hist->GetYaxis()->FindBin(std::abs(p4.eta()));
         double sf = trigger_hist->GetBinContent(bin_pt,bin_eta);

         std::cout << " Trigger (MuonSF function): " << sf << '\n';

         return sf;
     }

     template<typename LorentzVector>
     double GetIsoSF(const LorentzVector& p4) const
     {
         const Int_t bin_pt = iso_hist->GetXaxis()->FindBin(p4.pt());
         const Int_t bin_eta = iso_hist->GetYaxis()->FindBin(std::abs(p4.eta()));
         double sf = iso_hist->GetBinContent(bin_pt,bin_eta);

         std::cout << "Iso (MuonSF function) : " << sf << '\n';
         return sf;
     }

     template<typename LorentzVector>
     double GetIdSF(const LorentzVector& p4) const
     {
         const Int_t bin_pt = id_hist->GetXaxis()->FindBin(p4.pt());
         const Int_t bin_eta = id_hist->GetYaxis()->FindBin(std::abs(p4.eta()));
         double sf = id_hist->GetBinContent(bin_pt,bin_eta);
          std::cout << "Id (MuonSF function) : " << sf << '\n';
         return sf;
     }

     template<typename LorentzVector>
     double GetIdIsoSF(const LorentzVector& p4) const {
          std::cout << "GetIdSF(p4) * GetIsoSF(p4) (MuonSF function) : " << GetIdSF(p4) * GetIsoSF(p4) << '\n';
         return GetIdSF(p4) * GetIsoSF(p4); }

     static HistPtr LoadWeight(const std::string& weight_file_name, const std::string& hist_name)
     {
         auto file = root_ext::OpenRootFile(weight_file_name);
         return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
     }

     private:
     HistPtr id_hist, iso_hist, trigger_hist;
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
    LeptonWeights(const std::string& electron_idIsoInput, const std::string& electron_triggerInput2,
                  const std::string& muon_idIsoInput, const std::string& muon_triggerInput2,
                  const std::string& electron_idInput, const std::string& electron_IsoInput,
                  const std::string& muon_idInput, const std::string& muon_IsoInput, const std::string& muon_triggerInput) :
        electronSF(electron_idIsoInput, electron_triggerInput2),
        muonSF(muon_idIsoInput, muon_triggerInput2),
        electronSFPOG(std::make_shared<detail::ElectronScaleFactorPOG>(electron_idInput, electron_IsoInput)),
        muonSFPOG(std::make_shared<detail::MuonScaleFactorPOG>(muon_idInput, muon_IsoInput, muon_triggerInput))
    {
    }

    double GetIdIsoWeight(const Event& event) const
    {
        try {
            const Channel channel = static_cast<Channel>(event.channelId);

            if(channel == Channel::ETau) {
                if(!electronSFPOG || (std::abs(event.p4_2.eta()) < 2.1)){
                    // if(std::abs(event.p4_2.eta()) < 2.1){
                  // if(electronSFPOG == nullptr){
                        std::cout << "eta (Id Iso eTau): " <<  std::abs(event.p4_2.eta()) << '\n';
                        std::cout << "Et(Id Iso eTau) : " << event.p4_2.Et() << '\n';
                        std::cout << "SF (Id Iso eTau) for eta < 2.1: " <<  electronSF.GetIdIsoSF(event.p4_1) << '\n';
                        return electronSF.GetIdIsoSF(event.p4_1);
                   // }
                }

                     if ((std::abs(event.p4_2.eta()) >= 2.10) ){
                    // else{
                        std::cout << "eta (Id Iso eTau): " <<  std::abs(event.p4_2.eta()) << '\n';
                        std::cout << "Et(Id Iso eTau) : " << event.p4_2.Et() << '\n';
                        std::cout << "SF (Id Iso eTau) for eta > 2.1: " <<  electronSF.GetIdIsoSF(event.p4_1) << '\n';
                        return electronSFPOG->GetIdIsoSF(event.p4_1);
                  // }
                }
            }

            if(channel == Channel::MuTau) {
                // if(!electronSFPOG){
                    if(!muonSFPOG || (std::abs(event.p4_2.eta()) < 2.1)){

                        std::cout << "eta (Id Iso muTau): " <<  std::abs(event.p4_2.eta()) << '\n';
                        std::cout << "Et(Id Iso muTau) : " << event.p4_2.Et() << '\n';
                        std::cout << "SF (Id Iso muTau) for eta < 2.1: " <<  electronSF.GetIdIsoSF(event.p4_1) << '\n';
                        return muonSF.GetIdIsoSF(event.p4_1);
                   // }
                }
                // else if(muonSFPOG  != nullptr){
                     if ((std::abs(event.p4_2.eta()) >= 2.10) ){
                 // if (electronSFPOG){
                // else{
                        std::cout << "eta (Id Iso muTau): " <<  std::abs(event.p4_2.eta()) << '\n';
                        std::cout << "Et(Id Iso muTau) : " << event.p4_2.Et() << '\n';
                        std::cout << "SF (Id Iso muTau) for eta > 2.1: " <<  electronSF.GetIdIsoSF(event.p4_1) << '\n';
                        return muonSFPOG->GetIdIsoSF(event.p4_1);
                  // }
                }
            }
        } catch(std::runtime_error& e) {
            const EventIdentifier id(event);
            std::cerr << id << " ERROR: " << e.what() << std::endl;
        }
        return 1.0;
    }

    double GetTriggerWeight(const Event& event) const
    {
        try {
            const Channel channel = static_cast<Channel>(event.channelId);
            if(channel == Channel::ETau) {

                    if( !electronSFPOG || std::abs(event.p4_2.eta()) < 2.1){
                //    if(electronSFPOG == nullptr){
                        std::cout << "eta (Trigger eTau): " <<  std::abs(event.p4_2.eta()) << '\n';
                        std::cout << "Et(Trigger eTau) : " << event.p4_2.Et() << '\n';
                        std::cout << "SF (Trigger eTau) for eta  < 2.1: " <<  electronSF.GetIdIsoSF(event.p4_1) << '\n';
                        return electronSF.GetTriggerSF(event.p4_1);
                   }

                // else if (electronSFPOG  != nullptr) {
                     if ( std::abs(event.p4_2.eta()) >= 2.10 ){
                  // if (electronSFPOG){
                  // else{
                        std::cout << "eta (Trigger eTau): " <<  std::abs(event.p4_2.eta()) << '\n';
                        std::cout << "Et(Trigger eTau) : " << event.p4_2.Et() << '\n';
                        std::cout << "SF (Trigger eTau) for eta  > 2.1: " <<  electronSF.GetIdIsoSF(event.p4_1) << '\n';
                        return electronSFPOG->GetTriggerSF();
                   // }

                }
            }
            if(channel == Channel::MuTau){
                 // if(!electronSFPOG){
                    if( !muonSFPOG || std::abs(event.p4_2.eta()) < 2.1){
                //    if(muonSFPOG == nullptr){
                        std::cout << "eta (Trigger muTau): " <<  std::abs(event.p4_2.eta()) << '\n';
                        std::cout << "Et(Trigger muTau) : " << event.p4_2.Et() << '\n';
                        std::cout << "SF (Trigger muTau) for eta  < 2.1: " <<  electronSF.GetIdIsoSF(event.p4_1) << '\n';
                        return muonSF.GetTriggerSF(event.p4_1);
                    // }
                }

                      if ( std::abs(event.p4_2.eta()) >= 2.10 ){
                  // if (electronSFPOG){
                   // else{
                        std::cout << "eta (Trigger muTau): " <<  std::abs(event.p4_2.eta()) << '\n';
                        std::cout << "Et(Trigger muTau) : " << event.p4_2.Et() << '\n';
                        std::cout << "SF (Trigger muTau) for eta  > 2.1: " <<  electronSF.GetIdIsoSF(event.p4_1) << '\n';
                        return muonSFPOG->GetTriggerSF(event.p4_1);

                }

         }
        } catch(std::runtime_error& e) {
            const EventIdentifier id(event);
            std::cerr << id << " ERROR: " << e.what() << std::endl;
        }
        return 1.0;
    }

    virtual double Get(const Event& event) const override { return GetIdIsoWeight(event) * GetTriggerWeight(event); }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
    }

private:
    detail::LeptonScaleFactors electronSF, muonSF;
    std::shared_ptr<detail::ElectronScaleFactorPOG> electronSFPOG;
    std::shared_ptr<detail::MuonScaleFactorPOG> muonSFPOG;
};

} // namespace mc_corrections
} // namespace analysis




// #pragma once
//
// #include "AnalysisTools/Core/include/EventIdentifier.h"
// #include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
// #include "h-tautau/Analysis/include/AnalysisTypes.h"
// #include "WeightProvider.h"
//
// namespace analysis {
// namespace mc_corrections {
//
// namespace detail {
// class LeptonScaleFactors {
// public:
//     using SF = htt_utilities::ScaleFactor;
//
//     LeptonScaleFactors(const std::string& idIsoInput, const std::string& triggerInput) :
//         idIso(new SF()), trigger(new SF())
//     {
//         idIso->init_ScaleFactor(idIsoInput);
//         trigger->init_ScaleFactor(triggerInput);
//     }
//
//     template<typename LorentzVector>
//     double GetIdIsoSF(const LorentzVector& p4) const
//     {
//         return idIso->get_ScaleFactor(p4.pt(), p4.eta());
//     }
//
//     template<typename LorentzVector>
//     double GetTriggerSF(const LorentzVector& p4) const
//     {
//         return trigger->get_ScaleFactor(p4.pt(), p4.eta());
//     }
//
//     template<typename LorentzVector>
//     double GetTotalSF(const LorentzVector& p4) const { return GetIsoSF(p4) * GetTriggerSF(p4); }
//
// private:
//     std::shared_ptr<SF> idIso, trigger;
// };
//
// class MuonScaleFactorPOG {
// public:
//     using Hist = TH1;
//     using HistPtr = std::shared_ptr<Hist>;
//
//     MuonScaleFactorPOG(const std::string& idInput_B_F, const std::string& isoInput_B_F,
//                        const std::string& triggerInput_B_F, const std::string& idInput_G_H,
//                        const std::string& isoInput_G_H, const std::string& triggerInput_G_H) :
//         id_hist_B_F(LoadWeight(idInput_B_F,"MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")),
//         id_hist_G_H(LoadWeight(idInput_G_H,"MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")),
//         iso_hist_B_F(LoadWeight(isoInput_B_F,"LooseISO_TightID_pt_eta/pt_abseta_ratio")),
//         iso_hist_G_H(LoadWeight(isoInput_G_H,"LooseISO_TightID_pt_eta/pt_abseta_ratio")),
//         trigger_hist_B_F(LoadWeight(triggerInput_B_F,"IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio")),
//         trigger_hist_G_H(LoadWeight(triggerInput_G_H,"IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio")),
//         lumi_B_F(19.72), lumi_G_H(15.931)
//     {
// //        file_idInput_B_F = new TFile(idInput_B_F, "read");
//     }
//
//     template<typename LorentzVector>
//     double GetIdSF(const LorentzVector& p4) const
//     {
// //        TH2F* hist_B_F = (TH2F*)file_idInput_B_F->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
//         const Int_t bin_pt_B_F = id_hist_B_F->GetXaxis()->FindBin(p4.pt());
//         const Int_t bin_eta_B_F = id_hist_B_F->GetYaxis()->FindBin(std::abs(p4.eta()));
//         double sf_B_F = id_hist_B_F->GetBinContent(bin_pt_B_F,bin_eta_B_F);
//
//         const Int_t bin_pt_G_H = id_hist_G_H->GetXaxis()->FindBin(p4.pt());
//         const Int_t bin_eta_G_H = id_hist_G_H->GetYaxis()->FindBin(std::abs(p4.eta()));
//         double sf_G_H = id_hist_G_H->GetBinContent(bin_pt_G_H,bin_eta_G_H);
//
//         return ((sf_B_F * lumi_B_F) + (sf_G_H * lumi_G_H))/(lumi_B_F + lumi_G_H);
//     }
//
//     template<typename LorentzVector>
//     double GetIsoSF(const LorentzVector& p4) const
//     {
//         const Int_t bin_pt_B_F = iso_hist_B_F->GetXaxis()->FindBin(p4.pt());
//         const Int_t bin_eta_B_F = iso_hist_B_F->GetYaxis()->FindBin(std::abs(p4.eta()));
//         double sf_B_F = iso_hist_B_F->GetBinContent(bin_pt_B_F,bin_eta_B_F);
//
//         const Int_t bin_pt_G_H = iso_hist_G_H->GetXaxis()->FindBin(p4.pt());
//         const Int_t bin_eta_G_H = iso_hist_G_H->GetYaxis()->FindBin(std::abs(p4.eta()));
//         double sf_G_H = iso_hist_G_H->GetBinContent(bin_pt_G_H,bin_eta_G_H);
//
//         return ((sf_B_F * lumi_B_F) + (sf_G_H * lumi_G_H))/(lumi_B_F + lumi_G_H);
//     }
//
//     template<typename LorentzVector>
//     double GetTriggerSF(const LorentzVector& p4) const
//     {
//         const Int_t bin_pt_B_F = trigger_hist_B_F->GetXaxis()->FindBin(p4.pt());
//         const Int_t bin_eta_B_F = trigger_hist_B_F->GetYaxis()->FindBin(std::abs(p4.eta()));
//         double sf_B_F = trigger_hist_B_F->GetBinContent(bin_pt_B_F,bin_eta_B_F);
//
//         const Int_t bin_pt_G_H = trigger_hist_G_H->GetXaxis()->FindBin(p4.pt());
//         const Int_t bin_eta_G_H = trigger_hist_G_H->GetYaxis()->FindBin(std::abs(p4.eta()));
//         double sf_G_H = trigger_hist_G_H->GetBinContent(bin_pt_G_H,bin_eta_G_H);
//
//         return ((sf_B_F * lumi_B_F) + (sf_G_H * lumi_G_H))/(lumi_B_F + lumi_G_H);
//     }
//
//     template<typename LorentzVector>
//     double GetTotalSF(const LorentzVector& p4) const { return GetIdSF(p4) * GetIsoSF(p4) * GetTriggerSF(p4); }
//
//     static HistPtr LoadWeight(const std::string& weight_file_name, const std::string& hist_name)
//     {
//         auto file = root_ext::OpenRootFile(weight_file_name);
//         return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
//     }
//
// private:
//     HistPtr id_hist_B_F, id_hist_G_H, iso_hist_B_F, iso_hist_G_H, trigger_hist_B_F, trigger_hist_G_H;
//     double lumi_B_F, lumi_G_H;
//
// };
//
// } // namespace detail
//
// class LeptonWeights : public IWeightProvider {
// public:
//     using Event = ntuple::Event;
//
//     LeptonWeights(const std::string& electron_idIsoInput, const std::string& electron_triggerInput,
//                   const std::string& muon_idIsoInput, const std::string& muon_triggerInput) :
//         electronSF(electron_idIsoInput, electron_triggerInput),
//         muonSF(muon_idIsoInput, muon_triggerInput)
//     {
//     }
//
//     double GetIdIsoWeight(const Event& event) const
//     {
//         try {
//             const Channel channel = static_cast<Channel>(event.channelId);
//             if(channel == Channel::ETau) return electronSF.GetIdIsoSF(event.p4_1);
//             if(channel == Channel::MuTau) return muonSF.GetIdIsoSF(event.p4_1);
//         } catch(std::runtime_error& e) {
//             const EventIdentifier id(event);
//             std::cerr << id << " ERROR: " << e.what() << std::endl;
//         }
//         return 1.0;
//     }
//
//     double GetTriggerWeight(const Event& event) const
//     {
//         try {
//             const Channel channel = static_cast<Channel>(event.channelId);
//             if(channel == Channel::ETau) return electronSF.GetTriggerSF(event.p4_1);
//             if(channel == Channel::MuTau) return muonSF.GetTriggerSF(event.p4_1);
//         } catch(std::runtime_error& e) {
//             const EventIdentifier id(event);
//             std::cerr << id << " ERROR: " << e.what() << std::endl;
//         }
//         return 1.0;
//     }
//
//     virtual double Get(const Event& event) const override { return GetIdIsoWeight(event) * GetTriggerWeight(event); }
//
//     virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
//     {
//         throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
//     }
//
// private:
//     detail::LeptonScaleFactors electronSF, muonSF;
// //    detail::MuonScaleFactorPOG muonSF;
// };
//
// } // namespace mc_corrections
// } // namespace analysis
