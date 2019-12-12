/*! Definition of class EventCandidate.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Analysis/include/EventCandidate.h"

namespace analysis {

    EventCandidate::EventCandidate(const ntuple::Event& _event, UncertaintySource _uncertainty_source,
    UncertaintyScale _scale, analysis::Period _period, TauIdDiscriminator _tau_id_discriminator,
    TauIdDiscriminator _e_id_discriminator, DiscriminatorWP _tauVSeWP) : event(&_event),
    uncertainty_source(_uncertainty_source), scale(_scale), period(_period), tau_id_discriminator(_tau_id_discriminator),
    ele_id_discriminator(_e_id_discriminator), tauVSeWP(_tauVSeWP)  {}

    void EventCandidate::InitializeJecUncertainty(Period period, const std::string& working_path)
    {
        std::map<analysis::Period,std::string> file_uncertainty_sources = {
            {analysis::Period::Run2016,"h-tautau/McCorrections/data/2016/JES/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt"},
            {analysis::Period::Run2017,"h-tautau/McCorrections/data/2017/JES/Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt"},
            {analysis::Period::Run2018,"h-tautau/McCorrections/data/2018/JES/Autumn18_V8_MC_UncertaintySources_AK4PFchs.txt"}
        };
        if(!file_uncertainty_sources.count(period))
            throw exception("Period not found in file uncertainty source.");
        std::string full_path_source = tools::FullPath({working_path,file_uncertainty_sources.at(period)});
        jecUncertainties = std::make_shared<jec::JECUncertaintiesWrapper>(full_path_source);
    }

    const LepCollection& EventCandidate::GetLeptons()
    {
      if(!lepton_candidates) {
        CreateLeptons();
      }
      return *lepton_candidates;
    }

    const JetCollection& EventCandidate::GetJets()
    {
      if(!jet_candidates) {
          CreateJets();
      }
      return *jet_candidates;
    }

    const FatJetCollection& EventCandidate::GetFatJets()
    {
        if(!fatJets) {
            fatJets = std::shared_ptr<FatJetCollection>(new FatJetCollection());
            tuple_fatJets = std::make_shared<std::vector<ntuple::TupleFatJet>>();
            for(size_t n = 0; n < event->fatJets_p4.size(); ++n){
              tuple_fatJets->emplace_back(*event, n);
            }
            for(size_t n = 0; n < tuple_fatJets->size(); ++n){
              fatJets->emplace_back(tuple_fatJets->at(n));
            }
        }
        return *fatJets;
    }

    const MET& EventCandidate::GetMET()
    {
        if(!met) {
            CreateLeptons();
            if(jecUncertainties->JetUncertainties_withTotal().count(uncertainty_source))
                CreateJets();
        }
        return *met;
    }

    const ntuple::Event& EventCandidate::GetEvent() const
    {
        return *event;
    }

    UncertaintyScale EventCandidate::GetScale() const
    {
        return scale;
    }

    UncertaintySource EventCandidate::GetUncSource() const
    {
        return uncertainty_source;
    }

    void EventCandidate::CreateLeptons()
    {
        double shifted_met_px = 0;
        double shifted_met_py = 0;

        if(!met){
          tuple_met = std::make_shared<ntuple::TupleMet>(*event, MetType::PF);
          met = std::shared_ptr<MET>(new MET(*tuple_met, tuple_met->cov()));
        }


        lepton_candidates = std::make_shared<LepCollection>();
        tuple_leptons = std::make_shared<std::vector<ntuple::TupleLepton>>();
        for(size_t n = 0; n < event->lep_p4.size(); ++n){
          tuple_leptons->emplace_back(*event, n);
        }
        for(size_t n = 0; n < tuple_leptons->size(); ++n){
          lepton_candidates->emplace_back(tuple_leptons->at(n),tuple_leptons->at(n).iso());
        }
        for(size_t n = 0; n < tuple_leptons->size(); ++n) {
          auto tuple_lepton = tuple_leptons->at(n);
          LorentzVectorM lepton_p4(tuple_lepton.p4());
          LorentzVectorM corrected_lepton_p4(tuple_lepton.p4());
          if(tuple_lepton.leg_type() == analysis::LegType::tau){
              UncertaintyScale current_scale = uncertainty_source == UncertaintySource::TauES ? scale : UncertaintyScale::Central;
              double sf = TauESUncertainties::GetCorrectionFactor(period, tuple_lepton.decayMode(), tuple_lepton.gen_match(),
                                                                  uncertainty_source, current_scale,
                                                                  tuple_lepton.p4().pt(), tau_id_discriminator,
                                                                  ele_id_discriminator, tuple_lepton.p4().eta(),
                                                                  tauVSeWP);

              // static double GetCorrectionFactor(analysis::Period period, int decayMode, GenLeptonMatch genLeptonMatch,
              //                                   UncertaintySource unc_source, UncertaintyScale scale, double pt,
              //                                   TauIdDiscriminator tauVSjetDiscriminator, TauIdDiscriminator tauVSeDiscriminator,
              //                                   double eta, DiscriminatorWP tauVSeWP);

              // if(tuple_lepton.decayMode() == 0){
              //     double shifted_pt = lepton_p4.pt() * sf;
              //     corrected_lepton_p4 = LorentzVectorM(shifted_pt, lepton_p4.eta(), lepton_p4.phi(),lepton_p4.M());
              // }
              // else{
              corrected_lepton_p4 = lepton_p4 * sf;
              // }

              shifted_met_px += tuple_lepton.p4().px() - corrected_lepton_p4.px();
              shifted_met_py += tuple_lepton.p4().py() - corrected_lepton_p4.py();

          }
          lepton_candidates->at(n).SetMomentum(corrected_lepton_p4);
        }


        shifted_met_px += met->GetMomentum().px();
        shifted_met_py += met->GetMomentum().py();
        analysis::LorentzVectorXYZ shifted_met;
        double E = std::hypot(shifted_met_px,shifted_met_py);
        shifted_met.SetPxPyPzE(shifted_met_px,shifted_met_py,0,E);
        met->SetMomentum(shifted_met);


    }

    void EventCandidate::CreateJets()
    {
        if(jecUncertainties->JetUncertainties_withTotal().count(uncertainty_source)){
            tuple_met = std::make_shared<ntuple::TupleMet>(*event, MetType::PF);
            met = std::shared_ptr<MET>(new MET(*tuple_met, tuple_met->cov()));
        }

        jet_candidates = std::make_shared<JetCollection>();
        tuple_jets = std::make_shared<std::vector<ntuple::TupleJet>>();
        for(size_t n = 0; n < event->jets_p4.size(); ++n) {
          tuple_jets->emplace_back(*event, n);
        }
        for(size_t n = 0; n < tuple_jets->size(); ++n) {
          jet_candidates->emplace_back(tuple_jets->at(n));
        }

        if(jecUncertainties->JetUncertainties_withTotal().count(uncertainty_source)){
          const auto& other_jets_p4 = event->other_jets_p4;
          auto shifted_met_p4(met->GetMomentum());
          *jet_candidates = jecUncertainties->ApplyShift(*jet_candidates,uncertainty_source,scale,&other_jets_p4,&shifted_met_p4);
          met->SetMomentum(shifted_met_p4);
        }

    }


    std::shared_ptr<jec::JECUncertaintiesWrapper> EventCandidate::jecUncertainties;

} // namespace analysis
