/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/SummaryTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2017.h"
#include "h-tautau/Cuts/include/H_tautau_2016_baseline.h"
#include "h-tautau/Cuts/include/H_tautau_2017_baseline.h"
#include "h-tautau/Analysis/include/MetFilters.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/McCorrections/include/TauUncertainties.h"
#include "h-tautau/JetTools/include/JECUncertaintiesWrapper.h"

namespace analysis {

using analysis::UncertaintySource;
using analysis::UncertaintyScale;

class EventCandidate {
public:
  using LepCandidate = LeptonCandidate<ntuple::TupleLepton>;
  using LepCollection = std::vector<LepCandidate>;
  using JetCandidate = Candidate<ntuple::TupleJet>;
  using JetCollection = std::vector<JetCandidate>;
  using MET = MissingET<ntuple::TupleMet>;

  EventCandidate(const ntuple::Event& _event, UncertaintySource _uncertainty_source,
    UncertaintyScale _scale, analysis::Period _period) : event(_event), uncertainty_source(_uncertainty_source),
    scale(_scale), period(_period) {}

    static void InitializeJecUncertainty(const std::string& file_uncertainty_source)
    {
        if(!file_uncertainty_source.empty())
            jecUncertainties = std::make_shared<jec::JECUncertaintiesWrapper>(file_uncertainty_source);
    }

    const LepCollection& GetLeptons()
    {
      if(!lepton_candidates) {
        CreateLeptons();
      }
      return *lepton_candidates;
    }

  const JetCollection& GetJets(const std::string& file_uncertainty_source)
  {
      if(!jet_candidates) {
          CreateJets(file_uncertainty_source);
      }
      return *jet_candidates;
  }

  const MET& GetMET()
  {
      if(!met) {
          tuple_met = std::make_shared<ntuple::TupleMet>(event, MetType::PF);
          met = std::shared_ptr<MET>(new MET(*tuple_met, tuple_met->cov()));
      }
      return *met;
  }

  const ntuple::Event& GetEvent()
  {
      return *event;
  }

private:
  void CreateLeptons()
  {
      static const std::map<analysis::Period, std::map<int, double>> tau_correction_factor = {
        { analysis::Period::Run2016, { {0, analysis::uncertainties::tau_2016::sf_1prong},
                                       {1, analysis::uncertainties::tau_2016::sf_1prongPi0},
                                       {10, analysis::uncertainties::tau_2016::sf_3prong} }},
        { analysis::Period::Run2017, { {0, analysis::uncertainties::tau_2017::sf_1prong},
                                       {1, analysis::uncertainties::tau_2017::sf_1prongPi0},
                                       {10, analysis::uncertainties::tau_2017::sf_3prong} } },
       { analysis::Period::Run2018, { {0, analysis::uncertainties::tau_2017::sf_1prong},
                                      {1, analysis::uncertainties::tau_2017::sf_1prongPi0},
                                      {10, analysis::uncertainties::tau_2017::sf_3prong} } }
      };

      static const std::map<analysis::Period, double> tau_energyUncertainty = {
          { analysis::Period::Run2016, analysis::uncertainties::tau_2016::energyUncertainty},
          { analysis::Period::Run2017, analysis::uncertainties::tau_2017::energyUncertainty},
          { analysis::Period::Run2018, analysis::uncertainties::tau_2017::energyUncertainty}
      };
      double shifted_met_px = 0;
      double shifted_met_py = 0;
      bool met_shift_applied = false;
      tuple_met = std::make_shared<ntuple::TupleMet>(event, MetType::PF);
      met = std::shared_ptr<MET>(new MET(*tuple_met, tuple_met->cov()));

      lepton_candidates = std::make_shared<LepCollection>();
      tuple_leptons = std::make_shared<std::list<ntuple::TupleLepton>>();
      for(size_t n = 0; n < event.lep_p4.size(); ++n) {
          auto tuple_lepton = ntuple::TupleLepton(event, n); //forse meglio il candidate...
          if(tuple_lepton.leg_type() == analysis::LegType::tau && uncertainty_source == UncertaintySource::None){
              bool tau_es_set = false;
              double corr_factor = 1;
              double sf = 1;
              double tau_es_var = 1;
              bool tau_es_sf_set = false;
              if(tau_correction_factor.at(period).count(tuple_lepton.decayMode())){
                  corr_factor = tau_correction_factor.at(period).at(tau.decayMode());
                  sf = corr_factor;
                  tau_es_set = true;
                  if(scale != UncertaintyScale::Central) {
                      const int sign = static_cast<int>(scale);
                      tau_es_var = sign * tau_energyUncertainty.at(period);
                      sf = corr_factor + tau_es_var;
                      tau_es_sf_set = true;
                  }
              }

              if(tau_es_sf_set){
                  const auto shiftedMomentum_met = tuple_lepton.p4() * tau_es_var;
                  shifted_met_px -= shiftedMomentum_met.px();
                  shifted_met_py -= shiftedMomentum_met.py();
                  met_shift_applied = true;
              }
              if(tau_es_set){//modfiy
                  auto shiftedMomentum; //LorentzVector??
                  if(tuple_lepton.decayMode() == 0){
                      double shifted_pt = tuple_lepton.p4().pt() * sf;
                      //how????
                  }
                  else{
                      shiftedMomentum = tau.p4() * sf;
                  }
                  tauCandidate.SetMomentum(shiftedMomentum); ///how???
              }
          }

          if(met_shift_applied){
                shifted_met_px += met->GetMomentum().px();
                shifted_met_py += met->GetMomentum().py();
                analysis::LorentzVectorXYZ shifted_met;
                double E = std::hypot(shifted_met_px,shifted_met_py);
                shifted_met.SetPxPyPzE(shifted_met_px,shifted_met_py,0,E);
                met->SetMomentum(shifted_met); //correct???
          }


          tuple_leptons->push_back(tuple_lepton);
          lepton_candidates->push_back(LepCandidate(tuple_leptons->back()));
       }

  }

  void CreateJets(const std::string& file_uncertainty_source)
  {
      jet_candidates = std::make_shared<JetCollection>();
      tuple_jets = std::make_shared<std::list<ntuple::TupleJet>>();
      //shift
      if(uncertainty_source != UncertaintySource::None){
          InitializeJecUncertainty(file_uncertainty_source);
          const JetCollection& jets = shifted_event_info->GetJets();
          const auto& other_jets_p4 = event->other_jets_p4;
          auto shifted_met_p4(shifted_event_info->GetMET().GetMomentum());
          const JetCollection& corrected_jets = jecUncertainties.ApplyShift(jets,uncertainty_source,scale,&other_jets_p4,&shifted_met_p4);
          shifted_event_info->SetJets(corrected_jets); //jets = std::make_shared<JetCollection>(new_jets);
          shifted_event_info->SetMetMomentum(shifted_met_p4);
          //met = std::make_shared<MET>(*tuple_met, tuple_met->cov());
          //met->SetMomentum(new_met_p4);
      }

      //fill come correrarli??
      for(size_t n = 0; n < event.jets_p4.size(); ++n) {
          auto tuple_jet = ntuple::TupleJet(event, n);


          tuple_jets->push_back(tuple_jet);
          jet_candidates->push_back(JetCandidate(tuple_jets->back()));
      }
  }

  ntuple::Event event;
  UncertaintySource uncertainty_source;
  UncertaintyScale scale;
  analysis::Period period;
  std::shared_ptr<std::list<ntuple::TupleLepton>> tuple_leptons;
  std::shared_ptr<std::list<ntuple::TupleJet>> tuple_jets;
  std::shared_ptr<ntuple::TupleMet> tuple_met;
  std::vector<std::shared_ptr<LepCandidate>> lepton_candidates;
  std::vector<std::shared_ptr<JetCandidate>> jet_candidates;
  std::shared_ptr<MET> met;
  std::shared_ptr<jec::JECUncertaintiesWrapper> jecUncertainties;
};

} // namespace analysis
