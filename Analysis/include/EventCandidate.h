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
#include "h-tautau/Analysis/include/TauUncertainties.h"
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
    UncertaintyScale _scale, analysis::Period _period);

    static void InitializeJecUncertainty(const std::string& file_uncertainty_source)
    {
        jecUncertainties = std::make_shared<jec::JECUncertaintiesWrapper>(file_uncertainty_source);
    }

    LepCollection& GetLeptons();
    JetCollection& GetJets();
    MET& GetMET();
    const ntuple::Event& GetEvent() const;
    const UncertaintyScale& GetScale();

private:
    void CreateLeptons();
    void CreateJets();

  ntuple::Event event;
  UncertaintySource uncertainty_source;
  UncertaintyScale scale;
  analysis::Period period;
  std::shared_ptr<std::vector<ntuple::TupleLepton>> tuple_leptons;
  std::shared_ptr<std::vector<ntuple::TupleJet>> tuple_jets;
  std::shared_ptr<ntuple::TupleMet> tuple_met;
  std::shared_ptr<std::vector<LepCandidate>> lepton_candidates;
  std::shared_ptr<std::vector<JetCandidate>> jet_candidates;
  std::shared_ptr<MET> met;
  static std::shared_ptr<jec::JECUncertaintiesWrapper> jecUncertainties;
};

} // namespace analysis
