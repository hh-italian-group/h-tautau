/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/TauUncertainties.h"
#include "h-tautau/JetTools/include/JECUncertaintiesWrapper.h"
#include "AnalysisTools/Core/include/Tools.h"

namespace analysis {
    using LepCandidate = LeptonCandidate<ntuple::TupleLepton>;
    using LepCollection = std::vector<LepCandidate>;
    using JetCandidate = Candidate<ntuple::TupleJet>;
    using JetCollection = std::vector<JetCandidate>;
    using FatJetCandidate = Candidate<ntuple::TupleFatJet>;
    using FatJetCollection = std::vector<FatJetCandidate>;
    using MET = MissingET<ntuple::TupleMet>;

class EventCandidate {
public:

    EventCandidate(const ntuple::Event& _event, UncertaintySource _uncertainty_source,
                   UncertaintyScale _scale, Period _period, TauIdDiscriminator _tau_id_discriminator,
                   TauIdDiscriminator _ele_id_discriminator);

    EventCandidate(const EventCandidate& ) = default; //copy constructor
    EventCandidate(EventCandidate&& ) = default; // move constructor

    EventCandidate& operator= ( const EventCandidate& ) = default; //assignment

    static void InitializeJecUncertainty(Period period, const std::string& working_path);

    const LepCollection& GetLeptons();
    const JetCollection& GetJets();
    const FatJetCollection& GetFatJets();
    const MET& GetMET();
    const ntuple::Event& GetEvent() const;
    UncertaintyScale GetScale() const;
    UncertaintySource GetUncSource() const;

private:
    void CreateLeptons();
    void CreateJets();

    const ntuple::Event* event;
    UncertaintySource uncertainty_source;
    UncertaintyScale scale;
    analysis::Period period;
    std::shared_ptr<std::vector<ntuple::TupleLepton>> tuple_leptons;
    std::shared_ptr<std::vector<ntuple::TupleJet>> tuple_jets;
    std::shared_ptr<std::vector<ntuple::TupleFatJet>> tuple_fatJets;
    std::shared_ptr<FatJetCollection> fatJets;
    std::shared_ptr<ntuple::TupleMet> tuple_met;
    std::shared_ptr<std::vector<LepCandidate>> lepton_candidates;
    std::shared_ptr<std::vector<JetCandidate>> jet_candidates;
    std::shared_ptr<MET> met;
    static std::shared_ptr<jec::JECUncertaintiesWrapper> jecUncertainties;
    TauIdDiscriminator tau_id_discriminator;
    TauIdDiscriminator ele_id_discriminator;
    DiscriminatorWP tauVSeWP;
};

} // namespace analysis
