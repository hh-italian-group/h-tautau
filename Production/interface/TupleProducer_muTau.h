/*! Definition of an event tuple producer for the mu-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <memory>
#include <vector>
#include "TTree.h"
#include "Math/VectorUtil.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "h-tautau/Production/interface/BaseTupleProducer.h"

class TupleProducer_muTau: public BaseTupleProducer {
public:
    using SelectionResults = analysis::SelectionResults<MuonCandidate>;
    using HiggsCandidate = SelectionResults::HiggsCandidate;

public:
    TupleProducer_muTau(const edm::ParameterSet& iConfig) : BaseTupleProducer(iConfig, "muTau") {}

private:
    virtual void ProcessEvent(Cutter& cut) override;

    std::vector<MuonCandidate> CollectSignalMuons();
    std::vector<TauCandidate> CollectSignalTaus();

    void SelectSignalMuon(const MuonCandidate& muon, Cutter& cut) const;
    void SelectSignalTau(const TauCandidate& tau, Cutter& cut) const;

    void FillEventTuple(const SelectionResults& selection);
    void FillSyncTuple(const SelectionResults& selection);

};
