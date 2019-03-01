/*! Definition of an event tuple producer for the mu-mu channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Production/interface/BaseTupleProducer.h"

class TupleProducer_muMu: public BaseTupleProducer {
public:
    using SelectionResults = analysis::SelectionResults<MuonCandidate, MuonCandidate>;
    using SelectionResultsPtr = std::shared_ptr<SelectionResults>;
    using HiggsCandidate = SelectionResults::HiggsCandidate;

public:
    TupleProducer_muMu(const edm::ParameterSet& iConfig) : BaseTupleProducer(iConfig, analysis::Channel::MuMu) {}

private:
    virtual void ProcessEvent(Cutter& cut) override;
    void FillEventTuple(const SelectionResults& selection);

    std::vector<MuonCandidate> CollectSignalMuons();
    void SelectSignalMuon(const MuonCandidate& muon, Cutter& cut) const;
    void FillHiggsDaughtersIndexes(const SelectionResults& selection);

private:
    SelectionResultsPtr previous_selection;
};
