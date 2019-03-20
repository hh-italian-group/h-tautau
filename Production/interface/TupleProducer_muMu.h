/*! Definition of an event tuple producer for the mu-mu channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Production/interface/BaseTupleProducer.h"

class TupleProducer_muMu: public BaseTupleProducer {
public:
    using SelectionResultsBase = analysis::SelectionResultsBase;
    using SelectionResultsBasePtr = std::shared_ptr<SelectionResultsBase>;

public:
    TupleProducer_muMu(const edm::ParameterSet& iConfig) : BaseTupleProducer(iConfig, analysis::Channel::MuMu) {}

private:
    virtual void ProcessEvent(Cutter& cut) override;
    void FillEventTuple(const SelectionResultsBase& selection);

    std::vector<MuonCandidate> CollectSignalMuons();
    void SelectSignalMuon(const MuonCandidate& muon, Cutter& cut) const;

private:
    SelectionResultsBasePtr previous_selection;
};
