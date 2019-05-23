/*! Definition of an event tuple producer for the e-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <memory>
#include <vector>
#include "TTree.h"
#include "Math/VectorUtil.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "h-tautau/Production/interface/BaseTupleProducer.h"

class TupleProducer_eTau: public BaseTupleProducer {
public:
    using SelectionResultsBase = analysis::SelectionResultsBase;
    using SelectionResultsBasePtr = std::shared_ptr<SelectionResultsBase>;

public:
    TupleProducer_eTau(const edm::ParameterSet& iConfig) : BaseTupleProducer(iConfig, analysis::Channel::ETau) {}

private:
    virtual void ProcessEvent(Cutter& cut) override;
    void FillEventTuple(const SelectionResultsBase& selection);

    std::vector<ElectronCandidate> CollectSignalElectrons();
    std::vector<TauCandidate> CollectSignalTaus();

    void SelectSignalElectron(const ElectronCandidate& electron, Cutter& cut) const;
    void SelectSignalTau(const TauCandidate& tau, Cutter& cut) const;

};
