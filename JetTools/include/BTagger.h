/*! b-jet tagging.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"

namespace analysis {

enum class JetOrdering { NoOrdering, Pt, CSV, DeepCSV, DeepFlavour };
ENUM_NAMES(JetOrdering) = {
    { JetOrdering::NoOrdering, "NoOrdering" },
    { JetOrdering::Pt, "Pt" },
    { JetOrdering::CSV, "CSV" },
    { JetOrdering::DeepCSV, "DeepCSV" },
    { JetOrdering::DeepFlavour, "DeepFlavour" },
};

struct BTagger {
public:
    BTagger(Period _period, JetOrdering _ordering);

    double BTag(const ntuple::Event& event, size_t jet_index) const;
    double BTag(const ntuple::TupleJet& jet) const;
    bool Pass(const ntuple::Event& event, size_t jet_index, DiscriminatorWP wp = DiscriminatorWP::Medium) const;
    bool Pass(const ntuple::TupleJet& jet, DiscriminatorWP wp = DiscriminatorWP::Medium) const;

    double PtCut() const;
    double EtaCut() const;

private:
    Period period;
    JetOrdering ordering;
    const std::map<DiscriminatorWP,double>* cut;
};

}
