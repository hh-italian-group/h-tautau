/*! b-jet tagging.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"

namespace analysis {

enum class BTaggerKind { NoTagger, Pt, CSV, DeepCSV, DeepFlavour, HHbtag };
ENUM_NAMES(BTaggerKind) = {
    { BTaggerKind::NoTagger, "NoTagger" },
    { BTaggerKind::Pt, "Pt" },
    { BTaggerKind::CSV, "CSV" },
    { BTaggerKind::DeepCSV, "DeepCSV" },
    { BTaggerKind::DeepFlavour, "DeepFlavour" },
    { BTaggerKind::HHbtag, "HHbtag" },
};

struct BTagger {
public:
    BTagger(Period _period, BTaggerKind _tagger);

    double BTag(const ntuple::TupleJet& jet, bool use_base_tagger) const;
    bool Pass(const ntuple::TupleJet& jet, DiscriminatorWP wp) const;

    double PtCut() const;
    double EtaCut() const;

    Period GetPeriod() const;
    BTaggerKind GetTagger() const;
    BTaggerKind GetBaseTagger() const;

private:
    Period period;
    BTaggerKind tagger, base_tagger;
    const std::map<DiscriminatorWP, double>* cut;
};

}
