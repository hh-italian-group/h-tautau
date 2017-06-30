/*! Definition of loader of ntuple::Event.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
#include <bitset>
#include "EventTuple.h"

namespace ntuple {

class StorageMode {
public:
    static constexpr size_t NumberOfParts = 5;
    enum class EventPart { FirstTauIds = 0, SecondTauIds = 1, Jets = 2, FatJets = 3, GenInfo = 4 };

    explicit StorageMode(unsigned _mode = 0) : mode_bits(_mode) {}
    unsigned long Mode() const { return mode_bits.to_ulong(); }
    bool IsFull() const { return Mode() == 0; }
    bool IsMissing(EventPart part) const { return mode_bits[static_cast<size_t>(part)]; }
    bool IsPresent(EventPart part) const { return !IsMissing(part); }
    void SetPresence(EventPart part, bool presence) { mode_bits[static_cast<size_t>(part)] = !presence; }

private:
    std::bitset<NumberOfParts> mode_bits;
};

class EventLoader {
public:
    static StorageMode Load(Event& event, const Event* ref)
    {
        using EventPart = StorageMode::EventPart;

        StorageMode mode(event.storageMode);
        if(mode.IsFull()) return mode;
        if(!ref)
            throw analysis::exception("Can't load partially stored event without the reference.");
        if(event.run != ref->run || event.lumi != ref->lumi || event.evt != ref->evt)
            throw analysis::exception("Incompatible reverence event number.");
        StorageMode ref_mode(ref->storageMode);
        if(!ref_mode.IsFull())
            throw analysis::exception("Incomplete reference event. Ref event storage mode = %1%") % ref_mode.Mode();

        if(mode.IsMissing(EventPart::FirstTauIds)) {
            event.tauId_keys_1 = ref->tauId_keys_1;
            event.tauId_values_1 = ref->tauId_values_1;
        }

        if(mode.IsMissing(EventPart::SecondTauIds)) {
            event.tauId_keys_2 = ref->tauId_keys_2;
            event.tauId_values_2 = ref->tauId_values_2;
        }

        if(mode.IsMissing(EventPart::Jets)) {
            event.jets_p4 = ref->jets_p4;
            event.jets_csv = ref->jets_csv;
            event.jets_hadronFlavour = ref->jets_hadronFlavour;
            event.jets_mva = ref->jets_mva;
            event.jets_rawf = ref->jets_rawf;
        }

        if(mode.IsMissing(EventPart::FatJets)) {
            event.fatJets_p4 = ref->fatJets_p4;
            event.fatJets_csv = ref->fatJets_csv;
            event.fatJets_m_pruned = ref->fatJets_m_pruned;
            event.fatJets_m_softDrop = ref->fatJets_m_softDrop;
            event.fatJets_n_subjettiness_tau1 = ref->fatJets_n_subjettiness_tau1;
            event.fatJets_n_subjettiness_tau2 = ref->fatJets_n_subjettiness_tau2;
            event.fatJets_n_subjettiness_tau3 = ref->fatJets_n_subjettiness_tau3;

            event.subJets_p4 = ref->subJets_p4;
            event.subJets_csv = ref->subJets_csv;
            event.subJets_parentIndex = ref->subJets_parentIndex;
        }

        if(mode.IsMissing(EventPart::GenInfo)) {
            event.lhe_n_partons = ref->lhe_n_partons;
            event.lhe_n_b_partons = ref->lhe_n_b_partons;
            event.lhe_HT = ref->lhe_HT;
            event.genEventType = ref->genEventType;
            event.genEventWeight = ref->genEventWeight;
            event.genParticles_pdg = ref->genParticles_pdg;
            event.genParticles_p4 = ref->genParticles_p4;
            event.genJets_nTotal = ref->genJets_nTotal;
            event.genJets_p4 = ref->genJets_p4;
        }
        return mode;
    }
};

} // namespace nutple
