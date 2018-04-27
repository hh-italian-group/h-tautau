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

    static const StorageMode& Full() { static const StorageMode m(0); return m; }

    explicit StorageMode(unsigned _mode = 0) : mode_bits(_mode) {}
    unsigned long Mode() const { return mode_bits.to_ulong(); }
    bool IsFull() const { return Mode() == 0; }
    bool IsMissing(EventPart part) const { return mode_bits[static_cast<size_t>(part)]; }
    bool IsPresent(EventPart part) const { return !IsMissing(part); }
    void SetPresence(EventPart part, bool presence) { mode_bits[static_cast<size_t>(part)] = !presence; }

    bool operator==(const StorageMode& other) const { return Mode() == other.Mode(); }
    bool operator!=(const StorageMode& other) const { return Mode() != other.Mode(); }

private:
    std::bitset<NumberOfParts> mode_bits;
};

#define CP_BR(br) event.br = ref->br
#define RAW_ID(name, n) CP_BR(tauId_##name##_##n);

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
            throw analysis::exception("Incompatible reference event number.");
        StorageMode ref_mode(ref->storageMode);
        if(!ref_mode.IsFull())
            throw analysis::exception("Incomplete reference event. Ref event storage mode = %1%") % ref_mode.Mode();

        if(mode.IsMissing(EventPart::FirstTauIds)) {
            RAW_TAU_IDS(1)
        }

        if(mode.IsMissing(EventPart::SecondTauIds)) {
            RAW_TAU_IDS(2)
        }

        if(mode.IsMissing(EventPart::Jets)) {
            CP_BR(jets_p4);
            CP_BR(jets_csv);
            CP_BR(jets_rawf);
            CP_BR(jets_mva);
            CP_BR(jets_partonFlavour);
            CP_BR(jets_hadronFlavour);
        }

        if(mode.IsMissing(EventPart::FatJets)) {
            CP_BR(fatJets_p4);
            CP_BR(fatJets_csv);
            CP_BR(fatJets_m_pruned);
            CP_BR(fatJets_m_softDrop);
            CP_BR(fatJets_n_subjettiness_tau1);
            CP_BR(fatJets_n_subjettiness_tau2);
            CP_BR(fatJets_n_subjettiness_tau3);

            CP_BR(subJets_p4);
            CP_BR(subJets_csv);
            CP_BR(subJets_parentIndex);
        }

        if(mode.IsMissing(EventPart::GenInfo)) {
            CP_BR(genEventType);
            CP_BR(genEventWeight);
            CP_BR(lhe_n_partons);
            CP_BR(lhe_n_c_partons);
            CP_BR(lhe_n_b_partons);
            CP_BR(lhe_HT);
            CP_BR(lhe_H_m);
            CP_BR(lhe_hh_m);
            CP_BR(lhe_hh_cosTheta);
            CP_BR(genParticles_pdg);
            CP_BR(genParticles_p4);
            CP_BR(genParticles_nPromptElectrons);
            CP_BR(genParticles_nPromptMuons);
            CP_BR(genParticles_nPromptTaus);
            CP_BR(genJets_nTotal);
            CP_BR(jets_nTotal_partonFlavour_b);
            CP_BR(jets_nTotal_partonFlavour_c);
            CP_BR(jets_nTotal_hadronFlavour_b);
            CP_BR(jets_nTotal_hadronFlavour_c);
            CP_BR(genJets_p4);
            CP_BR(genJets_partonFlavour);
            CP_BR(genJets_hadronFlavour);
        }
        return mode;
    }
};

#undef CP_BR
#undef RAW_ID

} // namespace nutple
