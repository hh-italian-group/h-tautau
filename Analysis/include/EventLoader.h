/*! Definition of loader of ntuple::Event.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <bitset>
#include "h-tautau/Core/include/EventTuple.h"

namespace ntuple {

class StorageMode {
public:
    static constexpr size_t NumberOfParts = 6;
    enum class EventPart { FirstTauIds = 0, SecondTauIds = 1, Jets = 2, FatJets = 3, GenInfo = 4, OtherLeptons = 5 };

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

class EventLoader {
public:
    static StorageMode Load(Event& event, const Event* ref);
};

} // namespace nutple
