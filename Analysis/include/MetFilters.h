/*! Definition of interface to store MET filters results.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include<bitset>

namespace ntuple {

class MetFilters {
public:
    static constexpr size_t NumberOfFilters = 9;
    enum class Filter {
        PrimaryVertex = 0, BeamHalo = 1, HBHE_noise = 2, HBHEiso_noise = 3, ECAL_TP = 4, ee_badSC_noise = 5,
        badMuon = 6, badChargedHadron = 7, ecalBadCalib = 8
    };

    explicit MetFilters(unsigned filter_results = 0) : bits(filter_results) {}
    unsigned FilterResults() const { return static_cast<unsigned int>(bits.to_ulong()); }
    bool PassAll() const { return bits.all(); }
    bool Fail(Filter filter) const { return !bits[static_cast<size_t>(filter)]; }
    bool Pass(Filter filter) const { return !Fail(filter); }
    void SetResult(Filter filter, bool result) { bits[static_cast<size_t>(filter)] = result; }

private:
    std::bitset<NumberOfFilters> bits;
};

} // namespace nutple
