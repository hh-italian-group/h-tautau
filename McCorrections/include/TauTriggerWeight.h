/*! Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "WeightProvider.h"
#include "TauTriggerSFs2017.h"

namespace analysis {
namespace mc_corrections {
//namespace detail{

class TauTriggerWeight {
public:
    virtual double GetEfficiency(Channel channel, const LorentzVectorM& p4, GenLeptonMatch gen_match, int decay_mode,
                                 DiscriminatorWP iso_wp, bool isData) const = 0;

    virtual ~TauTriggerWeight() {}
};

class TauTriggerWeight2016  : public TauTriggerWeight {
public:
    using Event = ntuple::Event;
    using ptree = boost::property_tree::ptree;

    struct Key {
        bool isData{false};
        bool isGenuineTau{false};
        DiscriminatorWP isowp{DiscriminatorWP::Medium};
        int decayMode{0};

        static Key Parse(const std::string& key_string);
        bool operator < (const Key& other) const;
    };

    struct Parameters{
        double alpha, m_0, sigma, norm, n;

        double EvaluateEfficiency(double pt) const;
        static Parameters Parse(const ptree& entry);
    };

    TauTriggerWeight2016(const std::string& tauId_input);
    virtual double GetEfficiency(Channel channel, const LorentzVectorM& p4, GenLeptonMatch gen_match, int decay_mode,
                                 DiscriminatorWP iso_wp, bool isData) const override;

private:
    double EvaluateSF(double pt, GenLeptonMatch gen_match, int decay_mode, DiscriminatorWP iso_wp ) const;
    double EvaluateEfficiency(double pt, bool is_data, bool is_genuine, int decay_mode, DiscriminatorWP iso_wp) const;

    std::map<Key, Parameters> tauIdparam_map;
};

class TauTriggerWeight2017 : public TauTriggerWeight {
public:
    TauTriggerWeight2017(const std::string& tauTriggerInput, const std::string& tauTriggerInputOld,
                         std::string _tau_iso_wp);
    virtual double GetEfficiency(Channel channel, const LorentzVectorM& p4, GenLeptonMatch /*gen_match*/, int /*decay_mode*/,
                                 DiscriminatorWP /*iso_wp*/, bool isData) const override;

private:
    std::shared_ptr<TauTriggerSFs2017> tauSF;
};

} // namespace mc_corrections
} // namespace analysis
