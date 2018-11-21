/*! Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "WeightProvider.h"
#include "TextIO.h"
#include "TString.h"
#include "TauTriggerSFs2017.h"

namespace analysis {
namespace mc_corrections {
//namespace detail{

class TauTriggerWeight {
public:
    virtual double GetEfficiency(Channel channel, const LorentzVectorM& p4, GenMatch gen_match, int decay_mode,
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

        static Key Parse(const std::string& key_string)
        {
            Key key;
            static const std::map<std::string, bool> data_map = { {"data", true}, {"mc" , false } };
            static const std::map<std::string, bool> genuine_map = { {"genuine", true}, {"fake" , false } };

            std::vector<std::string> key_elements = SplitValueList(key_string, true, "_");
            if(key_elements.size() != 4)
                throw exception("More than 4 key params found");

            if(!data_map.count(key_elements.at(0)))
                throw exception("Invalid tauId key %1%." ) % key_string;
            key.isData = data_map.at(key_elements.at(0));
            if(!genuine_map.count(key_elements.at(1)))
                throw exception("Invalid tauId key %1%." ) % key_string;
            key.isGenuineTau = genuine_map.at(key_elements.at(1));
            if(key_elements.at(2).substr(key_elements.at(2).size() - 3) != "Iso")
                throw exception("Invalid tauId key %1%." ) % key_string;
            if(!TryParse(key_elements.at(2).substr(0, key_elements.at(2).size() - 3), key.isowp))
                throw exception("Invalid tauId key %1%." ) % key_string;
            if(key_elements.at(3).substr(0,2) != "dm")
                throw exception("Invalid tauId key %1%." ) % key_string;
            if(!TryParse(key_elements.at(3).substr(2), key.decayMode))
                throw exception("Invalid tauId key %1%." ) % key_string;

            return key;
        }

        bool operator < (const Key& other) const
        {
            if(isData != other.isData) return isData < other.isData;
            if(isGenuineTau != other.isGenuineTau) return isGenuineTau < other.isGenuineTau;
            if(isowp != other.isowp) return isowp < other.isowp;
            return decayMode < other.decayMode;
        }
    };

    struct Parameters{
        double alpha, m_0, sigma, norm, n;

        double EvaluateEfficiency(double pt) const
        {
            return pt < 1000 ? crystalball(pt, m_0, sigma, alpha, n, norm) : 1;
        }

        static Parameters Parse(const ptree& entry)
        {
            Parameters parameters;

            parameters.alpha = entry.get_child("alpha").get_value<double>();
            parameters.m_0 = entry.get_child("m_{0}").get_value<double>();
            parameters.sigma = entry.get_child("sigma").get_value<double>();
            parameters.norm = entry.get_child("norm").get_value<double>();
            parameters.n = entry.get_child("n").get_value<double>();

            return parameters;
        }
    };

    TauTriggerWeight2016(const std::string& tauId_input)
    {
        ptree property_tree;
        boost::property_tree::json_parser::read_json(tauId_input, property_tree);
        for (auto& type_entry : property_tree) {
            const Key key = Key::Parse(type_entry.first);
            tauIdparam_map[key] = Parameters::Parse(type_entry.second);
        }
    }

     virtual double GetEfficiency(Channel channel, const LorentzVectorM& p4, GenMatch gen_match, int decay_mode,
         DiscriminatorWP iso_wp, bool isData) const override
    {
        const bool is_genuine = gen_match == GenMatch::Tau;
        if(channel == Channel::TauTau)
            return EvaluateEfficiency(p4.pt(), isData, is_genuine, decay_mode, iso_wp);
        else
            throw exception ("channel %1% not supported") % channel;
    }

    private:
    double EvaluateSF(double pt, GenMatch gen_match, int decay_mode, DiscriminatorWP iso_wp ) const
    {
        const bool is_genuine = gen_match == GenMatch::Tau;
        const double eff_data = EvaluateEfficiency(pt, true, is_genuine, decay_mode, iso_wp);
        const double eff_mc = EvaluateEfficiency(pt, false, is_genuine, decay_mode, iso_wp);
        return eff_data / eff_mc;
    }

    double EvaluateEfficiency(double pt, bool is_data, bool is_genuine, int decay_mode, DiscriminatorWP iso_wp) const
    {
        const Key key{is_data, is_genuine, iso_wp, decay_mode};
        auto iter = tauIdparam_map.find(key);
        if(iter == tauIdparam_map.end()) {
            const std::string name = is_data ? "data" : "mc";
            throw exception("Key %1% not found for TauId weight.") % name;
        }
        return iter->second.EvaluateEfficiency(pt);
    }

    std::map<Key, Parameters> tauIdparam_map;
};

class TauTriggerWeight2017 : public TauTriggerWeight {
public:
    TauTriggerWeight2017(const std::string& tauTriggerInput, const std::string& tauTriggerInputOld, std::string _tau_iso_wp)
        : tauSF(std::make_shared<TauTriggerSFs2017>(tauTriggerInput, tauTriggerInputOld, _tau_iso_wp))
    {}

    virtual double GetEfficiency(Channel channel, const LorentzVectorM& p4, GenMatch /*gen_match*/, int /*decay_mode*/,
                                 DiscriminatorWP /*iso_wp*/, bool isData) const override
    {
        if(channel == Channel::ETau){
            return isData ? tauSF->getETauEfficiencyData(p4.pt(), p4.eta(), p4.phi(), TauTriggerSFs2017::kCentral)
                          : tauSF->getETauEfficiencyMC(p4.pt(), p4.eta(), p4.phi(), TauTriggerSFs2017::kCentral);
        }

        else if(channel == Channel::MuTau){
            return isData ? tauSF->getMuTauEfficiencyData(p4.pt(), p4.eta(), p4.phi(), TauTriggerSFs2017::kCentral)
                          : tauSF->getMuTauEfficiencyMC(p4.pt(), p4.eta(), p4.phi(), TauTriggerSFs2017::kCentral);
        }

        else if(channel == Channel::TauTau){
            return isData ? tauSF->getDiTauEfficiencyData(p4.pt(), p4.eta( ), p4.phi(),TauTriggerSFs2017::kCentral)
                          : tauSF->getDiTauEfficiencyMC(p4.pt(), p4.eta(), p4.phi(), TauTriggerSFs2017::kCentral);
        }

        throw exception ("channel not allowed");
    }
private:
    std::shared_ptr<TauTriggerSFs2017> tauSF;
};
}
}
