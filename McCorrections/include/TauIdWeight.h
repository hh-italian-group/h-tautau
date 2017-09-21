/*! Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <fstream>
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {


class TauIdWeight : public IWeightProvider {
public:
    using Event = ntuple::Event;

    TauIdWeight(const std::string& tauId_input, const std::string& _iso_type) : iso_type(_iso_type)
    {
        using namespace boost::property_tree;

        ptree pt;
        read_json(tauId_input, pt);

        for (auto& type : pt){
            std::cout << type.first << "\n";
            std::vector<double>& param_value_map = tauIdparam_map[type.first];
            for (auto& prop : type.second){
                if(param_value_map.count(prop.first))
                    throw exception("Param already found");
                param_value_map[prop.first] = prop.second;
                std::cout << prop.first << ": " << prop.second.get_value<std::string>() << "\n";
            }
        }

    }

    double GetWeight(const Event& event) const
    {
        const Channel channel = static_cast<Channel>(event.channelId);
        if(channel == Channel::TauTau) {
            std::string type_name_data = "data_" + iso_type + "Iso_dm" + event.tauId_values_1;
            std::string type_name_mc = "mc_" + iso_type + "Iso_dm" + event.tauId_values_1;
            auto iter_data = tauIdparam_map.find(type_name_data);
            if(iter != tauIdparam_map.end())
                return iter->second;
            throw exception("ttbar merge weight not found for genEventType = %1%.") % event.genEventType;
            return electronSF.GetIdIsoSF(event.p4_1);
        }
        return 1.0;
    }

    virtual double Get(const Event& event) const override { return GetWeight(event); }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
    }

private:
    std::map<std::string, std::vector<double>> tauIdparam_map;
    std::string iso_type;

};

} // namespace mc_corrections
} // namespace analysis
