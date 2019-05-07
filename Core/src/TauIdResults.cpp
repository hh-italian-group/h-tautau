/*! Definition of Tau ID discriminators.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Core/include/TauIdResults.h"
#include "AnalysisTools/Core/include/TextIO.h"

namespace analysis {

namespace tau_id {

    TauIdDescriptor::TauIdDescriptor(TauIdDiscriminator _discriminator, const std::string& _name_pattern, bool _has_raw,
                    const std::string& wp_list) :
        discriminator(_discriminator), name_pattern(_name_pattern), has_raw(_has_raw)
    {
        if(has_raw)
            raw_name = ToStringRaw();
        auto wp_names = SplitValueList(wp_list, false, ", \t", true);
        for(const auto& wp_name : wp_names) {
            const DiscriminatorWP wp = ::analysis::Parse<DiscriminatorWP>(wp_name);
            working_points[wp] = ToString(wp);
        }
    }

    std::string TauIdDescriptor::ToString(DiscriminatorWP wp) const
    {
        std::string name = name_pattern;
        boost::algorithm::replace_all(name, "{wp}", analysis::ToString(wp));
        boost::algorithm::replace_all(name, "{raw}", "");
        boost::algorithm::replace_all(name, "{Raw}", "");
        return name;
    }


    std::string TauIdDescriptor::ToStringRaw() const
    {
        std::string name = name_pattern;
        boost::algorithm::replace_all(name, "{wp}", "");
        boost::algorithm::replace_all(name, "{raw}", "raw");
        boost::algorithm::replace_all(name, "{Raw}", "Raw");
        return name;
    }


}

} // namespace analysis
