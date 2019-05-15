/*! TauId weight test.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "h-tautau/McCorrections/include/TauIdWeight.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"
#include <iostream>


struct Arguments {
    run::Argument<std::string> input_file{"input_file", "input file"};
    run::Argument<std::string> json_file{"json_file", "json file"};
    run::Argument<std::string> iso_type{"iso_type", "iso type"};
    run::Argument<analysis::Period> period{"period","period",analysis::Period::Run2017};
    run::Argument<analysis::SignalMode> mode{"mode", "Signal mode"};
};

namespace analysis {

class TauIdWeight_t {
public:
    using TauIdWeight = analysis::mc_corrections::TauIdWeight;

    TauIdWeight_t(const Arguments& _args) : args(_args), signalObjectSelector(args.mode())
    {
        if(args.period()==Period::Run2017) tauId_weight=std::make_shared<mc_corrections::TauIdWeight2017>();
        else if(args.period()==Period::Run2016) tauId_weight=std::make_shared<mc_corrections::TauIdWeight2016>();
    }

public:

    void Run()
    {
        auto inputFile = root_ext::OpenRootFile(args.input_file());
        auto eventTuple = ntuple::CreateEventTuple("tauTau",inputFile.get(),true,ntuple::TreeState::Full);

        for(const ntuple::Event& event : *eventTuple) {
            boost::optional<analysis::EventInfoBase> eventInfo = CreateEventInfo(event,signalObjectSelector);
            if(!eventInfo.is_initialized()) continue;
            GenLeptonMatch gen_match = eventInfo->GetLeg(2)->gen_match();
            std::cout << "TauID weight: " << tauId_weight->GetIdIsoSF(eventInfo->GetLeg(2)->p4(), gen_match, eventInfo->GetLeg(2)->decayMode(), DiscriminatorWP::Loose,DiscriminatorWP::Loose,DiscriminatorWP::Medium) << std::endl;
        }

    }


private:
    Arguments args;
    std::shared_ptr<TauIdWeight> tauId_weight;
    SignalObjectSelector signalObjectSelector;

};

} //namespace analysis

PROGRAM_MAIN(analysis::TauIdWeight_t, Arguments)
