#include <iostream>
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/GenParticle.h"
#include "AnalysisTools/Core/include/RootExt.h"

#include "AnalysisTools/Core/include/EventIdentifier.h"

struct Arguments {
    REQ_ARG(std::string, inputPath);
    REQ_ARG(analysis::Channel, channel);
    REQ_ARG(std::string, event_id);
};

namespace analysis {

class PrintGen {
public:
    PrintGen(const Arguments& _args): args(_args)
    {
        GenEvent::InitializeParticleDataTable("h-tautau/Analysis/config/pdg_name_type_charge.txt");
    }

    void Run()
    {
        auto file = root_ext::OpenRootFile(args.inputPath());
        auto tuple = ntuple::CreateEventTuple(ToString(args.channel()), file.get(), true, ntuple::TreeState::Full);

        for(const auto& event : *tuple) {
            GenEvent genEvent(event);
            const EventIdentifier EventId(event.run, event.lumi, event.evt);
            const EventIdentifier EventIdTest(args.event_id());
            if(EventId != EventIdTest) continue;
            genEvent.Print();

        }
    }
private:
    Arguments args;
};
} // namespace analysis
PROGRAM_MAIN(analysis::PrintGen, Arguments)
