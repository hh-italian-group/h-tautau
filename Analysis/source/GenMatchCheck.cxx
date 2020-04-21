/*! Analyzer to check genMatch information.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"

struct Arguments {
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(std::string, output_file);
    REQ_ARG(analysis::Period, period);
    REQ_ARG(analysis::SignalMode, mode);
};

namespace analysis {

class GenMatchCheckData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;

    TH1D_ENTRY(mass_bothGood, 50, 0, 300)
    TH1D_ENTRY(mass_bothBad, 50, 0, 300)
    TH1D_ENTRY(mass_1stGood_2ndBad, 50, 0, 300)
    TH1D_ENTRY(mass_1stBad_2ndGood, 50, 0, 300)
    TH1D_ENTRY(mass_allBad, 50, 0, 300)

};

class GenMatchCheck {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    GenMatchCheck(const Arguments& _args) : args(_args), output(root_ext::CreateRootFile(args.output_file())),
        anaData(output), signalObjectSelector(args.mode()) { }

    void Run()
    {
        static const std::map<Channel, std::pair<GenLeptonMatch,GenLeptonMatch>> genMatchInfo = {
            { Channel::ETau, { GenLeptonMatch::TauElectron,GenLeptonMatch::Tau } },
            { Channel::MuTau, { GenLeptonMatch::TauMuon,GenLeptonMatch::Tau } },
            { Channel::TauTau, { GenLeptonMatch::Tau,GenLeptonMatch::Tau } }
        };

        std::cout << boost::format("Processing input file '%1%' into output file '%2%'.\n")
                   % args.input_file() % args.output_file();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        auto originalTuple = ntuple::CreateEventTuple(args.tree_name(),originalFile.get(),true,ntuple::TreeState::Full);

        auto summaryTuple = ntuple::CreateSummaryTuple("summary", originalFile.get(), true, ntuple::TreeState::Full);
        summaryTuple->GetEntry(0);
        auto summaryInfo = std::make_shared<SummaryInfo>(summaryTuple->data(), Parse<Channel>(args.tree_name()));
        const Channel channel = Parse<Channel>(args.tree_name());
        const Long64_t n_entries = originalTuple->GetEntries();
        const BTagger bTagger(args.period(), BTaggerKind::DeepFlavour);
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
            originalTuple->GetEntry(current_entry);

            auto event = EventInfo::Create(originalTuple->data(), signalObjectSelector, bTagger,
                                           DiscriminatorWP::Medium, summaryInfo);
            if(!event) continue;
            // if(event->GetEnergyScale() != EventEnergyScale::Central) continue;
            if(!event->GetTriggerResults().AnyAcceptAndMatch()) continue;
            if((*event)->extraelec_veto || (*event)->extramuon_veto) continue;

            if(!genMatchInfo.count(channel))
                throw exception("Gen Match numbers not found for this channel.");
            std::pair<GenLeptonMatch,GenLeptonMatch> genMatchNumbers = genMatchInfo.at(channel);
            if(event->GetLeg(1)->gen_match() == genMatchNumbers.first && event->GetLeg(2)->gen_match() == genMatchNumbers.second)
                anaData.mass_bothGood().Fill((event->GetLeg(1)->p4() + event->GetLeg(2)->p4()).mass());
            if(event->GetLeg(1)->gen_match() != genMatchNumbers.first && event->GetLeg(2)->gen_match() != genMatchNumbers.second)
                anaData.mass_bothBad().Fill((event->GetLeg(1)->p4() + event->GetLeg(2)->p4()).mass());
            if(event->GetLeg(1)->gen_match() == genMatchNumbers.first && event->GetLeg(2)->gen_match() != genMatchNumbers.second)
                anaData.mass_1stGood_2ndBad().Fill((event->GetLeg(1)->p4() + event->GetLeg(2)->p4()).mass());
            if(event->GetLeg(1)->gen_match() != genMatchNumbers.first && event->GetLeg(2)->gen_match() == genMatchNumbers.second)
                anaData.mass_1stBad_2ndGood().Fill((event->GetLeg(1)->p4() + event->GetLeg(2)->p4()).mass());
            if(event->GetLeg(1)->gen_match() != genMatchNumbers.first || event->GetLeg(2)->gen_match() != genMatchNumbers.second)
                anaData.mass_allBad().Fill((event->GetLeg(1)->p4() + event->GetLeg(2)->p4()).mass());

        }

    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    GenMatchCheckData anaData;
    SignalObjectSelector signalObjectSelector;
};

} // namespace analysis

PROGRAM_MAIN(analysis::GenMatchCheck, Arguments)
