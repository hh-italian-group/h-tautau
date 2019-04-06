/*! Analyzer to check genMatch information.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"

struct Arguments {
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(std::string, output_file);
    REQ_ARG(analysis::Period, period);
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
        anaData(output) { }

    void Run()
    {
        static const std::map<Channel, std::pair<Float_t,Float_t>> genMatchInfo = {
            { Channel::ETau, { 3,5 } },
            { Channel::MuTau, { 4,5 } },
            { Channel::TauTau, { 5,5 } }
        };

        std::cout << boost::format("Processing input file '%1%' into output file '%2%'.\n")
                   % args.input_file() % args.output_file();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        auto originalTuple = ntuple::CreateEventTuple(args.tree_name(),originalFile.get(),true,ntuple::TreeState::Full);

        auto summaryTuple = ntuple::CreateSummaryTuple("summary", originalFile.get(), true, ntuple::TreeState::Full);
        summaryTuple->GetEntry(0);
        std::shared_ptr<SummaryInfo> summaryInfo(new SummaryInfo(summaryTuple->data()));
        const Channel channel = Parse<Channel>(args.tree_name());
        const Long64_t n_entries = originalTuple->GetEntries();
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
            originalTuple->GetEntry(current_entry);

            JetOrdering jet_ordering = args.period() == Period::Run2017 ? JetOrdering::DeepCSV : JetOrdering::CSV;

            boost::optional<analysis::EventInfoBase> event = CreateEventInfo(originalTuple->data(),TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,args.period(), jet_ordering, summaryInfo.get());
            if(!event.is_initialized()) continue;
            if(event->GetEnergyScale() != EventEnergyScale::Central) continue;
            if(!event->GetTriggerResults().AnyAcceptAndMatch()) continue;
            if((*event)->extraelec_veto || (*event)->extramuon_veto) continue;

            if(!genMatchInfo.count(channel))
                throw exception("Gen Match numbers not found for this channel.");
            std::pair<Float_t,Float_t> genMatchNumbers = genMatchInfo.at(channel);
            if(event->GetLeg(1)->gen_match() == static_cast<GenLeptonMatch>(genMatchNumbers.first) && event->GetLeg(2)->gen_match() == static_cast<GenLeptonMatch>(genMatchNumbers.second))
                anaData.mass_bothGood().Fill((event->GetLeg(1)->p4() + event->GetLeg(2)->p4()).mass());
            if(event->GetLeg(1)->gen_match() != static_cast<GenLeptonMatch>(genMatchNumbers.first) && event->GetLeg(2)->gen_match() != static_cast<GenLeptonMatch>(genMatchNumbers.second))
                anaData.mass_bothBad().Fill((event->GetLeg(1)->p4() + event->GetLeg(2)->p4()).mass());
            if(event->GetLeg(1)->gen_match() == static_cast<GenLeptonMatch>(genMatchNumbers.first) && event->GetLeg(2)->gen_match() != static_cast<GenLeptonMatch>(genMatchNumbers.second))
                anaData.mass_1stGood_2ndBad().Fill((event->GetLeg(1)->p4() + event->GetLeg(2)->p4()).mass());
            if(event->GetLeg(1)->gen_match() != static_cast<GenLeptonMatch>(genMatchNumbers.first) && event->GetLeg(2)->gen_match() == static_cast<GenLeptonMatch>(genMatchNumbers.second))
                anaData.mass_1stBad_2ndGood().Fill((event->GetLeg(1)->p4() + event->GetLeg(2)->p4()).mass());
            if(event->GetLeg(1)->gen_match() != static_cast<GenLeptonMatch>(genMatchNumbers.first) || event->GetLeg(2)->gen_match() != static_cast<GenLeptonMatch>(genMatchNumbers.second))
                anaData.mass_allBad().Fill((event->GetLeg(1)->p4() + event->GetLeg(2)->p4()).mass());

        }

    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    GenMatchCheckData anaData;
};

} // namespace analysis

PROGRAM_MAIN(analysis::GenMatchCheck, Arguments)
