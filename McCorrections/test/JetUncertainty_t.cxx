/*! Check SM weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/ConfigReader.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "McCorrections/include/JetCorrectionUncertainty.h"


struct Arguments {
    run::Argument<std::string> input_file{"input_file", "file to test"};
    run::Argument<std::string> unc_txt{"unc_txt", "uncertainty source txt"};
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> output_file{"output_file", "Output root file"};
};

namespace analysis {

class JetUncertaintyData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    TH1D_ENTRY(m_sv, 200, 0, 400)
    TH1D_ENTRY(m_bb, 200, 0, 400)
    TH1D_ENTRY(m_kinFit, 250, 200, 1200)
};


class JetUncertainty_t {
public:

    JetUncertainty_t(const Arguments& _args) :
        args(_args), output(root_ext::CreateRootFile(args.output_file())), anaData(output)
    {
        LoadInputs();
    }

    void Run() {}

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    JetUncertaintyData anaData;

    static const std::set<std::string>& GetEnabledBranches()
    {
        static const std::set<std::string> EnabledBranches_read = {
            "eventEnergyScale", "lhe_hh_m", "lhe_hh_cosTheta", "jets_p4", "SVfit_p4", "kinFit_m"
        };
        return EnabledBranches_read;
    }

    void LoadInputs()
    {
        auto inputFile = root_ext::OpenRootFile(args.input_file());
        ntuple::EventTuple eventTuple(args.tree_name(), inputFile.get(), true, {}, GetEnabledBranches());
        const Long64_t n_entries = eventTuple.GetEntries();

        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(args.unc_txt());

        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
            eventTuple.GetEntry(current_entry);
            const ntuple::Event& event = eventTuple.data();
            if (/*static_cast<EventEnergyScale>(event.eventEnergyScale) != analysis::EventEnergyScale::Central ||*/
                    event.jets_p4.size() < 2 )
                continue;

            std::vector<LorentzVectorE_Float> jets_corr;
            for(unsigned n = 0; n < event.jets_p4.size(); ++n) {
                LorentzVectorE_Float jet = event.jets_p4.at(n);

                jecUnc->setJetEta(jet.eta());
                jecUnc->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
                const double unc = jecUnc->getUncertainty(true);
                const int sign = -1;
                const double sf = (1.0 + (sign * unc));
                const auto shiftedMomentum = jet * sf;

                jets_corr.push_back(shiftedMomentum);
            }

            LorentzVectorE_Float jet_1 = jets_corr.at(0);
            LorentzVectorE_Float jet_2 = jets_corr.at(1);
            if (!(std::abs(jet_1.eta()) < 2.1 && std::abs(jet_2.eta()) < 2.1)) continue;

            LorentzVectorE_Float bb = jet_1 + jet_2;


            anaData.m_bb().Fill(bb.M());
            anaData.m_sv().Fill(event.SVfit_p4.M());
            anaData.m_kinFit().Fill(event.kinFit_m.at(0));

        } //end loop on entries

    }
};

} //namespace analysis

PROGRAM_MAIN(analysis::JetUncertainty_t, Arguments)
