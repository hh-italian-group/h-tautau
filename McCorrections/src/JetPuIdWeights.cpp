/*! Jet PU Id weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "h-tautau/McCorrections/include/JetPuIdWeights.h"
#include "AnalysisTools/Core/include/RootExt.h"


namespace analysis {
namespace mc_corrections {

JetPuIdWeights::JetPuIdWeights(const std::string& file_eff, const std::string& file_sf,
                               const std::string& file_mistag_eff, const std::string& file_mistag_sf,
                               const BTagger& _bTagger, Period _period) :
    bTagger(_bTagger), period(_period)
{
    static const std::map<Period ,std::string>  period_label = { { Period::Run2016, "2016" },
                                                                 { Period::Run2017, "2017" },
                                                                 { Period::Run2018, "2018" } };
    //Files can be found at: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID#Efficiencies_and_data_MC_scale_f
    const std::string eff = boost::str(boost::format("h2_eff_mc%1%_L") % period_label.at(period));
    const std::string sf = boost::str(boost::format("h2_eff_sf%1%_L") % period_label.at(period));
    const std::string eff_mistag = boost::str(boost::format("h2_mistag_mc%1%_L") % period_label.at(period));
    const std::string sf_mistag = boost::str(boost::format("h2_mistag_sf%1%_L") % period_label.at(period));

    eff_hist = std::shared_ptr<TH2F>(root_ext::ReadCloneObject<TH2F>(*root_ext::OpenRootFile(file_eff), eff, "", true));
    sf_hist = std::shared_ptr<TH2F>(root_ext::ReadCloneObject<TH2F>(*root_ext::OpenRootFile(file_sf), sf, "", true));
    eff_mistag_hist = std::shared_ptr<TH2F>(root_ext::ReadCloneObject<TH2F>(*root_ext::OpenRootFile(file_mistag_eff),
                                            eff_mistag, "", true));
    sf_mistag_hist = std::shared_ptr<TH2F>(root_ext::ReadCloneObject<TH2F>(*root_ext::OpenRootFile(file_mistag_sf),
                                           sf_mistag, "", true));
}

double JetPuIdWeights::GetEfficiency(std::shared_ptr<TH2F> hist, double pt, double eta) const
{
    int xBin = hist->GetXaxis()->FindFixBin(pt);
    xBin = std::min(hist->GetXaxis()->GetNbins(), std::max(1, xBin));
    int yBin = hist->GetYaxis()->FindFixBin(eta);
    yBin = std::min(hist->GetYaxis()->GetNbins(), std::max(1, yBin));
    return hist->GetBinContent(xBin, yBin);
}

double JetPuIdWeights::Get(EventInfo& eventInfo) const
{
    double MC = 1;
    double Data = 1;
    const auto sel_jets = SignalObjectSelector::CreateJetInfos(eventInfo.GetEventCandidate(), bTagger, false,
                                                               eventInfo.GetHttIndex(),
                                                               SignalObjectSelector::SelectedSignalJets());
    for(const auto& sel_jet_info : sel_jets) {
        const auto& jet = eventInfo.GetEventCandidate().GetJets().at(sel_jet_info.index);
        if(!(jet.GetMomentum().pt() < 50)) continue;


        double SF = 1;
        double eff = 1;
        auto index = eventInfo.FindGenMatch(jet);
        if(index.is_initialized() && index == sel_jet_info.index){ //jet from hard interaction, index of the closest gen jet, if found
            SF = GetEfficiency(sf_hist, jet.GetMomentum().pt(), std::abs(jet.GetMomentum().eta()));
            eff = GetEfficiency(eff_hist, jet.GetMomentum().pt(), std::abs(jet.GetMomentum().eta()));
        }
        else{ //jet from PileUp
            SF = GetEfficiency(sf_mistag_hist, jet.GetMomentum().pt(), std::abs(jet.GetMomentum().eta()));
            eff = GetEfficiency(eff_mistag_hist, jet.GetMomentum().pt(), std::abs(jet.GetMomentum().eta()));
        }

        DiscriminatorIdResults jet_pu_id(jet->GetPuId());
        bool jetPuIdOutcome = jet_pu_id.Passed(DiscriminatorWP::Loose);

        MC *= jetPuIdOutcome ? eff : 1 - eff;
        Data *= jetPuIdOutcome ? eff * SF : 1 - eff * SF;
    }

    return MC != 0 ? Data/MC : 0;
}

double JetPuIdWeights::Get(const ntuple::ExpressEvent& /*event*/) const
{
    throw exception("ExpressEvent is not supported in JetPuIdWeights::Get.");
}

} // namespace mc_corrections
} // namespace analysis
