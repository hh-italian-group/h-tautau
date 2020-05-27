/*! Jet PU Id weight.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "h-tautau/McCorrections/include/JetPuIdWeights.h"
#include "AnalysisTools/Core/include/RootExt.h"


namespace analysis {
namespace mc_corrections {

JetInfo::JetInfo(const JetCandidate& jet) :
    eff(0.), SF(0.)
{
    pt  = jet.GetMomentum().pt();
    eta = jet.GetMomentum().eta();
}

JetPuIdWeights::JetPuIdWeights(const std::string& file_eff, const std::string& file_sf, const BTagger& _bTagger,
                               Period _period) :
    bTagger(_bTagger), period(_period)
{
    static const std::map<Period ,std::string>  period_label = { { Period::Run2016, "2016" },
                                                                 { Period::Run2017, "2017" },
                                                                 { Period::Run2018, "2018" } };

    const std::string eff = boost::str(boost::format("h2_eff_mc%1%_L") % period_label.at(period));
    const std::string sf = boost::str(boost::format("h2_eff_sf%1%_L") % period_label.at(period));
    //Files can be found at: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID#Efficiencies_and_data_MC_scale_f
    eff_hist = std::shared_ptr<TH2F>(root_ext::ReadCloneObject<TH2F>(*root_ext::OpenRootFile(file_eff), eff, "", true));
    sf_hist = std::shared_ptr<TH2F>(root_ext::ReadCloneObject<TH2F>(*root_ext::OpenRootFile(file_sf), sf, "", true));
}

double JetPuIdWeights::GetEfficiency(std::shared_ptr<TH2F> hist, double pt, double eta) const
{
    int xBin = hist->GetXaxis()->FindFixBin(pt);
    xBin = std::min(hist->GetXaxis()->GetNbins(), std::max(1, xBin));
    int yBin = hist->GetYaxis()->FindFixBin(eta);
    yBin = std::min(hist->GetYaxis()->GetNbins(), std::max(1, yBin));
    return hist->GetBinContent(xBin, yBin);
}

void JetPuIdWeights::InitializeEff(JetInfo& jetInfo) const
{
    jetInfo.SF  = GetEfficiency(sf_hist, jetInfo.pt, std::abs(jetInfo.eta));
    jetInfo.eff = GetEfficiency(eff_hist, jetInfo.pt, std::abs(jetInfo.eta));
}

double JetPuIdWeights::Get(EventInfo& eventInfo) const
{
    std::vector<JetInfo> jetInfos;
    const auto sel_jets = SignalObjectSelector::CreateJetInfos(eventInfo.GetEventCandidate(), bTagger, false,
                                                               eventInfo.GetHttIndex(),
                                                               SignalObjectSelector::SelectedSignalJets());
    for(const auto& sel_jet_info : sel_jets) {
        const auto& jet = eventInfo.GetEventCandidate().GetJets().at(sel_jet_info.index);
        JetInfo jetInfo(jet);
        if(!(jetInfo.pt > bTagger.PtCut() && std::abs(jetInfo.eta) < bTagger.EtaCut())) continue;
        InitializeEff(jetInfo);
        DiscriminatorIdResults jet_pu_id(eventInfo.GetEventCandidate().GetJets().at(sel_jet_info.index)->GetPuId());
        jetInfo.jetPuIdOutcome = jetInfo.pt < 50 && !jet_pu_id.Passed(DiscriminatorWP::Loose);
        jetInfos.push_back(jetInfo);
    }

    return GetJetPuIdWeight(jetInfos);
}

double JetPuIdWeights::Get(const ntuple::ExpressEvent& /*event*/) const
{
    throw exception("ExpressEvent is not supported in JetPuIdWeights::Get.");
}

double JetPuIdWeights::GetJetPuIdWeight(const std::vector<JetInfo>& jetInfos)
{
    double MC = 1;
    double Data = 1;
    for (const auto& jetInfo : jetInfos) {
        MC *= jetInfo.jetPuIdOutcome ? jetInfo.eff : 1 - jetInfo.eff;
        Data *= jetInfo.jetPuIdOutcome ? jetInfo.eff * jetInfo.SF : 1 - jetInfo.eff * jetInfo.SF;
    }
    return MC != 0 ? Data/MC : 0;
}

} // namespace mc_corrections
} // namespace analysis
