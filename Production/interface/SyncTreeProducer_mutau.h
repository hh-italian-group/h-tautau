/*! Definition of a SyncTree producer for the mu-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

// system include files
#include <memory>
#include <vector>

// user include files

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "h-tautau/Production/interface/BaseEDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "TTree.h"
#include "Math/VectorUtil.h"


namespace analysis {

class SyncAnalyzerData_muTau : public  analysis::BaseEDAnalyzerData {
public:
    SyncAnalyzerData_muTau(std::shared_ptr<TFile> outputFile) : BaseEDAnalyzerData(outputFile) {}
    SyncAnalyzerData_muTau(const std::string& outputFileName) : BaseEDAnalyzerData(outputFileName) {}

    SELECTION_ENTRY(Selection)

    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)
    TH1D_ENTRY(Htautau_Mass, 60, 0.0, 300.0)
};

struct SelectionResultsV2_mutau : public SelectionResultsV2 {
    //finalState::bbMuTaujet muTau_MC;
    CandidateV2Ptr GetMuon() const { return higgs->GetDaughter(CandidateV2::Type::Muon); }
    CandidateV2Ptr GetTau() const { return higgs->GetDaughter(CandidateV2::Type::Tau); }

    virtual CandidateV2Ptr GetLeg(size_t leg_id) const override
    {
        if(leg_id == 1) return GetMuon();
        if(leg_id == 2) return GetTau();
        throw exception("Bad leg id = %1%") % leg_id;
    }

    //virtual const finalState::bbTauTau& GetFinalStateMC() const override { return muTau_MC; }

};

}



//
// class declaration
//

class SyncTreeProducer_mutau: public BaseEDAnalyzer {
   public:
      explicit SyncTreeProducer_mutau(const edm::ParameterSet&);
      ~SyncTreeProducer_mutau();

//      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	protected:
      virtual analysis::Channel ChannelId() const override;
   	  virtual void SelectSignalMuon(const analysis::CandidateV2Ptr& muon, analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override;
   	  virtual void SelectSignalTau(const analysis::CandidateV2Ptr& tau, analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override;
	  virtual void SelectJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override;
        virtual void SelectJetsTight(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override;
      virtual void SelectBJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      virtual analysis::SyncAnalyzerData_muTau& GetAnaData() override { return anaData; }
      virtual analysis::CandidateV2Ptr SelectHiggs(analysis::CandidateV2PtrVector& higgses) override;

      //double SVFit(const analysis::CandidateV2Ptr& higgs, const analysis::MissingETPtr& met);
      void FillSyncTree(const edm::Event& iEvent);
      void GetKinFitResults();
      void ProcessEvent(const edm::Event&, const edm::EventSetup&, const analysis::EventEnergyScale);



      // ----------member data --------------------------

      std::shared_ptr<ntuple::SyncTree> sync_tree;
      ntuple::SyncTree& syncTree;
      analysis::SyncAnalyzerData_muTau anaData;
      analysis::SelectionResultsV2_mutau selection;
      JetCorrectionUncertainty *jecUnc;
};
