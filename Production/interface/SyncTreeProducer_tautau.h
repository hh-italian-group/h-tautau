/*! Definition of a SyncTree producer for the tau-tau channel.
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

class SyncAnalyzerData_tauTau : public analysis::BaseEDAnalyzerData {
public:
    SyncAnalyzerData_tauTau(std::shared_ptr<TFile> outputFile) : BaseEDAnalyzerData(outputFile) {}
    SyncAnalyzerData_tauTau(const std::string& outputFileName) : BaseEDAnalyzerData(outputFileName) {}

    SELECTION_ENTRY(Selection)

    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)
    TH1D_ENTRY(Htautau_Mass, 60, 0.0, 300.0)
};


struct SelectionResultsV2_tautau : public SelectionResultsV2 {
    CandidateV2Ptr GetLeadingTau() const { return higgs->GetLeadingDaughter(CandidateV2::Type::Tau); }
    CandidateV2Ptr GetSubleadingTau() const { return higgs->GetSubleadingDaughter(CandidateV2::Type::Tau); }

    virtual CandidateV2Ptr GetLeg(size_t leg_id) const override
    {
        if(leg_id == 1) return GetLeadingTau();
        if(leg_id == 2) return GetSubleadingTau();
        throw exception("Bad leg id = %1%") % leg_id;
    }
};

}

//
// class declaration
//

class SyncTreeProducer_tautau: public BaseEDAnalyzer {
   public:
      explicit SyncTreeProducer_tautau(const edm::ParameterSet&);
      ~SyncTreeProducer_tautau();

  protected:
      virtual analysis::Channel ChannelId() const override;
      //virtual void SelectSignalElectron(const analysis::CandidateV2Ptr& electron, analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override;
      virtual void SelectSignalTau(const analysis::CandidateV2Ptr& tau, analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override;
      virtual void SelectJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override;
      virtual void SelectBJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut) override;



   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      virtual analysis::SyncAnalyzerData_tauTau& GetAnaData() override { return anaData; }
      virtual analysis::CandidateV2Ptr SelectHiggs(analysis::CandidateV2PtrVector& higgses) override;

      void FillSyncTree(const edm::Event& iEvent);


      // ----------member data --------------------------

      std::shared_ptr<ntuple::SyncTree> sync_tree;
      ntuple::SyncTree& syncTree;
      analysis::SyncAnalyzerData_tauTau anaData;
      analysis::SelectionResultsV2_tautau selection;
};
