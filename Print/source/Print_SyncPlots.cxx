/*! Print control plots that were selected to synchronize produced tree-toople.
Command line to compile standalone version:
  g++ -o Print_SyncPlots -std=c++0x -I$( root-config --incdir ) $( root-config --libs ) -D STANDALONE Print_SyncPlots.cxx

This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <set>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <memory>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TText.h>
#include <TLine.h>
#include <TPad.h>
#include <Rtypes.h>
#include <TError.h>

#include "AnalysisTools/Run/include/program_main.h"

struct EventId {
    unsigned runId;
    unsigned lumiBlock;
    unsigned eventId;
    static const EventId& Undef_event() {
        static const EventId undef_event;
        return undef_event;
    }

    EventId() : runId(std::numeric_limits<UInt_t>::max()), lumiBlock(std::numeric_limits<UInt_t>::max()),
                eventId(std::numeric_limits<UInt_t>::max()){}

    EventId(unsigned _runId, unsigned _lumiBlock, unsigned _eventId) : runId(_runId), lumiBlock(_lumiBlock),
                eventId(_eventId){}

    bool operator == (const EventId& other) const
    {
        return !(*this != other);
    }

    bool operator != (const EventId& other) const
    {
        return runId != other.runId || lumiBlock != other.lumiBlock || eventId != other.eventId;
    }

    bool operator < (const EventId& other) const
    {
        if(runId < other.runId) return true;
        if(runId > other.runId) return false;
        if(lumiBlock < other.lumiBlock) return true;
        if(lumiBlock > other.lumiBlock) return false;
        return eventId < other.eventId;
    }
};

std::ostream& operator <<(std::ostream& s, const EventId& event)
{
    s << "run = " << event.runId << ", lumi =" << event.lumiBlock << ", evt = " << event.eventId;
    return s;
}

struct Arguments {
    REQ_ARG(std::string, channelName);
    REQ_ARG(std::string, sampleName);
    REQ_ARG(std::string, myGroup);
    REQ_ARG(std::string, myRootFile);
    REQ_ARG(std::string, myTTreeName);
    REQ_ARG(std::string, group);
    REQ_ARG(std::string, groupRootFile);
    REQ_ARG(std::string, groupTTreeName);
};

class Print_SyncPlots {
public:
    typedef std::set<EventId> EventSet;
    typedef std::vector<EventId> EventVector;
    typedef std::map<EventId, size_t> EventToEntryMap;
    typedef std::pair<size_t, size_t> EntryPair;
    typedef std::map<EventId, EntryPair> EventToEntryPairMap;

    Print_SyncPlots(const Arguments& args)
        : channel(args.channelName()), sample(args.sampleName()), myGroup(args.myGroup()),
          myRootFile(args.myRootFile()), myTree(args.myTTreeName()),
          group(args.group()), groupRootFile(args.groupRootFile()), groupTree(args.groupTTreeName()),
          isFirstPage(true), isLastDraw(false)
    {
        std::cout << channel << " " << sample << std::endl;
        std::cout << myGroup << "  " << myRootFile << "  " << myTree << std::endl;
        std::cout << group << "  " << groupRootFile << "  " << groupTree << std::endl;

        Tmine = LoadTree(Fmine, myRootFile, myTree);
        Tother = LoadTree(Fother, groupRootFile, groupTree);
        my_eventES = CollectValues<Int_t>(Tmine, "eventEnergyScale");

        CollectEvents(Tmine, Tother);

        std::cout << "Mine: " << Fmine->GetName() << std::endl;
        std::cout << "Other: " << Fother->GetName() << std::endl;

        file_name = std::string("PlotsDiff_") + channel + "_" + sample + "_" + myGroup + "_" + group + ".pdf";

        gErrorIgnoreLevel = kWarning;
    }

    void Run()
    {
        using namespace std;

//        drawHistos("run", 200, 0, 200000);
//        drawHistos("lumi", 2000, 0, 2000);
//        drawHistos("evt", 2000, 0, 2000);

        drawHistos("npv", 50, 0, 50);
        drawHistos("npu", 50, 0, 50);
        //drawHistos("rho", 50, 0, 50);

        //event weights
        //drawHistos("mcweight", 40, 0, 2);
//        drawHistos("puweight", 40, 0, 2);
//        drawHistos("trigweight_1", 50, 0.90, 1.1);
//        drawHistos("trigweight_2", 50, 0.90, 1.1);
//        if (channel == "etau")
//            drawHistos("idweight_1", 50, 0.75, 1.);
//        if (channel == "mutau")
//            drawHistos("idweight_1", 50, 0.95, 1.);
//        //drawHistos("idweight_2", 50, 0.95, 1.);
//        if (channel == "etau" || channel == "mutau")
//            drawHistos("isoweight_1", 100, 0.90, 1.);
//        //drawHistos("isoweight_2", 100, 0.90, 1.);
//        drawHistos("fakeweight", 25, 0.5, 1.5);
//        drawHistos("etau_fakerate", 25, 0.5, 2.5);

////        drawHistos("effweight", 30, .8, 1.1);
//        drawHistos("weight", 30, .0, 1.2);
////        drawHistos("weight", 30, .8, 1.1);
////        drawHistos("embeddedWeight", 50, 0, 1);
//        drawHistos("decayModeWeight_1", 30, 0, 2.0);
//        drawHistos("decayModeWeight_2", 30, 0, 2.0);

        //SV Fit variables
        drawHistos("m_vis", 20, 0, 200);
        drawHistos("m_sv", 30, 0, 300);
        //drawHistos("pt_sv", 100, 0, 100);
        //drawHistos("eta_sv", 60, -3, 3);
        //drawHistos("phi_sv", 30, -3.5, 3.5);
//        drawHistos("m_sv_Up", 60, 0, 300);
//        drawHistos("m_sv_Down", 60, 0, 300);

//        //  First lepton :  muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau


        drawHistos("pt_1", 15, 0, 150);
        drawHistos("eta_1", 60, -3, 3);
        drawHistos("phi_1", 30, -3.5, 3.5);
////        if (channel == "etau" || channel == "mutau")
        drawHistos("iso_1", 40, 0., 0.12);
////        if (channel == "tautau")
////            drawHistos("iso_1", 100, 0., 1);
//        drawHistos("m_1", 20, 0, 0.4);
        drawHistos("q_1", 4, -2, 2);
////        if (channel == "etau")
////            drawHistos("mva_1", 30, 0.8, 1.1);
//        drawHistos("d0_1", 40, -0.1, 0.1);
//        drawHistos("dZ_1", 40, -0.5, 0.5);
////        drawHistos("passid_1", 2, -0.5, 1.5);
////        drawHistos("passiso_1", 2, -0.5, 1.5);
        drawHistos("pfmt_1", 50, 0, 200);
////        drawHistos("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", 40, 0, 20);
////        drawHistos("againstElectronMVA3raw_1", 2, -0.5, 1.5);
////        drawHistos("byIsolationMVA2raw_1", 300, -1, 1);
////        drawHistos("againstMuonLoose2_1", 2, -0.5, 1.5);
////        drawHistos("againstMuonMedium2_1", 2, -0.5, 1.5);
////        drawHistos("againstMuonTight2_1", 2, -0.5, 1.5);

//        //Second lepton :  hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT)  Tau for Tau Tau
        drawHistos("pt_2", 15, 0, 150);
        drawHistos("eta_2", 60, -3, 3);
        drawHistos("phi_2", 30, -3.5, 3.5);
        drawHistos("iso_2", 50, -1, 1);
//        drawHistos("m_2", 40, 0, 2);
        drawHistos("q_2", 4, -2, 2);
////        drawHistos("d0_2", 40, -0.1, 0.1);
////        drawHistos("dZ_2", 40, -0.5, 0.5);
////        drawHistos("mva_2", 20, -0.5, 1.5);
////        drawHistos("passid_2", 2, -0.5, 1.5);
////        drawHistos("passiso_2", 2, -0.5, 1.5);
        drawHistos("pfmt_2", 50, 0, 200);

        drawHistos("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", 40, 0, 2);
////        drawHistos("againstElectronMVA3raw_2", 2, -0.5, 1.5);
////        drawHistos("byIsolationMVA2raw_2", 40, 0.5, 1.2);
////        drawHistos("againstMuonLoose2_2", 2, -0.5, 1.5);
////        drawHistos("againstMuonMedium2_2", 2, -0.5, 1.5);
////        drawHistos("againstMuonTight2_2", 2, -0.5, 1.5);

//        //Dilepton
//        drawHistos("pt_tt", 100, 0, 200);

//        //Met related variables

        drawHistos("met", 20, 0, 200);
        drawHistos("metphi", 30, -3.5, 3.5);
        drawHistos("metcov00", 40, 0, 1000);
        drawHistos("metcov01", 40, 0, 1000);
        drawHistos("metcov10", 40, 0, 1000);
        drawHistos("metcov11", 40, 0, 1000);

////        drawHistos("mvamet", 20, 0, 200);
////        drawHistos("mvametphi", 30, -3.5, 3.5);
////        drawHistos("mvacov00", 40, 0, 1000);
////        drawHistos("mvacov01", 40, 0, 1000);
////        drawHistos("mvacov10", 40, 0, 1000);
////        drawHistos("mvacov11", 40, 0, 1000);

////        drawHistos("pzetavis", 300, 0, 300);
////        drawHistos("pzetamiss", 300, 0, 300);

//        std::cout << "after met" << std::endl;

//        // Jets
//        const std::vector<Int_t> my_njets = CollectValues<Int_t>(Tmine, "njets");
//        const std::vector<Int_t> other_njets = CollectValues<Int_t>(Tother, "njets");
////        const auto noJets_my = [&](size_t entry_id) -> bool {
////            return my_njets.at(entry_id) == 0;
////        };
//        const auto atLeast1jet_my = [&](size_t entry_id) -> bool {
//            return my_njets.at(entry_id) >= 1;
//        };
//        const auto atLeast2jets_my = [&](size_t entry_id) -> bool {
//            return my_njets.at(entry_id) >= 2;
//        };
////        const auto noJets_other = [&](size_t entry_id) -> bool {
////            return other_njets.at(entry_id) == 0;
////        };
//        const auto atLeast1jet_other = [&](size_t entry_id) -> bool {
//            return other_njets.at(entry_id) >= 1;
//        };
//        const auto atLeast2jets_other = [&](size_t entry_id) -> bool {
//            return other_njets.at(entry_id) >= 2;
//        };

//        //First Jet   : leading jet after applying Jet energy corrections (excluding hadronic Tau)
//        drawHistos("jpt_1", 50, 0, 300);
//        drawHistos("jeta_1", 50, -5, 5);
//        drawHistos("jphi_1", 30, -3.5, 3.5);
//        drawHistos("jptraw_1", 50, 0, 300);
//        //drawHistos("jptunc_1", 50, 0, 300, atLeast1jet_my, atLeast1jet_other, "njets>=1");
//        drawHistos("jmva_1", 40, 0.5, 1.1, atLeast1jet_my, atLeast1jet_other, "njets>=1");
////        drawHistos("jlrm_1", 120, -60, 60, atLeast1jet_my, atLeast1jet_other, "njets>=1");
////        drawHistos("jctm_1", 120, -60, 60, atLeast1jet_my, atLeast1jet_other, "njets>=1");
////        drawHistos("jpass_1" , 2, -0.5, 1.5, atLeast1jet_my, atLeast1jet_other, "njets>=1");

//            //Second Jet  : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
//        drawHistos("jpt_2", 50, 0, 300);
//        drawHistos("jeta_2", 50, -5, 5);
//        drawHistos("jphi_2", 30, -3.5, 3.5);
//        drawHistos("jptraw_2", 50, 0, 300);
////        drawHistos("jptunc_2", 50, 0, 300, atLeast2jets_my, atLeast2jets_other, "njets>=2");
//        drawHistos("jmva_2", 40, 0.5, 1.1, atLeast2jets_my, atLeast2jets_other, "njets>=2");
////        drawHistos("jlrm_2", 120, -60, 60, atLeast2jets_my, atLeast2jets_other, "njets>=2");
////        drawHistos("jctm_2", 120, -60, 60, atLeast2jets_my, atLeast2jets_other, "njets>=2");
////        drawHistos("jpass_2" , 2, -0.5, 1.5, atLeast2jets_my, atLeast2jets_other, "njets>=2");


//            //number of jets passing jet id ( pt > 30 )
        drawHistos("njets", 7, -0.5, 6.5);
        drawHistos("njetspt20", 7, -0.5, 6.5);


//        // b-jets
//        const std::vector<Int_t> my_nbjets = CollectValues<Int_t>(Tmine, "nbtag");
//        const std::vector<Int_t> other_nbjets = CollectValues<Int_t>(Tother, "nbtag");
//        const auto atLeast1bjet_my = [&](size_t entry_id) -> bool {
//            return my_nbjets.at(entry_id) >= 1;
//        };
//        const auto atLeast1bjet_other = [&](size_t entry_id) -> bool {
//            return other_nbjets.at(entry_id) >= 1;
//        };

//        const auto atLeast2bjets_my = [&](size_t entry_id) -> bool {
//            return my_nbjets.at(entry_id) >= 2;
//        };
//        const auto atLeast2bjets_other = [&](size_t entry_id) -> bool {
//            return other_nbjets.at(entry_id) >= 2;
//        };

//        const auto atLeast3bjets_my = [&](size_t entry_id) -> bool {
//            return my_nbjets.at(entry_id) >= 3;
//        };
//        const auto atLeast3bjets_other = [&](size_t entry_id) -> bool {
//            return other_nbjets.at(entry_id) >= 3;
//        };


//            //Candidate B Jets : leading jet (in CSV ordering) passing (pt > 20 + eta < 2.4)
//        drawHistos("bpt_1", 20, 0, 200);
//        drawHistos("beta_1", 50, -5, 5);
//        drawHistos("bphi_1", 30, -3.5, 3.5);
//        drawHistos("bcsv_1", 50, 0.5, 1.1);

//            //Candidate B Jets : subleading jet (in CSV ordering) passing (pt > 20 + eta < 2.4)
//        drawHistos("bpt_2", 20, 0, 200);
//        drawHistos("beta_2", 50, -5, 5);
//        drawHistos("bphi_2", 30, -3.5, 3.5);
//        drawHistos("bcsv_2", 50, 0.5, 1.1);
//        drawHistos("m_bb", 30, 0, 300, atLeast2bjets_my, atLeast2bjets_other, "nbtag>=2");
//        drawHistos("m_ttbb", 50, 0, 500, atLeast2bjets_my, atLeast2bjets_other, "nbtag>=2");

//            //Candidate B Jets : third jet (in CSV ordering) passing (pt > 20 + eta < 2.4)
//        drawHistos("bpt_3", 20, 0, 200, atLeast3bjets_my, atLeast3bjets_other, "nbtag>=3");
//        drawHistos("beta_3", 50, -5, 5, atLeast3bjets_my, atLeast3bjets_other, "nbtag>=3");
//        drawHistos("bphi_3", 30, -3.5, 3.5, atLeast3bjets_my, atLeast3bjets_other, "nbtag>=3");
//        drawHistos("bcsv_3", 50, 0.5, 1.1, atLeast3bjets_my, atLeast3bjets_other, "nbtag>=3");

//            //number of btags passing btag id (medium CSV WP) ( pt > 20 )
        drawHistos("dilepton_veto", 2, 0, 2);
        drawHistos("extraelec_veto", 2, 0, 2);
        drawHistos("extramuon_veto", 2, 0, 2);
        isLastDraw = true;
        drawHistos("nbtag", 5, -0.5, 4.5);
  //drawHistos("dilepton_veto", 2, 0, 2);
    }

private:
    static TTree* LoadTree(std::shared_ptr<TFile>& file, const std::string& fileName, const std::string& treeName)
    {
        file = std::shared_ptr<TFile>(new TFile(fileName.c_str(), "READ"));
        TTree* tree = (TTree*)file->Get(treeName.c_str());
        if(!tree) {
            std::ostringstream ss;
            ss << "File " << fileName << " is empty.";
            throw std::runtime_error(ss.str());

        }
        tree->SetBranchStatus("*", 0);
        return tree;
    }

    void drawHistos(const std::string& var, int nbins, float xmin, float xmax)
    {
        const auto SelectAll = [](size_t entry_id) -> bool { return true; };
        const auto SelectES = [&](size_t entry_id) -> bool { return my_eventES.at(entry_id) == 0; };
        drawHistos(var, nbins, xmin, xmax, SelectES, SelectAll, "All");
    }

    template<typename MySelector, typename OtherSelector>
    void drawHistos(const std::string& var, int nbins, float xmin, float xmax, const MySelector& my_selector,
                    const OtherSelector& other_selector, const std::string& selection_label)
    {
        try {

            std::shared_ptr<TH1F> Hmine_all(new TH1F(TString("Hmine") + var + "all","",nbins,xmin,xmax));

            std::shared_ptr<TH1F> Hother_all(new TH1F(TString("Hother")+var + "all","",nbins,xmin,xmax));

            std::shared_ptr<TH1F> Hmine_common(new TH1F(TString("Hmine")+var + "common","",nbins,xmin,xmax));

            std::shared_ptr<TH1F> Hother_common(new TH1F(TString("Hother")+var + "common","",nbins,xmin,xmax));

            std::shared_ptr<TH1F> Hmine_diff(new TH1F(TString("Hmine")+var + "diff","",nbins,xmin,xmax));

            std::shared_ptr<TH1F> Hother_diff(new TH1F(TString("Hother")+var + "diff","",nbins,xmin,xmax));

            std::shared_ptr<TH2F> Hmine_vs_other(new TH2F(TString("Hmine_vs_other") + var, "", nbins, xmin, xmax,
                                                          40, -1.0, 1.0));

            TBranch* myBranch = Tmine->GetBranch(var.c_str());
            if (!myBranch){
                std::ostringstream ss;
                ss << "My Branch '" << var << "' is not found.";
                throw std::runtime_error(ss.str());
            }
            TBranch* otherBranch = Tother->GetBranch(var.c_str());
            if (!otherBranch){
                std::ostringstream ss;
                ss << "Other Branch '" << var << "' is not found.";
                throw std::runtime_error(ss.str());
            }
            TClass *myClass, *otherClass;
            EDataType myType, otherType;
            myBranch->GetExpectedType(myClass, myType);
            otherBranch->GetExpectedType(otherClass, otherType);
            if(myClass || otherClass) {
                std::ostringstream ss;
                ss << "branches with complex objects are not supported for branch '" << var << "'.";
                throw std::runtime_error(ss.str());
            }
            if(myType == kFloat_t && otherType == kFloat_t)
                FillAllHistograms<Float_t, Float_t>(var, my_selector, other_selector, *Hmine_all, *Hother_all,
                                             *Hmine_common, *Hother_common, *Hmine_vs_other, *Hmine_diff, *Hother_diff);
            else if(myType == kDouble_t && otherType == kDouble_t)
                FillAllHistograms<Double_t, Double_t>(var, my_selector, other_selector, *Hmine_all, *Hother_all,
                                             *Hmine_common, *Hother_common, *Hmine_vs_other, *Hmine_diff, *Hother_diff);
            else if(myType == kDouble_t && otherType == kFloat_t)
                FillAllHistograms<Double_t, Float_t>(var, my_selector, other_selector, *Hmine_all, *Hother_all,
                                             *Hmine_common, *Hother_common, *Hmine_vs_other, *Hmine_diff, *Hother_diff);
            else if(myType == kFloat_t && otherType == kDouble_t)
                FillAllHistograms<Float_t, Double_t>(var, my_selector, other_selector, *Hmine_all, *Hother_all,
                                             *Hmine_common, *Hother_common, *Hmine_vs_other, *Hmine_diff, *Hother_diff);
            else if(myType == kInt_t && otherType == kInt_t)
                FillAllHistograms<Int_t, Int_t>(var, my_selector, other_selector, *Hmine_all, *Hother_all,
                                             *Hmine_common, *Hother_common, *Hmine_vs_other, *Hmine_diff, *Hother_diff);
            else if(myType == kInt_t && otherType == kDouble_t)
                FillAllHistograms<Int_t, Double_t>(var, my_selector, other_selector, *Hmine_all, *Hother_all,
                                             *Hmine_common, *Hother_common, *Hmine_vs_other, *Hmine_diff, *Hother_diff);
            else if(myType == kDouble_t && otherType == kInt_t)
                FillAllHistograms<Double_t, Int_t>(var, my_selector, other_selector, *Hmine_all, *Hother_all,
                                             *Hmine_common, *Hother_common, *Hmine_vs_other, *Hmine_diff, *Hother_diff);
            else if(myType == kBool_t && otherType == kBool_t)
                FillAllHistograms<Bool_t, Bool_t>(var, my_selector, other_selector, *Hmine_all, *Hother_all,
                                             *Hmine_common, *Hother_common, *Hmine_vs_other, *Hmine_diff, *Hother_diff);
            else if(myType == kBool_t && otherType == kChar_t)
                FillAllHistograms<Bool_t, Char_t>(var, my_selector, other_selector, *Hmine_all, *Hother_all,
                                             *Hmine_common, *Hother_common, *Hmine_vs_other, *Hmine_diff, *Hother_diff);
            else {
                std::ostringstream ss;
                ss << "Unknown branch type combination (" <<  myType << ", " << otherType
                   << ") for branch '" << var << "'.";
                throw std::runtime_error(ss.str());
            }
            DrawSuperimposedHistograms(Hmine_all, Hother_all, selection_label + " (all)", var);
            DrawSuperimposedHistograms(Hmine_common, Hother_common, selection_label + " (common)", var);
            DrawSuperimposedHistograms(Hmine_diff, Hother_diff, selection_label + " (diff)", var);
            Draw2DHistogram(Hmine_vs_other, selection_label, var);

        } catch(std::runtime_error& e){
            std::cerr << "WARNING: " << e.what() << std::endl;
        }
    }

    void DrawSuperimposedHistograms(std::shared_ptr<TH1F> Hmine, std::shared_ptr<TH1F> Hother,
                                    const std::string& selection_label, const std::string& var)
    {
        Hmine->SetTitle(selection_label.c_str());
        Hmine->GetYaxis()->SetTitle(selection_label.c_str());
        Hmine->GetXaxis()->SetTitle(var.c_str());
        Hmine->SetLineColor(1);
        Hmine->SetMarkerColor(1);
        Hmine->SetStats(0);

        Hother->GetYaxis()->SetTitle(selection_label.c_str());
        Hother->GetXaxis()->SetTitle(var.c_str());
        Hother->SetLineColor(2);
        Hother->SetMarkerColor(2);
        Hother->SetStats(0);

        TPad pad1("pad1","",0,0.2,1,1);
        TPad pad2("pad2","",0,0,1,0.2);

        pad1.cd();

        // Draw one histogram on top of the other
        if(Hmine->GetMaximum()>Hother->GetMaximum())
            Hmine->GetYaxis()->SetRangeUser(0,Hmine->GetMaximum()*1.1);
        else
            Hmine->GetYaxis()->SetRangeUser(0,Hother->GetMaximum()*1.1);
        Hmine->Draw("hist");
        Hother->Draw("histsame");
        DrawTextLabels(Hmine->Integral(0,Hmine->GetNbinsX()+1), Hother->Integral(0,Hother->GetNbinsX()+1));

        pad2.cd();

        // Draw the ratio of the historgrams
        std::unique_ptr<TH1F> HDiff((TH1F*)Hother->Clone("HDiff"));
        HDiff->Divide(Hmine.get());
        ///HDiff->GetYaxis()->SetRangeUser(0.9,1.1);
        HDiff->GetYaxis()->SetRangeUser(0.9,1.1);
        //HDiff->GetYaxis()->SetRangeUser(0.98,1.02);
        //HDiff->GetYaxis()->SetRangeUser(0.,2.0);
        HDiff->GetYaxis()->SetNdivisions(3);
        HDiff->GetYaxis()->SetLabelSize(0.1);
        HDiff->GetYaxis()->SetTitleSize(0.1);
        HDiff->GetYaxis()->SetTitleOffset(0.5);
        //HDiff->GetYaxis()->SetTitle(myGroup + " / " + group);
        HDiff->GetYaxis()->SetTitle("Ratio");
        HDiff->GetXaxis()->SetNdivisions(-1);
        HDiff->GetXaxis()->SetTitle("");
        HDiff->GetXaxis()->SetLabelSize(0.0001);
        HDiff->SetMarkerStyle(7);
        HDiff->SetMarkerColor(2);
        HDiff->Draw("histp");
        TLine line;
        line.DrawLine(HDiff->GetXaxis()->GetXmin(),1,HDiff->GetXaxis()->GetXmax(),1);

        canvas.Clear();
        pad1.Draw();
        pad2.Draw();

        PrintCanvas(var);
    }

    void Draw2DHistogram(std::shared_ptr<TH2F> Hmine_vs_other, const std::string& selection_label,
                         const std::string& var)
    {
        Hmine_vs_other->SetTitle(selection_label.c_str());
        const std::string my_name = var + "_mine";
        const std::string other_name = var + "_other";
        std::ostringstream y_name;
        y_name << "(" << other_name << " - " << my_name << ")/" << other_name;
        Hmine_vs_other->GetXaxis()->SetTitle(my_name.c_str());
        Hmine_vs_other->GetYaxis()->SetTitle(y_name.str().c_str());

        TPad pad1("pad1","", 0, 0, 1, 1);
        pad1.cd();
        Hmine_vs_other->Draw("colz");
        const size_t n_events = Hmine_vs_other->Integral(0, Hmine_vs_other->GetNbinsX() + 1,
                                                         0, Hmine_vs_other->GetNbinsY() + 1);
        DrawTextLabels(n_events, n_events);
        canvas.Clear();
        pad1.Draw();
        PrintCanvas(var, isLastDraw);
    }

    void DrawTextLabels(size_t n_events_mine, size_t n_events_other)
    {
        TText TXmine;
        TXmine.SetTextColor(1);
        TXmine.SetTextSize(.04);
        TText TXother;
        TXother.SetTextColor(2);
        TXother.SetTextSize(.04);

        //Print the integrals of the histograms a the top
        //TXmine.DrawTextNDC(.2,.965,myGroup+"_"+myRootFile+": "+(long)(Hmine->Integral(0,Hmine->GetNbinsX()+1)));
        //TXother.DrawTextNDC(.2,.93,group+"_"+groupRootFile+": "+(long)(Hother->Integral(0,Hother->GetNbinsX()+1)));
        TXmine.DrawTextNDC(.23,.84,myGroup+" : " + n_events_mine);
        TXother.DrawTextNDC(.53,.84,group+": " + n_events_other);
    }

    template<typename MyVarType, typename OtherVarType, typename MySelector, typename OtherSelector,
             typename Histogram, typename Histogram2D>
    void FillAllHistograms(const std::string& var, const MySelector& my_selector, const OtherSelector& other_selector,
                           Histogram& Hmine_all, Histogram& Hother_all,
                           Histogram& Hmine_common, Histogram& Hother_common, Histogram2D& Hmine_vs_other,
                           Histogram& Hmine_diff, Histogram& Hother_diff)
    {
        const std::vector<MyVarType> my_values = CollectValues<MyVarType>(Tmine, var);
        const std::vector<OtherVarType> other_values = CollectValues<OtherVarType>(Tother, var);
        FillCommonHistograms(var, my_values, other_values, my_selector, other_selector,
                             Hmine_common, Hother_common, Hmine_vs_other);
        FillInclusiveHistogram(my_values, my_selector, Hmine_all);
        FillInclusiveHistogram(other_values, other_selector, Hother_all);

        FillExclusiveHistogram(my_values, other_values, my_selector, other_selector, Hmine_diff, Hother_diff);
    }

    template<typename MyVarType, typename OtherVarType, typename MySelector, typename OtherSelector,
             typename Histogram, typename Histogram2D>
    void FillCommonHistograms(const std::string& var, const std::vector<MyVarType>& my_values,
                              const std::vector<OtherVarType>& other_values,
                              const MySelector& my_selector, const OtherSelector& other_selector,
                              Histogram& my_histogram, Histogram& other_histogram, Histogram2D& histogram2D)
    {
        std::cout << var << " bad events:\n";
        for(const auto& event_entry_pair : common_event_to_entry_pair_map) {
            const size_t my_entry = event_entry_pair.second.first;
            const size_t other_entry = event_entry_pair.second.second;
            if(!my_selector(my_entry) || !other_selector(other_entry))
                continue;
            const MyVarType& my_value = my_values.at(my_entry);
            const OtherVarType& other_value = other_values.at(other_entry);
            my_histogram.Fill(my_value);
            other_histogram.Fill(other_value);
            if(other_value) {
                const auto y_value = (other_value - my_value)/other_value;
                const auto diff = other_value - my_value;
                if (std::abs(y_value) >= 0.1 )
                    std::cout << "event = " << event_entry_pair.first << ", other - my = " << diff << std::endl;
                histogram2D.Fill(my_value, y_value);
            }
        }
    }

    template<typename VarType, typename Histogram, typename Selector>
    void FillInclusiveHistogram(const std::vector<VarType>& values, const Selector& selector, Histogram& histogram)
    {
        for(size_t n = 0; n < values.size(); ++n) {
            if(selector(n))
                histogram.Fill(values.at(n));
        }
    }

    template<typename MyVarType, typename OtherVarType, typename Histogram, typename MySelector, typename OtherSelector>
    void FillExclusiveHistogram(const std::vector<MyVarType>& my_values,
                                const std::vector<OtherVarType>& other_values,
                                const MySelector& my_selector, const OtherSelector& other_selector,
                                Histogram& my_histogram, Histogram& other_histogram)
    {
        for (const auto& event_entry_pair : common_event_to_entry_pair_map){
            const size_t my_entry = event_entry_pair.second.first;
            const size_t other_entry = event_entry_pair.second.second;
            if(my_selector(my_entry) && !other_selector(other_entry)){
                const MyVarType& my_value = my_values.at(my_entry);
                my_histogram.Fill(my_value);
            }
            if(!my_selector(my_entry) && other_selector(other_entry)){
                const OtherVarType& other_value = other_values.at(other_entry);
                other_histogram.Fill(other_value);
            }
        }
        for(const auto& event_entry : my_events_only_map) {
            if(my_selector(event_entry.second)) {
                const MyVarType& value = my_values.at(event_entry.second);
                my_histogram.Fill(value);
            }
        }
        for(const auto& event_entry : other_events_only_map) {
            if(other_selector(event_entry.second)) {
                const OtherVarType& value = other_values.at(event_entry.second);
                other_histogram.Fill(value);
            }
        }
    }

    void CollectEvents(TTree* my_tree, TTree* other_tree)
    {
        my_events = CollectEventIds(my_tree);
        const EventSet my_events_set(my_events.begin(), my_events.end());
        other_events = CollectEventIds(other_tree);
        const EventSet other_events_set(other_events.begin(), other_events.end());

        EventSet intersection, my_events_only, other_events_only;
        EventSetIntersection(my_events_set, other_events_set, intersection);
        EventSetDifference(my_events_set, other_events_set, my_events_only);
        EventSetDifference(other_events_set, my_events_set, other_events_only);

        EventToEntryMap my_event_to_entry_map, other_event_to_entry_map;
        FillEventToEntryMap(my_events, my_event_to_entry_map, true);
        FillEventToEntryMap(other_events, other_event_to_entry_map);

        std::cout << "Mine events" << std::endl;
        for(const auto& event_entry : my_event_to_entry_map) {
            if(intersection.count(event_entry.first)) {
                common_event_to_entry_pair_map[event_entry.first] =
                        std::pair<size_t, size_t>(event_entry.second, other_event_to_entry_map.at(event_entry.first));
            }
            if(my_events_only.count(event_entry.first)) {
                my_events_only_map[event_entry.first] = event_entry.second;
                std::cout << "eventId = " << event_entry.first << std::endl;
            }
        }

        std::cout << "Common events" << std::endl;
        for (const auto common_event_entry : common_event_to_entry_pair_map){
            std::cout << "eventId = " << common_event_entry.first << std::endl;
        }

        std::cout << "Other's events" << std::endl;
        for(const auto& event_entry : other_event_to_entry_map) {
            if(other_events_only.count(event_entry.first)){
                other_events_only_map[event_entry.first] = event_entry.second;
                std::cout << "eventId = " << event_entry.first << std::endl;
            }
        }

        std::cout << "# my events = " << my_events.size() << ", " << "# my unique events = " << my_events_set.size()
                  << "\n# other events = " << other_events.size()
                  << ", # other unique events = " << other_events_set.size()
                  << "\n# common events = " << intersection.size() << std::endl;
    }

    static void EventSetIntersection(const EventSet& first_set, const EventSet& second_set, EventSet& intersection_set)
    {
        const size_t max_intersection_size = std::max(first_set.size(), second_set.size());
        EventVector intersection_vector(max_intersection_size);
        const auto iter = std::set_intersection(first_set.begin(), first_set.end(),
                                                second_set.begin(), second_set.end(),
                                                intersection_vector.begin());
        intersection_vector.resize(iter - intersection_vector.begin());
        intersection_set.clear();
        intersection_set.insert(intersection_vector.begin(), intersection_vector.end());
    }

    static void EventSetDifference(const EventSet& first_set, const EventSet& second_set, EventSet& diff_set)
    {
        EventVector diff_vector(first_set.size());
        const auto iter = std::set_difference(first_set.begin(), first_set.end(),
                                              second_set.begin(), second_set.end(),
                                              diff_vector.begin());
        diff_vector.resize(iter - diff_vector.begin());
        diff_set.clear();
        diff_set.insert(diff_vector.begin(), diff_vector.end());
    }

    void FillEventToEntryMap(const EventVector& events, EventToEntryMap& event_to_entry_map, const bool filterES = false) const
    {
        for(size_t n = 0; n < events.size(); ++n) {
             if(filterES && my_eventES.at( n ) !=0) continue;
             event_to_entry_map[events[n]] = n;
        }
    }

    template<typename VarType>
    std::vector<VarType> CollectValues(TTree* tree, const std::string& name)
    {
        EnableBranch(tree, name, true);
        std::vector<VarType> result;
        VarType value;
        tree->SetBranchAddress(name.c_str(), &value);
        const Long64_t N = tree->GetEntries();
        for(Long64_t n = 0; n < N;++n) {
            if(tree->GetEntry(n) < 0)
                throw std::runtime_error("error while reading tree.");
            result.push_back(value);
        }
        EnableBranch(tree, name, false);
        return result;
    }

    template<typename VarType>
    std::vector<VarType> CollectValuesEx(TTree* tree, const std::string& name)
    {
        std::vector<VarType> result;

        TBranch* branch = tree->GetBranch(name.c_str());
        if (!branch){
            std::ostringstream ss;
            ss << "Branch '" << name << "' not found.";
            throw std::runtime_error(ss.str());
        }
        TClass *branch_class;
        EDataType branch_type;
        branch->GetExpectedType(branch_class, branch_type);
        if(branch_class) {
            std::ostringstream ss;
            ss << "branches with complex objects are not supported for branch '" << name << "'.";
            throw std::runtime_error(ss.str());
        }
        if(branch_type == kInt_t) {
            const std::vector<Int_t> original_result = CollectValues<Int_t>(tree, name);
            result.resize(original_result.size());
            std::copy(original_result.begin(), original_result.end(), result.begin());
        }
        if(branch_type == kUInt_t) {
            const std::vector<UInt_t> original_result = CollectValues<UInt_t>(tree, name);
            result.resize(original_result.size());
            std::copy(original_result.begin(), original_result.end(), result.begin());
        }
        if(branch_type == kULong64_t) {
            const std::vector<ULong64_t> original_result = CollectValues<ULong64_t>(tree, name);
            result.resize(original_result.size());
            std::copy(original_result.begin(), original_result.end(), result.begin());
        }
        if(branch_type == kFloat_t) {
            const std::vector<Float_t> original_result = CollectValues<Float_t>(tree, name);
            result.resize(original_result.size());
            std::copy(original_result.begin(), original_result.end(), result.begin());
        }
        if(branch_type == kDouble_t) {
            const std::vector<Double_t> original_result = CollectValues<Double_t>(tree, name);
            result.resize(original_result.size());
            std::copy(original_result.begin(), original_result.end(), result.begin());
        }
        if(branch_type == kChar_t) {
            const std::vector<Char_t> original_result = CollectValues<Char_t>(tree, name);
            result.resize(original_result.size());
            std::copy(original_result.begin(), original_result.end(), result.begin());
        }
        if(branch_type == kBool_t) {
            const std::vector<Bool_t> original_result = CollectValues<Bool_t>(tree, name);
            result.resize(original_result.size());
            std::copy(original_result.begin(), original_result.end(), result.begin());
        }

        return result;
    }


    EventVector CollectEventIds(TTree* tree)
    {
        const std::vector<UInt_t> run = CollectValuesEx<UInt_t>(tree, "run");
        const std::vector<UInt_t> lumi = CollectValuesEx<UInt_t>(tree, "lumi");
        const std::vector<UInt_t> evt = CollectValuesEx<UInt_t>(tree, "evt");

        EventVector events;
        for(size_t n = 0; n < evt.size(); ++n)
            events.push_back(EventId(run.at(n), lumi.at(n), evt.at(n)));
        return events;
    }

    void EnableBranch(TTree* tree, const std::string& name, bool enable)
    {
        UInt_t n_found = 0;
        tree->SetBranchStatus(name.c_str(), enable, &n_found);
        if(n_found != 1) {
            std::ostringstream ss;
            ss << "Branch '" << name << "' is not found.";
            throw std::runtime_error(ss.str());
        }
    }

    void PrintCanvas(const std::string& page_name, bool isLastPage = false)
    {
        std::ostringstream print_options, output_name;
        print_options << "Title:" << page_name;
        output_name << file_name;
        if(isFirstPage && !isLastPage)
            output_name << "(";
        else if(isLastPage && !isFirstPage)
            output_name << ")";
        isFirstPage = false;
        canvas.Print(output_name.str().c_str(), print_options.str().c_str());
    }

private:
    std::string channel, sample, myGroup, myRootFile, myTree, group, groupRootFile, groupTree;
    std::shared_ptr<TFile> Fmine, Fother;
    TTree *Tmine, *Tother;
    std::vector<Int_t> my_eventES;

    EventVector my_events, other_events;
    EventToEntryMap my_events_only_map, other_events_only_map;
    EventToEntryPairMap common_event_to_entry_pair_map;

    TCanvas canvas;
    std::string file_name;
    bool isFirstPage, isLastDraw;
};

PROGRAM_MAIN(Print_SyncPlots, Arguments)
