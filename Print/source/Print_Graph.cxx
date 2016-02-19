/*! Print TGraph.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <TTree.h>
#include <TGraph.h>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Print/include/RootPrintToPdf.h"
#include "Analysis/include/Htautau_TriggerEfficiency.h"

class Print_Graph {
public:
    typedef std::pair< std::string, std::string > FileNameTagPair;
    typedef std::vector< FileNameTagPair> FileNameTagVector;
    typedef std::shared_ptr<TFile> FilePtr;
    typedef std::vector< FilePtr > FileVector;

    template<typename ...Args>
    Print_Graph(const std::string& _outputFileName, const Args& ...args)
       : outputFileName(_outputFileName)
    {
        Initialize(args...);
        for(auto nametag : inputTags) {
            std::cout << nametag.first << std::endl;
            inputFiles.push_back(FilePtr(root_ext::OpenRootFile(nametag.first)));
        }
    }

    void Run()
    {
        std::shared_ptr<TCanvas> canvas(new TCanvas("","",50, 50, 700, 700));
        canvas->SetFillColor(0);
        canvas->SetBorderMode(0);
        canvas->SetFrameFillStyle(0);
        canvas->SetFrameLineColor(kWhite);
        canvas->SetFrameBorderMode(0);

        Int_t old_gErrorIgnoreLevel = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kWarning;
        canvas->Print((outputFileName + "(").c_str());
        gErrorIgnoreLevel = old_gErrorIgnoreLevel;

        gROOT->SetStyle("Plain");
        gStyle->SetOptStat(1111);
        gStyle->SetOptFit(111);


        canvas->cd();

        cms_tdr::setTDRStyle();

        std::shared_ptr<TPad> main_pad(root_ext::Adapter::NewPad(root_ext::Box(0.,0., 1, 1)));

        int W = 700;
        int H = 700;

        int H_ref = 700;
        int W_ref = 700;

        float T = 0.08*H_ref;
        float B = 0.14*H_ref;
        float L = 0.18*W_ref;
        float R = 0.05*W_ref;

        main_pad->SetLeftMargin( L/W );
        main_pad->SetRightMargin( R/W );
        main_pad->SetTopMargin( T/H );
        main_pad->SetBottomMargin( B/H );
        main_pad->SetTickx(0);
        main_pad->SetTicky(0);

        main_pad->Draw();
        main_pad->cd();


//        TH1D* Data_distr = root_ext::ReadCloneObject<TH1D>(*inputFiles.at(0), "pileup");
//        Data_distr->Scale( 1.0/ Data_distr->Integral() );
//        TH1D* MC_distr = root_ext::ReadCloneObject<TH1D>(*inputFiles.at(1), "pileup");
//        MC_distr->Scale( 1.0/ MC_distr->Integral() );

        TH1D* w_distr = root_ext::ReadCloneObject<TH1D>(*inputFiles.at(0), "weights");
//        Data_distr->Scale( 1.0/ Data_distr->Integral() );

        std::vector<double> x1, y1;
//        std::vector<double> x2, y2;
        for(int n = 1; n <= w_distr->GetNbinsX(); ++n){
            x1.push_back(w_distr->GetBinCenter(n));
            y1.push_back(w_distr->GetBinContent(n));
        }
//        for(int n = 1; n <= MC_distr->GetNbinsX(); ++n){
//            x2.push_back(MC_distr->GetBinCenter(n));
//            y2.push_back(MC_distr->GetBinContent(n));
//        }

        TGraph* gr = new TGraph(x1.size(), x1.data(), y1.data());
        gr->SetLineColor(kBlue); // kRed
        gr->SetMarkerColor(kBlue);
        gr->SetLineWidth(2);
        //gr->SetTitle("Option ACP example");
        gr->GetXaxis()->SetTitle("Number of interactions per bunch crossing");
        gr->GetYaxis()->SetTitle("Weight");
        //main_pad->SetLogy();
        gr->Draw("AC");

        gr->GetYaxis()->SetTitleOffset(1.45); //1.45
        gr->GetYaxis()->SetLabelSize(0.04);
        gr->GetYaxis()->SetTitleSize(0.04);
        gr->GetXaxis()->SetTitleOffset(1.05); //1.05
        gr->GetXaxis()->SetTitleSize(0.04);
        gr->GetXaxis()->SetLabelSize(0.04);
        gr->GetXaxis()->SetLabelOffset(0.015);
        gr->GetXaxis()->SetRangeUser(0, 60);

//        TGraph* gr2 = new TGraph(x2.size(), x2.data(), y2.data());
//        gr2->SetLineColor(kBlue);
//        gr2->SetMarkerColor(kBlue);
//        gr2->SetLineWidth(2);
        //gr2->Draw("same");

        std::shared_ptr<TLegend> legend(new TLegend (0.62, 0.55, 0.79, 0.65));
        legend->SetFillColor(0);
        legend->SetTextSize(0.035);
        legend->SetTextFont(42);
        legend->SetFillStyle (0);
        legend->SetFillColor (0);
        legend->SetBorderSize(0);
        legend->AddEntry(gr, "MC PU weight", "l");
        //legend->AddEntry(gr2, "Simulations", "l");
        legend->Draw("same");

        std::shared_ptr<TPaveText> text(new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC"));
        text->SetTextSize(0.05);
        text->SetTextFont(62);
        text->SetFillColor(0);
        text->SetBorderSize(0);
        text->SetMargin(0.01);
        text->SetTextAlign(12); // align left
        std::ostringstream ss_text;
        //ss_text << "e#tau_{h}";
        text->AddText(0.01,0.05, ss_text.str().c_str());
        text->Draw("same");

        cms_tdr::writeExtraText = true;
        cms_tdr::extraText = TString("Unpublished");
        cms_tdr::CMS_lumi(main_pad.get(), 2, 11);


        canvas->Draw();
        std::ostringstream print_options;
        print_options << "Title:" << "trala";
        std::cout << print_options.str() << std::endl;
        canvas->Print(outputFileName.c_str(), print_options.str().c_str());

        print_options << "Title:" << "title";
        canvas->Print(outputFileName.c_str(), print_options.str().c_str());

//        canvas->Clear();
//        canvas->Print(outputFileName.c_str());


        old_gErrorIgnoreLevel = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kWarning;
        canvas->Print((outputFileName+")").c_str());
        gErrorIgnoreLevel = old_gErrorIgnoreLevel;
    }

private:
    template<typename ...Args>
    void Initialize(const std::string& inputName, const Args& ...args)
    {
        const size_t split_index = inputName.find_first_of(':');
        const std::string fileName = inputName.substr(0, split_index);
        const std::string tagName = inputName.substr(split_index + 1);
        inputTags.push_back(FileNameTagPair(fileName, tagName));
        Initialize(args...);
    }

    void Initialize() {}

private:
    FileNameTagVector inputTags;
    FileVector inputFiles;
    std::string outputFileName;
};


