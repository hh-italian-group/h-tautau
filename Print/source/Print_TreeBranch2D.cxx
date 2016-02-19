/*! Print histogram for a tree branch with a specified name superimposing several files.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <TTree.h>
#include <TGraph.h>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Print/include/RootPrintToPdf.h"
#include "Analysis/include/Htautau_TriggerEfficiency.h"

class Print_TreeBranch2D {
public:

    Print_TreeBranch2D(const std::string& _outputFileName, const std::string& _inputFileName)
       : outputFileName(_outputFileName), inputFile(root_ext::OpenRootFile(_inputFileName))
    {
    }

    void Run()
    {
        std::shared_ptr<TCanvas> canvas(new TCanvas("","",50, 50, 700, 500));
        canvas->SetFillColor(0);
        canvas->SetBorderMode(0);
        canvas->SetFrameFillStyle(0);
        canvas->SetFrameLineColor(kWhite);
        canvas->SetFrameBorderMode(0);

        Int_t old_gErrorIgnoreLevel = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kWarning;
        canvas->Print((outputFileName + "[").c_str());
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

        using namespace analysis::Htautau_Summer13::trigger::Run2012ABCD::ETau;

        const double eta = 0.;
        const double pt_0 = 0, pt_step = 0.1;
        const Int_t n = 600;
        Double_t x[n], y[n], y2[n];
        for (Int_t i=0;i<n;i++) {
            const double pt = pt_0 + i*pt_step;
            x[i] = pt;
            y[i] = Data::tauEfficiency(pt, eta);
            y2[i] = MC::tauEfficiency(pt, eta);
        }

        TGraph* gr = new TGraph(n,x,y);
        gr->SetLineColor(kRed);
        gr->SetMarkerColor(kRed);
        gr->SetLineWidth(2);
        gr->SetTitle("Option ACP example");
        gr->GetXaxis()->SetTitle("P_{T}(#tau) [GeV]");
        gr->GetYaxis()->SetTitle("Efficiency");
        gr->Draw("AC");

        gr->GetYaxis()->SetTitleOffset(0.7); //1.45
        gr->GetYaxis()->SetLabelSize(0.04);
        gr->GetXaxis()->SetTitleOffset(1.00); //1.05
        gr->GetXaxis()->SetLabelSize(0.04);
        gr->GetXaxis()->SetLabelOffset(0.015);
        gr->GetXaxis()->SetRangeUser(0, 60);

        TGraph* gr2 = new TGraph(n,x,y2);
        gr2->SetLineColor(kBlue);
        gr2->SetMarkerColor(kBlue);
        gr2->SetLineWidth(2);
        gr2->Draw("same");

        std::shared_ptr<TLegend> legend(new TLegend (0.72, 0.55, 0.89, 0.65));
        legend->SetFillColor(0);
        legend->SetTextSize(0.04);
        legend->SetTextFont(42);
        legend->SetFillStyle (0);
        legend->SetFillColor (0);
        legend->SetBorderSize(0);
        legend->AddEntry(gr, "Data", "l");
        legend->AddEntry(gr2, "Simulations", "l");
        legend->Draw("same");

        std::shared_ptr<TPaveText> text(new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC"));
        text->SetTextSize(0.05);
        text->SetTextFont(62);
        text->SetFillColor(0);
        text->SetBorderSize(0);
        text->SetMargin(0.01);
        text->SetTextAlign(12); // align left
        std::ostringstream ss_text;
        ss_text << "e#tau_{h}";
        text->AddText(0.01,0.05, ss_text.str().c_str());
        text->Draw("same");

        cms_tdr::writeExtraText = true;
        cms_tdr::extraText = TString("Unpublished");
        cms_tdr::CMS_lumi(main_pad.get(), 12, 11);


        canvas->Draw();
        std::ostringstream print_options;
        print_options << "Title:" << "test";
        canvas->Print(outputFileName.c_str(), print_options.str().c_str());


        old_gErrorIgnoreLevel = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kWarning;
        canvas->Print((outputFileName+"]").c_str());
        gErrorIgnoreLevel = old_gErrorIgnoreLevel;
    }

private:
    std::string outputFileName;
    std::shared_ptr<TFile> inputFile;
};

