/*! Print TGraph.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <TTree.h>
#include <TGraph.h>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Print/include/RootPrintToPdf.h"
#include "Analysis/include/Htautau_TriggerEfficiency.h"

class Print_Graph_2 {
public:
    typedef std::pair< std::string, std::string > FileNameTagPair;
    typedef std::vector< FileNameTagPair> FileNameTagVector;
    typedef std::shared_ptr<TFile> FilePtr;
    typedef std::vector< FilePtr > FileVector;

    template<typename ...Args>
    Print_Graph_2(const std::string& _outputFileName, const Args& ...args)
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


        TH1D* distr1 = root_ext::ReadCloneObject<TH1D>(*inputFiles.at(0),
                                                       "2jets1btag/NoCuts/OS_Isolated/Central/EMBEDDED MuTau/m_sv");
//        distr1->Scale(1.0 / distr1->Integral());
        distr1->SetLineColor(kBlue);
        distr1->SetMarkerColor(kBlue);
        distr1->SetLineWidth(2);
//        distr1->GetXaxis()->SetTitle("M_{#tau#tau} (GeV), SVfit");
        //distr1->GetYaxis()->SetTitle("dN/(N dm_{#tau#tau}) (1/GeV)");
        //main_pad->SetLogy();
        distr1->SetMaximum(distr1->GetMaximum()*1.3);
        distr1->Draw();


        distr1->GetYaxis()->SetTitleOffset(1.55); //1.45
        distr1->GetYaxis()->SetLabelSize(0.04);
        distr1->GetYaxis()->SetTitleSize(0.04);
        distr1->GetXaxis()->SetTitleOffset(1.05); //1.05
        distr1->GetXaxis()->SetTitleSize(0.04);
        distr1->GetXaxis()->SetLabelSize(0.04);
        distr1->GetXaxis()->SetLabelOffset(0.015);
//        distr1->GetXaxis()->SetRangeUser(0, 60);

        TH1D* distr2 = root_ext::ReadCloneObject<TH1D>(*inputFiles.at(0),
                                                       "2jets1btag/NoCuts/OS_Isolated/Central/TT_EMBEDDED MuTau/m_sv");
//        distr2->Scale(1.0 / distr2->Integral());
        distr2->SetLineColor(kRed);
        distr2->SetMarkerColor(kRed);
        distr2->SetLineWidth(2);
        distr2->Draw("same");

        std::shared_ptr<TLegend> legend(new TLegend (0.55, 0.55, 0.75, 0.65));
        legend->SetFillColor(0);
        legend->SetTextSize(0.035);
        legend->SetTextFont(42);
        legend->SetFillStyle (0);
        legend->SetFillColor (0);
        legend->SetBorderSize(0);
//        legend->AddEntry(distr1, "h#rightarrow#tau#tau in H#rightarrowhh(m_{H}=300)", "l");
//        legend->AddEntry(distr2, "Z#rightarrow#tau#tau", "l");
        legend->AddEntry(distr1, "Z#rightarrow#tau#tau embedded", "l");
        legend->AddEntry(distr2, "t#bar{t} embedded", "l");
        legend->Draw("same");

        std::shared_ptr<TPaveText> text(new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC"));
        text->SetTextSize(0.04);
        text->SetTextFont(52);
        text->SetFillColor(0);
        text->SetBorderSize(0);
        text->SetMargin(0.01);
        text->SetTextAlign(12); // align left
        std::ostringstream ss_text;
        ss_text << "#mu#tau_{h}, 2jets-1btag";
        text->AddText(0.01,0.05, ss_text.str().c_str());
        text->Draw("same");

        cms_tdr::writeExtraText = true;
        cms_tdr::extraText = TString("Unpublished");
        cms_tdr::CMS_lumi(main_pad.get(), 2, 11);


        canvas->Draw();
        std::ostringstream print_options;
        print_options << "Title:" << "m_h";
        std::cout << print_options.str() << std::endl;
        canvas->Print(outputFileName.c_str(), print_options.str().c_str());

        old_gErrorIgnoreLevel = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kWarning;
        canvas->Print((outputFileName+"]").c_str());
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


