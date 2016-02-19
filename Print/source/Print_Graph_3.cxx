/*! Print TGraph.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <TTree.h>
#include <TGraph.h>

#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Print/include/RootPrintToPdf.h"
#include "Analysis/include/Htautau_TriggerEfficiency.h"

class Print_Graph_3 {
public:
    typedef std::pair< std::string, std::string > FileNameTagPair;
    typedef std::vector< FileNameTagPair> FileNameTagVector;
    typedef std::shared_ptr<TFile> FilePtr;
    typedef std::vector< FilePtr > FileVector;

    template<typename ...Args>
    Print_Graph_3(const std::string& _outputFileName, const Args& ...args)
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

        const std::string sep = "/";
        const std::string prefix = "2jets2btag/KinFitConvergedWithMassWindow/OS_Isolated";
//        const std::vector<std::string> scales = { "Central", "TauUp", "TauDown" };
        const std::vector<std::string> scales = { "Central", "JetUp", "JetDown" };
//        const std::string sample_name = "ggHhh300";
        const std::string sample_name = "TTbar";
        const std::string distr_name = "m_ttbb_kinfit";
//        const std::vector<std::string> legend_names = { "central #tau energy scale", "#tau energy scaled up by 1#sigma",
//                                                        "#tau energy scaled down by 1#sigma" };
        const std::vector<std::string> legend_names = { "central jets energy scale", "jets energy scaled up by 1#sigma",
                                                        "jets energy scaled down by 1#sigma" };
        std::vector<std::ostringstream> ss_paths(scales.size());

        for(size_t n = 0; n < ss_paths.size(); ++n)
            ss_paths.at(n) << prefix << sep << scales.at(n) << sep << sample_name << sep << distr_name;

        std::vector<TH1D*> distrs(ss_paths.size());

        distrs.at(0) = root_ext::ReadCloneObject<TH1D>(*inputFiles.at(0), ss_paths.at(0).str());
//        distr1->Scale(1.0 / distr1->Integral());
        distrs.at(0)->SetLineColor(kBlack);
        distrs.at(0)->SetMarkerColor(kBlack);
        distrs.at(0)->SetLineWidth(2);
//        distr1->GetXaxis()->SetTitle("M_{#tau#tau} (GeV), SVfit");
        //distr1->GetYaxis()->SetTitle("dN/(N dm_{#tau#tau}) (1/GeV)");
        //main_pad->SetLogy();
        distrs.at(0)->SetMaximum(distrs.at(0)->GetMaximum()*1.4);
        distrs.at(0)->Draw();

        distrs.at(0)->GetYaxis()->SetTitleOffset(1.55); //1.45
        distrs.at(0)->GetYaxis()->SetLabelSize(0.04);
        distrs.at(0)->GetYaxis()->SetTitleSize(0.04);
        distrs.at(0)->GetXaxis()->SetTitleOffset(1.2); //1.05
        distrs.at(0)->GetXaxis()->SetTitleSize(0.04);
        distrs.at(0)->GetXaxis()->SetLabelSize(0.04);
        distrs.at(0)->GetXaxis()->SetLabelOffset(0.015);
//        distr1->GetXaxis()->SetRangeUser(0, 60);

        distrs.at(1) = root_ext::ReadCloneObject<TH1D>(*inputFiles.at(0), ss_paths.at(1).str());
//        distr2->Scale(1.0 / distr2->Integral());
        distrs.at(1)->SetLineColor(kBlue);
        distrs.at(1)->SetMarkerColor(kBlue);
        distrs.at(1)->SetLineWidth(2);
        distrs.at(1)->Draw("same");

        distrs.at(2) = root_ext::ReadCloneObject<TH1D>(*inputFiles.at(0), ss_paths.at(2).str());
//        distr2->Scale(1.0 / distr2->Integral());
        distrs.at(2)->SetLineColor(kRed);
        distrs.at(2)->SetMarkerColor(kRed);
        distrs.at(2)->SetLineWidth(2);
        distrs.at(2)->Draw("same");

        std::shared_ptr<TLegend> legend(new TLegend (0.45, 0.55, 0.75, 0.65));
        legend->SetFillColor(0);
        legend->SetTextSize(0.035);
        legend->SetTextFont(42);
        legend->SetFillStyle (0);
        legend->SetFillColor (0);
        legend->SetBorderSize(0);
//        legend->AddEntry(distr1, "h#rightarrow#tau#tau in H#rightarrowhh(m_{H}=300)", "l");
//        legend->AddEntry(distr2, "Z#rightarrow#tau#tau", "l");
        for(size_t n = 0; n < distrs.size(); ++n)
            legend->AddEntry(distrs.at(n), legend_names.at(n).c_str(), "l");
        legend->Draw("same");

        std::shared_ptr<TPaveText> text(new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC"));
        text->SetTextSize(0.04);
        text->SetTextFont(52);
        text->SetFillColor(0);
        text->SetBorderSize(0);
        text->SetMargin(0.01);
        text->SetTextAlign(12); // align left
        std::ostringstream ss_text;
        ss_text << "#mu#tau_{h}, 2jets-2btag";
        text->AddText(0.01,0.05, ss_text.str().c_str());
        text->Draw("same");

        cms_tdr::writeExtraText = true;
        cms_tdr::extraText = TString("Unpublished");
        cms_tdr::CMS_lumi(main_pad.get(), 2, 11);


        canvas->Draw();
        std::ostringstream print_options;
        print_options << "Title:" << "title";
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



