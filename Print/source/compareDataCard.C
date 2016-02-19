//#include <iosstream>
//#include <iomanip>

void drawHisto(TCanvas* C,TString Group1,TFile* F1,TString Group2,TFile* F2,TString channel,TString category,TString sample,bool PrintDiff,int MaxBin){
  C->Clear();
  
  TH1F* H1 = (TH1F*) F1->Get(channel+"_"+category+"/"+sample);
  TH1F* H2 = (TH1F*) F2->Get(channel+"_"+category+"/"+sample);
  
  if(!H1||!H2||H1->Integral()!=H1->Integral()||H2->Integral()!=H2->Integral()){

    if(H1&&H1->Integral()==H1->Integral())H1->Draw("histpe");
    if(H2&&H2->Integral()==H2->Integral())H2->Draw("histpe");

    TText cat;
    cat.DrawTextNDC(.2,.5,channel+"_"+category+"/"+sample);
    if(!H1||H1->Integral()!=H1->Integral()) cat.DrawTextNDC(.2,.4,Group1+" is missing");
    if(!H2||H2->Integral()!=H2->Integral()) cat.DrawTextNDC(.2,.4,Group2+" is missing");

    if(PrintDiff)printf("| %.1f ",100);
    
    C->Print(TString(C->GetName())+".pdf");

    return;
  }
//   cout<<channel<<" "<<category<<" "<<sample<<endl;
//   cout<<H1->Integral()<<" "<<H2->Integral()<<endl;
//   TText cat;
//   cat.DrawTextNDC(.2,.5,channel+"_"+category+"/"+sample);
//   C->Print(TString(C->GetName())+".pdf");
//   return;

  //H1->SetTitle("");
  
  C->cd() ;
  TPad histoPad("histoPad","histoPad",0.,0.3,1.,1.,0,0) ;
  histoPad.Draw() ;
  histoPad.cd() ;

  H1->SetTitle(channel+"  "+category+"  "+sample);
  H1->GetXaxis()->SetTitle("mass");
  H1->GetYaxis()->SetTitle("");
  H1->GetYaxis()->SetRangeUser(0,(H1->GetMaximum()>H2->GetMaximum() ? H1->GetMaximum() : H2->GetMaximum())*1.2);
  H1->SetMarkerStyle(8);
  H1->SetMarkerSize(0.6);
  H1->SetLineWidth(2);
  H1->SetMarkerColor(1);
  H1->SetLineColor(1);
  H1->SetLineStyle(1);
  H1->SetFillColor(0);
  if (MaxBin>0) H1->GetXaxis()->SetRange(0,MaxBin);
  H1->Draw("hist");
  TH1F * H4 = H1->Clone() ;  
  H4->SetLineWidth(0);
  H4->SetFillColor(kGray);
  H4->SetFillStyle(3344);
  H4->Draw("e2same");
  H2->SetMarkerStyle(8);
  H2->SetMarkerSize(0.6);
  H2->SetMarkerColor(4);
  H2->SetLineWidth(2);
  H2->SetLineStyle(1);
  H2->SetLineColor(4);
  H2->SetFillColor(0);
  if (MaxBin>0) H2->GetXaxis()->SetRange(0,MaxBin);
  H2->Draw("hist e same");
  
  TText TXmine;
  TXmine.SetTextColor(1);
  TXmine.SetTextSize(.04);
  TText TXother;
  TXother.SetTextColor(H2->GetLineColor());
  TXother.SetTextSize(.04);
  TText TXdiff;
  TXdiff.SetTextColor(2);
  TXdiff.SetTextSize(.03);
  char yield1[100];
  const double H1_Integral = H1->Integral(1, H1->GetNbinsX());
  const double H2_Integral = H2->Integral(1, H2->GetNbinsX());
  sprintf(yield1,"%.2f",H1_Integral);
  char yield2[100];
  sprintf(yield2,"%.2f",H2_Integral);
  TXmine.DrawTextNDC(.55,.845,Group1+" : "+yield1);
  TXother.DrawTextNDC(.55,.81,Group2+" : "+yield2);
  char txt[100];
  //float diff=100*2*fabs(H1->Integral(1,H1->GetNbinsX())-H2->Integral(1,H2->GetNbinsX()))/(H1->Integral(1,H1->GetNbinsX())+H2->Integral(1,H2->GetNbinsX()));
  float diff=100*2*(H1_Integral-H2_Integral)/(H1_Integral+H2_Integral);
  //sprintf(txt,TString(Group1+"-"+Group2+" difference = %.1f",diff));
  sprintf(txt,"difference = %.1f",diff);
  TXdiff.DrawTextNDC(.25,.85,TString(txt)+"%");

  C->cd() ;
  TPad ratioPad("ratioPad","ratioPad",0.,0.3,1,.1,0,0) ;
  ratioPad.Draw() ;
  ratioPad.cd() ;
  ratioPad.SetGridy(1) ;
  
  TH1F * H3 = H1->Clone() ;  
  H3->Divide(H2) ;
  H3->SetMaximum(1.3) ;
  H3->SetMinimum(0.7) ;
  H3->SetTitle("") ;
  H3->GetXaxis()->SetTitle("") ;
  H3->GetXaxis()->SetLabelSize(0.08) ;
  H3->GetYaxis()->SetLabelSize(0.08) ;
  H3->GetYaxis()->SetTitleOffset(0.3) ;
  H3->GetYaxis()->SetTitleSize(0.12) ;
  H3->GetYaxis()->SetTitle(TString(Group1)+" / "+TString(Group2)) ;
  if (MaxBin>0) H3->GetXaxis()->SetRange(0,MaxBin);
  H3->Draw() ;

  //cout<<setiosflags(ios::fixed)<<precision(2);
  //cout<<"| "<<diff<<" "<<endl;
  if(PrintDiff)printf("| %.1f ",diff);
  C->Print(TString(C->GetName())+".pdf");
}


void drawCategory(TCanvas* C,TString Group1, TFile* F1,TString Group2, TFile* F2,TString channel, TString category,
                  bool PrintDiff,int MaxBin) {

    drawHisto(C,Group1,F1,Group2,F2,channel,category,"data_obs"                             ,1&&PrintDiff,MaxBin);
    drawHisto(C,Group1,F1,Group2,F2,channel,category,"ggHTohhTo2Tau2B300"                   ,1&&PrintDiff,MaxBin);
    drawHisto(C,Group1,F1,Group2,F2,channel,category,"ZTT"                                  ,1&&PrintDiff,MaxBin);
    if(channel != "tauTau") {
        drawHisto(C,Group1,F1,Group2,F2,channel,category,"ZL"                                   ,1&&PrintDiff,MaxBin);
        drawHisto(C,Group1,F1,Group2,F2,channel,category,"ZJ"                                   ,1&&PrintDiff,MaxBin);
    }
    drawHisto(C,Group1,F1,Group2,F2,channel,category,"ZLL"                                  ,1&&PrintDiff,MaxBin);
    if(channel != "tauTau" || category == "2jet0tag")
        drawHisto(C,Group1,F1,Group2,F2,channel,category,"W"                                    ,1&&PrintDiff,MaxBin);
    drawHisto(C,Group1,F1,Group2,F2,channel,category,"QCD"                                  ,1&&PrintDiff,MaxBin);
    drawHisto(C,Group1,F1,Group2,F2,channel,category,"TT"                                   ,1&&PrintDiff,MaxBin);
    drawHisto(C,Group1,F1,Group2,F2,channel,category,"VV"                                   ,1&&PrintDiff,MaxBin);
//  drawHisto(C,Group1,F1,Group2,F2,channel,category,"TT_emb"                                   ,1&&PrintDiff,MaxBin);
}


void compareDataCard(
		     TString channel = "tauTau",
		     TString year    = "2012",
		     TString Group1  = "CERN",
		     TString Path1   = "/afs/cern.ch/user/b/benitezj/public/datacards/2011/Aug4",
		     TString File1   = "muTauSM_svfitmass",
		     TString Group2  = "MIT",
		     TString Path2   = "/afs/cern.ch/user/b/bianchi/public/Roger/datacards2011v3",
		     TString File2   = "muTauSM", 
		     int MaxBin      = 0  
		     ){

  TFile F1(Path1+"/"+File1+".root","read");
  TFile F2(Path2+"/"+File2+".root","read");
  
  F1.ls();
  F2.ls();

  TString fname=channel+year+"_"+Group1+"_"+Group2+"_DataCardDiff";
  TCanvas C(fname,fname,700,700);
  
  
  
  C.Print(TString(C.GetName())+".pdf[");
  
  //cout<<" | 0jet_low Data | 0jet_low ZTT | 0jet_low W | 0jet_low QCD | 0jet_low TT | 0jet_low  ggH | 0jet_low  qqH | Boosted_high Data | Boosted_high ZTT | Boosted_high W | Boosted_high QCD | Boosted_high TT | Boosted_high ggH | Boosted_high qqH | VBF Data| VBF ZTT | VBF W | VBF QCD | VBF TT | VBF qqH | VBF ggH |"<<endl;

  //cout<<"| "<<Group1<<":"<<Group2<<TString(" | [[http://benitezj.web.cern.ch/benitezj/HTTSync/DataCards201252X/")<<fname<<".pdf][pdf]] ";

//   drawCategory(&C,Group1,&F1,Group2,&F2,channel,"inclusive",1,MaxBin       );
//   drawCategory(&C,Group1,&F1,Group2,&F2,channel,"nobtag"   ,1,MaxBin       );
//   drawCategory(&C,Group1,&F1,Group2,&F2,channel,"btag"     ,1,(int)MaxBin/2);

//  drawCategory(&C,Group1,&F1,Group2,&F2,channel,"inclusive"             ,1,MaxBin       );
  drawCategory(&C,Group1,&F1,Group2,&F2,channel,"2jet0tag" ,1,MaxBin       );
  drawCategory(&C,Group1,&F1,Group2,&F2,channel,"2jet1tag"   ,1,MaxBin       );
  drawCategory(&C,Group1,&F1,Group2,&F2,channel,"2jet2tag"                   ,1,(int)MaxBin/2);

  C.Print(TString(C.GetName())+".pdf]");
  
  cout<<"|"<<endl;

  gROOT->ProcessLine(".q");
}

