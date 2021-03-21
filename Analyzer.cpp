#define Analyzer_cxx
#include "Analyzer.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include<TH1F.h>
#include<TH2F.h>
#include<TF1.h>
#include <vector>
#include<TLorentzVector.h>
#include<TLegend.h>
#include<THStack.h>
#include<TGraph.h>
#include<TGraphPainter.h>
#include<TMath.h>
#include "TSystem.h"
#include "cConstants.h"
#include <TFile.h>
#include <cassert>
#include <cmath>
#include <iostream>
cConstantSpline::cConstantSpline(const TString& filename) : filename_(filename), spline_(nullptr) {}

void cConstantSpline::initspline(bool isDbkg) {
  if (!spline_) {
    TFile* f_ =TFile::Open(filename_);
    spline_.reset((TSpline3*)(f_->Get("sp_gr_varReco_Constant_Smooth")->Clone()));
    f_->Close();
  }
  assert(spline_.get());
}

double cConstantSpline::eval(double ZZMass, bool isDbkg) {
  initspline(isDbkg);
  return spline_->Eval(ZZMass);
}

namespace {
 
  cConstantSpline DVBF2jetsSpline("SmoothKDConstant_m4l_DjjVBF13TeV.root");
  cConstantSpline DVBF1jetSpline("SmoothKDConstant_m4l_DjVBF13TeV.root");
  cConstantSpline DZHhSpline("SmoothKDConstant_m4l_DjjZH13TeV.root");
  cConstantSpline DWHhSpline("SmoothKDConstant_m4l_DjjWH13TeV.root");
}

extern "C" float getDVBF2jetsConstant(float ZZMass){
  return DVBF2jetsSpline.eval(ZZMass, false);
}

extern "C" float getDVBF1jetConstant(float ZZMass){
  return DVBF1jetSpline.eval(ZZMass, false);
}

extern "C" float getDWHhConstant(float ZZMass){
  return DWHhSpline.eval(ZZMass, false);
}
extern "C" float getDZHhConstant(float ZZMass){
  return DZHhSpline.eval(ZZMass, false);
}

void Analyzer::Loop()
{
 
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      cout<<LepPt->at(0);
      // if (Cut(ientry) < 0) continue;
   }
}
void Analyzer::Fill_Histogram(TString s)
{
  TFile *f;  
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.  
      f = new TFile(s);
      TDirectory * dir = (TDirectory*)f->Get(s+":/ZZTree");
      dir->GetObject("candTree",tree);

   double w;
   Init(tree);
   h1 = (TH1F*)f->Get("ZZTree/Counters"); //dodano (radi citanja iz Countersa)
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
	  w=137000*xsec*overallEventWeight/h1->GetBinContent(40);
      if(s.Contains("ggH125"))
	  {
		  Histo_PFMET_ggH125->Fill(PFMET,w);
          Histo_nCleanedJetsPt30_ggH125->Fill(nCleanedJetsPt30,w);
          Histo_nCleanedJetsPt30BTagged_ggH125->Fill(nCleanedJetsPt30BTagged,w);
          Histo_ZZPt_ggH125->Fill(ZZPt,w);
          Histo_Z1Pt_ggH125->Fill(Z1Pt,w);
          Histo_Z2Pt_ggH125->Fill(Z2Pt,w);
          Histo_PhotonPt_ggH125->Fill(*max_element(PhotonPt->begin(),PhotonPt->end()),w); 
          Histo_nExtraLep_ggH125->Fill(nExtraLep,w);
          Histo_nExtraZ_ggH125->Fill(nExtraZ,w);
		  if(DiJetMass>0)
          Histo_DiJetMass_ggH125->Fill(DiJetMass,w);	  
	  }	  
      if(s.Contains("VBFH125"))
	  {
		  Histo_PFMET_VBFH125->Fill(PFMET,w);
          Histo_nCleanedJetsPt30_VBFH125->Fill(nCleanedJetsPt30,w);	
          Histo_nCleanedJetsPt30BTagged_VBFH125->Fill(nCleanedJetsPt30BTagged,w);
          Histo_ZZPt_VBFH125->Fill(ZZPt,w);
          Histo_Z1Pt_VBFH125->Fill(Z1Pt,w);
          Histo_Z2Pt_VBFH125->Fill(Z2Pt,w);	
          Histo_PhotonPt_VBFH125->Fill(*max_element(PhotonPt->begin(),PhotonPt->end()),w);
          Histo_nExtraLep_VBFH125->Fill(nExtraLep,w);
          Histo_nExtraZ_VBFH125->Fill(nExtraZ,w);
		  if(DiJetMass>0)
          Histo_DiJetMass_VBFH125->Fill(DiJetMass,w);		  
      }
	  if(s.Contains("ttH125"))
	  {
		  Histo_PFMET_ttH125->Fill(PFMET,w);
          Histo_nCleanedJetsPt30_ttH125->Fill(nCleanedJetsPt30,w);
          Histo_nCleanedJetsPt30BTagged_ttH125->Fill(nCleanedJetsPt30BTagged,w);
          Histo_ZZPt_ttH125->Fill(ZZPt,w);
          Histo_Z1Pt_ttH125->Fill(Z1Pt,w);
          Histo_Z2Pt_ttH125->Fill(Z2Pt,w);
          Histo_PhotonPt_ttH125->Fill(*max_element(PhotonPt->begin(),PhotonPt->end()),w);
          Histo_nExtraLep_ttH125->Fill(nExtraLep,w);
          Histo_nExtraZ_ttH125->Fill(nExtraZ,w);
		  if(DiJetMass>0)
          Histo_DiJetMass_ttH125->Fill(DiJetMass,w);		  
      }
	  if(s.Contains("ZZTo4lext1"))
	  {
		  Histo_PFMET_ZZTo4lext1->Fill(PFMET,w);
          Histo_nCleanedJetsPt30_ZZTo4lext1->Fill(nCleanedJetsPt30,w);
          Histo_nCleanedJetsPt30BTagged_ZZTo4lext1->Fill(nCleanedJetsPt30BTagged,w);
          Histo_ZZPt_ZZTo4lext1->Fill(ZZPt,w);
          Histo_Z1Pt_ZZTo4lext1->Fill(Z1Pt,w);
          Histo_Z2Pt_ZZTo4lext1->Fill(Z2Pt,w);
          Histo_PhotonPt_ZZTo4lext1->Fill(*max_element(PhotonPt->begin(),PhotonPt->end()),w);
          Histo_nExtraLep_ZZTo4lext1->Fill(nExtraLep,w);
          Histo_nExtraZ_ZZTo4lext1->Fill(nExtraZ,w);
		  if(DiJetMass>0)
          Histo_DiJetMass_ZZTo4lext1->Fill(DiJetMass,w);		  
      }
	  }
	   
}
void Analyzer :: Plot_Histogram()
{
	//PFMET
	gStyle->SetOptStat(0);
	TCanvas *PFMET_Canvas=new TCanvas("PFMET","PFMET",1600,800);
	TLegend* PFMET_Legend=new TLegend(0.75,0.75,0.9,0.9);
	Histo_PFMET_ggH125->Scale(1/Histo_PFMET_ggH125->Integral());
	Histo_PFMET_VBFH125->Scale(1/Histo_PFMET_VBFH125->Integral());
	Histo_PFMET_ttH125->Scale(1/Histo_PFMET_ttH125->Integral());
	Histo_PFMET_ZZTo4lext1->Scale(1/Histo_PFMET_ZZTo4lext1->Integral());
	PFMET_Legend->AddEntry(Histo_PFMET_ggH125,"ggH125","f");
	PFMET_Legend->AddEntry(Histo_PFMET_VBFH125,"VBFH125","f");
	PFMET_Legend->AddEntry(Histo_PFMET_ttH125,"ttH125","f");
	PFMET_Legend->AddEntry(Histo_PFMET_ZZTo4lext1,"qqZZ","f");
	Histo_PFMET_ggH125->SetLineColor(kBlue);
	Histo_PFMET_VBFH125->SetLineColor(kGreen);
	Histo_PFMET_ttH125->SetLineColor(kMagenta);
	Histo_PFMET_ZZTo4lext1->SetLineColor(kRed+2);
	Histo_PFMET_ggH125->SetTitle("Missing transverse momentum");
	Histo_PFMET_ggH125->GetXaxis()->SetTitle("Pt_{miss}/GeV");
	Histo_PFMET_ggH125->GetYaxis()->SetTitle("normalized to 1");
	Histo_PFMET_ggH125->Draw("HIST SAME");
	Histo_PFMET_VBFH125->Draw("HIST SAME");
	Histo_PFMET_ttH125->Draw("HIST SAME");
	Histo_PFMET_ZZTo4lext1->Draw("HIST SAME");
	PFMET_Legend->Draw();
	PFMET_Canvas->SaveAs("PFMET.pdf");
	
	//nCleanedJetsPt30
	TCanvas * nCleanedJetsPt30_Canvas=new TCanvas("nCleanedJetsPt30","nCleanedJetsPt30",1600,800);
	nCleanedJetsPt30_Canvas->cd();
	TLegend*  nCleanedJetsPt30_Legend=new TLegend(0.75,0.75,0.9,0.9);
	Histo_nCleanedJetsPt30_ggH125->Scale(1/Histo_nCleanedJetsPt30_ggH125->Integral());
	Histo_nCleanedJetsPt30_VBFH125->Scale(1/Histo_nCleanedJetsPt30_VBFH125->Integral());
	Histo_nCleanedJetsPt30_ttH125->Scale(1/Histo_nCleanedJetsPt30_ttH125->Integral());
	Histo_nCleanedJetsPt30_ZZTo4lext1->Scale(1/Histo_nCleanedJetsPt30_ZZTo4lext1->Integral());
	nCleanedJetsPt30_Legend->AddEntry(Histo_nCleanedJetsPt30_ggH125,"ggH125","f");
	nCleanedJetsPt30_Legend->AddEntry(Histo_nCleanedJetsPt30_VBFH125,"VBFH125","f");
	nCleanedJetsPt30_Legend->AddEntry(Histo_nCleanedJetsPt30_ttH125,"ttH125","f");
	nCleanedJetsPt30_Legend->AddEntry(Histo_nCleanedJetsPt30_ZZTo4lext1,"qqZZ","f");
	Histo_nCleanedJetsPt30_ggH125->SetLineColor(kBlue);
	Histo_nCleanedJetsPt30_VBFH125->SetLineColor(kGreen);
	Histo_nCleanedJetsPt30_ttH125->SetLineColor(kMagenta);
	Histo_nCleanedJetsPt30_ZZTo4lext1->SetLineColor(kRed+2);
	Histo_nCleanedJetsPt30_ZZTo4lext1->SetTitle("Jet number");
	Histo_nCleanedJetsPt30_ZZTo4lext1->GetXaxis()->SetTitle("jets");
	Histo_nCleanedJetsPt30_ZZTo4lext1->GetYaxis()->SetTitle("normalized to 1");
	Histo_nCleanedJetsPt30_ZZTo4lext1->Draw("HIST SAME");
	Histo_nCleanedJetsPt30_ggH125->Draw("HIST SAME");
	Histo_nCleanedJetsPt30_VBFH125->Draw("HIST SAME");
	Histo_nCleanedJetsPt30_ttH125->Draw("HIST SAME");
	nCleanedJetsPt30_Legend->Draw();
	nCleanedJetsPt30_Canvas->SaveAs("nCleanedJetsPt30.pdf");
	
	//nCleanedJetsPt30BTagged
	TCanvas * nCleanedJetsPt30BTagged_Canvas=new TCanvas("nCleanedJetsPt30BTagged","nCleanedJetsPt30BTagged",1600,800);
	nCleanedJetsPt30BTagged_Canvas->cd();
	TLegend*  nCleanedJetsPt30BTagged_Legend=new TLegend(0.75,0.75,0.9,0.9);
	Histo_nCleanedJetsPt30BTagged_ggH125->Scale(1/Histo_nCleanedJetsPt30BTagged_ggH125->Integral());
	Histo_nCleanedJetsPt30BTagged_VBFH125->Scale(1/Histo_nCleanedJetsPt30BTagged_VBFH125->Integral());
	Histo_nCleanedJetsPt30BTagged_ttH125->Scale(1/Histo_nCleanedJetsPt30BTagged_ttH125->Integral());
	Histo_nCleanedJetsPt30BTagged_ZZTo4lext1->Scale(1/Histo_nCleanedJetsPt30BTagged_ZZTo4lext1->Integral());
	nCleanedJetsPt30BTagged_Legend->AddEntry(Histo_nCleanedJetsPt30BTagged_ggH125,"ggH125","f");
	nCleanedJetsPt30BTagged_Legend->AddEntry(Histo_nCleanedJetsPt30BTagged_VBFH125,"VBFH125","f");
	nCleanedJetsPt30BTagged_Legend->AddEntry(Histo_nCleanedJetsPt30BTagged_ttH125,"ttH125","f");
	nCleanedJetsPt30BTagged_Legend->AddEntry(Histo_nCleanedJetsPt30BTagged_ZZTo4lext1,"qqZZ","f");
	Histo_nCleanedJetsPt30BTagged_ggH125->SetLineColor(kBlue);
	Histo_nCleanedJetsPt30BTagged_VBFH125->SetLineColor(kGreen);
	Histo_nCleanedJetsPt30BTagged_ttH125->SetLineColor(kMagenta);
	Histo_nCleanedJetsPt30BTagged_ZZTo4lext1->SetLineColor(kRed+2);
	Histo_nCleanedJetsPt30BTagged_ZZTo4lext1->SetTitle("Number of b-tagged jets");
	Histo_nCleanedJetsPt30BTagged_ZZTo4lext1->GetXaxis()->SetTitle("jets");
	Histo_nCleanedJetsPt30BTagged_ZZTo4lext1->GetYaxis()->SetTitle("normalized to 1");
	Histo_nCleanedJetsPt30BTagged_ZZTo4lext1->Draw("HIST SAME");
	Histo_nCleanedJetsPt30BTagged_ggH125->Draw("HIST SAME");
	Histo_nCleanedJetsPt30BTagged_VBFH125->Draw("HIST SAME");
	Histo_nCleanedJetsPt30BTagged_ttH125->Draw("HIST SAME");
	nCleanedJetsPt30BTagged_Legend->Draw();
	nCleanedJetsPt30BTagged_Canvas->SaveAs("nCleanedJetsPt30BTagged.pdf");
	
	//ZZPt
	TCanvas *ZZPt_Canvas=new TCanvas("ZZPt","ZZPt",1600,800);
	ZZPt_Canvas->cd();
	TLegend* ZZPt_Legend=new TLegend(0.75,0.75,0.9,0.9);
	Histo_ZZPt_ggH125->Scale(1/Histo_ZZPt_ggH125->Integral());
	Histo_ZZPt_VBFH125->Scale(1/Histo_ZZPt_VBFH125->Integral());
	Histo_ZZPt_ttH125->Scale(1/Histo_ZZPt_ttH125->Integral());
	Histo_ZZPt_ZZTo4lext1->Scale(1/Histo_ZZPt_ZZTo4lext1->Integral());
	ZZPt_Legend->AddEntry(Histo_ZZPt_ggH125,"ggH125","f");
	ZZPt_Legend->AddEntry(Histo_ZZPt_VBFH125,"VBFH125","f");
	ZZPt_Legend->AddEntry(Histo_ZZPt_ttH125,"ttH125","f");
	ZZPt_Legend->AddEntry(Histo_ZZPt_ZZTo4lext1,"qqZZ","f");
	Histo_ZZPt_ggH125->SetLineColor(kBlue);
	Histo_ZZPt_VBFH125->SetLineColor(kGreen);
	Histo_ZZPt_ttH125->SetLineColor(kMagenta);
	Histo_ZZPt_ZZTo4lext1->SetLineColor(kRed+2);
	Histo_ZZPt_ZZTo4lext1->SetTitle("Transverse momentum of 4 leptons");
	Histo_ZZPt_ZZTo4lext1->GetXaxis()->SetTitle("Pt_{ZZ}/GeV");
	Histo_ZZPt_ZZTo4lext1->GetYaxis()->SetTitle("normalized to 1");
	Histo_ZZPt_ZZTo4lext1->Draw("HIST SAME");
	Histo_ZZPt_ggH125->Draw("HIST SAME");
	Histo_ZZPt_VBFH125->Draw("HIST SAME");
	Histo_ZZPt_ttH125->Draw("HIST SAME");
	ZZPt_Legend->Draw();
	ZZPt_Canvas->SaveAs("ZZPt.pdf");

    //Z1Pt
	TCanvas *Z1Pt_Canvas=new TCanvas("Z1Pt","Z1Pt",1600,800);
	Z1Pt_Canvas->cd();
	TLegend* Z1Pt_Legend=new TLegend(0.75,0.75,0.9,0.9);
	Histo_Z1Pt_ggH125->Scale(1/Histo_Z1Pt_ggH125->Integral());
	Histo_Z1Pt_VBFH125->Scale(1/Histo_Z1Pt_VBFH125->Integral());
	Histo_Z1Pt_ttH125->Scale(1/Histo_Z1Pt_ttH125->Integral());
	Histo_Z1Pt_ZZTo4lext1->Scale(1/Histo_Z1Pt_ZZTo4lext1->Integral());
	Z1Pt_Legend->AddEntry(Histo_Z1Pt_ggH125,"ggH125","f");
	Z1Pt_Legend->AddEntry(Histo_Z1Pt_VBFH125,"VBFH125","f");
	Z1Pt_Legend->AddEntry(Histo_Z1Pt_ttH125,"ttH125","f");
	Z1Pt_Legend->AddEntry(Histo_Z1Pt_ZZTo4lext1,"qqZZ","f");
	Histo_Z1Pt_ggH125->SetLineColor(kBlue);
	Histo_Z1Pt_VBFH125->SetLineColor(kGreen);
	Histo_Z1Pt_ttH125->SetLineColor(kMagenta);
	Histo_Z1Pt_ZZTo4lext1->SetLineColor(kRed+2);
	Histo_Z1Pt_ggH125->SetTitle("Transverse momentum of Z1 boson");
	Histo_Z1Pt_ggH125->GetXaxis()->SetTitle("Pt_{Z1}/GeV");
	Histo_Z1Pt_ggH125->GetYaxis()->SetTitle("normalized to 1");
	Histo_Z1Pt_ggH125->Draw("HIST SAME");
	Histo_Z1Pt_ZZTo4lext1->Draw("HIST SAME");
	Histo_Z1Pt_VBFH125->Draw("HIST SAME");
	Histo_Z1Pt_ttH125->Draw("HIST SAME");
	Z1Pt_Legend->Draw();
	Z1Pt_Canvas->SaveAs("Z1Pt.pdf");

    //Z2Pt
	TCanvas *Z2Pt_Canvas=new TCanvas("Z2Pt","Z2Pt",1600,800);
	Z2Pt_Canvas->cd();
	TLegend* Z2Pt_Legend=new TLegend(0.75,0.75,0.9,0.9);
	Histo_Z2Pt_ggH125->Scale(1/Histo_Z2Pt_ggH125->Integral());
	Histo_Z2Pt_VBFH125->Scale(1/Histo_Z2Pt_VBFH125->Integral());
	Histo_Z2Pt_ttH125->Scale(1/Histo_Z2Pt_ttH125->Integral());
	Histo_Z2Pt_ZZTo4lext1->Scale(1/Histo_Z2Pt_ZZTo4lext1->Integral());
	Z2Pt_Legend->AddEntry(Histo_Z2Pt_ggH125,"ggH125","f");
	Z2Pt_Legend->AddEntry(Histo_Z2Pt_VBFH125,"VBFH125","f");
	Z2Pt_Legend->AddEntry(Histo_Z2Pt_ttH125,"ttH125","f");
	Z2Pt_Legend->AddEntry(Histo_Z2Pt_ZZTo4lext1,"qqZZ","f");
	Histo_Z2Pt_ggH125->SetLineColor(kBlue);
	Histo_Z2Pt_VBFH125->SetLineColor(kGreen);
	Histo_Z2Pt_ttH125->SetLineColor(kMagenta);
	Histo_Z2Pt_ZZTo4lext1->SetLineColor(kRed+2);
	Histo_Z2Pt_ggH125->SetTitle("Transverse momentum of Z2 boson");
	Histo_Z2Pt_ggH125->GetXaxis()->SetTitle("Pt_{Z2}/GeV");
	Histo_Z2Pt_ggH125->GetYaxis()->SetTitle("normalized to 1");
	Histo_Z2Pt_ggH125->Draw("HIST SAME");
	Histo_Z2Pt_ZZTo4lext1->Draw("HIST SAME");
	Histo_Z2Pt_VBFH125->Draw("HIST SAME");
	Histo_Z2Pt_ttH125->Draw("HIST SAME");
	Z2Pt_Legend->Draw();
	Z2Pt_Canvas->SaveAs("Z2Pt.pdf");

    //PhotonPt
	TCanvas *PhotonPt_Canvas=new TCanvas("PhotonPt","PhotonPt",1600,800);
	PhotonPt_Canvas->cd();
	TLegend* PhotonPt_Legend=new TLegend(0.75,0.75,0.9,0.9);
	Histo_PhotonPt_ggH125->Scale(1/Histo_PhotonPt_ggH125->Integral());
	Histo_PhotonPt_VBFH125->Scale(1/Histo_PhotonPt_VBFH125->Integral());
	Histo_PhotonPt_ttH125->Scale(1/Histo_PhotonPt_ttH125->Integral());
	Histo_PhotonPt_ZZTo4lext1->Scale(1/Histo_PhotonPt_ZZTo4lext1->Integral());
	PhotonPt_Legend->AddEntry(Histo_PhotonPt_ggH125,"ggH125","f");
	PhotonPt_Legend->AddEntry(Histo_PhotonPt_VBFH125,"VBFH125","f");
	PhotonPt_Legend->AddEntry(Histo_PhotonPt_ttH125,"ttH125","f");
	PhotonPt_Legend->AddEntry(Histo_PhotonPt_ZZTo4lext1,"qqZZ","f");
	Histo_PhotonPt_ggH125->SetLineColor(kBlue);
	Histo_PhotonPt_VBFH125->SetLineColor(kGreen);
	Histo_PhotonPt_ttH125->SetLineColor(kMagenta);
	Histo_PhotonPt_ZZTo4lext1->SetLineColor(kRed+2);
	Histo_PhotonPt_ggH125->SetTitle("Maximum transverse momentum of photon");
	Histo_PhotonPt_ggH125->GetXaxis()->SetTitle("Pt_{#gamma}/GeV");
	Histo_PhotonPt_ggH125->GetYaxis()->SetTitle("normalized to 1");
	Histo_PhotonPt_ggH125->Draw("HIST SAME");
	Histo_PhotonPt_ZZTo4lext1->Draw("HIST SAME");
	Histo_PhotonPt_VBFH125->Draw("HIST SAME");
	Histo_PhotonPt_ttH125->Draw("HIST SAME");
	PhotonPt_Legend->Draw();
	PhotonPt_Canvas->SaveAs("PhotonPt.pdf");

    //nExtraLep
	TCanvas *nExtraLep_Canvas=new TCanvas("nExtraLep","nExtraLep",1600,800);
	nExtraLep_Canvas->cd();
	TLegend* nExtraLep_Legend=new TLegend(0.75,0.75,0.9,0.9);
	Histo_nExtraLep_ggH125->Scale(1/Histo_nExtraLep_ggH125->Integral());
	Histo_nExtraLep_VBFH125->Scale(1/Histo_nExtraLep_VBFH125->Integral());
	Histo_nExtraLep_ttH125->Scale(1/Histo_nExtraLep_ttH125->Integral());
	Histo_nExtraLep_ZZTo4lext1->Scale(1/Histo_nExtraLep_ZZTo4lext1->Integral());
	//nExtraLep_Legend->AddEntry(Histo_nExtraLep_ggH125,"ggH125","f");  nema dodatnih leptona
	nExtraLep_Legend->AddEntry(Histo_nExtraLep_VBFH125,"VBFH125","f");
	nExtraLep_Legend->AddEntry(Histo_nExtraLep_ttH125,"ttH125","f");
	//nExtraLep_Legend->AddEntry(Histo_nExtraLep_ZZTo4lext1,"qqZZ","f"); nema dodatnih leptona
	Histo_nExtraLep_ggH125->SetLineColor(kBlue);
	Histo_nExtraLep_VBFH125->SetLineColor(kGreen);
	Histo_nExtraLep_ttH125->SetLineColor(kMagenta);
	Histo_nExtraLep_ZZTo4lext1->SetLineColor(kRed+2);
	Histo_nExtraLep_VBFH125->SetTitle("Number of extra leptons");
	Histo_nExtraLep_VBFH125->GetXaxis()->SetTitle("extra leptons");
	Histo_nExtraLep_VBFH125->GetYaxis()->SetTitle("normalized to 1");
	//Histo_nExtraLep_ggH125->Draw("HIST SAME");
	//Histo_nExtraLep_ZZTo4lext1->Draw("HIST SAME"); 
	Histo_nExtraLep_VBFH125->Draw("HIST SAME");
	Histo_nExtraLep_ttH125->Draw("HIST SAME");
	nExtraLep_Legend->Draw();
	nExtraLep_Canvas->SaveAs("nExtraLep.pdf");

    //nExtraZ
	TCanvas *nExtraZ_Canvas=new TCanvas("nExtraZ","nExtraZ",1600,800);
	nExtraZ_Canvas->cd();
	TLegend* nExtraZ_Legend=new TLegend(0.75,0.75,0.9,0.9);
	Histo_nExtraZ_ggH125->Scale(1/Histo_nExtraZ_ggH125->Integral());
	Histo_nExtraZ_VBFH125->Scale(1/Histo_nExtraZ_VBFH125->Integral());
	Histo_nExtraZ_ttH125->Scale(1/Histo_nExtraZ_ttH125->Integral());
	Histo_nExtraZ_ZZTo4lext1->Scale(1/Histo_nExtraZ_ZZTo4lext1->Integral());
	//nExtraZ_Legend->AddEntry(Histo_nExtraZ_ggH125,"ggH125","f");  nema dodatnih Z bozona
	nExtraZ_Legend->AddEntry(Histo_nExtraZ_VBFH125,"VBFH125","f");
	nExtraZ_Legend->AddEntry(Histo_nExtraZ_ttH125,"ttH125","f");
	//nExtraZ_Legend->AddEntry(Histo_nExtraZ_ZZTo4lext1,"qqZZ","f");  nema dodatnih Z bozona
	Histo_nExtraZ_ggH125->SetLineColor(kBlue);
	Histo_nExtraZ_VBFH125->SetLineColor(kGreen);
	Histo_nExtraZ_ttH125->SetLineColor(kMagenta);
	Histo_nExtraZ_ZZTo4lext1->SetLineColor(kRed+2);
	Histo_nExtraZ_VBFH125->SetTitle("Number of extra Z bosons");
	Histo_nExtraZ_VBFH125->GetXaxis()->SetTitle("extra Z bosons");
	Histo_nExtraZ_VBFH125->GetYaxis()->SetTitle("normalized to 1");
	//Histo_nExtraZ_ggH125->Draw("HIST SAME");
	//Histo_nExtraZ_ZZTo4lext1->Draw("HIST SAME");  
	Histo_nExtraZ_VBFH125->Draw("HIST SAME");
	Histo_nExtraZ_ttH125->Draw("HIST SAME");
	nExtraZ_Legend->Draw();
	nExtraZ_Canvas->SaveAs("nExtraZ.pdf");	
	
	//DiJetMass
	TCanvas *DiJetMass_Canvas=new TCanvas("DiJetMass","DiJetMass",1600,800);
	DiJetMass_Canvas->cd();
	TLegend* DiJetMass_Legend=new TLegend(0.75,0.75,0.9,0.9);
	Histo_DiJetMass_ggH125->Scale(1/Histo_DiJetMass_ggH125->Integral());
	Histo_DiJetMass_VBFH125->Scale(1/Histo_DiJetMass_VBFH125->Integral());
	Histo_DiJetMass_ttH125->Scale(1/Histo_DiJetMass_ttH125->Integral());
	Histo_DiJetMass_ZZTo4lext1->Scale(1/Histo_DiJetMass_ZZTo4lext1->Integral());
	DiJetMass_Legend->AddEntry(Histo_DiJetMass_ggH125,"ggH125","f");
	DiJetMass_Legend->AddEntry(Histo_DiJetMass_VBFH125,"VBFH125","f");
	DiJetMass_Legend->AddEntry(Histo_DiJetMass_ttH125,"ttH125","f");
	DiJetMass_Legend->AddEntry(Histo_DiJetMass_ZZTo4lext1,"qqZZ","f");
	Histo_DiJetMass_ggH125->SetLineColor(kBlue);
	Histo_DiJetMass_VBFH125->SetLineColor(kGreen);
	Histo_DiJetMass_ttH125->SetLineColor(kMagenta);
	Histo_DiJetMass_ZZTo4lext1->SetLineColor(kRed+2);
	Histo_DiJetMass_ZZTo4lext1->SetTitle("Invariant mass of 2 jets");
	Histo_DiJetMass_ZZTo4lext1->GetXaxis()->SetTitle("M/GeV");
	Histo_DiJetMass_ZZTo4lext1->GetYaxis()->SetTitle("normalized to 1");
	Histo_DiJetMass_ZZTo4lext1->Draw("HIST SAME");
	Histo_DiJetMass_ggH125->Draw("HIST SAME");
	Histo_DiJetMass_VBFH125->Draw("HIST SAME");
	Histo_DiJetMass_ttH125->Draw("HIST SAME");
	DiJetMass_Legend->Draw();
	DiJetMass_Canvas->SaveAs("DiJetMass.pdf");  
	//cout<<Histo_DiJetMass_ggH125->Integral()<<" "<<Histo_DiJetMass_VBFH125->Integral()<<" "<<Histo_DiJetMass_ttH125->Integral()<<" "<<Histo_DiJetMass_ZZTo4lext1->Integral()<<endl;
}

void Analyzer :: Categorize(TString s1)
{
	
	TFile *f;  
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.  
      f = new TFile(s1);
      TDirectory * dir = (TDirectory*)f->Get(s1+":/ZZTree");
      dir->GetObject("candTree",tree);
   Init(tree);
   h1 = (TH1F*)f->Get("ZZTree/Counters"); //dodano (radi citanja iz Countersa)
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   double w;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
	  w=137000*xsec*overallEventWeight/h1->GetBinContent(40);
	  float c_Mela2j = getDVBF2jetsConstant(ZZMass);
      float D_VBF2j=1./(1.+ c_Mela2j*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
	  float c_Mela1j = getDVBF1jetConstant(ZZMass);
      float D_VBF1j = 1./(1.+ c_Mela1j*p_JQCD_SIG_ghg2_1_JHUGen_JECNominal/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal));
	  float c_MelaWH = getDWHhConstant(ZZMass);
	  float D_WHh = 1./(1.+ c_MelaWH*(p_HadWH_mavjj_true_JECNominal*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_HadWH_mavjj_JECNominal*p_HadWH_SIG_ghw1_1_JHUGen_JECNominal));
	  float c_MelaZH = getDZHhConstant(ZZMass);
      float D_ZHh = 1./(1.+ c_MelaZH*(p_HadZH_mavjj_true_JECNominal*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_HadZH_mavjj_JECNominal*p_HadZH_SIG_ghz1_1_JHUGen_JECNominal));
	  if(s1.Contains("ggH125"))
	  {
	if(ZZMass>118 && ZZMass<130)
	{	
	if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && D_VBF2j>0.5 )
    Histo_ggH125->Fill(4.5,w);

    else if( nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && (D_WHh>0.5||D_ZHh>0.5) )
    Histo_ggH125->Fill(3.5,w); 

    else if( ( nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged_bTagSF==0 && (nExtraLep==1||nExtraZ>=1) )|| ( nCleanedJetsPt30==0 && nExtraLep>=1 ) )
    Histo_ggH125->Fill(2.5,w); 

    else if( nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF>=1 && nExtraLep ==0)
    Histo_ggH125->Fill(1.5,w); 
  
    else if( nExtraLep>=1 )
    Histo_ggH125->Fill(0.5,w);
	
    else if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>0.5 )
    Histo_ggH125->Fill(5.5,w);
 
    else
    Histo_ggH125->Fill(6.5,w);
	}
    }
	
	  if(s1.Contains("VBFH125"))
	  {
		  if(ZZMass>118 && ZZMass<130)
		  {
	if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && D_VBF2j>0.5 )
    Histo_VBFH125->Fill(4.5,w);

    else if( nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && (D_WHh>0.5||D_ZHh>0.5) )
    Histo_VBFH125->Fill(3.5,w); 

    else if( ( nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged_bTagSF==0 && (nExtraLep==1||nExtraZ>=1) )|| ( nCleanedJetsPt30==0 && nExtraLep>=1 ) )
    Histo_VBFH125->Fill(2.5,w); 

    else if( nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF>=1 && nExtraLep ==0)
    Histo_VBFH125->Fill(1.5,w); 
  
    else if( nExtraLep>=1 )
    Histo_VBFH125->Fill(0.5,w);
	
    else if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>0.5 )
    Histo_VBFH125->Fill(5.5,w);
 
    else
    Histo_VBFH125->Fill(6.5,w); 
		  }

	  }
	  if(s1.Contains("ttH125"))
	  {
	if(ZZMass>118 && ZZMass<130)
	{		
	if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && D_VBF2j>0.5 )
    Histo_ttH125->Fill(4.5,w);

    else if( nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && (D_WHh>0.5||D_ZHh>0.5) )
    Histo_ttH125->Fill(3.5,w); 

    else if( ( nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged_bTagSF==0 && (nExtraLep==1||nExtraZ>=1) )|| ( nCleanedJetsPt30==0 && nExtraLep>=1 ) )
    Histo_ttH125->Fill(2.5,w); 

    else if( nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF>=1 && nExtraLep ==0)
    Histo_ttH125->Fill(1.5,w); 
  
    else if( nExtraLep>=1 )
    Histo_ttH125->Fill(0.5,w);
	
    else if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>0.5 )
    Histo_ttH125->Fill(5.5,w);
 
    else
    Histo_ttH125->Fill(6.5,w); 
	}
	  }
	  
	  if(s1.Contains("ZZTo4lext1"))
	  {

    if(ZZMass>118 && ZZMass<130)
	{	
	if( nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && D_VBF2j>0.5 )
    Histo_ZZTo4lext1->Fill(4.5,w);

    else if( nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && (D_WHh>0.5||D_ZHh>0.5) )
    Histo_ZZTo4lext1->Fill(3.5,w); 

    else if( ( nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged_bTagSF==0 && (nExtraLep==1||nExtraZ>=1) )|| ( nCleanedJetsPt30==0 && nExtraLep>=1 ) )
    Histo_ZZTo4lext1->Fill(2.5,w); 

    else if( nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF>=1 && nExtraLep ==0)
    Histo_ZZTo4lext1->Fill(1.5,w); 
  
    else if( nExtraLep>=1 )
    Histo_ZZTo4lext1->Fill(0.5,w);
	
    else if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>0.5 )
    Histo_ZZTo4lext1->Fill(5.5,w);
 
    else
    Histo_ZZTo4lext1->Fill(6.5,w);
	}
	  }
  }
}

void Analyzer :: Categorize_Display()
{
	int i;
	gStyle->SetOptStat(0);
	const char *categories[7]={"ttH-leptonic tagged ","ttH-hadronic tagged ","VH-leptonic tagged ","VH-hadronic tagged ","VBF-2jet tagged "," VBF-1jet tagged"," Untagged"};
	Histo_ggH125->SetTitle("ggH");
	Histo_VBFH125->SetTitle("VBF");
	Histo_ttH125->SetTitle("ttH");
	Histo_ZZTo4lext1->SetTitle("qqZZ");
	for (i=1;i<=7;i++) 
	{
		Histo_ggH125->GetXaxis()->SetBinLabel(i,categories[i-1]);
		Histo_VBFH125->GetXaxis()->SetBinLabel(i,categories[i-1]);
		Histo_ttH125->GetXaxis()->SetBinLabel(i,categories[i-1]);
		Histo_ZZTo4lext1->GetXaxis()->SetBinLabel(i,categories[i-1]);
	}
	Histo_ggH125->SetFillColor(kBlue);
	Histo_VBFH125->SetFillColor(kGreen+2);
	Histo_ttH125->SetFillColor(kMagenta);
	Histo_ZZTo4lext1->SetFillColor(kRed);
	Histo_ggH125->GetYaxis()->SetTitle("Expected events");
	Histo_VBFH125->GetYaxis()->SetTitle("Expected events");
	Histo_ttH125->GetYaxis()->SetTitle("Expected events");
	Histo_ZZTo4lext1->GetYaxis()->SetTitle("Expected events");
	
	TCanvas *canvas_categories_1=new TCanvas("ggH125","ggH125",1400,800);
	canvas_categories_1->Divide(2,2);
	
	Histo_ggH125->GetYaxis()->SetTitleSize(0.05);
	Histo_ggH125->GetXaxis()->SetLabelSize(0.08);
    Histo_ggH125->GetYaxis()->SetLabelSize(0.04);
    Histo_ggH125->SetBarOffset(0.1);
	Histo_ggH125->SetBarWidth(0.8);
	canvas_categories_1->cd(1);
	gPad->SetLeftMargin(0.35);
	Histo_ggH125->Draw("hbar");
	
	Histo_VBFH125->GetYaxis()->SetTitleSize(0.05);
	Histo_VBFH125->GetXaxis()->SetLabelSize(0.08);
	Histo_VBFH125->GetYaxis()->SetLabelSize(0.04);
	Histo_VBFH125->SetBarOffset(0.1);
	Histo_VBFH125->SetBarWidth(0.8);
	canvas_categories_1->cd(2);
	gPad->SetLeftMargin(0.25);
	Histo_VBFH125->Draw("hbar");
	
	Histo_ttH125->GetYaxis()->SetTitleSize(0.05);
	Histo_ttH125->GetXaxis()->SetLabelSize(0.08);
	Histo_ttH125->GetYaxis()->SetLabelSize(0.04);
	Histo_ttH125->SetBarOffset(0.1);
	Histo_ttH125->SetBarWidth(0.8);
	canvas_categories_1->cd(3);
	gPad->SetLeftMargin(0.35);
	Histo_ttH125->Draw("hbar");
	
	Histo_ZZTo4lext1->GetYaxis()->SetTitleSize(0.05);
	Histo_ZZTo4lext1->GetXaxis()->SetLabelSize(0.08);
	Histo_ZZTo4lext1->GetYaxis()->SetLabelSize(0.04);
	Histo_ZZTo4lext1->SetBarOffset(0.1);
	Histo_ZZTo4lext1->SetBarWidth(0.8);
	canvas_categories_1->cd(4);
	gPad->SetLeftMargin(0.25);
	Histo_ZZTo4lext1->GetYaxis()->SetRangeUser(0., 80.);
	Histo_ZZTo4lext1->Draw("hbar");
	
	canvas_categories_1->SaveAs("Categorization.pdf");
	
	double VBF2j_exp_events=Histo_ggH125->Integral(5,5)+Histo_VBFH125->Integral(5,5)+Histo_ttH125->Integral(5,5)+Histo_ZZTo4lext1->Integral(5,5);
	double VHh_exp_events=Histo_ggH125->Integral(4,4)+Histo_VBFH125->Integral(4,4)+Histo_ttH125->Integral(4,4)+Histo_ZZTo4lext1->Integral(4,4);
	double VHl_exp_events=Histo_ggH125->Integral(3,3)+Histo_VBFH125->Integral(3,3)+Histo_ttH125->Integral(3,3)+Histo_ZZTo4lext1->Integral(3,3);
	double ttHh_exp_events=Histo_ggH125->Integral(2,2)+Histo_VBFH125->Integral(2,2)+Histo_ttH125->Integral(2,2)+Histo_ZZTo4lext1->Integral(2,2);
	double ttHl_exp_events=Histo_ggH125->Integral(1,1)+Histo_VBFH125->Integral(1,1)+Histo_ttH125->Integral(1,1)+Histo_ZZTo4lext1->Integral(1,1);
	double VBF1j_exp_events=Histo_ggH125->Integral(6,6)+Histo_VBFH125->Integral(6,6)+Histo_ttH125->Integral(6,6)+Histo_ZZTo4lext1->Integral(6,6);
	double Untagged_exp_events=Histo_ggH125->Integral(7,7)+Histo_VBFH125->Integral(7,7)+Histo_ttH125->Integral(7,7)+Histo_ZZTo4lext1->Integral(7,7);
	
	
	double exp_events[7]={ttHl_exp_events,ttHh_exp_events,VHl_exp_events,VHh_exp_events,VBF2j_exp_events,VBF1j_exp_events,Untagged_exp_events};
	
	
	for(i=1;i<=7;i++)
	{
	  Histo_ggH125->SetBinContent(i,Histo_ggH125->Integral(i,i)/exp_events[i-1]);	
	  Histo_VBFH125->SetBinContent(i,Histo_VBFH125->Integral(i,i)/exp_events[i-1]);
	  Histo_ttH125->SetBinContent(i,Histo_ttH125->Integral(i,i)/exp_events[i-1]);
	  Histo_ZZTo4lext1->SetBinContent(i,Histo_ZZTo4lext1->Integral(i,i)/exp_events[i-1]);
	}
	THStack *stack = new THStack("stack","Categorization (simulation) ; ;fraction");
	TCanvas *master_canvas=new TCanvas("master","master",1400,800);
	TLegend* legend=new TLegend(0.9,0.7,0.98,0.9);
	legend->AddEntry(Histo_ggH125,"ggH","f");
	legend->AddEntry(Histo_VBFH125,"VBFH","f");
	legend->AddEntry(Histo_ttH125,"ttH","f");
	legend->AddEntry(Histo_ZZTo4lext1,"qqZZ","f");
	gPad->SetLeftMargin(0.2);
	TLatex *t1=new TLatex(0.3,0.35,"0.62 exp. events");
	TLatex *t2=new TLatex(0.3,1.35,"1.10 exp. events");
	TLatex *t3=new TLatex(0.3,2.35,"1.16 exp. events");
	TLatex *t4=new TLatex(0.3,3.35,"10.17 exp. events");
	TLatex *t5=new TLatex(0.3,4.35,"12.28 exp. events");
	TLatex *t6=new TLatex(0.3,5.35,"26.77 exp. events");
	TLatex *t7=new TLatex(0.3,6.35,"235.09 exp. events");
	t1->SetTextColor(kWhite);
	t2->SetTextColor(kWhite);
	t1->SetTextColor(kWhite);
	t3->SetTextColor(kWhite);
	t4->SetTextColor(kWhite);
	t5->SetTextColor(kWhite);
	t6->SetTextColor(kWhite);
	t7->SetTextColor(kWhite);
	stack->Add(Histo_ggH125);
	stack->Add(Histo_VBFH125);
	stack->Add(Histo_ttH125);
	stack->Add(Histo_ZZTo4lext1);
	stack->Draw("hbar");
	legend->Draw();
	t1->Draw();
	t2->Draw();
	t3->Draw();
	t4->Draw();
	t5->Draw();
	t6->Draw();
	t7->Draw();
	master_canvas->SaveAs("Categorization_stack.pdf");
   //stack
}
 
void Analyzer :: TMVAMultiClass()
{
	// This loads the library
   TMVA::Tools::Instance();
 
   // to get access to the GUI and all tmva macros
   //
   //     TString tmva_dir(TString(gRootDir) + "/tmva");
   //     if(gSystem->Getenv("TMVASYS"))
   //        tmva_dir = TString(gSystem->Getenv("TMVASYS"));
   //     gROOT->SetMacroPath(tmva_dir + "/test/:" + gROOT->GetMacroPath() );
   //     gROOT->ProcessLine(".L TMVAMultiClassGui.C");
 
 
   //---------------------------------------------------------------
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;
   Use["MLP"]             = 1;
   Use["BDTG"]            = 1;
#ifdef R__HAS_TMVAGPU
   Use["DL_CPU"]          = 1;
   Use["DL_GPU"]          = 1;
#else
   Use["DL_CPU"]          = 1;
   Use["DL_GPU"]          = 0;
#endif
   Use["FDA_GA"]          = 0;
   Use["PDEFoam"]         = 1;
 
   //---------------------------------------------------------------
 
  
   // Create a new root output file.
   TString outfileName = "TMVAMulticlass.root";
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
 
   TMVA::Factory *factory = new TMVA::Factory( "TMVAMulticlass", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
 
   dataloader->AddVariable( "D_VBF2j:=1./(1.+ p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal)", 'F' );
   dataloader->AddVariable("D_VBF1j := 1./(1.+p_JQCD_SIG_ghg2_1_JHUGen_JECNominal/(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal))",'F');
   dataloader->AddVariable("D_WHh := 1./(1.+ (p_HadWH_mavjj_true_JECNominal*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_HadWH_mavjj_JECNominal*p_HadWH_SIG_ghw1_1_JHUGen_JECNominal))",'F');
   dataloader->AddVariable("D_ZHh := 1./(1.+ (p_HadZH_mavjj_true_JECNominal*p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal)/(p_HadZH_mavjj_JECNominal*p_HadZH_SIG_ghz1_1_JHUGen_JECNominal))",'F');
   
   dataloader->AddVariable("ZZMass",'F');
   dataloader->AddVariable( "nCleanedJetsPt30BTagged", 'F' );
   dataloader->AddVariable( "nExtraLep", 'F' );
   dataloader->AddVariable( "nCleanedJetsPt30", 'F' );
   
   
   TFile *input1(0);
   TString fname1 = "/home/public/data/2018_MC/ggH125/ZZ4lAnalysis.root";
   if (!gSystem->AccessPathName( fname1 )) {
      input1 = TFile::Open( fname1 ); // check if file in local directory exists
   }
   if (!input1) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   TFile *input2(0);
   TString fname2 = "/home/public/data/2018_MC/VBFH125/ZZ4lAnalysis.root";
   if (!gSystem->AccessPathName( fname2 )) {
      input2 = TFile::Open( fname2 ); // check if file in local directory exists
   }
   
   TFile *input3(0);
   TString fname3 = "/home/public/data/2018_MC/ttH125/ZZ4lAnalysis.root";
   if (!gSystem->AccessPathName( fname3 )) {
      input3 = TFile::Open( fname3 ); // check if file in local directory exists
   }
   TFile *input4(0);
   TString fname4 = "/home/public/data/2018_MC/ZZTo4lext1/ZZ4lAnalysis.root";
   if (!gSystem->AccessPathName( fname4 )) {
      input4 = TFile::Open( fname4 ); // check if file in local directory exists
   }
    
   TTree *Tree_ggH  = (TTree*)input1->Get("ZZTree/candTree");
   TTree *Tree_VBFH = (TTree*)input2->Get("ZZTree/candTree");
   TTree *Tree_ttH = (TTree*)input3->Get("ZZTree/candTree");
   TTree *Tree_qqZZ = (TTree*)input4->Get("ZZTree/candTree");
       
   
   gROOT->cd( outfileName+TString(":/") );
   dataloader->AddTree(Tree_ggH,"ggH");
   dataloader->AddTree(Tree_VBFH,"VBFH");
   dataloader->AddTree(Tree_ttH,"ttH");
   dataloader->AddTree(Tree_qqZZ,"qqZZ");
   
   dataloader->SetWeightExpression ("137000/28744188*xsec*overallEventWeight","ggH");
   dataloader->SetWeightExpression ("137000/1819984.75*xsec*overallEventWeight","VBFH");
   dataloader->SetWeightExpression ("137000/257544.9375*xsec*overallEventWeight","ttH");
   dataloader->SetWeightExpression ("137000/8398762*xsec*overallEventWeight","qqZZ");
   
    TCut cut="ZZMass>118 && ZZMass<130";
    dataloader->PrepareTrainingAndTestTree(cut, "SplitMode=Random:NormMode=NumEvents:!V" );
   
   //dataloader->PrepareTrainingAndTestTree(cut,"nTrain_ggH=1000:nTrain_VBFH=1000:nTrain_ttH=1000:nTrain_qqZZ=1000:SplitMode=Random:NormMode=NumEvents:!V" );
   factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
   
   /*if (Use["BDTG"]) // gradient boosted decision trees
      factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
   if (Use["MLP"]) // neural network
      factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:NCycles=1000:HiddenLayers=N+5,5:TestRate=5:EstimatorType=MSE");
   if (Use["FDA_GA"]) // functional discriminant with GA minimizer
      factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_GA", "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );
   if (Use["PDEFoam"]) // PDE-Foam approach
      factory->BookMethod( dataloader,  TMVA::Types::kPDEFoam, "PDEFoam", "!H:!V:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
 
 
   if (Use["DL_CPU"]) {
      TString layoutString("Layout=TANH|100,TANH|50,TANH|10,LINEAR");
      TString trainingStrategyString("TrainingStrategy=Optimizer=ADAM,LearningRate=1e-3,"
                                     "TestRepetitions=1,ConvergenceSteps=10,BatchSize=100");
      TString nnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                        "WeightInitialization=XAVIERUNIFORM:Architecture=GPU");
      nnOptions.Append(":");
      nnOptions.Append(layoutString);
      nnOptions.Append(":");
      nnOptions.Append(trainingStrategyString);
      factory->BookMethod(dataloader, TMVA::Types::kDL, "DL_CPU", nnOptions);
   }
   if (Use["DL_GPU"]) {
      TString layoutString("Layout=TANH|100,TANH|50,TANH|10,LINEAR");
      TString trainingStrategyString("TrainingStrategy=Optimizer=ADAM,LearningRate=1e-3,"
                                     "TestRepetitions=1,ConvergenceSteps=10,BatchSize=100");
      TString nnOptions("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                        "WeightInitialization=XAVIERUNIFORM:Architecture=GPU");
      nnOptions.Append(":");
      nnOptions.Append(layoutString);
      nnOptions.Append(":");
      nnOptions.Append(trainingStrategyString);
      factory->BookMethod(dataloader, TMVA::Types::kDL, "DL_GPU", nnOptions);
   }*/
 
 
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
 
   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
 
   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
 
   // --------------------------------------------------------------
 
   // Save the output
   outputFile->Close();
 
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAMulticlass is done!" << std::endl;
 
   delete factory;
   delete dataloader;
 
   // Launch the GUI for the root macros
   //if (!gROOT->IsBatch()) TMVAMultiClassGui( outfileName ); 
 
}	
 


