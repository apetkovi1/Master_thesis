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
}

extern "C" float getDVBF2jetsConstant(float ZZMass){
  return DVBF2jetsSpline.eval(ZZMass, false);
}

extern "C" float getDVBF1jetConstant(float ZZMass){
  return DVBF1jetSpline.eval(ZZMass, false);
}

void Analyzer::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Analyzer.C
//      root> Analyzer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
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
	  if(s1.Contains("ggH125"))
	  {
		if(nExtraLep==0 && (((nCleanedJetsPt30==2 || nCleanedJetsPt30==3)&& nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF==0))&& D_VBF2j>0.5)
        //VBF_2j_ggH125++*w;
        Histo_VBF_2j_ggH125->Fill(2,w);	
		if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>0.5 )
		Histo_VBF_1j_ggH125->Fill(2,w);
	  }
	  if(s1.Contains("VBFH125"))
	  {
		if(nExtraLep==0 && (((nCleanedJetsPt30==2 || nCleanedJetsPt30==3)&& nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF==0))&& D_VBF2j>0.5)
        //VBF_2j_VBFH125++*w;
        Histo_VBF_2j_VBFH125->Fill(2,w);
		if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>0.5 )
        Histo_VBF_1j_VBFH125->Fill(2,w);		
	  }
	  if(s1.Contains("ttH125"))
	  {
		if(nExtraLep==0 && (((nCleanedJetsPt30==2 || nCleanedJetsPt30==3)&& nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF==0))&& D_VBF2j>0.5)
        //VBF_2j_ttH125++*w;
        Histo_VBF_2j_ttH125->Fill(2,w);	
		if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>0.5 )
		Histo_VBF_1j_ttH125->Fill(2,w);
	  }
	  if(s1.Contains("ZZTo4lext1"))
	  {
		if(nExtraLep==0 && (((nCleanedJetsPt30==2 || nCleanedJetsPt30==3)&& nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF==0))&& D_VBF2j>0.5)
        Histo_VBF_2j_ZZTo4lext1->Fill(2,w);	
        if( nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>0.5 )
        Histo_VBF_1j_ZZTo4lext1->Fill(2,w);			
	  }
  }
}

void Analyzer :: Categorize_Display()
{
	//cout<<VBF_2j_ggH125<<" "<<VBF_2j_VBFH125<<" "<<VBF_2j_ttH125<<" "<<VBF_2j_ZZTo4lext1<<endl;
	cout<<Histo_VBF_2j_ggH125->Integral()<<" "<<Histo_VBF_2j_VBFH125->Integral()<<" "<<Histo_VBF_2j_ttH125->Integral()<<" "<<Histo_VBF_2j_ZZTo4lext1->Integral()<<endl;
	cout<<Histo_VBF_1j_ggH125->Integral()<<" "<<Histo_VBF_1j_VBFH125->Integral()<<" "<<Histo_VBF_1j_ttH125->Integral()<<" "<<Histo_VBF_1j_ZZTo4lext1->Integral()<<endl;
}
   
 


