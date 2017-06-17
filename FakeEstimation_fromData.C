#include <memory>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TH2.h>
#include <TH3.h>
#include <THStack.h>
#include <TCanvas.h>
#include "TTree.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TText.h"
#include "treeBase_data.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <utility> 
#include <TLorentzVector.h>
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TColor.h"
#include "TLine.h"
#include "TROOT.h"
#include <TStyle.h>
#include "TVector3.h"

using namespace std;

void prn(int n, int N, int d, TString& cf, TTree* ch) {

    if (cf != (ch->GetCurrentFile()->GetName())) {
        cout << "\n >> Current File: " << (ch->GetCurrentFile()->GetName()) << endl;
        cf = (ch->GetCurrentFile()->GetName());
    }

    if (n % d == 0) {
        //        fprintf(stdout, "\rProcessed events: \033[1;36;40m%6d of %6d\033[0m", n, N);
        fprintf(stdout, "\rProcessed events: %6d of %6d", n, N);
        fflush(stdout);
    }
    if (n == (N - 1)) {

        fprintf(stdout, "\n");
        fflush(stdout);
    }
}

void CanvasCreator(TH2D * h, string x_title, string y_title, string title, TString pathPNG) {
	
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptStat(0);

	TCanvas * canvas = new TCanvas(h->GetName(),h->GetName(),350,350,600,600);
        canvas->SetHighLightColor(2);
        canvas->Range(0,0,1,1);
        canvas->SetFillColor(0);
        canvas->SetBorderMode(0);
        canvas->SetBorderSize(2);
        canvas->SetTickx(1);
        canvas->SetTicky(1);
        canvas->SetLeftMargin(0.16);
        canvas->SetRightMargin(0.02);
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.13);
        canvas->SetFrameFillStyle(0);
        canvas->SetFrameBorderMode(0);

	h->SetTitle(title.c_str());
	h->SetTitleOffset(1.5);;
	h->GetXaxis()->SetTitle(x_title.c_str());
	h->GetXaxis()->SetTitleOffset(1.25);;
	h->GetYaxis()->SetTitle(y_title.c_str());
	h->GetYaxis()->SetTitleOffset(1.85);;
	h->Draw("COLTEXT");

	pathPNG.Append(h->GetName());
	pathPNG.Append(".png");
        canvas->SaveAs(pathPNG);
}
void CanvasCreator(TH1D * h, string x_title, string y_title, TString pathPNG) {
	
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	TCanvas * canvas = new TCanvas(h->GetName(),h->GetName(),350,350,500,500);
        canvas->SetHighLightColor(2);
        canvas->Range(0,0,1,1);
        canvas->SetFillColor(0);
        canvas->SetBorderMode(0);
        canvas->SetBorderSize(2);
        canvas->SetTickx(1);
        canvas->SetTicky(1);
        canvas->SetLeftMargin(0.16);
        canvas->SetRightMargin(0.02);
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.13);
        canvas->SetFrameFillStyle(0);
        canvas->SetFrameBorderMode(0);

	h->GetXaxis()->SetTitle(x_title.c_str());
	h->GetXaxis()->SetTitleOffset(1.25);;
	h->GetYaxis()->SetTitle(y_title.c_str());
	h->GetYaxis()->SetTitleOffset(1.85);;
	h->SetLineColor(kBlue);
	h->Draw();

	TString name = h->GetName();
	if (name.Contains("pt") && name.Contains("loose")) canvas->SetLogy();

	pathPNG.Append(h->GetName());
	pathPNG.Append(".png");
        canvas->SaveAs(pathPNG);
}

struct Sort{
  bool operator() (std::pair <TH1D*, std::pair <double, double > > comb1, std::pair <TH1D*, std::pair <double, double > > comb2)
        //{return ((comb1.first)->Integral() > (comb2.first)->Integral());}
        {return (((comb1.second).first) < ((comb2.second).first) && ((comb1.second).second) < ((comb2.second).second));}
};

void FakeEstimation_fromData() {

	TH1::SetDefaultSumw2();
	gStyle->SetOptTitle(0);

	double N_ll = 0.;
	double nevents = 0.;

	double N_0t = 0.;
	double N_1t = 0.;
	double N_2t = 0.;

	TH1D * met_0t_wjet = new TH1D("met_0t_wjet","met_0t_wjet",30,0.,300.);
	TH1D * met_1t_wjet = new TH1D("met_1t_wjet","met_1t_wjet",30,0.,300.);
	TH1D * met_2t_wjet = new TH1D("met_2t_wjet","met_2t_wjet",30,0.,300.);

	TH1D * met_0t_data = new TH1D("met_0t_data","met_0t_data",30,0.,300.);
	TH1D * met_1t_data = new TH1D("met_1t_data","met_1t_data",30,0.,300.);
	TH1D * met_2t_data = new TH1D("met_2t_data","met_2t_data",30,0.,300.);

	TH1D * n_0t = new TH1D("n_0t","n_0t",1,1.,2.);
	TH1D * n_1t = new TH1D("n_1t","n_1t",1,1.,2.);
	TH1D * n_2t = new TH1D("n_2t","n_2t",1,1.,2.);

	TString fileNames[8] = {"200317/DiTauData/Tau_Run2016B_03Feb2017_v2/H2TauTauTreeProducerTauTau/tree.root","200317/DiTauData/Tau_Run2016C_03Feb2017/H2TauTauTreeProducerTauTau/tree.root","200317/DiTauData/Tau_Run2016D_03Feb2017/H2TauTauTreeProducerTauTau/tree.root","200317/DiTauData/Tau_Run2016E_03Feb2017/H2TauTauTreeProducerTauTau/tree.root","200317/DiTauData/Tau_Run2016F_03Feb2017/H2TauTauTreeProducerTauTau/tree.root","200317/DiTauData/Tau_Run2016G_03Feb2017/H2TauTauTreeProducerTauTau/tree.root","200317/DiTauData/Tau_Run2016H_03Feb2017_v2/H2TauTauTreeProducerTauTau/tree.root","200317/DiTauData/Tau_Run2016H_03Feb2017_v3/H2TauTauTreeProducerTauTau/tree.root"};

	double IsoCut = 1.5;

	TH1D * prompt_rate = new TH1D("pr","pr",1,1.,2.);
	prompt_rate->SetBinContent(1,0.710724);
	prompt_rate->SetBinError(1,0.0083010);
	//prompt_rate->SetBinContent(1,0.639389);
	//prompt_rate->SetBinError(1,0.00479055);

	TH1D * fake_rate_wjet = new TH1D("fr_wjet","fr_wjet",1,1.,2.);
	fake_rate_wjet->SetBinContent(1,0.28725);
	fake_rate_wjet->SetBinError(1,0.0221622); // from W+jets
	//fake_rate_wjet->SetBinContent(1,0.13478);
	//fake_rate_wjet->SetBinError(1,0.00656193); // from W+jets

	TH1D * fake_rate_data = new TH1D("fr_data","fr_data",1,1.,2.);
	fake_rate_data->SetBinContent(1,0.276816);
	fake_rate_data->SetBinError(1,0.00069182); // from Data
	//fake_rate_data->SetBinContent(1,0.161887);
	//fake_rate_data->SetBinError(1,0.00059797); // from Data

	double p = prompt_rate->GetBinContent(1);
	double f_wjet = fake_rate_wjet->GetBinContent(1);
	double f_data = fake_rate_data->GetBinContent(1);

	for (int file = 0; file < 8; file++) {

	TFile* myFile = TFile::Open(fileNames[file]);
	TTree *myTree = (TTree*)(myFile->Get("tree"));

	treeBase_data tr(myTree);

	if (tr.fChain == 0) return;

	Long64_t nentries = tr.fChain->GetEntries();

	TString CurrentFile = "";

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        Long64_t ientry = tr.LoadTree(jentry);
        if (ientry < 0) break;
        tr.fChain->GetEntry(jentry);

	if (!(!tr.veto_dilepton && !tr.veto_thirdlepton && !tr.veto_otherlepton)) continue;
	if (!(tr.l1_againstMuon3>0.5 && tr.l1_againstElectronMVA6>0.5 && tr.l1_pt>40.)) continue;
	if (!(tr.l2_againstMuon3>0.5 && tr.l2_againstElectronMVA6>0.5 && tr.l2_pt>40.)) continue;
	if (!(tr.l1_charge != tr.l2_charge)) continue;
	
	double summt = tr.pfmet_mt1+tr.pfmet_mt2;
	
	if (!((tr.trigger_ditau35 && tr.trigger_matched_ditau35) || (tr.trigger_ditau35_combiso && tr.trigger_matched_ditau35_combiso))) continue;
	if (!(tr.Flag_HBHENoiseFilter && tr.Flag_HBHENoiseIsoFilter && tr.Flag_EcalDeadCellTriggerPrimitiveFilter && tr.Flag_goodVertices && tr.Flag_eeBadScFilter && tr.Flag_globalTightHalo2016Filter && tr.passBadMuonFilter && tr.passBadChargedHadronFilter && tr.badMuonMoriond2017 && tr.badCloneMuonMoriond2017)) continue;

	if (!(tr.pfmet_pt > 30)) continue;
	if (!(tr.mvis > 85 || tr.mvis < 55)) continue;
	if (!(tr.mt2 > 20)) continue;
	// +++ SR I +++
	if (!(tr.mt2 > 90)) continue;
	// +++ SR II +++
	//if (!(tr.mt2 < 90 && summt > 250 && tr.n_bjets == 0)) continue;

	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>1.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT>1.5) N_ll++;

	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5) nevents++;

	double weight_0t_wjet = f_wjet*f_wjet*p*p/pow(p-f_wjet,2);
	double weight_1t_wjet = f_wjet*f_wjet*p*(1-p)/pow(p-f_wjet,2);
	double weight_2t_wjet = f_wjet*f_wjet*(1-p)*(1-p)/pow(p-f_wjet,2);

	double weight_0t_data = f_data*f_data*p*p/pow(p-f_data,2);
	double weight_1t_data = f_data*f_data*p*(1-p)/pow(p-f_data,2);
	double weight_2t_data = f_data*f_data*(1-p)*(1-p)/pow(p-f_data,2);

	int both_loose_non_tight = 0;
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l1_byIsolationMVArun2v1DBoldDMwLT<4.5) both_loose_non_tight++;
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l2_byIsolationMVArun2v1DBoldDMwLT<4.5) both_loose_non_tight++;

	if (both_loose_non_tight == 2) {
		N_0t++;	
		n_0t->Fill(1.);
		met_0t_wjet->Fill(tr.pfmet_pt,weight_0t_wjet);
		met_0t_data->Fill(tr.pfmet_pt,weight_0t_data);
	}

	int one_loose_one_tight = 0;
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l1_byIsolationMVArun2v1DBoldDMwLT<4.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5) one_loose_one_tight++;
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l2_byIsolationMVArun2v1DBoldDMwLT<4.5 && tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5) one_loose_one_tight++;

	if (one_loose_one_tight == 1) {
		N_1t++;	
                n_1t->Fill(1.);
		met_1t_wjet->Fill(tr.pfmet_pt,weight_1t_wjet);
		met_1t_data->Fill(tr.pfmet_pt,weight_1t_data);
	}

	int both_tight = 0;
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5) both_tight++;
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5) both_tight++;

	if (both_tight == 2) {
		N_2t++;	
                n_2t->Fill(1.);
		met_2t_wjet->Fill(tr.pfmet_pt,weight_2t_wjet);
		met_2t_data->Fill(tr.pfmet_pt,weight_2t_data);
	}
}
}

	cout<<"N_ll: "<<N_ll<<endl;
	cout<<"N_0t: "<<N_0t<<endl;
	cout<<"N_1t: "<<N_1t<<endl;
	cout<<"N_2t: "<<N_2t<<endl;
	cout<<"sum: "<<N_0t+N_1t+N_2t<<endl;

	cout<<" Data Yields: "<<nevents<<endl;
/*
	TH1D * prompt_rate = new TH1D("pr","pr",1,1.,2.);
	prompt_rate->SetBinContent(1,0.639389);
	prompt_rate->SetBinError(1,0.00479055);

	TH1D * fake_rate = new TH1D("fr","fr",1,1.,2.);
	fake_rate->SetBinContent(1,0.13478);
	fake_rate->SetBinError(1,0.00656193);

	double p = prompt_rate->GetBinContent(1);
	double f = fake_rate->GetBinContent(1);
*/
	cout<<"--- fr from WJets -----------------------------------------------------------------------------------"<<endl;
	double single_fake_estimation = p*f_wjet*(-2.*f_wjet*p*N_0t + (f_wjet*(1-p)+p*(1-f_wjet))*N_1t -2.*(1-p)*(1-f_wjet)*N_2t)/pow(p-f_wjet,2);
	cout<<"))))))))))))))))))))))))  single fake estimation: "<<single_fake_estimation<<endl;

	double double_fake_estimation = f_wjet*f_wjet*(p*p*N_0t - p*(1-p)*N_1t + (1-p)*(1-p)*N_2t)/pow(p-f_wjet,2);
	cout<<"))))))))))))))))))))))))  double fake estimation: "<<double_fake_estimation<<endl;
	cout<<"))))))))))))))))))))))))  	contribution from n0t: "<<f_wjet*f_wjet*(p*p*N_0t)/pow(p-f_wjet,2)<<endl;
	cout<<"))))))))))))))))))))))))  	contribution from n1t: "<<f_wjet*f_wjet*(-1*p*(1-p)*N_1t)/pow(p-f_wjet,2)<<endl;
	cout<<"))))))))))))))))))))))))  	contribution from n2t: "<<f_wjet*f_wjet*((1-p)*(1-p)*N_2t)/pow(p-f_wjet,2)<<endl;
	
	double double_prompt_estimation = p*p*(f_wjet*f_wjet*N_0t - f_wjet*(1-f_wjet)*N_1t + (1-f_wjet)*(1-f_wjet)*N_2t)/pow(p-f_wjet,2);
	cout<<"))))))))))))))))))))))))  double prompt estimation: "<<double_prompt_estimation<<endl;
	
	cout<<"--- fr from QCD -----------------------------------------------------------------------------------"<<endl;
	single_fake_estimation = p*f_data*(-2.*f_data*p*N_0t + (f_data*(1-p)+p*(1-f_data))*N_1t -2.*(1-p)*(1-f_data)*N_2t)/pow(p-f_data,2);
	cout<<"))))))))))))))))))))))))  single fake estimation: "<<single_fake_estimation<<endl;

	double_fake_estimation = f_data*f_data*(p*p*N_0t - p*(1-p)*N_1t + (1-p)*(1-p)*N_2t)/pow(p-f_data,2);
	cout<<"))))))))))))))))))))))))  double fake estimation: "<<double_fake_estimation<<endl;
	
	double_prompt_estimation = p*p*(f_data*f_data*N_0t - f_data*(1-f_data)*N_1t + (1-f_data)*(1-f_data)*N_2t)/pow(p-f_data,2);
	cout<<"))))))))))))))))))))))))  double prompt estimation: "<<double_prompt_estimation<<endl;
	
	
        //TH1D * fake_rate = (TH1D *) fake_rate_wjet->Clone("fr_w_clone");
        TH1D * fake_rate = (TH1D *) fake_rate_data->Clone("fr_data_clone");
	TH1D * p_clone = (TH1D *) prompt_rate->Clone("p_clone");
	p_clone->Add(fake_rate,-1);

	TH1D * c0 = (TH1D *) prompt_rate->Clone("c0");
	c0->Multiply(prompt_rate);
	c0->Multiply(fake_rate);
	c0->Multiply(fake_rate);
	c0->Multiply(n_0t);
	c0->Divide(p_clone);
	c0->Divide(p_clone);

	TH1D * c1 = (TH1D *) prompt_rate->Clone("c1");
	c1->Multiply(fake_rate);
	c1->Multiply(fake_rate);

	TH1D * c2 = (TH1D *) prompt_rate->Clone("c2");
	c2->Multiply(fake_rate);
	c2->Multiply(fake_rate);
	c2->Multiply(prompt_rate);

	TH1D * c3 = (TH1D *) prompt_rate->Clone("c3");
	c3->Multiply(fake_rate);
	c3->Multiply(prompt_rate);

	c1->Add(c3);	
	c1->Add(c2,-2);	
	c1->Multiply(n_1t);
	c1->Divide(p_clone);
	c1->Divide(p_clone);

	TH1D * c4 = (TH1D *) prompt_rate->Clone("c4");
	c4->Multiply(fake_rate);

	TH1D * c5 = (TH1D *) prompt_rate->Clone("c5");
	c5->Multiply(fake_rate);
	c5->Multiply(fake_rate);

	c4->Add(c5,-1);
	c4->Add(c3,-1);
	c4->Add(c2);
	c4->Multiply(n_2t);
	c4->Divide(p_clone);
	c4->Divide(p_clone);

	c1->Add(c0,-2);
	c1->Add(c4,-2);

	cout<<"PF ++++++++++++++++++++++++++++++++++++++   "<<c1->GetBinContent(1)<<"   +-   "<<c1->GetBinError(1)<<endl;

// -------------------------------------------------  
//
//
	TH1D * clone_n0t = (TH1D *) n_0t->Clone("clone_n0t");
	TH1D * clone_n1t = (TH1D *) n_1t->Clone("clone_n1t");

	TH1D * c6 = (TH1D *) prompt_rate->Clone("c6");
	c6->Multiply(prompt_rate);

	TH1D * c7 = (TH1D *) prompt_rate->Clone("c7");

	clone_n0t->Add(n_1t);
	clone_n0t->Add(n_2t);
	c6->Multiply(clone_n0t);

	clone_n1t->Add(n_2t,+2);
	c7->Multiply(clone_n1t);

	c6->Add(c7,-1);
	c6->Add(n_2t);

	c6->Divide(p_clone);
        c6->Divide(p_clone);
	c6->Multiply(fake_rate);
	c6->Multiply(fake_rate);

	cout<<"FF ++++++++++++++++++++++++++++++++++++++   "<<c6->GetBinContent(1)<<"   +-   "<<c6->GetBinError(1)<<endl;


// -------------------------------------------------  
//
//
	TH1D * c8 = (TH1D *) fake_rate->Clone("c8");
	c8->Multiply(fake_rate);

	TH1D * c9 = (TH1D *) fake_rate->Clone("c9");

	c8->Multiply(clone_n0t);

	c9->Multiply(clone_n1t);

	c8->Add(c9,-1);
	c8->Add(n_2t);

	c8->Divide(p_clone);
        c8->Divide(p_clone);
	c8->Multiply(prompt_rate);
	c8->Multiply(prompt_rate);

	cout<<"PP ++++++++++++++++++++++++++++++++++++++   "<<c8->GetBinContent(1)<<"   +-   "<<c8->GetBinError(1)<<endl;


	//cout<<"FP+FF+PP ???   "<<c1->GetBinContent(1)/(f*p)+c6->GetBinContent(1)/(f*f)+c8->GetBinContent(1)/(p*p)<<endl;

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//
	cout<<"--- fr from WJets -----------------------------------------------------------------------------------"<<endl;
	met_0t_wjet->Add(met_1t_wjet,-1);
	met_0t_wjet->Add(met_2t_wjet);
	cout<<"   integral double fake: "<<met_0t_wjet->Integral()<<endl;
	cout<<"   integral double fake - including overflows: "<<met_0t_wjet->Integral(0,31)<<endl;
	
	cout<<"--- fr from QCD -----------------------------------------------------------------------------------"<<endl;
	met_0t_data->Add(met_1t_data,-1);
	met_0t_data->Add(met_2t_data);
	cout<<"   integral double fake: "<<met_0t_data->Integral()<<endl;
	cout<<"   integral double fake - including overflows: "<<met_0t_data->Integral(0,31)<<endl;

	THStack * myStack = new THStack("hs","");
	met_0t_wjet->SetLineColor(kBlue);
	met_0t_wjet->SetTitle("Fake Method (FR from wjets - 0.13)");
	myStack->Add(met_0t_wjet);
	
	met_0t_data->SetLineColor(kGreen);
	met_0t_data->SetTitle("Fake Method (FR from qcd - 0.16)");
	myStack->Add(met_0t_data);
	
   TH1F *QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8 = new TH1F("QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8","",30,0,300);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(1,654.6697);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(2,1421.451);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(3,1388.635);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(4,927.6854);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(5,505.3523);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(6,228.0307);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(7,95.2132);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(8,39.18341);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(9,16.34255);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(10,6.327331);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(11,2.12365);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(12,1.677342);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(13,0.5009792);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(14,0.78095);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(15,0.5511512);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(16,0.2312188);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(17,0.08023093);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(18,0.4843155);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(19,-0.1504541);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(20,0.1584545);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(21,0.3056902);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(22,0.02617667);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(23,0.005154244);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(24,0.4075724);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(25,0.3197701);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(26,0.20474);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(27,0.002132303);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(28,0.2811694);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(29,-0.090826);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(30,0.02573229);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinContent(31,0.2839665);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(1,6.929691);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(2,10.24464);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(3,10.17104);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(4,8.377146);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(5,6.247657);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(6,4.329233);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(7,2.881247);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(8,1.976311);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(9,1.375651);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(10,0.9800938);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(11,0.7767834);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(12,0.6219888);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(13,0.5284557);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(14,0.4929002);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(15,0.3668446);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(16,0.3359898);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(17,0.288951);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(18,0.3319417);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(19,0.2360256);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(20,0.2214898);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(21,0.2090605);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(22,0.1590649);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(23,0.1074392);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(24,0.2061217);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(25,0.2016908);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(26,0.1442044);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(27,0.1196818);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(28,0.1585169);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(29,0.06121306);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(30,0.07677593);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetBinError(31,0.2115993);

   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetLineColor(kRed);
   QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8->SetTitle("ABCD");
	
	myStack->Add(QCD_9edf23c6b4226cac4bed3656defe2bcf_ZZTo4L_pfmet_pt_ZZTo4L__8);
	TCanvas * canvas = new TCanvas("can","can",350,350,500,500);
	myStack->Draw("nostack");
	myStack->SetTitle("QCD shape");
	myStack->GetXaxis()->SetTitle("pf met");
	myStack->GetYaxis()->SetTitle("Events");
	canvas->BuildLegend();
	canvas->SaveAs("QCD_comparison.png");
	canvas->SetLogy();
	canvas->SaveAs("QCD_comparison_log.png");
}
