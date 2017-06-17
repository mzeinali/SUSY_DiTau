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

void QCD_fakeRate() {

	TH1::SetDefaultSumw2();


	TH1D * n_loose = new TH1D("n_loose","n_loose",1,1.,2.);
	TH1D * n_tight = new TH1D("n_tight","n_tight",1,1.,2.);

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

	TH1D * prompt_rate = new TH1D("pr","pr",1,1.,2.);
	prompt_rate->SetBinContent(1,0.639389);
	prompt_rate->SetBinError(1,0.00479055);

	TH1D * fake_rate_wjet = new TH1D("fr_wjet","fr_wjet",1,1.,2.);
	fake_rate_wjet->SetBinContent(1,0.13478);
	fake_rate_wjet->SetBinError(1,0.00656193); // from W+jets

	TH1D * fake_rate_data = new TH1D("fr_data","fr_data",1,1.,2.);
	fake_rate_data->SetBinContent(1,0.161887);
	fake_rate_data->SetBinError(1,0.00059797); // from Data

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

	//cout<<"Number of initial entries: "<<nentries<<endl;

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        Long64_t ientry = tr.LoadTree(jentry);
        if (ientry < 0) break;
        tr.fChain->GetEntry(jentry);

	if (!(!tr.veto_dilepton && !tr.veto_thirdlepton && !tr.veto_otherlepton)) continue;
	if (!(tr.l1_againstMuon3>0.5 && tr.l1_againstElectronMVA6>0.5 && tr.l1_pt>40. && abs(tr.l1_eta)<2.1)) continue;
	if (!(tr.l2_againstMuon3>0.5 && tr.l2_againstElectronMVA6>0.5 && tr.l2_pt>40. && abs(tr.l2_eta)<2.1)) continue;

	//if (!(tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5)) continue;
	//if (!(tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5)) continue;

	if (!((tr.trigger_ditau35 && tr.trigger_matched_ditau35) || (tr.trigger_ditau35_combiso && tr.trigger_matched_ditau35_combiso))) continue;
	if (!(tr.Flag_HBHENoiseFilter && tr.Flag_HBHENoiseIsoFilter && tr.Flag_EcalDeadCellTriggerPrimitiveFilter && tr.Flag_goodVertices && tr.Flag_eeBadScFilter && tr.Flag_globalTightHalo2016Filter && tr.passBadMuonFilter && tr.passBadChargedHadronFilter && tr.badMuonMoriond2017 && tr.badCloneMuonMoriond2017)) continue;

	//if (!(tr.l1_charge != tr.l2_charge)) continue; // check Yields
	if (!(tr.l1_charge == tr.l2_charge)) continue; // SS requirement, to make a pure QCD sample from data
	if (!(tr.pfmet_pt < 30)) continue; // pfmet requirement, to make a pure QCD sample from data, these 2 cuts result in a purity of 98%(94%) for loose(tight) tau selection
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>1.5 || tr.l2_byIsolationMVArun2v1DBoldDMwLT>1.5) n_loose->Fill(1.);
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5 || tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5) n_tight->Fill(1.);
	
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>1.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT>1.5) N_ll++;

	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5) nevents++;

	double weight_0t_wjet = f_wjet*f_wjet*p*p/pow(p-f_wjet,2);
	double weight_1t_wjet = f_wjet*f_wjet*p*(1-p)/pow(p-f_wjet,2);
	double weight_2t_wjet = f_wjet*f_wjet*(1-p)*(1-p)/pow(p-f_wjet,2);

	double weight_0t_data = f_data*f_data*p*p/pow(p-f_data,2);
	double weight_1t_data = f_data*f_data*p*(1-p)/pow(p-f_data,2);
	double weight_2t_data = f_data*f_data*(1-p)*(1-p)/pow(p-f_data,2);

	int both_loose_non_tight = 0;
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>1.5 && tr.l1_byIsolationMVArun2v1DBoldDMwLT<4.5) both_loose_non_tight++;
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>1.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT<4.5) both_loose_non_tight++;

	if (both_loose_non_tight == 2) {
		N_0t++;	
		n_0t->Fill(1.);
		met_0t_wjet->Fill(tr.pfmet_pt,weight_0t_wjet);
		met_0t_data->Fill(tr.pfmet_pt,weight_0t_data);
	}

	int one_loose_one_tight = 0;
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>1.5 && tr.l1_byIsolationMVArun2v1DBoldDMwLT<4.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5) one_loose_one_tight++;
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>1.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT<4.5 && tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5) one_loose_one_tight++;

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

	cout<<"*** n_loose->GetBinContent(1): "<<n_loose->GetBinContent(1)<<endl;
	cout<<"*** n_tight->GetBinContent(1): "<<n_tight->GetBinContent(1)<<endl;

	err = sqrt(pow(n_loose->GetBinError(1)/n_loose->GetBinContent(1),2) + pow(n_tight->GetBinError(1)/n_tight->GetBinContent(1),2)) * (n_tight->GetBinContent(1)/n_loose->GetBinContent(1));	
	cout<<"?????????? "<<err<<endl;

	n_tight->Divide(n_loose);
	cout<<"*** FR: "<<n_tight->GetBinContent(1)<<" +- "<<n_tight->GetBinError(1)<<endl;

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

	c6->Add(c1);
	cout<<"SF + FF +++++++++++++++++++++++++++++++++   "<<c6->GetBinContent(1)<<"   +-   "<<c6->GetBinError(1)<<endl;


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

}
