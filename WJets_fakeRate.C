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
#include "treeBase.h"
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

double getTauIDWeight(double pt, double eta, double dm) {
    return 0.95;
}

double getMuToTauWeightLoose(double eta) {
    auto aeta = std::abs(eta);
    if (aeta < 0.4)
        return 1.22;
    if (aeta < 0.8)
        return 1.12;
    if (aeta < 1.2)
        return 1.26;
    if (aeta < 1.7)
        return 1.22;
    if (aeta < 2.3)
        return 2.39;
    return 1.;
}

double getEToTauWeightVLoose(double eta) {
    auto aeta = std::abs(eta);
    if (aeta < 1.46) 
        return 1.21;
    if (aeta < 0.8)
        return 1.12;
    return 1.;
}

double getTauWeight(int gen_match, double pt, double eta, double dm) {
    if (gen_match == 5)
        return getTauIDWeight(pt, eta, dm);
    if (gen_match == 2 || gen_match == 4)
        return getMuToTauWeightLoose(eta);
    if (gen_match == 1 || gen_match == 3)
        return getEToTauWeightVLoose(eta);
    return 1.;
}   

void WJets_fakeRate() {

	TH1::SetDefaultSumw2();

	double N_ll = 0.;
	double nevents = 0.;
	double nevents_loose = 0.;
	double nevents_tight = 0.;
	double nevents_loose_perEvt = 0.;
	double nevents_tight_perEvt = 0.;

	double N_0t = 0.;
	double N_1t = 0.;
	double N_2t = 0.;

	TH1D * n_0t = new TH1D("n_0t","n_0t",1,1.,2.);
	TH1D * n_1t = new TH1D("n_1t","n_1t",1,1.,2.);
	TH1D * n_2t = new TH1D("n_2t","n_2t",1,1.,2.);


	TH2D * h_gen_match = new TH2D("gen_match","gen_match",6,1.,7.,6,1.,7.);

	TH2D * h_gen_match_loose = new TH2D("gen_match_loose","gen_match_loose",13,-6.,7.,13,-6.,7.);
	TH2D * h_gen_match_tight = new TH2D("gen_match_tight","gen_match_tight",13,-6.,7.,13,-6.,7.);

	//TH2I * h_isolation_loose = new TH2I("isolation_loose","isolation_loose",14,0,7,14,0,7);
	//TH2I * h_isolation_tight = new TH2I("isolation_tight","isolation_tight",14,0,7,14,0,7);

	TH1D * n_singleFake = new TH1D("n_singleFake","n_singleFake",1,1.,2.);

	TH1D * n_loose = new TH1D("n_loose","n_loose",1,1.,2.);
	TH1D * n_tight = new TH1D("n_tight","n_tight",1,1.,2.);

	TH1D * h_pt_loose = new TH1D("pt_loose","pt_loose",13,40.,300.);
	TH1D * h_pt_tight = new TH1D("pt_tight","pt_tight",13,40.,300.);

	TH1D * h_eta_loose = new TH1D("eta_loose","eta_loose",10,-2.5,2.5);
	TH1D * h_eta_tight = new TH1D("eta_tight","eta_tight",10,-2.5,2.5);

	TH1D * h_dm_loose = new TH1D("dm_loose","dm_loose",12,0.,12.);
	TH1D * h_dm_tight = new TH1D("dm_tight","dm_tight",12,0.,12.);

	TH1D * h_jetpt_loose = new TH1D("jetpt_loose","jetpt_loose",15,0.,300.);
	TH1D * h_jetpt_tight = new TH1D("jetpt_tight","jetpt_tight",15,0.,300.);

	TH1D * h_mt2_loose = new TH1D("mt2_loose","mt2_loose",15,0.,300.);
	TH1D * h_mt2_tight = new TH1D("mt2_tight","mt2_tight",15,0.,300.);

	TH1D * h_summt_loose = new TH1D("summt_loose","summt_loose",20,0.,600.);
	TH1D * h_summt_tight = new TH1D("summt_tight","summt_tight",20,0.,600.);

	TH1D * h_njets_loose = new TH1D("njets_loose","njets_loose",10,0.,10.);
	TH1D * h_njets_tight = new TH1D("njets_tight","njets_tight",10,0.,10.);

	TString fileNames[4] = {"reproduction/W1JetsToLNu_LO/H2TauTauTreeProducerTauTau/tree.root","reproduction/W2JetsToLNu_LO/H2TauTauTreeProducerTauTau/tree.root","reproduction/W3JetsToLNu_LO/H2TauTauTreeProducerTauTau/tree.root","reproduction/W4JetsToLNu_LO/H2TauTauTreeProducerTauTau/tree.root"};
	//TString fileNames[4] = {"200317/DiTauMC/W1JetsToLNu_LO/H2TauTauTreeProducerTauTau/tree.root","200317/DiTauMC/W2JetsToLNu_LO/H2TauTauTreeProducerTauTau/tree.root","200317/DiTauMC/W3JetsToLNu_LO/H2TauTauTreeProducerTauTau/tree.root","200317/DiTauMC/W4JetsToLNu_LO/H2TauTauTreeProducerTauTau/tree.root"};

	double lumi = 35900.;
	double X_sec[4] = {9644.5,3144.5,954.8,485.6};
	double nGenEvts[4] = {45367044.0,29878415.0,19798117.0,9170576.0};

	int gen_match = 6; // fake rate calculation
	double IsoCut = 1.5; // changed from 0.5 to 1.5 to check the results from VLoose to Loose ! on May 3.  
	TString pathPNG = "PlotsFakeRate_WNJets/";
	mkdir(pathPNG,0777);
	//mkdir(pathPNG.c_str(),0777);

	for (int file = 0; file < 4; file++) {

	TFile* myFile = TFile::Open(fileNames[file]);
	TTree *myTree = (TTree*)(myFile->Get("tree"));

	treeBase tr(myTree);

	if (tr.fChain == 0) return;

	Long64_t nentries = tr.fChain->GetEntries();

	TString CurrentFile = "";

	//cout<<"Number of initial entries: "<<nentries<<endl;

	double xsec_weight = lumi * X_sec[file] / nGenEvts[file];

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        Long64_t ientry = tr.LoadTree(jentry);
        if (ientry < 0) break;
        tr.fChain->GetEntry(jentry);

	if (!(!tr.veto_dilepton && !tr.veto_thirdlepton && !tr.veto_otherlepton)) continue;
	if (!(tr.l1_againstMuon3>0.5 && tr.l1_againstElectronMVA6>0.5 && tr.l1_pt>40. && abs(tr.l1_eta)<2.1)) continue;
	if (!(tr.l2_againstMuon3>0.5 && tr.l2_againstElectronMVA6>0.5 && tr.l2_pt>40. && abs(tr.l2_eta)<2.1)) continue;
	if (!(tr.l1_charge != tr.l2_charge)) continue;
	
	//if (!(tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5)) continue;
	//if (!(tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5)) continue;
	
	if (!((tr.trigger_ditau35 && tr.trigger_matched_ditau35) || (tr.trigger_ditau35_combiso && tr.trigger_matched_ditau35_combiso))) continue;
	if (!(tr.Flag_HBHENoiseFilter && tr.Flag_HBHENoiseIsoFilter && tr.Flag_EcalDeadCellTriggerPrimitiveFilter && tr.Flag_goodVertices && tr.Flag_eeBadScFilter && tr.Flag_globalTightHalo2016Filter && tr.passBadMuonFilter && tr.passBadChargedHadronFilter && tr.badMuonMoriond2017 && tr.badCloneMuonMoriond2017)) continue;

	double evt_weight = tr.weight * xsec_weight; // * getTauWeight(tr.l1_gen_match,tr.l1_pt,tr.l1_eta,tr.l1_decayMode) * getTauWeight(tr.l2_gen_match,tr.l2_pt,tr.l2_eta,tr.l2_decayMode); // The scale factors are not applied for MC closure test.

	double summt = tr.pfmet_mt1+tr.pfmet_mt2;

	double l1_isLoose = tr.l1_byIsolationMVArun2v1DBoldDMwLT>IsoCut ? +1. : -1.;
	double l2_isLoose = tr.l2_byIsolationMVArun2v1DBoldDMwLT>IsoCut ? +1. : -1.;

	double l1_isTight = tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5 ? +1. : -1.;
	double l2_isTight = tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5 ? +1. : -1.;

	h_gen_match_loose->Fill(tr.l1_gen_match*l1_isLoose,tr.l2_gen_match*l2_isLoose,evt_weight);
	h_gen_match_tight->Fill(tr.l1_gen_match*l1_isTight,tr.l2_gen_match*l2_isTight,evt_weight);

	bool alreadyFilled_looseFake = false;
	bool alreadyFilled_tightFake = false;

	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l1_gen_match == gen_match) {
		n_loose->Fill(1.,evt_weight);
		nevents_loose += evt_weight;
		h_pt_loose->Fill(tr.l1_pt,evt_weight);
		h_eta_loose->Fill(tr.l1_eta,evt_weight);
		h_dm_loose->Fill(tr.l1_decayMode,evt_weight);
		h_jetpt_loose->Fill(tr.l1_jet_pt,evt_weight);
		h_mt2_loose->Fill(tr.mt2,evt_weight);
		h_summt_loose->Fill(summt,evt_weight);
		h_njets_loose->Fill(tr.n_jets,evt_weight);
		alreadyFilled_looseFake = true;
		nevents_loose_perEvt += evt_weight;
	}
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l2_gen_match == gen_match) {
		n_loose->Fill(1.,evt_weight);
		nevents_loose += evt_weight;
		h_pt_loose->Fill(tr.l2_pt,evt_weight);
		h_eta_loose->Fill(tr.l2_eta,evt_weight);
		h_dm_loose->Fill(tr.l2_decayMode,evt_weight);
                h_jetpt_loose->Fill(tr.l2_jet_pt,evt_weight);
		if (!alreadyFilled_looseFake) {
			nevents_loose_perEvt += evt_weight;
			h_mt2_loose->Fill(tr.mt2,evt_weight);
			h_summt_loose->Fill(summt,evt_weight);
			h_njets_loose->Fill(tr.n_jets,evt_weight);
		}
	}
	
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5 && tr.l1_gen_match == gen_match) {
		n_tight->Fill(1.,evt_weight);
		nevents_tight += evt_weight;
		h_pt_tight->Fill(tr.l1_pt,evt_weight);
		h_eta_tight->Fill(tr.l1_eta,evt_weight);
		h_dm_tight->Fill(tr.l1_decayMode,evt_weight);
		h_jetpt_tight->Fill(tr.l1_jet_pt,evt_weight);
		h_mt2_tight->Fill(tr.mt2,evt_weight);
		h_summt_tight->Fill(summt,evt_weight);
		h_njets_tight->Fill(tr.n_jets,evt_weight);
		alreadyFilled_tightFake = true;
		nevents_tight_perEvt += evt_weight;
	}
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5 && tr.l2_gen_match == gen_match) {
		n_tight->Fill(1.,evt_weight);
		nevents_tight += evt_weight;
		h_pt_tight->Fill(tr.l2_pt,evt_weight);
		h_eta_tight->Fill(tr.l2_eta,evt_weight);
		h_dm_tight->Fill(tr.l2_decayMode,evt_weight);
                h_jetpt_tight->Fill(tr.l2_jet_pt,evt_weight);
		if (!alreadyFilled_tightFake) {
			nevents_tight_perEvt += evt_weight;
			h_mt2_tight->Fill(tr.mt2,evt_weight);
			h_summt_tight->Fill(summt,evt_weight);
			h_njets_tight->Fill(tr.n_jets,evt_weight);
		}
	}
	
	if (!(tr.pfmet_pt > 30)) continue;
	if (!(tr.mvis > 85 || tr.mvis < 55)) continue;
	if (!(tr.mt2 > 20)) continue;
	// +++ SR I +++
	//if (!(tr.mt2 > 90)) continue;
	// +++ SR II +++
	//if (!(tr.mt2 < 90 && summt > 250 && tr.n_bjets == 0)) continue;

	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l2_byIsolationMVArun2v1DBoldDMwLT>IsoCut) N_ll += evt_weight;

	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5) {
		nevents += evt_weight;
		h_gen_match->Fill(tr.l1_gen_match,tr.l2_gen_match,evt_weight);
		//if (tr.l1_gen_match+tr.l2_gen_match == 11) n_singleFake->Fill(1.,evt_weight);
		n_singleFake->Fill(1.,evt_weight);
	}

	int both_loose_non_tight = 0;
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l1_byIsolationMVArun2v1DBoldDMwLT<4.5) both_loose_non_tight++;
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l2_byIsolationMVArun2v1DBoldDMwLT<4.5) both_loose_non_tight++;

	if (both_loose_non_tight == 2) {
		N_0t += evt_weight;	
		n_0t->Fill(1.,evt_weight);
	}

	int one_loose_one_tight = 0;
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l1_byIsolationMVArun2v1DBoldDMwLT<4.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5) one_loose_one_tight++;
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l2_byIsolationMVArun2v1DBoldDMwLT<4.5 && tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5) one_loose_one_tight++;

	if (one_loose_one_tight == 1) {
		N_1t += evt_weight;	
                n_1t->Fill(1.,evt_weight);
	}

	int both_tight = 0;
	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5) both_tight++;
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5) both_tight++;

	if (both_tight == 2) {
		N_2t += evt_weight;	
                n_2t->Fill(1.,evt_weight);
	}

}
}
	cout<<"N_ll: "<<N_ll<<endl;
	cout<<"N_0t: "<<N_0t<<endl;
	cout<<"N_1t: "<<N_1t<<endl;
	cout<<"N_2t: "<<N_2t<<endl;
	cout<<"sum: "<<N_0t+N_1t+N_2t<<endl;

	TH1D * h_loose = new TH1D("h_loose","h_loose",1,0.,1.);
	TH1D * h_tight = new TH1D("h_tight","h_tight",1,0.,1.);
	//TH1D * h_check = new TH1D("h_c","h_c",1,1.,2.);
	//for (int i = 0 ; i < 9; i++) h_check->Fill(1);
	//cout<<"==================  "<<h_check->GetBinError(1)<<endl; // to check if the error assigned to a bin is sqrt(content), of course when there is not weight set 
	h_loose->SetBinContent(1,nevents_loose);
	h_loose->SetBinError(1,sqrt(nevents_loose));
	h_tight->SetBinContent(1,nevents_tight);
	h_tight->SetBinError(1,sqrt(nevents_tight));

	double err = sqrt(pow(h_loose->GetBinError(1)/h_loose->GetBinContent(1),2) + pow(h_tight->GetBinError(1)/h_tight->GetBinContent(1),2)) * (h_tight->GetBinContent(1)/h_loose->GetBinContent(1));	
	cout<<"?????????? "<<err<<endl;

	h_tight->Divide(h_loose);
	cout<<"*** FR: "<<h_tight->GetBinContent(1)<<" +- "<<h_tight->GetBinError(1)<<endl;

	//cout<<"*** yields: "<<nevents<<endl;
	cout<<"*** loose taus: "<<nevents_loose<<endl;
	cout<<"*** tight taus: "<<nevents_tight<<endl;
	cout<<"*** fake rate:  "<<nevents_tight/nevents_loose<<endl;
	cout<<"---------------------"<<endl;
	cout<<"*** n_loose->GetBinContent(1): "<<n_loose->GetBinContent(1)<<endl;
	cout<<"*** n_tight->GetBinContent(1): "<<n_tight->GetBinContent(1)<<endl;

	err = sqrt(pow(n_loose->GetBinError(1)/n_loose->GetBinContent(1),2) + pow(n_tight->GetBinError(1)/n_tight->GetBinContent(1),2)) * (n_tight->GetBinContent(1)/n_loose->GetBinContent(1));	
	cout<<"?????????? "<<err<<endl;

	n_tight->Divide(n_loose);
	cout<<"*** FR: "<<n_tight->GetBinContent(1)<<" +- "<<n_tight->GetBinError(1)<<endl;

	cout<<"))))))))))))))))))))))))  single fakes(MC): "<<nevents<<endl;
        int bin56 = h_gen_match->FindBin(5,6);
        int bin65 = h_gen_match->FindBin(6,5);
        cout<<"		bin56: "<<h_gen_match->GetBinContent(5,6)<<endl;
        cout<<"		bin65: "<<h_gen_match->GetBinContent(6,5)<<endl;
        cout<<"		sum: "<<h_gen_match->GetBinContent(5,6)+h_gen_match->GetBinContent(6,5)<<endl;
        cout<<"		n_singleFake: "<<n_singleFake->GetBinContent(1)<<" +- "<<n_singleFake->GetBinError(1)<<endl;

	TH1D * prompt_rate = new TH1D("pr","pr",1,1.,2.);
	prompt_rate->SetBinContent(1,0.710724);
	prompt_rate->SetBinError(1,0.0083010);
	//prompt_rate->SetBinContent(1,0.639389);
	//prompt_rate->SetBinError(1,0.00479055); // these are for VLoose/VTight
	double p = prompt_rate->GetBinContent(1);
	double f = n_tight->GetBinContent(1); // *** FR: 0.28725 +- 0.0221622
	double fake_estimation = p*f*(-2.*f*p*N_0t + (f*(1-p)+p*(1-f))*N_1t -2.*(1-p)*(1-f)*N_2t)/pow(p-f,2);
	cout<<"))))))))))))))))))))))))  fake estimation: "<<fake_estimation<<endl;
	
	TH1D * p_clone = (TH1D *) prompt_rate->Clone("p_clone");
	p_clone->Add(n_tight,-1);

	n_0t->Multiply(prompt_rate);
	n_0t->Multiply(prompt_rate);
	n_0t->Multiply(n_tight);
	n_0t->Multiply(n_tight);
	n_0t->Divide(p_clone);
	n_0t->Divide(p_clone);

	TH1D * c1 = (TH1D *) prompt_rate->Clone("c1");
	c1->Multiply(n_tight);
	c1->Multiply(n_tight);

	TH1D * c2 = (TH1D *) prompt_rate->Clone("c2");
	c2->Multiply(n_tight);
	c2->Multiply(n_tight);
	c2->Multiply(prompt_rate);

	TH1D * c3 = (TH1D *) prompt_rate->Clone("c3");
	c3->Multiply(n_tight);
	c3->Multiply(prompt_rate);

	c1->Add(c3);	
	c1->Add(c2,-2);	
	n_1t->Multiply(c1);
	n_1t->Divide(p_clone);
	n_1t->Divide(p_clone);

	TH1D * c4 = (TH1D *) prompt_rate->Clone("c4");
	c4->Multiply(n_tight);

	TH1D * c5 = (TH1D *) prompt_rate->Clone("c5");
	c5->Multiply(n_tight);
	c5->Multiply(n_tight);

	c4->Add(c5,-1);
	c4->Add(c3,-1);
	c4->Add(c2);
	n_2t->Multiply(c4);
	n_2t->Divide(p_clone);
	n_2t->Divide(p_clone);

	n_1t->Add(n_0t,-2);
	n_1t->Add(n_2t,-2);
	cout<<"++++++++++++++++++++++++++++++++++++++   "<<n_1t->GetBinContent(1)<<"   +-   "<<n_1t->GetBinError(1)<<endl;

////
	CanvasCreator(h_gen_match,"l1_gen_match","l2_gen_match","genmatch",pathPNG);		

	CanvasCreator(h_gen_match_loose,"l1_gen_match","l2_gen_match","negative: non-loose; positive: loose",pathPNG);		
	CanvasCreator(h_gen_match_tight,"l1_gen_match","l2_gen_match","negative: non-tight; positive: tight",pathPNG);		
////

	CanvasCreator(h_pt_loose,"p_{T}","Loose fake #taus",pathPNG);		
	CanvasCreator(h_pt_tight,"p_{T}","Tight fake #taus",pathPNG);		

	TH1D * h_FR_pt = (TH1D*) h_pt_tight->Clone("h_FR_pt");
	h_FR_pt->Divide(h_pt_loose);
	CanvasCreator(h_FR_pt,"p_{T}","Fake Rate",pathPNG);		
////

        CanvasCreator(h_eta_loose,"#eta","Loose fake #taus",pathPNG);     
        CanvasCreator(h_eta_tight,"#eta","Tight fake #taus",pathPNG);     

	TH1D * h_FR_eta = (TH1D*) h_eta_tight->Clone("h_FR_eta");
	h_FR_eta->Divide(h_eta_loose);
	CanvasCreator(h_FR_eta,"#eta","Fake Rate",pathPNG);		
////

        CanvasCreator(h_dm_loose,"decay mode","Loose fake #taus",pathPNG);     
        CanvasCreator(h_dm_tight,"decay mode","Tight fake #taus",pathPNG);     

	TH1D * h_FR_dm = (TH1D*) h_dm_tight->Clone("h_FR_dm");
	h_FR_dm->Divide(h_dm_loose);
	CanvasCreator(h_FR_dm,"decay mode","Fake Rate",pathPNG);		
////

	CanvasCreator(h_jetpt_loose,"associated jet p_{T}","Loose fake #taus",pathPNG);		
	CanvasCreator(h_jetpt_tight,"associated jet p_{T}","Tight fake #taus",pathPNG);		

	TH1D * h_FR_jetpt = (TH1D*) h_jetpt_tight->Clone("h_FR_jetpt");
	h_FR_jetpt->Divide(h_jetpt_loose);
	CanvasCreator(h_FR_jetpt,"associated jet p_{T}","Fake Rate",pathPNG);		
////

        CanvasCreator(h_mt2_loose,"M_{T2}","Events with loose fake #taus",pathPNG);     
        CanvasCreator(h_mt2_tight,"M_{T2}","Events with tight fake #taus",pathPNG);     

	TH1D * h_FR_mt2 = (TH1D*) h_mt2_tight->Clone("h_FR_mt2");
	h_FR_mt2->Divide(h_mt2_loose);
	CanvasCreator(h_FR_mt2,"M_{T2}","Fake Rate",pathPNG);		
////

        CanvasCreator(h_summt_loose,"M_{T}^{#tau1}+M_{T}^{#tau2}","Events with loose fake #taus",pathPNG);     
        CanvasCreator(h_summt_tight,"M_{T}^{#tau1}+M_{T}^{#tau2}","Events with tight fake #taus",pathPNG);     

	TH1D * h_FR_summt = (TH1D*) h_summt_tight->Clone("h_FR_summt");
	h_FR_summt->Divide(h_summt_loose);
	CanvasCreator(h_FR_summt,"M_{T}^{#tau1}+M_{T}^{#tau2}","Fake Rate",pathPNG);		
////

        CanvasCreator(h_njets_loose,"N_{jets}","Events with loose fake #taus",pathPNG);     
        CanvasCreator(h_njets_tight,"N_{jets}","Events with tight fake #taus",pathPNG);     

	TH1D * h_FR_njets = (TH1D*) h_njets_tight->Clone("h_FR_njets");
	h_FR_njets->Divide(h_njets_loose);
	CanvasCreator(h_FR_njets,"N_{jets}","Fake Rate",pathPNG);		
////

}
/*
)))))) hist.name:  W1Jets
)))))) hist.integral weighted:  112.743965149
)))))) hist.integral unweighted:  36.0663386583
         -->> entries found:  73.0
         w:  112.743965149  and error:  14.123093446
         u:  36.0663386583  and error:  4.51792050052
)))))) hist.name:  W2Jets
)))))) hist.integral weighted:  83.8979568481
)))))) hist.integral unweighted:  61.9288195968
         -->> entries found:  114.0
         w:  83.8979568481  and error:  8.40984408492
         u:  61.9288195968  and error:  6.20768063328
)))))) hist.name:  W3Jets
)))))) hist.integral weighted:  31.9809514284
)))))) hist.integral unweighted:  51.5176604986
         -->> entries found:  89.0
         w:  31.9809514284  and error:  3.59562746388
         u:  51.5176604986  and error:  5.79214514715
)))))) hist.name:  W4Jets
)))))) hist.integral weighted:  17.2988674045
)))))) hist.integral unweighted:  25.3007566929
         -->> entries found:  47.0
         w:  17.2988674045  and error:  2.65301851042
         u:  25.3007566929  and error:  3.88021789974
*/
