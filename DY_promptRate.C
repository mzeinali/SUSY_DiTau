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

void DY_promptRate() {

	TH1::SetDefaultSumw2();

	double nevents = 0.;
	double nevents_loose = 0.;
	double nevents_tight = 0.;
	double nevents_loose_perEvt = 0.;
	double nevents_tight_perEvt = 0.;

	TH2D * h_gen_match = new TH2D("gen_match","gen_match",7,0.,7.,7,0.,7.);

	TH2D * h_gen_match_loose = new TH2D("gen_match_loose","gen_match_loose",13,-6.,7.,13,-6.,7.);
	TH2D * h_gen_match_tight = new TH2D("gen_match_tight","gen_match_tight",13,-6.,7.,13,-6.,7.);

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

	TString fileName = "200317/DiTauMC/DYJetsToLL_M50_LO_ext2/H2TauTauTreeProducerTauTau/tree.root";
	//TString fileName = "200317/DiTauMC/DYJetsToLL_M50_LO_ext/H2TauTauTreeProducerTauTau/tree.root";

	double lumi = 35900.;
	double X_sec = 5765.4;
	double nGenEvts = 96658943.0 ;// for ext : 49144274.0; // for ext2 : 96658943.0;

	int gen_match = 5; // prompt rate calculation
	double IsoCut = 1.5; // changed on May 3 from 0.5 to 1.5 (VLoose to Loose)
	TString pathPNG = "PlotsPromptRate_DYJets/";
	mkdir(pathPNG,0777);

	TFile* myFile = TFile::Open(fileName);
	TTree *myTree = (TTree*)(myFile->Get("tree"));

	treeBase tr(myTree);

	if (tr.fChain == 0) return;

	Long64_t nentries = tr.fChain->GetEntries();

	TString CurrentFile = "";

	//cout<<"Number of initial entries: "<<nentries<<endl;

	double xsec_weight = lumi * X_sec / nGenEvts;

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
	//if (!(tr.l1_gen_match == 5 && tr.l2_gen_match == 5)) continue;
	
	double evt_weight = tr.weight * xsec_weight;
	nevents += evt_weight;

	double summt = tr.pfmet_mt1+tr.pfmet_mt2;

	double l1_isLoose = tr.l1_byIsolationMVArun2v1DBoldDMwLT>IsoCut ? +1. : -1.;
	double l2_isLoose = tr.l2_byIsolationMVArun2v1DBoldDMwLT>IsoCut ? +1. : -1.;

	double l1_isTight = tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5 ? +1. : -1.;
	double l2_isTight = tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5 ? +1. : -1.;

	h_gen_match_loose->Fill(tr.l1_gen_match*l1_isLoose,tr.l2_gen_match*l2_isLoose,evt_weight);
	h_gen_match_tight->Fill(tr.l1_gen_match*l1_isTight,tr.l2_gen_match*l2_isTight,evt_weight);

	bool alreadyFilled_loosePrompt = false;
	bool alreadyFilled_tightPrompt = false;

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
		alreadyFilled_loosePrompt = true;
		nevents_loose_perEvt += evt_weight;
	}
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>IsoCut && tr.l2_gen_match == gen_match) {
		n_loose->Fill(1.,evt_weight);
		nevents_loose += evt_weight;
		h_pt_loose->Fill(tr.l2_pt,evt_weight);
		h_eta_loose->Fill(tr.l2_eta,evt_weight);
		h_dm_loose->Fill(tr.l2_decayMode,evt_weight);
                h_jetpt_loose->Fill(tr.l2_jet_pt,evt_weight);
		if (!alreadyFilled_loosePrompt) {
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
		alreadyFilled_tightPrompt = true;
		nevents_tight_perEvt += evt_weight;
	}
	if (tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5 && tr.l2_gen_match == gen_match) {
		n_tight->Fill(1.,evt_weight);
		nevents_tight += evt_weight;
		h_pt_tight->Fill(tr.l2_pt,evt_weight);
		h_eta_tight->Fill(tr.l2_eta,evt_weight);
		h_dm_tight->Fill(tr.l2_decayMode,evt_weight);
                h_jetpt_tight->Fill(tr.l2_jet_pt,evt_weight);
		if (!alreadyFilled_tightPrompt) {
			nevents_tight_perEvt += evt_weight;
			h_mt2_tight->Fill(tr.mt2,evt_weight);
			h_summt_tight->Fill(summt,evt_weight);
			h_njets_tight->Fill(tr.n_jets,evt_weight);
		}
	}
	

	if (tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5 && tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5) {
		h_gen_match->Fill(tr.l1_gen_match,tr.l2_gen_match,evt_weight);
	}
}
cout<<"+++++++ nevents_loose "<<nevents_loose<<endl;
cout<<"+++++++ nevents_tight "<<nevents_tight<<endl;
	TH1D * h_loose = new TH1D("h_loose","h_loose",1,0.,1.);
	TH1D * h_tight = new TH1D("h_tight","h_tight",1,0.,1.);

	h_loose->SetBinContent(1,nevents_loose);
	h_loose->SetBinError(1,sqrt(nevents_loose));
	h_tight->SetBinContent(1,nevents_tight);
	h_tight->SetBinError(1,sqrt(nevents_tight));

	double err = sqrt(pow(h_loose->GetBinError(1)/h_loose->GetBinContent(1),2) + pow(h_tight->GetBinError(1)/h_tight->GetBinContent(1),2)) * (h_tight->GetBinContent(1)/h_loose->GetBinContent(1));	
	cout<<"?????????? "<<err<<endl;

	h_tight->Divide(h_loose);
	cout<<"*** PR: "<<h_tight->GetBinContent(1)<<" +- "<<h_tight->GetBinError(1)<<endl;
	cout<<"---------------------"<<endl;
	cout<<"*** n_loose->GetBinContent(1): "<<n_loose->GetBinContent(1)<<endl;
	cout<<"*** n_tight->GetBinContent(1): "<<n_tight->GetBinContent(1)<<endl;

	err = sqrt(pow(n_loose->GetBinError(1)/n_loose->GetBinContent(1),2) + pow(n_tight->GetBinError(1)/n_tight->GetBinContent(1),2)) * (n_tight->GetBinContent(1)/n_loose->GetBinContent(1));	
	cout<<"?????????? "<<err<<endl;

	n_tight->Divide(n_loose);
	cout<<"*** PR: "<<n_tight->GetBinContent(1)<<" +- "<<n_tight->GetBinError(1)<<endl;

	cout<<"*** yields: "<<nevents<<endl;

////
	CanvasCreator(h_gen_match,"l1_gen_match","l2_gen_match","genmatch",pathPNG);		

	CanvasCreator(h_gen_match_loose,"l1_gen_match","l2_gen_match","negative: non-loose; positive: loose",pathPNG);		
	CanvasCreator(h_gen_match_tight,"l1_gen_match","l2_gen_match","negative: non-tight; positive: tight",pathPNG);		
////

	CanvasCreator(h_pt_loose,"p_{T}","Loose prompt #taus",pathPNG);		
	CanvasCreator(h_pt_tight,"p_{T}","Tight prompt #taus",pathPNG);		

	TH1D * h_PR_pt = (TH1D*) h_pt_tight->Clone("h_PR_pt");
	h_PR_pt->Divide(h_pt_loose);
	CanvasCreator(h_PR_pt,"p_{T}","Prompt Rate",pathPNG);		
////

        CanvasCreator(h_eta_loose,"#eta","Loose prompt #taus",pathPNG);     
        CanvasCreator(h_eta_tight,"#eta","Tight prompt #taus",pathPNG);     

	TH1D * h_PR_eta = (TH1D*) h_eta_tight->Clone("h_PR_eta");
	h_PR_eta->Divide(h_eta_loose);
	CanvasCreator(h_PR_eta,"#eta","Prompt Rate",pathPNG);		
////

        CanvasCreator(h_dm_loose,"decay mode","Loose prompt #taus",pathPNG);     
        CanvasCreator(h_dm_tight,"decay mode","Tight prompt #taus",pathPNG);     

	TH1D * h_PR_dm = (TH1D*) h_dm_tight->Clone("h_PR_dm");
	h_PR_dm->Divide(h_dm_loose);
	CanvasCreator(h_PR_dm,"decay mode","Prompt Rate",pathPNG);		
////

	CanvasCreator(h_jetpt_loose,"associated jet p_{T}","Loose prompt #taus",pathPNG);		
	CanvasCreator(h_jetpt_tight,"associated jet p_{T}","Tight prompt #taus",pathPNG);		

	TH1D * h_PR_jetpt = (TH1D*) h_jetpt_tight->Clone("h_PR_jetpt");
	h_PR_jetpt->Divide(h_jetpt_loose);
	CanvasCreator(h_PR_jetpt,"associated jet p_{T}","Prompt Rate",pathPNG);		
////

        CanvasCreator(h_mt2_loose,"M_{T2}","Events with loose prompt #taus",pathPNG);     
        CanvasCreator(h_mt2_tight,"M_{T2}","Events with tight prompt #taus",pathPNG);     

	TH1D * h_PR_mt2 = (TH1D*) h_mt2_tight->Clone("h_PR_mt2");
	h_PR_mt2->Divide(h_mt2_loose);
	CanvasCreator(h_PR_mt2,"M_{T2}","Prompt Rate",pathPNG);		
////

        CanvasCreator(h_summt_loose,"M_{T}^{#tau1}+M_{T}^{#tau2}","Events with loose prompt #taus",pathPNG);     
        CanvasCreator(h_summt_tight,"M_{T}^{#tau1}+M_{T}^{#tau2}","Events with tight prompt #taus",pathPNG);     

	TH1D * h_PR_summt = (TH1D*) h_summt_tight->Clone("h_PR_summt");
	h_PR_summt->Divide(h_summt_loose);
	CanvasCreator(h_PR_summt,"M_{T}^{#tau1}+M_{T}^{#tau2}","Prompt Rate",pathPNG);		
////

        CanvasCreator(h_njets_loose,"N_{jets}","Events with loose prompt #taus",pathPNG);     
        CanvasCreator(h_njets_tight,"N_{jets}","Events with tight prompt #taus",pathPNG);     

	TH1D * h_PR_njets = (TH1D*) h_njets_tight->Clone("h_PR_njets");
	h_PR_njets->Divide(h_njets_loose);
	CanvasCreator(h_PR_njets,"N_{jets}","Prompt Rate",pathPNG);		
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
