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

void CanvasCreator(TH1D * h1, TH1D * h2, string m_name, string legend_name) {
	
	//gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);

	TCanvas * canvas = new TCanvas(m_name.c_str(),m_name.c_str(),632,112,500,502);
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

	h1->SetTitle(legend_name.c_str());
	h1->SetLineColor(2);
	h2->SetLineColor(4);
	h1->SetMaximum(20);
	h1->Draw();
	//h2->SetLineStyle(3);
        h2->Draw("same");

   TLegend * leg = new TLegend(0.6,0.5,0.9,0.7);
   leg->AddEntry(h1,"caloMET","f");
   leg->AddEntry(h2,"pfMET","f");
   leg->Draw();

	//canvas->SetLogy();
        canvas->SaveAs((m_name+legend_name+".png").c_str());
        canvas->SaveAs((m_name+legend_name+".C").c_str());
}

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

void CountSMS_SetLimit() {

	TFile * f = new TFile("/cmsdata2/mzeinali/SUSY_DiTau/CMSSW_8_0_13/src/SMS_TChipmStauSnu/H2TauTauTreeProducerTauTau/tree.root");
	TH3D * myWeights = (TH3D*)f->Get("SumGenWeightsSMS");

	TFile* myFile = TFile::Open("200317/DiTauMC/SMS_TChipmStauSnu/H2TauTauTreeProducerTauTau/tree.root");
	TTree *myTree = (TTree*)(myFile->Get("tree"));

	treeBase tr(myTree);

	if (tr.fChain == 0) return;

	Long64_t nentries = tr.fChain->GetEntries();

	TString CurrentFile = "";

	double nevents = 0;
	cout<<"Number of initial entries: "<<nentries<<endl;

	bool sr1 = false;
	bool sr2 = true;
	bool sr3 = false;
	
	TH1D * h_xsec = new TH1D("h_xsec","h_xsec",1000,0.,1000.);
	TH1D * h_PN_Bkg = new TH1D("h_PN_Bkg","h_PN_Bkg",1,0.,1.);
	// set the backgrounds here
	//if (sr1) h_PN_Bkg->SetBinContent(1,10.7532); // 10.7532 +- 5.36486
	//if (sr1) h_PN_Bkg->SetBinContent(1,4.83252); // 4.83252 +- 4.04178 , when l1_pt > 100. GeV applied in SR I.
	if (sr1) h_PN_Bkg->SetBinContent(1,5.92063); // 5.92063 +- 3.52785 , when l1_pt < 100. GeV applied in SR I.
	//if (sr2) h_PN_Bkg->SetBinContent(1,58.11); // 58.11 +- 9.20804
	//if (sr2) h_PN_Bkg->SetBinContent(1,39.4648); // 39.4648 +- 7.39612 , when dphi_tau1_tau2 > 2. added
	if (sr2) h_PN_Bkg->SetBinContent(1,14.2451); // 14.2451 +- 4.43295 , when sumMT cut changed to 300 GeV
	if (sr3) h_PN_Bkg->SetBinContent(1,3.);  

	//if (sr1) h_PN_Bkg->SetBinError(1,5.36486);  
	//if (sr1) h_PN_Bkg->SetBinError(1,4.04178);  
	if (sr1) h_PN_Bkg->SetBinError(1,3.52785);  
	//if (sr2) h_PN_Bkg->SetBinError(1,9.20804);  
	//if (sr2) h_PN_Bkg->SetBinError(1,7.39612);  
	if (sr2) h_PN_Bkg->SetBinError(1,4.43295);  
	if (sr3) h_PN_Bkg->SetBinError(1,0.1);  

	TH2D * h_Counting = new TH2D("h_Counting","h_Counting",1000,0.,1000.,1000,0.,1000.);
	TH2D * h_PN_MLSP_MChi = new TH2D("h_PN_MLSP_MChi","h_PN_MLSP_MChi",1000,0.,1000.,1000,0.,1000.);
	TH2D * C1C1_13TeV_NLONLL_LSP = new TH2D("C1C1_13TeV_NLONLL_LSP","C1C1_13TeV_NLONLL_LSP",1000,0.,1000.,1000,0.,1000.);
	TH2D * UP_C1C1_13TeV_NLONLL_LSP = new TH2D("UP_C1C1_13TeV_NLONLL_LSP","UP_C1C1_13TeV_NLONLL_LSP",1000,0.,1000.,1000,0.,1000.);
	TH2D * DOWN_C1C1_13TeV_NLONLL_LSP = new TH2D("DOWN_C1C1_13TeV_NLONLL_LSP","DOWN_C1C1_13TeV_NLONLL_LSP",1000,0.,1000.,1000,0.,1000.);

	ifstream fin("xsec_C1C1_wino.txt");

	double chiMass;
	double xsection;
	double error;
	vector<std::pair <double, double > > XSections;
	vector<std::pair <double, double > > XSectionErrors;

	while (fin >> chiMass >> xsection >> error) {
		h_xsec->Fill(chiMass,xsection);
		XSections.push_back(std::make_pair(chiMass,xsection));
		XSectionErrors.push_back(std::make_pair(chiMass,error));
	}


    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        Long64_t ientry = tr.LoadTree(jentry);
        if (ientry < 0) break;
        tr.fChain->GetEntry(jentry);

	if (!(!tr.veto_dilepton && !tr.veto_thirdlepton && !tr.veto_otherlepton)) continue;
	if (!(tr.l1_againstMuon3>0.5 && tr.l1_againstElectronMVA6>0.5 && tr.l1_pt>40. && abs(tr.l1_eta)<2.1)) continue;
	if (!(tr.l2_againstMuon3>0.5 && tr.l2_againstElectronMVA6>0.5 && tr.l2_pt>40. && abs(tr.l2_eta)<2.1)) continue;
	if (!(tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5)) continue;
	if (!(tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5)) continue;
	if (!(tr.l1_charge != tr.l2_charge)) continue;

	if (!(tr.Flag_HBHENoiseFilter && tr.Flag_HBHENoiseIsoFilter && tr.Flag_EcalDeadCellTriggerPrimitiveFilter && tr.Flag_goodVertices && tr.Flag_eeBadScFilter && tr.Flag_globalTightHalo2016Filter && tr.passBadMuonFilter && tr.passBadChargedHadronFilter && tr.badMuonMoriond2017 && tr.badCloneMuonMoriond2017)) continue;


	if (!(tr.pfmet_pt > 30)) continue;
	if (!(tr.mvis > 85 || tr.mvis < 55)) continue;
	if (!(tr.mt2 > 20)) continue;

	// ++++ SR I +++++++++++
	//if (sr1 && !(tr.mt2 > 90 && tr.l1_pt > 100.)) continue;
	if (sr1 && !(tr.mt2 > 90 && tr.l1_pt < 100.)) continue;

	double diff = tr.l1_phi - tr.l2_phi;	
	if (diff >= M_PI) diff -= 2. * M_PI;
	if (diff <  (-1. * M_PI)) diff += 2. * M_PI;

	// ++++ SR II +++++++++++
	//if (sr2 && !(tr.mt2 < 90 && tr.n_bjets == 0 && (tr.pfmet_mt1+tr.pfmet_mt2) > 300 && fabs(diff) > 2.)) continue; 
	if (sr2 && !(tr.mt2 < 90 && tr.n_bjets == 0 && (tr.pfmet_mt1+tr.pfmet_mt2) > 300)) continue; // on June 13, 250 chaned with 300 to check the effect of cutting on higher sumMT cuts
	
	// ++++ SR III +++++++++++
	if (sr3 && !(tr.mt2 < 90 && tr.n_bjets == 0 && (tr.pfmet_mt1+tr.pfmet_mt2) < 250 && (tr.pfmet_mt1+tr.pfmet_mt2) > 200)) continue;
	
	//if (!(tr.GenSusyMChargino == 400 && tr.GenSusyMNeutralino == 1)) continue;
	//if (!(tr.GenSusyMChargino == 200 && tr.GenSusyMNeutralino == 1)) continue;
	

	double weight_total = tr.weight * getTauWeight(tr.l1_gen_match,tr.l1_pt,tr.l1_eta,tr.l1_decayMode) * getTauWeight(tr.l2_gen_match,tr.l2_pt,tr.l2_eta,tr.l2_decayMode);
	nevents += weight_total;

	//if (tr.GenSusyMChargino != 0 || tr.GenSusyMNeutralino != 0)
		h_Counting->Fill(tr.GenSusyMChargino,tr.GenSusyMNeutralino,weight_total);
}

	for (int i = 1; i <= h_Counting->GetXaxis()->GetNbins(); i++) {
		double xsec = 0.;		
		double err = 0.;
		for (int s = 0; s < XSections.size(); s++) {
			if (XSections[s].first == i) {
				xsec = XSections[s].second;
				err = XSectionErrors[s].second;
			}
		}
		for (int j = 1; j <= h_Counting->GetYaxis()->GetNbins(); j++) {

			int global = myWeights->FindBin(i,j);
			double sumWeight = myWeights->GetBinContent(global);

			if (sumWeight== 0) continue; 	
			//if (i == 200 && j == 1) {
			//cout<<""<<endl;
			//cout<<">>> myWeights->FindBin(i,j): "<<global<<endl;
			//cout<<">>> myWeights->GetBinContent(global): "<<sumWeight<<endl;
			//cout<<">>> xsec: "<<xsec<<endl;

			global = C1C1_13TeV_NLONLL_LSP->FindBin(i,j);
			C1C1_13TeV_NLONLL_LSP->SetBinContent(global,xsec/1000.);		
			UP_C1C1_13TeV_NLONLL_LSP->SetBinContent(global,(xsec+err)/1000.);		
			DOWN_C1C1_13TeV_NLONLL_LSP->SetBinContent(global,(xsec-err)/1000.);		

		        global = h_Counting->FindBin(i,j);
			double content = h_Counting->GetBinContent(global);
			//cout<<">>> h_Counting->FindBin(i,j): "<<global<<endl;
			//cout<<">>> h_Counting->GetBinContent(global): "<<content<<endl;

			double yield = content * xsec * 35900. / (sumWeight * 1000.);

			h_PN_MLSP_MChi->SetBinContent(global,yield);
			//}	

		}
	}


	/*
	cout<<"Number of unweighted events (400,1): "<<nevents<<endl;

	int global = h_Counting->FindBin(200,1);
	double content = h_Counting->GetBinContent(global);

	cout<<"**** unweighted (200,1): "<<content<<endl;	
	*/
	int global = h_PN_MLSP_MChi->FindBin(200,100);
	double content = h_PN_MLSP_MChi->GetBinContent(global);

	cout<<"**** (200,100): "<<content<<endl;	

	global = h_PN_MLSP_MChi->FindBin(250,50);
	content = h_PN_MLSP_MChi->GetBinContent(global);

	cout<<"**** (250,50): "<<content<<endl;	

	global = h_PN_MLSP_MChi->FindBin(400,1);
	content = h_PN_MLSP_MChi->GetBinContent(global);

	cout<<"**** (400,1): "<<content<<endl;	
/*
	global = UP_C1C1_13TeV_NLONLL_LSP->FindBin(400,200);
	content = UP_C1C1_13TeV_NLONLL_LSP->GetBinContent(global);
	cout<<"**** UP xsec (400,200): "<<content<<endl;	

	global = DOWN_C1C1_13TeV_NLONLL_LSP->FindBin(400,200);
	content = DOWN_C1C1_13TeV_NLONLL_LSP->GetBinContent(global);
	cout<<"**** DOWN xsec (400,200): "<<content<<endl;	
*/
	TString postFix = "";
	if (sr1) postFix = "_SRI_tauptLT100";
	//else if (sr2) postFix = "_SRII_dphiAdded";
	else if (sr2) postFix = "_SRII_sumMT300";
	else if (sr3) postFix = "_SRIII";

	TString filename = "Yields";
	filename.Append(postFix);
	filename.Append(".root");

	TFile * f_yields = new TFile(filename,"recreate");

	h_xsec->Write();
	h_Counting->Write();
	h_PN_Bkg->Write();
	h_PN_MLSP_MChi->Write();

	f_yields->Close();

// this is something which is identical for various Signal Regions
/*
	filename = "referenceXSecs";
	filename.Append(postFix);
	filename.Append(".root");

	TFile * f_xsec = new TFile(filename,"recreate");

	C1C1_13TeV_NLONLL_LSP->Write();
	UP_C1C1_13TeV_NLONLL_LSP->Write();
	DOWN_C1C1_13TeV_NLONLL_LSP->Write();

	f_xsec->Close();
*/

}
