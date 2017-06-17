#include <memory>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TH2.h>
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

void SMS_Counting() {

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

	TH2D * h_Counting = new TH2D("h_Counting","h_Counting",1000,0.,1000.,1000,0.,1000.);
	TH2D * h_Yields = new TH2D("h_Yields","h_Yields",1000,0.,1000.,1000,0.,1000.);
	TH1D * h_xsec = new TH1D("h_xsec","h_xsec",1000,0.,1000.);

	ifstream fin("xsec_C1C1_wino.txt");

	double chiMass;
	double xsection;
	double error;
	vector<std::pair <double, double > > XSections;

	while (fin >> chiMass >> xsection >> error) {
		h_xsec->Fill(chiMass,xsection);
		XSections.push_back(std::make_pair(chiMass,xsection));
	}


    for (Long64_t jentry = 0; jentry < nentries; jentry++) {

        Long64_t ientry = tr.LoadTree(jentry);
        if (ientry < 0) break;
        tr.fChain->GetEntry(jentry);

	if (!(!tr.veto_dilepton && !tr.veto_thirdlepton && !tr.veto_otherlepton)) continue;
	if (!(tr.l1_againstMuon3>0.5 && tr.l1_againstElectronMVA6>0.5 && tr.l1_pt>40.)) continue;
	if (!(tr.l2_againstMuon3>0.5 && tr.l2_againstElectronMVA6>0.5 && tr.l2_pt>40.)) continue;
	if (!(tr.l1_byIsolationMVArun2v1DBoldDMwLT>4.5)) continue;
	if (!(tr.l2_byIsolationMVArun2v1DBoldDMwLT>4.5)) continue;
	if (!(tr.l1_charge != tr.l2_charge)) continue;
	if (!(tr.mvis>100 && tr.n_bjets==0 && tr.mt + tr.mt_leg2 > 150.)) continue;
/*
	if (!(tr.met_pt > 30)) continue;
	if (!(tr.mvis > 85 || tr.mvis < 55)) continue;
	if (!(tr.mt2 > 90 && tr.n_bjets == 0)) continue;
*/	
	if (!(tr.GenSusyMChargino == 400 && tr.GenSusyMNeutralino == 1)) continue;
	//if (!(tr.GenSusyMChargino == 200 && tr.GenSusyMNeutralino == 1)) continue;

	nevents += tr.weight;

	if (tr.GenSusyMChargino != 0 || tr.GenSusyMNeutralino != 0)
		h_Counting->Fill(tr.GenSusyMChargino,tr.GenSusyMNeutralino,tr.weight);
}

	for (int i = 1; i <= h_Counting->GetXaxis()->GetNbins(); i++) {
		double xsec= 0.;		
		for (int s = 0; s < XSections.size(); s++) {
			if (XSections[s].first == i) xsec = XSections[s].second;
		}
		for (int j = 1; j <= h_Counting->GetYaxis()->GetNbins(); j++) {

			int global = myWeights->FindBin(i,j);
			double sumWeight = myWeights->GetBinContent(global);

			if (sumWeight== 0) continue; 	
			//if (i == 400 && j == 1) {
			//cout<<""<<endl;
			//cout<<">>> myWeights->FindBin(i,j): "<<global<<endl;
			//cout<<">>> myWeights->GetBinContent(global): "<<sumWeight<<endl;
			//cout<<">>> xsec: "<<xsec<<endl;

		        global = h_Counting->FindBin(i,j);
			double content = h_Counting->GetBinContent(global);
			//cout<<">>> h_Counting->FindBin(i,j): "<<global<<endl;
			//cout<<">>> h_Counting->GetBinContent(global): "<<content<<endl;

			double yield = content * xsec * 35900. / (sumWeight * 1000.);

			h_Yields->SetBinContent(global,yield);
			//}	

		}
	}


	cout<<"Number of unweighted events (400,1): "<<nevents<<endl;

	int global = h_Counting->FindBin(400,1);
	double content = h_Counting->GetBinContent(global);

	cout<<"**** unweighted (400,1): "<<content<<endl;	

	global = h_Yields->FindBin(400,1);
	content = h_Yields->GetBinContent(global);

	cout<<">>> h_Yields->FindBin(i,j): "<<global<<endl;
	cout<<">>> h_Yields->GetBinContent(global): "<<content<<endl;

	cout<<"**** weighted (400,1): "<<content<<endl;	

	TFile * myfile = new TFile("file.root","recreate");

	h_xsec->Write();
	h_Counting->Write();
	h_Yields->Write();

	myfile->Close();

}
