//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun May 21 12:15:05 2017 by ROOT version 6.06/01
// from TTree tree/H2TauTauTreeProducerTauTau
// found on file: 200317/DiTauData/Tau_Run2016B_03Feb2017_v2/H2TauTauTreeProducerTauTau/tree.root
//////////////////////////////////////////////////////////

#ifndef treeBase_data_h
#define treeBase_data_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class treeBase_data {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   Int_t           bx;
   Int_t           orbit_number;
   Int_t           is_data;
   Double_t        nPU;
   Int_t           pass_leptons;
   Int_t           veto_dilepton;
   Int_t           veto_thirdlepton;
   Int_t           veto_otherlepton;
   Int_t           n_jets;
   Int_t           n_jets_puid;
   Int_t           n_jets_20;
   Int_t           n_jets_20_puid;
   Int_t           n_bjets;
   Int_t           n_jets_csvl;
   Int_t           n_vertices;
   Double_t        rho;
   Double_t        weight;
   Double_t        weight_vertex;
   Double_t        weight_njet;
   Double_t        weight_dy;
   Int_t           Flag_HBHENoiseFilter;
   Int_t           Flag_HBHENoiseIsoFilter;
   Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter;
   Int_t           Flag_goodVertices;
   Int_t           Flag_eeBadScFilter;
   Int_t           Flag_globalTightHalo2016Filter;
   Int_t           passBadMuonFilter;
   Int_t           passBadChargedHadronFilter;
   Int_t           badMuonMoriond2017;
   Int_t           badCloneMuonMoriond2017;
   Double_t        mvis;
   Double_t        dil_pt;
   Double_t        dil_eta;
   Double_t        dil_phi;
   Double_t        mt_total;
   Double_t        pzeta_met;
   Double_t        pzeta_vis;
   Double_t        pzeta_disc;
   Double_t        mt;
   Double_t        mt_leg2;
   Double_t        met_cov00;
   Double_t        met_cov10;
   Double_t        met_cov11;
   Double_t        met_phi;
   Double_t        met_pt;
   Double_t        pthiggs;
   Double_t        delta_r_l1_l2;
   Double_t        svfit_mass;
   Double_t        svfit_transverse_mass;
   Double_t        svfit_mass_error;
   Double_t        svfit_pt;
   Double_t        geninfo_mcweight;
   Int_t           geninfo_nup;
   Double_t        geninfo_htgen;
   Double_t        geninfo_invmass;
   Double_t        weight_gen;
   Double_t        genmet_pt;
   Double_t        genmet_phi;
   Double_t        vbf_mjj;
   Double_t        vbf_deta;
   Int_t           vbf_n_central20;
   Int_t           vbf_n_central;
   Double_t        vbf_jdphi;
   Double_t        vbf_dijetpt;
   Double_t        vbf_dijetphi;
   Double_t        vbf_dphidijethiggs;
   Double_t        vbf_mindetajetvis;
   Double_t        jet1_pt;
   Double_t        jet1_eta;
   Double_t        jet1_phi;
   Double_t        jet1_charge;
   Double_t        jet1_mass;
   Double_t        jet1_mva_pu;
   Double_t        jet1_id_pu;
   Double_t        jet1_flavour_parton;
   Double_t        jet1_csv;
   Double_t        jet1_genjet_pt;
   Double_t        jet2_pt;
   Double_t        jet2_eta;
   Double_t        jet2_phi;
   Double_t        jet2_charge;
   Double_t        jet2_mass;
   Double_t        jet2_mva_pu;
   Double_t        jet2_id_pu;
   Double_t        jet2_flavour_parton;
   Double_t        jet2_csv;
   Double_t        jet2_genjet_pt;
   Double_t        bjet1_pt;
   Double_t        bjet1_eta;
   Double_t        bjet1_phi;
   Double_t        bjet1_charge;
   Double_t        bjet1_mass;
   Double_t        bjet1_mva_pu;
   Double_t        bjet1_id_pu;
   Double_t        bjet1_flavour_parton;
   Double_t        bjet1_csv;
   Double_t        bjet1_genjet_pt;
   Double_t        bjet2_pt;
   Double_t        bjet2_eta;
   Double_t        bjet2_phi;
   Double_t        bjet2_charge;
   Double_t        bjet2_mass;
   Double_t        bjet2_mva_pu;
   Double_t        bjet2_id_pu;
   Double_t        bjet2_flavour_parton;
   Double_t        bjet2_csv;
   Double_t        bjet2_genjet_pt;
   Double_t        HT_allJets;
   Double_t        HT_jets;
   Double_t        HT_bJets;
   Double_t        HT_cleanJets;
   Double_t        HT_jets30;
   Double_t        HT_cleanJets30;
   Double_t        genboson_pt;
   Double_t        genboson_eta;
   Double_t        genboson_phi;
   Double_t        genboson_charge;
   Double_t        genboson_mass;
   Double_t        genboson_pdgId;
   Double_t        gen_top_1_pt;
   Double_t        gen_top_2_pt;
   Double_t        gen_top_weight;
   Double_t        puppimet_pt;
   Double_t        puppimet_phi;
   Double_t        puppimet_mt1;
   Double_t        puppimet_mt2;
   Double_t        puppimet_mttotal;
   Double_t        pfmet_pt;
   Double_t        pfmet_phi;
   Double_t        pfmet_mt1;
   Double_t        pfmet_mt2;
   Double_t        pfmet_mttotal;
   Double_t        l1_pt;
   Double_t        l1_eta;
   Double_t        l1_phi;
   Double_t        l1_charge;
   Double_t        l1_mass;
   Double_t        l1_jet_pt;
   Double_t        l1_jet_eta;
   Double_t        l1_jet_phi;
   Double_t        l1_jet_charge;
   Double_t        l1_jet_mass;
   Double_t        l1_dxy;
   Double_t        l1_dxy_error;
   Double_t        l1_dz;
   Double_t        l1_dz_error;
   Double_t        l1_weight;
   Double_t        l1_weight_trigger;
   Double_t        l1_weight_idiso;
   Int_t           l1_gen_match;
   Int_t           l1_decayMode;
   Double_t        l1_zImpact;
   Double_t        l1_dz_selfvertex;
   Double_t        l1_ptScale;
   Int_t           l1_againstElectronMVA6;
   Int_t           l1_againstMuon3;
   Double_t        l1_byCombinedIsolationDeltaBetaCorrRaw3Hits;
   Double_t        l1_byIsolationMVArun2v1DBoldDMwLTraw;
   Double_t        l1_byIsolationMVArun2v1DBnewDMwLTraw;
   Double_t        l1_byIsolationMVArun2v1DBdR03oldDMwLTraw;
   Int_t           l1_byCombinedIsolationDeltaBetaCorr3Hits;
   Int_t           l1_byIsolationMVArun2v1DBoldDMwLT;
   Int_t           l1_byIsolationMVArun2v1DBoldDMwLTNew;
   Int_t           l1_byIsolationMVArun2v1DBnewDMwLT;
   Int_t           l1_byIsolationMVArun2v1DBdR03oldDMwLT;
   Double_t        l1_chargedIsoPtSum;
   Double_t        l1_decayModeFinding;
   Double_t        l1_footprintCorrection;
   Double_t        l1_neutralIsoPtSum;
   Double_t        l1_puCorrPtSum;
   Double_t        l1_photonPtSumOutsideSignalCone;
   Double_t        l1_byTightIsolationMVArun2v1DBoldDMwLT;
   Double_t        l2_pt;
   Double_t        l2_eta;
   Double_t        l2_phi;
   Double_t        l2_charge;
   Double_t        l2_mass;
   Double_t        l2_jet_pt;
   Double_t        l2_jet_eta;
   Double_t        l2_jet_phi;
   Double_t        l2_jet_charge;
   Double_t        l2_jet_mass;
   Double_t        l2_dxy;
   Double_t        l2_dxy_error;
   Double_t        l2_dz;
   Double_t        l2_dz_error;
   Double_t        l2_weight;
   Double_t        l2_weight_trigger;
   Double_t        l2_weight_idiso;
   Int_t           l2_gen_match;
   Int_t           l2_decayMode;
   Double_t        l2_zImpact;
   Double_t        l2_dz_selfvertex;
   Double_t        l2_ptScale;
   Int_t           l2_againstElectronMVA6;
   Int_t           l2_againstMuon3;
   Double_t        l2_byCombinedIsolationDeltaBetaCorrRaw3Hits;
   Double_t        l2_byIsolationMVArun2v1DBoldDMwLTraw;
   Double_t        l2_byIsolationMVArun2v1DBnewDMwLTraw;
   Double_t        l2_byIsolationMVArun2v1DBdR03oldDMwLTraw;
   Int_t           l2_byCombinedIsolationDeltaBetaCorr3Hits;
   Int_t           l2_byIsolationMVArun2v1DBoldDMwLT;
   Int_t           l2_byIsolationMVArun2v1DBoldDMwLTNew;
   Int_t           l2_byIsolationMVArun2v1DBnewDMwLT;
   Int_t           l2_byIsolationMVArun2v1DBdR03oldDMwLT;
   Double_t        l2_chargedIsoPtSum;
   Double_t        l2_decayModeFinding;
   Double_t        l2_footprintCorrection;
   Double_t        l2_neutralIsoPtSum;
   Double_t        l2_puCorrPtSum;
   Double_t        l2_photonPtSumOutsideSignalCone;
   Double_t        l2_byTightIsolationMVArun2v1DBoldDMwLT;
   Double_t        l1_gen_pt;
   Double_t        l1_gen_eta;
   Double_t        l1_gen_phi;
   Double_t        l1_gen_charge;
   Double_t        l1_gen_mass;
   Double_t        l1_gen_pdgId;
   Double_t        l2_gen_pt;
   Double_t        l2_gen_eta;
   Double_t        l2_gen_phi;
   Double_t        l2_gen_charge;
   Double_t        l2_gen_mass;
   Double_t        l2_gen_pdgId;
   Double_t        l1_gen_vis_pt;
   Double_t        l1_gen_vis_eta;
   Double_t        l1_gen_vis_phi;
   Double_t        l1_gen_vis_charge;
   Double_t        l1_gen_vis_mass;
   Double_t        l2_gen_vis_pt;
   Double_t        l2_gen_vis_eta;
   Double_t        l2_gen_vis_phi;
   Double_t        l2_gen_vis_charge;
   Double_t        l2_gen_vis_mass;
   Int_t           l1_gen_decaymode;
   Int_t           l2_gen_decaymode;
   Double_t        l1_trigger_weight;
   Double_t        l1_trigger_weight_up;
   Double_t        l1_trigger_weight_down;
   Double_t        l2_trigger_weight;
   Double_t        l2_trigger_weight_up;
   Double_t        l2_trigger_weight_down;
   Double_t        mt2;
   Double_t        mt2_lep;
   Double_t        mt2_mvamet;
   Double_t        mt2_rawpfmet;
   Double_t        minDphiMETJets;
   Double_t        trigger_ditau35;
   Double_t        trigger_ditau35_combiso;
   Double_t        trigger_singletau140;
   Double_t        trigger_singletau120;
   Double_t        trigger_matched_ditau35;
   Double_t        trigger_matched_ditau35_combiso;
   Double_t        trigger_matched_singletau140;
   Double_t        trigger_matched_singletau120;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_orbit_number;   //!
   TBranch        *b_is_data;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_pass_leptons;   //!
   TBranch        *b_veto_dilepton;   //!
   TBranch        *b_veto_thirdlepton;   //!
   TBranch        *b_veto_otherlepton;   //!
   TBranch        *b_n_jets;   //!
   TBranch        *b_n_jets_puid;   //!
   TBranch        *b_n_jets_20;   //!
   TBranch        *b_n_jets_20_puid;   //!
   TBranch        *b_n_bjets;   //!
   TBranch        *b_n_jets_csvl;   //!
   TBranch        *b_n_vertices;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_weight_vertex;   //!
   TBranch        *b_weight_njet;   //!
   TBranch        *b_weight_dy;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_passBadMuonFilter;   //!
   TBranch        *b_passBadChargedHadronFilter;   //!
   TBranch        *b_badMuonMoriond2017;   //!
   TBranch        *b_badCloneMuonMoriond2017;   //!
   TBranch        *b_mvis;   //!
   TBranch        *b_dil_pt;   //!
   TBranch        *b_dil_eta;   //!
   TBranch        *b_dil_phi;   //!
   TBranch        *b_mt_total;   //!
   TBranch        *b_pzeta_met;   //!
   TBranch        *b_pzeta_vis;   //!
   TBranch        *b_pzeta_disc;   //!
   TBranch        *b_mt;   //!
   TBranch        *b_mt_leg2;   //!
   TBranch        *b_met_cov00;   //!
   TBranch        *b_met_cov10;   //!
   TBranch        *b_met_cov11;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_pthiggs;   //!
   TBranch        *b_delta_r_l1_l2;   //!
   TBranch        *b_svfit_mass;   //!
   TBranch        *b_svfit_transverse_mass;   //!
   TBranch        *b_svfit_mass_error;   //!
   TBranch        *b_svfit_pt;   //!
   TBranch        *b_geninfo_mcweight;   //!
   TBranch        *b_geninfo_nup;   //!
   TBranch        *b_geninfo_htgen;   //!
   TBranch        *b_geninfo_invmass;   //!
   TBranch        *b_weight_gen;   //!
   TBranch        *b_genmet_pt;   //!
   TBranch        *b_genmet_phi;   //!
   TBranch        *b_vbf_mjj;   //!
   TBranch        *b_vbf_deta;   //!
   TBranch        *b_vbf_n_central20;   //!
   TBranch        *b_vbf_n_central;   //!
   TBranch        *b_vbf_jdphi;   //!
   TBranch        *b_vbf_dijetpt;   //!
   TBranch        *b_vbf_dijetphi;   //!
   TBranch        *b_vbf_dphidijethiggs;   //!
   TBranch        *b_vbf_mindetajetvis;   //!
   TBranch        *b_jet1_pt;   //!
   TBranch        *b_jet1_eta;   //!
   TBranch        *b_jet1_phi;   //!
   TBranch        *b_jet1_charge;   //!
   TBranch        *b_jet1_mass;   //!
   TBranch        *b_jet1_mva_pu;   //!
   TBranch        *b_jet1_id_pu;   //!
   TBranch        *b_jet1_flavour_parton;   //!
   TBranch        *b_jet1_csv;   //!
   TBranch        *b_jet1_genjet_pt;   //!
   TBranch        *b_jet2_pt;   //!
   TBranch        *b_jet2_eta;   //!
   TBranch        *b_jet2_phi;   //!
   TBranch        *b_jet2_charge;   //!
   TBranch        *b_jet2_mass;   //!
   TBranch        *b_jet2_mva_pu;   //!
   TBranch        *b_jet2_id_pu;   //!
   TBranch        *b_jet2_flavour_parton;   //!
   TBranch        *b_jet2_csv;   //!
   TBranch        *b_jet2_genjet_pt;   //!
   TBranch        *b_bjet1_pt;   //!
   TBranch        *b_bjet1_eta;   //!
   TBranch        *b_bjet1_phi;   //!
   TBranch        *b_bjet1_charge;   //!
   TBranch        *b_bjet1_mass;   //!
   TBranch        *b_bjet1_mva_pu;   //!
   TBranch        *b_bjet1_id_pu;   //!
   TBranch        *b_bjet1_flavour_parton;   //!
   TBranch        *b_bjet1_csv;   //!
   TBranch        *b_bjet1_genjet_pt;   //!
   TBranch        *b_bjet2_pt;   //!
   TBranch        *b_bjet2_eta;   //!
   TBranch        *b_bjet2_phi;   //!
   TBranch        *b_bjet2_charge;   //!
   TBranch        *b_bjet2_mass;   //!
   TBranch        *b_bjet2_mva_pu;   //!
   TBranch        *b_bjet2_id_pu;   //!
   TBranch        *b_bjet2_flavour_parton;   //!
   TBranch        *b_bjet2_csv;   //!
   TBranch        *b_bjet2_genjet_pt;   //!
   TBranch        *b_HT_allJets;   //!
   TBranch        *b_HT_jets;   //!
   TBranch        *b_HT_bJets;   //!
   TBranch        *b_HT_cleanJets;   //!
   TBranch        *b_HT_jets30;   //!
   TBranch        *b_HT_cleanJets30;   //!
   TBranch        *b_genboson_pt;   //!
   TBranch        *b_genboson_eta;   //!
   TBranch        *b_genboson_phi;   //!
   TBranch        *b_genboson_charge;   //!
   TBranch        *b_genboson_mass;   //!
   TBranch        *b_genboson_pdgId;   //!
   TBranch        *b_gen_top_1_pt;   //!
   TBranch        *b_gen_top_2_pt;   //!
   TBranch        *b_gen_top_weight;   //!
   TBranch        *b_puppimet_pt;   //!
   TBranch        *b_puppimet_phi;   //!
   TBranch        *b_puppimet_mt1;   //!
   TBranch        *b_puppimet_mt2;   //!
   TBranch        *b_puppimet_mttotal;   //!
   TBranch        *b_pfmet_pt;   //!
   TBranch        *b_pfmet_phi;   //!
   TBranch        *b_pfmet_mt1;   //!
   TBranch        *b_pfmet_mt2;   //!
   TBranch        *b_pfmet_mttotal;   //!
   TBranch        *b_l1_pt;   //!
   TBranch        *b_l1_eta;   //!
   TBranch        *b_l1_phi;   //!
   TBranch        *b_l1_charge;   //!
   TBranch        *b_l1_mass;   //!
   TBranch        *b_l1_jet_pt;   //!
   TBranch        *b_l1_jet_eta;   //!
   TBranch        *b_l1_jet_phi;   //!
   TBranch        *b_l1_jet_charge;   //!
   TBranch        *b_l1_jet_mass;   //!
   TBranch        *b_l1_dxy;   //!
   TBranch        *b_l1_dxy_error;   //!
   TBranch        *b_l1_dz;   //!
   TBranch        *b_l1_dz_error;   //!
   TBranch        *b_l1_weight;   //!
   TBranch        *b_l1_weight_trigger;   //!
   TBranch        *b_l1_weight_idiso;   //!
   TBranch        *b_l1_gen_match;   //!
   TBranch        *b_l1_decayMode;   //!
   TBranch        *b_l1_zImpact;   //!
   TBranch        *b_l1_dz_selfvertex;   //!
   TBranch        *b_l1_ptScale;   //!
   TBranch        *b_l1_againstElectronMVA6;   //!
   TBranch        *b_l1_againstMuon3;   //!
   TBranch        *b_l1_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_l1_byIsolationMVArun2v1DBoldDMwLTraw;   //!
   TBranch        *b_l1_byIsolationMVArun2v1DBnewDMwLTraw;   //!
   TBranch        *b_l1_byIsolationMVArun2v1DBdR03oldDMwLTraw;   //!
   TBranch        *b_l1_byCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_l1_byIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_l1_byIsolationMVArun2v1DBoldDMwLTNew;   //!
   TBranch        *b_l1_byIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_l1_byIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_l1_chargedIsoPtSum;   //!
   TBranch        *b_l1_decayModeFinding;   //!
   TBranch        *b_l1_footprintCorrection;   //!
   TBranch        *b_l1_neutralIsoPtSum;   //!
   TBranch        *b_l1_puCorrPtSum;   //!
   TBranch        *b_l1_photonPtSumOutsideSignalCone;   //!
   TBranch        *b_l1_byTightIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_l2_pt;   //!
   TBranch        *b_l2_eta;   //!
   TBranch        *b_l2_phi;   //!
   TBranch        *b_l2_charge;   //!
   TBranch        *b_l2_mass;   //!
   TBranch        *b_l2_jet_pt;   //!
   TBranch        *b_l2_jet_eta;   //!
   TBranch        *b_l2_jet_phi;   //!
   TBranch        *b_l2_jet_charge;   //!
   TBranch        *b_l2_jet_mass;   //!
   TBranch        *b_l2_dxy;   //!
   TBranch        *b_l2_dxy_error;   //!
   TBranch        *b_l2_dz;   //!
   TBranch        *b_l2_dz_error;   //!
   TBranch        *b_l2_weight;   //!
   TBranch        *b_l2_weight_trigger;   //!
   TBranch        *b_l2_weight_idiso;   //!
   TBranch        *b_l2_gen_match;   //!
   TBranch        *b_l2_decayMode;   //!
   TBranch        *b_l2_zImpact;   //!
   TBranch        *b_l2_dz_selfvertex;   //!
   TBranch        *b_l2_ptScale;   //!
   TBranch        *b_l2_againstElectronMVA6;   //!
   TBranch        *b_l2_againstMuon3;   //!
   TBranch        *b_l2_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBoldDMwLTraw;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBnewDMwLTraw;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBdR03oldDMwLTraw;   //!
   TBranch        *b_l2_byCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBoldDMwLTNew;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_l2_chargedIsoPtSum;   //!
   TBranch        *b_l2_decayModeFinding;   //!
   TBranch        *b_l2_footprintCorrection;   //!
   TBranch        *b_l2_neutralIsoPtSum;   //!
   TBranch        *b_l2_puCorrPtSum;   //!
   TBranch        *b_l2_photonPtSumOutsideSignalCone;   //!
   TBranch        *b_l2_byTightIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_l1_gen_pt;   //!
   TBranch        *b_l1_gen_eta;   //!
   TBranch        *b_l1_gen_phi;   //!
   TBranch        *b_l1_gen_charge;   //!
   TBranch        *b_l1_gen_mass;   //!
   TBranch        *b_l1_gen_pdgId;   //!
   TBranch        *b_l2_gen_pt;   //!
   TBranch        *b_l2_gen_eta;   //!
   TBranch        *b_l2_gen_phi;   //!
   TBranch        *b_l2_gen_charge;   //!
   TBranch        *b_l2_gen_mass;   //!
   TBranch        *b_l2_gen_pdgId;   //!
   TBranch        *b_l1_gen_vis_pt;   //!
   TBranch        *b_l1_gen_vis_eta;   //!
   TBranch        *b_l1_gen_vis_phi;   //!
   TBranch        *b_l1_gen_vis_charge;   //!
   TBranch        *b_l1_gen_vis_mass;   //!
   TBranch        *b_l2_gen_vis_pt;   //!
   TBranch        *b_l2_gen_vis_eta;   //!
   TBranch        *b_l2_gen_vis_phi;   //!
   TBranch        *b_l2_gen_vis_charge;   //!
   TBranch        *b_l2_gen_vis_mass;   //!
   TBranch        *b_l1_gen_decaymode;   //!
   TBranch        *b_l2_gen_decaymode;   //!
   TBranch        *b_l1_trigger_weight;   //!
   TBranch        *b_l1_trigger_weight_up;   //!
   TBranch        *b_l1_trigger_weight_down;   //!
   TBranch        *b_l2_trigger_weight;   //!
   TBranch        *b_l2_trigger_weight_up;   //!
   TBranch        *b_l2_trigger_weight_down;   //!
   TBranch        *b_mt2;   //!
   TBranch        *b_mt2_lep;   //!
   TBranch        *b_mt2_mvamet;   //!
   TBranch        *b_mt2_rawpfmet;   //!
   TBranch        *b_minDphiMETJets;   //!
   TBranch        *b_trigger_ditau35;   //!
   TBranch        *b_trigger_ditau35_combiso;   //!
   TBranch        *b_trigger_singletau140;   //!
   TBranch        *b_trigger_singletau120;   //!
   TBranch        *b_trigger_matched_ditau35;   //!
   TBranch        *b_trigger_matched_ditau35_combiso;   //!
   TBranch        *b_trigger_matched_singletau140;   //!
   TBranch        *b_trigger_matched_singletau120;   //!

   treeBase_data(TTree *tree=0);
   virtual ~treeBase_data();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

treeBase_data::treeBase_data(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("200317/DiTauData/Tau_Run2016B_03Feb2017_v2/H2TauTauTreeProducerTauTau/tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("200317/DiTauData/Tau_Run2016B_03Feb2017_v2/H2TauTauTreeProducerTauTau/tree.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

treeBase_data::~treeBase_data()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t treeBase_data::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treeBase_data::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void treeBase_data::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("orbit_number", &orbit_number, &b_orbit_number);
   fChain->SetBranchAddress("is_data", &is_data, &b_is_data);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("pass_leptons", &pass_leptons, &b_pass_leptons);
   fChain->SetBranchAddress("veto_dilepton", &veto_dilepton, &b_veto_dilepton);
   fChain->SetBranchAddress("veto_thirdlepton", &veto_thirdlepton, &b_veto_thirdlepton);
   fChain->SetBranchAddress("veto_otherlepton", &veto_otherlepton, &b_veto_otherlepton);
   fChain->SetBranchAddress("n_jets", &n_jets, &b_n_jets);
   fChain->SetBranchAddress("n_jets_puid", &n_jets_puid, &b_n_jets_puid);
   fChain->SetBranchAddress("n_jets_20", &n_jets_20, &b_n_jets_20);
   fChain->SetBranchAddress("n_jets_20_puid", &n_jets_20_puid, &b_n_jets_20_puid);
   fChain->SetBranchAddress("n_bjets", &n_bjets, &b_n_bjets);
   fChain->SetBranchAddress("n_jets_csvl", &n_jets_csvl, &b_n_jets_csvl);
   fChain->SetBranchAddress("n_vertices", &n_vertices, &b_n_vertices);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("weight_vertex", &weight_vertex, &b_weight_vertex);
   fChain->SetBranchAddress("weight_njet", &weight_njet, &b_weight_njet);
   fChain->SetBranchAddress("weight_dy", &weight_dy, &b_weight_dy);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("passBadMuonFilter", &passBadMuonFilter, &b_passBadMuonFilter);
   fChain->SetBranchAddress("passBadChargedHadronFilter", &passBadChargedHadronFilter, &b_passBadChargedHadronFilter);
   fChain->SetBranchAddress("badMuonMoriond2017", &badMuonMoriond2017, &b_badMuonMoriond2017);
   fChain->SetBranchAddress("badCloneMuonMoriond2017", &badCloneMuonMoriond2017, &b_badCloneMuonMoriond2017);
   fChain->SetBranchAddress("mvis", &mvis, &b_mvis);
   fChain->SetBranchAddress("dil_pt", &dil_pt, &b_dil_pt);
   fChain->SetBranchAddress("dil_eta", &dil_eta, &b_dil_eta);
   fChain->SetBranchAddress("dil_phi", &dil_phi, &b_dil_phi);
   fChain->SetBranchAddress("mt_total", &mt_total, &b_mt_total);
   fChain->SetBranchAddress("pzeta_met", &pzeta_met, &b_pzeta_met);
   fChain->SetBranchAddress("pzeta_vis", &pzeta_vis, &b_pzeta_vis);
   fChain->SetBranchAddress("pzeta_disc", &pzeta_disc, &b_pzeta_disc);
   fChain->SetBranchAddress("mt", &mt, &b_mt);
   fChain->SetBranchAddress("mt_leg2", &mt_leg2, &b_mt_leg2);
   fChain->SetBranchAddress("met_cov00", &met_cov00, &b_met_cov00);
   fChain->SetBranchAddress("met_cov10", &met_cov10, &b_met_cov10);
   fChain->SetBranchAddress("met_cov11", &met_cov11, &b_met_cov11);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("pthiggs", &pthiggs, &b_pthiggs);
   fChain->SetBranchAddress("delta_r_l1_l2", &delta_r_l1_l2, &b_delta_r_l1_l2);
   fChain->SetBranchAddress("svfit_mass", &svfit_mass, &b_svfit_mass);
   fChain->SetBranchAddress("svfit_transverse_mass", &svfit_transverse_mass, &b_svfit_transverse_mass);
   fChain->SetBranchAddress("svfit_mass_error", &svfit_mass_error, &b_svfit_mass_error);
   fChain->SetBranchAddress("svfit_pt", &svfit_pt, &b_svfit_pt);
   fChain->SetBranchAddress("geninfo_mcweight", &geninfo_mcweight, &b_geninfo_mcweight);
   fChain->SetBranchAddress("geninfo_nup", &geninfo_nup, &b_geninfo_nup);
   fChain->SetBranchAddress("geninfo_htgen", &geninfo_htgen, &b_geninfo_htgen);
   fChain->SetBranchAddress("geninfo_invmass", &geninfo_invmass, &b_geninfo_invmass);
   fChain->SetBranchAddress("weight_gen", &weight_gen, &b_weight_gen);
   fChain->SetBranchAddress("genmet_pt", &genmet_pt, &b_genmet_pt);
   fChain->SetBranchAddress("genmet_phi", &genmet_phi, &b_genmet_phi);
   fChain->SetBranchAddress("vbf_mjj", &vbf_mjj, &b_vbf_mjj);
   fChain->SetBranchAddress("vbf_deta", &vbf_deta, &b_vbf_deta);
   fChain->SetBranchAddress("vbf_n_central20", &vbf_n_central20, &b_vbf_n_central20);
   fChain->SetBranchAddress("vbf_n_central", &vbf_n_central, &b_vbf_n_central);
   fChain->SetBranchAddress("vbf_jdphi", &vbf_jdphi, &b_vbf_jdphi);
   fChain->SetBranchAddress("vbf_dijetpt", &vbf_dijetpt, &b_vbf_dijetpt);
   fChain->SetBranchAddress("vbf_dijetphi", &vbf_dijetphi, &b_vbf_dijetphi);
   fChain->SetBranchAddress("vbf_dphidijethiggs", &vbf_dphidijethiggs, &b_vbf_dphidijethiggs);
   fChain->SetBranchAddress("vbf_mindetajetvis", &vbf_mindetajetvis, &b_vbf_mindetajetvis);
   fChain->SetBranchAddress("jet1_pt", &jet1_pt, &b_jet1_pt);
   fChain->SetBranchAddress("jet1_eta", &jet1_eta, &b_jet1_eta);
   fChain->SetBranchAddress("jet1_phi", &jet1_phi, &b_jet1_phi);
   fChain->SetBranchAddress("jet1_charge", &jet1_charge, &b_jet1_charge);
   fChain->SetBranchAddress("jet1_mass", &jet1_mass, &b_jet1_mass);
   fChain->SetBranchAddress("jet1_mva_pu", &jet1_mva_pu, &b_jet1_mva_pu);
   fChain->SetBranchAddress("jet1_id_pu", &jet1_id_pu, &b_jet1_id_pu);
   fChain->SetBranchAddress("jet1_flavour_parton", &jet1_flavour_parton, &b_jet1_flavour_parton);
   fChain->SetBranchAddress("jet1_csv", &jet1_csv, &b_jet1_csv);
   fChain->SetBranchAddress("jet1_genjet_pt", &jet1_genjet_pt, &b_jet1_genjet_pt);
   fChain->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   fChain->SetBranchAddress("jet2_eta", &jet2_eta, &b_jet2_eta);
   fChain->SetBranchAddress("jet2_phi", &jet2_phi, &b_jet2_phi);
   fChain->SetBranchAddress("jet2_charge", &jet2_charge, &b_jet2_charge);
   fChain->SetBranchAddress("jet2_mass", &jet2_mass, &b_jet2_mass);
   fChain->SetBranchAddress("jet2_mva_pu", &jet2_mva_pu, &b_jet2_mva_pu);
   fChain->SetBranchAddress("jet2_id_pu", &jet2_id_pu, &b_jet2_id_pu);
   fChain->SetBranchAddress("jet2_flavour_parton", &jet2_flavour_parton, &b_jet2_flavour_parton);
   fChain->SetBranchAddress("jet2_csv", &jet2_csv, &b_jet2_csv);
   fChain->SetBranchAddress("jet2_genjet_pt", &jet2_genjet_pt, &b_jet2_genjet_pt);
   fChain->SetBranchAddress("bjet1_pt", &bjet1_pt, &b_bjet1_pt);
   fChain->SetBranchAddress("bjet1_eta", &bjet1_eta, &b_bjet1_eta);
   fChain->SetBranchAddress("bjet1_phi", &bjet1_phi, &b_bjet1_phi);
   fChain->SetBranchAddress("bjet1_charge", &bjet1_charge, &b_bjet1_charge);
   fChain->SetBranchAddress("bjet1_mass", &bjet1_mass, &b_bjet1_mass);
   fChain->SetBranchAddress("bjet1_mva_pu", &bjet1_mva_pu, &b_bjet1_mva_pu);
   fChain->SetBranchAddress("bjet1_id_pu", &bjet1_id_pu, &b_bjet1_id_pu);
   fChain->SetBranchAddress("bjet1_flavour_parton", &bjet1_flavour_parton, &b_bjet1_flavour_parton);
   fChain->SetBranchAddress("bjet1_csv", &bjet1_csv, &b_bjet1_csv);
   fChain->SetBranchAddress("bjet1_genjet_pt", &bjet1_genjet_pt, &b_bjet1_genjet_pt);
   fChain->SetBranchAddress("bjet2_pt", &bjet2_pt, &b_bjet2_pt);
   fChain->SetBranchAddress("bjet2_eta", &bjet2_eta, &b_bjet2_eta);
   fChain->SetBranchAddress("bjet2_phi", &bjet2_phi, &b_bjet2_phi);
   fChain->SetBranchAddress("bjet2_charge", &bjet2_charge, &b_bjet2_charge);
   fChain->SetBranchAddress("bjet2_mass", &bjet2_mass, &b_bjet2_mass);
   fChain->SetBranchAddress("bjet2_mva_pu", &bjet2_mva_pu, &b_bjet2_mva_pu);
   fChain->SetBranchAddress("bjet2_id_pu", &bjet2_id_pu, &b_bjet2_id_pu);
   fChain->SetBranchAddress("bjet2_flavour_parton", &bjet2_flavour_parton, &b_bjet2_flavour_parton);
   fChain->SetBranchAddress("bjet2_csv", &bjet2_csv, &b_bjet2_csv);
   fChain->SetBranchAddress("bjet2_genjet_pt", &bjet2_genjet_pt, &b_bjet2_genjet_pt);
   fChain->SetBranchAddress("HT_allJets", &HT_allJets, &b_HT_allJets);
   fChain->SetBranchAddress("HT_jets", &HT_jets, &b_HT_jets);
   fChain->SetBranchAddress("HT_bJets", &HT_bJets, &b_HT_bJets);
   fChain->SetBranchAddress("HT_cleanJets", &HT_cleanJets, &b_HT_cleanJets);
   fChain->SetBranchAddress("HT_jets30", &HT_jets30, &b_HT_jets30);
   fChain->SetBranchAddress("HT_cleanJets30", &HT_cleanJets30, &b_HT_cleanJets30);
   fChain->SetBranchAddress("genboson_pt", &genboson_pt, &b_genboson_pt);
   fChain->SetBranchAddress("genboson_eta", &genboson_eta, &b_genboson_eta);
   fChain->SetBranchAddress("genboson_phi", &genboson_phi, &b_genboson_phi);
   fChain->SetBranchAddress("genboson_charge", &genboson_charge, &b_genboson_charge);
   fChain->SetBranchAddress("genboson_mass", &genboson_mass, &b_genboson_mass);
   fChain->SetBranchAddress("genboson_pdgId", &genboson_pdgId, &b_genboson_pdgId);
   fChain->SetBranchAddress("gen_top_1_pt", &gen_top_1_pt, &b_gen_top_1_pt);
   fChain->SetBranchAddress("gen_top_2_pt", &gen_top_2_pt, &b_gen_top_2_pt);
   fChain->SetBranchAddress("gen_top_weight", &gen_top_weight, &b_gen_top_weight);
   fChain->SetBranchAddress("puppimet_pt", &puppimet_pt, &b_puppimet_pt);
   fChain->SetBranchAddress("puppimet_phi", &puppimet_phi, &b_puppimet_phi);
   fChain->SetBranchAddress("puppimet_mt1", &puppimet_mt1, &b_puppimet_mt1);
   fChain->SetBranchAddress("puppimet_mt2", &puppimet_mt2, &b_puppimet_mt2);
   fChain->SetBranchAddress("puppimet_mttotal", &puppimet_mttotal, &b_puppimet_mttotal);
   fChain->SetBranchAddress("pfmet_pt", &pfmet_pt, &b_pfmet_pt);
   fChain->SetBranchAddress("pfmet_phi", &pfmet_phi, &b_pfmet_phi);
   fChain->SetBranchAddress("pfmet_mt1", &pfmet_mt1, &b_pfmet_mt1);
   fChain->SetBranchAddress("pfmet_mt2", &pfmet_mt2, &b_pfmet_mt2);
   fChain->SetBranchAddress("pfmet_mttotal", &pfmet_mttotal, &b_pfmet_mttotal);
   fChain->SetBranchAddress("l1_pt", &l1_pt, &b_l1_pt);
   fChain->SetBranchAddress("l1_eta", &l1_eta, &b_l1_eta);
   fChain->SetBranchAddress("l1_phi", &l1_phi, &b_l1_phi);
   fChain->SetBranchAddress("l1_charge", &l1_charge, &b_l1_charge);
   fChain->SetBranchAddress("l1_mass", &l1_mass, &b_l1_mass);
   fChain->SetBranchAddress("l1_jet_pt", &l1_jet_pt, &b_l1_jet_pt);
   fChain->SetBranchAddress("l1_jet_eta", &l1_jet_eta, &b_l1_jet_eta);
   fChain->SetBranchAddress("l1_jet_phi", &l1_jet_phi, &b_l1_jet_phi);
   fChain->SetBranchAddress("l1_jet_charge", &l1_jet_charge, &b_l1_jet_charge);
   fChain->SetBranchAddress("l1_jet_mass", &l1_jet_mass, &b_l1_jet_mass);
   fChain->SetBranchAddress("l1_dxy", &l1_dxy, &b_l1_dxy);
   fChain->SetBranchAddress("l1_dxy_error", &l1_dxy_error, &b_l1_dxy_error);
   fChain->SetBranchAddress("l1_dz", &l1_dz, &b_l1_dz);
   fChain->SetBranchAddress("l1_dz_error", &l1_dz_error, &b_l1_dz_error);
   fChain->SetBranchAddress("l1_weight", &l1_weight, &b_l1_weight);
   fChain->SetBranchAddress("l1_weight_trigger", &l1_weight_trigger, &b_l1_weight_trigger);
   fChain->SetBranchAddress("l1_weight_idiso", &l1_weight_idiso, &b_l1_weight_idiso);
   fChain->SetBranchAddress("l1_gen_match", &l1_gen_match, &b_l1_gen_match);
   fChain->SetBranchAddress("l1_decayMode", &l1_decayMode, &b_l1_decayMode);
   fChain->SetBranchAddress("l1_zImpact", &l1_zImpact, &b_l1_zImpact);
   fChain->SetBranchAddress("l1_dz_selfvertex", &l1_dz_selfvertex, &b_l1_dz_selfvertex);
   fChain->SetBranchAddress("l1_ptScale", &l1_ptScale, &b_l1_ptScale);
   fChain->SetBranchAddress("l1_againstElectronMVA6", &l1_againstElectronMVA6, &b_l1_againstElectronMVA6);
   fChain->SetBranchAddress("l1_againstMuon3", &l1_againstMuon3, &b_l1_againstMuon3);
   fChain->SetBranchAddress("l1_byCombinedIsolationDeltaBetaCorrRaw3Hits", &l1_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_l1_byCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("l1_byIsolationMVArun2v1DBoldDMwLTraw", &l1_byIsolationMVArun2v1DBoldDMwLTraw, &b_l1_byIsolationMVArun2v1DBoldDMwLTraw);
   fChain->SetBranchAddress("l1_byIsolationMVArun2v1DBnewDMwLTraw", &l1_byIsolationMVArun2v1DBnewDMwLTraw, &b_l1_byIsolationMVArun2v1DBnewDMwLTraw);
   fChain->SetBranchAddress("l1_byIsolationMVArun2v1DBdR03oldDMwLTraw", &l1_byIsolationMVArun2v1DBdR03oldDMwLTraw, &b_l1_byIsolationMVArun2v1DBdR03oldDMwLTraw);
   fChain->SetBranchAddress("l1_byCombinedIsolationDeltaBetaCorr3Hits", &l1_byCombinedIsolationDeltaBetaCorr3Hits, &b_l1_byCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("l1_byIsolationMVArun2v1DBoldDMwLT", &l1_byIsolationMVArun2v1DBoldDMwLT, &b_l1_byIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("l1_byIsolationMVArun2v1DBoldDMwLTNew", &l1_byIsolationMVArun2v1DBoldDMwLTNew, &b_l1_byIsolationMVArun2v1DBoldDMwLTNew);
   fChain->SetBranchAddress("l1_byIsolationMVArun2v1DBnewDMwLT", &l1_byIsolationMVArun2v1DBnewDMwLT, &b_l1_byIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("l1_byIsolationMVArun2v1DBdR03oldDMwLT", &l1_byIsolationMVArun2v1DBdR03oldDMwLT, &b_l1_byIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("l1_chargedIsoPtSum", &l1_chargedIsoPtSum, &b_l1_chargedIsoPtSum);
   fChain->SetBranchAddress("l1_decayModeFinding", &l1_decayModeFinding, &b_l1_decayModeFinding);
   fChain->SetBranchAddress("l1_footprintCorrection", &l1_footprintCorrection, &b_l1_footprintCorrection);
   fChain->SetBranchAddress("l1_neutralIsoPtSum", &l1_neutralIsoPtSum, &b_l1_neutralIsoPtSum);
   fChain->SetBranchAddress("l1_puCorrPtSum", &l1_puCorrPtSum, &b_l1_puCorrPtSum);
   fChain->SetBranchAddress("l1_photonPtSumOutsideSignalCone", &l1_photonPtSumOutsideSignalCone, &b_l1_photonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("l1_byTightIsolationMVArun2v1DBoldDMwLT", &l1_byTightIsolationMVArun2v1DBoldDMwLT, &b_l1_byTightIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("l2_pt", &l2_pt, &b_l2_pt);
   fChain->SetBranchAddress("l2_eta", &l2_eta, &b_l2_eta);
   fChain->SetBranchAddress("l2_phi", &l2_phi, &b_l2_phi);
   fChain->SetBranchAddress("l2_charge", &l2_charge, &b_l2_charge);
   fChain->SetBranchAddress("l2_mass", &l2_mass, &b_l2_mass);
   fChain->SetBranchAddress("l2_jet_pt", &l2_jet_pt, &b_l2_jet_pt);
   fChain->SetBranchAddress("l2_jet_eta", &l2_jet_eta, &b_l2_jet_eta);
   fChain->SetBranchAddress("l2_jet_phi", &l2_jet_phi, &b_l2_jet_phi);
   fChain->SetBranchAddress("l2_jet_charge", &l2_jet_charge, &b_l2_jet_charge);
   fChain->SetBranchAddress("l2_jet_mass", &l2_jet_mass, &b_l2_jet_mass);
   fChain->SetBranchAddress("l2_dxy", &l2_dxy, &b_l2_dxy);
   fChain->SetBranchAddress("l2_dxy_error", &l2_dxy_error, &b_l2_dxy_error);
   fChain->SetBranchAddress("l2_dz", &l2_dz, &b_l2_dz);
   fChain->SetBranchAddress("l2_dz_error", &l2_dz_error, &b_l2_dz_error);
   fChain->SetBranchAddress("l2_weight", &l2_weight, &b_l2_weight);
   fChain->SetBranchAddress("l2_weight_trigger", &l2_weight_trigger, &b_l2_weight_trigger);
   fChain->SetBranchAddress("l2_weight_idiso", &l2_weight_idiso, &b_l2_weight_idiso);
   fChain->SetBranchAddress("l2_gen_match", &l2_gen_match, &b_l2_gen_match);
   fChain->SetBranchAddress("l2_decayMode", &l2_decayMode, &b_l2_decayMode);
   fChain->SetBranchAddress("l2_zImpact", &l2_zImpact, &b_l2_zImpact);
   fChain->SetBranchAddress("l2_dz_selfvertex", &l2_dz_selfvertex, &b_l2_dz_selfvertex);
   fChain->SetBranchAddress("l2_ptScale", &l2_ptScale, &b_l2_ptScale);
   fChain->SetBranchAddress("l2_againstElectronMVA6", &l2_againstElectronMVA6, &b_l2_againstElectronMVA6);
   fChain->SetBranchAddress("l2_againstMuon3", &l2_againstMuon3, &b_l2_againstMuon3);
   fChain->SetBranchAddress("l2_byCombinedIsolationDeltaBetaCorrRaw3Hits", &l2_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_l2_byCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBoldDMwLTraw", &l2_byIsolationMVArun2v1DBoldDMwLTraw, &b_l2_byIsolationMVArun2v1DBoldDMwLTraw);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBnewDMwLTraw", &l2_byIsolationMVArun2v1DBnewDMwLTraw, &b_l2_byIsolationMVArun2v1DBnewDMwLTraw);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBdR03oldDMwLTraw", &l2_byIsolationMVArun2v1DBdR03oldDMwLTraw, &b_l2_byIsolationMVArun2v1DBdR03oldDMwLTraw);
   fChain->SetBranchAddress("l2_byCombinedIsolationDeltaBetaCorr3Hits", &l2_byCombinedIsolationDeltaBetaCorr3Hits, &b_l2_byCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBoldDMwLT", &l2_byIsolationMVArun2v1DBoldDMwLT, &b_l2_byIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBoldDMwLTNew", &l2_byIsolationMVArun2v1DBoldDMwLTNew, &b_l2_byIsolationMVArun2v1DBoldDMwLTNew);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBnewDMwLT", &l2_byIsolationMVArun2v1DBnewDMwLT, &b_l2_byIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBdR03oldDMwLT", &l2_byIsolationMVArun2v1DBdR03oldDMwLT, &b_l2_byIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("l2_chargedIsoPtSum", &l2_chargedIsoPtSum, &b_l2_chargedIsoPtSum);
   fChain->SetBranchAddress("l2_decayModeFinding", &l2_decayModeFinding, &b_l2_decayModeFinding);
   fChain->SetBranchAddress("l2_footprintCorrection", &l2_footprintCorrection, &b_l2_footprintCorrection);
   fChain->SetBranchAddress("l2_neutralIsoPtSum", &l2_neutralIsoPtSum, &b_l2_neutralIsoPtSum);
   fChain->SetBranchAddress("l2_puCorrPtSum", &l2_puCorrPtSum, &b_l2_puCorrPtSum);
   fChain->SetBranchAddress("l2_photonPtSumOutsideSignalCone", &l2_photonPtSumOutsideSignalCone, &b_l2_photonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("l2_byTightIsolationMVArun2v1DBoldDMwLT", &l2_byTightIsolationMVArun2v1DBoldDMwLT, &b_l2_byTightIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("l1_gen_pt", &l1_gen_pt, &b_l1_gen_pt);
   fChain->SetBranchAddress("l1_gen_eta", &l1_gen_eta, &b_l1_gen_eta);
   fChain->SetBranchAddress("l1_gen_phi", &l1_gen_phi, &b_l1_gen_phi);
   fChain->SetBranchAddress("l1_gen_charge", &l1_gen_charge, &b_l1_gen_charge);
   fChain->SetBranchAddress("l1_gen_mass", &l1_gen_mass, &b_l1_gen_mass);
   fChain->SetBranchAddress("l1_gen_pdgId", &l1_gen_pdgId, &b_l1_gen_pdgId);
   fChain->SetBranchAddress("l2_gen_pt", &l2_gen_pt, &b_l2_gen_pt);
   fChain->SetBranchAddress("l2_gen_eta", &l2_gen_eta, &b_l2_gen_eta);
   fChain->SetBranchAddress("l2_gen_phi", &l2_gen_phi, &b_l2_gen_phi);
   fChain->SetBranchAddress("l2_gen_charge", &l2_gen_charge, &b_l2_gen_charge);
   fChain->SetBranchAddress("l2_gen_mass", &l2_gen_mass, &b_l2_gen_mass);
   fChain->SetBranchAddress("l2_gen_pdgId", &l2_gen_pdgId, &b_l2_gen_pdgId);
   fChain->SetBranchAddress("l1_gen_vis_pt", &l1_gen_vis_pt, &b_l1_gen_vis_pt);
   fChain->SetBranchAddress("l1_gen_vis_eta", &l1_gen_vis_eta, &b_l1_gen_vis_eta);
   fChain->SetBranchAddress("l1_gen_vis_phi", &l1_gen_vis_phi, &b_l1_gen_vis_phi);
   fChain->SetBranchAddress("l1_gen_vis_charge", &l1_gen_vis_charge, &b_l1_gen_vis_charge);
   fChain->SetBranchAddress("l1_gen_vis_mass", &l1_gen_vis_mass, &b_l1_gen_vis_mass);
   fChain->SetBranchAddress("l2_gen_vis_pt", &l2_gen_vis_pt, &b_l2_gen_vis_pt);
   fChain->SetBranchAddress("l2_gen_vis_eta", &l2_gen_vis_eta, &b_l2_gen_vis_eta);
   fChain->SetBranchAddress("l2_gen_vis_phi", &l2_gen_vis_phi, &b_l2_gen_vis_phi);
   fChain->SetBranchAddress("l2_gen_vis_charge", &l2_gen_vis_charge, &b_l2_gen_vis_charge);
   fChain->SetBranchAddress("l2_gen_vis_mass", &l2_gen_vis_mass, &b_l2_gen_vis_mass);
   fChain->SetBranchAddress("l1_gen_decaymode", &l1_gen_decaymode, &b_l1_gen_decaymode);
   fChain->SetBranchAddress("l2_gen_decaymode", &l2_gen_decaymode, &b_l2_gen_decaymode);
   fChain->SetBranchAddress("l1_trigger_weight", &l1_trigger_weight, &b_l1_trigger_weight);
   fChain->SetBranchAddress("l1_trigger_weight_up", &l1_trigger_weight_up, &b_l1_trigger_weight_up);
   fChain->SetBranchAddress("l1_trigger_weight_down", &l1_trigger_weight_down, &b_l1_trigger_weight_down);
   fChain->SetBranchAddress("l2_trigger_weight", &l2_trigger_weight, &b_l2_trigger_weight);
   fChain->SetBranchAddress("l2_trigger_weight_up", &l2_trigger_weight_up, &b_l2_trigger_weight_up);
   fChain->SetBranchAddress("l2_trigger_weight_down", &l2_trigger_weight_down, &b_l2_trigger_weight_down);
   fChain->SetBranchAddress("mt2", &mt2, &b_mt2);
   fChain->SetBranchAddress("mt2_lep", &mt2_lep, &b_mt2_lep);
   fChain->SetBranchAddress("mt2_mvamet", &mt2_mvamet, &b_mt2_mvamet);
   fChain->SetBranchAddress("mt2_rawpfmet", &mt2_rawpfmet, &b_mt2_rawpfmet);
   fChain->SetBranchAddress("minDphiMETJets", &minDphiMETJets, &b_minDphiMETJets);
   fChain->SetBranchAddress("trigger_ditau35", &trigger_ditau35, &b_trigger_ditau35);
   fChain->SetBranchAddress("trigger_ditau35_combiso", &trigger_ditau35_combiso, &b_trigger_ditau35_combiso);
   fChain->SetBranchAddress("trigger_singletau140", &trigger_singletau140, &b_trigger_singletau140);
   fChain->SetBranchAddress("trigger_singletau120", &trigger_singletau120, &b_trigger_singletau120);
   fChain->SetBranchAddress("trigger_matched_ditau35", &trigger_matched_ditau35, &b_trigger_matched_ditau35);
   fChain->SetBranchAddress("trigger_matched_ditau35_combiso", &trigger_matched_ditau35_combiso, &b_trigger_matched_ditau35_combiso);
   fChain->SetBranchAddress("trigger_matched_singletau140", &trigger_matched_singletau140, &b_trigger_matched_singletau140);
   fChain->SetBranchAddress("trigger_matched_singletau120", &trigger_matched_singletau120, &b_trigger_matched_singletau120);
   Notify();
}

Bool_t treeBase_data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treeBase_data::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treeBase_data::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void treeBase_data::Loop()
{
//   In a ROOT session, you can do:
//      root> .L treeBase_data.C
//      root> treeBase_data t
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
      // if (Cut(ientry) < 0) continue;
   }
}
#endif // #ifdef treeBase_data_cxx
