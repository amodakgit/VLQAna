#include <iostream>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "Analysis/VLQAna/interface/MVASkim.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

using std::string;
using std::cout;
using std::endl;

MVASkim::MVASkim(const string& filename) {
  edm::Service<TFileService> fs;
  TFileDirectory _mvaFile = fs->mkdir(filename.c_str());
  _tree = _mvaFile.make<TTree>("RTree", "RTree");
  _tree->Branch("lepEta_Mu",             &_varList.lepEta_Mu,           "lepEta_Mu/F");
  _tree->Branch("lepPt_Mu",              &_varList.lepPt_Mu,            "lepPt_Mu/F");
  _tree->Branch("lepPhi_Mu",             &_varList.lepPhi_Mu,           "lepPhi_Mu/F");
  _tree->Branch("leadJetEta_Mu",         &_varList.leadJetEta_Mu,       "leadJetEta_Mu/F");
  _tree->Branch("leadJetPt_Mu",          &_varList.leadJetPt_Mu,        "leadJetPt_Mu/F");
  _tree->Branch("leadJetPhi_Mu",         &_varList.leadJetPhi_Mu,       "leadJetPhi_Mu/F");
  _tree->Branch("met_Mu",                &_varList.met_Mu,              "met_Mu/F");
  _tree->Branch("ST_Mu",                 &_varList.ST_Mu,               "ST_Mu/F");
  _tree->Branch("HT_Mu",                 &_varList.HT_Mu,               "HT_Mu/F");
  _tree->Branch("NBTags_Mu",             &_varList.NBTags_Mu,           "NBTags_Mu/F");
  _tree->Branch("NJets_Mu",              &_varList.NJets_Mu,            "NJets_Mu/F");
  _tree->Branch("Mass_Mu",               &_varList.Mass_Mu,             "Mass_Mu/F");
  _tree->Branch("DPHI_LepJet_Mu",        &_varList.DPHI_LepJet_Mu,      "DPHI_LepJet_Mu/F");
  _tree->Branch("DPHI_Metlep_Mu",        &_varList.DPHI_MetLep_Mu,      "DPHI_MetLep_Mu/F");
  _tree->Branch("FwdJetEta_Mu",          &_varList.FwdJetEta_Mu,        "FwdJetEta_Mu/F");
  _tree->Branch("FwdJetPt_Mu",           &_varList.FwdJetPt_Mu,         "FwdJetPt_Mu/F");
  _tree->Branch("FwdJetPhi_Mu",          &_varList.FwdJetPhi_Mu,        "FwdJetPhi_Mu/F");
  _tree->Branch("MT_Mu",                 &_varList.MT_Mu,               "MT_Mu/F");
  _tree->Branch("DR_LepCloJet_Mu",       &_varList.DR_LepCloJet_Mu,     "DR_LepCloJet_Mu/F");
  _tree->Branch("bVsW_ratio_Mu",         &_varList.bVsW_ratio_Mu,       "bVsW_ratio_Mu/F");
  _tree->Branch("Ext_Jet_TransPt_Mu",    &_varList.Ext_Jet_TransPt_Mu,  "Ext_Jet_TransPt_Mu/F");
  _tree->Branch("HTNoLead_Mu",           &_varList.HTNoLead_Mu,         "HTNoLead_Mu/F");
  _tree->Branch("Angle_MuZ_Jet_Mu",      &_varList.Angle_MuZ_Jet_Mu,    "Angle_MuZ_Jet_Mu/F");
  _tree->Branch("Angle_MuJet_Met_Mu",    &_varList.Angle_MuJet_Met_Mu,  "Angle_MuJet_Met_Mu/F");
  _tree->Branch("Evtwt_Mu",              &_varList.Evtwt_Mu,            "Evtwt_Mu/F");

  _tree->Branch("lepEta_Ele",             &_varList.lepEta_Ele,           "lepEta_Ele/F");
  _tree->Branch("lepPt_Ele",              &_varList.lepPt_Ele,            "lepPt_Ele/F");
  _tree->Branch("lepPhi_Ele",             &_varList.lepPhi_Ele,           "lepPhi_Ele/F");
  _tree->Branch("leadJetEta_Ele",         &_varList.leadJetEta_Ele,       "leadJetEta_Ele/F");
  _tree->Branch("leadJetPt_Ele",          &_varList.leadJetPt_Ele,        "leadJetPt_Ele/F");
  _tree->Branch("leadJetPhi_Ele",         &_varList.leadJetPhi_Ele,       "leadJetPhi_Ele/F");
  _tree->Branch("met_Ele",                &_varList.met_Ele,              "met_Ele/F");
  _tree->Branch("ST_Ele",                 &_varList.ST_Ele,               "ST_Ele/F");
  _tree->Branch("HT_Ele",                 &_varList.HT_Ele,               "HT_Ele/F");
  _tree->Branch("NBTags_Ele",             &_varList.NBTags_Ele,           "NBTags_Ele/F");
  _tree->Branch("NJets_Ele",              &_varList.NJets_Ele,            "NJets_Ele/F");
  _tree->Branch("Mass_Ele",               &_varList.Mass_Ele,             "Mass_Ele/F");
  _tree->Branch("DPHI_LepJet_Ele",        &_varList.DPHI_LepJet_Ele,      "DPHI_LepJet_Ele/F");
  _tree->Branch("DPHI_Metlep_Ele",        &_varList.DPHI_MetLep_Ele,      "DPHI_MetLep_Ele/F");
  _tree->Branch("FwdJetEta_Ele",          &_varList.FwdJetEta_Ele,        "FwdJetEta_Ele/F");
  _tree->Branch("FwdJetPt_Ele",           &_varList.FwdJetPt_Ele,         "FwdJetPt_Ele/F");
  _tree->Branch("FwdJetPhi_Ele",          &_varList.FwdJetPhi_Ele,        "FwdJetPhi_Ele/F");
  _tree->Branch("MT_Ele",                 &_varList.MT_Ele,               "MT_Ele/F");
  _tree->Branch("DR_LepCloJet_Ele",       &_varList.DR_LepCloJet_Ele,     "DR_LepCloJet_Ele/F");
  _tree->Branch("bVsW_ratio_Ele",         &_varList.bVsW_ratio_Ele,       "bVsW_ratio_Ele/F");
  _tree->Branch("Ext_Jet_TransPt_Ele",    &_varList.Ext_Jet_TransPt_Ele,  "Ext_Jet_TransPt_Ele/F");
  _tree->Branch("HTNoLead_Ele",           &_varList.HTNoLead_Ele,         "HTNoLead_Ele/F");
  _tree->Branch("Angle_MuZ_Jet_Ele",      &_varList.Angle_MuZ_Jet_Ele,    "Angle_MuZ_Jet_Ele/F");
  _tree->Branch("Angle_MuJet_Met_Ele",    &_varList.Angle_MuJet_Met_Ele,  "Angle_MuJet_Met_Ele/F");
  _tree->Branch("Evtwt_Ele",              &_varList.Evtwt_Ele,            "Evtwt_Ele/F");

  fs->file().ls();
}

MVASkim::~MVASkim() {
  //delete _mvaFile;
}

void MVASkim::fill(const TreeVariables& varList) {
  memcpy(&_varList, &varList, sizeof(varList));
  //fs->file().cd();
  //_mvaFile->cd();
  _tree->Fill();  
}

void MVASkim::close() {
  //_mvaFile->cd();
  //_tree->Print();
  //_tree->Write();
  //_mvaFile->Write();
  //_mvaFile->Close();
}

