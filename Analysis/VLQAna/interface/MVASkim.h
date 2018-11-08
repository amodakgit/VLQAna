#ifndef __MVASkim__h
#define __MVASkim__h

#include <fstream>
#include <string>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

class TTree;
class TFile;

typedef struct  
{
  float lepEta_Mu;
  float lepPt_Mu;
  float lepPhi_Mu;
  float leadJetEta_Mu;
  float leadJetPt_Mu;
  float leadJetPhi_Mu;
  float met_Mu;
  float ST_Mu;
  float HT_Mu;
  float NBTags_Mu;
  float NJets_Mu;
  float Mass_Mu;
  float DPHI_LepJet_Mu;
  float DPHI_MetLep_Mu;
  float FwdJetEta_Mu;
  float FwdJetPt_Mu;
  float FwdJetPhi_Mu;
  float MT_Mu;
  float DR_LepCloJet_Mu;
  float bVsW_ratio_Mu;
  float Ext_Jet_TransPt_Mu;
  float HTNoLead_Mu;
  float Angle_MuZ_Jet_Mu;
  float Angle_MuJet_Met_Mu;
  float Evtwt_Mu;

  float lepEta_Ele;
  float lepPt_Ele;
  float lepPhi_Ele;
  float leadJetEta_Ele;
  float leadJetPt_Ele;
  float leadJetPhi_Ele;
  float met_Ele;
  float ST_Ele;
  float HT_Ele;
  float NBTags_Ele;
  float NJets_Ele;
  float Mass_Ele;
  float DPHI_LepJet_Ele;
  float DPHI_MetLep_Ele;
  float FwdJetEta_Ele;
  float FwdJetPt_Ele;
  float FwdJetPhi_Ele;
  float MT_Ele;
  float DR_LepCloJet_Ele;
  float bVsW_ratio_Ele;
  float Ext_Jet_TransPt_Ele;
  float HTNoLead_Ele;
  float Angle_MuZ_Jet_Ele;
  float Angle_MuJet_Met_Ele;
  float Evtwt_Ele;
} TreeVariables;

class MVASkim {
    
 public:

  MVASkim(const std::string& filename);
  virtual ~MVASkim();

  void fill(const TreeVariables& varList);
  void close();

  //TFileDirectory* _mvaFile;
  TTree* _tree;

  TreeVariables _varList;
};
#endif
