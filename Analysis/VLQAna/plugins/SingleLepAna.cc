// -*- C++ -*-
//
// Package:    Analysis/VLQAna
// Class:      SingleLepAna
//
/**\class VLQAna SingleLepAna.cc Analysis/VLQAna/plugins/SingleLepAna.cc
 Description: [one line class summary]
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Devdatta Majumder
//         Created:  Fri, 27 Feb 2015 16:09:10 GMT
// Modified: Sadia Khalil
//           25 Mar 2016 17:11 CDT
//

#include <iostream>
#include <memory>
#include <vector>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "AnalysisDataFormats/BoostedObjects/interface/GenParticleWithDaughters.h"
#include "AnalysisDataFormats/BoostedObjects/interface/ResolvedVjj.h"

#include "Analysis/VLQAna/interface/Utilities.h"
#include "Analysis/VLQAna/interface/DileptonCandsProducer.h"
#include "Analysis/VLQAna/interface/CandidateFilter.h"
#include "Analysis/VLQAna/interface/MuonMaker.h"
#include "Analysis/VLQAna/interface/ElectronMaker.h"
#include "Analysis/VLQAna/interface/JetMaker.h"
#include "Analysis/VLQAna/interface/HT.h"
#include "Analysis/VLQAna/interface/ApplyLeptonIDSFs.h"
#include "Analysis/VLQAna/interface/ApplyLeptonTrigSFs.h"
#include "Analysis/VLQAna/interface/CandidateCleaner.h"
#include "Analysis/VLQAna/interface/METMaker.h"
#include "Analysis/VLQAna/interface/PickGenPart.h"
#include "Analysis/VLQAna/interface/JetID.h"
#include "Analysis/VLQAna/interface/MassReco.h"
#include "Analysis/VLQAna/interface/BTagSFUtils.h"

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TFitResult.h>

//
// class declaration
//

class SingleLepAna : public edm::EDFilter {
public:
    explicit SingleLepAna(const edm::ParameterSet&);
    ~SingleLepAna();
    int cutFlow(const edm::ParameterSet & param);
    int hltDecision(const edm::ParameterSet & param);
    void printing(std::string str){std::cout << str<<"\n";}
    
private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    void fillAdditionalPlots( vlq::ElectronCollection goodElectrons,double evtwt);
    double GetDYNLOCorr(const double dileppt);
    double ZptCorr(vlq::Candidate, double, double);
    double htCorr(double ht, double p0, double p1);
    bool solve_nu(const TLorentzVector &vlep, const TLorentzVector &vnu, double wmass, double& nuz1, double& nuz2);
    void adjust_e_for_mass(TLorentzVector& v, double mass);
    double getWjetHTRewSF(double ht, double ht_low_threshold);
    
    // ----------member data ---------------------------
    edm::EDGetTokenT<string>   t_evttype         ;
    edm::EDGetTokenT<int>      t_evtno           ;
    edm::EDGetTokenT<int>      t_lumi            ;
    edm::EDGetTokenT<int>      t_runno           ;
    edm::EDGetTokenT<double>   t_evtwtGen        ;
    edm::EDGetTokenT<double>   t_evtwtPV         ;
    edm::EDGetTokenT<double>   t_evtwtPVLow      ;
    edm::EDGetTokenT<double>   t_evtwtPVHigh     ;
    edm::EDGetTokenT<unsigned> t_npv             ;
    edm::EDGetTokenT<bool>     t_hltdecision_ele ;
    edm::EDGetTokenT<bool>     t_hltdecision_mu  ;
    edm::EDGetTokenT<double>   t_toppt           ;
    edm::ParameterSet GenHSelParams_             ;
    edm::ParameterSet genParams_                 ;
    const bool sys_                              ;
    bool btagsf_bcUp_                            ;
    bool btagsf_bcDown_                          ;
    bool btagsf_lUp_                             ;
    bool btagsf_lDown_                           ;
    const bool PileupUp_                         ;
    const bool PileupDown_                       ;
    const bool applyLeptonIDSFs_                 ;
    const bool applyLeptonTrigSFs_               ;
    const bool applyBTagSFs_                     ;
    const bool applyHtCorr_                      ;
    const bool applyTopPtCorr_                   ;
    const bool applyWjetHtCorr_                  ;
    const double jecShift_                       ;
    const double jerShift_                       ;
    const bool topPtRewUp_                       ;
    const bool topPtRewDown_                     ;
    const double wjetHtShift_                    ;
    ApplyLeptonIDSFs muidsf                      ;
    ApplyLeptonTrigSFs mutrigsf                  ;
    ApplyLeptonIDSFs elidsf                      ;
    ApplyLeptonTrigSFs eltrigsf                  ;
    METMaker metmaker                            ;
    MuonMaker muonmaker                          ;
    ElectronMaker electronmaker                  ;
    MuonMaker muonmakerloose                     ;
    MuonMaker muonmakerskim                      ;
    ElectronMaker electronmakerloose             ;
    JetMaker jetAK4maker                         ;
    JetMaker jetAK4BTaggedmaker                  ;
    JetMaker jetAK8maker                         ;
    JetMaker jetWTaggedmaker                     ;
    edm::Service<TFileService> fs                ;
    std::map<std::string, TH1D*> h1_             ;
    std::map<std::string, TH2D*> h2_             ;
    const std::string fnamebtagSF_               ;
    const std::string btageffmap_                ;
    std::unique_ptr<BTagSFUtils> btagsfutils_    ;
    vlq::ElectronCollection signalElectrons      ;
    vlq::MuonCollection signalMuons              ;
    vlq::MuonCollection skimMuons                ;
    vlq::JetCollection goodAK4Jets               ;
    vlq::JetCollection goodAK8Jets               ;
    vlq::MetCollection goodMet                   ;
    edm::ParameterSet signalSelection_           ;
    edm::ParameterSet wjetSelection_             ;
    edm::ParameterSet ttbarSelection_            ;
    edm::ParameterSet crSelection_               ;
};

using namespace std;

// constructors and destructor
SingleLepAna::SingleLepAna(const edm::ParameterSet& iConfig) :
t_evttype               (consumes<string> (iConfig.getParameter<edm::InputTag>("evttype"))),
t_evtno                 (consumes<int> (iConfig.getParameter<edm::InputTag>("evtno"))),
t_lumi                  (consumes<int> (iConfig.getParameter<edm::InputTag>("lumi"))),
t_runno                 (consumes<int> (iConfig.getParameter<edm::InputTag>("runno"))),
t_evtwtGen              (consumes<double> (iConfig.getParameter<edm::InputTag>("evtwtGen"))),
t_evtwtPV               (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPV"))),
t_evtwtPVLow            (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPVLow"))),
t_evtwtPVHigh           (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPVHigh"))),
t_npv                   (consumes<unsigned>(iConfig.getParameter<edm::InputTag>("npv"))),
t_hltdecision_ele       (consumes<bool>    (iConfig.getParameter<edm::InputTag>("hltdecisionEle"))),
t_hltdecision_mu        (consumes<bool>    (iConfig.getParameter<edm::InputTag>("hltdecisionMu"))),
t_toppt                 (consumes<double>  (iConfig.getParameter<edm::InputTag>("toppt"))),
GenHSelParams_          (iConfig.getParameter<edm::ParameterSet> ("GenHSelParams")),
sys_                    (iConfig.getParameter<bool>              ("sys")),
btagsf_bcUp_            (iConfig.getParameter<bool>              ("btagsf_bcUp")),
btagsf_bcDown_          (iConfig.getParameter<bool>              ("btagsf_bcDown")),
btagsf_lUp_             (iConfig.getParameter<bool>              ("btagsf_lUp")),
btagsf_lDown_           (iConfig.getParameter<bool>              ("btagsf_lDown")),
PileupUp_               (iConfig.getParameter<bool>              ("PileupUp")),
PileupDown_             (iConfig.getParameter<bool>              ("PileupDown")),
applyLeptonIDSFs_	(iConfig.getParameter<bool>              ("applyLeptonIDSFs")),
applyLeptonTrigSFs_     (iConfig.getParameter<bool>              ("applyLeptonTrigSFs")),
applyBTagSFs_           (iConfig.getParameter<bool>              ("applyBTagSFs")),
applyHtCorr_            (iConfig.getParameter<bool>              ("applyHtCorr")),
applyTopPtCorr_         (iConfig.getParameter<bool>              ("applyTopPtCorr")),
applyWjetHtCorr_        (iConfig.getParameter<bool>              ("applyWjetHtCorr")),
jecShift_               (iConfig.getParameter<double>            ("jecShift")),
jerShift_               (iConfig.getParameter<double>            ("jerShift")),
topPtRewUp_             (iConfig.getParameter<bool>              ("topPtRewUp")),
topPtRewDown_           (iConfig.getParameter<bool>              ("topPtRewDown")),
wjetHtShift_            (iConfig.getParameter<double>            ("wjetHtShift")),
muidsf                  (iConfig.getParameter<edm::ParameterSet> ("muIdSFsParams")),
mutrigsf                (iConfig.getParameter<edm::ParameterSet> ("muTrigSFsParams")),
elidsf                  (iConfig.getParameter<edm::ParameterSet> ("elIdSFsParams")),
eltrigsf                (iConfig.getParameter<edm::ParameterSet> ("elTrigSFsParams")),
metmaker                (iConfig.getParameter<edm::ParameterSet> ("metselParams"),consumesCollector()),
muonmaker               (iConfig.getParameter<edm::ParameterSet> ("muselParams"),consumesCollector()),
electronmaker           (iConfig.getParameter<edm::ParameterSet> ("elselParams"),consumesCollector()),
muonmakerloose          (iConfig.getParameter<edm::ParameterSet> ("musellooseParams"),consumesCollector()),
muonmakerskim           (iConfig.getParameter<edm::ParameterSet> ("muselskimParams"),consumesCollector()),
electronmakerloose      (iConfig.getParameter<edm::ParameterSet> ("elsellooseParams"),consumesCollector()),
jetAK4maker             (iConfig.getParameter<edm::ParameterSet> ("jetAK4selParams"),consumesCollector()),
jetAK4BTaggedmaker      (iConfig.getParameter<edm::ParameterSet> ("jetAK4BTaggedselParams"),consumesCollector()),
jetAK8maker             (iConfig.getParameter<edm::ParameterSet> ("jetAK8selParams"),consumesCollector()),
jetWTaggedmaker         (iConfig.getParameter<edm::ParameterSet> ("jetWTaggedselParams"),consumesCollector()),
fnamebtagSF_            (iConfig.getParameter<std::string>       ("fnamebtagSF")),
btageffmap_             (iConfig.getParameter<std::string>      ("btageffmap")),
//btagsfutils_            (new BTagSFUtils(fnamebtagSF_,BTagEntry::OP_MEDIUM,30., 670., 30., 670., 20., 1000.)),
btagsfutils_            (new BTagSFUtils(fnamebtagSF_,BTagEntry::OP_MEDIUM,20., 1000., 20., 1000., 20., 1000., btageffmap_)),
signalSelection_        (iConfig.getParameter<edm::ParameterSet> ("signalSelection")),
wjetSelection_          (iConfig.getParameter<edm::ParameterSet> ("wjetSelection")),
ttbarSelection_         (iConfig.getParameter<edm::ParameterSet> ("ttbarSelection")),
crSelection_            (iConfig.getParameter<edm::ParameterSet> ("crSelection"))
{
    cout <<"constructor finished" << endl;
}


SingleLepAna::~SingleLepAna() {}

bool SingleLepAna::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
    using namespace edm;
    
    // set all objects vectors to zero
    
    signalElectrons.clear();
    signalMuons.clear();
    skimMuons.clear();
    goodMet.clear();
    goodAK4Jets.clear();
    goodAK8Jets.clear();
    
    
    Handle<string>   h_evttype          ; evt.getByToken(t_evttype        , h_evttype     ) ;
    Handle<int>      h_evtno            ; evt.getByToken(t_evtno        , h_evtno     ) ;
    Handle<int>      h_lumi             ; evt.getByToken(t_lumi        , h_lumi     ) ;
    Handle<int>      h_runno             ; evt.getByToken(t_runno        , h_runno     ) ;
    Handle<double>   h_evtwtGen         ; evt.getByToken(t_evtwtGen       , h_evtwtGen    ) ;
    Handle<double>   h_evtwtPV          ; evt.getByToken(t_evtwtPV        , h_evtwtPV     ) ;
    Handle<double>   h_evtwtPVLow       ; evt.getByToken(t_evtwtPVLow     , h_evtwtPVLow  ) ;
    Handle<double>   h_evtwtPVHigh      ; evt.getByToken(t_evtwtPVHigh    , h_evtwtPVHigh ) ;
    Handle<unsigned> h_npv              ; evt.getByToken(t_npv            , h_npv         ) ;
    Handle<bool>     h_hltdecision_mu   ; evt.getByToken(t_hltdecision_mu , h_hltdecision_mu ) ;
    Handle<bool>     h_hltdecision_ele  ; evt.getByToken(t_hltdecision_ele, h_hltdecision_ele ) ;
    Handle<double>   h_toppt            ; evt.getByToken(t_toppt          , h_toppt       ) ;

    
    h1_["truePU"]->Fill(*h_evtwtPV.product());
    h1_["nPU"]   ->Fill(*h_npv.product());
    
    double evtwt = 1;
    
     if (PileupUp_)
     evtwt = (*h_evtwtGen.product()) * (*h_evtwtPVHigh.product()) ;
     else if (PileupDown_)
     evtwt = (*h_evtwtGen.product()) * (*h_evtwtPVLow.product()) ;
     else
     evtwt = (*h_evtwtGen.product()) * (*h_evtwtPV.product()) ;
    
    //// Single lepton analysis
    
    vlq::MuonCollection looseMuons;
    muonmakerloose (evt, looseMuons) ;
    muonmakerskim  (evt, skimMuons) ;
    muonmaker(evt, signalMuons) ;
    
    vlq::ElectronCollection looseElectrons;
    electronmakerloose(evt, looseElectrons) ;
    electronmaker(evt, signalElectrons) ;
    
    jetAK4maker(evt, goodAK4Jets);
    jetAK8maker(evt, goodAK8Jets);
    

    CandidateCleaner cleanjets(0.4);
    //cleanjets(goodAK4Jets, skimMuons);
    cleanjets(goodAK4Jets, signalMuons);
    cleanjets(goodAK4Jets, signalElectrons);
    
    vlq::JetCollection goodBTaggedAK4Jets;
    jetAK4BTaggedmaker(evt, goodBTaggedAK4Jets) ;
    
    cleanjets(goodBTaggedAK4Jets, signalMuons);
    cleanjets(goodBTaggedAK4Jets, signalElectrons);
    
    metmaker(evt, goodMet) ;
    
    HT htak4(goodAK4Jets) ;

    
    // single lepton final state
    
    // apply filters
    
    
    double btagsf(1) ;
    double btagsf_bcUp(1) ;
    double btagsf_bcDown(1) ;
    double btagsf_lUp(1) ;
    double btagsf_lDown(1) ;
    if ( applyBTagSFs_ && *h_evttype.product() != "EvtType_Data") {
        std::vector<double>csvs;
        std::vector<double>pts;
        std::vector<double>etas;
        std::vector<int>   flhads;
        
        for (vlq::Jet jet : goodAK4Jets) {
        //for (vlq::Jet jet : goodBTaggedAK4Jets) {
            csvs.push_back(jet.getCSV()) ;
            pts.push_back(jet.getPt()) ;
            etas.push_back(jet.getEta()) ;
            flhads.push_back(jet.getHadronFlavour()) ;
        }
        
        btagsfutils_->getBTagSFs(csvs, pts, etas, flhads, jetAK4BTaggedmaker.idxjetCSVDiscMin_, btagsf, btagsf_bcUp, btagsf_bcDown, btagsf_lUp, btagsf_lDown);
        if (btagsf_bcUp_)
            evtwt *= btagsf_bcUp;
        else if (btagsf_bcDown_)
            evtwt *= btagsf_bcDown;
        else  if (btagsf_lUp_)
            evtwt *= btagsf_lUp;
        else if (btagsf_lDown_)
            evtwt *= btagsf_lDown;
        else
            evtwt *= btagsf;
        
    }
    
    // apply top pt reweighting
    if ( applyTopPtCorr_ && *h_evttype.product() != "EvtType_Data") {
        
        double toppt_weight = *h_toppt.product();
        if(topPtRewUp_){
            evtwt *=  1.0;
        }
        else if(topPtRewDown_){
            evtwt *= toppt_weight * toppt_weight;
        }
        else{
            evtwt *= toppt_weight;
        }
    }
    
    // apply w+jet ht reweighting
    if ( applyWjetHtCorr_ &&  *h_evttype.product() != "EvtType_Data") {
        
        double wjet_rew = getWjetHTRewSF(htak4.getHT(), 700.);
        if(wjetHtShift_ == 1){
            evtwt *=  wjet_rew;
        }
        else if(wjetHtShift_ == 2){
            evtwt *= 1.0;
        }
        else if(wjetHtShift_ == -1){
            evtwt *= wjet_rew * wjet_rew;
        }
    }
    
    //unsigned int evtno = abs(*h_evtno.product());
    //unsigned int lumi = abs(*h_lumi.product());
    //unsigned int run = abs(*h_runno.product());

    if( (signalElectrons.size() + signalMuons.size()) == 1 ){
        
        string ch;
        TLorentzVector wbana_lep_p4;
        if(signalElectrons.size() == 1){
            wbana_lep_p4 = signalElectrons.at(0).getP4();
            ch = "_Ele";
            
            const bool hltdecision(*h_hltdecision_ele.product()) ;
            if ( !hltdecision && (*h_evttype.product()) != "EvtType_MC" ) return false;

            if((*h_evttype.product()) == "EvtType_MC" && applyLeptonIDSFs_){
                evtwt *= elidsf.IDSF(signalElectrons.at(0).getPt(),signalElectrons.at(0).getscEta());
                evtwt *= elidsf.RECOSF(signalElectrons.at(0).getPt(),signalElectrons.at(0).getscEta());
            }
            if((*h_evttype.product()) == "EvtType_MC" && applyLeptonTrigSFs_){
                evtwt *= eltrigsf(signalElectrons.at(0).getPt(),signalElectrons.at(0).getscEta());
            }
	}
	else if(signalMuons.size() == 1){
            wbana_lep_p4 = signalMuons.at(0).getP4();
            ch = "_Mu";

            const bool hltdecision(*h_hltdecision_mu.product()) ;
            // ignor muon events with pt > 500 GeV (keep this untill we know that the trigger works)
            //if(wbana_lep_p4.Pt() > 500) return false;

            if ( !hltdecision && (*h_evttype.product()) != "EvtType_MC") return false;

            if((*h_evttype.product()) == "EvtType_MC" && applyLeptonIDSFs_){
                evtwt *= muidsf.IDSF(signalMuons.at(0).getPt(),signalMuons.at(0).getEta());
                evtwt *= muidsf.IsoSF(signalMuons.at(0).getPt(),signalMuons.at(0).getEta());
            }
            if((*h_evttype.product()) == "EvtType_MC" && applyLeptonTrigSFs_){
                evtwt *= mutrigsf(signalMuons.at(0).getPt(),signalMuons.at(0).getEta());
            }
	}
	else{
	}
        
        double dr_ak4jet_lep = 999, fwd_ak4jet_eta = 0;
        TLorentzVector fwd_ak4jet_p4;
        
        
        for (unsigned int jet = 0; jet < goodAK4Jets.size(); jet++) {
            TLorentzVector jet_p4 = goodAK4Jets.at(jet).getP4();
            
            // calculate min dr between lepton and the closest jet
            double dr = jet_p4.DeltaR(wbana_lep_p4);
            if(dr < dr_ak4jet_lep && jet_p4.Pt() > 40){
                dr_ak4jet_lep = dr;
            }
            
            // find the forward jet
            if( fabs( jet_p4.Eta() ) > fwd_ak4jet_eta){
                fwd_ak4jet_eta = fabs( jet_p4.Eta() );
                fwd_ak4jet_p4  = jet_p4;
            }
        }
        
        
        // ak8 jets
        
        double dr_ak4_ak8 = 0.4;
        TLorentzVector leading_ak8jet_p4;
        for (unsigned int jet = 0; jet < goodAK8Jets.size() && goodAK4Jets.size() > 0; jet++) {
            TLorentzVector jet_p4 = goodAK8Jets.at(jet).getP4();
            
            // find the ak8 jet that matches to the leading ak4 jet
            
            double dr = jet_p4.DeltaR(goodAK4Jets.at(0).getP4());
            if(dr < dr_ak4_ak8){
                dr_ak4_ak8 = dr;
                leading_ak8jet_p4 = jet_p4;
            }
        }
        
        
        // analysis variables
        
        // define mt
        TLorentzVector met_p4 = goodMet.at(0).getP4();
        double mt = TMath::Sqrt( 2*wbana_lep_p4.Pt() * met_p4.Pt() * ( 1 - TMath::Cos(wbana_lep_p4.DeltaPhi(met_p4) ) ) );
        
        double st = 0, st_v2 = 0, st_nomet = 0, mass = 0, mass_v2 = 0, dphi = 999, dphi_metjet= 999, dphi_metlep = 999, angle_MuZ_Jet = 999, angle_MuJet_Met = 999;
        if(goodAK4Jets.size() > 0){
            
            // define st obtained using ak8
            st      = met_p4.Pt() + leading_ak8jet_p4.Pt() + wbana_lep_p4.Pt();

            // define st obtained using ak8
            st_nomet   = leading_ak8jet_p4.Pt() + wbana_lep_p4.Pt();
            
            // define st_v2 obtained using ak4
            st_v2   = met_p4.Pt() + goodAK4Jets.at(0).getP4().Pt() + wbana_lep_p4.Pt();

            // define mass obtained using ak8
            mass    = (met_p4 + leading_ak8jet_p4 + wbana_lep_p4).M();
            
            // define mass obtained using ak4
            mass_v2 = (met_p4 + goodAK4Jets.at(0).getP4() + wbana_lep_p4).M();
            
            // define dphi
            dphi = goodAK4Jets.at(0).getP4().DeltaPhi(wbana_lep_p4);

            // define met-jet dphi
            dphi_metjet = goodAK4Jets.at(0).getP4().DeltaPhi(met_p4);

            // define met-lep dphi
            dphi_metlep = wbana_lep_p4.DeltaPhi(met_p4);

            TVector3 Mu_TV3(wbana_lep_p4.Px(), wbana_lep_p4.Py(), wbana_lep_p4.Pz());
            TVector3 Jet_TV3(goodAK4Jets.at(0).getP4().Px(), goodAK4Jets.at(0).getP4().Py(), goodAK4Jets.at(0).getP4().Pz());
            TVector3 Z_TV3(0, 0, 1);
            TVector3 Met_TV3(met_p4.Px(), met_p4.Py(), 0);
            TVector3 muZ = Z_TV3.Cross(Mu_TV3);
            TVector3 muJet = Jet_TV3.Cross(Mu_TV3);
            angle_MuZ_Jet = muZ.Angle(Jet_TV3);
            angle_MuJet_Met = muJet.Angle(Met_TV3);
	    //std::cout << "angle_MuZ_Jet: " << angle_MuZ_Jet << ", angle_MuJet_Met: " << angle_MuJet_Met << std::endl;
        }
        

        //Histograms for b-tag efficiency map, after a preselection
        if (goodAK4Jets.size() > 0){
          for (vlq::Jet jet : goodAK4Jets) {
            if ( abs(jet.getHadronFlavour()) == 5) h2_["pt_eta_b_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
            else if ( abs(jet.getHadronFlavour()) == 4) h2_["pt_eta_c_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
            else if ( abs(jet.getHadronFlavour()) == 0) h2_["pt_eta_l_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
          }
          if ( goodBTaggedAK4Jets.size() > 0 ){
            for (vlq::Jet jet : goodBTaggedAK4Jets) {
              if ( abs(jet.getHadronFlavour()) == 5) h2_["pt_eta_b_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
              else if ( abs(jet.getHadronFlavour()) == 4) h2_["pt_eta_c_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
              else if ( abs(jet.getHadronFlavour()) == 0) h2_["pt_eta_l_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
            }
          }
        }
        
        // cutflows for the control and signal regions
        int cut_flow_ttbar   = cutFlow(ttbarSelection_);
        int cut_flow_cr      = cutFlow(crSelection_);
        int cut_flow_wjet    = cutFlow(wjetSelection_);
        int cut_flow_signal  = cutFlow(signalSelection_);
        
        // cutflow for the ttbar cr
        if(goodBTaggedAK4Jets.size() > 1){
            for(int cut = 0; cut < cut_flow_ttbar; cut++){
                
                string name = string("tt_")+to_string(cut) + ch;
                h1_[name+string("_MET")]    -> Fill(goodMet.at(0).getP4().Pt(), evtwt) ;
                h1_[name+string("_MT")]     -> Fill(mt, evtwt) ;
                h1_[name+string("_LepPt")]  -> Fill(wbana_lep_p4.Pt(), evtwt) ;
                h1_[name+string("_LepEta")] -> Fill(wbana_lep_p4.Eta(), evtwt) ;
                if(goodAK4Jets.size() > 0){
                    h1_[name+string("_JetPt")]  -> Fill(goodAK4Jets.at(0).getP4().Pt(), evtwt) ;
                    h1_[name+string("_JetEta")] -> Fill(goodAK4Jets.at(0).getP4().Eta(), evtwt) ;
                    h1_[name+string("_DR")]     -> Fill(dr_ak4jet_lep, evtwt) ;
                    h1_[name+string("_DPHI")]   -> Fill(dphi, evtwt) ;
                    h1_[name+string("_FwdJetPt")]   -> Fill(fwd_ak4jet_p4.Pt(), evtwt) ;
                    h1_[name+string("_FwdJetEta")]   -> Fill(fwd_ak4jet_p4.Eta(), evtwt) ;
                }
                
                h1_[name+string("_NJets")]    -> Fill(goodAK4Jets.size(), evtwt) ;
                h1_[name+string("_ST")]       -> Fill(st, evtwt) ;
                h1_[name+string("_STv2")]    -> Fill(st_v2, evtwt) ;
                h1_[name+string("_Mass")]     -> Fill(mass, evtwt) ;
                h1_[name+string("_Massv2")]  -> Fill(mass_v2, evtwt);
                h1_[name+string("_HT")]       -> Fill(htak4.getHT(), evtwt);
                h1_[name+string("_NBTags")]   -> Fill(goodBTaggedAK4Jets.size(), evtwt);
                h1_[name+string("_ST_nomet")]       -> Fill(st_nomet, evtwt) ;
                h1_[name+string("_DPHI_MetJet")]       -> Fill(dphi_metjet, evtwt) ;
                h1_[name+string("_DPHI_MetLep")]       -> Fill(dphi_metlep, evtwt) ;
                h1_[name+string("_Angle_MuZ_Jet")]       -> Fill(angle_MuZ_Jet, evtwt) ;
                h1_[name+string("_Angle_MuJet_Met")]       -> Fill(angle_MuJet_Met, evtwt) ;
            }
        }
        

        // cutflow for the alternative cr
        if(goodBTaggedAK4Jets.size() >= 1){
            for(int cut = 0; cut < cut_flow_cr; cut++){
                
                string name = string("cr_")+to_string(cut) + ch;
                h1_[name+string("_MET")]    -> Fill(goodMet.at(0).getP4().Pt(), evtwt) ;
                h1_[name+string("_MT")]     -> Fill(mt, evtwt) ;
                h1_[name+string("_LepPt")]  -> Fill(wbana_lep_p4.Pt(), evtwt) ;
                h1_[name+string("_LepEta")] -> Fill(wbana_lep_p4.Eta(), evtwt) ;
                if(goodAK4Jets.size() > 0){
                    h1_[name+string("_JetPt")]  -> Fill(goodAK4Jets.at(0).getP4().Pt(), evtwt) ;
                    h1_[name+string("_JetEta")] -> Fill(goodAK4Jets.at(0).getP4().Eta(), evtwt) ;
                    h1_[name+string("_DR")]     -> Fill(dr_ak4jet_lep, evtwt) ;
                    h1_[name+string("_DPHI")]   -> Fill(dphi, evtwt) ;
                    h1_[name+string("_FwdJetPt")]   -> Fill(fwd_ak4jet_p4.Pt(), evtwt) ;
                    h1_[name+string("_FwdJetEta")]   -> Fill(fwd_ak4jet_p4.Eta(), evtwt) ;
                }
                
                h1_[name+string("_NJets")]    -> Fill(goodAK4Jets.size(), evtwt) ;
                h1_[name+string("_ST")]       -> Fill(st, evtwt) ;
                h1_[name+string("_STv2")]    -> Fill(st_v2, evtwt) ;
                h1_[name+string("_Mass")]     -> Fill(mass, evtwt) ;
                h1_[name+string("_Massv2")]  -> Fill(mass_v2, evtwt);
                h1_[name+string("_HT")]       -> Fill(htak4.getHT(), evtwt);
                h1_[name+string("_NBTags")]   -> Fill(goodBTaggedAK4Jets.size(), evtwt);
                h1_[name+string("_ST_nomet")]       -> Fill(st_nomet, evtwt) ;
                h1_[name+string("_DPHI_MetJet")]       -> Fill(dphi_metjet, evtwt) ;
                h1_[name+string("_DPHI_MetLep")]       -> Fill(dphi_metlep, evtwt) ;
                h1_[name+string("_Angle_MuZ_Jet")]       -> Fill(angle_MuZ_Jet, evtwt) ;
                h1_[name+string("_Angle_MuJet_Met")]       -> Fill(angle_MuJet_Met, evtwt) ;
            }
        }
        
        
        // cutflow for the w+jet cr
        if(goodBTaggedAK4Jets.size() == 0 && htak4.getHT() > 400){
            for(int cut = 0; cut < cut_flow_wjet; cut++){
                
                string name = string("wjet_")+to_string(cut) + ch;
                
                h1_[name+string("_MET")]    -> Fill(goodMet.at(0).getP4().Pt(), evtwt) ;
                h1_[name+string("_MT")]     -> Fill(mt, evtwt) ;
                h1_[name+string("_LepPt")]  -> Fill(wbana_lep_p4.Pt(), evtwt) ;
                h1_[name+string("_LepEta")] -> Fill(wbana_lep_p4.Eta(), evtwt) ;
                
                if(goodAK4Jets.size() > 0){
                    h1_[name+string("_JetPt")]  -> Fill(goodAK4Jets.at(0).getP4().Pt(), evtwt) ;
                    h1_[name+string("_JetEta")] -> Fill(goodAK4Jets.at(0).getP4().Eta(), evtwt) ;
                    h1_[name+string("_DR")]     -> Fill(dr_ak4jet_lep, evtwt) ;
                    h1_[name+string("_DPHI")]   -> Fill(dphi, evtwt) ;
                    h1_[name+string("_FwdJetPt")]   -> Fill(fwd_ak4jet_p4.Pt(), evtwt) ;
                    h1_[name+string("_FwdJetEta")]   -> Fill(fwd_ak4jet_p4.Eta(), evtwt) ;
                    for (unsigned int i = 0; i < goodAK4Jets.size(); i++){
                      if (fabs(goodAK4Jets.at(i).getP4().Eta()) > 2.4){
                        h1_[name+string("_allFwdJetPt")]  -> Fill(goodAK4Jets.at(i).getP4().Pt(), evtwt) ;
                        h1_[name+string("_allFwdJetEta")]  -> Fill(goodAK4Jets.at(i).getP4().Eta(), evtwt) ;
                        h2_[name+string("_allFwdJetPtEta")]  -> Fill(goodAK4Jets.at(i).getP4().Pt(), goodAK4Jets.at(i).getP4().Eta(), evtwt) ;
                      }
                    }
                }
                
                h1_[name+string("_NJets")]    -> Fill(goodAK4Jets.size(), evtwt) ;
                h1_[name+string("_ST")]       -> Fill(st, evtwt) ;
                h1_[name+string("_STv2")]    -> Fill(st_v2, evtwt) ;
                h1_[name+string("_Mass")]     -> Fill(mass, evtwt) ;
                h1_[name+string("_Massv2")]  -> Fill(mass_v2, evtwt);
                h1_[name+string("_HT")]       -> Fill(htak4.getHT(), evtwt);
                h1_[name+string("_NBTags")]   -> Fill(goodBTaggedAK4Jets.size(), evtwt);
                h1_[name+string("_ST_nomet")]       -> Fill(st_nomet, evtwt) ;
                h1_[name+string("_DPHI_MetJet")]       -> Fill(dphi_metjet, evtwt) ;
                h1_[name+string("_DPHI_MetLep")]       -> Fill(dphi_metlep, evtwt) ;
                h1_[name+string("_Angle_MuZ_Jet")]       -> Fill(angle_MuZ_Jet, evtwt) ;
                h1_[name+string("_Angle_MuJet_Met")]       -> Fill(angle_MuJet_Met, evtwt) ;
            }
        }

        // cutflow for the signal region
        
        if(htak4.getHT() > 300){
            for(int cut = 0; cut < cut_flow_signal; cut++){
                
                string name = string("signal_")+to_string(cut) + ch;
                
                h1_[name+string("_MET")]    -> Fill(goodMet.at(0).getP4().Pt(), evtwt) ;
                h1_[name+string("_MT")]     -> Fill(mt, evtwt) ;
                h1_[name+string("_LepPt")]  -> Fill(wbana_lep_p4.Pt(), evtwt) ;
                h1_[name+string("_LepEta")] -> Fill(wbana_lep_p4.Eta(), evtwt) ;
                
                if(goodAK4Jets.size() > 0){
		    h1_[name+string("_JetCSV")]  ->Fill(goodAK4Jets.at(0).getCSV(), evtwt) ;
                    h1_[name+string("_JetPt")]  -> Fill(goodAK4Jets.at(0).getP4().Pt(), evtwt) ;
                    h1_[name+string("_JetEta")] -> Fill(goodAK4Jets.at(0).getP4().Eta(), evtwt) ;
                    h1_[name+string("_DR")]     -> Fill(dr_ak4jet_lep, evtwt) ;
                    h1_[name+string("_DPHI")]   -> Fill(dphi, evtwt) ;
                    h1_[name+string("_FwdJetPt")]   -> Fill(fwd_ak4jet_p4.Pt(), evtwt) ;
                    h1_[name+string("_FwdJetEta")]   -> Fill(fwd_ak4jet_p4.Eta(), evtwt) ;
                }
                
                h1_[name+string("_NJets")]    -> Fill(goodAK4Jets.size(), evtwt) ;
                h1_[name+string("_ST")]       -> Fill(st, evtwt) ;
                h1_[name+string("_STv2")]    -> Fill(st_v2, evtwt) ;
                h1_[name+string("_Mass")]     -> Fill(mass, evtwt) ;
                h1_[name+string("_Massv2")]  -> Fill(mass_v2, evtwt);
                h1_[name+string("_HT")]       -> Fill(htak4.getHT(), evtwt);
                h1_[name+string("_NBTags")]   -> Fill(goodBTaggedAK4Jets.size(), evtwt);
                h1_[name+string("_ST_nomet")]       -> Fill(st_nomet, evtwt) ;
                h1_[name+string("_DPHI_MetJet")]       -> Fill(dphi_metjet, evtwt) ;
                h1_[name+string("_DPHI_MetLep")]       -> Fill(dphi_metlep, evtwt) ;
                h1_[name+string("_Angle_MuZ_Jet")]       -> Fill(angle_MuZ_Jet, evtwt) ;
                h1_[name+string("_Angle_MuJet_Met")]       -> Fill(angle_MuJet_Met, evtwt) ;
                
                if(goodBTaggedAK4Jets.size() == 1){
                    h1_[name+string("_Mass1b")]     -> Fill(mass, evtwt) ;
                    h1_[name+string("_ST1b")]       -> Fill(st, evtwt) ;
                    if (htak4.getHT() > 300) h1_[name+string("_Mass1b")+string("_HT300")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 350) h1_[name+string("_Mass1b")+string("_HT350")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 400) h1_[name+string("_Mass1b")+string("_HT400")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 450) h1_[name+string("_Mass1b")+string("_HT450")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 500) h1_[name+string("_Mass1b")+string("_HT500")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 550) h1_[name+string("_Mass1b")+string("_HT550")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 600) h1_[name+string("_Mass1b")+string("_HT600")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 650) h1_[name+string("_Mass1b")+string("_HT650")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 700) h1_[name+string("_Mass1b")+string("_HT700")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 750) h1_[name+string("_Mass1b")+string("_HT750")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 800) h1_[name+string("_Mass1b")+string("_HT800")]     -> Fill(mass, evtwt) ;

                    if (st > 200) h1_[name+string("_Mass1b")+string("_ST200")]     -> Fill(mass, evtwt) ;
                    if (st > 250) h1_[name+string("_Mass1b")+string("_ST250")]     -> Fill(mass, evtwt) ;
                    if (st > 300) h1_[name+string("_Mass1b")+string("_ST300")]     -> Fill(mass, evtwt) ;
                    if (st > 350) h1_[name+string("_Mass1b")+string("_ST350")]     -> Fill(mass, evtwt) ;
                    if (st > 400) h1_[name+string("_Mass1b")+string("_ST400")]     -> Fill(mass, evtwt) ;
                    if (st > 450) h1_[name+string("_Mass1b")+string("_ST450")]     -> Fill(mass, evtwt) ;
                    if (st > 500) h1_[name+string("_Mass1b")+string("_ST500")]     -> Fill(mass, evtwt) ;
                    if (st > 550) h1_[name+string("_Mass1b")+string("_ST550")]     -> Fill(mass, evtwt) ;
                    if (st > 600) h1_[name+string("_Mass1b")+string("_ST600")]     -> Fill(mass, evtwt) ;
                    if (st > 650) h1_[name+string("_Mass1b")+string("_ST650")]     -> Fill(mass, evtwt) ;
                    if (st > 700) h1_[name+string("_Mass1b")+string("_ST700")]     -> Fill(mass, evtwt) ;
                    if (st > 750) h1_[name+string("_Mass1b")+string("_ST750")]     -> Fill(mass, evtwt) ;
                    if (st > 800) h1_[name+string("_Mass1b")+string("_ST800")]     -> Fill(mass, evtwt) ;
                    if (st > 850) h1_[name+string("_Mass1b")+string("_ST850")]     -> Fill(mass, evtwt) ;
                    if (st > 900) h1_[name+string("_Mass1b")+string("_ST900")]     -> Fill(mass, evtwt) ;
                    if (st > 950) h1_[name+string("_Mass1b")+string("_ST950")]     -> Fill(mass, evtwt) ;
                }
                if(goodBTaggedAK4Jets.size() > 1){
                    h1_[name+string("_Mass2b")]     -> Fill(mass, evtwt) ;
                    h1_[name+string("_ST2b")]       -> Fill(st, evtwt) ;
                    if (htak4.getHT() > 300) h1_[name+string("_Mass2b")+string("_HT300")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 350) h1_[name+string("_Mass2b")+string("_HT350")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 400) h1_[name+string("_Mass2b")+string("_HT400")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 450) h1_[name+string("_Mass2b")+string("_HT450")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 500) h1_[name+string("_Mass2b")+string("_HT500")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 550) h1_[name+string("_Mass2b")+string("_HT550")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 600) h1_[name+string("_Mass2b")+string("_HT600")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 650) h1_[name+string("_Mass2b")+string("_HT650")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 700) h1_[name+string("_Mass2b")+string("_HT700")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 750) h1_[name+string("_Mass2b")+string("_HT750")]     -> Fill(mass, evtwt) ;
                    if (htak4.getHT() > 800) h1_[name+string("_Mass2b")+string("_HT800")]     -> Fill(mass, evtwt) ;

                    if (st > 200) h1_[name+string("_Mass2b")+string("_ST200")]     -> Fill(mass, evtwt) ;
                    if (st > 250) h1_[name+string("_Mass2b")+string("_ST250")]     -> Fill(mass, evtwt) ;
                    if (st > 300) h1_[name+string("_Mass2b")+string("_ST300")]     -> Fill(mass, evtwt) ;
                    if (st > 350) h1_[name+string("_Mass2b")+string("_ST350")]     -> Fill(mass, evtwt) ;
                    if (st > 400) h1_[name+string("_Mass2b")+string("_ST400")]     -> Fill(mass, evtwt) ;
                    if (st > 450) h1_[name+string("_Mass2b")+string("_ST450")]     -> Fill(mass, evtwt) ;
                    if (st > 500) h1_[name+string("_Mass2b")+string("_ST500")]     -> Fill(mass, evtwt) ;
                    if (st > 550) h1_[name+string("_Mass2b")+string("_ST550")]     -> Fill(mass, evtwt) ;
                    if (st > 600) h1_[name+string("_Mass2b")+string("_ST600")]     -> Fill(mass, evtwt) ;
                    if (st > 650) h1_[name+string("_Mass2b")+string("_ST650")]     -> Fill(mass, evtwt) ;
                    if (st > 700) h1_[name+string("_Mass2b")+string("_ST700")]     -> Fill(mass, evtwt) ;
                    if (st > 750) h1_[name+string("_Mass2b")+string("_ST750")]     -> Fill(mass, evtwt) ;
                    if (st > 800) h1_[name+string("_Mass2b")+string("_ST800")]     -> Fill(mass, evtwt) ;
                    if (st > 850) h1_[name+string("_Mass2b")+string("_ST850")]     -> Fill(mass, evtwt) ;
                    if (st > 900) h1_[name+string("_Mass2b")+string("_ST900")]     -> Fill(mass, evtwt) ;
                    if (st > 950) h1_[name+string("_Mass2b")+string("_ST950")]     -> Fill(mass, evtwt) ;
                }
            }
        }
    }

    
    return true ;
}

// ------------ method called once each job just before starting event loop  ------------
void SingleLepAna::beginJob() {
    
    h1_["truePU"] = fs->make<TH1D>("truePU", "truePU", 51, -0.5, 50.5);
    h1_["nPU"] = fs->make<TH1D>("nPU", "nPU", 51, -0.5, 50.5);
    
    // b-tagging efficiency maps:
    h2_["pt_eta_b_all"] = fs->make<TH2D>("pt_eta_b_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_c_all"] = fs->make<TH2D>("pt_eta_c_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_l_all"] = fs->make<TH2D>("pt_eta_l_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 

    h2_["pt_eta_b_btagged"] = fs->make<TH2D>("pt_eta_b_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_c_btagged"] = fs->make<TH2D>("pt_eta_c_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_l_btagged"] = fs->make<TH2D>("pt_eta_l_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ;

    string full_name;
    
    for(int ch_i = 0; ch_i < 2; ch_i++){
        for(int cut = 0; cut < 11; cut++){
            
            
            string ch;
            if(ch_i == 0)
                ch = "_Mu";
            else
                ch = "_Ele";
            
            string name = string("signal_")+to_string(cut) + ch;
            
            full_name = name+string("_MET");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 1000);
            
            full_name = name+string("_MT");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 500);
            
            full_name = name+string("_LepPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 1000);
            
            full_name = name+string("_LepEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 120, -3, 3);
            
            full_name = name+string("_JetPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 200, 0, 2000);

 	    full_name = name+string("_JetCSV");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 150, -1.5, 1.5);
            
            full_name = name+string("_JetEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_FwdJetPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 500);
            
            full_name = name+string("_FwdJetEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);

            full_name = name+string("_DR");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
            full_name = name+string("_DPHI");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_NJets");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 20, -0.5, 19.5);
            
            full_name = name+string("_ST");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_STv2");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Massv2");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_ST1b");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_ST2b");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_HT");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_NBTags");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 5, -0.5, 4.5);
            
            full_name = name+string("_ST_nomet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_DPHI_MetJet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_DPHI_MetLep");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);

            full_name = name+string("_Angle_MuZ_Jet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
            full_name = name+string("_Angle_MuJet_Met");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
            full_name = name+string("_Mass1b")+string("_HT300");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT300");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_HT350");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT350");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_HT400");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT400");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_HT450");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT450");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_HT500");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT500");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_HT550");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT550");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_HT600");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT600");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_HT650");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT650");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_HT700");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT700");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_HT750");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT750");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_HT800");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_HT800");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            


            full_name = name+string("_Mass1b")+string("_ST200");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST200");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);

            full_name = name+string("_Mass1b")+string("_ST250");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST250");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);

            full_name = name+string("_Mass1b")+string("_ST300");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST300");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST350");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST350");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST400");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST400");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST450");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST450");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST500");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST500");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST550");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST550");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST600");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST600");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST650");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST650");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST700");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST700");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST750");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST750");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST800");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST800");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass1b")+string("_ST850");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST850");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);

            full_name = name+string("_Mass1b")+string("_ST900");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST900");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);

            full_name = name+string("_Mass1b")+string("_ST950");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass2b")+string("_ST950");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);

        }
    }
    
    
    for(int ch_i = 0; ch_i < 2; ch_i++){
        for(int cut = 0; cut < 11; cut++){
            
            
            string ch;
            if(ch_i == 0)
                ch = "_Mu";
            else
                ch = "_Ele";
            
            string name = string("tt_")+to_string(cut) + ch;
            
            full_name = name+string("_MET");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 1000);
            
            full_name = name+string("_MT");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 500);
            
            full_name = name+string("_LepPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 1000);
            
            full_name = name+string("_LepEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 120, -3, 3);
            
            full_name = name+string("_JetPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 200, 0, 2000);
            
            full_name = name+string("_JetEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_FwdJetPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 500);
            
            full_name = name+string("_FwdJetEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);

            
            full_name = name+string("_DR");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
            full_name = name+string("_DPHI");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_NJets");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 20, -0.5, 19.5);
            
            full_name = name+string("_ST");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_STv2");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Massv2");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_HT");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_NBTags");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 5, -0.5, 4.5);

            full_name = name+string("_ST_nomet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_DPHI_MetJet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_DPHI_MetLep");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);

            full_name = name+string("_Angle_MuZ_Jet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
            full_name = name+string("_Angle_MuJet_Met");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
        }
    }
    
    
    for(int ch_i = 0; ch_i < 2; ch_i++){
        for(int cut = 0; cut < 11; cut++){
            
            
            string ch;
            if(ch_i == 0)
                ch = "_Mu";
            else
                ch = "_Ele";
            
            string name = string("cr_")+to_string(cut) + ch;
            
            full_name = name+string("_MET");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 1000);
            
            full_name = name+string("_MT");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 500);
            
            full_name = name+string("_LepPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 1000);
            
            full_name = name+string("_LepEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 120, -3, 3);
            
            full_name = name+string("_JetPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 200, 0, 2000);
            
            full_name = name+string("_JetEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_FwdJetPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 500);
            
            full_name = name+string("_FwdJetEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);

            
            full_name = name+string("_DR");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
            full_name = name+string("_DPHI");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_NJets");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 20, -0.5, 19.5);
            
            full_name = name+string("_ST");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_STv2");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Massv2");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_HT");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_NBTags");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 5, -0.5, 4.5);

            full_name = name+string("_ST_nomet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_DPHI_MetJet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_DPHI_MetLep");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);

            full_name = name+string("_Angle_MuZ_Jet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
            full_name = name+string("_Angle_MuJet_Met");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
        }
    }
    
    
    
    for(int ch_i = 0; ch_i < 2; ch_i++){
        for(int cut = 0; cut < 11; cut++){
            
            
            string ch;
            if(ch_i == 0)
                ch = "_Mu";
            else
                ch = "_Ele";
            
            string name = string("wjet_")+to_string(cut) + ch;
            
            full_name = name+string("_MET");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 1000);
            
            full_name = name+string("_MT");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 500);
            
            full_name = name+string("_LepPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 1000);
            
            full_name = name+string("_LepEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 120, -3, 3);
            
            full_name = name+string("_JetPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 200, 0, 2000);
            
            full_name = name+string("_JetEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_allFwdJetPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 200, 0, 2000);
            
            full_name = name+string("_allFwdJetEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_allFwdJetPtEta");
            h2_[full_name]  =  fs->make<TH2D>(full_name.c_str(), full_name.c_str(), 200, 0, 2000, 100, -5, 5);
            
            full_name = name+string("_FwdJetPt");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 500);
            
            full_name = name+string("_FwdJetEta");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);

            
            full_name = name+string("_DR");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
            full_name = name+string("_DPHI");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_NJets");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 20, -0.5, 19.5);
            
            full_name = name+string("_ST");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_STv2");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Mass");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_Massv2");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_HT");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_NBTags");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 5, -0.5, 4.5);

            full_name = name+string("_ST_nomet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 300, 0, 3000);
            
            full_name = name+string("_DPHI_MetJet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);
            
            full_name = name+string("_DPHI_MetLep");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, -5, 5);

            full_name = name+string("_Angle_MuZ_Jet");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
            full_name = name+string("_Angle_MuJet_Met");
            h1_[full_name]  =  fs->make<TH1D>(full_name.c_str(), full_name.c_str(), 100, 0, 5);
            
        }
    }
    //LQ Histograms
    h1_["mass_LQ"]    = fs->make<TH1D>("mass_LQ", "massLQ", 3000, 0, 3000);
    h1_["mass_noveto_LQ"]  = fs->make<TH1D>("mass_noveto_LQ", "massLQ", 3000, 0, 3000);
    h1_["mass_skim_LQ"]  = fs->make<TH1D>("mass_skim_LQ", "massLQ", 3000, 0, 3000);
    h1_["leppt_LQ"]   = fs->make<TH1D>("leppt_LQ", "lepptLQ", 2000, 0, 2000);
    h1_["lepeta_LQ"]  = fs->make<TH1D>("lepeta_LQ", "lepetaLQ", 100, -5, 5);
    h1_["lepphi_LQ"]  = fs->make<TH1D>("lepphi_LQ", "lepphiLQ", 100, -5, 5);
    h1_["qpt_LQ"]     = fs->make<TH1D>("qpt_LQ", "qptLQ", 2000, 0, 2000);
    h1_["qeta_LQ"]    = fs->make<TH1D>("qeta_LQ", "qetaLQ", 100, -5, 5);
    h1_["qphi_LQ"]    = fs->make<TH1D>("qphi_LQ", "qphiLQ", 100, -5, 5);
    h1_["mt_LQ"]      = fs->make<TH1D>("mt_LQ", "mtLQ", 2000, 0, 2000);
    h1_["met_LQ"]     = fs->make<TH1D>("met_LQ", "metLQ", 1000, 0, 1000);
}




bool SingleLepAna::solve_nu(const TLorentzVector &vlep, const TLorentzVector &vnu, double wmass, double& nuz1, double& nuz2){
    //
    // Purpose: Solve for the neutrino longitudinal z-momentum that makes
    // the leptonic W have mass WMASS.
    //
    // Inputs:
    // ev - The event to solve.
    // wmass - The desired W mass.
    //
    // Outputs:
    // nuz1 - First solution (smaller absolute value).
    // nuz2 - Second solution.
    //
    // Returns:
    // True if there was a real solution. False if there were only
    // imaginary solutions. (In that case, we just set the imaginary
    // part to zero.)
    //
    TLorentzVector tmp;
    tmp = vlep;
    bool discrim_flag = true;
    
    //std::cout << "Before: " << vlep.Pt() << " " << vlep.Eta() << " " << vlep.Phi() << " " << vlep.M() << std::endl;
    
    
    double x = vlep.X()*vnu.X() + vlep.Y()*vnu.Y() + wmass*wmass/2;
    double a = vlep.Z()*vlep.Z() - vlep.E()*vlep.E();
    double b = 2*x*vlep.Z();
    double c = x*x - vnu.Perp2() * vlep.E()*vlep.E();
    double d = b*b - 4*a*c;
    if (d < 0){
        d = 0;
        discrim_flag = false;
    }
    
    nuz1 = (-b + sqrt(d))/2/a;
    nuz2 = (-b - sqrt(d))/2/a;
    if (abs(nuz1) > abs(nuz2))
        swap (nuz1, nuz2);
    
    //std::cout << "After: " << vlep.Pt() << " " << vlep.Eta() << " " << vlep.Phi() << " " << vlep.M() << std::endl;
    if (tmp != vlep){
        std::cout << "Before: " << tmp.Pt() << " " << tmp.Eta() << " " << tmp.Phi() << " " << tmp.M() << std::endl;
        std::cout << "After: " << vlep.Pt() << " " << vlep.Eta() << " " << vlep.Phi() << " " << vlep.M() << std::endl;
    }
    return discrim_flag;
}

void SingleLepAna::adjust_e_for_mass(TLorentzVector& v, double mass){
    //
    // Purpose: Adjust the energy component of V (leaving the 3-vector part
    // unchanged) so that it has mass MASS.
    //
    // Inputs:
    // v - The 4-vector to scale.
    // mass - The desired mass of the 4-vector.
    //
    // Outputs:
    // v - The scaled 4-vector.
    //
    v.SetE(sqrt(v.Vect().Mag2() + mass*mass));
}

int SingleLepAna::cutFlow(const edm::ParameterSet & pars){
    
    double lep_pt       = pars.getParameter<double>("lep_pt");
    double lep_eta      = pars.getParameter<double>("lep_eta");
    double met_et       = pars.getParameter<double>("met_et");
    double jet_pt       = pars.getParameter<double>("jet_pt");
    double jet_eta      = pars.getParameter<double>("jet_eta");
    double min_bdisc    = pars.getParameter<double>("min_bdisc");
    double max_bdisc    = pars.getParameter<double>("max_bdisc");
    double min_mt       = pars.getParameter<double>("min_mt");
    double max_mt       = pars.getParameter<double>("max_mt");
    double min_dr       = pars.getParameter<double>("min_dr");
    double max_dr       = pars.getParameter<double>("max_dr");
    double min_dphi     = pars.getParameter<double>("min_dphi");
    double max_dphi     = pars.getParameter<double>("max_dphi");
    double fwd_jet_eta  = pars.getParameter<double>("fwd_jet_eta");
    double min_st       = pars.getParameter<double>("min_st");
    
    
    TLorentzVector lep_p4;
    if(signalElectrons.size() == 1){
        lep_p4 = signalElectrons.at(0).getP4();
    }
    else if(signalMuons.size() == 1){
        lep_p4 = signalMuons.at(0).getP4();
    }
    else{
        
        cout<<"lepton selection bug" <<endl;
    }
    
    TLorentzVector met_p4 = goodMet.at(0).getP4();
    double mt = TMath::Sqrt( 2*lep_p4.Pt() * met_p4.Pt() * ( 1 - TMath::Cos(lep_p4.DeltaPhi(met_p4) ) ) );
    
    TLorentzVector leading_ak4jet_p4;
    double leading_ak4jet_bdisc = -99.,  dphi = -99.;
    
    if(goodAK4Jets.size() > 0){
        leading_ak4jet_p4    = goodAK4Jets.at(0).getP4();
        leading_ak4jet_bdisc = goodAK4Jets.at(0).getCSV();
        dphi = leading_ak4jet_p4.DeltaPhi(lep_p4);
    }
    
    double dr_ak4jet_lep = 999, fwd_ak4jet_eta = 0;
    TLorentzVector fwd_ak4jet_p4;
    
    double st = met_p4.Pt() + leading_ak4jet_p4.Pt() + lep_p4.Pt();
    
    for (unsigned int jet = 0; jet < goodAK4Jets.size(); jet++) {
        TLorentzVector jet_p4 = goodAK4Jets.at(jet).getP4();
        
        // calculate min dr between lepton and the closest jet
        double dr = jet_p4.DeltaR(lep_p4);
        if(dr < dr_ak4jet_lep && jet_p4.Pt() > 40){
            dr_ak4jet_lep = dr;
        }
        
        // find the forward jet
        if( fabs( jet_p4.Eta() ) > fwd_ak4jet_eta){
            fwd_ak4jet_eta = fabs( jet_p4.Eta() );
            fwd_ak4jet_p4  = jet_p4;
        }
    }
    
    
    
    if(lep_p4.Pt() < lep_pt         &&
       fabs( lep_p4.Eta() ) > lep_eta)            return 1;
    if(met_p4.Pt() < met_et)                      return 2;
    if(!(mt > min_mt && mt < max_mt) )            return 3;
    if(leading_ak4jet_p4.Pt() < jet_pt)           return 4;
    if(fabs( leading_ak4jet_p4.Eta() ) > jet_eta) return 5;
    if(!(leading_ak4jet_bdisc >= min_bdisc &&
         leading_ak4jet_bdisc <= max_bdisc ) )     return 6;
    if(!(fabs( dphi ) >= min_dphi &&
         fabs( dphi ) <= max_dphi) )              return 7;
    if(!(fabs( dr_ak4jet_lep ) >= min_dr &&
         fabs( dr_ak4jet_lep ) <= max_dr))        return 8;
    if( st < min_st)                              return 9;
    if(fabs(fwd_ak4jet_p4.Eta()) < fwd_jet_eta) return 10;
    
    return 11;
}

double SingleLepAna::getWjetHTRewSF(double ht, double ht_low_threshold){

    double scale_factor = 1.09189e+00 - ht * 2.43398e-04;
    return ht > ht_low_threshold ? scale_factor : 1.0;
}

void SingleLepAna::endJob() {
    
    return ;
}

DEFINE_FWK_MODULE(SingleLepAna);


  

