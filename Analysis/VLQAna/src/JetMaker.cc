#include "Analysis/VLQAna/interface/JetMaker.h"
#include "Analysis/VLQAna/interface/ApplyJER.h"
#include "Analysis/VLQAna/interface/Utilities.h"

#include <TRandom.h>

#define DEBUGMORE false
#define DEBUG false

using namespace std;
using namespace edm ;

JetMaker::JetMaker () {}

JetMaker::~JetMaker () {}

JetMaker::JetMaker (edm::ParameterSet const& iConfig, edm::ConsumesCollector && iC) :
JetIDParams_              (iConfig.getParameter<edm::ParameterSet>       ("JetIDParams")),
JetSubstrParams_          (iConfig.getParameter<edm::ParameterSet>       ("JetSubstrParams")),
SubjetParams_             (iConfig.getParameter<edm::ParameterSet>       ("SubjetParams")),
scaleJetP4_               (iConfig.getParameter<bool>                    ("scaleJetP4")),
scaledJetMass_            (iConfig.getParameter<double>                  ("scaledJetMass")),
jecShift_                 (iConfig.getParameter<double>                  ("jecShift")),
jerShift_                 (iConfig.getParameter<int>                     ("jerShift")),
newJECPayloadNames_       (iConfig.getParameter<std::vector<std::string>>("newJECPayloadNames")),
jecUncPayloadName_        (iConfig.getParameter<std::string>             ("jecUncPayloadName")),
jecAK8GroomedPayloadNames_(iConfig.getParameter<std::vector<std::string>>("jecAK8GroomedPayloadNames")),
doGroomedMassCorr_        (jecAK8GroomedPayloadNames_.size()>0 ? true : false),
jetID_(JetIDParams_, iC)
{
    
    std::string jettypeStr = iConfig.getParameter<std::string>("jettype") ;
    if      (jettypeStr == "AK4JET")             type_ = AK4JET;
    else if (jettypeStr == "AK8JET")             type_ = AK8JET;
    else edm::LogError("JetMaker::JetMaker") << " >>>> WrongJetType: " << type_<< " Check jet type !!!" ;
    
    t_npv               = iC.consumes<int>          (JetIDParams_.getParameter<edm::InputTag> ("npvLabel")              );
    t_rho               = iC.consumes<double>       (JetIDParams_.getParameter<edm::InputTag> ("rhoLabel")              );
    t_jetPt             = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetPtLabel")            );
    t_jetEta            = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetEtaLabel")           );
    t_jetPhi            = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetPhiLabel")           );
  //  t_jetMass           = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetMassLabel")          );
    t_jetEnergy         = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetEnergyLabel")        );
    t_jetHadronFlavour  = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetPartonFlavourLabel") );
    t_jetPartonFlavour  = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetHadronFlavourLabel") );
    t_jetCSV            = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetCSVLabel")           );
    t_jetCMVA           = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetCMVALabel")          );
    t_jetCvsB           = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetCvsBLabel")          );
    t_jetCvsL           = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetCvsLLabel")          );
 //   t_jetDoubleBAK8     = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetDoubleBAK8Label")    );
 //   t_jetDoubleBCA15    = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetDoubleBCA15Label")   );
    t_jetJEC            = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetJECLabel")           );
   // t_jetnumDaughters   = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetnumDaughtersLabel")  );
    t_jetY              = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetYLabel")             );
    t_jetArea           = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetAreaLabel")          );
    t_jetGenJetPt       = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetGenJetPtLabel")      );
    t_jetGenJetEta      = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetGenJetEtaLabel")     );
    t_jetGenJetPhi      = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetGenJetPhiLabel")     );
    t_jetGenJetE        = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetGenJetELabel")       );
    t_jetGenJetCharge   = iC.consumes<vector<float>>(JetIDParams_.getParameter<edm::InputTag> ("jetGenJetChargeLabel")  );
    
    if ( type_ == AK8JET ) {
        t_jettau1           = iC.consumes<vector<float>>(JetSubstrParams_.getParameter<edm::InputTag> ("jettau1Label")         );
        t_jettau2           = iC.consumes<vector<float>>(JetSubstrParams_.getParameter<edm::InputTag> ("jettau2Label")         );
        t_jettau3           = iC.consumes<vector<float>>(JetSubstrParams_.getParameter<edm::InputTag> ("jettau3Label")         );
        t_jetPrunedMass     = iC.consumes<vector<float>>(JetSubstrParams_.getParameter<edm::InputTag> ("jetPrunedMassLabel")   );
        t_jetTrimmedMass    = iC.consumes<vector<float>>(JetSubstrParams_.getParameter<edm::InputTag> ("jetTrimmedMassLabel")  );
        t_jetFilteredMass   = iC.consumes<vector<float>>(JetSubstrParams_.getParameter<edm::InputTag> ("jetFilteredMassLabel") );
        t_jetSoftDropMass   = iC.consumes<vector<float>>(JetSubstrParams_.getParameter<edm::InputTag> ("jetSoftDropMassLabel") );
        t_vjetsjIdx0        = iC.consumes<vector<float>>(JetSubstrParams_.getParameter<edm::InputTag> ("vjetsjIdx0Label")      );
        t_vjetsjIdx1        = iC.consumes<vector<float>>(JetSubstrParams_.getParameter<edm::InputTag> ("vjetsjIdx1Label")      );
        t_ak8sjHadFl        = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetHadronFlavourLabel"));
        t_ak8sjPt           = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetPtLabel")           );
        t_ak8sjEta          = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetEtaLabel")          );
        t_ak8sjPhi          = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetEtaLabel")          );
     //   t_ak8sjMass         = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetMassLabel")         );
        t_ak8sjCSV          = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetCSVLabel")          );
        t_ak8sjCMVA         = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetCMVALabel")         );
        t_ak8sjGenJetPt     = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetGenJetPtLabel")     );
        t_ak8sjGenJetEta    = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetGenJetEtaLabel")    );
        t_ak8sjGenJetPhi    = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetGenJetPhiLabel")    );
        t_ak8sjGenJetE      = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetGenJetELabel")      );
        t_ak8sjGenJetCharge = iC.consumes<vector<float>>(SubjetParams_   .getParameter<edm::InputTag> ("jetGenJetChargeLabel") );
        
        
        idxjettau2Bytau1Min_  = iConfig.getParameter<double>("jettau2Bytau1Min")  ;
        idxjettau2Bytau1Max_  = iConfig.getParameter<double>("jettau2Bytau1Max")  ;
        idxjettau3Bytau2Min_  = iConfig.getParameter<double>("jettau3Bytau2Min")  ;
        idxjettau3Bytau2Max_  = iConfig.getParameter<double>("jettau3Bytau2Max")  ;
        idxjetMassMin_        = iConfig.getParameter<double>("jetMassMin")  ;
        idxjetMassMax_        = iConfig.getParameter<double>("jetMassMax")  ;
        idxjetPrunedMassMin_  = iConfig.getParameter<double>("jetPrunedMassMin")  ;
        idxjetPrunedMassMax_  = iConfig.getParameter<double>("jetPrunedMassMax")  ;
        idxjetTrimmedMassMin_ = iConfig.getParameter<double>("jetTrimmedMassMin")  ;
        idxjetTrimmedMassMax_ = iConfig.getParameter<double>("jetTrimmedMassMax")  ;
        idxjetSoftDropMassMin_= iConfig.getParameter<double>("jetSoftDropMassMin")  ;
        idxjetSoftDropMassMax_= iConfig.getParameter<double>("jetSoftDropMassMax")  ;
        idxsjMassMin_         = iConfig.getParameter<double>("subjetMassMin")  ;
        idxsjCSVMin_          = iConfig.getParameter<double>("subjetCSVMin")  ;
        idxsjCSVMax_          = iConfig.getParameter<double>("subjetCSVMax")  ;
        idxsjHighestCSVMin_   = iConfig.getParameter<double>("subjetHighestCSVMin")  ;
        idxsjHighestCSVMax_   = iConfig.getParameter<double>("subjetHighestCSVMax")  ;
        btaggedcsvlOP_        = iConfig.getParameter<double>("btaggedcsvlOP")  ;
        btaggedcsvmOP_        = iConfig.getParameter<double>("btaggedcsvmOP")  ;
        btaggedcsvtOP_        = iConfig.getParameter<double>("btaggedcsvtOP")  ;
    }
    
    idxjetPtMin_ = iConfig.getParameter<double>("jetPtMin")  ;
    idxjetPtMax_ = iConfig.getParameter<double>("jetPtMax")  ;
    idxjetAbsEtaMax_ = iConfig.getParameter<double>("jetAbsEtaMax")  ;
    idxjetCSVDiscMin_ = iConfig.getParameter<double>("jetCSVDiscMin")  ;
    idxjetCSVDiscMax_ = iConfig.getParameter<double>("jetCSVDiscMax")  ;
    
    if (doGroomedMassCorr_) {
        jecAK8UncName_ = jecAK8GroomedPayloadNames_.back();
        std::vector<JetCorrectorParameters> vPar;
        for ( std::vector<std::string>::const_iterator it = jecAK8GroomedPayloadNames_.begin(); it != jecAK8GroomedPayloadNames_.end(); ++it) {
            JetCorrectorParameters pars(*it);
            vPar.push_back(pars) ;
        }
        ptr_jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
    }
    
    if (jecShift_ != 0 && !jecUncPayloadName_.empty() ) {
        ptr_jecUnc = boost::shared_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty(jecUncPayloadName_) ) ;
    }
    
    if ( newJECPayloadNames_.size() > 0 ) {
        std::vector<JetCorrectorParameters> vPar;
        for ( std::vector<std::string>::const_iterator it = newJECPayloadNames_.begin(); it != newJECPayloadNames_.end(); ++it) {
            JetCorrectorParameters pars(*it);
            vPar.push_back(pars) ;
        }
        ptr_newJEC_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
    }
    
}

void JetMaker::operator()(edm::Event& evt, vlq::JetCollection& jets) {
    
    Handle <int>            h_npv               ; evt.getByToken (t_npv                 , h_npv               );
    Handle <double>         h_rho               ; evt.getByToken (t_rho                 , h_rho               );
    Handle <vector<float>>  h_jetPt             ; evt.getByToken (t_jetPt               , h_jetPt             );
    Handle <vector<float>>  h_jetEta            ; evt.getByToken (t_jetEta              , h_jetEta            );
    Handle <vector<float>>  h_jetPhi            ; evt.getByToken (t_jetPhi              , h_jetPhi            );
  //  Handle <vector<float>>  h_jetMass           ; evt.getByToken (t_jetMass             , h_jetMass           );
    Handle <vector<float>>  h_jetEnergy         ; evt.getByToken (t_jetEnergy           , h_jetEnergy         );
    Handle <vector<float>>  h_jetHadronFlavour  ; evt.getByToken (t_jetHadronFlavour    , h_jetHadronFlavour  );
    Handle <vector<float>>  h_jetPartonFlavour  ; evt.getByToken (t_jetPartonFlavour    , h_jetPartonFlavour  );
    Handle <vector<float>>  h_jetCSV            ; evt.getByToken (t_jetCSV              , h_jetCSV            );
    Handle <vector<float>>  h_jetCMVA           ; evt.getByToken (t_jetCMVA             , h_jetCMVA           );
    Handle <vector<float>>  h_jetCvsB           ; evt.getByToken (t_jetCvsB             , h_jetCvsB           );
    Handle <vector<float>>  h_jetCvsL           ; evt.getByToken (t_jetCvsL             , h_jetCvsL           );
 //   Handle <vector<float>>  h_jetDoubleBAK8     ; evt.getByToken (t_jetDoubleBAK8       , h_jetDoubleBAK8     );
 //   Handle <vector<float>>  h_jetDoubleBCA15    ; evt.getByToken (t_jetDoubleBCA15      , h_jetDoubleBCA15    );
    Handle <vector<float>>  h_jetJEC            ; evt.getByToken (t_jetJEC              , h_jetJEC            );
 //   Handle <vector<float>>  h_jetnumDaughters   ; evt.getByToken (t_jetnumDaughters     , h_jetnumDaughters )  ;
    Handle <vector<float>>  h_jetY              ; evt.getByToken (t_jetY                , h_jetY              );
    Handle <vector<float>>  h_jetArea           ; evt.getByToken (t_jetArea             , h_jetArea           );
    Handle <vector<float>>  h_jetGenJetE        ; evt.getByToken (t_jetGenJetE          , h_jetGenJetE        );
    Handle <vector<float>>  h_jetGenJetEta      ; evt.getByToken (t_jetGenJetEta        , h_jetGenJetEta      );
    Handle <vector<float>>  h_jetGenJetPt       ; evt.getByToken (t_jetGenJetPt         , h_jetGenJetPt       );
    Handle <vector<float>>  h_jetGenJetPhi      ; evt.getByToken (t_jetGenJetPhi        , h_jetGenJetPhi      );
    
    if ((h_jetPt.product())->size() < 1) return ;
    
    const int npv(*h_npv) ;
    const double rho(*h_rho) ;
    
    for (unsigned ijet = 0; ijet < (h_jetPt.product())->size(); ++ijet) {
        
        double jetPt = (h_jetPt.product())->at(ijet) ;
        if (jetPt == 0.) continue;
        double jetAbsEta = abs((h_jetEta.product())->at(ijet)) ;
        if (jetAbsEta  > idxjetAbsEtaMax_) continue ;
        
        TLorentzVector  jetP4, uncorrJetP4, newJetP4;
        
        jetP4.SetPtEtaPhiE( (h_jetPt.product())->at(ijet),
                           (h_jetEta.product())->at(ijet),
                           (h_jetPhi.product())->at(ijet),
                           (h_jetEnergy.product())->at(ijet) ) ;
        
        uncorrJetP4 = jetP4 * (h_jetJEC.product())->at(ijet) ;
        
        double newJEC(1.0);
        if ( newJECPayloadNames_.size() > 0 ) {
            ptr_newJEC_->setJetPt ( uncorrJetP4.Pt()     );
            ptr_newJEC_->setJetEta( uncorrJetP4.Eta()    );
            ptr_newJEC_->setJetE  ( uncorrJetP4.Energy() );
            ptr_newJEC_->setJetA  ( (h_jetArea.product())->at(ijet) );
            ptr_newJEC_->setRho   ( rho ) ;
            ptr_newJEC_->setNPV   ( npv );
            newJEC = ptr_newJEC_->getCorrection();
            newJetP4 = uncorrJetP4*newJEC ;
        }
        else newJetP4 = jetP4 ;
        
#if DEBUGMORE
        cout
        << " \nold jec = " << 1./(h_jetJEC.product())->at(ijet)
        << " \nnew jec = " << newJEC
        << " \njet pt uncorrected   = " << uncorrJetP4.Pt()
        << " \njet mass uncorrected = " << uncorrJetP4.Mag()
        << " \njet pt before        = " << jetP4.Pt()
        << " \njet mass before      = " << jetP4.Mag()
        << " \njet pt newJec        = " << newJetP4.Pt()
        << " \njet mass newJEC      = " << newJetP4.Mag()
        << endl ;
#endif
        
        double ptsmear(1.) ;
        if ( jerShift_ != 0 ) {
            double pt_gen = (h_jetGenJetPt.product())->at(ijet) ;
            double pt_reco   = uncorrJetP4.Pt() ;
            double eta_reco = uncorrJetP4.Eta() ;
            double jerscalep4 = ApplyJERp4(eta_reco, jerShift_) ;
            if (pt_gen > 0.){ 
		ptsmear = std::max( 0.0, pt_gen + jerscalep4*(pt_reco - pt_gen) )/pt_reco;
	    }
            else if (jerscalep4 > 1. && pt_reco > 0.) {
               TRandom* rand = new TRandom();
               ptsmear = rand->Gaus(pt_reco, sqrt(jerscalep4*jerscalep4 - 1)*0.2)/pt_reco ; //// Assuming 20% JER
               delete rand; 
       	    }
            else ptsmear = 1.;

            newJetP4 *= ptsmear ;
#if DEBUG
            cout
            << " \n pt reco             = " << pt_reco
            << " \n pt gen              = " << pt_gen
            << " \n eta_reco            = " << eta_reco
            << " \n jerscalep4          = " << jerscalep4
            << " \n jer smear           = " << ptsmear
            << " \njet pt ptsmear       = " << newJetP4.Pt()
            << " \njet mass ptsmear     = " << newJetP4.Mag()
            << endl ;
#endif
        }
        
        double unc(0);
        if (jecShift_ != 0 ) {
            ptr_jecUnc->setJetEta( uncorrJetP4.Eta()    );
            ptr_jecUnc->setJetPt ( uncorrJetP4.Pt()     );
            unc = ptr_jecUnc->getUncertainty(true);
            newJetP4 *= (1 + jecShift_*unc) ;
        }
        
#if DEBUGMORE
        cout
        << " \njet pt jecshift      = " << newJetP4.Pt()
        << " \njet mass jecshift    = " << newJetP4.Mag()
        << endl ;
#endif
        
        
        if (jetID_(evt, ijet) == false) continue;
        
        
        double jetCSVDisc = (h_jetCSV.product())->at(ijet);
        if (jetPt      < idxjetPtMin_       ||
            jetPt      >  idxjetPtMax_      ||
            jetCSVDisc < idxjetCSVDiscMin_  ||
            jetCSVDisc >  idxjetCSVDiscMax_
            ) continue ;
        
        
        TLorentzVector genJetP4;
        genJetP4.SetPtEtaPhiE((h_jetGenJetPt.product())->at(ijet),
                              (h_jetGenJetEta.product())->at(ijet),
                              (h_jetGenJetPhi.product())->at(ijet),
                              (h_jetGenJetE.product())->at(ijet));
        
        //// Jet to put in the jet collection
        vlq::Jet jet ;
        
        if ( type_ == AK8JET ) {
            
            double jettau1         (-1000);
            double jettau2         (-1000);
            double jettau3         (-1000);
            double jetMass         (-1000);
            double jetPrunedMass   (-1000);
            double jetTrimmedMass  (-1000);
            double jetSoftDropMass (-1000);
            int    vjetsjIdx0      (-1000);
            int    vjetsjIdx1      (-1000);
            
            double vjetssj0HadFl (-1000) ;
            double vjetssj0Pt    (-1000) ;
            double vjetssj0Eta   (-1000) ;
            double vjetssj0Phi   (-1000) ;
            double vjetssj0Mass  (-1000) ;
            double vjetssj0CSV   (-1000) ;
            double vjetssj0CMVA  (-1000) ;
            
            double vjetssj1HadFl (-1000) ;
            double vjetssj1Pt    (-1000) ;
            double vjetssj1Eta   (-1000) ;
            double vjetssj1Phi   (-1000) ;
            double vjetssj1Mass  (-1000) ;
            double vjetssj1CSV   (-1000) ;
            double vjetssj1CMVA  (-1000) ;
            
            int nsubjetsbtaggedcsvl(0) ;
            TLorentzVector p4sj0, p4sj1;
            TLorentzVector p4gensj0, p4gensj1;
            
            Handle <vector<float>>  h_jettau1           ; evt.getByToken (t_jettau1             , h_jettau1        );
            Handle <vector<float>>  h_jettau2           ; evt.getByToken (t_jettau2             , h_jettau2        );
            Handle <vector<float>>  h_jettau3           ; evt.getByToken (t_jettau3             , h_jettau3        );
            Handle <vector<float>>  h_jetPrunedMass     ; evt.getByToken (t_jetPrunedMass       , h_jetPrunedMass  );
            Handle <vector<float>>  h_jetTrimmedMass    ; evt.getByToken (t_jetTrimmedMass      , h_jetTrimmedMass );
            Handle <vector<float>>  h_jetFilteredMass   ; evt.getByToken (t_jetFilteredMass     , h_jetFilteredMass);
            Handle <vector<float>>  h_jetSoftDropMass   ; evt.getByToken (t_jetSoftDropMass     , h_jetSoftDropMass);
            Handle <vector<float>>  h_ak8sjHadFl        ; evt.getByToken (t_ak8sjHadFl          , h_ak8sjHadFl     );
            Handle <vector<float>>  h_ak8sjPt           ; evt.getByToken (t_ak8sjPt             , h_ak8sjPt        );
            Handle <vector<float>>  h_ak8sjEta          ; evt.getByToken (t_ak8sjEta            , h_ak8sjEta       );
            Handle <vector<float>>  h_ak8sjPhi          ; evt.getByToken (t_ak8sjPhi            , h_ak8sjPhi       );
         //   Handle <vector<float>>  h_ak8sjMass         ; evt.getByToken (t_ak8sjMass           , h_ak8sjMass      );
            Handle <vector<float>>  h_ak8sjCSV          ; evt.getByToken (t_ak8sjCSV            , h_ak8sjCSV       );
            Handle <vector<float>>  h_ak8sjCMVA         ; evt.getByToken (t_ak8sjCMVA           , h_ak8sjCMVA      );
            Handle <vector<float>>  h_ak8sjGenJetPt     ; evt.getByToken (t_ak8sjGenJetPt       , h_ak8sjGenJetPt  );
            Handle <vector<float>>  h_ak8sjGenJetEta    ; evt.getByToken (t_ak8sjGenJetEta      , h_ak8sjGenJetEta );
            Handle <vector<float>>  h_ak8sjGenJetPhi    ; evt.getByToken (t_ak8sjGenJetPhi      , h_ak8sjGenJetPhi );
            Handle <vector<float>>  h_ak8sjGenJetE      ; evt.getByToken (t_ak8sjGenJetE        , h_ak8sjGenJetE   );
            Handle <vector<float>>  h_ak8sjGenJetCharge ; evt.getByToken (t_ak8sjGenJetCharge   , h_ak8sjGenJetE   );
            Handle <vector<float>>  h_vjetsjIdx0        ; evt.getByToken (t_vjetsjIdx0          , h_vjetsjIdx0     );
            Handle <vector<float>>  h_vjetsjIdx1        ; evt.getByToken (t_vjetsjIdx1          , h_vjetsjIdx1     );
            
            double massCorr(1.);
            if ( doGroomedMassCorr_ ) {
                ptr_jecAK8_->setJetPt ( uncorrJetP4.Pt()     );
                ptr_jecAK8_->setJetEta( uncorrJetP4.Eta()    );
                ptr_jecAK8_->setJetE  ( uncorrJetP4.Energy() );
                ptr_jecAK8_->setJetA  ( (h_jetArea.product())->at(ijet) );
                ptr_jecAK8_->setRho   ( rho );
                ptr_jecAK8_->setNPV   ( npv );
                massCorr = ptr_jecAK8_->getCorrection();
            }
            
            double masssmear(1.) ;
            if ( jerShift_ != 0 ) {
                double pt_gen = (h_jetGenJetPt.product())->at(ijet) ;
                double pt_reco   = uncorrJetP4.Pt() ;
                double jerscalemass = ApplyJERMass() ;
                masssmear = std::max( 0.0, pt_gen + jerscalemass*(pt_reco - pt_gen) )/pt_reco ;
            }
            
            newJetP4.SetVectM(newJetP4.Vect(), newJetP4.Mag()*massCorr*masssmear) ;
            
#if DEBUGMORE
            cout
            << " \njet pt massCorrSmear = " << newJetP4.Pt()
            << " \njet mass massCorrSmear " << newJetP4.Mag()
            << endl ;
#endif
            
            if (scaleJetP4_ ) newJetP4 *= scaledJetMass_/newJetP4.Mag() ;
            
#if DEBUGMORE
            cout
            << " \njet pt scaleJetP4    = " << newJetP4.Pt()
            << " \njet mass scaleJetP4  = " << newJetP4.Mag()
            << endl ;
#endif
            
            jettau1         = (h_jettau1.product())->at(ijet) ;
            jettau2         = (h_jettau2.product())->at(ijet) ;
            jettau3         = (h_jettau3.product())->at(ijet) ;
            jetMass         = newJetP4.Mag();
            jetPrunedMass   = masssmear * massCorr * (1 + jecShift_*unc) * (h_jetPrunedMass.product())->at(ijet) ;
            jetTrimmedMass  = masssmear * massCorr * (1 + jecShift_*unc) * (h_jetTrimmedMass.product())->at(ijet) ;
            jetSoftDropMass = masssmear * massCorr * (1 + jecShift_*unc) * (h_jetSoftDropMass.product())->at(ijet) ;
            vjetsjIdx0      = int((h_vjetsjIdx0.product())->at(ijet)) ;
            vjetsjIdx1      = int((h_vjetsjIdx1.product())->at(ijet)) ;
            if (vjetsjIdx0 >= 0 ) {
                vjetssj0HadFl   = (h_ak8sjHadFl.product())->at(vjetsjIdx0) ;
                vjetssj0Pt      = (h_ak8sjPt.product())->at(vjetsjIdx0) ;
                vjetssj0Eta     = (h_ak8sjEta.product())->at(vjetsjIdx0) ;
                vjetssj0Phi     = (h_ak8sjPhi.product())->at(vjetsjIdx0) ;
           //     vjetssj0Mass    = (h_ak8sjMass.product())->at(vjetsjIdx0) ;
                vjetssj0CSV     = (h_ak8sjCSV.product())->at(vjetsjIdx0) ;
                vjetssj0CMVA    = (h_ak8sjCMVA.product())->at(vjetsjIdx0) ;
                p4sj0.SetPtEtaPhiM(vjetssj0Pt, vjetssj0Eta, vjetssj0Phi, vjetssj0Mass) ;
                p4sj0 *= ptsmear * newJEC  * (1 + jecShift_*unc);
                p4gensj0.SetPtEtaPhiE((h_ak8sjGenJetPt.product())->at(vjetsjIdx0),
                                      (h_ak8sjGenJetEta.product())->at(vjetsjIdx0),
                                      (h_ak8sjGenJetPhi.product())->at(vjetsjIdx0),
                                      (h_ak8sjGenJetE.product())->at(vjetsjIdx0));
            }
            if (vjetsjIdx1 >= 0) {
                vjetssj1HadFl   = (h_ak8sjHadFl.product())->at(vjetsjIdx1) ;
                vjetssj1Pt      = (h_ak8sjPt.product())->at(vjetsjIdx1) ;
                vjetssj1Eta     = (h_ak8sjEta.product())->at(vjetsjIdx1) ;
                vjetssj1Phi     = (h_ak8sjPhi.product())->at(vjetsjIdx1) ;
             //   vjetssj1Mass    = (h_ak8sjMass.product())->at(vjetsjIdx1) ;
                vjetssj1CSV     = (h_ak8sjCSV.product())->at(vjetsjIdx1) ;
                vjetssj1CMVA    = (h_ak8sjCMVA.product())->at(vjetsjIdx1) ;
                p4sj1.SetPtEtaPhiM(vjetssj1Pt, vjetssj1Eta, vjetssj1Phi, vjetssj1Mass) ;
                p4sj1 *= ptsmear * newJEC  * (1 + jecShift_*unc);
                p4gensj1.SetPtEtaPhiE((h_ak8sjGenJetPt.product())->at(vjetsjIdx1),
                                      (h_ak8sjGenJetEta.product())->at(vjetsjIdx1),
                                      (h_ak8sjGenJetPhi.product())->at(vjetsjIdx1),
                                      (h_ak8sjGenJetE.product())->at(vjetsjIdx1));
            }
            
            if ( vjetssj0CSV >= btaggedcsvlOP_ ) ++nsubjetsbtaggedcsvl ;
            if ( vjetssj1CSV >= btaggedcsvlOP_ ) ++nsubjetsbtaggedcsvl ;
            
            double subjetHighestCSV = std::max(vjetssj0CSV,vjetssj1CSV) ;
            
            if ( abs(jettau2/jettau1) >= idxjettau2Bytau1Min_
                && abs(jettau2/jettau1) <  idxjettau2Bytau1Max_
                && abs(jettau3/jettau2) >= idxjettau3Bytau2Min_
                && abs(jettau3/jettau2) <  idxjettau3Bytau2Max_
                && abs(jetMass) >= idxjetMassMin_
                && abs(jetMass) <  idxjetMassMax_
                && abs(jetPrunedMass) >= idxjetPrunedMassMin_
                && abs(jetPrunedMass) <  idxjetPrunedMassMax_
                && abs(jetTrimmedMass) >=  idxjetTrimmedMassMin_
                && abs(jetTrimmedMass) < idxjetTrimmedMassMax_
                && abs(jetSoftDropMass) >=  idxjetSoftDropMassMin_
                && abs(jetSoftDropMass) < idxjetSoftDropMassMax_
                && p4sj0.Mag() >=  idxsjMassMin_
                && p4sj1.Mag() >= idxsjMassMin_
                && vjetssj0CSV >=  idxsjCSVMin_
                && vjetssj0CSV < idxsjCSVMax_
                && vjetssj1CSV >= idxsjCSVMin_
                && vjetssj1CSV <  idxsjCSVMax_
                && subjetHighestCSV >= idxsjHighestCSVMin_
                && subjetHighestCSV <  idxsjHighestCSVMax_
                ) {
                jet.setTau1         ( (h_jettau1.product())->at(ijet) ) ;
                jet.setTau2         ( (h_jettau2.product())->at(ijet) ) ;
                jet.setTau3         ( (h_jettau3.product())->at(ijet) ) ;
                jet.setFilteredMass ( ptsmear * massCorr * (1 + jecShift_*unc) * h_jetFilteredMass.product()->at(ijet) ) ;
                jet.setTrimmedMass  ( ptsmear * massCorr * (1 + jecShift_*unc) * h_jetTrimmedMass.product()->at(ijet) ) ;
                jet.setPrunedMass   ( ptsmear * massCorr * (1 + jecShift_*unc) * h_jetPrunedMass.product()->at(ijet) ) ;
                jet.setSoftDropMass ( ptsmear * massCorr * (1 + jecShift_*unc) * h_jetSoftDropMass.product()->at(ijet) ) ;
                jet.setHadronFlavourSubjet0   ( vjetssj0HadFl ) ;
                jet.setHadronFlavourSubjet1   ( vjetssj1HadFl ) ;
                jet.setPtSubjet0    ( p4sj0.Pt() ) ;
                jet.setPtSubjet1    ( p4sj1.Pt() ) ;
                jet.setEtaSubjet0   ( p4sj0.Eta() ) ;
                jet.setEtaSubjet1   ( p4sj1.Eta() ) ;
                jet.setCSVSubjet0   ( vjetssj0CSV ) ;
                jet.setCSVSubjet1   ( vjetssj1CSV ) ;
                jet.setCMVASubjet0  ( vjetssj0CMVA ) ;
                jet.setCMVASubjet1  ( vjetssj1CMVA ) ;
              //  jet.setNConsts      ( (h_jetnumDaughters.product())->at(ijet) ) ;
                jet.setGroomedMassCorr (masssmear * massCorr * (1 + jecShift_*unc) ) ;
                jet.setNSubjetsBTaggedCSVL (nsubjetsbtaggedcsvl) ;
            }
            else {
                continue ;
            }
            
        } //// if type_ == AK8JET
        
        jet.setP4           (newJetP4) ;
        jet.setIndex        (ijet)  ;
        jet.setPartonFlavour( (h_jetPartonFlavour.product())->at(ijet) ) ;
        jet.setHadronFlavour( (h_jetHadronFlavour.product())->at(ijet) ) ;
        jet.setCSV          ( (h_jetCSV.product())->at(ijet) ) ;
        jet.setCMVA         ( (h_jetCMVA.product())->at(ijet) ) ;
        jet.setCvsB         ( (h_jetCvsB.product())->at(ijet) ) ;
        jet.setCvsL         ( (h_jetCvsL.product())->at(ijet) ) ;
      //  jet.setDoubleBAK8   ( (h_jetDoubleBAK8.product())->at(ijet) ) ;
      //  jet.setDoubleBCA15  ( (h_jetDoubleBCA15.product())->at(ijet) ) ;
        jet.setGenJetP4     (genJetP4);
        jets.push_back(jet) ;
        
        
#if DEBUGMORE
        cout
        << " \njet pt finally       = " << jet.getPt()
        << " \njet mass finally     = " << jet.getMass()
        << endl ;
#endif
        
    } //// loop over jets
    
    std::sort(jets.begin(), jets.end(), Utilities::sortByPt<vlq::Jet>) ;
    
}