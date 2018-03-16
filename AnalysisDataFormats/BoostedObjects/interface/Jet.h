#ifndef ANALYSISDATAFORMATS_BOOSTEDOBJECTS_JET_HH
#define ANALYSISDATAFORMATS_BOOSTEDOBJECTS_JET_HH

#include "AnalysisDataFormats/BoostedObjects/interface/ResolvedVjj.h"

#include<iostream>

/*\ 
 * Simple Jet class to be used with B2G EDMNtuples
 * https://github.com/cmsb2g/B2GAnaFW
 */

namespace vlq{
  class Jet ;
  typedef std::vector<vlq::Jet> JetCollection ; 
  class Jet : public vlq::Candidate {
    friend class vlq::ResolvedVjj ; 
    public:
      Jet () : 
        index_(-1), 
        partonFlavour_(-1), 
        hadronFlavour_(-1), 
        csvDiscrimiator_(-1000), 
        cmvaDiscrimiator_(-1000), 
        CvsBDiscrimiator_(-1000), 
        CvsLDiscrimiator_(-1000), 
        doubleBAK8Discrimiator_(-1000), 
        doubleBCA15Discrimiator_(-1000), 
        nhf_ (-1000), 
        chf_ (-1000), 
        emf_ (-1000), 
        phf_ (-1000), 
        muf_ (-1000), 
        nconsts_(-1000),
        tau1_(-1000), 
        tau2_(-1000), 
        tau3_(-1000), 
        filteredMass_(-1000), 
        trimmedMass_(-1000), 
        prunedMass_(-1000), 
        softDropMass_(-1000), 
        nsubjets_  (-1000), 
        hadronFlavourSubjet0_(-1000), 
        hadronFlavourSubjet1_(-1000), 
        ptSubjet0_(-1000), 
        ptSubjet1_(-1000), 
        etaSubjet0_(-1000), 
        etaSubjet1_(-1000), 
        phiSubjet0_(-1000), 
        phiSubjet1_(-1000), 
        energySubjet0_(-1000), 
        energySubjet1_(-1000), 
        csvSubjet0_(-1000), 
        csvSubjet1_(-1000), 
        cmvaSubjet0_(-1000), 
        cmvaSubjet1_(-1000), 
        qSubjet0_(-1000), 
        qSubjet1_(-1000), 
        groomedMassCorr_(-1000), 
        nsubjetsBTaggedCSVL_(-1000), 
        isbtagged_(0), 
        istoptagged_(0), 
        ishtagged_(0), 
        iswtagged_(0),
        genJetP4_(-1000,-1000,-1000,-1000),
        genSubJet0P4_(-1000,-1000,-1000,-1000),
        genSubJet1P4_(-1000,-1000,-1000,-1000)
    {}
      Jet (const vlq::Jet& jet) :
        vlq::Candidate(jet),
        index_(jet.getIndex()), 
        partonFlavour_(jet.getPartonFlavour()),
        hadronFlavour_(jet.getHadronFlavour()), 
        csvDiscrimiator_(jet.getCSV()),
        cmvaDiscrimiator_(jet.getCMVA()),
        CvsBDiscrimiator_(jet.getCvsB()),
        CvsLDiscrimiator_(jet.getCvsL()),
        doubleBAK8Discrimiator_(jet.getDoubleBAK8()),
        doubleBCA15Discrimiator_(jet.getDoubleBCA15()),
        nhf_    (jet.getNHF     ()), 
        chf_    (jet.getCHF     ()), 
        emf_    (jet.getEMF     ()), 
        phf_    (jet.getPHF     ()), 
        muf_    (jet.getMUF     ()), 
        nconsts_(jet.getNConsts ()), 
        tau1_(jet.getTau1()),
        tau2_(jet.getTau2()),
        tau3_(jet.getTau3()),
        filteredMass_(jet.getFilteredMass()),
        trimmedMass_(jet.getTrimmedMass()),
        prunedMass_(jet.getPrunedMass()),
        softDropMass_(jet.getSoftDropMass()),
        nsubjets_  (jet.getNSubjets()),  
        hadronFlavourSubjet0_(jet.getHadronFlavourSubjet0()),  
        hadronFlavourSubjet1_(jet.getHadronFlavourSubjet1()),  
        ptSubjet0_(jet.getPtSubjet0()),  
        ptSubjet1_(jet.getPtSubjet1()),  
        etaSubjet0_(jet.getEtaSubjet0()),  
        etaSubjet1_(jet.getEtaSubjet1()),  
        phiSubjet0_(jet.getPhiSubjet0()),  
        phiSubjet1_(jet.getPhiSubjet1()),  
        energySubjet0_(jet.getEnergySubjet0()),  
        energySubjet1_(jet.getEnergySubjet1()),  
        csvSubjet0_(jet.getCSVSubjet0()),  
        csvSubjet1_(jet.getCSVSubjet1()),  
        cmvaSubjet0_(jet.getCMVASubjet0()),  
        cmvaSubjet1_(jet.getCMVASubjet1()),  
        qSubjet0_(jet.getQSubjet0()),  
        qSubjet1_(jet.getQSubjet1()),  
        groomedMassCorr_(jet.getGroomedMassCorr()),  
        nsubjetsBTaggedCSVL_(jet.getNSubjetsBTaggedCSVL()),  
        isbtagged_(jet.getIsbtagged()),
        istoptagged_(jet.getIstoptagged()),
        ishtagged_(jet.getIshtagged()),
        iswtagged_(jet.getIswtagged()),
        genJetP4_(jet.getGenJetP4()),
        genSubJet0P4_(jet.getGenJetP4()),
        genSubJet1P4_(jet.getGenJetP4())
    { } 
      ~Jet () {}

      int    getIndex                 () const { return index_ ; } 
      int    getPartonFlavour         () const { return partonFlavour_ ; } 
      int    getHadronFlavour         () const { return hadronFlavour_ ; } 
      double getCSV                   () const { return csvDiscrimiator_ ; }
      double getCMVA                  () const { return cmvaDiscrimiator_ ; }
      double getCvsB                  () const { return CvsBDiscrimiator_ ; }
      double getCvsL                  () const { return CvsLDiscrimiator_ ; }
      double getDoubleBAK8            () const { return doubleBAK8Discrimiator_ ; }
      double getDoubleBCA15           () const { return doubleBCA15Discrimiator_ ; }
      double getNHF                   () const { return nhf_     ;}
      double getCHF                   () const { return chf_     ;}
      double getEMF                   () const { return emf_     ;}
      double getPHF                   () const { return phf_     ;}
      double getMUF                   () const { return muf_     ;}
      int    getNConsts               () const { return nconsts_ ;}
      double getTau1                  () const { return tau1_                 ; }
      double getTau2                  () const { return tau2_                 ; }
      double getTau3                  () const { return tau3_                 ; }
      double getFilteredMass          () const { return filteredMass_         ; }
      double getTrimmedMass           () const { return  trimmedMass_         ; }
      double getPrunedMass            () const { return prunedMass_           ; }
      double getSoftDropMass          () const { return softDropMass_         ; }
      int    getNSubjets              () const { return nsubjets_             ; } 
      double getHadronFlavourSubjet0  () const { return hadronFlavourSubjet0_ ; } 
      double getHadronFlavourSubjet1  () const { return hadronFlavourSubjet1_ ; } 
      double getPtSubjet0             () const { return ptSubjet0_ ; } 
      double getPtSubjet1             () const { return ptSubjet1_ ; } 
      double getEtaSubjet0            () const { return etaSubjet0_ ; } 
      double getEtaSubjet1            () const { return etaSubjet1_ ; } 
      double getPhiSubjet0            () const { return phiSubjet0_ ; } 
      double getPhiSubjet1            () const { return phiSubjet1_ ; } 
      double getEnergySubjet0         () const { return energySubjet0_ ; } 
      double getEnergySubjet1         () const { return energySubjet1_ ; } 
      double getCSVSubjet0            () const { return csvSubjet0_ ; } 
      double getCSVSubjet1            () const { return csvSubjet1_ ; } 
      double getCMVASubjet0           () const { return cmvaSubjet0_ ; } 
      double getCMVASubjet1           () const { return cmvaSubjet1_ ; } 
      double getQSubjet0              () const { return qSubjet0_ ; } 
      double getQSubjet1              () const { return qSubjet1_ ; } 
      double getGroomedMassCorr       () const { return groomedMassCorr_ ; } 
      int    getNSubjetsBTaggedCSVL   () const { return nsubjetsBTaggedCSVL_ ; } 
      bool   getIsbtagged             () const { return isbtagged_ ; } 
      bool   getIstoptagged           () const { return istoptagged_ ; }
      bool   getIshtagged             () const { return ishtagged_ ; }
      bool   getIswtagged             () const { return iswtagged_ ; } 
      TLorentzVector getGenJetP4      () const { return genJetP4_  ; }
      TLorentzVector getGenSubJet0P4  () const { return genSubJet0P4_  ; }
      TLorentzVector getGenSubJet1P4  () const { return genSubJet1P4_  ; }

      void setIndex                 ( const int index                  ) {index_ = index ; } 
      void setPartonFlavour         ( const int partonflavour          ) { partonFlavour_ = partonflavour ; } 
      void setHadronFlavour         ( const int hadronflavour          ) { hadronFlavour_ = hadronflavour ; } 
      void setCSV                   ( const double csv                 ) {csvDiscrimiator_ = csv ; } 
      void setCMVA                  ( const double cmva                ) {cmvaDiscrimiator_ = cmva ; } 
      void setCvsB                  ( const double CvsB                ) {CvsBDiscrimiator_ = CvsB ; } 
      void setCvsL                  ( const double CvsL                ) {CvsLDiscrimiator_ = CvsL ; } 
      void setDoubleBAK8            ( const double doubleBAK8          ) {doubleBAK8Discrimiator_ = doubleBAK8 ; } 
      void setDoubleBCA15           ( const double doubleBCA15         ) {doubleBCA15Discrimiator_ = doubleBCA15 ; } 
      void setNHF                   ( const double nhf                 ) { nhf_  = nhf  ; }
      void setCHF                   ( const double chf                 ) { chf_  = chf  ; }
      void setEMF                   ( const double emf                 ) { emf_  = emf  ; }
      void setPHF                   ( const double phf                 ) { phf_  = phf  ; }
      void setMUF                   ( const double muf                 ) { muf_  = muf  ; }
      void setNConsts               ( const int    nconsts             ) { nconsts_  = nconsts  ; }
      void setTau1                  ( const double tau1                ) { tau1_ = tau1 ; }
      void setTau2                  ( const double tau2                ) { tau2_ = tau2 ; }
      void setTau3                  ( const double tau3                ) { tau3_ = tau3 ; }
      void setFilteredMass          ( const double mass                ) { filteredMass_ = mass ; }
      void setTrimmedMass           ( const double mass                ) { trimmedMass_ = mass ; }
      void setPrunedMass            ( const double mass                ) { prunedMass_ = mass ; }
      void setSoftDropMass          ( const double mass                ) { softDropMass_ = mass ; }
      void setNSubjets              ( const int    nsubjets            ) { nsubjets_ = nsubjets; }
      void setHadronFlavourSubjet0  ( const double hadronFlavourSubjet0) { hadronFlavourSubjet0_ = hadronFlavourSubjet0 ; }
      void setHadronFlavourSubjet1  ( const double hadronFlavourSubjet1) { hadronFlavourSubjet1_ = hadronFlavourSubjet1 ; }
      void setPtSubjet0             ( const double ptSubjet0           ) { ptSubjet0_ = ptSubjet0 ; }
      void setPtSubjet1             ( const double ptSubjet1           ) { ptSubjet1_ = ptSubjet1 ; }
      void setEtaSubjet0            ( const double etaSubjet0          ) { etaSubjet0_ = etaSubjet0 ; }
      void setEtaSubjet1            ( const double etaSubjet1          ) { etaSubjet1_ = etaSubjet1 ; }
      void setPhiSubjet0            ( const double phiSubjet0          ) { phiSubjet0_ = phiSubjet0 ; }
      void setPhiSubjet1            ( const double phiSubjet1          ) { phiSubjet1_ = phiSubjet1 ; }
      void setEnergySubjet0         ( const double energySubjet0       ) { energySubjet0_ = energySubjet0 ; }
      void setEnergySubjet1         ( const double energySubjet1       ) { energySubjet1_ = energySubjet1 ; }
      void setCSVSubjet0            ( const double csvSubjet0          ) { csvSubjet0_ = csvSubjet0 ; }
      void setCSVSubjet1            ( const double csvSubjet1          ) { csvSubjet1_ = csvSubjet1 ; }
      void setCMVASubjet0           ( const double cmvaSubjet0         ) { cmvaSubjet0_ = cmvaSubjet0 ; }
      void setCMVASubjet1           ( const double cmvaSubjet1         ) { cmvaSubjet1_ = cmvaSubjet1 ; }
      void setQSubjet0              ( const double qSubjet0            ) { qSubjet0_ = qSubjet0 ; }
      void setQSubjet1              ( const double qSubjet1            ) { qSubjet1_ = qSubjet1 ; }
      void setGroomedMassCorr       ( const double groomedMassCorr     ) { groomedMassCorr_ = groomedMassCorr ; }
      void setNSubjetsBTaggedCSVL   ( const int    nsubjetsBTaggedCSVL ) { nsubjetsBTaggedCSVL_ = nsubjetsBTaggedCSVL ; }
      void setIsbtagged             ( const bool   isbtagged           ) { isbtagged_ = isbtagged ; } 
      void setIstoptagged           ( const bool   istoptagged         ) { istoptagged_ = istoptagged ; }
      void setIshtagged             ( const bool   ishtagged           ) { ishtagged_ = ishtagged ; }
      void setIswtagged             ( const bool   iswtagged           ) { iswtagged_ = iswtagged ; } 
      void setGenJetP4              ( const TLorentzVector& p4         ) { genJetP4_ = p4 ; }
      void setGenSubJet0P4          ( const TLorentzVector& p4         ) { genSubJet0P4_ = p4 ; }
      void setGenSubJet1P4          ( const TLorentzVector& p4         ) { genSubJet1P4_ = p4 ; }

    protected:

      int index_ ; 
      int partonFlavour_ ;
      int hadronFlavour_ ; 
      double csvDiscrimiator_ ; 
      double cmvaDiscrimiator_ ; 
      double CvsBDiscrimiator_ ; 
      double CvsLDiscrimiator_ ; 
      double doubleBAK8Discrimiator_ ; 
      double doubleBCA15Discrimiator_ ; 

      double nhf_ ; 
      double chf_ ;
      double emf_ ;
      double phf_ ;
      double muf_ ;
      int    nconsts_ ; 

      double tau1_ ;
      double tau2_ ;
      double tau3_ ;
      double filteredMass_ ;
      double trimmedMass_ ;
      double prunedMass_ ; 
      double softDropMass_ ; 

      double nsubjets_ ; 
      double hadronFlavourSubjet0_ ; 
      double hadronFlavourSubjet1_ ; 
      double ptSubjet0_ ; 
      double ptSubjet1_ ; 
      double etaSubjet0_ ; 
      double etaSubjet1_ ; 
      double phiSubjet0_ ; 
      double phiSubjet1_ ; 
      double energySubjet0_ ; 
      double energySubjet1_ ; 
      double csvSubjet0_ ; 
      double csvSubjet1_ ; 
      double cmvaSubjet0_ ; 
      double cmvaSubjet1_ ; 
      double qSubjet0_ ; 
      double qSubjet1_ ; 

      double groomedMassCorr_ ; 
      int    nsubjetsBTaggedCSVL_ ; 

      bool isbtagged_ ; 
      bool istoptagged_ ;
      bool ishtagged_ ;
      bool iswtagged_ ; 

      TLorentzVector genJetP4_;
      TLorentzVector genSubJet0P4_;
      TLorentzVector genSubJet1P4_;
  }; 
}

#endif 
