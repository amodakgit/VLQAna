#ifndef ANALYSIS_VLQANA_APPLYLEPTONIDSFS_HH
#define ANALYSIS_VLQANA_APPLYLEPTONIDSFS_HH
#include <TH2.h>

class ApplyLeptonIDSFs {
  public:
    enum LEPTONIDTYPES_t {LOOSE, TIGHT} ; 
    ApplyLeptonIDSFs (edm::ParameterSet const& iConfig) : 
      sf_(1),
      zdecayMode_(iConfig.getParameter<std::string>("zdecayMode")),
      elidsfmap_(zdecayMode_ == "wel" ? iConfig.getParameter<std::string>("elidsfmap") : ""),
      elrecosfmap_(zdecayMode_ == "wel" ? iConfig.getParameter<std::string>("elrecosfmap") : "")
  {
    std::string lepidtypestr = iConfig.getParameter<std::string>("lepidtype") ; 
    if ( lepidtypestr == "LOOSE" ) type_ = LOOSE ;  
    else if ( lepidtypestr == "TIGHT" ) type_ = TIGHT ; 
    else edm::LogError("ApplyLeptonSF") << " >>>> WrongElectronIdType: " << type_<< " Check lepton id type !!!" ; 
  }
    ~ApplyLeptonIDSFs () {} 
    double IDSF (double pt, double eta){
      if (type_ == TIGHT && zdecayMode_ == "wel"){
        if(pt < 0.) pt = 10.1;
        if(pt > 500.) pt = 499.9;
        std::unique_ptr<TFile>f_idsfmap = std::unique_ptr<TFile>(new TFile(elidsfmap_.c_str())) ;
        std::unique_ptr<TH2>h2_effmap = std::unique_ptr<TH2>(dynamic_cast<TH2*>(f_idsfmap->Get("EGamma_SF2D"))) ;
        int binpt = h2_effmap->GetYaxis()->FindBin(pt);
        int bineta = h2_effmap->GetXaxis()->FindBin(eta);
        sf_ = h2_effmap->GetBinContent(bineta, binpt) ;
      }//end TIGHT and zelel

      if (type_ == TIGHT && zdecayMode_ == "wmu"){
        if(pt > 120.) pt = 120.;
        if (pt > 20 && pt <= 25){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.975305;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.986312;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.972437;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.986361;
        }
        else if (pt > 25 && pt <= 30){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.973942;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.986158;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.970745;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.983009;
        }
        else if (pt > 30 && pt <= 40){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.969471;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.987881;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.973607;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.984386;
        }
        else if (pt > 40 && pt <= 50){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.972574;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.990532;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.974868;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.9859;
        }
        else if (pt > 50 && pt <= 60){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.967713;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.985438;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.974566;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.981911;
        }
        else if (pt > 60 && pt <= 120){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.963148;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.988777;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.974918;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.992158;
        }
      }//TIGHT and zmumu
      return sf_ ; 
    }

    double IsoSF(double pt, double eta){
      if (type_ == TIGHT && zdecayMode_ == "wmu"){
        if(pt > 120.) pt = 120.;
        if (pt > 20 && pt <= 25){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.983849;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.991505;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.995189;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.98476;
        }
        else if (pt > 25 && pt <= 30){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.991877;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.995906;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 1.00043;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.993309;
        }
        else if (pt > 30 && pt <= 40){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.996227;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.997965;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.999432;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.993728;
        }
        else if (pt > 40 && pt <= 50){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.998324;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.998013;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.997695;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.995297;
        }
        else if (pt > 50 && pt <= 60){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 0.998603;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.998345;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.999144;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.996805;
        }
        else if (pt > 60 && pt <= 120){
          if (std::abs(eta) <= 2.4 && std::abs(eta) > 2.1)       sf_ = 1.00152;
          else if (std::abs(eta) <= 2.1 && std::abs(eta) > 1.2)  sf_ = 0.999213;
          else if (std::abs(eta) <= 1.2 && std::abs(eta) > 0.9)  sf_ = 0.999122;
          else if (std::abs(eta) <= 0.9 && std::abs(eta) > 0.0)  sf_ = 0.99878;
        }
      }//TIGHT and zmumu                                                                                                                                                                                                       
      return sf_ ;
    }
    double RECOSF (double pt, double eta){
      if (type_ == TIGHT && zdecayMode_ == "wel"){
        if(pt < 0.) pt = 10.1;
        if(pt > 500.) pt = 499.9;
        std::unique_ptr<TFile>f_recosfmap = std::unique_ptr<TFile>(new TFile(elrecosfmap_.c_str())) ;
        std::unique_ptr<TH2>h2_effmap = std::unique_ptr<TH2>(dynamic_cast<TH2*>(f_recosfmap->Get("EGamma_SF2D"))) ;
        int binpt = h2_effmap->GetYaxis()->FindBin(pt);
        int bineta = h2_effmap->GetXaxis()->FindBin(eta);
        sf_ = h2_effmap->GetBinContent(bineta, binpt) ;
      }//end TIGHT and zelel
      return sf_ ; 
    }
  private: 
    double sf_ ;
    const std::string zdecayMode_ ; 
    LEPTONIDTYPES_t type_ ; 
    const std::string elidsfmap_ ;
    const std::string elrecosfmap_ ;

};
#endif 
