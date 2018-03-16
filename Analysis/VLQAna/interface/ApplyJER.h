#ifndef ANALYSIS_INTERFACE_APPLYJER_H
#define ANALYSIS_INTERFACE_APPLYJER_H

double ApplyJERMass() {
  double massjerscale(1.12) ; 
  return massjerscale ; 
}

double ApplyJERp4 (double eta, int jerShift) {

  eta = std::abs(eta) ;
  if (eta >= 4.7) eta = 4.699 ;
  double jerscale(0) ; 

  if (eta >= 3.2 && eta < 5.0) {
    if (jerShift == 1) jerscale = 1.16  ;  
    else if (jerShift == 2) jerscale = 1.16 + 0.029 ;  
    else if (jerShift == -1) jerscale = 1.16 - 0.029 ;  
  }
  else if (eta >= 3.0 && eta < 3.2) {
    if (jerShift == 1) jerscale = 1.328  ;  
    else if (jerShift == 2) jerscale = 1.328 + 0.022;  
    else if (jerShift == -1) jerscale = 1.328 - 0.022;  
  }
  else if (eta >= 2.8 && eta < 3.0) {
    if (jerShift == 1) jerscale =  1.857 ;  
    else if (jerShift == 2) jerscale =  1.857 + 0.071;  
    else if (jerShift == -1) jerscale =  1.857 - 0.071;  
  }
  else if (eta >= 2.5 && eta < 2.8) {
    if (jerShift == 1) jerscale = 1.364  ;  
    else if (jerShift == 2) jerscale = 1.364 + 0.039;  
    else if (jerShift == -1) jerscale = 1.364 - 0.039;  
  }
  else if (eta >= 2.3 && eta < 2.5) {
    if (jerShift == 1) jerscale = 1.177 ;  
    else if (jerShift == 2) jerscale = 1.177 + 0.041;  
    else if (jerShift == -1) jerscale = 1.177 - 0.041;  
  }
  else if (eta >= 2.1 && eta < 2.3) {
    if (jerShift == 1) jerscale = 1.067 ;  
    else if (jerShift == 2) jerscale = 1.067 + 0.053;  
    else if (jerShift == -1) jerscale = 1.067 - 0.053 ;  
  }
  else if (eta >= 1.9 && eta < 2.1) {
    if (jerShift == 1) jerscale = 1.140  ;  
    else if (jerShift == 2) jerscale = 1.140 + 0.047;  
    else if (jerShift == -1) jerscale = 1.140 - 0.047;  
  }
  else if (eta >= 1.7 && eta < 1.9) {
    if (jerShift == 1) jerscale = 1.082  ;  
    else if (jerShift == 2) jerscale = 1.082 + 0.035;  
    else if (jerShift == -1) jerscale = 1.082 - 0.035 ;  
  }
  else if (eta >= 1.3 && eta < 1.7) {
    if (jerShift == 1) jerscale = 1.084  ;  
    else if (jerShift == 2) jerscale = 1.084 + 0.011;  
    else if (jerShift == -1) jerscale = 1.084 - 0.011 ;  
  }
  else if (eta >= 1.1 && eta < 1.3) {
    if (jerShift == 1) jerscale = 1.123  ;  
    else if (jerShift == 2) jerscale = 1.123 + 0.024;  
    else if (jerShift == -1) jerscale = 1.123 - 0.024;  
  }
  else if (eta >= 0.8 && eta < 1.1) {
    if (jerShift == 1) jerscale = 1.114  ;  
    else if (jerShift == 2) jerscale = 1.114 + 0.013;  
    else if (jerShift == -1) jerscale = 1.114 - 0.013;  
  }
  else if (eta >= 0.5 && eta < 0.8) {
    if (jerShift == 1) jerscale = 1.138  ;  
    else if (jerShift == 2) jerscale = 1.138 + 0.013;  
    else if (jerShift == -1) jerscale = 1.138 - 0.013;  
  }
  else if (eta >= 0.0 && eta < 0.5) {
    if (jerShift == 1) jerscale = 1.109  ;  
    else if (jerShift == 2) jerscale = 1.109 + 0.008;  
    else if (jerShift == -1) jerscale = 1.109 - 0.008;  
  }

  return jerscale ;

 // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
}

#endif 
