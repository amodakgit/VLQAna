#ifndef ANALYSISDATAFORMATS_BOOSTEDOBJECTS_GENPARTICLE_HH
#define ANALYSISDATAFORMATS_BOOSTEDOBJECTS_GENPARTICLE_HH

#include "AnalysisDataFormats/BoostedObjects/interface/Candidate.h"

namespace vlq {
  class GenParticle ;
  typedef std::vector<GenParticle> GenParticleCollection ; 
  class GenParticle : public Candidate {
    public: 
      GenParticle () : 
        pdgID_(-999999), 
        status_(-999999), 
        charge_(-999999), 
        mom0PdgID_(-999999), 
        mom1PdgID_(-999999), 
        dau0PdgID_(-999999), 
        dau1PdgID_(-999999), 
        mom0Status_(-999999), 
        mom1Status_(-999999), 
        dau0Status_(-999999), 
        dau1Status_(-999999)  
    {
    } 
      ~GenParticle (){} 
      GenParticle (const vlq::GenParticle& p ) :
        pdgID_(p.getPdgID()),
        status_(p.getStatus()), 
        charge_(p.getCharge()),
        mom0PdgID_(p.getMom0PdgID()),
        mom1PdgID_(p.getMom1PdgID()),
        dau0PdgID_(p.getDau0PdgID()),
        dau1PdgID_(p.getDau1PdgID()),
        mom0Status_(p.getMom0Status()),
        mom1Status_(p.getMom1Status()),
        dau0Status_(p.getDau0Status()),
        dau1Status_(p.getDau1Status()) 
    { 
      this->setP4( p.getP4() ) ; 
    } 

      int getPdgID () const { return pdgID_ ; }
      int getStatus () const { return status_ ; } 
      int getCharge () const {return charge_ ; }
      int getMom0PdgID () const { return mom0PdgID_ ; }
      int getMom1PdgID () const { return mom1PdgID_ ; }
      int getDau0PdgID () const { return dau0PdgID_ ; }
      int getDau1PdgID () const { return dau1PdgID_ ; }
      int getMom0Status () const { return mom0Status_ ; }
      int getMom1Status () const { return mom1Status_ ; }
      int getDau0Status () const { return dau0Status_ ; }
      int getDau1Status () const { return dau1Status_ ; }

      void setPdgID ( const int pdgID )  { pdgID_ = pdgID ; }
      void setStatus ( const int status ) { status_ = status ; } 
      void setCharge ( const int charge ) {charge_ = charge ; }
      void setMom0PdgID ( const int mom0PdgID ) { mom0PdgID_ = mom0PdgID ; }
      void setMom1PdgID ( const int mom1PdgID ) { mom1PdgID_ = mom1PdgID ; }
      void setDau0PdgID ( const int dau0PdgID ) { dau0PdgID_ = dau0PdgID ; }
      void setDau1PdgID ( const int dau1PdgID ) { dau1PdgID_ = dau1PdgID ; }
      void setMom0Status ( const int mom0Status ) { mom0Status_ = mom0Status ; }
      void setMom1Status ( const int mom1Status ) { mom1Status_ = mom1Status ; }
      void setDau0Status ( const int dau0Status ) { dau0Status_ = dau0Status ; }
      void setDau1Status ( const int dau1Status ) { dau1Status_ = dau1Status ; }

    protected:
      int pdgID_ ;
      int status_ ; 
      int charge_ ; 
      int mom0PdgID_ ;
      int mom1PdgID_ ;
      int dau0PdgID_ ;
      int dau1PdgID_ ;
      int mom0Status_ ;
      int mom1Status_ ;
      int dau0Status_ ;
      int dau1Status_ ;
  }; 
}
#endif
