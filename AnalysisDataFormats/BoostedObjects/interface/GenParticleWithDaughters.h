#ifndef ANALYSISDATAFORMATS_BOOSTEDOBJECTS_GENPARTICLEWITHDAUGHTERS_HH
#define ANALYSISDATAFORMATS_BOOSTEDOBJECTS_GENPARTICLEWITHDAUGHTERS_HH

#include "AnalysisDataFormats/BoostedObjects/interface/GenParticle.h" 

#include <TLorentzVector.h>

namespace vlq {
  class GenParticleWithDaughters ; 
  typedef std::vector<GenParticleWithDaughters> GenParticleWithDaughtersCollection ;  
  class GenParticleWithDaughters : public Candidate {
    friend class vlq::GenParticle ; 
    public:
      GenParticleWithDaughters () {}
      ~GenParticleWithDaughters () {} 
      GenParticleWithDaughters (const GenParticleWithDaughters& p) :
        moms(p.getMothers()),
        daus(p.getDaughters())
      {}
      void setMothers ( const vlq::GenParticleCollection& ps ) { 
        GenParticleCollection::const_iterator it ;
        for ( it = ps.begin(); it != ps.end(); ++it ) {
          GenParticle mom(*it) ; 
         moms.push_back(mom) ; 
        }
      } 
      void setDaughters ( const vlq::GenParticleCollection& ps ) { 
        GenParticleCollection::const_iterator it ;
        for ( it = ps.begin(); it != ps.end(); ++it ) {
          GenParticle dau(*it) ; 
          daus.push_back(dau) ; 
        }
      } 
      vlq::GenParticleCollection getMothers () const { return  moms; }
      vlq::GenParticleCollection getDaughters () const { return daus ; }
    private: 
      vlq::GenParticleCollection moms ; 
      vlq::GenParticleCollection daus ; 
  }; 

}

#endif 


