#include <iostream>

#include "Analysis/VLQAna/interface/PickGenPart.h"

using namespace edm;

PickGenPart::PickGenPart (const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC) :
  t_genparticles       (iC.consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  ids_                 (iConfig.getParameter<std::vector<int>>("ids")), 
  statuses_            (iConfig.getParameter<std::vector<int>>("statuses")), 
  checkstatus_         (iConfig.getParameter<bool>("checkstatus")), 
  mom0ids_             (iConfig.getParameter<std::vector<int>>("mom0ids")), 
  mom1ids_             (iConfig.getParameter<std::vector<int>>("mom1ids")), 
  checkmomid_          (iConfig.getParameter<bool>("checkmomid")),
  dau0ids_             (iConfig.getParameter<std::vector<int>>("dau0ids")), 
  dau1ids_             (iConfig.getParameter<std::vector<int>>("dau1ids")), 
  checkdauid_          (iConfig.getParameter<bool>("checkdauid")),
  debug_               (iConfig.getParameter<bool>("debug"))
{ 
}

PickGenPart::~PickGenPart () {}

const GenParticleCollection PickGenPart::operator() ( edm::Event& evt) {
    
    using namespace std;
    
    typedef edm::Handle <std::vector<float>> hfloat ;
    
    edm::Handle<std::vector<reco::GenParticle>> h_genparticles;
    evt.getByToken(t_genparticles, h_genparticles) ;
    
    genParts_.clear();
    
    for ( reco::GenParticle igen : *(h_genparticles.product())) {
        
        int pdgId(igen.pdgId()) ;
        if ( std::find(ids_.begin(), ids_.end(), pdgId) == ids_.end() ) continue ;
        
        int status(igen.status()) ;
        if ( checkstatus_ && std::find(statuses_.begin(), statuses_.end(), status) == statuses_.end() ) continue ;
        
        int mom0PdgId(0), mom1PdgId(0), dau0PdgId(0), dau1PdgId(0);
        int mom0Status(-1), mom1Status(-1), dau0Status(-1), dau1Status(-1);
        
        if (igen.numberOfMothers() > 0 && igen.mother(0) != nullptr) {
            mom0PdgId = igen.mother(0)->pdgId() ;
            mom0Status = igen.mother(0)->pdgId() ;
        }
        
        if ( checkmomid_ && std::find(mom0ids_.begin(), mom0ids_.end(), mom0PdgId) == mom0ids_.end() ) continue ;
        
        int charge(igen.charge()) ;
        
        if (igen.numberOfMothers() > 1 && igen.mother(1) != nullptr) {
            mom1PdgId = igen.mother(1)->pdgId() ;
            mom1Status = igen.mother(1)->pdgId() ;
        }
        
        if (igen.numberOfDaughters() > 0 && igen.daughter(0) != nullptr) {
            dau0PdgId = igen.daughter(0)->pdgId() ;
            dau0Status = igen.daughter(0)->pdgId() ;
        }
        
        if (igen.numberOfDaughters() > 1 && igen.daughter(1) != nullptr) {
            dau1PdgId = igen.daughter(1)->pdgId() ;
            dau1Status = igen.daughter(1)->pdgId() ;
        }
        
        TLorentzVector p4genpart ;
        p4genpart.SetPtEtaPhiM(igen.pt(), igen.eta(), igen.phi(), igen.mass());
        vlq::GenParticle genpart ;
        genpart.setP4( p4genpart ) ;
        genpart.setPdgID      ( pdgId      ) ;
        genpart.setStatus     ( status     ) ;
        genpart.setCharge     ( charge     ) ;
        genpart.setMom0PdgID  ( mom0PdgId  ) ;
        genpart.setMom1PdgID  ( mom1PdgId  ) ;
        genpart.setDau0PdgID  ( dau0PdgId  ) ;
        genpart.setDau1PdgID  ( dau1PdgId  ) ;
        genpart.setMom0Status ( mom0Status ) ;
        genpart.setMom1Status ( mom1Status ) ;
        genpart.setDau0Status ( dau0Status ) ;
        genpart.setDau1Status ( dau1Status ) ;
        genParts_.push_back(genpart) ;
        
    }

    return genParts_ ; 
}
