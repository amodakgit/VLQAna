import FWCore.ParameterSet.Config as cms

from Analysis.VLQAna.Muon_cfi import * 

defaultMuonMakerParameters = cms.PSet(
    defaultMuonParameters, 
    muidtype = cms.string("TIGHT"),
    muPtMin = cms.double(60),
    muPtMax = cms.double(10000),
    #muAbsEtaMax = cms.double(2.1),
    muAbsEtaMax = cms.double(2.4),
    muIsoMin = cms.double(0.0),
    muIsoMax = cms.double(0.15),
    )

skimMuonMakerParameters = cms.PSet(
    defaultMuonParameters, 
    muidtype = cms.string("TIGHT"),
    muPtMin = cms.double(60),
    muPtMax = cms.double(10000),
    muAbsEtaMax = cms.double(2.4),
    muIsoMin = cms.double(0.0),
    muIsoMax = cms.double(9999.9),
    )

looseMuonMakerParameters = cms.PSet(
    defaultMuonParameters, 
    muidtype = cms.string("LOOSE"),
    muPtMin = cms.double(15),
    muPtMax = cms.double(10000),
    muAbsEtaMax = cms.double(2.4),
    muIsoMin = cms.double(0.0),
    muIsoMax = cms.double(0.15),
    )
