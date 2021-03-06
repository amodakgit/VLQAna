import FWCore.ParameterSet.Config as cms

from Analysis.VLQAna.Electron_cfi import * 

defaultElectronMakerParameters = cms.PSet(
    defaultElectronParameters, 
    elidtype = cms.string("TIGHT"),
    elPtMin = cms.double(40),
    elPtMax = cms.double(10000),
    elAbsEtaMax = cms.double(2.1),
    )

looseElectronMakerParameters = cms.PSet(
    defaultElectronParameters, 
    elidtype = cms.string("LOOSE"),
    elPtMin = cms.double(15),
    elPtMax = cms.double(10000),
    elAbsEtaMax = cms.double(2.1),
    )
