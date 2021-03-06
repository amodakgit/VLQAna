import os, sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('isData', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Is data?"
                 )
options.register('lepID', 'TIGHT',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "lepton ID? Choose: 'TIGHT' or 'LOOSE'"
                 )
'''
options.register('newJECPayloadNames', 'Summer16_23Sep2016V4_MC',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "Select one of Summer16_23Sep2016BCDV4_DATA, Summer16_23Sep2016EFV4_DATA, Summer16_23Sep2016GV4_DATA, Summer16_23Sep2016HV4_DATA, Summer16_23Sep2016V4_MC" 
                )
'''
options.register('outFileName', 'singlelep.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('doPUReweightingOfficial', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Do pileup reweighting using official recipe"
                 )
options.register('applyLeptonIDSFs', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Apply lepton SFs to the MC"
                 )
options.register('applyLeptonTrigSFs', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Apply trigger SFs to the MC"
                )
options.register('applyBTagSFs', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Apply b-tagging SFs to the MC"
                 )
options.register('applyDYNLOCorr', False, ### Set to true only for DY process ### Only EWK NLO k-factor is applied
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Apply DY EWK k-factor to DY MC"
                 )
options.register('applyTopPtCorr', False, ### Set to true only for ttbar process
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Apply ttbar top pt correction to ttbar MC"
                 )
options.register('applyWjetHtCorr', False, ### Set to true only for ttbar process
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Apply w+jet ht correction to w+jet MC"
                 )
options.register('jecShift', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Apply w+jet ht correction to w+jet MC"
                 )
options.register('jerShift', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Apply w+jet ht correction to w+jet MC"
                 )

options.register('FileNames', 'bprime800',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Name of list of input files"
                 )
options.register('optimizeReco', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Optimize mass reconstruction"
                 )
options.register('applyHtCorr', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Optimize mass reconstruction"
                 )
options.register('doSkim', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Produce skim 1 or 0"
                 )
options.register('sys', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Do systematics"
                 )
options.register('short', False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "only signal"
                 )

options.setDefault('maxEvents', -1)
options.parseArguments()
print options
dataPath = ''

hltpaths_mu = [
	       "HLT_IsoMu24_v",
               "HLT_IsoTkMu24_v",
               #"HLT_Mu50_v",
               #"HLT_TkMu50_v",
               ]
hltpaths_ele = [
                "HLT_Ele27_eta2p1_WPLoose_Gsf_v",
                "HLT_Ele27_eta2p1_WPTight_Gsf_v",
		"HLT_Ele32_eta2p1_WPTight_Gsf_v",
                ]

if options.isData:
    options.optimizeReco = False
    options.applyLeptonIDSFs = False
    options.applyLeptonTrigSFs = False
    options.applyHtCorr = False
    options.applyWjetHtCorr = False
    options.doPUReweightingOfficial = False
    options.jecShift = False
    options.jerShift = False
    options.sys = False
else:
    print "applyWjetHtCorr ", options.applyWjetHtCorr
    print "isData ", options.isData

process = cms.Process("SingleLepAna")

from inputFiles_cfi import *

process.source = cms.Source(
                            "PoolSource",
                            fileNames = cms.untracked.vstring(
                                                              
#'/store/group/phys_b2g/B2GAnaFW_80X_V2p4/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_B2GAnaFW_80X_V2p4/161222_110143/0000/B2GEDMNtuple_1.root',
'file:/afs/cern.ch/work/a/amodak/public/VLQ/B2G/Production/CMSSW_8_0_26_patch1/src/Analysis/B2GAnaFW/test/B2GEDMNtuple.root',
                               )
                            )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Analysis.VLQAna.EventCleaner_cff")
process.evtcleaner.File_PUDistData      = cms.string(os.path.join(dataPath,'RunII2016Rereco_25ns_PUXsec69000nb.root'))
process.evtcleaner.File_PUDistDataLow   = cms.string(os.path.join(dataPath,'RunII2016Rereco_25ns_PUXsec65550nb.root'))
process.evtcleaner.File_PUDistDataHigh  = cms.string(os.path.join(dataPath,'RunII2016Rereco_25ns_PUXsec72450nb.root'))
process.evtcleaner.File_PUDistMC        = cms.string(os.path.join(dataPath,'PUDistMC_Summer2016_25ns_Moriond17MC_PoissonOOTPU.root'))

process.evtcleaner.isData = options.isData
process.evtcleaner.hltMuPaths  = cms.vstring (hltpaths_mu)
process.evtcleaner.hltElePaths = cms.vstring (hltpaths_ele)
process.evtcleaner.DoPUReweightingOfficial = cms.bool(options.doPUReweightingOfficial)

#process.evtcleaner.storeLHEWts = options.storeLHEWts

from Analysis.VLQAna.SingleLepAna_cfi import *



### Z candidate and jet selections
process.ana = ana.clone(
                        applyLeptonIDSFs  = cms.bool(options.applyLeptonIDSFs),
                        applyLeptonTrigSFs  = cms.bool(options.applyLeptonTrigSFs),
                        applyBTagSFs    = cms.bool(options.applyBTagSFs),
                        applyDYNLOCorr  = cms.bool(options.applyDYNLOCorr),
                        optimizeReco    = cms.bool(options.optimizeReco),
                        applyHtCorr     = cms.bool(options.applyHtCorr),
		        applyWjetHtCorr = cms.bool(options.applyWjetHtCorr),
			applyTopPtCorr  = cms.bool(options.applyTopPtCorr),	
                        doSkim          = cms.bool(options.doSkim),
                        sys             = cms.bool(options.sys),
                        short           = cms.bool(options.short),
                        )
process.ana.elselParams.elidtype = cms.string(options.lepID)
process.ana.muselParams.muidtype = cms.string(options.lepID)
#process.ana.muselParams.muIsoMax = cms.double(0.15)
#zdecaytype is automatic inside cfi
process.ana.muIdSFsParams.lepidtype = cms.string(options.lepID)
process.ana.elIdSFsParams.lepidtype = cms.string(options.lepID)
process.ana.elIdSFsParams.elidsfmap = cms.string(os.path.join(dataPath,"egammaIDEffi.txt_EGM2D.root"))
process.ana.elIdSFsParams.elrecosfmap = cms.string(os.path.join(dataPath,"egammaRECOEffi.txt_EGM2D.root"))
process.ana.elTrigSFsParams.eltrigsfmap = cms.string(os.path.join(dataPath,"SF_HLT_Ele32_eta2p1.root"))

#process.ana.BoostedZCandParams.ptMin = cms.double(150.)#not used in analysis
process.ana.jetAK8selParams.jetPtMin = cms.double(200)
process.ana.jetAK4BTaggedselParams.jetPtMin = cms.double(50)
process.ana.STMin = cms.double(1000.)
process.ana.vlqMass = cms.double(1200.) #M=1000
process.ana.bosonMass = cms.double(91.2) #Z
process.ana.recoPt = cms.double(0.)

process.ana.jetAK4selParams.jecShift = options.jecShift 
process.ana.jetAK4selParams.jerShift = options.jerShift 
process.ana.jetAK4BTaggedselParams.jecShift = options.jecShift 
process.ana.jetAK4BTaggedselParams.jerShift = options.jerShift 
process.ana.jetAK8selParams.jecShift = options.jecShift 
process.ana.jetAK8selParams.jerShift = options.jerShift
'''
#JEC Uncertainty
process.ana.jetAK4selParams.jecUncPayloadName = cms.string(os.path.join(dataPath,options.newJECPayloadNames+"_Uncertainty_AK4PFchs.txt"))
process.ana.jetAK4BTaggedselParams.jecUncPayloadName = cms.string(os.path.join(dataPath,options.newJECPayloadNames+"_Uncertainty_AK4PFchs.txt"))
process.ana.jetAK8selParams.jecUncPayloadName = cms.string(os.path.join(dataPath,options.newJECPayloadNames+"_Uncertainty_AK8PFchs.txt"))


#Mass Correction
masscorr = ['L2Relative', 'L3Absolute']
for i,c in enumerate(masscorr): masscorr[i] = os.path.join(dataPath,options.newJECPayloadNames+"_"+c+"_AK8PFchs.txt")
print "jec ak8 groom payload ", masscorr
process.ana.jetAK8selParams.jecAK8GroomedPayloadNames = cms.vstring(masscorr)


#JEC (ak4 & ak8)
if options.newJECPayloadNames != '':
  corrections = ['L1FastJet', 'L2Relative', 'L3Absolute']
  ak4chsCorr = []
  for c in corrections: ak4chsCorr.append(os.path.join(dataPath,options.newJECPayloadNames+"_"+c+"_AK4PFchs.txt"))
  ak8chsCorr = []
  for c in corrections: ak8chsCorr.append(os.path.join(dataPath,options.newJECPayloadNames+"_"+c+"_AK8PFchs.txt"))
  print "ak4chsCorr ", ak4chsCorr
  print "ak8chsCorr ", ak8chsCorr
  process.ana.jetAK4selParams.newJECPayloadNames = cms.vstring(ak4chsCorr)
  process.ana.jetAK4BTaggedselParams.newJECPayloadNames = cms.vstring(ak4chsCorr)
  process.ana.jetAK8selParams.newJECPayloadNames = cms.vstring(ak8chsCorr)
'''

if options.sys:
    process.anabcUp = process.ana.clone(
                                        sys = cms.bool(True),
                                        btagsf_bcUp = cms.bool(True),
                                        )
    process.anabcDown = process.ana.clone(
                                          sys = cms.bool(True),
                                          btagsf_bcDown = cms.bool(True),
                                          )
    process.analightUp = process.ana.clone(
                                           sys = cms.bool(True),
                                           btagsf_lUp = cms.bool(True),
                                           )
    process.analightDown = process.ana.clone(
                                             sys = cms.bool(True),
                                             btagsf_lDown = cms.bool(True),
                                             )
    process.anaJecUp = process.ana.clone(
                                         sys = cms.bool(True),
                                         jecShift = cms.double(1),
                                         )
    #process.anaJecUp.jetAK4selParams.jecShift = cms.double(1.)
    #process.anaJecUp.jetAK8selParams.jecShift = cms.double(1.)
    #process.anaJecUp.jetAK4BTaggedselParams.jecShift = cms.double(1.) 
    process.anaJecDown = process.ana.clone(
                                           sys = cms.bool(True),
                                           jecShift = cms.double(-1),
                                           )
    #process.anaJecDown.jetAK4selParams.jecShift = cms.double(-1.)
    #process.anaJecDown.jetAK8selParams.jecShift = cms.double(-1.)
    #process.anaJecDown.jetAK4BTaggedselParams.jecShift = cms.double(-1.) 
    process.anaJerUp = process.ana.clone(
                                         sys = cms.bool(True),
                                         jerShift = cms.double(2),
                                         )
    #process.anaJerUp.jetAK4selParams.jerShift = cms.int32(2)
    #process.anaJerUp.jetAK8selParams.jerShift = cms.int32(2)
    #process.anaJerUp.jetAK4BTaggedselParams.jerShift = cms.int32(2) 
    process.anaJerDown = process.ana.clone(
                                           sys = cms.bool(True),
                                           jerShift = cms.double(-1),
                                           )
    #process.anaJerDown.jetAK4selParams.jerShift = cms.int32(-1)
    #process.anaJerDown.jetAK8selParams.jerShift = cms.int32(-1)
    #process.anaJerDown.jetAK4BTaggedselParams.jerShift = cms.int32(-1) 
    process.anaPileupUp = process.ana.clone(
                                            sys = cms.bool(True),
                                            PileupUp = cms.bool(True),
                                            )
    process.anaPileupDown = process.ana.clone(
                                              sys = cms.bool(True),
                                              PileupDown = cms.bool(True),
                                              )
    process.anaTopPtRewUp = process.ana.clone(
                                              sys = cms.bool(True),
                                              topPtRewUp = cms.bool(True),
                                              )
    process.anaTopPtRewDown = process.ana.clone(
                                              sys = cms.bool(True),
                                              topPtRewDown = cms.bool(True),
                                              )
    process.anaWjetHtUp = process.ana.clone(
                                              sys = cms.bool(True),
                                              wjetHtShift = cms.double(2),
                                              )
    process.anaWjetHtDown = process.ana.clone(
                                              sys = cms.bool(True),
                                              wjetHtShift = cms.double(-1),
                                               )
    

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(
                                                         options.outFileName
                                                         )
                                   )

outCommand = ['keep *', 'drop *_evtcleaner_*_*', 'drop *_photons_*_*', 'drop *_photonjets_*_*', 'drop *_*Puppi_*_*', 'drop *_TriggerResults_*_*']

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(options.outFileName.split('.',1)[0]+'_skim.root'),
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
                               dropMetaData = cms.untracked.string('DROPPED'),#'type_label_instance_process'
                               outputCommands = cms.untracked.vstring(outCommand )
                               )

## Event counters
from Analysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone(isData=options.isData)
process.cleanedEvents = eventCounter.clone(isData=options.isData)
process.finalEvents = eventCounter.clone(isData=options.isData)

if options.sys:
    process.p = cms.Path(
                         process.allEvents
                         *process.evtcleaner
                         *cms.ignore(process.ana)
                         *cms.ignore(process.anabcUp)
                         *cms.ignore(process.anabcDown)
                         *cms.ignore(process.analightUp)
                         *cms.ignore(process.analightDown)
                         *cms.ignore(process.anaJecUp)
                         *cms.ignore(process.anaJecDown)
                         *cms.ignore(process.anaJerUp)
                         *cms.ignore(process.anaJerDown)
                         *cms.ignore(process.anaPileupUp)
                         *cms.ignore(process.anaPileupDown)
                         *cms.ignore(process.anaTopPtRewUp)
                         *cms.ignore(process.anaTopPtRewDown)
                         *cms.ignore(process.anaWjetHtUp)
                         *cms.ignore(process.anaWjetHtDown)
                         )
elif not options.sys and not options.isData:
    process.p = cms.Path(
                         process.allEvents
                         *process.evtcleaner
                         #*process.cleanedEvents
                         *cms.ignore(process.ana)
                         #*cms.ignore(process.anaH)
                         #* process.finalEvents
                         )
else:
    process.p = cms.Path(
                         process.allEvents
                         *process.evtcleaner
                         *cms.ignore(process.ana)
                         )

if options.doSkim:
    process.outpath = cms.EndPath(process.out)

