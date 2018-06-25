from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")

config.General.requestName = DUMMY_NAME
config.General.workArea = DUMMY_WORKDIR
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'SingleLepAna_cfg.py'
config.JobType.pyCfgParams = [DATa,  PUOFF, FILTERSIGNAL, WJetHTCOR, SIGNALTYPE, TopPTCOR]
config.JobType.inputFiles = [
'Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt',
'Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt',
'Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt',
'Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt',
'Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt',
'Summer16_23Sep2016V4_MC_L2L3Residual_AK8PFchs.txt',
'Summer16_23Sep2016V4_MC_L2L3Residual_AK4PFchs.txt',
'Summer16_23Sep2016V4_MC_L1FastJet_AK8PFchs.txt',
'Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt',
'Summer16_23Sep2016HV4_DATA_Uncertainty_AK8PFchs.txt',
'Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFchs.txt',
'Summer16_23Sep2016HV4_DATA_L3Absolute_AK8PFchs.txt',
'Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt',
'Summer16_23Sep2016HV4_DATA_L2Relative_AK8PFchs.txt',
'Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt',
'Summer16_23Sep2016HV4_DATA_L2L3Residual_AK8PFchs.txt',
'Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt',
'Summer16_23Sep2016HV4_DATA_L1FastJet_AK8PFchs.txt',
'Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt',
'Summer16_23Sep2016GV4_DATA_Uncertainty_AK8PFchs.txt',
'Summer16_23Sep2016GV4_DATA_Uncertainty_AK4PFchs.txt',
'Summer16_23Sep2016GV4_DATA_L3Absolute_AK8PFchs.txt',
'Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt',
'Summer16_23Sep2016GV4_DATA_L2Relative_AK8PFchs.txt',
'Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt',
'Summer16_23Sep2016GV4_DATA_L2L3Residual_AK8PFchs.txt',
'Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt',
'Summer16_23Sep2016GV4_DATA_L1FastJet_AK8PFchs.txt',
'Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt',
'Summer16_23Sep2016EFV4_DATA_Uncertainty_AK8PFchs.txt',
'Summer16_23Sep2016EFV4_DATA_Uncertainty_AK4PFchs.txt',
'Summer16_23Sep2016EFV4_DATA_L3Absolute_AK8PFchs.txt',
'Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt',
'Summer16_23Sep2016EFV4_DATA_L2Relative_AK8PFchs.txt',
'Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt',
'Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK8PFchs.txt',
'Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt',
'Summer16_23Sep2016EFV4_DATA_L1FastJet_AK8PFchs.txt',
'Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt',
'Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK8PFchs.txt',
'Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt',
'Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFchs.txt',
'Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt',
'Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt',
'Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt',
'Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt',
'Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt',
'Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK8PFchs.txt',
'Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt',
'RunII2016Rereco_25ns_PUXsec72450nb.root',
'RunII2016Rereco_25ns_PUXsec69000nb.root',
'RunII2016Rereco_25ns_PUXsec65550nb.root',
'PUDistMC_Summer2016_25ns_Moriond17MC_PoissonOOTPU.root',
'CSVv2_Moriond17_B_H.csv',
'egammaRECOEffi.txt_EGM2D.root',
'egammaIDEffi.txt_EGM2D.root',
'SF_HLT_Ele32_eta2p1.root',
'ttjets_bTagEff.root',
                             'btag-eff-subjet.root', 'CSVv2_ichep.csv']

config.section_("Data")
config.Data.inputDataset = DUMMY_DATASET
config.Data.inputDBS = 'phys03'
config.Data.splitting = DUMMY_BASE
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = DUMMY_NUMBER
config.Data.ignoreLocality = False
config.Data.publication = False
config.Data.outLFNDirBase = DUMMY_OUTPUT_PATH

config.section_("Site")
config.Site.storageSite = DUMMY_SITE
config.section_('User')
