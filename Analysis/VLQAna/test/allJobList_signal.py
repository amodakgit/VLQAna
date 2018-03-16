#! /bin/python

def job_list(isData):
    if isData:
        jobList = [
                   ['/SingleMuon/algomez-Run2016B-PromptReco-v2_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 'SingleMuon_B', '15', ''],
		   ['/SingleMuon/algomez-Run2016C-PromptReco-v2_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 'SingleMuon_C', '15', ''],
		   ['/SingleMuon/algomez-Run2016D-PromptReco-v2_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 'SingleMuon_D', '15', ''],
                   ['/SingleMuon/algomez-Run2016E-PromptReco-v2_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 'SingleMuon_E', '15', ''],
                   ['/SingleMuon/algomez-Run2016F-PromptReco-v1_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 'SingleMuon_F', '15', ''],
		   ['/SingleMuon/algomez-Run2016G-PromptReco-v1_B2GAnaFW_80X_V2p1-3e507461fd667ac6961fa4af5b123b09/USER', 'SingleMuon_G', '15', ''],
		   ['/SingleMuon/algomez-Run2016H-PromptReco-v2_B2GAnaFW_80X_V2p01p1-3e507461fd667ac6961fa4af5b123b09/USER', 'SingleMuon_H', '15', ''],
		   ['/SingleElectron/skhi-Run2016B-PromptReco-v2_B2GAnaFW_80X_V2p1-0f795150a28f873be7cbcf6c08b639fa/USER', 'SingleElectron_B', '15', ''],
		   ['/SingleElectron/skhi-Run2016C-PromptReco-v2_B2GAnaFW_80X_V2p1-0f795150a28f873be7cbcf6c08b639fa/USER', 'SingleElectron_C', '15', ''],
		   ['/SingleElectron/skhi-Run2016D-PromptReco-v2_B2GAnaFW_80X_V2p1-0f795150a28f873be7cbcf6c08b639fa/USER', 'SingleElectron_D', '15', ''],
		   ['/SingleElectron/skhi-Run2016E-PromptReco-v2_B2GAnaFW_80X_V2p1-0f795150a28f873be7cbcf6c08b639fa/USER', 'SingleElectron_E', '15', ''],
	           ['/SingleElectron/skhi-Run2016F-PromptReco-v1_B2GAnaFW_80X_V2p1-0f795150a28f873be7cbcf6c08b639fa/USER', 'SingleElectron_F', '15', ''],
	           ['/SingleElectron/skhi-Run2016G-PromptReco-v1_B2GAnaFW_80X_V2p1-0f795150a28f873be7cbcf6c08b639fa/USER', 'SingleElectron_G', '15', ''],
		   ['/SingleElectron/skhi-Run2016H-PromptReco-v2_B2GAnaFW_80X_V2p1-0f795150a28f873be7cbcf6c08b639fa/USER', 'SingleElectron_H', '15', ''], 
                   ]

        return jobList
    else:
        jobList = [
 		   ['/TprimeBToBW_M-700_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y700', '1',''],
                   ['/TprimeBToBW_M-800_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y800', '1',''],
	           ['/TprimeBToBW_M-900_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y900', '1',''],
                   ['/TprimeBToBW_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y1000', '1',''],
		   ['/TprimeBToBW_M-1200_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y1200', '1',''],
                   ['/TprimeBToBW_M-1300_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y1300', '1',''],
		   ['/TprimeBToBW_M-1400_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y1400', '1',''],
		   ['/TprimeBToBW_M-1500_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y1500', '1',''],
		   ['/TprimeBToBW_M-1600_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y1600', '1',''],
		   ['/TprimeBToBW_M-1700_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y700', '1',''],
		   ['/TprimeBToBW_M-1800_TuneCUETP8M1_13TeV-madgraph-pythia8/skhi-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'Y1800', '1',''],
                   ['/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/skhi-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'ST_t_top', '1',''],		  
                   ['/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/skhi-RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_B2GAnaFW_80X_V2p1-1fa4321daf986948ecc338ec9de72bc7/USER', 'ST_t_antitop', '1',''],
		   ['/TT_TuneCUETP8M1_13TeV-powheg-pythia8/grauco-B2GAnaFW_80X_V2p1-edbed0685401a5848e7d61871b3a63d8/USER', 'TTBar_Powheg', '20',''],
                   ['/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'TTBar_Powheg_NewTune', '20',''],
                   ['/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/grauco-B2GAnaFW_80X_V2p1-edbed0685401a5848e7d61871b3a63d8/USER', 'WJets_HT100_200', '5',''],
                   ['/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/grauco-B2GAnaFW_80X_V2p1-edbed0685401a5848e7d61871b3a63d8/USER', 'WJets_HT200_400', '5',''],
                   ['/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/grauco-B2GAnaFW_80X_V2p1-edbed0685401a5848e7d61871b3a63d8/USER', 'WJets_HT400_600', '5',''],
                   ['/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/grauco-B2GAnaFW_80X_V2p1-edbed0685401a5848e7d61871b3a63d8/USER', 'WJets_HT600_800', '5',''],
                   ['/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/grauco-B2GAnaFW_80X_V2p1-edbed0685401a5848e7d61871b3a63d8/USER', 'WJets_HT800_1200', '5',''],
                   ['/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/grauco-B2GAnaFW_80X_V2p1-edbed0685401a5848e7d61871b3a63d8/USER', 'WJets_HT1200_2500', '5',''],
                   ['/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/grauco-B2GAnaFW_80X_V2p1-edbed0685401a5848e7d61871b3a63d8/USER', 'WJets_HT2500_inf', '5',''],
                   ['/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-a300e501c1a543433113bdc094d47173/USER', 'DY_PT100_250', '5',''],
                   ['/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-a300e501c1a543433113bdc094d47173/USER', 'DY_PT250_400', '5',''],
                   ['/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-a300e501c1a543433113bdc094d47173/USER', 'DY_PT400_650', '5',''],
                   ['/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-a300e501c1a543433113bdc094d47173/USER', 'DY_PT650_inf', '5',''],
                   ['/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'ST_tw_top', '1',''],
                   ['/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/vorobiev-B2GAnaFW_RunIISpring16MiniAODv2_25ns_v80x_v2p1-24c95768e44153a28c3920000b3803cb/USER', 'ST_tw_antitop', '1',''],

		  ]
        return jobList



                   
                   
                   
                   
                   
                   
                   
                   
                   
