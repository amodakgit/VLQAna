Main Analyzer:
Analysis/VLQAna/plugins/SingleLepAna.cc

Configuration File:
Analysis/VLQAna/test/SingleLepAna_cfg.py

Job Configuration: which processes to run
Analysis/VLQAna/test/allJobList.py

Crab Config:
Analysis/VLQAna/test/crab_dummy_os2lana.py

Running in Crab Instruction:
run.py --isData 1 (Data)
run.py --isData 0 --applyTopPtCorr 1 (for TTJets)
run.py --isData 0 (all other mc)
