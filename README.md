2017 Recipie for MC Production

Step zero production
Help from Eliza

 cmsDriver.py Configuration/GenProduction/python/ThirteenTeV/doublej.py --fileout file:RSGravitonToZZ_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8_GEN-SIM.root -s GEN,SIM --mc --datatier GEN-SIM --beamspot Realistic25ns13TeVEarly2017Collision --conditions auto:phase1_2017_realistic --eventcontent RAWSIM --era Run2_2017 --python_filename RSGravitonToZZ_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8_GEN-SIM_cfg.py --no_exec -n 50


doublej.py file is in location
Fermilab:--->/uscms_data/d3/hacharya/MC13TeV_Resonance/MC_FromTexasA/ElizaHelp/CMSSW_9_4_0/src/Configuration/GenProduction/python/ThirteenTeV/
Step 0 work is done in Fermi lab---> /uscms_data/d3/hacharya/StefanHimalMCGeneration/

Step1 &2  should be done in cern...
Cern working directory is /afs/cern.ch/work/h/hacharya/ParallelJobCERN/
web:
https://twiki.cern.ch/twiki/bin/view/CMS/BPHMonteCarloContactInfo#Preparation_of_Requests_For_user

Step 0
>>cmsrel CMSSW_9_3_1
>>cd  CMSSW_9_3_1/src
>>copy configuration….whole path and cfi file inside python (or more inside)

>>scram b
>>cmsenv
>>voms-proxy-init --voms cms

>>cmsDriver.py Configuration/GenProduction/python/ThirteenTeV/doublej.py --fileout file:output_step0.root --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 93X_mc2017_realistic_v3 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN,SIM --nThreads 1 --geometry DB:Extended --era Run2_2017 --python_filename outputFileName_1_cfg.py  --customise Configuration/DataProcessing/Utils.addMonitoring -n 100


Step 1
2 steps to DIGI-RECO production under 9_4_X
>>cmsDriver.py step1 --fileout file:YourRootFile_step1.root --filein file:output_step0.root --pileup_input "dbs:/Neutrino_E-10_gun/RunIISummer17PrePremix-MCv2_correctPU_94X_mc2017_realistic_v9-v1/GEN-SIM-DIGI-RAW" --mc --eventcontent PREMIXRAW --datatier GEN-SIM-RAW --conditions 94X_mc2017_realistic_v10 --step DIGIPREMIX_S2,DATAMIX,L1,DIGI2RAW,HLT:2e34v40 --nThreads 8 --datamix PreMix --era Run2_2017 --python_filename yourConfigFile_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n -1
Step 2
>>cmsDriver.py step2 --filein file:YourRootFile_step1.root --fileout file:YourRootFile_step2.root --mc --eventcontent AODSIM runUnscheduled --datatier AODSIM --conditions 94X_mc2017_realistic_v10 --step RAW2DIGI,RECO,RECOSIM,EI --nThreads 8 --era Run2_2017 --python_filename YourConfigFile_2_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 1751 || exit $? ;


OR BY “Ozlem Ozcelik Ozludil” suggestion “STANDARD”
Step 0:
>>cmsDriver.py Configuration/GenProduction/python/ThirteenTeV/doublej.py --fileout file:output_step0.root --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 93X_mc2017_realistic_v3 --beamspot Realistic50ns13TeVCollision --step GEN,SIM --nThreads 1 --geometry DB:Extended --era Run2_2017 --python_filename outputFileName_1_cfg.py  --customise Configuration/DataProcessing/Utils.addMonitoring -n 100


In CMSSW_9_4_0/src
Step 1ls

cmsDriver.py --filein file:output_step0.root --fileout file:OutputFile_step1.root --pileup_input "dbs:/Neutrino_E-10_gun/RunIISummer17PrePremix-MCv2_correctPU_94X_mc2017_realistic_v9-v1/GEN-SIM-DIGI-RAW" --mc --eventcontent PREMIXRAW --datatier GEN-SIM-RAW --conditions 94X_mc2017_realistic_v10 --step DIGIPREMIX_S2,DATAMIX,L1,DIGI2RAW,HLT:2e34v40 --nThreads 8 --datamix PreMix --era Run2_2017 --python_filename yourConfigFile_cfg.py --customise Configuration/DataProcessing/Utils.addMonitoring  -n -1

If you do not want to execute cmsrun immediately  add --no_exec


Step 2
cmsDriver.py --filein file:OutputFile_step1.root --fileout file:OutputFile_step2.root --mc --eventcontent AODSIM runUnscheduled --datatier AODSIM --conditions 94X_mc2017_realistic_v10 --step RAW2DIGI,RECO,RECOSIM,EI --nThreads 8 --era Run2_2017 --python_filename YourConfigFile_2_cfg.py  --customise Configuration/DataProcessing/Utils.addMonitoring -n -1 



Now process from “FourMuAnalysis Program” after step 2




2016 Recipie for the MC Production

2016 working without pileup
In CMSSW_9_3_1
cmsenv
scram b

Step 0
cmsDriver.py Configuration/GenProduction/python/ThirteenTeV/doublej.py --fileout file:RSGravitonToZZ_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8_GEN-SIM.root -s GEN,SIM --mc --datatier GEN-SIM --beamspot Realistic50ns13TeVCollision --conditions auto:run2_mc --eventcontent RAWSIM --era Run2_2016 --python_filename RSGravitonToZZ_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8_GEN-SIM_cfg.py --no_exec -n 50


Step 1
cmsDriver.py -s DIGI,L1,DIGI2RAW,HLT:@relval2016 --datatier GEN-SIM-RAW --conditions auto:run2_mc --eventcontent RAWSIM --era Run2_2016 --filein file:RSGravitonToZZ_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8_GEN-SIM.root --fileout file:RSGravitonToZZ_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8_GEN-SIM-RAW.root --python_filename RSGravitonToZZ_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8_FULLSIM_GEN_SIM_RAW_cfg.py -n -1 --no_exec



Step 2

cmsDriver.py -s RAW2DIGI,L1Reco,RECO --datatier RECO --conditions auto:run2_mc --eventcontent AODSIM --era Run2_2016 --filein file:RSGravitonToZZ_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8_GEN-SIM-RAW.root --fileout file:RSGravitonToZZ_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8_GEN-SIM-RAW-RECO.root --python_filename RSGravitonToZZ_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8_FULLSIM_GEN_SIM_RAW_RECO_cfg.py -n -1 --no_exec

Does not  work  with FourMuon program to reconstruct 


Final 2016 MC recipe from Jordan Martins
Email : jordan.martins@cern.ch

Step 0 is done in fermi

Location : /uscms_data/d3/hacharya/StefanHimalMCGeneration/2016MC_WorkHelpFromJordan

>>scram p CMSSW CMSSW_7_1_42_patch1
>>cd src
>>eval `scram runtime -sh`  for lxplus OR  eval `scram runtime -csh` for fermi cluster
>>mkdir -p Configuration/GenProduction/python/
--->Our fragment is doublej.py
--->copy your fragment to Configuration/GenProduction/python/ 
---->go to CMSSW_7_1_42_patch1/src
>>cmsenv
>>scram b


Now do the step 0 step
>>cmsDriver.py Configuration/GenProduction/python/doublej_cfi.py --fileout file:GS.root --mc --eventcontent RAWSIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --conditions MCRUN2_71_V1::All --beamspot Realistic50ns13TeVCollision --step GEN,SIM --magField 38T_PostLS1 --python_filename GS_cfg.py --no_exec -n 100

Step 1

>>scram p CMSSW CMSSW_8_0_31
>>eval `scram runtime -sh`
>>cmsDriver.py step1 --filein file:GS.root --fileout file:step1.root  --pileup_input "dbs:/Neutrino_E-10_gun/RunIISpring15PrePremix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v2-v2/GEN-SIM-DIGI-RAW" --mc --eventcontent PREMIXRAW --datatier GEN-SIM-RAW --conditions 80X_mcRun2_asymptotic_2016_TrancheIV_v6 --step DIGIPREMIX_S2,DATAMIX,L1,DIGI2RAW,HLT:@frozen2016 --nThreads 8 --datamix PreMix --era Run2_2016 --python_filename step1_1_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n -1

Step 2

>>cmsDriver.py step2 --filein file:step1.root --fileout file:AOD.root --mc --eventcontent AODSIM --runUnscheduled --datatier AODSIM --conditions 80X_mcRun2_asymptotic_2016_TrancheIV_v6 --step RAW2DIGI,RECO,EI --nThreads 8 --era Run2_2016 --python_filename step2_2_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n -1

Now do Reconstruction as in data AOD

Works for four muon reconstruction in CMSSW_8_0_31 in cern to get final root files


2017 MC without Pileup

By Jordan:

export SCRAM_ARCH=slc6_amd64_gcc630
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram p CMSSW CMSSW_9_3_15_patch3
cd CMSSW_9_3_15_patch3/src
eval `scram runtime -sh`

- GS:

mkdir -p Configuration/GenProduction/python/

move your fragment to Configuration/GenProduction/python/

cmsDriver.py Configuration/GenProduction/python/GS_fragment.py --fileout file:GSout.root --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 93X_mc2017_realistic_v3 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN,SIM --geometry DB:Extended --era Run2_2017 --python_filename GS_fragment_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32(178571)" -n <number_of_evts>

cmsRun GS_fragment_cfg.py

- Step1 of DRNoPU:

cmsDriver.py step1 --filein GSout.root --fileout file:DR_step1.root --mc --eventcontent RAWSIM --pileup NoPileUp --datatier GEN-SIM-RAW --conditions 94X_mc2017_realistic_v10 --step DIGI,L1,DIGI2RAW,HLT:2e34v40 --nThreads 8 --geometry DB:Extended --era Run2_2017 --python_filename DR_1_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n <number_of_evts>

cmsRun DR_1_cfg.py

Step2 of DRNoPU:

cmsDriver.py step2 --filein file:DR_step1.root --fileout file:AOD.root --mc --eventcontent AODSIM --datatier AODSIM --conditions 94X_mc2017_realistic_v10 --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI --nThreads 8 --geometry DB:Extended --era Run2_2017 --python_filename DR_2_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n <number_of_evts>

cmsRun DR_2_cfg.py

(Also works in CMSSW_9_4_0 after step 0...)
