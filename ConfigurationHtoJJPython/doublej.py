import FWCore.ParameterSet.Config as cms

etafilter = cms.EDFilter("PythiaFilter",
	MaxEta = cms.untracked.double(9999.0),
	MinEta = cms.untracked.double(-9999.0),
	ParticleID = cms.untracked.int32(35)
)

from Configuration.Generator.PythiaUESettings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
	pythiaHepMCVerbosity = cms.untracked.bool(False),
	maxEventsToPrint = cms.untracked.int32(1),
	pythiaPylistVerbosity = cms.untracked.int32(1),
	displayPythiaCards = cms.untracked.bool(False),
	comEnergy = cms.double(13000.0),
	PythiaParameters = cms.PSet(
		pythiaUESettingsBlock,
		pythiaEtab = cms.vstring(
			'Higgs:useBSM = on',
			'HiggsBSM:gg2H2 = on',
			#'HiggsBSM:ffbar2H2 = on',
			'HiggsH2:coup2d = 10.0',
			'HiggsH2:coup2u = 10.0',
			'HiggsH2:coup2Z = 0.0',
			'HiggsH2:coup2W = 0.0',
			'HiggsA3:coup2H2Z = 0.0',
			'HiggsH2:coup2A3A3 = 0.0',
			'HiggsH2:coup2H1H1 = 0.0',
			'443:onMode = off',
			'443:onIfMatch 13 13',
			'333:onMode = off',
			'333:onIfMatch 13 13',
			#'223:onMode = off',
			#'223:onIfMatch 13 13',
			'553:onMode = off',
			'553:onIfMatch 13 13',
			############### For Floating Mass Distribution#######
			#'35:mMin = 15', 
			#'35:mMax = 25',
			#'35:m0   = 23.  ! Higgs mass',
			#'35:mWidth = 10 !Higgs Width',
			############# For Fixed Mass Distribution#############
			'35:mMin = 0',
			'35:mMax = 200',
			'35:m0   = 125.0',
			'35:mWidth = 0.00',
			'35:addChannel 1 1.00 100 443 443',
			# '35:addChannel 1 1.00 100 13 -13 553',
			# '35:addChannel 1 1.00 100 443 333',
			# '35:addChannel 1 1.00 100 443 223',
			# '35:addChannel 1 1.00 100 443 553',
			# '35:addChannel 1 1.00 100 333 553',
			# '35:addChannel 1 1.00 100 223 553',
			'35:onMode = off',
			# '35:onIfMatch 443 333'), ##Jpsi Phi decay channel!
			# '35:onIfMatch 13 -13 553'), ## Y(1S) mumu
			'35:onIfMatch 443 443'), ##Jpsi Jpsi decay channel!
			# '35:onIfMatch 443 223'), ##Jpsi Omega decay channel!
			# '35:onIfMatch 443 553'), ##Jpsi Upsilon decay channel!
			# '35:onIfMatch 333 553'), ##Phi Upsilon decay channel!
			# '35:onIfMatch 223 553'), ##Omega Upsilon decay channel!
			# This is a vector of ParameterSet names to be read, in this order
		parameterSets = cms.vstring(
			'pythiaEtab'
		)
	)
)


ProductionFilterSequence = cms.Sequence(generator+etafilter)
