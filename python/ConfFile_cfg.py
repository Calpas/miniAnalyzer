import FWCore.ParameterSet.Config as cms

process = cms.Process("h2taus")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

# MVA electron ID
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data forma
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
# add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)


runOnData=False

#########
#       #
# pfmet #
#       #
#########
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process, isData=runOnData )


###########
#         #
# mva met #
#         #
###########
# recorrect the jets
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('patJetsReapplyJEC'),
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')  # Do not forget 'L2L3Residual' on data!
)

# configure MVAMet
from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET
runMVAMET( process, jetCollectionPF = "patJetsReapplyJEC")
process.MVAMET.srcLeptons  = cms.VInputTag("slimmedMuons", "slimmedElectrons", "slimmedTaus")
process.MVAMET.requireOS   = cms.bool(False)

# For 8_0_20, MVA MET is pointing at the regular pfCandidates ...
process.tauMET.srcPFCands =  cms.InputTag("packedPFCandidates")

####################
#		   #
# tau gen particle #
#		   #
####################
process.analyzeTauSequence = cms.Sequence()
process.load('PhysicsTools/JetMCAlgos/TauGenJets_cfi')
process.tauGenJets = cms.EDProducer(
        "TauGenJetProducer",
        GenParticles =  cms.InputTag('prunedGenParticles'),
        includeNeutrinos = cms.bool( False ),
        verbose = cms.untracked.bool( False )
)
process.analyzeTauSequence += process.tauGenJets

process.load('PhysicsTools/JetMCAlgos/TauGenJetsDecayModeSelectorAllHadrons_cfi')
process.analyzeTauSequence += process.tauGenJetsSelectorAllHadrons

process.tauGenJetsSelectorElectronAndMuon = cms.EDFilter("TauGenJetDecayModeSelector",
	src    = cms.InputTag("tauGenJets"),
	select = cms.vstring('electron', 'muon'),
        filter = cms.bool(False),
)
process.analyzeTauSequence += process.tauGenJetsSelectorElectronAndMuon


##########################
#			 #
# configure miniAnalyzer #
#			 #
##########################
process.h2taus      = cms.EDAnalyzer('miniAnalyzer',
    genpruned       = cms.InputTag("prunedGenParticles"),
    gentaul         = cms.InputTag("tauGenJetsSelectorElectronAndMuon"),
    gentauh         = cms.InputTag("tauGenJetsSelectorAllHadrons"),
    rho             = cms.InputTag("fixedGridRhoFastjetAll"),
    genruninfo      = cms.InputTag("generator"),
    vertices        = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons           = cms.InputTag("slimmedMuons"),
    electrons       = cms.InputTag("slimmedElectrons"),
    taus            = cms.InputTag("slimmedTaus"),
    jets            = cms.InputTag("slimmedJets"),
    pu              = cms.InputTag("slimmedAddPileupInfo"),
    eleMediumIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
    eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
    mvaValuesMap    = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
    mets            = cms.InputTag("slimmedMETs", "", "h2taus"), 
    metspuppi       = cms.InputTag("slimmedMETsPuppi"),
    #bits            = cms.InputTag("TriggerResults","","HLT"),
    bits            = cms.InputTag("TriggerResults::HLT2"),
    objects         = cms.InputTag("selectedPatTrigger"),
    prescales       = cms.InputTag("patTrigger"),
    mutauFilterName = cms.vstring(
    "HLT_IsoMu18_v3",
    "HLT_IsoMu20_v4",
    "HLT_IsoMu22_v3",
    "HLT_IsoMu22_eta2p1_v2",
    "HLT_IsoMu24_v2",
    "HLT_IsoMu27_v4",
    "HLT_IsoTkMu18_v3",
    "HLT_IsoTkMu20_v5",
    "HLT_IsoTkMu22_eta2p1_v2",
    "HLT_IsoTkMu22_v3",
    "HLT_IsoTkMu24_v2",
    "HLT_IsoTkMu27_v4",
    "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v5",
    "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v5",
    "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v2",
    "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v2",
    "HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v2",
    ),
    etauFilterName = cms.vstring(
    "HLT_Ele23_WPLoose_Gsf_v4",
    "HLT_Ele24_eta2p1_WPLoose_Gsf_v2",
    "HLT_Ele25_WPTight_Gsf_v2",
    "HLT_Ele25_eta2p1_WPLoose_Gsf_v2",
    "HLT_Ele25_eta2p1_WPTight_Gsf_v2",
    "HLT_Ele27_WPLoose_Gsf_v2",
    "HLT_Ele27_WPTight_Gsf_v2",
    "HLT_Ele27_eta2p1_WPLoose_Gsf_v3",
    "HLT_Ele27_eta2p1_WPTight_Gsf_v3",
    "HLT_Ele32_eta2p1_WPTight_Gsf_v3",
    "HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v3",
    "HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v2",
    "HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v2",
    "HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v2",
    "HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v2",
    ),
    tautauFilterName = cms.vstring(
    "HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg_v2",
    "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v3",
    "HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v5",
    ),
    emuFilterName = cms.vstring(
    "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v4",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4",
    "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4",
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4",
    "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v2",
    "HLT_IsoMu18_v3",
    "HLT_IsoMu20_v4",
    "HLT_IsoMu22_v3",
    "HLT_IsoMu22_eta2p1_v2",
    "HLT_IsoMu24_v2",
    "HLT_IsoMu27_v4",
    "HLT_IsoTkMu18_v3",
    "HLT_IsoTkMu20_v5",
    "HLT_IsoTkMu22_eta2p1_v2",
    "HLT_IsoTkMu22_v3",
    "HLT_IsoTkMu24_v2",
    "HLT_IsoTkMu27_v4",
    "HLT_Ele23_WPLoose_Gsf_v4",
    "HLT_Ele24_eta2p1_WPLoose_Gsf_v2",
    "HLT_Ele25_WPTight_Gsf_v2",
    "HLT_Ele25_eta2p1_WPLoose_Gsf_v2",
    "HLT_Ele25_eta2p1_WPTight_Gsf_v2",
    "HLT_Ele27_WPLoose_Gsf_v2",
    "HLT_Ele27_WPTight_Gsf_v2",
    "HLT_Ele27_eta2p1_WPLoose_Gsf_v3",
    "HLT_Ele27_eta2p1_WPTight_Gsf_v3",
    "HLT_Ele32_eta2p1_WPTight_Gsf_v3",
    ),
 )


process.p = cms.Path( 
     process.fullPatMetSequence *
     process.analyzeTauSequence *
     process.egmGsfElectronIDSequence * 
     process.h2taus)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.TFileService = cms.Service("TFileService", fileName = cms.string('output.root'))

process.source = cms.Source("PoolSource",

    #eventsToProcess = cms.untracked.VEventRange('1:2786-1:2786'), # run:evt

    fileNames = cms.untracked.vstring(
      # DY
      #'file:/lstore/cms/calpas/h2taus/CMSSW_8_0_16/src/miniAOD/miniAnalyzer/miniAOD/DYJetsToLL_M-50_13TeV_PUSpring16RAWAODSIM_reHLT_80X_mcRun2.root',
      # SUSY160
      #'file:/lstore/cms/calpas/h2taus/CMSSW_8_0_16/src/miniAOD/miniAnalyzer/miniAOD/SUSYGluGluToHToTauTau_M-160_13TeV_RunIISpring16MiniAODv2_0.root',
      #'file:/lstore/cms/calpas/h2taus/CMSSW_8_0_16/src/miniAOD/miniAnalyzer/miniAOD/SUSYGluGluToHToTauTau_M-160_13TeV_RunIISpring16MiniAODv2_1.root',
      #'file:/lstore/cms/calpas/h2taus/CMSSW_8_0_16/src/miniAOD/miniAnalyzer/miniAOD/SUSYGluGluToHToTauTau_M-160_13TeV_RunIISpring16MiniAODv2_2.root',
      #'file:/lstore/cms/calpas/h2taus/CMSSW_8_0_16/src/miniAOD/miniAnalyzer/miniAOD/SUSYGluGluToHToTauTau_M-160_13TeV_RunIISpring16MiniAODv2_3.root',
      #'file:/lstore/cms/calpas/h2taus/CMSSW_8_0_16/src/miniAOD/miniAnalyzer/miniAOD/SUSYGluGluToHToTauTau_M-160_13TeV_RunIISpring16MiniAODv2_4.root',
      #'file:/lstore/cms/calpas/h2taus/CMSSW_8_0_16/src/miniAOD/miniAnalyzer/miniAOD/SUSYGluGluToHToTauTau_M-160_13TeV_RunIISpring16MiniAODv2_5.root',

#QCD
'/store/mc/RunIISpring16MiniAODv1/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/000E0956-5B0E-E611-8384-02163E011809.root'

# VBF125
#'file:/lstore/cms/calpas/h2taus/CMSSW_8_0_20/src/miniAOD/miniAnalyzer/miniAOD/VBFHToTauTau_M125_13TeV_MiniAODv2.root'
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/5066EAF6-4F3A-E611-8424-0025909077C6.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/0863B733-1A39-E611-AF47-0025905C53D8.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/0A6F49D3-2A39-E611-9788-0025905C96E8.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/0CEF531A-2A39-E611-BA2B-0025905D1D02.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/125513C8-3C39-E611-B006-0025905C54FE.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/223A3808-3339-E611-9128-0025905C96A4.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/242670BF-3C39-E611-A015-0025904C6622.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/32643904-3339-E611-A047-0025905C53B2.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/468A1AC1-3C39-E611-870D-0025905C3D3E.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/4E4E8192-2E39-E611-ACD3-0025905C54FE.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/5260BAEE-5539-E611-BE18-0025905C53F0.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/5668CF4C-4A39-E611-BBCD-00259090766E.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/56B2D2F5-5D39-E611-AF5F-001E6724862E.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/5E5B7EB5-0F39-E611-841B-0025904C66EC.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/6221201C-6039-E611-8384-001E6724862E.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/6ABB0A93-2539-E611-A0F5-0025904C66A6.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/74CBD530-6139-E611-88F9-001E67247C60.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/7AAAE021-3F39-E611-9115-0025905C2CE4.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/7AB54108-3339-E611-9242-0025905C5432.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/7E30AC57-4039-E611-9B3F-0025905C2CE4.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/8AB47F72-3539-E611-B67B-0025904C66EE.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/8E746820-3F39-E611-A632-0025905C96A4.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/96C1F39B-2E39-E611-8142-0025905C3D40.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/9CDFBBFB-1939-E611-8771-0025904C4F9C.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/A43A3EBF-3E39-E611-82EC-0025905C53AA.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/B85BE9A4-4139-E611-BB99-0025904E9010.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/C0C4CE9D-2E39-E611-943F-0025905C2CA4.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/C240773E-2639-E611-9C5C-0025905C96E8.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/E8384CB2-0F39-E611-BEF1-0025904C51F0.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/F4670E21-3F39-E611-9EA5-0025905C96EA.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/F6C63EBA-4139-E611-AE1C-0025905C43EA.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/FC331174-2E39-E611-87A6-0025905C96EA.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/80000/FC78A273-3539-E611-B45B-0025905D1CB2.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/2860B442-6D3A-E611-821C-00259090765E.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/5433B0FC-6C3A-E611-8FCA-002590907826.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/7C29D5C7-6D3A-E611-89F7-0025909082D2.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/AA81FA08-6D3A-E611-9187-002590909086.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/AE55BF2A-6D3A-E611-8CBC-00259090765E.root",
      #"/store/mc/RunIISpring16MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/F4E224EC-6D3A-E611-991D-00259090829E.root",
      ),
)



#from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
#process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
#  src     = cms.InputTag("patJetsReapplyJEC"),
#  levels  = ['L1FastJet', 'L2Relative', 'L3Absolute'],
#  payload = 'AK4PFchs' 
# ) 

#from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
#process.patJetsReapplyJEC = updatedPatJets.clone(
#  jetSource 	       = cms.InputTag("patJetsReapplyJEC"),
#  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
#  )

#if(isData):
#  process.patJetCorrFactorsReapplyJEC.levels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
