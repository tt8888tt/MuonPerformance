import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('SliceTestEfficiencyAnalysis',eras.Run2_2018,eras.run3_GEM)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '') #for rd
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v11', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.skipEvents = cms.untracked.uint32(0)

from glob import glob
process.source.fileNames.extend(
     ['root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/jlee/SingleMuon/Run2018C-v1/FULLRECOv1/AOD_605.root']
     #['root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/jlee/SingleMuon/Run2018C-v1/FULLRECOv1/AOD_%d.root'%i for i in range(600,605)]
   # ['file:'+f for f in glob('/xrootd/store/user/jlee/SingleMuon/Run2018C-v1/RECOv1/step3*.root')][:]
)

process.options = cms.untracked.PSet()
process.TFileService = cms.Service("TFileService",fileName = cms.string("histo2018.root"))

process.SliceTestEfficiencyAnalysis = cms.EDAnalyzer('SliceTestEfficiencyAnalysis',
    process.MuonServiceProxy,
    gemRecHits = cms.InputTag("gemRecHits"),
    muons = cms.InputTag("muons"),
    #muons = cms.InputTag("globalMuons"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    amcData = cms.InputTag("muonGEMDigis", "AMCdata", "reRECO"),
    gebStatusCol = cms.InputTag("muonGEMDigis", "gebStatus", "reRECO"),
    vfatStatusCol = cms.InputTag("muonGEMDigis", "vfatStatus", "reRECO"),

    minPt = cms.double(20),
    fiducialCut = cms.bool(True),
    amcBxCut = cms.bool(True),
)
process.p = cms.Path(process.SliceTestEfficiencyAnalysis)
