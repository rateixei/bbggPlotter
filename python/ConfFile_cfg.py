import FWCore.ParameterSet.Config as cms
from mAOD_GravFiles import *
from mAOD_RadFiles import *

process = cms.Process("bbggAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        GravFiles['260']
#		GravFiles['270']
    )
)

process.bbggAnalyzer = cms.EDAnalyzer('bbggPlotter'
)


process.p = cms.Path(process.bbggAnalyzer)
