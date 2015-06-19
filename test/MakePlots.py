import FWCore.ParameterSet.Config as cms
from microAOD_RadFiles import *
from microAOD_GravFiles import *

process = cms.Process("bbggAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 2000 )


## Available mass points:
# RadFiles: 320, 340, 350, 400, 600, 650, 700
# GravFiles: 260, 270, 280, 320, 350, 500, 550
NonResPhys14 = 'file:/afs/cern.ch/work/r/rateixei/work/DiHiggs/FLASHggPreSel/CMSSW_7_4_0_pre9/src/flashgg/MicroAOD/test/hhbbgg_hggVtx.root'

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
	RadFiles['700']
#        NonResPhys14
    )
)

process.load("flashgg.bbggPlotter.CfiFile_cfi")
process.bbggAnalyzer.OutFileName = cms.untracked.string('test_Rad700.root')


process.p = cms.Path(process.bbggAnalyzer)
