import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

from inputFiles_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'bVeto.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
    )
#####################################################
options.setDefault('maxEvents', -1000) 

options.parseArguments()

process = cms.Process("bVeto")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'), 
    ) 
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) ) # Leave it this way. 
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) # show cpu time
process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFilename) 
    )

from BpbH.BprimeTobHAnalysis.EventSelector_cfi import * 
from BpbH.BprimeTobHAnalysis.JMEUncertUntilParameters_cfi import * 

process.bVeto = cms.EDAnalyzer('BackgroundEstimationSkim',
    MaxEvents           = cms.int32(options.maxEvents),
    ReportEvery         = cms.int32(options.reportEvery),  
    InputTTree          = cms.string('skim/tree'),
    InputFiles          = cms.vstring(FileNames), 
	
    JetSelParams        = defaultJetSelectionParameters.clone(
		jetPtMin = cms.double(30),
    		jetCSVDiscMin   = cms.double(0.244),   
	),
    BJetSelParams       = defaultBJetSelectionParameters.clone(),
    HiggsJetSelParams   = defaultHiggsJetSelectionParameters.clone(
    		fatJetPrunedMassMin = cms.double(90),
    		fatJetPrunedMassMax = cms.double(140),
    		fatJetMassMin       = cms.double(-1),
    		fatJetMassMax       = cms.double(10000),
    		subjet1CSVDiscMin   = cms.double(0.679),
		subjet1CSVDiscMax   = cms.double(1.000),
    		subjet2CSVDiscMin   = cms.double(0.679),
    		subjet2CSVDiscMax   = cms.double(1.000),
	), 
    HTSelParams         = defaultHTSelectionParameters.clone(),
    EvtSelParams        = defaultEventSelectionParameters.clone(),

    #BuildMinTree        = cms.bool(False),
    BuildMinTree        = cms.bool(True),
    ) 
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.bVeto)

