import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

#from BackgroundEstimationSkimv1.BackgroundEstimationSkim.BpBpToBHBHinc.BprimeBprimeTobHbHinc_M_800_cfi import * 
#from BackgroundEstimationSkimv1.BackgroundEstimationSkim.Data.JetHT_Run2012BCD_cfi import * 
from inputFiles_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'bprimeTobH.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
    )
options.register('doHLTselection', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do HLT selection"
)
options.register('doGoodVertex', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do good vertex selection"
)
options.register('doPUReweighting', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do pileup reweighting"
)

options.setDefault('maxEvents', -1000) 

options.parseArguments()

process = cms.Process("BprimebH")

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

process.BprimebH = cms.EDAnalyzer('BackgroundEstimationSkim',
    MaxEvents           = cms.int32(options.maxEvents),
    ReportEvery         = cms.int32(options.reportEvery),  
    InputTTree          = cms.string('ntuple/tree'),
    InputFiles          = cms.vstring(FileNames), 

    HLTPaths            = defaultTriggerSelectionParameters.clone(), 
    DoHLTSelect     	= cms.bool(options.doHLTselection),
    DoGoodVtxSelect     = cms.bool(options.doGoodVertex),
    DoPUReweighting     = cms.bool(options.doPUReweighting),
    File_PUDistMC       = cms.string('pileup_Data_Summer12_53X_S10.root'),
    File_PUDistData     = cms.string('pileup_Data_Summer12_53X_S10.root'),
    Hist_PUDistMC       = cms.string('pileup_mc'),
    Hist_PUDistData     = cms.string('pileup_data'),
	
    HJetPtMin		= cms.double(150),
    HJetPtMax		= cms.double(100000),
    HJetAbsEtaMin	= cms.double(-1),
    HJetAbsEtaMax	= cms.double(2.4),
    JetPtMin		= cms.double(30),
    JetPtMax		= cms.double(100000),
    JetAbsEtaMin	= cms.double(-1),
    JetAbsEtaMax	= cms.double(2.4),
    numbJetMin		= cms.int32(0),
    numHiggsJetMin	= cms.int32(1),

    #BuildMinTree        = cms.bool(False),
    BuildMinTree        = cms.bool(True),
    ) 
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.BprimebH)

