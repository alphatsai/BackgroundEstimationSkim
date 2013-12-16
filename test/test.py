import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

#from BackgroundEstimationSkimv1.BackgroundEstimationSkim.BpBpToBHBHinc.BprimeBprimeTobHbHinc_M_800_cfi import * 
#from BackgroundEstimationSkimv1.BackgroundEstimationSkim.Data.JetHT_Run2012BCD_cfi import * 
from inputFiles_cfi import * 

FileNames = [
#'file:BprimeTobH_v1_10_1_HIn.root'
'file:/afs/cern.ch/user/j/jtsai/eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_10_1_HIn.root',
'file:/afs/cern.ch/user/j/jtsai/eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_1_1_Pf5.root',
'file:/afs/cern.ch/user/j/jtsai/eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_2_1_Y96.root',
'file:/afs/cern.ch/user/j/jtsai/eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_3_1_h6K.root',
'file:/afs/cern.ch/user/j/jtsai/eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_4_1_gI9.root',
'file:/afs/cern.ch/user/j/jtsai/eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_5_1_QTk.root',
'file:/afs/cern.ch/user/j/jtsai/eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_6_1_H6G.root',
'file:/afs/cern.ch/user/j/jtsai/eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_7_1_YrU.root',
'file:/afs/cern.ch/user/j/jtsai/eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_8_1_KTf.root',
'file:/afs/cern.ch/user/j/jtsai/eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_9_1_G4q.root'
]

options = VarParsing('python')

options.register('outFilename', 'test.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
    )
options.register('jetPtMin', 50.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum jet Pt"
    )
options.register('jetPtMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum jet Pt"
    )
options.register('bJetPtMin', 80.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum b jet Pt"
    )
options.register('fatJetPtMin', 300.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet Pt"
    )
options.register('fatJetPtMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet Pt"
    )
options.register('fatJetMassMin', 100.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet mass"
    )
options.register('fatJetMassMax', 150.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet mass"
    )
options.register('fatJetPrunedMassMin', 75.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet pruned mass"
    )
options.register('fatJetPrunedMassMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet pruned mass"
    )
options.register('fatJetTau2ByTau1Max', 0.5,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet tau2/tau1"
    )
options.register('subjet1CSVDiscMin', 0.679,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum subjet1 b discriminator"
    )
options.register('subjet1CSVDiscMax', 1.000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum subjet1 b discriminator"
    )
options.register('subjet2CSVDiscMin', 0.679,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum subjet2 b discriminator"
    )
options.register('subjet2CSVDiscMax', 1.000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum subjet2 b discriminator"
    )
options.register('hTMin', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum HT"
    )
options.register('hTMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum HT"
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

process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFilename) 
    )

from BpbH.BprimeTobH.TriggerSelector_cfi import * 
from BpbH.BprimeTobH.HiggsJetSelector_cfi import * 
from BpbH.BprimeTobH.HTSelector_cfi import * 
from BpbH.BprimeTobHAnalysis.EventSelector_cfi import * 

process.BprimebH = cms.EDAnalyzer('BackgroundEstimationSkim',
    MaxEvents           = cms.int32(options.maxEvents),
    ReportEvery         = cms.int32(options.reportEvery),  
    InputTTree          = cms.string('ntuple/tree'),
    InputFiles          = cms.vstring(FileNames), 
    #HLTPaths            = cms.vint32(3225,4136,4137,5089,5537,5538), 
    HLTPaths            = defaultTriggerSelectionParameters.clone(), 
    DoPUReweighting     = cms.bool(options.doPUReweighting),
    File_PUDistMC       = cms.string('pileup_Data_Summer12_53X_S10.root'),
    File_PUDistData     = cms.string('pileup_Data_Summer12_53X_S10.root'),
    Hist_PUDistMC       = cms.string('pileup_mc'),
    Hist_PUDistData     = cms.string('pileup_data'),
    JetSelParams        = defaultJetSelectionParameters.clone(
		jetPtMin = cms.double(30)
	), 
    BJetSelParams       = defaultBJetSelectionParameters.clone(
		jetCSVDiscMin = cms.double(0.244),	
		jetPtMin = cms.double(30)
	), #alpha 
    FatJetSelParams     = defaultFatJetSelectionParameters.clone(), 
    HiggsJetSelParams   = defaultHiggsJetSelectionParameters.clone(), 
    HTSelParams         = defaultHTSelectionParameters.clone(),
    EvtSelParams        = defaultEventSelectionParameters.clone(),
    BuildMinTree        = cms.bool(False),
    #BuildMinTree        = cms.bool(True),
    ) 

process.p = cms.Path(process.BprimebH)

