import FWCore.ParameterSet.Config as cms

from BpbH.BprimeTobH.TriggerSelector_cfi import *
from BpbH.BprimeTobH.JetSelector_cfi import *
from BpbH.BprimeTobH.BJetSelector_cfi import *
from BpbH.BprimeTobH.FatJetSelector_cfi import * 
from BpbH.BprimeTobH.HiggsJetSelector_cfi import * 
from BpbH.BprimeTobH.HTSelector_cfi import * 

defaultEventSelectionParameters = cms.PSet (
    TrigSelParams     = defaultTriggerSelectionParameters.clone(), 
    JetSelParams      = defaultJetSelectionParameters.clone(), 
    BJetSelParams     = defaultBJetSelectionParameters.clone(), 
    FatJetSelParams   = defaultFatJetSelectionParameters.clone(), 
    HiggsJetSelParams = defaultHiggsJetSelectionParameters.clone(), 
    HTSelParams       = defaultHTSelectionParameters.clone(),
    )
