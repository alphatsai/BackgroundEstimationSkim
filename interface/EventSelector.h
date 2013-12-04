#ifndef BPRIMETOBHANALYSIS_INTERFACE_EVENTSELECTOR_H
#define BPRIMETOBHANALYSIS_INTERFACE_EVENTSELECTOR_H

#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "BpbH/BprimeTobH/interface/TriggerSelector.h"
#include "BpbH/BprimeTobH/interface/VertexSelector.h"
#include "BpbH/BprimeTobH/interface/JetSelector.h"
#include "BpbH/BprimeTobH/interface/FatJetSelector.h"
#include "BpbH/BprimeTobH/interface/HTSelector.h"

class EventSelector {

  public : 

    enum SELECTYPES_t {EVTSEL, VTXSEL, JETSEL, BJETSEL, FATJETSEL, HIGGSJETSEL, HTSEL, LEPTONVETO} ; 

    EventSelector (edm::ParameterSet const& params, 
        const EvtInfoBranches &evt, const VertexInfoBranches &vtx, 
        const JetInfoBranches &jet, const JetInfoBranches &fatjet, const JetInfoBranches &subjet) ; 
    ~EventSelector() ; 

  private: 

    //// Configurables 
    //edm::ParameterSet evtSelParams_;
    const EvtInfoBranches* evtInfo_ ; 
    const VertexInfoBranches* vtxInfo_;
    const JetInfoBranches* jetInfo_;
    const JetInfoBranches* fatjetInfo_; 
    const JetInfoBranches* subjetInfo_; 

    //// Other members 
    edm::ParameterSet trigSelParams_;
    edm::ParameterSet vtxSelParams_;
    edm::ParameterSet jetSelParams_;
    edm::ParameterSet bjetSelParams_;
    edm::ParameterSet fatjetSelParams_;
    edm::ParameterSet higgsjetSelParams_; 

    const TriggerSelector* trigSelector_ ;
    const VertexSelector* vtxSelector_ ;
    const JetSelector* jetSelector_ ;
    const JetSelector* bjetSelector_ ;
    const FatJetSelector* fatjetSelector_ ; 
    const FatJetSelector* higgsjetSelector_ ; 
    const HTSelector* htSelector_ ; 

    bool passes_ ; 

} ; 

#endif 
