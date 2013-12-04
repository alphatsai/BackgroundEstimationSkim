#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"

EventSelector::EventSelector (const edm::ParameterSet& iConfig, const EvtInfoBranches &evt, const VertexInfoBranches &vtx,
    const JetInfoBranches &jet, const JetInfoBranches &fatjet, const JetInfoBranches &subjet) : 
  evtInfo_(&evt),
  vtxInfo_(&vtx),
  jetInfo_(&jet),
  fatjetInfo_(&fatjet),
  subjetInfo_(&subjet),
  passes_(0) 
{

  trigSelector_    = new TriggerSelector(iConfig.getParameter<edm::ParameterSet>("TrigSelParams")) ; 
  vtxSelector_     = new VertexSelector(vtx) ;
  jetSelector_     = new JetSelector(iConfig.getParameter<edm::ParameterSet>("JetSelParams")) ; 
  bjetSelector_    = new JetSelector(iConfig.getParameter<edm::ParameterSet>("BJetSelParams")) ; 
  fatjetSelector_  = new FatJetSelector(iConfig.getParameter<edm::ParameterSet>("FatJetSelParams")) ; 
  higgsjetSelector_= new FatJetSelector(iConfig.getParameter<edm::ParameterSet>("HiggsJetSelParams")) ; 
  htSelector_      = new HTSelector(iConfig.getParameter<edm::ParameterSet>("HTSelParams")) ;  

}

EventSelector::~EventSelector () {}


