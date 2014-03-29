// -*- C++ -*-
//
// Package:    BackgroundEstimationSkim
// Class:      BackgroundEstimationSkim
// 
/**\class BackgroundEstimationSkim BackgroundEstimationSkim.cc Bprime_kit/BackgroundEstimationSkim/src/BackgroundEstimationSkim.cc

Description: 
Analyzer class for Bprime -> b Higgs studies 
- National Taiwan University - 

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Jui-Fa Tsai
//         Created:  Tue Jul 16 19:48:47 CEST 2013
// Second Author:    Devdatta Majumder 
// $Id$
//
//

// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <vector>
#include <map>

// Root headers 
#include <TLorentzVector.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TEfficiency.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "../interface/format.h"
#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/TriggerBooking.h"
//#include "BpbH/BprimeTobH/interface/Njettiness.hh"
#include "BpbH/BprimeTobH/interface/Nsubjettiness.hh"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "BpbH/BprimeTobH/interface/TriggerSelector.h"
#include "BpbH/BprimeTobH/interface/VertexSelector.h"
#include "BpbH/BprimeTobH/interface/JetSelector.h"
#include "BpbH/BprimeTobH/interface/FatJetSelector.h"
//#include "BpbH/BprimeTobH/interface/HTSelector.h"
//#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"

#include "BpbH/BackgroundEstimationSkim/interface/reRegistGen.hh"
#include "BpbH/BackgroundEstimationSkim/interface/reRegistJet.hh"

///// Jet Correction
/*
#include "BpbH/BprimeTobHAnalysis/src/JMEUncertUtil.cc"
#include "BpbH/BprimeTobHAnalysis/src/BTagSFUtil.cc"
#include "BpbH/BprimeTobHAnalysis/src/ApplyBTagSF.cc"
#include "BpbH/BprimeTobHAnalysis/src/ApplyHiggsTagSF.cc"


#include "BpbH/BprimeTobHAnalysis/interface/JMEUncertUtil.h"
#include "BpbH/BprimeTobHAnalysis/interface/ApplyBTagSF.h"
//#include "BpbH/BackgroundEstimationSkim/interface/ApplyBTagSF.h"
#include "BpbH/BprimeTobHAnalysis/interface/ApplyHiggsTagSF.h"

#include "BpbH/BprimeTobHAnalysis/interface/HiggsBRscaleFactors.h" 
*/
//
// class declaration
//

class BackgroundEstimationSkim : public edm::EDAnalyzer{
	public:
		explicit BackgroundEstimationSkim(const edm::ParameterSet&);
		~BackgroundEstimationSkim();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
	
		edm::LumiReWeighting LumiWeights_; 

		// ----------member data ---------------------------

		//// Configurables 

		int                             maxEvents_; 
		const int                       reportEvery_; 
		const std::string               inputTTree_;
		const std::vector<std::string>  inputFiles_;
		const edm::ParameterSet         hltPaths_; 

		const double jetPtMin_ ; 
		const double jetPtMax_ ;
		const double bjetPtMin_ ; 
		const double bjetCSV_ ; 

		const edm::ParameterSet higgsJetSelParame_; 
		const edm::ParameterSet jetSelParams_; 
		const edm::ParameterSet bjetSelParams_; 
		const edm::ParameterSet evtSelParams_; 

		const edm::ParameterSet jmeParams_; 
		const double jesShift_;
		const double jerShift_; 
		const double SFbShift_;
		const double SFlShift_;

		const bool BuildMinTree_;

		TChain*            chain_;
		TTree*		   newtree;	

		GenInfoBranches    GenInfo;
		EvtInfoBranches    EvtInfo;
		VertexInfoBranches VtxInfo;
		JetInfoBranches    JetInfo;
		JetInfoBranches    FatJetInfo;
		JetInfoBranches    SubJetInfo;

		GenInfoBranches    newGenInfo;
		JetInfoBranches    newJetInfo;
		JetInfoBranches    newFatJetInfo;
		JetInfoBranches    newSubJetInfo;

		edm::Service<TFileService> fs; 
		
		// New branch

		long int _evtNo;
		int _lumiNo;
		int _runNo; 
		bool _mcFlag;
		double _PU;	
		double _evtwt; 

		long int evtNo_;
		int lumiNo_;
		int runNo_; 
		bool mcFlag_;
		double PU_;	
		double evtwt_; 

		TH1D* 	Evt_num;
		TH1D* 	cutFlow;
		TH1D* 	Higgs_num;
		TH1D* 	Higgs_pt;
		TH1D* 	Higgs_tau2Bytau1;
		TH1D* 	Higgs_mass;
		TH1D* 	Higgs_Prmass;
		TH1D* 	Higgs_subCSV;
		TH1D* 	bJet_num;
		TH1D* 	bJet_pt;
		TH1D* 	bJet_eta;
		TH1D* 	bJet_CSV;

};

//
// constructors and destructor
//
BackgroundEstimationSkim::BackgroundEstimationSkim(const edm::ParameterSet& iConfig) : 
	maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
	reportEvery_(iConfig.getParameter<int>("ReportEvery")),
	inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
	inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")),

	hltPaths_(iConfig.getParameter<edm::ParameterSet>("HLTPaths")),

	jetPtMin_(iConfig.getParameter<double>("JetPtMin")),
	jetPtMax_(iConfig.getParameter<double>("JetPtMax")),
  	bjetPtMin_(iConfig.getParameter<double>("BJetPtMin")),
  	bjetCSV_(iConfig.getParameter<double>("BJetCSV")),

	higgsJetSelParame_(iConfig.getParameter<edm::ParameterSet>("HiggsJetSelParams")),
	jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
	bjetSelParams_(iConfig.getParameter<edm::ParameterSet>("BJetSelParams")), 
	evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),

	jmeParams_(iConfig.getParameter<edm::ParameterSet>("JMEParams")),
	jesShift_(iConfig.getParameter<double>("JESShift")),
	jerShift_(iConfig.getParameter<double>("JERShift")),
	SFbShift_(iConfig.getParameter<double>("SFbShift")),
	SFlShift_(iConfig.getParameter<double>("SFlShift")),

	BuildMinTree_(iConfig.getParameter<bool>("BuildMinTree"))
{
	std::cout<<"Set Parameters !!!"<<std::endl;
} 



BackgroundEstimationSkim::~BackgroundEstimationSkim(){ 
	delete chain_;
}

// ------------ method called once each job just before starting event loop  ------------
void BackgroundEstimationSkim::beginJob(){ 
	chain_  = new TChain(inputTTree_.c_str());

	Evt_num		= fs->make<TH1D>("EvtInfo.Entries",	"", 1,   0, 1);
	Evt_num->GetXaxis()->SetBinLabel(1,"Entries");	
 
	for(unsigned i=0; i<inputFiles_.size(); ++i){
		chain_->Add(inputFiles_.at(i).c_str());
		TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");

		TH1D* _hevt = (TH1D*)f->Get("BprimebH/EvtInfo.Entries");
		int _events = _hevt->GetEntries();
		for( int _evt=0; _evt<_events; _evt++){ Evt_num->Fill(0); }

		f->Close();
	}

	GenInfo.Register(chain_);
	JetInfo.Register(chain_,"JetInfo");
	FatJetInfo.Register(chain_,"FatJetInfo");
	SubJetInfo.Register(chain_,"SubJetInfo");

	chain_->SetBranchAddress("EvtInfo.RunNo",	&_runNo);
	chain_->SetBranchAddress("EvtInfo.LumiNo",	&_lumiNo);
	chain_->SetBranchAddress("EvtInfo.EvtNo",	&_evtNo);
	chain_->SetBranchAddress("EvtInfo.McFlag",	&_mcFlag);
	chain_->SetBranchAddress("EvtInfo.PU",		&_PU);
	chain_->SetBranchAddress("EvtInfo.WeightEvt",	&_evtwt);

        if( BuildMinTree_ ) {
		newtree = fs->make<TTree>("tree", "") ;
		newtree->Branch("EvtInfo.EvtNo", 	&evtNo_, 	"EvtInfo.EvtNo/L"); // Store weight of Evt and PU for each event
		newtree->Branch("EvtInfo.LumiNo", 	&lumiNo_, 	"EvtInfo.LumiNo/I"); // Store weight of Evt and PU for each event
		newtree->Branch("EvtInfo.RunNo", 	&runNo_, 	"EvtInfo.RunNo/I"); // Store weight of Evt and PU for each event 
		newtree->Branch("EvtInfo.McFlag", 	&mcFlag_, 	"EvtInfo.McFlag/O"); // store weight of evt and pu for each event
		newtree->Branch("Evtinfo.PU", 		&PU_, 		"EvtInfo.PU/D"); // store weight of evt and pu for each event
		newtree->Branch("EvtInfo.WeightEvt", 	&evtwt_, 	"EvtInfo.WeightEvt/D"); // store weight of evt and pu for each event
		newGenInfo.RegisterTree(newtree);
		newJetInfo.RegisterTree(newtree,"JetInfo");
		newFatJetInfo.RegisterTree(newtree,"FatJetInfo");
		newSubJetInfo.RegisterTree(newtree,"SubJetInfo");
	}

	if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

	cutFlow		= fs->make<TH1D>("evtinfo.cutFlow",	"", 5,   0, 5); 
	Higgs_num	= fs->make<TH1D>("higgsjetinfo.Num",	"", 10,   0, 10); 
	Higgs_pt 	= fs->make<TH1D>("HiggsJetInfo.Pt",	"", 1500, 0, 1500);
	Higgs_tau2Bytau1= fs->make<TH1D>("HiggsJetInfo.Tau2ByTau1",	"", 10, 0, 1);
	Higgs_mass	= fs->make<TH1D>("HiggsJetInfo.Mass",	"", 3000, 0, 300);
	Higgs_Prmass	= fs->make<TH1D>("HiggsJetInfo.MassPruned",	"", 3000, 0, 300);
	Higgs_subCSV	= fs->make<TH1D>("HiggsJetInfo.SubCSV",	"", 110, 0, 1.1);
	bJet_num 	= fs->make<TH1D>("bJetInfo.Num",	"", 10, 0, 10); 
	bJet_pt		= fs->make<TH1D>("bJetInfo.Pt",		"", 1500, 0, 1500);
	bJet_eta	= fs->make<TH1D>("bJetInfo.Eta",	"", 600, -3, 3);
	bJet_CSV	= fs->make<TH1D>("bJetInfo.CSV", 	"", 100,  0, 1.);

	cutFlow->GetXaxis()->SetBinLabel(1,"All_Evt");	
	cutFlow->GetXaxis()->SetBinLabel(2,"bVeto_noWt");	

	return;  

}

// ------------ method called for each event  ------------
void BackgroundEstimationSkim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
	using namespace edm;
	using namespace std;

	if(  chain_ == 0) return;

	FatJetSelector fatjetSelHiggs(higgsJetSelParame_);
	JetSelector jetSelAK5(jetSelParams_); 
	pat::strbitset rethiggsjet = fatjetSelHiggs.getBitTemplate(); 
	pat::strbitset retjetidak5 = jetSelAK5.getBitTemplate(); 

	edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";
	
	//// Roop events ==================================================================================================	
	for(int entry=0; entry<maxEvents_; entry++){
		if( (entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_;
		chain_->GetEntry(entry);
		Evt_num->Fill(0);
		cutFlow->Fill(0);
		
		////  Higgs jets selection ================================================================================ 
		vector<TLorentzVector> higgsJets;
		for ( int i=0; i< FatJetInfo.Size; ++i ){
			rethiggsjet.set(false);
			if( fatjetSelHiggs( FatJetInfo, i, SubJetInfo, rethiggsjet)==0 ) continue; //higgs selection				
			TLorentzVector jet;
			jet.SetPtEtaPhiM(FatJetInfo.Pt[i], FatJetInfo.Eta[i], FatJetInfo.Phi[i], FatJetInfo.Mass[i]);
			higgsJets.push_back(jet);
			
			double tau2Bytau1 = FatJetInfo.tau2[i]/FatJetInfo.tau1[i];
			Higgs_pt->Fill(FatJetInfo.Pt[i]);
			Higgs_tau2Bytau1->Fill(tau2Bytau1);
			Higgs_mass->Fill(FatJetInfo.Mass[i]);
			Higgs_Prmass->Fill(FatJetInfo.MassPruned[i]);
			Higgs_subCSV->Fill(SubJetInfo.CombinedSVBJetTags[FatJetInfo.Jet_SubJet1Idx[i]]);
			Higgs_subCSV->Fill(SubJetInfo.CombinedSVBJetTags[FatJetInfo.Jet_SubJet2Idx[i]]);

		}
		Higgs_num->Fill(higgsJets.size());

		//// AK5 and bJet selection ================================================================================
		//// Preselection for AK5 Jet 
    		JetCollection myjets;
		for ( int i=0; i<JetInfo.Size; ++i ){ 
			retjetidak5.set(false);
			bool overlapWithCA8(false); 
			if( jetSelAK5(JetInfo, i, retjetidak5) == 0 ) continue; // AK5 pre-selection
		
			bJet_pt->Fill(JetInfo.Pt[i]);			
			bJet_CSV->Fill(JetInfo.CombinedSVBJetTags[i]);			

 			for ( unsigned int f=0; f<higgsJets.size(); ++f ){ // dR selection
				if( reco::deltaR(higgsJets[f].Eta(), higgsJets[f].Phi(), JetInfo.Eta[i], JetInfo.Phi[i])< 1.2 ){
					overlapWithCA8 = true; 
					break; 
				}
				else{
					overlapWithCA8 = false; 
				} 
			}
      			Jet thisjet(JetInfo, i);
      			if( !overlapWithCA8 ) myjets.push_back(thisjet) ;
			 
		}
		const int nbjets = myjets.size();
		bJet_num->Fill(nbjets);
		if( nbjets>0 ) continue;
		cutFlow->Fill(1);
			
		//// Store new tree, new branch with Jet correction  ==================================================================================================== 
		if( BuildMinTree_ ){
			evtNo_ 	= _evtNo;
			lumiNo_ = _lumiNo;
			runNo_ 	= _runNo; 
			mcFlag_ = _mcFlag;
			PU_ 	= _PU;	
			evtwt_ 	= _evtwt; 
			reRegistGen(GenInfo, newGenInfo); 	
			reRegistJet(JetInfo, newJetInfo);	
			reRegistJet(FatJetInfo, newFatJetInfo);	
			reRegistJet(SubJetInfo, newSubJetInfo);
			newtree->Fill();
		}
	} //// entry loop 

}

// ------------ method called once each job just after ending the event loop  ------------
void BackgroundEstimationSkim::endJob(){ 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BackgroundEstimationSkim::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BackgroundEstimationSkim);
