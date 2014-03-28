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

		const edm::ParameterSet higgsJetSelParame_; 
		const edm::ParameterSet jetSelParams_; 
		const edm::ParameterSet bjetSelParams_; 
		const edm::ParameterSet evtSelParams_; 

		const bool BuildMinTree_;

		TChain*            chain_;
		TTree*		   newtree;	

		edm::Service<TFileService> fs; 
		
		// branch

		JetInfoBranches    JetInfo;
		JetInfoBranches    FatJetInfo;
		JetInfoBranches    SubJetInfo;

		TH1D* 	h_cutflow;
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

	higgsJetSelParame_(iConfig.getParameter<edm::ParameterSet>("HiggsJetSelParams")),
	jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
	bjetSelParams_(iConfig.getParameter<edm::ParameterSet>("BJetSelParams")), 
	evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),

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
	
	h_cutflow = fs->make<TH1D>("h_cutflow" ,"Cut flow" ,6 ,0. ,6); 
	h_cutflow->GetXaxis()->SetBinLabel(1, "BeforeNtuplizer") ; 
	h_cutflow->GetXaxis()->SetBinLabel(2, "AfterNtuplizer") ; 
	h_cutflow->GetXaxis()->SetBinLabel(3, "BeforeSkim") ; 
	h_cutflow->GetXaxis()->SetBinLabel(4, "AfterSkim") ; 
	h_cutflow->GetXaxis()->SetBinLabel(5, "Before bVeto") ; 
	h_cutflow->GetXaxis()->SetBinLabel(6, "After bVeto") ; 

	for(unsigned i=0; i<inputFiles_.size(); ++i){
		chain_->Add(inputFiles_.at(i).c_str());
		TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
		TH1F* h_events = (TH1F*)f->Get("skim/h_cutflow") ; 
		f->Close();
		h_cutflow->Fill("BeforeNtuplizer", h_events->GetBinContent(1)) ; 
		h_cutflow->Fill("AfterNtuplizer", h_events->GetBinContent(2)) ;
		h_cutflow->Fill("BeforeSkim", h_events->GetBinContent(3)) ; 
		h_cutflow->Fill("AfterSkim", h_events->GetBinContent(4)) ;
	}

	JetInfo.Register(chain_,"JetInfo");
	FatJetInfo.Register(chain_,"FatJetInfo");
	SubJetInfo.Register(chain_,"SubJetInfo");

        if( BuildMinTree_ ) {
		fs->cd();
		newtree = fs->make<TTree>("tree", "") ;
		newtree = chain_->CloneTree(0);
	}

	if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

	Higgs_num	= fs->make<TH1D>("higgsjetinfo_Num",	"", 10,   0, 10); 
	Higgs_pt 	= fs->make<TH1D>("HiggsJetInfo_Pt",	"", 1500, 0, 1500);
	Higgs_tau2Bytau1= fs->make<TH1D>("HiggsJetInfo_Tau2ByTau1",	"", 10, 0, 1);
	Higgs_mass	= fs->make<TH1D>("HiggsJetInfo_Mass",	"", 3000, 0, 300);
	Higgs_Prmass	= fs->make<TH1D>("HiggsJetInfo_MassPruned",	"", 3000, 0, 300);
	Higgs_subCSV	= fs->make<TH1D>("HiggsJetInfo_SubCSV",	"", 110, 0, 1.1);
	bJet_num 	= fs->make<TH1D>("bJetInfo_Num",	"", 10, 0, 10); 
	bJet_pt		= fs->make<TH1D>("bJetInfo_Pt",		"", 1500, 0, 1500);
	bJet_eta	= fs->make<TH1D>("bJetInfo_Eta",	"", 600, -3, 3);
	bJet_CSV	= fs->make<TH1D>("bJetInfo_CSV", 	"", 100,  0, 1.);

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
		h_cutflow->Fill(4);
		
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
			 
			bJet_pt->Fill(JetInfo.Pt[i]);			
			bJet_eta->Fill(JetInfo.Eta[i]);			
			bJet_CSV->Fill(JetInfo.CombinedSVBJetTags[i]);			

		}
		const int nbjets = myjets.size();
		bJet_num->Fill(nbjets);
		if( nbjets>0 ) continue;
		h_cutflow->Fill(5);
			
		//// Store new tree, new branch with Jet correction  ==================================================================================================== 
		if( BuildMinTree_ ){
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
