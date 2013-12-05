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
// Original Author:  Eleni Petrakou,27 2-020,+41227674870,
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "../interface/format.h"
#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/TriggerBooking.h"
#include "BpbH/BprimeTobH/interface/Njettiness.hh"
#include "BpbH/BprimeTobH/interface/Nsubjettiness.hh"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "BpbH/BprimeTobH/interface/TriggerSelector.h"
#include "BpbH/BprimeTobH/interface/VertexSelector.h"
#include "BpbH/BprimeTobH/interface/JetSelector.h"
//#include "BpbH/BprimeTobH/interface/FatJetSelector.h"
//#include "BpbH/BprimeTobH/interface/HTSelector.h"
//#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"

#include "BpbH/BackgroundEstimationSkim/interface/reRegistGen.hh"
#include "BpbH/BackgroundEstimationSkim/interface/reRegistJet.hh"

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
	
		double dPhi(double p1,double p2);	
		double dR(double e1, double e2, double p1, double p2);	

		edm::LumiReWeighting LumiWeights_; 

		// ----------member data ---------------------------

		//// Configurables 

		int                             maxEvents_; 
		const int                       reportEvery_; 
		const std::string               inputTTree_;
		const std::vector<std::string>  inputFiles_;
		//const std::vector<string>       hltPaths_; 
		const edm::ParameterSet         hltPaths_; 
		const int                       doPUReweighting_;
		const std::string               file_PUDistMC_;
		const std::string               file_PUDistData_;
		const std::string               hist_PUDistMC_;
		const std::string               hist_PUDistData_;

		const edm::ParameterSet jetSelParams_; 
		const edm::ParameterSet bjetSelParams_; 
		const edm::ParameterSet evtSelParams_; 

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

		bool isData_; 
		double evtwt_; 
		double puweight_;
		float evtwtPu_; 

		TH1D* 	AK5_num;
		TH1D* 	AK5_pt;
		TH1D* 	AK5_CSV;
		TH1D* 	bJet_num;
		TH1D* 	bJet_pt;
		TH1D* 	bJet_CSV;
		TH1D* 	bJetVeto_num; 
		TH1D* 	bJetVetoMatchCA8_num;

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

	doPUReweighting_(iConfig.getParameter<bool>("DoPUReweighting")), 
	file_PUDistMC_(iConfig.getParameter<std::string>("File_PUDistMC")),
	file_PUDistData_(iConfig.getParameter<std::string>("File_PUDistData")),
	hist_PUDistMC_(iConfig.getParameter<std::string>("Hist_PUDistMC")),
	hist_PUDistData_(iConfig.getParameter<std::string>("Hist_PUDistData")),

	jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
	bjetSelParams_(iConfig.getParameter<edm::ParameterSet>("BJetSelParams")), 
	evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),
 
	isData_(0),
	evtwt_(1), 
	puweight_(1)  
{ 

	if( doPUReweighting_) LumiWeights_ = edm::LumiReWeighting(file_PUDistMC_, file_PUDistData_, hist_PUDistMC_, hist_PUDistData_);

}


BackgroundEstimationSkim::~BackgroundEstimationSkim(){ 
	delete chain_;
}

// ------------ method called once each job just before starting event loop  ------------
void BackgroundEstimationSkim::beginJob(){ 
	chain_  = new TChain(inputTTree_.c_str());
        newtree = fs->make<TTree>("tree", "") ; 

	for(unsigned i=0; i<inputFiles_.size(); ++i){
		chain_->Add(inputFiles_.at(i).c_str());
		TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
		f->Close();
	}

	EvtInfo.Register(chain_);
	VtxInfo.Register(chain_);
	GenInfo.Register(chain_);
	JetInfo.Register(chain_,"JetInfo");
	FatJetInfo.Register(chain_,"FatJetInfo");
	SubJetInfo.Register(chain_,"SubJetInfo");

	newtree->Branch("EvtInfo.WeightEvtPU", &evtwtPu_, "EvtInfo.WeightEvtPU/F"); // Store weight of Evt and PU for each event
	newGenInfo.RegisterTree(newtree);
	newJetInfo.RegisterTree(newtree,"JetInfo");
	newFatJetInfo.RegisterTree(newtree,"FatJetInfo");
	newSubJetInfo.RegisterTree(newtree,"SubJetInfo");

	if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

	AK5_num		= fs->make<TH1D>("AK5JetInfo.Num",	"", 10,   0, 10); 
	AK5_pt 		= fs->make<TH1D>("AK5JetInfo.Pt",	"", 1500, 0, 1500);
	AK5_CSV		= fs->make<TH1D>("AK5JetInfo.CSV", 	"", 100,  0, 1.);
	bJet_num 	= fs->make<TH1D>("bJetInfo.Num",	"", 10, 0, 10); 
	bJet_pt		= fs->make<TH1D>("bJetInfo.Pt",		"", 1500, 0, 1500);
	bJet_CSV	= fs->make<TH1D>("bJetInfo.CSV", 	"", 100,  0, 1.);
	bJetVeto_num 		= fs->make<TH1D>("bJetInfo.Num.Veto",		"", 10, 0, 10); 
	bJetVetoMatchCA8_num 	= fs->make<TH1D>("bJetInfo.NumMatchToCA8.Veto",	"", 10, 0, 10);

	return;  

}

double BackgroundEstimationSkim::dPhi(double p1, double p2){
        double dp = p1 - p2;
        if( fabs(dp+3.14159265358979323846*2.) < fabs(dp)) dp += 3.14159265358979323846*2.;
        else
                if( fabs(dp-3.14159265358979323846*2.) < fabs(dp)) dp -= 3.14159265358979323846*2.;
        return fabs(dp);
}
double BackgroundEstimationSkim::dR(double e1, double e2, double p1, double p2){
	return sqrt(pow(e1-e2,2)+pow(dPhi(p1,p2),2));
}

// ------------ method called for each event  ------------
void BackgroundEstimationSkim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
	using namespace edm;
	using namespace std;

	if(  chain_ == 0) return;

	JetSelector jetSelAK5(jetSelParams_); 
	JetSelector jetSelBJet(bjetSelParams_); 
	pat::strbitset retjetidak5 = jetSelAK5.getBitTemplate(); 
	pat::strbitset retjetidbjet = jetSelBJet.getBitTemplate(); 

	ofstream fout("Evt_NoJets.txt"); 
	if(  isData_ ){
		fout << "EvtInfo.RunNo " << " EvtInfo.LumiNo " << " EvtInfo.EvtNo " << std::endl;
	}

	edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";
	
	//// Roop events ==================================================================================================	
	for(int entry=0; entry<maxEvents_; entry++){
		if( (entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_; 
		chain_->GetEntry(entry);

		//// Event variables 
		std::vector<TLorentzVector>p4_Jets; 
		std::vector<TLorentzVector>p4_bJets; 

		bool passHLT(false); 

		int nGoodVtxs(0);
		int nAK5(0); 
		int nbjets(0); 
		int nbjetsNoCA8(0); 

		chain_->GetEntry(entry);

		isData_   = EvtInfo.McFlag ? 0 : 1; 
		if( !isData_ ) evtwt_    = EvtInfo.Weight; 
		if( doPUReweighting_ && !isData_ ) puweight_ = LumiWeights_.weight(EvtInfo.TrueIT[0]); 

		//// Vertex selection =====================================================================================
		VertexSelector vtxSel(VtxInfo); 
		nGoodVtxs = vtxSel.NGoodVtxs(); 
		if( nGoodVtxs < 1){ edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex "; continue; }

		//// Trigger selection =====================================================================================
		TriggerSelector trigSel(hltPaths_); 
		passHLT = trigSel.getTrigDecision(EvtInfo); 
		if( !passHLT ) continue; 

		//// Recall data no Jet =====================================================================================
		if( isData_ ){
			if( JetInfo.Size == 0 ) fout << EvtInfo.RunNo << " " << EvtInfo.LumiNo << " " << EvtInfo.EvtNo << std::endl; 
		}

		//// AK5 and bJet selection ================================================================================ 
		evtwt_ *= puweight_; 

		for ( int i=0; i<JetInfo.Size; ++i ){ 
			retjetidak5.set(false);
			retjetidbjet.set(false);
			bool overlapWithCA8(false); 

			if( jetSelAK5(JetInfo, i,retjetidak5) == 0 ) continue; 
			++nAK5; //AK5 Jet selection
			AK5_pt->Fill(JetInfo.Pt[i]); 
			AK5_CSV->Fill(JetInfo.CombinedSVBJetTags[i]);
 
			if( jetSelBJet(JetInfo, i,retjetidbjet) == 0 ) continue; 
			++nbjets; //b Jet selection
			bJet_pt->Fill(JetInfo.Pt[i]); 
			bJet_CSV->Fill(JetInfo.CombinedSVBJetTags[i]);

			for ( int f = 0; f < FatJetInfo.Size; ++f ){ // dR selection
				if( dR(FatJetInfo.Eta[f], JetInfo.Eta[i], FatJetInfo.Phi[f], JetInfo.Phi[i])< 1.2 ){
					overlapWithCA8 = true; 
					break; 
				}
				else{
					overlapWithCA8 = false; 
				} 
			} 
			if( overlapWithCA8 ) continue; 
			++nbjetsNoCA8; //// NO overlap with CA8  
		} //// AK5 jets END 

		//// Store mini tree and event infomation ================================================================== 
		AK5_num->Fill(nAK5);
		bJet_num->Fill(nbjets);

		if( nbjetsNoCA8>0 ) continue;
		bJetVeto_num->Fill(nbjetsNoCA8);
		bJetVetoMatchCA8_num->Fill(nbjets);
		evtwtPu_=evtwt_;
		reRegistGen(GenInfo,newGenInfo); 	
		reRegistJet(JetInfo,newJetInfo);	
		reRegistJet(FatJetInfo,newFatJetInfo);	
		reRegistJet(SubJetInfo,newSubJetInfo);	
		newtree->Fill();
		//cout<<entry<<" PASS!!!"<<endl;
	} //// entry loop 

	fout.close(); 

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
