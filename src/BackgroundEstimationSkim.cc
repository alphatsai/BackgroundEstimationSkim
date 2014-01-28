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

#include "BpbH/BprimeTobHAnalysis/interface/HiggsBRscaleFactors.h" 

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
		//const std::vector<string>       hltPaths_; 
		const edm::ParameterSet         hltPaths_; 
		const int                       doPUReweighting_;
		const std::string               file_PUDistMC_;
		const std::string               file_PUDistData_;
		const std::string               hist_PUDistMC_;
		const std::string               hist_PUDistData_;

		const edm::ParameterSet higgsJetSelParame_; 
		const edm::ParameterSet jetSelParams_; 
		const edm::ParameterSet bjetSelParams_; 
		const edm::ParameterSet evtSelParams_; 
		
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
		
		int evtNum_;
		bool isData_; 
		bool McFlag_; 
		double evtwt_; 
		double puweight_;
		double evtwtPu_; 

		TH1D* 	Evt_num;
		TH1D* 	cutFlow;
		TH1D* 	AK5_num;
		TH1D* 	AK5_pt;
		TH1D* 	AK5_eta;
		TH1D* 	AK5_CSV;
		TH1D* 	bJet_num;
		TH1D* 	bJet_pt;
		TH1D* 	bJet_eta;
		TH1D* 	bJet_CSV;
		TH1D* 	bJetVeto_pt;
		TH1D* 	bJetVeto_eta;
		TH1D* 	bJetVeto_CSV;
		TH1D* 	bJetVeto_num; 
		TH1D* 	bJetVeto_lead_pt;
		TH1D* 	bJetVeto_lead_eta;
		TH1D* 	bJetVeto_lead_CSV;
		TH1D* 	bJetVetoMatchCA8_num;
		TH1D* 	bJet1_pt;
		TH1D* 	bJet1_eta;
		TH1D* 	bJet1_CSV;
		TH1D* 	bJet1_lead_pt;
		TH1D* 	bJet1_lead_eta;
		TH1D* 	bJet1_lead_CSV;
		TH1D* 	bJet2_pt;
		TH1D* 	bJet2_eta;
		TH1D* 	bJet2_CSV;
		TH1D* 	bJet2_lead_pt;
		TH1D* 	bJet2_lead_eta;
		TH1D* 	bJet2_lead_CSV;

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

	higgsJetSelParame_(iConfig.getParameter<edm::ParameterSet>("HiggsJetSelParams")),
	jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
	bjetSelParams_(iConfig.getParameter<edm::ParameterSet>("BJetSelParams")), 
	evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),

	BuildMinTree_(iConfig.getParameter<bool>("BuildMinTree")), 

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

        if( BuildMinTree_ ) {
		newtree = fs->make<TTree>("tree", "") ; 
		newtree->Branch("EvtInfo.EvtNum", &evtNum_, "EvtInfo.EvtNum/I"); // Store weight of Evt and PU for each event
		newtree->Branch("EvtInfo.McFlag", &McFlag_, "EvtInfo.McFlag/O"); // Store weight of Evt and PU for each event
		newtree->Branch("EvtInfo.WeightEvtPU", &evtwtPu_, "EvtInfo.WeightEvtPU/D"); // Store weight of Evt and PU for each event
		newGenInfo.RegisterTree(newtree);
		newJetInfo.RegisterTree(newtree,"JetInfo");
		newFatJetInfo.RegisterTree(newtree,"FatJetInfo");
		newSubJetInfo.RegisterTree(newtree,"SubJetInfo");
	}

	if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

	Evt_num		= fs->make<TH1D>("EvtInfo.Entries",	"", 1,   0, 1); 
	cutFlow		= fs->make<TH1D>("EvtInfo.CutFlow",	"", 4,   0, 4); 
	AK5_num		= fs->make<TH1D>("AK5JetInfo.Num",	"", 10,   0, 10); 
	AK5_pt 		= fs->make<TH1D>("AK5JetInfo.Pt",	"", 1500, 0, 1500);
	AK5_eta 	= fs->make<TH1D>("AK5JetInfo.Eta",	"", 600, -3, 3);
	AK5_CSV		= fs->make<TH1D>("AK5JetInfo.CSV", 	"", 100,  0, 1.);
	bJet_num 	= fs->make<TH1D>("bJetInfo.Num",	"", 10, 0, 10); 
	bJet_pt		= fs->make<TH1D>("bJetInfo.Pt",		"", 1500, 0, 1500);
	bJet_eta	= fs->make<TH1D>("bJetInfo.Eta",	"", 600, -3, 3);
	bJet_CSV	= fs->make<TH1D>("bJetInfo.CSV", 	"", 100,  0, 1.);
	bJetVeto_num 	= fs->make<TH1D>("bJetInfo.Num.Veto",	"", 10, 0, 10); 
	bJetVeto_pt	= fs->make<TH1D>("bJetInfo.Pt.Veto",	"", 1500, 0, 1500);
	bJetVeto_eta	= fs->make<TH1D>("bJetInfo.Eta.Veto",	"", 600, -3, 3);
	bJetVeto_CSV	= fs->make<TH1D>("bJetInfo.CSV.Veto", 	"", 100,  0, 1.);
	bJetVeto_lead_pt	= fs->make<TH1D>("bJetInfo.lead.Pt.Veto",	"", 1500, 0, 1500);
	bJetVeto_lead_eta	= fs->make<TH1D>("bJetInfo.lead.Eta.Veto",	"", 600, -3, 3);
	bJetVeto_lead_CSV	= fs->make<TH1D>("bJetInfo.lead.CSV.Veto", 	"", 100,  0, 1.);
	bJetVetoMatchCA8_num 	= fs->make<TH1D>("bJetInfo.NumMatchToCA8.Veto",	"", 10, 0, 10);

	bJet1_pt	= fs->make<TH1D>("bJetInfo.Pt.1",	"", 1500, 0, 1500);
	bJet1_eta	= fs->make<TH1D>("bJetInfo.Eta.1",	"", 600, -3, 3);
	bJet1_CSV	= fs->make<TH1D>("bJetInfo.CSV.1", 	"", 100,  0, 1.);
	bJet1_lead_pt	= fs->make<TH1D>("bJetInfo.lead.Pt.1",	"", 1500, 0, 1500);
	bJet1_lead_eta	= fs->make<TH1D>("bJetInfo.lead.Eta.1",	"", 600, -3, 3);
	bJet1_lead_CSV	= fs->make<TH1D>("bJetInfo.lead.CSV.1", 	"", 100,  0, 1.);
	bJet2_pt	= fs->make<TH1D>("bJetInfo.Pt.2",	"", 1500, 0, 1500);
	bJet2_eta	= fs->make<TH1D>("bJetInfo.Eta.2",	"", 600, -3, 3);
	bJet2_CSV	= fs->make<TH1D>("bJetInfo.CSV.2", 	"", 100,  0, 1.);
	bJet2_lead_pt	= fs->make<TH1D>("bJetInfo.lead.Pt.2",	"", 1500, 0, 1500);
	bJet2_lead_eta	= fs->make<TH1D>("bJetInfo.lead.Eta.2",	"", 600, -3, 3);
	bJet2_lead_CSV	= fs->make<TH1D>("bJetInfo.lead.CSV.2", 	"", 100,  0, 1.);

	Evt_num->Sumw2();
	cutFlow->Sumw2();
	AK5_num->Sumw2();
	AK5_pt->Sumw2();
	AK5_eta->Sumw2();
	AK5_CSV->Sumw2();
	bJet_num->Sumw2();
	bJet_pt->Sumw2();
	bJet_eta->Sumw2();
	bJet_CSV->Sumw2();
	bJetVeto_pt->Sumw2();
	bJetVeto_eta->Sumw2();
	bJetVeto_CSV->Sumw2();
	bJetVeto_lead_pt->Sumw2();
	bJetVeto_lead_eta->Sumw2();
	bJetVeto_lead_CSV->Sumw2();
	bJetVeto_num->Sumw2();
	bJetVetoMatchCA8_num->Sumw2();
	
	bJet1_pt->Sumw2();
	bJet1_eta->Sumw2();
	bJet1_CSV->Sumw2();
	bJet1_lead_pt->Sumw2();
	bJet1_lead_eta->Sumw2();
	bJet1_lead_CSV->Sumw2();
	bJet2_pt->Sumw2();
	bJet2_eta->Sumw2();
	bJet2_CSV->Sumw2();
	bJet2_lead_pt->Sumw2();
	bJet2_lead_eta->Sumw2();
	bJet2_lead_CSV->Sumw2();

	Evt_num->GetXaxis()->SetBinLabel(1,"Entries");	
	cutFlow->GetXaxis()->SetBinLabel(1,"All_Evt");	
	cutFlow->GetXaxis()->SetBinLabel(2,"Trigger_Sel");	
	cutFlow->GetXaxis()->SetBinLabel(3,"Vertex_Sel");	
	cutFlow->GetXaxis()->SetBinLabel(4,"BJetVeto");	

	return;  

}

// ------------ method called for each event  ------------
void BackgroundEstimationSkim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
	using namespace edm;
	using namespace std;

	if(  chain_ == 0) return;

	FatJetSelector fatjetSelHiggs(higgsJetSelParame_);
	JetSelector jetSelAK5(jetSelParams_); 
	JetSelector jetSelBJet(bjetSelParams_); 
	pat::strbitset rethiggsjet = fatjetSelHiggs.getBitTemplate(); 
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
		Evt_num->Fill(0);

		//// Event variables 
		std::vector<TLorentzVector>p4_Jets; 
		std::vector<TLorentzVector>p4_bJets; 

		bool passHLT(false); 
		int nGoodVtxs(0);
		int nAK5(0); 
		int nbjets(0); 
		int nbjetsNoCA8(0); 

		isData_  = EvtInfo.McFlag ? 0 : 1; 
		if( !isData_ ) evtwt_    = EvtInfo.Weight; 
		if( doPUReweighting_ && !isData_ ) puweight_ = LumiWeights_.weight(EvtInfo.TrueIT[0]); 
		evtwt_ *= puweight_; 

		//// Higgs BR reweighting, for b'b'>bHbH sample
		if ( !isData_ ) {
			//int nbprimeBH(0), nbprimeBZ(0), nbprimeTW(0); 
			for (int igen=0; igen < GenInfo.Size; ++igen) {
				if ( GenInfo.Status[igen] == 3 && TMath::Abs(GenInfo.PdgID[igen]) == 25 && GenInfo.nDa[igen] >= 2 ) { //// H found
					int higgsDau0 = abs(GenInfo.Da0PdgID[igen]);
					int higgsDau1 = abs(GenInfo.Da1PdgID[igen]);
					int higgsDau = higgsDau0+higgsDau1; 
					if(higgsDau==10) evtwt_ *= HiggsBRscaleFactors::higgsBBSf;
					if(higgsDau==30) evtwt_ *= HiggsBRscaleFactors::higgsTauTauSf;
					if(higgsDau==26) evtwt_ *= HiggsBRscaleFactors::higgsMuMuSf;
					if(higgsDau==8)  evtwt_ *= HiggsBRscaleFactors::higgsCCSf;
					if(higgsDau==6)  evtwt_ *= HiggsBRscaleFactors::higgsSSSf;
					if(higgsDau==16) evtwt_ *= HiggsBRscaleFactors::higgsTTSf;
					if(higgsDau==42) evtwt_ *= HiggsBRscaleFactors::higgsGGSf;
					if(higgsDau==44) evtwt_ *= HiggsBRscaleFactors::higgsGammaGammaSf;
					if(higgsDau==45) evtwt_ *= HiggsBRscaleFactors::higgsZGammaSf;
					if(higgsDau==48) evtwt_ *= HiggsBRscaleFactors::higgsWWSf;
					if(higgsDau==46) evtwt_ *= HiggsBRscaleFactors::higgsZZSf; 
				}
			}
		}

		cutFlow->Fill(double(0),evtwt_);
		//// Trigger selection =====================================================================================
		TriggerSelector trigSel(hltPaths_); 
		passHLT = trigSel.getTrigDecision(EvtInfo); 
		if( !passHLT ) continue; 
		cutFlow->Fill(double(1),evtwt_);

		//// Vertex selection =====================================================================================
		VertexSelector vtxSel(VtxInfo); 
		nGoodVtxs = vtxSel.NGoodVtxs(); 
		if( nGoodVtxs < 1){ edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex "; continue; }
		cutFlow->Fill(double(2),evtwt_);

		//// Recall data no Jet =====================================================================================
		if( isData_ ){
			if( JetInfo.Size == 0 ) fout << EvtInfo.RunNo << " " << EvtInfo.LumiNo << " " << EvtInfo.EvtNo << std::endl; 
		}

		////  Higgs jets selection ================================================================================ 
		vector<TLorentzVector> higgsJets;
		for ( int i=0; i< FatJetInfo.Size; ++i ){ 
			rethiggsjet.set(false);
			if( fatjetSelHiggs(FatJetInfo, i, SubJetInfo, rethiggsjet)==0 ) continue; //higgs selection				
			TLorentzVector jet;
			jet.SetPtEtaPhiM(FatJetInfo.Pt[i], FatJetInfo.Eta[i], FatJetInfo.Phi[i], FatJetInfo.Mass[i]);
			higgsJets.push_back(jet);	
		}

		//// AK5 and bJet selection ================================================================================ 
		for ( int i=0; i<JetInfo.Size; ++i ){ 
			retjetidak5.set(false);
			retjetidbjet.set(false);
			bool overlapWithCA8(false); 

			if( jetSelAK5(JetInfo, i,retjetidak5) == 0 ) continue; 
			++nAK5; //AK5 Jet selection
			AK5_pt->Fill(double(JetInfo.Pt[i]),evtwt_); 
			AK5_eta->Fill(double(JetInfo.Eta[i]),evtwt_); 
			AK5_CSV->Fill(double(JetInfo.CombinedSVBJetTags[i]),evtwt_);
 
			if( jetSelBJet(JetInfo, i,retjetidbjet) == 0 ) continue; 
			++nbjets; //b Jet selection
			for ( unsigned int f=0; f<higgsJets.size(); ++f ){ // dR selection
				if( reco::deltaR(higgsJets[f].Eta(), higgsJets[f].Phi(), JetInfo.Eta[i], JetInfo.Phi[i])< 1.2 ){
					overlapWithCA8 = true; 
					break; 
				}
				else{
					overlapWithCA8 = false; 
				} 
			} 
			if( overlapWithCA8 ) continue; 
			++nbjetsNoCA8; //// NO overlap with CA8  
			bJet_pt->Fill(double(JetInfo.Pt[i]),evtwt_); 
			bJet_eta->Fill(double(JetInfo.Eta[i]),evtwt_); 
			bJet_CSV->Fill(double(JetInfo.CombinedSVBJetTags[i]),evtwt_);

		} //// AK5 jets END 

		//// Store event infomation and bJet veto ================================================================== 
		AK5_num->Fill(double(nAK5),evtwt_);
		bJet_num->Fill(double(nbjetsNoCA8),evtwt_);
		//bJet_num->Fill(double(nbjets),evtwt_);

		if(nbjetsNoCA8>=1){
			int lead=-1;
			int leadPt=-1;
			for( int i=0; i<JetInfo.Size; i++){
				bJet1_pt->Fill(double(JetInfo.Pt[i]),evtwt_);
				bJet1_eta->Fill(double(JetInfo.Eta[i]),evtwt_);
				bJet1_CSV->Fill(double(JetInfo.CombinedSVBJetTags[i]),evtwt_);
				if( leadPt<JetInfo.Pt[i] ){ 
					lead=i;
					leadPt=JetInfo.Pt[i];
				}
			}
			bJet1_lead_pt->Fill(double(JetInfo.Pt[lead]),evtwt_);
			bJet1_lead_eta->Fill(double(JetInfo.Eta[lead]),evtwt_);
			bJet1_lead_CSV->Fill(double(JetInfo.CombinedSVBJetTags[lead]),evtwt_);
				
		}
		if(nbjetsNoCA8>=2){
			int lead=-1;
			int leadPt=-1;
			for( int i=0; i<JetInfo.Size; i++){
				bJet2_pt->Fill(double(JetInfo.Pt[i]),evtwt_);
				bJet2_eta->Fill(double(JetInfo.Eta[i]),evtwt_);
				bJet2_CSV->Fill(double(JetInfo.CombinedSVBJetTags[i]),evtwt_);
				if( leadPt<JetInfo.Pt[i] ){ 
					lead=i;
					leadPt=JetInfo.Pt[i];
				}
			}	
			bJet2_lead_pt->Fill(double(JetInfo.Pt[lead]),evtwt_);
			bJet2_lead_eta->Fill(double(JetInfo.Eta[lead]),evtwt_);
			bJet2_lead_CSV->Fill(double(JetInfo.CombinedSVBJetTags[lead]),evtwt_);
		}


		if( nbjetsNoCA8>0 ) continue;  //Pass b-Jet veto
		int lead=-1;
		int leadPt=-1;
		bJetVeto_num->Fill(double(nbjetsNoCA8),evtwt_);
		bJetVetoMatchCA8_num->Fill(double(nbjets),evtwt_);
		cutFlow->Fill(double(3),evtwt_);
		for( int i=0; i<JetInfo.Size; i++){
			bJetVeto_pt->Fill(double(JetInfo.Pt[i]),evtwt_);
			bJetVeto_eta->Fill(double(JetInfo.Eta[i]),evtwt_);
			bJetVeto_CSV->Fill(double(JetInfo.CombinedSVBJetTags[i]),evtwt_);
			if( leadPt<JetInfo.Pt[i] ){ 
				lead=i;
				leadPt=JetInfo.Pt[i];
			}
		}
		bJetVeto_lead_pt->Fill(double(JetInfo.Pt[lead]),evtwt_);
		bJetVeto_lead_eta->Fill(double(JetInfo.Eta[lead]),evtwt_);
		bJetVeto_lead_CSV->Fill(double(JetInfo.CombinedSVBJetTags[lead]),evtwt_);

		//// Fill mini tree after bJetVeto =============================================================================================
		if( BuildMinTree_ ){
			if( isData_ ){
				McFlag_=0;
			}else{
				McFlag_=1;
			}
			evtwtPu_ = evtwt_;
			evtNum_	 = maxEvents_;
			reRegistGen(GenInfo,newGenInfo); 	
			reRegistJet(JetInfo,newJetInfo);	
			reRegistJet(FatJetInfo,newFatJetInfo);	
			reRegistJet(SubJetInfo,newSubJetInfo);	
			newtree->Fill();
		}
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
