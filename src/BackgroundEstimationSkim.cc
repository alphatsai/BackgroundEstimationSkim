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

		const double jetPtMin_ ; 
		const double jetPtMax_ ;
		const double bjetPtMin_ ; 
		const double bjetCSV_ ; 

		const edm::ParameterSet higgsJetSelParame_; 
		const edm::ParameterSet jetSelParams_; 
		const edm::ParameterSet bjetSelParams_; 
		const edm::ParameterSet evtSelParams_; 

		const bool selAK5_;
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
		
		bool isData_; 
		double evtwt_; 
		double puweight_;

		// New branch
		long int evtNum_;
		int runNum_;
		int lumiNum_;
		int totalEvts_;
		bool McFlag_; 
		double evtwtPu_; 
		double evtwt_new; 

		TH1D* 	Evt_num;
		TH1D* 	cutFlow;
		TH1D* 	cutFlow_unWt;
		TH1D* 	Higgs_num;
		TH1D* 	Higgs_pt;
		TH1D* 	Higgs_tau2Bytau1;
		TH1D* 	Higgs_mass;
		TH1D* 	AK5_num;
		TH1D* 	AK5_pt;
		TH1D* 	AK5_eta;
		TH1D* 	AK5_CSV;
		TH1D* 	bJet_num;
		TH1D* 	bJet_pt;
		TH1D* 	bJet_eta;
		TH1D* 	bJet_CSV;
		TH1D* 	bJet1_pt;
		TH1D* 	bJet1_eta;
		TH1D* 	bJet1_CSV;
		TH1D* 	bJet1_lead_pt;
		TH1D* 	bJet1_lead_eta;
		TH1D* 	bJet1_lead_CSV;

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

	jetPtMin_(iConfig.getParameter<double>("JetPtMin")),
	jetPtMax_(iConfig.getParameter<double>("JetPtMax")),
  	bjetPtMin_(iConfig.getParameter<double>("BJetPtMin")),
  	bjetCSV_(iConfig.getParameter<double>("BJetCSV")),

	higgsJetSelParame_(iConfig.getParameter<edm::ParameterSet>("HiggsJetSelParams")),
	jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
	evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),

	selAK5_(iConfig.getParameter<bool>("SelectAK5")),
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
		/*newtree->Branch("EvtInfo.EvtTotal", &totalEvts_, "EvtInfo.EvtTotal/I"); 
		newtree->Branch("EvtInfo.EvtNo", &evtNum_, "EvtInfo.EvtNo/L"); 
		newtree->Branch("EvtInfo.LumiNo", &lumiNum_, "EvtInfo.LumiNo/I"); 
		newtree->Branch("EvtInfo.RunNo", &runNum_, "EvtInfo.RunNo/I"); 
		newtree->Branch("EvtInfo.McFlag", &McFlag_, "EvtInfo.McFlag/O"); 
		newGenInfo.RegisterTree(newtree);
		newJetInfo.RegisterTree(newtree,"JetInfo");
		newFatJetInfo.RegisterTree(newtree,"FatJetInfo");
		newSubJetInfo.RegisterTree(newtree,"SubJetInfo");*/
		newtree = chain_->CloneTree(0);
		newtree->Branch("EvtInfo.EvtTotal", &totalEvts_, "EvtInfo.EvtTotal/I"); 
		newtree->Branch("EvtInfo.PU", &evtwtPu_, "EvtInfo.PU/D"); // Store weight of PU for each event
		newtree->Branch("EvtInfo.WeightEvt", &evtwt_new, "EvtInfo.WeightEvt/D"); // Store weight of Evt 
	}

	if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

	Evt_num		= fs->make<TH1D>("EvtInfo.Entries",	"", 1,   0, 1); 
	cutFlow		= fs->make<TH1D>("EvtInfo.CutFlow",	"", 5,   0, 5); 
	cutFlow_unWt	= fs->make<TH1D>("EvtInfo.CutFlow.UnWt","", 5,   0, 5); 
	Higgs_num	= fs->make<TH1D>("HiggsJetInfo.Num",	"", 10,   0, 10); 
	Higgs_pt 	= fs->make<TH1D>("HiggsJetInfo.Pt",	"", 1500, 0, 1500);
	Higgs_tau2Bytau1= fs->make<TH1D>("HiggsJetInfo.Tau2ByTau1",	"", 10, 0, 1);
	Higgs_mass	= fs->make<TH1D>("HiggsJetInfo.Mass",	"", 3000, 0, 300);
	AK5_num		= fs->make<TH1D>("AK5JetInfo.Num",	"", 10,   0, 10); 
	AK5_pt 		= fs->make<TH1D>("AK5JetInfo.Pt",	"", 1500, 0, 1500);
	AK5_eta 	= fs->make<TH1D>("AK5JetInfo.Eta",	"", 600, -3, 3);
	AK5_CSV		= fs->make<TH1D>("AK5JetInfo.CSV", 	"", 100,  0, 1.);
	bJet_num 	= fs->make<TH1D>("bJetInfo.Num",	"", 10, 0, 10); 
	bJet_pt		= fs->make<TH1D>("bJetInfo.Pt",		"", 1500, 0, 1500);
	bJet_eta	= fs->make<TH1D>("bJetInfo.Eta",	"", 600, -3, 3);
	bJet_CSV	= fs->make<TH1D>("bJetInfo.CSV", 	"", 100,  0, 1.);

	Evt_num->Sumw2();
	cutFlow->Sumw2();
	cutFlow_unWt->Sumw2();
	AK5_num->Sumw2();
	AK5_pt->Sumw2();
	AK5_eta->Sumw2();
	AK5_CSV->Sumw2();
	bJet_num->Sumw2();
	bJet_pt->Sumw2();
	bJet_eta->Sumw2();
	bJet_CSV->Sumw2();
	
	Evt_num->GetXaxis()->SetBinLabel(1,"Entries");	
	cutFlow->GetXaxis()->SetBinLabel(1,"All_Evt");	
	cutFlow->GetXaxis()->SetBinLabel(2,"Trigger_Sel");	
	cutFlow->GetXaxis()->SetBinLabel(3,"Vertex_Sel");	
	cutFlow->GetXaxis()->SetBinLabel(4,"Skim_HJet_Sel");	
	cutFlow->GetXaxis()->SetBinLabel(5,"Skim_BJet_Sel");	
	cutFlow_unWt->GetXaxis()->SetBinLabel(1,"All_Evt");	
	cutFlow_unWt->GetXaxis()->SetBinLabel(2,"Trigger_Sel");	
	cutFlow_unWt->GetXaxis()->SetBinLabel(3,"Vertex_Sel");	
	cutFlow_unWt->GetXaxis()->SetBinLabel(4,"Skim_HJet_Sel");	
	cutFlow_unWt->GetXaxis()->SetBinLabel(5,"Skim_BJet_Sel");	

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
		
		double evtwt_old(1);
		bool passHLT(false); 
		int nGoodVtxs(0);
		int nAK5(0); 
		int nHiggs(0); 

		isData_  = EvtInfo.McFlag ? 0 : 1; 
		if( !isData_ ) evtwt_old    = EvtInfo.Weight; 
		if( doPUReweighting_ && !isData_ ) puweight_ = LumiWeights_.weight(EvtInfo.TrueIT[0]); 

		//// Higgs BR reweighting, for b'b'>bHbH sample
		double br_(1);
		if ( !isData_ ) {
			//int nbprimeBH(0), nbprimeBZ(0), nbprimeTW(0); 
			for (int igen=0; igen < GenInfo.Size; ++igen) {
				if ( GenInfo.Status[igen] == 3 && TMath::Abs(GenInfo.PdgID[igen]) == 25 && GenInfo.nDa[igen] >= 2 ) { //// H found
					int higgsDau0 = abs(GenInfo.Da0PdgID[igen]);
					int higgsDau1 = abs(GenInfo.Da1PdgID[igen]);
					int higgsDau = higgsDau0+higgsDau1; 
					if(higgsDau==10) br_ *= HiggsBRscaleFactors::higgsBBSf;
					if(higgsDau==30) br_ *= HiggsBRscaleFactors::higgsTauTauSf;
					if(higgsDau==26) br_ *= HiggsBRscaleFactors::higgsMuMuSf;
					if(higgsDau==8)  br_ *= HiggsBRscaleFactors::higgsCCSf;
					if(higgsDau==6)  br_ *= HiggsBRscaleFactors::higgsSSSf;
					if(higgsDau==16) br_ *= HiggsBRscaleFactors::higgsTTSf;
					if(higgsDau==42) br_ *= HiggsBRscaleFactors::higgsGGSf;
					if(higgsDau==44) br_ *= HiggsBRscaleFactors::higgsGammaGammaSf;
					if(higgsDau==45) br_ *= HiggsBRscaleFactors::higgsZGammaSf;
					if(higgsDau==48) br_ *= HiggsBRscaleFactors::higgsWWSf;
					if(higgsDau==46) br_ *= HiggsBRscaleFactors::higgsZZSf; 
				}
			}
		}
		evtwt_ = evtwt_old*puweight_*br_; 

		//cout<<entry<<" "<<evtwt_old<<" "<<puweight_<<" "<<br_<<" "<<evtwt_old*puweight_*br_<<" "<<evtwt_<<endl; 	

		cutFlow_unWt->Fill(0);
		cutFlow->Fill("All_Evt", evtwt_);
		//// Trigger selection =====================================================================================
		TriggerSelector trigSel(hltPaths_); 
		passHLT = trigSel.getTrigDecision(EvtInfo); 
		if( !passHLT ) continue; 
		cutFlow_unWt->Fill(1);
		cutFlow->Fill("Trigger_Sel", evtwt_);

		//// Vertex selection =====================================================================================
		VertexSelector vtxSel(VtxInfo); 
		nGoodVtxs = vtxSel.NGoodVtxs(); 
		if( nGoodVtxs < 1){ edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex "; continue; }
		cutFlow_unWt->Fill(2);
		cutFlow->Fill("Vertex_Sel", evtwt_);

		//// Recall data no Jet =====================================================================================
		if( isData_ ){
			if( JetInfo.Size == 0 ) fout << EvtInfo.RunNo << " " << EvtInfo.LumiNo << " " << EvtInfo.EvtNo << std::endl; 
		}
		//cout<<entry<<" "<<EvtInfo.RunNo<<" "<<EvtInfo.LumiNo<<" "<<EvtInfo.EvtNo<<endl;
		
		////  Higgs jets selection ================================================================================ 
		for ( int i=0; i< FatJetInfo.Size; ++i ){
			rethiggsjet.set(false);
			if( fatjetSelHiggs( FatJetInfo, i, SubJetInfo, rethiggsjet)==0 ) continue; //higgs selection				
			nHiggs++;
		
			Higgs_pt->Fill(FatJetInfo.Pt[i]);
			double tau2Bytau1 = FatJetInfo.tau2[i]/FatJetInfo.tau1[i];
			Higgs_tau2Bytau1->Fill(tau2Bytau1);
			Higgs_mass->Fill(FatJetInfo.Mass[i]);
		}
		Higgs_num->Fill(nHiggs);

		if( nHiggs<1 ) continue;
		cutFlow_unWt->Fill(3);
		cutFlow->Fill("Skim_HJet_Sel", evtwt_);

		//// AK5 and bJet selection ================================================================================
		//// Preselection for AK5 Jet 
		for ( int i=0; i<JetInfo.Size; ++i ){ 
			retjetidak5.set(false);
			if( jetSelAK5(JetInfo, i, retjetidak5) == 0 ) continue; // AK5 pre-selection
			nAK5++;
			
			AK5_pt->Fill(JetInfo.Pt[i]);	
			AK5_eta->Fill(JetInfo.Eta[i]);
			AK5_CSV->Fill(JetInfo.CombinedSVBJetTags[i]);	
		}
		AK5_num->Fill(nAK5);
		if( selAK5_ && nAK5<2 ) continue;
		cutFlow_unWt->Fill(4);
		cutFlow->Fill("Skim_BJet_Sel", evtwt_);
		
		//// Store new tree, new branch with Jet correction  ==================================================================================================== 

		// Fill new tree, new branch
		if( BuildMinTree_ ){
			/*if( isData_ ){
				McFlag_=0;
			}else{
				McFlag_=1;
			}
			evtNum_    = EvtInfo.EvtNo;	
			lumiNum_   = EvtInfo.LumiNo;	
			runNum_    = EvtInfo.RunNo;	
			evtwt_new  = evtwt_old * br_;
			evtwtPu_   = puweight_;
			totalEvts_ = maxEvents_;
			reRegistGen(GenInfo, newGenInfo); 	
			reRegistJet(JetInfo, newJetInfo);	
			reRegistJet(FatJetInfo, newFatJetInfo);	
			reRegistJet(SubJetInfo, newSubJetInfo);*/
			evtwt_new  = evtwt_old * br_;
			evtwtPu_   = puweight_;
			totalEvts_ = maxEvents_;
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
