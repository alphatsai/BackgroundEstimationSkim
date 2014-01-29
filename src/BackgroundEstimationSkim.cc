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
#include "BpbH/BprimeTobHAnalysis/src/JMEUncertUtil.cc"
#include "BpbH/BprimeTobHAnalysis/src/BTagSFUtil.cc"
//#include "BpbH/BprimeTobHAnalysis/src/ApplyBTagSF.cc"
#include "BpbH/BprimeTobHAnalysis/src/ApplyHiggsTagSF.cc"

#include "BpbH/BprimeTobHAnalysis/interface/JMEUncertUtil.h"
//#include "BpbH/BprimeTobHAnalysis/interface/ApplyBTagSF.h"
#include "BpbH/BackgroundEstimationSkim/interface/ApplyBTagSF.h"
#include "BpbH/BprimeTobHAnalysis/interface/ApplyHiggsTagSF.h"

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
		
		bool isData_; 
		double evtwt_; 
		double puweight_;
		double higgsTagCorr_;

		// New branch
		int evtNum_;
		bool McFlag_; 
		bool CSVLSF_[128]; 
		bool CSVMSF_[128]; 
		bool CSVTSF_[128]; 
		double evtwtPu_; 
		double evtwtHiggsTagCorr_; 

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

	BuildMinTree_(iConfig.getParameter<bool>("BuildMinTree")), 

	isData_(0),
	evtwt_(1), 
	puweight_(1),
	higgsTagCorr_(1)  
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
		newtree->Branch("EvtInfo.WeightHiggsTagCorr", &evtwtHiggsTagCorr_, "EvtInfo.WeightHiggsTagCorr/D"); // Store weight of Evt and PU for each event
		newGenInfo.RegisterTree(newtree);
		newJetInfo.RegisterTree(newtree,"JetInfo");
		newtree->Branch("JetInfo.CSVLSF", &CSVLSF_[0], "JetInfo.CSVLSF[JetInfo.Size]/O"); // Store weight of Evt and PU for each event
		newtree->Branch("JetInfo.CSVMSF", &CSVMSF_[0], "JetInfo.CSVMSF[JetInfo.Size]/O"); // Store weight of Evt and PU for each event
		newtree->Branch("JetInfo.CSVTSF", &CSVTSF_[0], "JetInfo.CSVTSF[JetInfo.Size]/O"); // Store weight of Evt and PU for each event
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

    		JetCollection ak5jets_corr; 

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

			if ( !isData_ ) { //// Apply Higgs-tagging scale factor 
				int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[i];
				int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[i];
				ApplyHiggsTagSF* higgsTagSF = new ApplyHiggsTagSF(double(SubJetInfo.Pt[iSubJet1]), double(SubJetInfo.Pt[iSubJet2]), 
						double(SubJetInfo.Eta[iSubJet1]), double(SubJetInfo.Eta[iSubJet2]),
						SubJetInfo.GenFlavor[iSubJet1], SubJetInfo.GenFlavor[iSubJet2], 
						SubJetInfo.CombinedSVBJetTags[iSubJet1], SubJetInfo.CombinedSVBJetTags[iSubJet2]) ; 
				higgsTagCorr_ = higgsTagSF->GetHiggsTagSF() ;
				delete higgsTagSF ; 
			}
		}


		//// AK5 and bJet selection ================================================================================
		//// Preselection for AK5 Jet 
    		JetCollection myjets ;
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
		}


		//// Jet Correction: JER, JES
		JetCollection ak5jets_tmp; 
		if ( !isData_) {
			// Only AK5 jets not overlapping with Higgs jets 
			JMEUncertUtil* jmeUtil_jer = new JMEUncertUtil(jmeParams_, EvtInfo, myjets, "JER", jerShift_) ; 
			JetCollection ak5jets_jer = jmeUtil_jer->GetModifiedJetColl() ; 
			delete jmeUtil_jer ; 

			JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, EvtInfo, ak5jets_jer, "JES", jesShift_) ; 
			ak5jets_tmp = jmeUtil_jes->GetModifiedJetColl() ; 
			delete jmeUtil_jes ; 
		}
		else {
			// Only AK5 jets not overlapping with Higgs jets 
			JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, EvtInfo, myjets, "JES", jesShift_) ; 
			ak5jets_tmp = jmeUtil_jes->GetModifiedJetColl() ; 
			delete jmeUtil_jes ; 
		}

		//// AK5 Jet Selection
    		for (JetCollection::const_iterator ijet = ak5jets_tmp.begin(); ijet != ak5jets_tmp.end(); ++ijet) {
      			if ( ijet->Pt() < jetPtMin_ || ijet->Pt() > jetPtMax_) continue ;
			AK5_pt->Fill(double(ijet->Pt()),evtwt_); 
			AK5_eta->Fill(double(ijet->Eta()),evtwt_); 
			AK5_CSV->Fill(double(ijet->CombinedSVBJetTags()),evtwt_);
      			ak5jets_corr.push_back(*ijet);
			++nAK5; 
		}

 
		//// Jet Correction: BTag
		JetCollection ak5jets_bTag_corr;
		if ( !isData_ ){
			string algo;
			if( bjetCSV_ == 0.244 ){ 
				algo = "CSVL"; 
			}else if( bjetCSV_ == 0.679 ){
				algo = "CSVM"; 
			}else if( bjetCSV_ == 0.898 ){
				algo = "CSVT"; 
			}else{
				edm::LogError("ApplyBTagSF") << " Wrong b-tagging algo_ chosen: " << algo << ". Choose between CSVL:0.244, CSVM:0.679, or CSVT:0.898" ; 
			}
			ApplyBTagSF * btagsf =  new ApplyBTagSF(ak5jets_corr, bjetCSV_, algo, SFbShift_, SFlShift_) ;  
			ak5jets_bTag_corr =  btagsf->getBtaggedJetsWithSF () ; 
			delete btagsf ; 

		}else{
			ak5jets_bTag_corr = ak5jets_corr;
		}

		//// b Jet Selection
		JetCollection bjets_corr;
    		for (JetCollection::const_iterator ijet = ak5jets_bTag_corr.begin(); ijet != ak5jets_bTag_corr.end(); ++ijet) {
			if( ijet->Pt() < bjetPtMin_ || ijet->CombinedSVBJetTags() < bjetCSV_ ) continue;  
			bJet_pt->Fill(double(ijet->Pt()),evtwt_); 
			bJet_eta->Fill(double(ijet->Eta()),evtwt_); 
			bJet_CSV->Fill(double(ijet->CombinedSVBJetTags()),evtwt_);
			bjets_corr.push_back(*ijet);
			++nbjetsNoCA8; //// NO overlap with CA8  
		}  

		
		//// Fill histogram  ================================================================== 
		AK5_num->Fill(double(nAK5),evtwt_);
		bJet_num->Fill(double(nbjetsNoCA8),evtwt_);
		if(nbjetsNoCA8>=1){
			int ii=0;
			int lead=-1;
			int leadPt=-1;
    			for (JetCollection::const_iterator ijet = bjets_corr.begin(); ijet != bjets_corr.end(); ++ijet) {
				bJet1_pt->Fill(double(ijet->Pt()),evtwt_);
				bJet1_eta->Fill(double(ijet->Eta()),evtwt_);
				bJet1_CSV->Fill(double(ijet->CombinedSVBJetTags()),evtwt_);
				if( leadPt < ijet->Pt() ){ 
					lead = ii;
					leadPt = ijet->Pt();
				}
				ii++;	
			}
			bJet1_lead_pt->Fill(double(bjets_corr[lead].Pt()),evtwt_);
			bJet1_lead_eta->Fill(double(bjets_corr[lead].Eta()),evtwt_);
			bJet1_lead_CSV->Fill(double(bjets_corr[lead].CombinedSVBJetTags()),evtwt_);
				
		}
		if(nbjetsNoCA8>=2){
			int ii=0;
			int lead=-1;
			int leadPt=-1;
    			for (JetCollection::const_iterator ijet = bjets_corr.begin(); ijet != bjets_corr.end(); ++ijet) {
				bJet2_pt->Fill(double(ijet->Pt()),evtwt_);
				bJet2_eta->Fill(double(ijet->Eta()),evtwt_);
				bJet2_CSV->Fill(double(ijet->CombinedSVBJetTags()),evtwt_);
				if( leadPt<ijet->Pt() ){ 
					lead=ii;
					leadPt=ijet->Pt();
				}
				ii++;	
			}	
			bJet2_lead_pt->Fill(double(bjets_corr[lead].Pt()),evtwt_);
			bJet2_lead_eta->Fill(double(bjets_corr[lead].Eta()),evtwt_);
			bJet2_lead_CSV->Fill(double(bjets_corr[lead].CombinedSVBJetTags()),evtwt_);
		}

		//// Store new tree, new branch with Jet correction  ==================================================================================================== 
		if( nbjetsNoCA8>0 ) continue;  //Pass b-Jet veto
		bJetVeto_num->Fill(double(nbjetsNoCA8),evtwt_);
		//bJetVetoMatchCA8_num->Fill(double(nbjets),evtwt_);
		cutFlow->Fill(double(3),evtwt_);

		JetCollection allak5jets, allak5jets_corr;
		for ( int i=0; i<JetInfo.Size; ++i ){
      			Jet thisjet(JetInfo, i);
			allak5jets.push_back(thisjet);
		}
		//// 1. Jet Correction: JER, JES
		if ( !isData_) {
			// Only AK5 jets not overlapping with Higgs jets 
			JMEUncertUtil* jmeUtil_jer = new JMEUncertUtil(jmeParams_, EvtInfo, allak5jets, "JER", jerShift_) ; 
			JetCollection allak5jets_jer = jmeUtil_jer->GetModifiedJetColl() ; 
			delete jmeUtil_jer ; 

			JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, EvtInfo, allak5jets_jer, "JES", jesShift_) ; 
			allak5jets_corr = jmeUtil_jes->GetModifiedJetColl() ; 
			delete jmeUtil_jes ; 
		}else {
			// Only AK5 jets not overlapping with Higgs jets 
			JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, EvtInfo, allak5jets, "JES", jesShift_) ; 
			allak5jets_corr = jmeUtil_jes->GetModifiedJetColl() ; 
			delete jmeUtil_jes ; 

		}
		//// 2. Jet Correction: BTag, (For bVeto, don't need to apply) 
		ApplyBTagSF * btagCSVLsf =  new ApplyBTagSF(allak5jets_corr, 0.244, "CSVL", SFbShift_, SFlShift_) ;  
		ApplyBTagSF * btagCSVMsf =  new ApplyBTagSF(allak5jets_corr, 0.679, "CSVM", SFbShift_, SFlShift_) ;  
		ApplyBTagSF * btagCSVTsf =  new ApplyBTagSF(allak5jets_corr, 0.898, "CSVT", SFbShift_, SFlShift_) ;  
		const int allak5jets_num = allak5jets_corr.size();
		bool passBTagCSVL[allak5jets_num], passBTagCSVM[allak5jets_num], passBTagCSVT[allak5jets_num];	
		for( int ijet=0; ijet<allak5jets_num; ijet++){
			passBTagCSVL[ijet] =  btagCSVLsf->getBtaggedSF(ijet); 
			passBTagCSVM[ijet] =  btagCSVMsf->getBtaggedSF(ijet); 
			passBTagCSVT[ijet] =  btagCSVTsf->getBtaggedSF(ijet); 
		}
		delete btagCSVLsf;
		delete btagCSVMsf;
		delete btagCSVTsf;

		// Fill bVeto Histogram	
		int ii=0;
		int lead=-1;
		int leadPt=-1;
    		for (JetCollection::const_iterator ijet = allak5jets_corr.begin(); ijet != allak5jets_corr.end(); ++ijet) {
			bJetVeto_pt->Fill(double(ijet->Pt()),evtwt_);
			bJetVeto_eta->Fill(double(ijet->Eta()),evtwt_);
			bJetVeto_CSV->Fill(double(ijet->CombinedSVBJetTags()),evtwt_);
			if( leadPt<ijet->Pt() ){ 
				lead=ii;
				leadPt=ijet->Pt();
			}
			ii++;	
		}
		bJetVeto_lead_pt->Fill(double(allak5jets_corr[lead].Pt()),evtwt_);
		bJetVeto_lead_eta->Fill(double(allak5jets_corr[lead].Eta()),evtwt_);
		bJetVeto_lead_CSV->Fill(double(allak5jets_corr[lead].CombinedSVBJetTags()),evtwt_);


		// Fill new tree, new branch
		if( BuildMinTree_ ){
			cout<<entry<<endl;
			if( isData_ ){
				McFlag_=0;
			}else{
				McFlag_=1;
			}
			evtwtPu_ = evtwt_;
			evtNum_	 = maxEvents_;
			evtwtHiggsTagCorr_ = higgsTagCorr_;
			reRegistGen(GenInfo, newGenInfo); 	
			reRegistJet(allak5jets_corr, newJetInfo);	
			reRegistJet(FatJetInfo, newFatJetInfo);	
			reRegistJet(SubJetInfo, newSubJetInfo);
			for( int i=0; i<allak5jets_num; i++){
				CSVLSF_[i] = passBTagCSVL[i];
				CSVMSF_[i] = passBTagCSVM[i];
				CSVTSF_[i] = passBTagCSVT[i];
			}	
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
