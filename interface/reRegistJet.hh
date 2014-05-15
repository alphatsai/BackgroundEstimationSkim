#ifndef BACKGROUNDESTIMATION_INTERFACE_REREGISTJET_H
#define BACKGROUNDESTIMATION_INTERFACE_REREGISTJET_H

#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/JetCollection.h"
#include "BpbH/BprimeTobH/interface/Jet.h"

namespace{
	void reRegistJet( JetInfoBranches& OldJet, JetInfoBranches& NewJet){
		int size=0;
		for( int i=0; i<OldJet.Size; i++){
			size++;
			NewJet.Index[i] = OldJet.Index[i];
			NewJet.NTracks[i] = OldJet.NTracks[i];
			NewJet.Et[i] = OldJet.Et[i];
			NewJet.Pt[i] = OldJet.Pt[i];
			NewJet.Unc[i] = OldJet.Unc[i];
			NewJet.Eta[i] = OldJet.Eta[i];
			NewJet.Phi[i] = OldJet.Phi[i];
			NewJet.Energy[i] = OldJet.Energy[i];
			NewJet.Px[i] = OldJet.Px[i];
			NewJet.Py[i] = OldJet.Py[i];
			NewJet.Pz[i] = OldJet.Pz[i];
			NewJet.Mass[i] = OldJet.Mass[i];
			NewJet.Area[i] = OldJet.Area[i];
			NewJet.EtPruned[i] = OldJet.EtPruned[i];
			NewJet.PtPruned[i] = OldJet.PtPruned[i];
			NewJet.UncPruned[i] = OldJet.UncPruned[i];
			NewJet.EtaPruned[i] = OldJet.EtaPruned[i];
			NewJet.PhiPruned[i] = OldJet.PhiPruned[i];
			NewJet.EnergyPruned[i] = OldJet.EnergyPruned[i];
			NewJet.PxPruned[i] = OldJet.PxPruned[i];
			NewJet.PyPruned[i] = OldJet.PyPruned[i];
			NewJet.PzPruned[i] = OldJet.PzPruned[i];
			NewJet.MassPruned[i] = OldJet.MassPruned[i];
			NewJet.AreaPruned[i] = OldJet.AreaPruned[i];
			NewJet.tau1[i] = OldJet.tau1[i];
			NewJet.tau2[i] = OldJet.tau2[i];
			NewJet.tau3[i] = OldJet.tau3[i];
			NewJet.JetIDLOOSE[i] = OldJet.JetIDLOOSE[i];
			NewJet.JetIDTIGHT[i] = OldJet.JetIDTIGHT[i];
			NewJet.JetCharge[i] = OldJet.JetCharge[i];
			NewJet.QGTagsMLP[i] = OldJet.QGTagsMLP[i];
			NewJet.QGTagsLikelihood[i] = OldJet.QGTagsLikelihood[i];
			NewJet.NConstituents[i] = OldJet.NConstituents[i];
			NewJet.NCH[i] = OldJet.NCH[i];
			NewJet.CEF[i] = OldJet.CEF[i];
			NewJet.NHF[i] = OldJet.NHF[i];
			NewJet.NEF[i] = OldJet.NEF[i];
			NewJet.CHF[i] = OldJet.CHF[i];
			NewJet.PtCorrRaw[i] = OldJet.PtCorrRaw[i];
			NewJet.PtCorrL2[i] = OldJet.PtCorrL2[i];
			NewJet.PtCorrL3[i] = OldJet.PtCorrL3[i];
			NewJet.PtCorrL7g[i] = OldJet.PtCorrL7g[i];
			NewJet.PtCorrL7uds[i] = OldJet.PtCorrL7uds[i];
			NewJet.PtCorrL7c[i] = OldJet.PtCorrL7c[i];
			NewJet.PtCorrL7b[i] = OldJet.PtCorrL7b[i];
			NewJet.JetBProbBJetTags[i] = OldJet.JetBProbBJetTags[i];
			NewJet.JetProbBJetTags[i] = OldJet.JetProbBJetTags[i];
			NewJet.TrackCountHiPurBJetTags[i] = OldJet.TrackCountHiPurBJetTags[i];
			NewJet.CombinedSVBJetTags[i] = OldJet.CombinedSVBJetTags[i];
			NewJet.CombinedSVMVABJetTags[i] = OldJet.CombinedSVMVABJetTags[i];
			NewJet.SoftElecByIP3dBJetTags[i] = OldJet.SoftElecByIP3dBJetTags[i];
			NewJet.SoftElecByPtBJetTags[i] = OldJet.SoftElecByPtBJetTags[i];
			NewJet.SoftMuonBJetTags[i] = OldJet.SoftMuonBJetTags[i];
			NewJet.SoftMuonByIP3dBJetTags[i] = OldJet.SoftMuonByIP3dBJetTags[i];
			NewJet.SoftMuonByPtBJetTags[i] = OldJet.SoftMuonByPtBJetTags[i];
			NewJet.DoubleSVHighEffBJetTags[i] = OldJet.DoubleSVHighEffBJetTags[i];
			NewJet.GenJetPt[i] = OldJet.GenJetPt[i];
			NewJet.GenJetEta[i] = OldJet.GenJetEta[i];
			NewJet.GenJetPhi[i] = OldJet.GenJetPhi[i];
			NewJet.GenPt[i] = OldJet.GenPt[i];
			NewJet.GenEta[i] = OldJet.GenEta[i];
			NewJet.GenPhi[i] = OldJet.GenPhi[i];
			NewJet.GenPdgID[i] = OldJet.GenPdgID[i];
			NewJet.GenFlavor[i] = OldJet.GenFlavor[i];
			NewJet.GenMCTag[i] = OldJet.GenMCTag[i];
			NewJet.Jet_FatJetIdx[i] = OldJet.Jet_FatJetIdx[i];
			NewJet.Jet_SubJet1Idx[i] = OldJet.Jet_SubJet1Idx[i];
			NewJet.Jet_SubJet2Idx[i] = OldJet.Jet_SubJet2Idx[i];	
		}
		NewJet.Size=size;
	}

	void reRegistJet( JetCollection& OldJet, JetInfoBranches& NewJet){
		int size=0;
		for( JetCollection::const_iterator ijet = OldJet.begin(); ijet != OldJet.end(); ++ijet) {
			//	NewJet.Index[size] = ijet->Index();
			NewJet.NTracks[size] = ijet->NTracks();
			NewJet.Et[size] = ijet->Et();
			NewJet.Pt[size] = ijet->Pt();
			NewJet.Unc[size] = ijet->Unc();
			NewJet.Eta[size] = ijet->Eta();
			NewJet.Phi[size] = ijet->Phi();
			NewJet.Energy[size] = ijet->Energy();
			NewJet.Px[size] = ijet->Px();
			NewJet.Py[size] = ijet->Py();
			NewJet.Pz[size] = ijet->Pz();
			NewJet.Mass[size] = ijet->Mass();
			NewJet.Area[size] = ijet->Area();
			NewJet.EtPruned[size] = ijet->EtPruned();
			NewJet.PtPruned[size] = ijet->PtPruned();
			NewJet.UncPruned[size] = ijet->UncPruned();
			NewJet.EtaPruned[size] = ijet->EtaPruned();
			NewJet.PhiPruned[size] = ijet->PhiPruned();
			NewJet.EnergyPruned[size] = ijet->EnergyPruned();
			NewJet.PxPruned[size] = ijet->PxPruned();
			NewJet.PyPruned[size] = ijet->PyPruned();
			NewJet.PzPruned[size] = ijet->PzPruned();
			NewJet.MassPruned[size] = ijet->MassPruned();
			NewJet.AreaPruned[size] = ijet->AreaPruned();
			NewJet.tau1[size] = ijet->tau1();
			NewJet.tau2[size] = ijet->tau2();
			NewJet.tau3[size] = ijet->tau3();
			NewJet.JetIDLOOSE[size] = ijet->JetIDLOOSE();
			NewJet.JetIDTIGHT[size] = ijet->JetIDTIGHT();
			NewJet.JetCharge[size] = ijet->JetCharge();
			NewJet.QGTagsMLP[size] = ijet->QGTagsMLP();
			NewJet.QGTagsLikelihood[size] = ijet->QGTagsLikelihood();
			NewJet.NConstituents[size] = ijet->NConstituents();
			NewJet.NCH[size] = ijet->NCH();
			NewJet.CEF[size] = ijet->CEF();
			NewJet.NHF[size] = ijet->NHF();
			NewJet.NEF[size] = ijet->NEF();
			NewJet.CHF[size] = ijet->CHF();
			NewJet.PtCorrRaw[size] = ijet->PtCorrRaw();
			NewJet.PtCorrL2[size] = ijet->PtCorrL2();
			NewJet.PtCorrL3[size] = ijet->PtCorrL3();
			NewJet.PtCorrL7g[size] = ijet->PtCorrL7g();
			NewJet.PtCorrL7uds[size] = ijet->PtCorrL7uds();
			NewJet.PtCorrL7c[size] = ijet->PtCorrL7c();
			NewJet.PtCorrL7b[size] = ijet->PtCorrL7b();
			NewJet.JetBProbBJetTags[size] = ijet->JetBProbBJetTags();
			NewJet.JetProbBJetTags[size] = ijet->JetProbBJetTags();
			NewJet.TrackCountHiPurBJetTags[size] = ijet->TrackCountHiPurBJetTags();
			NewJet.CombinedSVBJetTags[size] = ijet->CombinedSVBJetTags();
			NewJet.CombinedSVMVABJetTags[size] = ijet->CombinedSVMVABJetTags();
			NewJet.SoftElecByIP3dBJetTags[size] = ijet->SoftElecByIP3dBJetTags();
			NewJet.SoftElecByPtBJetTags[size] = ijet->SoftElecByPtBJetTags();
			NewJet.SoftMuonBJetTags[size] = ijet->SoftMuonBJetTags();
			NewJet.SoftMuonByIP3dBJetTags[size] = ijet->SoftMuonByIP3dBJetTags();
			NewJet.SoftMuonByPtBJetTags[size] = ijet->SoftMuonByPtBJetTags();
			NewJet.DoubleSVHighEffBJetTags[size] = ijet->DoubleSVHighEffBJetTags();
			NewJet.GenJetPt[size] = ijet->GenJetPt();
			NewJet.GenJetEta[size] = ijet->GenJetEta();
			NewJet.GenJetPhi[size] = ijet->GenJetPhi();
			NewJet.GenPt[size] = ijet->GenPt();
			NewJet.GenEta[size] = ijet->GenEta();
			NewJet.GenPhi[size] = ijet->GenPhi();
			NewJet.GenPdgID[size] = ijet->GenPdgID();
			NewJet.GenFlavor[size] = ijet->GenFlavor();
			NewJet.GenMCTag[size] = ijet->GenMCTag();
			NewJet.Jet_FatJetIdx[size] = ijet->Jet_FatJetIdx();
			NewJet.Jet_SubJet1Idx[size] = ijet->Jet_SubJet1Idx();
			NewJet.Jet_SubJet2Idx[size] = ijet->Jet_SubJet2Idx();	
			size++;
		}
		NewJet.Size=size;
	}
}
#endif 
