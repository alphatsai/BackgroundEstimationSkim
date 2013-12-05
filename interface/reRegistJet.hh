#ifndef BACKGROUNDESTIMATION_INTERFACE_REREGISTJET_H
#define BACKGROUNDESTIMATION_INTERFACE_REREGISTJET_H

#include "BpbH/BprimeTobH/interface/format.h"

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
#endif 
