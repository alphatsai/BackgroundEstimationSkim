#ifndef BACKGROUNDESTIMATION_INTERFACE_REREGISTGEN_H
#define BACKGROUNDESTIMATION_INTERFACE_REREGISTGEN_H

#include "BpbH/BprimeTobH/interface/format.h"

void reRegistGen( GenInfoBranches& OldGen, GenInfoBranches& NewGen){
	int size=0;
	for( int i=0; i<OldGen.Size; i++){
		size++;
		NewGen.Pt[i] = OldGen.Pt[i];
		NewGen.Eta[i] = OldGen.Eta[i];
		NewGen.Phi[i] = OldGen.Phi[i];
		NewGen.Mass[i] = OldGen.Mass[i];
		NewGen.PdgID[i] = OldGen.PdgID[i];
		NewGen.Status[i] = OldGen.Status[i];
		NewGen.Mo0Pt[i] = OldGen.Mo0Pt[i];
		NewGen.Mo0Eta[i] = OldGen.Mo0Eta[i];
		NewGen.Mo0Phi[i] = OldGen.Mo0Phi[i];
		NewGen.Mo0Mass[i] = OldGen.Mo0Mass[i];
		NewGen.Mo0PdgID[i] = OldGen.Mo0PdgID[i];
		NewGen.Mo0Status[i] = OldGen.Mo0Status[i];
		NewGen.Mo1Pt[i] = OldGen.Mo1Pt[i];
		NewGen.Mo1Eta[i] = OldGen.Mo1Eta[i];
		NewGen.Mo1Phi[i] = OldGen.Mo1Phi[i];
		NewGen.Mo1Mass[i] = OldGen.Mo1Mass[i];
		NewGen.Mo1PdgID[i] = OldGen.Mo1PdgID[i];
		NewGen.Mo1Status[i] = OldGen.Mo1Status[i];
		NewGen.Da0Pt[i] = OldGen.Da0Pt[i];
		//std::cout<<OldGen.Da0Pt[i]<<" "<<NewGen.Da0Pt[i]<<std::endl;
		NewGen.Da0Eta[i] = OldGen.Da0Eta[i];
		NewGen.Da0Phi[i] = OldGen.Da0Phi[i];
		NewGen.Da0Mass[i] = OldGen.Da0Mass[i];
		NewGen.Da0PdgID[i] = OldGen.Da0PdgID[i];
		NewGen.Da0Status[i] = OldGen.Da0Status[i];
		NewGen.Da1Pt[i] = OldGen.Da1Pt[i];
		NewGen.Da1Eta[i] = OldGen.Da1Eta[i];
		NewGen.Da1Phi[i] = OldGen.Da1Phi[i];
		NewGen.Da1Mass[i] = OldGen.Da1Mass[i];
		NewGen.Da1PdgID[i] = OldGen.Da1PdgID[i];
		NewGen.Da1Status[i] = OldGen.Da1Status[i];
		NewGen.nMo[i] = OldGen.nMo[i];
		NewGen.nDa[i] = OldGen.nDa[i];
		NewGen.Mo1[i] = OldGen.Mo1[i];
		NewGen.Mo2[i] = OldGen.Mo2[i];
		NewGen.Da1[i] = OldGen.Da1[i];
		NewGen.Da2[i] = OldGen.Da2[i];
		NewGen.cQuark_pT[i] = OldGen.cQuark_pT[i];
		NewGen.cQuark_eta[i] = OldGen.cQuark_eta[i];
		NewGen.cQuark_phi[i] = OldGen.cQuark_phi[i];
		NewGen.cQuark_pdgID[i] = OldGen.cQuark_pdgID[i];
		NewGen.cQuark_status[i] = OldGen.cQuark_status[i];
		NewGen.cQuark_fromGSP[i] = OldGen.cQuark_fromGSP[i];
	}
	NewGen.Size=size;
	NewGen.ncQuarks=OldGen.ncQuarks;
}
#endif 
