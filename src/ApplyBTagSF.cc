#include "BpbH/BackgroundEstimationSkim/interface/ApplyBTagSF.h" //A1pha
#include "BpbH/BprimeTobHAnalysis/interface/SFb-pt_WITHttbar_payload_EPS13.h"
#include "BpbH/BprimeTobHAnalysis/interface/SFlightFuncs_EPS2013.h"

//ApplyBTagSF::ApplyBTagSF(JetCollection& jets, double bDisc, std::string algo_, std::string mode_) : 
ApplyBTagSF::ApplyBTagSF(JetCollection& jets, double bDisc, std::string algo, double SFbShift, double SFlShift) : 
  jets_(jets),
  bDisc_(bDisc),
  algo_(algo),
  //mode_(mode_),
  SFbShift_(SFbShift),
  SFlShift_(SFlShift) 
{
  btaggedJetsWithSF_.clear() ; 
}

ApplyBTagSF::~ApplyBTagSF() {
  btaggedJetsWithSF_.clear() ; 
  algo_.clear();
  //mode_.clear(); 
}

JetCollection ApplyBTagSF::getBtaggedJetsWithSF () { 

	int    jet_flavor;
	double jet_et;       
	double jet_phi;        
	double jet_eta;
	bool   isTagged = false;

	double BTagSF = 0.972;
	double BTageff = 0.705;
	double LightJetSF = 1.04;
	double LightJeteff = 0.01705;

	for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {

		jet_flavor = abs(ijet->GenFlavor()) ; 
		jet_et     = ijet->Et()        ; 
		jet_phi    = ijet->Eta()       ; 
		jet_eta    = ijet->Phi()       ; 
		isTagged   = ijet->CombinedSVBJetTags() > bDisc_ ; 

		int ptbin(-1) ; 
		for (int ii = 0; ii < int(sizeof(ptmin)/sizeof(ptmin[0])); ++ii) {  
			if ( jet_et >= ptmin[ii] && jet_et < ptmax[ii] ) {
				ptbin = ii ; 
			}
		}

		double errscale(1.);
		if ( jet_flavor == 4) errscale *= 2. ; 

		if ( jet_et < 20. )  jet_et = 20.;
		if ( jet_et > 800. ) jet_et = 800.;

		if (strcmp(algo_.c_str(), "CSVL") == 0) 
			BTagSF = SFb_CSVL->Eval(jet_et)*( 1 + (errscale*SFbShift_*SFb_CSVL_error[ptbin]) ) ; 
		else if (strcmp(algo_.c_str(), "CSVM") == 0) 
			BTagSF = SFb_CSVM->Eval(jet_et)*( 1 + (errscale*SFbShift_*SFb_CSVM_error[ptbin]) ) ;  
		else if (strcmp(algo_.c_str(), "CSVT") == 0) 
			BTagSF = SFb_CSVT->Eval(jet_et)*( 1 + (errscale*SFbShift_*SFb_CSVT_error[ptbin]) ) ;  
		else {
			edm::LogError("ApplyBTagSF") << " Wrong b-tagging algo_ chosen: " << algo_ << ". Choose between CSVL, CSVM, or CSVT." ; 
			return btaggedJetsWithSF_ ; 
		}

		if (strcmp(algo_.c_str(), "CSVL") == 0) {
			for (int ii = 0; ii < 4; ++ii) {
				if ( jet_eta >= SFlight_CSVL_etamin[ii] && jet_eta < SFlight_CSVL_etamax[ii] ) {
					LightJetSF = (GetSFlmean("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) + 
						SFlShift_ > 0 ? 
						SFlShift_*( (GetSFlmax("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) - 
								(GetSFlmean("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) ) : 
						SFlShift_*( (GetSFlmean("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) - 
								(GetSFlmin("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) ) ; 
				}
			}
		}
		else if (strcmp(algo_.c_str(), "CSVM") == 0) {
			for (int ii = 0; ii < 4; ++ii) {
				if ( jet_eta >= SFlight_CSVM_etamin[ii] && jet_eta < SFlight_CSVM_etamax[ii] ) {
					LightJetSF = (GetSFlmean("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) + 
						SFlShift_ > 0 ? 
						SFlShift_*( (GetSFlmax("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) - 
								(GetSFlmean("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) ) : 
						SFlShift_*( (GetSFlmean("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) - 
								(GetSFlmin("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) ) ; 
				}
			}
		}
		else if (algo_== "CSVT") {
			LightJetSF = (GetSFlmean("CSV","T",0.0, 2.4, "ABCD"))->Eval(jet_et) + 
				SFlShift_ > 0 ? 
				SFlShift_*( (GetSFlmax("CSV","T",0.0, 2.4, "ABCD"))->Eval(jet_et) - 
						(GetSFlmean("CSV","T",0.0, 2.4, "ABCD"))->Eval(jet_et) ) : 
				SFlShift_*( (GetSFlmean("CSV","T",0.0, 2.4, "ABCD"))->Eval(jet_et) - 
						(GetSFlmin("CSV","T",0.0, 2.4, "ABCD"))->Eval(jet_et) ) ; 
		}
		else {
			edm::LogError("ApplyBTagSF") << " Wrong b-tagging algo_ chosen: " << algo_ << ". Choose between CSVL, CSVM, or CSVT." ; 
			return btaggedJetsWithSF_ ; 
		}

		double phi = jet_phi;
		double sin_phi = sin(phi*1000000);
		double seed = abs(static_cast<int>(sin_phi*100000));

		//Initialize class
		BTagSFUtil* btsfutil = new BTagSFUtil( seed );

		btsfutil->modifyBTagsWithSF(isTagged, jet_flavor, BTagSF, BTageff, LightJetSF, LightJeteff);

		if (isTagged) btaggedJetsWithSF_.push_back(*ijet) ; 

		delete btsfutil ; 

	}

	return btaggedJetsWithSF_ ; 
}


// Alpha
bool ApplyBTagSF::getBtaggedSF( int ijet ) {
	int    jet_flavor;
	double jet_et;       
	double jet_phi;        
	double jet_eta;
	bool   isTagged = false;

	double BTagSF = 0.972;
	double BTageff = 0.705;
	double LightJetSF = 1.04;
	double LightJeteff = 0.01705;

	jet_flavor = abs(jets_[ijet].GenFlavor()); 
	jet_et     = jets_[ijet].Et(); 
	jet_phi    = jets_[ijet].Eta(); 
	jet_eta    = jets_[ijet].Phi(); 
	isTagged   = jets_[ijet].CombinedSVBJetTags() > bDisc_ ; 

	int ptbin(-1) ; 
	for (int ii = 0; ii < int(sizeof(ptmin)/sizeof(ptmin[0])); ++ii) {  
		if ( jet_et >= ptmin[ii] && jet_et < ptmax[ii] ) {
			ptbin = ii ; 
		}
	}

	double errscale(1.);
	if ( jet_flavor == 4) errscale *= 2. ; 

	if ( jet_et < 20. )  jet_et = 20.;
	if ( jet_et > 800. ) jet_et = 800.;

	if (strcmp(algo_.c_str(), "CSVL") == 0) 
		BTagSF = SFb_CSVL->Eval(jet_et)*( 1 + (errscale*SFbShift_*SFb_CSVL_error[ptbin]) ) ; 
	else if (strcmp(algo_.c_str(), "CSVM") == 0) 
		BTagSF = SFb_CSVM->Eval(jet_et)*( 1 + (errscale*SFbShift_*SFb_CSVM_error[ptbin]) ) ;  
	else if (strcmp(algo_.c_str(), "CSVT") == 0) 
		BTagSF = SFb_CSVT->Eval(jet_et)*( 1 + (errscale*SFbShift_*SFb_CSVT_error[ptbin]) ) ;  
	else {
		edm::LogError("ApplyBTagSF") << " Wrong b-tagging algo_ chosen: " << algo_ << ". Choose between CSVL, CSVM, or CSVT." ; 
		return 0; 
	}

	if (strcmp(algo_.c_str(), "CSVL") == 0) {
		for (int ii = 0; ii < 4; ++ii) {
			if ( jet_eta >= SFlight_CSVL_etamin[ii] && jet_eta < SFlight_CSVL_etamax[ii] ) {
				LightJetSF = (GetSFlmean("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) + 
					SFlShift_ > 0 ? 
					SFlShift_*( (GetSFlmax("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) - 
							(GetSFlmean("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) ) : 
					SFlShift_*( (GetSFlmean("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) - 
							(GetSFlmin("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) ) ; 
			}
		}
	}
	else if (strcmp(algo_.c_str(), "CSVM") == 0) {
		for (int ii = 0; ii < 4; ++ii) {
			if ( jet_eta >= SFlight_CSVM_etamin[ii] && jet_eta < SFlight_CSVM_etamax[ii] ) {
				LightJetSF = (GetSFlmean("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) + 
					SFlShift_ > 0 ? 
					SFlShift_*( (GetSFlmax("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) - 
							(GetSFlmean("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) ) : 
					SFlShift_*( (GetSFlmean("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) - 
							(GetSFlmin("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) ) ; 
			}
		}
	}
	else if (algo_== "CSVT") {
		LightJetSF = (GetSFlmean("CSV","T",0.0, 2.4, "ABCD"))->Eval(jet_et) + 
			SFlShift_ > 0 ? 
			SFlShift_*( (GetSFlmax("CSV","T",0.0, 2.4, "ABCD"))->Eval(jet_et) - 
					(GetSFlmean("CSV","T",0.0, 2.4, "ABCD"))->Eval(jet_et) ) : 
			SFlShift_*( (GetSFlmean("CSV","T",0.0, 2.4, "ABCD"))->Eval(jet_et) - 
					(GetSFlmin("CSV","T",0.0, 2.4, "ABCD"))->Eval(jet_et) ) ; 
	}
	else {
		edm::LogError("ApplyBTagSF") << " Wrong b-tagging algo_ chosen: " << algo_ << ". Choose between CSVL, CSVM, or CSVT." ; 
		return 0; 
	}

	double phi = jet_phi;
	double sin_phi = sin(phi*1000000);
	double seed = abs(static_cast<int>(sin_phi*100000));

	//Initialize class
	BTagSFUtil* btsfutil = new BTagSFUtil( seed );
	btsfutil->modifyBTagsWithSF(isTagged, jet_flavor, BTagSF, BTageff, LightJetSF, LightJeteff);
	delete btsfutil ; 

	return isTagged;

} 
