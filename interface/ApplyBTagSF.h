#ifndef BACKGROUNDESTIMATIONSKIM_INTERFACE_APPLYBTAGSF_H
#define BACKGROUNDESTIMATIONSKIM_INTERFACE_APPLYBTAGSF_H

#include "BpbH/BprimeTobHAnalysis/interface/BTagSFUtil.h"
#include "BpbH/BprimeTobH/interface/JetCollection.h" 

#include <string>

class ApplyBTagSF {

  public:
  //ApplyBTagSF(JetCollection&, double, std::string, std::string) ;
    ApplyBTagSF(JetCollection&, double, std::string, double, double) ; 
    ~ApplyBTagSF() ; 

    JetCollection getBtaggedJetsWithSF () ; 
    bool getBtaggedSF( int ijet ) ; 

  private:

    JetCollection jets_ ; 
    double bDisc_ ; 
    JetCollection btaggedJetsWithSF_ ; 

    std::string algo_ ;
    double SFbShift_ ; 
    double SFlShift_ ; 
    //std::string mode_ ; 

};

#endif 
