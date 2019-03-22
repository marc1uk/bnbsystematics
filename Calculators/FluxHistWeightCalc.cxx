#include "Base/WeightCalcCreator.h"
#include "Base/WeightCalc.h"

#include <iostream>

//#include "nutools/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandGaussQ.h"


#include "TFile.h"
#include "TH1F.h"

namespace evwgh {
  class FluxHistWeightCalc : public WeightCalc
  {
  public:
    FluxHistWeightCalc();
    void Configure(fhicl::ParameterSet const& pset);
    std::vector<std::vector<double> > GetWeight(bsim::Dk2Nu & e);
    
  private:    
    CLHEP::RandGaussQ *fGaussRandom;
    std::vector<double> fWeightArray;
    int fNmultisims;
    std::string fMode;
    std::string fGenieModuleLabel;

    //         pi+-,k+-,k0,mu+- 
    //         |  numu, numubar, nue, nuebar 
    //         |  |   50MeV bins
    //         |  |   |
    double fCV[4][4][200];
    double fRW[4][4][200];
    
    DECLARE_WEIGHTCALC(FluxHistWeightCalc)
  };
  FluxHistWeightCalc::FluxHistWeightCalc()
  {
  }

  void FluxHistWeightCalc::Configure(fhicl::ParameterSet const& p)
  {    
    //global config
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");

    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    //calc config
    fNmultisims = pset.get<int>("number_of_multisims");
    fMode       = pset.get<std::string>("mode");		
    std::string dataInput1 = pset.get< std::string >("cv_hist_file");
    std::string dataInput2 = pset.get< std::string >("rw_hist_file");

    cet::search_path sp("FW_SEARCH_PATH");
    std::string cvfile = sp.find_file(dataInput1);
    std::string rwfile = sp.find_file(dataInput2);
    
    std::string ptype[] = {"pi", "k", "k0", "mu"};
    std::string ntype[] = {"numu", "numubar", "nue", "nuebar"};

    TFile fcv(Form("%s",cvfile.c_str()));
    TFile frw(Form("%s",rwfile.c_str()));
    for (int iptyp=0;iptyp<4;iptyp++) {
      for (int intyp=0;intyp<4;intyp++) {
	for (int ibin=0;ibin<200;ibin++) {
	  fCV[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (fcv.Get(Form("h_%s_%s",ptype[iptyp].c_str(),ntype[intyp].c_str()))))->GetBinContent(ibin+1);
	  fRW[iptyp][intyp][ibin]=(dynamic_cast<TH1F*> (frw.Get(Form("h_%s_%s",ptype[iptyp].c_str(),ntype[intyp].c_str()))))->GetBinContent(ibin+1);
	}
      }
    }
    fcv.Close();
    frw.Close();

    CLHEP::HepRandomEngine* rng=new CLHEP::HepJamesRandom();
    fGaussRandom = new CLHEP::RandGaussQ(rng,0,1);
    fWeightArray.resize(fNmultisims);

    if (fMode.find("multisim") != std::string::npos )
      for (int i=0;i<fNmultisims;i++) fWeightArray[i]=fGaussRandom->shoot(rng,0,1.);
    else
      for (int i=0;i<fNmultisims;i++) fWeightArray[i]=1.;
  }

  std::vector<std::vector<double> > FluxHistWeightCalc::GetWeight(bsim::Dk2Nu & e)
  {
    //calculate weight(s) here 
    std::vector<std::vector<double> > weight;

    //within art it is possible to have multiple nu interactions in an event, but here it is always 1
    //leaving vector of vectors to minimize the diff in code
    weight.resize(1);
    for (unsigned int inu=0;inu<1;inu++) {
      weight[inu].resize(fNmultisims);
     
      int ptype=-9999;
      int ntype=-9999;
      int bin=-9999;
      
      if ( e.decay.ptype==211 || e.decay.ptype==-211 ) ptype = 0;
      else if ( e.decay.ptype==321 || e.decay.ptype==-321 ) ptype = 1;
      else if ( e.decay.ptype==130 ) ptype = 2;
      else if ( e.decay.ptype==13 || e.decay.ptype==-13 ) ptype = 3;
      else {
	throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ptype "<<e.decay.ptype<< std::endl;
      }
      
      if ( e.decay.ntype==14 ) ntype=0;
      else if ( e.decay.ntype==-14 ) ntype=1;
      else if ( e.decay.ntype==12 ) ntype=2;
      else if ( e.decay.ntype==-12 ) ntype=3;
      else {
	throw cet::exception(__FUNCTION__) << GetName()<<"::Unknown ptype "<<e.decay.ptype<< std::endl;
      }

      //nuray index should match the location of input reweighting histograms
      double enu=e.nuray[1].E;
      bin=enu/0.05;     
      for (int i=0;i<fNmultisims;i++) {
	double test = 1-(1-fRW[ptype][ntype][bin]/fCV[ptype][ntype][bin])*fWeightArray[i];
	
	// Guards against inifinite weights
	if(std::isfinite(test)){ weight[inu][i] = test;}
	else{weight[inu][i] = 1;}
      }
    }
    return weight;
  }
  REGISTER_WEIGHTCALC(FluxHistWeightCalc)
}
