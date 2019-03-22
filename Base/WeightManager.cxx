#ifndef WEIGHTMANAGER_CXX
#define WEIGHTMANAGER_CXX

#include <sstream>
#include <map>
#include <set>
#include "WeightManager.h"
#include "CLHEP/Random/Random.h"

namespace evwgh {

  WeightManager::WeightManager(const std::string name)
    : _name(name)
  {
    _configured = false;
  }

  const std::string& WeightManager::Name() const
  { return _name; }

  size_t WeightManager::Configure(fhicl::ParameterSet const & p) 
  {
   
    // Get list of weight functions
    std::vector<std::string> rw_func = p.get<std::vector<std::string>>("weight_functions");

    // Loop over all the functions and register them
    for (auto ifunc=rw_func.begin(); ifunc!=rw_func.end(); ifunc++) 
    {
      fhicl::ParameterSet const &ps_func = p.get<fhicl::ParameterSet> (*ifunc);
      std::string func_type = ps_func.get<std::string>("type");
  
      WeightCalc* wcalc=WeightCalcFactory::Create(func_type+"WeightCalc");
      if ( wcalc == NULL ) {
	std::cerr<< "Function " << *ifunc << " requested in fcl file has not been registered!" << std::endl;
	throw std::exception();
      } 
      if ( fWeightCalcMap.find(*ifunc)!= fWeightCalcMap.end() ) {
	std::cerr<< "Function " << *ifunc << " has been requested multiple times in fcl file!" << std::endl;
	throw std::exception();
      }
  
      //mf::LogInfo("") << "Configuring weight calculator " << *ifunc;
  
      // Create random engine for each rw function (name=*ifunc) (and seed it with random_seed set in the fcl)

      wcalc->SetName(*ifunc);
      wcalc->Configure(p);
      Weight_t* winfo=new Weight_t();
      winfo->fWeightCalcType=func_type;
      winfo->fWeightCalc=wcalc;
      winfo->fNmultisims=ps_func.get<int>("number_of_multisims");
  
      std::pair<std::string, Weight_t*> pwc(*ifunc,winfo);
      fWeightCalcMap.insert(pwc);
    }

    _configured = true;

    return fWeightCalcMap.size();
  }


 
  // 
  // CORE FUNCTION
  //
  MCEventWeight WeightManager::Run(bsim::Dk2Nu & e, const int inu)
  {     

    if (!_configured) {
      std::cerr<< "WeightManager was not configured!"<<std::endl;
      throw std::exception();
    }
    //
    // Loop over all functions ang calculate weights
    //
    MCEventWeight mcwgh;
    for (auto it = fWeightCalcMap.begin() ;it != fWeightCalcMap.end(); it++) {

      auto const & weights = it->second->GetWeight(e);
      if(weights.size() == 0){
        std::vector<double> empty;
        std::pair<std::string, std::vector <double> > p("empty",empty);
        mcwgh.fWeight.insert(p);
      } 
      else{
        std::pair<std::string, std::vector<double> > 
          p(it->first+"_"+it->second->fWeightCalcType,
            weights[inu]);
        mcwgh.fWeight.insert(p);
      }
    }

    return mcwgh;
  }



  void WeightManager::PrintConfig() {
    
    return; 
  }

}

#endif
