//g++ -o rwgh rwgh.cxxcd -L lib/ -lBNBSysBase -lBNBSysCalc -I${FHICLCPP_INC} -I${CETLIB_EXCEPT_INC} -I${CETLIB_INC} -I${CLHEP_INC} -I${ROOTSYS}/include -I./ `root-config --glibs` -L${FHICLCPP_LIB} -lfhiclcpp -L${CLHEP_LIB_DIR} -lCLHEP -I${BOOST_INC} -L${BOOST_LIB} -lprogram_options -lfilesystem
#include <iostream>
#include <glob.h>
#include <utility>

#include "Base/WeightManager.h"
#include "Base/MCEventWeight.h"

#include "cetlib/filepath_maker.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/readWeightLocations.h"
#include "dk2nu/tree/calcLocationWeights.h"


namespace fs=boost::filesystem ;
using namespace boost::program_options;

int to_pdg(int id);

int main(int ac, char* av[])
{
  std::string searchpath;
  std::string outputfn;
  std::string fhiclfile;
  options_description opt("Options");
  opt.add_options()
    ("help,h", "Print help message")
    ("input,i",value<std::string>(&searchpath)->required(),"Path pattern for input files. Put it in quotes or escape *.")
    ("fhicl,f",value<std::string>(&fhiclfile)->required(),"Configuration fhicl file.")
    ("output,o",value<std::string>(&outputfn)->default_value("hist.root"),"Output file name.");

  variables_map vm;
  
  try {
    store(parse_command_line(ac,av,opt),vm);
    notify(vm);
    if (vm.count("help")) {
      std::cerr<<opt<<std::endl;
      return 1;
    } 
  } catch (error& e) {
    std::cerr << e.what()<<std::endl<<std::endl;
    std::cerr << opt <<std::endl;
    return 1;
  }
  
  fs::path fhiclfp(fhiclfile);
  cet::filepath_first_absolute_or_lookup_with_dot fpm(fhiclfp.parent_path().string());
  fhicl::ParameterSet pm;
  fhicl::make_ParameterSet(fhiclfp.filename().string(),fpm,pm);

  glob_t glob_result;
  std::cout<<"Searching "<<searchpath<<std::endl;
  glob(searchpath.c_str(),GLOB_TILDE,NULL,&glob_result);
  std::vector<std::string> filelist;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    filelist.push_back(std::string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  if (filelist.size()>0) 
    std::cout <<"Found "<<filelist.size()<<" files"<<std::endl; 
  else {
    std::cerr<<"Failed to find any input files"<<std::endl;
    return -1;
  }

  evwgh::WeightManager wm;  
  wm.Configure(pm.get<fhicl::ParameterSet>("microboone_eventweight"));

  TChain* bnb=new TChain("h101");
  for (auto ifile : filelist) {
    bnb->Add(ifile.c_str());
  }
  
  float beamwgt;
  int ntp;
  int npart;
  int id[20];
  float ini_pos[20][3];
  float ini_mom[20][3];
  float ini_eng[20];
  float ini_t[20];
  float fin_mom[20][3];
  float fin_pol[20][3];
  
  bnb->SetBranchAddress("beamwgt",&beamwgt);
  bnb->SetBranchAddress("ntp",&ntp);
  bnb->SetBranchAddress("npart",&npart);
  bnb->SetBranchAddress("id",id);
  bnb->SetBranchAddress("ini_pos",ini_pos);
  bnb->SetBranchAddress("ini_mom",ini_mom);
  bnb->SetBranchAddress("ini_eng",ini_eng);
  bnb->SetBranchAddress("ini_t",ini_t);
  bnb->SetBranchAddress("fin_mom",fin_mom);
  bnb->SetBranchAddress("fin_pol",fin_pol);

  bsim::Dk2Nu* dk2nu= new bsim::Dk2Nu;
  
  std::string nu[]={"nue","nuebar","numu","numubar"};
  TH1D* hFluxCV[4];
  for (int i=0;i<4;i++) {
    hFluxCV[i]=new TH1D(nu[i].c_str(),";Energy (GeV);",200,0,10);
  }
  TH1D* hFluxMS[4][1000];
  for (int i=0;i<4;i++) {
    for (int j=0;j<1000;j++) 
      hFluxMS[i][j]=new TH1D(Form("%s_%i",nu[i].c_str(),j),";Energy (GeV);",200,0,10);
  }
  std::map<std::string, std::vector<std::vector<TH1D*> > > hSysMap;
  Long64_t ientry=0;
  while (bnb->GetEntry(ientry)) {
    if (ientry%100==0) 
      std::cout<<"On ientry "<<ientry<<std::endl;
    ientry++;
    //id(0)     = 4 (neutrino)
    //id(1)     = 5,6,8,9,10,11,12 (nu parent)
    //..
    //id(npart-1) = 14 (primary proton)
    //fill only relevant parts of dk2nu (used in calculators)
    //this is not actually tgtexit, it first particle produced in p+Be (and used with that assumption in Calculators)
    //should switch to using ancestor tree in Calculators
    dk2nu->tgtexit.tvx=ini_pos[npart-2][0];
    dk2nu->tgtexit.tvy=ini_pos[npart-2][1];
    dk2nu->tgtexit.tvz=ini_pos[npart-2][2];
    dk2nu->tgtexit.tpx=ini_mom[npart-2][0];
    dk2nu->tgtexit.tpy=ini_mom[npart-2][1];
    dk2nu->tgtexit.tpz=ini_mom[npart-2][2];
    dk2nu->tgtexit.tptype=to_pdg(id[npart-2]);
    dk2nu->tgtexit.tgen=npart;

    dk2nu->decay.vx       = ini_pos[0][0] ;
    dk2nu->decay.vy       = ini_pos[0][1] ;
    dk2nu->decay.vz       = ini_pos[0][2] ;
    dk2nu->decay.pdpx     = fin_mom[1][0] ;
    dk2nu->decay.pdpy     = fin_mom[1][1] ;
    dk2nu->decay.pdpz     = fin_mom[1][2] ;
    
    dk2nu->decay.ppdxdz   = ini_mom[1][0]/ini_mom[1][2] ;
    dk2nu->decay.ppdydz   = ini_mom[1][1]/ini_mom[1][2] ;
    dk2nu->decay.pppz     = ini_mom[1][2] ;
    dk2nu->decay.ppenergy = ini_eng[1];
    dk2nu->decay.ptype    = to_pdg(id[1]);

    //if parent is a muon fill muon parent info
    if ( id[1] == 5 ||
	 id[1] == 6) {
      dk2nu->decay.mupare  = ini_eng[2];
      dk2nu->decay.muparpx = fin_mom[2][0];
      dk2nu->decay.muparpy = fin_mom[2][1];
      dk2nu->decay.muparpz = fin_mom[2][2];
    } else {
      dk2nu->decay.mupare  = -9999.;
      dk2nu->decay.muparpx = -9999.;
      dk2nu->decay.muparpy = -9999.;
      dk2nu->decay.muparpz = -9999.;
    }
    
    //NEUTRINO info:
    if ( ntp == 1 )
      dk2nu->decay.ntype = 12;
    else if ( ntp == 2 )
      dk2nu->decay.ntype = -12;
    else if ( ntp == 3 )
      dk2nu->decay.ntype = 14;
    else if ( ntp == 4 )
      dk2nu->decay.ntype = -14;
    else
      std::cerr<<"Neutrino type not recognized!!! ntp = "<<ntp<<" ; ientry = "<<ientry<<std::endl;
    
    //Get the neutrino energy in the parent decay cm
    double parent_mass=sqrt(dk2nu->decay.ppenergy*dk2nu->decay.ppenergy-
			    dk2nu->decay.pppz*dk2nu->decay.pppz*(dk2nu->decay.ppdxdz*dk2nu->decay.ppdxdz +
							       dk2nu->decay.ppdydz*dk2nu->decay.ppdydz +
							       1.));
    double parent_energy = sqrt(dk2nu->decay.pdpx*dk2nu->decay.pdpx +
				dk2nu->decay.pdpy*dk2nu->decay.pdpy +
				dk2nu->decay.pdpz*dk2nu->decay.pdpz + parent_mass*parent_mass);
    
    double gamma         = parent_energy / parent_mass;
    double beta[3];
    beta[0] = dk2nu->decay.pdpx/parent_energy;
    beta[1] = dk2nu->decay.pdpy/parent_energy;
    beta[2] = dk2nu->decay.pdpz/parent_energy;
    double partial = gamma * ( beta[0] * ini_mom[0][0] + beta[1] * ini_mom[0][1] + beta[2]*ini_mom[0][2] );
    dk2nu->decay.necm = gamma * ini_eng[0] - partial;
    
    dk2nu->decay.nimpwt = beamwgt;
    
    //Fill Ndecay (check parent type, neutrino type and if it is a 2 or 3 body decay)
    if (id[1] == 10 && ntp == 1) 
      dk2nu->decay.ndecay = 1;
    else if (id[1] == 10 && ntp == 2) 
      dk2nu->decay.ndecay = 2;
    else if (id[1] == 10 && ntp == 3) 
      dk2nu->decay.ndecay = 3;
    else if (id[1] == 10 && ntp == 4) 
      dk2nu->decay.ndecay = 4;
    else if (id[1] == 11 && ntp == 3) {
      //check if it is a two or three body decay
      if (fabs((parent_mass*parent_mass-0.105658389*0.105658389)/(2.*parent_mass)-dk2nu->decay.necm)/dk2nu->decay.necm <= 0.001)
	//two body decay (numu + mu+)
	dk2nu->decay.ndecay = 5;
      else
	//three body decay (numu + pi0 + mu+)
	dk2nu->decay.ndecay = 7;
    } else if (id[1] == 11 && ntp == 1) 
      dk2nu->decay.ndecay = 6;
    else if (id[1] == 12 && ntp == 4) {
      if (fabs((parent_mass*parent_mass-0.105658389*0.105658389)/(2.*parent_mass)-dk2nu->decay.necm)/dk2nu->decay.necm <= 0.001)
	//two body decay (numu + mu+)
	dk2nu->decay.ndecay = 8;
      else
	//three body decay (numu + pi0 + mu+)
	dk2nu->decay.ndecay = 10;
    } else if (id[1] == 12 && ntp == 2) 
      dk2nu->decay.ndecay = 9;
    else if (id[1] == 5 ) 
      dk2nu->decay.ndecay = 11;
    else if (id[1] == 6 ) 
      dk2nu->decay.ndecay = 12;
    else if (id[1] == 8 ) 
      dk2nu->decay.ndecay = 13;
    else if (id[1] == 9 ) 
      dk2nu->decay.ndecay = 14;
    
    //random neutrino:
    bsim::NuRay anuray_rnd(ini_mom[0][0], ini_mom[0][1], ini_mom[0][2], ini_eng[0], 1);
    dk2nu->nuray.push_back(anuray_rnd);
   
    TVector3 xyzDk(dk2nu->decay.vx,dk2nu->decay.vy,dk2nu->decay.vz);  // origin of decay
    double enu_xy,wgt_xy;
     
    //neutrino at MiniBooNE (energy info used by histogram reweighting calculators
    TVector3 xyzMiniBooNE(0.,   189.614,   54134.);  // << miniboone is 540m downstream of the target
    bsim::calcEnuWgt(dk2nu,xyzMiniBooNE,enu_xy,wgt_xy);
    TVector3 p3_mb = enu_xy * (xyzMiniBooNE - xyzDk).Unit();
    bsim::NuRay anuray_mb(p3_mb.x(), p3_mb.y(), p3_mb.z(), enu_xy, wgt_xy);
    dk2nu->nuray.push_back(anuray_mb);
    
//    //neutrino energy at MicroBooNE:
//    TVector3 xyzMicroBooNE(0., 0., 46336.3525);
//    bsim::calcEnuWgt(dk2nu,xyzMicroBooNE,enu_xy,wgt_xy);
//    TVector3 p3_ub = enu_xy * (xyzMicroBooNE - xyzDk).Unit();
//    bsim::NuRay anuray_ub(p3_ub.x(), p3_ub.y(), p3_ub.z(), enu_xy, wgt_xy);
//    dk2nu->nuray.push_back(anuray_ub);
    
    //neutrino energy at ANNIE:
    TVector3 xyzANNIE(0., 0., 10000.0);
    // "Based on accurate survey data, the distance between the center of the beryllium target
    // and the center of the SciBar detector is taken to be 99.9 m, with the SciBooNE detector
    // located on beam axis within a tolerance of a few centimeters." (DOI: 10.1103/PhysRevD.78.112004)
    bsim::calcEnuWgt(dk2nu,xyzANNIE,enu_xy,wgt_xy);
    TVector3 p3_an = enu_xy * (xyzANNIE - xyzDk).Unit();
    bsim::NuRay anuray_an(p3_an.x(), p3_an.y(), p3_an.z(), enu_xy, wgt_xy);
    dk2nu->nuray.push_back(anuray_an);

    std::vector<double> totwgh(1000,1.);
    evwgh::MCEventWeight wght=wm.Run(*dk2nu,0);
    double rwgh=1;
    for (auto w:wght.fWeight) 
      if (w.second.size()==1) {
	rwgh*=w.second[0]; //assume this is reweight
      }
    for (auto w:wght.fWeight) {
      //check if it is multisim
      if (w.second.size()>1) {
	if (hSysMap.find(w.first)==hSysMap.end()) {
	  std::vector<std::vector<TH1D*> > hsys(4);
	  for (int ii=0;ii<4;ii++) {
	    hsys[ii].resize(1000);
	    for (int jj=0;jj<1000;jj++) {
	      hsys[ii][jj]=new TH1D(Form("%s_%s_%i",nu[ii].c_str(),w.first.c_str(),jj),"",200,0,10);
	    }

	  }
	  hSysMap.insert(std::make_pair(w.first,hsys));
	}
	for (unsigned int ims=0;ims<w.second.size();ims++) {
	  totwgh[ims]*=w.second[ims];	  
	  hSysMap[w.first][ntp-1][ims]->Fill(enu_xy,beamwgt*wgt_xy/3.14159*w.second[ims]*rwgh);
	}
      }
    }

    hFluxCV[ntp-1]->Fill(enu_xy,beamwgt*wgt_xy/3.14159*rwgh);
    for (int ims=0;ims<1000;ims++) {
      hFluxMS[ntp-1][ims]->Fill(enu_xy,beamwgt*wgt_xy/3.14159*totwgh[ims]*rwgh);
    }
    dk2nu->nuray.clear();
  }
  
  TFile fout(outputfn.c_str(),"RECREATE");
  for (int inu=0;inu<4;inu++) {
    for (int ims=0;ims<1000;ims++) {
      hFluxMS[inu][ims]->Write();
      for (auto sys : hSysMap) {
	sys.second[inu][ims]->Write();
      }
    }
  }
  for (int inu=0;inu<4;inu++) 
    hFluxCV[inu]->Write();
  fout.Close();
  
  std::cout<<ientry<<std::endl;
}

int to_pdg(int id)
{
  int pdg=0;
  if (id==8)
    pdg    = 211; 
  else if (id==9)
    pdg    = -211; 
  else if (id==5)
    pdg    = -13; 
  else if (id==6)
    pdg    = 13; 
  else if (id==11)
    pdg    = 321; 
  else if (id==12)
    pdg    = -321; 
  else if (id==10)
    pdg    = 130; 
  else if (id==14)
    pdg    = 2212;
  else if (id==15)
    pdg    = -2212;
  else if (id==13)
    pdg    = 2112;
  else if (id==25)
    pdg    = -2112;

  else if (id==18)
    pdg    = 3122;
  
  else if (id==19)
    pdg    = 3222;
  else if (id==20)
    pdg    = 3212;
  else if (id==21)
    pdg    = 3112;

  else if (id==16)
    pdg    = 310;
  
  else
    std::cerr<<"Parent type not recognized!!! ptype = "<<id<<std::endl;

  return pdg;

}
