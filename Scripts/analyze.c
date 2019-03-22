TMatrixD* getMatrix(TH1D* cv, std::vector<TH1D*> ms);
TH1D* getErrHist(std::string name, TH1D* cv, TMatrixD* mat);

void analyze(std::string fname,std::string ntype="numu")
{
  const int NSYS=11;
  const int NMS=1000;
  std::string sys[]={"piplus_PrimaryHadronSWCentralSplineVariation","kplus_PrimaryHadronFeynmanScaling","kzero_PrimaryHadronSanfordWang","nucleoninexsec_FluxUnisim","nucleonqexsec_FluxUnisim","nucleontotxsec_FluxUnisim","piontotxsec_FluxUnisim","pionqexsec_FluxUnisim","pioninexsec_FluxUnisim","expskin_FluxUnisim","horncurrent_FluxUnisim"};
  int hcol[]={kBlue, kRed, kMagenta, kGreen+2};
  TH1D* hMS; 
  TFile fin(fname.c_str());
  TH1D* hCV=(TH1D*)fin.Get(ntype.c_str());
  hCV->SetDirectory(0);
  std::vector<TH1D*> hVecMS;//total MS
  
  std::map<std::string,std::vector<TH1D*> > hSysMapMS;
  for (int isys=0;isys<NSYS;isys++) {
    std::vector<TH1D* > vec;
    hSysMapMS[sys[isys]]=vec;
  }
  for (int ims=0;ims<NMS;ims++) {
    hMS=(TH1D*)fin.Get(Form("%s_%i",ntype.c_str(),ims));
    hMS->SetDirectory(0);
    hVecMS.push_back(hMS);
    for (int isys=0;isys<NSYS;isys++) {
      hMS=(TH1D*)fin.Get(Form("%s_%s_%i",ntype.c_str(),sys[isys].c_str(),ims));
      hMS->SetDirectory(0);
      hSysMapMS[sys[isys]].push_back(hMS);
    }
  }
  fin.Close();
  
  std::map<std::string,TMatrixD* > hErrMap;
  for (auto sm: hSysMapMS) {
    cout <<"Build matrix for "<<sm.first<<endl;
    hErrMap[sm.first]=getMatrix(hCV, sm.second);
  }
  
  //fill histograms other than uni
  std::map<std::string,TH1D* > hErrHist;
  TMatrixD* unimat=new TMatrixD(hCV->GetNbinsX(),hCV->GetNbinsX());
  for (auto sm: hErrMap) {
    if (sm.first.find("FluxUnisim")==std::string::npos) {
      hErrHist[sm.first]=getErrHist(sm.first+"_err",hCV,hErrMap[sm.first]);
    } else {
      (*unimat)+=(*sm.second);
    }
  }
 
  std::cout<<"Build total error matrix"<<std::endl;
  TMatrixD* toterr=getMatrix(hCV,hVecMS);
  TH1D* hTotErr=getErrHist("hTotErr",hCV,toterr);
  for (int i=0;i<hCV->GetNbinsX();i++) 
    hCV->SetBinError(i+1,sqrt((*toterr)(i,i)));
  hErrHist["BeamUni"]=getErrHist("hUniErr",hCV,unimat);

  TCanvas* c1=new TCanvas("c1","",550,500);
  hCV->Draw("histe");

  TCanvas* c2=new TCanvas("c2","",550,500);
  hTotErr->SetLineColor(kBlack);
  hTotErr->SetLineWidth(2);
  hTotErr->SetStats(0);
  hTotErr->Draw("hist");
  TLegend* leg=new TLegend(0.5,0.4,0.88,0.88);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hTotErr,"Total","L");
  int ii=0;
  for (auto hm:hErrHist) {
    hm.second->SetLineColor(hcol[ii%4]);
    hm.second->Draw("samehist");
    leg->AddEntry(hm.second,hm.first.c_str(),"L");
    ii++;
  }
  leg->Draw();
}
TMatrixD* getMatrix(TH1D* cv, std::vector<TH1D*> ms)
{
  int nbin=cv->GetNbinsX();
  int nms=ms.size();
  
  double mat[nbin][nbin];
  for (int i=0;i<nbin;i++) {
    for (int j=0;j<nbin;j++) {
      mat[i][j]=0;
      for (int ims=0;ims<nms;ims++) {
	if (cv->GetBinContent(i+1)>0&&cv->GetBinContent(j+1)>0)
	  mat[i][j]+=(cv->GetBinContent(i+1)-ms[ims]->GetBinContent(i+1))*(cv->GetBinContent(j+1)-ms[ims]->GetBinContent(j+1))/(double(nms));
      }
    }
  }

  TMatrixD* newmat=new TMatrixD(nbin,nbin);
  for (int i=0;i<nbin;i++) {
    for (int j=0;j<nbin;j++) {
      (*newmat)(i,j)=mat[i][j];
    }
  }

  return newmat;
}

TH1D* getErrHist(std::string name, TH1D* cv, TMatrixD* mat)
{
  TH1D* h=(TH1D*)cv->Clone(name.c_str());
  h->Reset();
  for (int i=0;i<cv->GetNbinsX();i++) {
    if (cv->GetBinContent(i+1)>0)
      h->SetBinContent(i+1,sqrt((*mat)(i,i))/cv->GetBinContent(i+1));
  }

  return h;
}
