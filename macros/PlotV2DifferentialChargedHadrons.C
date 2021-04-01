#include "PlotV2LYZ.C"
#include "PlotV2EtaSubEventPlane.C"
#include "PlotV2ThreeEtaSubEventPlane.C"
#include "PlotV2FHCalEventPlane.C"
#include "PlotV2ScalarProduct.C"
#include "PlotV2QCumulant.C"
// #include "DrawTGraph.C"

vector<TGraphErrors*>* PlotV2DifferentialChargedHadrons(TString inputFirstRunFileName = "FirstRun_UrQMD_7.7_Reco.root", TString inputSecondRunFileName = "SecondRun_UrQMD_7.7_Reco.root")
{
  bool ETASUBEVENTPLANE_1 = 0;
  bool ETASUBEVENTPLANE_2 = 1;
  bool THREEETASUBEVENTPLANE_1 = 0;
  bool THREEETASUBEVENTPLANE_2 = 1;
  bool FHCALEVENTPLANE_1 = 0;
  bool FHCALEVENTPLANE_2 = 1;
  bool LYZ_SUM_1 = 0;
  bool LYZ_SUM_2 = 0;
  bool LYZ_SUM_PRODUCT_1 = 0;
  bool LYZ_SUM_PRODUCT_2 = 0;
  bool SCALARPRODUCT_1 = 0;
  bool SCALARPRODUCT_2 = 1;
  bool QCUMULANT = 0;
  bool HIGHORDERQCUMULANT = 0;
  bool LYZEP = 0;
  Double_t maxpt = 3.6;    // max pt for differential flow
  Double_t minpt = 0.;     // min pt for differential flow
  Double_t maxptRF = 3.;   // max pt for reference flow
  Double_t minptRF = 0.2;  // min pt for reference flow
  Double_t eta_cut = 1.5;  // pseudorapidity acceptance window for flow measurements
  Double_t eta_gap = 0.05; // +-0.05, eta-gap between 2 eta sub-event of two-particle cumulants method with eta-gap
  const int ratioToMethod = 4;
  const double J1rootJ0 = 0.519147;
  double X[npt];
  for (int ipt=0; ipt<npt; ipt++)
  {
    X[ipt] = (pTBin[ipt] + pTBin[ipt+1]) / 2.;
  }
  const double errX[npt] = {0.};
  bool bUseProduct = 1;
  Int_t nmethod = 9;
  TString title[]={"v_{2}{#Psi_{2,TPC}}","v_{2}^{SP}{Q_{2,TPC}}","v_{2}{2}","v_{2}{4}","v_{2}{2,#eta-gap}","v_{2}{LYZ, Sum}","v_{2}{LYZ, Prod.}","v_{2}{#Psi_{1,FHCal}}","v_{2}{#Psi_{2,3-sub}}"};
  const int markerStyle[] = {24,22,27,21,20,25,28,26,23};
  const float markerSize = 1.3;
  TGraphErrors *graph[ncent][nmethod];
  TFile *firun1 = new TFile(inputFirstRunFileName.Data(),"read");
  TFile *firun2 = new TFile(inputSecondRunFileName.Data(),"read");
  auto *prV2TPCEP3D = (TProfile3D*) firun2->Get("prV2EtaSubEventPlane");
  auto *prV2SP3D = (TProfile3D*) firun2->Get("prV2ScalarProduct");
  auto *prV2FHCalEP3D = (TProfile3D*) firun2->Get("prV2FHCalEventPlane");
  auto *prV2ThreeSubEP3D = (TProfile3D*) firun2->Get("prV2ThreeEtaSub");
  for (int i = 0; i < ncent; i++)
  {
    TProfile *prV2TPCEP = PlotV2EtaSubDifferentialVersusPt(prV2TPCEP3D,bin_cent[i],bin_cent[i+1]-1,eta_cut);
    TProfile *prV2SP = PlotV2ScalarProductDifferentialVersusPt(prV2SP3D,bin_cent[i],bin_cent[i+1]-1,eta_cut);
    TProfile *prV2FHCalEP = PlotV2FHCalEPDifferentialVersusPt(prV2FHCalEP3D,bin_cent[i],bin_cent[i+1]-1,eta_cut);
    TProfile *prV2ThreeEtaSub = PlotV2ThreeEtaSubEPDifferentialVersusPt(prV2ThreeSubEP3D,bin_cent[i],bin_cent[i+1]-1,eta_cut);
    graph[i][0] = Converter(prV2TPCEP);
    graph[i][1] = Converter(prV2SP);
    graph[i][7] = Converter(prV2FHCalEP);
    graph[i][8] = Converter(prV2ThreeEtaSub);
  }
  
  // QCumulant
  const int npid_plotQC = 12;
  TProfile *pCorrelator2EtaGap = (TProfile*)firun1->Get("pCorrelator2EtaGap");
  TProfile *pCorrelator2 = (TProfile*)firun1->Get("pCorrelator2");
  TProfile *pCorrelator4 = (TProfile*)firun1->Get("pCorrelator4");
  TProfile2D *pReducedCorrelator2EtaGap[npid_plotQC]; // <<2'>> (with eta-gap)
  TProfile2D *pReducedCorrelator2[npid_plotQC]; // <<2'>>
  TProfile2D *pReducedCorrelator4[npid_plotQC]; // <<4'>>
  TProfile2D *pCov22RedEtaGap[npid_plotQC];
  TProfile *pCov24 = (TProfile*)firun1->Get("pCov24");
  TProfile2D *pCov22Red[npid_plotQC];
  TProfile2D *pCov24Red[npid_plotQC];
  TProfile2D *pCov42Red[npid_plotQC];
  TProfile2D *pCov44Red[npid_plotQC];
  TProfile2D *pCov2Red4Red[npid_plotQC];

  for (int i=0; i<npid_plotQC-4; i++)
  {
    pReducedCorrelator2EtaGap[i] = (TProfile2D*)firun1->Get(Form("pReducedCorrelator2EtaGap_pid%i",i));
    pReducedCorrelator2[i]   = (TProfile2D*) firun1->Get(Form("pReducedCorrelator2_pid%i",i));
    pReducedCorrelator4[i]   = (TProfile2D*) firun1->Get(Form("pReducedCorrelator4_pid%i",i));
    pCov22RedEtaGap[i]       = (TProfile2D*) firun1->Get(Form("pCov22RedEtaGap_pid%i",i));
    pCov22Red[i]             = (TProfile2D*) firun1->Get(Form("pCov22Red_pid%i",i));
    pCov24Red[i]             = (TProfile2D*) firun1->Get(Form("pCov24Red_pid%i",i));
    pCov42Red[i]             = (TProfile2D*) firun1->Get(Form("pCov42Red_pid%i",i));
    pCov44Red[i]             = (TProfile2D*) firun1->Get(Form("pCov44Red_pid%i",i));
    pCov2Red4Red[i]          = (TProfile2D*) firun1->Get(Form("pCov2Red4Red_pid%i",i));
  }

  TProfile *pReducedCorrelator2EtaGap_cent[npid_plotQC][ncent];
  TProfile *pReducedCorrelator2_cent[npid_plotQC][ncent];
  TProfile *pReducedCorrelator4_cent[npid_plotQC][ncent];
  TProfile *pCov22RedEtaGap_cent[npid_plotQC][ncent];
  TProfile *pCov22Red_cent[npid_plotQC][ncent];
  TProfile *pCov24Red_cent[npid_plotQC][ncent];
  TProfile *pCov42Red_cent[npid_plotQC][ncent];
  TProfile *pCov44Red_cent[npid_plotQC][ncent];
  TProfile *pCov2Red4Red_cent[npid_plotQC][ncent];
  for (int i = 0; i < npid_plotQC-4; i++)
  {
    for (int c = 0; c < ncent; c++)
    {
      pReducedCorrelator2EtaGap_cent[i][c] = (TProfile*)pReducedCorrelator2EtaGap[i]->ProfileX(Form("%s_cent%i",pReducedCorrelator2EtaGap[i]->GetName(),c), c+1, c+1);
      pReducedCorrelator2_cent[i][c] = (TProfile*)pReducedCorrelator2[i]->ProfileX(Form("%s_cent%i",pReducedCorrelator2[i]->GetName(),c), c+1, c+1);
      pReducedCorrelator4_cent[i][c] = (TProfile*)pReducedCorrelator4[i]->ProfileX(Form("%s_cent%i",pReducedCorrelator4[i]->GetName(),c), c+1, c+1);
      pCov22RedEtaGap_cent[i][c] = (TProfile*)pCov22RedEtaGap[i]->ProfileX(Form("%s_cent%i",pCov22RedEtaGap[i]->GetName(),c), c+1, c+1);
      pCov22Red_cent[i][c] = (TProfile*)pCov22Red[i]->ProfileX(Form("%s_cent%i",pCov22Red[i]->GetName(),c), c+1, c+1);
      pCov24Red_cent[i][c] = (TProfile*)pCov24Red[i]->ProfileX(Form("%s_cent%i",pCov24Red[i]->GetName(),c), c+1, c+1); 
      pCov42Red_cent[i][c] = (TProfile*)pCov42Red[i]->ProfileX(Form("%s_cent%i",pCov42Red[i]->GetName(),c), c+1, c+1); 
      pCov44Red_cent[i][c] = (TProfile*)pCov44Red[i]->ProfileX(Form("%s_cent%i",pCov44Red[i]->GetName(),c), c+1, c+1); 
      pCov2Red4Red_cent[i][c] = (TProfile*)pCov2Red4Red[i]->ProfileX(Form("%s_cent%i",pCov2Red4Red[i]->GetName(),c), c+1, c+1); 
    }
  }
  for (int i=8;i<npid_plotQC;i++)
  {
    for (int c = 0; c < ncent; c++)
    {
      pReducedCorrelator2EtaGap_cent[i][c] = (TProfile*)  pReducedCorrelator2EtaGap_cent[i-8][c] -> Clone();
      pReducedCorrelator2_cent[i][c] = (TProfile*)  pReducedCorrelator2_cent[i-8][c] -> Clone();
      pReducedCorrelator4_cent[i][c] = (TProfile*)  pReducedCorrelator4_cent[i-8][c] -> Clone();
      pCov22RedEtaGap_cent[i][c] = (TProfile*)  pCov22RedEtaGap_cent[i-8][c] -> Clone();
      pCov22Red_cent[i][c] = (TProfile*)  pCov22Red_cent[i-8][c] -> Clone();
      pCov24Red_cent[i][c] = (TProfile*)  pCov24Red_cent[i-8][c] -> Clone();
      pCov42Red_cent[i][c] = (TProfile*)  pCov42Red_cent[i-8][c] -> Clone();
      pCov44Red_cent[i][c] = (TProfile*)  pCov44Red_cent[i-8][c] -> Clone();
      pCov2Red4Red_cent[i][c] = (TProfile*)  pCov2Red4Red_cent[i-8][c] -> Clone();

      pReducedCorrelator2EtaGap_cent[i][c] -> Add(pReducedCorrelator2EtaGap_cent[i-4][c]);
      pReducedCorrelator2_cent[i][c] -> Add(pReducedCorrelator2_cent[i-4][c]);
      pReducedCorrelator4_cent[i][c] -> Add(pReducedCorrelator4_cent[i-4][c]);
      pCov22RedEtaGap_cent[i][c] -> Add(pCov22RedEtaGap_cent[i-4][c]);
      pCov22Red_cent[i][c] -> Add(pCov22Red_cent[i-4][c]);
      pCov24Red_cent[i][c] -> Add(pCov24Red_cent[i-4][c]);
      pCov42Red_cent[i][c] -> Add(pCov42Red_cent[i-4][c]);
      pCov44Red_cent[i][c] -> Add(pCov44Red_cent[i-4][c]);
      pCov2Red4Red_cent[i][c] -> Add(pCov2Red4Red_cent[i-4][c]);
    }
  }
  double v2Dif[3][ncent][npid_plotQC][npt], v2eDif[3][ncent][npid_plotQC][npt];
  for (int icent=0; icent<ncent; icent++){ // loop over centrality classes
    // 2QC
    term cor2 = term(pCorrelator2,icent);
    double v22 = sqrt(cor2.mVal);
    double ev22 = sqrt(1./(4.*cor2.mVal)*cor2.mMSE);
    // 4QC
    term cor4 = term(pCorrelator4,icent);
    double cov24 = Covariance(pCov24,pCorrelator2,pCorrelator4,icent,icent,icent);
    double v24 = pow(2*pow(cor2.mVal,2)-cor4.mVal,0.25);
    double ev24 = sqrt( 1./pow(v24,6)*(cor2.mVal*cor2.mVal*cor2.mMSE+1./16*cor4.mMSE-0.5*cor2.mVal*cov24) );
    // 2QC Gapped
    term cor2Gap = term(pCorrelator2EtaGap,icent);
    double v22Gap = sqrt(cor2Gap.mVal);
    double ev22Gap = sqrt(1./(4.*cor2Gap.mVal)*cor2Gap.mMSE);
    for (int id=0;id<npid_plotQC;id++){ // Differential flow calculation
      for(int ipt=0; ipt<npt; ipt++){ // loop for all pT bin

        // v22
        term cor2red = term(pReducedCorrelator2_cent[id][icent],ipt);
        double v22Dif = cor2red.mVal/v22;
        double cov22prime = Covariance(pCov22Red_cent[id][icent],pCorrelator2,pReducedCorrelator2_cent[id][icent],ipt,icent,ipt);
        double ev22Dif = sqrt(0.25*pow(cor2.mVal,-3)*(pow(cor2red.mVal,2)*cor2.mMSE
                            + 4*pow(cor2.mVal,2)*cor2red.mMSE - 4*cor2.mVal*cor2red.mVal*cov22prime));
        v2Dif[0][icent][id][ipt] = v22Dif;
        v2eDif[0][icent][id][ipt] = ev22Dif;
        
        // v24
        term cor4red = term(pReducedCorrelator4_cent[id][icent],ipt);
        double cov24prime = Covariance(pCov24Red_cent[id][icent],pCorrelator2,pReducedCorrelator4_cent[id][icent],ipt,icent,ipt);
        double cov42prime = Covariance(pCov42Red_cent[id][icent],pCorrelator4,pReducedCorrelator2_cent[id][icent],ipt,icent,ipt);
        double cov44prime = Covariance(pCov44Red_cent[id][icent],pCorrelator4,pReducedCorrelator4_cent[id][icent],ipt,icent,ipt);
        double cov2prime4prime = Covariance(pCov2Red4Red_cent[id][icent],pReducedCorrelator2_cent[id][icent],pReducedCorrelator4_cent[id][icent],ipt,ipt,ipt);
        double v24Dif = (2.*cor2.mVal*cor2red.mVal-cor4red.mVal)*pow(v24,-3);
        double ev24Dif = sqrt( pow(v24,-14)
            * (pow(2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal,2.)
            * cor2.mMSE
            + 9./16*pow(2.*cor2.mVal*cor2red.mVal-cor4red.mVal,2.)*cor4.mMSE
            + 4*pow(cor2.mVal,2)*pow(v24,8)*cor2red.mMSE
            + pow(v24,8)*cor4red.mMSE
            - 1.5*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
            * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
            * cov24
            - 4*cor2.mVal*pow(v24,4)
            * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
            * cov22prime
            + 2*pow(v24,4)
            * (2*cor2.mVal*cor2.mVal*cor2red.mVal-3*cor2.mVal*cor4red.mVal+2*cor4.mVal*cor2red.mVal)
            * cov24prime
            + 3*cor2.mVal*pow(v24,4)*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
            * cov42prime
            - 1.5*pow(v24,4)*(2*cor2.mVal*cor2red.mVal-cor4red.mVal)
            * cov44prime
            - 4*cor2.mVal*pow(v24,8)*cov2prime4prime));
        v2Dif[1][icent][id][ipt] = v24Dif;
        v2eDif[1][icent][id][ipt] = ev24Dif;
        // v22 Gapped
        term cor2redGap = term(pReducedCorrelator2EtaGap_cent[id][icent],ipt);
        double v22DifGap = cor2redGap.mVal/v22Gap;
        double cov22primeGap = Covariance(pCov22RedEtaGap_cent[id][icent],pCorrelator2EtaGap,pReducedCorrelator2EtaGap_cent[id][icent],ipt,icent,ipt);
        double ev22DifGap = sqrt(0.25*pow(cor2Gap.mVal,-3)*(pow(cor2redGap.mVal,2)*cor2Gap.mMSE
                            + 4*pow(cor2Gap.mVal,2)*cor2redGap.mMSE - 4*cor2Gap.mVal*cor2redGap.mVal*cov22primeGap));
        v2Dif[2][icent][id][ipt] = v22DifGap;
        v2eDif[2][icent][id][ipt] = ev22DifGap;

      } // end of loop for all pT bin
    } // end of loop over PID
    graph[icent][2] = new TGraphErrors(npt, X, v2Dif[0][icent][8], errX, v2eDif[0][icent][8]); // v2{2}
    graph[icent][3] = new TGraphErrors(npt, X, v2Dif[1][icent][8], errX, v2eDif[1][icent][8]); // v2{4}
    graph[icent][4] = new TGraphErrors(npt, X, v2Dif[2][icent][8], errX, v2eDif[2][icent][8]); // v2{2,eta-gap}

  } // end of loop over centrality classes

  
  // LYZ
  if (LYZ_SUM_1 || LYZ_SUM_2 || LYZ_SUM_PRODUCT_1 || LYZ_SUM_PRODUCT_2)
  {
    double theta[thetabins];
    for (int thetabin = 0; thetabin < thetabins; ++thetabin)
    {
      theta[thetabin] = thetabin * TMath::Pi() / (2.0 * thetabins);
    }
    TProfile *prReGthetaSum[ncent][thetabins];
    TProfile *prImGthetaSum[ncent][thetabins];
    TProfile *prReGthetaProduct[ncent][thetabins];
    TProfile *prImGthetaProduct[ncent][thetabins];
    TProfile *prRefMult = (TProfile*) firun1->Get("prRefMult");
    TProfile *prQ2x = (TProfile*) firun1->Get("prQ2x");
    TProfile *prQ2y = (TProfile*) firun1->Get("prQ2y");
    TProfile *prQ2ModSq = (TProfile*) firun1->Get("prQ2ModSq");
    // TFile *fi_fixedGFS = new TFile("FirstRun_fixedGFSum.root","read");
    for (int i = 0; i < ncent; ++i)
    {
      for (int j = 0; j < thetabins; ++j)
      {
        prReGthetaSum[i][j] = (TProfile*) firun1->Get(Form("prReGthetaSum_mult%d_theta%d", i, j));
        prImGthetaSum[i][j] = (TProfile*) firun1->Get(Form("prImGthetaSum_mult%d_theta%d", i, j));
        if (bUseProduct)
        {
          prReGthetaProduct[i][j] = (TProfile*) firun1->Get(Form("prReGthetaProduct_mult%d_theta%d", i, j));
          prImGthetaProduct[i][j] = (TProfile*) firun1->Get(Form("prImGthetaProduct_mult%d_theta%d", i, j));
        }
      }
    }
    // Differential flow
    TProfile *prReDenom[thetabins];
    TProfile *prImDenom[thetabins];
    TProfile *prReNumer[thetabins][ncent];
    TProfile *prImNumer[thetabins][ncent];

    TProfile *prReDenomPro[thetabins];
    TProfile *prImDenomPro[thetabins];
    TProfile *prReNumerPro[thetabins][ncent];
    TProfile *prImNumerPro[thetabins][ncent];

    for (int i = 0; i < thetabins; i++)
    {
      prReDenom[i] = (TProfile*) firun2->Get(Form("prReDenom_theta%i",i));
      prImDenom[i] = (TProfile*) firun2->Get(Form("prImDenom_theta%i",i));

      for (int j = 0; j < ncent; j++)
      {
        prReNumer[i][j] = (TProfile*) firun2->Get(Form("prReNumer_theta%i_cent%i", i, j));
        prImNumer[i][j] = (TProfile*) firun2->Get(Form("prImNumer_theta%i_cent%i", i, j));
      }
      if (bUseProduct)
      {
        prReDenomPro[i] = (TProfile*) firun2->Get(Form("prReDenomPro_theta%i", i));
        prImDenomPro[i] = (TProfile*) firun2->Get(Form("prImDenomPro_theta%i", i));

        for (int j = 0; j < ncent; j++)
        {
          prReNumerPro[i][j] = (TProfile*) firun2->Get(Form("prReNumerPro_theta%i_cent%i", i, j));
          prImNumerPro[i][j] = (TProfile*) firun2->Get(Form("prImNumerPro_theta%i_cent%i", i, j));
        }    
      }   
    }
    TProfile *prMultPOI[ncent];
    for (int ic = 0; ic < ncent; ic++)
    {
      prMultPOI[ic] = (TProfile*) firun2->Get(Form("prMultPOI_cent%i",ic));
    }

    double dChi2[ncent]={0.}, dChi2Pro[ncent] = {0.};
    double v2LYZInt[ncent]={0.}, v2eLYZInt[ncent]={0.};
    double v2LYZIntPro[ncent]={0.}, v2eLYZIntPro[ncent]={0.};
    double dVtheta[ncent][thetabins] = {{0.}}, dVthetaPro[ncent][thetabins] = {{0.}};

    // Sum
    for (int ic = 0; ic < ncent; ic++)
    {
      double refmult = prRefMult->GetBinContent(ic+1);
      double v2int = 0., v2eint = 0., v2theta[thetabins] = {0.};
      int thetacount = 0;

      for (int it = 0; it < thetabins; it++)
      {
        TH1F *hGtheta = FillHistGtheta(prReGthetaSum[ic][it], prImGthetaSum[ic][it]);
        float r0theta = GetR0(hGtheta);
        if (r0theta!=0) 
        {
          v2int += rootJ0 / r0theta;
          v2theta[it] = rootJ0 / r0theta;
          thetacount++;
        }
      }
      if (thetacount!=0) v2int /= (float)thetacount; // refmult
      else {v2int = 0.;}
      
      float modQ2sqmean = prQ2ModSq->GetBinContent(ic+1);
      float Q2xmean = prQ2x->GetBinContent(ic+1);
      float Q2ymean = prQ2y->GetBinContent(ic+1);
      float chi2 = v2int/sqrt(modQ2sqmean-Q2xmean*Q2xmean-Q2ymean*Q2ymean-pow(v2int,2));
      
      float temp=0.;
      for(int it=0; it<thetabins; it++)
      {
        double arg = theta[it];
        temp+=exp(sqr(rootJ0/chi2)*cos(arg)/2.)*
          BesselJ0(2.*rootJ0*sin(arg/2.))+
          exp(-sqr(rootJ0/chi2)*cos(arg)/2.)*
          BesselJ0(2.*rootJ0*cos(arg/2.));
      }
      float neve = prRefMult->GetBinEntries(ic+1);
      float err2mean = v2int*sqrt(temp/2./neve/thetabins)/rootJ0/J1rootJ0;

      for (int it = 0; it < thetabins; it++) dVtheta[ic][it] = v2theta[it];
      v2LYZInt[ic] = v2int;
      dChi2[ic] = chi2;
      v2eLYZInt[ic] = err2mean;

    } // end of V2RP calculation
    

    // Product
    if (bUseProduct) {
      for (int ic = 0; ic < ncent; ic++)
      {
        double refmult = prRefMult->GetBinContent(ic+1);
        double v2int = 0., v2eint = 0., v2theta[thetabins] = {0.};
        int thetacount = 0;

        for (int it = 0; it < thetabins; it++)
        {
          TH1F *hGtheta = FillHistGtheta(prReGthetaProduct[ic][it], prImGthetaProduct[ic][it]);
          float r0theta = GetR0(hGtheta);
          if (r0theta!=0) 
          {
            v2int += rootJ0 / r0theta;
            v2theta[it] = rootJ0 / r0theta;
            thetacount++;
          }
        }
        if (thetacount!=0) v2int /= (float)thetacount*refmult;
        else {v2int = 0.;}
        
        float modQ2sqmean = prQ2ModSq->GetBinContent(ic+1);
        float Q2xmean = prQ2x->GetBinContent(ic+1);
        float Q2ymean = prQ2y->GetBinContent(ic+1);
        float chi2 = v2int/sqrt(modQ2sqmean-Q2xmean*Q2xmean-Q2ymean*Q2ymean-pow(v2int,2));
        
        float temp=0.;
        for(int it=0; it<thetabins; it++)
        {
          double arg = theta[it];
          temp+=exp(sqr(rootJ0/chi2)*cos(arg)/2.)*
            BesselJ0(2.*rootJ0*sin(arg/2.))+
            exp(-sqr(rootJ0/chi2)*cos(arg)/2.)*
            BesselJ0(2.*rootJ0*cos(arg/2.));
        }
        float neve = prRefMult->GetBinEntries(ic+1);
        float err2mean = v2int*sqrt(temp/2./neve/thetabins)/rootJ0/J1rootJ0;

        for (int it = 0; it < thetabins; it++) dVthetaPro[ic][it] = v2theta[it];
        v2LYZIntPro[ic] = v2int;
        dChi2Pro[ic] = chi2;
        v2eLYZIntPro[ic] = err2mean;

      } // end of V2RP calculation    
    } 

    // Differential v2 LYZ Sum
    double v2diff[ncent][npt]={0.};
    double v2diffe[ncent][npt]={0.};
    double v2diff_check[ncent][npt][thetabins]={{{0.}}};
    for (int ic = 0; ic < ncent; ic++)
    {
      double thetacount = 0.;
      for (int thetabin = 0; thetabin < thetabins; thetabin++)
      {
        if (dVtheta[ic][thetabin] != 0)
        {
          thetacount++;
          double re = prReDenom[thetabin]->GetBinContent(ic+1);
          double im = prImDenom[thetabin]->GetBinContent(ic+1);
          TComplex cDenominator = TComplex(re, im);
          if (cDenominator.Rho()==0) {
            cerr<<"WARNING: modulus of cDenominator is zero" << "Cent:" << bin_cent[ic] << "-"<< bin_cent[ic+1] <<"%, Theta=" << theta[thetabin] <<endl;
          }
          
          for (int ipt = 0; ipt < npt; ipt++)
          {
            double reNum = prReNumer[thetabin][ic]->GetBinContent(ipt+1);
            double imNum = prImNumer[thetabin][ic]->GetBinContent(ipt+1);
            TComplex cNumeratorPOI = TComplex(reNum, imNum);
            if (cDenominator.Rho()!=0) {
              double reRatio = (cNumeratorPOI/cDenominator).Re();
              double dVetaPOI = reRatio * dVtheta[ic][thetabin];
              v2diff[ic][ipt] += dVetaPOI;
              v2diff_check[ic][ipt][thetabin] = dVetaPOI;
            }
          }
        }
      }
      double neve = prReDenom[0]->GetBinEntries(ic+1);
      /* Computation of statistical error bars on the average estimates */
      double temp = 0.;
      double arg[thetabins];
      for(int k1=0; k1<thetabins; k1++)
      {
        // float arg=((float) it)*TMath::Pi()/(thetabins-1.);
        arg[k1] = theta[k1];

        /* Loop over the theta angles, to compute the statistical error */
        temp += (exp(sqr(rootJ0/dChi2[ic])*cos(arg[k1])/2.)*
        BesselJ0(2.*rootJ0*sin(arg[k1]/2.)) -
        exp(-sqr(rootJ0/dChi2[ic])*cos(arg[k1])/2.)*
        BesselJ0(2.*rootJ0*cos(arg[k1]/2.)))*cos(arg[k1]);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {    
        v2diff[ic][ipt] /= thetacount;
        double rpmult = prMultPOI[ic]->GetBinContent(ipt+1);
        v2diffe[ic][ipt] = sqrt(temp/rpmult/neve/thetabins)/2./J1rootJ0;
      }
    } // end of Diff LYZ Sum

    // Differential v2 LYZ Pro
    double v2diffPro[ncent][npt]={0.};
    double v2diffePro[ncent][npt]={0.};
    double v2diff_checkPro[ncent][npt][thetabins]={{{0.}}};
    if (bUseProduct){
    for (int ic = 0; ic < ncent; ic++)
    {
      double thetacount = 0.;
      for (int thetabin = 0; thetabin < thetabins; thetabin++)
      {
        if (dVthetaPro[ic][thetabin] != 0)
        {
          thetacount++;
          double re = prReDenomPro[thetabin]->GetBinContent(ic+1);
          double im = prImDenomPro[thetabin]->GetBinContent(ic+1);
          TComplex cDenominator = TComplex(re, im);
          if (cDenominator.Rho()==0) {
            cerr<<"WARNING: modulus of cDenominator (Product) is zero" << "Cent:" << bin_cent[ic] << "-"<< bin_cent[ic+1] <<"%, Theta=" << theta[thetabin] <<endl;
          }
          
          for (int ipt = 0; ipt < npt; ipt++)
          {
            double reNum = prReNumerPro[thetabin][ic]->GetBinContent(ipt+1);
            double imNum = prImNumerPro[thetabin][ic]->GetBinContent(ipt+1);
            TComplex cNumeratorPOI = TComplex(reNum, imNum);
            if (cDenominator.Rho()!=0) {
              double reRatio = (cNumeratorPOI/cDenominator).Re();
              double dVetaPOI = reRatio * dVthetaPro[ic][thetabin];
              v2diffPro[ic][ipt] += dVetaPOI;
              v2diff_checkPro[ic][ipt][thetabin] = dVetaPOI;
            }
          }
        }
      }
      double neve = prReDenomPro[0]->GetBinEntries(ic+1);
      /* Computation of statistical error bars on the average estimates */
      double temp = 0.;
      double arg[thetabins];
      for(int k1=0; k1<thetabins; k1++)
      {
        // float arg=((float) it)*TMath::Pi()/(thetabins-1.);
        arg[k1] = theta[k1];

        /* Loop over the theta angles, to compute the statistical error */
        temp += (exp(sqr(rootJ0/dChi2Pro[ic])*cos(arg[k1])/2.)*
        BesselJ0(2.*rootJ0*sin(arg[k1]/2.)) -
        exp(-sqr(rootJ0/dChi2Pro[ic])*cos(arg[k1])/2.)*
        BesselJ0(2.*rootJ0*cos(arg[k1]/2.)))*cos(arg[k1]);
      }
      for (int ipt = 0; ipt < npt; ipt++)
      {    
        v2diffPro[ic][ipt] /= thetacount;
        double rpmult = prMultPOI[ic]->GetBinContent(ipt+1);
        v2diffePro[ic][ipt] = sqrt(temp/rpmult/neve/thetabins)/2./J1rootJ0; // 
      }
    }} // end of Diff LYZ Product

    for (int ic = 0; ic < ncent; ic++)
    {
      graph[ic][5] = new TGraphErrors(npt, X, v2diff[ic], errX, v2diffe[ic]);
      graph[ic][6] = new TGraphErrors(npt, X, v2diffPro[ic], errX, v2diffePro[ic]);
    }
  }
  else
  {
    for (int ic = 0; ic < ncent; ic++)
    {
      graph[ic][5] = new TGraphErrors(npt, X, errX, errX, errX);
      graph[ic][6] = new TGraphErrors(npt, X, errX, errX, errX);
    }
  }
  
  for (int ic = 0; ic < ncent; ic++)
  {
    for (int i=0; i<nmethod; i++)
    {
    // graph[ic][i]->RemovePoint(0);
    graph[ic][i]->SetTitle(title[i].Data());
    graph[ic][i]->SetMarkerStyle(markerStyle[i]);
    graph[ic][i]->SetMarkerSize(markerSize);
    graph[ic][i]->GetXaxis()->SetTitle("p_{T}, GeV/c");
    graph[ic][i]->GetYaxis()->SetTitle("v_{2}");
    graph[ic][i]->SetDrawOption("P PLC PMC");
    }
  }

  vector<TGraphErrors*> vGr[ncent];
  for (int ic = 0; ic < ncent; ic++)
  {
    vGr[ic].push_back(graph[ic][ratioToMethod]);
    // for (int i=0; i<nmethod; i++)
    // {
    //   if (i==ratioToMethod) continue;
    //   vGr[ic].push_back(graph[ic][i]);
    // }
    vGr[ic].push_back(graph[ic][0]);
    vGr[ic].push_back(graph[ic][1]);
    // vGr[ic].push_back(graph[ic][3]);
    vGr[ic].push_back(graph[ic][7]);
    vGr[ic].push_back(graph[ic][8]);

  TCanvas *can = (TCanvas*)DrawTGraph(vGr[ic],Form("%.0f-%.0f%%",bin_cent[ic],bin_cent[ic+1]),0.79, 1.21, minpt, 2.8, -0.005, 0.25,
                                      // 0.65, 0.05, 0.9, 0.5,
                                      0.2, 0.45, 0.4, 0.88,
                                      "GEANT4, UrQMD, Au+Au at #sqrt{s_{NN}}=7.7GeV", Form("Ch. hadrons, |#eta|<%1.1f",eta_cut),1,Form("Ratio to %s",title[ratioToMethod].Data()));
  can->SetName(Form("%.0f-%.0f%%",bin_cent[ic],bin_cent[ic+1]));
  can->SaveAs(Form("DiffFlow_UrQMD_7.7_%.0f-%.0f%%.png",bin_cent[ic],bin_cent[ic+1]));
  }


  return vGr;
}