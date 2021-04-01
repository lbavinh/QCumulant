#define PLOTV2VSPTQCUMULANT
#include "DrawTGraph.C"

double Covariance(TProfile *const &hcovXY, TProfile *const &hX, TProfile *const &hY, Int_t binXY=0, Int_t binX=0, Int_t binY=0){

  double mSumWXY = hcovXY->GetBinEntries(binXY+1);
  double sumWX = hX->GetBinEntries(binX+1);
  double sumWY = hY->GetBinEntries(binY+1);

  double meanXY = hcovXY -> GetBinContent(binXY+1);
  double meanX = hX -> GetBinContent(binX+1);
  double meanY = hY -> GetBinContent(binY+1);
  double mVal = (meanXY-meanX*meanY)/(sumWX*sumWY/mSumWXY-1.); // Cov(x,y)/(sumWX*sumWY/sumWXY)
  return mVal;
}

struct term{ // structure for "Mean squared error of MEAN" calculation, using unbiased estimator for the root of variance
  term(){
    mVal = 0;
    mMSE = 0;
  }
  term(TProfile *const &pr, Int_t bin=0){

    double Neff = pr -> GetBinEffectiveEntries(bin+1);
    mVal = pr -> GetBinContent(bin+1);
    pr -> SetErrorOption("s");
    double stdevW = pr -> GetBinError(bin+1);
    double S2 = stdevW*stdevW/(1-1./Neff);
    mMSE = S2/Neff;
  };
 public: 
  double mVal; // weithted mean value
  double mMSE; // Mean squared error of mean, https://en.wikipedia.org/wiki/Mean_squared_error

};

void PlotV2vsPtQCumulant(TString inFileName = "test.root")
{
  TFile *inFile = new TFile(inFileName.Data(),"read");
  // TFile *outFile = new TFile("graphs_v2.root","recreate"); // Save TGraphErrors
  TString outDirName = "pics";
  
  // Flags
  bool saveAsPNG = false;
  int excludeMethod = -1; // not including i-th method in v2 plotting, where i=0,1,2,3 correspond v22,v24,v2eta-sub,v22eta-gap, respectively
  int drawDifferentialFlowTill = 3; // Draw v2 vs pT (10% centrality cut) till: 0: no drawing; 1: till 10%; 2: till 20%; etc.
  // Constants
  const int npid_plotQC = 12; // CH+, pion+, kaon+, proton, CH-, pion-, kaon-, antiproton, CH, pions, kaons, protons+antiproton
  const std::vector<TString> pidNames = {"hadron_pos", "pion_pos", "kaon_pos", "proton", "hadron_neg", "pion_neg", "kaon_neg", "proton_bar", "hadron", "pion", "kaon","proton_antiproton"};
  const std::vector<TString> pidFancyNames = {"h^{+}", "#pi^{+}", "K^{+}", "p", "h^{-}", "#pi^{-}", "K^{-}", "#bar{p}", "h^{#pm}","#pi^{#pm}","K^{#pm}","p(#bar{p})"};

  const int nmethod = 3; // 2QC, 4QC, EP, 2QC-gapped

  const int npt = 16; // 0.5 - 3.6 GeV/c - number of pT bins
  const double bin_pT[npt+1]={0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.2,3.6};
  const int binMinPtRFP = 1;  // 0.2 GeV 
  const int binMaxPtRFP = 14; // 2.8 GeV
  const double minptRFP = 0.2;
  const double maxptRFP = 3.0;

  const double maxpt = 2.5; // for v2 vs pt plotting
  const double minpt = 0.;  // for v2 vs pt plotting

  const int ncent = 9; // 0-80 %
  const double bin_cent[ncent] = {2.5,7.5,15,25,35,45,55,65,75};
  const double bin_centE[ncent] = {0};
  const vector<pair<int,int>> centrality = {{0,5},{5,10},{10,20},{20,30},{30,40},{40,50},{50,60},{60,70},{70,80}};
  const float eta_gap = 0.05;

  const double mincent = 0.;  // for v2 vs centrality
  const double maxcent = 60.; // for v2 vs centrality

  const double minV2int = -0.005; // for v2 vs centrality plotting
  const double maxV2int = 0.1; // for v2 vs centrality plotting
  const double minV2dif = -0.01; // for v2 vs pt plotting
  const double maxV2dif = 0.2; // for v2 vs pt plotting


  vector <Double_t> coordinateLeg = {0.18,0.63,0.45,0.889};
  vector<pair<Double_t,Double_t>> rangeRatio = {{0.84,1.16},{0.84,1.16}};
  vector<pair<Double_t,Double_t>> rangeRatioRF ={{0.65,1.11},{0.65,1.11}};
  int marker[]={21,20,22}; // 2QC, 4QC, 2QC-gapped

  // Input hist

  TProfile *pCorrelator2EtaGap = (TProfile*)inFile->Get("pCorrelator2EtaGap");
  TProfile *pCorrelator2 = (TProfile*)inFile->Get("pCorrelator2");
  TProfile *pCorrelator4 = (TProfile*)inFile->Get("pCorrelator4");
  TProfile2D *pReducedCorrelator2EtaGap[npid_plotQC]; // <<2'>> (with eta-gap)
  TProfile2D *pReducedCorrelator2[npid_plotQC]; // <<2'>>
  TProfile2D *pReducedCorrelator4[npid_plotQC]; // <<4'>>
  TProfile2D *pCov22RedEtaGap[npid_plotQC];
  TProfile *pCov24 = (TProfile*)inFile->Get("pCov24");
  TProfile2D *pCov22Red[npid_plotQC];
  TProfile2D *pCov24Red[npid_plotQC];
  TProfile2D *pCov42Red[npid_plotQC];
  TProfile2D *pCov44Red[npid_plotQC];
  TProfile2D *pCov2Red4Red[npid_plotQC];

  for (int i=0; i<npid_plotQC-4; i++)
  {
    pReducedCorrelator2EtaGap[i] = (TProfile2D*)inFile->Get(Form("pReducedCorrelator2EtaGap_pid%i",i));
    pReducedCorrelator2[i]   = (TProfile2D*) inFile->Get(Form("pReducedCorrelator2_pid%i",i));
    pReducedCorrelator4[i]   = (TProfile2D*) inFile->Get(Form("pReducedCorrelator4_pid%i",i));
    pCov22RedEtaGap[i]       = (TProfile2D*) inFile->Get(Form("pCov22RedEtaGap_pid%i",i));
    pCov22Red[i]             = (TProfile2D*) inFile->Get(Form("pCov22Red_pid%i",i));
    pCov24Red[i]             = (TProfile2D*) inFile->Get(Form("pCov24Red_pid%i",i));
    pCov42Red[i]             = (TProfile2D*) inFile->Get(Form("pCov42Red_pid%i",i));
    pCov44Red[i]             = (TProfile2D*) inFile->Get(Form("pCov44Red_pid%i",i));
    pCov2Red4Red[i]          = (TProfile2D*) inFile->Get(Form("pCov2Red4Red_pid%i",i));
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
  // OUTPUT
  TGraphErrors *grDifFl[nmethod][ncent][npid_plotQC];    // v2(pt)
  
  // Filling pT bin
  double pt[npt];
  double ept[npt]={0}; // error bin pT = 0.0
  for (int ipt=0; ipt<npt; ipt++){
    pt[ipt] = ( bin_pT[ipt] + bin_pT[ipt+1] ) / 2.;
  }
  
  double v2Dif[nmethod][ncent][npid_plotQC][npt], v2eDif[nmethod][ncent][npid_plotQC][npt];
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
        if (id==8 && icent==3) cout << v22Dif <<" ";
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
        // if (id==8 && icent==3) cout << v24Dif <<" ";
        // if (id==8 && icent==3) cout << ev24Dif <<" ";
        // v22 Gapped
        term cor2redGap = term(pReducedCorrelator2EtaGap_cent[id][icent],ipt);
        double v22DifGap = cor2redGap.mVal/v22Gap;
        double cov22primeGap = Covariance(pCov22RedEtaGap_cent[id][icent],pCorrelator2EtaGap,pReducedCorrelator2EtaGap_cent[id][icent],ipt,icent,ipt);
        double ev22DifGap = sqrt(0.25*pow(cor2Gap.mVal,-3)*(pow(cor2redGap.mVal,2)*cor2Gap.mMSE
                            + 4*pow(cor2Gap.mVal,2)*cor2redGap.mMSE - 4*cor2Gap.mVal*cor2redGap.mVal*cov22primeGap));
        v2Dif[2][icent][id][ipt] = v22DifGap;
        v2eDif[2][icent][id][ipt] = ev22DifGap;
        // if (id==8 && icent==3) cout << v22DifGap <<" ";
        // if (id==8 && icent==3) cout << ev22DifGap <<" ";
      } // end of loop for all pT bin
      for (int i=0; i<nmethod; i++){
        grDifFl[i][icent][id] = new TGraphErrors(npt,pt,v2Dif[i][icent][id],ept,v2eDif[i][icent][id]);
        grDifFl[i][icent][id] -> SetMarkerStyle(marker[i]);
        grDifFl[i][icent][id] -> SetMarkerSize(1.5);
        grDifFl[i][icent][id] -> SetDrawOption("P");
      }
    } // end of loop over PID
  } // end of loop over centrality classes
  cout << endl;
  const char *grTitleDF[nmethod]={"v_{2}{2};p_{T} [GeV/c];v_{2}",
                                  "v_{2}{4};p_{T} [GeV/c];v_{2}",
                                  "v_{2}{2,#eta-gap};p_{T} [GeV/c];v_{2}"};
  // outFile -> cd();
  for (int imeth=0; imeth<nmethod; imeth++){
    for (int id=0;id<npid_plotQC;id++){
      for (int icent=0;icent<ncent;icent++){
        grDifFl[imeth][icent][id] -> SetTitle(grTitleDF[imeth]);
        // grDifFl[imeth][icent][id] -> Write(Form("gr_cent%i_%i_%i",icent,imeth,id));
      }
    }
  }
  TString level = (TString) Form("UrQMD, Au+Au at #sqrt{s_{NN}}=7.7 GeV");
  if (saveAsPNG) gSystem->Exec(Form("mkdir -p ./%s/",outDirName.Data()));
  TCanvas *cV2PT[ncent][npid_plotQC];
  char hname[800];
  for (int icent=0; icent<drawDifferentialFlowTill; icent++){
    for (int id=0;id<npid_plotQC;id++){
      std::vector<TGraphErrors*> vgrv2pt;
      vgrv2pt.push_back(grDifFl[2][icent][id]); // v2{gapped 2QC}
      for (int i=0; i<nmethod-1; i++){
        if (i==excludeMethod) continue;
        vgrv2pt.push_back(grDifFl[i][icent][id]);
      }
      sprintf(hname,"%s, %i-%i%%",pidFancyNames.at(id).Data(),centrality.at(icent).first,centrality.at(icent).second);
      cV2PT[icent][id] = (TCanvas*) DrawTGraph(vgrv2pt,"",rangeRatio.at(0).first, rangeRatio.at(0).second, minpt, maxpt, minV2dif, maxV2dif,
                                               coordinateLeg.at(0), coordinateLeg.at(1), coordinateLeg.at(2), coordinateLeg.at(3),
                                               level.Data(), hname);
      cV2PT[icent][id] -> SetName(hname);
      if (saveAsPNG) cV2PT[icent][id] -> SaveAs(Form("./%s/DifferentialFlow_Centrality%i-%i_%s.png",outDirName.Data(),centrality.at(icent).first,centrality.at(icent).second,pidNames.at(id).Data()));
    }
  }
  inFile->Close();
  // outFile->Close();

}
