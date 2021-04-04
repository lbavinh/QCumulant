#ifndef FLOWANALYSISWITHQCUMULANT_H
#define FLOWANALYSISWITHQCUMULANT_H
#include <iostream>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TComplex.h>
#include "QVector.h"

using std::cerr;
using std::cout;
using std::endl;

class FlowAnalysisWithQCumulant
{
public:
  FlowAnalysisWithQCumulant();
  virtual ~FlowAnalysisWithQCumulant();
  void Init();
  void Zero(); // Reset variables for new event loop
  void ProcessFirstTrackLoopRP(const Double_t &eta, const Double_t &phi);
  void ProcessFirstTrackLoopPOI(const Int_t &ipt, const Double_t &eta, const Double_t &phi, const Int_t &pid, const Double_t &charge);
  void ProcessEventAfterFirstTrackLoop(const Int_t &icent);
  void SetHarmonic(Int_t i) { this->fHarmonic = i; }
  void SetEtaGap(Double_t d) { this->fEtaGap = d; }
  void SaveHist();

  TComplex Qstar(const TComplex &Q);
  Double_t CalCor22(const TComplex &Q2, const Double_t &M, const Double_t &w2);
  Double_t CalCor24(const TComplex &Q2, const TComplex &Q4, const Double_t &M, const Double_t &w4);
  Double_t CalRedCor22(const TComplex &Q2, const TComplex &p2, const Double_t &M,
                     const Double_t &mp, const Double_t &mq, const Double_t &wred2);
  Double_t CalRedCor24(const TComplex &Q2, const TComplex &Q4, const TComplex &p2, const TComplex &q2,
                     const TComplex &q4, const Double_t &M, const Double_t &mp, const Double_t &mq, const Double_t &wred4);

private:
  Int_t fHarmonic;
  Double_t fEtaGap;
  // 2,QC & 4,QC without eta-gap
  Double_t Qx2, Qy2, Qx4, Qy4;
  TComplex Q2, Q4;
  Double_t px2[npt][npid], py2[npt][npid];
  TComplex p2[npt][npid], q2[npt][npid], q4[npt][npid];
  Double_t qx2[npt][npid], qy2[npt][npid], qx4[npt][npid], qy4[npt][npid];
  Double_t M = 0.;
  Double_t mq[npt][npid], mp[npt][npid];
  Double_t redCor22[npt][npid], redCor24[npt][npid];
  Double_t w2, w4;
  Double_t wred2[npt][npid], wred4[npt][npid];
  Double_t cor22, cor24;
  // Int_t icent;

  // 2,QC with eta-gap
  Double_t Qx2Gap[neta], Qy2Gap[neta];
  Double_t px2Gap[neta][npt][npid], py2Gap[neta][npt][npid];
  TComplex Q2Gap[neta], p2Gap[neta][npt][npid];
  Double_t MGap[neta];
  Double_t mpGap[neta][npt][npid];
  Double_t w2Gap;
  Double_t wred2Gap[neta][npt][npid];
  Double_t cor22Gap;
  Double_t redCor22Gap[neta][npt][npid];

  TProfile *pCorrelator2EtaGap;
  TProfile *pCorrelator2;
  TProfile *pCorrelator4;
  TProfile2D *pReducedCorrelator2EtaGap[npid];
  TProfile2D *pReducedCorrelator2[npid];
  TProfile2D *pReducedCorrelator4[npid];

  TProfile *pCov24;               // <2>*<4>
  TProfile2D *pCov22Red[npid];    // <2>*<2'>
  TProfile2D *pCov24Red[npid];    // <2>*<4'>
  TProfile2D *pCov42Red[npid];    // <4>*<2'>
  TProfile2D *pCov44Red[npid];    // <4>*<4'>
  TProfile2D *pCov2Red4Red[npid]; // <2'>*<4'>
  TProfile2D *pCov22RedEtaGap[npid];

  ClassDef(FlowAnalysisWithQCumulant, 0);
};

#endif

