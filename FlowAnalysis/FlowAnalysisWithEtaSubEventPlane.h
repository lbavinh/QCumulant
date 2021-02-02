#ifndef FLOWANALYSISWITHETASUBEVENTPLANE_H
#define FLOWANALYSISWITHETASUBEVENTPLANE_H
#include <iostream>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>
#include <TH2F.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile3D.h>
#include <TDatabasePDG.h>
#include <TString.h>

#include "QVector.h"
#include "../constants.C"
using std::cerr;
using std::cout;
using std::endl;

class FlowAnalysisWithEtaSubEventPlane
{
public:
  FlowAnalysisWithEtaSubEventPlane(bool bFirstRun, TString inputFileFromFirstRun);
  virtual ~FlowAnalysisWithEtaSubEventPlane();
  void Zero(); // Reset variables for new event loop
  void ProcessFirstTrackLoop(const double &eta, const double &phi, const double &pt);
  void ProcessEventAfterFirstTrackLoop(const double &dCent);
  void ProcessSecondTrackLoop(const double &eta, const double &phi, const double &pt, const double &dCent);
  void GetRes(TString inputFileFromFirstRun);
  void SaveHist();
private:
  bool fFirstRun;
  bool fMultCut;
  double fPsi_L;
  double fPsi_R;
  QVector Qvector_L;
  QVector Qvector_R;
  double fRes2[ncent];

  TProfile *fPrRes;
  TProfile3D *fPrV2EtaSubEventPlane;
  ClassDef(FlowAnalysisWithEtaSubEventPlane,0);

};

#endif