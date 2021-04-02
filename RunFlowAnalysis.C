// C++ headers
#include <iostream>
#include <fstream>
// ROOT headers
#include <TStopwatch.h>
#include <TChain.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TDatabasePDG.h>
#include <TComplex.h>
#include <TString.h>
#include <TEnv.h>

// PicoDst headers
#include <PicoDstMCEvent.h>
#include <PicoDstRecoEvent.h>
#include <PicoDstMCTrack.h>
#include <PicoDstRecoTrack.h>
#include <PicoDstFHCal.h>

#include <IReader.h>
#include <PicoDstReader.h>
// Flow method headers
#include <QVector.h>
#include <FlowAnalysisWithEtaSubEventPlane.h>
#include <FlowAnalysisWithThreeEtaSubEventPlane.h>
#include <FlowAnalysisWithFHCalEventPlane.h>
#include <FlowAnalysisWithLeeYangZeros.h>
#include <FlowAnalysisWithScalarProduct.h>
#include <FlowAnalysisWithQCumulant.h>
#include <FlowAnalysisWithHighOrderQCumulant.h>
#include <FlowAnalysisWithLeeYangZerosEventPlane.h>

#include "utilities.C"

using std::cout;
using std::cerr;
using std::endl;
// Flags for flow methods, where *_1 is the first run over the data and *_2 is the second run (need to move these flags to config file!!)
Bool_t ETASUBEVENTPLANE_1 = 0;        // Eta-sub EP (first run)
Bool_t ETASUBEVENTPLANE_2 = 1;        // Eta-sub EP (second run)
Bool_t THREEETASUBEVENTPLANE_1 = 0;   // 3 eta-sub method (first run)
Bool_t THREEETASUBEVENTPLANE_2 = 0;   // 3 eta-sub method (second run)
Bool_t FHCALEVENTPLANE_1 = 0;         // FHCal EP (w.r.t. 1-st order harmonic) (first run)
Bool_t FHCALEVENTPLANE_2 = 0;         // FHCal EP (w.r.t. 1-st order harmonic) (second run)
Bool_t LYZ_SUM_1 = 0;                 // Lee-Yang Zeros using sum generating function (first run)
Bool_t LYZ_SUM_2 = 0;                 // Lee-Yang Zeros using sum generating function (second run)
Bool_t LYZ_SUM_PRODUCT_1 = 0;         // Lee-Yang Zeros using product generating function (first run) (integrated with sum GF at the moment, will be separated soon)
Bool_t LYZ_SUM_PRODUCT_2 = 0;         // Lee-Yang Zeros using product generating function (second run) (integrated with sum GF at the moment, will be separated soon)
Bool_t SCALARPRODUCT_1 = 0;           // Scalar product using eta-sub method (first run)
Bool_t SCALARPRODUCT_2 = 0;           // Scalar product using eta-sub method (second run)
Bool_t QCUMULANT = 0;                 // Q-Cumulants: 2- and 4-particle cumulants obtained by both standard and subevent methods 
Bool_t HIGHORDERQCUMULANT = 0;        // Q-Cumulants: 2- up to 8-particle cumulants using recursive algorithm
Bool_t LYZEP = 0;                     // one needs to run LYZ_SUM_1 & 2 (or LYZ_SUM_PRODUCT_1 & 2) before set this flag to kTRUE
Bool_t readMCTracks = 0; // 0 - read reco tracks, 1 - read MC tracks
// Kinetic cuts by default if not using config file
Double_t maxpt = 3.6;     // max pt for differential flow
Double_t minpt = 0.;      // min pt for differential flow
Double_t maxptRF = 3.;    // max pt for reference flow
Double_t minptRF = 0.2;   // min pt for reference flow
Double_t eta_cut = 1.5;   // pseudorapidity acceptance window for flow measurements 
Double_t eta_gap = 0.05;  // +-0.05, eta-gap between 2 eta sub-event
Int_t Nhits_cut = 16;     // minimum nhits of reconstructed tracks
Double_t DCAcut = 0.5;
Double_t pid_probability = 0.9;
Long64_t Nevents = -1;

Int_t debug = 1;

std::string format = "picodst";

void readConfig(const TString& _strFileName)
{
  if (_strFileName.Length() == 0)
  {
    return;
  }

  TEnv env(_strFileName);

  Nhits_cut = env.GetValue("Nhits_cut", 0);

  maxpt = env.GetValue("maxpt", 0.);
  minpt = env.GetValue("minpt", 0.);
  maxptRF = env.GetValue("maxptRF", 0.);
  minptRF = env.GetValue("minptRF", 0.);
  eta_cut = env.GetValue("eta_cut", 0.);
  eta_gap = env.GetValue("eta_gap", 0.);
  DCAcut = env.GetValue("DCAcut", 0.);
  pid_probability = env.GetValue("pid_probability", 0.);

  debug = env.GetValue("debug", 0);
  Nevents = env.GetValue("Nevents", 0);

  format = env.GetValue("format", "");
}

Bool_t trackCut(PicoDstRecoTrack *const &recoTrack)
{
  if (!recoTrack) { return false; }
  Double_t pt = recoTrack->GetPt();
  Double_t eta = recoTrack->GetEta();
  if (pt < minpt || pt > maxpt || fabs(eta) > eta_cut)     return false;
  if (fabs(recoTrack->GetDCAx()) > DCAcut)        return false; // DCAx cut
  if (fabs(recoTrack->GetDCAy()) > DCAcut)        return false; // DCAy cut
  if (fabs(recoTrack->GetDCAz()) > DCAcut)        return false; // DCAz cut
  if (recoTrack->GetNhits() < Nhits_cut)          return false; // TPC hits cut    
  return true;                                   
}

Bool_t trackCut(PicoDstMCTrack *const &mcTrack)
{
  if (!mcTrack) { return false; }
  Double_t pt = mcTrack->GetPt();
  Double_t eta = mcTrack->GetEta();
  auto particle = (TParticlePDG*) TDatabasePDG::Instance()->GetParticle(mcTrack->GetPdg());
  if (!particle) { return false; }
  Double_t charge = 1./3.*particle->Charge();
  if (pt < minpt || pt > maxpt || fabs(eta) > eta_cut || charge == 0) { return false; }
  return true;
}

Int_t findBin(Double_t pt)
{
  Int_t ipt = -1;
  for (Int_t j = 0; j < npt; j++) { if (pt >= pTBin[j] && pt < pTBin[j + 1]) ipt = j; }
  return ipt;
}

Int_t findId(PicoDstRecoTrack *recoTrack)
{
  Int_t fId = -1;

  Double_t charge = recoTrack->GetCharge();
  if (recoTrack->GetTofFlag() != 0 && recoTrack->GetTofFlag() != 4)
  {
    if (recoTrack->GetPidProbPion() > pid_probability && charge > 0)   fId = 1; // pion+
    if (recoTrack->GetPidProbKaon() > pid_probability && charge > 0)   fId = 2; // kaon+
    if (recoTrack->GetPidProbProton() > pid_probability && charge > 0) fId = 3; // proton
    if (recoTrack->GetPidProbPion() > pid_probability && charge < 0)   fId = 5; // pion-
    if (recoTrack->GetPidProbKaon() > pid_probability && charge < 0)   fId = 6; // kaon-
    if (recoTrack->GetPidProbProton() > pid_probability && charge < 0) fId = 7; // antiproton
  }

  return fId;
}

void RunFlowAnalysis(TString inputFileName, TString outputFileName, TString configFileName = "")
{
  TStopwatch timer;
  timer.Start();

  if (configFileName.Length() > 0)
  {
    readConfig(configFileName);
  }

  if (debug)
  {
    cout << "Nevents = " << Nevents << endl;
    cout << "Nhits_cut = " << Nhits_cut << endl;

    cout << "maxpt = " << maxpt << endl;
    cout << "minpt = " << minpt << endl;
    cout << "maxptRF = " << maxptRF << endl;
    cout << "minptRF = " << minptRF << endl;
    cout << "eta_cut = " << eta_cut << endl;
    cout << "eta_gap = " << eta_gap << endl;
    cout << "DCAcut = " << DCAcut << endl;
    cout << "pid_probability = " << pid_probability << endl;
    cout << "format = " << format << endl;
  }

  if ((LYZ_SUM_1 && LYZ_SUM_PRODUCT_1) || (LYZ_SUM_2 && LYZ_SUM_PRODUCT_2))
  {
    cerr << "Both LYZ_SUM_1(2) and LYZ_SUM_PRODUCT_1(2) are TRUE. Set one of them to FASLE to run flow analysis." << endl;
    return;
  }

  // Configure input information
  TChain *chain = initChain(inputFileName, format.c_str());

  PicoDstMCEvent *mcEvent = nullptr;

  IReader* reader = nullptr;
  if (format == "picodst")
  {
    reader = new PicoDstReader();
  }
  if (!reader)
  {
    cerr << "No valid format is set!" << endl;
    return;
  }

  reader->Init(chain);

  // Configure output information
  TFile *fo = new TFile(outputFileName.Data(), "recreate");
  
  FlowAnalysisWithEtaSubEventPlane  *flowEtaSub  = nullptr; // Eta-sub Event Plane
  FlowAnalysisWithThreeEtaSubEventPlane *flowThreeEtaSub = nullptr; // 3-Eta-sub Event Plane
  FlowAnalysisWithFHCalEventPlane   *flowFHCalEP = nullptr; // FHCal Event Plane
  FlowAnalysisWithLeeYangZeros      *flowLYZ     = nullptr; // Lee Yang Zeros
  FlowAnalysisWithScalarProduct     *flowSP      = nullptr; // Scalar Product
  FlowAnalysisWithQCumulant         *flowQC      = nullptr; // Q-Cumulant
  FlowAnalysisWithHighOrderQCumulant *flowHighQC = nullptr;
  FlowAnalysisWithLeeYangZerosEventPlane *flowLYZEP = nullptr;

  if (ETASUBEVENTPLANE_1) {
    flowEtaSub = new FlowAnalysisWithEtaSubEventPlane();
    flowEtaSub->SetFirstRun(true);
    flowEtaSub->SetHarmonic(3);
    flowEtaSub->SetEtaGap(eta_gap);
    flowEtaSub->Init();
  }
  if (ETASUBEVENTPLANE_2) {
    flowEtaSub = new FlowAnalysisWithEtaSubEventPlane();
    flowEtaSub->SetFirstRun(false);
    flowEtaSub->SetHarmonic(3);
    flowEtaSub->SetEtaGap(eta_gap);
    flowEtaSub->SetDebugFlag(debug);
    flowEtaSub->SetInputFileFromFirstRun("FirstRun.root"); // need to be improve!!!
    flowEtaSub->Init();
  }
  if (THREEETASUBEVENTPLANE_1) {
    flowThreeEtaSub = new FlowAnalysisWithThreeEtaSubEventPlane();
    flowThreeEtaSub->SetFirstRun(true);
    flowThreeEtaSub->SetEtaGap(eta_gap);
    flowThreeEtaSub->Init();
  }
  if (THREEETASUBEVENTPLANE_2) {
    flowThreeEtaSub = new FlowAnalysisWithThreeEtaSubEventPlane();
    flowThreeEtaSub->SetFirstRun(false);
    flowThreeEtaSub->SetEtaGap(eta_gap);
    flowThreeEtaSub->SetDebugFlag(debug);
    flowThreeEtaSub->SetInputFileFromFirstRun("FirstRun.root");
    flowThreeEtaSub->Init();
  }
  if (FHCALEVENTPLANE_1) {
    flowFHCalEP = new FlowAnalysisWithFHCalEventPlane();
    flowFHCalEP->SetFirstRun(true);
    flowFHCalEP->SetEtaGap(eta_gap);
    flowFHCalEP->Init();
  }
  if (FHCALEVENTPLANE_2) {
    flowFHCalEP = new FlowAnalysisWithFHCalEventPlane();
    flowFHCalEP->SetFirstRun(false);
    flowFHCalEP->SetEtaGap(eta_gap);
    flowFHCalEP->SetDebugFlag(debug);
    flowFHCalEP->SetInputFileFromFirstRun("FirstRun.root");
    flowFHCalEP->Init();
  }

  if (LYZ_SUM_1) {
    flowLYZ = new FlowAnalysisWithLeeYangZeros();
    flowLYZ->SetFirstRun(true);
    flowLYZ->Init();
  }
  if (LYZ_SUM_2) {
    flowLYZ = new FlowAnalysisWithLeeYangZeros();
    // flowLYZ->SetDebugFlag(true);
    flowLYZ->SetFirstRun(false);
    flowLYZ->SetInputFileFromFirstRun("FirstRun.root"); // need to be improve!!!
    flowLYZ->Init();
  }

  if (LYZ_SUM_PRODUCT_1) {
    flowLYZ = new FlowAnalysisWithLeeYangZeros();
    flowLYZ->SetFirstRun(true);
    flowLYZ->SetUseProduct(true);
    flowLYZ->Init();
  }
  if (LYZ_SUM_PRODUCT_2) {
    flowLYZ = new FlowAnalysisWithLeeYangZeros();
    flowLYZ->SetDebugFlag(debug);
    flowLYZ->SetFirstRun(false);
    flowLYZ->SetUseProduct(true);
    flowLYZ->SetInputFileFromFirstRun("FirstRun.root"); // need to be improve!!!
    flowLYZ->Init();
  }
  QVector *Q2 = new QVector(); // for LYZ only, need to be improve to be also feasible to QC.
  // Start event loop
  if (SCALARPRODUCT_1) {
    flowSP = new FlowAnalysisWithScalarProduct();
    flowSP->SetFirstRun(true);
    flowSP->SetEtaGap(eta_gap);
    flowSP->Init();
  }
  if (SCALARPRODUCT_2) {
    flowSP = new FlowAnalysisWithScalarProduct();
    flowSP->SetFirstRun(false);
    flowSP->SetEtaGap(eta_gap);
    flowSP->SetInputFileFromFirstRun("FirstRun.root"); // need to be improve!!!
    flowSP->Init();
  }
  if (QCUMULANT) {
    flowQC = new FlowAnalysisWithQCumulant();
    flowQC->SetEtaGap(eta_gap);
    flowQC->Init();
  }

  if (HIGHORDERQCUMULANT) {
    flowHighQC = new FlowAnalysisWithHighOrderQCumulant();
    flowHighQC->Init();
  }
  if (LYZEP) {
    // if (inputHistogramFileName=="" || inputHistFromLYZSecondRun=="")
    // {
    //   cerr << "Input files with needed histograms for Lee Yang Zeros Event Plane aren't set" << endl;
    //   return;
    // }
    flowLYZEP = new FlowAnalysisWithLeeYangZerosEventPlane();
    // flowLYZEP->SetInputFileFromFirstAndSecondRun(inputHistogramFileName, inputHistFromLYZSecondRun);
    flowLYZEP->SetInputFileFromFirstAndSecondRun("FirstRun.root", "SecondRun.root");
    flowLYZEP->Init();
  }
  Int_t icent, reco_mult, ipt, fId;
  Double_t bimp, cent, pt, eta, phi, charge, energy;
  Long64_t chain_size = chain->GetEntries();
  Long64_t n_entries = (Nevents < chain_size && Nevents > 0) ? Nevents : chain_size;
  cout << "Hi Master, let's do some physics together..." << endl;
  for (Int_t iEv = 0; iEv < n_entries; iEv++)
  {
    if (iEv % 10000 == 0)
      std::cout << "Event [" << iEv << "/" << n_entries << "]" << std::endl;
    // chain->GetEntry(iEv);
    mcEvent = reader->ReadMcEvent(iEv);
 
    // Read MC event
    bimp = mcEvent->GetB();
    cent = CentB(bimp);
    if (cent == -1)
      continue;
    icent = GetCentBin(cent);
    
    if (readMCTracks) reco_mult = reader->GetMcTrackSize();
    else reco_mult = reader->GetRecoTrackSize();

    if (ETASUBEVENTPLANE_1 || ETASUBEVENTPLANE_2) flowEtaSub->Zero();
    if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->Zero();
    if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2) flowFHCalEP->Zero();
    if (LYZ_SUM_1 || LYZ_SUM_2 || LYZ_SUM_PRODUCT_1 || LYZ_SUM_PRODUCT_2) flowLYZ->Zero();
    if (SCALARPRODUCT_1 || SCALARPRODUCT_2) flowSP->Zero();
    Q2->Zero();
    if (QCUMULANT) flowQC->Zero();
    if (HIGHORDERQCUMULANT) flowHighQC->Zero();
    if (LYZEP) flowLYZEP->Zero();
    
    if ((FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2 || THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) && !readMCTracks)
    {
      Int_t Nmodules = reader->GetNFHCalModules();
      for (Int_t iModule = 0; iModule < Nmodules; iModule++)
      {
        auto module = (PicoDstFHCal *) reader->ReadFHCalModule(iModule);
        energy = module->GetEnergy();
        phi = GetFHCalPhi(iModule);
        eta = (iModule < 45) ? -3. : 3.; // Left & Right FHCal
        if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2) flowFHCalEP->ProcessFirstTrackLoop(eta, phi, energy);
        if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->ProcessFirstTrackLoopFHCal(eta, phi, energy);
      }
    }
    
    for (Int_t iTrk = 0; iTrk < reco_mult; iTrk++)
    { // Track loop
      if (readMCTracks)
      {
        auto mcTrack = (PicoDstMCTrack *) reader->ReadMcTrack(iTrk);
        
        if (!mcTrack) { cout << "skip mcTrack" << endl; continue;}
        pt = mcTrack->GetPt();
        eta = mcTrack->GetEta();
        phi = mcTrack->GetPhi();
        if ((FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2 || THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) && readMCTracks)
        {
          if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2) flowFHCalEP->ProcessFirstTrackLoop(eta, phi, pt);
          if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->ProcessFirstTrackLoopFHCal(eta, phi, pt);          
        }
        if (!trackCut(mcTrack)) { continue; } // TPC cut
        auto particle = (TParticlePDG*) TDatabasePDG::Instance()->GetParticle(mcTrack->GetPdg());   
        charge = 1./3.*particle->Charge();
        fId = findId(mcTrack);
      }
      else
      {
        auto recoTrack = (PicoDstRecoTrack *) reader->ReadRecoTrack(iTrk);
        if (!trackCut(recoTrack))
        {
          continue;
        }
        pt = recoTrack->GetPt();
        eta = recoTrack->GetEta();
        phi = recoTrack->GetPhi();
        charge = recoTrack->GetCharge();
        fId = findId(recoTrack);
      }
      ipt = findBin(pt);
      if (pt > minptRF && pt < maxptRF)
      { // Reference Flow pt cut
        if (ETASUBEVENTPLANE_1 || ETASUBEVENTPLANE_2) flowEtaSub->ProcessFirstTrackLoop(eta, phi, pt);
        if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->ProcessFirstTrackLoopTPC(eta, phi, pt);
        if (SCALARPRODUCT_1 || SCALARPRODUCT_2) flowSP->ProcessFirstTrackLoop(eta, phi);
        if (LYZ_SUM_1 || LYZ_SUM_2 || LYZ_SUM_PRODUCT_1 || LYZ_SUM_PRODUCT_2) flowLYZ->ProcessFirstTrackLoop(phi, pt, icent);
        Q2->CalQVector(phi, 1.);
        if (QCUMULANT) flowQC->ProcessFirstTrackLoopRP(eta, phi);
        if (HIGHORDERQCUMULANT) flowHighQC->ProcessFirstTrackLoopRP(phi);
      }

      if (QCUMULANT) flowQC->ProcessFirstTrackLoopPOI(ipt, eta, phi, fId, charge);
    } // end of track loop
    
    if (ETASUBEVENTPLANE_1 || ETASUBEVENTPLANE_2) flowEtaSub->ProcessEventAfterFirstTrackLoop(cent);
    if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->ProcessEventAfterFirstTrackLoop(cent);    
    if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2) flowFHCalEP->ProcessEventAfterFirstTrackLoop(cent);
    if (SCALARPRODUCT_1 || SCALARPRODUCT_2) flowSP->ProcessEventAfterFirstTrackLoop(cent);
    if (LYZ_SUM_1 || LYZ_SUM_2 || LYZ_SUM_PRODUCT_1 || LYZ_SUM_PRODUCT_2) flowLYZ->ProcessEventAfterFirstTrackLoop(Q2, icent);
    if (QCUMULANT) flowQC->ProcessEventAfterFirstTrackLoop(icent);
    if (HIGHORDERQCUMULANT) flowHighQC->ProcessEventAfterFirstTrackLoop(icent);
    if (LYZEP) flowLYZEP->ProcessEventAfterFirstTrackLoop(Q2, icent);
    if (ETASUBEVENTPLANE_2 || FHCALEVENTPLANE_2 || THREEETASUBEVENTPLANE_2 || LYZ_SUM_2 || LYZ_SUM_PRODUCT_2 || SCALARPRODUCT_2 || LYZEP)
    {
      for (Int_t iTrk = 0; iTrk < reco_mult; iTrk++)
      { // 2nd Track loop
        if (readMCTracks)
        {
          auto mcTrack = (PicoDstMCTrack *) reader->ReadMcTrack(iTrk);
          pt = mcTrack->GetPt();
          eta = mcTrack->GetEta();
          phi = mcTrack->GetPhi();
          if (!trackCut(mcTrack)) { continue; }
        }
        else 
        {
          auto recoTrack = (PicoDstRecoTrack *) reader->ReadRecoTrack(iTrk); 
          pt = recoTrack->GetPt();
          eta = recoTrack->GetEta();
          phi = recoTrack->GetPhi();
          charge = recoTrack->GetCharge();
          if (!trackCut(recoTrack)) { continue; }
        }


        if (ETASUBEVENTPLANE_2) flowEtaSub->ProcessSecondTrackLoop(eta, phi, pt, cent);
        if (THREEETASUBEVENTPLANE_2) flowThreeEtaSub->ProcessSecondTrackLoop(eta, phi, pt, cent);
        if (FHCALEVENTPLANE_2) flowFHCalEP->ProcessSecondTrackLoop(eta, phi, pt, cent);
        if (SCALARPRODUCT_2) flowSP->ProcessSecondTrackLoop(eta, phi, pt, cent);
        if (LYZ_SUM_2 || LYZ_SUM_PRODUCT_2) flowLYZ->ProcessSecondTrackLoop(phi, pt, icent);
        if (LYZEP) flowLYZEP->ProcessSecondTrackLoop(eta, phi, pt, cent);
      }
    }
    
  } // end event loop
  
  // Writing output
  fo->cd();
  if (ETASUBEVENTPLANE_1 || ETASUBEVENTPLANE_2) flowEtaSub->SaveHist();
  if (THREEETASUBEVENTPLANE_1 || THREEETASUBEVENTPLANE_2) flowThreeEtaSub->SaveHist();
  if (FHCALEVENTPLANE_1 || FHCALEVENTPLANE_2) flowFHCalEP->SaveHist();
  if (SCALARPRODUCT_1 || SCALARPRODUCT_2) flowSP->SaveHist();
  if (LYZ_SUM_1 || LYZ_SUM_2 || LYZ_SUM_PRODUCT_1 || LYZ_SUM_PRODUCT_2) flowLYZ->SaveHist();
  if (QCUMULANT) flowQC->SaveHist();
  if (HIGHORDERQCUMULANT) flowHighQC->SaveHist();
  if (LYZEP) flowLYZEP->SaveHist();
  fo->Close();

  timer.Stop();
  timer.Print();
}