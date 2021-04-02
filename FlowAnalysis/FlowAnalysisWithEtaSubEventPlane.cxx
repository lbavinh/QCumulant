#include <FlowAnalysisWithEtaSubEventPlane.h>
ClassImp(FlowAnalysisWithEtaSubEventPlane);

FlowAnalysisWithEtaSubEventPlane::FlowAnalysisWithEtaSubEventPlane() :
  fFirstRun(kTRUE),
  fMultCut(kTRUE),
  fDebug(kFALSE),
  fHarmonic(2),
  fPsi_L(0.),
  fPsi_R(0.),
  fQvector_L(nullptr),
  fQvector_R(nullptr),
  fRes2(),
  fEtaGap(0.),
  fstrInputFileFromFirstRun(""),
  fPrRes(nullptr),
  fPrV2EtaSubEventPlane(nullptr)
{
}

FlowAnalysisWithEtaSubEventPlane::~FlowAnalysisWithEtaSubEventPlane()
{
}

void FlowAnalysisWithEtaSubEventPlane::Init()
{
  fPrRes = new TProfile("prRes", "EP resolution", ncent, &bin_cent[0]);
  fQvector_L = new QVector(fHarmonic);
  fQvector_R = new QVector(fHarmonic);
  if (!fFirstRun) 
  {
    fPrV2EtaSubEventPlane = new TProfile3D("prV2EtaSubEventPlane", "", ncent, &bin_cent[0], npt, &pTBin[0], netaBin, &etaBin[0]);
    GetRes();
  }  
}

void FlowAnalysisWithEtaSubEventPlane::Zero()
{
  fPsi_L = 0.;
  fPsi_R = 0.;
  fQvector_L->Zero();
  fQvector_R->Zero();
}

void FlowAnalysisWithEtaSubEventPlane::ProcessFirstTrackLoop(const Double_t &eta, const Double_t &phi, const Double_t &pt)
{
  if (eta < - fEtaGap)
  {
    fQvector_L->CalQVector(phi, pt);
  }
  if (eta > fEtaGap)
  {
    fQvector_R->CalQVector(phi, pt);
  }
}

void FlowAnalysisWithEtaSubEventPlane::ProcessEventAfterFirstTrackLoop(const Double_t &dCent)
{
  if (fQvector_L->GetMult() > mult_EP_cut && fQvector_R->GetMult() > mult_EP_cut)
  {
    fMultCut = kFALSE;
    fQvector_L->WeightQVector();
    fQvector_R->WeightQVector();
    fPsi_L = TMath::ATan2(fQvector_L->Y(), fQvector_L->X())/fHarmonic;
    fPsi_R = TMath::ATan2(fQvector_R->Y(), fQvector_R->X())/fHarmonic;
    if (fFirstRun) fPrRes->Fill(dCent, TMath::Cos( fHarmonic * (fPsi_L - fPsi_R) ));
  }
  else
  {
    fMultCut = kTRUE;
  }
  
}

void FlowAnalysisWithEtaSubEventPlane::GetRes()
{
  if (!fFirstRun)
  {
    if (fstrInputFileFromFirstRun == "") 
    { cerr << "Warning: in FlowAnalysisWithEtaSubEventPlane::GetRes() fstrInputFileFromFirstRun="" " << endl;}
    TFile *fi = new TFile(fstrInputFileFromFirstRun.Data(), "read");
    fPrRes = (TProfile*)fi->Get("prRes");
    for (int ic = 0; ic < ncent; ic++)
    {
      fRes2[ic] = TMath::Sqrt(fPrRes->GetBinContent(ic+1));
    }
    if (fDebug)
    {
      cout << "TPC EP Resolution w.r.t. " << fHarmonic << "-th harmonic (2-eta-sub):" << endl;
      for (Int_t ic = 0; ic < ncent; ic++)
      {
        cout << fRes2[ic] <<", ";
      }
      cout << endl;
    }
  }
}

void FlowAnalysisWithEtaSubEventPlane::ProcessSecondTrackLoop(const Double_t &eta, const Double_t &phi, const Double_t &pt, const Double_t &dCent)
{
  if (!fMultCut && !fFirstRun)
  {
    Double_t v2EtaSubEventPlane = -999.0;
    if (eta < -fEtaGap)
    {
      v2EtaSubEventPlane = TMath::Cos( fHarmonic * (phi - fPsi_R) );
    }
    else if (eta > fEtaGap)
    {
      v2EtaSubEventPlane = TMath::Cos( fHarmonic * (phi - fPsi_L) );
    }
    else { return; }
    int icent = fPrRes->FindBin(dCent) - 1;
    if (fRes2[icent] != 0)
    { 
      v2EtaSubEventPlane /= fRes2[icent];
      fPrV2EtaSubEventPlane->Fill(dCent, pt, eta, v2EtaSubEventPlane);
    }
  }
}

void FlowAnalysisWithEtaSubEventPlane::SaveHist()
{
  fPrRes->Write();
  if (!fFirstRun) fPrV2EtaSubEventPlane->Write();
}

void FlowAnalysisWithEtaSubEventPlane::SaveHist(TDirectoryFile *const &outputDir)
{
  outputDir->Add(fPrRes);
  if (!fFirstRun) outputDir->Add(fPrV2EtaSubEventPlane);
  outputDir->Write();
}