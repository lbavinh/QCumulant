#include <FlowAnalysisWithEtaSubEventPlane.h>
ClassImp(FlowAnalysisWithEtaSubEventPlane);

FlowAnalysisWithEtaSubEventPlane::FlowAnalysisWithEtaSubEventPlane(bool bFirstRun, TString inputFileFromFirstRun) :fFirstRun(bFirstRun)
{
  fPrRes = new TProfile("prRes", "EP resolution", ncent, &bin_cent[0]);
  if (!fFirstRun) 
  {
    fPrV2EtaSubEventPlane = new TProfile3D("prV2EtaSubEventPlane", "", ncent, &bin_cent[0], npt, &pTBin[0], netaBin, &etaBin[0]);
    GetRes(inputFileFromFirstRun);
  }
}

FlowAnalysisWithEtaSubEventPlane::~FlowAnalysisWithEtaSubEventPlane()
{
}

void FlowAnalysisWithEtaSubEventPlane::Zero()
{
  fPsi_L = 0.;
  fPsi_R = 0.;
  Qvector_L.Zero();
  Qvector_R.Zero();
}

void FlowAnalysisWithEtaSubEventPlane::ProcessFirstTrackLoop(const double &eta, const double &phi, const double &pt)
{
  if (eta < - eta_gapEP)
  {
    Qvector_L.CalQVector(phi, pt);
  }
  if (eta > eta_gapEP)
  {
    Qvector_R.CalQVector(phi, pt);
  }
}

void FlowAnalysisWithEtaSubEventPlane::ProcessEventAfterFirstTrackLoop(const double &dCent)
{
  if (Qvector_L.GetMult() > mult_EP_cut && Qvector_R.GetMult() > mult_EP_cut)
  {
    fMultCut = false;
    Qvector_L.WeightQVector();
    Qvector_R.WeightQVector();
    fPsi_L = 0.5 * TMath::ATan2(Qvector_L.Y(), Qvector_L.X());
    fPsi_R = 0.5 * TMath::ATan2(Qvector_R.Y(), Qvector_R.X());
    fPrRes->Fill(dCent, TMath::Cos( 2.0 * (fPsi_L - fPsi_R) ));
  }
  else
  {
    fMultCut = true;
  }
  
}

void FlowAnalysisWithEtaSubEventPlane::GetRes(TString inputFileFromFirstRun)
{
  if (!fFirstRun)
  {
    if (!inputFileFromFirstRun) { cerr << "Warning: inputFileFromFirstRun=NULL" << endl;}
    TFile *fi = new TFile(inputFileFromFirstRun.Data(), "read");
    fPrRes = (TProfile*)fi->Get("prRes");
    for (int ic = 0; ic < ncent; ic++)
    {
      fRes2[ic] = TMath::Sqrt(fPrRes->GetBinContent(ic+1));
    }
  }
}

void FlowAnalysisWithEtaSubEventPlane::ProcessSecondTrackLoop(const double &eta, const double &phi, const double &pt, const double &dCent)
{
  if (!fMultCut && !fFirstRun)
  {
    double v2EtaSubEventPlane = -999.0;
    if (eta < -eta_gapEP)
    {
      v2EtaSubEventPlane = TMath::Cos( 2.0 * (phi - fPsi_R) );
    }
    else if (eta > eta_gapEP)
    {
      v2EtaSubEventPlane = TMath::Cos( 2.0 * (phi - fPsi_L) );
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