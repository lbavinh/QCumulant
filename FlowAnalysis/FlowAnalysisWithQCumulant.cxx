/**
 * Elliptic flow v2 measurements using Q-Cumulant
 * proposed by Ante Bilandzic in https://arxiv.org/abs/1010.0233
 * coded by Vinh Ba Luong (lbavinh@gmail.com)
 * 25/11/2020
 */
#include <FlowAnalysisWithQCumulant.h>
ClassImp(FlowAnalysisWithQCumulant);
FlowAnalysisWithQCumulant::FlowAnalysisWithQCumulant() :
  fHarmonic(2),
  fEtaGap(0),
  pCorrelator2EtaGap(nullptr),
  pCorrelator2(nullptr),
  pCorrelator4(nullptr),
  pCov24(nullptr)
{
  for (Int_t id = 0; id < npid; id++){
    pReducedCorrelator2EtaGap[id] = nullptr;
    pReducedCorrelator2[id] = nullptr;
    pReducedCorrelator4[id] = nullptr;
    pCov22Red[id] = nullptr;
    pCov24Red[id] = nullptr;
    pCov42Red[id] = nullptr;
    pCov44Red[id] = nullptr;
    pCov2Red4Red[id] = nullptr;
    pCov22RedEtaGap[id] = nullptr;
  }
  Zero();
}

FlowAnalysisWithQCumulant::~FlowAnalysisWithQCumulant()
{
}

void FlowAnalysisWithQCumulant::Init()
{
  pCorrelator2EtaGap = new TProfile("pCorrelator2EtaGap", "2nd order correlator with eta-gap, TPC", ncent, 0, ncent); // <<2>> (with eta-gap)
  pCorrelator2EtaGap->Sumw2();
  pCorrelator2 = new TProfile("pCorrelator2", "2nd order correlator", ncent, 0, ncent); // <<2>>
  pCorrelator2->Sumw2();
  pCorrelator4 = new TProfile("pCorrelator4", "4th order correlator", ncent, 0, ncent); // <<4>>
  pCorrelator4->Sumw2();

  for (Int_t i = 0; i < npid; i++)
  {
    pReducedCorrelator2EtaGap[i] = new TProfile2D(Form("pReducedCorrelator2EtaGap_pid%i", i), Form("Reduced 2nd order correlator with eta-gap of %s (TPC", pidNames.at(i).Data()), npt, 0, npt, ncent, 0, ncent);
    pReducedCorrelator2EtaGap[i]->Sumw2();
    pReducedCorrelator2[i] = new TProfile2D(Form("pReducedCorrelator2_pid%i", i), Form("Reduced 2nd order correlator of %s", pidNames.at(i).Data()), npt, 0, npt, ncent, 0, ncent);
    pReducedCorrelator2[i]->Sumw2();
    pReducedCorrelator4[i] = new TProfile2D(Form("pReducedCorrelator4_pid%i", i), Form("Reduced 4th order correlator of %s", pidNames.at(i).Data()), npt, 0, npt, ncent, 0, ncent);
    pReducedCorrelator4[i]->Sumw2();
  }

  pCov24 = new TProfile("pCov24", "Covariance(<2>,<4>)", ncent, 0, ncent); // <2>*<4>
  pCov24->Sumw2();

  for (Int_t i = 0; i < npid; i++)
  {
    pCov22RedEtaGap[i] = new TProfile2D(Form("pCov22RedEtaGap_pid%i", i), Form("Covariance(<2>,<2'>) with eta-gap of %s (TPC)", pidNames.at(i).Data()), npt, 0, npt, ncent, 0, ncent);
    pCov22RedEtaGap[i]->Sumw2();
    pCov22Red[i] = new TProfile2D(Form("pCov22Red_pid%i", i), Form("Covariance(<2>,<2'>) of %s", pidNames.at(i).Data()), npt, 0, npt, ncent, 0, ncent);
    pCov22Red[i]->Sumw2();
    pCov24Red[i] = new TProfile2D(Form("pCov24Red_pid%i", i), Form("Covariance(<2>,<4'>) of %s", pidNames.at(i).Data()), npt, 0, npt, ncent, 0, ncent);
    pCov24Red[i]->Sumw2();
    pCov42Red[i] = new TProfile2D(Form("pCov42Red_pid%i", i), Form("Covariance(<4>,<2'>) of %s", pidNames.at(i).Data()), npt, 0, npt, ncent, 0, ncent);
    pCov42Red[i]->Sumw2();
    pCov44Red[i] = new TProfile2D(Form("pCov44Red_pid%i", i), Form("Covariance(<4>,<4'>) of %s", pidNames.at(i).Data()), npt, 0, npt, ncent, 0, ncent);
    pCov44Red[i]->Sumw2();
    pCov2Red4Red[i] = new TProfile2D(Form("pCov2Red4Red_pid%i", i), Form("Covariance(<4'>,<2'>) of %s", pidNames.at(i).Data()), npt, 0, npt, ncent, 0, ncent);
    pCov2Red4Red[i]->Sumw2();
  }
}

void FlowAnalysisWithQCumulant::Zero()
{
  Qx2 = 0.;
  Qy2 = 0.;
  Qx4 = 0.;
  Qy4 = 0.;
  Q2 = TComplex(0., 0.);
  Q4 = TComplex(0., 0.);
  M = 0.;
  w2 = 0.;
  w4 = 0.;
  cor22 = 0.;
  cor24 = 0.;
  for (Int_t ipt = 0; ipt < npt; ipt++)
  {
    for (Int_t ipid = 0; ipid < npid; ipid++)
    {
      px2[ipt][ipid] = 0.;
      py2[ipt][ipid] = 0.;
      qx2[ipt][ipid] = 0.;
      qy2[ipt][ipid] = 0.;
      qx4[ipt][ipid] = 0.;
      qy4[ipt][ipid] = 0.;

      mq[ipt][ipid] = 0.;
      mp[ipt][ipid] = 0.;
      redCor22[ipt][ipid] = 0.;
      redCor24[ipt][ipid] = 0.;
      wred2[ipt][ipid] = 0.;
      wred4[ipt][ipid] = 0.;
      p2[ipt][ipid] = TComplex(0., 0.);
      q2[ipt][ipid] = TComplex(0., 0.);
      q4[ipt][ipid] = TComplex(0., 0.);
    }
  }

  // 2,QC with eta-gap

  w2Gap = 0.;
  cor22Gap = 0.;
  for (Int_t ieta = 0; ieta < neta; ieta++)
  {
    Qx2Gap[ieta] = 0.;
    Qy2Gap[ieta] = 0.;
    MGap[ieta] = 0.;
    Q2Gap[ieta] = TComplex(0., 0.);
    for (Int_t ipt = 0; ipt < npt; ipt++)
    {
      for (Int_t ipid = 0; ipid < npid; ipid++)
      {
        px2Gap[ieta][ipt][ipid] = 0;
        py2Gap[ieta][ipt][ipid] = 0;
        mpGap[ieta][ipt][ipid] = 0;
        wred2Gap[ieta][ipt][ipid] = 0;
        redCor22Gap[ieta][ipt][ipid] = 0;
        p2Gap[ieta][ipt][ipid] = TComplex(0., 0.);
      }
    }
  }
}

void FlowAnalysisWithQCumulant::ProcessFirstTrackLoopRP(const Double_t &eta, const Double_t &phi)
{
  Double_t cos2nphi = TMath::Cos(2.*fHarmonic*phi);
  Double_t sin2nphi = TMath::Sin(2.*fHarmonic*phi);
  Double_t cosnphi = TMath::Cos(fHarmonic*phi);
  Double_t sinnphi = TMath::Sin(fHarmonic*phi);

  // 2,4-QC
  Qx2 += cosnphi;
  Qy2 += sinnphi;
  Qx4 += cos2nphi;
  Qy4 += sin2nphi;
  M++;

  // 2-QC, eta-gapped
  if (eta < -fEtaGap)
  { // Left TPC subevent selection
    Qx2Gap[0] += cosnphi;
    Qy2Gap[0] += sinnphi;
    MGap[0]++;
  }
  if (eta > fEtaGap)
  { // Right TPC subevent selection
    Qx2Gap[1] += cosnphi;
    Qy2Gap[1] += sinnphi;
    MGap[1]++;
  }

}

void FlowAnalysisWithQCumulant::ProcessFirstTrackLoopPOI(const Int_t &ipt, const Double_t &eta, const Double_t &phi, const Int_t &pid, const Double_t &charge)
{
  Double_t cos2nphi = TMath::Cos(2.*fHarmonic*phi);
  Double_t sin2nphi = TMath::Sin(2.*fHarmonic*phi);
  Double_t cosnphi = TMath::Cos(fHarmonic*phi);
  Double_t sinnphi = TMath::Sin(fHarmonic*phi);

  // Differential Flow of 2,4-QC
  if (charge > 0)
  {
    px2[ipt][0] += cosnphi;
    py2[ipt][0] += sinnphi;
    mp[ipt][0]++;

    qx2[ipt][0] += cosnphi;
    qy2[ipt][0] += sinnphi;
    qx4[ipt][0] += cos2nphi;
    qy4[ipt][0] += sin2nphi;
    mq[ipt][0]++;
  }
  if (charge < 0)
  {
    px2[ipt][4] += cosnphi;
    py2[ipt][4] += sinnphi;
    mp[ipt][4]++;

    qx2[ipt][4] += cosnphi;
    qy2[ipt][4] += sinnphi;
    qx4[ipt][4] += cos2nphi;
    qy4[ipt][4] += sin2nphi;
    mq[ipt][4]++;
  }
  if (pid > 0)
  {
    px2[ipt][pid] += cosnphi;
    py2[ipt][pid] += sinnphi;
    mp[ipt][pid]++;

    qx2[ipt][pid] += cosnphi;
    qy2[ipt][pid] += sinnphi;
    qx4[ipt][pid] += cos2nphi;
    qy4[ipt][pid] += sin2nphi;
    mq[ipt][pid]++;
  }

  // Differential Flow of 2-QC, eta-gapped
  if (eta < -fEtaGap)
  { // Left TPC subevent selection
    if (charge > 0)
    { // Inclusive positively charged hadrons
      px2Gap[1][ipt][0] += cosnphi;
      py2Gap[1][ipt][0] += sinnphi;
      mpGap[1][ipt][0]++;
    }
    if (charge < 0)
    { // Inclusive negatively charged hadrons
      px2Gap[1][ipt][4] += cosnphi;
      py2Gap[1][ipt][4] += sinnphi;
      mpGap[1][ipt][4]++;
    }
    if (pid > 0)
    { // Identified charged hadrons
      px2Gap[1][ipt][pid] += cosnphi;
      py2Gap[1][ipt][pid] += sinnphi;
      mpGap[1][ipt][pid]++;
    }
  } // end of Left TPC subevent selection
  if (eta > fEtaGap)
  { // Right TPC subevent selection
    if (charge > 0)
    { // Inclusive positively charged hadrons
      px2Gap[0][ipt][0] += cosnphi;
      py2Gap[0][ipt][0] += sinnphi;
      mpGap[0][ipt][0]++;
    }
    if (charge < 0)
    { // Inclusive negatively charged hadrons
      px2Gap[0][ipt][4] += cosnphi;
      py2Gap[0][ipt][4] += sinnphi;
      mpGap[0][ipt][4]++;
    }
    if (pid > 0)
    { // Identified charged hadrons
      px2Gap[0][ipt][pid] += cosnphi;
      py2Gap[0][ipt][pid] += sinnphi;
      mpGap[0][ipt][pid]++;
    }
  } // end of Right TPC subevent selection

}

void FlowAnalysisWithQCumulant::ProcessEventAfterFirstTrackLoop(const Int_t &icent)
{
  // 2,QC & 4,QC without eta-gap
  Q2 = TComplex(Qx2, Qy2);
  w2 = M * (M - 1); // w(<2>)
  Q4 = TComplex(Qx4, Qy4);
  w4 = M * (M - 1) * (M - 2) * (M - 3); // w(<4>)
  if (w2 != 0 && w4 != 0)
  {
    cor22 = CalCor22(Q2, M, w2);                   // <2>
    cor24 = CalCor24(Q2, Q4, M, w4);               // <4>
    pCorrelator2->Fill(0.5 + icent, cor22, w2); // <<2>>
    pCorrelator4->Fill(0.5 + icent, cor24, w4); // <<4>>
    // TProfile for covariance calculation in statistic error
    pCov24->Fill(0.5 + icent, cor22 * cor24, w2 * w4); // <2>*<4>
    for (Int_t ipt = 0; ipt < npt; ipt++)
    {
      for (Int_t id = 0; id < npid; id++)
      {
        wred2[ipt][id] = mp[ipt][id] * M - mq[ipt][id];                           // w(<2'>)
        wred4[ipt][id] = (mp[ipt][id] * M - 3 * mq[ipt][id]) * (M - 1) * (M - 2); // w(<4'>)
        if (mp[ipt][id] == 0 || wred2[ipt][id] == 0 || wred4[ipt][id] == 0)
          continue;
        p2[ipt][id] = TComplex(px2[ipt][id], py2[ipt][id]);
        q2[ipt][id] = TComplex(qx2[ipt][id], qy2[ipt][id]);
        q4[ipt][id] = TComplex(qx4[ipt][id], qy4[ipt][id]);
        redCor22[ipt][id] = CalRedCor22(Q2, p2[ipt][id], M, mp[ipt][id], mq[ipt][id], wred2[ipt][id]); // <2'>
        pReducedCorrelator2[id]->Fill(0.5 + ipt, 0.5 + icent, redCor22[ipt][id], wred2[ipt][id]);
        redCor24[ipt][id] = CalRedCor24(Q2, Q4, p2[ipt][id], q2[ipt][id], q4[ipt][id], M, mp[ipt][id], mq[ipt][id], wred4[ipt][id]); // <4'>
        pReducedCorrelator4[id]->Fill(0.5 + ipt, 0.5 + icent, redCor24[ipt][id], wred4[ipt][id]);                                 // <<4'>>

        // TProfile for covariance calculation in statistic error
        pCov22Red[id]->Fill(0.5 + ipt, 0.5 + icent, cor22 * redCor22[ipt][id], w2 * wred2[ipt][id]); // <2>*<2'>
        pCov24Red[id]->Fill(0.5 + ipt, 0.5 + icent, cor22 * redCor24[ipt][id], w2 * wred4[ipt][id]);
        pCov42Red[id]->Fill(0.5 + ipt, 0.5 + icent, cor24 * redCor22[ipt][id], w4 * wred2[ipt][id]);
        pCov44Red[id]->Fill(0.5 + ipt, 0.5 + icent, cor24 * redCor24[ipt][id], w4 * wred4[ipt][id]);
        pCov2Red4Red[id]->Fill(0.5 + ipt, 0.5 + icent, redCor22[ipt][id] * redCor24[ipt][id], wred2[ipt][id] * wred4[ipt][id]);
      }
    }
  }

  // 2-QC, eta-gapped: multi-particle correlation calculation
  if (MGap[0] != 0 && MGap[1] != 0)
  {
    for (Int_t ieta = 0; ieta < neta; ieta++)
    {
      Q2Gap[ieta] = TComplex(Qx2Gap[ieta], Qy2Gap[ieta]);
    }
    w2Gap = MGap[0] * MGap[1];
    cor22Gap = CalRedCor22(Q2Gap[0], Q2Gap[1], MGap[0], MGap[1], 0., w2Gap); // <2>
    pCorrelator2EtaGap->Fill(0.5 + icent, cor22Gap, w2Gap);

    for (Int_t ieta = 0; ieta < neta; ieta++)
    {
      for (Int_t ipt = 0; ipt < npt; ipt++)
      { // <2'>
        for (Int_t id = 0; id < npid; id++)
        {
          if (mpGap[ieta][ipt][id] == 0)
            continue;
          p2Gap[ieta][ipt][id] = TComplex(px2Gap[ieta][ipt][id], py2Gap[ieta][ipt][id]);
          wred2Gap[ieta][ipt][id] = mpGap[ieta][ipt][id] * MGap[ieta];
          redCor22Gap[ieta][ipt][id] = CalRedCor22(Q2Gap[ieta], p2Gap[ieta][ipt][id], MGap[ieta], mpGap[ieta][ipt][id], 0., wred2Gap[ieta][ipt][id]); // <2'>
          pReducedCorrelator2EtaGap[id]->Fill(0.5 + ipt, 0.5 + icent, redCor22Gap[ieta][ipt][id], wred2Gap[ieta][ipt][id]);
          // TProfile for covariance calculation in statistic error
          pCov22RedEtaGap[id]->Fill(0.5 + ipt, 0.5 + icent, cor22Gap * redCor22Gap[ieta][ipt][id], w2Gap * wred2Gap[ieta][ipt][id]);
        }
      }
    }
  }
}

void FlowAnalysisWithQCumulant::SaveHist()
{

  pCorrelator2EtaGap->Write();
  pCorrelator2->Write();
  pCorrelator4->Write();
  for (Int_t id = 0; id < npid; id++)
  {
    pReducedCorrelator2EtaGap[id]->Write();
    pReducedCorrelator2[id]->Write();
    pReducedCorrelator4[id]->Write();
  }

  pCov24->Write();
  for (Int_t id = 0; id < npid; id++)
  {
    pCov22RedEtaGap[id]->Write();
    pCov22Red[id]->Write();
    pCov24Red[id]->Write();
    pCov42Red[id]->Write();
    pCov44Red[id]->Write();
    pCov2Red4Red[id]->Write();
  }
}

TComplex FlowAnalysisWithQCumulant::Qstar(const TComplex &Q)
{
  TComplex QStar = TComplex::Conjugate(Q);
  return QStar;
}

Double_t FlowAnalysisWithQCumulant::CalCor22(const TComplex &Q2, const Double_t &M, const Double_t &w2)
{
  // single-event average 2-particle azimuthal correlation <2>
  Double_t Q2Square = Q2.Rho2();
  Double_t coor22 = Q2Square - M;
  return coor22 / w2;
}

Double_t FlowAnalysisWithQCumulant::CalCor24(const TComplex &Q2, const TComplex &Q4, const Double_t &M, const Double_t &w4)
{
  // single-event average 4-particle azimuthal correlation <4>

  TComplex Q2Star = Qstar(Q2);
  TComplex Q4Star = Qstar(Q4);

  Double_t Q2Square = Q2.Rho2();
  Double_t Q4Square = Q4.Rho2();
  Double_t ReQQQ = (Q4 * Q2Star * Q2Star).Re();

  Double_t coor24 = (Q2Square * Q2Square + Q4Square - 2 * ReQQQ - 4 * (M - 2) * Q2Square + 2 * M * (M - 3));

  return coor24 / w4;
}

Double_t FlowAnalysisWithQCumulant::CalRedCor22(const TComplex &Q2, const TComplex &p2, const Double_t &M, const Double_t &mp,
                                              const Double_t &mq, const Double_t &wred2)
{

  // Calculate the average reduced single-event 2-particle correlations
  TComplex Q2Star = TComplex::Conjugate(Q2);
  Double_t coor22 = (p2 * Q2Star - mq).Re();

  return coor22 / wred2;
}

Double_t FlowAnalysisWithQCumulant::CalRedCor24(const TComplex &Q2, const TComplex &Q4, const TComplex &p2, const TComplex &q2,
                                              const TComplex &q4, const Double_t &M, const Double_t &mp, const Double_t &mq, const Double_t &wred4)
{

  // Calculate the average reduced single-event 2-particle correlations
  TComplex Q2Star = TComplex::Conjugate(Q2);
  TComplex Q4Star = TComplex::Conjugate(Q4);
  TComplex q2Star = TComplex::Conjugate(q2);
  Double_t Q2Square = Q2.Rho2();
  TComplex coorc = p2 * Q2 * Q2Star * Q2Star - q4 * Q2Star * Q2Star - p2 * Q2 * Q4Star - 2.0 * M * p2 * Q2Star - 2.0 * mq * Q2Square + 7.0 * q2 * Q2Star - Q2 * q2Star + q4 * Q4Star + 2.0 * p2 * Q2Star + 2.0 * mq * M - 6.0 * mq;
  Double_t coor24 = coorc.Re();
  return coor24 / wred4;
}