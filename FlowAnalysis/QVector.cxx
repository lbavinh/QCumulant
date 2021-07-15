#include <QVector.h>

ClassImp(QVector);

QVector::QVector() : fNHarmonic(2), fQx(0), fQy(0), fMult(0), fWeight(0) 
{
}

QVector::QVector(Double_t n) : fNHarmonic(n), fQx(0), fQy(0), fMult(0), fWeight(0) 
{
}

QVector::~QVector()
{
}

void QVector::Zero()
{
  fQx = 0.;
  fQy = 0.;
  fMult = 0.;
  fWeight = 0.;
}

void QVector::CalQVector(const Double_t &phi, const Double_t &weight)
{
  fQx += weight * TMath::Cos( fNHarmonic * phi );
  fQy += weight * TMath::Sin( fNHarmonic * phi );
  fWeight += weight;
  fMult++;
}

void QVector::WeightQVector()
{
  if (fMult == 0) { cout << "Warning! fMult==0!" << endl;}
  if (fMult != 0)
  {
    if (fWeight<0) fWeight *= -1.;
    fQx /= fWeight;
    fQy /= fWeight;
  }
}