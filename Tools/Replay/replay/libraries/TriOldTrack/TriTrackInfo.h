#ifndef ROOT_TriTrackInfo
#define ROOT_TriTrackInfo

//////////////////////////////////////////////////////////////////////////
//
// TriTrackInfo
//
//////////////////////////////////////////////////////////////////////////

#include "TVector3.h"
class TriTrack;
class TriSpectrometer;

class TriTrackInfo {
public:
  TriTrackInfo() : fP(kBig), fDp(kBig), fX(kBig), fY(kBig), fTheta(kBig), 
    fPhi(kBig), fPvect(kBig,kBig,kBig), fOK(0), fSpectro(0) {}

  TriTrackInfo( Double_t p, Double_t dp, Double_t x, Double_t y, Double_t th,
		Double_t ph, Double_t px, Double_t py, Double_t pz ) :
    fP(p), fDp(dp), fX(x), fY(y), fTheta(th), fPhi(ph), fOK(1), fSpectro(0)
  { fPvect.SetXYZ(px,py,pz); }

  TriTrackInfo( Double_t p, Double_t dp, Double_t x, Double_t y, Double_t th,
		Double_t ph, const TVector3& pvect ) :
    fP(p), fDp(dp), fX(x), fY(y), fTheta(th), fPhi(ph), fPvect(pvect), fOK(1),
    fSpectro(0)
  {}

  TriTrackInfo( const TriTrackInfo& t ) :
    fP(t.fP), fDp(t.fDp), fX(t.fX), fY(t.fY), fTheta(t.fTheta), 
    fPhi(t.fPhi), fPvect(t.fPvect), fOK(t.fOK), fSpectro(t.fSpectro) {}
  TriTrackInfo& operator=( const TriTrackInfo& );
  TriTrackInfo& operator=( const TriTrack& );
  virtual ~TriTrackInfo() {}

  void      Clear( Option_t* opt="" );

  Bool_t    IsOK()  const { return fOK; }
  Double_t  GetPx() const { return fPvect.X(); }
  Double_t  GetPy() const { return fPvect.Y(); }
  Double_t  GetPz() const { return fPvect.Z(); }

  Double_t  GetP()     const { return fP; }
  Double_t  GetDp()    const { return fDp; }
  Double_t  GetX()     const { return fX; }
  Double_t  GetY()     const { return fY; }
  Double_t  GetTheta() const { return fTheta; }
  Double_t  GetPhi()   const { return fPhi; }
  const TVector3& GetPvect() const { return fPvect; }

  TriSpectrometer* GetSpectrometer() const { return fSpectro; }

  void Set( Double_t p, Double_t dp, Double_t x, Double_t y, Double_t th,
	    Double_t ph, Double_t px, Double_t py, Double_t pz )
  { 
    fP = p; fDp = dp; fX = x; fY = y; fTheta = th; fPhi = ph; 
    fPvect.SetXYZ(px,py,pz); fOK = 1; 
  }
  void Set( Double_t p, Double_t dp, Double_t x, Double_t y, Double_t th,
	    Double_t ph, const TVector3& pvect )
  { 
    fP = p; fDp = dp; fX = x; fY = y; fTheta = th; fPhi = ph; fPvect = pvect; 
    fOK = 1;
  }
  void SetP( Double_t p );
  void SetSpectrometer( TriSpectrometer* s ) { fSpectro = s; }

protected:
  Double_t  fP;         // Track momentum (GeV)
  Double_t  fDp;        // TRANSPORT delta
  Double_t  fX;         // TRANSPORT x at target (m)
  Double_t  fY;         // TRANSPORT y at target (m)
  Double_t  fTheta;     // Tangent of TRANSPORT theta at target
  Double_t  fPhi;       // Tangent of TRANSPORT phi at target
  TVector3  fPvect;     // Momentum vector in lab (GeV)
  Int_t     fOK;        // Data ok (0:no 1:yes)

  TriSpectrometer* fSpectro; //! Spectrometer creating this track

private:
  static const Double_t kBig;

  ClassDef(TriTrackInfo,1)  // Track information
};

//_____________________________________________________________________________
inline
void TriTrackInfo::Clear( Option_t* )
{ 
  fPhi = fTheta = fY = fX = fDp = fP = kBig; fPvect.SetXYZ(kBig,kBig,kBig); 
  fOK = 0; 
}

//_____________________________________________________________________________
inline
void TriTrackInfo::SetP( Double_t p )
{
  // Set momentum to new value. Changes p, dp, and pvect

  Double_t p0 = fP / (1.0+fDp);
  fP = p;
  fPvect.SetMag(p);
  if( p0 != 0.0 )
    fDp = fP/p0 - 1.0;
}

#endif

