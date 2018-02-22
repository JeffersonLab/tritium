#ifndef ROOT_TriFadcBPM
#define ROOT_TriFadcBPM

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TriFadcBPM                                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaBeamDet.h"
#include "TVectorT.h"

class TriFadcBPM : public THaBeamDet {

public:
  TriFadcBPM( const char* name, const char* description = "",
		   THaApparatus* a = NULL );
  virtual ~TriFadcBPM();

  virtual void       Clear( Option_t* ="" );
  virtual Int_t      Decode( const THaEvData& );
  virtual Int_t      Process();

  virtual TVector3 GetPosition()  const { return fPosition; }
  virtual TVector3 GetDirection()  const { return fDirection; }

  Double_t GetRawSignal0() {return fRawSignal(0);}
  Double_t GetRawSignal1() {return fRawSignal(1);}
  Double_t GetRawSignal2() {return fRawSignal(2);}
  Double_t GetRawSignal3() {return fRawSignal(3);}

  Double_t GetRotPosX() {return fRotPos(0); }
  Double_t GetRotPosY() {return fRotPos(1); }


protected:

  virtual Int_t  ReadDatabase( const TDatime& date );
  virtual Int_t  DefineVariables( EMode mode = kDefine );

  //  TriFadcBPM() {}
  //  TriFadcBPM( const TriFadcBPM& ) {}
  //  TriFadcBPM& operator=( const TriFadcBPM& ) { return *this; }


  TVectorD  fRawSignal;     // induced signal of the antennas
  TVectorD  fPedestals;
  TVectorD  fCorSignal;     // pedestal subtracted signal

  TVectorD  fRotPos;        // position in the BPM system, arbitrary units
  
  TMatrix fRot2HCSPos;     
  
  TVector3  fOffset;
  TVector3  fPosition;   // Beam position at the BPM (meters)
  TVector3  fDirection;  // Beam direction at the BPM
                         // always points along z-axis

  Int_t fNfired;
  Double_t fCalibRot;

  ClassDef(TriFadcBPM,0)   // Generic BPM class
};

////////////////////////////////////////////////////////////////////////////////

#endif
