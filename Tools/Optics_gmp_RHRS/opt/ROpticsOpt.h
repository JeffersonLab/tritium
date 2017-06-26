///////////////////////////////////////////////////////////////////////////////
//
// ROpticsOpt
//
// HRS optics matrix optimization class
// Based on THaVDC
//
// Units used:
//   For X, Y, and Z coordinates of track    -  meters
//   For Theta and Phi angles of track       -  tan(angle)
//   For Momentums, Masses                   -  GeV, GeV/c^2
//
// Author: Jin Huang <jinhuang@jlab.org>
//
// Modification:
//   Jun 25, 2010 Updated for APEX optics calibration
//   Aug 01, 2013 Updated for G2P optics calibration (Chao Gu)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_ROpticsOpt
#define ROOT_ROpticsOpt

#include <vector>

#include "TRotation.h"
#include "TMath.h"

#include "THaTrackingDetector.h"
#include "THaString.h"

///////////////////////////////////////////////////////////////////////////////
// Debug Definitions
// place this section below any other head files
///////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG_LEVEL
#undef DEBUG_LEVEL
#endif
//     DEBUG_LEVEL;
//     =0      or not define: no debug, full speed
//     >=1     enable debug extra warning (suggested setting)
//     >=2     above + enable debug assert
//     >=3     above + enable debug extra info
//     >=4     above + massive info (in a for or while)
//     >=5     Decode dump
#define DEBUG_LEVEL 3
#include "DebugDef.h"
#include "InputR.h"

///////////////////////////////////////////////////////////////////////////////
class TCanvas;
class THaTrack;
class TClonesArray;
class TTree;
class TVector3;

class THaMatrixElement;

class ROpticsOpt : public THaTrackingDetector {
public:
    ROpticsOpt(const char* name = "Optimizer", const char* description = "Optimizer for HRS Optics", THaApparatus* a = NULL);
    virtual ~ROpticsOpt();

    ///////////////////////////////////////////////////////////////////////////
    // Database input/output
    ///////////////////////////////////////////////////////////////////////////
    TString OldComments;
    Int_t LoadDataBase(TString DataBaseName); // Database file -> Memory
    Int_t SaveDataBase(TString DataBaseName); // Memory -> Database file

    virtual void Print(const Option_t* opt) const;

    UInt_t Matrix2Array(Double_t Array[], Bool_t FreeParaFlag[] = NULL) // fCurrentMatrixElems -> Array
    {
        assert(fCurrentMatrixElems);
        return Matrix2Array(Array, (*fCurrentMatrixElems), FreeParaFlag);
    }
    UInt_t Matrix2Array(Double_t Array[], const std::vector<THaMatrixElement> &Matrix, Bool_t FreeParaFlag[] = NULL);

    UInt_t Array2Matrix(const Double_t Array[]) // Array -> fCurrentMatrixElems
    {
        assert(fCurrentMatrixElems);
        return Array2Matrix(Array, (*fCurrentMatrixElems));
    }
    UInt_t Array2Matrix(const Double_t Array[], std::vector<THaMatrixElement> &Matrix);

    ///////////////////////////////////////////////////////////////////////////
    // Data storage
    ///////////////////////////////////////////////////////////////////////////
    
    enum {
        MaxNEventData = 50, MaxNRawData = 2000000, kNUM_PRECOMP_POW = 10, kMaxDataGroup = 180 * 5 * 5
    };
    
    
    
    UInt_t LoadRawData(TString DataFileName, UInt_t NLoad = MaxNRawData, UInt_t MaxDataPerGroup = (UInt_t) - 1); // load data to Rawdata[]

    
    //typedef struct {
    struct EventData {
        Double_t Data[MaxNEventData]; // [CommonIdx]
        Double_t powers[kNUM_PRECOMP_POW][5]; // {(x), th, y, ph, abs(th) }
    };
    //} EventData;
    EventData fRawData[MaxNRawData]; // [fNRawData]
    UInt_t fNRawData;
    UInt_t fNCalibData; // for dp calib only

    enum CommonIdx {
        kCutID = 0, // cut ID in order of tree2ascii cut file
        kX = 1, // R.tr.r_x
        kTh = 2, // R.tr.r_th
        kY = 3, // R.tr.r_y
        kPhi = 4, // R.tr.r_ph
        kurb_e = 5, // Beam energy
        kBeamX = 6, // urb.x or rb.x
        kBeamY = 7, // urb.y or rb.y
	kBeamVZ = 8, //R.tr.vz
    };

    enum ExtraSieveIdx {
        kRealTh = 30, // real target th from survey
        kRealPhi, // real target ph from survey
        kRealTgX, // real target x from survey, beam
        kRealTgY, // real target y from survey, beam
        kRealThMatrix, // expected target th before extended target corrections
        kCalcTh, // calculated th from matrix
        kCalcPh, // calculated ph from matrix
        kSieveX,
        kSieveY,
        kSieveZ,
        kBeamZ
    };

    enum ExtraVertexIdx {
        kRealReactZ = 10, //expected ReactZ
        kCalcTgY, //calculated Tg_y
        kCalcReactZ //calculated ReactZ
    };

    enum ExtraDpIdx {
        kExtraDataFlag = 15, //Whether this event is for optimization; 0=used for optimization, 1=for plotting only
        kKineID, //Delta Scan Kinematics ID
        kCentralp, //Central Momentum
        kElossDp, //Radiation Loss for this event in unit of dp
        kScatterAngle, //Scattering Angle
        kDpKinOffsets, //=dp-dp_kin, based on sieve hole survey
        kRealDpKin, //expected dp_kin, before radiation correction
        kRealDpKinMatrix, //real dp kin before extended target corrections
        kCalcDpKinMatrix, //calculated dp kin before extended target corrections
        kCalcDpKin, //calculated dp_kin, before radiation correction
        kRealDpKinExcitations,/*Do not append more index*/ //first index of expected dp_kins for all excitation states
	kElossTgBefore,
	kElossTgAfter,
	kTravelLengthBefore,
	kTravelLengthAfter
    };

    ///////////////////////////////////////////////////////////////////////////
    // Optimization related Commands
    ///////////////////////////////////////////////////////////////////////////
    const TVector3 GetSieveHoleTCS(UInt_t Col, UInt_t Row);
    const TVector3 GetSieveHoleCorrectionTCS(UInt_t nfoil, UInt_t Col, UInt_t Row);

    Double_t TravelLength(TVector3 ReactionVertex, TVector3 MomDirectionHCS);
    //    Double_t ElossTarget(TVector3 ReactionVertex, TVector3 MomDirectionHCS);
    Double_t ElossTargetBefore(TVector3 ReactionVertex, TVector3 MomDirectionHCS);
    Double_t ElossTargetAfter(TVector3 ReactionVertex, TVector3 MomDirectionHCS);
    void PrepareSieve(void);
    TCanvas* CheckSieve(Int_t PlotFoilID = 0);
    TCanvas* CheckSieveAccu(Int_t PlotFoilID = 0);
    Double_t SumSquareDTh(void);
    Double_t SumSquareDPhi(void);

    Double_t fArbitaryVertexShift[100]; // compensate bias due to event selections, array of [FoilID]
    void PrepareVertex(void);
    TCanvas* CheckVertex(void);
    Double_t SumSquareDTgY();

    Double_t fArbitaryDpKinShift[100]; // compensate bias due to dp event selections, array of [KineID]
    void PrepareDp(void);
    TCanvas* CheckDp(void);
    TCanvas* CheckDpGlobal(void);
    Double_t SumSquareDp(Bool_t IncludeExtraData = kFALSE);

    TRotation fTCSInHCS; // transformations vector from TCS to HCS
    TVector3 fPointingOffset; // Optical point in lab coordinate system

    inline Double_t ScatMom(Double_t DM, Double_t Ma, Double_t P0, Double_t Theta)
    {
        // Calc Scattered Electron Momentum
        // Assuming Electron with P0, Scattering on a fix target of Mass Ma Assuming
        // recoil system is a resonance of DM and scattering angle is Theta.
        return (-DM * DM - 2 * Ma * DM + 2 * Ma * P0) / (2 * (Ma + P0 - P0 * TMath::Cos(Theta)));
    }

    ///////////////////////////////////////////////////////////////////////////
    // declarations for target vertex reconstruction
    ///////////////////////////////////////////////////////////////////////////

    enum ECoordTypes {
        kTransport, kRotatingTransport
    };

    enum EFPMatrixElemTags {
        T000 = 0, Y000, P000
    };

    enum {
        kPORDER = 7
    };

    friend class THaMatrixElement;
    std::vector<THaMatrixElement> * fCurrentMatrixElems;
    // initial matrix elements
    std::vector<THaMatrixElement> fTMatrixElems;
    std::vector<THaMatrixElement> fDMatrixElems;
    std::vector<THaMatrixElement> fPMatrixElems;
    std::vector<THaMatrixElement> fPTAMatrixElems; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fYMatrixElems;
    std::vector<THaMatrixElement> fYTAMatrixElems; // involves abs(theta_fp)
    std::vector<THaMatrixElement> fFPMatrixElems; // matrix elements used in

    std::vector<THaMatrixElement> fRMatrixElems; // Path-length corrections (meters)

    void CalcMatrix(const double x, std::vector<THaMatrixElement> &matrix);
    // Double_t DoPoly(const int n, const std::vector<double> &a, const double x);
    // Double_t PolyInv(const double x1, const double x2, const double xacc, const double y, const int norder, const std::vector<double> &a);
    Double_t CalcTargetVar(const std::vector<THaMatrixElement> &matrix, const double powers[][5]);

    ///////////////////////////////////////////////////////////////////////////
    // Inherited from THaTrackingDetector
    ///////////////////////////////////////////////////////////////////////////

    virtual Int_t Decode(const THaEvData&)
    {
        return 0;
    }

    virtual Int_t CoarseTrack(TClonesArray&)
    {
        return 0;
    }

    virtual Int_t FineTrack(TClonesArray&)
    {
        return 0;
    }

    virtual EStatus Init(const TDatime&)
    {
        return fStatus = kOK;
    };

private:
    ClassDef(ROpticsOpt, 0) // HRS Optics Optimizer
};

///////////////////////////////////////////////////////////////////////////////
// class for storing matrix element data
///////////////////////////////////////////////////////////////////////////////

class THaMatrixElement {
public:

    THaMatrixElement() : iszero(true), pw(3), order(0), v(0), poly(ROpticsOpt::kPORDER), OptOrder(0)
    {
    }
    bool match(const THaMatrixElement& rhs) const;

    bool iszero; // whether the element is zero
    std::vector<int> pw; // exponents of matrix element
    //   e.g. D100 = { 1, 0, 0 }
    int order;
    double v; // its computed value
    std::vector<double> poly; // the associated polynomial

    void SkimPoly(); //reduce order to highest non-zero poly

    UInt_t OptOrder; //order optimize to
};
#endif
