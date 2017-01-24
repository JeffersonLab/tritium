#define NOCASE 1

/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <iomanip>
/*}}}*/
/*ROOT Includes{{{*/
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TNamed.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TError.h>
#include <TVirtualFitter.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TCut.h>
#include <TMultiGraph.h>
#include <TCutG.h>
/*}}}*/

//*ANALYZER Includes{{{*/
//#include <THaAnalysisObject.h>
//#include <THaApparatus.h>
//#include <THaHRS.h>
//#include <THaShower.h>
//#include <THaCherenkov.h>
//#include <THaScintillator.h>
//#include <THaUnRasteredBeam.h>
//#include <THaRasteredBeam.h>
//#include <THaDecData.h>
//#include <THaScalerGroup.h>
//#include <THaReactionPoint.h>
//#include <THaElectronKine.h>
//#include <THaGoldenTrack.h>
//#include <THaExtTarCor.h>
//#include <THaNormAna.h>
//#include <THaAnalyzer.h>
//#include <THaRun.h>
//#include <THaEvent.h>
//*}}}*/
#include "XEM_Target.h"

using namespace std;
  
  
/*Constants{{*/
//const Double_t Na=TMath::Na();           //Avogadro's Number
//const Double_t Qe=TMath::Qe();
const Double_t GeVToMeV = 1000.0; 
const Double_t MeVToGeV = 0.001;
const Double_t PION_MASS = 0.135;          //GeV Mass of Pion
const Double_t E_MASS = 0.510999*MeVToGeV; //Electron Mass
const Double_t P_MASS = 0.938272;          //Proton Mass
const Double_t N_MASS = 0.939565;          //Neutron Mass
const Double_t D_MASS = 1.8756;            //Deutrilium mass
const Double_t PM_SQ = P_MASS*P_MASS;      //:->
const Double_t NM_SQ = N_MASS*N_MASS;      //:->
const Double_t PI_SQ = PION_MASS*PION_MASS;//:->
const Double_t UToGeV = 0.931494028;       //Convert from atomic mass unit to GeV unit
const Double_t Alpha=1/137.036;
const Double_t META = 0.547;
const Double_t HC_SQ = 0.3893857;          //mB/sR, Look like is different with MEV2SR_TO_NBARNSR, FIX_HERE
const Double_t HBARC_XEM=197.33;           //hbar*c=197.33 MeV^2 fm^2
//const Double_t FM2TONBARN=1e-7;          //1 fm^2=1e7 nbarn
//const Double_t MEV2SR_TO_NBARNSR=HBARC*HBARC*FM2TONBARN; //default cross section 1/(MeV^2*sr) to nbarn/sr
  
class XEM_VAR2
{
 public:
  XEM_VAR2()
    {
      SetValue();
    };
  XEM_VAR2(XEM_VAR2 const&){};
  XEM_VAR2& operator=(XEM_VAR2 const&){return *this;};
  XEM_VAR2& operator*=(const Double_t& aNum)
    {
      this->First *= aNum;
      this->Second *= aNum;
      return *this;
    };

  Double_t First; // value of first variable
  Double_t Second;// value of second variable

  void SetValue
    (
     const Double_t& aFirst=0,
     const Double_t& aSecond=0
     )
  {
    First=aFirst;
    Second=aSecond;
  };
    
  void Print()
  {
    Printf("%-*s%f",20,"First Value:  ",First);
    Printf("%-*s%f",20,"Second Value:  ",Second);
  };
};

class XEM_VAR3
{
 public:
  XEM_VAR3()
    {
      SetValue();
    };
  XEM_VAR3(XEM_VAR3 const&){};
  XEM_VAR3& operator=(XEM_VAR3 const&){return *this;};
  XEM_VAR3& operator*=(const Double_t& aNum)
    {
      this->First *= aNum;
      this->Second *= aNum;
      this->Third *= aNum;
      return *this;
    };

  Double_t First; // value of first variable
  Double_t Second;// value of second variable
  Double_t Third;// value of third variable
  
  void SetValue
    (
     const Double_t& aFirst=0,
     const Double_t& aSecond=0,
     const Double_t& aThird=0
     )
  {
    First=aFirst;
    Second=aSecond;
    Third=aThird;
  };
    
  void Print()
  {
    Printf("%-*s%f",20,"First Value:  ",First);
    Printf("%-*s%f",20,"Second Value:  ",Second);
    Printf("%-*s%f",20,"Third Value:  ",Third);
  };
};

class XEM_VAR4
{
 public:
  XEM_VAR4()
    {
      SetValue();
    };
  XEM_VAR4(XEM_VAR4 const&){};
  XEM_VAR4& operator=(XEM_VAR4 const&){return *this;};
  XEM_VAR4& operator*=(const Double_t& aNum)
    {
      this->First *= aNum;
      this->Second *= aNum;
      this->Third *= aNum;
      this->Forth *= aNum;
      return *this;
    };

  Double_t First; // value of first variable
  Double_t Second;// value of second variable
  Double_t Third;// value of third variable
  Double_t Forth; // value of forth variable

  void SetValue
    (
     const Double_t& aFirst=0,
     const Double_t& aSecond=0,
     const Double_t& aThird=0,
     const Double_t& aForth=0
     )
  {
    First=aFirst;
    Second=aSecond;
    Third=aThird;
    Forth=aForth;
  };
    
  void Print()
  {
    Printf("%-*s%f",20,"First Value:  ",First);
    Printf("%-*s%f",20,"Second Value:  ",Second);
    Printf("%-*s%f",20,"Third Value:  ",Third);
    Printf("%-*s%f",20,"Forth Value:  ",Forth);

  };
};


class XEM_SIG
{
 public:
  XEM_SIG()
    {
      SetValue();
    };
  XEM_SIG(XEM_SIG const&){};
  XEM_SIG& operator=(XEM_SIG const&){return *this;};
  XEM_SIG& operator*=(const Double_t& aNum)
    {
      this->Value *= aNum;
      this->Factor *= aNum;
      this->Plus *=aNum;
      return *this;
    };

  Double_t Value;//value of this variable
  Double_t Factor;//correction factor
  Double_t Plus; //Additional output if needed, most of time does not need this
  
  void SetValue
    (
     const Double_t& aValue=0,
     const Double_t& aFactor=0,
     const Double_t& aPlus=0
     )
  {
    Value=aValue;
    Factor=aFactor;
    Plus=aPlus;
  };

  void Print()
  {
    Printf("%-*s%e",20,"Value:  ",Value);
    Printf("%-*s%f",20,"Factor:  ",Factor);
    Printf("%-*s%f",20,"Plus:  ",Plus);
  };
};

class XEM_XS
{
 public:
  XEM_XS()
    {
      SetValue();
    };
  XEM_XS(XEM_XS const&){};
  XEM_XS& operator=(XEM_XS const&){return *this;};
  XEM_XS& operator*=(const Double_t& aNum)
    {
      this->QE *= aNum;
      this->DIS *= aNum;
      /* this->Born *= aNum; */
      /* this->Mott *= aNum; */
      this->All *=aNum;
      this->Factor *=aNum;
      return *this;
    };

  Double_t QE;   // Quasi-Elastic Cross Section
  Double_t DIS;  // Dis Cross Section
  /* Double_t Born; // Born Cross Section */
  /* Double_t Mott; // Mott Cross Section */
  Double_t All;  // Sum of all Cross Section
  Double_t Factor; //Correction Factor
    
  void SetValue
    (
     const Double_t& aQE=0,
     const Double_t& aDIS=0,
     /* const Double_t& aBorn=0, */
     /* const Double_t& aMott=0, */
     const Double_t& aAll=0,
     const Double_t& aFactor=0.
     )
  {
    QE = aQE;
    DIS = aDIS;
    /* Born = aBorn; */
    /* Mott = aMott; */
    All = aAll; 
    Factor = aFactor;
  };

  void Print()
  {
    Printf("%-*s%e",20,"QE:  ",QE);
    Printf("%-*s%e",20,"DIS:  ",DIS);
    /* Printf("%-*s%e",20,"Born:  ",Born); */
    /* Printf("%-*s%e",20,"Mott:  ",Mott); */
    Printf("%-*s%e",20,"All:  ",All);
    Printf("%-*s%f",20,"Factor:  ",Factor);
  };
};
  


inline Double_t gGet_Max(double a, double b){
  if(a>=b)
    return a;
  else if(a<b)
    return b;
  else{
    cerr <<"****** Can not make MAX(a,b) right! Check a and b !!!"<<endl;
    return 0;
  } 
}


//Subroutines in XEM_SigCal.h
inline Double_t gGet_EMC_Func(Int_t aA, Double_t aXbj);
inline Double_t gGet_EMC_Func_Slac(Int_t aA, Double_t aXbj);
inline Double_t gGet_Dis_HighX_Cor(Int_t aA, Double_t aXbj);
 
//Subroutines in XEM_Fy2Sig.h
inline Double_t gGet_Y(double vE0, double vEp, double vTheta, XEM_TGT* vTarget);
inline void gGet_Sig_Bar_DF(double vE0,double vEp,double vTheta, double vY,double vPt,XEM_VAR2* sig_pn );
inline Double_t gGet_Fy(double vY, XEM_TGT* vTarget);
inline Double_t gGet_He3_Fy(double vY, XEM_TGT* vTarget);
inline void gGet_FF(int vFlag, double vQsq, double* vFF);
inline Double_t gGet_Tail_Corr(Int_t kA,Double_t kX);

//Subroutines in XEM_F1F2.h
inline Double_t gGet_Resmod316(double aQsq, double aWsq, int aflag);
//inline XEM_VAR4* gGet_Christy806(double aQsq, double aWsq);
inline void gGet_Christy806(double aQsq, double aWsq,XEM_VAR4 *sig);
inline Double_t gGet_NMC(double aXbj, double aQsq);
inline Double_t gGet_Resmod_Hack(double aQsq,double aWsq, double *aXval);
inline Double_t gGet_Resd(double aQsq, double aWsq, double *aXval);
inline void gGet_Pind(double aQsq, double aWsq,XEM_VAR4* sig);
 
