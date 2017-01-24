#define NOCASE 1

/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
/*}}}*/
/*ROOT Includes{{{*/
#include <TSystem.h>
#include <TString.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TDatime.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
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
///*ANALYZER Includes{{{*/
////#include <THaAnalysisObject.h>
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
///*}}}*/
using namespace std;

/*Customization{{{*/
const TString XGT2_DIR="/work/halla/e08014/disk1/yez/HRS_Cali/XGT2"; //this is for inputfiles,outfiles,logfiles
const TString ROOTFILES_DIR="/work/halla/e08014/disk1/Rootfiles";//this is for rootfiles
const Double_t GCPMT_For_1Photon=100;

/*}}}*/

Bool_t IsBackground=kFALSE;

/*Constants and enum{{{*/
const Double_t Pi=TMath::Pi();
const Double_t Na=TMath::Na();//Avogadro's Number
const Double_t Qe=TMath::Qe();
const Double_t DegToRad=TMath::DegToRad();
const Double_t RadToDeg=TMath::RadToDeg();
const Double_t PROTON_MASS=0.938272; //GeV Mass of proton
const Double_t Alpha=1/137.036;
const Int_t    DQ_Name_Length=60;//for DataQuality Output Name Length
//Form("%-*s%-f",DQ_Name_Length,"Name",Value);
const Double_t HBARC=197.33; //hbar*c=197.33 MeV^2 fm^2
const Double_t FM2TONBARN=1e-7; //1 fm^2=1e7 nbarn
const Double_t MEV2SR_TO_NBARNSR=HBARC*HBARC*FM2TONBARN; //default cross section 1/(MeV^2*sr) to nbarn/sr

///////////////////////////////////
/// PreScale Factors
///////////////////////////////////

Int_t gGet_PS(const Int_t& RunNo,const TString& aArm, Int_t *aPar )
{
  TString filename_ps = ROOTFILES_DIR+Form("e08014_%s_%d.root",aArm,RunNo);
  TFile *file_ps = new TFile(filename_ps);
  TArrayI prescales = Run_Data->GetParameters()->GetPrescales();
  for(int i=0;i<8;i++){
    aPar[i] = prescales[i];
  }
  filename_ps->Close();  return &aPar;  
}

///////////////////////////////////
/// Energy Loss
///////////////////////////////////

/*inline Double_t gCal_Eloss(const Double_t& aE,const Int_t& aZ,const Double_t& aA,const Double_t& aM,const Double_t& aI,const Double_t& aDensity,const Double_t& aThickness){{{*/
inline Double_t gCal_Eloss(const Double_t& aE,const Int_t& aZ,const Double_t& aA,const Double_t& aM,const Double_t& aI,const Double_t& aDensity,const Double_t& aThickness)
{
  //aE,aM,aI GeV, aA kg/mol, aDensity g/cm^3
  Double_t lE=aE*1000; //GeV to MeV
  Double_t lM=aM*1000; //GeV to MeV
  Double_t lI=aI*1000; //GeV to MeV
  Double_t lA=aA*1000; //kg/mol to g/mol
  Double_t lThickness=aThickness/10; //kg/m^2 to g/cm^2

  Double_t K=0.307075; //MeV cm^2/g
  Double_t z=1;//electron
  Double_t me=0.511; //MeV
  Double_t beta,gamma,omega_p,Tmax,delta;

  beta=1-me*me/(lE*lE);
  gamma=1/sqrt(1-beta*beta);

  omega_p=28.816*sqrt(aDensity*aZ/lA/1000)*1e-6; //MeV
  Tmax=2*me*beta*beta*gamma*gamma/(1+2*gamma*me/lM+me*me/(lM*lM));
  delta=log(omega_p/lI)+log(beta*gamma)-0.5;
  return (K*z*z*aZ/lA/(beta*beta)*( 0.5*log(2*me*beta*beta*gamma*gamma*Tmax/lI/lI)-beta*beta-delta ))*lThickness*1e-3;//GeV
}
/*}}}*/




/*class XGT2_VAR{{{*/
class XGT2_VAR
{
 public:
  XGT2_VAR()
    {
      SetValue();
    };
  XGT2_VAR(XGT2_VAR const&){};
  XGT2_VAR& operator=(XGT2_VAR const&){};
  XGT2_VAR& operator*=(const Double_t& aNum)
    {
      this->Value *= aNum;
      this->Stat_Err *= aNum;
      this->Sys_Err *= aNum;
      return *this;
    };
  
  Double_t Value;//value of this variable
  Double_t Stat_Err;//statistical error of this variable
  Double_t Sys_Err;//systematic error of this variable
  
  void SetValue
    (
     const Double_t& aValue=0,
     const Double_t& aStat_Err=0,
     const Double_t& aSys_Err=0
     )
  {
    Value=aValue;
    Stat_Err=aStat_Err;
    Sys_Err=aSys_Err;
  };
  
  void Print()
  {
    Printf("%-*s%e",20,"Value:",Value);
    Printf("%-*s%e",20,"Stat_Err:",Stat_Err);
    Printf("%-*s%e",20,"Sys_Err:",Sys_Err);
  }
};
/*}}}*/


//////////////////////////
// Logger
/////////////////////////

/*class XGT2_Logger{{{*/
class XGT2_Logger
{
 public:
  /*static XGT2_Logger* Instance(){{{*/
  static XGT2_Logger* Instance()
  {
    // Only allow one instance of class to be generated.
    if (!m_pInstance)
      {
	Printf("Creating XGT2_Logger System.");
	m_pInstance = new XGT2_Logger;
      }
    return m_pInstance;
  }
  /*}}}*/
		/*static XGT2_Logger* Destroy(){{{*/
  static XGT2_Logger* Destroy()
  {
    // Only allow one instance of class to be generated.
    if ( m_pInstance)
      delete m_pInstance;
    Printf("Deleting XGT2_Logger System.");
    m_pInstance=NULL;
  }
  /*}}}*/
  /*const TString GetLogFileName(){{{*/
  const TString GetLogFileName()
  {
    return LogFile_Name;
  }
  /*}}}*/
  
  void SetLogFileName(const char* aLogFile_Name)
  {
    LogFile_Name=aLogFile_Name;
  }
  
  
  /*void WriteToLogFile(const char* aSrcFile, const char* aSrcLine,const char* aLog){{{*/
  void WriteToLogFile(const char* aSrcFile, int aSrcLine,const char* aLog)
  {
    ofstream outfile(LogFile_Name.Data(),ios::app);
    if ( !outfile.good() )
      {
	Error(__FILE__,Form("At line %d. %s cannot be appended.",__LINE__,LogFile_Name.Data()));
	return;
      }
    TDatime stamp;
    outfile<<Form("%s File-\"%s\" Line-%d: %s",stamp.AsString(),aSrcFile,aSrcLine,aLog)<<endl;
    outfile.close();
  }
  /*}}}*/
  
  
 private:
  XGT2_Logger(){};  // Private so that it can  not be called
  XGT2_Logger(XGT2_Logger const&){};             // copy constructor is private
  XGT2_Logger& operator=(XGT2_Logger const&){};  // assignment operator is private
  static XGT2_Logger* m_pInstance;
  TString LogFile_Name;
}
;/*}}}*/

XGT2_Logger* XGT2_Logger::m_pInstance = NULL;
