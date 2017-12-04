#define NOCASE 1
/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
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
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
//#include <TMatrix.h>
/*}}}*/
//*analyzer Includes{{{*/
//#include <THaAnalysisObject.h>
//#include <THaApparatus.h>
//#include <THaHRS.h>
//#include <THaShower.h>
//#include <THaCherenkov.h>
//#include <THaScintillator.h>
#include <THaUnRasteredBeam.h>
//#include <THaRasteredBeam.h>
//#include <THaDecData.h>
//#include <THaScalerGroup.h>
//#include <THaReactionPoint.h>
//#include <THaElectronKine.h>
//#include <THaGoldenTrack.h>
//#include <THaExtTarCor.h>
#include <THaNormAna.h>
#include <THaAnalyzer.h>
#include <THaRun.h>
#include <THaEvent.h>

//*}}}*/
using namespace std;

/*Customization{{{*/
const TString TRI_DIR="/work/halla/triton//yez/SRCs/Cross_Sections/"; //this is for inputfiles,outfiles,logfiles
const TString ROOTFILES_DIR="/volatile/halla/triton/Tritium_Rootfiles/";     //this is for rootfiles
const TString ROOTFILES_NAME = "tritium"; //such as tritium.root
const TString SAMCFILES_DIR="/volatile/halla/tritium/yez/Monte_Calo/rootfiles/";    //this is for SAMC Rootfiles
const TString XEMC_Table_Path="/work/halla/tritium/yez/Monte_Calo/Cross_Section_Tables/"; 
TString D_Name="D_N_";

ofstream outlog; 
Double_t Target_Length = 0.0;
Double_t Target_Length_Default = 0.25; // Default length for tritium target cells are 25cm;
/*}}}*/

Bool_t IsBackground=kFALSE;
/*Constants and enum{{{*/
const Double_t Pi=TMath::Pi();
const Double_t Na=TMath::Na();             //Avogadro's Number
const Double_t Qe=TMath::Qe();
const Int_t    DQ_Name_Length=60;         //for DataQuality Output Name Length
const Double_t CM2ToNB = 1e33;//Convert units of cm^2 into nbarn, 1 Barn=1e28*M^2=1e24*CM^2, 1 Barn = 1e9 nBarn   
//All Energy Unit are in MeV, all Length Units are in cm!
const Double_t GeVToMeV=1000.00;
const Double_t MeVToGeV=0.001;
const Double_t CM2M=0.01;
const Double_t PROTON_MASS=0.938272; //GeV Mass of proton

enum E_INPUT_VARS_Def
{
	k_INPUT_RUNNO,
	k_INPUT_ARM,
	k_INPUT_SPEC_SET_MOM,                   //GeV
	k_INPUT_TARGET_NAME,
	k_INPUT_TARGET_THICKNESS,              //density*Length kg/m^2
	k_INPUT_TARGET_THICKNESS_STAT_ERR,     //density*Length kg/m^2
	k_INPUT_TARGET_A,                      //Atomic mass kg/mol
	k_INPUT_TARGET_D,                      //target density kg/m^3
	k_INPUT_ISEXTGT,                       //1=extended target 0=not
	k_INPUT_MAIN_TRIGGER,
	k_INPUT_S_TREE_NAME,
	k_INPUT_RUNNOCHAIN
};
/*}}}*/

const TString Left_Trigger_Name = "LeftT_N_";//will replace by LeftT1, LeftT2 ...
const TString Right_Trigger_Name = "RightT_N_";//will replace by RightT1, RightT2 ...

//////////////////////////////////////////
// BCM Constant, needed by TRI_Beam.h
//////////////////////////////////////////
//Need to update the names of these variables
const double Clock_Freq = 103700;
const TString LEFT_CLOCK = "LeftLclock";
const TString RIGHT_CLOCK= "RightRclock";
const TString LEFT_BCM_U = "Leftu1";
const TString LEFT_BCM_D = "Leftd1";
const TString RIGHT_BCM_U= "Rightu1";
const TString RIGHT_BCM_D= "Rightd1";

/*BCM{{{*/
// Note: Only Left U1,D1 are calibrated, but D1 had problem during run# 3835~3888.
//       Right BCM looks like having problems.
const Double_t u_Coeff[2][3]=
  {
    //[2]
    //0=Left 1=Right
    //[3]
    //0=1, 1=3, 2=10 u1,u3,u10
//Alpha: U1    U3    U10         //Updated in 07/20/2012
    {2070.0,    0.0, 0.0},            //u Coeff to calculate Charge for left arm
    {1035.0, 3178.0, 0.0},            //u Coeff to calculate Charge for right arm
  };
//                 Offset        D1     D3    D10        //Updated in 07/20/2012
const Double_t u_Offset[2][3]={{400.0,   0.0, 0.0},    //the offsets of the u1,u3,u10 for left arm
			       {200.0, 174.0, 0.0},};  //the offsets of the u1,u3,u10 for right arm
const Double_t d_Coeff[2][3]=
  {
    //[2]
    //0=Left 1=Right
    //[3]
    //0=1, 1=3, 2=10 u1,u3,u10
   
//Alpha: D1    D3    D10        //Updated in 07/20/2012
    {1263.0,    0.0, 0.0},      //d Coeff to calculate Charge for left arm
    {1263.0, 3902.0, 0.0},      //d Coeff to calculate Charge for right arm
  };
//                 Offset        D1     D3  D10       //Updated in 07/20/2012
const Double_t d_Offset[2][3]={{20.0,  0.0, 0.0},     //the offsets of the d1,d3,d10 for left arm
			       {37.0, 22.0, 0.0},};   //the offsets of the d1,d3,d10 for right arm
/*}}}*/
