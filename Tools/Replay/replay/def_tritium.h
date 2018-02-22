//////////////////////////////////////////////////////////////////////////
//
// Bigbite lib replay script definition header
//
// https://hallaweb.jlab.org/wiki/index.php/ReplayCore
//
// It is suggested to use this file set replay environmental variable for 
// different machine. So, they can use same replay script
//	
//////////////////////////////////////////////////////////////////////////
//	
//	Author : Jin Huang (jinhuang@jlab.org)    Aug 2007
//	Modify History:
//
//      27 Mar 08 ole  Change ROOT file names for DetL and PhysL replays
//                     from _det_nnnn to _det_L_nnnn etc.
//
//////////////////////////////////////////////////////////////////////////
//
//  To Do:
//  - Raw File Format - change to tritium
//  - File path names
//
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_REPLAY_DEF
#define ROOT_REPLAY_DEF


//////////////////////////////////////////////////////////////////
// extra structures
//////////////////////////////////////////////////////////////////

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

//*ANALYZER Includes{{{*/
#include <THaAnalysisObject.h>
#include <THaApparatus.h>
#include <THaHRS.h>
#include <THaVDC.h>
#include <THaShower.h>
#include <THaCherenkov.h>
#include <THaScintillator.h>
#include <THaUnRasteredBeam.h>
#include <THaRasteredBeam.h>
#include <THaRaster.h>
#include <THaIdealBeam.h>
#include <THaDecData.h>
#include <THaReactionPoint.h>
#include <THaElectronKine.h>
#include <THaGoldenTrack.h>
#include <THaExtTarCor.h>
#include <THaAnalyzer.h>
#include <THaRun.h>
#include <THaBPM.h>
#include <THaEvent.h>
//*}}}*/


#include "TObject.h" //trick to avoid error: `Int_t' does not name a type
#include "TString.h"

//Used for setting paths. %s is necessary so that it can be used to Form the paths.
//This allows changing the directory in fewer places for ease of portability.
const char* REPLAY_DIR_PREFIX = "./%s";
const char* ROOTFILE_DIR_PREFIX = "./%s";

typedef struct _sReplaySetUp
{
    Int_t DefReplayNum; //defaut replay event num
    char* OutFileFormat;//output file format
    char* OutDefineFile;//out define
    char* CutDefineFile;//cut define
}sReplaySetUp;


//////////////////////////////////////////////////////////////////
// common defines
//////////////////////////////////////////////////////////////////
static const char* PATHS[] = {
  "/cache/halla/triton/raw",
  "/adaq1/data1",
   0
};

static const char* RAW_DATA_FORMAT="%s/triton_%d.dat.%d";
//static const char* RAW_DATA_FORMAT="%s/dvcs14_%d.dat.%d";

TString STD_REPLAY_OUTPUT_DIR= TString::Format(ROOTFILE_DIR_PREFIX,"Rootfiles");
TString CUSTOM_REPLAY_OUTPUT_DIR= TString::Format(ROOTFILE_DIR_PREFIX,"ScratchROOTfiles");
TString SUMMARY_PHYSICS_FORMAT= TString::Format(REPLAY_DIR_PREFIX,"summaryfiles/summaryphy_%d.log");

//Previously used:
//static const char* STD_REPLAY_OUTPUT_DIR=Form(ROOTFILE_DIR_PREFIX,"Rootfiles");
//static const char* CUSTOM_REPLAY_OUTPUT_DIR=Form(ROOTFILE_DIR_PREFIX,"ScratchROOTfiles");
//static const char* SUMMARY_PHYSICS_FORMAT=Form(REPLAY_DIR_PREFIX,"summaryfiles/summaryphy_%d.log");

//not used since Nov 15, 2008
static const Int_t ANA_MARK_INTERVAL=1000;

// #include "ReplayCore.C"  

// DO NOT DELETE!
#endif//#define ROOT_REPLAY_DEF

