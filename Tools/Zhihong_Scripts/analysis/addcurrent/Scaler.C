//////////////////////////////////////////////////////////////
// Extract trigger events in between two consecutive scaler
// reading, so we can cut out events we don't want and still
// get the correct total values of reading in one run
//
// Zhihong Ye, 08/23/2012


/*Include{{{*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <errno.h>
#include <sstream>
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
/*}}}*/

using namespace::std;
using namespace::TMath;

char* gInfile;
char* gOutfile;
int gLast=0;
int gFirst=0;
int getargs(int argc, char**argv);

int main(int argc,char** argv) 
{
  int err = getargs(argc,argv);
  if(err<0){
    cerr<<"****** I can not recognize the parameters!"<<endl;
    return -1;
  }

  /*Definition{{{*/
  //In Scaler tree
  double ls_t1c, ls_t2c,ls_t3c,ls_t4c,ls_t5c,ls_t6c,ls_t7c,ls_t8c;
  double rs_t1c, rs_t2c,rs_t3c,rs_t4c,rs_t5c,rs_t6c,rs_t7c,rs_t8c;
  double ls_clkcount, rs_clkcount;
  //In T tree
  double left_t1c,left_t2c,left_t3c,left_t4c,left_t5c,left_t6c,left_t7c,left_t8c;
  double right_t1c,right_t2c,right_t3c,right_t4c,right_t5c,right_t6c,right_t7c,right_t8c;
  double left_clkcount,right_clkcount, right_LT, left_LT,right_LT_err, left_LT_err;
  //Variables for store values of previous events
  double prev_t1c, prev_t2c, prev_t3c, prev_t4c, prev_t5c, prev_t6c, prev_t7c, prev_t8c;

  //New variables to add, real time counts
  double rt_t1c, rt_t2c, rt_t3c, rt_t4c, rt_t5c, rt_t6c, rt_t7c, rt_t8c;
  double left_rt_t1c, left_rt_t2c, left_rt_t3c, left_rt_t4c, left_rt_t5c, left_rt_t6c, left_rt_t7c, left_rt_t8c;
  double right_rt_t1c,right_rt_t2c,right_rt_t3c,right_rt_t4c,right_rt_t5c,right_rt_t6c,right_rt_t7c,right_rt_t8c;

  TString infilename; infilename=gInfile;
  // cerr<<" Which Runlist ?(Run_List.txt) "; cin >> infilename;
  ifstream infile(infilename);

  TString ROOTFILE_PATH = "/work/halla/e08014/disk1/Rootfiles/";
  TString Log_Path="./Log/";
  int Counter=0;
  TString filename;
  while(!(infile.eof())){
    infile >> filename;
    if(filename=="END")
      break;
    Counter++;
    if(Counter>=gFirst&&Counter<=gLast)
      {
	cerr<<"Reading in Rootfile --"<<ROOTFILE_PATH+filename<<endl;
	TString LogFile = filename;

	TFile *f1 = new TFile(ROOTFILE_PATH+filename,"update");
	// //Obtain PS
	// THaRun *aRun = (THaRun*) f1->Get("Run_Data"); 
	// TArrayI aps = aRun->GetParameters()->GetPrescales();
	// Int_t PS[8];
	// for(int i=0;i<8;i++){
	//   PS[i] = aps[i];
	//   cerr << Form(" PS%d = %d ",i+1,PS[i]);
	// }
	// cerr<<endl;

	TTree  *T = (TTree*)f1->GetObjectChecked("T","Tree");
	TTree *LS = (TTree*)f1->GetObjectChecked("LEFT","Tree");
	TTree *RS = (TTree*)f1->GetObjectChecked("RIGHT","Tree");

	/*Set Branch Address{{{*/	
	LS->SetBranchAddress("t1c",&ls_t1c);
	LS->SetBranchAddress("t2c",&ls_t2c);
	LS->SetBranchAddress("t3c",&ls_t3c);
	LS->SetBranchAddress("t4c",&ls_t4c);
	//    LS->SetBranchAddress("t5c",&ls_t5c);
	LS->SetBranchAddress("t6c",&ls_t6c);
	LS->SetBranchAddress("t7c",&ls_t7c);	
	LS->SetBranchAddress("t8c",&ls_t8c);
	LS->SetBranchAddress("clkcount",&ls_clkcount);

	RS->SetBranchAddress("t1c",&rs_t1c);
	RS->SetBranchAddress("t2c",&rs_t2c);
	RS->SetBranchAddress("t3c",&rs_t3c);
	RS->SetBranchAddress("t4c",&rs_t4c);
	//    RS->SetBranchAddress("t5c",&rs_t5c);
	RS->SetBranchAddress("t6c",&rs_t6c);
	RS->SetBranchAddress("t7c",&rs_t7c);	
	RS->SetBranchAddress("t8c",&rs_t8c);
	RS->SetBranchAddress("clkcount",&rs_clkcount);

	T->SetBranchAddress("left_t1c",&left_t1c);
	T->SetBranchAddress("left_t2c",&left_t2c);
	T->SetBranchAddress("left_t3c",&left_t3c);
	T->SetBranchAddress("left_t4c",&left_t4c);
	//    T->SetBranchAddress("left_t5c",&left_t5c);
	T->SetBranchAddress("left_t6c",&left_t6c);
	T->SetBranchAddress("left_t7c",&left_t7c);
	T->SetBranchAddress("left_t8c",&left_t8c);
	T->SetBranchAddress("left_clkcount",  &left_clkcount);

	T->SetBranchAddress("right_t1c",&right_t1c);
	T->SetBranchAddress("right_t2c",&right_t2c);
	T->SetBranchAddress("right_t3c",&right_t3c);
	T->SetBranchAddress("right_t4c",&right_t4c);
	//    T->SetBranchAddress("right_t5c",&right_t5c);
	T->SetBranchAddress("right_t6c",&right_t6c);
	T->SetBranchAddress("right_t7c",&right_t7c);
	T->SetBranchAddress("right_t8c",&right_t8c);
	T->SetBranchAddress("right_clkcount",  &right_clkcount);

	//LEFT Scaler Current
	TBranch* LS_t1c = LS->Branch("rt_t1c",  &rt_t1c,  "rt_t1c/D");
	TBranch* LS_t2c = LS->Branch("rt_t2c",  &rt_t2c,  "rt_t2c/D");
	TBranch* LS_t3c = LS->Branch("rt_t3c",  &rt_t3c,  "rt_t3c/D");
	TBranch* LS_t4c = LS->Branch("rt_t4c",  &rt_t4c,  "rt_t4c/D");
	//    TBranch* LS_t5c = LS->Branch("rt_t5c",  &rt_t5c,  "rt_t5c/D");
	TBranch* LS_t6c = LS->Branch("rt_t6c",  &rt_t6c,  "rt_t6c/D");
	TBranch* LS_t7c = LS->Branch("rt_t7c",  &rt_t7c,  "rt_t7c/D");
	TBranch* LS_t8c = LS->Branch("rt_t8c",  &rt_t8c,  "rt_t8c/D");

	//RIGHT Scaler Current
	TBranch* RS_t1c = RS->Branch("rt_t1c",  &rt_t1c,  "rt_t1c/D");
	TBranch* RS_t2c = RS->Branch("rt_t2c",  &rt_t2c,  "rt_t2c/D");
	TBranch* RS_t3c = RS->Branch("rt_t3c",  &rt_t3c,  "rt_t3c/D");
	TBranch* RS_t4c = RS->Branch("rt_t4c",  &rt_t4c,  "rt_t4c/D");
	//    TBranch* RS_t5c = RS->Branch("rt_t5c",  &rt_t5c,  "rt_t5c/D");
	TBranch* RS_t6c = RS->Branch("rt_t6c",  &rt_t6c,  "rt_t6c/D");
	TBranch* RS_t7c = RS->Branch("rt_t7c",  &rt_t7c,  "rt_t7c/D");
	TBranch* RS_t8c = RS->Branch("rt_t8c",  &rt_t8c,  "rt_t8c/D");

	//LEFT Scaler current in T
	TBranch* TL_t1c = T->Branch("left_rt_t1c",  &left_rt_t1c,  "left_rt_t1c/D");
	TBranch* TL_t2c = T->Branch("left_rt_t2c",  &left_rt_t2c,  "left_rt_t2c/D");
	TBranch* TL_t3c = T->Branch("left_rt_t3c",  &left_rt_t3c,  "left_rt_t3c/D");
	TBranch* TL_t4c = T->Branch("left_rt_t4c",  &left_rt_t4c,  "left_rt_t4c/D");
	//    TBranch* TL_t5c = T->Branch("left_rt_t5c",  &left_rt_t5c,  "left_rt_t5c/D");
	TBranch* TL_t6c = T->Branch("left_rt_t6c",  &left_rt_t6c,  "left_rt_t6c/D");
	TBranch* TL_t7c = T->Branch("left_rt_t7c",  &left_rt_t7c,  "left_rt_t7c/D");
	TBranch* TL_t8c = T->Branch("left_rt_t8c",  &left_rt_t8c,  "left_rt_t8c/D");
	// TBranch* TL_LT = T->Branch("left_LT",  &left_LT,  "left_LT/D");
	// TBranch* TL_LT_Err = T->Branch("left_LT_err",  &left_LT_err,  "left_LT_err/D");

	//RIGHT Scaler current in T
	TBranch* TR_t1c = T->Branch("right_rt_t1c",  &right_rt_t1c,  "right_rt_t1c/D");
	TBranch* TR_t2c = T->Branch("right_rt_t2c",  &right_rt_t2c,  "right_rt_t2c/D");
	TBranch* TR_t3c = T->Branch("right_rt_t3c",  &right_rt_t3c,  "right_rt_t3c/D");
	TBranch* TR_t4c = T->Branch("right_rt_t4c",  &right_rt_t4c,  "right_rt_t4c/D");
	//    TBranch* TR_t5c = T->Branch("right_rt_t5c",  &right_rt_t5c,  "right_rt_t5c/D");
	TBranch* TR_t6c = T->Branch("right_rt_t6c",  &right_rt_t6c,  "right_rt_t6c/D");
	TBranch* TR_t7c = T->Branch("right_rt_t7c",  &right_rt_t7c,  "right_rt_t7c/D");
	TBranch* TR_t8c = T->Branch("right_rt_t8c",  &right_rt_t8c,  "right_rt_t8c/D");
	// TBranch* TR_LT = T->Branch("right_LT",  &right_LT,  "right_LT/D");
	// TBranch* TR_LT_Err = T->Branch("right_LT_err",  &right_LT_err,  "right_LT_err/D");
	/*}}}*/

	/////////////////////////////////////////////////////////////
	//Now calculate and fill the events in Scaler LEFT tree
	//////////////////////////////////////////////////////////////
	/*LEFT_Tree{{{*/
	int NL = LS->GetEntries();
	//Save the values of the first events
	LS->GetEntry(0);
	prev_t1c = ls_t1c; prev_t2c = ls_t2c; prev_t3c = ls_t3c; prev_t4c = ls_t4c;
	//prev_t5c = ls_t5c; 
	prev_t6c = ls_t6c; prev_t7c = ls_t7c; prev_t8c = ls_t8c;
	rt_t1c = 0.0;  rt_t2c = 0.0;  rt_t3c = 0.0;  rt_t4c = 0.0;
	//rt_t5c = 0.0;  
	rt_t6c = 0.0;  rt_t7c = 0.0;  rt_t8c = 0.0;
	LS_t1c->Fill(); LS_t2c->Fill(); LS_t3c->Fill(); LS_t4c->Fill(); 
	//LS_t5c->Fill(); 
	LS_t6c->Fill(); LS_t7c->Fill(); LS_t8c->Fill(); 

	//Then skip the first reading
	for(int iL=1;iL<NL;iL++){
	  LS->GetEntry(iL);

	  rt_t1c=ls_t1c-prev_t1c; rt_t2c=ls_t2c-prev_t2c; rt_t3c=ls_t3c-prev_t3c; rt_t4c=ls_t4c-prev_t4c;
	  //rt_t5c=ls_t5c-prev_t5c; 
	  rt_t6c=ls_t6c-prev_t6c; rt_t7c=ls_t7c-prev_t7c; rt_t8c=ls_t8c-prev_t8c;
	  LS_t1c->Fill(); LS_t2c->Fill(); LS_t3c->Fill(); LS_t4c->Fill(); 
	  //LS_t5c->Fill(); 
	  LS_t6c->Fill(); LS_t7c->Fill(); LS_t8c->Fill(); 

	  prev_t1c = ls_t1c; prev_t2c = ls_t2c; prev_t3c = ls_t3c; prev_t4c = ls_t4c;
	  //  prev_t5c = ls_t5c;
	  prev_t6c = ls_t6c; prev_t7c = ls_t7c; prev_t8c = ls_t8c;
	}//  for(int iL=1;iL<NL;iL++){
	LS->Write("",TObject::kOverwrite); 
	/*}}}*/

	/////////////////////////////////////////////////////////////
	//Now calculate and fill the events in Scaler Right tree
	//////////////////////////////////////////////////////////////
	/*RIGHT_Tree{{{*/
	int NR = RS->GetEntries();
	//Save the values of the first events
	RS->GetEntry(0);
	prev_t1c = rs_t1c; prev_t2c = rs_t2c; prev_t3c = rs_t3c; prev_t4c = rs_t4c;
	//prev_t5c = rs_t5c;
	prev_t6c = rs_t6c; prev_t7c = rs_t7c; prev_t8c = rs_t8c;
	rt_t1c = 0.0;  rt_t2c = 0.0;  rt_t3c = 0.0;  rt_t4c = 0.0;
	//rt_t5c = 0.0;  
	rt_t6c = 0.0;  rt_t7c = 0.0;  rt_t8c = 0.0;
	RS_t1c->Fill(); RS_t2c->Fill(); RS_t3c->Fill(); RS_t4c->Fill(); 
	//RS_t5c->Fill();
	RS_t6c->Fill(); RS_t7c->Fill(); RS_t8c->Fill(); 

	//Then skip the first reading
	for(int iR=1;iR<NR;iR++){
	  RS->GetEntry(iR);

	  rt_t1c=rs_t1c-prev_t1c; rt_t2c=rs_t2c-prev_t2c; rt_t3c=rs_t3c-prev_t3c; rt_t4c=rs_t4c-prev_t4c;
	  // rt_t5c=rs_t5c-prev_t5c; 
	  rt_t6c=rs_t6c-prev_t6c; rt_t7c=rs_t7c-prev_t7c; rt_t8c=rs_t8c-prev_t8c;
	  RS_t1c->Fill(); RS_t2c->Fill(); RS_t3c->Fill(); RS_t4c->Fill(); 
	  //RS_t5c->Fill();
	  RS_t6c->Fill(); RS_t7c->Fill(); RS_t8c->Fill(); 

	  prev_t1c = rs_t1c; prev_t2c = rs_t2c; prev_t3c = rs_t3c; prev_t4c = rs_t4c;
	  // prev_t5c = rs_t5c; 
	  prev_t6c = rs_t6c; prev_t7c = rs_t7c; prev_t8c = rs_t8c;
	}//  for(int iR=1;iR<NR;iR++){
	RS->Write("",TObject::kOverwrite); 
	/*}}}*/

	/////////////////////////////////////////////////////////////
	//Now calculate and fill the events in Scaler T tree
	//////////////////////////////////////////////////////////////
	/*T_Tree{{{*/
	double temp_t1c, temp_t2c, temp_t3c, temp_t4c, temp_t5c, temp_t6c, temp_t7c, temp_t8c;
	double temp_LT,temp_LT_err;
	double prevL_t1c,prevL_t2c,prevL_t3c,prevL_t4c,prevL_t5c,prevL_t6c,prevL_t7c,prevL_t8c;
	double prevR_t1c,prevR_t2c,prevR_t3c,prevR_t4c,prevR_t5c,prevR_t6c,prevR_t7c,prevR_t8c;
	double clk_prevL, clk_prevR;
	int N_prevL =0, N_prevR=0;
	int NT = T->GetEntries();

	TString Left_DAQ_Cut, Right_DAQ_Cut;
	Int_t Left_Trigger, Right_Trigger;

	//Save the first event:
	T->GetEntry(0);
	prevL_t1c = left_t1c; prevL_t2c = left_t2c; prevL_t3c = left_t3c; prevL_t4c = left_t4c;
	//prevL_t5c = left_t5c;
	prevL_t6c = left_t6c; prevL_t7c = left_t7c; prevL_t8c = left_t8c;
	left_rt_t1c = 0.0;  left_rt_t2c = 0.0;  left_rt_t3c = 0.0;  left_rt_t4c = 0.0;
	//left_rt_t5c = 0.0;  
	left_rt_t6c = 0.0;  left_rt_t7c = 0.0;  left_rt_t8c = 0.0;
	TL_t1c->Fill(); TL_t2c->Fill(); TL_t3c->Fill(); TL_t4c->Fill(); 
	//TL_t5c->Fill(); 
	TL_t6c->Fill(); TL_t7c->Fill(); TL_t8c->Fill(); 

	prevR_t1c = right_t1c; prevR_t2c = right_t2c; prevR_t3c = right_t3c; prevR_t4c = right_t4c;
	//prevR_t5c = right_t5c;
	prevR_t6c = right_t6c; prevR_t7c = right_t7c; prevR_t8c = right_t8c;
	right_rt_t1c = 0.0;  right_rt_t2c = 0.0;  right_rt_t3c = 0.0;  right_rt_t4c = 0.0;
	//right_rt_t5c = 0.0; 
	right_rt_t6c = 0.0;  right_rt_t7c = 0.0;  right_rt_t8c = 0.0;
	TR_t1c->Fill(); TR_t2c->Fill(); TR_t3c->Fill(); TR_t4c->Fill(); 
	//TR_t5c->Fill();
	TR_t6c->Fill(); TR_t7c->Fill(); TR_t8c->Fill(); 

	clk_prevL = left_clkcount; clk_prevR = right_clkcount;

	for(int iT=1;iT<NT;iT++){
	  T->GetEntry(iT);

	  if(left_clkcount>clk_prevL){
	    temp_t1c = left_t1c-prevL_t1c; 	temp_t2c = left_t2c-prevL_t2c; 	temp_t3c = left_t3c-prevL_t3c; 	temp_t4c = left_t4c-prevL_t4c; 	
	    //	temp_t5c = left_t5c-prevL_t5c;
	    temp_t6c = left_t6c-prevL_t6c; 	temp_t7c = left_t7c-prevL_t7c; 	temp_t8c = left_t8c-prevL_t8c; 
	    // Left_DAQ_Cut = Form("DBB.evtypebits>>%d&1 && left_clkcount>=%f && left_clkcount<%f", 3, clk_prevL,left_clkcount);
	    // Left_Trigger = T->GetEntries(Left_DAQ_Cut.Data());
	    // temp_LT = Left_Trigger*PS[2]/temp_t3c;
	    // temp_LT_err = temp_LT * (sqrt(1.0/Left_Trigger+1.0/temp_t3c));

	    // if(!(iT/1000))
	    //   cerr<<Form("Left_DAQ = %d, Left_t3c = %f, LT = %f", Left_Trigger, temp_t3c, temp_LT)<<endl;

	    //Now loop back to the previous scaler reading time and fill all current values
	    for(int i=N_prevL;i<iT;i++){
	      T->GetEntry(i);
	      left_rt_t1c = temp_t1c;left_rt_t2c = temp_t2c;left_rt_t3c = temp_t3c;left_rt_t4c = temp_t4c;
	      // left_rt_t5c = temp_t5c;
	      left_rt_t6c = temp_t6c;left_rt_t7c = temp_t7c;left_rt_t8c = temp_t8c;
	      //	  left_LT = temp_LT; left_LT_err = temp_LT_err;

	      TL_t1c->Fill(); TL_t2c->Fill(); TL_t3c->Fill(); TL_t4c->Fill(); 
	      // TL_t5c->Fill(); 
	      TL_t6c->Fill(); TL_t7c->Fill(); TL_t8c->Fill(); 
	      //  TL_LT->Fill(); TL_LT_Err->Fill();
	    }

	    //Store values for the next scaler reading
	    prevL_t1c = left_t1c; prevL_t2c = left_t2c; prevL_t3c = left_t3c; prevL_t4c = left_t4c;
	    //prevL_t5c = left_t5c; 
	    prevL_t6c = left_t6c; prevL_t7c = left_t7c; prevL_t8c = left_t8c;
	    clk_prevL = left_clkcount;
	    N_prevL = iT;
	  }

	  if(right_clkcount>clk_prevR){
	    temp_t1c = right_t1c-prevR_t1c;	temp_t2c = right_t2c-prevR_t2c;	temp_t3c = right_t3c-prevR_t3c;	temp_t4c = right_t4c-prevR_t4c; 	
	    //temp_t5c = right_t5c-prevR_t5c;	
	    temp_t6c = right_t6c-prevR_t6c;	temp_t7c = right_t7c-prevR_t7c;	temp_t8c = right_t8c-prevR_t8c; 

	    // Right_DAQ_Cut = Form("DBB.evtypebits>>%d&1 && right_clkcount>=%f && right_clkcount<%f", 1, clk_prevR,right_clkcount);
	    // Right_Trigger = T->GetEntries(Right_DAQ_Cut.Data());
	    // temp_LT = Right_Trigger*PS[0]/temp_t1c;
	    // temp_LT_err = temp_LT * (sqrt(1.0/Right_Trigger+1.0/temp_t1c));

	    //	cerr<<Form("Right_DAQ = %d, Right_t1c = %f, LT = %f", Right_Trigger, temp_t1c,temp_LT)<<endl;

	    //Now loop back to the previous scaler reading time and fill all current values
	    for(int i=N_prevR;i<iT;i++){
	      T->GetEntry(i);
	      right_rt_t1c = temp_t1c;right_rt_t2c = temp_t2c;right_rt_t3c = temp_t3c;right_rt_t4c = temp_t4c;
	      //right_rt_t5c = temp_t5c;
	      right_rt_t6c = temp_t6c;right_rt_t7c = temp_t7c;right_rt_t8c = temp_t8c;
	      //  right_LT = temp_LT; right_LT_err = temp_LT_err;

	      TR_t1c->Fill(); TR_t2c->Fill(); TR_t3c->Fill(); TR_t4c->Fill(); 
	      //TR_t5c->Fill();
	      TR_t6c->Fill(); TR_t7c->Fill(); TR_t8c->Fill(); 
	      //	  TR_LT->Fill();  TR_LT_Err->Fill();
	    }
	    //Store values for the next scaler reading
	    prevR_t1c = right_t1c; prevR_t2c = right_t2c; prevR_t3c = right_t3c; prevR_t4c = right_t4c;
	    //prevR_t5c = right_t5c; 
	    prevR_t6c = right_t6c; prevR_t7c = right_t7c; prevR_t8c = right_t8c;
	    clk_prevR = right_clkcount;
	    N_prevR = iT;
	  }
	}//  for(int iT=1;iT<NT;iT++){
	T->Write("",TObject::kOverwrite);
	/*}}}*/

	//Save and close the rootfiles
	f1->Close();
      }
  }// while(!(infile.eof()))

  infile.close();
}//END


int getargs(int argc,char** argv)
{
  char* argptr;
  TString aPar;
  bool noStop;  
  int err=-1;

  for(int i=0;i<argc;i++){
    argptr = argv[i];
    if(*argptr=='-'){
      argptr++;
      switch (*argptr){
      case 'I':
	if(*(++argptr))
	  gInfile = argptr;
	err=0;
	break;
      case 'O':
	if(*(++argptr))
	  gOutfile = argptr;
	err=0;
	break;
      case 'F':
	if(*(++argptr)){
	  aPar=argptr;
	  gFirst = aPar.Atoi();
	}
	err=0;
	break;
      case 'L':
	if(*(++argptr)){
	  aPar=argptr;
	  gLast = aPar.Atoi();
	}
	err=0;
	break;
      case 'h':
	cerr <<"================================================"<<endl;
	cerr <<"Option: " <<endl;
	cerr <<"-I[Rootfile_List.txt] List of rootfiles" <<endl;
	cerr <<"-O Outputfile" <<endl;
	cerr <<"-F First Line" <<endl;
	cerr <<"-L Lat Line"<<endl;
	cerr <<endl<<"================================================"<<endl;
	noStop=false;
	err=-1;
	goto OUT;
	break;
      default:
	cerr <<"================================================"<<endl;
	cerr <<"Option: " <<endl;
	cerr <<"-I[Rootfile_List.txt] List of rootfiles" <<endl;
	cerr <<"-O Outputfile" <<endl;
	cerr <<"-F First Line" <<endl;
	cerr <<"-L Lat Line"<<endl;
	cerr <<endl<<"================================================"<<endl;
	err=-1;
	break;
      }

    }
    noStop=true;
  OUT:if(!noStop){break;}
    continue;
  }
  return err; 
}
		
