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

#include "BCM.h"

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

  ifstream infile(gInfile);
  TString ROOTFILE_PATH = "/work/halla/e08014/disk1/Rootfiles/";
  TString Log_Path ="./Log/";
  /*Definition{{{*/
  //In Scaler tree
  TString filename;
  double ls_bcm_d1c,ls_bcm_d3c,ls_bcm_u1c,ls_bcm_u3c;
  double rs_bcm_d1c,rs_bcm_d3c,rs_bcm_u1c,rs_bcm_u3c;
  double ls_clkcount, rs_clkcount;
  //In T tree
  double left_bcm_d1c,left_bcm_d3c,left_bcm_u1c,left_bcm_u3c;
  double right_bcm_d1c,right_bcm_d3c,right_bcm_u1c,right_bcm_u3c;
  double left_clkcount,right_clkcount;

  double left_current, right_current;
  double clk_prev=0.0, u1c_prev=0.0, d1c_prev=0.0,u3c_prev=0.0, d3c_prev=0.0;
  double time = 0.0,right_time=0.0,left_time = 0.0;
  //Current
  double I_u1=0.0, I_d1=0.0, I_u3=0.0, I_d3=0.0;
  double IR_u1=0.0, IR_d1=0.0, IR_u3=0.0, IR_d3=0.0;
  double IL_u1=0.0, IL_d1=0.0, IL_u3=0.0, IL_d3=0.0;
  //Charge 
  double CL_u1=0.0, CL_d1=0.0, CL_u3=0.0, CL_d3=0.0;
  double CR_u1=0.0, CR_d1=0.0, CR_u3=0.0, CR_d3=0.0;
  double C_u1=0.0, C_d1=0.0, C_u3=0.0, C_d3=0.0;

  /*}}}*/

  /*BCM calibration constants defined in BCM.h {{{*/
  //LFET
  double offset_Lu1 = u_Offset[0][0];
  double offset_Ld1 = d_Offset[0][0];
  double alpha_Lu1 = u_Coeff[0][0];
  double alpha_Ld1 = d_Coeff[0][0];
  double offset_Lu3 = u_Offset[0][1];
  double offset_Ld3 = d_Offset[0][1];
  double alpha_Lu3 = u_Coeff[0][1];
  double alpha_Ld3 = d_Coeff[0][1];
  //RIGHT
  double offset_Ru1 = u_Offset[1][0];
  double offset_Rd1 = d_Offset[1][0];
  double alpha_Ru1 = u_Coeff[1][0];
  double alpha_Rd1 = d_Coeff[1][0];
  double offset_Ru3 = u_Offset[1][1];
  double offset_Rd3 = d_Offset[1][1];
  double alpha_Ru3 = u_Coeff[1][1];
  double alpha_Rd3 = d_Coeff[1][1];
  /*}}}*/
  int Counter =0;
  while(!(infile.eof())){
    infile >> filename;
    Counter++;
    if(Counter>=gFirst&&Counter<=gLast)
      { 
	cerr<<"Reading in Rootfile --"<<ROOTFILE_PATH+filename<<endl;
	TString LogFile = filename;
	LogFile.ReplaceAll(".root",".log");
	LogFile = Log_Path+LogFile;
	ofstream outfile(LogFile);
	cerr<<"Saving the log in "<<LogFile<<endl;	
	TFile *f1 = new TFile(ROOTFILE_PATH+filename,"update");
	TTree  *T = (TTree*)f1->GetObjectChecked("T","Tree");
	TTree *LS = (TTree*)f1->GetObjectChecked("LEFT","Tree");
	TTree *RS = (TTree*)f1->GetObjectChecked("RIGHT","Tree");

	/*Set Branch Address{{{*/	
	LS->SetBranchAddress("bcm_d1c",&ls_bcm_d1c);
	LS->SetBranchAddress("bcm_d3c",&ls_bcm_d3c);
	LS->SetBranchAddress("bcm_u1c",&ls_bcm_u1c);
	LS->SetBranchAddress("bcm_u3c",&ls_bcm_u3c);
	LS->SetBranchAddress("clkcount",&ls_clkcount);

	RS->SetBranchAddress("bcm_d1c",&rs_bcm_d1c);
	RS->SetBranchAddress("bcm_d3c",&rs_bcm_d3c);
	RS->SetBranchAddress("bcm_u1c",&rs_bcm_u1c);
	RS->SetBranchAddress("bcm_u3c",&rs_bcm_u3c);
	RS->SetBranchAddress("clkcount",&rs_clkcount);

	T->SetBranchAddress("left_clkcount", &left_clkcount);
	T->SetBranchAddress("left_bcm_d1c",  &left_bcm_d1c);
	T->SetBranchAddress("left_bcm_d3c",  &left_bcm_d3c);
	T->SetBranchAddress("left_bcm_u1c",  &left_bcm_u1c);
	T->SetBranchAddress("left_bcm_u3c",  &left_bcm_u3c);

	T->SetBranchAddress("right_clkcount", &right_clkcount);
	T->SetBranchAddress("right_bcm_d1c",  &right_bcm_d1c);
	T->SetBranchAddress("right_bcm_d3c",  &right_bcm_d3c);
	T->SetBranchAddress("right_bcm_u1c",  &right_bcm_u1c);
	T->SetBranchAddress("right_bcm_u3c",  &right_bcm_u3c);

	//LEFT Scaler Current
	TBranch* I_LS_D1 = LS->Branch("I_d1",  &I_d1,  "I_d1/D");
	TBranch* I_LS_D3 = LS->Branch("I_d3",  &I_d3,  "I_d3/D");
	TBranch* I_LS_U1 = LS->Branch("I_u1",  &I_u1,  "I_u1/D");
	TBranch* I_LS_U3 = LS->Branch("I_u3",  &I_u3,  "I_u3/D");
	TBranch* I_LS_Current = LS->Branch("left_current",  &left_current,  "left_current/D");
	//RIGHT Scaler Current
	TBranch* I_RS_D1 = RS->Branch("I_d1",  &I_d1,  "I_d1/D");
	TBranch* I_RS_D3 = RS->Branch("I_d3",  &I_d3,  "I_d3/D");
	TBranch* I_RS_U1 = RS->Branch("I_u1",  &I_u1,  "I_u1/D");
	TBranch* I_RS_U3 = RS->Branch("I_u3",  &I_u3,  "I_u3/D");
	TBranch* I_RS_Current = RS->Branch("right_current", &right_current, "right_current/D");
	//LEFT Scaler current in T
	TBranch* I_LT_D1 = T->Branch("IL_d1",  &IL_d1,  "IL_d1/D");
	TBranch* I_LT_D3 = T->Branch("IL_d3",  &IL_d3,  "IL_d3/D");
	TBranch* I_LT_U1 = T->Branch("IL_u1",  &IL_u1,  "IL_u1/D");
	TBranch* I_LT_U3 = T->Branch("IL_u3",  &IL_u3,  "IL_u3/D");
	TBranch* I_LT_Current = T->Branch("left_current",   &left_current,  "left_current/D");
	//RIGHT Scaler current in T
	TBranch* I_RT_D1 = T->Branch("IR_d1",  &IR_d1,  "IR_d1/D");
	TBranch* I_RT_D3 = T->Branch("IR_d3",  &IR_d3,  "IR_d3/D");
	TBranch* I_RT_U1 = T->Branch("IR_u1",  &IR_u1,  "IR_u1/D");
	TBranch* I_RT_U3 = T->Branch("IR_u3",  &IR_u3,  "IR_u3/D");
	TBranch* I_RT_Current = T->Branch("right_current",  &right_current, "right_current/D");

	//LEFT Scaler Charge
	TBranch* C_LS_D1 = LS->Branch("C_d1",  &C_d1,  "C_d1/D");
	TBranch* C_LS_D3 = LS->Branch("C_d3",  &C_d3,  "C_d3/D");
	TBranch* C_LS_U1 = LS->Branch("C_u1",  &C_u1,  "C_u1/D");
	TBranch* C_LS_U3 = LS->Branch("C_u3",  &C_u3,  "C_u3/D");
	TBranch* C_LS_Time = LS->Branch("time",  &time,  "time/D");
	//RIGHT Scaler Charge
	TBranch* C_RS_D1 = RS->Branch("C_d1",  &C_d1,  "C_d1/D");
	TBranch* C_RS_D3 = RS->Branch("C_d3",  &C_d3,  "C_d3/D");
	TBranch* C_RS_U1 = RS->Branch("C_u1",  &C_u1,  "C_u1/D");
	TBranch* C_RS_U3 = RS->Branch("C_u3",  &C_u3,  "C_u3/D");
	TBranch* C_RS_Time = RS->Branch("time", &time, "time/D");
	//LEFT Scaler Charge in T
	TBranch* C_LT_D1 = T->Branch("CL_d1",  &CL_d1,  "CL_d1/D");
	TBranch* C_LT_D3 = T->Branch("CL_d3",  &CL_d3,  "CL_d3/D");
	TBranch* C_LT_U1 = T->Branch("CL_u1",  &CL_u1,  "CL_u1/D");
	TBranch* C_LT_U3 = T->Branch("CL_u3",  &CL_u3,  "CL_u3/D");
	TBranch* C_LT_Time = T->Branch("left_time",   &left_time,  "left_time/D");
	//RIGHT Scaler Chrarge in T
	TBranch* C_RT_D1 = T->Branch("CR_d1",  &CR_d1,  "CR_d1/D");
	TBranch* C_RT_D3 = T->Branch("CR_d3",  &CR_d3,  "CR_d3/D");
	TBranch* C_RT_U1 = T->Branch("CR_u1",  &CR_u1,  "CR_u1/D");
	TBranch* C_RT_U3 = T->Branch("CR_u3",  &CR_u3,  "CR_u3/D");
	TBranch* C_RT_Time = T->Branch("right_time",  &right_time, "right_time/D");
	/*}}}*/

	/////////////////////////////////////////////////////////////
	//Now calculate and fill the current in Scaler LEFT tree
	//////////////////////////////////////////////////////////////
	/*LEFT_Tree{{{*/
	int NL = LS->GetEntries();
	//Save the values of the first events
	LS->GetEntry(0);
	clk_prev = ls_clkcount;
	d1c_prev = ls_bcm_d1c;	u1c_prev = ls_bcm_u1c;
	d3c_prev = ls_bcm_d3c;	u3c_prev = ls_bcm_u3c;
	//Always fill the current of first event as zero
	left_current=0.0;
	I_LS_Current->Fill(); 

	//Then skip the first reading
	for(int iL=1;iL<NL;iL++){
	  LS->GetEntry(iL);
	  time = (ls_clkcount - clk_prev)/1024.0;
	  I_u1 = (((ls_bcm_u1c-u1c_prev)/time)-offset_Lu1)/alpha_Lu1;
	  I_d1 = (((ls_bcm_d1c-d1c_prev)/time)-offset_Ld1)/alpha_Ld1; 
	  I_u3 = (((ls_bcm_u3c-u3c_prev)/time)-offset_Lu3)/alpha_Lu3;
	  I_d3 = (((ls_bcm_d3c-d3c_prev)/time)-offset_Ld3)/alpha_Ld3; 
	  //Take the average of four values.
	  //left_current = (I_u1+I_u3+I_d1+I_d3)/4.0; 
	  left_current = I_d1; //Left scaler has problem, d1 might be reliable 
	  if(left_current<0||left_current>200.0)
	    left_current=-1000.0;

	  C_u1 = I_u1*time; C_u3 = I_u3*time;	C_d1 = I_d1*time; C_d3 = I_d3*time;
	  I_LS_D1->Fill(); I_LS_D3->Fill(); I_LS_U1->Fill(); I_LS_U3->Fill();	I_LS_Current->Fill();
	  C_LS_D1->Fill(); C_LS_D3->Fill(); C_LS_U1->Fill(); C_LS_U3->Fill();	C_LS_Time->Fill();

	  outfile<<Form(" >>> LEFT Tree: %d events have been proceeded!",iL)<<endl;
	  outfile<< Form("     clk_prev=%f, ls_clk=%f, left_current=%f", clk_prev,ls_clkcount,left_current)<<endl;

	  //Store the values for the next event
	  clk_prev = ls_clkcount;
	  d1c_prev = ls_bcm_d1c;	u1c_prev = ls_bcm_u1c;
	  d3c_prev = ls_bcm_d3c;	u3c_prev = ls_bcm_u3c;
	}//  for(int iL=1;iL<NL;iL++){
	/*}}}*/

	/////////////////////////////////////////////////////////////
	//Now calculate and fill the current in Scaler Right tree
	//////////////////////////////////////////////////////////////
	/*RIGHT_Tree{{{*/
	int NR = RS->GetEntries();
	//Save the values of the first events
	RS->GetEntry(0);
	clk_prev = rs_clkcount;
	d1c_prev = rs_bcm_d1c;	u1c_prev = rs_bcm_u1c;
	d3c_prev = rs_bcm_d3c;	u3c_prev = rs_bcm_u3c;
	//Always fill the current of first event as zero
	right_current=0.0;
	I_RS_Current->Fill(); 

	//Then skip the first reading
	for(int iR=1;iR<NR;iR++){
	  RS->GetEntry(iR);
	  time = (rs_clkcount - clk_prev)/1024.0;
	  I_u1 = (((rs_bcm_u1c-u1c_prev)/time)-offset_Ru1)/alpha_Ru1;
	  I_d1 = (((rs_bcm_d1c-d1c_prev)/time)-offset_Rd1)/alpha_Rd1; 
	  I_u3 = (((rs_bcm_u3c-u3c_prev)/time)-offset_Ru3)/alpha_Ru3;
	  I_d3 = (((rs_bcm_d3c-d3c_prev)/time)-offset_Rd3)/alpha_Rd3; 
	  //Take the average of four values.
	  right_current = (I_u1+I_u3+I_d1+I_d3)/4.0; 

	  if(right_current<0||right_current>200.0)
	    right_current=-1000.0;

	  C_u1 = I_u1*time; C_u3 = I_u3*time;	C_d1 = I_d1*time; C_d3 = I_d3*time;
	  I_RS_D1->Fill(); I_RS_D3->Fill(); I_RS_U1->Fill(); I_RS_U3->Fill();	I_RS_Current->Fill();
	  C_RS_D1->Fill(); C_RS_D3->Fill(); C_RS_U1->Fill(); C_RS_U3->Fill();	C_RS_Time->Fill();

	  outfile<<Form(" >>> RIGHT Tree: %d events have been proceeded!",iR)<<endl;
	  outfile<<Form("     clk_prev=%f, rs_clk=%f, right_current=%f", clk_prev,rs_clkcount,right_current)<<endl;

	  //Store the values for the next event
	  clk_prev = rs_clkcount;	
	  d1c_prev = rs_bcm_d1c;	u1c_prev = rs_bcm_u1c;
	  d3c_prev = rs_bcm_d3c;	u3c_prev = rs_bcm_u3c;
	}//  for(int iR=1;iR<NR;iR++){
	/*}}}*/

	/////////////////////////////////////////////////////////////
	//Now calculate and fill the current in Scaler T tree
	//////////////////////////////////////////////////////////////
	/*T_Tree{{{*/
	double clk_prevR=0.0, u1c_prevR=0.0, d1c_prevR=0.0,u3c_prevR=0.0, d3c_prevR=0.0;
	double clk_prevL=0.0, u1c_prevL=0.0, d1c_prevL=0.0,u3c_prevL=0.0, d3c_prevL=0.0;
	double timeL=0.0, timeR=0.0;
	int N_prevL =0, N_prevR=0;
	double current=0.0;
	int NT = T->GetEntries();
	left_current=0.0; right_current=0.0;

	for(int iT=1;iT<NT;iT++){
	  T->GetEntry(iT);
	  //      right_clkcount+=4*1024;

	  if(left_clkcount>clk_prevL){
	    timeL = (left_clkcount - clk_prevL)/1024.0;
	    IL_u1 = (((left_bcm_u1c-u1c_prevL)/timeL)-offset_Ru1)/alpha_Ru1;
	    IL_d1 = (((left_bcm_d1c-d1c_prevL)/timeL)-offset_Rd1)/alpha_Rd1; 
	    IL_u3 = (((left_bcm_u3c-u3c_prevL)/timeL)-offset_Ru3)/alpha_Ru3;
	    IL_d3 = (((left_bcm_d3c-d3c_prevL)/timeL)-offset_Rd3)/alpha_Rd3; 
	    //Take the average of four values.
	    //current = (IL_u1+IL_u3+IL_d1+IL_d3)/4.0; 
	    current = IL_d1; //Left Sclaer has problems, d1 may be working better. 

	    if(current<0||current>200.0)
	      current=-1000.0;

	    outfile<< Form(">>> T for  Left:  clk_prev=%12.2f, left_clk=%12.2f, left_current=%12.2f, right_current=%12.2f",
			   clk_prevL,left_clkcount,current,right_current)<<endl;

	    //Now loop back to the previous scaler reading time and fill all current values
	    for(int i=N_prevL;i<iT;i++){
	      T->GetEntry(i);
	      left_current =current;
	      left_time = timeL;
	      CL_u1 = IL_u1*timeL; CL_u3 = IL_u3*timeL;	CL_d1 = IL_d1*timeL; CL_d3 = IL_d3*timeL;
	      I_LT_D1->Fill(); I_LT_D3->Fill(); I_LT_U1->Fill(); I_LT_U3->Fill();	I_LT_Current->Fill();
	      C_LT_D1->Fill(); C_LT_D3->Fill(); C_LT_U1->Fill(); C_LT_U3->Fill();	C_LT_Time->Fill();
	    }
	    //Store values for the next scaler reading
	    u1c_prevL = left_bcm_u1c;	u3c_prevL = left_bcm_u3c;
	    d1c_prevL = left_bcm_d1c;	d3c_prevL = left_bcm_d3c;
	    clk_prevL = left_clkcount;	N_prevL = iT;
	  }

	  if(right_clkcount>clk_prevR){
	    timeR = (right_clkcount - clk_prevR)/1024.0;
	    IR_u1 = (((right_bcm_u1c-u1c_prevR)/timeR)-offset_Ru1)/alpha_Ru1;
	    IR_d1 = (((right_bcm_d1c-d1c_prevR)/timeR)-offset_Rd1)/alpha_Rd1; 
	    IR_u3 = (((right_bcm_u3c-u3c_prevR)/timeR)-offset_Ru3)/alpha_Ru3;
	    IR_d3 = (((right_bcm_d3c-d3c_prevR)/timeR)-offset_Rd3)/alpha_Rd3; 
	    //Take the average of four values.
	    //current = (IR_u1+IR_u3+IR_d1+IR_d3)/4.0; 
	    current = (IR_u1+IR_d1)/2.0; 
	    if(current<0||current>200.0)
	      current=-1000.0;
	    outfile<< Form(">>> T for Right: clk_prev=%12.2f, right_clk=%12.2f, left_current=%12.2f, right_current=%12.2f", 
			   clk_prevR,right_clkcount,left_current,current)<<endl;

	    //Now loop back to the previous scaler reading time and fill all current values
	    for(int j=N_prevR;j<iT;j++){
	      T->GetEntry(j);
	      right_current = current;
	      if(!(j%10000))
		cerr<<Form("---> Left Current = %6.2f, Right Current = %6.2f", left_current, right_current)<<endl;
	      right_time = timeR;
	      CR_u1 = IR_u1*time; CR_u3 = IR_u3*time;	CR_d1 = IR_d1*time; CR_d3 = IR_d3*time;
	      I_RT_D1->Fill(); I_RT_D3->Fill(); I_RT_U1->Fill(); I_RT_U3->Fill();	I_RT_Current->Fill();
	      C_RT_D1->Fill(); C_RT_D3->Fill(); C_RT_U1->Fill(); C_RT_U3->Fill();	C_RT_Time->Fill();
	    }
	    //Store values for the next scaler reading
	    u1c_prevR = right_bcm_u1c;	u3c_prevR = right_bcm_u3c;
	    d1c_prevR = right_bcm_d1c;	d3c_prevR = right_bcm_d3c;
	    clk_prevR = right_clkcount;	N_prevR = iT;
	  }
	}//  for(int iT=1;iT<NT;iT++){
	/*}}}*/

	//Save and close the rootfiles
	LS->Write("",TObject::kOverwrite); RS->Write("",TObject::kOverwrite); T->Write("",TObject::kOverwrite);
	f1->Close();
	outfile.close();
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
