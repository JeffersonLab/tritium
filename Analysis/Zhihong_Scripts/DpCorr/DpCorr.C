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
  //int err = getargs(argc,argv);
 // if(err<0){
 //   cerr<<"****** I can not recognize the parameters!"<<endl;
 //   return -1;
//  }
 
  gFirst = 1;
  gLast =2;

  /*Definition{{{*/
  //ifstream infile(gInfile);
  ifstream infile("test.txt");
  TString ROOTFILE_PATH = "/work/halla/e08014/disk1/Rootfiles/";
  TString Log_Path ="./Log/";
  TString filename;
  //In T tree
  double R_tr_x, R_tr_y, R_tr_th, R_tr_ph, dpR_Old, EpR_Old,Nu_Old,Q2_Old, Xbj_Old, ExtdpR_Old; 
  double E0, SinSQ, Ep0, dpR, EpR, Q2R, NuR, XbjR, ExtdpR, ExtEpR, ExtQ2R, ExtNuR, ExtXbjR;
  /*}}}*/
  
  /*Dp Corr Const{{{*/
  //Addtional Correction for Delta on HRS-R, since the D-Term in optics was not calibrated.
  //The parameters were fitted from SAMC simulation data.
  //--- Z. Ye, 02/18/2013
  //===============================================================================
  //--X_fp:   
  const double X0=3.034150e-04,  X1=-6.629199e-03, X2=-8.535028e-04;
  //--Theta_fp: 
  const double T0=-1.060184e-04, T1=1.061409e-01,  T2=6.751897e-01;
  //--Phi_fp:   
  const double P0=6.492564e-05,  P1=2.171117e-05,  P2=-3.002779e-01;
  //--Y_fp:   
  const double Y0=2.933949e-05,  Y1=-5.402877e-04, Y2=-4.617217e-01, Y3=7.100115e-01, Y4=5.226948e+02;
  //--Delta:   
  const double D0=2.668049e-05,  D1=-4.939113e-03, D2=-5.875315e-02,  D3=2.054553e+00;
  //===============================================================================
  /*}}}*/
  
  int Counter =0;
  while(!(infile.eof())){
	  infile >> filename;
	  Counter++;
	  if(Counter>=gFirst&&Counter<=gLast){
		  cerr<<"Reading in Rootfile --"<<ROOTFILE_PATH+filename<<endl;
		  TFile *f1 = new TFile(ROOTFILE_PATH+filename,"update");
		  TTree  *T = (TTree*)f1->GetObjectChecked("T","Tree");

		  /*Set Branch Address{{{*/	
		  T->SetBranchAddress("R.tr.r_x", &R_tr_x);
		  T->SetBranchAddress("R.tr.r_y", &R_tr_y);
		  T->SetBranchAddress("R.tr.r_th", &R_tr_th);
		  T->SetBranchAddress("R.tr.r_ph", &R_tr_ph);
		  T->SetBranchAddress("R.gold.p", &EpR_Old);
		  T->SetBranchAddress("R.gold.dp", &dpR_Old);
		  T->SetBranchAddress("ExTgtR.dp", &ExtdpR_Old);
		  T->SetBranchAddress("EKR.nu",    &Nu_Old);
		  T->SetBranchAddress("EKR.Q2",    &Q2_Old);
		  T->SetBranchAddress("EKR.x_bj",  &Xbj_Old);

		  //Add New Branches
		  TBranch* DP = T->Branch("dpR",  &dpR,  "dpR/D");
		  TBranch* EP = T->Branch("EpR",  &EpR,  "EpR/D");
		  TBranch* Q2 = T->Branch("Q2R",  &Q2R,  "Q2R/D");
		  TBranch* NU = T->Branch("NuR",  &NuR,  "NuR/D");
		  TBranch* XBJ = T->Branch("XbjR",  &XbjR,  "XbjR/D");
		  TBranch* ExtDP = T->Branch("ExtdpR",  &ExtdpR,  "ExtdpR/D");
		  TBranch* ExtEP = T->Branch("ExtEpR",  &ExtEpR,  "ExtEpR/D");
		  TBranch* ExtQ2 = T->Branch("ExtQ2R",  &ExtQ2R,  "ExtQ2R/D");
		  TBranch* ExtNU = T->Branch("ExtNuR",  &ExtNuR,  "ExtNuR/D");
		  TBranch* ExtXBJ= T->Branch("ExtXbjR", &ExtXbjR, "ExtXbjR/D");
		  /*}}}*/

		  /////////////////////////////////////////////////////////////
		  //Now calculate and fill the current in Scaler LEFT tree
		  //////////////////////////////////////////////////////////////
		  /*Add_Branch{{{*/
		  int N = T->GetEntries();
		  for(int i=0;i<N;i++){
			  T->GetEntry(i);

			  if(EpR_Old>0.5&&EpR<5.0){
				  Ep0 = EpR_Old/(1.+dpR_Old);
				  E0 = Nu_Old + EpR_Old;
				  SinSQ = Q2_Old/(4.0*E0*EpR);

				  dpR = dpR_Old - (X0 + X1*pow(R_tr_x, 1) + X2*pow(R_tr_x, 2))
					  - (Y0 + Y1*pow(R_tr_y, 1) + Y2*pow(R_tr_y, 2) + Y3*pow(R_tr_y, 3) + Y4*pow(R_tr_y, 4))
					  - (T0 + T1*pow(R_tr_th,1) + T2*pow(R_tr_th,2))
					  - (P0 + P1*pow(R_tr_ph,1) + P2*pow(R_tr_ph,2))
					  - (D0 + D1*pow(dpR_Old,1) + D2*pow(dpR_Old,2) + D3*pow(dpR_Old,3));

				  ExtdpR = ExtdpR_Old - (X0 + X1*pow(R_tr_x, 1) + X2*pow(R_tr_x, 2))
					  - (Y0 + Y1*pow(R_tr_y, 1) + Y2*pow(R_tr_y, 2) + Y3*pow(R_tr_y, 3) + Y4*pow(R_tr_y, 4))
					  - (T0 + T1*pow(R_tr_th,1) + T2*pow(R_tr_th,2))
					  - (P0 + P1*pow(R_tr_ph,1) + P2*pow(R_tr_ph,2))
					  - (D0 + D1*pow(ExtdpR_Old,1) + D2*pow(ExtdpR_Old,2) + D3*pow(ExtdpR_Old,3));
				  EpR = Ep0*(1.+dpR); 			   
				  NuR = E0 - EpR;
				  Q2R = 4.0*E0*EpR*SinSQ;
				  XbjR = Q2R/(2.0*0.938272013*NuR);

				  ExtEpR = Ep0*(1.+ExtdpR); 			   
				  ExtNuR = E0 - ExtEpR;
				  ExtQ2R = 4.0*E0*ExtEpR*SinSQ;
				  ExtXbjR= ExtQ2R/(2.0*0.938272013*ExtNuR);
			  }
			  else{
				  dpR = -1.0; EpR = -1.0; NuR = -1.0; Q2R =-1.0; XbjR = -1.0;
				  ExtdpR = -1.0; ExtEpR = -1.0; ExtNuR = -1.0; ExtQ2R =-1.0; ExtXbjR = -1.0;
			  }

			  DP->Fill(); EP->Fill(); NU->Fill(); Q2->Fill(); XBJ->Fill();
			  ExtDP->Fill(); ExtEP->Fill(); ExtNU->Fill(); ExtQ2->Fill(); ExtXBJ->Fill();

			  if(!(i%50000))
				  cerr<<"--- Event# have been processed: "<<i<<endl;
		  }//  for(int i=1;i<N;i++){
		  /*}}}*/

		  //Save and close the rootfiles
		  T->Write("",TObject::kOverwrite);
		  delete DP;
		  delete EP;
		  delete Q2;
		  delete NU;
		  delete XBJ;
		  delete ExtDP;
		  delete ExtEP;
		  delete ExtQ2;
		  delete ExtNU;
		  delete ExtXBJ;
		  f1->Close();
		  cerr<<"=== C ==="<<endl;
	  }
	  }
	  infile.close();
  }//END

/*getargs{{{*/
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
/*}}}*/
