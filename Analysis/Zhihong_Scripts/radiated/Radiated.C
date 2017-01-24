/*Include{{{*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <errno.h>
#include <sstream>
//#include <iterator>
//from root
#include <TMath.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMatrix.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
/*}}}*/
using namespace::std;
using namespace::TMath;

#include "XEMC.h"

TString Kin, Tgt;
int getargs(int argc, char**argv);

int main(int argc,char** argv){
  // TString Tgt, Kin;
  // cerr<<"--- Target (H2,He3,He4,C12,Ca40,Ca48) = "; cin >> Tgt;
  // cerr<<"--- Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
  
  int err = getargs(argc, argv);

  if(err<0){
    cerr<<"****** I can not recognize the parameters!"<<endl;
    return -1;
  }

  double Ep, Theta, P0;
  if(Kin=="3.1"){
    Theta = 21.00; P0 = 2.905;
  }
  else if(Kin=="3.2"){
    Theta = 21.00; P0 = 3.055;
  }
  else if(Kin=="4.1"){
    Theta = 23.00; P0 = 2.855;
  }
  else if(Kin=="4.2"){
    Theta = 23.00; P0 = 3.035;
  }
  else if(Kin=="5.0"){
    Theta = 25.00; P0 = 2.505;
  }
  else if(Kin=="5.05"){
    Theta = 25.00; P0 = 2.650;
  }
  else if(Kin=="5.1"){
    Theta = 25.00; P0 = 2.795;
  }
  else if(Kin=="5.2"){
    Theta = 25.00; P0 = 2.995;
  }
  else if(Kin=="6.5"){
    Theta = 28.00;P0 = 2.845;
  }
  else
    cerr<<"*** ERROR, unknown kinematics setting!"<<endl;
  
  XEMC* test = new XEMC();
  test->Init(Form("./input/%s_Input.dat",Tgt.Data()));
   
  double xs_rad = -1.0, xs_born=-1.0, E0 =3.356, RctPt = 0.0;;
  
  TString fKin = "Kin"+Kin;
  ofstream outfile1(Form("./output/%s_%s_Rad.dat",Tgt.Data(),fKin.Data()));

  outfile1 << Form("%8s %8s %12s %12s","Ep","RctPt","xs_born","xs_rad")<<endl;

  const int Bin=20;

  // XS at 10 locations for Ep-3%
  RctPt = 0.0;
  Ep = P0 * (1 - 0.03);
  cerr<<"--- Working on Ep -3% ..."<<endl;
  for(int i=0; i<Bin;i++){
    RctPt = i*20./Bin - 9.5;

    test->Process(E0,Ep,Theta,RctPt);

    xs_rad = test->XS_Rad();
    xs_born = test->XS_Born();
    
   	outfile1 << Form("%8.4f %8.4f  %12.6e  %12.6e", Ep, RctPt, xs_born, xs_rad)<<endl;
    cerr <<Form("Ep=%f, Rct=%f, RadC = %f", Ep,RctPt,xs_rad/xs_born)
	 <<endl;
  }  
  // XS at 10 locations for Ep-0%
  RctPt = 0.0;
  Ep = P0 * (1 - 0.0);
  cerr<<"--- Working on Ep -0% ..."<<endl;
  for(int i=0; i<Bin;i++){
    RctPt = i*20./Bin - 9.5;

    test->Process(E0,Ep,Theta,RctPt);

    xs_rad = test->XS_Rad();
    xs_born = test->XS_Born();
     
   	outfile1 << Form("%8.4f %8.4f  %12.6e  %12.6e", Ep, RctPt, xs_born, xs_rad)<<endl;
    cerr <<Form("Ep=%f, Rct=%f, RadC = %f", Ep,RctPt,xs_rad/xs_born)
	 <<endl;
  }  
  // XS at 10 locations for Ep+2%
  RctPt = 0.0;
  Ep = P0 * (1 + 0.02);
  cerr<<"--- Working on Ep +2% ..."<<endl;
  for(int i=0; i<Bin;i++){
    RctPt = i*20./Bin - 9.5;

    test->Process(E0,Ep,Theta,RctPt);

    xs_rad = test->XS_Rad();
    xs_born = test->XS_Born();
    
   	outfile1 << Form("%8.4f %8.4f  %12.6e  %12.6e", Ep, RctPt, xs_born, xs_rad)<<endl;
    cerr <<Form("Ep=%f, Rct=%f, RadC = %f", Ep,RctPt,xs_rad/xs_born)
	 <<endl;
  }  


  delete test;
  outfile1.close();  
  return 0;
}


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
      case 'T':
	if(*(++argptr))
	  Tgt = argptr;
	err=0;
	break;
      case 'K':
	if(*(++argptr))
	  Kin = argptr;
	err=0;
	break;
      case 'h':
	cerr <<"================================================"<<endl;
	cerr <<"Option: " <<endl;
	cerr <<"-T Target Name (H2,He3,He4,C12,Ca40,Ca48)" <<endl;
	cerr <<"-K Kinematics Setting (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5)"<<endl;
	cerr <<endl<<"================================================"<<endl;
	noStop=false;
	err=-1;
	goto OUT;
	break;
      default:	
	cerr <<"================================================"<<endl;
	cerr <<"Option: " <<endl;
	cerr <<"-T Target Name (H2,He3,He4,C12,Ca40,Ca48)" <<endl;
	cerr <<"-K Kinematics Setting (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5)"<<endl;
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
