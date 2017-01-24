////////////////////////////////////////////////
// Generate Cross Section Look-up Table
//  --Zhihong Ye, 07/03/2012
//
// 10/02/2012, Change
////////////////////////////////////////////////
/*Include{{{*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <cmath>
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

const double DEGtoRAD = 3.1415926 / 180.0;
int gAngle;
int gMaxBin, gMinBin;	
char* gTarget;
int getargs(int argc, char**argv);

int main(int argc,char** argv){

	int err = getargs(argc,argv);
	if(err>=0){
		TString Target;
		Target = gTarget;

		Double_t E0 = 3.356; //GeV
		Double_t Ep_Min = 0.0, Ep_Max=0;   //GeV
		Double_t Theta = (double)(gAngle);
        Double_t Q2 = 0.0, SinSQ = 0.0, Xbj = 0.0;
		Int_t A=1; 
		Int_t Z=0;

		TString Target_Input = Form("./input/%s_Input.dat",Target.Data());	
		TString Table_File = Form("./%s/XS_%s_T%d_%d_%d.table",Target.Data(),Target.Data(),gAngle,gMinBin,gMaxBin);
        TString Log_File = Table_File;
		Log_File.ReplaceAll(".table",".log");
		ofstream outlog(Log_File.Data());
		/*Set Kin{{{*/
		const Double_t Ep_Bin_Size = 0.01; //1 MeV per bin
		if(gAngle==21)  {
			Ep_Min = 2.730; Ep_Max = 3.240; }//51 Bins 
		else if(gAngle==23) {   
			Ep_Min = 2.680; Ep_Max = 3.220; } //54 Bins
		else if(gAngle==25) {   
			Ep_Min = 2.350; Ep_Max = 3.180; } //83 Bins
		else if(gAngle==28) {   
			Ep_Min = 2.680; Ep_Max = 3.020; } //34 Bins
		else {
			outlog<< "I don't understand the Kinematics setting!"<<endl;}
		/*}}}*/

		/*Set Target{{{*/
		if(Target == "H2") {
			A = 2; Z = 1;}
		else if(Target == "He3") {
			A = 3; Z = 2;}
		else if(Target == "He4") {
			A = 4; Z = 2;}
		else if(Target == "C12") { 
			A = 12; Z = 6;}
		else if(Target == "Al") { //The Aluminium Window for the target cell, density use the average of up and down-stream.
			A = 27; Z = 13;}
		else if(Target == "Dummy") { 
			A = 27; Z = 13;}
		else if(Target == "Ca40") {
			A = 40; Z = 20;}
		else if(Target == "Ca48") {
			A = 48; Z = 20;}
		else{
			outlog<<"I don't understand the Target!"<<endl;}
		/*}}}*/

		Double_t Factor = 1.01;
		if(A<=4)
			Factor = 1.1;

		/*Gen_XS_Table{{{*/
		//Define Acceptace Range and Bin-Size
		Double_t Delta_Theta = 0.07*TMath::RadToDeg();      Int_t Bin_Theta = 20;   Double_t Size_Theta = Delta_Theta/Bin_Theta;

		//Define Output file name
		ofstream outfile(Table_File.Data());
        if(gMinBin==0)
			outfile << Form("%5s %3s %6s %10s %6s %10s %6s %16s %16s %16s", "A","Z","Scatt_Ang","Ep","EBin","Theta","TBin","xs_born","xs_rad","Tail")<<endl;

		//Define a event to calculate radiated cross section
		XEMC* Event = new XEMC(); 
		Event->Init(Target_Input.Data());
		Int_t err = -1;
		Double_t mom, theta, cs_rad,cs_born,cs_qe,cs_dis,tail;
		for(int j=gMinBin;  j<gMaxBin;   j++){
			mom = Ep_Min + j * Ep_Bin_Size;
			if(mom>Ep_Max){
				outlog<<Form("Out of Momentum Setting Range at Bin=%d !!!",j)<<endl;
				break;
			}
			for(int i=0;  i<Bin_Theta; i++){
				theta = Theta - Delta_Theta/2. + i * Size_Theta; 

				SinSQ = pow(sin(theta*DEGtoRAD/2.),2);
				Q2 = 4.0 * E0 * mom * SinSQ;
				Xbj = Q2 / (2.0*(E0-mom));

				if(Xbj<=Factor*A){//If Xbj > A, the XS basically is zero but I want to see a larger xbj due to the He3 problems.
					err = Event->Process(E0,mom,theta,A,Z,0.0);	
					if(err>=0){
						cs_rad = Event->XS_Rad();
						cs_born = Event->XS_Born();
						cs_qe = Event->XS_QE();
						cs_dis = Event->XS_DIS();
						tail = Event->Elastic_Tail();
						outfile << Form("%5d %3d %6.3f %10.7f %6d %10.7f %6d %16.6e %16.6e %16.6e", A,Z,Theta,mom,j,theta,i,cs_born,cs_rad,tail)<<endl;
						//if(!i%10)
						//	outlog << Form("%5d %3d %6.3f %10.7f %6d %10.7f %6d %16.6e %16.6e %16.6e", A,Z,Theta,mom,j,theta,i,cs_born,cs_rad,tail)<<endl;
					}
				}
				else{
					cs_rad = 0.0;
					cs_born = 0.0;
					cs_qe = 0.0;
					cs_dis = 0.0;
					tail = 0.0;
				}
			}
		} 
		delete Event;
		outfile.close();
		/*}}}*/

		outlog << Form("Cross Section Table for %s at Angle=%d, Ep Bin from %d to %d has been generated!",Target.Data(), gAngle,gMinBin, gMaxBin)<<endl;

	}
	return 0;
}

/*int getargs(int argc,char** argv){{{*/
int getargs(int argc,char** argv){
	char* argptr;
	bool noStop;  
	int err=-1;

	for(int i=0;i<argc;i++){
		argptr = argv[i];
		if(*argptr=='-'){
			argptr++;
			switch (*argptr){
				case 'A':
					if(*(++argptr))
						gAngle = atoi(argptr);
					err=0;
					break;
				case 'T':
					if(*(++argptr))
						gTarget = argptr;
					err=0;
					break;
				case 'L':
					if(*(++argptr))
						gMaxBin = atoi(argptr);
					err=0;
					break;
				case 'F':
					if(*(++argptr))
						gMinBin = atoi(argptr);
					err=0;
					break;
				case 'h':
					cerr <<"================================================"<<endl;
					cerr <<"Option: " <<endl;
					cerr <<"-A[21,23,25,28]  Scattering Angle" <<endl;
					cerr <<"-T[H2,He3,He4,C12,Ca40,Ca48,...]  Target Name" <<endl;
					cerr <<"-F[1] First Bin" <<endl;
					cerr <<"-L[100] Last Bin" <<endl;
					cerr <<endl<<"================================================"<<endl;
					noStop=false;
					err=-1;
					goto OUT;
					break;
				default:
					cerr <<"A code to generate cross section tables for X>2 Experiment: " << argptr << endl;
					cerr <<endl<<"================================================"<<endl;
					cerr <<"Option: " <<endl;
					cerr <<"-A[21,23,25,28]  Scattering Angle" <<endl;
					cerr <<"-T[H2,He3,He4,C12,Ca40,Ca48,...]  Target Name" <<endl;
					cerr <<"-F[1] First Bin" <<endl;
					cerr <<"-L[100] Last Bin" <<endl;
					cerr <<"             Zhihong Ye 10/02/2012" <<endl;
					cerr <<"================================================"<<endl;
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
/**/
