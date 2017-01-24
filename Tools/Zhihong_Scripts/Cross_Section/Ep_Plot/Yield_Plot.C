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
//#include <iterator>
//from root
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMatrix.h>
#include <TFile.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TTree.h>
#include <TStopwatch.h>
/*}}}*/
using namespace::std;
using namespace::TMath;
#include "./XGT2_Main.h"
#include "/work/halla/e08014/disk1/yez/Ep_Cross_Section/SRC/XGT2_XEMC.h"

char* gFile1;
char* gFile2;
char* gFile3;
char* gFile4;
char* gTgt;
char* gCom;
int getargs(int argc, char** argv);

/*Get_Target_Name{{{*/
TString gGet_Target_Name(const TString& aName){
	TString aTarget_Name="MM";

	if(aName.Contains("H2"))
		aTarget_Name="H2";
	else if(aName.Contains("He3"))
		aTarget_Name="He3";
	else if(aName.Contains("He4"))
		aTarget_Name="He4";
	else if(aName.Contains("C12"))
		aTarget_Name="C12";
	else if(aName.Contains("Ca40"))
		aTarget_Name="Ca40";
	else if(aName.Contains("Ca48"))
		aTarget_Name="Ca48";
	else
		cerr<<"*** I don't know the target name!!!"<<endl;
	return aTarget_Name;
}
/*}}}*/

/*Get_Kin_Name{{{*/
TString gGet_Kin_Name(const TString& aKin){
	TString aKin_Name="MM";
	if(aKin.Contains("3.1"))
		aKin_Name="3.1";
	else if(aKin.Contains("3.2"))
		aKin_Name="3.2";
	else if(aKin.Contains("4.1"))
		aKin_Name="4.1";
	else if(aKin.Contains("4.2"))
		aKin_Name="4.2";
	else if(aKin.Contains("5.0") && !aKin.Contains("5.05") )
		aKin_Name="5.0";
	else if(aKin.Contains("5.05"))
		aKin_Name="5.05";
	else if(aKin.Contains("5.1"))
		aKin_Name="5.1";
	else if(aKin.Contains("5.2"))
		aKin_Name="5.2";
	else if(aKin.Contains("6.5"))
		aKin_Name="6.5";
	else
		cerr<<"*** I don't know the kinematics name!!!"<<endl;
	return aKin_Name;
}
/*}}}*/
TString aFileName1, aFileName2, aFileName3,aFileName4;
TString Target, Comment;

/*Yield Plot{{{*/
void Yield_Plot()
{  
	/*ROOT Define{{{*/
	// Fill color
	gStyle->SetStatColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetPadColor(0);
	gStyle->SetFrameFillColor(0);

	// Border mode
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetPalette(1,0);

	// Font
	gStyle->SetTextFont(32);
	gStyle->SetLabelFont(32, "xyz");
	gStyle->SetTitleFont(22, "xyz");
	gStyle->SetStatFont(32);
	gStyle->SetTitleFillColor(0);

	// Opt
	//    gStyle->SetOptDate(0);
	//    gStyle->SetOptTitle(1);
	gStyle->SetOptStat(0);
	//    gStyle->SetOptFit(1);
	//    // Margin
	// //    gStyle->SetPadLeftMargin(0.1);
	// //    gStyle->SetPadRightMargin(0.1);
	// //    gStyle->SetPadTopMargin(0.1);
	// //    gStyle->SetPadBottomMargin(0.1);


	//void HallA_style() {
	gROOT->SetStyle("Plain");
	gStyle->SetPaperSize(TStyle::kUSLetter);
	gStyle->SetPaperSize(18,22);
	gStyle->SetOptFit(1111);
	gStyle->SetPalette(1);
	gStyle->SetNdivisions(505);

	gStyle->SetCanvasColor(0);
	gStyle->SetPadTopMargin(.05);
	gStyle->SetPadLeftMargin(.15);
	gStyle->SetPadRightMargin(.1);
	gStyle->SetPadBottomMargin(.15);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetLabelFont(32,"X");
	gStyle->SetLabelFont(32,"Y");

	// prepare gStyle to be useful
	//   1 = solid
	//   2 = long dash (30 10)
gStyle->SetTitleFillColor(-1); 
	//   3 = dotted (4 8)
	//   4 = dash-dot (15 12 4 12)
	//   5 = short dash ( 15 15 )
	//   6 = dash-dot-dot 
	gStyle->SetLineStyleString(1,"[]");
	gStyle->SetLineStyleString(2,"[30 10]");
	gStyle->SetLineStyleString(3,"[4 8]");
	gStyle->SetLineStyleString(4,"[15 12 4 12]");
	gStyle->SetLineStyleString(5,"[15 15]");
	gStyle->SetLineStyleString(6,"[15 12 4 12 4 12]");
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetNdivisions(505,"Y");
	gStyle->SetOptDate(0);
	gStyle->SetDateY(.98);
	gStyle->SetStripDecimals(kFALSE);

	gStyle->SetHistFillColor(0); 
	gStyle->SetFillColor(0); 
	gStyle->SetFillStyle(4000); 
	gStyle->SetStatStyle(0); 
	gStyle->SetTitleStyle(0); 
	gStyle->SetCanvasBorderSize(0); 
	gStyle->SetFrameBorderSize(0); 
	gStyle->SetLegendBorderSize(0); 
	gStyle->SetStatBorderSize(0); 
	gStyle->SetTitleBorderSize(0); 

	/*}}}*/
	
	/*Define{{{*/
	double E0 = 3.356;
	//	TString Target = "C12";
	//	cerr<<"Target (H2,He3,He4,C12,Ca40,Ca48) = "; cin >> Target;
	double Theta = 21.00;
	TString Kin = "3.1";

	Double_t p,nu,Q2,xbj,angle,Ei,Ep,tmp,Livetime,Ntg,Ni,Nf_EX,Nf_AL,Nf_SAMC,Nf_SAMC_Gen,BinCor,XS_MC;
	Double_t XS,XS_Stat_Err,XS_Sys_Err,Yield_MC,Yield_EX,Yield_AL,Y_MC_Sys_Err,Y_EX_Sys_Err,Y_AL_Sys_Err,Y_MC_Stat_Err,Y_EX_Stat_Err,Y_AL_Stat_Err;
	TString comment;
	Double_t Min=1e50, Max = -1.0;
	Double_t Min_X=1e50, Max_X = -1.0;
	Double_t Min_Yield=1e50, Max_Yield=-1.0;
	Double_t Min_R=1e50, Max_R=-1.0;
	Double_t pXS1[200],pXS_Err1[200],pXbj1[200],pXbj_Err1[200],pXS_MC1[200],pYield_MC1[200],pYield_EX1[200],pYield_AL1[200],pYield_AL_Err1[200],pYield_MC_Err1[200],pYield_EX_Err1[200];
	Double_t pXS2[200],pXS_Err2[200],pXbj2[200],pXbj_Err2[200],pXS_MC2[200],pYield_MC2[200],pYield_EX2[200],pYield_AL2[200],pYield_AL_Err2[200],pYield_MC_Err2[200],pYield_EX_Err2[200];
	Double_t pXS3[200],pXS_Err3[200],pXbj3[200],pXbj_Err3[200],pXS_MC3[200],pYield_MC3[200],pYield_EX3[200],pYield_AL3[200],pYield_AL_Err3[200],pYield_MC_Err3[200],pYield_EX_Err3[200];
	Double_t pXS4[200],pXS_Err4[200],pXbj4[200],pXbj_Err4[200],pXS_MC4[200],pYield_MC4[200],pYield_EX4[200],pYield_AL4[200],pYield_AL_Err4[200],pYield_MC_Err4[200],pYield_EX_Err4[200];

	Double_t pY_R1[200],pY_R2[200],pY_R_Err1[200],pY_R_Err2[200];
	Double_t pY_R3[200],pY_R4[200],pY_R_Err3[200],pY_R_Err4[200];
	/*}}}*/

	/*First{{{*/
	TString infile1 = aFileName1;
	cerr<<"==== Reading 1st file from "<<infile1<<endl;
	ifstream run_1; run_1.open(infile1);

	Int_t Count1 =0;
	// Target = gGet_Target_Name(infile1.Data());
	Kin = gGet_Kin_Name(infile1.Data());

	TString Arm1="X";
	if(infile1.Contains("L_"))
		Arm1="L";
	else if(infile1.Contains("R_"))
		Arm1="R";


	if(Kin=="MM"){
		cerr<<"1st: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
	}
	if(Kin=="3.1"||Kin=="3.2")
		Theta = 21.00; 
	else if(Kin=="4.1"||Kin=="4.2")
		Theta = 23.00;
	else if(Kin=="5.05"||Kin=="5.0"||Kin=="5.1"||Kin=="5.2")
		Theta =  25.00;
	else if(Kin=="6.5")
		Theta = 28.00;
	else
		cerr<<"*** ERROR, unknown kinematics setting!"<<endl;

	TString Kin_Name1 = Form("Kin%s",Kin.Data());

	Get_XS* xs1 = new Get_XS();
	xs1->Load_Table(Target.Data(),Kin_Name1.Data());

	comment.ReadLine(run_1);
	while ( run_1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
	{
		//	if(XS_Stat_Err<0.0){
		//		Y_EX_Stat_Err = Yield_EX * sqrt(1./Nf_EX+1./Ni);
		//		XS_Stat_Err = XS * sqrt( pow(Y_EX_Stat_Err/Yield_EX,2)+pow(Y_MC_Stat_Err/Yield_MC,2));
		//	}
		//		XS*=1e33; XS_Stat_Err*=1e33;  XS_Sys_Err *= 1e33; //From cm^2/MeV/sr to nb/MeV/sr
		//		Yield_EX*=1e33; Y_EX_Stat_Err*=1e33; Y_EX_Sys_Err*=1e33;
		if ( angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
		{
			pXS1[Count1] = XS;      pXS_Err1[Count1] = XS_Stat_Err;
			pXbj1[Count1] = xbj;    pXbj_Err1[Count1] = 0.0;

			pYield_MC1[Count1]=Yield_MC; 
			pYield_EX1[Count1]=Yield_EX;
			pYield_AL1[Count1]=Yield_AL;
			pYield_MC_Err1[Count1]=Y_MC_Stat_Err; 
			pYield_EX_Err1[Count1]=Y_EX_Stat_Err;
			pYield_AL_Err1[Count1]=Y_AL_Stat_Err;
			pY_R1[Count1]=Yield_EX/Yield_MC;
			pY_R_Err1[Count1]=pY_R1[Count1] * sqrt( pow(Y_EX_Stat_Err/Yield_EX,2) + pow(Y_MC_Stat_Err/Yield_MC,2) );

			pXS_MC1[Count1] = xs1->gGet_XS_Born(E0,Ep,Theta);
			//pXS_MC1[Count1] = xs1->gGet_XS_Born_ByKin(E0,Ep/1000.0,Theta);

			//cerr << Form("   --Xbj = %f, XS_EX = %e, XS_MC=%e ", xbj, XS, pXS_MC1[Count1] ) <<endl;
			if(Min_Yield>=Yield_MC&&Yield_MC>=1e3)
				Min_Yield = Yield_MC;
			if(Max_Yield<=Yield_MC)
				Max_Yield = Yield_MC;
			if(Min_Yield>=Yield_EX&&Yield_EX>=1e3)
				Min_Yield = Yield_EX;
			if(Max_Yield<=Yield_EX)
				Max_Yield = Yield_EX;

			if(Min_X>=xbj)
				Min_X = xbj;	
			if(Max_X<=xbj)
				Max_X = xbj;

			if(Min>=XS)
				Min = XS;	
			if(Max<=XS)
				Max = XS;
			if(Min>= pXS_MC1[Count1])
				Min =  pXS_MC1[Count1];	
			if(Max<= pXS_MC1[Count1])
				Max =  pXS_MC1[Count1];

			if(Min_R>=pY_R1[Count1])
				Min_R = pY_R1[Count1];
			if(Max_R<=pY_R1[Count1])
				Max_R = pY_R1[Count1];

			Count1++;
		}
		if ( !run_1.good() )
			break;
	}
	run_1.close();
	delete xs1;
	cerr << "  Total XS Points for 1st = " << Count1 <<endl;
	/*}}}*/

	/*Second{{{*/
	TString Kin_Name2 ="";
	Bool_t bMore = kTRUE;
	TString infile2 = aFileName2;
	cerr<<"==== Reading 2nd file from "<<infile2<<endl;
	TString Arm2="X";
	if(infile2 =="x")
		bMore = kFALSE;
	if(infile1.Contains("L_"))
		Arm2="L";
	else if(infile1.Contains("R_"))
		Arm2="R";

	Int_t Count2=0;
	if(bMore){
		ifstream run_2; run_2.open(infile2);
		//	Target = gGet_Target_Name(infile2.Data());
		Kin = gGet_Kin_Name(infile2.Data());
		if(Kin=="MM"){
			cerr<<"2nd: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
		}
		if(Kin=="3.1"||Kin=="3.2")
			Theta = 21.00;
		else if(Kin=="4.1"||Kin=="4.2")
			Theta = 23.00;
		else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
			Theta =  25.00;
		else if(Kin=="6.5")
			Theta = 28.00;
		else
			cerr<<"*** ERROR, unknown kinematics setting!"<<endl;
		Kin_Name2 = Form("Kin%s",Kin.Data());
		Get_XS* xs2 = new Get_XS();
		xs2->Load_Table(Target.Data(),Kin_Name2.Data());
		//xs2->Load_Table_ByKin(Target.Data(),Kin_Name2.Data());

		comment.ReadLine(run_2);
		while ( run_2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{

			//		if(XS_Stat_Err<0.0){
			//		  Y_EX_Stat_Err = Yield_EX * sqrt(1./Nf_EX+1./Ni);
			//		  XS_Stat_Err = XS * sqrt( pow(Y_EX_Stat_Err/Yield_EX,2)+pow(Y_MC_Stat_Err/Yield_MC,2));
			//		}
			//			XS*=1e33; XS_Stat_Err*=1e33;  XS_Sys_Err *= 1e33; //From cm^2/MeV/sr to nb/MeV/sr
			//            Yield_EX*=1e33; Y_EX_Stat_Err*=1e33; Y_EX_Sys_Err*=1e33;
			if ( angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				pXS2[Count2] = XS;      pXS_Err2[Count2] = XS_Stat_Err;
				pXbj2[Count2] = xbj;    pXbj_Err2[Count2] = 0.0;
				pYield_MC2[Count2]=Yield_MC; 
				pYield_EX2[Count2]=Yield_EX;
				pYield_AL2[Count2]=Yield_AL;
				pYield_MC_Err2[Count2]=Y_MC_Stat_Err; 
				pYield_EX_Err2[Count2]=Y_EX_Stat_Err;
				pYield_AL_Err2[Count2]=Y_AL_Stat_Err;
				pY_R2[Count2]=Yield_EX/Yield_MC;
				pY_R_Err2[Count2]=pY_R2[Count2] * sqrt( pow(Y_EX_Stat_Err/Yield_EX,2) + pow(Y_MC_Stat_Err/Yield_MC,2) );

				pXS_MC2[Count2] = xs2->gGet_XS_Born(E0,Ep,Theta);
				//pXS_MC2[Count2] = xs2->gGet_XS_Born_ByKin(E0,Ep/1000.0,Theta);
			//	cerr << Form("   --Xbj = %f, XS_EX = %e, XS_MC=%e ", xbj, XS, pXS_MC2[Count2] ) <<endl;

				if(Min_Yield>=Yield_MC&&Yield_MC>1e3)
					Min_Yield = Yield_MC;
				if(Max_Yield<=Yield_MC)
					Max_Yield = Yield_MC;
				if(Min_Yield>=Yield_EX&&Yield_EX>1e3)
					Min_Yield = Yield_EX;
				if(Max_Yield<=Yield_EX)
					Max_Yield = Yield_EX;
				if(Min>=XS)
					Min = XS;	
				if(Max<=XS)
					Max = XS;
				if(Min_X>=xbj)
					Min_X = xbj;	
				if(Max_X<=xbj)
					Max_X = xbj;

				if(Min>= pXS_MC2[Count2])
					Min =  pXS_MC2[Count2];	
				if(Max<= pXS_MC2[Count2])
					Max =  pXS_MC2[Count2];

				if(Min_R>=pY_R2[Count2])
					Min_R = pY_R2[Count2];
				if(Max_R<=pY_R2[Count2])
					Max_R = pY_R2[Count2];

				Count2++;
			}
			if ( !run_2.good() )
				break;
		}
		run_2.close();
		delete xs2;
		cerr << "  Total XS Points for 2nd  = " << Count2 <<endl;
	}
	/*}}}*/

	/*Third{{{*/
	Bool_t bMore1 = kTRUE;
	TString infile3 = aFileName3;
//		if(infile3 =="x")
//			bMore1 = kFALSE;
	TString Arm3="X";
	if(infile3.Contains("L_"))
		Arm3="L";
	else if(infile3.Contains("R_"))
		Arm3="R";

	Int_t Count3=0;
	TString Kin_Name3; 


	if(bMore1){
	Kin = gGet_Kin_Name(infile3.Data());

	if(Kin=="MM"){
		cerr<<"3rd: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
	}
	if(Kin=="3.1"||Kin=="3.2")
		Theta = 21.00; 
	else if(Kin=="4.1"||Kin=="4.2")
		Theta = 23.00;
	else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
		Theta =  25.00;
	else if(Kin=="6.5")
		Theta = 28.00;
	else
		cerr<<"*** ERROR, unknown kinematics setting!"<<endl;

	 Kin_Name3 = Form("Kin%s",Kin.Data());

		ifstream run_3; run_3.open(infile3);

		// Target = gGet_Target_Name(infile3.Data());

		cerr<<"==== Reading 3rd file from "<<infile3<<endl;

		Get_XS* xs3 = new Get_XS();
		xs3->Load_Table(Target.Data(),Kin_Name3.Data());
		//xs3->Load_Table_ByKin(Target.Data(),Kin_Name3.Data());

		comment.ReadLine(run_3);
		while ( run_3>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				pXS3[Count3] = XS;      pXS_Err3[Count3] = XS_Stat_Err;
				pXbj3[Count3] = xbj;    pXbj_Err3[Count3] = 0.0;

				pYield_MC3[Count3]=Yield_MC; 
				pYield_EX3[Count3]=Yield_EX;
				pYield_AL3[Count3]=Yield_AL;
				pYield_MC_Err3[Count3]=Y_MC_Stat_Err; 
				pYield_EX_Err3[Count3]=Y_EX_Stat_Err;
				pYield_AL_Err3[Count3]=Y_AL_Stat_Err;
				pY_R3[Count3]=Yield_EX/Yield_MC;
				pY_R_Err3[Count3]=pY_R3[Count3] * sqrt( pow(Y_EX_Stat_Err/Yield_EX,2) + pow(Y_MC_Stat_Err/Yield_MC,2) );

				pXS_MC3[Count3] = xs3->gGet_XS_Born(E0,Ep,Theta);
				//pXS_MC3[Count3] = xs3->gGet_XS_Born_ByKin(E0,Ep/1000.0,Theta);

			//	cerr << Form("   --Xbj = %f, XS_EX = %e, XS_MC=%e ", xbj, XS, pXS_MC3[Count3] ) <<endl;
				if(Min_Yield>=Yield_MC&&Yield_MC>=1e3)
					Min_Yield = Yield_MC;
				if(Max_Yield<=Yield_MC)
					Max_Yield = Yield_MC;
				if(Min_Yield>=Yield_EX&&Yield_EX>=1e3)
					Min_Yield = Yield_EX;
				if(Max_Yield<=Yield_EX)
					Max_Yield = Yield_EX;

				if(Min_X>=xbj)
					Min_X = xbj;	
				if(Max_X<=xbj)
					Max_X = xbj;

				if(Min>=XS)
					Min = XS;	
				if(Max<=XS)
					Max = XS;
				if(Min>= pXS_MC3[Count3])
					Min =  pXS_MC1[Count3];	
				if(Max<= pXS_MC3[Count3])
					Max =  pXS_MC3[Count3];

				if(Min_R>=pY_R3[Count3])
					Min_R = pY_R3[Count3];
				if(Max_R<=pY_R3[Count3])
					Max_R = pY_R3[Count3];

				Count3++;
			}
			if ( !run_3.good() )
				break;
		}
		run_3.close();
		delete xs3;
		cerr << "  Total XS Points for 3rd = " << Count3 <<endl;
	}
	/*}}}*/

	/*Forth{{{*/
	TString Kin_Name4 ="";
	Bool_t bMore2= kTRUE;

	TString infile4 = aFileName4;
	cerr<<"==== Reading 4th file from "<<infile4<<endl;
	TString Arm4="X";
	if(infile4 =="x")
		bMore2= kFALSE;
	if(infile4.Contains("L_"))
		Arm4="L";
	else if(infile4.Contains("R_"))
		Arm4="R";

	Int_t Count4=0;
	if(bMore2){
		ifstream run_4; run_4.open(infile4);
		//	Target = gGet_Target_Name(infile4.Data());
		Kin = gGet_Kin_Name(infile4.Data());
		if(Kin=="MM"){
			cerr<<"4th: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
		}
		if(Kin=="3.1"||Kin=="3.2")
			Theta = 21.00;
		else if(Kin=="4.1"||Kin=="4.2")
			Theta = 23.00;
		else if(Kin=="5.0"||Kin=="5.05"||Kin=="5.1"||Kin=="5.2")
			Theta =  25.00;
		else if(Kin=="6.5")
			Theta = 28.00;
		else
			cerr<<"*** ERROR, unknown kinematics setting!"<<endl;
		Kin_Name4 = Form("Kin%s",Kin.Data());
		Get_XS* xs4 = new Get_XS();
		xs4->Load_Table(Target.Data(),Kin_Name4.Data());
		//xs4->Load_Table_ByKin(Target.Data(),Kin_Name4.Data());
		comment.ReadLine(run_4);
		while ( run_4>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{

			//		if(XS_Stat_Err<0.0){
			//		  Y_EX_Stat_Err = Yield_EX * sqrt(1./Nf_EX+1./Ni);
			//		  XS_Stat_Err = XS * sqrt( pow(Y_EX_Stat_Err/Yield_EX,2)+pow(Y_MC_Stat_Err/Yield_MC,2));
			//		}
			//			XS*=1e33; XS_Stat_Err*=1e33;  XS_Sys_Err *= 1e33; //From cm^2/MeV/sr to nb/MeV/sr
			//            Yield_EX*=1e33; Y_EX_Stat_Err*=1e33; Y_EX_Sys_Err*=1e33;
			if ( angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				pXS4[Count4] = XS;      pXS_Err4[Count4] = XS_Stat_Err;
				pXbj4[Count4] = xbj;    pXbj_Err4[Count4] = 0.0;
				pYield_MC4[Count4]=Yield_MC; 
				pYield_EX4[Count4]=Yield_EX;
				pYield_AL4[Count4]=Yield_AL;
				pYield_MC_Err4[Count4]=Y_MC_Stat_Err; 
				pYield_EX_Err4[Count4]=Y_EX_Stat_Err;
				pYield_AL_Err4[Count4]=Y_AL_Stat_Err;
				pY_R4[Count4]=Yield_EX/Yield_MC;
				pY_R_Err4[Count4]=pY_R4[Count4] * sqrt( pow(Y_EX_Stat_Err/Yield_EX,2) + pow(Y_MC_Stat_Err/Yield_MC,2) );

				pXS_MC4[Count4] = xs4->gGet_XS_Born(E0,Ep,Theta);
				//pXS_MC4[Count4] = xs4->gGet_XS_Born_ByKin(E0,Ep/1000.0,Theta);
			//	cerr << Form("   --Xbj = %f, XS_EX = %e, XS_MC=%e ", xbj, XS, pXS_MC4[Count4] ) <<endl;

				if(Min_Yield>=Yield_MC&&Yield_MC>1e3)
					Min_Yield = Yield_MC;
				if(Max_Yield<=Yield_MC)
					Max_Yield = Yield_MC;
				if(Min_Yield>=Yield_EX&&Yield_EX>1e3)
					Min_Yield = Yield_EX;
				if(Max_Yield<=Yield_EX)
					Max_Yield = Yield_EX;
				if(Min>=XS)
					Min = XS;	
				if(Max<=XS)
					Max = XS;
				if(Min_X>=xbj)
					Min_X = xbj;	
				if(Max_X<=xbj)
					Max_X = xbj;

				if(Min>= pXS_MC4[Count4])
					Min =  pXS_MC4[Count4];	
				if(Max<= pXS_MC4[Count4])
					Max =  pXS_MC4[Count4];

				if(Min_R>=pY_R4[Count4])
					Min_R = pY_R4[Count4];
				if(Max_R<=pY_R4[Count4])
					Max_R = pY_R4[Count4];

				Count4++;
			}
			if ( !run_4.good() )
				break;
		}
		run_4.close();
		delete xs4;
		cerr << "  Total XS Points for 4th  = " << Count4 <<endl;
	}
	/*}}}*/

	/*XS{{{*/

   Double_t Xbj_Max = 2.2;
   Double_t Xbj_Min = 0.7;
  
   if(Target=="H2")
	   Xbj_Max = 2.2;
   else if(Target=="He3")
	   Xbj_Max = 3.2;
   else
	   Xbj_Max = 4.2;

	if(Kin_Name1.Contains("5.0"))
		Xbj_Max = 1.8;

	TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	//TH2F* xs = new TH2F("xs","",1000,Xbj_Min,Xbj_Max,1000,Min/2.,Max*2.0);
	TH2F* xs = new TH2F("xs","",1000,Xbj_Min,Xbj_Max,1000,Min/1.1,Max*1.1);
	xs->SetStats(kFALSE);
	xs->SetXTitle("x_{bj}");
	xs->GetXaxis()->CenterTitle(1);
	xs->GetXaxis()->SetTitleFont(32);
	xs->GetXaxis()->SetTitleSize(0.06);
	xs->SetYTitle("#frac{d#sigma}{dE'd#Omega} (nb/MeV/sr)");
	xs->GetYaxis()->CenterTitle(1);
	xs->GetYaxis()->SetTitleFont(32);
	xs->GetYaxis()->SetTitleSize(0.04);
	xs->Draw(); 
	//gPad->SetLogy(1);
	gPad->SetLogy(1);

	TGraphErrors *pl_x1 = new TGraphErrors(Count1,pXbj1,pXS1,pXbj_Err1,pXS_Err1);
	pl_x1-> SetMarkerColor(2);
	pl_x1-> SetMarkerSize(1.0);
	pl_x1-> SetMarkerStyle(24);
	pl_x1-> GetXaxis()->SetTitle("x_{bj}");
	pl_x1-> GetYaxis()->SetTitle("#sigma");

	TGraphErrors *pl_x11 = new TGraphErrors(Count1,pXbj1,pXS_MC1,pXbj_Err1,pXbj_Err1);
	pl_x11-> SetMarkerColor(46);
	pl_x11-> SetMarkerSize(0.8);
	pl_x11-> SetLineColor(46);
	pl_x11-> SetLineWidth(2.0);
	pl_x11-> SetMarkerStyle(20);
	pl_x11-> GetXaxis()->SetTitle("x_{bj}");
	pl_x11-> GetYaxis()->SetTitle("#sigma_{mc}");

	TGraphErrors *pl_x3 = new TGraphErrors(Count3,pXbj3,pXS3,pXbj_Err3,pXS_Err3);
	TGraphErrors *pl_x13 = new TGraphErrors(Count3,pXbj3,pXS_MC3,pXbj_Err3,pXbj_Err3);
	if(bMore1){
		// pl_x1-> SetTitle(Form("C12,P_{0}=2.905 GeV/c, #theta=21^{o}"));
		pl_x3-> SetMarkerColor(4);
		pl_x3-> SetMarkerSize(1.0);
		pl_x3-> SetMarkerStyle(28);
		pl_x3-> GetXaxis()->SetTitle("x_{bj}");
		pl_x3-> GetYaxis()->SetTitle("#sigma");

		// pl_x1-> SetTitle(Form("C12,P_{0}=2.905 GeV/c, #theta=21^{o}"));
		pl_x13-> SetMarkerColor(9);
		pl_x13-> SetMarkerSize(0.8);
		pl_x13-> SetLineColor(9);
		pl_x13-> SetLineWidth(2.0);
		pl_x13-> SetMarkerStyle(20);
		pl_x13-> GetXaxis()->SetTitle("x_{bj}");
		pl_x13-> GetYaxis()->SetTitle("#sigma_{mc}");
	}
	/*}}}*/

	/*Yield{{{*/
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	//TH2F* y = new TH2F("y","",1000,Xbj_Min,Xbj_Max,1000,Min_Yield/2.,Max_Yield*2.0);
	TH2F* y = new TH2F("y","",1000,Xbj_Min,Xbj_Max,1000,Min_Yield/1.1,Max_Yield*1.1);
	y->SetStats(kFALSE);
	y->SetXTitle("x_{bj}");
	y->GetXaxis()->CenterTitle(1);
	y->GetXaxis()->SetTitleFont(32);
	y->GetXaxis()->SetTitleSize(0.06);
	y->SetYTitle("Yield");
	y->GetYaxis()->CenterTitle(1);
	y->GetYaxis()->SetTitleFont(32);
	y->GetYaxis()->SetTitleSize(0.04);
	y->Draw(); 
	//gPad->SetLogy(1);
	gPad->SetLogy(1);

	TGraphErrors *pl_y1 = new TGraphErrors(Count1,pXbj1,pYield_EX1,pXbj_Err1,pYield_EX_Err1);
	pl_y1-> SetMarkerColor(2);
	pl_y1-> SetMarkerSize(1.0);
	pl_y1-> SetMarkerStyle(20);
	pl_y1-> GetXaxis()->SetTitle("x_{bj}");
	pl_y1-> GetYaxis()->SetTitle("Yield");

	TGraphErrors *pl_y21 = new TGraphErrors(Count1,pXbj1,pYield_AL1,pXbj_Err1,pYield_AL_Err1);
	pl_y21-> SetMarkerColor(45);
	pl_y21-> SetMarkerSize(1.0);
	pl_y21-> SetMarkerStyle(24);
	pl_y21-> GetXaxis()->SetTitle("x_{bj}");
	pl_y21-> GetYaxis()->SetTitle("Yield");

	TGraphErrors *pl_y11 = new TGraphErrors(Count1,pXbj1,pYield_MC1,pXbj_Err1,pYield_MC_Err1);
	pl_y11-> SetMarkerColor(46);
	pl_y11-> SetMarkerSize(0.6);
	pl_y11-> SetLineColor(46);
	pl_y11-> SetLineWidth(2.0);
	pl_y11-> SetMarkerStyle(20);
	pl_y11-> GetXaxis()->SetTitle("x_{bj}");
	pl_y11-> GetYaxis()->SetTitle("Yield");

	TGraphErrors *pl_y3 = new TGraphErrors(Count3,pXbj3,pYield_EX3,pXbj_Err3,pYield_EX_Err3);
	TGraphErrors *pl_y13 = new TGraphErrors(Count3,pXbj3,pYield_MC3,pXbj_Err3,pYield_MC_Err3);
	TGraphErrors *pl_y23 = new TGraphErrors(Count3,pXbj3,pYield_AL3,pXbj_Err3,pYield_AL_Err3);
	if(bMore1){
		pl_y3-> SetMarkerColor(4);
		pl_y3-> SetMarkerSize(1.1);
		pl_y3-> SetMarkerStyle(29);
		pl_y3-> GetXaxis()->SetTitle("x_{bj}");
		pl_y3-> GetYaxis()->SetTitle("Yield");

		pl_y23-> SetMarkerColor(38);
		pl_y23-> SetMarkerSize(1.1);
		pl_y23-> SetMarkerStyle(30);
		pl_y23-> GetXaxis()->SetTitle("x_{bj}");
		pl_y23-> GetYaxis()->SetTitle("Yield");

		pl_y13-> SetMarkerColor(9);
		pl_y13-> SetMarkerSize(0.6);
		pl_y13-> SetLineColor(9);
		pl_y13-> SetLineWidth(2.0);
		pl_y13-> SetMarkerStyle(20);
		pl_y13-> GetXaxis()->SetTitle("x_{bj}");
		pl_y13-> GetYaxis()->SetTitle("Yield");
	}
	/*}}}*/ 

	/*Yield  Ratio{{{*/
	TCanvas *c3 = new TCanvas("c3","c3",1200,800);
	TH2F* r = new TH2F("r","",1000,Xbj_Min,Xbj_Max,1000,0.0, 2.0);
	r->SetStats(kFALSE);
	r->SetXTitle("x_{bj}");
	r->GetXaxis()->CenterTitle(1);
	r->GetXaxis()->SetTitleFont(32);
	r->GetXaxis()->SetTitleSize(0.06);
	r->SetYTitle("Yield_{EX}/Yield_{MC}");
	r->GetYaxis()->CenterTitle(1);
	r->GetYaxis()->SetTitleFont(32);
	r->GetYaxis()->SetTitleSize(0.04);
	r->Draw(); 

	TGraphErrors *pl_r1 = new TGraphErrors(Count1,pXbj1,pY_R1,pXbj_Err1,pY_R_Err1);
	TGraphErrors *pl_r3 = new TGraphErrors(Count3,pXbj3,pY_R3,pXbj_Err3,pY_R_Err3);
	if(bMore1){
		pl_r1-> SetMarkerColor(2);
		pl_r1-> SetMarkerSize(1.0);
		pl_r1-> SetMarkerStyle(24);
		pl_r1-> GetXaxis()->SetTitle("x_{bj}");
		pl_r1-> GetYaxis()->SetTitle("Yield_{EX}/Yield_{MC}");

		pl_r3-> SetMarkerColor(4);
		pl_r3-> SetMarkerSize(1.0);
		pl_r3-> SetMarkerStyle(28);
		pl_r3-> GetXaxis()->SetTitle("x_{bj}");
		pl_r3-> GetYaxis()->SetTitle("Yield_{EX}/Yield_{MC}");
	}
	/*}}}*/

	/*Third{{{*/
	TGraphErrors *pl_x2 = new TGraphErrors(Count2,pXbj2,pXS2,pXbj_Err2,pXS_Err2);
	TGraphErrors *pl_x12 = new TGraphErrors(Count2,pXbj2,pXS_MC2,pXbj_Err2,pXbj_Err2);
	TGraphErrors *pl_y2 = new TGraphErrors(Count2,pXbj2,pYield_EX2,pXbj_Err2,pYield_EX_Err2);
	TGraphErrors *pl_y12 = new TGraphErrors(Count2,pXbj2,pYield_MC2,pXbj_Err2,pYield_MC_Err2);
	TGraphErrors *pl_y22 = new TGraphErrors(Count2,pXbj2,pYield_AL2,pXbj_Err2,pYield_AL_Err2);
	TGraphErrors *pl_r2 = new TGraphErrors(Count2,pXbj2,pY_R2,pXbj_Err2,pY_R_Err2);

	if(bMore){
		pl_x2-> SetMarkerColor(6);
		pl_x2-> SetMarkerSize(1.0);
		pl_x2-> SetMarkerStyle(25);
		pl_x2-> GetXaxis()->SetTitle("x_{bj}");
		pl_x2-> GetYaxis()->SetTitle("#sigma");

		pl_x12-> SetMarkerColor(41);
		pl_x12-> SetMarkerSize(0.8);
		pl_x12-> SetLineColor(6);
		pl_x12-> SetLineWidth(2.0);
		pl_x12-> SetMarkerStyle(20);
		pl_x12-> GetXaxis()->SetTitle("x_{bj}");
		pl_x12-> GetYaxis()->SetTitle("#sigma_{mc}");

		pl_y2-> SetMarkerColor(6);
		pl_y2-> SetMarkerSize(1.0);
		pl_y2-> SetMarkerStyle(21);
		pl_y2-> GetXaxis()->SetTitle("x_{bj}");
		pl_y2-> GetYaxis()->SetTitle("Yield");

		pl_y22-> SetMarkerColor(28);
		pl_y22-> SetMarkerSize(1.0);
		pl_y22-> SetMarkerStyle(25);
		pl_y22-> GetXaxis()->SetTitle("x_{bj}");
		pl_y22-> GetYaxis()->SetTitle("Yield");

		pl_y12-> SetMarkerColor(41);
		pl_y12-> SetMarkerSize(0.6);
		pl_y12-> SetLineColor(41);
		pl_y12-> SetLineWidth(2.0);
		pl_y12-> SetMarkerStyle(20);
		pl_y12-> GetXaxis()->SetTitle("x_{bj}");
		pl_y12-> GetYaxis()->SetTitle("Yield");

		pl_r2-> SetMarkerColor(6);
		pl_r2-> SetMarkerSize(1.0);
		pl_r2-> SetMarkerStyle(25);
		pl_r2-> GetXaxis()->SetTitle("x_{bj}");
		pl_r2-> GetYaxis()->SetTitle("Yield_{EX}/Yield_{MC}");
	}
	/*}}}*/

	/*Forth{{{*/
	TGraphErrors *pl_x4 = new TGraphErrors(Count4,pXbj4,pXS4,pXbj_Err4,pXS_Err4);
	TGraphErrors *pl_x14 = new TGraphErrors(Count4,pXbj4,pXS_MC4,pXbj_Err4,pXbj_Err4);
	TGraphErrors *pl_y4 = new TGraphErrors(Count4,pXbj4,pYield_EX4,pXbj_Err4,pYield_EX_Err4);
	TGraphErrors *pl_y14 = new TGraphErrors(Count4,pXbj4,pYield_MC4,pXbj_Err4,pYield_MC_Err4);
	TGraphErrors *pl_y24 = new TGraphErrors(Count4,pXbj4,pYield_AL4,pXbj_Err4,pYield_AL_Err4);
	TGraphErrors *pl_r4 = new TGraphErrors(Count4,pXbj4,pY_R4,pXbj_Err4,pY_R_Err4);

	if(bMore2){
		//pl_x4-> SetTitle(form("c12,p_{0}=3.055 gev/c, #theta=21^{o}"));
		pl_x4-> SetMarkerColor(8);
		pl_x4-> SetMarkerSize(1.0);
		pl_x4-> SetMarkerStyle(26);
		pl_x4-> GetXaxis()->SetTitle("x_{bj}");
		pl_x4-> GetYaxis()->SetTitle("#sigma");

		//pl_x14-> SetTitle(form("c12,p_{0}=3.055 gev/c, #theta=21^{o}"));
		pl_x14-> SetMarkerColor(3);
		pl_x14-> SetMarkerSize(0.8);
		pl_x14-> SetLineColor(3);
		pl_x14-> SetLineWidth(2.0);
		pl_x14-> SetMarkerStyle(20);
		pl_x14-> GetXaxis()->SetTitle("x_{bj}");
		pl_x14-> GetYaxis()->SetTitle("#sigma_{mc}");

		pl_y4-> SetMarkerColor(8);
		pl_y4-> SetMarkerSize(1.0);
		pl_y4-> SetMarkerStyle(22);
		pl_y4-> GetXaxis()->SetTitle("x_{bj}");
		pl_y4-> GetYaxis()->SetTitle("Yield");

		pl_y24-> SetMarkerColor(30);
		pl_y24-> SetMarkerSize(1.0);
		pl_y24-> SetMarkerStyle(26);
		pl_y24-> GetXaxis()->SetTitle("x_{bj}");
		pl_y24-> GetYaxis()->SetTitle("Yield");

		pl_y14-> SetMarkerColor(3);
		pl_y14-> SetMarkerSize(0.6);
		pl_y14-> SetLineColor(3);
		pl_y14-> SetLineWidth(2.0);
		pl_y14-> SetMarkerStyle(20);
		pl_y14-> GetXaxis()->SetTitle("x_{bj}");
		pl_y14-> GetYaxis()->SetTitle("Yield");

		pl_r4-> SetMarkerColor(8);
		pl_r4-> SetMarkerSize(1.0);
		pl_r4-> SetMarkerStyle(26);
		pl_r4-> GetXaxis()->SetTitle("x_{bj}");
		pl_r4-> GetYaxis()->SetTitle("Yield_{EX}/Yield_{MC}");
	}
	/*}}}*/

	/*Plot{{{*/

	/*XS{{{*/
	c1->cd();
	//gPad->SetLogy(1);
	gPad->SetLogy(1);

	TLegend *l1 = new TLegend(0.65,0.59,0.80, 0.94,Form("%s Target",Target.Data()));
	l1->SetBorderSize(0);
	l1->SetTextSize(0.03);
	l1->SetTextFont(22);
	l1->AddEntry(pl_x1, Form("XS_{EX}:  %s (HRS-%s)",Kin_Name1.Data(), Arm1.Data()),"p");
	pl_x1->Draw("P");

	if(bMore1){
		l1->AddEntry(pl_x3, Form("XS_{EX}:  %s (HRS-%s)",Kin_Name3.Data(), Arm3.Data()),"p");
		pl_x3->Draw("P");
	}

	l1->AddEntry(pl_x11,Form("XS_{MC}:  %s",Kin_Name1.Data()),"pl");
	pl_x11->Draw("L");

	if(bMore){
		l1->AddEntry(pl_x2, Form("XS_{EX}:  %s (HRS-%s)",Kin_Name2.Data(), Arm2.Data()),"p");
		pl_x2->Draw("P");
	}
	if(bMore2){
		l1->AddEntry(pl_x4, Form("XS_{EX}:  %s (HRS-%s)",Kin_Name4.Data(), Arm4.Data()),"p");
		pl_x4->Draw("P");
	}
	if(bMore){
		l1->AddEntry(pl_x12,Form("XS_{MC}:  %s",Kin_Name2.Data()),"pl");
		pl_x12->Draw("L");
	}

	l1->Draw();
	c1->Print(Form("%s_XS_Kin%s_%s_xbj.pdf",Target.Data(),Kin.Data(),Comment.Data()));
	c1->Print(Form("%s_XS_Kin%s_%s_xbj.png",Target.Data(),Kin.Data(),Comment.Data()));
 	/*}}}*/

	/*Yield{{{*/
	c2->cd();
	//gPad->SetLogy(1);
	gPad->SetLogy(1);

	TLegend *l2 = new TLegend(0.65,0.59,0.90, 0.94,Form("%s Target",Target.Data()));
	l2->SetBorderSize(0);
	l2->SetTextSize(0.03);
	l2->SetTextFont(22);

	l2->AddEntry(pl_y1, Form("Y_{EX}:  %s (HRS-%s)",Kin_Name1.Data(), Arm1.Data()),"p");
	pl_y1->Draw("P");
	l2->AddEntry(pl_y11,Form("Y_{MC}^{%s}:  %s",Arm1.Data(),Kin_Name1.Data()),"pl");
	pl_y11->Draw("Lp");

	if(bMore1){
		l2->AddEntry(pl_y3, Form("Y_{EX}:  %s (HRS-%s)",Kin_Name3.Data(), Arm3.Data()),"p");
		pl_y3->Draw("P");
		l2->AddEntry(pl_y13,Form("Y_{MC}^{%s}:  %s",Arm3.Data(),Kin_Name3.Data()),"pl");
		pl_y13->Draw("Lp");
	}

	if(bMore){
		l2->AddEntry(pl_y2, Form("Y_{EX}:  %s (HRS-%s)",Kin_Name2.Data(), Arm2.Data()),"p");
		pl_y2->Draw("P");
		l2->AddEntry(pl_y12,Form("Y_{MC}^{%s}:  %s",Arm2.Data(),Kin_Name2.Data()),"pl");
		pl_y12->Draw("Lp");	
	}
	if(bMore2){
		l2->AddEntry(pl_y4, Form("Y_{EX}:  %s (HRS-%s)",Kin_Name4.Data(), Arm4.Data()),"p");
		pl_y4->Draw("P");
		l2->AddEntry(pl_y14,Form("Y_{MC}^{%s}:  %s",Arm4.Data(),Kin_Name4.Data()),"pl");
		pl_y14->Draw("Lp");	
	}

	if(Target=="H2"||Target=="He3"||Target=="He4"){
		l2->AddEntry(pl_y21, Form("Y_{AL}:  %s (HRS-%s)",Kin_Name1.Data(), Arm1.Data()),"p");
		pl_y21->Draw("P");

		if(bMore1){
			l2->AddEntry(pl_y23, Form("Y_{AL}:  %s (HRS-%s)",Kin_Name3.Data(), Arm3.Data()),"p");
			pl_y23->Draw("P");
		}
		if(bMore){
			l2->AddEntry(pl_y22, Form("Y_{AL}:  %s (HRS-%s)",Kin_Name2.Data(), Arm2.Data()),"p");
			pl_y22->Draw("P");
		}
		if(bMore2){
			l2->AddEntry(pl_y24, Form("Y_{AL}:  %s (HRS-%s)",Kin_Name4.Data(), Arm4.Data()),"p");
			pl_y24->Draw("P");
		}
	}


	l2->Draw();
	c2->Print(Form("%s_Yield_Kin%s_%s_xbj.pdf",Target.Data(),Kin.Data(),Comment.Data()));
	c2->Print(Form("%s_Yield_Kin%s_%s_xbj.png",Target.Data(),Kin.Data(),Comment.Data()));
	/*}}}*/

	/*Yield Ratio{{{*/
	c3->cd();
	gPad->SetLogy(0);

	TMultiGraph *mg3 = new TMultiGraph();
	TLegend *l3 = new TLegend(0.65,0.59,0.90, 0.94,Form("%s Target",Target.Data()));
	l3->SetBorderSize(0);
	l3->SetTextSize(0.025);
	l3->SetTextFont(22);

	l3->AddEntry(pl_r1, Form("Y_{EX}/Y_{MC}:  %s (HRS-%s)",Kin_Name1.Data(), Arm1.Data()),"p");
	mg3->Add(pl_r1);

	if(bMore){
		l3->AddEntry(pl_r2, Form("Y_{EX}/Y_{MC}:  %s (HRS-%s)",Kin_Name2.Data(), Arm2.Data()),"p");
		mg3->Add(pl_r2);
	}

	if(bMore1){	
		l3->AddEntry(pl_r3, Form("Y_{EX}/Y_{MC}:  %s (HRS-%s)",Kin_Name3.Data(), Arm3.Data()),"p");
		mg3->Add(pl_r3);
	}
	if(bMore2){
		l3->AddEntry(pl_r4, Form("Y_{EX}/Y_{MC}:  %s (HRS-%s)",Kin_Name4.Data(), Arm4.Data()),"p");
		mg3->Add(pl_r4);
	}
	mg3->Draw("P");
	l3->Draw();
	c3->Print(Form("%s_YR_Kin%s_%s_xbj.pdf",Target.Data(),Kin.Data(),Comment.Data()));
	c3->Print(Form("%s_YR_Kin%s_%s_xbj.png",Target.Data(),Kin.Data(),Comment.Data()));
//	c3->Print(Form("%s_YR_Kin%s_xbj.png",Target.Data(),Kin.Data()));
	/*}}}*/

	/*}}}*/
}
/*}}}*/


int main(int argc, char** argv){



	int err = getargs(argc,argv);

	if(err>=0){
		aFileName1 = gFile1;
		aFileName2 = gFile2;
		aFileName3 = gFile3;
		aFileName4 = gFile4;
		Target = gTgt;
		Comment = gCom;
		Yield_Plot(); 
	}
	else
		cerr<<"*** ERROR, I can not recoginize the parameters!"<<endl;

	return 0;
}

/*int getargs(int argc,char** argv){{{*/
int getargs(int argc,char** argv)
{
	char* argptr;
	bool noStop;  
	int err=-1;

	for(int i=0;i<argc;i++){
		argptr = argv[i];
		if(*argptr=='-'){
			argptr++;
			switch (*argptr){
				case '1':
					if(*(++argptr))
						gFile1 = argptr;
					err=0;
					break;
				case '2':
					if(*(++argptr))
						gFile2 = argptr;
					err=0;
					break;
				case '3':
					if(*(++argptr))
						gFile3 = argptr;
					err=0;
					break;
				case '4':
					if(*(++argptr))
						gFile4 = argptr;
					err=0;
					break;
				case 'T':
					if(*(++argptr))
						gTgt = argptr;
					err=0;
					break;
				case 'C':
					if(*(++argptr))
						gCom = argptr;
					err=0;
					break;
				case 'h':
					cerr <<"================================================"<<endl;
					cerr <<"Option: " <<endl;
					cerr <<"-N[*.dat] Input File,N=1,2,3,4, *.dat =x for nothing" <<endl;
					cerr <<"-T[H2,He3,He4,C12,Ca40,Ca48] Target Name" <<endl;
					cerr <<"-C[] Comment" <<endl;
					cerr <<endl<<"================================================"<<endl;
					noStop=false;
					err=-1;
					goto OUT;
					break;
				default:
					cerr <<"================================================"<<endl;
					cerr <<"Option: " <<endl;
					cerr <<"-N[*.dat] Input File,N=1,2,3,4, *.dat =x for nothing" <<endl;
					cerr <<"-T[H2,He3,He4,C12,Ca40,Ca48] Target Name" <<endl;
					cerr <<"-C[] Comment" <<endl;
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
