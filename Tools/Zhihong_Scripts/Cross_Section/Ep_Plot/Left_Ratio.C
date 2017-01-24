#include "XGT2_Main.h"
#include "/work/halla/e08014/disk1/yez/Ep_Cross_Section/SRC/XGT2_XEMC.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"

//const double He4_Norm =1.0/0.864;
const double He4_Norm =1.0;

/*Get_Target_Name{{{*/
TString gGet_Target_Name(const TString& aName){
	TString aTarget_Name="MM";

	if(aName.Contains("H2")||aName.Contains("D2"))
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
	else if(aKin.Contains("5.0") && !aKin.Contains("5.05"))
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

/*Cuts on the edge of Xbj{{{*/
// Cut_L_Min Cut_L_Max Cut_R_Min Cut_R_Max
const double E0 = 3.356;
const double Xbj_Cut_21_1[4] = {1.25, 2.15,1.25,2.15}; 
const double Xbj_Cut_21_2[4] = {1.65, 3.00,1.65,3.00}; 
const double Xbj_Cut_23_1[4] = {1.30, 2.15,1.30,2.15}; 
const double Xbj_Cut_23_2[4] = {1.85, 3.00,1.85,3.00};
const double Xbj_Cut_QE_1[4] = {0.85, 1.20,0.85,1.20}; 
const double Xbj_Cut_QE_2[4] = {1.05, 1.55,1.05,1.55}; 
const double Xbj_Cut_25_1[4] = {1.35, 2.10,1.35,2.10};
const double Xbj_Cut_25_2[4] = {2.00, 3.00,2.00,3.00}; 
const double Xbj_Cut_28_1[4] = {1.80, 3.00,1.80,3.00}; 
/*
const double Xbj_Cut_21_1[4] = {1.25, 1.75,1.25,1.75}; 
const double Xbj_Cut_21_2[4] = {2.10, 3.00,2.10,3.00}; 
const double Xbj_Cut_23_1[4] = {1.45, 2.05,1.45,2.05}; 
const double Xbj_Cut_23_2[4] = {1.95, 3.00,1.95,3.00};
const double Xbj_Cut_QE_1[4] = {0.85, 1.15,0.85,1.15}; 
const double Xbj_Cut_QE_2[4] = {1.05, 1.50,1.05,1.50}; 
const double Xbj_Cut_25_1[4] = {1.40, 2.10,1.40,2.10};
const double Xbj_Cut_25_2[4] = {2.10, 3.00,2.10,3.00}; 
const double Xbj_Cut_28_1[4] = {1.90, 3.00,1.90,3.00}; 
*/
double Ep_Cut_21_1[4] = {2.80, 3.00,2.80,3.00}; 
double Ep_Cut_21_2[4] = {2.94, 3.18,2.94,3.18}; 
double Ep_Cut_23_1[4] = {2.74, 2.94,2.74,2.94}; 
double Ep_Cut_23_2[4] = {2.92, 3.16,2.92,3.16};
double Ep_Cut_QE_1[4] = {2.42, 2.58,2.42,2.58}; 
double Ep_Cut_QE_2[4] = {2.56, 2.74,2.56,2.74}; 
double Ep_Cut_25_1[4] = {2.70, 2.88,2.70,2.88};
double Ep_Cut_25_2[4] = {2.88, 3.12,2.88,3.12}; 
double Ep_Cut_28_1[4] = {2.74, 2.98,2.74,2.98};

inline void gGet_Ep_Edge(){
	double Theta = 0.0, SinSQ = 0.0;
	for(int i=0;i<4;i++){
		Theta = 21.00;
		SinSQ = pow(sin(Theta*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2
		Ep_Cut_21_1[i]=E0/(1.0+4.0*E0*SinSQ/(2.0*PROTON_MASS*Xbj_Cut_21_1[i])); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
		Ep_Cut_21_2[i]=E0/(1.0+4.0*E0*SinSQ/(2.0*PROTON_MASS*Xbj_Cut_21_2[i])); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
		Theta = 23.00;
		SinSQ = pow(sin(Theta*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2
		Ep_Cut_23_1[i]=E0/(1.0+4.0*E0*SinSQ/(2.0*PROTON_MASS*Xbj_Cut_23_1[i])); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
		Ep_Cut_23_2[i]=E0/(1.0+4.0*E0*SinSQ/(2.0*PROTON_MASS*Xbj_Cut_23_2[i])); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
		Theta = 25.00;
		SinSQ = pow(sin(Theta*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2
		Ep_Cut_QE_1[i]=E0/(1.0+4.0*E0*SinSQ/(2.0*PROTON_MASS*Xbj_Cut_QE_1[i])); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
		Ep_Cut_QE_2[i]=E0/(1.0+4.0*E0*SinSQ/(2.0*PROTON_MASS*Xbj_Cut_QE_2[i])); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
		Ep_Cut_25_1[i]=E0/(1.0+4.0*E0*SinSQ/(2.0*PROTON_MASS*Xbj_Cut_25_1[i])); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
		Ep_Cut_25_2[i]=E0/(1.0+4.0*E0*SinSQ/(2.0*PROTON_MASS*Xbj_Cut_25_2[i])); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
		Theta = 28.00;
		SinSQ = pow(sin(Theta*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2
		Ep_Cut_28_1[i]=E0/(1.0+4.0*E0*SinSQ/(2.0*PROTON_MASS*Xbj_Cut_28_1[i])); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
	}
}
/*}}}*/

/*Get He3 RadCorr{{{*/
inline double gGet_RadCorr_He3(const double aXbj){
  double aRadCor = 1-0.1*(aXbj-2.45)+0.025*pow((aXbj-2.45),2)-0.72*pow((aXbj-2.45),3);
  //return 1.0/aRadCor;
  return 1.0;
}
/*}}}*/

void Merge(const TString& aInputFileName, double* aXbj, double* aR_EX, double* aR_EX_Err, double* aR_MC, int* aN);
void E02019(const TString& aInputFileName, double aTheta, double* aXbj, double* aR_EX, double* aR_EX_Err, int* aN);
void E02019_A_He3(const TString& aTarget, const double aTheta,double* aXbj, double* aR_EX, double* aR_EX_Err_Up,double* aR_EX_Err_Down, int* aN);
void E02019_He4_He3(const TString& aTarget, const double aTheta,double* aXbj, double* aR_EX, double* aR_EX_Err_Up,double* aR_EX_Err_Down, int* aN);
void CLAS_He4_He3( double* aXbj_CS, double* aR_EX_CS, double* aR_EX_Err_CS, int *aN_CS);
TString Target0,Target,Com0,Com;
/*Comment{{{*/
/*Format in aFileName:
  Target 1:
  line#1 --> The first  setting of each angle on  L-HRS, like L_Kin3.1 ...
  line#2 --> The second setting of each angle on  L-HRS, like L_Kin3.2 ...
  line#3 --> The first  setting of each angle on  R-HRS, like R_Kin3.1 ...
  line#4 --> The second setting of each angle on  R-HRS, like R_Kin3.2 ...
  Target 2:  
  line#5 --> The first  setting of each angle on  L-HRS, like L_Kin3.1 ...
  line#6 --> The second setting of each angle on  L-HRS, like L_Kin3.2 ...

  line#7 --> The first  setting of each angle on  R-HRS, like R_Kin3.1 ...
  line#8 --> The second setting of each angle on  R-HRS, like R_Kin3.2 ...
  */
/*}}}*/

/*main{{{*/
int main()
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
	//    gStyle->SetOptStat(1);
	gStyle->SetOptFit(1);
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

	cerr<<" 1 vs 2 (He3 or He4 or C12) for 1 = "; cin >> Target0;
	cerr<<" 1 vs 2 (H2  or He3 or He4) for 2 = "; cin >> Target;
	cerr<<" Com for 1 = "; cin >> Com0;
	cerr<<" Com for 2 = "; cin >> Com;

	gGet_Ep_Edge();//Get Ep Range from Xbj edge
	int A1 = 0, A2 = 0, Z1 = 0, Z2=0;
	/*Get Target Info 1{{{*/
	if(Target0 == "H2"||Target == "D2"){
		A1 = 2; Z1 = 1;  
	}
	else if(Target0 == "He3"){
		A1 = 3; Z1 = 2;  
	}
	else if(Target0 == "He4"){
		A1 = 4; Z1 = 2; 
	}
	else if(Target0 == "Be"){
		A1 = 9; Z1 = 4; 
	}
	else if(Target0 == "C12"){
		A1 = 12; Z1 = 6; 
	}
	else if(Target0 == "Al"){
		A1 = 27; Z1 = 13; 
	}
	else if(Target0 == "Ca40"){
		A1 = 40; Z1 = 20;
	}
	else if(Target0 == "Ca48"){
		A1 = 48; Z1 = 20; 
	}
	else
		cerr<<" **** ERROR, I don't know the target!!!"<<endl;
	/*}}}*/
	/*Get Target Info 2{{{*/
	if(Target == "H2"||Target == "D2"){
		A2 = 2; Z2 = 1;  
	}
	else if(Target == "He3"){
		A2 = 3; Z2 = 2;  
	}
	else if(Target == "He4"){
		A2 = 4; Z2 = 2; 
	}
	else if(Target == "Be"){
		A2 = 9; Z2 = 4; 
	}
	else if(Target == "C12"){
		A2 = 12; Z2 = 6; 
	}
	else if(Target == "Al"){
		A2 = 27; Z2 = 13; 
	}
	else if(Target == "Ca40"){
		A2 = 40; Z2 = 20;
	}
	else if(Target == "Ca48"){
		A2 = 48; Z2 = 20; 
	}
	else
		cerr<<" **** ERROR, I don't know the target!!!"<<endl;
	/*}}}*/

	/*Get Ratio{{{*/
	//bool bMC = kTRUE;
	bool bMC = kFALSE;
	double Theta =21;
	double Zero[100] = {100*0.0};

	double xbj_21[100], R_EX_21[100], R_EX_Err_21[100], R_MC_21[100];
	int N_21 = 0;
	Merge(Form("./input/R21_%s_%s.rt",Target0.Data(),Target.Data()), xbj_21, R_EX_21, R_EX_Err_21, R_MC_21, &N_21);
	double Q2_21_D = 4.0*E0*Ep_Cut_21_1[0]*pow(TMath::Sin(21.*3.1415926/180.0*0.5),2);
	double Q2_21_U = 4.0*E0*Ep_Cut_21_2[3]*pow(TMath::Sin(21.*3.1415926/180.0*0.5),2);
	
	double xbj_23[100], R_EX_23[100], R_EX_Err_23[100], R_MC_23[100];
	int N_23 = 0;
	if(!(Target=="H2"))
		Merge(Form("./input/R23_%s_%s.rt",Target0.Data(),Target.Data()), xbj_23, R_EX_23, R_EX_Err_23, R_MC_23, &N_23);
	double Q2_23_D = 4.0*E0*Ep_Cut_23_1[0]*pow(TMath::Sin(23.*3.1415926/180.0*0.5),2);
	double Q2_23_U = 4.0*E0*Ep_Cut_23_2[3]*pow(TMath::Sin(23.*3.1415926/180.0*0.5),2);
	
	double xbj_QE[100], R_EX_QE[100], R_EX_Err_QE[100], R_MC_QE[100];
	int N_QE = 0;
	if(!(Target=="H2"))
		Merge(Form("./input/RQE_%s_%s.rt",Target0.Data(),Target.Data()), xbj_QE, R_EX_QE, R_EX_Err_QE, R_MC_QE, &N_QE);
	double Q2_25_D = 4.0*E0*Ep_Cut_QE_1[0]*pow(TMath::Sin(25.*3.1415926/180.0*0.5),2);
	
	double xbj_25[100], R_EX_25[100], R_EX_Err_25[100], R_MC_25[100];
	int N_25 = 0;
	Merge(Form("./input/R25_%s_%s.rt",Target0.Data(),Target.Data()), xbj_25, R_EX_25, R_EX_Err_25, R_MC_25, &N_25);
	double Q2_25_U = 4.0*E0*Ep_Cut_25_2[3]*pow(TMath::Sin(25.*3.1415926/180.0*0.5),2);

	double xbj_28[100], R_EX_28[100], R_EX_Err_28[100], R_MC_28[100];
	int N_28 = 0;
	if(Target0=="C12" && Target=="He3"){
		Merge(Form("./input/R28_%s_%s.rt",Target0.Data(),Target.Data()), xbj_28, R_EX_28, R_EX_Err_28, R_MC_28, &N_28);
	}
	double Q2_28_D = 4.0*E0*Ep_Cut_28_1[0]*pow(TMath::Sin(28.*3.1415926/180.0*0.5),2);
	double Q2_28_U = 4.0*E0*Ep_Cut_28_1[3]*pow(TMath::Sin(28.*3.1415926/180.0*0.5),2);
	
	//Add E02019 Data
	double xbj_NF[1000], R_EX_NF[1000], R_EX_Err_NF[1000],R_EX_Err_Up[1000],R_EX_Err_Down[1000];
	//Only fit one angle 18,22,26,32,40,50, 0->for all:
	//if(!(Target0=="C12"&&Target=="He4")){
	int N_NF = 0;  Double_t Theta_NF=18;

	if(Target=="He3" && !(Target0=="He4"))
		E02019_A_He3(Target0.Data(), Theta_NF,xbj_NF, R_EX_NF, R_EX_Err_Up,  R_EX_Err_Down, &N_NF);
	if(Target=="He3" && Target0=="He4")
		E02019_He4_He3(Target0.Data(), Theta_NF,xbj_NF, R_EX_NF, R_EX_Err_Up,  R_EX_Err_Down, &N_NF);
	
	if(0){
		Theta_NF = 18; //0->all
		E02019(Form("./input/RNF_%s_%s.rt",Target0.Data(),Target.Data()), Theta_NF,xbj_NF, R_EX_NF, R_EX_Err_NF, &N_NF);
		Theta_NF = 22; //0->all
		E02019(Form("./input/RNF_%s_%s.rt",Target0.Data(),Target.Data()), Theta_NF,xbj_NF, R_EX_NF, R_EX_Err_NF, &N_NF);
		Theta_NF = 26; //0->all
		E02019(Form("./input/RNF_%s_%s.rt",Target0.Data(),Target.Data()), Theta_NF,xbj_NF, R_EX_NF, R_EX_Err_NF, &N_NF);
	}
	cerr<<Form("@@@@@@@@@@ N_21 = %d, N_23 = %d, N_25 = %d, N_QE = %d, N_28 = %d, N_NF = %d", N_21,N_23,N_25,N_QE,N_28,N_NF)<<endl;

	//Add CLAS Data
	double xbj_CS[1000], R_EX_CS[1000], R_EX_Err_CS[1000];
	int N_CS = 0;  
	if(Target=="He3" && Target0=="He4")
		CLAS_He4_He3( xbj_CS, R_EX_CS, R_EX_Err_CS, &N_CS);
	cerr<<Form("@@@@@@@@@@ N_CLAS = %d", N_CS)<<endl;
/*}}}*/

	/*Save Ratio Values into files{{{*/
/*    ofstream outf_21(Form("Ratio_21_%s_%s_%s_%s.dat", Target0.Data(),Target.Data(),Com0.Data(),Com.Data() ));
    for(int i = 0;i<N_21;i++)
		outf_21 <<Form("%5d %10.8f, %12.8f, %12.8f",21, xbj_21[i], R_EX_21[i],R_EX_Err_21[i])<<endl;
    outf_21.close();
 
	if(!(Target=="H2")){
		ofstream outf_23(Form("Ratio_23_%s_%s_%s_%s.dat", Target0.Data(),Target.Data(),Com0.Data(),Com.Data() ));
		for(int i = 0;i<N_23;i++)
			outf_23 <<Form("%5d %10.8f, %12.8f, %12.8f",23, xbj_23[i], R_EX_23[i],R_EX_Err_23[i])<<endl;
		outf_23.close();
	}

	ofstream outf_25(Form("Ratio_25_%s_%s_%s_%s.dat", Target0.Data(),Target.Data(),Com0.Data(),Com.Data() ));
	if(!(Target=="H2")){
		for(int i = 0;i<N_QE;i++)
			outf_25 <<Form("%5d %10.8f, %12.8f, %12.8f",25, xbj_QE[i], R_EX_QE[i],R_EX_Err_QE[i])<<endl;
	}
	for(int i = 0;i<N_25;i++)
		outf_25 <<Form("%5d %10.8f, %12.8f, %12.8f",25, xbj_25[i], R_EX_25[i],R_EX_Err_25[i])<<endl;
    outf_25.close();

	if(Target0=="C12"&&Target=="He3"){
		ofstream outf_28(Form("Ratio_28_%s_%s_%s_%s.dat", Target0.Data(),Target.Data(),Com0.Data(),Com.Data() ));
		for(int i = 0;i<N_28;i++)
			outf_28 <<Form("%5d %10.8f, %12.8f, %12.8f",28, xbj_28[i], R_EX_28[i],R_EX_Err_28[i])<<endl;
		outf_28.close();
	}
*/
	/*}}}*/

	/*Plot{{{*/
	cerr<<endl<<"@@@@@@@@ Now starting to make some plots ..."<<endl;

	gStyle->SetOptStat(0);

	double RMAX = -1000,RMIN=0.4;
	for(int i=0; i<N_21; i++)
		RMAX = RMAX?R_EX_21[i]:RMAX<R_EX_21[i];
	for(int i=0; i<N_23; i++)
		RMAX = RMAX?R_EX_23[i]:RMAX<R_EX_23[i];
	for(int i=0; i<N_25; i++)
		RMAX = RMAX?R_EX_25[i]:RMAX<R_EX_25[i];
	if(Target=="He3"){
		for(int i=0; i<N_28; i++)
			RMAX = RMAX?R_EX_28[i]:RMAX<R_EX_28[i];
	}
	for(int i=0; i<N_21; i++)
		RMAX = RMAX?R_MC_21[i]:RMAX<R_MC_21[i];
	for(int i=0; i<N_23; i++)
		RMAX = RMAX?R_MC_23[i]:RMAX<R_MC_23[i];
	for(int i=0; i<N_25; i++)
		RMAX = RMAX?R_MC_25[i]:RMAX<R_MC_25[i];
	if(Target=="He3"){
		for(int i=0; i<N_28; i++)
			RMAX = RMAX?R_MC_28[i]:RMAX<R_MC_28[i];
	}

	double XMax = 3.2;
	if(Target=="He4")
		XMax=3.2;
	if(Target=="He3")
		XMax=3.2;
	if(Target=="H2")
		XMax=2.2;

	if(Target=="He3" &&Target0=="C12")
		RMAX = 20.0;
	if(Target=="He3" &&Target0=="He4")
		RMAX = 8.0;
	if(Target=="He4" &&Target0=="C12"){
		RMIN = 0.8; RMAX = 6.0;
	}
	if(Target=="Ca40" &&Target0=="Ca48"){
		RMAX = 1.25; RMIN=0.75;
	}
	if(Target=="C12" &&Target0=="Ca40")
		RMAX = 1.5;
	if(Target=="C12" &&Target0=="Ca48")
		RMAX = 1.5;
	if(Target=="H2" &&Target0=="He3")
		RMAX = 3.0;
	if(Target=="H2" &&Target0=="He4")
		RMAX = 5.1;
	if(Target=="H2" &&Target0=="C12"){
		RMIN=0.1; RMAX = 7.5;
	}
	if(Target=="H2" &&Target0=="Ca40"){
		RMIN=0.1; RMAX = 9.5;
	}
	if(Target=="H2" &&Target0=="Ca48"){
		RMIN=0.1; RMAX = 9.5;
	}

	TCanvas *c3 = new TCanvas("c3","c3",1200,700);
	c3->cd();
	TH2F *h3 = new TH2F("h3","",300,0.7,XMax,300,RMIN,RMAX);
	h3->SetStats(kFALSE);
	h3->SetXTitle("x");
	h3->GetXaxis()->CenterTitle(1);
	h3->GetXaxis()->SetTitleFont(22);
	h3->GetXaxis()->SetTitleSize(0.08);
	h3->GetXaxis()->SetTitleOffset(0.7);
	h3->GetYaxis()->SetTitle(Form("#frac{#sigma_{%s}}{#sigma_{%s}} #times #frac{%d}{%d}",Target0.Data(),Target.Data(),A2,A1));
	h3->GetYaxis()->CenterTitle(1);
	h3->GetYaxis()->SetTitleFont(22);
	h3->GetYaxis()->SetTitleSize(0.07);
	h3->GetYaxis()->SetTitleOffset(0.8);
	h3->GetYaxis()->CenterTitle();
	gPad->SetLogy(0);
	h3->Draw();

	if(Target=="He3" &&Target0=="C12")
		//RMAX = 10.0;
		RMAX = 5.0;
	if(Target=="He3" &&Target0=="He4")
		RMAX = 8.0;
	if(Target=="He4" &&Target0=="C12")
		RMAX = 5.5;
	if(Target=="Ca40" &&Target0=="Ca48"){
		RMIN = 0.8;		RMAX = 1.25;
	}
	if(Target=="C12" &&Target0=="Ca40"){
		RMAX = 1.5; RMIN=0.7;
	}
	if(Target=="C12" &&Target0=="Ca48"){
		RMAX = 1.5; RMIN=0.7;
	}
	if(Target=="H2" &&Target0=="C12"){
		 RMIN=0.1;
	}

	TCanvas *c4 = new TCanvas("c4","c4",1200,700);
	c4->cd();
	TH2F *h4 = new TH2F("h4","",300,0.7,XMax,300,RMIN,RMAX);
	h4->SetStats(kFALSE);
	h4->SetXTitle("x");
	h4->GetXaxis()->CenterTitle(1);
	h4->GetXaxis()->SetTitleFont(22);
	h4->GetXaxis()->SetTitleSize(0.08);
	h4->GetXaxis()->SetTitleOffset(0.7);
	h4->GetYaxis()->SetTitle(Form("#frac{#sigma_{%s}}{#sigma_{%s}} #times #frac{%d}{%d}",Target0.Data(),Target.Data(),A2,A1));
	h4->GetYaxis()->CenterTitle(1);
	h4->GetYaxis()->SetTitleFont(22);
	h4->GetYaxis()->SetTitleSize(0.07);
	h4->GetYaxis()->SetTitleOffset(0.8);
	h4->GetYaxis()->CenterTitle();
	gPad->SetLogy(0);
	h4->Draw();

	TLegend *l1 = new TLegend(0.17,0.60,0.50, 0.94,Form("#sigma_{%s} / #sigma_{%s} Ratio: ",Target0.Data(),Target.Data()));
	l1->SetTextSize(0.043);
	l1->SetTextFont(32);

	TGraphErrors *ex21 = new TGraphErrors(N_21,xbj_21, R_EX_21, Zero,R_EX_Err_21);
	ex21->SetMarkerStyle(26);
	ex21->SetMarkerSize(1.5);
	ex21->SetMarkerColor(2);
	ex21->SetLineColor(2);

	TGraphErrors *ex23 = new TGraphErrors(N_23,xbj_23, R_EX_23, Zero,R_EX_Err_23);
	ex23->SetMarkerStyle(25);
	ex23->SetMarkerSize(1.5);
	ex23->SetMarkerColor(4);
	ex23->SetLineColor(4);

	TGraphErrors *exQE = new TGraphErrors(N_QE,xbj_QE, R_EX_QE, Zero,R_EX_Err_QE);
	exQE->SetMarkerStyle(24);
	exQE->SetMarkerSize(1.5);
	exQE->SetMarkerColor(1);
	exQE->SetLineColor(1);

	TGraphErrors *ex25 = new TGraphErrors(N_25,xbj_25, R_EX_25, Zero,R_EX_Err_25);
	ex25->SetMarkerStyle(24);
	ex25->SetMarkerSize(1.5);
	ex25->SetMarkerColor(1);
	ex25->SetLineColor(1);

	TGraphErrors *ex28 = new TGraphErrors(N_28,xbj_28, R_EX_28, Zero,R_EX_Err_28);
	ex28->SetMarkerStyle(27);
	ex28->SetMarkerSize(1.9);
	ex28->SetMarkerColor(6);
	ex28->SetLineColor(6);

	//if(!(Target0=="C12"&&Target=="He4")){
	if(Target=="He3"){
		//TGraphErrors *exNF = new TGraphErrors(N_NF,xbj_NF, R_EX_NF, Zero,R_EX_Err_NF);
		TGraphAsymmErrors *exNF = new TGraphAsymmErrors(N_NF,xbj_NF, R_EX_NF, Zero,Zero,R_EX_Err_Down,R_EX_Err_Up);
		exNF->SetMarkerStyle(30);
		exNF->SetMarkerSize(1.8);
		exNF->SetMarkerColor(8);
		exNF->SetLineColor(8);
		c3->cd(); exNF->Draw("P");
		c4->cd(); exNF->Draw("P");
		//l1->AddEntry(exNF,Form("E02019 Data (%2.0f^{0})",Theta_NF),"p");
		l1->AddEntry(exNF,Form("E02019 Data (18^{0},Q^{2}=2.7 GeV/c^{2})",Theta_NF),"p");
	}

	//CLAS Data
    cerr<<"---- Plotting CLAS Data"<<endl;
	if(Target0=="He4"&&Target=="He3"){
		TGraphErrors *exCS = new TGraphErrors(N_CS,xbj_CS, R_EX_CS, Zero,R_EX_Err_CS);
		exCS->SetMarkerStyle(20);
		exCS->SetMarkerSize(1.40);
		exCS->SetMarkerColor(7);
		exCS->SetLineColor(7);
		c3->cd(); exCS->Draw("P");
		c4->cd(); exCS->Draw("P");
		l1->AddEntry(exCS,Form("CLAS Data"),"p");
	}

	Theta = 28.0;
	if(Target0=="C12" && Target=="He3"){
		c3->cd(); 	ex28->Draw("P");
		c4->cd(); 	ex28->Draw("P");
		//ex28->Fit("pol0","","",1.7, 2.0);
		l1->AddEntry(ex28, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_28_D,Q2_28_U),"p");
	}

	Theta = 25.0;
	c3->cd(); ex25->Draw("P");
	c4->cd(); ex25->Draw("P");
	l1->AddEntry(ex25, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_25_D,Q2_25_U),"p");
//	ex25->Fit("pol0","","",1.5, 1.8);
	if(!(Target=="H2")){
		c3->cd(); exQE->Draw("P");
		c4->cd(); exQE->Draw("P");
	}
	Theta = 23.0;
	if(!(Target=="H2")){
		c3->cd();	ex23->Draw("P");
		c4->cd();	ex23->Draw("P");
		l1->AddEntry(ex23, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_23_D,Q2_23_U),"p");
	//	ex23->Fit("pol0","","",1.7, 2.0);
	}
	Theta = 21.0;
	c3->cd();	ex21->Draw("P");
	c4->cd();	ex21->Draw("P");
	l1->AddEntry(ex21, Form("%2.0f^{o}: %4.3f<Q^{2}<%4.3f (GeV/c^{2})",Theta,Q2_21_D,Q2_21_U),"p");
//	ex21->Fit("pol0","","",1.75, 2.0);

    gStyle->SetOptFit(0);
    int N_ALL = N_21+N_23+N_25+N_28;
	double R_EX_ALL[N_ALL], R_EX_Err_ALL[N_ALL], xbj_ALL[N_ALL];
	int count = 0;
	for(int i=0;i<N_21;i++){
         R_EX_ALL[count] = R_EX_21[i];
         R_EX_Err_ALL[count] = R_EX_Err_21[i];
         xbj_ALL[count] = xbj_21[i];
		 count++;
	}
	for(int i=0;i<N_23;i++){
         R_EX_ALL[count] = R_EX_23[i];
         R_EX_Err_ALL[count] = R_EX_Err_23[i];
         xbj_ALL[count] = xbj_23[i];
		 count++;
	}
	for(int i=0;i<N_25;i++){
         R_EX_ALL[count] = R_EX_25[i];
         R_EX_Err_ALL[count] = R_EX_Err_25[i];
         xbj_ALL[count] = xbj_25[i];
		 count++;
	}
	for(int i=0;i<N_28;i++){
         R_EX_ALL[count] = R_EX_28[i];
         R_EX_Err_ALL[count] = R_EX_Err_28[i];
         xbj_ALL[count] = xbj_28[i];
		 count++;
	}
	
	if(Target=="Ca40" &&Target0=="Ca48"){
        TLine *bb = new TLine(0.7,1.0,XMax,1.0);
		bb->SetLineStyle(9);
		bb->SetLineColor(7);
		bb->SetLineWidth(2);
        TLine *b2 = new TLine(0.7,0.92,XMax,0.92);
		b2->SetLineStyle(9);
		b2->SetLineColor(8);
		b2->SetLineWidth(2);
		c3->cd(); bb->Draw();b2->Draw();
		c4->cd(); bb->Draw();b2->Draw();
	}

   	TCanvas *c5 = new TCanvas("c5","c5",1200,700);
	c5->cd();
	TH2F *h5 = new TH2F("h5","",300,0.7,XMax,300,RMIN,RMAX);
	h5->SetStats(kFALSE);
	h5->SetXTitle("x");
	h5->GetXaxis()->CenterTitle(1);
	h5->GetXaxis()->SetTitleFont(32);
	h5->GetXaxis()->SetTitleSize(0.08);
	h5->GetXaxis()->SetTitleOffset(0.7);
	h5->GetYaxis()->SetTitle(Form("#frac{#sigma_{%s}}{#sigma_{%s}} #times #frac{%d}{%d}",Target0.Data(),Target.Data(),A2,A1));
	h5->GetYaxis()->CenterTitle(1);
	h5->GetYaxis()->SetTitleFont(32);
	h5->GetYaxis()->SetTitleSize(0.07);
	h5->GetYaxis()->SetTitleOffset(0.6);
	h5->GetYaxis()->CenterTitle();
	gPad->SetLogy(0);
	h5->Draw();

	TGraphErrors *exall = new TGraphErrors(N_ALL,xbj_ALL, R_EX_ALL, Zero,R_EX_Err_ALL);
	exall->SetMarkerStyle(24);
	exall->SetMarkerSize(1.3);
	exall->SetMarkerColor(2);
	exall->SetLineColor(2);
	c5->cd();	exall->Draw("P");
	exall->Fit("pol0","","",1.5, 1.8);

	if(bMC){	
		TGraphErrors *mc21 = new TGraphErrors(N_21,xbj_21, R_MC_21, Zero, Zero);
		mc21->SetLineWidth(2);
		mc21->SetLineColor(2);
		c3->cd(); mc21->Draw("l");
		c4->cd(); mc21->Draw("l");
		l1->AddEntry(mc21, Form("R_{MC} at %d^{0}",21),"l");

		TGraphErrors *mc23 = new TGraphErrors(N_23,xbj_23, R_MC_23, Zero, Zero);
		mc23->SetLineWidth(2);
		mc23->SetLineColor(4);
		c3->cd();	mc23->Draw("l");
		c4->cd();	mc23->Draw("l");
		l1->AddEntry(mc23, Form("R_{MC} at %d^{0}",23),"l");

		TGraphErrors *mcQE = new TGraphErrors(N_QE,xbj_QE, R_MC_QE, Zero, Zero);
		mcQE->SetLineWidth(2);
		mcQE->SetLineColor(6);
		c3->cd();	mcQE->Draw("l");
		c4->cd();	mcQE->Draw("l");
		l1->AddEntry(mcQE, Form("R_{MC} at %d^{0} (QE) ",25),"l");

		TGraphErrors *mc25 = new TGraphErrors(N_25,xbj_25, R_MC_25, Zero, Zero);
		mc25->SetLineWidth(2);
		mc25->SetLineColor(6);
		c3->cd();	mc25->Draw("l");
		c4->cd();	mc25->Draw("l");
		l1->AddEntry(mc25, Form("R_{MC} at %d^{0}",25),"l");

		TGraphErrors *mc28 = new TGraphErrors(N_28,xbj_28, R_MC_28, Zero, Zero);
		mc28->SetLineWidth(2);
		mc28->SetLineColor(3);
		if(Target=="C12"&&Target=="He3"){
			c3->cd();	mc28->Draw("l");
			c4->cd();	mc28->Draw("l");
			l1->AddEntry(mc28, Form("R_{MC} at %d^{0}",28),"l");
		}
	}

	c3->cd(); l1->Draw();
	c4->cd(); l1->Draw();
	if(bMC){
		c3->Print(Form("./left/%s_%s_XS_Ratio_MC_%s_%s_L.png",Target0.Data(),Target.Data(),Com0.Data(),Com.Data()));
		c3->Print(Form("./left/%s_%s_XS_Ratio_MC_%s_%s_L.pdf",Target0.Data(),Target.Data(),Com0.Data(),Com.Data()));
	}
	else{
		c3->Print(Form("./left/%s_%s_XS_Ratio_%s_%s_L.png",Target0.Data(),Target.Data(),Com0.Data(),Com.Data()));
		c3->Print(Form("./left/%s_%s_XS_Ratio_%s_%s_L.pdf",Target0.Data(),Target.Data(),Com0.Data(),Com.Data()));
	}
	if(bMC){
		c4->Print(Form("./left/%s_%s_XS_Ratio_MC_%s_%s_Zoom_L.png",Target0.Data(),Target.Data(),Com0.Data(),Com.Data()));
		c4->Print(Form("./left/%s_%s_XS_Ratio_MC_%s_%s_Zoom_L.pdf",Target0.Data(),Target.Data(),Com0.Data(),Com.Data()));
	}
	else{
		c4->Print(Form("./left/%s_%s_XS_Ratio_%s_%s_Zoom_L.png",Target0.Data(),Target.Data(),Com0.Data(),Com.Data()));
		c4->Print(Form("./left/%s_%s_XS_Ratio_%s_%s_Zoom_L.pdf",Target0.Data(),Target.Data(),Com0.Data(),Com.Data()));
	}

	/*}}}*/
	return 0;
}
/*}}}*/

/*gGet_Xbj_Cut{{{*/
void gGet_Xbj_Cut(const TString &aKinName, const TString & aArm, double* aXbj_Min, double* aXbj_Max){

	if(aKinName.Contains("Kin3.1")){
		if(aArm=="L"){
			aXbj_Min[0] = Xbj_Cut_21_1[0];		
			aXbj_Max[0] = Xbj_Cut_21_1[1];		
		}
		else if(aArm=="R"){
			aXbj_Min[0] = Xbj_Cut_21_1[2];		
			aXbj_Max[0] = Xbj_Cut_21_1[3];		
		}
	}
	else if(aKinName.Contains("Kin3.2")){
		if(aArm=="L"){
			aXbj_Min[0] = Xbj_Cut_21_2[0];		
			aXbj_Max[0] = Xbj_Cut_21_2[1];		
		}
		else if(aArm=="R"){
			aXbj_Min[0] = Xbj_Cut_21_2[2];		
			aXbj_Max[0] = Xbj_Cut_21_2[3];		
		}
	}

	else if(aKinName.Contains("Kin4.1")){
		if(aArm=="L"){
			aXbj_Min[0] = Xbj_Cut_23_1[0];		
			aXbj_Max[0] = Xbj_Cut_23_1[1];		
		}
		else if(aArm=="R"){
			aXbj_Min[0] = Xbj_Cut_23_1[2];		
			aXbj_Max[0] = Xbj_Cut_23_1[3];		
		}
	}
	else if(aKinName.Contains("Kin4.2")){
		if(aArm=="L"){
			aXbj_Min[0] = Xbj_Cut_23_2[0];		
			aXbj_Max[0] = Xbj_Cut_23_2[1];		
		}
		else if(aArm=="R"){
			aXbj_Min[0] = Xbj_Cut_23_2[2];		
			aXbj_Max[0] = Xbj_Cut_23_2[3];		
		}
	}

	else if(aKinName.Contains("Kin5.0") && ! aKinName.Contains("Kin5.05")){
		if(aArm=="L"){
			aXbj_Min[0] = Xbj_Cut_QE_1[0];		
			aXbj_Max[0] = Xbj_Cut_QE_1[1];		
		}
		else if(aArm=="R"){
			aXbj_Min[0] = Xbj_Cut_QE_1[2];		
			aXbj_Max[0] = Xbj_Cut_QE_1[3];		
		}
	}
	else if(aKinName.Contains("Kin5.05")){
		if(aArm=="L"){
			aXbj_Min[0] = Xbj_Cut_QE_2[0];		
			aXbj_Max[0] = Xbj_Cut_QE_2[1];		
		}
		else if(aArm=="R"){
			aXbj_Min[0] = Xbj_Cut_QE_2[2];		
			aXbj_Max[0] = Xbj_Cut_QE_2[3];		
		}
	}

	else if(aKinName.Contains("Kin5.1")){
		if(aArm=="L"){
			aXbj_Min[0] = Xbj_Cut_25_1[0];		
			aXbj_Max[0] = Xbj_Cut_25_1[1];		
		}
		else if(aArm=="R"){
			aXbj_Min[0] = Xbj_Cut_25_1[2];		
			aXbj_Max[0] = Xbj_Cut_25_1[3];		
		}
	}
	else if(aKinName.Contains("Kin5.2")){
		if(aArm=="L"){
			aXbj_Min[0] = Xbj_Cut_25_2[0];		
			aXbj_Max[0] = Xbj_Cut_25_2[1];		
		}
		else if(aArm=="R"){
			aXbj_Min[0] = Xbj_Cut_25_2[2];		
			aXbj_Max[0] = Xbj_Cut_25_2[3];		
		}
	}

	else if(aKinName.Contains("Kin6.5")){
		if(aArm=="L"){
			aXbj_Min[0] = Xbj_Cut_28_1[0];		
			aXbj_Max[0] = Xbj_Cut_28_1[1];		
		}
		else if(aArm=="R"){
			aXbj_Min[0] = Xbj_Cut_28_1[2];		
			aXbj_Max[0] = Xbj_Cut_28_1[3];		
		}
	}

	else{
		aXbj_Min[0]=-1.0;
		aXbj_Max[0]=10.0;
	}

}
/*}}}*/

/*Merge{{{*/
void Merge(const TString& aInputFileName, double* aXbj, double* aR_EX, double* aR_EX_Err, double* aR_MC, int* aN){

	/*Pre-define{{{*/	
	int aline = 0;
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/";

	double scale = 1.0;
	//First Target:
	const TString aName_T1L1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T1L2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T1R1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T1R2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	//Second Target:
	const TString aName_T2L1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T2L2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T2R1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T2R2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const double factor = gGet_InputFile_Var(aInputFileName, aline).Atof(); 

	double Theta = 21.00;
	TString Target = "C12";
	TString Kin = "3.1";
	double Xbj_Min = 1000.0, Xbj_Max = -1000.0;
	Double_t nu,Q2,xbj,angle,Ep,Livetime,Ntg,Ni,Nf_EX,Nf_SAMC,Nf_SAMC_Gen,BinCor;
	Double_t XS,XS_Stat_Err,XS_Sys_Err,XS_MC,Yield_MC,Yield_EX,Y_MC_Sys_Err,Y_EX_Sys_Err,Y_MC_Stat_Err,Y_EX_Stat_Err;
	Double_t Nf_AL, Yield_AL, Y_AL_Stat_Err, Y_AL_Sys_Err;
	TString comment;
	Double_t aRadCorr = 0.0;
	/*}}}*/

	/*First Target{{{*/
	cerr<<endl<<"------------------- First Target ----------------------------------------------------------"<<endl<<endl;

	/*Define{{{*/
	Double_t XS_T1L1[200], XS_Err_T1L1[200], Ep_T1L1[200], Xbj_T1L1[200],  Xbj_Err_T1L1[200],XS_MC_T1L1[200];
	Double_t XS_T1L2[200], XS_Err_T1L2[200], Ep_T1L2[200], Xbj_T1L2[200],  Xbj_Err_T1L2[200],XS_MC_T1L2[200];
	Double_t XS_T1R1[200], XS_Err_T1R1[200], Ep_T1R1[200], Xbj_T1R1[200],  Xbj_Err_T1R1[200],XS_MC_T1R1[200];
	Double_t XS_T1R2[200], XS_Err_T1R2[200], Ep_T1R2[200], Xbj_T1R2[200],  Xbj_Err_T1R2[200],XS_MC_T1R2[200];
	/*}}}*/

	/*First{{{*/
	Target = gGet_Target_Name(aName_T1L1.Data());

	TString infile_T1L1 = aInputFilePath +Target+"/"+Com0+"/"+ aName_T1L1 +"_"+ Com0;
	cerr<<"==== Reading 1st file from "<<infile_T1L1<<endl;
	ifstream run_T1L1; run_T1L1.open(infile_T1L1);

	Int_t Count_T1L1 =0;
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;

	Kin = gGet_Kin_Name(infile_T1L1.Data());

	if(Kin=="MM"){
		cerr<<"1st: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
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

	TString Kin0 = Kin;
	TString Arm_T1L1="X";
	if(aName_T1L1.Contains("L_"))
		Arm_T1L1="L";
	else if(aName_T1L1.Contains("R_"))
		Arm_T1L1="R";

	TString KinName_T1L1 = Form("Kin%s",Kin.Data());
	gGet_Xbj_Cut(KinName_T1L1.Data(), Arm_T1L1.Data(),&Xbj_Min, &Xbj_Max);
	cerr<<Form("@@@  For %s, Xbj_Min = %f, Xbj_Max = %f", KinName_T1L1.Data(), Xbj_Min,Xbj_Max)<<endl;

	Get_XS* xs_T1L1 = new Get_XS();
	xs_T1L1->Load_Table(Target.Data(),KinName_T1L1.Data());

	comment.ReadLine(run_T1L1);
	while ( run_T1L1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
	{
		if ( xbj>=Xbj_Min && xbj<=Xbj_Max && angle<180 &&xbj>0.5&&xbj<=3.00 && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
		{
			Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
			Ep_T1L1[Count_T1L1] = Ep; 
			XS_T1L1[Count_T1L1] = XS * scale; 
			XS_Err_T1L1[Count_T1L1] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
			//if(Target.Contains("Ca"))
			if(1)
				XS_Err_T1L1[Count_T1L1] = pow(XS_Stat_Err,1);
			Xbj_T1L1[Count_T1L1] = xbj;    Xbj_Err_T1L1[Count_T1L1] = 0.0;
			XS_MC_T1L1[Count_T1L1] = xs_T1L1->gGet_XS_Born(E0,Ep,Theta);

			cerr<<Form("%s at %s: Ep = %5.3f, xbj=%5.3f, XS_EX=%e, XS_MC=%e",Target.Data(),
				   	KinName_T1L1.Data(), Ep_T1L1[Count_T1L1], Xbj_T1L1[Count_T1L1], 
					XS_T1L1[Count_T1L1],XS_MC_T1L1[Count_T1L1])<<endl;

			Count_T1L1++;
		}
		if ( !run_T1L1.good() )
			break;
	}
	run_T1L1.close();
	delete xs_T1L1;
	cerr << "  Total XS Points for 1st = " << Count_T1L1 <<endl;
	/*}}}*/

	/*Second{{{*/
	TString KinName_T1L2 ="";
	Bool_t bMore_T1L2 = kTRUE;
	TString infile_T1L2 = aInputFilePath +Target+"/"+Com0+"/"+  aName_T1L2 +"_"+ Com0;
	cerr<<"==== Reading 2nd file from "<<infile_T1L2<<endl;
	TString Arm_T1L2="X";
	if(aName_T1L2 =="0")
		bMore_T1L2 = kFALSE;
	if(aName_T1L2.Contains("L_"))
		Arm_T1L2="L";
	else if(aName_T1L2.Contains("R_"))
		bMore_T1L2 = kFALSE; //Arm_T1L2="R";

	Int_t Count_T1L2=0;
	if(bMore_T1L2){
		ifstream run_T1L2; run_T1L2.open(infile_T1L2);
		Target = gGet_Target_Name(infile_T1L2.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

	Kin = gGet_Kin_Name(infile_T1L2.Data());
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
		KinName_T1L2 = Form("Kin%s",Kin.Data());
		gGet_Xbj_Cut(KinName_T1L2.Data(), Arm_T1L2.Data(),&Xbj_Min, &Xbj_Max);
		cerr<<Form("@@@  For %s, Xbj_Min = %f, Xbj_Max = %f", KinName_T1L2.Data(), Xbj_Min,Xbj_Max)<<endl;

		Get_XS* xs_T1L2 = new Get_XS();
		xs_T1L2->Load_Table(Target.Data(),KinName_T1L2.Data());

		comment.ReadLine(run_T1L2);
		while ( run_T1L2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( xbj>=Xbj_Min && xbj<=Xbj_Max && angle<180 &&xbj>0.5 &&xbj<=3.00  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				Ep_T1L2[Count_T1L2] = Ep; 
				XS_T1L2[Count_T1L2] = XS * scale;    
				XS_Err_T1L2[Count_T1L2] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
			//	if(Target.Contains("Ca"))
				if(1)
						XS_Err_T1L2[Count_T1L2] = pow(XS_Stat_Err,2);
				Xbj_T1L2[Count_T1L2] = xbj;    Xbj_Err_T1L2[Count_T1L2] = 0.0;
				XS_MC_T1L2[Count_T1L2] = xs_T1L2->gGet_XS_Born(E0,Ep,Theta);

				cerr<<Form("%s at %s: Ep = %5.3f, xbj=%5.3f, XS_EX=%e, XS_MC=%e",Target.Data(),
				   	KinName_T1L2.Data(), Ep_T1L2[Count_T1L2], Xbj_T1L2[Count_T1L2], 
					XS_T1L2[Count_T1L2],XS_MC_T1L2[Count_T1L2])<<endl;

				Count_T1L2++;
			}
			if ( !run_T1L2.good() )
				break;
		}
		run_T1L2.close();
		delete xs_T1L2;
		cerr << "  Total XS Points for 2nd  = " << Count_T1L2 <<endl;
	}
	/*}}}*/

	/*Third{{{*/
	TString KinName_T1R1 ="";
	Bool_t bMore_T1R1 = kTRUE;
	TString infile_T1R1 =  aInputFilePath +Target+"/"+Com0+"/"+ aName_T1R1 +"_"+ Com0;
	cerr<<"==== Reading 3rd file from "<<infile_T1R1<<endl;
	TString Arm_T1R1="X";
	if(aName_T1R1 =="0")
		bMore_T1R1 = kFALSE;
	if(aName_T1R1.Contains("L_"))
		Arm_T1R1="L";
	else if(aName_T1R1.Contains("R_"))
		bMore_T1R1 = kFALSE; //Arm_T1R1="R";

	Int_t Count_T1R1=0;
	if(bMore_T1R1){
		ifstream run_T1R1; run_T1R1.open(infile_T1R1);
		Target = gGet_Target_Name(infile_T1R1.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

	Kin = gGet_Kin_Name(infile_T1R1.Data());
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
		KinName_T1R1 = Form("Kin%s",Kin.Data());
		gGet_Xbj_Cut(KinName_T1R1.Data(), Arm_T1R1.Data(),&Xbj_Min, &Xbj_Max);
		cerr<<Form("@@@  For %s, Xbj_Min = %f, Xbj_Max = %f", KinName_T1R1.Data(), Xbj_Min,Xbj_Max)<<endl;

		Get_XS* xs_T1R1 = new Get_XS();
		xs_T1R1->Load_Table(Target.Data(),KinName_T1R1.Data());

		comment.ReadLine(run_T1R1);
		while ( run_T1R1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( xbj>=Xbj_Min && xbj<=Xbj_Max && angle<180 &&xbj>0.5 &&xbj<=3.00  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				Ep_T1R1[Count_T1R1] = Ep; 
				XS_T1R1[Count_T1R1] = XS * scale; 
				XS_Err_T1R1[Count_T1R1] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
			//	if(Target.Contains("Ca"))
				if(1)
					XS_Err_T1R1[Count_T1R1] = pow(XS_Stat_Err,1);
				Xbj_T1R1[Count_T1R1] = xbj;    Xbj_Err_T1R1[Count_T1R1] = 0.0;
				XS_MC_T1R1[Count_T1R1] = xs_T1R1->gGet_XS_Born(E0,Ep,Theta);

				cerr<<Form("%s at %s: Ep = %5.3f, xbj=%5.3f, XS_EX=%e, XS_MC=%e",Target.Data(),
				   	KinName_T1R1.Data(), Ep_T1R1[Count_T1R1], Xbj_T1R1[Count_T1R1], 
					XS_T1R1[Count_T1R1],XS_MC_T1R1[Count_T1R1])<<endl;

				Count_T1R1++;
			}
			if ( !run_T1R1.good() )
				break;
		}
		run_T1R1.close();
		delete xs_T1R1;
		cerr << "  Total XS Points for 2nd  = " << Count_T1R1 <<endl;
	}
	/*}}}*/

	/*Forth{{{*/
	TString KinName_T1R2 ="";
	Bool_t bMore_T1R2 = kTRUE;
	TString infile_T1R2 = aInputFilePath +Target+"/"+Com0+"/"+  aName_T1R2 +"_"+ Com0;
	cerr<<"==== Reading 4th file from "<<infile_T1R2<<endl;
	TString Arm_T1R2="X";
	if(aName_T1R2 =="0")
		bMore_T1R2 = kFALSE;
	if(aName_T1R2.Contains("L_"))
		Arm_T1R2="L";
	else if(aName_T1R2.Contains("R_"))
		bMore_T1R2 = kFALSE; //Arm_T1R2="R";

	Int_t Count_T1R2=0;
	if(bMore_T1R2){
		ifstream run_T1R2; run_T1R2.open(infile_T1R2);
		Target = gGet_Target_Name(infile_T1R2.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

	Kin = gGet_Kin_Name(infile_T1R2.Data());
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
		KinName_T1R2 = Form("Kin%s",Kin.Data());
		gGet_Xbj_Cut(KinName_T1R2.Data(), Arm_T1R2.Data(),&Xbj_Min, &Xbj_Max);
		cerr<<Form("@@@  For %s, Xbj_Min = %f, Xbj_Max = %f", KinName_T1R2.Data(), Xbj_Min,Xbj_Max)<<endl;

		Get_XS* xs_T1R2 = new Get_XS();
		xs_T1R2->Load_Table(Target.Data(),KinName_T1R2.Data());

		comment.ReadLine(run_T1R2);
		while ( run_T1R2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( xbj>=Xbj_Min && xbj<=Xbj_Max && angle<180 &&xbj>0.5  &&xbj<=3.00 && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				Ep_T1R2[Count_T1R2] = Ep; 
				XS_T1R2[Count_T1R2] = XS * scale;
				XS_Err_T1R2[Count_T1R2] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
			//	if(Target.Contains("Ca"))
				if(1)
					XS_Err_T1R2[Count_T1R2] = pow(XS_Stat_Err,1);
				Xbj_T1R2[Count_T1R2] = xbj;    Xbj_Err_T1R2[Count_T1R2] = 0.0;
				XS_MC_T1R2[Count_T1R2] = xs_T1R2->gGet_XS_Born(E0,Ep,Theta);

				cerr<<Form("%s at %s: Ep = %5.3f, xbj=%5.3f, XS_EX=%e, XS_MC=%e",Target.Data(),
				   	KinName_T1R2.Data(), Ep_T1R2[Count_T1R2], Xbj_T1R2[Count_T1R2], 
					XS_T1R2[Count_T1R2],XS_MC_T1R2[Count_T1R2])<<endl;

				Count_T1R2++;
			}
			if ( !run_T1R2.good() )
				break;
		}
		run_T1R2.close();
		delete xs_T1R2;
		cerr << "  Total XS Points for 2nd  = " << Count_T1R2 <<endl;
	}
	/*}}}*/

	/*Merge{{{*/
	//Add First Two Input
	Double_t Ep_T1L[200],Xbj_T1L[200], XS_T1L[200], XS_Err_T1L[200], XS_MC_T1L[200];
	Double_t XS_Sum_T1L = 0.0, XS_Err_Sum_T1L = 0.0; Int_t fSum_T1L = 0;
	Int_t Count_T1L = 0, J_T1L = 0;
	for(int i=0;i<Count_T1L1;i++){
		XS_Sum_T1L = XS_T1L1[i];
		XS_Err_Sum_T1L = pow(XS_Err_T1L1[i]/XS_T1L1[i],2);
		fSum_T1L = 1;
		if(bMore_T1L2){
			for(int j=0;j<Count_T1L2;j++){
				//Find the verlapping bin for two setting, such as Kin3.1 and Kin3.2
				if(abs(Xbj_T1L1[i]-Xbj_T1L2[j])< 0.0001){
				//if(abs(Ep_T1L1[i]-Ep_T1L2[j])< 0.0001){
					XS_Sum_T1L += XS_T1L2[j];
					XS_Err_Sum_T1L +=	pow(XS_Err_T1L2[j]/XS_T1L2[j],2);
					fSum_T1L++; J_T1L = j;
					break; //No two xbj bin within one setting
				}
			}
		}
		Ep_T1L[i] = Ep_T1L1[i];
		Xbj_T1L[i] = Xbj_T1L1[i];
		XS_MC_T1L[i] = XS_MC_T1L1[i];
		XS_T1L[i] = XS_Sum_T1L / fSum_T1L;
		XS_Err_T1L[i] = XS_T1L[i] * sqrt(XS_Err_Sum_T1L) / fSum_T1L; 
		Count_T1L++;
	}
	for(int k=J_T1L+1;k<Count_T1L2;k++){
		Ep_T1L[Count_T1L] = Ep_T1L2[k];
		Xbj_T1L[Count_T1L] = Xbj_T1L2[k];
		XS_T1L[Count_T1L] = XS_T1L2[k];
		XS_Err_T1L[Count_T1L] = XS_Err_T1L2[k];
		XS_MC_T1L[Count_T1L] = XS_MC_T1L2[k];
		Count_T1L++;	
	}
	//Add Last Two Input
	Double_t Ep_T1R[200],Xbj_T1R[200], XS_T1R[200], XS_Err_T1R[200], XS_MC_T1R[200];
	Double_t XS_Sum_T1R = 0.0, XS_Err_Sum_T1R = 0.0; Int_t fSum_T1R = 0;
	Int_t Count_T1R = 0, J_T1R = 0;
	if(bMore_T1R1){
		for(int i=0;i<Count_T1R1;i++){
			XS_Sum_T1R = XS_T1R1[i];
			XS_Err_Sum_T1R = pow(XS_Err_T1R1[i]/XS_T1R1[i],2);
			fSum_T1R = 1;
			if(bMore_T1R2){
				for(int j=0;j<Count_T1R2;j++){
					//Find the verlapping bin for two setting, such as Kin3.1 and Kin3.2
					if(abs(Xbj_T1R1[i]-Xbj_T1R2[j])< 0.0001){
					//if(abs(Ep_T1R1[i]-Ep_T1R2[j])< 0.0001){
						XS_Sum_T1R += XS_T1R2[j];
						XS_Err_Sum_T1R +=	pow(XS_Err_T1R2[j]/XS_T1R2[j],2);
						fSum_T1R++; J_T1R = j;
						break; //No two xbj bin within one setting
					}
				}
			}
			Ep_T1R[i] = Ep_T1R1[i];
			Xbj_T1R[i] = Xbj_T1R1[i];
			XS_MC_T1R[i] = XS_MC_T1R1[i];
			XS_T1R[i] = XS_Sum_T1R / fSum_T1R;
			XS_Err_T1R[i] = XS_T1R[i] * sqrt(XS_Err_Sum_T1R) / fSum_T1R; 
			Count_T1R++;
		}
		for(int k=J_T1R+1;k<Count_T1R2;k++){
			Ep_T1R[Count_T1R] = Ep_T1R2[k];
			Xbj_T1R[Count_T1R] = Xbj_T1R2[k];
			XS_T1R[Count_T1R] = XS_T1R2[k];
			XS_Err_T1R[Count_T1R] = XS_Err_T1R2[k];
			XS_MC_T1L[Count_T1R] = XS_MC_T1R2[k];
			Count_T1R++;	
		}
	}

	vector<double> vEp1,vXbj1, vXbj_Err1,vXS1, vXS_Err1, vXS_MC1;
	double XS_MC_T1, XS_T1, XS_Err_T1;
	double X_Min1 = 100.0, X_Max1 = -100.0, Y_Min1 = 10000.0, Y_Max1 = -10000.0;
	int N1 = 0, J1R= 0;
	bool bFound_T1L = kFALSE;
	if(bMore_T1R1){
		for(int i=0;i<Count_T1L; i++){
			bFound_T1L = kFALSE;	
			for(int j=0;j<Count_T1R; j++){
				if( abs(Xbj_T1L[i] - Xbj_T1R[j])<0.0001 ){
				//if( abs(Ep_T1L[i] - Ep_T1R[j])<0.0001 ){
					vEp1.push_back(Ep_T1L[i]);	
					vXbj1.push_back(Xbj_T1L[i]);	
					bFound_T1L = kTRUE;	
					vXbj_Err1.push_back(0.0);	
					XS_T1 = (XS_T1L[i]+XS_T1R[j])/2.0;
					XS_Err_T1 = XS_T1 * sqrt( pow(XS_Err_T1L[i]/XS_T1L[i],2) + pow(XS_Err_T1R[j]/XS_T1R[j],2) )/2.0;
					//XS_MC_T1 = 	(XS_MC_T1L[i]+XS_MC_T1R[j])/2.0;
					XS_MC_T1 = 	XS_MC_T1L[i];
					vXS1.push_back(XS_T1);
					vXS_Err1.push_back(XS_Err_T1);
					vXS_MC1.push_back( XS_MC_T1 );
					N1++; J1R = j;
					if(Xbj_T1L[i] < X_Min1)
						X_Min1 = Xbj_T1L[i];
					if(Xbj_T1L[i] > X_Max1)
						X_Max1 = Xbj_T1L[i];
					if(XS_T1 < Y_Min1)
						Y_Min1 = XS_T1;
					if(XS_T1 > Y_Max1)
						Y_Max1 = XS_T1;
					if(XS_MC_T1 < Y_Min1)
						Y_Min1 = XS_MC_T1;
					if(XS_T1 > Y_Max1)
						Y_Max1 = XS_MC_T1;
				}
			}
			if(!bFound_T1L){
				vEp1.push_back(Ep_T1L[i]);	
				vXbj1.push_back(Xbj_T1L[i]);	
				vXbj_Err1.push_back(0.0);	
				XS_T1 = XS_T1L[i];
				XS_Err_T1 = XS_Err_T1L[i];
				XS_MC_T1 = 	XS_MC_T1L[i];
				vXS1.push_back(XS_T1);
				vXS_Err1.push_back(XS_Err_T1);
				vXS_MC1.push_back( XS_MC_T1 );
				N1++;
				if(Xbj_T1L[i] < X_Min1)
					X_Min1 = Xbj_T1L[i];
				if(Xbj_T1L[i] > X_Max1)
					X_Max1 = Xbj_T1L[i];
				if(XS_T1 < Y_Min1)
					Y_Min1 = XS_T1;
				if(XS_T1 > Y_Max1)
					Y_Max1 = XS_T1;
				if(XS_MC_T1 < Y_Min1)
					Y_Min1 = XS_MC_T1;
				if(XS_T1 > Y_Max1)
					Y_Max1 = XS_MC_T1;
			}
		}

		for(int k =J1R+1; k<Count_T1R; k++){
			vEp1.push_back(Ep_T1R[k]);	
			vXbj1.push_back(Xbj_T1R[k]);	
			vXbj_Err1.push_back(0.0);	
			XS_T1 = XS_T1R[k];
			XS_Err_T1 = XS_Err_T1R[k];
			XS_MC_T1 = 	XS_MC_T1R[k];
			vXS1.push_back(XS_T1);
			vXS_Err1.push_back(XS_Err_T1);
			vXS_MC1.push_back( XS_MC_T1 );
			N1++;
			if(Xbj_T1R[k] < X_Min1)
				X_Min1 = Xbj_T1R[k];
			if(Xbj_T1R[k] > X_Max1)
				X_Max1 = Xbj_T1R[k];
			if(XS_T1 < Y_Min1)
				Y_Min1 = XS_T1;
			if(XS_T1 > Y_Max1)
				Y_Max1 = XS_T1;
			if(XS_MC_T1 < Y_Min1)
				Y_Min1 = XS_MC_T1;
			if(XS_T1 > Y_Max1)
				Y_Max1 = XS_MC_T1;
		}
	}
	else{
		for(int i=0;i<Count_T1L; i++){
			vEp1.push_back(Ep_T1L[i]);	
			vXbj1.push_back(Xbj_T1L[i]);	
			vXbj_Err1.push_back(0.0);	
			XS_T1 = XS_T1L[i];
			XS_Err_T1 = XS_Err_T1L[i];
			XS_MC_T1 = 	XS_MC_T1L[i];
			vXS1.push_back(XS_T1);
			vXS_Err1.push_back(XS_Err_T1);
			vXS_MC1.push_back( XS_MC_T1 );
			N1++;
			if(Xbj_T1L[i] < X_Min1)
				X_Min1 = Xbj_T1L[i];
			if(Xbj_T1L[i] > X_Max1)
				X_Max1 = Xbj_T1L[i];
			if(XS_T1 < Y_Min1)
				Y_Min1 = XS_T1;
			if(XS_T1 > Y_Max1)
				Y_Max1 = XS_T1;
			if(XS_MC_T1 < Y_Min1)
				Y_Min1 = XS_MC_T1;
			if(XS_T1 > Y_Max1)
				Y_Max1 = XS_MC_T1;
		}
	}
  
	/*RadCorr for He3 at x>2.45{{{*/
	if(Target=="He3"){
		for(int i=0;i<vXbj1.size(); i++){
			if(vXbj1[i]>=2.45){	
				aRadCorr = gGet_RadCorr_He3(vXbj1[i]);
				vXS1[i]*=aRadCorr;
				vXS_Err1[i]*=aRadCorr;
			}
		}
	}
	/*}}}*/

	/*}}}*/

	/*Plot{{{*/
	TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	gPad->SetLogy(1);
	TH2F* xs1 = new TH2F("xs1","",1000,0.7,3.2,1000,Y_Min1/2.,Y_Max1*2.0);
	xs1->SetStats(kFALSE);
	xs1->SetXTitle("x");
	xs1->GetXaxis()->CenterTitle(1);
	xs1->GetXaxis()->SetTitleFont(32);
	xs1->GetXaxis()->SetTitleSize(0.06);
	xs1->SetYTitle("#frac{d#sigma}{dE'd#Omega} (nb/MeV/sr)");
	xs1->GetYaxis()->CenterTitle(1);
	xs1->GetYaxis()->SetTitleFont(32);
	xs1->GetYaxis()->SetTitleSize(0.04);
	xs1->Draw(); 

	TGraphErrors *pl_x1 = new TGraphErrors(N1,&vXbj1[0],&vXS1[0],&vXbj_Err1[0],&vXS_Err1[0]);
	pl_x1-> SetMarkerColor(2);
	pl_x1-> SetMarkerSize(1.0);
	pl_x1-> SetMarkerStyle(22);
	pl_x1-> GetXaxis()->SetTitle("x");
	pl_x1-> GetYaxis()->SetTitle("#sigma");
	pl_x1->Draw("P");

	TGraphErrors *pl_x11 = new TGraphErrors(N1,&vXbj1[0],&vXS_MC1[0],&vXbj_Err1[0],&vXbj_Err1[0]);
	pl_x11-> SetMarkerColor(0);
	pl_x11-> SetMarkerSize(0.0);
	pl_x11-> SetLineColor(4);
	pl_x11-> SetLineWidth(2);
	pl_x11-> SetMarkerStyle(20);
	pl_x11-> GetXaxis()->SetTitle("x");
	pl_x11-> GetYaxis()->SetTitle("#sigma_{mc}");
	pl_x11->Draw("L");

	TLegend *l1 = new TLegend(0.65,0.60,0.90, 0.90,Form("%s Target",Target.Data()));
	l1->SetTextSize(0.03);
	l1->SetTextFont(22);
	l1->AddEntry(pl_x1,  Form("XS_{EX} at %3.1f^{0}",Theta),"p");
	pl_x1->Draw("P");
	l1->AddEntry(pl_x11, Form("XS_{MC} at %3.1f^{0}",Theta),"l");
	pl_x11->Draw("L");
	l1->Draw();

/*	if(aInputFileName.Contains("QE")){
		c1->Print(Form("XS_%s_%2.0f_QE_L_%s.png", Target.Data(), Theta,Com0.Data()));
		c1->Print(Form("XS_%s_%2.0f_QE_L_%s.pdf", Target.Data(), Theta,Com0.Data()));
	}
	else{
		c1->Print(Form("XS_%s_%2.0f_%s_L.png", Target.Data(), Theta,Com0.Data()));
		c1->Print(Form("XS_%s_%2.0f_%s_L.pdf", Target.Data(), Theta,Com0.Data()));
	}*/
	/*}}}*/
	/*}}}*/

	/*Second Target{{{*/
	cerr<<endl<<"------------------ Second Target ----------------------------------------------------------"<<endl<<endl;

	/*Define{{{*/
	Double_t XS_T2L1[200], XS_Err_T2L1[200], Ep_T2L1[200], Xbj_T2L1[200],  Xbj_Err_T2L1[200],XS_MC_T2L1[200];
	Double_t XS_T2L2[200], XS_Err_T2L2[200], Ep_T2L2[200], Xbj_T2L2[200],  Xbj_Err_T2L2[200],XS_MC_T2L2[200];
	Double_t XS_T2R1[200], XS_Err_T2R1[200], Ep_T2R1[200], Xbj_T2R1[200],  Xbj_Err_T2R1[200],XS_MC_T2R1[200];
	Double_t XS_T2R2[200], XS_Err_T2R2[200], Ep_T2R2[200], Xbj_T2R2[200],  Xbj_Err_T2R2[200],XS_MC_T2R2[200];
	/*}}}*/

	/*First{{{*/
	Target = gGet_Target_Name(aName_T2L1.Data());
	
	TString infile_T2L1 = aInputFilePath +Target+"/"+Com+"/"+  aName_T2L1 +"_"+ Com;
	cerr<<"==== Reading 1st file from "<<infile_T2L1<<endl;
	ifstream run_T2L1; run_T2L1.open(infile_T2L1);

	Int_t Count_T2L1 =0;
	if(Target=="He4")
		scale = He4_Norm;
	else
		scale = 1.0;

	Kin = gGet_Kin_Name(infile_T2L1.Data());

	TString Arm_T2L1="X";
	if(aName_T2L1.Contains("L_"))
		Arm_T2L1="L";
	else if(aName_T2L1.Contains("R_"))
		Arm_T2L1="R";


	if(Kin=="MM"){
		cerr<<"1st: Kin (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5) = "; cin >> Kin;
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

	TString KinName_T2L1 = Form("Kin%s",Kin.Data());
	gGet_Xbj_Cut(KinName_T2L1.Data(), Arm_T2L1.Data(),&Xbj_Min, &Xbj_Max);
	cerr<<Form("@@@  For %s, Xbj_Min = %f, Xbj_Max = %f", KinName_T2L1.Data(), Xbj_Min,Xbj_Max)<<endl;

	Get_XS* xs_T2L1 = new Get_XS();
	xs_T2L1->Load_Table(Target.Data(),KinName_T2L1.Data());

	comment.ReadLine(run_T2L1);
	while ( run_T2L1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
	{
		if ( xbj>=Xbj_Min && xbj<=Xbj_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
		{
			Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
			Ep_T2L1[Count_T2L1]=Ep;
			XS_T2L1[Count_T2L1] = XS * scale;
			XS_Err_T2L1[Count_T2L1] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
			//if(Target.Contains("Ca"))
			if(1)
				XS_Err_T2L1[Count_T2L1] = pow(XS_Stat_Err,1);
			Xbj_T2L1[Count_T2L1] = xbj;    Xbj_Err_T2L1[Count_T2L1] = 0.0;
			XS_MC_T2L1[Count_T2L1] = xs_T2L1->gGet_XS_Born(E0,Ep,Theta);

			cerr<<Form("%s at %s: Ep = %5.3f, xbj=%5.3f, XS_EX=%e, XS_MC=%e",Target.Data(),
					KinName_T2L1.Data(), Ep_T2L1[Count_T2L1], Xbj_T2L1[Count_T2L1], 
					XS_T2L1[Count_T2L1],XS_MC_T2L1[Count_T2L1])<<endl;

			Count_T2L1++;
		}
		if ( !run_T2L1.good() )
			break;
	}
	run_T2L1.close();
	delete xs_T2L1;
	cerr << "  Total XS Points for 1st = " << Count_T2L1 <<endl;
	/*}}}*/

	/*Second{{{*/
	TString KinName_T2L2 ="";
	Bool_t bMore_T2L2 = kTRUE;
	TString infile_T2L2 = aInputFilePath +Target+"/"+Com+"/"+  aName_T2L2 +"_"+ Com;
	cerr<<"==== Reading 2nd file from "<<infile_T2L2<<endl;
	TString Arm_T2L2="X";
	if(aName_T2L2 =="0")
		bMore_T2L2 = kFALSE;
	if(aName_T2L2.Contains("L_"))
		Arm_T2L2="L";
	else if(aName_T2L2.Contains("R_"))
		bMore_T2L2 = kFALSE; //Arm_T2L2="R";

	Int_t Count_T2L2=0;
	if(bMore_T2L2){
		ifstream run_T2L2; run_T2L2.open(infile_T2L2);
		Target = gGet_Target_Name(infile_T2L2.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

		Kin = gGet_Kin_Name(infile_T2L2.Data());
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
		KinName_T2L2 = Form("Kin%s",Kin.Data());
		gGet_Xbj_Cut(KinName_T2L2.Data(), Arm_T2L2.Data(),&Xbj_Min, &Xbj_Max);
		cerr<<Form("@@@  For %s, Xbj_Min = %f, Xbj_Max = %f", KinName_T2L2.Data(), Xbj_Min,Xbj_Max)<<endl;

		Get_XS* xs_T2L2 = new Get_XS();
		xs_T2L2->Load_Table(Target.Data(),KinName_T2L2.Data());

		comment.ReadLine(run_T2L2);
		while ( run_T2L2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( xbj>=Xbj_Min && xbj<=Xbj_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				Ep_T2L2[Count_T2L2]=Ep;
				XS_T2L2[Count_T2L2] = XS * scale; 
				XS_Err_T2L2[Count_T2L2] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
			//	if(Target.Contains("Ca"))
				if(1)
					XS_Err_T2L2[Count_T2L2] = pow(XS_Stat_Err,1);
				Xbj_T2L2[Count_T2L2] = xbj;    Xbj_Err_T2L2[Count_T2L2] = 0.0;
				XS_MC_T2L2[Count_T2L2] = xs_T2L2->gGet_XS_Born(E0,Ep,Theta);

				cerr<<Form("%s at %s: Ep = %5.3f, xbj=%5.3f, XS_EX=%e, XS_MC=%e",Target.Data(),
						KinName_T2L2.Data(), Ep_T2L2[Count_T2L2], Xbj_T2L2[Count_T2L2], 
						XS_T2L2[Count_T2L2],XS_MC_T2L2[Count_T2L2])<<endl;

				Count_T2L2++;
			}
			if ( !run_T2L2.good() )
				break;
		}
		run_T2L2.close();
		delete xs_T2L2;
		cerr << "  Total XS Points for 2nd  = " << Count_T2L2 <<endl;
	}
	/*}}}*/

	/*Third{{{*/
	TString KinName_T2R1 ="";
	Bool_t bMore_T2R1 = kTRUE;
	TString infile_T2R1 = aInputFilePath +Target+"/"+Com+"/"+  aName_T2R1 +"_"+ Com;
	cerr<<"==== Reading 3rd file from "<<infile_T2R1<<endl;
	TString Arm_T2R1="X";
	if(aName_T2R1 =="0")
		bMore_T2R1 = kFALSE;
	if(aName_T2R1.Contains("L_"))
		Arm_T2R1="L";
	else if(aName_T2R1.Contains("R_"))
		bMore_T2R1 = kFALSE; //Arm_T2R1="R";

	Int_t Count_T2R1=0;
	if(bMore_T2R1){
		ifstream run_T2R1; run_T2R1.open(infile_T2R1);
		Target = gGet_Target_Name(infile_T2R1.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

		Kin = gGet_Kin_Name(infile_T2R1.Data());
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
		KinName_T2R1 = Form("Kin%s",Kin.Data());
		gGet_Xbj_Cut(KinName_T2R1.Data(), Arm_T2R1.Data(),&Xbj_Min, &Xbj_Max);
		cerr<<Form("@@@  For %s, Xbj_Min = %f, Xbj_Max = %f", KinName_T2R1.Data(), Xbj_Min,Xbj_Max)<<endl;

		Get_XS* xs_T2R1 = new Get_XS();
		xs_T2R1->Load_Table(Target.Data(),KinName_T2R1.Data());

		comment.ReadLine(run_T2R1);
		while ( run_T2R1>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( xbj>=Xbj_Min && xbj<=Xbj_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				Ep_T2R1[Count_T2R1]=Ep;
				XS_T2R1[Count_T2R1] = XS* scale;   
				XS_Err_T2R1[Count_T2R1] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2));
			//	if(Target.Contains("Ca"))
				if(1)
					XS_Err_T2R1[Count_T2R1] = pow(XS_Stat_Err,1);
				Xbj_T2R1[Count_T2R1] = xbj;    Xbj_Err_T2R1[Count_T2R1] = 0.0;
				XS_MC_T2R1[Count_T2R1] = xs_T2R1->gGet_XS_Born(E0,Ep,Theta);

				cerr<<Form("%s at %s: Ep = %5.3f, xbj=%5.3f, XS_EX=%e, XS_MC=%e",Target.Data(),
						KinName_T2R1.Data(), Ep_T2R1[Count_T2R1], Xbj_T2R1[Count_T2R1], 
						XS_T2R1[Count_T2R1],XS_MC_T2R1[Count_T2R1])<<endl;

				Count_T2R1++;
			}
			if ( !run_T2R1.good() )
				break;
		}
		run_T2R1.close();
		delete xs_T2R1;
		cerr << "  Total XS Points for 2nd  = " << Count_T2R1 <<endl;
	}
	/*}}}*/

	/*Forth{{{*/
	TString KinName_T2R2 ="";
	Bool_t bMore_T2R2 = kTRUE;
	TString infile_T2R2 = aInputFilePath +Target+"/"+Com+"/"+ aName_T2R2 +"_"+ Com;
	cerr<<"==== Reading 4th file from "<<infile_T2R2<<endl;
	TString Arm_T2R2="X";
	if(aName_T2R2 =="0")
		bMore_T2R2 = kFALSE;
	if(aName_T2R2.Contains("L_"))
		Arm_T2R2="L";
	else if(aName_T2R2.Contains("R_"))
		bMore_T2R2 = kFALSE;  //Arm_T2R2="R";

	Int_t Count_T2R2=0;
	if(bMore_T2R2){
		ifstream run_T2R2; run_T2R2.open(infile_T2R2);
		Target = gGet_Target_Name(infile_T2R2.Data());
		if(Target=="He4")
		scale = He4_Norm;
		else
			scale = 1.0;

		Kin = gGet_Kin_Name(infile_T2R2.Data());
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
		KinName_T2R2 = Form("Kin%s",Kin.Data());
		gGet_Xbj_Cut(KinName_T2R2.Data(), Arm_T2R2.Data(),&Xbj_Min, &Xbj_Max);
		cerr<<Form("@@@  For %s, Xbj_Min = %f, Xbj_Max = %f", KinName_T2R2.Data(), Xbj_Min,Xbj_Max)<<endl;

		Get_XS* xs_T2R2 = new Get_XS();
		xs_T2R2->Load_Table(Target.Data(),KinName_T2R2.Data());

		comment.ReadLine(run_T2R2);
		while ( run_T2R2>>XS>>XS_Stat_Err>>XS_Sys_Err>>XS_MC>>xbj>>angle>>Ep>>nu>>Q2>>BinCor>>Livetime>>Ntg>>Ni>>Nf_EX>>Nf_AL>>Nf_SAMC>>Nf_SAMC_Gen>>Yield_EX>>Y_EX_Stat_Err>>Y_EX_Sys_Err>>Yield_MC>>Y_MC_Stat_Err>>Y_MC_Sys_Err>>Yield_AL>>Y_AL_Stat_Err>>Y_AL_Sys_Err )
		{
			if ( xbj>=Xbj_Min && xbj<=Xbj_Max && angle<180 &&xbj>0.5  && XS>0.0 && sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,2))/XS <0.2)
			{
				Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
				Ep_T2R2[Count_T2R2]=Ep;
				XS_T2R2[Count_T2R2] = XS * scale;     
				XS_Err_T2R2[Count_T2R2] = sqrt(pow(XS_Sys_Err,2)+pow(XS_Stat_Err,1));
			//	if(Target.Contains("Ca"))
				if(1)
					XS_Err_T2R2[Count_T2R2] = pow(XS_Stat_Err,1);
				Xbj_T2R2[Count_T2R2] = xbj;    Xbj_Err_T2R2[Count_T2R2] = 0.0;
				XS_MC_T2R2[Count_T2R2] = xs_T2R2->gGet_XS_Born(E0,Ep,Theta);

				cerr<<Form("%s at %s: Ep = %5.3f, xbj=%5.3f, XS_EX=%e, XS_MC=%e",Target.Data(),
						KinName_T2R2.Data(), Ep_T2R2[Count_T2R2], Xbj_T2R2[Count_T2R2], 
						XS_T2R2[Count_T2R2],XS_MC_T2R2[Count_T2R2])<<endl;

				Count_T2R2++;
			}
			if ( !run_T2R2.good() )
				break;
		}
		run_T2R2.close();
		delete xs_T2R2;
		cerr << "  Total XS Points for 2nd  = " << Count_T2R2 <<endl;
	}
	/*}}}*/

	/*Merge{{{*/
	//Add First Two Input
	Double_t Ep_T2L[200],Xbj_T2L[200], XS_T2L[200], XS_Err_T2L[200], XS_MC_T2L[200];
	Double_t XS_Sum_T2L = 0.0, XS_Err_Sum_T2L = 0.0; Int_t fSum_T2L = 0;
	Int_t Count_T2L = 0, J_T2L = 0;
	for(int i=0;i<Count_T2L1;i++){
		XS_Sum_T2L = XS_T2L1[i];
		XS_Err_Sum_T2L = pow(XS_Err_T2L1[i]/XS_T2L1[i],2);
		fSum_T2L = 1;
		if(bMore_T2L2){
			for(int j=0;j<Count_T2L2;j++){
				//Find the verlapping bin for two setting, such as Kin3.1 and Kin3.2
				if(abs(Xbj_T2L1[i]-Xbj_T2L2[j])< 0.0001){
				//if(abs(Ep_T2L1[i]-Ep_T2L2[j])< 0.0001){
					XS_Sum_T2L += XS_T2L2[j];
					XS_Err_Sum_T2L +=	pow(XS_Err_T2L2[j]/XS_T2L2[j],2);
					fSum_T2L++; J_T2L = j;
					break; //No two xbj bin within one setting
				}
			}
		}
		Ep_T2L[i] = Ep_T2L1[i];
		Xbj_T2L[i] = Xbj_T2L1[i];
		XS_MC_T2L[i] = XS_MC_T2L1[i];
		XS_T2L[i] = XS_Sum_T2L / fSum_T2L;
		XS_Err_T2L[i] = XS_T2L[i] * sqrt(XS_Err_Sum_T2L) / fSum_T2L; 
		Count_T2L++;
	}
	for(int k=J_T2L+1;k<Count_T2L2;k++){
		Ep_T2L[Count_T2L] = Ep_T2L2[k];
		Xbj_T2L[Count_T2L] = Xbj_T2L2[k];
		XS_T2L[Count_T2L] = XS_T2L2[k];
		XS_Err_T2L[Count_T2L] = XS_Err_T2L2[k];
		XS_MC_T2L[Count_T2L] = XS_MC_T2L2[k];
		Count_T2L++;	
	}

	//Add Last Two Input
	Double_t Ep_T2R[200], Xbj_T2R[200], XS_T2R[200], XS_Err_T2R[200], XS_MC_T2R[200];
	Double_t XS_Sum_T2R = 0.0, XS_Err_Sum_T2R = 0.0; Int_t fSum_T2R = 0;
	Int_t Count_T2R = 0, J_T2R = 0;
	if(bMore_T2R1){
		for(int i=0;i<Count_T2R1;i++){
			XS_Sum_T2R = XS_T2R1[i];
			XS_Err_Sum_T2R = pow(XS_Err_T2R1[i]/XS_T2R1[i],2);
			fSum_T2R = 1;
			if(bMore_T2R2){
				for(int j=0;j<Count_T2R2;j++){
					//Find the verlapping bin for two setting, such as Kin3.1 and Kin3.2
					if(abs(Xbj_T2R1[i]-Xbj_T2R2[j])< 0.0001){
					//if(abs(Ep_T2R1[i]-Ep_T2R2[j])< 0.0001){
						XS_Sum_T2R += XS_T2R2[j];
						XS_Err_Sum_T2R +=	pow(XS_Err_T2R2[j]/XS_T2R2[j],2);
						fSum_T2R++; J_T2R = j;
						break; //No two xbj bin within one setting
					}
				}
			}
			Ep_T2R[i] = Ep_T2R1[i];
			Xbj_T2R[i] = Xbj_T2R1[i];
			XS_MC_T2R[i] = XS_MC_T2R1[i];
			XS_T2R[i] = XS_Sum_T2R / fSum_T2R;
			XS_Err_T2R[i] = XS_T2R[i] * sqrt(XS_Err_Sum_T2R) / fSum_T2R; 
			Count_T2R++;
		}
		for(int k=J_T2R+1;k<Count_T2R2;k++){
			Ep_T2R[Count_T2R] = Ep_T2R2[k];
			Xbj_T2R[Count_T2R] = Xbj_T2R2[k];
			XS_T2R[Count_T2R] = XS_T2R2[k];
			XS_MC_T2R[Count_T2R] = XS_MC_T2R2[k];
			XS_Err_T2R[Count_T2R] = XS_Err_T2R2[k];
			Count_T2R++;	
		}
	}

	vector<double> vEp2, vXbj2, vXbj_Err2, vXS2, vXS_Err2, vXS_MC2;
	double XS_MC_T2, XS_T2,XS_Err_T2;
	double X_Min2 = 100.0, X_Max2 = -100.0, Y_Min2 = 10000.0, Y_Max2 = -10000.0;
	int N2 = 0, J2R = 0;
	bool bFound_T2L = kFALSE;
	if(bMore_T2R1){
		for(int i=0;i<Count_T2L; i++){
			bFound_T2L = kFALSE;
			for(int j=0;j<Count_T2R; j++){
				if( abs(Xbj_T2L[i] - Xbj_T2R[j])<0.0001 ){
				//if( abs(Ep_T2L[i] - Ep_T2R[j])<0.0001 ){
					bFound_T2L = kTRUE;
					vEp2.push_back(Ep_T2L[i]);	
					vXbj2.push_back(Xbj_T2L[i]);	
					vXbj_Err2.push_back(0.0);	
					XS_T2 = (XS_T2L[i]+XS_T2R[j])/2.0;
					XS_Err_T2 = XS_T2 * sqrt( pow(XS_Err_T2L[i]/XS_T2L[i],2) + pow(XS_Err_T2R[j]/XS_T2R[j],2) )/2.0;
					XS_MC_T2 = (XS_MC_T2L[i]+XS_MC_T2R[j])/2.0;
					vXS2.push_back(XS_T2);
					vXS_Err2.push_back(XS_Err_T2);
					vXS_MC2.push_back( XS_MC_T2 );
					N2++; J2R = j;

					if(Xbj_T2L[i] < X_Min2)
						X_Min2 = Xbj_T2L[i];
					if(Xbj_T2L[i] > X_Max2)
						X_Max2 = Xbj_T2L[i];
					if(XS_T2 < Y_Min2)
						Y_Min2 = XS_T2;
					if(XS_T2 > Y_Max2)
						Y_Max2 = XS_T2;
					if(XS_MC_T2 < Y_Min2)
						Y_Min2 = XS_MC_T2;
					if(XS_T2 > Y_Max2)
						Y_Max2 = XS_MC_T2;
				}
			}
			if(!bFound_T2L){
				vEp2.push_back(Ep_T2L[i]);	
				vXbj2.push_back(Xbj_T2L[i]);	
				vXbj_Err2.push_back(0.0);	
				XS_T2 = XS_T2L[i];
				XS_Err_T2 = XS_Err_T2L[i];
				XS_MC_T2 = XS_MC_T2L[i];
				vXS2.push_back(XS_T2);
				vXS_Err2.push_back(XS_Err_T2);
				vXS_MC2.push_back( XS_MC_T2 );
				N2++;

				if(Xbj_T2L[i] < X_Min2)
					X_Min2 = Xbj_T2L[i];
				if(Xbj_T2L[i] > X_Max2)
					X_Max2 = Xbj_T2L[i];
				if(XS_T2 < Y_Min2)
					Y_Min2 = XS_T2;
				if(XS_T2 > Y_Max2)
					Y_Max2 = XS_T2;
				if(XS_MC_T2 < Y_Min2)
					Y_Min2 = XS_MC_T2;
				if(XS_T2 > Y_Max2)
					Y_Max2 = XS_MC_T2;
			}
		}
		for(int k = J2R+1; k<Count_T2R; k++){
			vEp2.push_back(Ep_T2R[k]);	
			vXbj2.push_back(Xbj_T2R[k]);	
			vXbj_Err2.push_back(0.0);	
			XS_T2 = XS_T2R[k];
			XS_Err_T2 = XS_Err_T2R[k];
			XS_MC_T2 = XS_MC_T2R[k];
			vXS2.push_back(XS_T2);
			vXS_Err2.push_back(XS_Err_T2);
			vXS_MC2.push_back( XS_MC_T2 );
			N2++;

			if(Xbj_T2R[k] < X_Min2)
				X_Min2 = Xbj_T2R[k];
			if(Xbj_T2R[k] > X_Max2)
				X_Max2 = Xbj_T2R[k];
			if(XS_T2 < Y_Min2)
				Y_Min2 = XS_T2;
			if(XS_T2 > Y_Max2)
				Y_Max2 = XS_T2;
			if(XS_MC_T2 < Y_Min2)
				Y_Min2 = XS_MC_T2;
			if(XS_T2 > Y_Max2)
				Y_Max2 = XS_MC_T2;
		}
	}
	else{
		for(int i=0;i<Count_T2L; i++){
			vEp2.push_back(Ep_T2L[i]);	
			vXbj2.push_back(Xbj_T2L[i]);	
			vXbj_Err2.push_back(0.0);	
			XS_T2 = XS_T2L[i];
			XS_Err_T2 = XS_Err_T2L[i];
			XS_MC_T2 = XS_MC_T2L[i];
			vXS2.push_back(XS_T2);
			vXS_Err2.push_back(XS_Err_T2);
			vXS_MC2.push_back( XS_MC_T2 );
			N2++;

			if(Xbj_T2L[i] < X_Min2)
				X_Min2 = Xbj_T2L[i];
			if(Xbj_T2L[i] > X_Max2)
				X_Max2 = Xbj_T2L[i];
			if(XS_T2 < Y_Min2)
				Y_Min2 = XS_T2;
			if(XS_T2 > Y_Max2)
				Y_Max2 = XS_T2;
			if(XS_MC_T2 < Y_Min2)
				Y_Min2 = XS_MC_T2;
			if(XS_T2 > Y_Max2)
				Y_Max2 = XS_MC_T2;
		}
	}

	/*RadCorr for He3 at x>2.45{{{*/
	if(Target=="He3"){
		for(int i=0;i<vXbj2.size(); i++){
			if(vXbj2[i]>=2.45){	
				aRadCorr = gGet_RadCorr_He3(vXbj2[i]);
				cerr<<Form("++++++ xbj = %f, radcor = %f", vXbj2[i], aRadCorr)<<endl;
				vXS2[i]*=aRadCorr;
				vXS_Err2[i]*=aRadCorr;
			}
		}
	}
	/*}}}*/
	
	/*}}}*/

	/*Plot{{{*/
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	gPad->SetLogy(1);
	TH2F* xs2 = new TH2F("xs2","",1000,0.7,3.2,1000,Y_Min2/2.,Y_Max2*2.0);
	xs2->SetStats(kFALSE);
	xs2->SetXTitle("x");
	xs2->GetXaxis()->CenterTitle(1);
	xs2->GetXaxis()->SetTitleFont(32);
	xs2->GetXaxis()->SetTitleSize(0.06);
	xs2->SetYTitle("#frac{d#sigma}{dE'd#Omega} (nb/MeV/sr)");
	xs2->GetYaxis()->CenterTitle(1);
	xs2->GetYaxis()->SetTitleFont(32);
	xs2->GetYaxis()->SetTitleSize(0.04);
	xs2->Draw(); 

	TGraphErrors *pl_x2 = new TGraphErrors(N2,&vXbj2[0],&vXS2[0],&vXbj_Err2[0],&vXS_Err2[0]);
	pl_x2-> SetMarkerColor(2);
	pl_x2-> SetMarkerSize(1.0);
	pl_x2-> SetMarkerStyle(22);
	pl_x2-> GetXaxis()->SetTitle("x");
	pl_x2-> GetYaxis()->SetTitle("#sigma");
	pl_x2-> Draw("P");

	TGraphErrors *pl_x21 = new TGraphErrors(N2,&vXbj2[0],&vXS_MC2[0],&vXbj_Err2[0],&vXbj_Err2[0]);
	pl_x21-> SetMarkerColor(0);
	pl_x21-> SetMarkerSize(0.0);
	pl_x21-> SetLineColor(4);
	pl_x21-> SetLineWidth(2);
	pl_x21-> SetMarkerStyle(20);
	pl_x21-> GetXaxis()->SetTitle("x");
	pl_x21-> GetYaxis()->SetTitle("#sigma_{mc}");
	pl_x21-> Draw("L");

	TLegend *l2 = new TLegend(0.65,0.60,0.90, 0.90,Form("%s Target",Target.Data()));
	l2->SetTextSize(0.03);
	l2->SetTextFont(22);
	l2->AddEntry(pl_x2,  Form("XS_{EX} at %3.1f^{0}",Theta),"p");
	pl_x2->Draw("P");
	l2->AddEntry(pl_x21, Form("XS_{MC} at %3.1f^{0}",Theta),"l");
	pl_x21->Draw("L");
	l2->Draw();

/*	if(aInputFileName.Contains("QE")){
		c2->Print(Form("XS_%s_%2.0f_QE_L_%s.png", Target.Data(), Theta,Com.Data()));
		c2->Print(Form("XS_%s_%2.0f_QE_L_%s.pdf", Target.Data(), Theta,Com.Data()));
	}
	else{
		c2->Print(Form("XS_%s_%2.0f_%s_L.png", Target.Data(), Theta,Com.Data()));
		c2->Print(Form("XS_%s_%2.0f_%s_L.pdf", Target.Data(), Theta,Com.Data()));
	}*/
	/*}}}*/
	/*}}}*/

	/*Ratio{{{*/  
	ofstream outfile(Form("./left/XS_%s_%s_%s_%s_%s.dat",Kin0.Data(),Target0.Data(),Target.Data(),Com0.Data(),Com.Data()));
	outfile << Form("%7s %14s %14s %14s %14s %14s %14s", "Xbj","Ratio","XS1","XS2","Ratio_Err","XS_Err1","XS_Err2" )<<endl;

	Int_t Count=0;
	Double_t x1,x2,ep1,ep2,ratio,ratio_mc,err, max;
	vector<Double_t> vEp, vXbj, vXbj_Err,vR, vR_Err,vR_MC;
    
	Double_t X_Min = 10, X_Max = -10, R_Min = 1000.0, R_Max = -1000.0;
	int i,j; max=-1.0;
	for( i=0;i<N1;i++){
		x1 = vXbj1[i];
		ep1 = vEp1[i];
		for( j=0;j<N2;j++){
			x2 = vXbj2[j];
			ep2 = vEp2[j];
			if(x1>=0.9&&x2>=0.9&&abs(x1-x2)<0.0001){
			//if(x1>=0.9&&x2>=0.9&&abs(ep1-ep2)<0.0001){
				ratio = factor*vXS1[i]/vXS2[j];
				err = ratio * sqrt( pow(vXS_Err2[j]/vXS2[j],2) + pow(vXS_Err1[i]/vXS1[i],2) );
				ratio_mc = factor*vXS_MC1[i]/vXS_MC2[j];
				if(ratio>1e-30 && ratio<1000.&& (err/ratio)<1.0){
					vEp.push_back(ep1);
					vXbj.push_back(x1);
					vXbj_Err.push_back(0.0);
					vR.push_back(ratio);
					vR_MC.push_back(ratio_mc);
					vR_Err.push_back(err);
					cerr <<Form("-----> Found one: xbj1=%f, xbj2=%f, XS1=%f, XS2=%f,Ratio=%f ",x1,x2,vXS1[i],vXS2[j],factor*vXS1[i]/vXS2[j])<<endl;

					if(R_Min > ratio)
						R_Min = ratio;
					if(R_Min > ratio_mc)
						R_Min = ratio_mc;
					if(R_Max < ratio)
						R_Max = ratio;
					//				if(R_Max < ratio_mc)
					//					R_Max = ratio_mc;
					if(X_Min > x1)
						X_Min = x1;
					if(X_Max < x1)
						X_Max = x1;
					Count++;
					outfile << Form("%7.5f %14e %14e %14e %14e %14e %14e", x1,ratio, vXS1[i], vXS2[j],err, vXS_Err1[i], vXS_Err2[j] )<<endl;
				}
			}
		}
	}

	//Return the value
	for(int k=0; k<Count; k++){
		aXbj[k] = vXbj[k];
		aR_EX[k] = vR[k];
		aR_EX_Err[k] = vR_Err[k];
		aR_MC[k] = vR_MC[k];
	}
	aN[0] = Count;
	cerr<<"******************** Factor = "<< factor<<endl;
	outfile.close();
	/*}}}*/
}
/*}}}*/

/*E02-019 Data{{{*/
void E02019(const TString& aInputFileName, const double aTheta,double* aXbj, double* aR_EX, double* aR_EX_Err, int* aN){
	int aline = 0;
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/Ratio/";

	ofstream outfile(Form("./left/E02019_XS_%s_%s.dat",Target0.Data(),Target.Data()));
	outfile << Form("%7s %14s %14s %14s %14s %14s %14s", "Xbj","Ratio","XS1","XS2","Ratio_Err","XS_Err1","XS_Err2" )<<endl;

	const TString aName_T1 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const TString aName_T2 = gGet_InputFile_Var(aInputFileName, aline); aline++;
	const double factor = gGet_InputFile_Var(aInputFileName, aline).Atof(); 
	cerr<<Form("==== Reading E02-019 Data from %s and %s ", aName_T1.Data(), aName_T2.Data())<<endl;

	TString Ascii_File1 = Form("%s/e02019/%s",aInputFilePath.Data(),aName_T1.Data());
	cerr<<Form(" === Reading E02019 data from %s", Ascii_File1.Data())<<endl;
	ifstream infile1(Ascii_File1);
	TString Ascii_File2 = Form("%s/e02019/%s",aInputFilePath.Data(),aName_T2.Data());
	cerr<<Form(" === Reading E02019 data from %s", Ascii_File2.Data())<<endl;
	ifstream infile2(Ascii_File2);

	Double_t xbj, Q2, nA, nTheta, nP0, nEp, nXS_Rad, nXS_Born,nXS_CC,nFrac_Err;  
	const double nE0 = 5.776;
	/*First Target:{{{*/
	vector<double> vXbj1, vXS1, vXS_Err1,vP1,vTheta1;
	int N1=0;
	while(!(infile1.eof())){
		infile1 >> nA >> nTheta >> nP0 >> nEp >> nXS_Rad >> nXS_Born >> nXS_CC >> nFrac_Err;
		if(abs(nTheta-aTheta)<1.0|| (aTheta-0.0)<1){
		//if(abs(nTheta-aTheta)<1.0|| (nTheta)<aTheta){
			Q2 = 4.0 * nE0 * nEp * pow(sin(nTheta*3.1415926/180./2.0),2);
			xbj = Q2/2.0/0.938272/(nE0-nEp);
			//	cerr<<Form("=== E0=%f, Ep=%f, Theta=%f",nE0,nEp,nTheta)<<endl;
			N1++;
			vXbj1.push_back(xbj);
			vXS1.push_back(nXS_Born);
			vXS_Err1.push_back(nXS_Born*nFrac_Err);
			vTheta1.push_back(nTheta);
			vP1.push_back(nP0);
		}
	}
	/*}}}*/

	/* Second  Target:{{{*/
	vector<double> vXbj2, vXS2, vXS_Err2,vP2,vTheta2;
	int N2=0;
	while(!(infile2.eof())){
		infile2 >> nA >> nTheta >> nP0 >> nEp >> nXS_Rad >> nXS_Born >> nXS_CC >> nFrac_Err;

		if(abs(nTheta-aTheta)<1.0||(aTheta-0.0)<1){
		//if(abs(nTheta-aTheta)<1.0|| (nTheta)<aTheta){
			Q2 = 4.0 * nE0 * nEp * pow(sin(nTheta*3.1415926/180./2.0),2);
			xbj = Q2/2.0/0.938272/(nE0-nEp);
			//	cerr<<Form("===  E0=%f, Ep=%f, Theta=%f",nE0,nEp,nTheta)<<endl;
			N2++;
			vXbj2.push_back(xbj);
			vXS2.push_back(nXS_Born);
			vXS_Err2.push_back(nXS_Born*nFrac_Err);
			vTheta2.push_back(nTheta);
			vP2.push_back(nP0);
		}
	}
	/*}}}*/

	/*Ratio{{{*/  
	Int_t Count=0;
	Double_t x1,x2,ratio,err;
	vector<Double_t> vXbj, vXbj_Err,vR, vR_Err;

	for(int i=0;i<N1;i++){
		x1 = vXbj1[i]; 
		for(int j=0;j<N2;j++){
			x2 = vXbj2[j];
			if(abs(x1-x2)<0.0001 && abs(vP1[i]-vP2[j])<0.1&& abs(vTheta1[i]-vTheta2[j])<1.0){
				ratio = factor*vXS1[i]/vXS2[j];
				err = ratio * sqrt( pow(vXS_Err2[j]/vXS2[j],2) + pow(vXS_Err1[i]/vXS1[i],2) );
				if(ratio>1e-30 && ratio<100.&&(err/ratio)<0.5){
					vXbj.push_back(x1);
					vXbj_Err.push_back(0.0);
					vR.push_back(ratio);
					vR_Err.push_back(err);
					cerr <<Form("-----> Found one: xbj1=%f, xbj2=%f, XS1=%f, XS2=%f,Ratio=%f ",x1,x2,vXS1[i],vXS2[j],factor*vXS1[i]/vXS2[j])<<endl;
					Count++;
					outfile << Form("%7.5f %14e %14e %14e %14e %14e %14e", x1,ratio, vXS1[i], vXS2[j],err, vXS_Err1[i], vXS_Err2[j] )<<endl;
				}
			}
		}
	}

	//Return the value
	for(int k=0; k<Count; k++){
		aXbj[aN[0]+k] = vXbj[k];
		aR_EX[aN[0]+k] = vR[k];
		aR_EX_Err[aN[0]+k] = vR_Err[k];
	}
	aN[0] += Count;
	cerr<<Form("******************** Count = %d, Factor =%f ",Count, factor)<<endl;
	outfile.close();
	/*}}}*/

}
/*}}}*/

/*gGetTargetAZ{{{*/
void gGetTargetAZ(const TString& aTarget, int* aA, int* aZ){
   if(aTarget=="H2"||aTarget=="D2"){
     aA[0] = 2; aZ[0] = 1;
   }
   else if(aTarget=="He3"){
     aA[0] = 3; aZ[0] = 2;
   }
   else if(aTarget=="He4"){
     aA[0] = 4; aZ[0] = 2;
   }
   else if(aTarget=="C12"){
     aA[0] = 12; aZ[0] = 6;
   }
   else if(aTarget=="Ca40"){
     aA[0] = 40; aZ[0] = 20;
   }
   else if(aTarget=="Ca48"){
     aA[0] = 48; aZ[0] = 20;
   }
   else{
	   cerr<<" **** I don't know this target:  "<<aTarget.Data()<<endl;
	   aA[0]=0; aZ[0]=0;
   } 
}
/*}}}*/

/*E02-019 A/He3{{{*/
void E02019_A_He3(const TString& aTarget, const double aTheta,double* aXbj, double* aR_EX, double* aR_EX_Err_Up,double* aR_EX_Err_Down, int* aN){
	int aline = 0;
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/Ratio/";
	TString aInputFile = aInputFilePath + "Nadia_A_He3.dat";
	ifstream outfile(aInputFile);
    Double_t nTheta,nXbj,nRatio,nErr1,nErr2;
	Double_t aErr_Up,aErr_Down;
    int aA,aZ,nA;
    gGetTargetAZ(aTarget.Data(), &aA,&aZ);
    
	int aNNN = 0;
	while(!outfile.eof()){
		outfile >> nA >> nTheta >> nXbj >> nRatio >> nErr1 >> nErr2;
		if(abs(nA-aA)<1.0&&(abs(aTheta-nTheta)<1.0 ||(aTheta-0.0)<1.0)&&nRatio>=1e-33){
			aXbj[aNNN] = nXbj;
			aR_EX[aNNN] = nRatio;
			aErr_Up = abs(nErr1 - nRatio);
			aErr_Down = abs(nRatio - nErr2);	
			aR_EX_Err_Up[aNNN] = aErr_Up;
			aR_EX_Err_Down[aNNN] = aErr_Down;
			aNNN++;
		}
	}
	aN[0]=aNNN;
	outfile.close();
}
/*}}}*/

/*E02-019 A/He3{{{*/
void E02019_He4_He3(const TString& aTarget, const double aTheta,double* aXbj, double* aR_EX, double* aR_EX_Err_Up,double* aR_EX_Err_Down, int* aN){
	int aline = 0;
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/Ratio/";
	TString aInputFile = aInputFilePath + "Nadia_A_He3.dat";
	ifstream outfile(aInputFile);
    Double_t nTheta,nXbj,nRatio,nErr,nErr1,nErr2;
	Double_t aErr_Up,aErr_Down;
    int aA,aZ,nA;
    gGetTargetAZ(aTarget.Data(), &aA,&aZ);
    double aXbj_Max = 1.55;
	int aNNN = 0;
	while(!outfile.eof()){
		outfile >> nA >> nTheta >> nXbj >> nRatio >> nErr1 >> nErr2;
		if(abs(nA-aA)<1.0&&(abs(aTheta-nTheta)<1.0 ||(aTheta-0.0)<1.0)&&nRatio>=1e-33&&nXbj<=aXbj_Max){
			aXbj[aNNN] = nXbj;
			aR_EX[aNNN] = nRatio;
			aErr_Up = abs(nErr1 - nRatio);
			aErr_Down = abs(nRatio - nErr2);	
			aR_EX_Err_Up[aNNN] = aErr_Up;
			aR_EX_Err_Down[aNNN] = aErr_Down;
			aNNN++;
		}
	}

	TString aInputFile1 = aInputFilePath + "Nadia_He4_He3.dat";
	ifstream outfile1(aInputFile1);
	while(!outfile1.eof()){
		outfile1 >> nXbj >> nRatio >> nErr;
			aXbj[aNNN] = nXbj;
			aR_EX[aNNN] = 1./nRatio;
			aErr_Up = 1./(nRatio-nErr);
			aErr_Down = 1./(nRatio + nErr);	
			aR_EX_Err_Up[aNNN] = aErr_Up - 1./nRatio;
			aR_EX_Err_Down[aNNN] = 1./nRatio - aErr_Down;
			aNNN++;
	}
	aN[0]=aNNN;
	outfile.close();
	outfile1.close();
}
/*}}}*/

/*CLAS_He4_He3{{{*/
void CLAS_He4_He3( double* aXbj_CS, double* aR_EX_CS, double* aR_EX_Err_CS, int *aN_CS){
	TString aInputFilePath = "/work/halla/e08014/disk1/yez/Ep_Cross_Section/results/Ratio/";
	TString aInputFile = aInputFilePath + "CLAS_He4_He3.dat";
	ifstream infile(aInputFile);

	double cXbj, cR, cR_Err;
	int aN = 0;

	while(!(infile.eof())){
      infile >> cXbj >> cR >> cR_Err;	
      aXbj_CS[aN] = cXbj;
	  aR_EX_CS[aN] = cR/1.15;
	  aR_EX_Err_CS[aN] = cR_Err/1.15;
	  aN++;	
	}
	aN_CS[0] = aN;
	infile.close();
}
/*}}}*/
