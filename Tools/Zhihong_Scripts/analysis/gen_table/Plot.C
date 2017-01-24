#include "./XGT2_Main.h"
#include "./XGT2_XEMC.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
/*Cuts on the edge of Ep{{{*/
// Cut_L_Min Cut_L_Max Cut_R_Min Cut_R_Max
const double E0 = 3.356;
const double Ep_Cut_21_1[4] = {2.80, 3.00,2.80,3.00}; 
const double Ep_Cut_21_2[4] = {2.94, 3.18,2.94,3.18}; 
const double Ep_Cut_23_1[4] = {2.74, 2.94,2.74,2.94}; 
const double Ep_Cut_23_2[4] = {2.92, 3.16,2.92,3.16};
const double Ep_Cut_QE_1[4] = {2.42, 2.58,2.42,2.58}; 
const double Ep_Cut_QE_2[4] = {2.56, 2.74,2.56,2.74}; 
const double Ep_Cut_25_1[4] = {2.70, 2.88,2.70,2.88};
const double Ep_Cut_25_2[4] = {2.88, 3.12,2.88,3.12}; 
const double Ep_Cut_28_1[4] = {2.74, 2.98,2.74,2.98}; 
/*}}}*/

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

/*gGet_Ep_Cut{{{*/
void gGet_Ep_Cut(const TString &aKinName, const TString & aArm, double* aEp_Min, double* aEp_Max, double* aTheta){

	if(aKinName.Contains("Kin3.1")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_21_1[0];		
			aEp_Max[0] = Ep_Cut_21_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_21_1[2];		
			aEp_Max[0] = Ep_Cut_21_1[3];		
		}
		aTheta[0] = 21.00;
	}
	else if(aKinName.Contains("Kin3.2")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_21_2[0];		
			aEp_Max[0] = Ep_Cut_21_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_21_2[2];		
			aEp_Max[0] = Ep_Cut_21_2[3];		
		}
		aTheta[0] = 21.00;
	}

	else if(aKinName.Contains("Kin4.1")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_23_1[0];		
			aEp_Max[0] = Ep_Cut_23_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_23_1[2];		
			aEp_Max[0] = Ep_Cut_23_1[3];		
		}
		aTheta[0] = 23.00;
	}
	else if(aKinName.Contains("Kin4.2")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_23_2[0];		
			aEp_Max[0] = Ep_Cut_23_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_23_2[2];		
			aEp_Max[0] = Ep_Cut_23_2[3];		
		}
		aTheta[0] = 23.00;
	}

	else if(aKinName.Contains("Kin5.0") && ! aKinName.Contains("Kin5.05")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_QE_1[0];		
			aEp_Max[0] = Ep_Cut_QE_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_QE_1[2];		
			aEp_Max[0] = Ep_Cut_QE_1[3];		
		}
		aTheta[0] = 25.00;
	}
	else if(aKinName.Contains("Kin5.05")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_QE_2[0];		
			aEp_Max[0] = Ep_Cut_QE_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_QE_2[2];		
			aEp_Max[0] = Ep_Cut_QE_2[3];		
		}
		aTheta[0] = 25.00;
	}

	else if(aKinName.Contains("Kin5.1")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_25_1[0];		
			aEp_Max[0] = Ep_Cut_25_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_25_1[2];		
			aEp_Max[0] = Ep_Cut_25_1[3];		
		}
		aTheta[0] = 25.00;
	}
	else if(aKinName.Contains("Kin5.2")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_25_2[0];		
			aEp_Max[0] = Ep_Cut_25_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_25_2[2];		
			aEp_Max[0] = Ep_Cut_25_2[3];		
		}
		aTheta[0] = 25.00;
	}

	else if(aKinName.Contains("Kin6.5")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_28_1[0];		
			aEp_Max[0] = Ep_Cut_28_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_28_1[2];		
			aEp_Max[0] = Ep_Cut_28_1[3];		
		}
		aTheta[0] = 28.00;
	}

	else{
		aEp_Min[0]= 1000.0;
		aEp_Max[0]=-1000.0;
		aTheta[0] = 0.00;
	}

}
/*}}}*/


//void Plot(){
int main(){

	double Ep_Min=0.0, Ep_Max=0.0;
	double Step = 0.01;
    const double E0 = 3.356;
	double Ep = 0.0, Theta = 0.0, Q2 = 0, xbj = 0.0;
	double Y_Min=1000000.0, Y_Max=-1000000.0;
	double X_Min=1000000.0, X_Max=-1000000.0;

	TString Target1="He3"; 
  	cerr<<"Target 1 (H2,He3,He4,C12,Ca40,Ca48)= "; cin >> Target1;
	TString Kin1 = "3.1";
   	cerr<<"Kin 1 (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5)= "; cin >> Kin1;
	Kin1 = "Kin"+Kin1;
	Get_XS* xs1 = new Get_XS();
	xs1->Load_Table(Target1.Data(),Kin1.Data());

	TString Target2="He3"; 
  	cerr<<"Target 2 (H2,He3,He4,C12,Ca40,Ca48)= "; cin >> Target2;
	TString Kin2 = "3.1";
   	cerr<<"Kin 2 (3.1,3.2,4.1,4.2,5.0,5.05,5.1,5.2,6.5)= "; cin >> Kin2;
	Kin2 = "Kin"+Kin2;
	Get_XS* xs2 = new Get_XS();
	xs2->Load_Table(Target2.Data(),Kin2.Data());

	/*Read In{{{*/
	double XS1[100], Xbj1[100],XS_Err1[100], Xbj_Err1[100];
	int i=0;
	gGet_Ep_Cut(Kin1.Data(),"L",&Ep_Min,&Ep_Max, &Theta);
	Ep = 0.0, xbj = 0.0;;
	while(xbj<=3.2){
		Ep = Ep_Min + i*Step;
		Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
		xbj = Q2/2.0/0.938272/(E0-Ep);
		Xbj1[i] = xbj;
		XS1[i] = xs1->gGet_XS_Born(E0,Ep,Theta);
        Xbj_Err1[i] = 0.0; XS_Err1[i] = 0.0;
        if(Y_Min > XS1[i])
			Y_Min = XS1[i];
        if(Y_Max < XS1[i])
			Y_Max = XS1[i];
		if(X_Min > Xbj1[i])
			X_Min = Xbj1[i];
        if(X_Max < Xbj1[i])
			X_Max = Xbj1[i];
		i++;
	}
	const int N1 = i;

	double XS2[100], Xbj2[100],XS_Err2[100], Xbj_Err2[100];
	gGet_Ep_Cut(Kin2.Data(),"L",&Ep_Min,&Ep_Max, &Theta);
	int j=0;
	Ep = 0.0, xbj = 0.0;;
	while(xbj<=3.2){
		Ep = Ep_Min + j*Step;
		Q2 = 4.0 * E0 * Ep * pow(sin(Theta*3.1415926/180./2.0),2);
		xbj = Q2/2.0/0.938272/(E0-Ep);
		Xbj2[j] = xbj;
		XS2[j] = xs2->gGet_XS_Born(E0,Ep,Theta);
        Xbj_Err2[j] = 0.0; XS_Err2[j] = 0.0;
		if(Y_Min > XS2[j])
			Y_Min = XS2[j];
        if(Y_Max < XS2[j])
			Y_Max = XS2[j];
		if(X_Min > Xbj2[j])
			X_Min = Xbj2[j];
        if(X_Max < Xbj2[j])
			X_Max = Xbj2[j];
			j++;
	}
	const int N2 = j;

	delete xs1;
	delete xs2;
    /*}}}*/

	/*Draw{{{*/
	TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	gPad->SetLogy(1);
	TH2F* xs = new TH2F("xs1","",1000,X_Min-0.1,X_Max+0.1,1000,Y_Max*1e-6,Y_Max*2.0);
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

	TGraphErrors *pl_x1 = new TGraphErrors(N1,&Xbj1[0],&XS1[0],&Xbj_Err1[0],&XS_Err1[0]);
	pl_x1-> SetMarkerColor(2);
	pl_x1-> SetMarkerSize(1.5);
	pl_x1-> SetMarkerStyle(22);
	pl_x1-> GetXaxis()->SetTitle("x_{bj}");
	pl_x1-> GetYaxis()->SetTitle("#sigma");
	pl_x1->Draw("P");

	TGraphErrors *pl_x2 = new TGraphErrors(N2,&Xbj2[0],&XS2[0],&Xbj_Err2[0],&XS_Err2[0]);
	pl_x2-> SetMarkerColor(4);
	pl_x2-> SetMarkerSize(1.5);
	pl_x2-> SetMarkerStyle(24);
	pl_x2-> GetXaxis()->SetTitle("x_{bj}");
	pl_x2-> GetYaxis()->SetTitle("#sigma");
	pl_x2->Draw("P");

	TLegend *l1 = new TLegend(0.65,0.60,0.90, 0.90,Form("%s Target",Target1.Data()));
	l1->SetTextSize(0.03);
	l1->SetTextFont(22);
	l1->AddEntry(pl_x1,  Form("XS_{EX} at %3.1f^{0}",Theta),"p");
	l1->AddEntry(pl_x2,  Form("XS_{EX} at %3.1f^{0} (Hacked)",Theta),"p");
	l1->Draw();

	c1->Print(Form("%s_%s_XS.pdf",Target1.Data(), Kin1.Data()));
	c1->Print(Form("%s_%s_XS.png",Target1.Data(), Kin1.Data()));
    /*}}}*/
}
