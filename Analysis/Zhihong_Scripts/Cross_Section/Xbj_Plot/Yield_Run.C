
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
#include <TLine.h>
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
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TLegend.h>
#include <TRandom3.h>
//#include <TMatrix.h>
/*}}}*/

using namespace::std;
using namespace::TMath;
const TString XGT2_DIR = "/work/halla/e08014/disk1/yez/Xbj_Cross_Section/";

/*inline vector<TString> Tokens(TString aline,TString aDelim=" "){{{*/
inline vector<TString> Tokens(TString aline,TString aDelim=" ")
{
	Int_t i;
	TObjArray* InObjArray;
	TObjString* os;
	TString s;
	vector<TString> OutStringVec;
	OutStringVec.clear();

	InObjArray=aline.Tokenize(aDelim);
	for ( i=0; i<InObjArray->GetEntriesFast(); i++ )
	{
		os=(TObjString*)InObjArray->At(i);
		s=os->GetString();
		OutStringVec.push_back(s);
	}
	return OutStringVec;
}
/*}}}*/

/*inline vector<Int_t> gGet_RunNoChain(const Int_t& aRunNo,const TString& aString){{{*/
inline vector<Int_t> gGet_RunNoChain(const Int_t& aRunNo,const TString& aString)
{
	vector<Int_t> output;
	if(aRunNo>1000)
		output.push_back(aRunNo);
	vector<TString> input=Tokens(aString);
	unsigned int i;
	for ( i=0; i<input.size(); i++ )
	{
		if ( isdigit(*input[i].Data()) || input[i].BeginsWith("-") )
		{
			vector<TString> tmp=Tokens(input[i],'-');
			if ( tmp.size()<=2 )
			{
				// Int_t FirstRunNo;
				Int_t LastRunNo;
				Int_t thisrunno;
				size_t index=input[i].First('-');
				if ( index==0 )
				{
					thisrunno=aRunNo;
					LastRunNo=tmp[0].Atoi();
				}
				else if ( index<input[i].Length() )
				{
					thisrunno=tmp[0].Atoi();
					LastRunNo=tmp[1].Atoi();
				}
				else
				{
					thisrunno=tmp[0].Atoi();
					LastRunNo=thisrunno;
				}
				if ( thisrunno>LastRunNo )
				{
					thisrunno ^= LastRunNo;
					LastRunNo ^= thisrunno;
					thisrunno ^= LastRunNo;
				}
				while ( thisrunno<=LastRunNo )
				{
					if ( thisrunno!=aRunNo )
						output.push_back(thisrunno);
					thisrunno++;
				}
			}
		}
	}
	for ( i=0; i<output.size(); i++ )
	{
		if ( i==0 && output.size()>1 ){
//			cerr<<"      RunNoChain is "<<output[i];
		}
		else if ( i==0 && i==(output.size()-1) ){
//			cerr<<"      RunNoChain is "<<output[i]<<endl;
		}
		else if ( i==(output.size()-1) ){
//			cerr<<" "<<output[i]<<endl;
		}
		else{
	//		cerr<<" "<<output[i];
		}
	}
	return output;
}
/*}}}*/

/*inline TString gGet_InputFile_Var(const TString& aFile_Name,const Int_t& aLine){{{*/
inline TString gGet_InputFile_Var(const TString& aFile_Name,const Int_t& aLine)
{
	ifstream file(aFile_Name.Data());
	TString content;
	Int_t i=0;
	TString outstring;
	while ( content.ReadLine(file) )
	{
		if ( i==aLine )
		{
			outstring=content;
			break;
		}
		i++;
		if ( !file.good() )
			break;
	}
	file.close();
	if ( outstring=="\"\"" )
		outstring="";
	return outstring;
}
/*}}}*/

char* gTarget, *gKin, *gArm, *gCom;
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
				case 'T':
					if(*(++argptr))
						gTarget = argptr;
					err=0;
					break;
				case 'K':
					if(*(++argptr))
						gKin = argptr;
					err=0;
					break;
				case 'A':
					if(*(++argptr))
						gArm = argptr;
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
					cerr <<"-T[He3,He4,C12,Ca40,Ca48] Target" <<endl;
					cerr <<"-K[3.1,3.2] Kinematic"<<endl;
					cerr <<"-A[L, R] HRS Arm"<<endl;
					cerr <<"-C[DpTh40] Comment, such as Dp=4%"<<endl;
					cerr <<endl<<"================================================"<<endl;
					noStop=false;
					err=-1;
					goto OUT;
					break;
				default:
					cerr <<"================================================"<<endl;
					cerr <<"Option: " <<endl;
					cerr <<"-T[He3,He4,C12,Ca40,Ca48] Target" <<endl;
					cerr <<"-K[3.1,3.2] Kinematic"<<endl;
					cerr <<"-A[L, R] HRS Arm"<<endl;
					cerr <<"-C[DpTh40] Comment, such as Dp=4%"<<endl;
					cerr <<endl<<"================================================"<<endl;
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

//void Yield_Run(){
int main(int argc,char** argv){
	int pserr = getargs(argc,argv);

     TString Target = gTarget;
	 TString aKin = gKin; TString Kin = Form("Kin%s",aKin.Data());
	 TString Arm = gArm;
	 TString Com = gCom;

/*	
	//Here are where you need to change
	TString Target = "He3"; cerr<<"--- Target (H2,He3,He4,C12,Ca40 and Ca48) = "; cin >> Target;
	if(Target!="H2"&&Target!="He3"&&Target!="He4"&&Target!="C12"&&Target!="Ca40"&&Target!="Ca48") exit(-1);
	TString Kin = "3.1"; cerr<<"--- Kin (e.g. 3.1) = "; cin >> Kin; Kin = "Kin"+Kin;
    if(Kin!="Kin3.1"&&Kin!="Kin3.2"&&Kin!="Kin4.1"&&Kin!="Kin4.2"&&Kin!="Kin5.0"&&Kin!="Kin5.05"&&Kin!="Kin5.1"&&Kin!="Kin5.2"&&Kin!="Kin6.5")  exit(-2);
	TString Arm = "L";      cerr<<"--- Arm (L or R) = "; cin >> Arm;
	TString Com = "DpTh40"; cerr<<"--- Comment = "; cin >> Com;
*/
	if(Arm!="L"&&Arm!="R") exit(-3);

	/*Read In RunNo{{{*/
	TString sTarget="", sKin = "";
	if(Target.Contains("H")){
		sTarget = Target; sTarget.ReplaceAll("H","h");
	}
	if(Target.Contains("C")){
		sTarget = Target; sTarget.ReplaceAll("C","c");
	} 
	sKin=Kin; sKin.ReplaceAll("Kin","kin");

	TString Input_Parameter = Form("%s/xs_input/%sHRS/xs_%s_%s.dat",XGT2_DIR.Data(),Arm.Data(),sTarget.Data(),sKin.Data());
	ifstream input_par(Input_Parameter.Data());
	TString one_line;

	int RunNo_Temp = gGet_InputFile_Var(Input_Parameter,0).Atoi();
	vector<int> vRunNoChain = gGet_RunNoChain(RunNo_Temp, gGet_InputFile_Var(Input_Parameter,15));  

	const int N=vRunNoChain.size();
	int RunNo[N];
	for(int i=0;i<N;i++)
		RunNo[i] = vRunNoChain[i];
	/*}}}*/

	TString InputFileName = Form("%s/results/%s/%s/%s_%s_%s_Yield_%s.root",
			XGT2_DIR.Data(),Target.Data(),Com.Data(),Arm.Data(),Kin.Data(),Target.Data(),Com.Data());
	cerr<<" &&& Reading in RootFile --> "<<InputFileName.Data()<<endl;

	double Ep, xbj;
	double Y_EX[N], Y_EX_Stat[N], Y_EX_Sys[N], Y_MC[N],Y_MC_Stat[N],Y_MC_Sys[N]; 

	TFile *f1 = new TFile(InputFileName.Data(),"r");
	if(f1->IsZombie())	exit(-4);
	TTree *T = (TTree*) gDirectory->Get("T");

	const int N_Evt = T->GetEntries();

	T->SetBranchAddress("Ep",  &Ep);
	T->SetBranchAddress("xbj", &xbj);
	for(int i=0;i<N;i++){
		T->SetBranchAddress(Form("Y_EX_%d",     RunNo[i]), &Y_EX[i]);
		T->SetBranchAddress(Form("Y_EX_Stat_%d",RunNo[i]), &Y_EX_Stat[i]);
		T->SetBranchAddress(Form("Y_EX_Sys_%d" ,RunNo[i]), &Y_EX_Sys[i]);
		T->SetBranchAddress(Form("Y_MC_%d",     RunNo[i]), &Y_MC[i]);
		T->SetBranchAddress(Form("Y_MC_Stat_%d",RunNo[i]), &Y_MC_Stat[i]);
		T->SetBranchAddress(Form("Y_MC_Sys_%d" ,RunNo[i]), &Y_MC_Sys[i]);
	}

	double Y_EX_Sum[N], Y_EX_Sum_Err[N], Y_MC_Sum[N], Y_MC_Sum_Err[N];
	for(int i=0;i<N;i++){
		Y_EX_Sum[i] = 0; 
		Y_MC_Sum[i] = 0; 
		Y_EX_Sum_Err[i] = 0; 
		Y_MC_Sum_Err[i] = 0; 
	}

	double Ep_Bin[N_Evt], Xbj_Bin[N_Evt];
	double Y_R[N][N_Evt], Y_R_Err[N][N_Evt];
	double Ep_Max = -1000.0, Ep_Min = 1000.0,Xbj_Max = -1000.0, Xbj_Min = 1000.0, Y_EX_Err=0.0, Y_MC_Err = 0.0;
	for(int j=0;j<N_Evt;j++){
		T->GetEntry(j);
		Ep_Bin[j] = Ep;
		if(Ep>Ep_Max) Ep_Max = Ep;
		if(Ep<Ep_Min) Ep_Min = Ep;

		Xbj_Bin[j] = Ep;
		if(xbj>Xbj_Max) Xbj_Max = xbj;
		if(xbj<Xbj_Min) Xbj_Min = xbj;


		for(int i=0;i<N;i++){
			Y_EX_Err = sqrt( pow(Y_EX_Stat[i],2)+pow(Y_EX_Sys[i],2) ); 
			Y_MC_Err = sqrt( pow(Y_MC_Stat[i],2)+pow(Y_MC_Sys[i],2) );
			if(Y_EX[i]>1e-30&&Y_MC[i]>1e-30){

				if(Y_EX[i]>1e-30&&Y_MC[i]>1e-30){
					Y_R[i][j] = Y_EX[i] / Y_MC[i];
					Y_R_Err[i][j] = Y_R[i][j]*sqrt( pow(Y_EX_Err/Y_EX[i],2)+pow(Y_MC_Err/Y_MC[i],2) );
				}
				else{
					Y_R[i][j] = 0.0;
					Y_R_Err[i][j] = 0.0;
				}
				cerr<<Form("Run#%d: Ep = %f, Y_R = %e (%e)", RunNo[i],Ep_Bin[j], Y_R[i][j],Y_R_Err[i][j])<<endl;
				if(Y_R[i][j]>1e-3&&Y_R_Err[i][j]>1e-10){
					Y_EX_Sum[i] += Y_EX[i]; 
					Y_MC_Sum[i] += Y_MC[i]; 
					Y_EX_Sum_Err[i] += Y_EX_Err*Y_EX_Err; 
					Y_MC_Sum_Err[i] += Y_MC_Err*Y_MC_Err;
				}	
			}
		}
	}

	for(int i=0;i<N;i++){
		Y_EX_Sum_Err[i] = sqrt(Y_EX_Sum_Err[i]); 
		Y_MC_Sum_Err[i] = sqrt(Y_MC_Sum_Err[i]); 
	}
	
	double Y_Min =1000., Y_Max = -1000.,Y_R_Avg =0.0;
	double Y_R_Sum[N][1], Y_R_Sum_Err[N][1];
	int count = 0;
	for(int i=0;i<N;i++){
		if(Y_EX_Sum[i]>1e-30&&Y_MC_Sum[i]>1e-30){
			Y_R_Sum[i][0] = Y_EX_Sum[i] / Y_MC_Sum[i];
			Y_R_Sum_Err[i][0] = Y_R_Sum[i][0]*sqrt( pow(Y_EX_Sum_Err[i]/Y_EX_Sum[i],2)+pow(Y_MC_Sum_Err[i]/Y_MC_Sum[i],2) );
			cerr<<Form("Run#%d: Y_R = %e (%e)", RunNo[i],Y_R_Sum[i][0],Y_R_Sum_Err[i][0])<<endl;
			if(Y_Min>(Y_R_Sum[i][0]+Y_R_Sum_Err[i][0]))
				Y_Min =Y_R_Sum[i][0]+Y_R_Sum_Err[i][0]; 
			if(Y_Max<(Y_R_Sum[i][0]+Y_R_Sum_Err[i][0]))
				Y_Max =Y_R_Sum[i][0]+Y_R_Sum_Err[i][0]; 

			Y_R_Avg += Y_R_Sum[i][0]; count++;
		}else{
			Y_R_Sum[i][0] = 0.0;
			Y_R_Sum_Err[i][0] = 0.0;
		}
	}
    Y_R_Avg /= count;
	cerr<<"------- Y_R_Avg = "<<Y_R_Avg<<endl;

	/*Yield  Ratio{{{*/
    double Y_Min_Out = 0.75, Y_Max_Out = 1.25;
    if(Target=="He4"){
		Y_Min_Out = 0.55; Y_Max_Out = 1.05;
	}
	TCanvas *c3 = new TCanvas("c3","c3",1920,1080);
	//c3->Divide(1,2);
	c3->cd(1);
	TH2F* r = new TH2F("r","",1000,-1,N+2,1000,Y_Min_Out, Y_Max_Out);
	r->SetStats(kFALSE);
	r->SetXTitle("RunNo");
	r->GetXaxis()->CenterTitle(1);
	r->GetXaxis()->SetTitleFont(32);
	r->GetXaxis()->SetTitleSize(0.06);
	r->GetXaxis()->SetTitleOffset(0.6);
	r->SetYTitle("Yield_{EX}/Yield_{MC}");
	r->GetYaxis()->CenterTitle(1);
	r->GetYaxis()->SetTitleFont(32);
	r->GetYaxis()->SetTitleSize(0.04);
	r->Draw(); 
	

	double RN[N][1],RN_All[N];
	double Y_R_All[N], Y_R_All_Err[N];
	for(int i=0;i<N;i++){
		RN[i][0]=i+1;
		RN_All[i]=i+1;
        Y_R_All[i] = Y_R_Sum[i][0];
        Y_R_All_Err[i] = Y_R_Sum_Err[i][0];
	}

	gStyle->SetOptFit(1);
	TGraphErrors* pl_r_all = new TGraphErrors(N,RN_All,Y_R_All,0,Y_R_All_Err);
	pl_r_all->SetLineColor(5);
	pl_r_all->Draw("P");
	pl_r_all->Fit("pol0","Q");
	TF1* fit = pl_r_all->GetFunction("pol0");
	fit->SetLineColor(5);
    fit->SetLineColor(7);
    fit->SetLineWidth(1.5);
	fit->SetLineStyle(9);
	Y_R_Avg = fit->GetParameter(1);   
    double Chi = fit->GetChisquare();
	int N_DOF = fit->GetNdim();

	cerr<<Form("%s  %s  HRS-%s %10d %10.3f %10.3f/%d", Target.Data(), Kin.Data(), Arm.Data(), N, Y_R_Avg, Chi,N_DOF)<<endl;

	TGraphErrors *pl_r[N];
	int N_Color = 1, N_Style = 20;
	for(int i=0;i<N;i++){
		pl_r[i] = new TGraphErrors(1,&RN[i][0],&Y_R_Sum[i][0],0,&Y_R_Sum_Err[i][0]);
		pl_r[i]-> SetMarkerSize(1.5);
		N_Color++;
		if(N_Color>9&&N_Color<37)	N_Color = 37;
		if(N_Color>49)	N_Color = 1;
		if(N_Color==2) N_Color = 28;
		if(Y_R_Sum[i][0]<0.1) N_Color = 2;
		if(N_Color==5)		
			pl_r[i]-> SetMarkerColor(1);
        else
			pl_r[i]-> SetMarkerColor(N_Color);
		N_Style++;
		if(N_Style>30)	N_Style = 2;
		if(N_Style>5&&N_Style<20) N_Style = 20;
		pl_r[i]-> SetMarkerStyle(N_Style);
		pl_r[i]->Draw("P");
	}

	TLatex t1; N_Color = 1;
	t1.DrawLatex(0.2*(1),Y_Max_Out*1.01,Form("%s@ %s on HRS-%s(%d Runs)", Target.Data(),Kin.Data(),Arm.Data(),N));
	for(int i=0;i<N;i++){
        t1.SetTextAngle(90);
		if(N>10)
		t1.SetTextSize(0.02);
		else
			t1.SetTextSize(0.03);
		N_Color++;
		if(N_Color>9&&N_Color<37)	N_Color = 37;
		if(N_Color>49)	N_Color = 1;
		if(N_Color==2) N_Color = 28;
		if(Y_R_Sum[i][0]<0.1) N_Color = 2;
		if(N_Color==5) 
			t1.SetTextColor(1);
		else
			t1.SetTextColor(N_Color);
		//if(Y_R_Sum[i][0]<1e-3){
		if(Y_R_Sum[i][0]<Y_Min_Out || Y_R_Sum[i][0]>Y_Max_Out){
			t1.SetTextColor(2);
			t1.DrawLatex(i+1,Y_R_Avg*0.91,Form("#%d (Bad?)", RunNo[i]));
		}
		else if((Y_R_Sum[i][0]+Y_R_Sum_Err[i][0])>Y_Max_Out)
			t1.DrawLatex(i+1,(Y_R_Sum[i][0]-Y_R_Sum_Err[i][0])*0.95,Form("#%d", RunNo[i]));
		else
			t1.DrawLatex(i+1,(Y_R_Sum[i][0]+Y_R_Sum_Err[i][0])*1.01,Form("#%d", RunNo[i]));
	}

	TLine *line = new TLine(0, Y_R_Avg, N+1, Y_R_Avg);
    line->SetLineColor(7);
    line->SetLineWidth(1.5);
	line->SetLineStyle(9);
	line->Draw("same");
/*
	c3->cd(2);
	TH2F* r1= new TH2F("r1","",1000,-2,N+2,1000,Y_Min_Out, Y_Max_Out);
	r1->SetStats(kFALSE);
	r1->SetXTitle("RunNo");
	r1->GetXaxis()->CenterTitle(1);
	r1->GetXaxis()->SetTitleFont(32);
	r1->GetXaxis()->SetTitleSize(0.06);
	r1->SetYTitle("Yield_{EX}/Yield_{MC}");
	r1->GetYaxis()->CenterTitle(1);
	r1->GetYaxis()->SetTitleFont(32);
	r1->GetYaxis()->SetTitleSize(0.04);
	r1->Draw();
		
	TGraphErrors *pl_r1[N];
	N_Color = 1; N_Style = 20;
	for(int i=0;i<N;i++){
		pl_r1[i] = new TGraphErrors(1,&RN[i][0],&Y_R_Sum[i][0],0,&Y_R_Sum_Err[i][0]);
		pl_r1[i]-> SetMarkerSize(1.5);
		N_Color++;
		if(N_Color>9&&N_Color<37)	N_Color = 37;
		if(N_Color>49)	N_Color = 1;
		if(N_Color==2) N_Color = 28;
		if(Y_R_Sum[i][0]<0.1) N_Color = 2;
		if(N_Color==5)		
			pl_r1[i]-> SetMarkerColor(1);
        else
			pl_r1[i]-> SetMarkerColor(N_Color);
		N_Style++;
		if(N_Style>30)	N_Style = 2;
		if(N_Style>5&&N_Style<20) N_Style = 20;
		pl_r1[i]-> SetMarkerStyle(N_Style);
		pl_r1[i]->Draw("P");
	}

	TLatex t2; N_Color = 1;
	t2.DrawLatex(0.3*N,Y_Max*0.98,Form("%s@ %s on HRS-%s", Target.Data(),Kin.Data(),Arm.Data()));
	for(int i=0;i<N;i++){
        t2.SetTextAngle(90);
		if(N>10)
		t2.SetTextSize(0.02);
		else
			t2.SetTextSize(0.03);
		N_Color++;
		if(N_Color>9&&N_Color<37)	N_Color = 37;
		if(N_Color>49)	N_Color = 1;
		if(N_Color==2) N_Color = 28;
		if(Y_R_Sum[i][0]<0.1) N_Color = 2;
		if(N_Color==5) 
			t2.SetTextColor(1);
		else
			t2.SetTextColor(N_Color);
		//if(Y_R_Sum[i][0]<1e-3){
		if(Y_R_Sum[i][0]<Y_Min || Y_R_Sum[i][0]>Y_Max){
			t2.SetTextColor(2);
			t2.DrawLatex(i+1,(Y_R_Sum[i-1][0]+Y_R_Sum_Err[i-1][0])*0.91,Form("#%d (Bad?)", RunNo[i]));
		}
		else if((Y_R_Sum[i][0]+Y_R_Sum_Err[i][0])>Y_Max)
			t2.DrawLatex(i+1,(Y_R_Sum[i][0]-Y_R_Sum_Err[i][0])*0.95,Form("#%d", RunNo[i]));
		else
			t2.DrawLatex(i+1,(Y_R_Sum[i][0]+Y_R_Sum_Err[i][0])*1.01,Form("#%d", RunNo[i]));
	}
	line->Draw("same");

*/

	c3->Print(Form("YS_%s_%s_%s_%s.pdf",Target.Data(),Arm.Data(),Kin.Data(),Com.Data()));
	c3->Print(Form("YS_%s_%s_%s_%s.png",Target.Data(),Arm.Data(),Kin.Data(),Com.Data()));

	f1->Close();

}
