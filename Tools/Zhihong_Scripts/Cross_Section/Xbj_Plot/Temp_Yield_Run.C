
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
			cerr<<"      RunNoChain is "<<output[i];
		}
		else if ( i==0 && i==(output.size()-1) ){
			cerr<<"      RunNoChain is "<<output[i]<<endl;
		}
		else if ( i==(output.size()-1) ){
			cerr<<" "<<output[i]<<endl;
		}
		else{
			cerr<<" "<<output[i];
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

{
	//Here are where you need to change
	TString Target = "He3";
	TString Kin = "Kin3.1"; 
	TString Arm = "L";     
	TString Com = "DpTh40";

	if(Arm!="L"&&Arm!="R") exit(-3);

	/*Read In RunNo{{{*/
	TString sTarget="he3", sKin = "kin3.1";

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

	double Ep, Xbj;
	double Y_EX[N], Y_EX_Stat[N], Y_EX_Sys[N], Y_MC[N],Y_MC_Stat[N],Y_MC_Sys[N]; 

	TFile *f1 = new TFile(InputFileName.Data(),"r");
	if(f1->IsZombie())	exit(-4);
	TTree *T = (TTree*) gDirectory->Get("T");

	const int N_Evt = T->GetEntries();

	T->SetBranchAddress("Ep",  &Ep);
	T->SetBranchAddress("Xbj", &Xbj);
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
		if(Xbj>Xbj_Max) Xbj_Max = Xbj;
		if(Xbj<Xbj_Min) Xbj_Min = Xbj;


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
	double Y_R_Sum[N], Y_R_Sum_Err[N];
	int count = 0;
	cerr<<endl<<"============================================"<<endl;
	for(int i=0;i<N;i++){
		if(Y_EX_Sum[i]>1e-30&&Y_MC_Sum[i]>1e-30){
			Y_R_Sum[i] = Y_EX_Sum[i] / Y_MC_Sum[i];
			Y_R_Sum_Err[i] = Y_R_Sum[i]*sqrt( pow(Y_EX_Sum_Err[i]/Y_EX_Sum[i],2)+pow(Y_MC_Sum_Err[i]/Y_MC_Sum[i],2) );
			cerr<<Form("Run#%d: Y_R = %e (%e)", RunNo[i],Y_R_Sum[i],Y_R_Sum_Err[i])<<endl;
			if(Y_Min>(Y_R_Sum[i]+Y_R_Sum_Err[i]))
				Y_Min =Y_R_Sum[i]+Y_R_Sum_Err[i]; 
			if(Y_Max<(Y_R_Sum[i]+Y_R_Sum_Err[i]))
				Y_Max =Y_R_Sum[i]+Y_R_Sum_Err[i]; 

			Y_R_Avg += Y_R_Sum[i]; count++;
		}else{
			Y_R_Sum[i] = 0.0;
			Y_R_Sum_Err[i] = 0.0;
		}
	}
	cerr<<"============================================"<<endl;
    Y_R_Avg /= count;
	cerr<<"------- Y_R_Avg = "<<Y_R_Avg<<endl;

	/*Yield  Ratio{{{*/
    double Y_Min_Out = 0.5, Y_Max_Out = 1.5;
//    if(Target=="He4")
//		Y_Min = 0.65; Y_Max = 1.05;
	TCanvas *c1 = new TCanvas("c1","c1",1920,1080);
	TH2F* r_all = new TH2F("r_all","",1000,-2,N+2,1000,Y_Min_Out, Y_Max_Out);
	r_all->SetStats(kFALSE);
	r_all->SetXTitle("RunNo");
	r_all->GetXaxis()->CenterTitle(1);
	r_all->GetXaxis()->SetTitleFont(32);
	r_all->GetXaxis()->SetTitleSize(0.06);
	r_all->SetYTitle("Yield_{EX}/Yield_{MC}");
	r_all->GetYaxis()->CenterTitle(1);
	r_all->GetYaxis()->SetTitleFont(32);
	r_all->GetYaxis()->SetTitleSize(0.04);
	r_all->Draw(); 
	
	double RN[N];
	for(int i=0;i<N;i++)
		RN[i]=i+1;
	
	TGraphErrors* pl_r_all = new TGraphErrors(N,RN,Y_R_Sum,0,Y_R_Sum_Err);
	pl_r_all->Draw("P");
	pl_r_all->Fit("pol0");
//	double Y_R_Fit = pl_r_all->GetFunction("pol0")->GetParamter(0);

//	cerr<<" Y-R-Fit = "<<Y_R_Fit<<endl;

	//	f1->Close();

}
