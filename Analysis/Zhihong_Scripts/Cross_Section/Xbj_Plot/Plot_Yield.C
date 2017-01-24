
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

//void Plot_Yield(){
int main(){
	//Here are where you need to change
	//const TString Com = "DpTh40";
	TString Target = "He3"; cerr<<"--- Target (H2,He3,He4,C12,Ca40 and Ca48) = "; cin >> Target;
	TString Kin = "Kin3.1"; cerr<<"--- Kin (e.g. 3.1) = "; cin >> Kin; Kin="Kin"+Kin;
	TString Arm = "L";      cerr<<"--- Arm (L or R) = "; cin >> Arm;
	TString Com = "J02";    cerr<<"--- Comment = "; cin >> Com;

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

	double Ep,xbj;
	double Y_EX[N], Y_EX_Stat[N], Y_EX_Sys[N], Y_MC[N],Y_MC_Stat[N],Y_MC_Sys[N]; 

	TFile *f1 = new TFile(InputFileName.Data(),"r");
	if(f1->IsZombie()){
	   cerr<<"*** I cannot fine ROOTFile: "<<InputFileName.Data()<<" ***"<<endl;
	   exit(-1);
	}
		
	TTree *T = (TTree*) gDirectory->Get("T");

	const int N_Evt = T->GetEntries();

	T->SetBranchAddress("Ep", &Ep);
	T->SetBranchAddress("xbj", &xbj);
	for(int i=0;i<N;i++){
		T->SetBranchAddress(Form("Y_EX_%d",     RunNo[i]), &Y_EX[i]);
		T->SetBranchAddress(Form("Y_EX_Stat_%d",RunNo[i]), &Y_EX_Stat[i]);
		T->SetBranchAddress(Form("Y_EX_Sys_%d" ,RunNo[i]), &Y_EX_Sys[i]);
		T->SetBranchAddress(Form("Y_MC_%d",     RunNo[i]), &Y_MC[i]);
		T->SetBranchAddress(Form("Y_MC_Stat_%d",RunNo[i]), &Y_MC_Stat[i]);
		T->SetBranchAddress(Form("Y_MC_Sys_%d" ,RunNo[i]), &Y_MC_Sys[i]);
	}

	double Y_R[N][N_Evt], Y_R_Err[N][N_Evt];
	double Ep_Bin[N_Evt],xbj_Bin[N_Evt];
	double Ep_Max = -1000.0, Ep_Min = 1000.0, xbj_Max = -1000.0, xbj_Min = 1000.0, Y_EX_Err=0.0, Y_MC_Err = 0.0;
	for(int j=0;j<N_Evt;j++){
		T->GetEntry(j);
		Ep_Bin[j] = Ep;
		if(Ep>Ep_Max) Ep_Max = Ep;
		if(Ep<Ep_Min) Ep_Min = Ep;

		xbj_Bin[j] = xbj;
		if(xbj>xbj_Max) xbj_Max = xbj;
		if(xbj<xbj_Min) xbj_Min = xbj;


		for(int i=0;i<N;i++){
			//cerr<<Form("Run#%d: Ep = %f, Y_R = %e (%e)", RunNo[i],Ep_Bin[j], Y_EX[i],Y_EX_Stat[i])<<endl;
			Y_EX_Err = sqrt( pow(Y_EX_Stat[i],2)+pow(Y_EX_Sys[i],2) ); 
			Y_MC_Err = sqrt( pow(Y_MC_Stat[i],2)+pow(Y_MC_Sys[i],2) );
			if(Y_EX[i]>1e-30&&Y_MC[i]>1e-30){
				Y_R[i][j] = Y_EX[i] / Y_MC[i];
				Y_R_Err[i][j] = Y_R[i][j]*sqrt( pow(Y_EX_Err/Y_EX[i],2)+pow(Y_MC_Err/Y_MC[i],2) );
			}
			else{
				Y_R[i][j] = 0.0;
				Y_R_Err[i][j] = 0.0;
			}
		}
	}

	/*Yield  Ratio{{{*/
	double Y_Min =0, Y_Max = 2.0;
    if(N>10)
		Y_Min = -1.0;
	TCanvas *c3 = new TCanvas("c3","c3",1200,800);
	TH2F* r = new TH2F("r","",1000,xbj_Min-0.2,xbj_Max+0.2,1000,Y_Min, Y_Max);
	r->SetStats(kFALSE);
	r->SetXTitle("X");
	r->GetXaxis()->CenterTitle(1);
	r->GetXaxis()->SetTitleFont(32);
	r->GetXaxis()->SetTitleSize(0.06);
	r->GetXaxis()->SetTitleOffset(0.7);
	r->SetYTitle("Yield_{EX}/Yield_{MC}");
	r->GetYaxis()->CenterTitle(1);
	r->GetYaxis()->SetTitleFont(32);
	r->GetYaxis()->SetTitleSize(0.05);
	r->GetYaxis()->SetTitleOffset(0.7);
	r->Draw(); 

	TGraphErrors *pl_r[N];
	int N_Color = 1, N_Style = 20;
	for(int i=0;i<N;i++){
		pl_r[i] = new TGraphErrors(N_Evt,&xbj_Bin[0],&Y_R[i][0],0,&Y_R_Err[i][0]);
		pl_r[i]-> SetMarkerSize(1.5);
		N_Color++;
		if (N_Color==5) N_Color++;
		if(N_Color>9&&N_Color<37)	N_Color = 37;
		if(N_Color>49)	N_Color = 1;
		pl_r[i]-> SetMarkerColor(N_Color);
		N_Style++;
		if(N_Style>30)	N_Style = 2;
		if(N_Style>5&&N_Style<20) N_Style = 20;
		pl_r[i]-> SetMarkerStyle(N_Style);
	}

	TMultiGraph *mg3 = new TMultiGraph();
	TLegend *l3 = new TLegend(0.35,0.17,0.60, 0.17+N*0.03,Form("%s Target",Target.Data()));
	l3->SetBorderSize(0);
	l3->SetTextSize(0.025);
	if(N>10)
		l3->SetTextSize(0.015);
	l3->SetTextFont(22);

	for(int i=0;i<N;i++){
		l3->AddEntry(pl_r[i], Form(" %s (HRS-%s) #%d",Kin.Data(), Arm.Data(), RunNo[i]),"p");
		mg3->Add(pl_r[i]);
	}
	mg3->Draw("P");
	l3->Draw();
	/*}}}*/
	
	c3->Print(Form("%s_YR_%s_%s_%s.pdf",Target.Data(),Arm.Data(),Kin.Data(),Com.Data()));
	c3->Print(Form("%s_YR_%s_%s_%s.png",Target.Data(),Arm.Data(),Kin.Data(),Com.Data()));

//	f1->Close();
}
