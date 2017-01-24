////////////////////////////////////////
// VDC Tracking Efficiency Study
// For x>2 experiment
// -- Zhihong Ye 07/08/2011
////////////////////////////////////////
/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
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
#include <TRandom3.h>
//#include <TMatrix.h>
/*}}}*/

using namespace std;

TString ROOTFILES_DIR = "/work/halla/e08014/disk1/Root_0904";
TString ROOTFILES_NAME = "e08014_R4";

Bool_t gIsExist(int kRunno);
Bool_t gIsBad(int kRunno);
/*inline TChain* gAddTree(const Int_t aRunNo,const TString& aArm,TString aTreeName="T"){{{*/
inline TChain* gAddTree(const Int_t aRunNo, TString aTreeName="T")
{
	TChain* aTree=new TChain(aTreeName.Data(),"Main Tree");
	Int_t index=0;
	TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),aRunNo);
	if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
	{
		aTree->Add(Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo));
	}
	else{  
		cerr<<Form("       *ERROR*, no rootfile founded for Run#%d",aRunNo)<<endl;
		return 0;
	}
	index++;
	File_Form=Form("%s_%d_%d.root",ROOTFILES_NAME.Data(),aRunNo,index);
	while ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
	{
		aTree->Add(Form("%s/%s_%d_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo,index));
		index++;
		File_Form=Form("%s_%d_%d.root",ROOTFILES_NAME.Data(),aRunNo,index);
	}
	aTree->LoadTree(0);

	return aTree;
}
/*}}}*/

//void R_vdc_eff(){
int main(){
	gSystem->Exec("cp -p Bad_R.txt Bad_R1.txt");
	gSystem->Exec("cp -p R_VDC_Eff_All.dat R_VDC_Eff_All.dat.old");
	
	Int_t  runnum;    TString filename;
	ofstream myfile1; myfile1.open("R_VDC_Eff_All1.dat");
	ofstream badfile; badfile.open("Bad_R.txt");
	ifstream runfile("Run_ChartR.txt"); 

	const int NUMR = 1000;
	TString TriggerCut1 = "((DBB.evtypebits>>1)&1)";
	TString electron_cut = "R.cer.asum_c>700.&&R.cer.asum_c<2100.&&(R.sh.asum_c+R.ps.asum_c)>700.0";
	TString cut_t1 =TriggerCut1 + "&&"+electron_cut;

	TString cut_n0,cut_n1,cut_n2,cut_n3,cut_n4;
	cut_n1 = "&&R.tr.n==1";  cut_n2 = "&&R.tr.n==2";
	cut_n3 = "&&R.tr.n==3";  cut_n4 = "&&R.tr.n==4"; cut_n0 = "&&R.tr.n==0";
	//  TString cut_nn = "R.vdc.u1.nhit==5 &&R.vdc.u2.nhit==5 &&R.vdc.v1.nhit==5 && R.vdc.v2.nhit==5 &&";
	//  TString cut_nn = "R.vdc.u1.nhit>=4 &&R.vdc.u2.nhit>=4 &&R.vdc.v1.nhit>=4 && R.vdc.v2.nhit>=4 &&";
	//  cut_nn += "R.vdc.u1.nhit<=6 &&R.vdc.u2.nhit<=6 &&R.vdc.v1.nhit<=6 && R.vdc.v2.nhit<=6 &&";
	TString cut_nn = "R.s1.nthit==1&&R.s2.nthit==1&&";

	double rnum[NUMR]; 
	double T1track0[NUMR], T1track1[NUMR], T1track2[NUMR], T1track3[NUMR], T1track4[NUMR];
	double track0,track1,track2,track3,track4;

	int jj;  double evt,evt0,evt1,evt2,evt3,evt4;
	int t1r,t3r,t6r,t7r,tgt;
	TString Kin[NUMR],Target[NUMR],kine,target,comment;
	double Mom[NUMR],Angle[NUMR];
	double mom,angle,kin;

	TFile *f1 = new TFile("R_VDC_Eff_All.root","UPDATE");
	TTree *vdc = new TTree("vdc","VDC Efficiency");

	vdc->Branch("runnum",&runnum, "runnum/I");
	vdc->Branch("mom",&mom, "mom/D");
	vdc->Branch("angle",&angle, "angle/D");
	vdc->Branch("t1r",&t1r, "t1r/I");
	vdc->Branch("track0", &track0,"track0/D");
	vdc->Branch("track1", &track1,"track1/D");
	vdc->Branch("track2", &track2,"track2/D");
	vdc->Branch("track3", &track3,"track3/D");
	vdc->Branch("track4", &track4,"track4/D");

	Bool_t skip;
	Bool_t IsExt = kFALSE;
	while(!(runfile.eof())){

		runfile >> runnum >> kine >> mom >> angle >> target >>comment; 
		IsExt = gIsExist(runnum);

		if(!IsExt){

			cerr << "Working on Run#" << runnum 
				<<" for target "<< target<<" at "<<kine
			<<" with Momentum "<<mom<<"..."
			<<endl;

		skip = kFALSE;
		if(target=="H2")
			tgt=2;
		else if(target=="He3")
			tgt=3;
		else if(target=="He4")
			tgt=4;
		else if(target=="C12")
			tgt=12;
		else if(target=="Ca40")
			tgt=40;
		else if(target=="Ca48")
			tgt=48;
		else if(target=="Dummy_20cm")
			tgt=2720;
		else if(target=="Dummy_10cm")
			tgt=2710;
		else if(target=="Dummy_4cm")
			tgt=2704;
		else if(target=="Optics")
			tgt=55;
		else
			skip = kTRUE;

		if(skip)
			cerr << "Not for HRS-R, skipping ..."<<endl;
		else{
			TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),runnum);
			if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) ){
				rnum[jj] =  runnum;
				Mom[jj] = mom;
				Target[jj] = target;
				Kin[jj] = kin;
				Angle[jj] = angle;

				TChain *T = gAddTree(runnum,"T");

				evt  = T->GetEntries(cut_nn+cut_t1);
				evt0 = T->GetEntries(cut_nn+cut_t1+cut_n0);
				evt1 = T->GetEntries(cut_nn+cut_t1+cut_n1);
				evt2 = T->GetEntries(cut_nn+cut_t1+cut_n2);
				evt3 = T->GetEntries(cut_nn+cut_t1+cut_n3);
				evt4 = T->GetEntries(cut_nn+cut_t1+cut_n4);

				T1track0[jj] = 100.0*evt0/evt; track0 = 100.0*evt0/evt;
				T1track1[jj] = 100.0*evt1/evt; track1 = 100.0*evt1/evt;
				T1track2[jj] = 100.0*evt2/evt; track2 = 100.0*evt2/evt;
				T1track3[jj] = 100.0*evt3/evt; track3 = 100.0*evt3/evt;
				T1track4[jj] = 100.0*evt4/evt; track4 = 100.0*evt4/evt;

				cerr <<"Tracking Eff for run n=1,2,3,4 is "
					<< T1track0[jj] <<" "
					<< T1track1[jj] <<" "
					<< T1track2[jj] <<" "
					<< T1track3[jj] <<" "
					<< T1track4[jj] <<endl;

				myfile1 << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(1) << runnum << " ";
				myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << Mom[jj] << " ";
				myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << 100.*evt0/evt << " ";
				myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << 100.*evt1/evt << " ";
				myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << 100.*evt2/evt << " ";
				myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << 100.*evt3/evt << " ";
				myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << 100.*evt4/evt << endl;

				jj++; T->Delete();
				f1->cd(); vdc->Fill();
			}
			else{
				cerr<<" **** No Rootfile found for Run#="<< runnum<<endl;
				badfile << runnum <<endl;
			}
		}
	}
	}
	f1->cd(); vdc->Write("",TObject::kOverwrite);
	f1->Close(); runfile.close(); myfile1.close(); badfile.close();
}


Bool_t gIsExist(int kRunno){
	ifstream infile("R_VDC_Eff_All.dat");
	int runno = 0;
	double aa;

	bool IsExt = kFALSE;

	while(!infile.eof()){
		infile >> runno >> aa >> aa >> aa >> aa >> aa >> aa;

		if(runno == kRunno){
			IsExt = kTRUE;
			break;
		}
	}

	return IsExt;
}

Bool_t gIsBad(int kRunno){
	ifstream infile("Bad_R1.txt");
	int runno = 0;
	bool IsExt = kFALSE;

	while(!infile.eof()){
		infile >> runno;

		if(runno == kRunno){
			IsExt = kTRUE;
			break;
		}
	}

	return IsExt;
}
