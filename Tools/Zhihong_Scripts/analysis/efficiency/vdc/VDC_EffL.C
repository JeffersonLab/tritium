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

TString ROOTFILES_DIR = "/work/halla/e08014/disk1/Rootfiles";
TString ROOTFILES_NAME = "e08014";

Bool_t gIsExist(int kRunno);
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

int main(){  
	gSystem->Exec("cp -p Bad_L.txt Bad_L1.txt");
	gSystem->Exec("cp -p L_VDC_Eff_All.dat L_VDC_Eff_All.dat.old");
	
	Int_t  runnum;    TString filename;
	ofstream myfile1; myfile1.open("L_VDC_Eff_All1.dat");
	ofstream badfile; badfile.open("Bad_L.txt");
	ifstream runfile("Run_ChartL.txt");
	const int NUMR = 1000;

	TString TriggerCut3 = "((DBB.evtypebits>>3)&1)";
	TString electron_cut = "L.cer.asum_c>700.&&L.cer.asum_c<2100.&&(L.prl1.asum_c+L.prl2.asum_c)>700.";
	TString cut_t3 = TriggerCut3+ "&&" +electron_cut;

	TString cut_n0,cut_n1,cut_n2,cut_n3,cut_n4;
	cut_n0 = "&&L.tr.n==0";  cut_n1 = "&&L.tr.n==1";  cut_n2 = "&&L.tr.n==2";
	cut_n3 = "&&L.tr.n==3";  cut_n4 = "&&L.tr.n==4";
	//  TString cut_nn = "L.vdc.u1.nhit==5 &&L.vdc.u2.nhit==5 &&L.vdc.v1.nhit==5 && L.vdc.v2.nhit==5 &&";
	//  TString cut_nn = "L.vdc.u1.nhit>=5 &&L.vdc.u2.nhit>=5 &&L.vdc.v1.nhit>=5 && L.vdc.v2.nhit>=5 &&";
	//  cut_nn += "L.vdc.u1.nhit<=7 &&L.vdc.u2.nhit<=7 &&L.vdc.v1.nhit<=7 && L.vdc.v2.nhit<=7 &&";
	TString cut_nn = "L.s1.nthit==1&&L.s2.nthit==1&&";

	double rnum[NUMR]; 
	double T3track0[NUMR], T3track1[NUMR], T3track2[NUMR], T3track3[NUMR], T3track4[NUMR];
	double track0,track1,track2,track3,track4;

	int jj;  double evt,evt0,evt1,evt2,evt3,evt4;
	int t1r,t3r,t6r,t7r,tgt;
	TString Kin[NUMR],Target[NUMR],kine,target,comment;
	double Mom[NUMR],Angle[NUMR];
	double mom,angle,kin;

	TFile *f1 = new TFile("L_VDC_Eff_All.root","UPDATE");
	TTree *vdc = new TTree("vdc","VDC Efficiency");

	vdc->Branch("runnum",&runnum, "runnum/I");
	vdc->Branch("mom",&mom, "mom/D");
	vdc->Branch("angle",&angle, "angle/D");
	vdc->Branch("t3r",&t3r, "t3r/I");
	vdc->Branch("track0", &track0,"track0/D");
	vdc->Branch("track1", &track1,"track1/D");
	vdc->Branch("track2", &track2,"track2/D");
	vdc->Branch("track3", &track3,"track3/D");
	vdc->Branch("track4", &track4,"track4/D");

	Bool_t skip;
	Bool_t IsExt;
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
			else if(target=="He3"||target=="He3R")
				tgt=3;
			else if(target=="He4")
				tgt=4;
			else if(target=="C12"||target=="C12R")
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
			else if(target=="BeO")
				tgt=9;
			else
				skip = kTRUE;

			if(skip)
				cerr << "Not for HRS-L, skipping ..."<<endl;
			else{
				TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),runnum);
				if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) ){
					rnum[jj] =  runnum;
					Mom[jj] = mom;
					Target[jj] = target;
					Kin[jj] = kine;
					Angle[jj] = angle;

					TChain *T = gAddTree(runnum,"T");

					evt  = T->GetEntries(cut_nn+cut_t3);
					evt0 = T->GetEntries(cut_nn+cut_t3+cut_n0);
					evt1 = T->GetEntries(cut_nn+cut_t3+cut_n1);
					evt2 = T->GetEntries(cut_nn+cut_t3+cut_n2);
					evt3 = T->GetEntries(cut_nn+cut_t3+cut_n3);
					evt4 = T->GetEntries(cut_nn+cut_t3+cut_n4);

					T3track0[jj] = 100.0*evt0/evt; track0 = 100.0*evt0/evt;
					T3track1[jj] = 100.0*evt1/evt; track1 = 100.0*evt1/evt;
					T3track2[jj] = 100.0*evt2/evt; track2 = 100.0*evt2/evt;
					T3track3[jj] = 100.0*evt3/evt; track3 = 100.0*evt3/evt;
					T3track4[jj] = 100.0*evt4/evt; track4 = 100.0*evt4/evt;

					double track0_err = 100 * sqrt( abs(1/evt0 - 1/evt)  );
					double track1_err = 100 * sqrt( abs(1/evt1 - 1/evt)  );
					double track2_err = 100 * sqrt( abs(1/evt2 - 1/evt)  );
					double track3_err = 100 * sqrt( abs(1/evt3 - 1/evt)  );
					double track4_err = 100 * sqrt( abs(1/evt4 - 1/evt)  );

					cerr <<"Tracking Eff for run n=1,2,3,4 is "
						<< T3track0[jj] <<" "
						<< T3track1[jj] <<" "
						<< T3track2[jj] <<" "
						<< T3track3[jj] <<" "
						<< T3track4[jj] <<endl;

					myfile1 << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(0) << runnum << "   ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << Mom[jj] << " ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << 100.*evt0/evt << " ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << track0_err << " ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << 100.*evt1/evt << " ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << track1_err << " ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << 100.*evt2/evt << " ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << track2_err << " ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << 100.*evt3/evt << " ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << track3_err << " ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << 100.*evt4/evt << " ";
					myfile1 << setiosflags(ios::left) << setw(7) << setiosflags(ios::fixed) << setprecision(4) << track4_err << endl;

					jj++; T->Delete();
					f1->cd(); vdc->Fill();
				}
				else{
					badfile <<runnum<<endl;
                    cerr<<" **** No Rootfile found for Run#="<< runnum<<endl;
				}
			}
		}
	}
	f1->cd(); vdc->Write("",TObject::kOverwrite);
	f1->Close(); runfile.close();myfile1.close();badfile.close();

}


Bool_t gIsExist(int kRunno){
	ifstream infile("L_VDC_Eff_All.dat");
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
	ifstream infile("Bad_L1.txt");
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
