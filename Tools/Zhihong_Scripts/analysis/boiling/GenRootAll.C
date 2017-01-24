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
#include <TRandom3.h>
//#include <TMatrix.h>
/*}}}*/

void GenRootAll(){
	gSystem->Exec("ls -l *.dat");
	TString Infile_Name = ""; cerr<<" --- Input File = "; cin >> Infile_Name;
	ifstream infile(Infile_Name);
	
	TString RootFile_Name = Infile_Name;
	RootFile_Name.ReplaceAll(".dat","_All.root");
	TFile *f1 = new TFile(RootFile_Name.Data(),"recreate");
	TTree *T  = new TTree("T", "A New Tree");

	int runsize = 0;
	infile >> runsize; cerr<<" ---  RunList = "<<runsize<<", ";
	const int Size = runsize; 
	double VZ, I, Y, Y_Err;
	int Bin, RunNo, RunChain[Size];
	for(int i=0;i<Size;i++){
		infile >> RunChain[i];
	    cerr<<" Run# = "<<RunChain[i]<< " ";
	}
    cerr<<endl;

	T->Branch("runsize",  &runsize, "runsize/I");
	T->Branch("RunNo",  &RunNo, "RunNo/I");
	T->Branch("Bin", &Bin, "Bin/I");
	T->Branch("VZ",  &VZ, "VZ/D");
	T->Branch("I",  &I, "I/D");
	T->Branch("Y",  &Y, "Y/D");
	T->Branch("Y_Err",  &Y_Err, "Y_Err/D");

	Bin=0;
    while(!(infile.eof())){
	    infile >> VZ;
		for(int i=0; i<Size; i++){
          infile >> RunNo >> I >> Y >> Y_Err;
          cerr <<" @@@ "<<VZ<<"  "<< RunNo <<"  "<< I <<"  "<< Y <<"  "<< Y_Err<<endl;

		  T->Fill();
		}
		Bin++;
	}

    f1->cd(); T->Write(); f1->Close();
}
