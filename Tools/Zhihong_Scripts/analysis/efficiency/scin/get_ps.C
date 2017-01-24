#include <TMath.h>

void get_ps(Int_t nrun=2883)
{
  gROOT->Reset();

  TString gDir = "/w/halla-2/e08014/Rootfiles/";
  Int_t  run_number;
  TString file1 = gDir + Form("e08014_less_%d.root",3590);
  TFile *f1 = new TFile(file1);
  if (f1->IsZombie())
    { 
      cout << "Sorry dude, file with run number \"" << nrun << "\" could not be opened." << endl; 
      cout << "Check to make sure that it exists and is really a root file." << endl;
      exit(1);
    }
  else
    {
      TArrayI prescales = Run_Data->GetParameters()->GetPrescales();
      
      cout << "Run: " << nrun << endl;
      for( int i = 0; i<prescales->GetSize(); ++i ) 
	cout << "prescale[" << i+1 << "] = " << prescales[i] << endl; 
    }

}

