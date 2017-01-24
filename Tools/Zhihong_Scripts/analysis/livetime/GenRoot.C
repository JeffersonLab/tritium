{
  TString infile;
  cerr<<"What is the input file?  "; cin >> infile;
  ifstream inf(infile.Data());
  TString com; com.ReadLine(inf);

  TString rootfile = infile.ReplaceAll(".table",".root");

  Double_t lt_avg, lt_avg_err, lt, lt_err, rate, Imax, Imean;
  int runno;

  TFile *f1 = new TFile(rootfile,"recreate");
  TTree* T = new TTree("T","A New Tree");

  T->Branch("lt_avg", &lt_avg, "lt_avg/D");
  T->Branch("lt_avg_err", &lt_avg_err, "lt_avg_err/D");
  T->Branch("lt", &lt , "lt/D");
  T->Branch("lt_err", &lt_err, "lt_err/D");
  T->Branch("Imax", &Imax, "Imax/D");
  T->Branch("rate", &rate , "rate/D");
  T->Branch("Imean", &Imean, "Imean/D");
  T->Branch("runno", &runno, "runno/I");

  while(!(inf.eof())){
    inf >> runno >> lt_avg >> lt_avg_err >> lt >> lt_err >> rate >> Imax >> Imean;
    T->Fill();
  }

  f1->cd(); T->Write(); f1->Close();
  inf.close();
}
