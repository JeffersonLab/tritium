void plot(){

  TCanvas *cc = new TCanvas("cc","Trigger 1",1200,800);
  cc->Divide(1,2);cc->cd(1);
  
  ifstream file1("L_Scin_Eff_a.txt");
  const int NUMR = 2000;
  float rnum[NUMR] ; 
  float Teff[NUMR];

  int i=0;
  while(!(file1.eof())){
    file1 >> rnum[i] >> Teff[i];
    i++;
  }
 
 const int TTT1 = i;
  gStyle->SetOptStat(0);

  TH2F *hh = new TH2F("hh","",820,3500,4320,820,99.0,100.1);
  hh->GetXaxis()->SetTitle("Run Number");
  hh->GetXaxis()->CenterTitle();
  hh->GetYaxis()->SetTitle("HRS-L T3 Trigger Efficiency (%)");
  hh->GetYaxis()->CenterTitle();
  hh->Draw();
  
  TGraph *t3 = new TGraph(TTT1,rnum,Teff);
  t3->SetMarkerStyle(20);
  t3->SetMarkerColor(2);
  t3->Draw("P");
 
  cc->cd(2);
  ifstream file2("R_Scin_Eff_a.txt");
  const int NUMR2 = 2000;
  float rnum2[NUMR2] ; 
  float Teff2[NUMR2];

  int j=0;
  while(!(file2.eof())){
    file2 >> rnum2[j] >> Teff2[j];
    j++;
  }
 
  const int TTT2 = j;
  gStyle->SetOptStat(0);

  TH2F *hf = new TH2F("hf","",820,3500,4320,820,99.0,100.1);
  hf->GetXaxis()->SetTitle("Run Number");
  hf->GetXaxis()->CenterTitle();
  hf->GetYaxis()->SetTitle("HRS-R T1 Trigger Efficiency (%)");
  hf->GetYaxis()->CenterTitle();
  hf->Draw();
  
  TGraph *t1 = new TGraph(TTT2,rnum2,Teff2);
  t1->SetMarkerStyle(20);
  t1->SetMarkerColor(4);
  t1->Draw("P");
 
}
