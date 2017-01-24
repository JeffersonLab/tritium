void plot(){

  gStyle->SetOptStat(0);
  TCanvas *cc = new TCanvas("cc","cc",1200,800);
  TH2F *hh = new TH2F("hh","",810,3500,4310,200,90.0,101.0);

 ifstream infile("L_Det_Eff.dat");

  Float_t run_number=0, comment;

  Float_t cer_t1, cer_t6,calo_t1,calo_t6;
  Float_t cert1[700], cert6[700],calot1[700],calot6[700];
  Float_t runn[700];

  int i =0;
  while(!infile.eof()){
    infile >> run_number >> cer_t1 >> cer_t6 >> calo_t1 >> calo_t6 >> comment;

    

    runn[i] = run_number;
    cert1[i] = cer_t1;
    cert6[i] = cer_t6;
    calot1[i] = calo_t1;
    calot6[i] = calo_t6;
    
    i++;
}
  const int count = i;

 
  hh->GetXaxis()->SetTitle("HRS-R Cerenkov Detection Efficiencies");
  hh->GetXaxis()->CenterTitle();
  hh->GetYaxis()->SetTitle("Detection Efficiency (%)");
  hh->GetYaxis()->CenterTitle();
  hh->Draw();

  TGraph *gr1 = new TGraph(count,runn,cert1);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(4);
  gr1->Draw("P");
  
  TGraph *gr2 = new TGraph(count,runn,cert6);
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerColor(6);
  gr2->Draw("P");
  
  TLegend *l1 = new TLegend(0.4,0.4,0.6,0.6);
  l1->AddEntry(gr1,"Trigger T1","p");
  l1->AddEntry(gr2,"Trigger T6","p");
  l1->Draw();

  TCanvas *cf = new TCanvas("cf","cf",1200,800);
  TH2F *hf = new TH2F("hf","",810,3500,4310,200,80.0,110.0);
  hf->GetXaxis()->SetTitle("HRS-R Calorimeter Detection Efficiencies");
  hf->GetXaxis()->CenterTitle();
  hf->GetYaxis()->SetTitle("Detection Efficiency (%)");
  hf->GetYaxis()->CenterTitle();
  hf->Draw();

  TGraph *gr3 = new TGraph(count,runn,calot1);
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerColor(4);
  gr3->Draw("P");
  
  TGraph *gr4 = new TGraph(count,runn,calot6);
  gr4->SetMarkerStyle(21);
  gr4->SetMarkerColor(6);
  gr4->Draw("P");
  
  TLegend *l2 = new TLegend(0.4,0.4,0.6,0.6);
  l2->AddEntry(gr3,"Trigger T1","p");
  l2->AddEntry(gr4,"Trigger T6","p");
  l2->Draw();


}
