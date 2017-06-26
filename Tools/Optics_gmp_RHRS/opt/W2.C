void W2(){
  Double_t mean=0;
  Double_t sigma=0,ymax;
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(3,2);

  c1->cd(1);
  //TFile *f0 = new TFile("/Users/tong/Tritium_analysis/optics/rootfiles/right_gmp_22006.root");
  TFile *f0 = new TFile("/Users/tong/Loop120116/right_gmp_22006.root");
  T->Draw("EKRxe.W2>>h1(100,0.8,1)");
  TF1 *h1fit = new TF1("h1fit","gaus",0.87+0.005,0.89+0.005);
  h1->Fit("h1fit","R");
  mean = h1fit->GetParameter(1);
  sigma = h1fit->GetParameter(2);
  cout<<"mean:"<<mean<<endl;
  TPaveText *t1 = new TPaveText(0.52,0.5,0.9,0.7,"NDC");
  t1->AddText(Form("P_{gaus} %.4f",mean));
  t1->AddText(Form("#sigma = %.4f",sigma));
  t1->SetFillColor(0);
  t1->SetShadowColor(0);
  t1->Draw();
  c1->cd(1)->Update();
  ymax=c1->cd(1)->GetUymax();
  TLine *line1 = new TLine(mean-0.01,0,mean-0.01,ymax);
  TLine *line2 = new TLine(mean+0.01,0,mean+0.01,ymax);
  line1->Draw("same");
  line2->Draw("same");

  c1->cd(2);
  //TFile *f0 = new TFile("/Users/tong/Tritium_analysis/optics/rootfiles/right_gmp_22007.root");
  TFile *f0 = new TFile("/Users/tong/Loop120116/right_gmp_22007.root");
  T->Draw("EKRxe.W2>>h2(100,0.8,1)");
  TF1 *h2fit = new TF1("h2fit","gaus",0.87+0.005,0.89+0.005);
  h2->Fit("h2fit","R");
  mean = h2fit->GetParameter(1);
  sigma = h2fit->GetParameter(2);
  cout<<"mean:"<<mean<<endl;
  TPaveText *t2 = new TPaveText(0.52,0.5,0.9,0.7,"NDC");
  t2->AddText(Form("P_{gaus} %.4f",mean));
  t2->AddText(Form("#sigma = %.4f",sigma));
  t2->SetFillColor(0);
  t2->SetShadowColor(0);
  t2->Draw();
  c1->cd(2)->Update();
  ymax=c1->cd(2)->GetUymax();
  TLine *line1 = new TLine(mean-0.01,0,mean-0.01,ymax);
  TLine *line2 = new TLine(mean+0.01,0,mean+0.01,ymax);
  line1->Draw("same");
  line2->Draw("same");

  c1->cd(3);
  //TFile *f0 = new TFile("/Users/tong/Tritium_analysis/optics/rootfiles/right_gmp_22008.root");
  TFile *f0 = new TFile("/Users/tong/Loop120116/right_gmp_22008.root");
  T->Draw("EKRxe.W2>>h3(100,0.8,1)");
  TF1 *h3fit = new TF1("h3fit","gaus",0.87+0.005,0.89+0.005);
  h3->Fit("h3fit","R");
  mean = h3fit->GetParameter(1);
  sigma = h3fit->GetParameter(2);
  cout<<"mean:"<<mean<<endl;
  TPaveText *t3 = new TPaveText(0.52,0.5,0.9,0.7,"NDC");
  t3->AddText(Form("P_{gaus} %.4f",mean));
  t3->AddText(Form("#sigma = %.4f",sigma));
  t3->SetFillColor(0);
  t3->SetShadowColor(0);
  t3->Draw();
  c1->cd(3)->Update();
  ymax=c1->cd(3)->GetUymax();
  TLine *line1 = new TLine(mean-0.01,0,mean-0.01,ymax);
  TLine *line2 = new TLine(mean+0.01,0,mean+0.01,ymax);
  line1->Draw("same");
  line2->Draw("same");

  c1->cd(4);
 //TFile *f0 = new TFile("/Users/tong/Tritium_analysis/optics/rootfiles/right_gmp_22011.root");
  TFile *f0 = new TFile("/Users/tong/Loop120116/right_gmp_22011.root");
  T->Draw("EKRxe.W2>>h4(100,0.8,1)");
  TF1 *h4fit = new TF1("h4fit","gaus",0.87+0.005,0.89+0.005);
  h4->Fit("h4fit","R");
  mean = h4fit->GetParameter(1);
  sigma = h4fit->GetParameter(2);
  cout<<"mean:"<<mean<<endl;
  TPaveText *t4 = new TPaveText(0.52,0.5,0.9,0.7,"NDC");
  t4->AddText(Form("P_{gaus} %.4f",mean));
  t4->AddText(Form("#sigma = %.4f",sigma));
  t4->SetFillColor(0);
  t4->SetShadowColor(0);
  t4->Draw();
  c1->cd(4)->Update();
  ymax=c1->cd(4)->GetUymax();
  TLine *line1 = new TLine(mean-0.01,0,mean-0.01,ymax);
  TLine *line2 = new TLine(mean+0.01,0,mean+0.01,ymax);
  line1->Draw("same");
  line2->Draw("same");

  c1->cd(5);
  //TFile *f0 = new TFile("/Users/tong/Tritium_analysis/optics/rootfiles/right_gmp_22013.root");
  TFile *f0 = new TFile("/Users/tong/Loop120116/right_gmp_22013.root");
  T->Draw("EKRxe.W2>>h5(100,0.8,1)");
  TF1 *h5fit = new TF1("h5fit","gaus",0.87+0.005,0.89+0.005);
  h5->Fit("h5fit","R");
  mean = h5fit->GetParameter(1);
  sigma = h5fit->GetParameter(2);
  cout<<"mean:"<<mean<<endl;
  TPaveText *t5 = new TPaveText(0.52,0.5,0.9,0.7,"NDC");
  t5->AddText(Form("P_{gaus} %.4f",mean));
  t5->AddText(Form("#sigma = %.4f",sigma));
  t5->SetFillColor(0);
  t5->SetShadowColor(0);
  t5->Draw();
  c1->cd(5)->Update();
  ymax=c1->cd(5)->GetUymax();
  TLine *line1 = new TLine(mean-0.01,0,mean-0.01,ymax);
  TLine *line2 = new TLine(mean+0.01,0,mean+0.01,ymax);
  line1->Draw("same");
  line2->Draw("same");






}
