//////////////////////////////////////
//     F(y) Parameters Fitting      //
//     TMinuit Minimization         //
//    ---  Zhihong Ye 12/24/2012    //
//////////////////////////////////////
#include "fynk_minuit2.hh"
//Maximium Value of Y to fit:
Double_t Ymax = 0.09;
double yaxis_min = 1e-08;
//TString Target_Name =" ";
XEM_TGT* Target_Old = new XEM_TGT();  
XEM_TGT* Target_New = new XEM_TGT();  
int Angle = 0;

//Scaling Function F(y), from N. Fomin's Thesis (143P)
/*DoGraphFit{{{*/
//void Fy_Plot()
  int main()
{
  cerr<<"Target (H2,He3,He4,C12) = "; cin >> Target_Name;

  gStyle->SetOptStat(0);
  int n = 0; 
  TCanvas *c1=new TCanvas("c1","",1400,800);
  c1->cd(); gPad->SetLogy(1);


  TLegend *ll = new TLegend(0.6,0.3,0.85, 0.55,Form("%s Target",Target_Name.Data()));
  ll->SetTextSize(0.02);
  
  Angle = 18;
  LoadData();
  n=vY.size();

  if (yMin-0.1>yaxis_min)
    yaxis_min = yMin-0.1;
   TH2F *h1 = new TH2F("h1",Form("%s: F(y) Fitting",Target_Name.Data()),1000, xMin-0.1,xMax+0.1,1000,1e-6,20.);
  h1->SetXTitle("y");
  h1->SetYTitle("F(y)");
  h1->Draw();

  TGraphErrors *gr18=new TGraphErrors(n,&vY[0],&vFy[0],&vZero[0],&vFyErr[0]);
  gr18->SetMarkerStyle(24);
  gr18->SetMarkerSize(1.0);
  gr18->SetMarkerColor(2);
  gr18->Draw("P");
  ll->AddEntry(gr18,"E02-019 Data at 18^{o}","p");

  Angle = 22;
  LoadData();
  n=vY.size();
  TGraphErrors *gr22=new TGraphErrors(n,&vY[0],&vFy[0],&vZero[0],&vFyErr[0]);
  gr22->SetMarkerStyle(25);
  gr22->SetMarkerColor(4);
  gr22->SetMarkerSize(1.0);
  gr22->Draw("P");
  ll->AddEntry(gr22,"E02-019 Data at 22^{o}","p");

  Angle = 26;
  LoadData();
  n=vY.size();
  TGraphErrors *gr26=new TGraphErrors(n,&vY[0],&vFy[0],&vZero[0],&vFyErr[0]);
  gr26->SetMarkerStyle(26);
  gr26->SetMarkerSize(1.0);
  gr26->SetMarkerColor(6);
  gr26->Draw("P");
  ll->AddEntry(gr26,"E02-019 Data at 26^{o}","p");

  Angle = 32;
  LoadData();
  n=vY.size();
  TGraphErrors *gr32=new TGraphErrors(n,&vY[0],&vFy[0],&vZero[0],&vFyErr[0]);
  gr32->SetMarkerStyle(27);
  gr32->SetMarkerColor(7);
  gr32->Draw("P");
  ll->AddEntry(gr32,"E02-019 Data at 32^{o}","p");

  Angle = 40;
  LoadData();
  n=vY.size();
  TGraphErrors *gr40=new TGraphErrors(n,&vY[0],&vFy[0],&vZero[0],&vFyErr[0]);
  gr40->SetMarkerStyle(28);
  gr40->SetMarkerColor(8);
  gr40->SetMarkerSize(1.0);
  gr40->Draw("P");
  ll->AddEntry(gr40,"E02-019 Data at 40^{o}","p");

  Angle = 50;
  LoadData();
  n=vY.size();
  TGraphErrors *gr50=new TGraphErrors(n,&vY[0],&vFy[0],&vZero[0],&vFyErr[0]);
  gr50->SetMarkerStyle(30);
  gr50->SetMarkerSize(1.0);
  gr50->SetMarkerColor(9);
  gr50->Draw("P");
  ll->AddEntry(gr50,"E02-019 Data at 50^{o}","p");


  TF1 *thisF1=new TF1("thisF1",func,xMin-0.001,xMax+0.001,5);
  thisF1->SetParameters(Target_Old->f0,Target_Old->B, Target_Old->alpha, Target_Old->a, Target_Old->b);
  thisF1->SetParNames("f0","B","alpha","a","b");
  thisF1->SetLineColor(2);
  thisF1->Draw("same");
  ll->AddEntry(thisF1,"F(y) Function with Nadia's Fit","l");

  TF1 *thisF2=new TF1("thisF2",func,xMin-0.001,xMax+0.001,5);
  thisF2->SetParameters(Target_New->f0,Target_New->B, Target_New->alpha, Target_New->a, Target_New->b);
  thisF2->SetParNames("f0","B","alpha","a","b");
  thisF2->SetLineColor(4);
  thisF2->Draw("same");
  ll->AddEntry(thisF2,"F(y) Function with Zhihong's Fit","l");

  ll->Draw();

  c1->Print(Form("%s.png",Target_Name.Data()));
  c1->Print(Form("%s.pdf",Target_Name.Data()));
  
  // c1->Print(Form("%s_%d.png",Target_Name.Data(),Angle));
  // c1->Print(Form("%s_%d.pdf",Target_Name.Data(),Angle));


  delete gr18;  
  delete gr22;  
  delete gr26;  
  delete gr32;  
  delete gr40;  
  delete gr50;  
  delete h1;
  delete thisF1;
  delete thisF2;
  delete c1;
}
/*}}}*/

/*func{{{*/
Double_t func(Double_t *x, Double_t *par)
{
  //Set Parameters
  XEM_TGT* kTarget = new XEM_TGT();
  kTarget->GetValueAZ(Target->A,Target->Z,1);
  Double_t kY = x[0];
  kTarget->f0 = par[0]; 
  kTarget->B = par[1]; 
  kTarget->alpha = par[2]; 
  kTarget->a = par[3]; 
  kTarget->b = par[4]; 
  
  double  kFy_Fit = gGet_Fy(kY, kTarget);
  delete kTarget;
  return kFy_Fit;
 }
/*}}}*/

/*LoadData{{{*/
void LoadData(){
  TString Root_File = Form("%s_XEMC.root",Target_Name.Data());
  TChain *T = new TChain("T");
  T->Add(Root_File);
  int nTotal = T->GetEntries();

  Double_t xbj,Ep,Theta, xs_rad, xs_born,xs_dis,xs_qe, nXS_Rad, nXS_Born,nXS_CC,nFrac_Err;
  Int_t A,Z;

  T->SetBranchAddress("xbj", &xbj);
  T->SetBranchAddress("Ep",  &Ep);
  T->SetBranchAddress("Theta", &Theta);
  T->SetBranchAddress("xs_rad", &xs_rad);
  T->SetBranchAddress("xs_born",&xs_born);
  T->SetBranchAddress("xs_dis", &xs_dis);
  T->SetBranchAddress("xs_qe", &xs_qe);
  T->SetBranchAddress("nXS_Rad", &nXS_Rad);
  T->SetBranchAddress("nXS_Born", &nXS_Born);
  T->SetBranchAddress("nXS_CC", &nXS_CC);
  T->SetBranchAddress("nFrac_Err", &nFrac_Err);
  T->SetBranchAddress("A", &A);
  T->SetBranchAddress("Z", &Z);

  Double_t nXS_QE;
  XEM_VAR2* aFy = new XEM_VAR2();
  Double_t  aFy_EX = 0.0, aFy_EX_Err = 0.0, aY = 0.0;
  vEp.clear(); vTheta.clear(); vXbj.clear(); vZero.clear();
  vY.clear();  vXS.clear(); vXSErr.clear(); vFy.clear(); vFyErr.clear(); 
  for(int i=0;i<nTotal; i++){
    T->GetEntry(i);
    if(i==0){
      Target->GetValueAZ(A,Z,1);
      Target->Print();
      
      cerr<<Form("Initializing Old parameters for target %s, A=%d, Z=%d",Target_Name.Data(),A,Z)<<endl;
      Target_Old->GetValueAZ(A,Z,2);
      Target_Old->Print();

      cerr<<Form("Initializing New parameters for target %s, A=%d, Z=%d",Target_Name.Data(),A,Z)<<endl;
      Target_New->GetValueAZ(A,Z,1);
      Target_New->Print();

    }
    if(abs(Angle-Theta)<1e-5 || (Angle-0.0)<1e-5){
      //F(y)
      nXS_QE = (nXS_Born - xs_dis);
      //nXS_QE = xs_qe;
      //nXS_QE = (xs_born - xs_dis);
      aFy = gCal_Sig2Fy(E0, Ep, Theta * pi /180., nXS_QE, Target);
      aY = aFy->First;
      aFy_EX = aFy->Second;
      aFy_EX_Err =  aFy_EX * nFrac_Err;
      
      //    if(aY>-1e4 &&aY<=Ymax+0.01 ){
      //    if(aY>-1e4){
      if(aY>-2.0 && aY < 0.2){
	vEp.push_back(Ep);
	vTheta.push_back(Theta);
	vXbj.push_back(xbj);
	vZero.push_back(0.0);
	
	vXS.push_back(nXS_QE);
	vXSErr.push_back(nXS_QE*nFrac_Err);
	
	vFy.push_back(aFy_EX);
	vFyErr.push_back(aFy_EX_Err);
	vY.push_back(aY);
	
	if(xMin >= aY)
	  xMin = aY;
	if(xMax <=aY)
	  xMax = aY;
	if(yMin >= aFy_EX)
	  yMin = aFy_EX;
	if(yMax <= aFy_EX)
	  yMax = aFy_EX;
      }
    }
  }
  yMin = 1e-30;

  delete aFy;
  T->Delete();
}
 /*}}}*/

/*LoadDataAscii{{{*/
void LoadDataAscii(){
  TString Ascii_File = Form("data/%s.dat",Target_Name.Data());
  ifstream infile(Ascii_File);

  Int_t A,Z;

 /*Get Target Info{{{*/
  if(Target_Name == "H2"){
    A = 2; Z = 1;  
  }
  else if(Target_Name == "He3"){
    A = 3; Z = 2;  
  }
  else if(Target_Name == "He4"){
    A = 4; Z = 2; 
  }
  else if(Target_Name == "Be"){
    A = 9; Z = 4; 
  }
  else if(Target_Name == "C12"){
    A = 12; Z = 6; 
  }
  else if(Target_Name == "Al"){
    A = 27; Z = 13; 
  }
  else if(Target_Name == "Ca40"){
    A = 40; Z = 20;
  }
  else if(Target_Name == "Ca48"){
    A = 48; Z = 20; 
  }
  else if(Target_Name == "Cu"){
    A = 64; Z = 29;
  }
  else if(Target_Name == "Au"){
    A = 197; Z = 79; 
  }
  else
    cerr<<" **** ERROR, I don't know the target!!!"<<endl;
  /*}}}*/
  cerr<<Form("Initializing parameters for target %s, A=%d, Z=%d",
	     Target_Name.Data(),A,Z)<<endl;
  Target->GetValueAZ(A,Z,1);
  Target->Print();
  
  Double_t nA, nTheta, nP0, nEp, nXS_Rad, nXS_Born,nXS_CC,nFrac_Err, nXS_QE;  
  Double_t xs_dis = -1.,xs_qe = -1.,xs_born = -1.,xbj,Q2;
  XEM_VAR2* aFy = new XEM_VAR2();
  Double_t  aFy_EX = 0.0, aFy_EX_Err = 0.0, aY = 0.0;
  int i=0;
  vEp.clear(); vTheta.clear(); vXbj.clear(); vZero.clear();
  vY.clear();  vXS.clear(); vXSErr.clear(); vFy.clear(); vFyErr.clear(); 
  while(!(infile.eof())){
    infile >> nA >> nTheta >> nP0 >> nEp >> nXS_Rad >> nXS_Born >> nXS_CC >> nFrac_Err;
    
    Q2 = 4.0 * E0 * nEp * pow(sin(nTheta*3.1415926/180./2.0),2);
    xbj = Q2/2.0/0.938272/(E0-nEp);

    xs_dis = XEMC_Born(E0*GeVToMeV, nEp*GeVToMeV, nTheta, A, Z, 3); //3->For DIS_F1F2IN09
    //   xs_dis = XEMC_Born(E0*GeVToMeV, nEp*GeVToMeV, nTheta, A, Z, 4); //4->For DIS_F1F2IN06  
    //xs_qe = XEMC_Born(E0*GeVToMeV, nEp*GeVToMeV, nTheta, A, Z, 2); //2->For QE_XEM  
    // xs_born = XEMC_Born(E0*GeVToMeV, nEp*GeVToMeV, nTheta, A, Z, 1); //1->For DIS_F1F2IN09 + QE
    // nXS_QE = xs_qe;
    //nXS_QE = xs_born - xs_dis;
    nXS_QE = nXS_Born - xs_dis;

    if(abs(Angle-nTheta)<1e-5 || (Angle-0.0)<1e-5){
      //  	   	    if(1){
    //F(y)
    aFy = gCal_Sig2Fy(E0, nEp, nTheta * pi /180., nXS_QE, Target);
    aY = aFy->First;
    aFy_EX = aFy->Second;
    aFy_EX_Err =  aFy_EX * nFrac_Err;

    //    if(aY>-1e4 &&aY<=Ymax+0.01 ){
    if(aY>-1e4 && aY>=-0.8){
    vEp.push_back(nEp);
    vTheta.push_back(nTheta);
    vXbj.push_back(xbj);
    vZero.push_back(0.0);
    nXS_QE = (nXS_Born - xs_dis);
    vXS.push_back(nXS_QE);
    vXSErr.push_back(nXS_QE*nFrac_Err);

    vFy.push_back(aFy_EX);
    vFyErr.push_back(aFy_EX_Err);
    vY.push_back(aY);

    if(xMin >= aY)
      xMin = aY;
    if(xMax <=aY)
      xMax = aY;
    if(yMin >= aFy_EX)
      yMin = aFy_EX;
    if(yMax <= aFy_EX)
      yMax = aFy_EX;
    }
    }
    i++;
  }
  delete aFy;
 }
/*}}}*/

/*output{{{*/
//Output the results into a file
void output(Double_t *outPar)
{ 
  TString name =  Form("%s_%d_Par.dat",Target_Name.Data(),Angle);

	ofstream out(name,ios::out);

	cerr<<endl<<"The Fitting Result is: "<<endl;

	cerr <<setiosflags(ios::fixed)
		<<"f0    = "
		<<setw(14)<<setprecision(8)<<outPar[0]<<endl;
	cerr <<setiosflags(ios::fixed)
		<<"bigB  = "
		<<setw(14)<<setprecision(8)<<outPar[1]<<endl;
	cerr <<setiosflags(ios::fixed)
		<<"alpha = "
		<<setw(14)<<setprecision(8)<<outPar[2]<<endl;

	cerr <<setiosflags(ios::fixed)
		<<"a     = "
		<<setw(14)<<setprecision(8)<<outPar[3]<<endl;
	cerr <<setiosflags(ios::fixed)
		<<"b     = "
		<<setw(14)<<setprecision(8)<<outPar[4]<<endl;

	out <<"======================================================================================="<<endl;
	out <<setiosflags(ios::fixed);
	out <<"   F(y) = (f0-bigB)*(alpha*alpha*Exp(-a*a*y*y)/(alpha*alpha+y*y)) + bigB*Exp(-b*|y|) "<<endl;
	out <<"----------------------------------------------------------------------------------------"<<endl;
	out <<" Where for Angle =  "<<Angle<<": "<<endl;
	out <<setiosflags(ios::fixed)
		<<"         f0    = "
		<<setw(14)<<setprecision(8)<<outPar[0]<<endl;
	out <<setiosflags(ios::fixed)
		<<"         bigB  = "
		<<setw(14)<<setprecision(8)<<outPar[1]<<endl;
	out <<setiosflags(ios::fixed)
		<<"         alpha = "
		<<setw(14)<<setprecision(8)<<outPar[2]<<endl;

	out <<setiosflags(ios::fixed)
		<<"         a     = "
		<<setw(14)<<setprecision(8)<<outPar[3]<<endl;
	out <<setiosflags(ios::fixed)
		<<"         b     = "
		<<setw(14)<<setprecision(8)<<outPar[4]<<endl;
	out <<" Note: Some values have to be scaled: alpha = alpha/1000.0, y = y*1000.0."<<endl;
	out <<"========================================================================================"<<endl;
}
/*}}}*/
