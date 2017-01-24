/////////////////////////////////////
//   Dummy Contaimination Study    //
//       For X>2 Experiment        //
// --  Zhihong Ye, 09/30/2011      //
/////////////////////////////////////

#include <iostream>
#include <iomanip>
using namespace std;

Double_t dgaus(Double_t* x, Double_t* par);

void Get_Dummy_Contaim(){

  /* Cuts {{{*/
  TString GeneralCut = "L.tr.n==1 && abs(L.tr.x)<0.9 && abs(L.tr.y)<0.5 && abs(L.tr.th)<0.8 && abs(L.tr.ph)<0.5&&L.tr.p>2.0&&L.tr.p<4.0";
  GeneralCut += "&& abs(exL.ph)<0.02 && abs(exL.th)<0.04 && abs(exL.dp)<0.06 && abs(rpl.z)<0.2";
  TString TriggerCut = "((DBB.evtypebits>>3)&1)&&DBB.edtpr[0]==0";
  TString electrons = "(L.prl1.e+L.prl2.e)/L.tr.p/1000.0>=0.5&&L.prl2.e>=100.0&&L.cer.asum_c>=50.0";
  TString cute1 = TriggerCut+"&&" + GeneralCut + "&&"+electrons;
  /*}}}*/

  /*BCM for L-HRS{{{*/
  Double_t alpha_u = 2070, offset_u = 400;
  Double_t alpha_d = 1263, offset_d = 20;
  /*}}}*/

  Double_t i_u1,i_d1,q_u1,q_d1, time,charge,diff,livetime,count1,count2,yield1,yield2;
  Int_t PS1=1, PS2=1.,PS3=1.,PS4=1.,PS6=1,PS7=1.;
  Double_t Left_LiveTime = -1.0;  Double_t Right_LiveTime = -1.0;
  Double_t Right_Raw_Trigger, Left_Raw_Trigger;
  Int_t count,i,j,k,curr,bcm_u1c,bcm_d1c;
  Double_t Yield_Dum,Err_Dum,Yield_He3,Err_He3;
  Double_t Yield_He4,Err_He4,Yield_H2,Err_H2;
 
  Double_t Density_Dum = 2.81; //g/cc;
  Double_t RP_Cut = 0.15;
  
  Bool_t Dummy = kFALSE; Bool_t He3 = kTRUE;
  Bool_t He4 = kFALSE; Bool_t H2 = kFALSE;

  ///////////////////////
  //  Dummy 
  ///////////////////////
  
  if(Dummy){
    cerr<<endl<< "Woring on Dummy20cm ..."<<endl;

    ofstream f_dum("Yield_Dum.txt");	
    Double_t Thick_Dum_U = 0.1581; //cm
    Double_t Thick_Dum_D = 0.1585; //cm

    TFile *dummy = new TFile("Dummy20cm_L_4100.root","r");
    TTree *T_Dum = (TTree*)gDirectory->Get("T");
 
    TArrayI prescales_dum = Run_Data->GetParameters()->GetPrescales();
    PS1 = prescales_dum[0]; PS2 = prescales_dum[1]; PS3 = prescales_dum[2];
    PS4 = prescales_dum[3]; PS6 = prescales_dum[5]; PS7 = prescales_dum[6];
  
    TChain *cT_Dum = new TChain("T");
    cT_Dum->Add("Dummy20cm_L_4100.root");
    cT_Dum->Add("Dummy20cm_L_4100_1.root");
  
    /*Get Livetime{{{*/
    Int_t DAQ_T3 = (Int_t) cT_Dum->GetEntries("(DBB.evtypebits>>3)&1");
    Left_Raw_Trigger = cT_Dum->GetMaximum("left_t3c") - cT_Dum->GetMinimum("left_t3c");
    Left_LiveTime  = (Double_t) 1.0*PS3*DAQ_T3 /Left_Raw_Trigger;
    /*}}}*/

    /*Get Charge {{{*/
    bcm_u1c = cT_Dum->GetMaximum("left_bcm_u1c") - cT_Dum->GetMinimum("left_bcm_u1c");
    bcm_d1c = cT_Dum->GetMaximum("left_bcm_d1c") - cT_Dum->GetMinimum("left_bcm_d1c");
    time =  (cT_Dum->GetMaximum("left_clkcount") - cT_Dum->GetMinimum("left_clkcount"))/1024;
    i_u1 = ((bcm_u1c/time)-offset_u)/alpha_u;
    i_d1 = ((bcm_d1c/time)-offset_d)/alpha_d;
    q_u1 = i_u1*time;
    q_d1 = i_d1*time;
    charge = q_u1;
    diff = 100*(q_u1 - q_d1)/((q_u1+q_d1)/2.0);
    /*}}}*/
  
    for(i=0;i<80;i++){
      RP_Cut = 0.04+ i*0.001;
      count1 = PS3*cT_Dum->GetEntries(TCut(cute1+Form("&&(rpl.z+0.00664)<0.0&&abs(rpl.z+0.00664)<%f",RP_Cut)));
      yield1 = count1/charge/Left_LiveTime/Thick_Dum_U;
      count2 = PS3*cT_Dum->GetEntries(TCut(cute1+Form("&&(rpl.z+0.00664)>0.0&&abs(rpl.z+0.00664)<%f",RP_Cut)));
      yield2 = count2/charge/Left_LiveTime/Thick_Dum_D;
    
      Yield_Dum = yield1 + yield2;
      Err_Dum = 1.0/sqrt(count1+count2)/charge/Left_LiveTime;
    
      f_dum <<setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<< RP_Cut <<"   ";
      f_dum <<setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<<Yield_Dum<<"   ";
      f_dum <<setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<<Err_Dum<<endl;
    }    

    cerr<<"=============================================="<<endl;
    cerr<<" -- The L-HRS Live Time is: "<<Left_LiveTime<<endl;
    cerr<< " --- Charge = "<< charge << " Events = "<< count <<endl 
        << " --- Yield = "<< Yield_Dum <<", Error = "<<Err_Dum <<endl;
    cerr<<"=============================================="<<endl;

    cT_Dum->Delete();  T_Dum->Delete(); dummy->Close();
  }

  /////////////////////////
  // He3
  /////////////////////////
 
  if(He3){
    cerr <<endl<< "Woring on He3 ..."<<endl;
    ofstream f_he3("Yield_He3.txt");

    Double_t Thick_He3_U = 0.0274; //cm
    Double_t Thick_He3_D = 0.0328; //cm
  
    TFile *he3 = new TFile("He3_L_4085.root","r");
    TTree *T_He3 = (TTree*)gDirectory->Get("T");
 
    TArrayI prescales_he3 = Run_Data->GetParameters()->GetPrescales();
    PS1 = prescales_he3[0]; PS2 = prescales_he3[1]; PS3 = prescales_he3[2];
    PS4 = prescales_he3[3]; PS6 = prescales_he3[5]; PS7 = prescales_he3[6];
  
    TChain *cT_He3 = new TChain("T");
    cT_He3->Add("He3_L_4085.root");
    cT_He3->Add("He3_L_4085_1.root");
    cT_He3->Add("He3_L_4085_2.root");
  
    /*Get Livetime{{{*/
    Int_t DAQ_T3 = (Int_t) cT_He3->GetEntries("(DBB.evtypebits>>3)&1");
    Left_Raw_Trigger = cT_He3->GetMaximum("left_t3c") - cT_He3->GetMinimum("left_t3c");
    Left_LiveTime  = (Double_t) 1.0*PS3*DAQ_T3 /Left_Raw_Trigger;
    /*}}}*/

    /*Get Charge {{{*/
    bcm_u1c = cT_He3->GetMaximum("left_bcm_u1c") - cT_He3->GetMinimum("left_bcm_u1c");
    bcm_d1c = cT_He3->GetMaximum("left_bcm_d1c") - cT_He3->GetMinimum("left_bcm_d1c");
    time =  (cT_He3->GetMaximum("left_clkcount") - cT_He3->GetMinimum("left_clkcount"))/1024;
    i_u1 = ((bcm_u1c/time)-offset_u)/alpha_u;
    i_d1 = ((bcm_d1c/time)-offset_d)/alpha_d;
    q_u1 = i_u1*time;
    q_d1 = i_d1*time;
    charge = q_u1;
    diff = 100*(q_u1 - q_d1)/((q_u1+q_d1)/2.0);
    /*}}}*/
   
    for(i=0;i<80;i++){
      RP_Cut = 0.04+ i*0.001;
      count1 = PS3*cT_He3->GetEntries(TCut(cute1+Form("&&(rpl.z)<0.0&&abs(rpl.z)<%f",RP_Cut)));
      yield1 = count1/charge/Left_LiveTime/Thick_He3_U;
      count2 = PS3*cT_He3->GetEntries(TCut(cute1+Form("&&(rpl.z)>0.0&&abs(rpl.z)<%f",RP_Cut)));
      yield2 = count2/charge/Left_LiveTime/Thick_He3_D;
   
      Yield_He3 = yield1 +yield2;
      Err_He3 = 1.0/sqrt(count1+count2)/charge/Left_LiveTime;
    
      f_he3 <<setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<< RP_Cut <<"   ";
      f_he3 <<setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<<Yield_He3<<"   ";
      f_he3 <<setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<<Err_He3<<endl;
    }
    cerr<<"=============================================="<<endl;
    cerr<<" -- The L-HRS Live Time is: "<<Left_LiveTime<<endl;
    cerr<< " --- Charge = "<< charge << " Events = "<< count <<endl 
        << " --- Yield = "<< Yield_He3 <<", Error = "<<Err_He3 <<endl;
    cerr<<"=============================================="<<endl;
     
    cT_He3->Delete(); T_He3->Delete(); he3->Close();
  }

  /////////////////////////
  // He4
  /////////////////////////

  if(He4){
    cerr <<endl<< "Woring on He4 ..."<<endl;
    ofstream f_he4("Yield_He4.txt"); 
    
    Double_t Thick_He4_U = 0.0274; //cm
    Double_t Thick_He4_D = 0.0328; //cm
  
    TFile *he4 = new TFile("He4_L_3674.root","r");
    TTree *T_He4 = (TTree*)gDirectory->Get("T");
 
    TArrayI prescales_he4 = Run_Data->GetParameters()->GetPrescales();
    PS1 = prescales_he4[0]; PS2 = prescales_he4[1]; PS3 = prescales_he4[2];
    PS4 = prescales_he4[3]; PS6 = prescales_he4[5]; PS7 = prescales_he4[6];
  
    TChain *cT_He4 = new TChain("T");
    cT_He4->Add("He4_L_3674.root");
    cT_He4->Add("He4_L_3674_1.root");
    //  cT_He4->Add("He4_L_3674_2.root");
  
    /*Get Livetime{{{*/
    Int_t DAQ_T3 = (Int_t) cT_He4->GetEntries("(DBB.evtypebits>>3)&1");
    Left_Raw_Trigger = cT_He4->GetMaximum("left_t3c") - cT_He4->GetMinimum("left_t3c");
    Left_LiveTime  = (Double_t) 1.0*PS3*DAQ_T3 /Left_Raw_Trigger;
    /*}}}*/

    /*Get Charge {{{*/
    bcm_u1c = cT_He4->GetMaximum("left_bcm_u1c") - cT_He4->GetMinimum("left_bcm_u1c");
    bcm_d1c = cT_He4->GetMaximum("left_bcm_d1c") - cT_He4->GetMinimum("left_bcm_d1c");
    time =  (cT_He4->GetMaximum("left_clkcount") - cT_He4->GetMinimum("left_clkcount"))/1024;
    i_u1 = ((bcm_u1c/time)-offset_u)/alpha_u;
    i_d1 = ((bcm_d1c/time)-offset_d)/alpha_d;
    q_u1 = i_u1*time;
    q_d1 = i_d1*time;
    charge = q_u1;
    diff = 100*(q_u1 - q_d1)/((q_u1+q_d1)/2.0);
    /*}}}*/
   
    for(i=0;i<80;i++){
      RP_Cut = 0.04+ i*0.001;
      count1 = PS3*cT_He4->GetEntries(TCut(cute1+Form("&&(rpl.z)<0.0&&abs(rpl.z)<%f",RP_Cut)));
      yield1 = count1/charge/Left_LiveTime/Thick_He4_U;
      count2 = PS3*cT_He4->GetEntries(TCut(cute1+Form("&&(rpl.z)>0.0&&abs(rpl.z)<%f",RP_Cut)));
      yield2 = count2/charge/Left_LiveTime/Thick_He4_D;
   
      Yield_He4 = yield1 +yield2;
      Err_He4 = 1.0/sqrt(count1+count2)/charge/Left_LiveTime;
    
      f_he4 << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<< RP_Cut <<"   ";
      f_he4 << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<<Yield_He4<<"   ";
      f_he4 << setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<<Err_He4<<endl;
    }
    cerr<<"=============================================="<<endl;
    cerr<<" -- The L-HRS Live Time is: "<<Left_LiveTime<<endl;
    cerr<< " --- Charge = "<< charge << " Events = "<< count <<endl 
        << " --- Yield = "<< Yield_He4 <<", Error = "<<Err_He4 <<endl;
    cerr<<"=============================================="<<endl;
     
    cT_He4->Delete(); T_He4->Delete(); he4->Close();
  }

  /////////////////////////
  // H2
  /////////////////////////
  
  if(H2){
    cerr <<endl<< "Woring on H2 ..."<<endl;
    ofstream f_h2("Yield_H2.txt");
    
    Double_t Thick_H2_U = 0.0274; //cm
    Double_t Thick_H2_D = 0.0328; //cm
  
    TFile *h2 = new TFile("H2_L_3682.root","r");
    TTree *T_H2 = (TTree*)gDirectory->Get("T");
 
    TArrayI prescales_h2 = Run_Data->GetParameters()->GetPrescales();
    PS1 = prescales_h2[0]; PS2 = prescales_h2[1]; PS3 = prescales_h2[2];
    PS4 = prescales_h2[3]; PS6 = prescales_h2[5]; PS7 = prescales_h2[6];
  
    TChain *cT_H2 = new TChain("T");
    cT_H2->Add("H2_L_3682.root");
    cT_H2->Add("H2_L_3682_1.root");
    // cT_H2->Add("H2_L_3862_2.root");
  
    /*Get Livetime{{{*/
    Int_t DAQ_T3 = (Int_t) cT_H2->GetEntries("(DBB.evtypebits>>3)&1");
    Left_Raw_Trigger = cT_H2->GetMaximum("left_t3c") - cT_H2->GetMinimum("left_t3c");
    Left_LiveTime  = (Double_t) 1.0*PS3*DAQ_T3 /Left_Raw_Trigger;
    /*}}}*/

    /*Get Charge {{{*/
    bcm_u1c = cT_H2->GetMaximum("left_bcm_u1c") - cT_H2->GetMinimum("left_bcm_u1c");
    bcm_d1c = cT_H2->GetMaximum("left_bcm_d1c") - cT_H2->GetMinimum("left_bcm_d1c");
    time =  (cT_H2->GetMaximum("left_clkcount") - cT_H2->GetMinimum("left_clkcount"))/1024;
    i_u1 = ((bcm_u1c/time)-offset_u)/alpha_u;
    i_d1 = ((bcm_d1c/time)-offset_d)/alpha_d;
    q_u1 = i_u1*time;
    q_d1 = i_d1*time;
    charge = q_u1;
    diff = 100*(q_u1 - q_d1)/((q_u1+q_d1)/2.0);
    /*}}}*/
   
    for(i=0;i<80;i++){
      RP_Cut = 0.04+ i*0.001;
      count1 = PS3*cT_H2->GetEntries(TCut(cute1+Form("&&(rpl.z)<0.0&&abs(rpl.z)<%f",RP_Cut)));
      yield1 = count1/charge/Left_LiveTime/Thick_H2_U;
      count2 = PS3*cT_H2->GetEntries(TCut(cute1+Form("&&(rpl.z)>0.0&&abs(rpl.z)<%f",RP_Cut)));
      yield2 = count2/charge/Left_LiveTime/Thick_H2_D;
   
      Yield_H2 = yield1 +yield2;
      Err_H2 = 1.0/sqrt(count1+count2)/charge/Left_LiveTime;
    
      f_h2 <<setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<< RP_Cut <<"   ";
      f_h2 <<setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<<Yield_H2<<"   ";
      f_h2 <<setiosflags(ios::left) << setw(5) << setiosflags(ios::fixed) << setprecision(3)<<Err_H2<<endl;
    }
    cerr<<"=============================================="<<endl;
    cerr<<" -- The L-HRS Live Time is: "<<Left_LiveTime<<endl;
    cerr<< " --- Charge = "<< charge << " Events = "<< count <<endl 
        << " --- Yield = "<< Yield_H2 <<", Error = "<<Err_H2 <<endl;
    cerr<<"=============================================="<<endl;
     
    cT_H2->Delete(); T_H2->Delete(); h2->Close();
  }

}

void plot(){

  float yield_dum[80],err_dum[80],rpcut_dum[80],yield_he3[80],err_he3[80],rpcut_he3[80],ratio[80];
  float errc1[80]={10*0.0};
  ifstream infile1("Yield_Dum.txt");
 
  int  i=0;
  for(i=0;i<80;i++){
    infile1 >> rpcut_dum[i] >> yield_dum[i] >> err_dum[i];
    cerr <<"Cuts = "<<rpcut_dum[i]<<" Yield = "<<yield_dum[i]<<endl;
  }
  const int TTT1 = 80;  infile1.close();
 
  ifstream infile2("Yield_He3.txt");
  for(i=0;i<80;i++){
    infile2 >> rpcut_he3[i] >> yield_he3[i] >> err_he3[i];
    cerr <<"Cuts = "<<rpcut_he3[i]<<" Yield = "<<yield_he3[i]<<endl;

    ratio[i] = yield_dum[i]/yield_he3[i]*100.0;
  }
  const int TTT2 = 80;infile2.close();
 
  gStyle->SetOptStat(0);
  TCanvas *cc = new TCanvas("cc","cc",1200,800);
  TH2F *hh = new TH2F("hh","",200,0.055,0.125,200,0.0,40.0);
  hh->SetTitle("ReactPoint.Z Cut vs Yield(HRS-L)");
  hh->GetXaxis()->SetTitle("abs(ReactPoint.Z)");
  hh->GetXaxis()->CenterTitle();
  hh->GetYaxis()->SetTitle("Dummy Contaim(%)");
  hh->GetYaxis()->CenterTitle();
  hh->Draw();
  TGraph *gr1 = new TGraph(TTT1,rpcut_dum,ratio);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerColor(4);
  gr1->Draw("P");
  
}



Float_t dgaus(Double_t* x, Double_t* par){

  Double_t t1, t2, mean1,mean2, sigma1, sigma2;
  t1 = par[0];       t1 = par[3];
  mean1  = par[1];    mean2 = par[4];
  sigma1 = par[2];   sigma2 = par[5];
  

  Double_t dg = t1* TMath::Gaus(x[0],mean1,sigma1);
  dg += t2* TMath::Gaus(x[0],mean2,sigma2);
  //  dg += par[6]+par[7]*x[0];//+par[8]*x[0]*x[0];
  
  return dg;
}

Double_t pgaus(Double_t *x, Double_t *par) 
{
  
  Double_t gs = par[0]*TMath::Gaus(x[0],par[1],par[2]);
  Double_t pl = par[3]+par[4]*x[0]+par[5]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0];
  Double_t ep = par[7]*exp(-(x[0]-par[1])/par[8]);
  Double_t sum = gs+pl+ep;
  return(sum);
}

/*
  Double_t convolution(Double_t *x, Double_t *par) 
  {
  Double_t fland,sum,xlow,xupp,step,i,xx; 
  Double_t np,sc,invsq2pi;
  Double_t shift,mean;
  //  Fit parameters:
  //     par[0]= Weight of Total Fucntion
  //     par[1]= Tao of Exponential function
  //     par[2]= Sigma of Gaussian in Convolution Function 
  //     par[3]= Shift of Convolution Function Peak
  //     par[4]= Ratio of Gaussian vs Convlolution


  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  np= 200.0;            //  number of convolution steps
  sc= 8.0;              // convolution extends to +- sc Gaussian sigmas
  sum= 0.0;
  shift = -par[3];
  //  mean = 0.19;          //Mean value is from the fitting of prompt
  mean = par[4];
  x[0] = x[0]-mean;
  
  //Range of Convolution integral
  // Because the Range of Exponential Function is (0 ~ Inf)
  xlow = 0.0;
  xupp = x[0] + sc*par[2];
  step = (xupp - xlow)/np;
   
  for(i =1.0; i<=np/2;i++)
  {
  xx = xlow + (i-0.5)*step-shift;
  fland = TMath::Gaus(xx,x[0],par[2]);
  sum += fland*exp(-xx/par[1]);
      
  xx = xupp - (i-0.5)*step-shift;
  fland = TMath::Gaus(xx,x[0],par[2]);
  sum+= fland*exp(-xx/par[1]);
  }

  //  return (par[0]*sum*step*invsq2pi/(par[2]*par[1]*exp(par[3]/par[1])));
 
  // BackGround
  Double_t poly = par[5]+par[6]*x[0]+par[7]*x[0]*x[0]+par[8]*x[0]*x[0]*x[0];
  return (par[0]*sum*step*invsq2pi/(par[2]*par[1]*exp(par[3]/par[1]))+poly);
  }
*/



/*
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);c1->cd();
  TH1F *hdum_vz = new TH1F("hdum_vz", "Dummy Vertex Z", 300,-0.2,0.2);
  T1->Draw("rpl.z>>hdum_vz",TCut(cute1));
  
  Double_t par[9];
  // TF1 *g1  = new TF1("g1","gaus",-0.15,0.0);
  // TF1 *g2  = new TF1("g2","gaus", 0.0,0.15);
  // TF1 *g3  = new TF1("g3","pol2",-0.15,0.15);
  // TF1 *total = new TF1("total","gaus(0)+gaus(3)+pol2(6)",-0.15,0.15);
  //  TF1 *dgfit = new TF1("fit1",dgaus,-0.2,0.2,8);
  //  dgfit->SetParNames("Peak_L","Mean_L","Sigma_L","Peak_R","Mean_R","Sigma_R","pol1","pol2");
  //  dgfit->SetParameters(90000,-0.1,0.02,50000,0.1,0.02,10.0,10.0);
  //  hdum_vz->Fit("fit1","","",-0.15,0.15);
  
  TF1 *g1  = new TF1("g1","landau",-0.15,0.0);
  TF1 *g2  = new TF1("g2","landau", 0.0,0.15);
  TF1 *g3  = new TF1("g3","pol2",-0.15,0.15);
  TF1 *total = new TF1("total","landau(0)+landau(3)+pol2(6)",-0.15,0.15);
  total->SetLineColor(2);
  hdum_vz->Fit(g1,"R");
  hdum_vz->Fit(g2,"R+");
  hdum_vz->Fit(g3,"R+");

  g1->GetParameters(&par[0]);
  g2->GetParameters(&par[3]);
  g3->GetParameters(&par[6]);
  total->SetParameters(par);
  total->SetParNames("Peak_L","Mean_L","Sigma_L","Peak_R","Mean_R","Sigma_R","pol1","pol2","pol3");
  TCanvas *c2 = new TCanvas("c2","c2",1200,800); c2->cd();
  TH1F *hdum; hdum = hdum_vz;
  hdum->Fit(total,"R");
*/

//  cT_Dum->Delete();  T_Dum->Delete(); dummy->Close();
  
