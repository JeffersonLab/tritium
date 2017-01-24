//////////////////////////////////////
//     F(y) Parameters Fitting      //
//     TMinuit Minimization         //
//    ---  Zhihong Ye 12/24/2012    //
//////////////////////////////////////
#include "fynk_minuit.hh"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"



Int_t main()
{ 
  gStyle->SetOptFit(1);  
  gStyle->SetOptStat(0);

  cerr<<"Target (H2,He3,He4,C12) = "; cin >> Target_Name;
  LoadData();
  DoMinuit();
  return 0;
}

/*DoMunuit{{{*/
void DoMinuit()
{
  cout << "      *************************************************" << endl; 
  cout << "      *          Minimization Fitting for             *" << endl;
  cout << "      *              F(y) Function                    *" << endl;
  cout << "      *             Z. Ye 09/20/2010                  *" << endl;
  cout << "      *************************************************" << endl; 
  cout << endl;
  gSystem->Load("libMinuit");
  static const Int_t iParaNum=5;

  TMinuit *pMinuit = new TMinuit(iParaNum);  //initialize TMinuit with a maximum of iParaNum params
  pMinuit->SetFCN(myFcn);

  Double_t arglist[10];
  Int_t ierflg = 0;

  /// set print level
  /*  SET PRIntout  <level>
      Sets the print level, determining how much output will be
      produced. Allowed values and their meanings are displayed
      after a SHOw PRInt command, and are currently <level>=:
      [-1]  no output except from SHOW commands
      [0]  minimum output
      [1]  default value, normal output
      [2]  additional output giving intermediate results.
      [3]  maximum output, showing progress of minimizations. */
  arglist[0] = 1;
  pMinuit->mnexcm("SET PRIntout",arglist,1,ierflg);

  /*
    SET NOWarnings
    Supresses Minuit warning messages.
    SET WARnings
    Instructs Minuit to output warning messages when suspicious
    conditions arise which may indicate unreliable results.
    This is the default.
  */
  arglist[0] = 1;
  pMinuit->mnexcm("SET NOWarnings",arglist,0,ierflg);

  //Set Name of Parameters
  Char_t *namePar[iParaNum]={"f0","bigB","alpha","a","b"};
 
  Double_t inPar[iParaNum]={Target->f0,   
			    Target->B,   
			    Target->alpha,
			    Target->a,
			    Target->b};	
 
  //Set Stepsize of Parameters
  Double_t step[iParaNum]={ 0.0000001,0.0000001,0.0001,0.0000001,0.0000001};
 
  //Set Min of Parameters, value==0 for No-Bound
  // Double_t minVal[iParaNum]={ 0.0001,0.0001,100.0,0.0001,0.0001};
  Double_t minVal[iParaNum]={ 0.0,0.0,0.0,0.0,0.0};

  //Set Max of Parameters, value==0 for No-Bound
  // Double_t maxVal[iParaNum]={ 0.01,0.01,300.0,0.01,0.01};
  Double_t maxVal[iParaNum]={ 0.0,0.0,0.0,0.0,0.0};

  //Initializing Parameters
  for(int ii=0;ii<iParaNum;ii++)
    {
      pMinuit->DefineParameter(ii,namePar[ii], 
			       inPar[ii], step[ii], minVal[ii], maxVal[ii]);
    }
  //Fix parameters
    arglist[0] = 3;   pMinuit->mnexcm("FIX ", arglist ,1,ierflg);

   arglist[0] = 1;
   pMinuit->mnexcm("SET ERR", arglist,1,ierflg);

  /*
  //SET LIMits  0->[parno]  1->[lolim]  2->[uplim]
  //	  arglist[0] = 0;	   //this means set all parameters with the same limit
  //	  arglist[0] = 1;	   //this means set 1st parameters with the specified limit
  //	  arglist[0] = 2;	   //this means set 2nd parameters with the sepecified limit
  arglist[0] = 0;
  arglist[1] = 0.;
  arglist[2] = 0.;
  pMinuit->mnexcm("SET LIMits",arglist,3,ierflg);
  */

  //Set Iteration Times and Use Migrad to start minimization
  arglist[0] = 500;
  arglist[1] = 0.01;
  pMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  //Get Results

  for(int ii = 0; ii < iParaNum; ii++)
    {
      pMinuit->GetParameter(ii,outPar[ii],err[ii]);    
    }
  //Put the results into a file
  cerr<<"Before:"<<endl;
  DoPlot(inPar,"Before");
  output(inPar);

  cerr<<"After:"<<endl;
  DoPlot(outPar,"After");
  output(outPar);

  //	LoadData();
}
/*}}}*/

// Calculate F(y), then get the Chi-Square Difference between Expertment Data,
// TMinuite will find a set of parameters to make this Chi-Square minimized
// Chi_Square = (F(y)- F_exp)**2
/*myFcn{{{*/
/*
  Instead of calculating Fy, calculate XS directly and direclty compare with expeirment XS
*/
void myFcn(Int_t& npar, Double_t* deriv, Double_t& f, Double_t *par, Int_t flag)
{
  
  Double_t aXS_EX = 0.0, aXS_EX_Err = 0.0;
  Double_t aXS_Fit = 0.0, aEp =0.0,aTheta = 0.0;
  Double_t aDelta = 0.0, aSum=0.0;
  //Set Parameters
  Target->f0 = par[0]; 
  Target->B = par[1]; 
  Target->alpha = par[2]; 
  Target->a = par[3]; 
  Target->b = par[4]; 
  
  for ( unsigned int ii=0; ii< vXbj.size(); ii++ )
    {
      aXS_EX = vXS[ii];
      aXS_EX_Err = vXSErr[ii];
      aEp = vEp[ii];
      aTheta = vTheta[ii];
      //F(y)
      aXS_Fit = FyToSig(Target, aEp, aTheta);

      aDelta = (aXS_EX - aXS_Fit) / aXS_EX_Err;
      aSum+=aDelta*aDelta;

    }
  
  f =  sqrt(aSum)/vXbj.size();
  //f = aSum;
  
  cerr << " >>> Chisq = " << f <<endl;
}
/*}}}*/

//Scaling Function F(y), from N. Fomin's Thesis (143P)

/*FyToSig{{{*/
Double_t FyToSig(XEM_TGT* kTarget, double kEp, double kTheta){

 //Get DeltaE value 
  Double_t DeltaE = 0.0;
  bool kCCor = kTRUE;
  if(kCCor){
    if(kTarget->A==1) //  !Hydrogen
      DeltaE = 0.0;
    else if(kTarget->A==2) // !Deuterium
      DeltaE = 0.0;
    else if(kTarget->A==3) // !Helium-3
      DeltaE = 0.00085;
    else if(kTarget->A==4) // !Helium-4
      DeltaE = 0.0010;
    else if(kTarget->A==9) // !Beryllium
      DeltaE = 0.001875;
    else if(kTarget->A==12) // !Carbon
      DeltaE = 0.00292;
    else if(kTarget->A==27) // !Aluminum
      DeltaE = 0.0061;
    else if(kTarget->A==64) // !Copper
      DeltaE = 0.0102;
    else if(kTarget->A==197) // !Gold
      DeltaE = 0.0199;
    else{
      cerr<<"**** Target not in the list! Coulumb Correction give up~ "<<endl;
      DeltaE = 0.0;
    }     
  }
  
  Double_t FF1 = (E0 - DeltaE)/E0;
  Double_t kE0_cc = E0 + DeltaE;
  Double_t kEp_cc = kEp + DeltaE;
 
  Double_t kXS = gCal_Fy2Sig(kE0_cc, kEp_cc, kTheta*pi/180.,kTarget) * FF1*FF1;
  //cerr<<"XS = "<<kXS<<endl;
  return kXS;

}
/*}}}*/


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
  for(int i=0;i<nTotal; i++){
    T->GetEntry(i);
    if(i==0){
      cerr<<Form("Initializing parameters for target %s, A=%d, Z=%d",Target_Name.Data(),A,Z)<<endl;
      Target->GetValueAZ(A,Z);
      Target->Print();
    }
    
    if(abs(Angle-Theta)<1e-5 && Ep>3.2 && xbj<2.0){

    vEp.push_back(Ep);
    vTheta.push_back(Theta);
    vXbj.push_back(xbj);
    vZero.push_back(0.0);
    nXS_QE = (nXS_Born - xs_dis);
    vXS.push_back(nXS_QE);
    vXSErr.push_back(nXS_QE*nFrac_Err);
    
    if(xMin >= Ep)
      xMin = Ep;
    if(xMax <=Ep)
      xMax = Ep;
    if(yMin >= nXS_QE)
      yMin = nXS_QE;
    if(yMax <= nXS_QE)
      yMax = nXS_QE;
    }
  }
  T->Delete();
}

/*DoFit{{{*/
void DoPlot(Double_t *par, const TString& kCom)
{
  int n=vXbj.size();

  TCanvas *c1=new TCanvas("c1","",800,600);
  c1->cd(1); gPad->SetLogy(1);
  TH2F *h1 = new TH2F("h1",Form("XS Fitting (%s)",kCom.Data()),1000, xMin-0.1,xMax+0.1,1000,yMin-0.01,yMax+0.01);
  h1->SetXTitle("Ep");
  h1->SetYTitle("#sigma");
  h1->Draw();
  TGraphErrors *gr=new TGraphErrors(n,&vEp[0],&vXS[0],&vZero[0],&vXSErr[0]);
  gr->SetMarkerStyle(21);
  gr->SetMarkerColor(4);
  gr->Draw("P");
  
  TF1 *thisF=new TF1("thisF",func,xMin-0.001,xMax+0.001,5);
  thisF->SetParameters(par[0],par[1],par[2],par[3],par[4]);
  // thisF->SetParameters(Target->f0,Target->B,Target->alpha,Target->a,Target->b);
  thisF->SetParNames("f0","B","alpha","a","b");
  thisF->SetLineColor(2);
  thisF->Draw("same");
     
  c1->Print(Form("fit_%s.png",kCom.Data()));

  delete gr;
  delete c1;
}
/*}}}*/

/*func{{{*/
Double_t func(Double_t *x, Double_t *par)
{
  //Set Parameters
  XEM_TGT* kTarget = new XEM_TGT();
  kTarget->GetValueAZ(Target->A,Target->Z);
  Double_t kEp = x[0];
  kTarget->f0 = par[0]; 
  kTarget->B = par[1]; 
  kTarget->alpha = par[2]; 
  kTarget->a = par[3]; 
  kTarget->b = par[4]; 
  
  Double_t kXS = FyToSig(kTarget, kEp,Angle);
  delete kTarget;
  return kXS;
 }
/*}}}*/

/*output{{{*/
//Output the results into a file
void output(Double_t *outPar)
{ 
  TString name =  Form("%s_%2.0f_Par.dat",Target_Name.Data(),Angle);

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
