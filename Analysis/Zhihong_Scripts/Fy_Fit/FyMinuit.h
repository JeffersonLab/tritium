/*Headers{{{*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <string>
#include <cstdio>

#include "math.h"
#include "stdio.h"
#include <stdlib.h>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TObject.h"
#include "TSystem.h"
#include <TVirtualFitter.h>
#include <TMinuit.h>
//#include <TFitterMinuit.h>
#include <TApplication.h>
#include <vector>
using namespace std;

/*}}}*/
#include "XEMC_Born.h"
//#include "XEM_Main.h"
#include "XEM_Sig2Fy.h"

Double_t FyToSig(XEM_TGT* kTarget, double kEp, double kTheta);
void myFcn(Int_t& npar, Double_t* deriv, Double_t& f, Double_t *par, Int_t flag);
Double_t func(Double_t *x, Double_t *par);
void DoPlot(Double_t* par,const TString& kCom);
void output(Double_t *outPar, const int aChoice);
void LoadData();
void LoadDataXGT2();
void LoadDataAscii();
void DoMinuit(double *par);
void DoGraphFit(double *par);
void DoNkFit(double* par);
void DoRatio(double* par_old, double* par_new, const int aChoice);
void DoRatio_Angle(double* par_old, double* par_new, const int aChoice);

static const Double_t pi = 3.1415926;
//Input Data:
const double E0 = 5.766;
TString Target_Name = "C12";
Double_t outPar[5], err[5];

vector <double> vEp;
vector <double> vE0;
vector <double> vZero;
vector <double> vTheta; 
vector <double> vXS; 
vector <double> vXbj; 
vector <double> vXSErr; 
vector <double> vFy,vFy_XGT2,vFy_NF;
vector <double> vFyErr,vFyErr_XGT2,vFyErr_NF;
vector <double> vY,vY_XGT2,vY_NF, vTheta_XGT2, vTheta_NF;

XEM_TGT* Target = new XEM_TGT();  
XEM_SIG* Fy2Sig = new XEM_SIG();
Double_t xMin = 1000., xMax = -1000.;
Double_t yMin = 1000., yMax = -1000.;
Double_t xbjMin = 1000., xbjMax = -1000.;
