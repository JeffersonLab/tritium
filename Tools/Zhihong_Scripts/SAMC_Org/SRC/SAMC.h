#ifndef SAMC_PHYSICAL_CONSTANT_H
#define SAMC_PHYSICAL_CONSTANT_H
/*Include{{{*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <errno.h>
#include <sstream>
//#include <iterator>
//from root
#include <TMath.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMatrix.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
/*}}}*/

#define CHAR_LEN 1024

#define NELEMENTS 4 //for RValue
#define LINE_SLOPE 0
#define LINE_INTERSECTION 1
#define LINE_SIGN 2


using namespace::std;
using namespace::TMath;

static const double ALPHA=1/137.036; //fine-structure constant
static const double PI=3.14159265358979323846;
static const double AP=ALPHA/PI;
static const double ELECTRON_MASS=0.510998918; //MeV
static const double M2=ELECTRON_MASS*ELECTRON_MASS; //MeV
static const double AMU=931.494043; //MeV
static const double PROTON_MASS=938.272029; //MeV
static const double DELTA_MASS=1232.; //MeV
static const double PI_MASS=139.57; //MeV
static const double MU_P=2.7928456; //MeV/T proton magnetic momentum
static const double MU_N=-1.91304184; //MeV/T neutron magnetic momentum
static const double HBARC=197.326968; //hbar*c MeV fm
static const double MEV2SR_TO_NBARNSR=3.89379*1e11; //1/(MeV^2 sr)--->nbarn/sr
static const double K_P=1.7927; //anomalous magnetic moment of Proton, used for Proton MeV fm
static const double K_N=-1.91348; //anomalous magnetic moment of Netutron, used for Neutron MeV fm
static const double K_HE=-4.2; //anomalous magnetic moment of He3, used for He3 and He4 MeV fm

//default unit
//MeV cm g rad
//cross section output unit = nbarn/sr
//
bool IsDebug=false; //output to debugfile
bool IsMyModel=false; //true: Use my own one-body model false: use old
bool IsMultiScat=true; //enable multi-scattering
bool IsEnergyLoss=false; //enable Energy Loss
int  Which_Beam_Profile=1; //0=Norm 1=Gaus 2=gaus+norm=triangle
int  Which_Kin=0; //0=phase space 1=elastic 2=quasi-elastic
double HRS_L=116; //Length Of HRS ( Left or Right,here _L=Length) (cm)
double Q1_Exit_Eff_L=40; //Q1 Effective Length  (cm)
double D_Entrance_Eff_L=400; //D Entrance Effective Length  (cm)
double D_Exit_Eff_L=700; //D Exit Effective Length  (cm)
double Q3_Entrance_Eff_L=775; //Q3 Entrance Effective Length  (cm)
double Q3_Exit_Eff_L=865; //Q3 Exit Effective Length  (cm)
double FP_Eff_L=1130; //FP Effective Length  (cm)
double Q1_Radius=14.92; //Q1 Radius  (cm)
double D_X_Radius=40; //Dipole X dir Radius  (cm)
double D_Y_L=12.5; //Dipole Y Length when x=0  (cm)
double Q3_Entrance_Radius=30; //Q3 Entrance Radius  (cm)
double Q3_Exit_Radius=30; //Q3 Exit Radius  (cm)
double VDC_Res_x=0.013;	//resolution of Wire chamer x  (cm)
double VDC_Res_y=0.013;	//resolution of Wire chamer y  (cm)
double VDC_Res_th=0.3;	//resolution of Wire chamer th (mr)
double VDC_Res_ph=0.3;	//resolution of Wire chamer ph (mr)
double delta_dp=0.14;	//d(dp) dp full width for generator
double delta_th=0.18;	//d(th) th full width for generator(tan(th))
double delta_ph=0.09;	//d(ph) th full width for generator(tan(ph))
double gaus_x_sigma=0;	//if Which_Beam_Profile==1/2 sigma of x beam for generator(cm)
double gaus_y_sigma=0;	//if Which_Beam_Profile==1/2 sigma of y beam for generator(cm)
double raster_x_size=0;	//raster x full size for generator(cm)
double raster_y_size=0;	//raster y full size for generator(cm)
double beam_x_center=0;	//beam x center for generator(cm)
double beam_y_center=0;	//beam y center for generator(cm)
double z0=0;	//target center for generator(cm)
double T_L=0;	//target length for generator(cm)
double T_H=0;	//target height for generator(cm)
double Fermi_Moment=0;	//Fermi_Moment,fermi momentum of target (MeV)
double NIE=0;	//NIE,nucleon interaction energy (MeV)
double DEL_SEP=0;	//DEL_SEP,delta separation energy (MeV)

bool IsRadCor=true; //enable Radiative correction
bool IsQFS_Q2dep=false; //enable Q2 dependence when calculating Q.E. Peak
/*Only for Tail{{{*/
int Form_Factor_Id=17; //which form factor used
bool  IsInternal_Bre=true; //enable/disable internal bremsstrahlung
bool  IsExternal_Bre=true; //enable/disable external bremsstrahlung
bool  IsMulti_Photon_Cor=true; //enable/disable multiple photon correction
bool  IsAppro=false; //use sigma_p instead of sigma_ex (approximate exact internal bremsstrahlung)
double DeltaE; //default=5 MeV only for radiative correction
/*}}}*/
double D_x=0; //X offset(in TCS) between TCS and HCS
double D_y=0; //Y offset(in TCS) between TCS and HCS
double E0;//=incident beam energy for generator (MeV)
double P0;//=HRS Setting Momentum for generator(MeV)

string File_Name; //input file name

/*for rvalue{{{*/
string RfunDB_FileName; //Rfun DB input file name
int    fNCuts=-1;
int**  fXY=NULL;
double**  fLineProperty=NULL;
/*}}}*/

#endif
