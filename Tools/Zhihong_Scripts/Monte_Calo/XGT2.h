/*
Variables used by SAMC, simply extracted from SAMCEvent.h
Zhihong Ye, 2011/08/08
*/


#ifndef XGT2_h
#define XGT2_h

//C++ headers
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <cmath>
#include "unistd.h"
#include "math.h"
#include "stdio.h"
using namespace std;

//ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TMath.h"
#include "TObject.h"
#include "TStorage.h"
#include "TChain.h"


int Id; //Event Id
int IsPassed; //if pass through the magnet, 1=true, 0=false(not pass)
int IsQualified; //if rec var is in the range of gen, 1=true, 0=false(not

/*only for RadCor{{{*/
double mp_cor; //Multiple Photon Correction
double tail; //radiative tail
double tail_err; //error for radiative tail
double R; //Phys.Rev.D 12 1884 A83
 
double cs_qfs;//quasi-free scattering cross section
double cs_del;//Delta
double cs_x;  //DIS
double cs_r1500; //Resonace 1500 MeV
double cs_r1700; //Resonace 1700 MeV
double cs_2N; //Dip region
double cs_q;//sum of all the cross section above
 
double cs_p; //sigma_p (peak approximation for internal radation)
double cs_b_bt; //the part before target in sigma_b
double cs_b_at; //the part after target in sigma_b
double cs_b; //sigma_b (external bremsstrahlung ionization)
double cs_ex; //sigma_ex (extract internal radation)
double cs_qt;
/*}}}*/

/*Member Data from random{{{*/
double beam_x; //cm
double beam_y; //cm
double reactz_gen; // interaction z position cm
double E_s;//=incident beam energy (MeV)
double E_p;//=incident beam energy (MeV)
//gen=generator
double th_tg_gen; //theta target
double ph_tg_gen; //phi target
double dp_gen; //dp at target
/*}}}*/
  
/*Member Data derived from variables above{{{*/
double x_tg_gen; //
double y_tg_gen; //
/*}}}*/
  
/*Member Data for refinement{{{*/
double dp_ref; //dp for John.LeRose matrix
double x_tg_ref; //x tg for John.LeRose matrix
double y_tg_ref; //y tg for John.LeRose matrix
double th_tg_ref; //th tg for John.LeRose matrix
double ph_tg_ref; //ph tg for John.LeRose matrix
/*}}}*/
  
/*Member Data for focal plane {{{*/
double x_fp; //
double y_fp; //
double th_fp; //
double ph_fp; //
double th_fp_no_ort; //th_fp without TXFIT orthogonalization
double q1ex_x,q1ex_y;//0:x 1:y
double dent_x,dent_y;
double dext_x,dext_y;
double q3en_x,q3en_y;
double q3ex_x,q3ex_y;
  
int IsPassedQ1Ex;
int IsPassedDipoleEn;
int IsPassedDipoleEx;
int IsPassedQ3En;
int IsPassedQ3Ex;

/* bool IsPassedQ1Ex; */
/* bool IsPassedDipoleEn; */
/* bool IsPassedDipoleEx; */
/* bool IsPassedQ3En; */
/* bool IsPassedQ3Ex; */
/*}}}*/
  
/*Member Data for reconstruct target varibles {{{*/
double x_tg_rec; //cm
double y_tg_rec; //cm
double th_tg_rec; //
double ph_tg_rec; //
double dp_rec;
double Angle_rec;
double reactz_rec;
double rvalue;

double cs_M;//mott cross section
double cs_Final;//final cross section
/*}}}*/

/***********************************************************/
/**               Set BranchAddress                       **/
/***********************************************************/

void SetTree(TTree *tree){

  tree->SetBranchAddress("Id",&Id); //Event Id
  tree->SetBranchAddress("IsPassed",&IsPassed); //if pass through the magnet, 1=true, 0=false(not pass)
  tree->SetBranchAddress("IsQualified",&IsQualified); //if rec var is in the range of gen, 1=true, 0=false(not


  /*only for RadCor{{{*/
  tree->SetBranchAddress("mp_cor",&mp_cor   ); //Multiple Photon Correction
  tree->SetBranchAddress("tail",&tail   ); //radiative tail
  tree->SetBranchAddress("R",&R   ); //Phys.Rev.D 12 1884 A83
 
  tree->SetBranchAddress("cs_qfs",&cs_qfs   );//quasi-free scattering cross section
  tree->SetBranchAddress("cs_del",&cs_del   );//Delta
  tree->SetBranchAddress("cs_x",&cs_x   );  //DIS
  tree->SetBranchAddress("cs_r1500",&cs_r1500   ); //Resonace 1500 MeV
  tree->SetBranchAddress("cs_r1700",&cs_r1700   ); //Resonace 1700 MeV
  tree->SetBranchAddress("cs_2N",&cs_2N   ); //Dip region
  tree->SetBranchAddress("cs_q",&cs_q   );//sum of all the cross section above
 
  tree->SetBranchAddress("cs_p",&cs_p   ); //sigma_p (peak approximation for internal radation)
  tree->SetBranchAddress("cs_b_bt",&cs_b_bt   ); //the part before target in sigma_b
  tree->SetBranchAddress("cs_b_at",&cs_b_at   ); //the part after target in sigma_b
  tree->SetBranchAddress("cs_b",&cs_b   ); //sigma_b (external bremsstrahlung ionization)
  tree->SetBranchAddress("cs_ex",&cs_ex   ); //sigma_ex (extract internal radation)
  tree->SetBranchAddress("cs_qt",&cs_qt   );
  /*}}}*/
  

  /*Member Data from random{{{*/
  tree->SetBranchAddress("beam_x",&beam_x    ); //cm
  tree->SetBranchAddress("beam_y",&beam_y    ); //cm
  tree->SetBranchAddress("reactz_gen",&reactz_gen    ); // interaction z position cm
  tree->SetBranchAddress("E_s",&E_s    );//=incident beam energy (MeV)
  tree->SetBranchAddress("E_p",&E_p    );//=incident beam energy (MeV)
  //gen=generator
  tree->SetBranchAddress("th_tg_gen",&th_tg_gen    ); //theta target
  tree->SetBranchAddress("ph_tg_gen",&ph_tg_gen    ); //phi target
  tree->SetBranchAddress("dp_gen",&dp_gen    ); //dp at target
  /*}}}*/
  
  /*Member Data derived from variables above{{{*/
  tree->SetBranchAddress("x_tg_gen",&x_tg_gen    ); //
  tree->SetBranchAddress("y_tg_gen",&y_tg_gen    ); //
  /*}}}*/
  
  /*Member Data for refinement{{{*/
  tree->SetBranchAddress("dp_ref",&dp_ref    ); //dp for John.LeRose matrix
  tree->SetBranchAddress("x_tg_ref",&x_tg_ref    ); //x tg for John.LeRose matrix
  tree->SetBranchAddress("y_tg_ref",&y_tg_ref    ); //y tg for John.LeRose matrix
  tree->SetBranchAddress("th_tg_ref",&th_tg_ref    ); //th tg for John.LeRose matrix
  tree->SetBranchAddress("ph_tg_ref",&ph_tg_ref    ); //ph tg for John.LeRose matrix
  /*}}}*/
  
  /*Member Data for focal plane {{{*/
  tree->SetBranchAddress("x_fp",&x_fp    ); //
  tree->SetBranchAddress("y_fp",&y_fp    ); //
  tree->SetBranchAddress("th_fp",&th_fp    ); //
  tree->SetBranchAddress("ph_fp",&ph_fp    ); //
  tree->SetBranchAddress("th_fp_no_ort",&th_fp_no_ort    ); //th_fp without TXFIT orthogonalization
 
  tree->SetBranchAddress("q1ex_x",&q1ex_x    );//0:x 1:y
  tree->SetBranchAddress("dent_x",&dent_x    );
  tree->SetBranchAddress("dext_x",&dext_x    );
  tree->SetBranchAddress("q3en_x",&q3en_x    );
  tree->SetBranchAddress("q3ex_x",&q3ex_x    );
 
  tree->SetBranchAddress("q1ex_y",&q1ex_y    );//0:x 1:y
  tree->SetBranchAddress("dent_y",&dent_y    );
  tree->SetBranchAddress("dext_y",&dext_y    );
  tree->SetBranchAddress("q3en_y",&q3en_y    );
  tree->SetBranchAddress("q3ex_y",&q3ex_y    );
  
  tree->SetBranchAddress("IsPassedQ1Ex",&IsPassedQ1Ex    );
  tree->SetBranchAddress("IsPassedDipoleEn",&IsPassedDipoleEn    );
  tree->SetBranchAddress("IsPassedDipoleEx",&IsPassedDipoleEx    );
  tree->SetBranchAddress("IsPassedQ3En",&IsPassedQ3En    );
  tree->SetBranchAddress("IsPassedQ3Ex",&IsPassedQ3Ex    );

  /*}}}*/
  
  /*Member Data for reconstruct target varibles {{{*/
  tree->SetBranchAddress("x_tg_rec",&x_tg_rec    ); //cm
  tree->SetBranchAddress("y_tg_rec",&y_tg_rec    ); //cm
  tree->SetBranchAddress("th_tg_rec",&th_tg_rec    ); //
  tree->SetBranchAddress("ph_tg_rec",&ph_tg_rec    ); //
  tree->SetBranchAddress("dp_rec",&dp_rec    );
  tree->SetBranchAddress("reactz_rec",&reactz_rec    );
  tree->SetBranchAddress("Angle_rec",&Angle_rec    );
  tree->SetBranchAddress("rvalue",&rvalue    );

  tree->SetBranchAddress("cs_M",&cs_M    );//mott cross section
  tree->SetBranchAddress("cs_Final",&cs_Final    );//final cross section
/*}}}*/
}

#endif
