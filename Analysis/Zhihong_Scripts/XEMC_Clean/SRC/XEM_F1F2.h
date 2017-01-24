#include "XEM_Constant.h"
extern "C"
{
 void f1f2in09_(const double*,const double*, const double*, const double*, double*, double*, double*);

}

inline void gCal_F1F2In06(int kA, int kZ, double kQsq, double kWsq,XEM_VAR3* f1f2);

inline void gCal_F1F2(const int kA, const int kZ, const double kQsq, const double kWsq,XEM_VAR3* f1f2)
{
/*
//--------------------------------------------------------------------
// This is the updated version of F1F2 subroutine. For more info:
//    http://arxiv.org/abs/1203.2262
// Original Code locates at:
//    https:/userweb.jlab.org/~bosted/fits.html
//
// Directly call the FORTRAN subroutine instead of converting into C++
//
//    --- Zhihong Ye, 12/03/2012
//--------------------------------------------------------------------
*/

  f1f2->First=0.0;
  f1f2->Second=0.0;
  f1f2->Third=0.0;

// gCal_F1F2In06( kA, kZ, kQsq, kWsq, f1f2);
 
  //=======================================
  // Updated F1F2 in Fortran Subroutines
  //=======================================
  double kF1 = -1.1, kF2 = -1.2, kRC = -1.3;
  double dZ = (double) (kZ); double dA = (double) (kA);
  f1f2in09_(&dZ, &dA, &kQsq, &kWsq, &kF1, &kF2, &kRC);
// DEBUG -- For D2, Can not get F1 and F2 while Wsq < 0, Z. Ye, 12/05/2012
	  if(kA==2 && kZ==1 && kWsq<-1e-33){
		 kF1 = 0.0; kF2 = 0.0; kRC = 0.0;
//	 cerr<<Form("++++ Z=%d, A=%d, Q2=%f, W2=%f, F1=%f, F2=%f, RC=%f", kZ, kA, kQsq, kWsq, kF1, kF2, kRC)<<endl;
	  }

  f1f2->First = kF1;
  f1f2->Second = kF2;
  f1f2->Third=kRC;
  return;
}

inline void gCal_F1F2In06(int kA, int kZ, double kQsq, double kWsq,XEM_VAR3* f1f2)
{
/*
//--------------------------------------------------------------------
// Fit to inelastic cross sections for A(e,e')X
// valid for all W<3 GeV and all Q2<10 GeV2
// 
// Inputs: Z, A (real*8) are Z and A of nucleus 
//         (use Z=0., A=1. to get free neutron)
//         Qsq (real*8) is 4-vector momentum transfer squared (positive in
//                     chosen metric)
//         Wsq (real*8) is invarinat mass squared of final state calculated
//                     assuming electron scattered from a free proton
//                 
// outputs: F1, F2 (real*8) are structure functions per nucleus
// Version of 10/20/2006 P. Bosted
//--------------------------------------------------------------------
*/

  f1f2->First=0.0;
  f1f2->Second=0.0;
  f1f2->Third=0.0;
  XEM_VAR4* temp = new XEM_VAR4();

  Double_t  nu = (kWsq - PM_SQ + kQsq)/2.0/P_MASS;
  Double_t Qv = sqrt(nu*nu + kQsq);
  
  Double_t Es,kf,pf,dW2DPF,dW2DES,Wsqp;
  Double_t Rc=0.,sigt=0.,sigl=0.,F1p=0.,F1d=0.,W1=0.0,W2=0.0;  
  Double_t Rcp=0.,siglp=0.,sigtp=0.,F1pp=0.,F1dp=0.;
  Double_t kXval[50];
  for(int i=0;i<50;i++)
    kXval[i] = cXvald0[i];

  // Cross section for proton or neutron
  if(kA < 2){ 
    gGet_Christy806(kQsq,kWsq,temp);
    F1p = temp->First;    Rc = temp->Second;
    sigl = temp->Third;   sigt = temp->Forth;
  
    // If neutron, subtract proton from deuteron. 
    // Factor of two to convert from per nucleon to per deuteron
    if(kZ < 0.5){ 
      F1d = gGet_Resmod_Hack(kQsq,kWsq,kXval);
      F1p = F1d * 2.0 - F1p;
    }
      
    W1 = F1p/P_MASS ;
    W2 = W1 * (1.0+ Rc)/(1.0 + nu*nu/kQsq);
  }// if(kA < 2){ 
  
  //For deuteron
  else if(kA == 2){
    // get Fermi-smeared R from Erics proton fit
    gGet_Pind(kQsq, kWsq,temp);
    F1p = temp->First;    Rc = temp->Second;
    sigl = temp->Third;   sigt = temp->Forth;
  
    // get fit to F1 in deuteron, per nucleon
    F1d=gGet_Resd(kQsq, kWsq,kXval);
      // convert to W1 per deuteron
    W1 = F1d/P_MASS*2.0;
    W2 = W1 * (1. + Rc) / (1.0 + nu*nu /kQsq);
  }//else if(kA == 2){
  
  //For nuclei
  else if(kA>2){
    sigt = 0.;    sigl = 0.;
    F1d = 0.;     F1p = 0.;
    // Modifed to use Superscaling from Sick, Donnelly, Maieron,
    // nucl-th/0109032
    if(kA==2){       kf=0.085;      Es=0.0022;    }
    if(kA==3){       kf=0.180;      Es=0.010;     }
    if(kA==4){       kf=0.200;      Es=0.015;     }
    if(kA>4 ){       kf=0.165;      Es=0.015;     }    
    if(kA>7 ){       kf=0.228;      Es=0.020;     } 
    if(kA>16){       kf=0.230;      Es=0.025;     }
    if(kA>25){       kf=0.236;      Es=0.018;     }
    if(kA>38){       kf=0.241;      Es=0.028;     }
    if(kA>55){       kf=0.241;      Es=0.023;     } 
    if(kA>60){       kf=0.245;      Es=0.028;     }
    
    //adjust pf to give right width based on kf
    pf = 0.5*kf; 
    // assume this is 2*pf*Qv
    dW2DPF = 2.0*Qv;
    dW2DES = 2.0*(nu + P_MASS); 
    for(int ism=0; ism<15;ism++){
      Wsqp = kWsq + cXX[ism]*pf*dW2DPF-Es*dW2DES;
      if(Wsqp > 1.159){
	gGet_Christy806(kQsq,Wsqp,temp);
	F1pp = temp->First; 
	Rcp = temp->Second;
	siglp = temp->Third; 
	sigtp = temp->Forth;

   	F1dp = gGet_Resmod_Hack(kQsq,Wsqp,kXval);
     	F1d = F1d + F1dp * cFy[ism];
	F1p = F1p + F1pp * cFy[ism];
	sigt = sigt + sigtp * cFy[ism];
	sigl = sigl + siglp * cFy[ism];
      }
    }
  

    Rc = 0.;
    if(sigt > 0.){
      Rc = sigl / sigt;
      W1 = (2.0*kZ*F1d+(kA-2.0*kZ)*(2.0*F1d-F1p))/P_MASS; 
      W2 = W1*(1.0+Rc)/(1.0+nu*nu /kQsq); 
    }
  }// else if(kA>2){
  else{
    cerr<<"****** ERROR, Wrong Target, Cannot calculate F1F2!"<<endl;
    return;
  }
  delete temp;
  f1f2->First = P_MASS*W1;
  f1f2->Second = nu*W2;
  f1f2->Third=Rc;
  return;
  }

//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////
// Christy fit to proton
//inline XEM_VAR4* gGet_Christy806(double aQsq, double aWsq)
inline void gGet_Christy806(double aQsq, double aWsq, XEM_VAR4 *sig)
{
  //  XEM_VAR4* sig = new XEM_VAR4();
  Double_t aR=0.0,aF1=0.0,aFL=0.0,aNu=0.0,aSigL=0.0,aSigT=0.0,aW1p=0.0,aW2p=0.0;
  if(aWsq > 1.155){
    Double_t xb = aQsq/(aWsq+aQsq-PM_SQ);
    if(xb>0.0){
      aSigT = gGet_Resmod316(aQsq,aWsq,1);
      aSigL = gGet_Resmod316(aQsq,aWsq,2);

      aF1 = aSigT*(aWsq-PM_SQ)/8.0/TMath::Pi()/TMath::Pi()/Alpha/(HC_SQ*GeVToMeV);
      aFL = aSigL*2.0*xb*(aWsq-PM_SQ)/8.0/TMath::Pi()/TMath::Pi()/Alpha/(HC_SQ*GeVToMeV);
      aR = aFL/(2.0*xb*aF1);
      //Why we need this?! FIX_HERE
      aNu = aQsq/2.0/P_MASS/xb;
      aW1p = aF1/P_MASS;
      aW2p = aW1p/(1.0+aNu*aNu/aQsq)*(1.0+aR);
    }
  }
  sig->First = aF1;
  sig->Second = aR;
  sig->Third = aSigL;
  sig->Forth = aSigT;

  return;
}//inline XEM_VAR4* gGet_Christy806(double aQsq, double aWsq)

//////////////////////////////////////////////////////////////////////
//  Version 031606  -  Author:  M.E. Christy                        //
//  Fit form is empirical.  Please do not try to interpret physics  //
//  from it.  This routine needs the parameter files f1parms.dat    //
//  and fLparms.dat.  Units are ub/Sr/Gev.                          //
//////////////////////////////////////////////////////////////////////
inline Double_t gGet_Resmod316(double aQsq, double aWsq, int aflag)
{
  //Double_t* cXval=new double[50];
  double* cXval=new double[50];
  
  for(int i=0;i<50;i++){
    if(aflag==1)
      cXval[i] = cXval1[i];
    else if(aflag==2)
      cXval[i] = cXvalL[i];
  }

  Double_t aW = sqrt(aWsq);
  Double_t Wdiff = aW - (P_MASS+PION_MASS);
 
  //Single Pion Branching Retias
  Double_t Pion_BR[2][7] = {{1.0,0.5,0.65,0.65,0.4,0.65,0.6},
			    {0.0,0.0,0.0,0.0,0.0,0.0,0.0}};
  if(aflag==2){
    Pion_BR[0][1]=cXval[48]; 
    Pion_BR[0][5]=cXval[47];
  }
  for(int k=0;k<7;k++)
    Pion_BR[1][k] = 1.0 - Pion_BR[0][k];  

  Double_t ang[7]={1.,0.,2.,3.,0.,1.,3.};
  Double_t aX0[7]={0.165,0.165,0.165,0.165,0.165,0.165,0.165};
  aX0[3]=0.6;

  double Q2_Low; int lmax=1;
  bool lowq2 = kFALSE;
  if(aflag==1)
    Q2_Low =0.0;
  else
    Q2_Low =0.1;
  if(aQsq<Q2_Low){
    lmax=2;
    lowq2 = kTRUE;
  }
  Double_t Q2=aQsq;// A new Q2 to replace Qsq if necessary
  Double_t dQ2 = 0.05;
  Double_t sig_res=0.0,sig_nr=0.0,sig_temp=0.0,sig=0.0,slope=0.0;
  Double_t dIP,dIP2,xb,xpr,k,kcm,epicm,ppicm,ppi2cm,epi2cm,eetacm,petacm;
  Double_t kr[7],kcmr[7],epicmr[7],ppicmr[7],ppi2cmr[7],epi2cmr[7],eetacmr[7],petacmr[7],pgam[7];
  Double_t sigr[7],aMass[7],aWidth[7],aIntWidth[7],height[7];
  Double_t pwid[2][7],rescoef[6][4],nr_coef[3][4];
  Int_t aNum=0;

  for(int iLoop=1;iLoop<=lmax;iLoop++){
    if(iLoop==1&&lowq2)
      Q2 = Q2_Low;
    else if(iLoop==2&&lowq2)
      Q2 = Q2_Low + dQ2;
 
    //Dipole Parameterization
    dIP = 1.0/pow((1.0+Q2/0.71),2);
    dIP2 = dIP*dIP;
    xb = Q2/(Q2+aWsq-PM_SQ);
    xpr = 1.0+(aWsq-(P_MASS+PION_MASS)*(P_MASS+PION_MASS))/(Q2+cXval[49]);
    xpr = 1.0/xpr;
    
    k = (aWsq -PM_SQ)/2.0/P_MASS;
    kcm = (aWsq-PM_SQ)/2.0/aW;
    epicm = (aWsq + PI_SQ -PM_SQ )/2.0/aW;
    ppicm = sqrt(gGet_Max(0.0,(epicm*epicm - PI_SQ)));
    epi2cm = (aWsq + 4.0*PI_SQ -PM_SQ )/2./aW;
    ppi2cm = sqrt(gGet_Max(0.0,(epi2cm*epi2cm-4.0*PI_SQ)));
    eetacm = (aWsq + META*META -PM_SQ)/2.0/aW;
    petacm = sqrt(gGet_Max(0.0,(eetacm*eetacm - META*META)));

    aNum = 0;
    for(int i=0;i<6;i++){
      aMass[i] = cXval[aNum];	 //  Read in resonance masses
      aNum++;
    }
    for(int i=0;i<6;i++){
      aIntWidth[i] = cXval[aNum];//Read in resonance widths    
      aWidth[i] = aIntWidth[i];
      aNum++;
    }
    
    if(aflag==2){//  Put in 4th resonance region  !!!
      aMass[6] = cXval[40];
      aIntWidth[6] = cXval[41];
      aWidth[6] = aIntWidth[6];
    }
    else{
      aMass[6] = cXval[46];
      aIntWidth[6] = cXval[47];
      aWidth[6] = aIntWidth[6] ;
    }

    for(int i=0;i<7;i++){
      kr[i] = (aMass[i]*aMass[i]-PM_SQ)/2.0/P_MASS;
      kcmr[i] = (aMass[i]*aMass[i]-PM_SQ)/2.0/aMass[i];
      epicmr[i] = (aMass[i]*aMass[i] + PI_SQ -PM_SQ)/2.0/aMass[i];
      ppicmr[i] = sqrt(gGet_Max(0.0,(epicmr[i]*epicmr[i]-PI_SQ)));
      epi2cmr[i] = (aMass[i]*aMass[i]+4.0*PI_SQ-PM_SQ)/2.0/aMass[i];
      ppi2cmr[i] = sqrt(gGet_Max(0.0,(epi2cmr[i]*epi2cmr[i] - pow((2.0*PION_MASS),2))));
      eetacmr[i] = (aMass[i]*aMass[i] + META*META -PM_SQ)/2.0/aMass[i];
      petacmr[i] =  sqrt(gGet_Max(0.0,(eetacmr[i]*eetacmr[i] - META*META)));

      //Calculate partial widths   
      pwid[0][i] = aIntWidth[i]
	         * pow((ppicm/ppicmr[i]),(2.*ang[i]+1.))
	         * pow(((ppicmr[i]*ppicmr[i]+aX0[i]*aX0[i])/(ppicm*ppicm+aX0[i]*aX0[i])),ang[i]);
      pwid[1][i] = aIntWidth[i]
	         * pow((ppi2cm/ppi2cmr[i]),(2.0*ang[i]+1.))
	         * pow(((ppi2cmr[i]*ppi2cmr[i]+aX0[i]*aX0[i])/(ppi2cm*ppi2cm+aX0[i]*aX0[i])),ang[i]);
   
      if(i==1)
	pwid[1][i] = aIntWidth[i]
	           * pow((petacm/petacmr[i]),(2.0*ang[i]+1.))
	           * pow(((petacmr[i]*petacmr[i]+aX0[i]*aX0[i])/(petacm*petacm+aX0[i]*aX0[i])),ang[i]);
      
      pgam[i] = pow((kcm/kcmr[i]),2)
	      * (kcmr[i]*kcmr[i]+aX0[i]*aX0[i])/(kcm*kcm+aX0[i]*aX0[i]);
      
      pgam[i] *= aIntWidth[i];
      aWidth[i] = Pion_BR[0][i]*pwid[0][i]+Pion_BR[1][i]*pwid[1][i];
    }// for(int i=0;i<7;i++), 
    // End resonance kinematics and Widths calculations   

    // Begin resonance Q^2 dependence calculations //
    for(int i=0;i<6;i++){ 
      for(int j=0;j<4;j++){
	rescoef[i][j]=cXval[aNum];
	aNum++;
      }
      
      if(aflag==1){
	if(i==5) //FIX_HERE, double check the formula
	  height[i] = rescoef[i][0]
	    /pow((1.0+ Q2/rescoef[i][1]),(rescoef[i][2]+rescoef[i][3]*Q2));
	
	height[i] = rescoef[i][0]
	  *(1.+rescoef[i][1]*Q2/(1.+rescoef[i][2]*Q2))
	  /pow((1.+Q2/0.71),rescoef[i][3]);
      }
      else{
	height[i] = rescoef[i][0]*Q2/(1.+rescoef[i][3]*Q2)
	  /pow((1.+Q2/rescoef[i][1]),rescoef[i][2]);// Galster Shape             
      }
    }// for(int i=0;i<6;i++), 
    // End resonance Q^2 dependence calculations   //

    for(int i=0;i<3;i++){//  Non-Res coefficients 
      for(int j=0;j<4;j++){
	nr_coef[i][j]=cXval[aNum];
	  aNum++;
      }
    }

    if(aflag==2)// !!!  4th resonance region  !!!
      height[6] = cXval[42]*Q2/(1.0+cXval[43]*Q2)
	/pow((1.+Q2/cXval[44]),cXval[45]);//     !!!  Galster Shape  !!!
    else
      height[6] = cXval[48]*dIP2; 
   
// Calculate Breit-Wigners for the 3 resonance regions   //
    sig_res = 0.0;
    for(int i=0;i<7;i++){
      sigr[i] = Pion_BR[0][i]*pwid[0][i]+Pion_BR[1][i]*pwid[1][i];
      sigr[i] *= pgam[i]
	/(pow((aWsq - aMass[i]*aMass[i]),2) + pow((aMass[i]*aWidth[i]),2));
      sigr[i] *= height[i]*kr[i]/k*kcmr[i]/kcm/aIntWidth[i];
      sig_res = sig_res + sigr[i];              
    }    
    // Finish resonances / start non-res background calculation   //
    sig_nr = 0.;
    if(aflag==1){
      for(int i=0;i<2;i++)
	sig_nr = sig_nr+nr_coef[i][0]*pow(Wdiff,((2.*(i+1)+1.)/2.))
	       / pow((Q2+nr_coef[i][1]),(nr_coef[i][2]+nr_coef[i][3]*Q2+cXval[44+i]*Q2*Q2));
      sig_nr = sig_nr*xpr;
    }
    else if(aflag==2){
      for(int i=0;i<1;i++){
	sig_nr = sig_nr +nr_coef[i][0]*pow(Wdiff,(i+1)/2.)
	       * pow((1.-xpr),2.)*Q2/(1.+nr_coef[i][2]*Q2)
	       / pow((1.+ Q2/nr_coef[i][1]),nr_coef[i][3])/aWsq;
      }
    }
    sig = sig_res + sig_nr;
    
    if(iLoop==1) 
      sig_temp = sig;
  } // for(int iLoop=1;iLoop<lmax;iLoop++)

  Q2 = aQsq;
  if(lowq2){
    if(aflag==1){
      slope = (sig_temp - sig)/dQ2;
      sig = sig_temp + slope*(Q2_Low-Q2);
    }
    else{
      slope = sig/Q2_Low;
      sig = sig - slope*(Q2_Low-Q2);     
    }
  }
  delete[]  cXval;
  return sig;
}//inline Double_t gGet_Resmod316(double aQsq, double aWsq, int aflag)

////////////////////////////////////////////////////////////////
// Calculate proton with Fermi smearing of a deuteron 
////////////////////////////////////////////////////////////////
inline void gGet_Pind(double aQsq, double aWsq, XEM_VAR4* sig)
{
  Double_t nu = (aWsq-PM_SQ+aQsq)/2.0/P_MASS;
  Double_t Qv = sqrt(nu*nu + aQsq);
  XEM_VAR4* temp = new XEM_VAR4();
    
  Double_t W2p=-1000.0,sigl=-1000.0,siglp=-1000.0,sigt=-1000.0,sigtp=-1000.0;
  Double_t R=-1000.0,Rp=-1000.0,F1=-1000.0,F1p=-1000.0;
  Double_t pz=-1000.0;
// Do fast 20 bins if above threshold
  if(aWsq>1.16){
    for(int ism=0;ism<20;ism++){  
      // try with energy term zero. Fix sign of qv * pz
      W2p = pow((D_MASS+nu-P_MASS),2)-
	Qv*Qv + 2.0*Qv*cAvpz[ism]-cAvp2[ism];
      if(W2p>1.155){ 
	gGet_Christy806(aQsq,W2p,temp);
	F1p = temp->First; Rp = temp->Second;
	siglp = temp->Third; sigtp = temp->Forth;
	
	sigt += sigtp*cFyd[ism]/10.;
	sigl += siglp*cFyd[ism]/10.;
	F1   += F1p*cFyd[ism]/10.;
      }
    }
  }
  else{
    for(int ism=0;ism<200;ism++){  
      pz = -1.0+0.01*(ism+0.5);
      // Need avp2f term to get right behavior x>1// 
      W2p = pow((D_MASS+nu-sqrt(PM_SQ + cAvp2f[ism])),2)
	- Qv*Qv+2.0*Qv*pz-cAvp2f[ism];
      if(W2p>1.155){
	gGet_Christy806(aQsq,W2p,temp);
	F1p = temp->First; Rp = temp->Second;
	siglp = temp->Third; sigtp = temp->Forth;
	sigt += sigtp*cFydf[ism]/100.;
	sigl += siglp*cFydf[ism]/100.;
	F1   += F1p*cFydf[ism]/100.;
      }
    }
  }
  
  if(!(sigt==0.)) 
    R = sigl / sigt;
  
  sig->First = F1;
  sig->Second = R;
  sig->Third = sigl;
  sig->Forth = sigt;
  
  delete temp;
  return;
}//inline XEM_VAR4* gGet_Pind(double aQsq, double aWsq)
      
////////////////////////////////////////////////////////////////
// Calculate dueteron F1 by Fermi smearing of proton plus neutron 
////////////////////////////////////////////////////////////////
inline Double_t gGet_Resd(double aQsq, double aWsq,double *aXval)
{
  double* cXval= new double[50];
  for(int i=0;i<50;i++)
    cXval[i] = aXval[i];
  
  Double_t nu = (aWsq - PM_SQ + aQsq)/2./P_MASS;
  Double_t Qv = sqrt(nu*nu + aQsq);
  Double_t F1 = 0.0;
  Double_t pz,W2p,xbj,sigp;

  // Do fast 20 bins if abvoe threshold
  for(int ism=0;ism<200;ism++){
    pz = -1.0+0.01*(ism+0.5);
    // Need avp2f term to get right behavior x>1//
    W2p = pow((D_MASS + nu - sqrt(PM_SQ + cAvp2f[ism])),2) 
      -Qv*Qv + 2.0*Qv*pz - cAvp2f[ism];
    xbj = aQsq/(W2p+aQsq-PM_SQ);
    
    if(W2p>1.155){
      sigp = gGet_Resmod_Hack(aQsq,W2p,cXval); //sigp is really F1deut/2.0
      F1 += sigp*cFydf[ism]/100.; 
    }
  }
  delete[] cXval;
  return F1;
}//inline Double_t gGet_Resd(double aQsq, double aWsq)


////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////
inline Double_t gGet_Resmod_Hack(double aQsq,double aWsq,double* aXval)
{

  double* cXval= new double[50];
  for(int i=0;i<50;i++)
    cXval[i] = aXval[i];
  
  Double_t sig = 0.;
  if(aWsq<1.07327*1.07327 ||  aWsq>25 || aQsq < 0.0 || aQsq > 10.01)
    return sig;
  
  XEM_VAR4* temp = new XEM_VAR4();

  //Single Pion Branching Retias
  Double_t Pion_BR[2][7] = {{1.0,0.5,0.65,0.65,0.4,0.65,0.6},
			    {0.0,0.5,0.35,0.35,0.6,0.35,0.4}};
  // angular momenta
  //      ang(1) = 1.       //////  P33(1232)
  //      ang(2) = 0.       //////  S11(1535)
  //      ang(3) = 2.       //////  D13(1520)
  //      ang(4) = 3.       //////  F15(1680)
  //      ang(5) = 0.       //////  S15(1650)
  //      ang(6) = 1.       //////  P11(1440) roper   
  //      ang(7) = 3.       //////  ? 4th resonance region
  const Double_t ang[7]={1.,0.,2.,3.,0.,1.,3.};
  const Double_t aX0[7]={0.165,0.165,0.165,0.6,0.165,0.165,0.165};
  
  // remember w2
  // Double_t W2_Save = aWsq;
  Double_t Q2 = aQsq;

  Double_t aMass[7],aIntWidth[7];
  // uses xvals of 1-12, 47, and 48, move masses, wdiths into local variables
  // pyb changed to be fixed
  Int_t aNum = 0;
  for(int i=0;i<6;i++){  
    aMass[i] = cXval0[aNum];
    aNum++;
  }
  for(int i=0;i<6;i++){  
    aIntWidth[i] = cXval0[aNum];
    aNum++;
  }
  aMass[6] = cXval[46];
  aIntWidth[6] = cXval[47];
  
  // precalculate w-dependent quantites in 0.1 MeV bins
  Double_t w,w2,wdiff,wr,k,kcm,epicm,ppicm,epi2cm,ppi2cm,eetacm,petacm,dIP,dIP2;
  Double_t sig_res,sig_nr,xpr,F1,xb,nu;
  Double_t kr[7],kcmr[7],epicmr[7],ppicmr[7],epi2cmr[7],ppi2cmr[7],eetacmr[7],petacmr[7],pgam[7],width[7];
  Double_t height[7],sig_rsv[7];
  Double_t pwid[2][7],rescoef[6][4],nr_coeff[6][4];
  Double_t sigr[7][5001];

  for(int i=0;i<7;i++){

    for(int j=0;j<5001;j++){
      sigr[i][j]=0.0;
    }
  }

  for(int iw=1073;iw<=5000;iw++){
    w = (iw+0.5)*MeVToGeV;
    w2 = w*w;
    wdiff = w-(P_MASS+PION_MASS);
    wr = wdiff/w;

    // Calculate kinematics needed for threshold Relativistic B-W 
    k = (w2-PM_SQ)/2.0/P_MASS;
    kcm = (w2-PM_SQ)/ 2./ w;
    epicm = (w2+PI_SQ -PM_SQ ) / 2. / w;
    ppicm = sqrt(gGet_Max(0.0,(epicm*epicm - PI_SQ)));
    epi2cm = (w2+4.*PI_SQ -PM_SQ )/2.0/w;
    ppi2cm = sqrt(gGet_Max(0.0,(epi2cm*epi2cm-4.0*PI_SQ)));
    eetacm = (w2+META*META-PM_SQ)/2.0/w;
    petacm =  sqrt(gGet_Max(0.0,(eetacm*eetacm-META*META)));
 
    for(int i=0;i<7;i++){
      kr[i] = (aMass[i]*aMass[i]-PM_SQ)/2./P_MASS;
      kcmr[i] = (aMass[i]*aMass[i]-PM_SQ)/2./aMass[i];
      epicmr[i] = (aMass[i]*aMass[i] + PI_SQ -PM_SQ )/2./aMass[i];
      ppicmr[i] = sqrt(gGet_Max(0.0,(epicmr[i]*epicmr[i]-PI_SQ)));
      epi2cmr[i] = (aMass[i]*aMass[i]+4.0*PI_SQ -PM_SQ )/2./aMass[i];
      ppi2cmr[i] = sqrt(gGet_Max(0.0,(epi2cmr[i]*epi2cmr[i]-4.*PI_SQ)));
      eetacmr[i] = (aMass[i]*aMass[i] + META*META -PM_SQ )/2./aMass[i];
      petacmr[i] =  sqrt(gGet_Max(0.0,(eetacmr[i]*eetacmr[i]-META*META)));
      
      // Calculate partial widths
      pwid[0][i] = aIntWidth[i]
	*pow((ppicm/ppicmr[i]),(2.*ang[i]+1.))
	*pow((ppicmr[i]*ppicmr[i]+aX0[i]*aX0[i])/(ppicm*ppicm+aX0[i]*aX0[i]),ang[i]);

      if(!(i==1)) 
	pwid[1][i] = aIntWidth[i]
	  *pow((ppi2cm/ppi2cmr[i]),(2.*ang[i]+1.))
	  *pow((ppi2cmr[i]*ppi2cmr[i]+aX0[i]*aX0[i])/(ppi2cm*ppi2cm+aX0[i]*aX0[i]),ang[i]);
      else
	pwid[1][i] = aIntWidth[1]
	  *pow((petacm/petacmr[i]),(2.*ang[i]+1.))
	  *pow(((petacmr[i]*petacmr[i]+aX0[i]*aX0[i])/(petacm*petacm+aX0[i]*aX0[i])),ang[i]);
      
      pgam[i] = pow((kcm/kcmr[i]),2)
	*(kcmr[i]*kcmr[i]+aX0[i]*aX0[i])/(kcm*kcm+aX0[i]*aX0[i]);
      pgam[i] = aIntWidth[i]*pgam[i];
      width[i] = Pion_BR[0][i]*pwid[0][i]+Pion_BR[1][i]*pwid[1][i];
      sigr[i][iw] = width[i]*pgam[i]
	/(pow((w2-aMass[i]*aMass[i]),2)+ pow((aMass[i]*width[i]),2)) 
	*kr[i]/k*kcmr[i]/kcm/aIntWidth[i];
    } // loop on i
  }	      // loop on iw

  // resonance height coefficients. xvals of 13-36
  if(!aNum==12)
    cerr <<"****** ERROR,A, aNum is incorrect"<<endl;
  
  for(int i=0;i<6;i++){
    for(int j=0;j<4;j++){
      rescoef[i][j]=cXval[aNum];
      aNum++;
    }
  }
  //  Non-Res coefficients xvals of 37-44
  if(!aNum==36)
    cerr <<"****** ERROR,B, aNum is incorrect"<<endl;

  for(int i=0;i<6;i++){
    for(int j=0;j<4;j++){
      nr_coeff[i][j]=cXval[aNum];
      aNum++;
    }
  }
  
  // Begin resonance Q^2 dependence calculations   //
  // uses xvals 49
  for(int i=0;i<6;i++){
    height[i] = rescoef[i][0]
      *(1.0+rescoef[i][1]*Q2/(1.0+rescoef[i][2]*Q2))
      /pow((1.+ Q2/0.71),rescoef[i][3]);
  }
  
  dIP = 1./pow((1.+Q2/0.71),2);  
  dIP2 = dIP*dIP;
  height[6] = cXval[48]*dIP2; 
  Int_t IW_Mean = (Int_t)(1000.*sqrt(aWsq));
  sig_res = 0.0;
  for(int i=0;i<7;i++){
    sig_rsv[i] =  height[i]*sigr[i][IW_Mean];
    sig_res += sig_rsv[i];
  }

// Begin non-resonant part uses xvals 45, 46, 50
// Depends on both W2 and Q2 so can't easily precalculate
  sig_nr = 0.;
  xpr = 1.+(aWsq-pow((P_MASS+PION_MASS),2))/(Q2+cXval[49]);
  xpr = 1./xpr;
  w = sqrt(aWsq);
  wdiff = w - (P_MASS + PION_MASS);
  
  for(int i=0;i<2;i++)
    sig_nr += nr_coeff[i][0]*pow(wdiff,((2*(i+1)+1)/2.0))
      /pow((Q2+nr_coeff[i][1]),(nr_coeff[i][2]+nr_coeff[i][3]*Q2+cXval[44+i]*Q2*Q2));
  sig_nr = sig_nr * xpr;
  sig = sig_res + sig_nr; 
       
  // changed to use F1 instead of sigt
  F1 = sig * (aWsq-PM_SQ)/8.0/TMath::Pi()/TMath::Pi()/Alpha/(HC_SQ*GeVToMeV);
  sig = F1;
  xb = Q2/(aWsq+Q2-PM_SQ);
  nu = Q2/P_MASS/2.0/xb;

// Now include Aji's hack to fix f2n/f2p
// if neutron then do a low x tweak, because peter's f2n is not well constrained at low x
// we don't want to go to higher x values since we know nmc f2n/f2p is for "bound"
// these limits are chosen, so that we suffer only  minimal damage in transition                            
  Double_t xbj1=0.38;                 //x<xbj1 --> use tweaked f2n
  Double_t xbj2=0.42;                 //x1<=x<x2 --> smooth transition from above to peters
  Double_t f1_pet=F1;                 //x>=x2  --> peters's f2n
  Double_t F1p,Rc,W1n,W2n,F2,f2n_fix=0.0,f1n_fix=0.0,f1_fix=0.0,frac=0.0;

  if (xb<xbj2){
    Double_t f2nf2p_nmc = gGet_NMC(xb,Q2);
    //    temp = gGet_Christy806(Q2,aWsq);
    gGet_Christy806(Q2,aWsq,temp);
    F1p = temp->First;
    Rc = temp->Second;

    W1n = F1p / P_MASS;
    W2n = W1n*(1.0+ Rc)/(1.+nu*nu/Q2);
    F2 = nu * W2n; 
    f2n_fix=F2*f2nf2p_nmc;  //eric f2p times nmc f2n/f2p
    f1n_fix=f2n_fix*(P_MASS/nu)*(1.+nu*nu/Q2)/(1.0+ Rc);
    f1_fix = (f1n_fix+F1p)/2.0;
  }
           
  if (xb<xbj1)
    F1=f1_fix;
  else if ((xb>=xbj1)&&(xb<xbj2)){  
    frac = (xb-xbj1)/(xbj2-xbj1);
    F1 = f1_pet*frac + f1_fix*(1.-frac);
  }
  else if(xb>=xbj2)
    F1=f1_pet;
  sig = F1;
  
  delete[] cXval;
  delete temp;
  return sig;
}//inline Double_t gGet_Resmod_Hack(double aQsq,double aWsq,double* aXval)


// NMC "bound" f2n/f2p 
// ref: Nuc Phy B 371 (1992) 3
inline Double_t gGet_NMC(double aXbj, double aQsq)
{
  Double_t a_x = 0.979
    -(1.692*aXbj)
    +(2.797*aXbj*aXbj)
    -(4.313*aXbj*aXbj*aXbj)
    +(3.075*aXbj*aXbj*aXbj*aXbj);
  
  Double_t  b_x = -(0.171*aXbj)
    +(0.244*aXbj*aXbj);

  Double_t f2nf2p=pow((aQsq/20.),b_x);
  f2nf2p *= a_x*(1.+(aXbj*aXbj/aQsq));
  return f2nf2p;
}//inline Double_t gGet_Resmod_Hack(double aQsq,double aWsq)
