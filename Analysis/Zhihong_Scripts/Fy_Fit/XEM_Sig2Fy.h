inline XEM_VAR2* gCal_Sig2Fy( double aE0,double aEp, double aTheta, double aXS_QE, XEM_TGT* aTarget)
{
  /*Comment Region{{{*/
  /*
    This code is to extract the F(y) and y-scaling from known Cross Sections. It is
    the reverse procedure of Fy2Sig.  
     -- Zhihong Ye, 12/26/2012
   */

  /* Original Code for Fy2Sig
    Subroutine to compute sigma (nb/sr/MeV) from Y-scaling model. This subroutine
    is derived from the subroutine FTOSIG by R. McKeown. It uses the same
    basic method, but calculates the scaling variable Y from a different
    formula, in which the separation energy is used to calculate the mass
    of the recoiling (A-1) system. F_TO_SIG also has an argument list
    different from FTOSIG, and hence is NOT interchangeable with FTOSIG.

    INPUT ARGUMENTS:
    E1:         (R*8) - Incident energy in GeV
    E2:         (R*8) - Outgoing energy in GeV
    THETA:      (R*8) - Scattering angle in rad
    A:          (R*8) - 'A' of target nucleus
    Z:          (R*8) - 'Z' of target nucleus
    M1:         (R*8) - Mass of target nucleus in GeV/c2.
    MR:         (R*8) - Mass of recoiling nucleon in GeV/c2.
    ES:         (R*8) - Separation energy in GeV.
    OUTPUT ARGUMENTS:
    SIG:        (R*8) - Calculated cross section
    Y:          (R*8) - Scaling variable Y in 1/GeV.
    J. Arrington, 3/5/98
    AUX(1) = RESOL parameter for DEEPSIG smearing (not used).
    AUX(2) = ES (Separation energy)
    AUX(3) = F(0)  parameter for F(y) function
    AUX(4) = BigB  parameter for F(y)
    AUX(5) = a     parameter for F(y)
    AUX(6) = b     parameter for F(y)
    AUX(7) = alpha parameter for F(y)
  */
  /*}}}*/  
 //Get DeltaE value 
  Double_t DeltaE = 0.0;
  bool kCCor = kFALSE;
  if(kCCor){
    if(aTarget->A==1) //  !Hydrogen
      DeltaE = 0.0;
    else if(aTarget->A==2) // !Deuterium
      DeltaE = 0.0;
    else if(aTarget->A==3) // !Helium-3
      DeltaE = 0.0; //DeltaE = 0.00085;
    else if(aTarget->A==4) // !Helium-4
      DeltaE = 0.0; // DeltaE = 0.0010;
    else if(aTarget->A==9) // !Beryllium
      DeltaE = 0.001875;
    else if(aTarget->A==12) // !Carbon
      DeltaE = 0.00292;
	else if(aTarget->A==27) // !Aluminum
		DeltaE = 0.0061;
	else if(aTarget->A==40) // !Cacium 40
		DeltaE = 0.00743;
	else if(aTarget->A==48) // !Cacium 48
		DeltaE = 0.00699;
	else if(aTarget->A==64) // !Copper
		DeltaE = 0.0102;
    else if(aTarget->A==197) // !Gold
      DeltaE = 0.0199;
    else{
      cerr<<"**** Target not in the list! Coulumb Correction give up~ "<<endl;
      DeltaE = 0.0;
    }     
  }
  Double_t FF1 = (aE0 - DeltaE)/aE0;
  aE0 += DeltaE;
  aEp += DeltaE;

  Double_t CS_Theta = cos(aTheta/2.0);
  Double_t SN_Theta = sin(aTheta/2.0);
  //Double_t TN_Theta = tan(aTheta/2.0);
  Double_t SN_SQ = SN_Theta*SN_Theta;
  //Double_t CS_SQ = CS_Theta*CS_Theta;
  
  Double_t Qsq = 4.0*aE0*aEp*SN_SQ;
  Double_t nu = aE0 - aEp;
  Double_t Qv2 = Qsq + nu*nu;
  Double_t Qv = sqrt(Qv2);
  Double_t x_local = Qsq/2.0/P_MASS/nu;
  //Double_t y = nu/aEp;
  //Double_t Wsq = -Qsq + P_MASS*P_MASS + 2.0*P_MASS*nu;
  //Double_t Elastic_Peak = aE0/(1.0+2.0*aE0*SN_SQ/aTarget_Mass);

  //Calculate Mott Cross Section
  Double_t sig_mott = CS_Theta/(274.0*aE0*SN_SQ);
  sig_mott = pow(sig_mott*0.019733,2);

  XEM_VAR2* aFyy = new XEM_VAR2();
  Double_t Fy = -11111.0;
  Double_t y_cal = gGet_Y(aE0,aEp,aTheta,aTarget);
 
  if(y_cal<-1e5){ 
    aFyy->First = y_cal;
    aFyy->Second = Fy;
    return aFyy;
  }
  Double_t fact=0.,sig_p=0.,sig_n=0.;
  Double_t dwdy = Qv/sqrt(PM_SQ+Qv2+y_cal*y_cal+2.0*y_cal*Qv);
  
   if(aTarget->A > 1.5){
    XEM_VAR2* sig_np = new XEM_VAR2();
    gGet_Sig_Bar_DF(aE0,aEp,aTheta,y_cal,0.0,sig_np);
    sig_p = sig_np->First/1000.0;
    sig_n = sig_np->Second/1000.0;
    fact = (aTarget->Z*sig_p+(aTarget->A-aTarget->Z)*sig_n)/dwdy;
    delete sig_np;
   }
  else 
    fact=0.0;

  Fy = aXS_QE / (fact*1e6);
  
  
  /*John's Correction{{{*/
  /*
	 Double_t johns_fac=1.0;
	 if (aTarget->A==3)
	 johns_fac=gGet_Max(1.,1.+y_cal*1.4*2.5);
	 else if (aTarget->A==4)
	 johns_fac=gGet_Max(1.,1.+y_cal*1.4*1.75);
	 else if ((aTarget->A==9)||(aTarget->A==12))  
	 johns_fac=gGet_Max(1.,1.+y_cal*1.4*2.5);
	 else if ((aTarget->A==56)||(aTarget->A==64))
	 johns_fac=gGet_Max(1.,1.+y_cal*1.4*3);
	 else if (aTarget->A==197)
	 johns_fac=gGet_Max(1.,1.+y_cal*1.4*4);

  //  if(!aTarget->A==2)
  //    cerr <<"***** Oh Crap!!! Target has been changed! Find it out"<<endl;

  Fy /=johns_fac;
  */
  /*}}}*/

  /*Nadia's Tail Correction{{{*/
  /*
  if(x_local>2.0)
    x_local=2.0;
  
  Double_t x_low=0.0,x_high=2.0;
  if(aTarget->A>2){
    if((aTarget->A==64)||(aTarget->A==4)||(aTarget->A==197)){
      x_low=1.2; x_high=1.4;
    }
    else{
      x_low=1.4; x_high=1.6;
    }
  }
  else if(aTarget->A==2){
    x_low=1.4;  x_high=1.45;
  }
 
  Double_t corfact=1.0;
  Double_t my_frac=1.0;
  if(x_local>=x_low){ 
    corfact=  gGet_Tail_Corr(aTarget->A,x_local);
    if(x_local<x_high){
      my_frac=(x_local-x_low)/(x_high-x_low);	      
      corfact=my_frac*corfact+(1.0-my_frac);
    }
    if(corfact>0.0&&corfact<10.0)//Just to make sure corfact does not blow up
       Fy/=corfact;
  }
*/
  /*}}}End of Tail Correction*/

  aFyy->First = y_cal;
  aFyy->Second = Fy/(FF1*FF1);

  return aFyy;
 }
