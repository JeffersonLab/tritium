inline void gCal_Fy2Sig( double aE0,double aEp, double aTheta, XEM_TGT* aTarget, XEM_SIG* fy2sig)
{
  /*Info{{{*/
  /*
    Subroutine to compute sigma (nb/sr/MeV) from Y-scaling model. This subroutine
    is derived from the subroutine FTOSIG by R. McKeown. It uses the same
    basic method, but calculates the scaling variable Y from a different
    formula, in which the separation energy is used to calculate the mass
    of the recoiling (A-1) system. F_TO_SIG also has an argument list
    different from FTOSIG, and hence is NOT interchangeable with FTOSIG.

    INPUT ARGUMENTS:
    E1:         (R*8) - Incident energy in GeV
    E2:         (R*8) - Outgoing energy in GeV
    THETA:      (R*8) - Scattering angle in Rad 
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
  fy2sig->Value=0.0;
  fy2sig->Factor=0.0;

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

  Double_t y_cal = gGet_Y(aE0,aEp,aTheta,aTarget);
  
  if(y_cal<-1e5){
    fy2sig->Value=0.0;
    fy2sig->Factor=0.0;
    return;
  }
  Double_t fact=0.,sig_p=0.,sig_n=0.,sig_out=0.;
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

  const Double_t Xbj_Max_He3 = 2.85; 
  //if(aTarget->A==3 &&aTarget->Z==2 && x_local >= Xbj_Max_He3)
  if(aTarget->A==3 &&aTarget->Z==2) // New twist for He3 only, -- Z. Ye 05/19/2014
	  sig_out = fact*gGet_He3_Fy(y_cal,aTarget)*1e6;//Add a flat tail for He3 at y<0.80
  else
	  sig_out = fact*gGet_Fy(y_cal,aTarget)*1e6;

  /*Old Correction not in used{{{*/
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
  
//  sig_out*=johns_fac;
  
  if(x_local>2.0)
    x_local=2.0;
  
  Double_t x_low=0.0,x_high=2.0;
  if(aTarget->A>2){
    if((aTarget->A==64)||(aTarget->A==4)||(aTarget->A==197)){
      x_low=1.2;
      x_high=1.4;
    }
    else{
      x_low=1.4;
      x_high=1.6;
    }
  }
  else  if(aTarget->A==2){
    x_low=1.4;
    x_high=1.45;
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
      sig_out*=corfact;
  }
*/
  /*}}}*/

  fy2sig->Value = sig_out;
  fy2sig->Factor = fact;
  return;
 }

//Calculate Y-Scaling variable
inline Double_t gGet_Y(double vE0, double vEp, double vTheta, XEM_TGT* vTarget)
{
  Double_t nu = vE0 - vEp;
  Double_t SN_Theta = sin(vTheta/2.0);
  Double_t SN_SQ = SN_Theta*SN_Theta;
  Double_t Qsq = 4.0*vE0*vEp*SN_SQ; //4-momentum
  Double_t Qv2 = Qsq + nu*nu; //3-Momentum
  Double_t Qv = sqrt(Qv2);

  //Solving y
  Double_t Mass_A = vTarget->Mass;

  Double_t y_cal=0.0;
  Double_t Mrec = P_MASS;
  Double_t Mass_A_1 = Mass_A + vTarget->ESep - Mrec; //Mass of (A-1) system
  Double_t temp, coeff_a, coeff_b, coeff_c,root2;
  
  temp = Qv2 + Mrec*Mrec - Mass_A_1*Mass_A_1 - (Mass_A+nu)*(Mass_A+nu);
  coeff_a = 4.0*(Qv2-(Mass_A+nu)*(Mass_A+nu));
  coeff_b = 4.0*Qv*temp;
  coeff_c = temp*temp-4.0*Mass_A_1*Mass_A_1*(Mass_A+nu)*(Mass_A+nu);
  root2 = coeff_b*coeff_b - 4.0*coeff_a*coeff_c;

  if(root2<0.0||coeff_a==0)
    return -1e6;  //No Solution
  else
    y_cal = -(coeff_b+sqrt(root2))/(2.0*coeff_a);

  return y_cal;    
}

inline Double_t gGet_Fy(double vY, XEM_TGT* vTarget)
{
  Double_t fy,y,f0,B,alpha,a,b;
  y=1000.0*vY;
  f0 = vTarget->f0;
  B = vTarget->B;
  a = vTarget->a;
  b = vTarget->b;
  alpha = vTarget->alpha;

  if(vTarget->A==2)
    fy = (f0-B)*alpha*alpha*exp(-(a*y)*(a*y))/(alpha*alpha+y*y)+B*exp(-b*abs(y));
  else
    fy = (f0-B)*alpha*alpha*exp(-(a*y)*(a*y))/(alpha*alpha+y*y)+B*exp(-(b*y)*(b*y));

  fy*=1000.0;
 
  return fy;
}

inline Double_t gGet_He3_Fy( double vY, XEM_TGT* vTarget){
  //This is a special function for He3 only to prevent the He3 Fy falls to zero quickly when x->3
  // -- Z. Ye, 05/19/2014
  Double_t fy,y,f0,B,alpha,a,b;
  f0 = vTarget->f0;
  B = vTarget->B;
  a = vTarget->a;
  b = vTarget->b;
  alpha = vTarget->alpha;
 
  const double Y_Max_He3 = -0.8; //-0.8 roughly corresponds to xbj->3
  if(vY>Y_Max_He3)
	  y=1000.0*vY;//If it is still not at the tail, use the right distribution
  else
	  y=1000.0*Y_Max_He3; //Add a flat line after y<=-0.8 to prevend the Fy falls into zero quickly

  fy = 1000.0*(f0-B)*alpha*alpha*exp(-(a*y)*(a*y))/(alpha*alpha+y*y)+B*exp(-(b*y)*(b*y));

  return fy;
}

inline void gGet_Sig_Bar_DF(double vE0,double vEp,double vTheta, double vY,double vPt,XEM_VAR2* sig_pn)
{
  Double_t CS_Theta = cos(vTheta/2.0);
  Double_t SN_Theta = sin(vTheta/2.0);
  Double_t TN_Theta = tan(vTheta/2.0);
  Double_t SN_SQ = SN_Theta*SN_Theta;
  // Double_t CS_SQ = CS_Theta*CS_Theta;
  Double_t TN_SQ = TN_Theta*TN_Theta;
  
  Double_t Qsq = 4.0*vE0*vEp*SN_SQ; //4-momentum
  Double_t nu = vE0 - vEp;
  Double_t Tau =Qsq/(4.0*P_MASS*P_MASS);
  Double_t Qv2 = Qsq + nu*nu; //3-Momentum
  Double_t Qv = sqrt(Qv2);
  Double_t vYp = vY*vY + vPt*vPt;

  //Calculate Mott Cross Section, Final Energy of Struck nulceon
  Double_t sig_mott =  HC_SQ*pow(Alpha*CS_Theta,2)/pow(2.0*vE0*SN_SQ,2);
  Double_t En_Final = sqrt(PM_SQ+Qv2+vYp+2.0*Qv*vY); 
  Double_t E_Bar = sqrt(vYp+PM_SQ);
  Double_t Qsq_Bar = Qv2 - pow((En_Final-E_Bar),2);
  Double_t Tau_Bar = Qsq_Bar/(4.0*PM_SQ);

  Double_t GFF[4]={0.0,0.,0.,0.}; //GEP,GEN,GMP,GMN
  Double_t gep,gen,gmp,gmn,tmp;

  //NOTE: It need to be discussed which form factors we need to use. There are many options.
  gGet_FF(11,Qsq,GFF);
  tmp = GFF[2];
  gGet_FF(0,Qsq,GFF);
  GFF[2] = tmp;

  gep=GFF[0]; gen=GFF[1];
  gmp=GFF[2]; gmn=GFF[3];

  Double_t f1p = (gep + Tau*gmp)/(1 + Tau);
  Double_t f1n = (gen + Tau*gmn)/(1 + Tau);
  Double_t f2p = (gmp - gep)/(1 + Tau);
  Double_t f2n = (gmn - gen)/(1 + Tau);
	
  // DeForest cross sections.
  Double_t Qrat = Qsq/Qv2;
  Double_t t1 = Qsq_Bar*TN_SQ/2.0 + 0.25*Qrat*(Qsq_Bar - Qsq);
  Double_t t2 = 0.25*Qrat*Qrat*(E_Bar + En_Final)*(E_Bar + En_Final) + (Qrat/2 + TN_SQ)*vPt*vPt;
  Double_t sig_p = t1*pow((f1p + f2p),2) + t2*(f1p*f1p + Tau_Bar*f2p*f2p);
  Double_t sig_n = t1*pow((f1n + f2n),2) + t2*(f1n*f1n + Tau_Bar*f2n*f2n);
   
  sig_p = sig_mott*sig_p/(En_Final*E_Bar);
  sig_n = sig_mott*sig_n/(En_Final*E_Bar);
 
  sig_pn->First = sig_p;
  sig_pn->Second = sig_n;
  return;
}

inline void gGet_FF(int vFlag, double vQsq, double* vFF)
{

/* Converted from nform_xem.f
C  ----------------------------------------------------------
C
C   CALCULATE NUCLEON FORM FACTORS
C
C   IG = 1 - DIPOLE WITH FORM FACTOR SCALING AND GEN=0.0
C      = 2 - IJL FIVE PARAMETER MODEL FIT
C        3 - GEP AND GMP FROM IJL, GMN=DIPOLE, GEN=GALSTER
C            WE CALL THIS THE "BEST FIT" NUCLEON FORM FACTORS
C        4 - BEST FIT EXCEPT GEN = 0.0
C        5 - BLATNIK + ZOVKO VDM FIT
C        6 - JANNSENS 1966 STANDARD FIT
C        7 - DIPOLE + F1N = 0.0
C        8 - GES = 0.0,  GMS = 1.0
C        9 - GES = 1.0,  GMS = 0.0
C       10 - HOHLER1 - PROTON AND NEUTRON FROM FIT 8.2
C       11 - HOHLER2 - PROTON FROM FIT 5.3, NEUTRON FROM 8.2
C       12 - GARI + KRUMPELMANN, Z PHYS. A322,689(1985)
C       13 - KORNER + KURODA, PHYS. REV. D16,2165(1977)
C       14 - JRA fit for proton, Kelly (PRC 70 068202) for neutron
C   QQG = INPUT Q SQUARED (GEV**2)
C
C  ------------------------------------------------------------
*/
  
//IJL PARAMETERS FOR 5 PARAMETER DIPOLE FIT (IN GEV UNITS), PHYS LETT. 43B, 191(1973)
  Double_t GAM=0.25, BR=0.672, BW=1.102, BF=0.112,AF=-0.052;
  Double_t RMN2=0.8817,RMW2=0.6146, RMF2=1.0384, RMR2=0.5852, GAMR=0.112 , PI=3.14159;
  Double_t RMPI= 0.139,  RMPI2=0.019321;
  //PARAMETERS FOR BLATNIK AND ZOVKO VDM FIT,ACTA PHYSICA AUSTRIACA 39, 62(1974)
  Double_t TRO=0.585, TROP =1.30, TROPP=2.10, TFI=1.039, TOM=0.614 , TOMP=1.40;
  //FITTED PARAMETERS
  Double_t RMUS=-0.060, RMUV=1.853, BS=-0.91, BV=-1.10;
  //PARAMETERS FOR GARI AND KRUMPELMANN, Z.PHYS.A322,689(1985)
  Double_t RMRHO2=0.6022, CRHO=0.377, RKRHO=6.62, RKS=-0.12, RMOMG2=0.6147, COMG=0.411, RKOMG=0.163, RKV=3.706;
  Double_t RLAM12=0.632, RLAM22=5.153, RLAMQ2=0.0841;
  //PARAMETERS FOR KORNER AND KURODA
  Double_t VRH1=0.593, VRH2=1.593, VRH3=2.593, VOM1=0.614, VOM2=1.614, VOM3=2.614;
  Double_t RMUP=2.792782, RMUN=-1.913148;
  
  //CONVERT TO FM**-2
  Double_t QQG = vQsq;
  Double_t QQ  = QQG/pow(0.197328,2);
  Double_t TAU = QQG/(4.0*RMN2);
 
  Double_t GEP=0.,GEN=0.,GMP=0.,GMN=0.;
  Double_t GT,T1,T2,ALPH,TOP,BOT,RHO,F1S,F1V,F2S,F2V,GD;
  Double_t RS,RV,F1E,F2E,F1M,F2M;
  Double_t F1,F2,F3,GES,GMS,GEV,GMV;
  Double_t F1RHO,F2RHO,F1P,F2P;
  Double_t QQP,C1,C2,C3,C4,F2VK,F2SK,F1N,F2N;
  // DIPOLE
  if(vFlag==0){
    GEP = 1./pow((1.0 + QQG/0.71),2);
    GEN = 0.0;
    GMP = RMUP*GEP;
    GMN = RMUN*GEP;
  }

  //IJL 5 PARAMTER JOB
  if(vFlag==1||vFlag==2||vFlag==3){ 
    GT = 0.5/pow((1. + GAM*QQG),2);
    T1 = sqrt(QQG + 4.0*RMPI2);
    T2 = sqrt(QQG);
    ALPH = 2.*T1*log((T1 + T2)/(2.0*RMPI))/(T2*PI);
    TOP = RMR2 + 8.0*GAMR*RMPI/PI;
    BOT = RMR2 + QQG + (4.0*RMPI2 + QQG)*GAMR*ALPH/RMPI;
    RHO = TOP/BOT;
    F1S = GT*((1.0 - BW - BF) + BW/(1.0 + QQG/RMW2)
		       + BF/(1.0 + QQG/RMF2));
    F1V = GT*((1.0 - BR) + BR*RHO);
    F2S = GT*((-0.12 - AF)/(1.0 + QQG/RMW2)
		       + AF/(1. + QQG/RMF2));
    F2V = GT*(3.706*RHO);
    
    GEP = F1V + F1S - TAU*(F2V + F2S);
    GEN = F1S - F1V - TAU*(F2S - F2V);
    GMP = F1V + F1S + F2V + F2S;
    GMN = F1S - F1V + F2S - F2V;

    if(vFlag==1||vFlag==3){
      GD = 1./pow((1.+QQG/.71),2);
      GMN = RMUN*GD;
      GEN = -RMUN*TAU*GD/(1.0 + 5.6*TAU);;
      if(vFlag==1)
	GEN = 0.0;
    }
  }

  //BLATNIK AND ZOVKO
  if(vFlag == 4){
    RS = 1./((1.+QQG/TOM)*(1.+QQG/TFI)*(1.+QQG/TOMP));
    RV = 1./((1.+QQG/TRO)*(1.+QQG/TROP)*(1.+QQG/TROPP));
    F1E = (0.5-TAU*(RMUS+2.*RMN2*BS))*RS;
    F2E = (0.5-TAU*(RMUV+2.*RMN2*BV))*RV;
    F1M = (0.5+RMUS-0.5*BS*QQG)*RS;
    F2M = (0.5+RMUV-0.5*BV*QQG)*RV;
    GEP = F1E + F2E;
    GMP = F1M + F2M;
    GEN = F1E - F2E;
    GMN = F1M - F2M;
  }

  //JANNSSENS
  if(vFlag==5){
    F1 = 1. + QQ/15.7;
    F2 = 1. + QQ/26.7;
    F3 = 1. + QQ/8.19;
    GES = 0.5*(2.5/F1 - 1.6/F2 + 0.10);
    GMS = 0.44*(3.33/F1 - 2.77/F2 + 0.44);
    GEV = 0.5*(1.16/F3 - 0.16);
    GMV = 2.353*(1.11/F3 - 0.11);
    GEP = GES + GEV;
    GMP = GMS + GMV;
    GEN = GES - GEV;
    GMN = GMS - GMV;
  }

  //DIPOLE + F1N = 0.0
  if(vFlag==6){
    GEP = 1./pow((1. + QQG/0.71),2);
    GEN = -RMUN*TAU*GEP;
    GMP = RMUP*GEP;
    GMN = RMUN*GEP;
  }

  if(vFlag==7){
    GEP = 0.0;
    GEN = 0.0;
    GMP = 1.0;
    GMN = 0.0;
  }

  if(vFlag==8){
    GEP = 1.0;
    GEN = 0.0;
    GMP = 0.0;
    GMN = 0.0;
  }

  //HOHLER1 AND HOHLER2
  if(vFlag==9||vFlag==10){
    F1RHO = 0.5*(0.955 + 0.090/pow((1. + QQG/0.355),2))/
      (1. + QQG/0.536);
    F2RHO = 0.5*(5.335 + 0.962/(1. + QQG/0.268))/
      (1. + QQG/0.603);
    F1S = 0.71/(0.6129 + QQG) - 0.64/(1.0404 + QQG)
      - 0.13/(3.240 + QQG);
    F1V = F1RHO + 0.05/(1.464 + QQG) - 0.52/(6.0025 + QQG)
      + 0.28/(8.7025 + QQG);
    F2S = -0.11/(0.6129 + QQG) + 0.13/(1.0404 + QQG)
      - 0.02/(3.240 + QQG);
    F2V = F2RHO - 1.99/(1.464 + QQG) + 0.20/(6.0025 + QQG)
      + 0.19/(8.7025 + QQG);
    
    GEP = F1V + F1S - TAU*(F2V + F2S);
    GEN = F1S - F1V - TAU*(F2S - F2V);
    GMP = F1V + F1S + F2V + F2S;
    GMN = F1S - F1V + F2S - F2V;
 
    if(vFlag==9){
      //   HOHLER2 - USE PROTON FIT 5.3
      F1P = F1RHO + 0.67/(0.6129 + QQG) - 0.39/(0.9216 + QQG)
	-0.54/(2.7556 + QQG);
      F2P = F2RHO + 0.04/(0.6129 + QQG) - 1.88/(1.2996 + QQG)
	+ 0.24/(10.1761 + QQG);
      GEP = F1P - TAU*F2P;
      GMP = F1P + F2P;
    }
  }

  // GARI AND KRUMPELMANN
  if(vFlag==11){
    QQP = QQG*log(((RLAM22+QQG)/RLAMQ2))/log(RLAM22/RLAMQ2);
    C1 = RLAM12/(RLAM12 + QQP);
    C2 = RLAM22/(RLAM22 + QQP);
    F1 = C1*C2;
    F2 = F1*C2;
    C3 = RMRHO2/(RMRHO2 + QQG);
    C4 = RMOMG2/(RMOMG2 + QQG);
    F1V = (C3*CRHO + (1-CRHO))*F1;
    F2VK = (C3*CRHO*RKRHO + (RKV-CRHO*RKRHO))*F2;
    F1S = (C4*COMG + (1-COMG))*F1;
    F2SK = (C4*COMG*RKOMG + (RKS-COMG*RKOMG))*F2;
    F1P = 0.5*(F1S + F1V);
    F1N = 0.5*(F1S - F1V);
    F2P = 0.5*(F2SK + F2VK);
    F2N = 0.5*(F2SK - F2VK);
    GEP = F1P - TAU*F2P;
    GMP = F1P + F2P;
    GEN = F1N - TAU*F2N;
    GMN = F1N + F2N;
  }

  // KORNER AND KURODA
  if(vFlag==12){
    F1S = (1/(1+QQG/VOM1))*(1/(1+QQG/VOM2));
    F1V = (1/(1+QQG/VRH1))*(1/(1+QQG/VRH2));
    F2S = F1S*(1/(1+QQG/VOM3));
    F2V = F1V*(1/(1+QQG/VRH3));
    F1P = 0.5*F1S + 0.5*F1V;
    F1N = 0.5*F1S - 0.5*F1V;
    F2P = (RMUP-1)*(-0.0335*F2S + 1.0335*F2V);
    F2N =    -RMUN*(-0.0335*F2S - 1.0335*F2V);
    GEP = F1P - TAU*F2P;
    GMP = F1P + F2P;
    GEN = F1N - TAU*F2N;
    GMN = F1N + F2N;
  }

  // JRA for proton, KELLY for neutron
  Double_t p2,p4,p6,p8,p10,p12,q2,q4,q6,q8,q10,q12;
  if(vFlag==13){
    p2 = 3.19;
    p4 = 1.355;
    p6 = 0.151;
    p8 = -0.114*1e-01;
    p10 = 0.533*1e-03;
    p12 = -0.900*1e-05;
    q2 = 3.226;
    q4 = 1.508;
    q6 = -0.3773;
    q8 = 0.611;
    q10 =-0.1853;
    q12 = 0.1596*1e-01;
    
    Double_t Q  = sqrt(QQG);
    GMP=RMUP/(1. + p2*pow(Q,2) + p4*pow(Q,4) + p6*pow(Q,6) + p8*pow(Q,8) 
	      + p10*pow(Q,10) + p12*pow(Q,12));
        
    GEP=1./(1. + q2*pow(Q,2) + q4*pow(Q,4) + q6*pow(Q,6) + q8*pow(Q,8)
	    +q10*pow(Q,10) + q12*pow(Q,12));

    Double_t a1 = 2.33;
    Double_t b1=14.72;
    Double_t b2=24.20;
    Double_t b3=84.10;
    Double_t aa = 1.70;
    Double_t bb=3.30;
    Double_t Ldip=0.71;
    
    GMN = RMUN*(1.0 + a1*TAU)/(1.0+b1*TAU+b2*pow(TAU,2)+b3*pow(TAU,3));
    GEN = aa*TAU/(1.0+bb*TAU)*1.0/pow((1.0+QQG/Ldip),2);
  }

  vFF[0] = GEP; 
  vFF[1] = GEN;
  vFF[2] = GMP;
  vFF[3] = GMN;
  return;
}//inline void gGet_FF(double vFlag, double vQsq, double* vFF)


inline Double_t gGet_Tail_Corr(Int_t kA,Double_t kX)
{
   Double_t aa=1.;
   Double_t bb=0.;
   Double_t cc=0.;
   Double_t dd=0.;
   Double_t ee=0.;
   Double_t ff=0.;
   
      if(kA==2){
        aa = 1.72816025139459;
        bb = 2.53114253906281;
        cc = -2.72505067059703;
        dd = -1.58637090989274;
        ee = -16.3972900673533;
      }          
      else if(kA==3){
        aa = 0.552199789237622;
        bb = 0.8;
        cc = 0.06864880408328;
        dd = -0.320972192919132;
        ee = 0;
      }        
      else if(kA==4){
        bb = 0.466102123780537;
        cc = 0.0156369553828335;
        dd = -0.122243059123825;
        aa = 0.682462282515971;
      }
      else if(kA==9){
        bb = 0.463011918692135;
        cc = 0.0125252624899601;
        dd = -0.101843839213646;
        aa = 0.674455752091906;
      }
      else if (kA==12){
        bb = 0.222436834975864;
        cc = 0.00769270345172033;
        dd = -0.060282702596254;
        aa = 0.840262866196151;
      }
      else if((kA==40)||(kA==48)||kA==27){ //FIX_HERE
        bb = 0.041323394008416;
        cc = 0.00447016532137148;
        dd = -0.0303635977582275;
        aa = 1.00675406673173;
      }
      else if((kA==63)||(kA==64)||(kA==56)){
        bb = 0.041323394008416;
        cc = 0.00447016532137148;
        dd = -0.0303635977582275;
        aa = 1.00675406673173;
      }
      else if (kA==197){
        bb = 0.0667337559531751;
        cc = 0.00448892579200859;
        dd = -0.0334460588480325;
        aa = 0.981274819686673;
      }
      else{
        cerr<<Form(" *** Tail Correction, can not find the target A=%d",kA)<<endl;
      }
     Double_t tail_cor=aa*exp(bb*kX)+cc*pow(kX,6)+dd*pow(kX,4)+ee*pow(kX,2)+ff;

     return tail_cor;
}//inline Double_t gGet_Tail_Corr(Int_t kA,Double_t kX)

