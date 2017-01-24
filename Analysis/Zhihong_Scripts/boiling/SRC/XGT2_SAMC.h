#include "XGT2_XEMC.h"

/*inline Double_t gGet_Rad_Corr_Bump(){{{*/
inline Double_t gGet_Rad_Corr_Bump(const TString& aTarget, const Double_t aRctPt){
  double aRad_Corr_Fact = 1.0;

  return aRad_Corr_Fact;
}
/*}}}*/

/*inline Double_t gGet_Delta_Correct(const TString& aArm,const TString& aKin, const Double_t& aMom){{{*/
inline Double_t gGet_Delta_Correct(const TString& aArm,const TString& aKin, const Double_t& aDelta)
{
  /*
     Due to the RQ3 field was scaled down, and the Delta-Matrix was not calibrated, while in the SAMC, 
	 all Delta, Phi and Theta are all updated in the HRS-R transportation functions, so the Yield of real data
	 and simulation data don't agree with each other, so an addtional correction on Yield due to the un-updated D marix is applied.
	 However, the correct way is to apply the correction on Delta in the REAL data instead. It is just a study.
     --- Z. Ye, 02/13/2013
  */
	Double_t aDelta_Corr = 1.0;
    Double_t aA = 0., aB = 0., aC = 0.;


  if(aArm == "L"){
	  //No Correction for HRS-L
	  aDelta_Corr = 1.0;
     
  }
       
  if(aArm == "R"){

	  if(aKin == "Kin3.1"){
		  aA = 0.9137; aB = 4.018; aC = 76.84;
	  }
	  else if(aKin == "Kin3.2"){
		  aA = 0.9054; aB = 6.265; aC = 65.54;
	  }
	  else if(aKin == "Kin4.1"){
		  aA = 0.9087; aB = 3.509; aC = 53.57;
	  }
	  else if(aKin == "Kin4.2"){
		  aA = 0.9091; aB = 5.933; aC = 85.99;
	  }
	  else if(aKin == "Kin5.0"){
		  aA = 0.9095; aB = 1.620; aC = 61.99;
	  }
	  else if(aKin == "Kin5.05"){
		  aA = 0.9095; aB = 1.620; aC = 61.99;
	  }
	  else if(aKin == "Kin5.1"){
		  aA = 0.9052; aB = 3.581; aC = 50.65;
	  }
	  else if(aKin == "Kin5.2"){
		  aA = 0.9114; aB = 5.992; aC = 70.39;
	  }
	  else if(aKin == "Kin6.5"){
		  aA = 0.9104; aB = 4.453; aC = 52.75;
	  }
	  else{
		  cerr<<"   ****ERROR, in Delta_Correction, Kine is not recognized!" <<endl;
		  aA = 1.0; aB = 0.0; aC = 0.0;
	  }

	  aDelta_Corr = 1.0/(aA + aB * aDelta + aC * aDelta * aDelta);
	  //In case the correction goes crazy, don't do the correction. It shouldn't! 
	  if(aDelta_Corr<0.1 || aDelta_Corr > 2.0){
		  cerr<<"   ****ERROR, in Delta_Correction, Wrong Correction Factor!" <<endl;
		  aDelta_Corr = 1.0;
	  }
	  //aDelta_Corr *= 1.33; //1.33 is a noralization factor
  }
  return aDelta_Corr;
}
/*}}}*/

/*Int_t gGet_Nf_SAMC(TChain* aTree){{{*/
Int_t gGet_Nf_SAMC(TChain* aTree,const Double_t aBin,const Double_t aBinCut, const Double_t aAccCut[9],const TString &aTarget, const TString &aKin,const Double_t aEp)
{
 /* Bin Cuts:
      |Xbj_tg -aBin| < aBinCut;
    Acc Cuts:
      |x_fp|       < aAccCut[0];
      |y_fp|       < aAccCut[1];
      |th_fp|      < aAccCut[2];
      |ph_fp|      < aAccCut[3];
      |dp_rec|     < aAccCut[4];
      |y_tg_rec|   < aAccCut[5];
      |th_tg_rec|  < aAccCut[6];
      |ph_tg_rec|  < aAccCut[7];
      |reactz_rec| < aAccCut[8];
    Weight:
      w[i] = aCS_Center / aCS[i];
 */

 Int_t aNf = 0;   
 Int_t aNevt = aTree->GetEntries();

 Double_t x_fp,y_fp,th_fp,ph_fp,dp_rec,y_tg_rec,th_tg_rec,ph_tg_rec,reactz_rec,Angle_rec;
 Double_t E_s,E_p, Theta, Qsq,Xbj;
 Bool_t aIsPassed; 

 aTree->SetBranchAddress("x_fp",         &x_fp);
 aTree->SetBranchAddress("y_fp",         &y_fp);
 aTree->SetBranchAddress("th_fp",        &th_fp);
 aTree->SetBranchAddress("ph_fp",        &ph_fp);
 aTree->SetBranchAddress("dp_rec",       &dp_rec);
 aTree->SetBranchAddress("y_tg_rec",     &y_tg_rec);
 aTree->SetBranchAddress("th_tg_rec",    &th_tg_rec);
 aTree->SetBranchAddress("ph_tg_rec",    &ph_tg_rec);
 aTree->SetBranchAddress("reactz_rec",   &reactz_rec);
 aTree->SetBranchAddress("Angle_rec",    &Angle_rec);
 aTree->SetBranchAddress("E_s",          &E_s);
 aTree->SetBranchAddress("E_p",          &E_p);
 aTree->SetBranchAddress("Theta",          &Theta);
 aTree->SetBranchAddress("Qsq",          &Qsq);
 aTree->SetBranchAddress("Xbj",          &Xbj);
 // aTree->SetBranchAddress("cs_Final",     &cs_Final);
 aTree->SetBranchAddress("IsPassed",     &aIsPassed);

 //Get Average Cross Section at the Center Angle
 Double_t aCS_Center = 0.0;
 Get_XS* aXS = new Get_XS();
 if(bXEMC_Model_New){
	 aXS->Load_Table(aTarget.Data(),aKin.Data());
	 aCS_Center = aXS->gGet_XS_Center(aEp) * gGet_Rad_Corr_Bump(aTarget.Data(),0.0); //CS for Ep with center angle
 }
 else{
	 aXS->Load_Table_ByKin(aTarget.Data(),aKin.Data());
	 aCS_Center = aXS->gGet_XS_Center_ByKin(aEp) * gGet_Rad_Corr_Bump(aTarget.Data(),0.0); //CS for Ep with center angle
 }

 cerr<<" aCS_Center = "<< aCS_Center<<endl;
 Double_t aCS = 0.0, aSum = 0.0, aWeight = 1.0;
 Int_t aCount = 0;

 for(int i=0;i<aNevt;i++){
  aTree->GetEntry(i);  
  if(fabs(x_fp)       < aAccCut[0] &&
     fabs(y_fp)       < aAccCut[1] &&
     fabs(th_fp)      < aAccCut[2] &&
     fabs(ph_fp)      < aAccCut[3] &&
     fabs(dp_rec)     < aAccCut[4] &&
     fabs(y_tg_rec)   < aAccCut[5] &&
     fabs(th_tg_rec)  < aAccCut[6] &&
     fabs(ph_tg_rec)  < aAccCut[7] &&
     fabs(reactz_rec) < aAccCut[8] &&
     fabs(Xbj    - aBin)< aBinCut &&
     aIsPassed ){

//    aCS_Center = aXS->gGet_XS_Center(E_p); //CS for Ep with center angle
	  if(bXEMC_Model_New)
		  aCS = aXS->gGet_XS(E_s,E_p,Theta); //CS for this event
	  else 
		  aCS = aXS->gGet_XS_ByKin(E_s,E_p,Theta); //CS for this event
	  //Calculate the value of weight for each event in the cuts:
     aWeight = aCS / aCS_Center * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec);
    aCount++;
	if(aCS>1e-33)
		aSum += aWeight;
  }
 }
 outlog <<"      SAMC Events in this bin before correction = " << aCount <<endl;
 cerr <<"      SAMC Events in this bin before correction = " << aCount <<endl;
 outlog <<"      SAMC Events in this bin after Correction  = "<<aSum<<endl<<endl;
 cerr <<"      SAMC Events in this bin after Correction  = "<<aSum<<endl<<endl;
 aNf = (Int_t) aSum;
 
 delete aXS;
 return aNf ; 
}
/*}}}*/

/*XGT2_VAR* gGet_XS_SAMC_SUM(TChain* aTree){{{*/
XGT2_VAR* gGet_XS_SAMC_SUM(TChain* aTree,const Double_t aBin,const Double_t aBinCut, const Double_t aAccCut[9], const TString &aArm, const TString &aTarget, const TString &aKin, const Double_t aP0, double *aBinCor)
{
 /* Bin Cuts:
      |Dp_tg -aBin| < aCut;
    Acc Cuts:
      |x_fp|       < aAccCut[0];
      |y_fp|       < aAccCut[1];
      |th_fp|      < aAccCut[2];
      |ph_fp|      < aAccCut[3];
      |dp_rec|     < aAccCut[4];
      |y_tg_rec|   < aAccCut[5];
      |th_tg_rec|  < aAccCut[6];
      |ph_tg_rec|  < aAccCut[7];
      |reactz_rec| < aAccCut[8];
    Weight:
      w[i] = aCS_Center / aCS[i];
 */
 XGT2_VAR* sum = new XGT2_VAR();

 Int_t aNevt = aTree->GetEntries();

 Double_t x_fp,y_fp,th_fp,ph_fp,dp_rec,y_tg_rec,th_tg_rec,ph_tg_rec,reactz_rec,Angle_rec;
 Double_t E_s,E_p, Theta, Qsq,Xbj;
 Bool_t aIsPassed; 

 aTree->SetBranchAddress("x_fp",         &x_fp);
 aTree->SetBranchAddress("y_fp",         &y_fp);
 aTree->SetBranchAddress("th_fp",        &th_fp);
 aTree->SetBranchAddress("ph_fp",        &ph_fp);
 aTree->SetBranchAddress("dp_rec",       &dp_rec);
 aTree->SetBranchAddress("y_tg_rec",     &y_tg_rec);
 aTree->SetBranchAddress("th_tg_rec",    &th_tg_rec);
 aTree->SetBranchAddress("ph_tg_rec",    &ph_tg_rec);
 aTree->SetBranchAddress("reactz_rec",   &reactz_rec);
 aTree->SetBranchAddress("Angle_rec",    &Angle_rec);
 aTree->SetBranchAddress("E_s",          &E_s);
 aTree->SetBranchAddress("E_p",          &E_p);
 aTree->SetBranchAddress("Theta",          &Theta);
 aTree->SetBranchAddress("Qsq",          &Qsq);
 aTree->SetBranchAddress("Xbj",          &Xbj);
 // aTree->SetBranchAddress("cs_Final",     &cs_Final);
 aTree->SetBranchAddress("IsPassed",     &aIsPassed);

 //Get Average Cross Section at the Center Angle
 Get_XS* aXS = new Get_XS();
 if(bXEMC_Model_New)
	 aXS->Load_Table(aTarget.Data(),aKin.Data());
 else
	 aXS->Load_Table_ByKin(aTarget.Data(),aKin.Data());
 
 Double_t aCS = 0.0, aCS_Center = 0.0, aSum = 0.0, aSum_Center=0.0;
 Int_t aCount = 0;

 for(int i=0;i<aNevt;i++){
  aTree->GetEntry(i);  
  if(fabs(x_fp)       < aAccCut[0] &&
     fabs(y_fp)       < aAccCut[1] &&
     fabs(th_fp)      < aAccCut[2] &&
     fabs(ph_fp)      < aAccCut[3] &&
     fabs(dp_rec)     < aAccCut[4] &&
     fabs(y_tg_rec)   < aAccCut[5] &&
     fabs(th_tg_rec)  < aAccCut[6] &&
     fabs(ph_tg_rec)  < aAccCut[7] &&
     fabs(reactz_rec) < aAccCut[8] &&
     fabs(E_p - aBin)< aBinCut &&
     aIsPassed ){
    
	  if(bXEMC_Model_New){
         aCS = aXS->gGet_XS(E_s,E_p,Theta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
         aCS_Center = aXS->gGet_XS(E_s,aBin,Theta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
	  }
	  else{
	     aCS = aXS->gGet_XS_ByKin(E_s,E_p,Theta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
	     aCS_Center = aXS->gGet_XS_ByKin(E_s,aBin,Theta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
	  }
	  aCount++;
    
	if(!(aCount%1000))
       cerr<<Form(" --- Reading SAMC Event = %d, XS = %e", aCount, aCS)<<endl;
	if(aCS>1e-33)
		aSum += aCS;
	if(aCS_Center>1e-33)
		aSum_Center += aCS_Center;
  }
 }
 outlog <<"      SAMC XS Sum in this bin = " << aSum <<endl;
 cerr   <<"      SAMC XS Sum in this bin = " << aSum <<endl;
/*
 //Apply Delta Correction on SAMC Yield! A temp. method --Z. Ye 02/13/2013
 Double_t aDelta = (aBin - aP0)/aP0;
 Double_t aDelta_Corr = gGet_Delta_Correct(aArm.Data(), aKin.Data(), aDelta);
 cerr<<Form("      Delta = %f, Correction = %f", aDelta, aDelta_Corr)<<endl;
 aSum *= aDelta_Corr;
*/
//Bin Center Correction from XS(Ep_Bin, Theta_Avg) to XS(Ep_Bin, Theta_Center) -- Z. Ye, 03/04/2012
 aBinCor[0] = aXS->gGet_XS_Center(aBin)/(aSum_Center/aCount); 

 sum->Value = aSum;
 if(aCount>0)
	 sum->Stat_Err = aSum * 1.0 / sqrt(aCount);
 else
	 sum->Stat_Err = 0.0;
 sum->Sys_Err = 0.0;

 delete aXS;
 return sum ; 
}
/*}}}*/

/*XGT2_VAR* gGet_XS_SAMC_SUM_Xbj(TChain* aTree){{{*/
XGT2_VAR* gGet_XS_SAMC_SUM_Xbj(TChain* aTree,const Double_t aBin,const Double_t aBinCut, const Double_t aAccCut[9], const TString &aArm, const TString &aTarget, const TString &aKin, const Double_t aP0, double* aBinCor)
{
 /* Bin Cuts:
      |Xbj -aBin| < aCut;
    Acc Cuts:
      |x_fp|       < aAccCut[0];
      |y_fp|       < aAccCut[1];
      |th_fp|      < aAccCut[2];
      |ph_fp|      < aAccCut[3];
      |dp_rec|     < aAccCut[4];
      |y_tg_rec|   < aAccCut[5];
      |th_tg_rec|  < aAccCut[6];
      |ph_tg_rec|  < aAccCut[7];
      |reactz_rec| < aAccCut[8];
    Weight:
      w[i] = aCS_Center / aCS[i];
 */
 XGT2_VAR* sum = new XGT2_VAR();

 Int_t aNevt = aTree->GetEntries();

 Double_t x_fp,y_fp,th_fp,ph_fp,dp_rec,y_tg_rec,th_tg_rec,ph_tg_rec,reactz_rec,Angle_rec;
 Double_t E_s,E_p, Theta, Qsq,Xbj;
 Bool_t aIsPassed; 

 aTree->SetBranchAddress("x_fp",         &x_fp);
 aTree->SetBranchAddress("y_fp",         &y_fp);
 aTree->SetBranchAddress("th_fp",        &th_fp);
 aTree->SetBranchAddress("ph_fp",        &ph_fp);
 aTree->SetBranchAddress("dp_rec",       &dp_rec);
 aTree->SetBranchAddress("y_tg_rec",     &y_tg_rec);
 aTree->SetBranchAddress("th_tg_rec",    &th_tg_rec);
 aTree->SetBranchAddress("ph_tg_rec",    &ph_tg_rec);
 aTree->SetBranchAddress("reactz_rec",   &reactz_rec);
 aTree->SetBranchAddress("Angle_rec",    &Angle_rec);
 aTree->SetBranchAddress("E_s",          &E_s);
 aTree->SetBranchAddress("E_p",          &E_p);
 aTree->SetBranchAddress("Theta",        &Theta);
 aTree->SetBranchAddress("Qsq",          &Qsq);
 aTree->SetBranchAddress("Xbj",          &Xbj);
 // aTree->SetBranchAddress("cs_Final",     &cs_Final);
 aTree->SetBranchAddress("IsPassed",     &aIsPassed);

 //Get Average Cross Section at the Center Angle
 Get_XS* aXS = new Get_XS();
 if(bXEMC_Model_New)
	 aXS->Load_Table(aTarget.Data(),aKin.Data());
 else
	 aXS->Load_Table_ByKin(aTarget.Data(),aKin.Data());
 
 Double_t aSinSQ,aEp_Bin, aCS = 0.0, aCS_Center = 0.0, aSum = 0.0, aSum_Center=0.0;
 Int_t aCount = 0;

 for(int i=0;i<aNevt;i++){
  aTree->GetEntry(i);  
  if(fabs(x_fp)       < aAccCut[0] &&
     fabs(y_fp)       < aAccCut[1] &&
     fabs(th_fp)      < aAccCut[2] &&
     fabs(ph_fp)      < aAccCut[3] &&
     fabs(dp_rec)     < aAccCut[4] &&
     fabs(y_tg_rec)   < aAccCut[5] &&
     fabs(th_tg_rec)  < aAccCut[6] &&
     fabs(ph_tg_rec)  < aAccCut[7] &&
     fabs(reactz_rec) < aAccCut[8] &&
     fabs(Xbj - aBin)< aBinCut &&
     aIsPassed ){

	  aSinSQ = pow(sin(Theta*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2
	  aEp_Bin = E_s/(1.0+4.0*E_s*aSinSQ/(2.0*PROTON_MASS*aBin)); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 

	  if(bXEMC_Model_New){
         aCS = aXS->gGet_XS(E_s,E_p,Theta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
         aCS_Center = aXS->gGet_XS(E_s,aEp_Bin,Theta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
	  }
	  else{
	     aCS = aXS->gGet_XS_ByKin(E_s,E_p,Theta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
	     aCS_Center = aXS->gGet_XS_ByKin(E_s,aEp_Bin,Theta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
	  }
	  aCount++;
    
	if(!(aCount%1000))
       cerr<<Form(" --- Reading SAMC Event = %d, XS = %e", aCount, aCS)<<endl;
	if(aCS>1e-33)
		aSum += aCS;
	if(aCS_Center>1e-33)
		aSum_Center += aCS_Center;
  }
 }
 outlog <<"      SAMC XS Sum in this bin = " << aSum <<endl;
 cerr   <<"      SAMC XS Sum in this bin = " << aSum <<endl;
/*
 //Apply Delta Correction on SAMC Yield! A temp. method --Z. Ye 02/13/2013
 Double_t aDelta = (aBin - aP0)/aP0;
 Double_t aDelta_Corr = gGet_Delta_Correct(aArm.Data(), aKin.Data(), aDelta);
 cerr<<Form("      Delta = %f, Correction = %f", aDelta, aDelta_Corr)<<endl;
 aSum *= aDelta_Corr;
*/
//Bin Center Correction from XS(Ep_Bin, Theta_Avg) to XS(Ep_Bin, Theta_Center) -- Z. Ye, 03/04/2012
 aBinCor[0] = aXS->gGet_XS_Center(aEp_Bin)/(aSum_Center/aCount); 

 sum->Value = aSum;
 if(aCount>0)
	 sum->Stat_Err = aSum * 1.0 / sqrt(aCount);
 else
	 sum->Stat_Err = 0.0;
 sum->Sys_Err = 0.0;

 delete aXS;
 return sum ; 
}
/*}}}*/

/*inline TString gGet_XGT2_Cut_SAMC(const TString& aArm,const TString& aTarget, const Double_t& aMom){{{*/
inline TString gGet_XGT2_Cut_SAMC(const TString& aArm,const TString& aTarget, const Double_t& aMom)
{
  TString Acc_Cut_Fp, Acc_Cut_Tg, Calo_Cuts, Cer_Cuts, Trigger_Cut;
  Double_t aX_Fp=0,aY_Fp=0,aTheta_Fp=0,aPhi_Fp=0;      // Focal Plane Cuts
  Double_t aDp_Tg=0,aY_Tg=0,aTheta_Tg=0,aPhi_Tg=0,aRctPt_Z=0;     // Target Plane Cuts
  
  Double_t gMom, gX_Fp, gY_Fp, gTh_Fp, gPh_Fp, gY_Tg, gTh_Tg, gPh_Tg, gDp_Tg,gRctPt_Z;
  TString gKin,gTgt;
  
  if(aArm == "L"){
      //Search the Cut Table
    ifstream InCut("Acc_Left_Table.dat");
    TString comment;
    comment.ReadLine(InCut);
    comment.ReadLine(InCut);
    while(!(InCut.eof())){
      InCut >> gMom >> gKin >> gTgt
	    >> gX_Fp >> gY_Fp >> gTh_Fp >> gPh_Fp
	    >> gDp_Tg >> gY_Tg >> gTh_Tg >> gPh_Tg >> gRctPt_Z;
//       outlog << "Mom = "<< gMom <<"  Old Mom = "<<aMom<<endl;
//       outlog << "Tgt = "<< gTgt <<"  Old Tgt = "<<aTarget.Data()<<endl;	

    if(aTarget.Data() == gTgt && (gMom - aMom) < 0.02 && (gMom - aMom) >-0.02){
	aX_Fp = gX_Fp;   aY_Fp = gY_Fp; aTheta_Fp = gTh_Fp; aPhi_Fp = gPh_Fp; 
	aDp_Tg = gDp_Tg; aY_Tg = gY_Tg; aTheta_Tg = gTh_Tg; aPhi_Tg = gPh_Tg; aRctPt_Z = gRctPt_Z;
	outlog <<"   === Applying Accpetance Cuts on SAMC Data: "<<endl
	     <<"   === Focal:   X_Fp = "<<aX_Fp <<", Y_Fp = "<<aY_Fp<<", Theta_Fp = "<<aTheta_Fp<<", Phi_Fp = "<<aPhi_Fp<<endl
	     <<"   === Target: Dp_Tg = "<<aDp_Tg<<", Y_Tg = "<<aY_Tg<<", Theta_Tg = "<<aTheta_Tg<<", Phi_Tg = "<<aPhi_Tg<<", RctPt_Z = "<<aRctPt_Z<<endl;
	cerr <<"   === Applying Accpetance Cuts on SAMC Data: "<<endl
	     <<"   === Focal:   X_Fp = "<<aX_Fp <<", Y_Fp = "<<aY_Fp<<", Theta_Fp = "<<aTheta_Fp<<", Phi_Fp = "<<aPhi_Fp<<endl
	     <<"   === Target: Dp_Tg = "<<aDp_Tg<<", Y_Tg = "<<aY_Tg<<", Theta_Tg = "<<aTheta_Tg<<", Phi_Tg = "<<aPhi_Tg<<", RctPt_Z = "<<aRctPt_Z<<endl;
	break;
      }
    }
  }
       
  if(aArm == "R"){
    //Search the Cut Table
    ifstream InCut("Acc_Right_Table.dat");
    TString comment;
    comment.ReadLine(InCut);
    comment.ReadLine(InCut);
    while(!(InCut.eof())){
      InCut >> gMom >> gKin >> gTgt
	    >> gX_Fp >> gY_Fp >> gTh_Fp >> gPh_Fp
	    >> gDp_Tg >> gY_Tg >> gTh_Tg >> gPh_Tg >> gRctPt_Z;
     
      if(aTarget.Data() == gTgt.Data() && (gMom - aMom) < 0.002 && (gMom - aMom) >- 0.002){
	aX_Fp = gX_Fp;   aY_Fp = gY_Fp; aTheta_Fp = gTh_Fp; aPhi_Fp = gPh_Fp; 
	aDp_Tg = gDp_Tg; aY_Tg = gY_Tg; aTheta_Tg = gTh_Tg; aPhi_Tg = gPh_Tg; aRctPt_Z = gRctPt_Z;
	
	outlog <<"   ^^^ Applying Accpetance Cuts on SAMC Data: "<<endl
	     <<"   ^^^^Focal:   X_Fp = "<<aX_Fp <<" Y_Fp"<<aY_Fp<<" Theta_Fp"<<aTheta_Fp<<"Phi_Fp"<<aPhi_Fp<<endl
	     <<"   ^^^^Target: Dp_Tg = "<<aDp_Tg<<" Y_Tg"<<aY_Tg<<" Theta_Tg"<<aTheta_Tg<<"Phi_Tg"<<aPhi_Tg<<", RctPt_Z = "<<aRctPt_Z<<endl;
	cerr <<"   ^^^ Applying Accpetance Cuts on SAMC Data: "<<endl
	     <<"   ^^^^Focal:   X_Fp = "<<aX_Fp <<" Y_Fp"<<aY_Fp<<" Theta_Fp"<<aTheta_Fp<<"Phi_Fp"<<aPhi_Fp<<endl
	     <<"   ^^^^Target: Dp_Tg = "<<aDp_Tg<<" Y_Tg"<<aY_Tg<<" Theta_Tg"<<aTheta_Tg<<"Phi_Tg"<<aPhi_Tg<<", RctPt_Z = "<<aRctPt_Z<<endl;
	break;
      }
    }
  }

    //Focal Plane Cuts
    Acc_Cut_Fp = Form("fabs(x_fp)<=%4.3f && fabs(y_fp)<=%4.3f && fabs(th_fp)<=%4.3f && fabs(ph_fp)<=%4.3f"
		      , aX_Fp, aY_Fp, aTheta_Fp, aPhi_Fp);
    //Target Plane Cuts
    Acc_Cut_Tg = Form("fabs(dp_rec)<=%4.3f && fabs(y_tg_rec)<=%4.3f && fabs(th_tg_rec)<=%4.3f && fabs(ph_tg_rec)<=%4.3f && fabs(reactz_rec)<%4.3f&&IsPassed",
		       aDp_Tg, aY_Tg, aTheta_Tg, aPhi_Tg, aRctPt_Z);  

  TString aCuts = Acc_Cut_Fp+"&&"+ Acc_Cut_Tg;
  return aCuts;
}
/*}}}*/

/*inline TChain* gGet_SAMC_Chain(const TString& aSAMC_File){{{*/
inline TChain* gGet_SAMC_Chain(const TString& aSAMC_File)
{
  TChain *aTree = new TChain("SAMC");
  TString aSAMC_Chain = SAMCFILES_DIR+aSAMC_File;
  aTree->Add(aSAMC_Chain.Data());
  for(int i=1;i<10;i++){
    aSAMC_Chain = aSAMC_Chain.ReplaceAll(Form("%d_xem.root",i-1),Form("%d_xem.root",i));
    outlog <<"       Add SAMC Files -> "<<aSAMC_Chain.Data()<<endl;
    cerr <<"       Add SAMC Files -> "<<aSAMC_Chain.Data()<<endl;
    aTree->Add(aSAMC_Chain.Data());
  }
  return aTree;
}/*}}}*/

/*inline Double_t gGet_XS_Born_XEMC(const TString& aTarget, const TString& aKin, const Double_t aE0, const Double_t aEp, const Double_t aTheta){{{*/
inline Double_t gGet_XS_Born_XEMC(const TString& aTarget, const TString& aKin, const Double_t aE0, const Double_t aEp, const Double_t aTheta)
{
	//Get Born Cross Section at the Center Angle
	Get_XS* aXS = new Get_XS();
	Double_t aCS_Born=0.0;
	if(bXEMC_Model_New){
		aXS->Load_Table(aTarget.Data(),aKin.Data());
		aCS_Born = aXS->gGet_XS_Born(aE0,aEp,aTheta); //CS for this event
	}
	else{
		aXS->Load_Table_ByKin(aTarget.Data(),aKin.Data());
		aCS_Born = aXS->gGet_XS_Born_ByKin(aE0,aEp,aTheta); //CS for this event
	}

	delete aXS;
	return aCS_Born ; 
}
/*}}}*/

/*inline Double_t gGet_XS_Rad_XEMC(const TString& aTarget, const TString& aKin, const Double_t aE0, const Double_t aEp, const Double_t aTheta){{{*/
inline Double_t gGet_XS_Rad_XEMC(const TString& aTarget, const TString& aKin, const Double_t aE0, const Double_t aEp, const Double_t aTheta)
{
	//Get Born Cross Section at the Center Angle
	Get_XS* aXS = new Get_XS();
	Double_t aCS_Rad = 0.0;
	if(bXEMC_Model_New){
		aXS->Load_Table(aTarget.Data(),aKin.Data());
		aCS_Rad = aXS->gGet_XS(aE0,aEp,aTheta); //CS for this event
	}
	else{
		aXS->Load_Table_ByKin(aTarget.Data(),aKin.Data());
		aCS_Rad = aXS->gGet_XS_ByKin(aE0,aEp,aTheta); //CS for this event
	}

	delete aXS;
	return aCS_Rad ; 
}
/*}}}*/
