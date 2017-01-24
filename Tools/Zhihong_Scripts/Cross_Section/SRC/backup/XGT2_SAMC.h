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

		if(aKin.Contains("Kin3.1")){
			aA = 0.9137; aB = 4.018; aC = 76.84;
		}
		else if(aKin.Contains("Kin3.2")){
			aA = 0.9054; aB = 6.265; aC = 65.54;
		}
		else if(aKin.Contains("Kin4.1")){
			aA = 0.9087; aB = 3.509; aC = 53.57;
		}
		else if(aKin.Contains("Kin4.2")){
			aA = 0.9091; aB = 5.933; aC = 85.99;
		}
		else if(aKin.Contains("Kin5.0")&&!(aKin.Contains("Kin5.0"))){
			aA = 0.9095; aB = 1.620; aC = 61.99;
		}
		else if(aKin.Contains("Kin5.05")){
			aA = 0.9095; aB = 1.620; aC = 61.99;
		}
		else if(aKin.Contains("Kin5.1")){
			aA = 0.9052; aB = 3.581; aC = 50.65;
		}
		else if(aKin.Contains("Kin5.2")){
			aA = 0.9114; aB = 5.992; aC = 70.39;
		}
		else if(aKin.Contains("Kin6.5")){
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
Int_t gGet_Nf_SAMC(TChain* aTree,const Double_t aBin,const Double_t aBinCut, const Double_t aAccCut[9],const TString &aTarget, const TString &aKin,const Double_t aEp, Get_XS* aXS)
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
	Int_t aIsPassed; 

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
    //Get_XS* aXS = new Get_XS();
    //aXS->Load_Table(aTarget.Data(),aKin.Data());
    aCS_Center = aXS->gGet_XS_Center(aEp) * gGet_Rad_Corr_Bump(aTarget.Data(),0.0); //CS for Ep with center angle

	cerr<<" aCS_Center = "<< aCS_Center<<endl;
	Double_t aCS = 0.0, aSum = 0.0, aWeight = 1.0;
	Int_t aCount = 0;

	for(int i=0;i<aNevt;i++){
		aTree->GetEntry(i);  
		if(fabs(x_fp)       <= aAccCut[0] &&
				fabs(y_fp)       <= aAccCut[1] &&
				fabs(th_fp)      <= aAccCut[2] &&
				fabs(ph_fp)      <= aAccCut[3] &&
				fabs(dp_rec)     <= aAccCut[4] &&
				fabs(y_tg_rec)   <= aAccCut[5] &&
				fabs(th_tg_rec)  <= aAccCut[6] &&
				fabs(ph_tg_rec)  <= aAccCut[7] &&
				fabs(reactz_rec) <= aAccCut[8] &&
				fabs(Xbj - aBin) <= aBinCut &&
				aIsPassed==1 ){

			//    aCS_Center = aXS->gGet_XS_Center(E_p); //CS for Ep with center angle
			aCS = aXS->gGet_XS(E_s,E_p,Theta); //CS for this event
			//Calculate the value of weight for each event in the cuts:
			aWeight = aCS / aCS_Center * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec);
			aCount++;
			if(aCS>1e-33)
				aSum += aWeight;
		}
	}
	outlog <<"      SAMC Events in this bin before correction = " << aCount <<endl;
	cerr   <<"      SAMC Events in this bin before correction = " << aCount <<endl;
	outlog <<"      SAMC Events in this bin after Correction  = " <<aSum<<endl<<endl;
	cerr   <<"      SAMC Events in this bin after Correction  = " <<aSum<<endl<<endl;
	aNf = (Int_t) aSum;

	//delete aXS;
	return aNf ; 
}
/*}}}*/

/*XGT2_VAR* gGet_XS_SAMC_SUM(TChain* aTree){{{*/
XGT2_VAR* gGet_XS_SAMC_SUM(TChain* aTree,const Double_t aBin,const Double_t aBinCut, const Double_t aAccCut[9], const Double_t aVZ_Cut[2], const TString &aArm, const TString &aTarget, const TString &aKin, const Double_t aP0, const Double_t aAngle, const TString& aBin_Variable, double *aBinCor, Get_XS* aXS)
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
	Int_t aIsPassed; 

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
	aTree->SetBranchAddress("E_s",          &E_s);//GeV
	aTree->SetBranchAddress("E_p",          &E_p);//GeV
	aTree->SetBranchAddress("Theta",          &Theta);
	aTree->SetBranchAddress("Qsq",          &Qsq);
	aTree->SetBranchAddress("Xbj",          &Xbj);
	// aTree->SetBranchAddress("cs_Final",     &cs_Final);
	aTree->SetBranchAddress("IsPassed",     &aIsPassed);

    ////Get Average Cross Section at the Center Angle
	//Get_XS* aXS = new Get_XS();
	//aXS->Load_Table(aTarget.Data(),aKin.Data());

	Double_t kTheta = 0.0, aCS = 0.0, aCS_Center = 0.0, aSum = 0.0, aSum_Center=0.0;
	Int_t aCount = 0;
	
	/*Target Density{{{*/
	const double aVZ_Max = 0.1020, aVZ_Min = -0.1020;
	const double aStep = (aVZ_Max-aVZ_Min)/51.;//0.004
	double aVZ[52], aRho[52], aRho1[52],aRho2[52],aRho3[52];
	for(int i=0;i<52;i++){
		aVZ[i] = 0.0; aRho[i] = 1.0; aRho2[i] = 1.0; aRho3[i] = 1.0;
	}
	double aRho_Sum = 0;
	TString aDummy="X";

	if(aTarget=="H2"||aTarget=="He3"||aTarget=="He4"){
		TString aDensity_Profile = Form("/work/halla/e08014/disk1/yez/Ep_Cross_Section/Target/%s.rho",aTarget.Data());
		cerr<<Form("\t Reading Density Profile for %s from %s",aTarget.Data(), aDensity_Profile.Data())<<endl;
		ifstream density(aDensity_Profile.Data()); 
		density >> aDummy >> aDummy >> aDummy >> aDummy >> aDummy;
		for(int i=0;i<52;i++){
			density >> aVZ[i] >> aRho[i] >> aRho1[i] >> aRho2[i] >> aRho3[i];
			aRho_Sum += aRho3[i];
		}
		density.close();
		aRho_Sum /=52;
	}
	else
		aRho_Sum = 1.0;
	/*}}}*/
    
	/*Ep Bining{{{*/
	int aVZ_Bin = 0; double aWeight = 0.0;
	if(aBin_Variable=="Ep"){
		cerr<<"      --- Calculating SAMC Events in Ep binning"<<endl;
		for(int i=0;i<aNevt;i++){
			aTree->GetEntry(i);  
			//Put the real P0 and Angle values from the real data
				double E_p_new = aP0 * (1+dp_rec);
				double aSinA = sin(aAngle*TMath::DegToRad());
				double aCosA = cos(aAngle*TMath::DegToRad());
				double Angle_new = acos( (aCosA-ph_tg_rec*aSinA)/sqrt(1.+pow(th_tg_rec,2)+pow(ph_tg_rec,2))); 
			if(fabs(x_fp)       < aAccCut[0] &&
					fabs(y_fp)       < aAccCut[1] &&
					fabs(th_fp)      < aAccCut[2] &&
					fabs(ph_fp)      < aAccCut[3] &&
					fabs(dp_rec)     < aAccCut[4] &&
					fabs(y_tg_rec)   < aAccCut[5] &&
					fabs(th_tg_rec)  < aAccCut[6] &&
					fabs(ph_tg_rec)  < aAccCut[7] &&
					fabs(reactz_rec) < aAccCut[8] &&
					fabs(E_p_new - aBin)< aBinCut &&
					reactz_rec >= aVZ_Cut[0] &&
					reactz_rec < aVZ_Cut[1] &&
					//fabs(E_p - aBin)< aBinCut &&
					aIsPassed==1 ){

				//kTheta =  Angle_rec / TMath::DegToRad(); 
				kTheta =  Angle_new / TMath::DegToRad(); 
				//aCS = aXS->gGet_XS(E_s,E_p,kTheta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
				aCS = aXS->gGet_XS(E_s,E_p_new,kTheta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
				aCS_Center = aXS->gGet_XS(E_s,aBin,kTheta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event

				if(aTarget=="H2"||aTarget=="He3"||aTarget=="He4"){
					aVZ_Bin = (reactz_rec-aVZ_Min)/aStep;
					aWeight = aRho3[aVZ_Bin]/aRho_Sum; 
				}
				else
					aWeight = 1.0;

				if(aCS>1e-33){
					aCount++;
					aSum += aCS*aWeight;
				}
				if(!(aCount%1000))
					cerr<<Form(" --- Reading SAMC Event = %d, XS = %e", aCount, aCS)<<endl;
				if(aCS_Center>1e-33)
					aSum_Center += aCS_Center*aWeight;
			}
		}
		//Bin Center Correction from XS(Ep_Bin, Theta_Avg) to XS(Ep_Bin, Theta_Center) -- Z. Ye, 03/04/2012
		aBinCor[0] = aXS->gGet_XS_Center(aBin)/(aSum/aCount)*aWeight; 
		//aBinCor[0] = aXS->gGet_XS_Center(aBin)/(aSum_Center/aCount); 
	}
    /*End of Ep Bining}}}*/
    
    /*Xbj Bining{{{*/
	else if(aBin_Variable=="xbj"){
	    Double_t aSinSQ = 0.0, aEp_Bin = 0.0;
		for(int i=0;i<aNevt;i++){
			aTree->GetEntry(i);  
			//Put the real P0 and Angle values from the real data
			double E_p_new = aP0 * (1+dp_rec);
			double aSinA = sin(aAngle*TMath::DegToRad());
			double aCosA = cos(aAngle*TMath::DegToRad());
			double Angle_new = acos( (aCosA-ph_tg_rec*aSinA)/sqrt(1.+pow(th_tg_rec,2)+pow(ph_tg_rec,2))); 
			double Q2_new = 4.0*E_s*E_p_new*pow( sin(Angle_new/2.), 2);
			double xbj_new = Q2_new/2./TARGET_MASS_A/(E_s-E_p_new);

			if(fabs(x_fp)       < aAccCut[0] &&
					fabs(y_fp)       < aAccCut[1] &&
					fabs(th_fp)      < aAccCut[2] &&
					fabs(ph_fp)      < aAccCut[3] &&
					fabs(dp_rec)     < aAccCut[4] &&
					fabs(y_tg_rec)   < aAccCut[5] &&
					fabs(th_tg_rec)  < aAccCut[6] &&
					fabs(ph_tg_rec)  < aAccCut[7] &&
					fabs(reactz_rec) < aAccCut[8] &&
					fabs(xbj_new - aBin)< aBinCut &&
					//fabs(Xbj - aBin)< aBinCut &&
					aIsPassed==1 ){

				aVZ_Bin = (reactz_rec-aVZ_Min)/aStep;
				aWeight = aRho3[aVZ_Bin]/aRho_Sum; 

				aSinSQ = pow(sin(Angle_new/2.0),2); //Define Sin(Theta/2.)^2
				aEp_Bin = E_s/(1.0+4.0*E_s*aSinSQ/(2.0*TARGET_MASS_A*aBin)); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 

				//kTheta =  Angle_rec / TMath::DegToRad(); 
				kTheta =  Angle_new / TMath::DegToRad(); 
				//aCS = aXS->gGet_XS(E_s,E_p,kTheta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
				aCS = aXS->gGet_XS(E_s,E_p_new,kTheta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
				aCS_Center = aXS->gGet_XS(E_s,aEp_Bin,kTheta) * gGet_Rad_Corr_Bump(aTarget.Data(),reactz_rec); //CS for this event
				
				if(aCS>1e-33){
					aSum += aCS*aWeight;
					aCount++;
				}
				if(!(aCount%1000))
					cerr<<Form(" --- Reading SAMC Event = %d, XS = %e", aCount, aCS)<<endl;
				if(aCS_Center>1e-33)
					aSum_Center += aCS_Center*aWeight;
			}
		}
		Double_t aTheta0 = 0.0;
		/*Get Theta0{{{*/
		if(aKin.Contains("3.1")||aKin.Contains("3.2"))
			aTheta0 = 21.0;
		else if(aKin.Contains("4.1")||aKin.Contains("4.2"))
			aTheta0 = 23.0;
		else if(aKin.Contains("5.0")||aKin.Contains("5.1")||aKin.Contains("5.05")||aKin.Contains("5.2"))
			aTheta0 = 25.0;
		else if(aKin.Contains("6.5"))
			aTheta0 = 28.0;
		else 
			cerr<<"***ERROR***, in XGT2_SAMC.h, ln398, unknown Kinematic Setting!"<<endl;
		/*}}}*/
		aTheta0 = aAngle;//Use the real data scattering angle

		aSinSQ = pow(sin(aTheta0*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2
        aEp_Bin = E_s/(1.0+4.0*E_s*aSinSQ/(2.0*TARGET_MASS_A*aBin)); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 

		//Bin Center Correction from XS(Ep_Bin, Theta_Avg) to XS(Ep_Bin, Theta_Center) -- Z. Ye, 03/04/2012
		aBinCor[0] = aXS->gGet_XS_Center(aEp_Bin)/(aSum/aCount)*aWeight; 
		//aBinCor[0] = aXS->gGet_XS_Center(aEp_Bin)/(aSum_Center/aCount); 
	}
	else{
		outlog <<"      *** ERROR, I am not able to calculate the SAMC SUM" <<endl;
		cerr   <<"      *** ERROR, I am not able to calculate the SAMC SUM" <<endl;
	}
    /*End of Xbj Bining}}}*/

	outlog <<"      SAMC XS Sum in this bin = " << aSum <<endl;
	cerr   <<"      SAMC XS Sum in this bin = " << aSum <<endl;

	sum->Value = aSum;
	if(aCount>0)//Swape Stat Err and Sys Err, Z. Ye 06/12/2014
		sum->Sys_Err = aSum * 1.0 / sqrt(aCount);
	else
		sum->Sys_Err = 0.0;
	sum->Stat_Err = 0.0;

	//delete aXS;
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
	Acc_Cut_Fp = Form("abs(x_fp)<=%4.3f && abs(y_fp)<=%4.3f && abs(th_fp)<=%4.3f && abs(ph_fp)<=%4.3f"
			, aX_Fp, aY_Fp, aTheta_Fp, aPhi_Fp);
	//Target Plane Cuts
	Acc_Cut_Tg = Form("abs(dp_rec)<=%4.3f && abs(y_tg_rec)<=%4.3f && abs(th_tg_rec)<=%4.3f && abs(ph_tg_rec)<=%4.3f && abs(reactz_rec)<%4.3f&&IsPassed",
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
inline Double_t gGet_XS_Born_XEMC(const TString& aTarget, const TString& aKin, const Double_t aE0, const Double_t aEp, const Double_t aTheta, Get_XS* aXS)
{
	//Get Born Cross Section at the Center Angle
	Double_t aCS_Born=0.0;
	//Get_XS* aXS = new Get_XS();
	//aXS->Load_Table(aTarget.Data(),aKin.Data());
	aCS_Born = aXS->gGet_XS_Born(aE0,aEp,aTheta); //CS for this event

	//delete aXS;
	return aCS_Born ; 
}
/*}}}*/

/*inline Double_t gGet_XS_Rad_XEMC(const TString& aTarget, const TString& aKin, const Double_t aE0, const Double_t aEp, const Double_t aTheta){{{*/
inline Double_t gGet_XS_Rad_XEMC(const TString& aTarget, const TString& aKin, const Double_t aE0, const Double_t aEp, const Double_t aTheta)
{
	//Get Born Cross Section at the Center Angle
	Get_XS* aXS = new Get_XS();
	Double_t aCS_Rad = 0.0;
	aXS->Load_Table(aTarget.Data(),aKin.Data());
	aCS_Rad = aXS->gGet_XS(aE0,aEp,aTheta); //CS for this event

	delete aXS;
	return aCS_Rad ; 
}
/*}}}*/
