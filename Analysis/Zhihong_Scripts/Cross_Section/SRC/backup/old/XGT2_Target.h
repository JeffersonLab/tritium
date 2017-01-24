/* 
Update: Add Boiling Effects along Z_react, in 01/09/2013
*/
/*Boiling Correction{{{*/
double Boiling_Corr(const int aTarget_A, const double aVZ, const double aCurrent, const double aRho){
	TString aTarget_Name;
	if(aTarget_A==2)
		aTarget_Name = "H2";
	else if(aTarget_A==3)
		aTarget_Name = "He3";
	else if(aTarget_A==4)
		aTarget_Name = "He4";
	else
		cerr<<"   ***ERROR, in XGT2_Target.h, Line=14"<<endl;

	TString Boiling_File=Form("%s/Target/%s_Boiling.out",XGT2_DIR.Data(), aTarget_Name.Data());
	const int tBin = 60; //Bin Number used in Boiling Study, all equal to 60
	const double step = 24.0/tBin; //+- 12cm used on VZ cut
	double *tVZ = new double[tBin];
	double *tSlop = new double[tBin];

	ifstream infile(Boiling_File); 
	TString aCom; infile >> aCom >> aCom; //Reading the first line

	for(int i=0;i<tBin;i++){
		infile >> tVZ[i] >> tSlop[i];
		tVZ[i]*=100.0; //Convert from meter to cm
	}
	infile.close();

	int bin = (int) ((aVZ + 12.0)/step);
	double slop = ( tSlop[bin+1] - (tSlop[bin+1]-tSlop[bin])/(tVZ[bin+1]-tVZ[bin]) * (tVZ[bin+1]-aVZ) );
	double rho = aRho * (1 + slop*aCurrent);

//	cerr<<Form("Bin=%d, Slop=%f, VZ=%f/%f, I=%f, Rho0=%f, Rho=%f",bin,slop,aVZ,tVZ[bin], aCurrent, aRho, rho)<<endl;

	delete tVZ; delete tSlop;
	return rho;
}

/*}}}*/

/*Correction-Inverse{{{*/
double Boiling_Corr_Inverse(const int aTarget_A, const double aVZ, const double aCurrent, const double aRho){
	TString aTarget_Name;
	if(aTarget_A==2)
		aTarget_Name = "H2";
	else if(aTarget_A==3)
		aTarget_Name = "He3";
	else if(aTarget_A==4)
		aTarget_Name = "He4";
	else
		cerr<<"   ***ERROR, in XGT2_Target.h, Line=14"<<endl;	
	
	TString Boiling_File=Form("%s/Target/%s_Boiling.out",XGT2_DIR.Data(), aTarget_Name.Data());
	const int tBin = 60; //Bin Number used in Boiling Study, all equal to 60
	const double step = 24.0/tBin; //+- 0.12 used on VZ cut
	double *tVZ = new double[tBin];
	double *tSlop = new double[tBin];

	ifstream infile(Boiling_File); 
	TString aCom; infile >> aCom >> aCom; //Reading the first line

	for(int i=0;i<tBin;i++){
		infile >> tVZ[i] >> tSlop[i];
		tVZ[i]*=100.0; //Convert from meter to cm
	}
	infile.close();

	int bin = (int) ((aVZ + 12.)/step);
	double slop = ( tSlop[bin+1] - (tSlop[bin+1]-tSlop[bin])/(tVZ[bin+1]-tVZ[bin]) * (tVZ[bin+1]-aVZ) );
	double rho = aRho / (1 + slop*aCurrent);

	//   cerr<<Form("Bin=%d, Slop=%f, VZ=%f/%f, I=%f, Rho0=%f, Rho=%f",bin,slop,aVZ,tVZ[bin], aCurrent, aRho, rho)<<endl;

	delete tVZ; delete tSlop;
	return rho;
}

/*}}}*/

/*inline void gGet_Cryo_Density(  const Int_t aRunNo, const Int_t aTarget_A, const Double_t aCurrent, const TString& aArm){{{*/
inline void gGet_Cryo_Density(  const Int_t aRunNo, const Int_t aTarget_A, const Double_t aCurrent, const TString& aArm, const Double_t aTarget_Length_Cut, double* aTarget_Density, double* aTarget_Density_Err){

	//Boiling Factor in each VZ bin
	const double aTarget_Length = 20.00;//cm
   	//Symmetric Cut, e.g., cut=7.5 cm, then we only count the target luminosity on these range
	//const double aLength_Corrected = 2.0*aTarget_Length_Cut/aTarget_Length;	
	const double aLength_Corrected = 1.0;
    const int aBin_Size =200;
    double aVZ[aBin_Size],aRho[aBin_Size],aRho_New[aBin_Size],aBF[aBin_Size];
    double aBF_Err[aBin_Size],aRho_Err[aBin_Size],aRho_New_Err[aBin_Size];

	TString aTarget_Density_Profile, aCom;
	ifstream aInfile;
	int aBin = 0;
	double aDensity_Old = 0.0, aDensity_New =0.0, aDensity_New_Err =0.0;
	//H2
	if(aTarget_A == 2){
		//Boiling_Factor = 0.188; aDensity = 3.352 * aLength_Corrected;
		aDensity_Old = 3.3811 * aLength_Corrected; //g/cm2
        aTarget_Density_Profile = "H2_new.rho";
	}	
	//He3
	else if(aTarget_A == 3){ 
		if(aRunNo<3753)
			aDensity_Old = 0.4206 * aLength_Corrected; //g/cm2
		else 
			//aDensity_Old = 0.55926 * aLength_Corrected; //g/cm2
			aDensity_Old = 0.5832 * aLength_Corrected; //g/cm2
		aTarget_Density_Profile = "He3_new.rho";
	}
	//He4
	else if(aTarget_A == 4){ 
		aDensity_Old = 0.64901 * aLength_Corrected;//g/cm2
        aTarget_Density_Profile = "He4_new.rho";
	}
	else
		cerr<<"****** ERROR, in XGT2_Target.h, gGet_Cryo_Density, Line#117, unknown Target_Name!!!"<<endl;

	aInfile.open(Form("%s/Target/%s",XGT2_DIR.Data(),aTarget_Density_Profile.Data()));
	aInfile >> aCom >> aCom >> aCom >> aCom >> aCom;
	for(int i=0;i<aBin_Size;i++){
		aInfile >> aVZ[i] >> aRho[i] >> aRho_Err[i] >> aBF[i] >> aBF_Err[i]; //aVZ[i] is in cm here
		if(abs(aVZ[i])<=aTarget_Length_Cut){ //Cut on Target Length, it is symmetric now.
			//aRho_New[i] = Boiling_Corr(aTarget_A, aVZ[i], aCurrent, aRho[i]); //Older File with BF only
			aRho_New[i] = aRho[i] * (1 + aBF[i]*aCurrent);

			aRho_New[i]*= aDensity_Old; //aRho is a normalized density distribution, with the sum equal to one
			aRho_New_Err[i] = aRho_New_Err[i] * sqrt( pow(aRho_Err[i]/aRho[i],2)+pow(aBF_Err[i]/aBF[i],2));
			aDensity_New+=aRho_New[i];
			aDensity_New_Err+=pow(aRho_New_Err[i]/aRho_New[i],2);
		}
	}
	aInfile.close();
    aDensity_New_Err = aDensity_New*sqrt(aDensity_New_Err);

	aTarget_Density[0]=aDensity_New;
	aTarget_Density_Err[0]=aDensity_New_Err;
	
//	cerr<<Form("      ---2: I = %5.2f, Target Thickness, old = %6.4f, Corrected = %6.4f", aCurrent, aDensity_Old, aDensity_New)<<endl;
}
/*}}}*/

/*inline double gGet_Density( const Int_t aRunNo, const Int_t aTarget_A, const double aCurrent, const TString& aArm){{{*/
inline double gGet_Density( const Int_t aRunNo, const Int_t aTarget_A, const Double_t aCurrent, const TString& aArm, const Double_t aTarget_Length_Cut){
	//Symmetric Cut, e.g., cut=7.5 cm, then we only count the target luminosity on these range
	//const double aLength_Corrected = aTarget_Length_Cut/10.0; 
	const double aLength_Corrected = 1.0;

  //Boiling Factor are calculated at 100uA
  //Values can be updated here.
  Double_t Boiling_Factor = 0.0, aDensity = 0.0;
  if(aTarget_A == 2){ //H2
	  //Boiling_Factor = 0.188; aDensity = 3.352;
	 aDensity = 3.352 * aLength_Corrected;
   	 if(aArm=="R"){
		  Boiling_Factor = 0.41179; //HRS-R
	  }
	  else {
		  Boiling_Factor = 0.46394;//HRS-L
	  }
  }	
  else if(aTarget_A == 3){ //He3
	  if(aArm=="R"){
		  Boiling_Factor = 0.25440; //HRS-R
	  }
	  else{ 
		  Boiling_Factor = 0.26573; //HRS-L
	  }
	  if(aRunNo<3753)
		  aDensity = 0.426 * aLength_Corrected;
	  else 
		  aDensity = 0.592 * aLength_Corrected;
  }
  else if(aTarget_A == 4){ //He4
	 aDensity = 0.648 * aLength_Corrected;
   	 if(aArm=="R"){
		  Boiling_Factor = 0.284774; //HRS-R
	  }
	  else{ 
		  Boiling_Factor = 0.30141; //HRS-L
	  }
  }
  else if(aTarget_A == 12){ //C12
	  Boiling_Factor = 0.0; aDensity = 0.8918;
  }
  else if(aTarget_A == 40){ //Ca40
	  Boiling_Factor = 0.0; aDensity = 0.8890;
  }
   else if(aTarget_A == 48){ //Ca48
    Boiling_Factor = 0.0; aDensity = 0.8190;
  }
  else
    Boiling_Factor = 0.0;
 
// Boiling_Factor = gGet_Boiling_Factor(aTarget_A, aArm.Data());
  double aDensity_Real = aDensity * (1.0 - Boiling_Factor * aCurrent /100.00); 
 // cerr<<Form("      ---1: I = %5.2f, Target Thickness, old = %6.4f, Corrected = %6.4f", aCurrent, aDensity, aDensity_Real)<<endl;
  return aDensity_Real;
}
/*}}}*/

/*inline XGT2_VAR* gCal_Ntg(const vector<Int_t>& aRunNoChain, const Double_t aTarget_Thickness, const Int_t aTarget_A){{{*/
inline XGT2_VAR* gCal_Ntg(const vector<Int_t>& aRunNoChain,const TString& aArm,const Int_t aTarget_A, const Double_t aTarget_Length_Cut, const Double_t aTarget_Thickness,const Double_t aTarget_Thickness_Stat_Err, double* aNtgChain){
  XGT2_VAR* aNtg = new XGT2_VAR();
  Double_t aNtg_Sum = 0.0, aNtg_Err=0.0;
  for (unsigned int i=0; i<aRunNoChain.size(); i++ ){
      Int_t aRunNo=aRunNoChain[i];
      TChain *aT_Tree = gAddTree(aRunNo,"T"); 
	  Double_t aCurrent = gGet_Current(aT_Tree, aArm.Data());
	  Double_t aTarget_Thickness_Corrected = 0.0;
	  Double_t aTarget_Thickness_Corrected_Err = 0.0;
	  //Cryogenic Targets --> LH2, He3,He4
	  if(aTarget_A == 2 || aTarget_A == 3 || aTarget_A ==4){
		//    aTarget_Thickness_Corrected = gGet_Density( aRunNo, aTarget_A, aCurrent, aArm.Data(), aTarget_Length_Cut); //Average Boiling Effect
		 //gGet_Cryo_Density( aRunNo, aTarget_A, aCurrent, aArm.Data(), aTarget_Length_Cut, &aTarget_Thickness_Corrected, &aTarget_Thickness_Corrected_Err); //Z-Dependence Boiling Effect
		  gGet_Cryo_Density( aRunNo, aTarget_A, aCurrent, aArm.Data(), 10.0, &aTarget_Thickness_Corrected, &aTarget_Thickness_Corrected_Err); //Z-Dependence Boiling Effect
		  Double_t aNtg_Corrected = aTarget_Thickness_Corrected * Na / aTarget_A ; 
		  aNtgChain[i] = aNtg_Corrected;
	  }
	  else{
		  aNtgChain[i] = aTarget_Thickness * Na / aTarget_A ;
	  }

	  outlog<<Form("      For Run#=%d, Current = %6.2f, rho= %8.6e, Ntg = %4.3e",aRunNo, aCurrent, aTarget_Thickness_Corrected, aNtgChain[i])<<endl;
	  cerr  <<Form("      For Run#=%d, Current = %6.2f, rho= %8.6e, Ntg = %4.3e",aRunNo, aCurrent, aTarget_Thickness_Corrected, aNtgChain[i])<<endl;
      aNtg_Sum += aNtgChain[i];
      aNtg_Err += pow((aTarget_Thickness_Stat_Err/aTarget_Thickness),2);
      delete aT_Tree;
  }
  //Take the average value of all runs
  if(aRunNoChain.size()>1e-16){
	  aNtg->Value    = aNtg_Sum/aRunNoChain.size(); 
	  aNtg->Stat_Err = 0.0; 
	  aNtg->Sys_Err  = aNtg->Value * sqrt(aNtg_Err/aRunNoChain.size()); 
  }
  else
	  aNtg->SetValue(-1.0,0.0,0.0);

  return aNtg;
}
/*}}}*/
