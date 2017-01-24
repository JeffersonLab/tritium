/* 
Update: Add Boiling Effects along Z_react, in 01/09/2013
*/

/*inline double gGet_DensityDistribution(const double aVZ,const int aTarget_A){{{*/
inline double gGet_DensityDistribution(const double aVZ, const int aTarget_A){

	//Density Distribution: Polynomials Functions of VZ in meter
	const double aPol_H2[7] = {0.9508, -4.9160, 74.34, 1603.00, -2.900e04, -2.191e05, 3.513e06}; 
    const double aPol_He3[7]= {0.8690, -2.8580, 214.5,-2212.00, -9.511e04,  7.468e05, 1.617e07};
    const double aPol_He4[7]= {0.8826, -7.4000, 238.1, 1105.00, -1.158e05,  2.604e04, 1.698e07};
	double kVZ = aVZ/100.00; //cm -> m
    //This is just a distribution but not an absolute value of density
	double aDensity = 0.0;
	if(aTarget_A == 2)
		aDensity =  aPol_H2[0] +  aPol_H2[1]*pow(kVZ,1) +  aPol_H2[2]*pow(kVZ,2) +  aPol_H2[3]*pow(kVZ,3) +  aPol_H2[4]*pow(kVZ,4) +  aPol_H2[5]*pow(kVZ,5) +  aPol_H2[6]*pow(kVZ,6);
	else if(aTarget_A == 3)
		aDensity = aPol_He3[0] + aPol_He3[1]*pow(kVZ,1) + aPol_He3[2]*pow(kVZ,2) + aPol_He3[3]*pow(kVZ,3) + aPol_He3[4]*pow(kVZ,4) + aPol_He3[5]*pow(kVZ,5) + aPol_He3[6]*pow(kVZ,6);
	else if(aTarget_A == 4)
		aDensity = aPol_He4[0] + aPol_He4[1]*pow(kVZ,1) + aPol_He4[2]*pow(kVZ,2) + aPol_He4[3]*pow(kVZ,3) + aPol_He4[4]*pow(kVZ,4) + aPol_He4[5]*pow(kVZ,5) + aPol_He4[6]*pow(kVZ,6);
    else
		aDensity = 1.0;

	return aDensity;
}
/*}}}*/

/*inline double gGet_DensityByBinZ( const double aDensity, const Int_t aRunNo, const Int_t aTarget_A, const Double_t aCurrent, const TString& aArm){{{*/
inline double gGet_DensityByBinZ( const double aTarget_Length_Cut, const Int_t aRunNo, const Int_t aTarget_A, const Double_t aCurrent, const TString& aArm){

	//Boiling Factor in each VZ bin
	const double aTarget_Length = 20.00;//cm
	const double aBin_Size = 1.00;//cm
	//(cm):                       -9.5  -8.5 -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5
	const double aBF_H2_L[20]  = { 9.1,  9.1, 9.1, 12.4, 13.9, 12.9, 12.8, 16.6, 22.8, 26.2, 26.4, 25.4, 27.1, 26.5, 28.4, 28.0, 27.9, 28.1, 28.1, 28.1};
	const double aBF_H2_R[20]  = { 7.7,  7.7, 7.7, 11.1, 12.6, 11.6, 11.5, 15.3, 21.6, 25.2, 25.3, 24.3, 26.0, 25.5, 27.3, 26.9, 26.9, 27.1, 27.1, 27.1};
	//(cm):                       -9.5  -8.5 -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5
	const double aBF_He3_L[20] = {20.9, 20.9, 20.9, 22.2, 22.1, 23.4, 30.3, 37.2, 35.9, 36.1, 35.6, 36.6, 37.1, 37.5, 38.6, 37.6, 38.2, 36.7, 36.7, 36.7};
	const double aBF_He3_R[20] = {20.7, 20.7, 20.7, 22.0, 22.0, 23.3, 30.2, 37.1, 35.8, 35.9, 35.5, 36.5, 37.0, 37.4, 38.5, 37.5, 38.1, 36.5, 36.5, 36.5};
	//(cm):                       -9.5  -8.5 -7.5  -6.5  -5.5  -4.5  -3.5  -2.5  -1.5  -0.5   0.5   1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5
	const double aBF_He4_L[20] = {17.5, 17.5, 17.5, 17.3, 17.4, 17.5, 22.0, 32.0, 36.7, 36.5, 34.6, 35.3, 35.5, 35.0, 36.0, 35.6, 36.0, 34.0, 34.0, 34.0};
	const double aBF_He4_R[20] = {17.2, 17.2, 17.2, 17.0, 17.1, 17.2, 21.7, 31.8, 36.5, 36.2, 34.4, 35.1, 35.3, 34.8, 35.7, 35.4, 35.8, 33.8, 33.8, 33.8};

	double aBin = 0.0, aDensity = 0.0, aDensity_Old = 0.0, aDensity_Sum =0.0;
	//H2
	if(aTarget_A == 2){
		//Boiling_Factor = 0.188; aDensity = 3.352;
		aDensity_Old = 3.352; //g/cm2
		aDensity = aDensity_Old / aTarget_Length; //g/cm3
		if(aArm=="L"){
			aDensity_Sum =0.0;
			for(int i=0;i<20;i++){
				aBin = -9.5 + i*aBin_Size;//cm
				aDensity *= gGet_DensityDistribution(aBin, aTarget_A); //g/cm3
				aDensity *= (1.0 - aBF_H2_L[i] * aCurrent /100.00); //g/cm3
				aDensity_Sum += aDensity * aBin_Size; //g/cm2	
			}
		}
		else if(aArm=="R"){
			aDensity_Sum =0.0;
			for(int i=0;i<20;i++){
				aBin = -9.5 + i*aBin_Size;//cm
				aDensity *= gGet_DensityDistribution(aBin, aTarget_A); //g/cm3
				aDensity *= (1.0 - aBF_H2_R[i] * aCurrent /100.00); //g/cm3
				aDensity_Sum += aDensity * aBin_Size; //g/cm2	
			}
		}
		else
			cerr<<"****** ERROR, in XGT2_Target.h, gGet_DensitybyBinZ, Line#68, unknown HRS Arm!!!"<<endl;
	}	
	//He3
	else if(aTarget_A == 3){ 
		if(aRunNo<3753)
			aDensity_Old = 0.426; //g/cm2
		else 
			aDensity_Old = 0.592; //g/cm2
		aDensity = aDensity_Old / aTarget_Length; //g/cm3
		
		if(aArm=="L"){
			aDensity_Sum =0.0;
			for(int i=0;i<20;i++){
				aBin = -9.5 + i*aBin_Size;//cm
				aDensity *= gGet_DensityDistribution(aBin, aTarget_A); //g/cm3
				aDensity *= (1.0 - aBF_He3_L[i] * aCurrent /100.00); //g/cm3
				aDensity_Sum += aDensity * aBin_Size; //g/cm2	
			}
		}
		else if(aArm=="R"){
			aDensity_Sum =0.0;
			for(int i=0;i<20;i++){
				aBin = -9.5 + i*aBin_Size;//cm
				aDensity *= gGet_DensityDistribution(aBin, aTarget_A); //g/cm3
				aDensity *= (1.0 - aBF_He3_R[i] * aCurrent /100.00); //g/cm3
				aDensity_Sum += aDensity * aBin_Size; //g/cm2	
			}
		}
		else
			cerr<<"****** ERROR, in XGT2_Target.h, gGet_DensitybyBinZ, Line#96, unknown HRS Arm!!!"<<endl;
	}
	//He4
	else if(aTarget_A == 4){ 
		aDensity_Old = 0.648;//g/cm2
		aDensity = aDensity_Old / aTarget_Length; //g/cm3
		if(aArm=="L"){
			aDensity_Sum =0.0;
			for(int i=0;i<20;i++){
				aBin = -9.5 + i*aBin_Size;//cm
				aDensity *= gGet_DensityDistribution(aBin, aTarget_A); //g/cm3
				aDensity *= (1.0 - aBF_He4_L[i] * aCurrent /100.00); //g/cm3
				aDensity_Sum += aDensity * aBin_Size; //g/cm2	
			}
		}
		else if(aArm=="R"){
			aDensity_Sum =0.0;
			for(int i=0;i<20;i++){
				aBin = -9.5 + i*aBin_Size;//cm
				aDensity *= gGet_DensityDistribution(aBin, aTarget_A); //g/cm3
				aDensity *= (1.0 - aBF_He4_R[i] * aCurrent /100.00); //g/cm3
				aDensity_Sum += aDensity * aBin_Size; //g/cm2	
			}
		}
		else
			cerr<<"****** ERROR, in XGT2_Target.h, gGet_DensitybyBinZ, Line#120, unknown HRS Arm!!!"<<endl;
	}
    //Add other targets in case I call this by mistake
	else if(aTarget_A == 12){ //C12
		aDensity = 0.8918;//g/cm2
		aDensity_Sum = aDensity; 
	}
	else if(aTarget_A == 40){ //Ca40
		aDensity = 0.8890;//g/cm2
		aDensity_Sum = aDensity; 
	}
	else if(aTarget_A == 48){ //Ca48
		aDensity = 0.8190;//g/cm2
		aDensity_Sum = aDensity; 
	}
	else
		cerr<<"****** ERROR, in XGT2_Target.h, gGet_DensitybyBinZ, Line#136, unknown Target_Name!!!"<<endl;

	cerr<<Form("      --- Target Thickness, old = %6.4f, Corrected = %6.4f", aDensity_Old, aDensity_Sum)<<endl;
	return aDensity_Sum;
}

/*}}}*/

/*inline double gGet_Density( const double aTarget_Length_Cut, const Int_t aRunNo, const Int_t aTarget_A, const double aCurrent, const TString& aArm){{{*/
inline double gGet_Density( const double aTarget_Length_Cut, const Int_t aRunNo, const Int_t aTarget_A, const Double_t aCurrent, const TString& aArm){
  
//const double aLength_Corrected = 2.0*aTarget_Length_Cut / 0.20; //0.20meter is the total length of target cell 
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
  return aDensity_Real;
}
/*}}}*/

/*inline XGT2_VAR* gCal_Ntg(const vector<Int_t>& aRunNoChain, const Double_t aTarget_Thickness, const Int_t aTarget_A){{{*/
inline XGT2_VAR* gCal_Ntg(const vector<Int_t>& aRunNoChain,const TString& aArm, const TString& aTarget, const Double_t aTarget_Length_Cut, double* aNtgChain){
  int aTarget_A = 0;
	if(aTarget=="H2")
      aTarget_A = 2;
	else if(aTarget=="He3")
      aTarget_A = 3;
	else if(aTarget=="He4")
      aTarget_A = 4;
	else if(aTarget=="C12")
      aTarget_A =12;
	else if(aTarget=="Ca40")
		aTarget_A = 40;
	else if(aTarget=="Ca48")
		aTarget_A =48;
    else
		aTarget_A = -1;

	XGT2_VAR* aNtg = new XGT2_VAR();
  Double_t aNtg_Sum = 0.0, aNtg_Err=0.0;
  for (unsigned int i=0; i<aRunNoChain.size(); i++ ){
      Int_t aRunNo=aRunNoChain[i];
      TChain *aT_Tree = gAddTree(aRunNo,"T"); 
      Double_t aCurrent = gGet_Current(aT_Tree, aArm.Data());

      //Cryogenic Targets --> LH2, He3,He4
	  Double_t aTarget_Thickness_Corrected = gGet_Density(aTarget_Length_Cut, aRunNo, aTarget_A, aCurrent, aArm.Data()); //Average Boiling Effect
	  Double_t aNtg_Corrected = aTarget_Thickness_Corrected * Na / aTarget_A ; 
	  aNtgChain[i] = aNtg_Corrected;

	  outlog<<Form("      For Run#=%d, Current = %6.2f, Ntg = %4.3e",aRunNo, aCurrent, aNtgChain[i])<<endl;
	  cerr<<Form("      For Run#=%d, Current = %6.2f, Ntg = %4.3e",aRunNo, aCurrent, aNtgChain[i])<<endl;
      aNtg_Sum += aNtgChain[i];
      //aNtg_Err += pow((aTarget_Thickness_Stat_Err/aTarget_Thickness),2);
      aNtg_Err += 0.0;
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
