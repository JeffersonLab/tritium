/*inline void gGet_Cryo_Density(  const Int_t aRunNo, const Int_t aTarget_A, const Double_t aCurrent){{{*/
inline void gGet_Cryo_Density(  const Int_t aRunNo, const Int_t aTarget_A, const Double_t aCurrent, double* aTarget_Density, double* aTarget_Density_Err){

	//Boiling Factor in each VZ bin
	const double aTarget_Length = 20.00;//cm
	const double aBin_Width = 1.0;
    const int aBin_Size =200;
    double aVZ[aBin_Size],aRho[aBin_Size],aRho_Boil[aBin_Size],aBF[aBin_Size];
    double aBF_Err[aBin_Size],aRho_Err[aBin_Size],aRho_Boil_Err[aBin_Size];

	TString aTarget_Density_Profile, aCom, aOutput_FileName;
	ifstream aInfile;
	int aBin = 0;
	double aDensity_Zero = 0.0, aDensity_Boil =0.0, aDensity_Boil_Err =0.0;
	//H2
	if(aTarget_A == 2){
		aDensity_Zero = 3.3811 * aBin_Width; //g/cm2
        aTarget_Density_Profile = "H2_new.rho";
        aOutput_FileName = Form("H2_I%d.rho", (int)(aCurrent));
	}	
	//He3
	else if(aTarget_A == 3){ 
		if(aRunNo<3753)
			aDensity_Zero = 0.4206 * aBin_Width; //g/cm2
		else 
			//aDensity_Zero = 0.55926 * aBin_Width; //g/cm2
			aDensity_Zero = 0.5832 * aBin_Width; //g/cm2
		aTarget_Density_Profile = "He3_new.rho";
        aOutput_FileName = Form("He3_I%d.rho", (int)(aCurrent));
	}
	//He4
	else if(aTarget_A == 4){ 
		aDensity_Zero = 0.64901 * aBin_Width;//g/cm2
        aTarget_Density_Profile = "He4_new.rho";
        aOutput_FileName = Form("He4_I%d.rho", (int)(aCurrent));
	}
	else
		cerr<<"****** ERROR, in XGT2_Target.h, gGet_Cryo_Density, Line#117, unknown Target_Name!!!"<<endl;

    ofstream aOutfile(aOutput_FileName.Data());

	aInfile.open(Form("%s/Target/%s",XGT2_DIR.Data(),aTarget_Density_Profile.Data()));
	aInfile >> aCom >> aCom >> aCom >> aCom >> aCom;
	for(int i=0;i<aBin_Size;i++){
		//     Vertex-Z-Bin  Density_Norm   Err      Boiling_Factor
		aInfile >> aVZ[i] >> aRho[i] >> aRho_Err[i] >> aBF[i] >> aBF_Err[i]; //aVZ[i] is in cm here
		aRho_Boil[i] = aRho[i] * (1 + aBF[i]*aCurrent);//Adding the boiling effect

		aRho_Boil[i]*= aDensity_Zero; //aRho is a normalized density distribution, with the sum equal to one
		aRho_Boil_Err[i] = aRho_Boil_Err[i] * sqrt( pow(aRho_Err[i]/aRho[i],2)+pow(aBF_Err[i]/aBF[i],2));
		aDensity_Boil+=aRho_Boil[i];
		aDensity_Boil_Err+=pow(aRho_Boil_Err[i]/aRho_Boil[i],2);
		aOutfile << aVZ[i] << aRho_Boil[i] << aRho_Boil_Err[i]<<endl; //aVZ[i] is in cm here
	}
	aInfile.close();
	aOutfile.close();

    aDensity_Boil_Err = aDensity_Boil*sqrt(aDensity_Boil_Err);

	aTarget_Density[0]=aDensity_Boil;
	aTarget_Density_Err[0]=aDensity_Boil_Err;
	
//	cerr<<Form("      ---2: I = %5.2f, Target Thickness, old = %6.4f, Corrected = %6.4f", aCurrent, aDensity_Zero, aDensity_Boil)<<endl;
}
/*}}}*/

int main(){

 TString Target = "He3";
 int Target_A = 3;
 int RunNo = 4299;
 double Current = 120; //uA
 double Density = 0.0, Density_Err = 0.; //g/cm2

 gGet_Cryo_Density(RunNo,Target_A,Current, &Density, &Density_Err);

 return 0;
}
