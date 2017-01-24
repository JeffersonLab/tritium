/*Correction{{{*/
double Correction(const TString& aTarget, const double aVZ, const double aCurrent, const double aRho){
	//TString Boiling_File=Form("%s_Boiling.out",aTarget.Data());
	TString Boiling_File=Form("%s_Boiling.out",aTarget.Data());
	const int tBin = 60; //Bin Number used in Boiling Study, all equal to 60
	const double step = 0.24/tBin; //+- 0.12 used on VZ cut
	double *tVZ = new double[tBin];
	double *tSlop = new double[tBin];

	ifstream infile(Boiling_File); 
	TString aCom; infile >> aCom >> aCom; //Reading the first line

	for(int i=0;i<tBin;i++){
		infile >> tVZ[i] >> tSlop[i];
	}
	infile.close();

	int bin = (int) ((aVZ + 0.12)/step);
	double slop = ( tSlop[bin+1] - (tSlop[bin+1]-tSlop[bin])/(tVZ[bin+1]-tVZ[bin]) * (tVZ[bin+1]-aVZ) );
	double rho = aRho * (1 + slop*aCurrent);

	cerr<<Form("Bin=%d, Slop=%f, VZ=%f/%f, I=%f, Rho0=%f, Rho=%f",bin,slop,aVZ,tVZ[bin], aCurrent, aRho, rho)<<endl;

	delete tVZ; delete tSlop;
	return rho;
}

/*}}}*/

/*Correction-Inverse{{{*/
double Correction_Inverse(const TString& aTarget, const double aVZ, const double aCurrent, const double aRho){
	TString Boiling_File=Form("%s_Boiling.out",aTarget.Data());
	const int tBin = 60; //Bin Number used in Boiling Study, all equal to 60
	const double step = 0.24/tBin; //+- 0.12 used on VZ cut
	double *tVZ = new double[tBin];
	double *tSlop = new double[tBin];

	ifstream infile(Boiling_File); 
	TString aCom; infile >> aCom >> aCom; //Reading the first line

	for(int i=0;i<tBin;i++){
		infile >> tVZ[i] >> tSlop[i];
	}
	infile.close();

	int bin = (int) ((aVZ + 0.12)/step);
	double slop = ( tSlop[bin+1] - (tSlop[bin+1]-tSlop[bin])/(tVZ[bin+1]-tVZ[bin]) * (tVZ[bin+1]-aVZ) );
	double rho = aRho / (1 + slop*aCurrent);

	//   cerr<<Form("Bin=%d, Slop=%f, VZ=%f/%f, I=%f, Rho0=%f, Rho=%f",bin,slop,aVZ,tVZ[bin], aCurrent, aRho, rho)<<endl;

	delete tVZ; delete tSlop;
	return rho;
}

/*}}}*/

void easy(){
 TString Target ; cerr<<" --- Target: "; cin >> Target;
 TString Input_File = Form("%s.rho",Target.Data()) ;
 TString Output_File = Input_File;
 Output_File.ReplaceAll(".rho","_new.rho");
 double rho0=0.0, I0=0;

 if(Target=="H2"){
	 rho0 = 2.65219; I0 = 45;
 }
 else if(Target=="He3"){
	 rho0 = 0.29017;I0 = 120;
 }
 else if(Target=="He4"){
	 rho0 = 0.46245;I0 = 95;
 }
 else
	 cerr<<"*** ERROR ***"<<endl;

 ofstream outf(Output_File);

 double tVZ[52], tRho[52], tRho_Err[52], tBF[52], tBF_Err[52];
 double aVZ[200], aRho[200], aRho_Err[200], aBF[200], aBF_Err[200];

 ifstream inf;inf.open(Input_File);
 TString aaa;
 inf >> aaa >> aaa >> aaa >> aaa >> aaa;
 for(int i=0;i<52;i++){
   inf >> tVZ[i] >> tRho[i]>> tRho_Err[i]>> tBF[i]>> tBF_Err[i];
   cerr<<"VZ = "<<tVZ[i]<<"  rho="<<tRho[i]<<endl;
 }
 inf.close();

 double sum=0;
 double step = 0.20/200;
 int bin = 0;
 for(int i=0;i<200;i++){
   aVZ[i] = -0.10 + step*(i+0.5); 
 
   bin = (int)((aVZ[i] + 0.102)/0.004);
   aRho[i] = ( tRho[bin+1] - (tRho[bin+1]-tRho[bin])/(tVZ[bin+1]-tVZ[bin]) * (tVZ[bin+1]-aVZ[i]) );
   aRho_Err[i] = 0.5* aRho[i] * sqrt( pow(tRho_Err[bin+1]/tRho[bin+1],2) + pow(tRho_Err[bin]/tRho[bin],2) );
   aBF[i] = ( tBF[bin+1] - (tBF[bin+1]-tBF[bin])/(tVZ[bin+1]-tVZ[bin]) * (tVZ[bin+1]-aVZ[i]) );
   aBF_Err[i] = 0.5* aBF[i] * sqrt( pow(tBF_Err[bin+1]/tBF[bin+1],2) + pow(tBF_Err[bin]/tBF[bin],2) );
   sum+=aRho[i];
 }

    double rho_sum = 0;
	outf<<Form("%12s %12s %12s %12s %12s", "VZ","Y0","Y0_Err","BF","BF_Err")<<endl; 
	for(int i=0;i<200;i++){
		aRho[i]/=sum;
		aRho_Err[i]/=sum;
//		aRho[i]*= rho0/20.0;
//		aRho_Err[i]*= rho0/20.0;

		rho_sum +=aRho[i];
		aVZ[i]*=100;
		outf<<Form("  %12.8f %12.8f %12.8f %12.8f %12.8f",aVZ[i],aRho[i], aRho_Err[i], aBF[i], aBF_Err[i])<<endl; 
	}
	outf.close();

	cerr<<" --- Rho_0="<<rho0<< ", Rho_Sum = "<< rho_sum<<endl;
}
