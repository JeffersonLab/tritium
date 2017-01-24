#define NOCASE 1

/*Double_t Get_Mott_CrossSection(const Double_t& aAngle,const Double_t& aEi){{{*/
/*Double_t Get_Mott_CrossSection(const Double_t& aAngle,const Double_t& aEi) 
{
  return ALPHA*ALPHA*pow(cos(aAngle/2*TMath::DegToRad()),2)/4/aEi/aEi/pow(sin(aAngle/2*TMath::DegToRad()),4);// unit=1/(MeV^2*sr)
}*/
/*}}}*/

/*inline Double_t gCal_Eloss(const Double_t& aE,const Int_t& aZ,const Double_t& aA,const Double_t& aM,const Double_t& aI,const Double_t& aDensity,const Double_t& aThickness){{{*/
inline Double_t gCal_Eloss(const Double_t& aE,const Int_t& aZ,const Double_t& aA,const Double_t& aM,const Double_t& aI,const Double_t& aDensity,const Double_t& aThickness)
{
  //aE,aM,aI GeV, aA kg/mol, aDensity g/cm^3
  Double_t lE=aE*1000; //GeV to MeV
  Double_t lM=aM*1000; //GeV to MeV
  Double_t lI=aI*1000; //GeV to MeV
  Double_t lA=aA*1000; //kg/mol to g/mol
  Double_t lThickness=aThickness/10; //kg/m^2 to g/cm^2

  Double_t K=0.307075; //MeV cm^2/g
  Double_t z=1;//electron
  Double_t me=0.511; //MeV
  Double_t beta,gamma,omega_p,Tmax,delta;

  beta=1-me*me/(lE*lE);
  gamma=1/sqrt(1-beta*beta);

  omega_p=28.816*sqrt(aDensity*aZ/lA/1000)*1e-6; //MeV
  Tmax=2*me*beta*beta*gamma*gamma/(1+2*gamma*me/lM+me*me/(lM*lM));
  delta=log(omega_p/lI)+log(beta*gamma)-0.5;
  return (K*z*z*aZ/lA/(beta*beta)*( 0.5*log(2*me*beta*beta*gamma*gamma*Tmax/lI/lI)-beta*beta-delta ))*lThickness*1e-3;//GeV
}
/*}}}*/

/*inline XGT2_VAR* gCal_P0(TTree* aTree,const TString& aArm,const Double_t& aSpec_Mom,Double_t& oDipole){{{*/
//Fix Here,ZYE
inline XGT2_VAR* gCal_P0(TTree* aTree,const TString& aArm,Double_t& oDipole)
{
  //return GeV, oDipole is kG.
  Int_t nentries=(Int_t)aTree->GetEntries();
  Double_t dipole;
  Double_t P0=0;

  aTree->SetBranchAddress(Form("Hac%s_D1_NMR_SIG",aArm.Data()), &dipole); // NMR field
   
  Double_t gamma[4];
  if ( aArm=="L" )
    {
      gamma[0]=0;
      gamma[1]=270.2;
      gamma[2]=0;
      gamma[3]=-1.6e-03;
    }
  else
    {
      gamma[0]=0;
      gamma[1]=269.8;
      gamma[2]=0;
      gamma[3]=-1.6e-03;
    }
    
  int i,j;
  Double_t thisp0;
  Int_t NoOfValidEntries=0;
  oDipole=0;
  for(i=0; i<nentries; i++)
    {
      aTree->GetEntry(i);
      if(fabs(dipole)<10)
	{
	  thisp0=0;
	  for ( j=0; j<4; j++ )
	    thisp0 += gamma[j]*pow(dipole*10,j);
	  P0 += thisp0;
	  oDipole+=dipole*10;
	  NoOfValidEntries++;
	}
    }
  if ( NoOfValidEntries!=0 )
    {
      P0 /= NoOfValidEntries;
      oDipole /= NoOfValidEntries;
    }
  else
    {
      P0=0;
      oDipole=0;
    }

  XGT2_VAR* output=new XGT2_VAR();
  output->SetValue(P0/1000,P0/1000*1e-4,0);
 
  return output;//MeV to GeV
}
/*}}}*/

/*inline void gLoad_EpicsTable(const int* aRunNo, const TString& aArm, double* aE0, double* aDEE, double* aMom, double* aAngle){{{*/
inline void gLoad_EpicsTable(const int aRunNo, const TString& aArm, double* aE0, double* aDEE, double* aMom, double* aAngle){
  ifstream aInfile(Form("%s/Table/EPICS_Run_List.table",XGT2_DIR.Data()));
  TString aCom; aCom.ReadLine(aInfile);

  double aTgt_P, aI, aAngle_Left, aAngle_Right, aMom_Left, aMom_Right, aE0_List, aDEE_List;
  int aRunNo_List = 0;
  TString aTgt_Name;
  //bool bAngle_Correct = kTRUE;
  bool bAngle_Correct = kFALSE;

  while(!aInfile.eof()){
	  aInfile >> aRunNo_List >> aTgt_Name >> aTgt_P >> aI >> aE0_List >> aDEE_List >> aMom_Left >> aAngle_Left >> aMom_Right >> aAngle_Right;

	  if(aRunNo == aRunNo_List){
	     *aE0= aE0_List;
		 *aDEE = aDEE_List;

		 if(aArm=="L"){
			 *aAngle = aAngle_Left;
			 *aMom = aMom_Left;
             
			 /*Offsets using JinGe's Optics Matrix{{{*/
			 
			 if((aRunNo >= 3565 && aRunNo <= 3656) && fabs(aAngle_Left - 25.00)<0.5)
				 //*aAngle = aAngle_Left;
				 *aAngle = aAngle_Left + (24.98 - 25.00);//Using Dien's new pointing study
			 else if((aRunNo >= 3657 && aRunNo <= 3683) && fabs(aAngle_Left - 21.00)<0.5)
				 //*aAngle = aAngle_Left + (21.03 - 21.00);
				 *aAngle = aAngle_Left + (20.98 - 21.00);//Using Dien's new pointing study
			 else if((aRunNo >= 3684 && aRunNo <= 3708) && fabs(aAngle_Left - 23.00)<0.5)
				 //*aAngle = aAngle_Left + (23.00 - 23.00);
				 *aAngle = aAngle_Left + (22.98 - 23.00);//Using Dien's new pointing study
			 else if((aRunNo >= 3735 && aRunNo <= 3891) && fabs(aAngle_Left - 25.00)<0.5)
				 //*aAngle = aAngle_Left + (24.99 - 25.00);
				 *aAngle = aAngle_Left + (24.98 - 25.00);//Using Dien's new pointing study
			 else if((aRunNo >= 3892 && aRunNo <= 3916) && fabs(aAngle_Left - 21.00)<0.5)
				 //*aAngle = aAngle_Left + (21.03 - 21.00);
				 *aAngle = aAngle_Left + (24.98 - 25.00);//Using Dien's new pointing study
			 else if((aRunNo >= 3917 && aRunNo <= 4071) && fabs(aAngle_Left - 28.00)<0.5)
				 //*aAngle = aAngle_Left + (27.98 - 28.00); //Or 28.03 from Multi-C, 28.00 from C12
				 *aAngle = aAngle_Left + (27.98 - 28.00);//Using Dien's new pointing study
			 else if((aRunNo >= 4073 && aRunNo <= 4103) && fabs(aAngle_Left - 21.00)<0.5)
				 //*aAngle = aAngle_Left + (21.04 - 21.00);                                        
				 *aAngle = aAngle_Left + (20.99 - 21.00);//Using Dien's new pointing study
			 else if((aRunNo >= 4112 && aRunNo <= 4179) && fabs(aAngle_Left - 23.00)<0.5)
				 //*aAngle = aAngle_Left + (23.00 - 23.00);                                        
				 *aAngle = aAngle_Left + (22.98 - 23.00);//Using Dien's new pointing study
			 else if((aRunNo >= 4181 && aRunNo <= 4241) && fabs(aAngle_Left - 25.00)<0.5)
				 //*aAngle = aAngle_Left + (24.98 - 25.00);                                        
				 *aAngle = aAngle_Left + (24.97 - 25.00);//Using Dien's new pointing study
			 else if((aRunNo >= 4242 && aRunNo <= 4250) && fabs(aAngle_Left - 21.00)<0.5)
				 //*aAngle = aAngle_Left + (21.02 - 21.00);                                        
				 *aAngle = aAngle_Left + (20.99 - 21.00);//Using Dien's new pointing study
			 else if((aRunNo >= 4251 && aRunNo <= 4299) && fabs(aAngle_Left - 28.00)<0.5)
				 //*aAngle = aAngle_Left + (27.98 - 28.00);                                        
				 *aAngle = aAngle_Left + (27.98 - 28.00);//Using Dien's new pointing study
			 else{
				 *aAngle = aAngle_Left;
				 cerr <<Form("      *** ERROR - In XGT2_Physics.h, RunNo#%d is not in the table!!!",aRunNo)<<endl;
				 outlog <<Form("      *** ERROR - In XGT2_Physics.h, RunNo#%d is not in the table!!!",aRunNo)<<endl;
			 }
			 
			 /*}}}*/
             
			 /*New Optics Matrix{{{*/
			 /*
			 if((aRunNo >= 3565 && aRunNo <= 3656) && fabs(aAngle_Left - 25.00)<0.5)
				 *aAngle = aAngle_Left;
			 else if((aRunNo >= 3657 && aRunNo <= 3683) && fabs(aAngle_Left - 21.00)<0.5)
				 *aAngle = aAngle_Left + (21.01 - 21.00);
			 else if((aRunNo >= 3684 && aRunNo <= 3708) && fabs(aAngle_Left - 23.00)<0.5)
				 *aAngle = aAngle_Left + (23.00 - 23.00);
			 else if((aRunNo >= 3735 && aRunNo <= 3891) && fabs(aAngle_Left - 25.00)<0.5)
				 *aAngle = aAngle_Left + (25.00 - 25.00);
			 //21 data are R-Only
			 else if((aRunNo >= 3917 && aRunNo <= 4071) && fabs(aAngle_Left - 28.00)<0.5)
				 *aAngle = aAngle_Left + (28.03 - 28.00); //Or 28.03 from Multi-C, 28.00 from C12
			 else if((aRunNo >= 4073 && aRunNo <= 4103) && fabs(aAngle_Left - 21.00)<0.5)
				 *aAngle = aAngle_Left + (21.01 - 21.00);                                        
			 else if((aRunNo >= 4112 && aRunNo <= 4179) && fabs(aAngle_Left - 23.00)<0.5)
				 *aAngle = aAngle_Left + (23.00 - 23.00);                                        
			 else if((aRunNo >= 4181 && aRunNo <= 4241) && fabs(aAngle_Left - 25.00)<0.5)
				 *aAngle = aAngle_Left + (24.98 - 25.00);                                        
			 else if((aRunNo >= 4242 && aRunNo <= 4250) && fabs(aAngle_Left - 21.00)<0.5)
				 *aAngle = aAngle_Left + (21.02 - 21.00);                                        
			 else if((aRunNo >= 4251 && aRunNo <= 4299) && fabs(aAngle_Left - 28.00)<0.5)
				 *aAngle = aAngle_Left + (28.00 - 28.00);                                        
			 else{
				 *aAngle = aAngle_Left;
				 cerr <<Form("      *** ERROR - In XGT2_Physics.h, RunNo#%d is not in the table!!!",aRunNo)<<endl;
				 outlog <<Form("      *** ERROR - In XGT2_Physics.h, RunNo#%d is not in the table!!!",aRunNo)<<endl;
			 }
			 */
			 /*}}}*/

			 /*Old Correction{{{*/
			 /*
				if(bAngle_Correct){
				if(fabs(aAngle_Left-21.00)<0.5)
			  *aAngle = aAngle_Left + (21.214002 - 21.00); //From Pointing Study

			  if(fabs(aAngle_Left-23.00)<0.5)
			  *aAngle = aAngle_Left + (23.212207 - 23.00); //From Pointing Study

			  if(fabs(aAngle_Left-25.00)<0.5)
			  *aAngle = aAngle_Left + (25.210281 - 25.00); //From Pointing Study

			  if(fabs(aAngle_Left-28.00)<0.5)
			  *aAngle = aAngle_Right + (28.207145-28.00); //From Pointing Study
			  }
			  */
			 /*}}}*/
		 }

		 if(aArm=="R"){
			 *aAngle = aAngle_Right;
			 *aMom = aMom_Right;

			 if((aRunNo >= 3565 && aRunNo <= 3656) && fabs(aAngle_Right - 25.00)<0.5)
				 *aAngle = aAngle_Right;
			 else if((aRunNo >= 3657 && aRunNo <= 3683) && fabs(aAngle_Right - 21.00)<0.5)
				 *aAngle = aAngle_Right + (21.04 - 21.00);
			 else if((aRunNo >= 3684 && aRunNo <= 3708) && fabs(aAngle_Right - 23.00)<0.5)
				 *aAngle = aAngle_Right + (23.01 - 23.00);
			 else if((aRunNo >= 3735 && aRunNo <= 3891) && fabs(aAngle_Right - 25.00)<0.5)
				 *aAngle = aAngle_Right + (25.00 - 25.00);
			 else if((aRunNo >= 3892 && aRunNo <= 3916) && fabs(aAngle_Right - 21.00)<0.5)
				 *aAngle = aAngle_Right + (21.03 - 21.00);
			 else if((aRunNo >= 3917 && aRunNo <= 4071) && fabs(aAngle_Right - 28.00)<0.5)
				 *aAngle = aAngle_Right + (27.99 - 28.00); 
			 else if((aRunNo >= 4073 && aRunNo <= 4101) && fabs(aAngle_Right - 28.00)<0.5)
				 *aAngle = aAngle_Right + (27.99 - 28.00);                                        
			 else if((aRunNo >= 4112 && aRunNo <= 4179) && fabs(aAngle_Right - 23.00)<0.5)
				 *aAngle = aAngle_Right + (23.04 - 23.00);                                        
			 else if((aRunNo >= 4181 && aRunNo <= 4241) && fabs(aAngle_Right - 25.00)<0.5)
				 *aAngle = aAngle_Right + (25.00 - 25.00);                                        
			 else if((aRunNo >= 4242 && aRunNo <= 4250) && fabs(aAngle_Right - 21.00)<0.5)
				 *aAngle = aAngle_Right + (21.03 - 21.00);                                        
			 else if((aRunNo >= 4251 && aRunNo <= 4299) && fabs(aAngle_Right - 28.00)<0.5)
				 *aAngle = aAngle_Right + (27.99 - 28.00);                                        
			 else{
				 *aAngle = aAngle_Right;
				 cerr   <<Form("      *** ERROR - In XGT2_Physics.h, RunNo#%d is not in the table!!!",aRunNo)<<endl;
				 outlog <<Form("      *** ERROR - In XGT2_Physics.h, RunNo#%d is not in the table!!!",aRunNo)<<endl;
			 }
  			 /*Old Correction{{{*/
		   /*
			 if(bAngle_Correct){
				 if(fabs(aAngle_Right-21.00)<0.5)
					 *aAngle = aAngle_Right + (20.971107-21.00); //From Pointing Study

				 if(fabs(aAngle_Right-23.00)<0.5)
					 *aAngle = aAngle_Right + (22.969339-23.00); //From Pointing Study

				 if(fabs(aAngle_Right-25.00)<0.5)
					 *aAngle = aAngle_Right + (24.967436-25.00); //From Pointing Study

				 if(fabs(aAngle_Right-28.00)<0.5)
					 *aAngle = aAngle_Right + (27.964330 - 28.00); //From Pointing Study
			 }
			   */
		   /*}}}*/
		 }
		 break;
	  }
  }
  aInfile.close();
}
/*}}}*/

/*inline double gLoad_P0Table(const int* aRunNo, const TString& aArm){{{*/
inline double gLoad_P0Table(const int aRunNo, const TString& aArm, double* aMom_Set){

	ifstream aInfileL(Form("%s/Table/E08014_Mom_Left.table",XGT2_DIR.Data()));
	ifstream aInfileR(Form("%s/Table/E08014_Mom_Right.table",XGT2_DIR.Data()));

	TString aCom;
	int aRunNo_List = 0;
	double aP0_List = 0.0, aP0_Set = 0.0;
	TString aTarget_List, aKin_List;
	double aP0 = -1.00;

	if(aArm=="L"){
		aCom.ReadLine(aInfileL);
		while(!aInfileL.eof()){
			aInfileL >> aRunNo_List >> aP0_List >> aP0_Set >> aTarget_List >> aKin_List;
			if(aRunNo == aRunNo_List){
				aP0 = aP0_List;
				*aMom_Set = aP0_Set;
				if(aP0<1.0 || fabs(aP0_List-aP0_Set)>0.05 ) //To avoid bad calculation of field from Dipole
					aP0 = aP0_Set;
				break;
			}
		}
		aInfileL.close();
	}

	if(aArm=="R"){
		aCom.ReadLine(aInfileR);
		while(!aInfileR.eof()){
			aInfileR >> aRunNo_List >> aP0_List >> aP0_Set >> aTarget_List >> aKin_List;
			if(aRunNo == aRunNo_List){
				aP0 = aP0_List;
				*aMom_Set = aP0_Set;
				if(aP0<1.0 || fabs(aP0_List-aP0_Set)>0.05)//To avoid bad calculation of field from Dipole
					aP0 = aP0_Set;
				break;
			}
		}
		aInfileR.close();
	}

	return aP0;
}
/*}}}*/

/*inline void gCheck_P0(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aE_Tree_Name,double* aP0Chain){{{*/
inline void gCheck_P0(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aE_Tree_Name,double* aP0Chain)
{

	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		Int_t aRunNo=aRunNoChain[i];
		/*
		   TChain *aE_Tree = gAddTree(aRunNo,aE_Tree_Name);
		   aE_Tree->Add(Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo));
		   Double_t aDipole;
		   XGT2_VAR* aP0=gCal_P0(aE_Tree,aArm,aDipole);
		   aP0Chain[i] = aP0->Value;
		   */

		double aE0, aDEE, aAngle, aMom,aMom_Set;
		gLoad_EpicsTable(aRunNo, aArm.Data(), &aE0, &aDEE, &aMom, &aAngle);
		aP0Chain[i] = aMom;

		double aP0 = gLoad_P0Table(aRunNo, aArm.Data(), &aMom_Set);
		if(fabs(aP0-aMom_Set)>0.005)
			aP0Chain[i] = aMom;
        else
			aP0Chain[i] = aP0;	

		outlog << Form("      For Run#%d, P0_C = %6.4f, P0_T = %6.4f",aRunNo,aP0, aMom)<<endl; 
		cerr   << Form("      For Run#%d, P0_C = %6.4f, P0_T = %6.4f",aRunNo,aP0, aMom)<<endl; 

		//		   delete aE_Tree;delete aP0;
	}
}
/*}}}*/

/*inline void gCheck_Angle(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aT_Tree_Name, double* aAngleChain){{{*/
inline void gCheck_Angle(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aT_Tree_Name, double* aAngleChain)
{
	Int_t aRunNo;
	double aE0, aDEE, aAngle, aMom;
	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		aRunNo=aRunNoChain[i];
		gLoad_EpicsTable(aRunNo, aArm.Data(), &aE0, &aDEE, &aMom, &aAngle);

		if(aArm=="L") 
			aAngleChain[i] = aAngle; 
		else if(aArm=="R") 
			aAngleChain[i] = -aAngle;//Note: Right Arm angle has to be negative!!! 

		outlog << Form("      For Run#%d, Angle = %6.4f ",aRunNo,aAngle)<<endl; 
		cerr << Form("      For Run#%d, Angle = %6.4f ",aRunNo,aAngle)<<endl; 
	}
}
/*}}}*/

/*inline void gCheck_E0(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aT_Tree_Name, double* aAngleChain){{{*/
inline void gCheck_E0(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aT_Tree_Name, double* aE0Chain, double* aDEEChain)
{
	Int_t aRunNo;
	double aE0, aDEE, aAngle, aMom;
	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		aRunNo=aRunNoChain[i];
		gLoad_EpicsTable(aRunNo, aArm.Data(), &aE0, &aDEE, &aMom, &aAngle);

		aE0Chain[i] = 3.356; 
		//aE0Chain[i] = aE0/1000.; //From MeV to GeV
		aDEEChain[i] = aDEE/1000.;//From MeV to GeV

		cerr << Form("      For Run#%d, E0 = %6.4f, dEE = %6.4f",aRunNo,aE0Chain[i],aDEEChain[i])<<endl; 
		outlog << Form("      For Run#%d, E0 = %6.4f, dEE = %6.4f",aRunNo,aE0Chain[i],aDEEChain[i])<<endl; 
	}

}
/*}}}*/

/*inline Double_t gGet_Angle_Corrected(const vector<int>& aRunNoChain, const TString& aArm, double aAngle){{{*/
inline Double_t gGet_Angle_Corrected(const vector<int>& aRunNoChain, const TString& aArm, double aAngle){

	double aAngle_Corrected = 0.0;

	if(aArm=="L"){
		if(fabs(aAngle-21.00)<0.5)
			aAngle_Corrected = 21.214002; //From Pointing Study

		if(fabs(aAngle-23.00)<0.5)
			aAngle_Corrected = 23.212207; //From Pointing Study

		if(fabs(aAngle-25.00)<0.5)
			aAngle_Corrected = 25.210281; //From Pointing Study

		if(fabs(aAngle-28.00)<0.5)
			aAngle_Corrected = 28.207145; //From Pointing Study
	}
	if(aArm=="R"){

		if(fabs(aAngle-21.00)<0.5)
			aAngle_Corrected = 20.971107; //From Pointing Study

		if(fabs(aAngle-23.00)<0.5)
			aAngle_Corrected = 22.969339; //From Pointing Study

		if(fabs(aAngle-25.00)<0.5)
			aAngle_Corrected = 24.967436; //From Pointing Study

		if(fabs(aAngle-28.00)<0.5)
			aAngle_Corrected = 27.964330; //From Pointing Study
	}

	return aAngle_Corrected;
}
/*}}}*/
