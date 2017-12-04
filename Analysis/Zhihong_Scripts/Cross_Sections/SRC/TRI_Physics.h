#define NOCASE 1

/*Double_t Get_Mott_CrossSection(const Double_t& aAngle,const Double_t& aEi){{{*/
/*Double_t Get_Mott_CrossSection(const Double_t& aAngle,const Double_t& aEi) 
{
  return ALPHA*ALPHA*pow(cos(aAngle/2*TMath::DegToRad()),2)/4/aEi/aEi/pow(sin(aAngle/2*TMath::DegToRad()),4);// unit=1/(MeV^2*sr)
}*/
/*}}}*/

/*inline TRI_VAR* gCal_P0(TTree* aTree,const TString& aArm,const Double_t& aSpec_Mom,Double_t& oDipole){{{*/
//Fix Here,ZYE
inline TRI_VAR* gCal_P0(TTree* aTree,const TString& aArm,Double_t& oDipole)
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

  TRI_VAR* output=new TRI_VAR();
  output->SetValue(P0/1000,P0/1000*1e-4,0);
 
  return output;//MeV to GeV
}
/*}}}*/

/*inline void gCheck_P0(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aE_Tree_Name,double* aP0Chain){{{*/
inline void gCheck_P0(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aE_Tree_Name,double* aP0Chain)
{

	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		Int_t aRunNo=aRunNoChain[i];
        
        TChain *aE_Tree = gAddTree(aRunNo,aE_Tree_Name);
        aE_Tree->Add(Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo));
        Double_t aDipole;
        TRI_VAR* aP0=gCal_P0(aE_Tree,aArm,aDipole);
        aP0Chain[i] = aP0->Value;

		outlog << Form("      For Run#%d, P0 = %6.4f",aRunNo,aP0->Value)<<endl; 
		cerr   << Form("      For Run#%d, P0 = %6.4f",aRunNo,aP0->Value)<<endl; 

        delete aE_Tree;
        delete aP0;
    }
}
/*}}}*/

/*inline void gCheck_E0(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aE_Tree_Name, double* aAngleChain){{{*/
inline void gCheck_E0(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aE_Tree_Name, double* aE0Chain, double* aDE0Chain)
{
	Int_t aRunNo;
	double aE0, aDE0;
	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		aRunNo=aRunNoChain[i];
        TChain *aE_Tree = gAddTree(aRunNo,aE_Tree_Name);

        Int_t nentries=(Int_t)aE_Tree->GetEntries();
        aE_Tree->SetBranchAddress("HALLA:p", &aE0); //Beam Energy
        aE_Tree->SetBranchAddress("HALLA:dpp", &aDE0); //Beam Energy Spread

        ////Take the last EPCIS value
        //aE_Tree->GetEntry(nentries-1);
        //aE0Chain[i] = aE0;
        //aDE0Chain[i] = aDE0;

        //Take the average values of the entire run, maybe an issue if beam energy changes dramatically during one run
        aE0Chain[i] = 0.0;
        aDE0Chain[i] = 0.0;
        for(int j=0;j<nentries;j++){
            aE_Tree->GetEntry(j);
            aE0Chain[i] += aE0;
            aDE0Chain[i] += aDE0;
        }
        aE0Chain[i] /= nentries;
        aDE0Chain[i] /= nentries;

        cerr << Form("      For Run#%d, E0 = %6.4f, dEE = %6.4f",aRunNo,aE0Chain[i],aDE0Chain[i])<<endl; 
        outlog << Form("      For Run#%d, E0 = %6.4f, dEE = %6.4f",aRunNo,aE0Chain[i],aDE0Chain[i])<<endl; 
	}

}
/*}}}*/
