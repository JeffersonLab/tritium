/*inline TRI_VAR* gGetLiveTime(const Int_t aRunNo, const TString& aArm){{{{*/
inline TRI_VAR* gGetLiveTime(const Int_t aRunNo, const TString& aArm){
    TRI_VAR* aLiveTime = new TRI_VAR();
    
    TString aScaler_Name="", aRaw_Branch = "";
    int aMain_Trigger=0;
    if(aArm=="L"){
        aScaler_Name = LEFT_SCALER_TREE;
        aMain_Trigger = MAIN_TRIGGER_LEFT;
        aRaw_Branch = Form("%s%d", LEFT_TRIGGER_SCALER.Data(), aMain_Trigger);
    
    }else if(aArm=="R"){
        aScaler_Name = RIGHT_SCALER_TREE;
        aMain_Trigger = MAIN_TRIGGER_RIGHT;
        aRaw_Branch = Form("%s%d", RIGHT_TRIGGER_SCALER.Data(), aMain_Trigger);
    
    }else{
        cerr<<Form("*** ERROR, I don't know the Arm=%s **** ", aArm.Data())<<endl;
        aLiveTime->Value = 0.0;
        aLiveTime->Sys_Err = 0.0;
        aLiveTime->Stat_Err = 0.0;
    }

    TString aDAQ_Cut = Form("D%s.evtypebits>>%d&1",aArm.Data(), aMain_Trigger);
    TTree* aT_Tree = gAddTree(aRunNo, "T");
   	Int_t aDAQ_Trigger=(Int_t)aT_Tree->GetEntries(aDAQ_Cut.Data());
    //	DAQ_Trigger *=aPS; //Corrected if PS is larger than one
    aT_Tree->SetNotify(0);
    
    TTree* aS_Tree = gAddTree(aRunNo, aScaler_Name.Data());
	Double_t aRaw_Trigger;
	aS_Tree->SetBranchAddress(aRaw_Branch.Data(),&aRaw_Trigger);
	Int_t aScaler_Nentries=aS_Tree->GetEntries();
	aS_Tree->GetEntry(aScaler_Nentries-1);

	if ( aRaw_Trigger>1e-16 && aDAQ_Trigger>1e-16 )
	{
		aLiveTime->Value=(Double_t) aDAQ_Trigger/aRaw_Trigger;
		aLiveTime->Sys_Err=aLiveTime->Value*sqrt(fabs(1./aDAQ_Trigger-1./aRaw_Trigger)); //Z. Ye, fix the "+" into "-" since DAQ_Trigger and RAW_Trigger are correlated. 01/15/2015
		aLiveTime->Stat_Err= 0.0;
	}
	return aLiveTime;
}
/*}}}*/

/*inline void gGetGoodEvents(const vector<Int_t>& aRunNoChain, const TString& aArm, const TString& aCut_e, const TString& aCut_pi, const TString& aCut, double* aNeleChain, double* aNpiChain, double* aNtotChain ){{{{*/
inline void gGetGoodEvents(const vector<Int_t>& aRunNoChain, const TString& aArm, const TString& aCut_e, const TString& aCut_pi, const TString& aCut_tot, double* aNeleChain, double* aNpiChain, double* aNtotChain, double* aLiveTimeChain, double* aLiveTimeErrChain){
    Int_t aN_ele=0, aN_pi=0, aN_tot=0; 
    for (unsigned int i=0; i<aRunNoChain.size(); i++ ){
        Int_t aRunNo=aRunNoChain[i];
        TChain *aTree = gAddTree(aRunNo,"T"); 

        aN_ele = aTree->GetEntries(TCut(aCut_e));
        aN_pi  = aTree->GetEntries(TCut(aCut_pi));
        aN_tot = aTree->GetEntries(TCut(aCut_tot));

        TRI_VAR* aLT = gGetLiveTime(aRunNo, aArm.Data());
        aLiveTimeChain[i] = aLT->Value;
        aLiveTimeErrChain[i] = aLT->Sys_Err;

        aNeleChain[i]= aN_ele/aLiveTimeChain[i];
        aNpiChain[i] = aN_pi/aLiveTimeChain[i];
        aNtotChain[i]= aN_tot/aLiveTimeChain[i];

        cout<<Form("&&&  For Run#%d, LiveTime=%4.3f, Nele=%e, Npi=%e, Ntot=%e", 
                aRunNo, aLiveTimeChain[i], aNeleChain[i], aNpiChain[i], aNtotChain[i])<<endl;

        delete aTree;
        delete aLT;
    } // for ( Int_t i=0; i<aRunNoChain.size(); i++ )
    cout<<endl;
}
/*}}}*/
