inline TRI_VAR* gGetLiveTime(const Int_t aRunNo, const TString& aArm){
    /*Define these quantities here if you haven't done elsewhere{{{*/
    const int MAIN_TRIGGER_LEFT = 1;
    const int MAIN_TRIGGER_RIGHT = 4;
    const TString LEFT_SCALER_TREE = "evLeft";
    const TString RIGHT_SCALER_TREE= "evRight";
    const TString LEFT_TRIGGER_SCALER = "evLeftT";
    const TString RIGHT_TRIGGER_SCALER = "evRightT";
    /*}}}*/

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
