/////////////////////////////////
// PreScale Factors
/////////////////////////////////
/*inline void gGet_PS(const Int_t& RunNo,const TString& aArm, Int_t *aPar ){{{*/
inline void gGet_PS(const Int_t& aRunNo,Int_t aPar[8] )
{
	TString filename_ps = (Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo));
	if ( gSystem->FindFile(ROOTFILES_DIR.Data(),filename_ps) )
	{	 
		TFile *file_ps = new TFile(filename_ps,"r"); 
		THaRun *aRun = (THaRun*) file_ps->Get("Run_Data"); 
		TArrayI ps = aRun->GetParameters()->GetPrescales();

		for(int i=0;i<8;i++){
			aPar[i] = ps[i];
		}

		file_ps->Close();
	}
	else{
		outlog <<"File Cannot be found "<< filename_ps <<endl;
		cerr <<"File Cannot be found "<< filename_ps <<endl;
	}
	//return &aPar;  
}
/*}}}*/

/*inline void gCheck_PS(const vector<Int_t>& aRunNoChain,const TString& aArm){{{*/
inline void gCheck_PS(const vector<Int_t>& aRunNoChain, int* kPS)
{
	Int_t aPS[8];

	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		Int_t aRunNo=aRunNoChain[i];

		gGet_PS(aRunNo, aPS);   //The values of last run in the chain would be given and we need to make sure PS for all runs are the same.

		outlog <<Form("      For Run#%d,",aRunNo); 
		cerr <<Form("      For Run#%d,",aRunNo); 
		for(int j=0;j<8;j++){
			kPS[i*8+j] = aPS[j];
			outlog << Form(" PS%d = %d ",j+1,kPS[i*8+j]);
			cerr << Form(" PS%d = %d ",j+1,kPS[i*8+j]);
		}
		outlog<<endl;
		cerr<<endl;
	}
}
/*}}}*/

/*inline TRI_VAR* gCal_LiveTime_Avg(TTree* aT_Tree,const TString& aDAQ_Cut,TTree* aS_Tree,const TString& aRaw_Branch){{{*/
inline TRI_VAR* gCal_LiveTime_Avg(TTree* aT_Tree,const Int_t aRunNo,const TString& aDAQ_Cut,TTree* aS_Tree,const TString& aRaw_Branch, Int_t aPS)
{
	Int_t DAQ_Trigger=(Int_t)aT_Tree->GetEntries(aDAQ_Cut.Data());
	DAQ_Trigger *=aPS; //Corrected if PS is larger than one
	aT_Tree->SetNotify(0);
	Double_t Raw_Trigger;
	aS_Tree->SetBranchAddress(aRaw_Branch.Data(),&Raw_Trigger);
	Int_t scalerT_nentries=aS_Tree->GetEntries();
	aS_Tree->GetEntry(scalerT_nentries-1);

	TRI_VAR* LiveTime=new TRI_VAR();
	if ( Raw_Trigger>1e-16 && DAQ_Trigger>1e-16 )
	{
		LiveTime->Value=(Double_t) DAQ_Trigger/Raw_Trigger;
		LiveTime->Sys_Err=LiveTime->Value*sqrt(fabs(1./DAQ_Trigger-1./Raw_Trigger)); //Z. Ye, fix the "+" into "-" since DAQ_Trigger and RAW_Trigger are correlated. 01/15/2015
		LiveTime->Stat_Err= 0.0;
	}
	return LiveTime;
}
/*}}}*/

/*inline TRI_VAR* gGet_LiveTime_Avg(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aCuts){{{*/
inline TRI_VAR* gGet_LiveTime_Avg(const vector<Int_t>& aRunNoChain,const TString& aArm, const TString& aS_Tree_Name, const int aMain_Trigger, const Int_t* aPSChain)
{
	Int_t aPS[8];
	Int_t index=0;
	Double_t aLiveTime_Avg=0.0,aLiveTime_Avg_Err=0.0;

	TRI_VAR* aLiveTime=new TRI_VAR();

    TString aTrigger_Name = "";
	if(aArm =="L"){
        aTrigger_Name=Left_Trigger_Name;
    }
	else if(aArm =="R"){
        aTrigger_Name=Right_Trigger_Name;
	}
	else {
		aLiveTime->Value = -1.0;
		aLiveTime->Stat_Err =  -1.0;
		return aLiveTime;
	}

    aTrigger_Name.ReplaceAll("_N_", Form("%d", aMain_Trigger));//eg., LefT_N_, to LeftT1, or right_T_N_ to right_T1
    for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		Int_t aRunNo=aRunNoChain[i];

		TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),aRunNo);
		if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
		{         
			// 	  //Prescale Factor
			// 	  gGet_PS(aRunNo, aPS);
			for(int j=0;j<8;j++)
				aPS[j] = aPSChain[i*8+j];

			TChain *aT_Tree = gAddTree(aRunNo,"T");
			TChain *aS_Tree = gAddTree(aRunNo,aS_Tree_Name.Data());

			//LiveTime
			TString aDAQ_Cut = Form("(%s.evtypebits&(1<<%d))>0",D_Name.Data(),aMain_Trigger);
			TRI_VAR* aLT=gCal_LiveTime_Avg(aT_Tree,aRunNo,aDAQ_Cut.Data(),aS_Tree,aTrigger_Name.Data(),aPS[aMain_Trigger-1]);

			aLiveTime_Avg += aLT->Value;
			aLiveTime_Avg_Err += aLT->Sys_Err;

			index++;
			delete aS_Tree;
			delete aT_Tree;
			delete aLT;
		}
		else{
			aLiveTime->Value = -2.0;
			aLiveTime->Sys_Err =  -2.0;
		}
		aLiveTime->Stat_Err =  0.0;

	}

	if(index!=0){
		aLiveTime->Value = aLiveTime_Avg/index;
		aLiveTime->Sys_Err =  aLiveTime_Avg_Err/index;
		aLiveTime->Stat_Err = 0.0;
	}
	else{
		aLiveTime->Value = -3.0;
		aLiveTime->Sys_Err =  -3.0;
		aLiveTime->Stat_Err =  -3.0;
	}

	return aLiveTime;
}
/*}}}*/

/*inline void gCheck_LiveTime(const vector<Int_t>& aRunNoChain,const TString& aArm){{{*/
inline void gCheck_LiveTime(const vector<Int_t>& aRunNoChain,const TString& aArm, const TString& aS_Tree_Name, const int aMain_Trigger, double *aLiveTimeChain,Double_t* aLiveTimeErrorChain, const int* aPSChain)
{
	Int_t aPS[8];
    TString aTrigger_Name = "";
	if(aArm =="L"){ //LEFT Arm BCM has problem
        aTrigger_Name=Left_Trigger_Name;
	}
	else if(aArm =="R"){
        aTrigger_Name=Right_Trigger_Name;
	}
	else{ 
		outlog <<"   *************** Wrong ARM! Stop now!!!"<<endl;
		cerr <<"   *************** Wrong ARM! Stop now!!!"<<endl;
    }

    aTrigger_Name.ReplaceAll("_N_", Form("%d", aMain_Trigger));//eg., LefT_N_, to LeftT1, or right_T_N_ to right_T1
	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		Int_t aRunNo=aRunNoChain[i];

		TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),aRunNo);
		if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
		{
			//Prescale Factor
			//gGet_PS(aRunNo, aPS);
			for(int j=0;j<8;j++)
				aPS[j] = aPSChain[i*8+j];

			TChain *aT_Tree = gAddTree(aRunNo,"T");
			TChain *aS_Tree = gAddTree(aRunNo,aS_Tree_Name.Data());

			//LiveTime
			TString aDAQ_Cut = Form("(%s.evtypebits&(1<<%d))>0",D_Name.Data(),aMain_Trigger);
			TRI_VAR* aLT=gCal_LiveTime_Avg(aT_Tree,aRunNo,aDAQ_Cut.Data(),aS_Tree,aTrigger_Name.Data(),aPS[aMain_Trigger-1]);

			if(aLT->Value<=1.0||aLT->Value<0.3){
				outlog <<Form("      For Run#%d, the LiveTime = %5.4f",aRunNo,aLT->Value)<<endl;
				cerr <<Form("      For Run#%d, the LiveTime = %5.4f",aRunNo,aLT->Value)<<endl;
			}
			else{
				outlog <<Form("***ERROR, For Run#%d, the LiveTime = %5.4f",aRunNo,aLT->Value)<<endl;
				cerr <<Form("***ERROR, For Run#%d, the LiveTime = %5.4f",aRunNo,aLT->Value)<<endl;
            }
			aLiveTimeChain[i] = aLT->Value;
			aLiveTimeErrorChain[i] = aLT->Sys_Err;

			delete aS_Tree;
			delete aT_Tree;
			delete aLT;
		}
		else{
			outlog <<"    XXXXXX Rootfile does not exist. Check please!!!"<<endl;
			cerr <<"    XXXXXX Rootfile does not exist. Check please!!!"<<endl;
		}

	}

}
/*}}}*/
