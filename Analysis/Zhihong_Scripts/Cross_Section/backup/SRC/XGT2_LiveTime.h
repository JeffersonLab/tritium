/////////////////////////////////
// PreScale Factors
/////////////////////////////////
//Swape Sys. and Stat. Errors, 12/19/2013
//
//
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

/*inline XGT2_VAR* gCal_LiveTime(TTree* aT_Tree,const TString& aDAQ_Cut,TTree* aS_Tree,const TString& aRaw_Branch){{{*/
inline XGT2_VAR* gCal_LiveTime(TTree* aT_Tree,const Int_t aRunNo,const TString& aDAQ_Cut,TTree* aS_Tree,const TString& aRaw_Branch,const TString& aI_Name, Int_t aPS)
{
	Int_t DAQ_Trigger=(Int_t)aT_Tree->GetEntries(aDAQ_Cut.Data());
	DAQ_Trigger *=aPS; //Corrected if PS is larger than one
	aT_Tree->SetNotify(0);
	Double_t Raw_Trigger,Raw_Trigger_Mean,rms;

	// Run "Scaler" to add branch rt_t?s for #counts in each scaler time windows, and the total counts after beam trip cut will just be adding up those counts
	Int_t Raw_Trigger_Count = aS_Tree -> GetEntries(Form("%s>=%f",aI_Name.Data(),Beam_Trip_Cut));
	gGet_Leaf_Mean_and_RMS(aS_Tree,aRaw_Branch.Data(), Form("%s>=%f",aI_Name.Data(),Beam_Trip_Cut),Raw_Trigger_Mean,rms); 
	Raw_Trigger = Raw_Trigger_Mean * Raw_Trigger_Count;

	//  outlog<<Form("DAQ = %d, RAW = %f",DAQ_Trigger,Raw_Trigger)<<endl;

	XGT2_VAR* LiveTime=new XGT2_VAR();
	if ( Raw_Trigger>1e-16 && DAQ_Trigger>1e-16)
	{
		LiveTime->Value=(Double_t) DAQ_Trigger/Raw_Trigger;
		LiveTime->Sys_Err=LiveTime->Value*sqrt(fabs(1./DAQ_Trigger-1./Raw_Trigger));
		LiveTime->Stat_Err= 0.0;
		//LiveTime->Stat_Err=LiveTime->Value*sqrt(1./DAQ_Trigger+1./Raw_Trigger/Raw_Trigger);
	}
	return LiveTime;
}
/*}}}*/

/*inline XGT2_VAR* gCal_LiveTime_Avg(TTree* aT_Tree,const TString& aDAQ_Cut,TTree* aS_Tree,const TString& aRaw_Branch){{{*/
inline XGT2_VAR* gCal_LiveTime_Avg(TTree* aT_Tree,const Int_t aRunNo,const TString& aDAQ_Cut,TTree* aS_Tree,const TString& aRaw_Branch, Int_t aPS)
{
	Int_t DAQ_Trigger=(Int_t)aT_Tree->GetEntries(aDAQ_Cut.Data());
	DAQ_Trigger *=aPS; //Corrected if PS is larger than one
	aT_Tree->SetNotify(0);
	Double_t Raw_Trigger;
	aS_Tree->SetBranchAddress(aRaw_Branch.Data(),&Raw_Trigger);
	Int_t scalerT_nentries=aS_Tree->GetEntries();
	aS_Tree->GetEntry(scalerT_nentries-1);

	XGT2_VAR* LiveTime=new XGT2_VAR();
	if ( Raw_Trigger>1e-16 && DAQ_Trigger>1e-16 )
	{
		LiveTime->Value=(Double_t) DAQ_Trigger/Raw_Trigger;
		LiveTime->Sys_Err=LiveTime->Value*sqrt(fabs(1./DAQ_Trigger-1./Raw_Trigger)); //Z. Ye, fix the "+" into "-" since DAQ_Trigger and RAW_Trigger are correlated. 01/15/2015
		LiveTime->Stat_Err= 0.0;
		//LiveTime->Stat_Err=LiveTime->Value*sqrt(1./DAQ_Trigger+1./Raw_Trigger/Raw_Trigger);
	}
	return LiveTime;
}
/*}}}*/

/*inline XGT2_VAR* gGet_LiveTime_Avg(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aCuts){{{*/
inline XGT2_VAR* gGet_LiveTime_Avg(const vector<Int_t>& aRunNoChain,const TString& aArm, const Int_t* aPSChain, const Double_t aBeam_Trip_Cut)
{
	Int_t aPS[8], aMain_Trigger=0; 
	Int_t index=0;
	Double_t aLiveTime_Avg=0.0,aLiveTime_Avg_Err=0.0;
	TString aS_Tree_Name, aI_Name;

	XGT2_VAR* aLiveTime=new XGT2_VAR();

	if(aArm =="L"){
		aMain_Trigger = 3;
		aS_Tree_Name = "LEFT"; 
		aI_Name = "left_current";
		//aS_Tree_Name = "RIGHT"; 
		//aI_Name = "right_current";
	}
	else if(aArm =="R"){
		aMain_Trigger = 1;
		aS_Tree_Name = "RIGHT";    
		aI_Name = "right_current";//LEFT BCM had problem for x>2 experiment
	}
	else {
		aLiveTime->Value = -1.0;
		aLiveTime->Stat_Err =  -1.0;
		return aLiveTime;
	}

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
			TString aDAQ_Cut = Form("(DBB.evtypebits&(1<<%d))>0 && %s>%f",aMain_Trigger, aI_Name.Data(), aBeam_Trip_Cut);
			XGT2_VAR* aLT=gCal_LiveTime_Avg(aT_Tree,aRunNo,aDAQ_Cut.Data(),aS_Tree,Form("t%dc",aMain_Trigger),aPS[aMain_Trigger-1]);
			//XGT2_VAR* aLT=gCal_LiveTime(aT_Tree,aRunNo,aDAQ_Cut.Data(),aS_Tree,Form("rt_t%dc",aMain_Trigger),aI_Name.Data(),aPS[aMain_Trigger-1]);

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

///////////////////////////////////
// Live Time
///////////////////////////////////
/*inline void gLoad_LiveTime_Chain(const vector<Int_t>& aRunNoChain,const TString& aArm){{{*/
inline void gLoad_LiveTime_Chain(const vector<Int_t>& aRunNoChain,const TString& aArm, double *aLiveTimeChain,Double_t* aLiveTimeStatErrorChain,Double_t* aLiveTimeSysErrorChain, const Int_t aFlag)
{ 
	Double_t kLT_Avg, kLT_Avg_Err, kLT, kLT_Err, kI_Max, kI_Mean, kRate;
	Int_t kRunNo,kPS;
	ifstream table;
	Double_t aLT, aLT_Stat_Err,aLT_Sys_Err;

	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		if(aArm=="R")
			table.open(Form("%s/Table/Right_LiveTime.table",XGT2_DIR.Data()));
		else if(aArm=="L")
			table.open(Form("%s/Table/Left_LiveTime.table",XGT2_DIR.Data()));
		TString aCom; aCom.ReadLine(table);

		Int_t aRunNo=aRunNoChain[i];
		aLT=-1.0;  aLT_Stat_Err=0.0; aLT_Sys_Err=0.0;
		while(!(table.eof())){
			table >> kRunNo >> kPS >> kLT_Avg >> kLT_Avg_Err >> kLT >> kLT_Err >> kRate >> kI_Max >> kI_Mean;
			//outlog<<Form("%d %f %f %f %f %f %f",kRunNo,kLT_Avg,kLT_Avg_Err,kLT,kLT_Err,kI_Max,kI_Mean)<<endl;
			if(kRunNo==aRunNo){
				if(aFlag==1){
					aLT= kLT;
					aLT_Sys_Err= kLT_Err;
					aLT_Stat_Err= 0.00;
					if(aLT>1.5||aLT<1e-33){
						aLT= kLT_Avg;
						aLT_Sys_Err= kLT_Avg_Err;
						aLT_Stat_Err = 0.00;
					}
				}
				else if(aFlag==2){
					aLT= kLT_Avg;
					aLT_Sys_Err= kLT_Avg_Err;
					aLT_Stat_Err = 0.00;
				}
				else{
					aLT= -1.00;
					aLT_Stat_Err= 0.00;
					aLT_Sys_Err = 0.00;
				}
				//In case LEFT BCM go crazy, use the average livetime
				break;
			}
		}//	while(!(table.eof()))
		table.close();

		if(aLT<=1.2&&aLT>0.4){
			outlog <<Form("      For Run#%d, the LiveTime = %5.4f",aRunNo,aLT)<<endl;
			cerr <<Form("      For Run#%d, the LiveTime = %5.4f",aRunNo,aLT)<<endl;
		}
		else{ 
            Int_t aPS[8];
			TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),aRunNo);
			if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
			{         
				//Prescale Factor
				gGet_PS(aRunNo, aPS);
								//LiveTime
				Int_t aMain_Trigger=0;
				TString aS_Tree_Name = "FAKE"; 
				if(aArm=="L"){
                  aMain_Trigger=3; aS_Tree_Name = "LEFT";
				}
				else if(aArm=="R"){
				  aMain_Trigger=1; aS_Tree_Name = "RIGHT";
				}
				else
					cerr<<"***ERROR, Wrong Arm in Livetime Calculation!"<<endl;

				TChain *aT_Tree = gAddTree(aRunNo,"T");
				TChain *aS_Tree = gAddTree(aRunNo,aS_Tree_Name.Data());

				TString aDAQ_Cut = Form("(DBB.evtypebits&(1<<%d))>0",aMain_Trigger);
				XGT2_VAR* aLiveTime=gCal_LiveTime_Avg(aT_Tree,aRunNo,aDAQ_Cut.Data(),aS_Tree,Form("t%dc",aMain_Trigger),aPS[aMain_Trigger-1]);
				//aLT = 1.00; aLT_Stat_Err=0.0; aLT_Sys_Err=0.01;
				aLT = aLiveTime->Value; aLT_Stat_Err=aLiveTime->Stat_Err; aLT_Sys_Err=aLiveTime->Sys_Err;
				outlog <<Form("***ERROR, Can not find LT for Run#%d, re-calculate the LiveTime = %5.4f!",aRunNo,aLT)<<endl;
				  cerr <<Form("***ERROR, Can not find LT for Run#%d, re-calculate the LiveTime = %5.4f!",aRunNo,aLT)<<endl;
			}
		}
		aLiveTimeChain[i] = aLT;
		aLiveTimeStatErrorChain[i] = aLT_Stat_Err;
		aLiveTimeSysErrorChain[i] = aLT_Sys_Err;
	}//for ( unsigned int i=0; i<aRunNoChain.size(); i++ )

}
/*}}}*/

/*inline void gCheck_LiveTime(const vector<Int_t>& aRunNoChain,const TString& aArm){{{*/
inline void gCheck_LiveTime(const vector<Int_t>& aRunNoChain,const TString& aArm, double *aLiveTimeChain,Double_t* aLiveTimeErrorChain, const int* aPSChain, const Double_t aBeam_Trip_Cut)
{
	Int_t aPS[8], aMain_Trigger=0; 
	TString aS_Tree_Name, aI_Name;

	if(aArm =="L"){ //LEFT Arm BCM has problem
		aMain_Trigger = 3;
		//aS_Tree_Name = "LEFT"; 
		//aI_Name = "left_current";
		aS_Tree_Name = "RIGHT"; 
		aI_Name = "right_current";
	}
	else if(aArm =="R"){
		aMain_Trigger = 1;
		aS_Tree_Name = "RIGHT";
		aI_Name = "right_current";
	}
	else 
		outlog <<"   XXXXXXXXXXXXXXX Wrong ARM! Stop now!!!"<<endl;
		cerr <<"   XXXXXXXXXXXXXXX Wrong ARM! Stop now!!!"<<endl;

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
			TString aDAQ_Cut = Form("(DBB.evtypebits&(1<<%d))>0 && %s>%f",aMain_Trigger, aI_Name.Data(), aBeam_Trip_Cut);
			XGT2_VAR* aLT=gCal_LiveTime(aT_Tree,aRunNo,aDAQ_Cut.Data(),aS_Tree,Form("rt_t%dc",aMain_Trigger),aI_Name.Data(),aPS[aMain_Trigger-1]);

			if(aLT->Value<0.1){
				if(aArm =="L")
					aS_Tree_Name = "LEFT"; 
				delete aS_Tree;  TChain *aS_Tree = gAddTree(aRunNo,aS_Tree_Name.Data());
				TString aDAQ_Cut = Form("(DBB.evtypebits&(1<<%d))>0",aMain_Trigger);
				aLT=gCal_LiveTime_Avg(aT_Tree,aRunNo,aDAQ_Cut.Data(),aS_Tree,Form("t%dc",aMain_Trigger),aPS[aMain_Trigger-1]);
			}

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
