/*inline TString gGet_Delta_Correct(const TString& aArm, const bool aIsExtTgt){{{*/
inline TString gGet_Delta_Correct(const TString& aArm, const bool aIsExtTgt){

	//Addtional Correction for Delta on HRS-R, since the D-Term in optics was not calibrated.
	//The parameters were fitted from SAMC simulation data.
	//--- Z. Ye, 02/18/2013
	//===============================================================================
    //--X_fp:   
	const double X0=3.034150e-04,  X1=-6.629199e-03, X2=-8.535028e-04;
    //--Theta_fp: 
	const double T0=-1.060184e-04, T1=1.061409e-01,  T2=6.751897e-01;
    //--Phi_fp:   
	const double P0=6.492564e-05,  P1=2.171117e-05,  P2=-3.002779e-01;
     //--Y_fp:   
	const double Y0=2.933949e-05,  Y1=-5.402877e-04, Y2=-4.617217e-01, Y3=7.100115e-01, Y4=5.226948e+02;
   //--Delta:   
	const double D0=2.668049e-05,  D1=-4.939113e-03, D2=-5.875315e-02,  D3=2.054553e+00;
    //===============================================================================

	TString aEX_Dp, aFP_Th, aFP_Ph, aFP_X, aFP_Y;
	if ( aIsExtTgt ){
		aEX_Dp= Form("ExTgt%s.dp",aArm.Data());
	}
	else{
		aEX_Dp =Form("%s.gold.dp",aArm.Data());
	}
	aFP_X = Form("%s.tr.r_x",aArm.Data());
	aFP_Y = Form("%s.tr.r_y",aArm.Data());
	aFP_Th = Form("%s.tr.r_th",aArm.Data());
	aFP_Ph = Form("%s.tr.r_ph",aArm.Data());

	TString aDp_Corrected = Form("(%s -(%e+%e*pow(%s,1)+%e*pow(%s,2)+%e*pow(%s,3))-(%e+%e*pow(%s,1)+%e*pow(%s,2))-(%e+%e*pow(%s,1)+%e*pow(%s,2))-(%e+%e*pow(%s,1)+%e*pow(%s,2))-(%e+%e*pow(%s,1)+%e*pow(%s,2)+%e*pow(%s,3)+%e*pow(%s,4)))",	
									    aEX_Dp.Data(),D0, D1, aEX_Dp.Data(), D2, aEX_Dp.Data(), D3, aEX_Dp.Data(),
										X0, X1, aFP_X.Data(),  X2, aFP_X.Data(),
										T0, T1, aFP_Th.Data(), T2, aFP_Th.Data(),
										P0, P1, aFP_Ph.Data(), P2, aFP_Ph.Data(),
										Y0, Y1, aFP_Y.Data(),  Y2, aFP_Y.Data(), Y3, aFP_Y.Data(), Y4, aFP_Y.Data() );
   return aDp_Corrected.Data();
}
/*}}}*/

/*XGT2_VAR* gGet_Nf_EX(const vector<Int_t>& aRunNoChain,const TString& aArm, const Double_t aEp_Bin, const Double_t aEp_Size, const TString& aElectronCuts, const Double_t* aLiveTimeChain,const Double_t* aLiveTimeStatErrChain,const Double_t* aLiveTimeSysErrChain, const Int_t* aPSChain, const Double_t* aE0Chain, const Double_t* aP0Chain, const Double_t* aAngleChain, double* aNf_EX_Chain,double* aNf_EX_Stat,double* aNf_EX_Sys, const bool aIsExtTgt){{{*/
XGT2_VAR* gGet_Nf_EX(const vector<Int_t>& aRunNoChain,const TString& aArm, const TString& aElectronCuts, const double aVZ, const double aVZ_Step, const Double_t* aLiveTimeChain,const Double_t* aLiveTimeStatErrChain,const Double_t* aLiveTimeSysErrChain, const Int_t* aPSChain, double* aNf_EX_Chain,double* aNf_EX_Stat,double* aNf_EX_Sys){
    XGT2_VAR* nf = new XGT2_VAR();

	Int_t aPS[8],aNf[aRunNoChain.size()];
	Double_t aNf_Stat_Err[aRunNoChain.size()], aNf_Sys_Err[aRunNoChain.size()];
	Int_t aMain_Trigger; 
	Double_t aNf_Stat_Err_Total, aNf_Sys_Err_Total,aNf_Total;
	TString aS_Tree_Name;
	aNf_Total = 0.; aNf_Stat_Err_Total =0.0; aNf_Sys_Err_Total=0.0;

	TString aEX_VZ, aEX_Dp, aEX_Th, aEX_Ph;
	TString aEX_Current ="", aEX_Cuts = "";
	Int_t aRunNo= 0 ;
	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		aRunNo=aRunNoChain[i];

		TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),aRunNo);
		if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
		{
			/////////////////////////////////////////////
			// Define Acceptance Variables
			////////////////////////////////////////////
			/*Define Acceptance Variables{{{*/
			//Apply Delta Correction in the data -- Z. Ye 02/18/2013
			aEX_VZ =Form("RctPt%s.z",aArm.Data());
			aEX_Dp = gGet_Delta_Correct(aArm.Data(), 1);
			if(aArm=="R")
				aEX_Dp= Form("ExTgt%s.dp",aArm.Data());
			aEX_Th = Form("ExTgt%s.th",aArm.Data());
			aEX_Ph = Form("ExTgt%s.ph",aArm.Data());

			//Current Branch for cutting out Beam Trip events
			// must run "AddCurrent" first to add those branches in rootfiles
			if(aArm=="R"||aArm=="L")
				aEX_Current = "right_current"; //LEFT BCM had problem for x>2 experiment
			else
				aEX_Current = "left_current";

			aEX_Cuts=Form("abs(%s-%f)<=%f&&%s>=%f",
					aEX_VZ.Data(),
					aVZ,
					aVZ_Step,
					aEX_Current.Data(),
					Beam_Trip_Cut
					);
			aEX_Cuts+="&&"+aElectronCuts;

			if(aArm =="L"){
				aMain_Trigger = 3;
				aS_Tree_Name = "LEFT"; 
			}
			else if(aArm =="R"){
				aMain_Trigger = 1;
				aS_Tree_Name = "RIGHT"; 
			}
			else
				outlog << Form("   *** ERROR, I cann't understand the name %s ARM", aArm.Data())<<endl;

			for(int j=0;j<8;j++)
				aPS[j] = aPSChain[i*8+j];
			
			TChain *aT_Tree = gAddTree(aRunNo,"T");

			aNf[i] = aPS[aMain_Trigger-1]*aT_Tree->GetEntries(aEX_Cuts.Data());
	        //Ignore those runs that only contribute events less or equal to 0 counts otherwise the stat_err will blow up,2012/10/08
			if(aNf[i]>1e-16&&aLiveTimeChain[i]>1e-16){
                aNf_EX_Chain[i] = aNf[i]/aLiveTimeChain[i];
				aNf_Stat_Err[i]= aNf_EX_Chain[i] * sqrt(1.0/aNf[i] + pow(aLiveTimeStatErrChain[i]/aLiveTimeChain[i] ,2));
				aNf_Sys_Err[i]=aLiveTimeSysErrChain[i];

				aNf_Total+= aNf_EX_Chain[i];
				aNf_Stat_Err_Total += pow(aLiveTimeStatErrChain[i]/aLiveTimeChain[i],2);
				aNf_Sys_Err_Total  += pow(aNf_Sys_Err[i] /aNf_EX_Chain[i],2);

				aNf_EX_Stat[i] = aNf_Stat_Err[i];
				aNf_EX_Sys[i] = aNf_Sys_Err[i];
				outlog << Form("   === RunNo=%d, PS%d=%d, LiveTime/Err=%4.3f / %4.3f, Events/Err=%7f / %e",
						aRunNo,aMain_Trigger,aPS[aMain_Trigger-1],aLiveTimeChain[i], aLiveTimeStatErrChain[i],aNf_EX_Chain[i],aNf_EX_Stat[i])<<endl;
				cerr   << Form("   === RunNo=%d, PS%d=%d, LiveTime/Err=%4.3f / %4.3f, Events/Err=%7f / %e",
						aRunNo,aMain_Trigger,aPS[aMain_Trigger-1],aLiveTimeChain[i], aLiveTimeStatErrChain[i],aNf_EX_Chain[i],aNf_EX_Stat[i])<<endl;
			}
			delete aT_Tree;
		}
		else{
			outlog << Form("   **** Cann't find Rootfile %s in %s", File_Form.Data(),ROOTFILES_DIR.Data())<<endl;
			cerr   << Form("   **** Cann't find Rootfile %s in %s", File_Form.Data(),ROOTFILES_DIR.Data())<<endl;
		}
	}   
	nf->Value    = aNf_Total;
	nf->Stat_Err = aNf_Total * sqrt(1.0/aNf_Total + aNf_Stat_Err_Total);
	nf->Sys_Err  = aNf_Total * sqrt(aNf_Sys_Err_Total);
    return nf;
}
/*}}}*/
