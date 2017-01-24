//#include "XGT2_Physics.h"
inline void gCheck_E0(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aT_Tree_Name, double* aE0Chain, double* aDEEChain);
inline void gCheck_Angle(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aT_Tree_Name, double* aAngleChain);
inline void gCheck_P0(const vector<Int_t>& aRunNoChain,const TString& aArm,const TString& aE_Tree_Name,double* aP0Chain);

/*Cuts on the edge of Ep{{{*/
// Cut_L_Min Cut_L_Max Cut_R_Min Cut_R_Max
const double E0 = 3.356;
const double Ep_Cut_21_1[4] = {2.80, 3.00,2.80,3.00}; 
const double Ep_Cut_21_2[4] = {2.94, 3.18,2.94,3.18}; 
const double Ep_Cut_23_1[4] = {2.74, 2.94,2.74,2.94}; 
const double Ep_Cut_23_2[4] = {2.92, 3.16,2.92,3.16};
const double Ep_Cut_QE_1[4] = {2.42, 2.58,2.42,2.58}; 
const double Ep_Cut_QE_2[4] = {2.56, 2.74,2.56,2.74}; 
const double Ep_Cut_25_1[4] = {2.70, 2.88,2.70,2.88};
const double Ep_Cut_25_2[4] = {2.88, 3.12,2.88,3.12}; 
const double Ep_Cut_28_1[4] = {2.74, 2.98,2.74,2.98}; 
/*}}}*/

/*gGet_Ep_Cut{{{*/
void gGet_Ep_Cut(const TString &aKinName, const TString & aArm, double* aEp_Min, double* aEp_Max){

	if(aKinName.Contains("Kin3.1")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_21_1[0];		
			aEp_Max[0] = Ep_Cut_21_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_21_1[2];		
			aEp_Max[0] = Ep_Cut_21_1[3];		
		}
	}
	else if(aKinName.Contains("Kin3.2")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_21_2[0];		
			aEp_Max[0] = Ep_Cut_21_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_21_2[2];		
			aEp_Max[0] = Ep_Cut_21_2[3];		
		}
	}

	else if(aKinName.Contains("Kin4.1")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_23_1[0];		
			aEp_Max[0] = Ep_Cut_23_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_23_1[2];		
			aEp_Max[0] = Ep_Cut_23_1[3];		
		}
	}
	else if(aKinName.Contains("Kin4.2")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_23_2[0];		
			aEp_Max[0] = Ep_Cut_23_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_23_2[2];		
			aEp_Max[0] = Ep_Cut_23_2[3];		
		}
	}

	else if(aKinName.Contains("Kin5.0") && ! aKinName.Contains("Kin5.05")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_QE_1[0];		
			aEp_Max[0] = Ep_Cut_QE_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_QE_1[2];		
			aEp_Max[0] = Ep_Cut_QE_1[3];		
		}
	}
	else if(aKinName.Contains("Kin5.05")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_QE_2[0];		
			aEp_Max[0] = Ep_Cut_QE_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_QE_2[2];		
			aEp_Max[0] = Ep_Cut_QE_2[3];		
		}
	}

	else if(aKinName.Contains("Kin5.1")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_25_1[0];		
			aEp_Max[0] = Ep_Cut_25_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_25_1[2];		
			aEp_Max[0] = Ep_Cut_25_1[3];		
		}
	}
	else if(aKinName.Contains("Kin5.2")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_25_2[0];		
			aEp_Max[0] = Ep_Cut_25_2[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_25_2[2];		
			aEp_Max[0] = Ep_Cut_25_2[3];		
		}
	}

	else if(aKinName.Contains("Kin6.5")){
		if(aArm=="L"){
			aEp_Min[0] = Ep_Cut_28_1[0];		
			aEp_Max[0] = Ep_Cut_28_1[1];		
		}
		else if(aArm=="R"){
			aEp_Min[0] = Ep_Cut_28_1[2];		
			aEp_Max[0] = Ep_Cut_28_1[3];		
		}
	}

	else{
		aEp_Min[0]=-1.0;
		aEp_Max[0]=10.0;
	}

}
/*}}}*/

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
    const double Dp_Offset = 0.0015; //Add offset based on the Acceptance Study
	//===============================================================================
    if(!(aArm=="R")){
	   cerr<<"*** ERROR !!! ***, Wrong Arm, must be HRS-R!"<<endl;
	   outlog<<"*** ERROR !!! ***, Wrong Arm, must be HRS-R!"<<endl;
	   return "XXX";
	}

	TString aEX_Dp, aFP_Th, aFP_Ph, aFP_X, aFP_Y;
	if ( aIsExtTgt ){
		aEX_Dp= "ExTgtR.dp";
	}
	else{
		aEX_Dp ="R.gold.dp";
	}
	aFP_X =  "R.tr.r_x";
	aFP_Y =  "R.tr.r_y";
	aFP_Th = "R.tr.r_th";
	aFP_Ph = "R.tr.r_ph";

	TString aDp_Corrected = Form("(%s -(%e+%e*pow(%s,1)+%e*pow(%s,2)+%e*pow(%s,3))-(%e+%e*pow(%s,1)+%e*pow(%s,2))-(%e+%e*pow(%s,1)+%e*pow(%s,2))-(%e+%e*pow(%s,1)+%e*pow(%s,2))-(%e+%e*pow(%s,1)+%e*pow(%s,2)+%e*pow(%s,3)+%e*pow(%s,4))-%e)",	
									    aEX_Dp.Data(),D0, D1, aEX_Dp.Data(), D2, aEX_Dp.Data(), D3, aEX_Dp.Data(),
										X0, X1, aFP_X.Data(),  X2, aFP_X.Data(),
										T0, T1, aFP_Th.Data(), T2, aFP_Th.Data(),
										P0, P1, aFP_Ph.Data(), P2, aFP_Ph.Data(),
										Y0, Y1, aFP_Y.Data(),  Y2, aFP_Y.Data(), Y3, aFP_Y.Data(), Y4, aFP_Y.Data(),
										Dp_Offset);
   return aDp_Corrected.Data();
}
/*}}}*/

/*XGT2_VAR* gGet_Nf_EX(const vector<Int_t>& aRunNoChain,const TString& aArm, const TString& aKin, const Double_t aBin, const Double_t aBin_Size, const TString& aElectronCuts,const TString& aBin_Variable, const Double_t* aLiveTimeChain,const Double_t* aLiveTimeStatErrChain,const Double_t* aLiveTimeSysErrChain, const Int_t* aPSChain, const Double_t* aE0Chain, const Double_t* aP0Chain, const Double_t* aAngleChain, double* aNf_EX_Chain,double* aNf_EX_Stat,double* aNf_EX_Sys, const bool aIsExtTgt){{{*/
XGT2_VAR* gGet_Nf_EX(const vector<Int_t>& aRunNoChain,const TString& aArm, const TString& aKin, const TString& aTarget,const Double_t aBin, const Double_t aBin_Size, const TString& aElectronCuts,const TString aBin_Variable, const Double_t* aLiveTimeChain,const Double_t* aLiveTimeStatErrChain,const Double_t* aLiveTimeSysErrChain, const Int_t* aPSChain, const Double_t* aE0Chain, const Double_t* aP0Chain, const Double_t* aAngleChain, const bool aIsExtTgt, double* aNf_EX_Chain,double* aNf_EX_Stat,double* aNf_EX_Sys){
    XGT2_VAR* nf = new XGT2_VAR();

	Int_t aPS[8],aNf[aRunNoChain.size()];
	Double_t aNf_Stat_Err[aRunNoChain.size()], aNf_Sys_Err[aRunNoChain.size()];
	Int_t aMain_Trigger; 
	Double_t aNf_Stat_Err_Total, aNf_Sys_Err_Total,aNf_Total;
	TString aS_Tree_Name;
	aNf_Total = 0.; aNf_Stat_Err_Total =0.0; aNf_Sys_Err_Total=0.0;

	Double_t aSinA =0.0, aCosA = 0.0, aSinSQ = 0.0;
	TString aEX_Dp, aEX_Th, aEX_Ph, aEX_SA, aEX_SinSQ, aEX_Ep, aEX_Nu, aEX_Q2, aEX_Xbj;
	TString aEX_Current ="", aEX_Cuts = "";
	Double_t aAngle=0.0, aE0 = 0.0, aP0 = 0.0;
	Int_t aRunNo= 0 ;

	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		aRunNo=aRunNoChain[i];

		TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),aRunNo);
		if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
		{
			aE0 = aE0Chain[i];
			aP0 = aP0Chain[i];
			aAngle = aAngleChain[i];	
			if(aArm=="R")
				aAngle = - fabs(aAngle); //Force to be negative on HRS-R, sometimes I forgot the sign
			if(aArm=="L")
				aAngle =   fabs(aAngle); //Force to be positive on HRS-L, sometimes I forgot the sign

			cerr<<Form("  === Run# = %d, E0 = %7.5f, P0 = %7.5f, Angle = %7.5f", aRunNo, aE0, aP0, aAngle)<<endl;
			outlog<<Form("  === Run# = %d, E0 = %7.5f, P0 = %7.5f, Angle = %7.5f", aRunNo, aE0, aP0, aAngle)<<endl;

			/////////////////////////////////////////////
			// Define Acceptance Variables
			////////////////////////////////////////////
			aSinA = sin(aAngle*TMath::DegToRad());
			aCosA = cos(aAngle*TMath::DegToRad());
			aSinSQ = pow(sin(aAngle*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2

			/*Define Acceptance Variables{{{*/
			if ( aIsExtTgt ){
				aEX_Dp= Form("ExTgt%s.dp",aArm.Data());
                //Apply Delta Correction in the data -- Z. Ye 02/18/2013
                if(aArm=="R")
	             aEX_Dp = gGet_Delta_Correct(aArm.Data(), aIsExtTgt);

				aEX_Th = Form("ExTgt%s.th",aArm.Data());
                if(aArm=="L")
					aEX_Ph = Form("(ExTgt%s.ph-0.002)",aArm.Data());//Add a Phi_tg offset here, ZYe 08/07/2014
				else 
					aEX_Ph = Form("(ExTgt%s.ph+0.002)",aArm.Data());//Add a Phi_tg offset here, ZYe 08/07/2014
				aEX_SA = Form("acos((%f-%s*%f)/sqrt(1.0+pow(%s,2)+pow(%s,2)))",aCosA,aEX_Ph.Data(),aSinA,aEX_Th.Data(),aEX_Ph.Data());
				aEX_SinSQ = Form("pow(sin(%s/2.0),2)",aEX_SA.Data()); //Define Sin(Theta/2.)^2

				aEX_Ep =Form("(%7.5f*(1.0+%s))",aP0,aEX_Dp.Data());
				aEX_Nu= Form("(%7.5f-%7.5f*(1.0+%s))",aE0,aP0,aEX_Dp.Data());
				aEX_Q2 = Form("4.0*%7.5f*(%7.5f*(1.0+%s))*%s",aE0,aP0,aEX_Dp.Data(),aEX_SinSQ.Data());
				aEX_Xbj=Form("%s/(2.0*%8.6f*%s)",aEX_Q2.Data(),PROTON_MASS,aEX_Nu.Data());
				//aEX_Xbj=Form("ExtEK%s.x_bj",aArm.Data());
			}
			else{
				aEX_Dp =Form("%s.gold.dp",aArm.Data());
                //Apply Delta Correction in the data -- Z. Ye 02/18/2013
				if(aArm=="R")
	             aEX_Dp = gGet_Delta_Correct(aArm,  aIsExtTgt);

				aEX_Th = Form("%s.gold.th",aArm.Data());
				if(aArm=="L")
					aEX_Ph = Form("(%s.gold.ph-0.002)",aArm.Data());//Add a Phi_tg offset here, ZYe 08/07/2014
				else
					aEX_Ph = Form("(%s.gold.ph+0.002)",aArm.Data());//Add a Phi_tg offset here, ZYe 08/07/2014
				aEX_SA = Form("acos((%f-%s*%f)/(sqrt(1.0+pow(%s,2)+pow(%s,2))))",aCosA,aEX_Ph.Data(),aSinA,aEX_Th.Data(),aEX_Ph.Data());
				aEX_SinSQ = Form("pow(sin(%s/2.0),2)",aEX_SA.Data()); //Define Sin(Theta/2.)^2

				aEX_Ep =Form("(%7.5f*(1.0+%s))",aP0,aEX_Dp.Data());
				aEX_Nu= Form("(%7.5f-%7.5f*(1.0+%s))",aE0,aP0,aEX_Dp.Data());
				aEX_Q2 = Form("4.0*%7.5f*(%7.5f*(1.0+%s))*%s",aE0,aP0,aEX_Dp.Data(),aEX_SinSQ.Data());
				aEX_Xbj=Form("%s/(2.0*%8.6f*%s)",aEX_Q2.Data(),PROTON_MASS,aEX_Nu.Data());
				//aEX_Xbj=Form("EK%s.x_bj",aArm.Data());
			}

			/*			cerr<<"============================================="<<endl;
						cerr<<Form("EX_Th = %s,  EX_Ph = %s, EX_SA = %s, EX_SinSQ = %s", aEX_Th.Data(),aEX_Ph.Data(), aEX_SA.Data(), aEX_SinSQ.Data() )<<endl;
					cerr<<"============================================="<<endl;
					cerr<<"============================================="<<endl;
					cerr<<Form("EX_Dp = %s", aEX_Dp.Data() )<<endl;
					cerr<<"============================================="<<endl;

					outlog<<"============================================="<<endl;
					outlog<<Form("EX_Th = %s,  EX_Ph = %s, EX_SA = %s, EX_SinSQ = %s", aEX_Th.Data(),aEX_Ph.Data(), aEX_SA.Data(), aEX_SinSQ.Data() )<<endl;
					outlog<<"============================================="<<endl;
         */
			//Current Branch for cutting out Beam Trip events
			// must run "AddCurrent" first to add those branches in rootfiles
			if(aArm=="R"||aArm=="L")
				aEX_Current = "right_current"; //LEFT BCM had problem for x>2 experiment
			else
				aEX_Current = "left_current";

			if(aBin_Variable=="Ep"){//Add a xbj<2.80 cut for He3 to exclude the elastic peak
				aEX_Cuts=Form("abs(%s-%f)<=%f&&%s>=%f",
						aEX_Ep.Data(),
						aBin,
						aBin_Size,
						aEX_Current.Data(),
						Beam_Trip_Cut
						);
		//	  if(aTarget=="He3")
		//		  aEX_Cuts+=Form("&&%s=<%f",aEX_Xbj.Data(), 2.8);
			}
			else if(aBin_Variable=="xbj"){
				double aEp_Min =1000.0, aEp_Max = -10000.0;
				gGet_Ep_Cut(aKin.Data(), aArm.Data(), &aEp_Min, &aEp_Max);
				aEX_Cuts=Form("%s>=%f&&%s<%f&&abs(%s-%f)<=%f&&%s>=%f",
						aEX_Ep.Data(),
						aEp_Min,
						aEX_Ep.Data(),
						aEp_Max,
						aEX_Xbj.Data(),
						aBin,
						aBin_Size,
						aEX_Current.Data(),
						Beam_Trip_Cut
						);
			}
			else
				cerr<<"*** ERROR***, in XGT2_Data.h, Ln166, unknown binning vairable!!"<<endl;
			aEX_Cuts+="&&"+aElectronCuts;
		    /*
			cerr<<"============================================="<<endl;
			cerr<<"    Cuts = "<< aEX_Cuts.Data()<<endl;
			cerr<<"============================================="<<endl;
			outlog<<"============================================="<<endl;
			outlog<<"    Cuts = "<< aEX_Cuts.Data()<<endl;
			outlog<<"============================================="<<endl;
		    */
			/*}}}*/

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
//Swape Stat Err and Sys Err, Z. Ye 06/12/2014

			for(int j=0;j<8;j++)
				aPS[j] = aPSChain[i*8+j];
			
			TChain *aT_Tree = gAddTree(aRunNo,"T");

			aNf[i] = aPS[aMain_Trigger-1]*aT_Tree->GetEntries(aEX_Cuts.Data());
	        //Ignore those runs that only contribute events less or equal to 0 counts otherwise the stat_err will blow up,2012/10/08
			if(aNf[i]>1e-16&&aLiveTimeChain[i]>1e-16){
                aNf_EX_Chain[i] = aNf[i]/aLiveTimeChain[i];
				//aNf_Stat_Err[i]= aNf_EX_Chain[i] * sqrt(1.0/aNf[i] + pow(aLiveTimeStatErrChain[i]/aLiveTimeChain[i] ,2));
				aNf_Stat_Err[i]=aNf_EX_Chain[i] * sqrt(1.0/aNf[i]);//LiveTime Stat is always zero, ZYe 06/12/2014
				aNf_Sys_Err[i] =aNf_EX_Chain[i] * aLiveTimeSysErrChain[i]/aLiveTimeChain[i];

				aNf_Total+= aNf_EX_Chain[i];
				aNf_Stat_Err_Total += pow(aLiveTimeStatErrChain[i]/aLiveTimeChain[i],2);//Zero basically
				aNf_Stat_Err_Total += pow(aNf_Stat_Err[i] /aNf_EX_Chain[i],2);
				aNf_Sys_Err_Total  += pow(aNf_Sys_Err[i]  /aNf_EX_Chain[i],2);

				aNf_EX_Stat[i] = aNf_Stat_Err[i];
				aNf_EX_Sys[i] = aNf_Sys_Err[i];
				outlog << Form("                 PS%d=%d, LiveTime/Err=%4.3f / %4.3f, Events/Err=%7.1f / %7.1f",
						aMain_Trigger,aPS[aMain_Trigger-1],aLiveTimeChain[i], aLiveTimeStatErrChain[i],aNf_EX_Chain[i],aNf_EX_Stat[i])<<endl;
				cerr   << Form("                 PS%d=%d, LiveTime/Err=%4.3f / %4.3f, Events/Err=%7.1f / %7.1f",
						aMain_Trigger,aPS[aMain_Trigger-1],aLiveTimeChain[i], aLiveTimeStatErrChain[i],aNf_EX_Chain[i],aNf_EX_Stat[i])<<endl;
			}
			delete aT_Tree;
		}
		else{
			outlog << Form("   **** Cann't find Rootfile %s in %s", File_Form.Data(),ROOTFILES_DIR.Data())<<endl;
			cerr   << Form("   **** Cann't find Rootfile %s in %s", File_Form.Data(),ROOTFILES_DIR.Data())<<endl;
		}
	}   
	nf->Value    = aNf_Total;
	nf->Stat_Err = aNf_Total * sqrt(1.0/aNf_Total + aNf_Stat_Err_Total);//FIX_HERE, not sure now ZYe 06/12/2014
	nf->Sys_Err  = aNf_Total * sqrt(aNf_Sys_Err_Total);
    return nf;
}
/*}}}*/

/*inline Double_t gGet_AlDummy_Factor(const TString& aKin){{{*/
inline Double_t gGet_AlDummy_Factor(const TString& aKin){

	const double Alumimium_Density = 2.708; //g/cm3
	const double Al_Thickness = 0.02998; //cm cell-down:0.02718, -up:0.03277, take the average
	const double Dummy_Thickness = 0.42800; //cm dummy-down:0.1577, -up:0.1584, take the average

    double RadCor_Ratio = 0.0; //XS_rad(Al)/XS_rad(Dummy)
	cerr<<"      Get Al-Dummy Correction Factor for "<<aKin.Data()<<" ..."<<endl;
	if(aKin.Contains("Kin3.1"))
		RadCor_Ratio = 1.1148;
	else if(aKin.Contains("Kin3.2"))
		RadCor_Ratio = 1.1239;
	else if(aKin.Contains("Kin4.1"))
		RadCor_Ratio = 1.1128;
	else if(aKin.Contains("Kin4.2"))
		RadCor_Ratio = 1.1249;
	else if(aKin.Contains("Kin5.0")&&!(aKin.Contains("Kin5.05")))
		RadCor_Ratio = 1.0741;
	else if(aKin.Contains("Kin5.05"))
		RadCor_Ratio = 1.0993;
	else if(aKin.Contains("Kin5.1"))
		RadCor_Ratio = 1.1096;
	else if(aKin.Contains("Kin5.2"))
		RadCor_Ratio = 1.1240;
	else if(aKin.Contains("Kin6.5"))
		RadCor_Ratio = 1.1381;
    else{
        cerr<<" **** ERROR, I can not find the Kin in gGet_AlDummy_Factor () in XGT2_Physics.h "<<endl;
		RadCor_Ratio = 1.115;
	}

	double aFactor = Al_Thickness / Dummy_Thickness * RadCor_Ratio;

	return aFactor;
}
/*}}}*/

/*XGT2_VAR* gGet_Nf_AL(const TString& aKin,const TString& aArm,const Double_t aBin,const Double_t aBin_Size,const TString& aElectronCuts,const bool aIsExtTgt){{{*/
XGT2_VAR* gGet_Nf_AL(const TString& aKin, const TString& aArm, const Double_t aBin, const Double_t aBin_Size, const TString& aElectronCuts, const TString& aBin_Variable,const bool aIsExtTgt){

	/*SetUp Run List{{{*/
	TString aDumInputFileName = Form("%s/xs_input/%sHRS/dum_%s.dat", XGT2_DIR.Data(), aArm.Data(), aKin.Data());
	ifstream aDumList(aDumInputFileName.Data());
	Int_t aDumRunNo = 0;
	vector<Int_t> aRunNoChain;
    cerr<<"      --- Dummy RunNoChain = ";    outlog<<"      --- Dummy RunNoChain = ";
	while(!(aDumList.eof())){
		aDumList >> aDumRunNo;
	    aRunNoChain.push_back(aDumRunNo);
		cerr<<aDumRunNo<<" ";     outlog<<aDumRunNo<<" ";
	}
	cerr<<endl; 	outlog<<endl;
    const int aChain_Size = aRunNoChain.size();

	Double_t* aLiveTimeChain = new double[aChain_Size];
	Double_t* aLiveTimeStatErrChain = new double[aChain_Size];
	Double_t* aLiveTimeSysErrChain = new double[aChain_Size];
	gLoad_LiveTime_Chain(aRunNoChain,aArm,aLiveTimeChain,aLiveTimeStatErrChain,aLiveTimeSysErrChain,2);//1->Beam Trip Cut, 2->NONE Cuts

	Int_t* aPSChain = new int[aChain_Size*8];
	gCheck_PS(aRunNoChain,aPSChain);

	Double_t* aE0Chain = new double[aChain_Size];//GeV
	Double_t* adEEChain = new double[aChain_Size];//GeV
	gCheck_E0(aRunNoChain,aArm,"E",aE0Chain, adEEChain); //Check P0 for each run from the Dipole readings. 

	Double_t* aP0Chain = new double[aChain_Size];//GeV/c
	gCheck_P0(aRunNoChain,aArm,"E",aP0Chain); //Check P0 for each run from the Dipole readings. 

	Double_t* aAngleChain = new double[aChain_Size];//Degree
	gCheck_Angle(aRunNoChain,aArm,"E",aAngleChain); //Check P0 for each run from the Dipole readings. 
	/*}}}*/

	XGT2_VAR* nf = new XGT2_VAR();

	Int_t aPS[8],aNf[aRunNoChain.size()];
	//Double_t aNf_Stat_Err[aRunNoChain.size()], aNf_Sys_Err[aRunNoChain.size()];
	Int_t aMain_Trigger; 
	Double_t aNf_Stat_Err_Total, aNf_Sys_Err_Total,aNf_Total;
	TString aS_Tree_Name;
	aNf_Total = 0.; aNf_Stat_Err_Total =0.0; aNf_Sys_Err_Total=0.0;

	Double_t aSinA =0.0, aCosA = 0.0, aSinSQ = 0.0;
	TString aEX_Dp, aEX_Th, aEX_Ph, aEX_SA, aEX_SinSQ, aEX_Ep, aEX_Nu, aEX_Q2, aEX_Xbj;
	TString aEX_Current ="", aEX_Cuts = "";
	Double_t aAngle=0.0, aE0 = 0.0, aP0 = 0.0;
	Int_t aRunNo= 0 ;

	Double_t* aNf_AL_Chain = new double[aChain_Size];
	Double_t* aNf_AL_StatChain = new double[aChain_Size];
	Double_t* aNf_AL_SysChain = new double[aChain_Size];

	//Consider the difference of the Al-Windows and Dummy-Foils thickness and radiative effect
	const double aFact_AlDummy = gGet_AlDummy_Factor(aKin.Data());

	for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
	{
		aRunNo=aRunNoChain[i];

		TString File_Form=Form("%s_%d.root",ROOTFILES_NAME.Data(),aRunNo);
		if ( gSystem->FindFile(ROOTFILES_DIR.Data(),File_Form) )
		{
			aE0 = aE0Chain[i];
			aP0 = aP0Chain[i];
			aAngle = aAngleChain[i];	
			if(aArm=="R")
				aAngle = - fabs(aAngle); //Force to be negative on HRS-R, sometimes I forgot the sign
			if(aArm=="L")
				aAngle =   fabs(aAngle); //Force to be positive on HRS-L, sometimes I forgot the sign

			cerr<<Form("  === Dummy Run# = %d, E0 = %7.5f, P0 = %7.5f, Angle = %7.5f", aRunNo, aE0, aP0, aAngle)<<endl;
			outlog<<Form("  === Dummy Run# = %d, E0 = %7.5f, P0 = %7.5f, Angle = %7.5f", aRunNo, aE0, aP0, aAngle)<<endl;

			/////////////////////////////////////////////
			// Define Acceptance Variables
			////////////////////////////////////////////
			aSinA = sin(aAngle*TMath::DegToRad());
			aCosA = cos(aAngle*TMath::DegToRad());
			aSinSQ = pow(sin(aAngle*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2

			/*Define Acceptance Variables{{{*/
			if ( aIsExtTgt ){
				aEX_Dp= Form("ExTgt%s.dp",aArm.Data());
                //Apply Delta Correction in the data -- Z. Ye 02/18/2013
                if(aArm=="R")
	             aEX_Dp = gGet_Delta_Correct(aArm.Data(), aIsExtTgt);

				aEX_Th = Form("ExTgt%s.th",aArm.Data());
                if(aArm=="L")
					aEX_Ph = Form("(ExTgt%s.ph-0.002)",aArm.Data());//Add an Phi_TG_Offset here, for L_HRS only for now
				else
					aEX_Ph = Form("(ExTgt%s.ph+0.002)",aArm.Data());//Add an Phi_TG_Offset here, for L_HRS only for now
				aEX_SA = Form("acos((%f-%s*%f)/sqrt(1.0+pow(%s,2)+pow(%s,2)))",aCosA,aEX_Ph.Data(),aSinA,aEX_Th.Data(),aEX_Ph.Data());
				aEX_SinSQ = Form("pow(sin(%s/2.0),2)",aEX_SA.Data()); //Define Sin(Theta/2.)^2

				aEX_Ep =Form("(%7.5f*(1.0+%s))",aP0,aEX_Dp.Data());
				aEX_Nu= Form("(%7.5f-%7.5f*(1.0+%s))",aE0,aP0,aEX_Dp.Data());
				aEX_Q2 = Form("4.0*%7.5f*(%7.5f*(1.0+%s))*%s",aE0,aP0,aEX_Dp.Data(),aEX_SinSQ.Data());
				aEX_Xbj=Form("%s/(2.0*%8.6f*%s)",aEX_Q2.Data(),PROTON_MASS,aEX_Nu.Data());
				//aEX_Xbj=Form("ExtEK%s.x_bj",aArm.Data());
			}
			else{
				aEX_Dp =Form("%s.gold.dp",aArm.Data());
                //Apply Delta Correction in the data -- Z. Ye 02/18/2013
				if(aArm=="R")
	             aEX_Dp = gGet_Delta_Correct(aArm,  aIsExtTgt);

				aEX_Th = Form("%s.gold.th",aArm.Data());
				if(aArm=="L")
					aEX_Ph = Form("(%s.gold.ph-0.002)",aArm.Data());//Add an Phi_TG_Offset here, for L_HRS only for now
				else
					aEX_Ph = Form("(%s.gold.ph+0.002)",aArm.Data());//Add an Phi_TG_Offset here, for L_HRS only for now
				aEX_SA = Form("acos((%f-%s*%f)/(sqrt(1.0+pow(%s,2)+pow(%s,2))))",aCosA,aEX_Ph.Data(),aSinA,aEX_Th.Data(),aEX_Ph.Data());
				aEX_SinSQ = Form("pow(sin(%s/2.0),2)",aEX_SA.Data()); //Define Sin(Theta/2.)^2

				aEX_Ep =Form("(%7.5f*(1.0+%s))",aP0,aEX_Dp.Data());
				aEX_Nu= Form("(%7.5f-%7.5f*(1.0+%s))",aE0,aP0,aEX_Dp.Data());
				aEX_Q2 = Form("4.0*%7.5f*(%7.5f*(1.0+%s))*%s",aE0,aP0,aEX_Dp.Data(),aEX_SinSQ.Data());
				aEX_Xbj=Form("%s/(2.0*%8.6f*%s)",aEX_Q2.Data(),PROTON_MASS,aEX_Nu.Data());
				//aEX_Xbj=Form("EK%s.x_bj",aArm.Data());
			}
			//Current Branch for cutting out Beam Trip events
			// must run "AddCurrent" first to add those branches in rootfiles
			if(aArm=="R"||aArm=="L")
				aEX_Current = "right_current"; //LEFT BCM had problem for x>2 experiment
			else
				aEX_Current = "left_current";

			if(aBin_Variable=="Ep")
				aEX_Cuts=Form("abs(%s-%f)<=%f&&%s>=%f",
						aEX_Ep.Data(),
						aBin,
						aBin_Size,
						aEX_Current.Data(),
						Beam_Trip_Cut
						);

			else if(aBin_Variable=="xbj")
				aEX_Cuts=Form("abs(%s-%f)<=%f&&%s>=%f",
						aEX_Xbj.Data(),
						aBin,
						aBin_Size,
						aEX_Current.Data(),
						Beam_Trip_Cut
						);
			else
				cerr<<"*** ERROR***, in XGT2_Data.h, Ln166, unknown binning vairable!!"<<endl;
			aEX_Cuts+="&&"+aElectronCuts;
		    /*
			cerr<<"============================================="<<endl;
			cerr<<"    Cuts = "<< aEX_Cuts.Data()<<endl;
			cerr<<"============================================="<<endl;
			outlog<<"============================================="<<endl;
			outlog<<"    Cuts = "<< aEX_Cuts.Data()<<endl;
			outlog<<"============================================="<<endl;
		    */
			/*}}}*/

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

			aNf[i] = aFact_AlDummy * aPS[aMain_Trigger-1]*aT_Tree->GetEntries(aEX_Cuts.Data());
	        //Ignore those runs that only contribute events less or equal to 0 counts otherwise the stat_err will blow up,2012/10/08
			if(aNf[i]>1e-16&&aLiveTimeChain[i]>1e-16){
                aNf_AL_Chain[i] = aNf[i]/aLiveTimeChain[i];
				aNf_AL_StatChain[i]= aNf_AL_Chain[i] * sqrt(1.0/aNf[i] + pow(aLiveTimeStatErrChain[i]/aLiveTimeChain[i] ,2));
				aNf_AL_SysChain[i]=aLiveTimeSysErrChain[i];

				aNf_Total+= aNf_AL_Chain[i];
				aNf_Stat_Err_Total += pow(aLiveTimeStatErrChain[i]/aLiveTimeChain[i],2);
				aNf_Sys_Err_Total  += pow(aNf_AL_SysChain[i] /aNf_AL_Chain[i],2);

				outlog << Form("                       PS%d=%d, LiveTime/Err=%4.3f / %4.3f, Events/Err=%7f / %7.1f",
						aMain_Trigger,aPS[aMain_Trigger-1],aLiveTimeChain[i], aLiveTimeStatErrChain[i],aNf_AL_Chain[i],aNf_AL_StatChain[i])<<endl;
				cerr   << Form("                       PS%d=%d, LiveTime/Err=%4.3f / %4.3f, Events/Err=%7f / %7.1f",
						aMain_Trigger,aPS[aMain_Trigger-1],aLiveTimeChain[i], aLiveTimeStatErrChain[i],aNf_AL_Chain[i],aNf_AL_StatChain[i])<<endl;
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
   
	delete aLiveTimeChain;
	delete aLiveTimeStatErrChain;
	delete aLiveTimeSysErrChain;
	delete aPSChain;
	delete aE0Chain;
	delete adEEChain;
	delete aP0Chain;
	delete aAngleChain;
    delete aNf_AL_Chain;
    delete aNf_AL_StatChain;
    delete aNf_AL_SysChain;

	return nf;
}
/*}}}*/

