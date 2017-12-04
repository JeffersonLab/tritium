
/*TRI_VAR* gGet_Nf_EX(const vector<Int_t>& aRunNoChain,const TString& aArm, const TString& aKin, const Double_t aBin, const Double_t aBin_Size, const TString& aElectronCuts,const TString& aBin_Variable, const Double_t* aLiveTimeChain,const Double_t* aLiveTimeStatErrChain,const Double_t* aLiveTimeSysErrChain, const Int_t* aPSChain, const Double_t* aE0Chain, const Double_t* aP0Chain, const Double_t* aAngleChain, double* aNf_EX_Chain,double* aNf_EX_Stat,double* aNf_EX_Sys, const bool aIsExtTgt){{{*/
TRI_VAR* gGet_Nf_EX(const vector<Int_t>& aRunNoChain,const TString& aArm, const TString& aKin, const TString& aTarget,const int aMain_Trigger, const Double_t aBin, const Double_t aBin_Size, const TString& aElectronCuts, const double *aVZ_Cut, const TString aBin_Variable, const Double_t* aLiveTimeChain,const Double_t* aLiveTimeStatErrChain,const Double_t* aLiveTimeSysErrChain, const Int_t* aPSChain, const Double_t* aE0Chain, const Double_t* aP0Chain, const Double_t* aAngleChain, const bool aIsExtTgt, double* aNf_EX_Chain,double* aNf_EX_Stat,double* aNf_EX_Sys){
    TRI_VAR* nf = new TRI_VAR();

	Int_t aPS[8],aNf[aRunNoChain.size()];
	Double_t aNf_Stat_Err[aRunNoChain.size()], aNf_Sys_Err[aRunNoChain.size()];
	Double_t aNf_Stat_Err_Total, aNf_Sys_Err_Total,aNf_Total;
	aNf_Total = 0.; aNf_Stat_Err_Total =0.0; aNf_Sys_Err_Total=0.0;

	Double_t aSinA =0.0, aCosA = 0.0, aSinSQ = 0.0;
	TString aEX_Dp, aEX_Th, aEX_Ph, aEX_SA, aEX_SinSQ, aEX_Ep, aEX_Nu, aEX_Q2, aEX_Xbj;
	TString aEX_Current ="", aEX_Cuts = "";
	Double_t aAngle=0.0, aE0 = 0.0, aP0 = 0.0;
	Int_t aRunNo= 0 ;

    Double_t TARGET_MASS_A = PROTON_MASS;

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
				aEX_Dp= Form("ex%s.dp",aArm.Data());
				aEX_Th = Form("ex%s.th",aArm.Data());
                aEX_Ph = Form("(ex%s.ph)",aArm.Data());
                aEX_SA = Form("acos((%f-%s*%f)/sqrt(1.0+pow(%s,2)+pow(%s,2)))",aCosA,aEX_Ph.Data(),aSinA,aEX_Th.Data(),aEX_Ph.Data());
				aEX_SinSQ = Form("pow(sin(%s/2.0),2)",aEX_SA.Data()); //Define Sin(Theta/2.)^2

				aEX_Ep =Form("(%7.5f*(1.0+%s))",aP0,aEX_Dp.Data());
				aEX_Nu= Form("(%7.5f-%7.5f*(1.0+%s))",aE0,aP0,aEX_Dp.Data());
				aEX_Q2 = Form("4.0*%7.5f*(%7.5f*(1.0+%s))*%s",aE0,aP0,aEX_Dp.Data(),aEX_SinSQ.Data());
				aEX_Xbj=Form("%s/(2.0*%8.6f*%s)",aEX_Q2.Data(),TARGET_MASS_A,aEX_Nu.Data());
			}
			else{
				aEX_Dp =Form("%s.gold.dp",aArm.Data());
				aEX_Th = Form("%s.gold.th",aArm.Data());
				aEX_Ph = Form("(%s.gold.ph)",aArm.Data());
				aEX_SA = Form("acos((%f-%s*%f)/(sqrt(1.0+pow(%s,2)+pow(%s,2))))",aCosA,aEX_Ph.Data(),aSinA,aEX_Th.Data(),aEX_Ph.Data());
				aEX_SinSQ = Form("pow(sin(%s/2.0),2)",aEX_SA.Data()); //Define Sin(Theta/2.)^2

				aEX_Ep =Form("(%7.5f*(1.0+%s))",aP0,aEX_Dp.Data());
				aEX_Nu= Form("(%7.5f-%7.5f*(1.0+%s))",aE0,aP0,aEX_Dp.Data());
				aEX_Q2 = Form("4.0*%7.5f*(%7.5f*(1.0+%s))*%s",aE0,aP0,aEX_Dp.Data(),aEX_SinSQ.Data());
				aEX_Xbj=Form("%s/(2.0*%8.6f*%s)",aEX_Q2.Data(),TARGET_MASS_A,aEX_Nu.Data());
			}
            //TString aEX_VZ = Form("rp%s.z", aArm.Data());
            TString aEX_VZ = "";
            if(aArm=="R") aEX_VZ = "rpr.z";
            if(aArm=="L") aEX_VZ = "rpl.z";

            
			if(aBin_Variable=="Ep"){
				aEX_Cuts=Form("abs(%s-%f)<=%f&&(%s>=%f && %s<%f)",
						aEX_Ep.Data(),
						aBin,
						aBin_Size,
                        aEX_VZ.Data(),
                        aVZ_Cut[0],
                        aEX_VZ.Data(),
                        aVZ_Cut[1]
						);
			}else if(aBin_Variable=="xbj")
				aEX_Cuts=Form("abs(%s-%f)<=%f&&(%s>=%f && %s<%f)",
						aEX_Xbj.Data(),
						aBin,
						aBin_Size,
						aEX_VZ.Data(),
                        aVZ_Cut[0],
                        aEX_VZ.Data(),
                        aVZ_Cut[1]
                        );
            else
                cerr<<"*** ERROR***, in TRI_Data.h, Ln166, unknown binning vairable!!"<<endl;
			
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
	nf->Stat_Err = aNf_Total * sqrt(1.0/aNf_Total);
	nf->Sys_Err  = aNf_Total * sqrt(aNf_Sys_Err_Total);
    return nf;
}
/*}}}*/
