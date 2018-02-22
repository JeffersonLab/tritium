#define NOCASE 1

//[3]-> BCM-Linear, BCMx3, BCMx10
//[2]-> Left, Right
const double d_Coeff[2] = { 3.264, 3.253 }; //, * 1e-4
const double u_Coeff[2] = { 3.954, 3.942 };//, * 1e-4
const double d_Offset[2] = {0.1055, 0.04668 };
const double u_Offset[2] = {0.1375, 0.07459 };

/*Double_t gGet_Rate_Avg(aRunNo,const TString& aArm){{{*/
Double_t gGet_Rate_Avg(Int_t aRunNo,const TString& aArm){
    //Note in use now
    Double_t aRate=-1000.0;

    TString aScaler_Tree_Name = LEFT_SCALER_TREE;
    TString aTRIGGER_RATE = Form("%s%d_r", LEFT_TRIGGER_SCALER.Data(), MAIN_TRIGGER_LEFT);
    if(aArm=="R"){
        aScaler_Tree_Name = RIGHT_SCALER_TREE;
        aTRIGGER_RATE = Form("%s%d_r", RIGHT_TRIGGER_SCALER.Data(), MAIN_TRIGGER_RIGHT);
    }
    TChain *aScaler_Tree = gAddTree(aRunNo,aScaler_Tree_Name.Data()); 

    if(aArm =="L"){
        Double_t offset_d = d_Offset[0];
        Double_t alpha_d = d_Coeff[0] * 1e-4;
        
        Double_t rate = 0.0, rms=0.0;
        Double_t I_max_d = aScaler_Tree->GetMaximum(LEFT_BCM_DR.Data()) * alpha_d;
        gGet_Leaf_Mean_and_RMS(aScaler_Tree, aTRIGGER_RATE.Data(), Form("%s*%f>%f", LEFT_BCM_DR.Data(), alpha_d, 0.5*I_max_d), &rate, &rms);
        
       aRate = rate;
    }
    else if(aArm == "R"){ 
        Double_t offset_d = d_Offset[1];
        Double_t alpha_d = d_Coeff[1] * 1e-4;

        Double_t rate = 0.0, rms=0.0;
        Double_t I_max_d = aScaler_Tree->GetMaximum(RIGHT_BCM_DR.Data()) * alpha_d;
        gGet_Leaf_Mean_and_RMS(aScaler_Tree, aTRIGGER_RATE.Data(), Form("%s*%f>%f", RIGHT_BCM_DR.Data(), alpha_d, 0.5*I_max_d), &rate, &rms);
       
        aRate = rate;
    }
    else
        aRate= -1000.0;

    delete aScaler_Tree;
    return aRate;
}
/*}}}*/

/*Double_t gGet_Current_Avg(aRunNo,const TString& aArm){{{*/
Double_t gGet_Current_Avg(Int_t aRunNo,const TString& aArm){
    //Note in use now
    Double_t aCurrent=-1000.0;

    TString aScaler_Tree_Name = LEFT_SCALER_TREE;
    if(aArm=="R")
        aScaler_Tree_Name = RIGHT_SCALER_TREE;
    TChain *aScaler_Tree = gAddTree(aRunNo,aScaler_Tree_Name.Data()); 

    if(aArm =="L"){
        Double_t offset_u = u_Offset[0];
        Double_t offset_d = d_Offset[0];
        Double_t alpha_u = u_Coeff[0] * 1e-4;
        Double_t alpha_d = d_Coeff[0] * 1e-4;
        
        Double_t bcm_ur = 0.0, bcm_dr=0.0, rms=0.0;
        Double_t I_max_u = aScaler_Tree->GetMaximum(LEFT_BCM_UR.Data()) * alpha_u;
        Double_t I_max_d = aScaler_Tree->GetMaximum(LEFT_BCM_DR.Data()) * alpha_d;
        gGet_Leaf_Mean_and_RMS(aScaler_Tree, LEFT_BCM_UR.Data(), Form("%s*%f>%f", LEFT_BCM_UR.Data(), alpha_u, 0.5*I_max_u), &bcm_ur, &rms);
        gGet_Leaf_Mean_and_RMS(aScaler_Tree, LEFT_BCM_DR.Data(), Form("%s*%f>%f", LEFT_BCM_DR.Data(), alpha_d, 0.5*I_max_d), &bcm_dr, &rms);
        
        //gGet_Leaf_Mean_and_RMS(aScaler_Tree, LEFT_BCM_UR.Data(), Form("%s>%f", LEFT_BCM_UR.Data(), BCM_TRIP_CUT), &bcm_ur, &rms);
        //gGet_Leaf_Mean_and_RMS(aScaler_Tree, LEFT_BCM_DR.Data(), Form("%s>%f", LEFT_BCM_DR.Data(), BCM_TRIP_CUT), &bcm_dr, &rms);

        Double_t i_u = bcm_ur * alpha_u;
        Double_t i_d = bcm_dr * alpha_d;

        aCurrent = (i_d+i_d)/2.0;
    }
    else if(aArm == "R"){ //Don't use Left Arm bcm scaler
        Double_t offset_u = u_Offset[1] * 1e-4;
        Double_t offset_d = d_Offset[1] * 1e-4;
        Double_t alpha_u = u_Coeff[1];
        Double_t alpha_d = d_Coeff[1];

        Double_t bcm_ur = 0.0, bcm_dr=0.0, rms=0.0;
        Double_t I_max_u = aScaler_Tree->GetMaximum(RIGHT_BCM_UR.Data()) * alpha_u;
        Double_t I_max_d = aScaler_Tree->GetMaximum(RIGHT_BCM_DR.Data()) * alpha_d;
        gGet_Leaf_Mean_and_RMS(aScaler_Tree, RIGHT_BCM_UR.Data(), Form("%s*%f>%f", RIGHT_BCM_UR.Data(), alpha_u, 0.5*I_max_u), &bcm_ur, &rms);
        gGet_Leaf_Mean_and_RMS(aScaler_Tree, RIGHT_BCM_DR.Data(), Form("%s*%f>%f", RIGHT_BCM_DR.Data(), alpha_d, 0.5*I_max_d), &bcm_dr, &rms);
        
        //gGet_Leaf_Mean_and_RMS(aScaler_Tree, RIGHT_BCM_UR.Data(), Form("%s>%f", RIGHT_BCM_UR.Data(), BCM_TRIP_CUT), &bcm_ur, &rms);
        //gGet_Leaf_Mean_and_RMS(aScaler_Tree, RIGHT_BCM_DR.Data(), Form("%s>%f", RIGHT_BCM_DR.Data(), BCM_TRIP_CUT), &bcm_dr, &rms);
       
        Double_t i_u = bcm_ur * alpha_u;
        Double_t i_d = bcm_dr * alpha_d;

        aCurrent = (i_d+i_d)/2.0;
    }
    else
        aCurrent= -1000.0;

    delete aScaler_Tree;
    return aCurrent;
}
/*}}}*/

/*Double_t gGet_BeamCharge_Avg(TTree *aScaler_Tree,const TString& aArm){{{*/
Double_t gGet_BeamCharge_Avg(TTree *aScaler_Tree,const TString& aArm, Double_t *aRuntime){
    //Note in use now
    Double_t aCharge=-1000.0;
    if(aArm =="L"){
        Double_t offset_u = u_Offset[0];
        Double_t offset_d = d_Offset[0];
        Double_t alpha_u = u_Coeff[0] * 1e-4;
        Double_t alpha_d = d_Coeff[0] * 1e-4;

        Double_t bcm_uc = aScaler_Tree->GetMaximum(LEFT_BCM_U.Data()) - aScaler_Tree->GetMinimum(LEFT_BCM_U.Data());
        Double_t bcm_dc = aScaler_Tree->GetMaximum(LEFT_BCM_D.Data()) - aScaler_Tree->GetMinimum(LEFT_BCM_D.Data());
        Double_t time =  aScaler_Tree->GetMaximum(LEFT_CLOCK.Data())/Clock_Freq - aScaler_Tree->GetMinimum(LEFT_CLOCK.Data())/Clock_Freq;
        Double_t i_u = bcm_uc*alpha_u;//+offset_u*time;
        Double_t i_d = bcm_dc*alpha_d;//+offset_d*time;

        aCharge = (i_d+i_d)/2.0;//* (4.69e4/1.57e7); // either upstream or downstream scaler has problems.
        *aRuntime = time;
    }
    else if(aArm == "R"){ //Don't use Left Arm bcm scaler
        Double_t offset_u = u_Offset[1];
        Double_t offset_d = d_Offset[1];
        Double_t alpha_u = u_Coeff[1] * 1e-4;
        Double_t alpha_d = d_Coeff[1] * 1e-4;

        Double_t bcm_uc = aScaler_Tree->GetMaximum(RIGHT_BCM_U.Data()) - aScaler_Tree->GetMinimum(RIGHT_BCM_U.Data());
        Double_t bcm_dc = aScaler_Tree->GetMaximum(RIGHT_BCM_D.Data()) - aScaler_Tree->GetMinimum(RIGHT_BCM_D.Data());
        Double_t time =  aScaler_Tree->GetMaximum(RIGHT_CLOCK.Data())/Clock_Freq - aScaler_Tree->GetMinimum(RIGHT_CLOCK.Data())/Clock_Freq;
        Double_t i_u = bcm_uc*alpha_u;//+offset_u*time;
        Double_t i_d = bcm_dc*alpha_d;//+offset_d*time;

        *aRuntime = time;
        aCharge = (i_d+i_d)/2.0;//* (4.69e4/1.57e7); // either upstream or downstream scaler has problems.
    }
    else
        aCharge= -1000.0;

    return aCharge;
}
/*}}}*/

/*inline TRI_VAR* gCal_Ne(const  vector<Int_t>& aRunNoChain,const TString& aArm,Double_t& oRuntime){{{*/
inline TRI_VAR* gCal_Ne(const vector<Int_t>& aRunNoChain,const TString& aArm,Double_t* aRuntime, double* aNeChain, double* aCurrentChain, double* aRateChain)
{ 
    TString aScaler_Tree_Name = LEFT_SCALER_TREE;
    if(aArm=="R")
        aScaler_Tree_Name = RIGHT_SCALER_TREE;

  //return uC
  Double_t oCharge_Sum=0.0,oCharge_Stat_Sum=0.0,oCharge_Sys_Sum=0.0,oRuntime_Sum=0.0, oCharge=0.0,oRuntime=0.0;     
  for (unsigned int i=0; i<aRunNoChain.size(); i++ ){
      Int_t aRunNo=aRunNoChain[i];
	  TChain *aScaler_Tree = gAddTree(aRunNo,aScaler_Tree_Name.Data()); 
               
      oCharge= gGet_BeamCharge_Avg(aScaler_Tree, aArm.Data(), &oRuntime);
      oCharge_Sum += oCharge;
	  oCharge_Stat_Sum += 1.0/oCharge;
	  oCharge_Sys_Sum  += 0.0;
      oRuntime_Sum += oRuntime;

      aNeChain[i] = oCharge/(Qe*1e6);

      aCurrentChain[i] = gGet_Current_Avg(aRunNo, aArm.Data());
      aRateChain[i] = gGet_Rate_Avg(aRunNo, aArm.Data());
        
      cout<<Form("--- Charge = %e, time =%e mins, Rate=%5.3f,  I = %4.3f", oCharge,oRuntime/60., aRateChain[i], aCurrentChain[i])<<endl;

      delete aScaler_Tree;
    } // for ( Int_t i=0; i<aRunNoChain.size(); i++ )
  
  *aRuntime = oRuntime_Sum;
  
  TRI_VAR* output=new TRI_VAR();
  output->Value=oCharge_Sum/(Qe*1e6);
  output->Stat_Err=0;
  output->Sys_Err=(output->Value * sqrt(oCharge_Stat_Sum/2.0)); //A Factor of 1/2 is from the average of (u,d)
	
  return output;
}
/*}}}*/
