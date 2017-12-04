#define NOCASE 1
//////////////////////////////////////
//  Beam Current and Charge
//////////////////////////////////////

/*inline gGet_Current_Avg(TTree* aTree, const TString& aArm){{{*/
inline double gGet_Current_Avg(TTree* aT_Tree, const TString& aArm){
  Double_t aCurrent;
  if(aArm =="L"){
    Double_t offset_u = u_Offset[0][0];
    Double_t offset_d = d_Offset[0][0];
    Double_t alpha_u = u_Coeff[0][0];
    Double_t alpha_d = d_Coeff[0][0];
    
    Double_t bcm_u1c = aT_Tree->GetMaximum(LEFT_BCM_U.Data()) - aT_Tree->GetMinimum(LEFT_BCM_U.Data());
    Double_t bcm_d1c = aT_Tree->GetMaximum(LEFT_BCM_D.Data()) - aT_Tree->GetMinimum(LEFT_BCM_D.Data());
    Double_t time =  aT_Tree->GetMaximum(LEFT_CLOCK.Data())/Clock_Freq - aT_Tree->GetMinimum(LEFT_CLOCK.Data())/Clock_Freq;
    Double_t i_u1 = ((bcm_u1c/time)-offset_u)/alpha_u;
    Double_t i_d1 = ((bcm_d1c/time)-offset_d)/alpha_d;
    
    aCurrent = (i_d1+i_u1)/2.0; // either upstream or downstream scaler has problems.
  }
  else if(aArm == "R"){ //Don't use Left Arm bcm scaler
    Double_t offset_u = u_Offset[1][0];
    Double_t offset_d = d_Offset[1][0];
    Double_t alpha_u = u_Coeff[0][1];
    Double_t alpha_d = d_Coeff[0][1];
    
    Double_t bcm_u1c = aT_Tree->GetMaximum(RIGHT_BCM_U.Data()) - aT_Tree->GetMinimum(RIGHT_BCM_U.Data());
    Double_t bcm_d1c = aT_Tree->GetMaximum(RIGHT_BCM_D.Data()) - aT_Tree->GetMinimum(RIGHT_BCM_D.Data());
    Double_t time =  aT_Tree->GetMaximum(RIGHT_CLOCK.Data())/Clock_Freq - aT_Tree->GetMinimum(RIGHT_CLOCK.Data())/Clock_Freq;
    Double_t i_u1 = ((bcm_u1c/time)-offset_u)/alpha_u;
    Double_t i_d1 = ((bcm_d1c/time)-offset_d)/alpha_d; 

    aCurrent = (i_d1+i_u1)/2.0;
  }
  else
    aCurrent = -1000.0;
  
  return aCurrent;
}
/*}}}*/

/*Double_t gGet_BeamCharge_Avg(TTree *aT_Tree,const TString& aArm){{{*/
Double_t gGet_BeamCharge_Avg(TTree *aT_Tree,const TString& aArm){
    //Note in use now
    Double_t charge=-1000.0, time = -1000.00;
  
  if(aArm =="L")
    time =  aT_Tree->GetMaximum(LEFT_CLOCK.Data())/Clock_Freq - aT_Tree->GetMinimum(LEFT_CLOCK.Data())/Clock_Freq;
  else if(aArm=="R")
    time =  aT_Tree->GetMaximum(RIGHT_CLOCK.Data())/Clock_Freq - aT_Tree->GetMinimum(RIGHT_CLOCK.Data())/Clock_Freq;
  else
    time = 0.0;

  charge = time * gGet_Current_Avg(aT_Tree, aArm.Data());
  return charge;
}
/*}}}*/

/*inline TRI_VAR* gCal_Ne(const  vector<Int_t>& aRunNoChain,const TString& aArm,Double_t& oRuntime){{{*/
inline TRI_VAR* gCal_Ne(const vector<Int_t>& aRunNoChain,const TString& aArm, const TString& aArm_Scaler,Double_t& oRuntime, double* aNeChain)
{
  //return uC
  Double_t oCharge_Sum=0.0,oCharge_Stat_Sum=0.0,oCharge_Sys_Sum=0.0,oRuntime_Sum=0.0, aCharge=0.0,aRuntime=0.0;     
  for (unsigned int i=0; i<aRunNoChain.size(); i++ ){
      Int_t aRunNo=aRunNoChain[i];
	  TChain *aT_Tree = gAddTree(aRunNo,aArm_Scaler.Data()); 
               
      TString aCLOCK="";
      if(aArm=="L") aCLOCK = LEFT_CLOCK;
      if(aArm=="R") aCLOCK = RIGHT_CLOCK;

      aT_Tree->SetBranchAddress(aCLOCK.Data(), &aRuntime);
      Int_t nentries=aT_Tree->GetEntries();
      aT_Tree->GetEntry(nentries-1);
      aRuntime /= Clock_Freq;
      oRuntime_Sum += aRuntime;
      
      aCharge= gGet_BeamCharge_Avg(aT_Tree, aArm.Data());
      oCharge_Sum += aCharge;
	  oCharge_Stat_Sum += 1.0/aCharge;
	  oCharge_Sys_Sum  += 0.0;

      aNeChain[i] = aCharge/(Qe*1e6);
	  outlog <<"   For Run#"<<aRunNo << ", Run Time = "<<aRuntime <<", Ne = "<<aCharge/(Qe*1e6)<<endl;
      cerr <<"   For Run#"<<aRunNo << ", Run Time = "<<aRuntime <<", Ne = "<<aCharge/(Qe*1e6)<<endl;
      delete aT_Tree;
    } // for ( Int_t i=0; i<aRunNoChain.size(); i++ )
  
  oRuntime = oRuntime_Sum;
  
  TRI_VAR* output=new TRI_VAR();
  output->Value=oCharge_Sum/(Qe*1e6);
  output->Stat_Err=0;
  output->Sys_Err=(output->Value * sqrt(oCharge_Stat_Sum/2.0)); //A Factor of 1/2 is from the average of (u1,d1)
	
  return output;
}
/*}}}*/
