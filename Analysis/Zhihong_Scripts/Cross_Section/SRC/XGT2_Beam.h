#define NOCASE 1
//////////////////////////////////////
//  Beam Current and Charge
//////////////////////////////////////
/*inline gGet_Current(TTree* aTree, const TString& aArm){{{*/
inline double gGet_Current(TTree* aTree, const TString& aArm){
  //Note: Those branches are not in the rootfiles unless "AddCurrent" has been run for each rootfile to add those variables 
  Double_t aCurrent;
 
  if(aArm =="X"){ ////In X>2, we only use RIGHT scaler
    Double_t aIL_d1, aIL_d3,aIL_u1,aIL_u3,aIL_Avg,rms;    
    gGet_Leaf_Mean_and_RMS(aTree,"IL_u1",Form("IL_u1>=%f&&IL_u1<200.0",Beam_Trip_Cut),aIL_u1,rms); //Get Current from Upstream BCM x3
    gGet_Leaf_Mean_and_RMS(aTree,"IL_u3",Form("IL_u3>=%f&&IL_u3<200.0",Beam_Trip_Cut),aIL_u3,rms); //Get Current from Upstream BCM x3
    gGet_Leaf_Mean_and_RMS(aTree,"IL_d1",Form("IL_d1>=%f&&IL_d1<200.0",Beam_Trip_Cut),aIL_d1,rms); //Get Current from Downstream BCM x1
    gGet_Leaf_Mean_and_RMS(aTree,"IL_d3",Form("IL_d3>=%f&&IL_d3<200.0",Beam_Trip_Cut),aIL_d3,rms); //Get Current from Downstream BCM x3
    gGet_Leaf_Mean_and_RMS(aTree,"left_current",Form("left_current>=%f&&left_current<200.0", Beam_Trip_Cut),aIL_Avg,rms); //Get Current from average of all BCM 
    aCurrent = aIL_Avg; // either upstream or downstream scaler has problems.
  }
//  else if(aArm == "R")
  else if(aArm == "R"||aArm == "L"){ //In X>2, we only use RIGHT scaler
    Double_t aIR_d1, aIR_d3,aIR_u1,aIR_u3,aIR_Avg,rms;    
    gGet_Leaf_Mean_and_RMS(aTree,"IR_u1",Form("IR_u1>=%f&&IR_u1<200.0",Beam_Trip_Cut),aIR_u1,rms); //Get Current from Upstream BCM x3
    gGet_Leaf_Mean_and_RMS(aTree,"IR_u3",Form("IR_u3>=%f&&IR_u3<200.0",Beam_Trip_Cut),aIR_u3,rms); //Get Current from Upstream BCM x3
    gGet_Leaf_Mean_and_RMS(aTree,"IR_d1",Form("IR_d1>=%f&&IR_d1<200.0",Beam_Trip_Cut),aIR_d1,rms); //Get Current from Downstream BCM x1
    gGet_Leaf_Mean_and_RMS(aTree,"IR_d3",Form("IR_d3>=%f&&IR_d3<200.0",Beam_Trip_Cut),aIR_d3,rms); //Get Current from Downstream BCM x3
    gGet_Leaf_Mean_and_RMS(aTree,"right_current",Form("right_current>=%f&&right_current<200.0", Beam_Trip_Cut),aIR_Avg,rms); //Get Current from average of all BCM 
    aCurrent = aIR_Avg; 
  }
  else
    aCurrent = -1000.0;
  
  return aCurrent;
}
/*}}}*/

/*inline gCheck_Current(TTree* aTree, const TString& aCuts){{{*/
inline void gCheck_Current(const vector<Int_t>& aRunNoChain){
  double aCurrent_Left,aCurrent_Right,rms;
  
  for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
    {
      Int_t aRunNo=aRunNoChain[i];
      
      TChain *aTree = gAddTree(aRunNo, "T");
      gGet_Leaf_Mean_and_RMS(aTree,"left_current",  Form("left_current>=%f",Beam_Trip_Cut),aCurrent_Left,rms); //Get Current from Upstream BCM
      gGet_Leaf_Mean_and_RMS(aTree,"right_current", Form("right_current>=%f",Beam_Trip_Cut),aCurrent_Right,rms); //Get Current from Downstream BCM

      outlog<<Form("      For Run = %d, Current_Left = %6.2f, Current_Right = %6.2f",aRunNo, aCurrent_Left, aCurrent_Right)<<endl;
      cerr<<Form("      For Run = %d, Current_Left = %6.2f, Current_Right = %6.2f",aRunNo, aCurrent_Left, aCurrent_Right)<<endl;
      delete aTree;
    }
}
/*}}}*/

/*Double_t gGet_BeamCharge(TTree *aT_Tree,const TString& aArm){{{*/
Double_t gGet_BeamCharge(TTree *aT_Tree,const TString& aArm){
  //Note: Those branches are not in the rootfiles unless "AddCurrent" has been run for each rootfile to add those variables 
  Double_t charge=-1000.0;
  Int_t aEntries = aT_Tree->GetEntries(Form("right_current>=%f",Beam_Trip_Cut)); 
  Double_t aC_u1,aC_u3,aC_d1,aC_d3,rms;
  
  if(aArm =="L"||aArm=="R"){ //Only use BCM info from right scaler
      gGet_Leaf_Mean_and_RMS(aT_Tree,"C_u1",  Form("right_current>=%f",Beam_Trip_Cut),aC_u1,rms); //Get Charge from Upstream BCM
      gGet_Leaf_Mean_and_RMS(aT_Tree,"C_u3",  Form("right_current>=%f",Beam_Trip_Cut),aC_u3,rms); //Get Charge from Upstream BCM
      gGet_Leaf_Mean_and_RMS(aT_Tree,"C_d1",  Form("right_current>=%f",Beam_Trip_Cut),aC_d1,rms); //Get Charge from Upstream BCM
      gGet_Leaf_Mean_and_RMS(aT_Tree,"C_d3",  Form("right_current>=%f",Beam_Trip_Cut),aC_d3,rms); //Get Charge from Upstream BCM
      //charge = (aC_u1+aC_u3+aC_d1+aC_d3)/4.0*aEntries; //Sum = Mean*Entries    
      charge = (aC_u1+aC_d1)/2.0*aEntries; //Sum = Mean*Entries   Only use U1 and D1 since we run at very high current, 12/19/2013 
  }
  return charge;
}
/*}}}*/

/*Double_t gGet_BeamCharge_Avg(TTree *aT_Tree,const TString& aArm){{{*/
Double_t gGet_BeamCharge_Avg(TTree *aT_Tree,const TString& aArm){
    //Note in use now
    Double_t charge=-1000.0, time = -1000.00;
  
  if(aArm =="L")
    time =  aT_Tree->GetMaximum("left_clkcount")/1024 - aT_Tree->GetMinimum("left_clkcount")/1024;
  else if(aArm=="R")
    time =  aT_Tree->GetMaximum("right_clkcount")/1024 - aT_Tree->GetMinimum("right_clkcount")/1024;
  else
    time = 0.0;

  charge = time * gGet_Current_Avg(aT_Tree, aArm.Data());
  return charge;
}
/*}}}*/

/*inline gGet_Current_Avg(TTree* aTree, const TString& aArm){{{*/
inline double gGet_Current_Avg(TTree* aT_Tree, const TString& aArm){
  Double_t aCurrent;
  if(aArm =="X"){
    Double_t offset_u = u_Offset[0][0];
    Double_t offset_d = d_Offset[0][0];
    Double_t alpha_u = u_Coeff[0][0];
    Double_t alpha_d = d_Coeff[0][0];
    
    Double_t bcm_u1c = aT_Tree->GetMaximum("left_bcm_u1c") - aT_Tree->GetMinimum("left_bcm_u1c");
    Double_t bcm_d1c = aT_Tree->GetMaximum("left_bcm_d1c") - aT_Tree->GetMinimum("left_bcm_d1c");
    Double_t time =  aT_Tree->GetMaximum("left_clkcount")/1024 - aT_Tree->GetMinimum("left_clkcount")/1024;
    Double_t i_u1 = ((bcm_u1c/time)-offset_u)/alpha_u;
    Double_t i_d1 = ((bcm_d1c/time)-offset_d)/alpha_d;
    
    aCurrent = (i_d1+i_u1)/2.0; // either upstream or downstream scaler has problems.
  }
  else if(aArm == "R"||aArm=="L"){ //Don't use Left Arm bcm scaler
    Double_t offset_u = u_Offset[1][0];
    Double_t offset_d = d_Offset[1][0];
    Double_t alpha_u = u_Coeff[0][1];
    Double_t alpha_d = d_Coeff[0][1];
    
    Double_t bcm_u1c = aT_Tree->GetMaximum("right_bcm_u1c") - aT_Tree->GetMinimum("right_bcm_u1c");
    Double_t bcm_d1c = aT_Tree->GetMaximum("right_bcm_d1c") - aT_Tree->GetMinimum("right_bcm_d1c");
    Double_t  time =  aT_Tree->GetMaximum("right_clkcount")/1024 - aT_Tree->GetMinimum("right_clkcount")/1024;
    Double_t  i_u1 = ((bcm_u1c/time)-offset_u)/alpha_u;
    Double_t  i_d1 = ((bcm_d1c/time)-offset_d)/alpha_d; 

    aCurrent = (i_d1+i_u1)/2.0;
  }
  else
    aCurrent = -1000.0;
  
  return aCurrent;
}
/*}}}*/

/*inline XGT2_VAR* gCal_Ne(const  vector<Int_t>& aRunNoChain,const TString& aArm,Double_t& oRuntime){{{*/
inline XGT2_VAR* gCal_Ne(const vector<Int_t>& aRunNoChain,const TString& aArm,Double_t& oRuntime, double* aNeChain)
{
  //return uC
  Double_t oCharge_Sum=0.0,oCharge_Stat_Sum=0.0,oCharge_Sys_Sum=0.0,oRuntime_Sum=0.0, aCharge=0.0,aRuntime=0.0;     
  for (unsigned int i=0; i<aRunNoChain.size(); i++ ){
      Int_t aRunNo=aRunNoChain[i];
      //The other option to get the charge is to get from "T"  tree, which should be the same since I sync. them well.
      //Again, here we only use scaler info from Right, but not left
	  TChain *aT_Tree = gAddTree(aRunNo,"RIGHT"); 
               
      aT_Tree->SetBranchAddress("clkcount", &aRuntime);
      Int_t nentries=aT_Tree->GetEntries(Form("right_current>=%f",Beam_Trip_Cut));
      aT_Tree->GetEntry(nentries-1);
      aRuntime /= 1024.;
      oRuntime_Sum += aRuntime;
      
      aCharge= gGet_BeamCharge(aT_Tree, aArm.Data());
      oCharge_Sum += aCharge;
	  oCharge_Stat_Sum += 1.0/aCharge;
	  oCharge_Sys_Sum  += 0.0;

      aNeChain[i] = aCharge/(Qe*1e6);
	  outlog <<"   For Run#"<<aRunNo << ", Run Time = "<<aRuntime <<", Ne = "<<aCharge/(Qe*1e6)<<endl;
      cerr <<"   For Run#"<<aRunNo << ", Run Time = "<<aRuntime <<", Ne = "<<aCharge/(Qe*1e6)<<endl;
      delete aT_Tree;
    } // for ( Int_t i=0; i<aRunNoChain.size(); i++ )
  
  oRuntime = oRuntime_Sum;
  
  XGT2_VAR* output=new XGT2_VAR();
  output->Value=oCharge_Sum/(Qe*1e6);
  output->Stat_Err=0;
  output->Sys_Err=(output->Value * sqrt(oCharge_Stat_Sum/2.0)); //A Factor of 1/2 is from the average of (u1,d1)
	
  return output;
}
/*}}}*/
