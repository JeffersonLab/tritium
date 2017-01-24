#define NOCASE 1
//////////////////////////////////////////
// BCM Constant
//////////////////////////////////////////

/*BCM{{{*/
// Note: Only Left U1,D1 are calibrated, but D1 had problem during run# 3835~3888.
//       Right BCM looks like having problems.
const Double_t u_Coeff[2][3]=
  {
    //[2]
    //0=Left 1=Right
    //[3]
    //0=1, 1=3, 2=10 u1,u3,u10
//Alpha: U1    U3    U10         //Updated in 07/20/2012
    {2070.0,    0.0, 0.0},            //u Coeff to calculate Charge for left arm
    {1035.0, 3178.0, 0.0},            //u Coeff to calculate Charge for right arm
  };
//                 Offset        D1     D3    D10        //Updated in 07/20/2012
const Double_t u_Offset[2][3]={{400.0,   0.0, 0.0},    //the offsets of the u1,u3,u10 for left arm
			       {200.0, 174.0, 0.0},};  //the offsets of the u1,u3,u10 for right arm
const Double_t d_Coeff[2][3]=
  {
    //[2]
    //0=Left 1=Right
    //[3]
    //0=1, 1=3, 2=10 u1,u3,u10
   
//Alpha: D1    D3    D10        //Updated in 07/20/2012
    {1263.0,    0.0, 0.0},      //d Coeff to calculate Charge for left arm
    {1263.0, 3902.0, 0.0},      //d Coeff to calculate Charge for right arm
  };
//                 Offset        D1     D3  D10       //Updated in 07/20/2012
const Double_t d_Offset[2][3]={{20.0,  0.0, 0.0},     //the offsets of the d1,d3,d10 for left arm
			       {37.0, 22.0, 0.0},};   //the offsets of the d1,d3,d10 for right arm
/*}}}*/

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
   // gGet_Leaf_Mean_and_RMS(aTree,"IL_u3",Form("IL_u3>=%f&&IL_u3<200.0",Beam_Trip_Cut),aIL_u3,rms); //Get Current from Upstream BCM x3
    gGet_Leaf_Mean_and_RMS(aTree,"IL_d1",Form("IL_d1>=%f&&IL_d1<200.0",Beam_Trip_Cut),aIL_d1,rms); //Get Current from Downstream BCM x1
    //gGet_Leaf_Mean_and_RMS(aTree,"IL_d3",Form("IL_d3>=%f&&IL_d3<200.0",Beam_Trip_Cut),aIL_d3,rms); //Get Current from Downstream BCM x3
    gGet_Leaf_Mean_and_RMS(aTree,"left_current",Form("left_current>=%f&&left_current<200.0", Beam_Trip_Cut),aIL_Avg,rms); //Get Current from average of all BCM 
    aCurrent = aIL_Avg; // either upstream or downstream scaler has problems.
  }
//  else if(aArm == "R")
  else if(aArm == "R"||aArm == "L"){ //In X>2, we only use RIGHT scaler
    Double_t aIR_d1, aIR_d3,aIR_u1,aIR_u3,aIR_Avg,rms;    
    gGet_Leaf_Mean_and_RMS(aTree,"IR_u1",Form("IR_u1>=%f&&IR_u1<200.0",Beam_Trip_Cut),aIR_u1,rms); //Get Current from Upstream BCM x3
    //gGet_Leaf_Mean_and_RMS(aTree,"IR_u3",Form("IR_u3>=%f&&IR_u3<200.0",Beam_Trip_Cut),aIR_u3,rms); //Get Current from Upstream BCM x3
    gGet_Leaf_Mean_and_RMS(aTree,"IR_d1",Form("IR_d1>=%f&&IR_d1<200.0",Beam_Trip_Cut),aIR_d1,rms); //Get Current from Downstream BCM x1
    //gGet_Leaf_Mean_and_RMS(aTree,"IR_d3",Form("IR_d3>=%f&&IR_d3<200.0",Beam_Trip_Cut),aIR_d3,rms); //Get Current from Downstream BCM x3
  //  gGet_Leaf_Mean_and_RMS(aTree,"right_current",Form("right_current>=%f&&right_current<200.0", Beam_Trip_Cut),aIR_Avg,rms); //Get Current from average of all BCM 
    aCurrent = (aIR_u1+aIR_d1)/2.0; 
  }
  else
    aCurrent = -1000.0;
  
  return aCurrent;
}
/*}}}*/

/*Double_t gGet_BeamCharge(TTree *aT_Tree,const TString& aArm){{{*/
Double_t gGet_BeamCharge(TTree *aT_Tree,const TString& aArm){
  Double_t charge=-1000.0;
  Int_t aEntries = aT_Tree->GetEntries(Form("right_current>=%f",Beam_Trip_Cut)); 
  Double_t aC_u1,aC_u3,aC_d1,aC_d3,rms;
  
  if(aArm =="L"||aArm=="R"){ //Only use BCM info from right scaler
      gGet_Leaf_Mean_and_RMS(aT_Tree,"C_u1",  Form("right_current>=%f",Beam_Trip_Cut),aC_u1,rms); //Get Charge from Upstream BCM
//      gGet_Leaf_Mean_and_RMS(aT_Tree,"C_u3",  Form("right_current>=%f",Beam_Trip_Cut),aC_u3,rms); //Get Charge from Upstream BCM
      gGet_Leaf_Mean_and_RMS(aT_Tree,"C_d1",  Form("right_current>=%f",Beam_Trip_Cut),aC_d1,rms); //Get Charge from Upstream BCM
  //    gGet_Leaf_Mean_and_RMS(aT_Tree,"C_d3",  Form("right_current>=%f",Beam_Trip_Cut),aC_d3,rms); //Get Charge from Upstream BCM
      //charge = (aC_u1+aC_u3+aC_d1+aC_d3)/4.0*aEntries; //Sum = Mean*Entries    
      charge = (aC_u1+aC_d1)/2.0*aEntries; //Sum = Mean*Entries //ZYe, 07/18/2013, remove D3 and U3
  }
  return charge;
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
  
  //Double_t U1_Stat_Err=0.1;      //Fix Me
  //Double_t Offset_Stat_Err=0.1;  //Fix Me
  //Double_t Runtime_Stat_Err=1;   //Fix Me
  //Double_t Coeff_Stat_Err=7.5;   //Fix Me
   
  XGT2_VAR* output=new XGT2_VAR();
  output->Value=oCharge_Sum/(Qe*1e6);
  //output->Stat_Err=sqrt(pow(U1_Stat_Err,2)+pow(Offset_Stat_Err*oRuntime,2)+pow(Runtime_Stat_Err*offset_u,2)+pow(Coeff_Stat_Err*oCharge_Sum,2))/alpha_u/(Qe*1e6);
  //output->Stat_Err=sqrt(output->Value)/2.0; //A Factor of 1/2 is from the average of (u1,u3,d1,d3)
  output->Stat_Err=(output->Value * sqrt(oCharge_Stat_Sum))/2.0; //A Factor of 1/2 is from the average of (u1,u3,d1,d3)
  output->Sys_Err=0;
	
  return output;
}
/*}}}*/

/*inline gCheck_Current(TTree* aTree, const TString& aCuts){{{*/
inline void gCheck_Current(const vector<Int_t>& aRunNoChain, const TString aArm, double* aCurrent){
  double aCurrent_Left,aCurrent_Right,rms;
  
  for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
    {
      Int_t aRunNo=aRunNoChain[i];
      
      TChain *aTree = gAddTree(aRunNo, "T");
      gGet_Leaf_Mean_and_RMS(aTree,"left_current",  Form("left_current>=%f",Beam_Trip_Cut),aCurrent_Left,rms); //Get Current from Upstream BCM
      gGet_Leaf_Mean_and_RMS(aTree,"right_current", Form("right_current>=%f",Beam_Trip_Cut),aCurrent_Right,rms); //Get Current from Downstream BCM
	  if(aArm=="R"){
		  aCurrent[i] = aCurrent_Right;
	  }
	  else if(aArm=="L"){
		  //aCurrent[i] = aCurrent_Left;
		  aCurrent[i] = aCurrent_Right; //Left Scaler has problem
	  }
      else{
	   cerr<<"*** ERROR ***, I can not recognize the ARM = "<< aArm.Data()<<endl;
		  aCurrent[i] = -1.0;
	  }
      outlog<<Form("      For Run = %d, Current_Left = %6.2f, Current_Right = %6.2f",aRunNo, aCurrent_Left, aCurrent_Right)<<endl;
      cerr<<Form("      For Run = %d, Current_Left = %6.2f, Current_Right = %6.2f",aRunNo, aCurrent_Left, aCurrent_Right)<<endl;
      delete aTree;
    }
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

/*Double_t gGet_BeamCharge_Avg(TTree *aT_Tree,const TString& aArm){{{*/
Double_t gGet_BeamCharge_Avg(TTree *aT_Tree,const TString& aArm){
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

/*inline XGT2_VAR* gCal_Ne_Avg(c onst vector<Int_t>& aRunNoChain,const TString& aArm,Double_t& oU1,Double_t& oRuntime,Double_t& oCharge){{{*/
//Fix Here,ZYE
inline XGT2_VAR* gCal_Ne_Avg(const vector<Int_t>& aRunNoChain,const TString& aArm,Double_t& oRuntime)
{
  //return uC
  TString ArmName;
  Int_t ArmIndex;
  Double_t offset_u,offset_d,alpha_u, alpha_d;
  Double_t oCharge,oU1;
  if ( aArm=="R" )
    {
      ArmName="right";
      ArmIndex=1;
      offset_u = u_Offset[1][0];
      offset_d = d_Offset[1][0];
      alpha_u = u_Coeff[0][1];
      alpha_d = d_Coeff[0][1];
    }
  else
    {
      //have to put if else, otherwise ROOT will have error, no idea
      ArmName="left";
      ArmIndex=0;
      offset_u = u_Offset[0][0];
      offset_d = d_Offset[0][0];
      alpha_u = u_Coeff[0][0];
      alpha_d = d_Coeff[0][0];
    }
          
   Double_t oCharge_Sum=0.0,oRuntime_Sum=0.0, aCharge=0.0,aRuntime=0.0;     
   for ( unsigned int i=0; i<aRunNoChain.size(); i++ )
    {
      Int_t aRunNo=aRunNoChain[i];
     
      TString File_Form=Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo);
	  TChain *aT_Tree = gAddTree(aRunNo,"T"); 
      
	  aT_Tree->SetBranchAddress(Form("%s_bcm_u1c",ArmName.Data()), &oU1);
      aT_Tree->SetBranchAddress(Form("%s_clkcount",ArmName.Data()), &aRuntime);
      Int_t nentries=aT_Tree->GetEntries();
      aT_Tree->GetEntry(nentries-1);
      aRuntime /= 1024.;
      oRuntime_Sum += aRuntime;
      
      aCharge=(oU1-offset_u*aRuntime)/alpha_u;//uC
      oCharge_Sum += aCharge;
      //outlog << "oU1 = "<<oU1<<", offset_u = "<<offset_u<<endl;
      outlog <<"   For Run#"<<aRunNo << ", Run Time = "<<aRuntime <<", Ne = "<<aCharge/(Qe*1e6)<<endl;
      cerr <<"   For Run#"<<aRunNo << ", Run Time = "<<aRuntime <<", Ne = "<<aCharge/(Qe*1e6)<<endl;
      delete aT_Tree;
    } // for ( Int_t i=0; i<aRunNoChain.size(); i++ )
  
  oCharge = oCharge_Sum;
  oRuntime = oRuntime_Sum;
  
  Double_t U1_Stat_Err=0.1;      //Fix Me
  Double_t Offset_Stat_Err=0.1;  //Fix Me
  Double_t Runtime_Stat_Err=1;   //Fix Me
  Double_t Coeff_Stat_Err=7.5;   //Fix Me
  Double_t Ne=oCharge/(Qe*1e6);  //Fix Me
  
  XGT2_VAR* output=new XGT2_VAR();
  output->Value=Ne;
  output->Stat_Err=sqrt(pow(U1_Stat_Err,2)+pow(Offset_Stat_Err*oRuntime,2)+pow(Runtime_Stat_Err*offset_u,2)+pow(Coeff_Stat_Err*oCharge,2))/alpha_u/(Qe*1e6);
  //output->Stat_Err=sqrt(output->Value);
  output->Sys_Err=0;
	
  return output;
}
 
/*inline XGT2_VAR* gCal_Ne_Single(const Int_t aRunNo,const TString& aArm,Double_t& oRuntime){{{*/
inline XGT2_VAR* gCal_Ne_Single(const Int_t aRunNo,const TString& aArm,Double_t& oRuntime)
{
  //return uC
  TString ArmName;
  Double_t offset_u,offset_d,alpha_u, alpha_d;
  //if ( aArm=="R")
  if ( aArm=="R" || aArm=="L"){ //Use Right Arm Scaler since Left arm are totally bad.
    ArmName="right";
    offset_u = u_Offset[1][0];
    offset_d = d_Offset[1][0];
    alpha_u = u_Coeff[0][1];
    alpha_d = d_Coeff[0][1];
  }
  else{ //Don't use here, since Left arm are bad!
    ArmName="left";
    offset_u = u_Offset[0][0];
    offset_d = d_Offset[0][0];
    alpha_u = u_Coeff[0][0];
    alpha_d = d_Coeff[0][0];
  }
             
  Double_t aCharge=0.0,aRuntime=0.0;     
     
  TString File_Form=Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),aRunNo);
  TChain *aT_Tree = gAddTree(aRunNo,"T"); 
       
  aT_Tree->SetBranchAddress(Form("%s_clkcount",ArmName.Data()), &aRuntime);
  Int_t nentries=aT_Tree->GetEntries();
  aT_Tree->GetEntry(nentries-1);
  aRuntime /= 1024.;
  aCharge= gGet_BeamCharge(aT_Tree, aArm.Data());//uC
  Double_t Ne=aCharge/(Qe*1e6);  //Fix Me
  outlog <<"   For Run#"<<aRunNo << ", Run Time = "<<aRuntime <<", Ne = "<<Ne<<endl;
  cerr <<"   For Run#"<<aRunNo << ", Run Time = "<<aRuntime <<", Ne = "<<Ne<<endl;
  delete aT_Tree;
   
  //Double_t U1_Stat_Err=0.1;      //Fix Me
  //Double_t Offset_Stat_Err=0.1;  //Fix Me
  //Double_t Runtime_Stat_Err=1;   //Fix Me
  //Double_t Coeff_Stat_Err=7.5;   //Fix Me
   
  XGT2_VAR* output=new XGT2_VAR();
  output->Value=Ne;
  //output->Stat_Err=sqrt(pow(U1_Stat_Err,2)+pow(Offset_Stat_Err*oRuntime,2)+pow(Runtime_Stat_Err*offset_u,2)+pow(Coeff_Stat_Err*aCharge,2))/alpha_u/(Qe*1e6);
  //output->Stat_Err=sqrt(output->Value);
  output->Sys_Err=0;
        
  return output;
}
/*}}}*/

/*inline Bool_t gIsRaster(TTree* aTree,const TString& aRB_Name,const TString& aURB_Name,const TString& aBeamDir,Double_t& oCenter,Double_t& oWidth){{{*/
inline Bool_t gIsRaster(TTree* aTree,const TString& aRB_Name,const TString& aURB_Name,const TString& aBeamDir,Double_t& oCenter,Double_t& oWidth)
{
  Double_t rb_rawcur_mean[2];
  Double_t rb_rawcur_rms[2];
  gGet_Leaf_Mean_and_RMS(aTree,Form("%s.Raster.rawcur.x",aRB_Name.Data()),"",rb_rawcur_mean[0],rb_rawcur_rms[0]);
  gGet_Leaf_Mean_and_RMS(aTree,Form("%s.Raster.rawcur.y",aRB_Name.Data()),"",rb_rawcur_mean[1],rb_rawcur_rms[1]);

  if ( rb_rawcur_rms[0]>DQ_Name_Length && rb_rawcur_rms[1]>DQ_Name_Length )
    {
      gGet_Leaf_Mean_and_RMS(aTree,Form("%s.%s",aRB_Name.Data(),aBeamDir.Data()),"",oCenter,oWidth);
      oWidth*=2;
      return kTRUE;
    }
  else
    {
      gGet_Leaf_Mean_and_RMS(aTree,Form("%s.%s",aURB_Name.Data(),aBeamDir.Data()),"",oCenter,oWidth);
      oWidth*=2;
      return kFALSE;
    }
}
/*}}}*/

