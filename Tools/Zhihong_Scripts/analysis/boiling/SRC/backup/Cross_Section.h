/******************************************/
/*     X>2 Cross Section Calculation      */
/*         Zhihong Ye, 10/14/2011         */
/*   --Modified from H Yao's CSR code     */
/******************************************/
//#include "XGT2_Main.h"
#include "XGT2_LiveTime.h"
#include "XGT2_Beam.h"
#include "XGT2_Cuts.h"
#include "XGT2_Check.h"
#include "XGT2_Physics.h"

/* Formart of Input files:
   RunNo      --> Int_t,   First RunNumber, usally it would be the best runtime
   Arm        --> TString, "L" for HRS-L, "R" for HRS-R
   P          --> Float_t, Center Momentum of HRS, (GeV/c)
   Kin        --> TString, Which Kinematic Settings, Kin3.1, Kin3.2 ...
   Target     --> TString, Which Target, H2,He3,He4,C12,Ca40,Ca48 ...
   Tgt Thick  --> Float_t, Target Thickness, g/cm2
   Thick Err  --> Float_t, Target Thickness Sys Error
   A          --> Int_t,   Target A, 
   IsExtTgt   --> Int_t,   Is Extected Target? 1(Yes): H2,He3,He4, 0(No): C12,Ca40,Ca48
   Trigger    --> Int_t,   Main Trigger, 1,2,6 for HRS-R, 3,4,7 for HRS-L, usally we use 1 and 3 for R and L
   Scaler     --> TString, Scaler Tree, "LEFT" for HRS-L, "RIGHT" for HRS-R
   AccBinFile --> TString, AccBin File for Acc Correction coefficiences & Bin-Correction Coefficiencies, in ~/Bin
   RunChain   --> Int_t,   Run Numbers for Runs in the same setting, besides RunNo
   bCheck     --> Int_t,   Check Data and Cuts before doing Calculating, will plot PID,FP,TG distributions, 1/0 - Yes/No
   bCal_XS    --> Int_t,   Calculate Cross Section or not. Of cause we usally choose 1(Yes). 0(No) is for testing.
*/

/*int CrossSection(const TString& aFileName,const TString& aLogFileName){{{*/
int Cross_Section(const TString& aFileName,const TString& aLogFileName)
{
  //Build a logger
  XGT2_Logger* log=XGT2_Logger::Instance();
  log->SetLogFileName(aLogFileName);
  Int_t i=0,j=0;
  
  /////////////////////////////
  // Initialize Input values
  /////////////////////////////
  /*Read InputFile{{{*/
  cerr << endl << "@@@@@@ Initialize Input values ..."<<endl;
  Int_t RunNo=-1;
  Int_t aline=0;
  RunNo=gGet_InputFile_Var(aFileName,aline).Atoi(); aline++;
  const TString Arm=gGet_InputFile_Var(aFileName,aline);aline++;                cerr << "   ARM = " << Arm <<endl;
  const Double_t P=gGet_InputFile_Var(aFileName,aline).Atof();aline++;          cerr << "   P = " << P <<endl;
  const TString Kin=gGet_InputFile_Var(aFileName,aline);aline++;                cerr << "   Kin = " << Kin<<endl;
  const TString Target_Name=gGet_InputFile_Var(aFileName,aline);aline++;        cerr << "   Target Name = " << Target_Name<<endl;
  const Double_t Target_Thickness=gGet_InputFile_Var(aFileName,aline).Atof();aline++;  //check the unit! g/cm^2
  const Double_t Target_Thickness_Stat_Err=gGet_InputFile_Var(aFileName,aline).Atof();aline++; //check the unit! g/cm^2
  const Double_t Target_A=gGet_InputFile_Var(aFileName,aline).Atof();aline++;   cerr << "   Target A = " << Target_A<<endl;
  const Bool_t IsExtTgt=gGet_InputFile_Var(aFileName,aline).Atoi();aline++;     cerr << "   Extended Target? " << IsExtTgt<<endl;
  const Int_t Main_Trigger=gGet_InputFile_Var(aFileName,aline).Atoi();aline++;  cerr << "   Main Trigger " << Main_Trigger<<endl;
  const TString S_Tree_Name=gGet_InputFile_Var(aFileName,aline);aline++;        cerr << "   Scaler Tree " << S_Tree_Name<<endl;
  const TString AccBin_File=gGet_InputFile_Var(aFileName,aline);aline++;        cerr << "   AccBin File = " << AccBin_File<<endl;
  const TString T_Tree_Name="T";
  const TString E_Tree_Name="E";
  const TString DECAT_Name="DBB";
  RunNo=gAskRunNo(RunNo);
  if ( RunNo<0 ){
    log->WriteToLogFile(__FILE__,__LINE__,Form("Error No %d is input.",RunNo));
    log->Destroy();
    return -1;
  }
  const vector<Int_t> RunNoChain=gGet_RunNoChain(RunNo,gGet_InputFile_Var(aFileName,aline));aline++;	
  gCheckDirectory(Form("%s/results/%s",XGT2_DIR.Data(),gGetScriptname(__FILE__).Data()));
  const Bool_t bCheck_Data=gGet_InputFile_Var(aFileName,aline).Atoi();aline++;     cerr << "   Check Data? " << bCheck_Data<<endl;
  //const Bool_t bCal_XS = kFALSE;
  const Bool_t bCal_XS = gGet_InputFile_Var(aFileName,aline).Atoi();aline++;     cerr << "   Calculate Cross Section? " << bCal_XS<<endl; 
  /*}}}end of reading input*/
  
  ////////////////////////////////
  // Target Luminocity
  ////////////////////////////////
  /*Get Ntg{{{*/
  cerr << endl << "@@@@@@ Calculating Target Luminocity ..."<<endl;
  XGT2_VAR* Ntg=new XGT2_VAR();
  Ntg->Value=Target_Thickness/Target_A*Na; // number of target particles in unit area /cm^2
  Ntg->Stat_Err=Target_Thickness_Stat_Err*Na/Target_A;
  cerr <<"     Ntg = "<< Ntg->Value<<endl;
  log->WriteToLogFile(__FILE__,__LINE__,Form("     Ntg = %f",Ntg->Value));
  /*}}}end of Get Ntg*/
  
  ///////////////////////////////////
  // Number of Incident particles
  ///////////////////////////////////
  /*Get Ni{{{*/
  cerr << endl << "@@@@@@ Calculating Incident Particle ..."<<endl;	
  Double_t u1,runtime,charge;
  const XGT2_VAR* Ni=gCal_Ni(RunNoChain,Arm,u1,runtime,charge);
  cerr <<"     Ni = "<< Ni->Value<<endl;	
  log->WriteToLogFile(__FILE__,__LINE__,Form("     Ni = %f",Ni->Value));
  /*}}}end of Get Ni*/
  
  //////////////////////////
  //Check PreScaler
  //////////////////////////
  /*Get PS{{{*/
  cerr << endl << "@@@@@@ Checking PreScaler Factors for every run ..." <<endl;	
  gCheck_PS(RunNoChain);
  /*}}}end of PS*/
  
  //////////////////////////
  //Check LiveTime
  //////////////////////////
  /*LiveTime{{{*/
  cerr  << endl<< "@@@@@@ Checking Live Time for every run ..." <<endl; //For multi runs, the calculation results is incorrect if we chain all runs together!
  gCheck_LiveTime(RunNoChain,Arm);
  
  //Livetime Correction has been done in gGet_Nf(), so we only evaluate the Average value for all runs
  XGT2_VAR* Livetime = gGet_LiveTime_Avg(RunNoChain,Arm); 
  cerr << "      Average Lifetime = " << Livetime->Value<<endl;
  
  /*}}}end of LiveTime*/
  
  //////////////////////////////
  // Set Cuts
  //////////////////////////////
  /*Set Cuts{{{*/
  cerr << endl << "@@@@@@ Obtaining Cuts ..." <<endl;
  TString Cuts,ElectronCuts;
  
  cerr <<"     Getting Electron Cuts for "<<Target_Name <<" Run = "<<RunNo <<" on "<<Arm<<" Arm:"<<endl; 
  ElectronCuts=gGet_XGT2_Cut(RunNo,Arm,Target_Name,P);
  cerr <<"--------------------------------------------------------------------------------------"<<endl;
  cerr <<"     --> "<<ElectronCuts.Data()<<endl;
  cerr <<"--------------------------------------------------------------------------------------"<<endl;
  /*}}}end of Cuts*/
  
  //////////////
  // Chain Runs
  //////////////
  cerr << endl << "@@@@@@ Chain Runs ..."<<endl;
  TChain* T_Tree=gGetTree(RunNoChain,Arm,T_Tree_Name);
  
  //////////////////////////
  //Check Data and Cuts
  //////////////////////////
  /*Check Data and Cuts{{{*/
  if(bCheck_Data){
    cerr << endl << "@@@@@@ Checking Data and Cuts ..." <<endl;
    gCheck_Data_Cuts(RunNoChain,Arm,ElectronCuts.Data(),IsExtTgt);
  }
  /*}}}end of Check*/
  
  /*----------------------------------------------------------------------------------------------------------------------*/	
  ////////////////////////////////////
  // Calculate Cross Section
  ////////////////////////////////////
  cerr <<endl<< "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl;
  cerr       << "oO0Oo              --- Calculating Cross Section ---              oO0Oo" <<endl;
  cerr       << "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl<<endl;
  /////////////////////
  // Define variables
  /////////////////////
  XGT2_VAR* Total_E_Eff=new XGT2_VAR();
  XGT2_VAR* Total_Pion_Rej=new XGT2_VAR();
  Total_E_Eff->SetValue(1,0,0); //For X>2 Experiment, PID efficiencies are all above 99%, so we ignore the values and put them into sysmetic errors
  Total_Pion_Rej->SetValue(1,0,0);
  
  Double_t q3m, epsilon, xbj, qsq, scattering_angle, Ei, Ep, MC;
  XGT2_VAR* Cross_Section=new XGT2_VAR();
  XGT2_VAR* Nf=new XGT2_VAR();  //Total Number of real experiment events after all cuts
  XGT2_VAR* Acc=new XGT2_VAR(); //Acceptance Correction, dOmega_eff = dOmega*Acc->Value
  XGT2_VAR* Bin=new XGT2_VAR(); //Bin Center Correction, dSigma_Center = dSigma*Bin->Value
  
  Double_t rms;//tmp use
  TString outfile_name=Form("%s/results/%s/%s_%d.dat",XGT2_DIR.Data(),Target_Name.Data(),Kin.Data(),RunNo);
  cerr<<"      Save Output file to "<<outfile_name.Data()<<endl;
  ofstream outfile(outfile_name.Data());
  TString outstring=Form("%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s",
			 16,"XS(/MeV\xc2\xb2/sr)",
			 16,"XS_Stat_Err",
			 16,"XS_Sys_Err",
			 16,"thetatg",
			 16,"phitg",
			 16,"ytg(m)",
			 16,"dp",
			 16,"Q2",
			 16,"xbj",
			 16,"q3m",
			 16,"epsilon",
			 16,"angle(\xc2\xb0",
			 16,"Ei(MeV)",
			 16,"Ep(MeV)",
			 16,"Mott_C(/MeV\xc2\xb2/sr)",
			 16,"Acc_Cor",
			 16,"Bin_Cor",
			 16,"LiveTime",
			 16,"Ntg",
			 16,"Ni",
			 16,"Nf");
  outfile<<outstring<<endl;
			 
  TChain *T_Tree_Single = new TChain(T_Tree_Name);
  TString File_=Form("e08014_%d.root",RunNo);
  T_Tree_Single->Add(Form("%s/e08014_%d.root",ROOTFILES_DIR.Data(),RunNo));
			 
  ///////////////////////////
  //Get Beam Energy
  ///////////////////////////
  gGet_Leaf_Mean_and_RMS(T_Tree,"HALLA_p","HALLA_p>=0",Ei,rms);
  // gGet_Leaf_Mean_and_RMS(T_Tree_Single,"HALLA_p","HALLA_p>=0",Ei,rms);
  // gGet_Leaf_Mean_and_RMS(T_Tree,"HALLA_p",ElectronCuts,Ei,rms);
  if ( fabs(Ei)<1e-30 ){
    log->WriteToLogFile(__FILE__,__LINE__,Form("Error No Ei, %d ROOT File is too small?",RunNo));
    log->Destroy();      return -5;
  }
  cerr <<endl<< "      Beam Energy = " << Ei <<endl<<endl;
			 
  /////////////////////////////////////////
  // Scattered Electrons Momentum
  /////////////////////////////////////////
  /*Get P0{{{*/
  // cerr << endl << "@@@@@@ Calculating Center Momentum for every run ..."<<endl;	
  gCheck_P0(RunNoChain,Arm,E_Tree_Name);
			 
  //TChain* E_Tree=gGetTree(RunNoChain,Arm,E_Tree_Name);    //Read the EPICS values of last runs in the chain, somehow does not give the right ones.
  TChain* E_Tree=new TChain(E_Tree_Name);  //Aussming the first run get the most correct EPICS value
  E_Tree->Add(Form("%s/e08014_%d.root",ROOTFILES_DIR.Data(),RunNo));
  Double_t Dipole;
  XGT2_VAR* P0=gCal_P0(E_Tree,Arm,Dipole);
  if ( fabs(P0->Value)<1e-30 ) {
    log->WriteToLogFile(__FILE__,__LINE__,Form("Error No P0, %d ROOT File is too small?",RunNo));
    log->Destroy();
    return -3;
  }
  else if(fabs(P0->Value - P)>0.01) { //10MeV diviation would be large
    //    log->WriteToLogFile(__FILE__,__LINE__,"Error Wrong P0, Check Dipole field!");
    //     log->Destroy();
    //     return -4;
    P0->Value = P*1000; //If field is wrong for all runs, use the setting value
    P0->Stat_Err = 0.0;
    cerr <<"****** Cannot obtain correct P0 from Dipole field, use the setting value, P0 = "<<P0->Value<<endl;
			   
  }
  else {
    *P0 *= 1000;//GeV to MeV
    cerr << "      Using the value of first run, P0 = "<<P0->Value<<endl; 
  }

  log->WriteToLogFile(__FILE__,__LINE__,Form("     P0 = %6.2f",P0->Value));
  delete E_Tree;
  /*}}}end of Get P0*/

  ///////////////////////////
  // Get Scattering Angle
  ///////////////////////////
  // gGet_Leaf_Mean_and_RMS(T_Tree_Single,Form("EK%s.angle",Arm.Data()),ElectronCuts,scattering_angle,rms);
  if ( IsExtTgt )
    gGet_Leaf_Mean_and_RMS(T_Tree_Single,Form("ExtEK%s.angle",Arm.Data()),ElectronCuts,scattering_angle,rms);
  else
    gGet_Leaf_Mean_and_RMS(T_Tree_Single,Form("EK%s.angle",Arm.Data()),ElectronCuts,scattering_angle,rms);

  cerr <<endl<< "      Scattering Angle = " << scattering_angle*RadToDeg<<endl<<endl;; 
  delete T_Tree_Single; 

  /////////////////////////////////////////////
  // Define Acceptance Variables
  ////////////////////////////////////////////
  TString Data_Tg_Theta,Data_Tg_Phi,Data_Tg_Y,Data_Tg_Dp,Data_Q2,Data_Xbj;
  if ( IsExtTgt ){
    Data_Tg_Theta= Form("ExTgt%s.th",Arm.Data());
    Data_Tg_Phi= Form("ExTgt%s.ph",Arm.Data());
    Data_Tg_Y= Form("ExTgt%s.y",Arm.Data());
    Data_Tg_Dp= Form("ExTgt%s.dp",Arm.Data());
    Data_Q2= Form("ExtEK%s.Q2",Arm.Data());
    Data_Xbj=Form("ExtEK%s.x_bj",Arm.Data());
  }
  else{
    Data_Tg_Theta =Form("%s.gold.th",Arm.Data());
    Data_Tg_Phi =Form("%s.gold.ph",Arm.Data());
    Data_Tg_Y =Form("%s.gold.y",Arm.Data());
    Data_Tg_Dp =Form("%s.gold.dp",Arm.Data());
    Data_Q2 =Form("EK%s.Q2",Arm.Data());
    Data_Xbj=Form("EK%s.x_bj",Arm.Data());
  }

  //////////////////////////////
  // Get Bins
  //////////////////////////////
  /*Get Bins{{{*/
  Double_t AccBin[6][3], Acc_var_Value[6];
  //  gGet_Acc_Bin(RunNo,Arm,AccBin);
  //  gGet_Acc_Bin_ByKin(Arm.Data(),Kin.Data(),Target_Name.Data(),AccBin);
  //  gGet_Acc_Bin_ByName(AccBin_File.Data(),AccBin);
  gGet_Acc_Bin_Xbj(AccBin_File.Data(),AccBin);

  for ( i=0; i<6; i++ ){
    for ( j=0; j<3; j++){
      Printf("      AccBin[%d][%d]=%f",i,j,AccBin[i][j]);
    }
  }
  cerr<<endl;

  Int_t TotalBin = 1;    Int_t NumOfAccBin[6];
  for ( i=0; i<6; i++ ){
    NumOfAccBin[i]=TMath::Nint((AccBin[i][1]-AccBin[i][0])/AccBin[i][2]);
    if(i==0)
      cerr <<"      Bins for Theta_Tg is ";
    else if(i==1)
      cerr <<"      Bins for Phi_Tg is ";
    else if(i==2)
      cerr <<"      Bins for Y_Tg is ";
    else if(i==3)
      cerr <<"      Bins for Dp_Tg is ";
    else if(i==4)
      cerr <<"      Bins for Q2 is ";
    else if(i==5)
      cerr <<"      Bins for X_Bj is ";
    else
      cerr << "      Wrong!!!";
    cerr << NumOfAccBin[i]<<endl;
    TotalBin *= NumOfAccBin[i];
  }
  cerr <<endl<< "      Total Number of Bins = " << TotalBin<<endl<<endl;
  /*}}}end of Get Bins*/

  ///////////////////////////////////////////
  //Get Acceptance & Bin Corrections
  ///////////////////////////////////////////
  /*Get Corrections{{{*/
  TString Acc_Input_FileName=Form("%s/Acceptance/%s_%s_%s_Acceptance.dat",XGT2_DIR.Data(),Arm.Data(),Kin.Data(),Target_Name.Data());
  ifstream Acc_Input_File(Acc_Input_FileName.Data());
  Printf("  +++ Reading Acceptance File from %s.",Acc_Input_FileName.Data());
  if ( !Acc_Input_File.good() ){
    log->WriteToLogFile(__FILE__,__LINE__,
			Form("Error: %s Acceptance File is not found for Run %d",Acc_Input_FileName.Data(),RunNo));
    log->Destroy();
    return -4;
  }
  
  TString comment;
  comment.ReadLine(Acc_Input_File);
  Double_t acc_var[6],Acc_value,AccStatError,AccSysError,Bin_value,BinStatError,BinSysError;;
  Double_t thisbin[6];  //[6] 0:Theta 1:Phi 2:Y 3:Dp 4:Q2 5:Xbj
  Int_t NumOfAccInput=0, NoCal=0;;
  map<Double_t,Double_t> map_err;
  /*}}}end of Get Corrections*/

  /*Calculate and Save{{{*/
  if(bCal_XS){
    Ep=P0->Value;//ave of Ep
    MC=Get_Mott_CrossSection(scattering_angle*RadToDeg,Ei);
  
    if ( Total_Pion_Rej->Value>0&&Total_E_Eff->Value>0&&Ni->Value>0&&Ntg->Value>0&&Livetime->Value>0&&P0->Value>0 ){
    
      //Theta_Tg
      Acc_var_Value[0]=AccBin[0][0]+AccBin[0][2]/2;
      while ( Acc_var_Value[0]<AccBin[0][1] ){
	thisbin[0]=0;
	gGet_Leaf_Mean_and_RMS(T_Tree,Data_Tg_Theta,
			       Form("abs(%s-%f)<=%f",Data_Tg_Theta.Data(),Acc_var_Value[0],AccBin[0][2]/2),thisbin[0],rms);
	if ( thisbin[0]==0 )
	  thisbin[0]=Acc_var_Value[0];
	//Phi_Tg
	Acc_var_Value[1]=AccBin[1][0]+AccBin[1][2]/2;
	while ( Acc_var_Value[1]<AccBin[1][1] ){
	  thisbin[1]=0;
	  gGet_Leaf_Mean_and_RMS(T_Tree,Data_Tg_Phi,
				 Form("abs(%s-%f)<=%f",Data_Tg_Phi.Data(),Acc_var_Value[1],AccBin[1][2]/2),thisbin[1],rms);
	  if ( thisbin[1]==0 )
	    thisbin[1]=Acc_var_Value[1];
	  //Y_Tg
	  Acc_var_Value[2]=AccBin[2][0]+AccBin[2][2]/2;
	  while ( Acc_var_Value[2]<AccBin[2][1] ){
	    thisbin[2]=0;
	    gGet_Leaf_Mean_and_RMS(T_Tree,Data_Tg_Y,
				   Form("abs(%s-%f)<=%f",Data_Tg_Y.Data(),Acc_var_Value[2],AccBin[2][2]/2),thisbin[2],rms);
	    if ( thisbin[2]==0 )
	      thisbin[2]=Acc_var_Value[2];
	    //Dp_Tg
	    Acc_var_Value[3]=AccBin[3][0]+AccBin[3][2]/2;
	    while ( Acc_var_Value[3]<AccBin[3][1] ){
	      thisbin[3]=0;
	      gGet_Leaf_Mean_and_RMS(T_Tree,Data_Tg_Dp,
				     Form("abs(%s-%f)<=%f",Data_Tg_Dp.Data(),Acc_var_Value[3],AccBin[3][2]/2),thisbin[3],rms);
	      if ( thisbin[3]==0 )
		thisbin[3]=Acc_var_Value[3];
	      //Q2
	      Acc_var_Value[4]=AccBin[4][0]+AccBin[4][2]/2;
	      while ( Acc_var_Value[4]<AccBin[4][1] ){
		thisbin[4]=0;
		gGet_Leaf_Mean_and_RMS(T_Tree,Data_Q2,
				       Form("abs(%s-%f)<=%f",Data_Q2.Data(),Acc_var_Value[4],AccBin[4][2]/2),thisbin[4],rms);
		if ( thisbin[4]==0 )
		  thisbin[4]=Acc_var_Value[4];
		//X_Bj
		Acc_var_Value[5]=AccBin[5][0]+AccBin[5][2]/2;
		while ( Acc_var_Value[5]<AccBin[5][1] ){
		  thisbin[5]=0;
		  gGet_Leaf_Mean_and_RMS(T_Tree,Data_Xbj,
					 Form("abs(%s-%f)<=%f",Data_Xbj.Data(),Acc_var_Value[5],AccBin[5][2]/2),thisbin[5],rms);
		  if ( thisbin[5]==0 )
		    thisbin[5]=Acc_var_Value[5];
	  
		  if ( NoCal%5==0 ){
		    Printf(Form(" +++ Calculating %d Points ...", NoCal));
		  }
		  
		  cerr <<endl<< Form("      Reading Correction Factors of #%d Bins ...", NoCal) <<endl;
                  Acc_Input_File>>acc_var[0]>>acc_var[1]>>acc_var[2]>>acc_var[3]>>acc_var[4]>>acc_var[5]
                                >>Acc_value>>AccStatError>>AccSysError
                                >>Bin_value>>BinStatError>>BinSysError;    
                  cerr <<Form("       -- Theta = %f, Phi = %f, Y = %f, Dp = %f, Qsq = %f, Xbj = %f",acc_var[0],acc_var[1],acc_var[2],acc_var[3],acc_var[4],acc_var[5]) <<endl;
                  cerr <<Form("       -- Acc = %e, Acc_Stat_Err = %e, Acc_Sys_Err=%e", Acc_value,AccStatError,AccSysError)<<endl;
                  cerr <<Form("       -- Bin = %e, Bin_Stat_Err = %e, Bin_Sys_Err=%e", Bin_value,BinStatError,BinSysError)<<endl<<endl;
                  
		  
// 		  Cuts=Form("abs(%s-%f)<=%f&&abs(%s-%f)<=%f&&abs(%s-%f)<=%f&&abs(%s-%f)<=%f&&abs(%s-%f)<=%f&&abs(%s-%f)<=%f",
// 			    Data_Tg_Theta.Data(),
// 			    Acc_var_Value[0],
// 			    AccBin[0][2]/2,
// 			    Data_Tg_Phi.Data(),
// 			    Acc_var_Value[1],
// 			    AccBin[1][2]/2,
// 			    Data_Tg_Y.Data(),
// 			    Acc_var_Value[2],
// 			    AccBin[2][2]/2,
// 			    Data_Tg_Dp.Data(),
// 			    Acc_var_Value[3],
// 			    AccBin[3][2]/2,
// 			    Data_Q2.Data(),
// 			    Acc_var_Value[4],
// 			    AccBin[4][2]/2,
// 			    Data_Xbj.Data(),
// 			    Acc_var_Value[5],
// 			    AccBin[5][2]/2);

                  //Only bin on Theta, Phi,Y,and Dp, so I remove Xbj and Qsq bin. If we bin on xbj or Qsq instead of Dp, we can change here --Z. Ye, 01/09/2012
                  Cuts=Form("abs(%s-%f)<=%f&&abs(%s-%f)<=%f&&abs(%s-%f)<=%f&&abs(%s-%f)<=%f",
                            Data_Tg_Theta.Data(),
                            Acc_var_Value[0],
                            AccBin[0][2]/2,
                            Data_Tg_Phi.Data(),
                            Acc_var_Value[1],
                            AccBin[1][2]/2,
                            Data_Tg_Y.Data(),
                            Acc_var_Value[2],
                            AccBin[2][2]/2,
                            Data_Tg_Dp.Data(),
                            Acc_var_Value[3],
                            AccBin[3][2]/2);
		  Cuts+="&&"+ElectronCuts;
                  cerr << Cuts.Data()<<endl<<endl;
                 
		  if ( IsExtTgt ){
		    gGet_Leaf_Mean_and_RMS(T_Tree,Form("ExtEK%s.q3m",Arm.Data()),Cuts,q3m,rms);
		    gGet_Leaf_Mean_and_RMS(T_Tree,Form("ExtEK%s.epsilon",Arm.Data()),Cuts,epsilon,rms);
                    gGet_Leaf_Mean_and_RMS(T_Tree,Form("ExtEK%s.x_bj",Arm.Data()),Cuts,xbj,rms);
                    gGet_Leaf_Mean_and_RMS(T_Tree,Form("ExtEK%s.Q2",Arm.Data()),Cuts,qsq,rms);
		  }
		  else{
		    gGet_Leaf_Mean_and_RMS(T_Tree,Form("EK%s.q3m",Arm.Data()),Cuts,q3m,rms);
		    gGet_Leaf_Mean_and_RMS(T_Tree,Form("EK%s.epsilon",Arm.Data()),Cuts,epsilon,rms);
                    gGet_Leaf_Mean_and_RMS(T_Tree,Form("EK%s.x_bj",Arm.Data()),Cuts,xbj,rms);
                    gGet_Leaf_Mean_and_RMS(T_Tree,Form("EK%s.Q2",Arm.Data()),Cuts,qsq,rms);
                  }
																			  
		  //ROOT 5.18 bug: have to add this after GetEntries(sel) fixed in 5.20
		  T_Tree->SetNotify(0);
 
                  //DEBUG
		  // Printf("DEBUG --- thisbin[0]=%f,thisbin[1]=%f,thisbin[2]=%f,thisbin[3]=%f,theta=%f,phi=%f,ytg=%f,dp=%f,acc=%e",
		  //   Acc_var_Value[0],Acc_var_Value[1],Acc_var_Value[2],Acc_var_Value[3],acc_var[0],acc_var[1],acc_var[2],acc_var[3],Acc_value);
																			  
		  /*Calculate total real events from data after cuts{{{*/
		  //method 1, asuming all runs have same livetime and ps, which might not be true
		  //Nf->Value = PS[Main_Trigger-1]*T_Tree->GetEntries(Cuts.Data());
		  //method 2, correct livetime and ps for each runs, then sum all events
		  Nf->Value = gGet_Nf(RunNoChain,Arm,Cuts.Data());
																			  
		  //Nf->Value=T_Tree->GetEntries(ElectronCuts.Data());  
		  cerr << "      Total Number of Events = " << Nf->Value<<endl;
		  Nf->Stat_Err=sqrt(Nf->Value);
		  /*}}}*/
																			  
		  if ( Nf->Value>0 ){
		    Acc->SetValue(Acc_value,AccStatError,0);
		    Bin->SetValue(Bin_value,BinStatError,0);
		    // 		    Acc->SetValue(1.0,0.0,0); //Temp disable those corrections
		    // 		    Bin->SetValue(1.0,0.0,0);
																			    
		    if ( Acc->Value>0 ){
		      //    cerr <<"============================================================================================"<<endl;
		      //    Printf("Nf=%f,Total_Pion_Rej=%f,Total_E_Eff=%f,Acc=%f,Livetime=%f,Ni=%f,Ntg=%f,P0=%f,AccBin[3][2]=%f,AccBin[0][2]=%f,AccBin[1][2]=%f",
		      //    Nf->Value,Total_Pion_Rej->Value,Total_E_Eff->Value,Acc->Value,Livetime->Value,Ni->Value,Ntg->Value,P0->Value,AccBin[3][2],AccBin[0][2],AccBin[1][2]);
		      //    Printf("Nf_err=%f,Total_Pion_Rej_err=%f,Total_E_Eff_err=%f,Acc_err=%f,Livetime_err=%f,Ni_err=%f, Ntg_err=%f,P0_err=%f",
		      //    Nf->Stat_Err,Total_Pion_Rej->Stat_Err,Total_E_Eff->Stat_Err,Acc->Stat_Err,Livetime->Stat_Err,Ni->Stat_Err,Ntg->Stat_Err,P0->Stat_Err);
		      //    Printf("Nf_rel_err=%f,Total_Pion_Rej_rel_err=%f,Total_E_Eff_rel_err=%f,Acc_rel_err=%f,Livetime_rel_err=%f,Ni_rel_err=%f,Ntg_rel_err=%f,P0_rel_err=%f",
		      //    Nf->Stat_Err/Nf->Value,Total_Pion_Rej->Stat_Err/Total_Pion_Rej->Value,Total_E_Eff->Stat_Err/Total_E_Eff->Value,Acc->Stat_Err/Acc->Value,
		      //    Livetime->Stat_Err/Livetime->Value,Ni->Stat_Err/Ni->Value,Ntg->Stat_Err/Ntg->Value,P0->Stat_Err/P0->Value);
																			      
		      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		      //Since I only bin on X,Q2 or W2, but not dE,dTheta,dPhi, the cross section now is not differential but a integral value
		      //Fix Me Here: Dp, Theta,Phi also need correction factors even we do not bin on them, since we still have to make cuts on them.
		      // Cross_Section->Value = 
		      //	Nf->Value*Total_Pion_Rej->Value*Bin->Value*Acc->Value/(Total_E_Eff->Value*Ntg->Value*Ni->Value);
		      // unit=1
                      
                      // It is a differential cross section -- dXS/(dTheta*dPhi*dP. However, the actual definition of XS depends on what varialbes we bin on.
                      Cross_Section->Value = 
                        Nf->Value*Total_Pion_Rej->Value*Bin->Value/(Total_E_Eff->Value*Ntg->Value*Ni->Value*(Acc->Value*AccBin[0][2]*AccBin[1][2]*P0->Value*AccBin[3][2]));
		      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/					  
		      //    Printf("Proportional: Nf=%e,Pion_Rej=%e    inverse proportional: Total_E_Eff=%e,Acc=%e,Livetime=%e,Ni=%e,
		      //    Ntg=%e,P0=%e,AccBin[3][2]=%f,AccBin[0][2]=%f,AccBin[1][2]=%f",
		      //    Nf->Value,Total_Pion_Rej->Value,Total_E_Eff->Value,Acc->Value,Livetime->Value,Ni->Value,
		      //    Ntg->Value,P0->Value,AccBin[3][2],AccBin[0][2],AccBin[1][2]);
																			      
		      map_err.clear();
		      map_err.insert(make_pair(Nf->Value,Nf->Stat_Err));
		      map_err.insert(make_pair(Ni->Value,Ni->Stat_Err));
		      map_err.insert(make_pair(Ntg->Value,Ntg->Stat_Err));
		      map_err.insert(make_pair(Total_Pion_Rej->Value,Total_Pion_Rej->Stat_Err));
		      map_err.insert(make_pair(Total_E_Eff->Value,Total_E_Eff->Stat_Err));
		      map_err.insert(make_pair(Livetime->Value,Livetime->Stat_Err));
		      map_err.insert(make_pair(Acc->Value,Acc->Stat_Err));
		      map_err.insert(make_pair(Bin->Value,Bin->Stat_Err));
		      map_err.insert(make_pair(P0->Value,P0->Stat_Err));
		      Cross_Section->Stat_Err=gCal_Simple_Error(map_err,Cross_Section->Value);
																			      
		      //Fix me: Err
		      //Cross_Section->Stat_Err = 
		      // Cross_Section->Value*sqrt(pow(Nf->Stat_Err/Nf->Value,2)
		      //+pow(Total_Pion_Rej->Stat_Err/Total_Pion_Rej->Value,2)+pow(Total_E_Eff->Stat_Err/Total_E_Eff->Value,2)
		      //+pow(Acc->Stat_Err/Acc->Value,2)+pow(Livetime->Stat_Err/Livetime->Value,2)+pow(Ni->Stat_Err/Ni->Value,2)
		      //+pow(Ntg->Stat_Err/Ntg->Value,2)+pow(P0->Stat_Err/P0->Value,2));
																			      
		      Cross_Section->Sys_Err=0;
		      outstring=Form("%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e",
				     16,Cross_Section->Value,
				     16,Cross_Section->Stat_Err,
				     16,Cross_Section->Sys_Err,
				     16,thisbin[0],  //Theta_Tg
				     16,thisbin[1],  //Phi_Tg
				     16,thisbin[2],  //Y_Tg
				     16,thisbin[3],  //Dp_Tg
				     16,qsq,  //Q2
				     16,xbj,  //Xbj
				     16,q3m,
				     16,epsilon,
				     16,scattering_angle*RadToDeg,
				     16,Ei,
				     16,Ep,
				     16,MC,
				     16,Acc->Value,
				     16,Bin->Value,
				     16,Livetime->Value,
				     16,Ntg->Value,
				     16,Ni->Value,
				     16,Nf->Value);
		      outfile<<outstring<<endl;
		      Printf("      --> Cross_Section=%e,Xbj=%3.2f,MC=%e",Cross_Section->Value,thisbin[5],MC);
		      Printf("      --> Q2=%e,q3m=%e,epsilon=%e",thisbin[4],q3m,epsilon);    
		    }
		  }						      
		  NoCal++;
		  NumOfAccInput++;
		  Acc_var_Value[5]+=AccBin[5][2]; //Xbj
		}
		Acc_var_Value[4]+=AccBin[4][2];  //Q2
	      }
	      Acc_var_Value[3]+=AccBin[3][2]; //Dp_Tg
	    }
	    Acc_var_Value[2]+=AccBin[2][2];  //Y_Tg
	  }
	  Acc_var_Value[1]+=AccBin[1][2]; //Phi_Tg
	}
	Acc_var_Value[0]+=AccBin[0][2]; //Theta_Tg
      }
    }//if ( Total_Pion_Rej->Value>0&&Total_E_Eff->Value>0&&Ni->Value>0&&Ntg->Value>0&&Livetime->Value>0&&P0->Value>0 )
  }// If(bCal_XS)
  /*}}}*/

  outfile.close();
  delete T_Tree;
  log->WriteToLogFile(__FILE__,__LINE__,Form("%d Success. Result is saved at %s",RunNo,outfile_name.Data()));
  log->Destroy();
  return 0;
}
/*}}}*/
