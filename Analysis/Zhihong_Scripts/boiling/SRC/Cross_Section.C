/******************************************/
/*     X>2 Cross Section Calculation      */
/*         Zhihong Ye, 10/14/2011         */
/*         --Bin by Bin Correction        */
/******************************************/
/*Include{{{*/
#include "./SRC/XGT2_Main.h"
#include "./SRC/XGT2_LiveTime.h"
#include "./SRC/XGT2_Beam.h"
#include "./SRC/XGT2_Cuts.h"
#include "./SRC/XGT2_Data.h"
#include "./SRC/XGT2_Target.h"
//#include "./SRC/XGT2_Check.h"
#include "./SRC/XGT2_Physics.h"
#include "./SRC/XGT2_SAMC.h"
/*}}}*/

/* Formart of Input files:{{{*/
/* RunNo      --> Int_t,   First RunNumber, usally it would be the best runtime
   Arm        --> TString, "L" for HRS-L, "R" for HRS-R
   P          --> Float_t, Center Momentum of HRS, (GeV/c)
   Angle      --> Float_t, Center Scattering Angle of HRS, (Degree)
   Kin        --> TString, Which Kinematic Settings, Kin3.1, Kin3.2 ...
   Target     --> TString, Which Target, H2,He3,He4,C12,Ca40,Ca48 ...
   Tgt Thick  --> Float_t, Target Thickness, g/cm2
   Thick Err  --> Float_t, Target Thickness Sys Error
   A          --> Int_t,   Target A, 
   Z          --> Int_t,   Target Z, 
   IsExtTgt   --> Int_t,   Is Extected Target? 1(Yes): H2,He3,He4, 0(No): C12,Ca40,Ca48
   Trigger    --> Int_t,   Main Trigger, 1,2,6 for HRS-R, 3,4,7 for HRS-L, usally we use 1 and 3 for R and L
   Scaler     --> TString, Scaler Tree, "LEFT" for HRS-L, "RIGHT" for HRS-R
   AccBinFile --> TString, AccBin File for Acc Correction coefficiences & Bin-Correction Coefficiencies, in ~/Bin
   SAMC_File  --> TString, SAMC Monte Calo Simulation Data in this Kinematics
   RunChain   --> Int_t,   Run Numbers for Runs in the same setting, besides RunNo
   Beam_Trip  --> Float_t, Beam Trip Cut Points, like Current > Beam_Trip
   bCheck     --> Int_t,   Check Data and Cuts before doing Calculating, will plot PID,FP,TG distributions, 1/0 - Yes/No
   bCal_XS    --> Int_t,   Calculate Cross Section or not. Of cause we usally choose 1(Yes). 0(No) is for testing. */
/*}}}*/

/*int CrossSection(const TString& aInputFile,const TString& aLogFile){{{*/
int CrossSection(const TString& aInputFile,const TString& aLogFile){
  
  TString aInputFileName = "./xs_input/"+aInputFile;
  TString aLogFileName = "./log/"+aLogFile;
  outlog.open(aLogFile); 
  //Build a logger
  XGT2_Logger* log=XGT2_Logger::Instance();
  log->SetLogFileName(aLogFileName);
  Int_t i=0,j=0;

  /////////////////////////////
  // Initialize Input values
  /////////////////////////////
  //Set Value of Beam Energy
  Double_t Ebeam = 3.356; //GeV
  /*Read InputFile{{{*/
  outlog << endl << Form("@@@@@@ Initialize Input values from %s ...",aInputFileName.Data())<<endl;
  cerr   << endl << Form("@@@@@@ Initialize Input values from %s ...",aInputFileName.Data())<<endl;
  Int_t RunNo=-1;
  Int_t aline=0;
  RunNo=gGet_InputFile_Var(aInputFileName,aline).Atoi(); aline++;
  const TString Arm=gGet_InputFile_Var(aInputFileName,aline);aline++;                outlog << "   ARM = " << Arm <<endl;
  const Double_t P=gGet_InputFile_Var(aInputFileName,aline).Atof();aline++;          outlog << "   P0 = " << P <<endl;
  const Double_t Angle=gGet_InputFile_Var(aInputFileName,aline).Atof();aline++;      outlog << "   Angle = " << Angle <<endl;
  const TString Kin=gGet_InputFile_Var(aInputFileName,aline);aline++;                outlog << "   Kin = " << Kin<<endl;
  const TString Target_Name=gGet_InputFile_Var(aInputFileName,aline);aline++;        outlog << "   Target Name = " << Target_Name<<endl;
  Double_t Target_Thickness=gGet_InputFile_Var(aInputFileName,aline).Atof();aline++;  outlog << Form("   Target Thickness = %f g/cm^2",Target_Thickness)<<endl; //check the unit! g/cm^2
  const Double_t Target_Thickness_Stat_Err=gGet_InputFile_Var(aInputFileName,aline).Atof();aline++; //check the unit! g/cm^2
  const Int_t Target_A=gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++;      outlog << "   Target A = " << Target_A<<endl;
  const Int_t Target_Z=gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++;      outlog << "   Target Z = " << Target_Z<<endl;
  const Bool_t IsExtTgt=gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++;     outlog << "   Extended Target? " << IsExtTgt<<endl;
  const Int_t Main_Trigger=gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++;  outlog << "   Main Trigger " << Main_Trigger<<endl;
  const TString S_Tree_Name=gGet_InputFile_Var(aInputFileName,aline);aline++;        outlog << "   Scaler Tree " << S_Tree_Name<<endl;
  const TString AccBin_File=gGet_InputFile_Var(aInputFileName,aline);aline++;        outlog << "   AccBin File = " << AccBin_File<<endl;
  const TString SAMC_File=gGet_InputFile_Var(aInputFileName,aline);aline++;          outlog << "   SAMC File = " << SAMC_File<<endl;
  const TString T_Tree_Name="T";
  const TString E_Tree_Name="E";
  const TString DECAT_Name="DBB";
  RunNo=gAskRunNo(RunNo);
  if ( RunNo<0 ){
    log->WriteToLogFile(__FILE__,__LINE__,Form("Error No %d is input.",RunNo));
    log->Destroy();
    return -1;
  }
  const vector<Int_t> RunNoChain=gGet_RunNoChain(RunNo,gGet_InputFile_Var(aInputFileName,aline));aline++;	
  const Int_t Chain_Size = RunNoChain.size();
  gCheckDirectory(Form("%s/results/%s",XGT2_DIR.Data(),gGetScriptname(__FILE__).Data()));
  Beam_Trip_Cut = gGet_InputFile_Var(aInputFileName,aline).Atof();aline++; outlog << Form("   Beam Tip Cut = %6.2f uA",Beam_Trip_Cut)<<endl;
  const Bool_t bCheck_Data= gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++;   outlog << "   Check Data? " << bCheck_Data<<endl;
  //const Bool_t bCal_XS = kFALSE;
  const Bool_t bCal_XS = gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++;     outlog << "   Calculate Cross Section? " << bCal_XS<<endl; 
  /*}}}end of reading input*/
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // End of Initializing Input values
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
  //////////////////
  // Set Cuts
  //////////////////
  /*Set Cuts{{{*/
  outlog << endl << "@@@@@@ Obtaining Cuts ..." <<endl;
  cerr   << endl << "@@@@@@ Obtaining Cuts ..." <<endl;
  TString EX_Cuts,ElectronCuts;

  outlog <<"     Getting Electron Cuts for "<<Target_Name <<" Run = "<<RunNo <<" on "<<Arm<<" Arm:"<<endl; 
  cerr   <<"     Getting Electron Cuts for "<<Target_Name <<" Run = "<<RunNo <<" on "<<Arm<<" Arm:"<<endl; 
  ElectronCuts=gGet_XGT2_Cut(RunNo,Arm,Target_Name.Data(),Kin.Data(),P,Main_Trigger); //Target_Length is also given in this subroutine
  outlog <<"--------------------------------------------------------------------------------------"<<endl;
  cerr <<"--------------------------------------------------------------------------------------"<<endl;
  outlog <<"     --> "<<ElectronCuts.Data()<<endl;
  cerr   <<"     --> "<<ElectronCuts.Data()<<endl;
  outlog <<"--------------------------------------------------------------------------------------"<<endl;
  cerr   <<"--------------------------------------------------------------------------------------"<<endl;
  /*}}}end of Cuts*/

  ////////////////////////////////
  // Target Luminocity
  ////////////////////////////////
  /*Get Ntg{{{*/
  outlog << endl << "@@@@@@ Calculating Target Luminocity ..."<<endl;
  cerr   << endl << "@@@@@@ Calculating Target Luminocity ..."<<endl;
        
  if(IsExtTgt){
    outlog<< endl << "       Adjust Target Length for Extanded Targets --- "
	<< Form("Length after cut for %s is: %5.3f cm.", Target_Name.Data(), Target_Length)
	<<endl;
    cerr  << endl << "       Adjust Target Length for Extanded Targets --- "
	<< Form("Length after cut for %s is: %5.3f cm.", Target_Name.Data(), Target_Length)
	<<endl;
    Target_Thickness *= Target_Length/Target_Length_Default;
  }
        
  XGT2_VAR* Ntg= gCal_Ntg(RunNoChain,Arm.Data(),Target_A,Target_Thickness,Target_Thickness_Stat_Err);
  //      Not corrected by boiling effect
  //      XGT2_VAR* Ntg = new XGT2_VAR();
  // 	Ntg->Value=Target_Thickness/Target_A*Na; // number of target particles in unit area /cm^2
  // 	Ntg->Stat_Err=Target_Thickness_Stat_Err*Na/Target_A;
  outlog <<Form("     Avg Ntg = %4.3e; Before Boiling Corrected = %4.3e",Ntg->Value,Target_Thickness/Target_A*Na)<<endl;
  cerr   <<Form("     Avg Ntg = %4.3e; Before Boiling Corrected = %4.3e",Ntg->Value,Target_Thickness/Target_A*Na)<<endl;
  log->WriteToLogFile(__FILE__,__LINE__,Form("     Ntg = %f, No Boiling Corrected = %f",Ntg->Value,Target_Thickness/Target_A*Na));
  /*}}}end of Get Ntg*/

  ///////////////////////////////////
  // Number of Incident particles
  ///////////////////////////////////
  /*Get Ne{{{*/
  outlog << endl << "@@@@@@ Calculating Incident Particle ..."<<endl;	
  cerr   << endl << "@@@@@@ Calculating Incident Particle ..."<<endl;	
  Double_t runtime;
  const XGT2_VAR* Ne=gCal_Ne(RunNoChain,Arm,runtime);
  outlog <<Form("     Ne = %4.3e", Ne->Value)<<endl;	
  cerr   <<Form("     Ne = %4.3e", Ne->Value)<<endl;	
  log->WriteToLogFile(__FILE__,__LINE__,Form("     Ne = %f, Runtime = %f",Ne->Value, runtime));
  /*}}}end of Get Ne*/

  //////////////////////////
  // Check PreScalerR
  //////////////////////////
  /*Get PS{{{*/
  outlog << endl << "@@@@@@ Checking and Storing PreScaler Factors for every run ..." <<endl;	
  cerr   << endl << "@@@@@@ Checking and Storing PreScaler Factors for every run ..." <<endl;	
  Int_t* PSChain = new int[Chain_Size*8];
  gCheck_PS(RunNoChain,PSChain);
  /*}}}end of PS*/

  //////////////////////////
  // Check LiveTime
  //////////////////////////
  /*LiveTime{{{*/
  outlog  << endl<< "@@@@@@ Checking and Storing Live Time for every run ..." <<endl; //For multi runs, the calculation results is incorrect if we chain all runs together!
  cerr    << endl<< "@@@@@@ Checking and Storing Live Time for every run ..." <<endl; //For multi runs, the calculation results is incorrect if we chain all runs together!
  Double_t* LiveTimeChain = new double[Chain_Size];
  Double_t* LiveTimeErrorChain = new double[Chain_Size]; //Store the Error of Livetime. Will use it for total Error_Bar calculation;
  for(int i=0;i<Chain_Size;i++){ 
    LiveTimeChain[i] = -1.0; //Initialize the livetime
    LiveTimeErrorChain[i] = -1.0;
  }
  
  gLoad_LiveTime_Chain(RunNoChain,Arm,LiveTimeChain,LiveTimeErrorChain,1);//1->Beam Trip Cut, 2->NONE Cuts
  XGT2_VAR* Livetime = new XGT2_VAR();
  Livetime->Value = 0.0; Livetime->Stat_Err = 0.0;
  for(int i=0;i<Chain_Size;i++){ 
    Livetime->Value += LiveTimeChain[i];
    Livetime->Stat_Err += LiveTimeErrorChain[i];
  }
  Livetime->Value/=Chain_Size; Livetime->Stat_Err/=Chain_Size;
  
  //Check livetime for each run, and store all values for future use
  //gCheck_LiveTime(RunNoChain,Arm, LiveTimeChain,LiveTimeErrorChain,PSChain,Beam_Trip_Cut);  
  //Livetime Correction has been done in gGet_Nf(), so we only evaluate the Average value for all runs
 // XGT2_VAR* Livetime = gGet_LiveTime_Avg(RunNoChain,Arm,PSChain,Beam_Trip_Cut); 
  
  

  outlog << "      Average Lifetime = " << Livetime->Value<<endl;
  cerr   << "      Average Lifetime = " << Livetime->Value<<endl;

  /*}}}end of LiveTime*/

  //////////////////////
  // Check Current
  //////////////////////
  gCheck_Current(RunNoChain);

  ///////////////////////////////
  //Check Central Momentum 
  ///////////////////////////////  
  outlog  << endl<< "@@@@@@ Checking Central Momentum for each run ..." <<endl; //
  cerr    << endl<< "@@@@@@ Checking Central Momentum for each run ..." <<endl; //
  gCheck_P0(RunNoChain,Arm,E_Tree_Name);  //FIX_HERE: Check P0 from Quadropole too!


  ///////////////////////////////
  // Check Scattering Angle
  ///////////////////////////////
  /*Check Scattering Angle{{{*/
  // outlog  << endl<< "@@@@@@ Checking Scattering Angle for each run ..." <<endl; //
  // gCheck_Angle(RunNoChain,Arm,T_Tree_Name,ElectronCuts.Data(),IsExtTgt);
  // //FIX_HERE: Add this feature!
  /*}}}*/

  //////////////////////////
  // Check Data and Cuts
  //////////////////////////
  /*Check Data and Cuts{{{*/
  //if(bCheck_Data){
  //  outlog << endl << "@@@@@@ Checking Data and Cuts ..." <<endl;
  //  gCheck_Data_Cuts(RunNoChain,Arm,ElectronCuts.Data(),IsExtTgt);
  //}
  /*}}}end of Check*/

  //////////////
  // Chain Runs
  //////////////
  outlog << endl << "@@@@@@ Chain Runs ..."<<endl;
  cerr   << endl << "@@@@@@ Chain Runs ..."<<endl;
  TChain* T_Tree=gGetTree(RunNoChain,Arm,T_Tree_Name);

  ////////////////////////////////////
  // Calculate Cross Section
  ////////////////////////////////////
  outlog <<endl<< "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl;
  outlog       << "oO0Oo              --- Calculating Cross Section ---              oO0Oo" <<endl;
  outlog       << "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl<<endl;

  cerr <<endl<< "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl;
  cerr       << "oO0Oo              --- Calculating Cross Section ---              oO0Oo" <<endl;
  cerr       << "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl<<endl;

  /////////////////////
  // Define variables
  /////////////////////
  /*Define variables{{{*/
  XGT2_VAR* Total_E_Eff=new XGT2_VAR();
  XGT2_VAR* Total_Pion_Rej=new XGT2_VAR();
  //For X>2 Experiment, PID efficiencies are all above 99%, so we ignore the values and put them into sysmetic errors
  Total_E_Eff->SetValue(1,0,0); 
  Total_Pion_Rej->SetValue(1,0,0);

  Double_t nu=-1000., xbj=-1000., qsq=-1000., Ep=-1000.0; 
  XGT2_VAR* Cross_Section=new XGT2_VAR();
  XGT2_VAR* Nf_EX=new XGT2_VAR();  //Total Number of real events in each bin after all cuts
  XGT2_VAR* Nf_SAMC=new XGT2_VAR();  //Total Number of simulation events in each bin after all cuts

  //Double_t rms;//tmp use
  TString outfile_name=Form("%s/results/%s/%s_%s_%d.dat",XGT2_DIR.Data(),Target_Name.Data(),Arm.Data(),Kin.Data(),RunNo);
   outlog<<"      Backup and Save Output file to "<<outfile_name.Data()<<endl;
   cerr  <<"      Backup and Save Output file to "<<outfile_name.Data()<<endl;
  gSystem->Exec("cp -vpf "+outfile_name+" "+outfile_name+".old");
  gSystem->Exec("chmod -w "+outfile_name+".old");
  gSystem->Exec("chmod +w "+outfile_name);
   
  ofstream outfile(outfile_name.Data());
  TString outstring=Form("%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s",
			 16,"XS(cm\xc2\xb2/MeV/sr)",
			 16,"XS_Stat_Err",
			 16,"XS_Sys_Err",
			 16,"Xbj",
			 16,"Angle",
			 16,"Ep(MeV/c)",
			 16,"Nu",
			 16,"Q2",
			 16,"Ei(MeV)",
			 16,"P0(MeV/c)",
			 16,"LiveTime",
			 16,"Ntg",
			 16,"Ne",
			 16,"Nf_EX",
			 16,"Nf_SAMC",
			 16,"Nf_SAMC_Gen");
  outfile<<outstring<<endl;

  Ebeam*=1000.0; //GeV to MeV;
  Double_t P_Center = P*1000.0; //GeV to MeV;
  Double_t SinA = sin(Angle*TMath::DegToRad());
  Double_t CosA = cos(Angle*TMath::DegToRad());
  Double_t SinSQ = pow(sin(Angle*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2
  /*}}}*/
	
  ///////////////////////////
  //Get Beam Energy
  ///////////////////////////
  TChain *T_Tree_Single = new TChain(T_Tree_Name);
  T_Tree_Single->Add(Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),RunNo));
        
  /*Get Beam Energy{{{*/
  // Double Ei;
  // gGet_Leaf_Mean_and_RMS(T_Tree,"HALLA_p","HALLA_p>=0",Ei,rms);
  // // gGet_Leaf_Mean_and_RMS(T_Tree_Single,"HALLA_p","HALLA_p>=0",Ei,rms);
  // // gGet_Leaf_Mean_and_RMS(T_Tree,"HALLA_p",ElectronCuts,Ei,rms);
  // if ( fabs(Ei)<1e-30 ){
  //   log->WriteToLogFile(__FILE__,__LINE__,Form("Error No Ei, %d ROOT File is too small?",RunNo));
  //   log->Destroy();      return -5;
  // }
  // outlog <<endl<< "++++++ Beam Energy = " << Ei <<endl<<endl;
  /*}}}*/

  /////////////////////////////////////////
  // Scattered Electrons Momentum
  /////////////////////////////////////////
  /*Get P0{{{*/
  //TChain* E_Tree=gGetTree(RunNoChain,Arm,E_Tree_Name);    //Read the EPICS values of last runs in the chain, somehow does not give the right ones.
  TChain* E_Tree=new TChain(E_Tree_Name);  //Aussming the first run get the most correct EPICS value
  E_Tree->Add(Form("%s/%s_%d.root",ROOTFILES_DIR.Data(),ROOTFILES_NAME.Data(),RunNo));
  Double_t Dipole;
  XGT2_VAR* P0=gCal_P0(E_Tree,Arm,Dipole);
  if ( fabs(P0->Value)<1e-30 ) {
    log->WriteToLogFile(__FILE__,__LINE__,Form("Error No P0, %d ROOT File is too small?",RunNo));
    log->Destroy();
    return -3;
  }
  else if(fabs(P0->Value - P)>0.01) { //10MeV diviation would be large
    P0->Value = P*1000; //If field is wrong for all runs, use the setting value
    P0->Stat_Err = 0.0;
    outlog <<"****** Cannot obtain correct P0 from Dipole field, use the setting value, P0 = "<<P0->Value<<endl;
    cerr   <<"****** Cannot obtain correct P0 from Dipole field, use the setting value, P0 = "<<P0->Value<<endl;
  }
  else {
    *P0 *= 1000;//GeV to MeV
    outlog << "++++++ Using the value of first run, P0 = "<<P0->Value<<endl; 
    cerr   << "++++++ Using the value of first run, P0 = "<<P0->Value<<endl; 
  }

  log->WriteToLogFile(__FILE__,__LINE__,Form("     P0 = %6.2f",P0->Value));
  delete E_Tree;
  /*}}}end of Get P0*/

  // ///////////////////////////
  // // Get Scattering Angle
  // ///////////////////////////
  /*Get Scattering Angle{{{*/
  // if ( IsExtTgt )
  //   gGet_Leaf_Mean_and_RMS(T_Tree_Single,Form("ExtEK%s.angle",Arm.Data()),ElectronCuts,scattering_angle,rms);
  // else
  //   gGet_Leaf_Mean_and_RMS(T_Tree_Single,Form("EK%s.angle",Arm.Data()),ElectronCuts,scattering_angle,rms);
  // outlog <<endl<< "++++++ Scattering Angle = " << scattering_angle*RadToDeg<<endl<<endl;; 
  /*}}}*/
  delete T_Tree_Single; 

  /////////////////////////////////////////////
  // Define Acceptance Variables
  ////////////////////////////////////////////
  /*Define Acceptance Variables{{{*/
  //For Real Data
  TString EX_Dp,EX_Nu,EX_Q2,EX_Xbj;
  if ( IsExtTgt ){
    EX_Dp = Form("ExTgt%s.dp",Arm.Data());
    EX_Th = Form("ExTgt%s.th",Arm.Data());
    EX_Ph = Form("ExTgt%s.ph",Arm.Data());
	EX_SA = From("acos((%f-%s*%f)/(sqrt(1.0+pow(%s,2)+pow(%s,2))))",CosA,EX_Ph.Data(),SinA,EX_Th.Data(),EX_Ph.Data());
 	EX_Nu = Form("ExtEK%s.nu",Arm.Data());
    EX_Q2 = Form("ExtEK%s.Q2",Arm.Data());
    EX_Xbj= Form("ExtEK%s.x_bj",Arm.Data());
  }
  else{
    EX_Dp = Form("%s.gold.dp",Arm.Data());
	Ex_Th = Form("%s.gold.th",Arm.Data());
	Ex_Ph = Form("%s.gold.ph",Arm.Data());
	EX_SA = From("acos((%f-%s*%f)/(sqrt(1.0+pow(%s,2)+pow(%s,2))))/%f",CosA,EX_Ph.Data(),SinA,EX_Th.Data(),EX_Ph.Data(),TMath::DegToRad());
    EX_Nu = Form("EK%s.nu",Arm.Data());
    EX_Q2 = Form("EK%s.Q2",Arm.Data());
    EX_Xbj= Form("EK%s.x_bj",Arm.Data());
  }
  //Current Branch for cutting out Beam Trip events
  // must run "AddCurrent" first to add those branches in rootfiles

  TString EX_Current =" ";
  if(Arm=="R"||Arm=="L")
    EX_Current = "right_current"; //LEFT BCM had problem for x>2 experiment
  else
    EX_Current = "left_current";
 
  /*}}}*/

  //////////////////////////////
  // Get Bins
  //////////////////////////////
  /*Get Bins{{{*/
  Double_t AccBin[4][3], Acc_var_Value[4];
  outlog <<"++++++ Getting Acceptance Binning from File "<<AccBin_File.Data()<<endl;
   cerr  <<"++++++ Getting Acceptance Binning from File "<<AccBin_File.Data()<<endl;
  gGet_Acc_Bin(AccBin_File.Data(),AccBin);
  for ( i=0; i<4; i++ ){
    for ( j=0; j<3; j++){
      outlog << Form("      AccBin[%d][%d]=%f",i,j,AccBin[i][j])<<endl;;
      cerr   << Form("      AccBin[%d][%d]=%f",i,j,AccBin[i][j])<<endl;;
    }
  }
  outlog<<endl;
  cerr<<endl;

  Int_t TotalBin = 1;    Int_t NumOfAccBin[4];
  for ( i=0; i<4; i++ ){
    NumOfAccBin[i]=TMath::Nint((AccBin[i][1]-AccBin[i][0])/AccBin[i][2]);
    if(i==0){
	  outlog <<"      Bins for Dp is ";
      cerr   <<"      Bins for Dp is ";
	}
	else if(i==1){
      outlog <<"      Bins for Nu is ";
      cerr   <<"      Bins for Nu is ";
	}
    else if(i==2){
      outlog <<"      Bins for Q2 is ";
      cerr   <<"      Bins for Q2 is ";
	}
    else if(i==3){
      outlog <<"      Bins for Xbj is ";
      cerr   <<"      Bins for Xbj is ";
	}
    else{
      outlog << "      Wrong!!!";
      cerr   << "      Wrong!!!";
	}
    outlog << NumOfAccBin[i]<<endl;
    cerr   << NumOfAccBin[i]<<endl;
    TotalBin *= NumOfAccBin[i];
  }
  outlog <<endl<< "      Total Number of Bins = " << TotalBin<<endl<<endl;
  cerr   <<endl<< "      Total Number of Bins = " << TotalBin<<endl<<endl;
  /*}}}end of Get Bins*/

  // Double_t thisbin[4];  //[6] 0:Dp 1:Nu 2:Q2 3:Xbj
  Int_t NumOfAccInput=0, NoCal=0;;
  map<Double_t,Double_t> map_err;

  /////////////
  // SAMC
  /////////////
  /*Get Info From SAMC Simulation{{{*/
  outlog <<endl<<"++++++ Info in SAMC:"<<endl;
  cerr   <<endl<<"++++++ Info in SAMC:"<<endl;
  Double_t SAMC_dDp    = 0.12;    outlog << Form("       Full Range of Dp    in SAMC = %f PCT",  SAMC_dDp)<<endl;     // MeV/c
  Double_t SAMC_dTheta = 0.18;    outlog << Form("       Full Range of Theta in SAMC = %f Rad",  SAMC_dTheta)<<endl;  // Rad  
  Double_t SAMC_dPhi   = 0.09;    outlog << Form("       Full Range of Phi   in SAMC = %f Rad",  SAMC_dPhi)<<endl;    // Rad
  Double_t SAMC_dOmega = SAMC_dTheta*SAMC_dPhi;
  Double_t SAMC_BinCut[4]; //0->Dp, 1->Nu, 2->Qsq, 3->Xbj

  //Chain Data
  TChain* SAMC_Tree=new TChain("SAMC"); 
  TString SAMC_RootFile = SAMCFILES_DIR+SAMC_File;
  SAMC_Tree->Add(SAMC_RootFile.Data());
  outlog<<Form("       Added SAMC Rootfile %s.",SAMC_RootFile.Data())<<endl;
  cerr  <<Form("       Added SAMC Rootfile %s.",SAMC_RootFile.Data())<<endl;
  XGT2_VAR* Nf_SAMC_Gen_All=new XGT2_VAR();  //Total Number of simulation events in the file, without any cuts 
  Nf_SAMC_Gen_All->Value = SAMC_Tree->GetEntries();
  Nf_SAMC_Gen_All->Stat_Err = sqrt(Nf_SAMC_Gen_All->Value);
  outlog << Form("       Total Events generated in SAMC = %f ",Nf_SAMC_Gen_All->Value)<<endl;
  cerr   << Form("       Total Events generated in SAMC = %f ",Nf_SAMC_Gen_All->Value)<<endl;
  outlog<<endl;
  cerr  <<endl;

  //Bin Variable3s For SAMC Data
  TString SAMC_Dp,SAMC_Nu,SAMC_Q2,SAMC_Xbj;
  SAMC_Dp = "dp_rec"; 
  SAMC_Nu = "(E_s-E_p)";
  SAMC_Q2 = "Qsq";
  SAMC_Xbj= "Xbj";

  //Get SAMC Accptance Cuts. Read the values from the same Acc File as for Data
  Double_t SAMC_AccCut[9]; //0->x_fp,1->y_fp,2->th_fp,3->ph_fp,4->dp_tg,5->y_tg,6->th_tg,7->ph_tg,8->reactz
  gGet_XGT2_Acc_Cut(Arm, Target_Name, P, SAMC_AccCut);

  /*SAMC Cuts, Not in use now*/
  // outlog <<"     Getting SAMC Cuts for "<<Target_Name <<" on "<<Arm<<" Arm (should be the same as cuts in EX data):"<<endl; 
  // TString SAMC_Cuts = gGet_XGT2_Cut_SAMC(Arm,Target_Name,P);
  // outlog <<"-------- Acceptance Cuts -------------------------------------------------------------"<<endl;
  // outlog <<"     --> "<<SAMC_Cuts.Data()<<endl;
  // outlog <<"--------------------------------------------------------------------------------------"<<endl;
	
  /*}}}*/

  ///////////////////////////////////
  // Calculate XS and Save
  ///////////////////////////////////
  /*Calculate XS and Save{{{*/
  if(bCal_XS){
    if ( Total_Pion_Rej->Value>0&&Total_E_Eff->Value>0&&Ne->Value>0&&Ntg->Value>0&&Livetime->Value>0&&P0->Value>0 ){
      /*Dp{{{*/
      Acc_var_Value[0]=AccBin[0][0]+AccBin[0][2]/2;
      while ( Acc_var_Value[0]<=AccBin[0][1] ){
	/*Nu{{{*/
	Acc_var_Value[1]=AccBin[1][0]+AccBin[1][2]/2;
	while ( Acc_var_Value[1]<=AccBin[1][1] ){
	  /*Q2{{{*/
	  Acc_var_Value[2]=AccBin[2][0]+AccBin[2][2]/2;
	  while ( Acc_var_Value[2]<=AccBin[2][1] ){
	    /*Xbj{{{*/
	    Acc_var_Value[3]=AccBin[3][0]+AccBin[3][2]/2;
	    while ( Acc_var_Value[3]<=AccBin[3][1] ){
	      if ( NoCal%5==0 ){
		outlog<<Form(" +++ Calculating %d Points ...", NoCal)<<endl;
		cerr  <<Form(" +++ Calculating %d Points ...", NoCal)<<endl;
	      }
	      /*EX_Cuts{{{*/
	      EX_Cuts=Form("abs(%s-%f)<=%f&&abs(%s-%f)<=%f&&abs(%s-%f)<=%f&&abs(%s-%f)<=%f&&%s>=%f",
			   EX_Dp.Data(),
			   Acc_var_Value[0],
			   AccBin[0][2]/2,
			   EX_Nu.Data(),
			   Acc_var_Value[1],
			   AccBin[1][2]/2,
			   EX_Q2.Data(),
			   Acc_var_Value[2],
			   AccBin[2][2]/2,
			   EX_Xbj.Data(),
			   Acc_var_Value[3],
			   AccBin[3][2]/2,
			   EX_Current.Data(),
			   Beam_Trip_Cut
			   );
	      EX_Cuts+="&&"+ElectronCuts;
	      outlog <<"      Real Data Cuts in this bin: "<<endl
	       <<"      Real Data Cuts in this bin: "<<endl
		   <<"      ------------------------------------------------"<<endl
		   <<EX_Cuts.Data()<<endl
		   <<"      ------------------------------------------------"<<endl<<endl;
           cerr  <<"      Real Data Cuts in this bin: "<<endl
	       <<"      Real Data Cuts in this bin: "<<endl
		   <<"      ------------------------------------------------"<<endl
		   <<EX_Cuts.Data()<<endl
		   <<"      ------------------------------------------------"<<endl<<endl;
	      /*}}}*/
	      xbj = Acc_var_Value[3]; //We are binning on Xbj so use its bin-center values
	      Ep=Ebeam/(1.0+4.0*Ebeam*SinSQ/(2.0*PROTON_MASS*xbj)); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
	      qsq = 4.0*Ebeam*Ep*SinSQ;
	      nu = qsq/(2.0*PROTON_MASS*xbj);

	      outlog<<Form("******** xbj=%7.4f, Ep = %7.2f MeV/c",xbj,Ep)<<endl;
	      cerr  <<Form("******** xbj=%7.4f, Ep = %7.2f MeV/c",xbj,Ep)<<endl;

	      //ROOT 5.18 bug: have to add this after GetEntries(sel) fixed in 5.20
	      T_Tree->SetNotify(0);

	      /*Get_Nf_Bin{{{*/
	      //Calculate Number of Events in each bin, For SAMC data, weighted by cross section, for Real data, do only correction.
	      //1->Dp, 2->Nu, 3->Nu, 4->Xbj
	      SAMC_BinCut[0] =  AccBin[0][2]/2; SAMC_BinCut[1] =  AccBin[1][2]/2; SAMC_BinCut[2] =  AccBin[2][2]/2; SAMC_BinCut[3] =  AccBin[3][2]/2;
	      Nf_SAMC->Value = gGet_Nf_SAMC(SAMC_Tree,Acc_var_Value,SAMC_BinCut, SAMC_AccCut,Target_Name.Data(),Kin.Data(),Ep);
	      Nf_SAMC->Stat_Err=sqrt(Nf_SAMC->Value);
	      outlog << "      Total Number of SAMC Events = " << Nf_SAMC->Value<<endl<<endl;
	      cerr   << "      Total Number of SAMC Events = " << Nf_SAMC->Value<<endl<<endl;

	      Nf_EX->Value = gGet_Nf_EX(RunNoChain,Arm,EX_Cuts.Data(),LiveTimeChain, PSChain);
	      Nf_EX->Stat_Err=sqrt(Nf_EX->Value);
	      outlog << "      Total Number of EX Events = " << Nf_EX->Value<<endl;
	      cerr   << "      Total Number of EX Events = " << Nf_EX->Value<<endl;
	      /*}}}*/

	      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////																		  
		  if ( Nf_SAMC->Value>0 && Nf_EX->Value>0){
			  //Note: Live_Time is already corrected when calculating Nf_EX for each run#, Ntg is already in cm^-2			    
			  Cross_Section->Value = Nf_EX->Value*Total_Pion_Rej->Value/(Total_E_Eff->Value*Ntg->Value*Ne->Value)
				  *Nf_SAMC_Gen_All->Value/(Nf_SAMC->Value*SAMC_dOmega*(P_Center*SAMC_dDp));
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		map_err.clear();
		map_err.insert(make_pair(Nf_EX->Value,Nf_EX->Stat_Err));
		map_err.insert(make_pair(Nf_SAMC->Value,Nf_SAMC->Stat_Err));
		map_err.insert(make_pair(Nf_SAMC_Gen_All->Value,Nf_SAMC_Gen_All->Stat_Err));
		map_err.insert(make_pair(Ne->Value,Ne->Stat_Err));
		map_err.insert(make_pair(Ntg->Value,Ntg->Stat_Err));
		//		      map_err.insert(make_pair(Total_Pion_Rej->Value,Total_Pion_Rej->Stat_Err));
		//		      map_err.insert(make_pair(Total_E_Eff->Value,Total_E_Eff->Stat_Err));
		//		      map_err.insert(make_pair(Livetime->Value,Livetime->Stat_Err));
		//		      map_err.insert(make_pair(P0->Value,P0->Stat_Err));
		Cross_Section->Stat_Err=gCal_Simple_Error(map_err,Cross_Section->Value);

		Cross_Section->Sys_Err=0;
		outstring=Form("%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e",
			       16,Cross_Section->Value,
			       16,Cross_Section->Stat_Err,
			       16,Cross_Section->Sys_Err,
			       16,xbj, //thisbin[3],  //Xbj
			       16,Angle, //Deg
			       16,Ep, 
			       16,nu, //thisbin[1],  //Nu
			       16,qsq, //thisbin[2],  //Q2
			       16,Ebeam, //MeV
			       16,P_Center, //MeV/c
			       16,Livetime->Value,
			       16,Ntg->Value,
			       16,Ne->Value,
			       16,Nf_EX->Value,
			       16,Nf_SAMC->Value,
			       16,Nf_SAMC_Gen_All->Value);
		outfile<<outstring<<endl;
		outlog <<endl
		     <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		outlog << Form("      --> Xbj=%4.2f, Cross_Section=%e",Acc_var_Value[3], Cross_Section->Value)<<endl;
		outlog <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;   
	    cerr <<endl
		     <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		cerr << Form("      --> Xbj=%4.2f, Cross_Section=%e",Acc_var_Value[3], Cross_Section->Value)<<endl;
		cerr <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;   
	      }	// if ( Nf_SAMC->Value>0 && Nf_EX->Value>0){					      
	      NoCal++;
	      NumOfAccInput++;

	      Acc_var_Value[3]+=AccBin[3][2]; 
	    }/*}}}Xbj*/
	    Acc_var_Value[2]+=AccBin[2][2];
	  }  /*}}}Q2*/
	  Acc_var_Value[1]+=AccBin[1][2];
	} /*}}}Nu*/
	Acc_var_Value[0]+=AccBin[0][2];
      } /*}}}Dp*/
    }//if ( Total_Pion_Rej->Value>0&&Total_E_Eff->Value>0&&Ne->Value>0&&Ntg->Value>0&&Livetime->Value>0&&P0->Value>0 )
  }// If(bCal_XS)
  /*}}}*/

  /*free memory{{{*/
  outfile.close();
  outlog.close();
  delete Total_E_Eff;
  delete Total_Pion_Rej;
  delete Cross_Section;
  delete Nf_SAMC;
  delete Nf_EX;
  delete Nf_SAMC_Gen_All;
  delete LiveTimeChain;
  delete LiveTimeErrorChain;
  delete PSChain;
  delete T_Tree;
  /*}}}*/
  gSystem->Exec("chmod  -w "+outfile_name);
  log->WriteToLogFile(__FILE__,__LINE__,Form("%d Success. Result is saved at %s",RunNo,outfile_name.Data()));
  log->Destroy();
  return 0;
}
/*}}}*/

char* gInput;
char* gLog;
int getargs(int argc, char**argv);

/*int main{{{*/
int main(int argc, char**argv){
  int pserr = getargs(argc,argv);
  if(pserr<0){
    outlog<<"****** I can not recognize the parameters!"<<endl;
    cerr<<"****** I can not recognize the parameters!"<<endl;
    return -1;
  }

  TString InputFileName = "xs_input/xs_c12_kin3.1.dat"; //outlog <<"Input File Name: "<<endl; cin >> InputFileName; 
  InputFileName = gInput;
  TString LogFileName = "log/xs_c12_kin3.1.log";  // outlog <<"Log File Name: "<<endl; cin >> LogFileName;
  LogFileName = gLog;

  int err = CrossSection(InputFileName.Data(), LogFileName.Data());

  if(err>0)
    outlog<<"@@@@@@@@@ Code runs successfully!"<<endl;
    cerr  <<"@@@@@@@@@ Code runs successfully!"<<endl;
   
}
/*}}}*/
/*int getargs(int argc,char** argv){{{*/
int getargs(int argc,char** argv)
{
  char* argptr;
  bool noStop;  
  int err=-1;

  for(int i=0;i<argc;i++){
    argptr = argv[i];
    if(*argptr=='-'){
      argptr++;
      switch (*argptr){
      case 'I':
	if(*(++argptr))
	  gInput = argptr;
	err=0;
	break;
      case 'L':
	if(*(++argptr))
	  gLog = argptr;
	err=0;
	break;
      case 'h':
	cerr <<"================================================"<<endl;
	cerr <<"Option: " <<endl;
	cerr <<"-I[xs_c12_kin3.1.dat] Input File" <<endl;
	cerr <<"-L[c12_kin3.1.log] Log File"<<endl;
	cerr <<endl<<"================================================"<<endl;
	noStop=false;
	err=-1;
	goto OUT;
	break;
      default:
	cerr <<"================================================"<<endl;
	cerr <<"Option: " <<endl;
	cerr <<"-I[xs_c12_kin3.1.dat] Input File" <<endl;
	cerr <<"-L[c12_kin3.1.log] Log File"<<endl;
	cerr <<endl<<"================================================"<<endl;
	err=-1;
	break;
      }

    }
    noStop=true;
  OUT:if(!noStop){break;}
    continue;
  }
  return err; 
}
/*}}}*/		
