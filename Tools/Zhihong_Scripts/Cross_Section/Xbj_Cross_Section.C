/******************************************/
/*     X>2 Cross Section Calculation      */
/*         Zhihong Ye, 10/14/2011         */
/*         --Bin by Bin Correction        */
/******************************************/
/*Include{{{*/
#include "./SRC/XGT2_Main.h"
#include "./SRC/XGT2_Tools.h"
#include "./SRC/XGT2_LiveTime.h"
#include "./SRC/XGT2_Beam.h"
#include "./SRC/XGT2_Data.h"
#include "./SRC/XGT2_Cuts.h"
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
int CrossSection(const TString& aInputFile,const TString& aLogFile,const TString& aCommentName){

	TString aInputFileName = "./xs_input/"+aInputFile;
	TString aLogFileName = "./log/"+aLogFile+aCommentName;
	outlog.open("./log/"+aLogFile+aCommentName); 
	//Build a logger
	XGT2_Logger* log=XGT2_Logger::Instance();
	log->SetLogFileName(aLogFileName);
	Int_t j=0;

	/////////////////////////////
	// Initialize Input values
	/////////////////////////////
	//Set Value of Beam Energy
	const Double_t Ebeam = 3.356; //GeV
	/*Read InputFile{{{*/
	outlog << endl << Form("@@@@@@ Initialize Input values from %s ...",aInputFileName.Data())<<endl;
	cout   << endl << Form("@@@@@@ Initialize Input values from %s ...",aInputFileName.Data())<<endl;
	Int_t RunNo=-1;
	Int_t aline=0;
	RunNo=gGet_InputFile_Var(aInputFileName,aline).Atoi(); aline++;
	const TString Arm=gGet_InputFile_Var(aInputFileName,aline);aline++;                outlog << "   ARM = " << Arm <<endl;
	const Double_t P=gGet_InputFile_Var(aInputFileName,aline).Atof();aline++;          outlog << "   P0 = " << P <<endl;
	Double_t Angle=gGet_InputFile_Var(aInputFileName,aline).Atof();aline++;      outlog << "   Angle = " << Angle <<endl;
	const TString Kin=gGet_InputFile_Var(aInputFileName,aline);aline++;                outlog << "   Kin = " << Kin<<endl;
	const TString Target_Name=gGet_InputFile_Var(aInputFileName,aline);aline++;        outlog << "   Target Name = " << Target_Name<<endl;
	Double_t Target_Thickness=gGet_InputFile_Var(aInputFileName,aline).Atof();aline++; outlog << Form("   Target Thickness = %f g/cm^2",Target_Thickness)<<endl; //check the unit! g/cm^2
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
	int* RunChain = new int[Chain_Size];
	for(int i=0;i<Chain_Size;i++)
		RunChain[i] = RunNoChain[i];

	gCheckDirectory(Form("%s/results/%s",XGT2_DIR.Data(),gGetScriptname(__FILE__).Data()));
	Beam_Trip_Cut = gGet_InputFile_Var(aInputFileName,aline).Atof();aline++;   outlog << Form("   Beam Tip Cut = %6.2f uA",Beam_Trip_Cut)<<endl;
	const Bool_t bCheck_Data= gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++; outlog << "   Check Data? " << bCheck_Data<<endl;
	const Bool_t bCal_XS = gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++;    outlog << "   Calculate Cross Section? " << bCal_XS<<endl; 
	/*}}}end of reading input*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// End of Initializing Input values
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////
	// Number of Incident particles
	///////////////////////////////////
	/*Get Ne{{{*/
	outlog << endl << "@@@@@@ Calculating Incident Particle ..."<<endl;	
	cout   << endl << "@@@@@@ Calculating Incident Particle ..."<<endl;	
	Double_t runtime;
	Double_t* NeChain = new double[Chain_Size];
	const XGT2_VAR* Ne=gCal_Ne(RunNoChain,Arm,runtime, NeChain);
	outlog <<Form("     Ne = %4.3e", Ne->Value)<<endl;	
	cout   <<Form("     Ne = %4.3e", Ne->Value)<<endl;	
	log->WriteToLogFile(__FILE__,__LINE__,Form("     Ne = %f, Runtime = %f",Ne->Value, runtime));
	/*}}}end of Get Ne*/

	//////////////////////////
	// Check PreScalerR
	//////////////////////////
	/*Get PS{{{*/
	outlog << endl << "@@@@@@ Checking and Storing PreScaler Factors for every run ..." <<endl;	
	cout   << endl << "@@@@@@ Checking and Storing PreScaler Factors for every run ..." <<endl;	
	Int_t* PSChain = new int[Chain_Size*8];
	gCheck_PS(RunNoChain,PSChain);
	/*}}}end of PS*/

	//////////////////////////
	// LiveTime
	//////////////////////////
	/*LiveTime{{{*/
	outlog  << endl<< "@@@@@@ Checking and Storing Live Time for every run ..." <<endl; 
	cout    << endl<< "@@@@@@ Checking and Storing Live Time for every run ..." <<endl;
	Double_t* LiveTimeChain = new double[Chain_Size];
	Double_t* LiveTimeStatErrChain = new double[Chain_Size]; //Store the Stat Error of Livetime for total Error_Bar calculation;
	Double_t* LiveTimeSysErrChain = new double[Chain_Size];  //Store the Stat Error of Livetime for total Error_Bar calculation;
	for(int i=0;i<Chain_Size;i++){ 
		LiveTimeChain[i] = -1.0; //Initialize the livetime
		LiveTimeStatErrChain[i] = -1.0;
		LiveTimeSysErrChain[i] = -1.0;
	}
	//Obtain Livetime for each run from tables
	gLoad_LiveTime_Chain(RunNoChain,Arm,LiveTimeChain,LiveTimeStatErrChain,LiveTimeSysErrChain,2);//1->Beam Trip Cut, 2->NONE Cuts
	//Define a Average livetime for all runs, just for reference and checking used. LiveTime is corrected run by run.
	XGT2_VAR* Livetime = new XGT2_VAR();
	Livetime->Value = 0.0; Livetime->Stat_Err = 0.0; Livetime->Sys_Err = 0.0;
	for(int i=0;i<Chain_Size;i++){ 
		Livetime->Value    += LiveTimeChain[i];
		Livetime->Stat_Err += pow(LiveTimeStatErrChain[i]/LiveTimeChain[i],2);
		Livetime->Sys_Err  += pow(LiveTimeSysErrChain[i]/LiveTimeChain[i],2);
	}
	Livetime->Value/=Chain_Size; 
	Livetime->Stat_Err=Livetime->Value*sqrt(Livetime->Stat_Err);
	Livetime->Sys_Err=Livetime->Value*sqrt(Livetime->Sys_Err);

	outlog << "      Average Lifetime = " << Livetime->Value<<endl;
	cout   << "      Average Lifetime = " << Livetime->Value<<endl;
	/*}}}end of LiveTime*/

	//////////////////////
	// Check Current
	//////////////////////
	//Check and Print out values of currents from each Arm for each run.
	outlog  << endl<< "@@@@@@ Checking Currents for each run ..." <<endl; 
	gCheck_Current(RunNoChain);

	//////////////////////////
	// Check Data and Cuts
	//////////////////////////
	/*Check Data and Cuts{{{*/
	//	if(bCheck_Data){
	//	  outlog << endl << "@@@@@@ Checking Data and Cuts ..." <<endl;
	//	  gCheck_Data_Cuts(RunNoChain,Arm,ElectronCuts.Data(),IsExtTgt);
	//	}
	/*}}}end of Check*/

	//////////////
	// Chain Runs
	//////////////
	outlog << endl << "@@@@@@ Chain Runs ..."<<endl;
	cout   << endl << "@@@@@@ Chain Runs ..."<<endl;
	TChain* T_Tree=gGetTree(RunNoChain,Arm,T_Tree_Name);

	////////////////////////////////////
	// Calculate Cross Section
	////////////////////////////////////
	outlog <<endl<< "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl;
	outlog       << "oO0Oo              --- Calculating Cross Section ---              oO0Oo" <<endl;
	outlog       << "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl<<endl;

	cout <<endl<< "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl;
	cout       << "oO0Oo              --- Calculating Cross Section ---              oO0Oo" <<endl;
	cout       << "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl<<endl;

	/*Kinematics Settings{{{*/
	////////////////////////////////////////////
	//Check Beam Energy and Store all P0
	/////////////////////////////////////////////  
	outlog  << endl<< "@@@@@@ Checking Beam Energy for each run ..." <<endl; 
	cout    << endl<< "@@@@@@ Checking Beam Energy for each run ..." <<endl; 
	Double_t* E0Chain = new double[Chain_Size]; //GeV
	Double_t* dEEChain = new double[Chain_Size];//GeV
	gCheck_E0(RunNoChain,Arm,E_Tree_Name,E0Chain, dEEChain); //Check P0 for each run from the Dipole readings. 

	////////////////////////////////////////////
	//Check Central Momentum and Store all P0
	/////////////////////////////////////////////  
	outlog  << endl<< "@@@@@@ Checking Central Momentum for each run ..." <<endl; 
	cout    << endl<< "@@@@@@ Checking Central Momentum for each run ..." <<endl; 
	Double_t* P0Chain = new double[Chain_Size];//GeV/c
	gCheck_P0(RunNoChain,Arm,E_Tree_Name,P0Chain); //Check P0 for each run from the Dipole readings. 
	Double_t P0_Avg = 0.0;
	for(int i=0;i<Chain_Size;i++)
		P0_Avg+=P0Chain[i];
	P0_Avg/=Chain_Size;

	////////////////////////////////////////////
	//Check Scattering Angle and Store all values
	/////////////////////////////////////////////  
	outlog  << endl<< "@@@@@@ Checking Scattering Angle for each run ..." <<endl; 
	cout    << endl<< "@@@@@@ Checking Scattering Angle for each run ..." <<endl; 
	Double_t* AngleChain = new double[Chain_Size];//Degree
	gCheck_Angle(RunNoChain,Arm,E_Tree_Name,AngleChain); //Check P0 for each run from the Dipole readings. 
	Double_t Angle_Avg = 0.0;
	for(int i=0;i<Chain_Size;i++)
		Angle_Avg+=AngleChain[i];
	Angle_Avg/=Chain_Size;
	/*}}}*/

	//////////////////
	// Set Cuts
	//////////////////
	/*Set Cuts{{{*/
	outlog << endl << "@@@@@@ Obtaining Cuts ..." <<endl;
	cout   << endl << "@@@@@@ Obtaining Cuts ..." <<endl;
	TString ElectronCuts;
	Double_t Target_Length_Cut = 0.0; //in meter, Symmetry Cuts, like |VZ|<Target_Length_Cut
	TString AccCut_FileName = Form("%s/Acceptance/Acc_%s_Table.dat_%s",XGT2_DIR.Data(), Arm.Data(),aCommentName.Data());
	////If there is not a specific acceptance file, use the default one
	if (!(gSystem->FindFile(Form("%s/Acceptance/",XGT2_DIR.Data()),AccCut_FileName))){
	        cout<<endl<<"******* No acceptance file is found. Use the default one!!"<<endl;
                AccCut_FileName = Form("%s/Acceptance/Acc_%s_Table.dat",XGT2_DIR.Data(), Arm.Data());
        }
        cout<<endl<<Form("        The Acceptance file is: %s", AccCut_FileName.Data())<<endl;

	outlog <<"     Getting Electron Cuts for "<<Target_Name <<" Run = "<<RunNo <<" on "<<Arm<<" Arm, from "<< AccCut_FileName.Data() <<endl; 
	cout   <<"     Getting Electron Cuts for "<<Target_Name <<" Run = "<<RunNo <<" on "<<Arm<<" Arm, from "<< AccCut_FileName.Data() <<endl; 
	ElectronCuts = gGet_XGT2_Cut(RunNo,Arm,Target_Name.Data(),Kin.Data(),Angle_Avg,P,Main_Trigger,AccCut_FileName,IsExtTgt, &Target_Length_Cut); 
	outlog <<"--------------------------------------------------------------------------------------"<<endl;
	cout <<"--------------------------------------------------------------------------------------"<<endl;
	outlog <<"     --> "<<ElectronCuts.Data()<<endl;
	cout   <<"     --> "<<ElectronCuts.Data()<<endl;
	outlog <<"--------------------------------------------------------------------------------------"<<endl;
	cout   <<"--------------------------------------------------------------------------------------"<<endl;
	
	TString ElectronCuts_AL="1";
	if(Target_Name == "H2"||Target_Name == "D2"||Target_Name =="He3" || Target_Name =="He4"){
		ElectronCuts_AL = gGet_XGT2_Cut(RunNo,Arm,"Dummy",Kin.Data(),Angle_Avg,P,Main_Trigger,AccCut_FileName,IsExtTgt, &Target_Length_Cut); 
		outlog <<"-------------------Dummy Cuts:--------------------------------------------------------"<<endl;
		cout <<"----------Dummy Cuts:-------------------------------------------------------------------"<<endl;
		outlog <<"     --> "<<ElectronCuts_AL.Data()<<endl;
		cout   <<"     --> "<<ElectronCuts_AL.Data()<<endl;
		outlog <<"--------------------------------------------------------------------------------------"<<endl;
		cout   <<"--------------------------------------------------------------------------------------"<<endl;
	}

	/*}}}end of Cuts*/

	/////////////////////
	// Define variables
	/////////////////////
	/*Define variables{{{*/
	XGT2_VAR* Total_E_Eff=new XGT2_VAR();
	XGT2_VAR* Total_Pion_Rej=new XGT2_VAR();
	//For X>2 Experiment, PID efficiencies are all above 99%, so we ignore the values and put them into sysmetic errors
	Total_E_Eff->SetValue(1,0,0.01); 
	Total_Pion_Rej->SetValue(1,0,0.01);
	Double_t nu=-1000., xbj=-1000., bin_size = -1000.0, qsq=-1000., Ep=-1000.0; 
	Double_t XS_Born_XEMC = 0.0;
	Double_t BinCor = 0.0;
	XGT2_VAR* Cross_Section=new XGT2_VAR();
	XGT2_VAR* XS_Avg=new XGT2_VAR();
	XGT2_VAR* Nf_EX=new XGT2_VAR();  //Total Number of real events in each bin after all cuts
	XGT2_VAR* Nf_AL=new XGT2_VAR();  //Total Number of real events for Dummy in each bin after all cuts
	XGT2_VAR* Yield_EX=new XGT2_VAR();    //Experimental Yield in this Bin
	XGT2_VAR* Yield_AL=new XGT2_VAR();    //Experimental Dummy Yield in this Bin
	XGT2_VAR* Yield_MC=new XGT2_VAR();  //SAMC Yield in this Bin
	XGT2_VAR* Yield_EX_Avg=new XGT2_VAR();    //Experimental Yield in this Bin
	XGT2_VAR* Yield_MC_Avg=new XGT2_VAR();  //SAMC Yield in this Bin
	XGT2_VAR* XS_SAMC_SUM=new XGT2_VAR(); //Sum of cross section values for all SAMC events in one bin

	//Double_t rms;//tmp use

	TString OutFile_Dir = Form("%s/results/%s/%s",
			XGT2_DIR.Data(),Target_Name.Data(),aCommentName.Data()); 
	if(!gSystem->OpenDirectory(OutFile_Dir.Data()))
		gSystem->mkdir(OutFile_Dir);	
	TString outfile_name=Form("%s/results/%s/%s/%s_%s_%s_Yield.dat_%s",
			XGT2_DIR.Data(),Target_Name.Data(),aCommentName.Data(),Arm.Data(),Kin.Data(),Target_Name.Data(),aCommentName.Data());
	outlog<<"      Backup and Save Output file to "<<outfile_name.Data()<<endl;
	cout  <<"      Backup and Save Output file to "<<outfile_name.Data()<<endl;
	gSystem->Exec("cp -vpf "+outfile_name+" "+outfile_name+".old");
	gSystem->Exec("chmod -w "+outfile_name+".old");
	gSystem->Exec("chmod +w "+outfile_name);

	ofstream outfile(outfile_name.Data());
	TString outstring=Form("%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s",
			16,"XS(nb/MeV/sr)",
			16,"XS_Stat_Err",
			16,"XS_Sys_Err",
			16,"XS_MC",
			16,"Xbj",
			16,"Angle",
			16,"Ep(GeV/c)",
			16,"Nu",
			16,"Q2",
			16,"BinCor",
			16,"LiveTime",
			16,"Ntg",
			16,"Ne",
			16,"Nf_EX",
			16,"Nf_AL",
			16,"XS_SAMC_SUM",
			16,"Nf_SAMC_Gen",
			16,"Yield_EX",
			16,"Y_EX_Stat_Err",
			16,"Y_EX_Sys_Err",
			16,"Yield_MC",
			16,"Y_MC_Stat_Err",
			16,"Y_MC_Sys_Err",
			16,"Yield_AL",
			16,"Y_AL_Stat_Err",
			16,"Y_AL_Sys_Err"
					);
	outfile<<outstring<<endl;

	//Double_t SinSQ = pow(sin(Angle_Avg*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2
	//I should use the setting Angle instead the real value since I calcualte the cross section with this angle
	Double_t SinSQ = pow(sin(Angle*TMath::DegToRad()/2.0),2); //Define Sin(Theta/2.)^2
	/*}}}*/

	////////////////////////////////
	// Target Luminocity
	////////////////////////////////
	/*Get Ntg{{{*/
	outlog << endl << "@@@@@@ Calculating Target Luminocity ..."<<endl;
	cout   << endl << "@@@@@@ Calculating Target Luminocity ..."<<endl;

	XGT2_VAR* Ntg = new XGT2_VAR();
	Double_t* NtgChain = new double[Chain_Size];
	if(IsExtTgt){
		Target_Length_Cut*=100.0; //From meter to centi-meter
		cout<<Form("      Cutting on Extended Long Target <= %5.4f cm.", Target_Length_Cut)<<endl;
		//Boiling Effect for long cryo-targets will be corrected in this subroutine
		Ntg= gCal_Ntg(RunNoChain,Arm.Data(),Target_A,Target_Length_Cut, Target_Thickness,Target_Thickness_Stat_Err,NtgChain);
	}
	else{
		//Foil targets are not corrected by boiling effect
		Ntg->Value=Target_Thickness/Target_A*Na; // number of target particles in unit area /cm^2
		Ntg->Stat_Err=Target_Thickness_Stat_Err*Na/Target_A;
		for(int i=0;i<Chain_Size;i++)
			NtgChain[i]=Ntg->Value;
	}
	outlog <<Form("     Avg Ntg = %4.3e; Before Boiling Corrected = %4.3e",Ntg->Value,Target_Thickness/Target_A*Na)<<endl;
	cout   <<Form("     Avg Ntg = %4.3e; Before Boiling Corrected = %4.3e",Ntg->Value,Target_Thickness/Target_A*Na)<<endl;
	log->WriteToLogFile(__FILE__,__LINE__,Form("     Ntg = %f, No Boiling Corrected = %f",Ntg->Value,Target_Thickness/Target_A*Na));
	/*}}}end of Get Ntg*/

	//////////////////////////////
	// Get Bins
	//////////////////////////////
	/*Get Bins{{{*/
	const TString Bin_Variable = "xbj";//Bin on Xbj
	//const TString Bin_Variable = "Ep";//Bin on Ep
	Double_t AccBin[3], Acc_var_Value;
	outlog <<"++++++ Getting Acceptance Binning from File "<<AccBin_File.Data()<<endl;
	cout   <<"++++++ Getting Acceptance Binning from File "<<AccBin_File.Data()<<endl;
	gGet_Acc_Bin(AccBin_File.Data(),AccBin,Target_Name.Data(),Bin_Variable.Data());
	for ( j=0; j<3; j++){
		outlog << Form("      AccBin[%d]=%f",j,AccBin[j])<<endl;
		cout   << Form("      AccBin[%d]=%f",j,AccBin[j])<<endl;
	}
	if(Target_Name.Contains("He3")&&AccBin[1]>3.0)
		AccBin[1]=3.2; //Heavier target x->4.0, check few more out of range points for error study
	if(Target_Name.Contains("H2")&&AccBin[1]>2.0)
		AccBin[1]=2.2; //Heavier target x->4.0, check few more out of range points for error study

	Int_t TotalBin = TMath::Nint((AccBin[1]-AccBin[0])/AccBin[2]);

	outlog <<endl<< "      Total Number of Bins = " << TotalBin<<endl<<endl;
	cout   <<endl<< "      Total Number of Bins = " << TotalBin<<endl<<endl;
	/*}}}*/

	Int_t NumOfAccInput=0, NoCal=0;
	map<Double_t,Double_t> map_err;

	/////////////
	// SAMC
	/////////////	
	/*Get Info From SAMC Simulation{{{*/
	outlog <<endl<<"++++++ Info in SAMC:"<<endl;
	cout   <<endl<<"++++++ Info in SAMC:"<<endl;
	Double_t SAMC_dDp    = 0.12;    outlog << Form("       Full Range of Dp    in SAMC = %f PCT",  SAMC_dDp)<<endl;     //% 
	Double_t SAMC_dTheta = 0.18;    outlog << Form("       Full Range of Theta in SAMC = %f Rad",  SAMC_dTheta)<<endl;  // Rad  
	Double_t SAMC_dPhi   = 0.09;    outlog << Form("       Full Range of Phi   in SAMC = %f Rad",  SAMC_dPhi)<<endl;    // Rad
	Double_t SAMC_dOmega = SAMC_dTheta*SAMC_dPhi;
	Double_t SAMC_BinCut; 

	//Chain Data
	TChain* SAMC_Tree=new TChain("SAMC"); 
	TString SAMC_RootFile = SAMCFILES_DIR+SAMC_File;
	SAMC_Tree->Add(SAMC_RootFile.Data());
	outlog<<Form("       Added SAMC Rootfile %s.",SAMC_RootFile.Data())<<endl;
	cout  <<Form("       Added SAMC Rootfile %s.",SAMC_RootFile.Data())<<endl;
	XGT2_VAR* Nf_SAMC_Gen_All=new XGT2_VAR();  //Total Number of simulation events in the file, without any cuts
        Nf_SAMC_Gen_All->Value = SAMC_Tree->GetEntries();
        Nf_SAMC_Gen_All->Stat_Err = 0.0;
	Nf_SAMC_Gen_All->Sys_Err = sqrt(Nf_SAMC_Gen_All->Value);//Set it to be the Sys Error, 01/13/2015
	outlog << Form("       Total Events generated in SAMC = %f ",Nf_SAMC_Gen_All->Value)<<endl;
	cout   << Form("       Total Events generated in SAMC = %f ",Nf_SAMC_Gen_All->Value)<<endl;
	outlog<<endl;
	cout  <<endl;

	//Get SAMC Accptance Cuts. Read the values from the same Acc File as for Data
	Double_t SAMC_AccCut[9]; //0->x_fp,1->y_fp,2->th_fp,3->ph_fp,4->dp_tg,5->y_tg,6->th_tg,7->ph_tg,8->reactz
	gGet_XGT2_Acc_Cut(Arm, Target_Name, P, SAMC_AccCut,AccCut_FileName.Data());
	/*}}}*/

	////////////////
	// Set Tree
	////////////////
	/*Set Tree{{{*/
	TString Root_Out = Form("%s/results/%s/%s/%s_%s_%s_Yield_%s.root",
			XGT2_DIR.Data(),Target_Name.Data(),aCommentName.Data(),Arm.Data(),Kin.Data(),Target_Name.Data(),aCommentName.Data());
	gSystem->Exec("cp -vpf "+Root_Out+" "+Root_Out+".old");
	gSystem->Exec("chmod -w "+Root_Out+".old");
	gSystem->Exec("chmod +w "+Root_Out);

	Double_t Nf_EX_Chain[Chain_Size];	
	Double_t Nf_EX_StatChain[Chain_Size];	
	Double_t Nf_EX_SysChain[Chain_Size];	
	Double_t Y_EX_Chain[Chain_Size];	
	Double_t Y_EX_StatChain[Chain_Size];	
	Double_t Y_EX_SysChain[Chain_Size];	
	Double_t Y_MC_Chain[Chain_Size];//For Cryo-Target with different boilng effect for each run.	
	Double_t Y_MC_StatChain[Chain_Size];	
	Double_t Y_MC_SysChain[Chain_Size];	
	Double_t XS_Chain[Chain_Size];	
	Double_t XS_StatChain[Chain_Size];	
	Double_t XS_SysChain[Chain_Size];	

	TFile *fxs = new TFile(Root_Out.Data(),"recreate");
	TTree *T = new TTree("T","A Tree for XS");

	T->Branch("Ntg",             &Ntg->Value,             "Ntg/D");
	T->Branch("Nf_EX",           &Nf_EX->Value,           "Nf_EX/D");
	T->Branch("Nf_EX_Stat",      &Nf_EX->Stat_Err,        "Nf_EX_Stat/D");
	T->Branch("Nf_EX_Sys",       &Nf_EX->Sys_Err,         "Nf_EX_Sys/D");
	T->Branch("Y_EX",            &Yield_EX->Value,        "Y_EX/D");
	T->Branch("Y_EX_Stat",       &Yield_EX->Stat_Err,     "Y_EX_Stat/D");
	T->Branch("Y_EX_Sys",        &Yield_EX->Sys_Err,      "Y_EX_Sys/D");
	T->Branch("Y_AL",            &Yield_AL->Value,        "Y_AL/D");
	T->Branch("Y_AL_Stat",       &Yield_AL->Stat_Err,     "Y_AL_Stat/D");
	T->Branch("Y_AL_Sys",        &Yield_AL->Sys_Err,      "Y_AL_Sys/D");
        T->Branch("Y_EX_Avg",        &Yield_EX_Avg->Value,    "Y_EX_Avg/D");
        T->Branch("Y_EX_AvgStat",    &Yield_EX_Avg->Stat_Err, "Y_EX_AvgStat/D");
	T->Branch("Y_EX_AvgSys",     &Yield_EX_Avg->Sys_Err,  "Y_EX_AvgSys/D");
	T->Branch("XS",              &Cross_Section->Value,   "XS/D");
	T->Branch("XS_Stat",         &Cross_Section->Stat_Err,"XS_Stat/D");
	T->Branch("XS_Sys",          &Cross_Section->Sys_Err, "XS_Sys/D");
	T->Branch("XS_Avg",          &XS_Avg->Value,   "XS_Avg/D");
	T->Branch("XS_Avg_Stat",     &XS_Avg->Stat_Err,"XS_Avg_Stat/D");
	T->Branch("XS_Avg_Sys",      &XS_Avg->Sys_Err, "XS_Avg_Sys/D");

	T->Branch("Nf_MC_All",       &Nf_SAMC_Gen_All->Value, "Nf_MC_All/D");
	T->Branch("Y_MC",            &Yield_MC->Value,      "Y_MC/D");
	T->Branch("Y_MC_Stat",       &Yield_MC->Stat_Err,   "Y_MC_Stat/D");
	T->Branch("Y_MC_Sys",        &Yield_MC->Sys_Err,    "Y_MC_Sys/D");
	T->Branch("Y_MC_Avg",        &Yield_MC_Avg->Value,      "Y_MC_Avg/D");
	T->Branch("Y_MC_AvgStat",    &Yield_MC_Avg->Stat_Err,   "Y_MC_AvgStat/D");
	T->Branch("Y_MC_AvgSys",     &Yield_MC_Avg->Sys_Err,    "Y_MC_AvgSys/D");
	T->Branch("XS_MC",           &XS_SAMC_SUM->Value,     "XS_MC_SUM/D");
	T->Branch("XS_MC_Stat",      &XS_SAMC_SUM->Stat_Err,  "XS_MC_Stat/D");
	T->Branch("XS_MC_Sys",       &XS_SAMC_SUM->Sys_Err,   "XS_MC_Sys/D");
	T->Branch("XS_Born_MC",      &XS_Born_XEMC,           "XS_Born_MC/D");

	T->Branch(    "Nf_EX_Chain",     Nf_EX_Chain,     Form("Nf_EX_Chain[%d]/D",Chain_Size));
	T->Branch("Nf_EX_StatChain", Nf_EX_StatChain, Form("Nf_EX_StatChain[%d]/D",Chain_Size));
	T->Branch( "Nf_EX_SysChain",  Nf_EX_SysChain,  Form("Nf_EX_SysChain[%d]/D",Chain_Size));
	T->Branch(    "Y_EX_Chain",      Y_EX_Chain,      Form("Y_EX_Chain[%d]/D",Chain_Size));
	T->Branch("Y_EX_StatChain",  Y_EX_StatChain,  Form("Y_EX_StatChain[%d]/D",Chain_Size));
	T->Branch( "Y_EX_SysChain",   Y_EX_SysChain,   Form("Y_EX_SysChain[%d]/D",Chain_Size));
	T->Branch(    "Y_MC_Chain",      Y_MC_Chain,      Form("Y_MC_Chain[%d]/D",Chain_Size));
	T->Branch("Y_MC_StatChain",  Y_MC_StatChain,  Form("Y_MC_StatChain[%d]/D",Chain_Size));
	T->Branch( "Y_MC_SysChain",   Y_MC_SysChain,   Form("Y_MC_SysChain[%d]/D",Chain_Size));
	T->Branch(    "XS_Chain",        XS_Chain,        Form("XS_Chain[%d]/D",Chain_Size));
	T->Branch("XS_StatChain",    XS_StatChain,    Form("XS_StatChain[%d]/D",Chain_Size));
	T->Branch( "XS_SysChain",     XS_SysChain,     Form("XS_SysChain[%d]/D",Chain_Size));

	for(int i=0; i<Chain_Size; i++){
		T->Branch(Form("Run_%d",    RunChain[i]), &RunChain[i],            Form("Run_%d/I",    RunChain[i]));
		T->Branch(Form("Ne_%d",     RunChain[i]), &NeChain[i],             Form("Ne_%d/D",     RunChain[i]));
		T->Branch(Form("Ntg_%d",     RunChain[i]),&NtgChain[i],            Form("Ntg_%d/D",    RunChain[i]));
		T->Branch(Form("E0_%d",     RunChain[i]), &E0Chain[i],             Form("E0_%d/D",     RunChain[i]));
		T->Branch(Form("dEE_%d",    RunChain[i]), &dEEChain[i],            Form("dEE_%d/D",    RunChain[i]));
		T->Branch(Form("P0_%d",     RunChain[i]), &P0Chain[i],             Form("P0_%d/D",     RunChain[i]));
		T->Branch(Form("Angle_%d",  RunChain[i]), &AngleChain[i],          Form("Angle_%d/D",  RunChain[i]));
		T->Branch(Form("LT_%d",     RunChain[i]), &LiveTimeChain[i],       Form("LT_%d/D",     RunChain[i]));
		T->Branch(Form("LT_Stat_%d",RunChain[i]), &LiveTimeStatErrChain[i],Form("LT_Stat_%d/D",RunChain[i]));
		T->Branch(Form("LT_Sys_%d", RunChain[i]), &LiveTimeSysErrChain[i], Form("LT_Sys_%d/D", RunChain[i]));

		T->Branch(Form("Nf_EX_%d",     RunChain[i]), &Nf_EX_Chain[i],     Form("Nf_EX_%d/D",     RunChain[i]));
		T->Branch(Form("Nf_EX_Stat_%d",RunChain[i]), &Nf_EX_StatChain[i], Form("Nf_EX_Stat_%d/D",RunChain[i]));
		T->Branch(Form("Nf_EX_Sys_%d", RunChain[i]), &Nf_EX_SysChain[i],  Form("Nf_EX_Sys_%d/D", RunChain[i]));
		T->Branch(Form("Y_EX_%d",      RunChain[i]), &Y_EX_Chain[i],      Form("Y_EX_%d/D",      RunChain[i]));
		T->Branch(Form("Y_EX_Stat_%d", RunChain[i]), &Y_EX_StatChain[i],  Form("Y_EX_Stat_%d/D", RunChain[i]));
		T->Branch(Form("Y_EX_Sys_%d",  RunChain[i]), &Y_EX_SysChain[i],   Form("Y_EX_Sys_%d/D",  RunChain[i]));
		T->Branch(Form("Y_MC_%d",      RunChain[i]), &Y_MC_Chain[i],      Form("Y_MC_%d/D",      RunChain[i]));
		T->Branch(Form("Y_MC_Stat_%d", RunChain[i]), &Y_MC_StatChain[i],  Form("Y_MC_Stat_%d/D", RunChain[i]));
		T->Branch(Form("Y_MC_Sys_%d",  RunChain[i]), &Y_MC_SysChain[i],   Form("Y_MC_Sys_%d/D",  RunChain[i]));
		T->Branch(Form("XS_%d",        RunChain[i]), &XS_Chain[i],        Form("XS_%d/D",        RunChain[i]));
		T->Branch(Form("XS_Stat_%d",   RunChain[i]), &XS_StatChain[i],    Form("XS_Stat_%d/D",   RunChain[i]));
		T->Branch(Form("XS_Sys_%d",    RunChain[i]), &XS_SysChain[i],     Form("XS_Sys_%d/D",    RunChain[i]));
	}

	int ChainSize = Chain_Size;
	T->Branch("ChainSize", &ChainSize, "ChainSize/I");

	T->Branch("xbj", &xbj, "xbj/D");
	T->Branch("nu",  &nu,  "nu/D");
	T->Branch("qsq", &qsq, "qsq/D");
	T->Branch("Ep",  &Ep,  "Ep/D");
	//T->Branch("P",&P,"P/D");
	//T->Branch("Ebeam",&Ebeam,"Ebeam/D");
	T->Branch("Angle",&Angle,"Angle/D");
	T->Branch("BinCor",&BinCor,"BinCor/D");
	/*}}}*/

	///////////////////////////////////
	// Calculate XS and Save
	///////////////////////////////////
	const Double_t P_MeV = P0_Avg*GeVToMeV; //From GeV to MeV during XS calculation
	double Jump = 1;
	/*Calculate XS and Save{{{*/
	if ( Total_Pion_Rej->Value>0&&Total_E_Eff->Value>0&&Ne->Value>0&&Ntg->Value>0){
		/*Xbj{{{*/
		if(xbj>=1.6&&xbj<2.4)
			Jump=3.0;
		else if(xbj>=2.4)
			Jump=4.0;
		else
			Jump=1.5;
		bin_size =	Jump*AccBin[2]; //the Xbj bin size 
		Acc_var_Value=AccBin[0]+bin_size/2.; //--> Xbj

		//Cut on the bin that specified by Bin_Type
		while ( Acc_var_Value<=AccBin[1] ){
			if ( NoCal%5==0 ){
				outlog<<Form(" +++ Calculating %d Points ...", NoCal)<<endl;
				cout  <<Form(" +++ Calculating %d Points ...", NoCal)<<endl;
			}

			xbj = Acc_var_Value; //We are binning on Xbj so use its bin-center values
			Ep=Ebeam/(1.0+4.0*Ebeam*SinSQ/(2.0*PROTON_MASS*xbj)); //use xbj = Q2/2Mp/(Ei-Ep),Q2 = 4Ei*Ep*Sin_SQ 
			qsq = 4.0*Ebeam*Ep*SinSQ;
			nu = qsq/(2.0*PROTON_MASS*xbj);

			outlog<<Form("******** xbj=%7.4f, Ep = %7.5f GeV/c",xbj,Ep)<<endl;
			cout  <<Form("******** xbj=%7.4f, Ep = %7.5f GeV/c",xbj,Ep)<<endl;

			//ROOT 5.18 bug: have to add this after GetEntries(sel) fixed in 5.20
			T_Tree->SetNotify(0);
			//COME_HERE -- Just a mark to quick access here.

			/*Get_Nf_Bin{{{*/ 
			//Calculate Number of Events in each bin, For SAMC data, weighted by cross section, for Real data, do only correction.
			Nf_EX = gGet_Nf_EX(RunNoChain,Arm,Kin,Target_Name,xbj, bin_size/2.,ElectronCuts.Data(),Bin_Variable.Data(),LiveTimeChain,LiveTimeStatErrChain,LiveTimeSysErrChain, PSChain,E0Chain, P0Chain, AngleChain,IsExtTgt, Nf_EX_Chain,Nf_EX_StatChain,Nf_EX_SysChain);
			//Nf_EX = gGet_Nf_EX_Xbj(RunNoChain,Arm,Kin,Target_Name,xbj, bin_size/2.,ElectronCuts.Data(),Bin_Variable.Data(),LiveTimeChain,LiveTimeStatErrChain,LiveTimeSysErrChain, PSChain,E0Chain, P0Chain, AngleChain,IsExtTgt, Nf_EX_Chain,Nf_EX_StatChain,Nf_EX_SysChain);

			if(Target_Name == "H2"||Target_Name == "D2"||Target_Name =="He3" || Target_Name =="He4")
				Nf_AL = gGet_Nf_AL(Kin, Arm,xbj, bin_size/2.,ElectronCuts_AL.Data(),Bin_Variable.Data(),IsExtTgt);
			else{
			   Nf_AL->Value   = 1e-33;
			   Nf_AL->Sys_Err = 1e-33;
			   Nf_AL->Stat_Err= 1e-33;
			}

			Nf_AL->Value   = 1e-33;
			Nf_AL->Sys_Err = 1e-33;
			Nf_AL->Stat_Err= 1e-33;

			outlog << Form("      Total EX Events = %5.0f, AL Event = %5.0f", Nf_EX->Value, Nf_AL->Value)<<endl;
			cout   << Form("      Total EX Events = %5.0f, AL Event = %5.0f", Nf_EX->Value, Nf_AL->Value)<<endl;
			/*}}}*/

			/*Yield{{{*/
			/*For Experiment Data{{{*/
			if(Nf_EX->Value>5. && Ne->Value>5. && Total_Pion_Rej->Value > 1e-33 && Total_E_Eff->Value > 1e-33){
				Yield_EX->Value = Nf_EX->Value*Total_Pion_Rej->Value/(Ne->Value*Total_E_Eff->Value) * CM2ToNB; //Eq(4.9) in page118 of Nadia's Thesis
				Yield_EX->Stat_Err = Yield_EX->Value 
					* sqrt( pow( Nf_EX->Stat_Err / Nf_EX->Value, 2) 
							+ pow( Ne->Stat_Err / Ne->Value, 2) 
							+ pow( Total_Pion_Rej->Stat_Err / Total_Pion_Rej->Value, 2) 
							+ pow( Total_E_Eff->Stat_Err / Total_E_Eff->Value, 2)) ; 
				Yield_EX->Sys_Err = Yield_EX->Value 
					* sqrt( pow( Nf_EX->Sys_Err / Nf_EX->Value, 2) 
							+ pow( Ne->Sys_Err / Ne->Value, 2) 
							+ pow( Total_Pion_Rej->Sys_Err / Total_Pion_Rej->Value, 2) 
							+ pow( Total_E_Eff->Sys_Err / Total_E_Eff->Value, 2)) ; 

				if(Nf_AL->Value > 5.0){
					Yield_AL->Value = Nf_AL->Value*Total_Pion_Rej->Value/(Ne->Value*Total_E_Eff->Value) * CM2ToNB; //Eq(4.9) in page118 of Nadia's Thesis
					Yield_AL->Stat_Err = Yield_AL->Value 
						* sqrt( pow( Nf_AL->Stat_Err / Nf_AL->Value, 2) 
								+ pow( Ne->Stat_Err / Ne->Value, 2) 
								+ pow( Total_Pion_Rej->Stat_Err / Total_Pion_Rej->Value, 2) 
								+ pow( Total_E_Eff->Stat_Err / Total_E_Eff->Value, 2)) ; 
					Yield_AL->Sys_Err = Yield_AL->Value 
						* sqrt( pow( Nf_AL->Sys_Err / Nf_AL->Value, 2) 
								+ pow( Ne->Sys_Err / Ne->Value, 2) 
								+ pow( Total_Pion_Rej->Sys_Err / Total_Pion_Rej->Value, 2) 
								+ pow( Total_E_Eff->Sys_Err / Total_E_Eff->Value, 2)) ; 
				}
				else{
					Yield_AL->Value   = 1e-30;	
					Yield_AL->Stat_Err= 1e-30;	
					Yield_AL->Sys_Err = 1e-30;	
				}

				Yield_EX_Avg->Value   = 1e-30;
				Yield_EX_Avg->Stat_Err= 1e-30;
				Yield_EX_Avg->Sys_Err = 1e-30;
				int N_Y_EX = 0;
				for(int i=0;i<Chain_Size;i++){
					if(Nf_EX_Chain[i]>1){
						Y_EX_Chain[i] = Nf_EX_Chain[i]*Total_Pion_Rej->Value/(NeChain[i]*Total_E_Eff->Value) * CM2ToNB;
						Y_EX_StatChain[i] = Y_EX_Chain[i] 
							* sqrt( pow( Nf_EX_StatChain[i] / Nf_EX_Chain[i], 2) 
									+ 1.0/ NeChain[i] 
									+ pow( Total_Pion_Rej->Stat_Err / Total_Pion_Rej->Value, 2) 
									+ pow( Total_E_Eff->Stat_Err / Total_E_Eff->Value, 2)) ; 
						Y_EX_SysChain[i] = Y_EX_Chain[i] 
							* sqrt( pow( Nf_EX_SysChain[i] / Nf_EX_Chain[i], 2) 
									+ pow( Total_Pion_Rej->Sys_Err / Total_Pion_Rej->Value, 2) 
									+ pow( Total_E_Eff->Sys_Err / Total_E_Eff->Value, 2)) ; 
					}
					else{
						Y_EX_Chain[i]    = 1e-30;
						Y_EX_StatChain[i]= 1e-30;
						Y_EX_SysChain[i] = 1e-30;
					}

					if( Nf_EX_Chain[i]>1 && Y_EX_Chain[i]>1. && (Y_EX_StatChain[i]/Y_EX_Chain[i])<1.0 ){
						N_Y_EX ++;
						Yield_EX_Avg->Value    += Y_EX_Chain[i];
						Yield_EX_Avg->Stat_Err += pow( Y_EX_StatChain[i], 2 ); 
						Yield_EX_Avg->Sys_Err  += pow( Y_EX_SysChain[i] , 2 ); 
					}
					cout  <<Form("   ---Run#%d, Yield_EX = %7.2e", RunChain[i],Y_EX_Chain[i])<<endl;
					outlog<<Form("   ---Run#%d, Yield_EX = %7.2e", RunChain[i],Y_EX_Chain[i])<<endl;
				}
				if(N_Y_EX>0){
					Yield_EX_Avg->Value   /= N_Y_EX;
					Yield_EX_Avg->Stat_Err = sqrt( Yield_EX_Avg->Stat_Err/ N_Y_EX);
					Yield_EX_Avg->Sys_Err  = sqrt( Yield_EX_Avg->Sys_Err / N_Y_EX);
				}
			}
			else
				Yield_EX->SetValue(-1.0,1e-30,1e-30);
            /*}}}Y_EX*/

			/*For SAMC Data{{{*/
			Yield_MC->SetValue(-1.0,1e-30,1e-30);
			if(bCal_XS){ //Only calculate Y_MC when calculating CS
				if(Nf_EX->Value>1e-33 && Ne->Value>1e-33)
					XS_SAMC_SUM = gGet_XS_SAMC_SUM(SAMC_Tree,xbj,bin_size/2.,SAMC_AccCut,Arm.Data(), Target_Name.Data(),Kin.Data(),P0_Avg,Angle_Avg,Bin_Variable.Data(),&BinCor);
				else
					XS_SAMC_SUM->SetValue(-1.0,1e-30,1e-30);

				if(Ntg->Value>1. && XS_SAMC_SUM->Value>1e-33 && Nf_SAMC_Gen_All->Value>1.){	

					cout<<Form("-------------------------------------------------")<<endl;
					outlog<<Form("-------------------------------------------------")<<endl;
					cout<<Form(" P=%e, XS_Sum=%7.2e, ", P_MeV, XS_SAMC_SUM->Value )<<endl;
				    outlog<<Form(" P=%e, XS_Sum=%7.2e, ", P_MeV, XS_SAMC_SUM->Value )<<endl;
					cout<<Form("-------------------------------------------------")<<endl;
					outlog<<Form("-------------------------------------------------")<<endl;
					
					Yield_MC->Value =  Ntg->Value*XS_SAMC_SUM->Value*(SAMC_dOmega*(P_MeV*SAMC_dDp)/Nf_SAMC_Gen_All->Value); //Eq(4.18) in page 130 of Nadia's Thesis,set Det Effi to one and quote its sys err later
					Yield_MC->Stat_Err = Yield_MC->Value 
						* sqrt( pow( Ntg->Stat_Err / Ntg->Value, 2)
								+ pow( XS_SAMC_SUM->Stat_Err / XS_SAMC_SUM->Value, 2) 
								+ pow( Nf_SAMC_Gen_All->Stat_Err / Nf_SAMC_Gen_All->Value, 2) );
					Yield_MC->Sys_Err = Yield_MC->Value 
						* sqrt( pow( Ntg->Sys_Err / Ntg->Value, 2)
								+ pow( XS_SAMC_SUM->Sys_Err / XS_SAMC_SUM->Value, 2) 
								+ pow( Nf_SAMC_Gen_All->Sys_Err / Nf_SAMC_Gen_All->Value, 2) );

					Yield_MC_Avg->Value   = 0.0;
					Yield_MC_Avg->Stat_Err= 0.0;
					Yield_MC_Avg->Sys_Err = 0.0;
					int N_Y_MC = 0;
					for(int i=0;i<Chain_Size;i++){
						Y_MC_Chain[i] = NtgChain[i] *XS_SAMC_SUM->Value*(SAMC_dOmega*(P_MeV*SAMC_dDp)/Nf_SAMC_Gen_All->Value);                
						Y_MC_StatChain[i] = Y_MC_Chain[i]
							* sqrt( 1.0 / NtgChain[i]
									+ pow( XS_SAMC_SUM->Stat_Err / XS_SAMC_SUM->Value, 2) 
									+ pow( Nf_SAMC_Gen_All->Stat_Err / Nf_SAMC_Gen_All->Value, 2) );
						Y_MC_SysChain[i] = Y_MC_Chain[i]
							* sqrt(   pow( Ntg->Sys_Err / Ntg->Value, 2)
									+ pow( XS_SAMC_SUM->Sys_Err / XS_SAMC_SUM->Value, 2) 
									+ pow( Nf_SAMC_Gen_All->Sys_Err / Nf_SAMC_Gen_All->Value, 2) );

						if( Y_MC_Chain[i]>1. && (Y_MC_StatChain[i]/Y_MC_Chain[i])<1.0 ){
							N_Y_MC ++;
							Yield_MC_Avg->Value   += Y_MC_Chain[i];
							Yield_MC_Avg->Stat_Err+= pow( Y_MC_StatChain[i], 2 ); 
							Yield_MC_Avg->Sys_Err += pow( Y_MC_SysChain[i] , 2 ); 
						}
						cout<<Form("   ---Run#%d, Yield_MC = %7.2e", RunChain[i],Y_MC_Chain[i])<<endl;
						outlog<<Form("   ---Run#%d, Yield_MC = %7.2e", RunChain[i],Y_MC_Chain[i])<<endl;
					}
					if(N_Y_MC>0){
						Yield_MC_Avg->Value   /= N_Y_MC;
						Yield_MC_Avg->Stat_Err = sqrt( Yield_MC_Avg->Stat_Err / N_Y_MC);
						Yield_MC_Avg->Sys_Err  = sqrt( Yield_MC_Avg->Sys_Err  / N_Y_MC);
					}
				}
				else
					Yield_MC->SetValue(-1.0,1e-30,1e-30);
			}//if(bCal_XS) in Y_MC
			/*Y_MC}}}*/

			cout   << Form( "       Yield_EX=%7.2e(%7.2e), Yield_MC=%7.2e(%7.2e), Yield_AL=%7.2e(%7.2e)", 
					Yield_EX->Value,Yield_EX->Stat_Err, Yield_MC->Value,Yield_AL->Stat_Err, Yield_AL->Value,Yield_AL->Stat_Err)<<endl;
			outlog << Form( "       Yield_EX=%7.2e(%7.2e), Yield_MC=%7.2e(%7.2e), Yield_AL=%7.2e(%7.2e)",
					Yield_EX->Value,Yield_EX->Stat_Err, Yield_MC->Value,Yield_AL->Stat_Err, Yield_AL->Value,Yield_AL->Stat_Err)<<endl;
			/*}}}*/

			/*Calculate XS and Store{{{*/																	  
			Cross_Section->SetValue(-1.0,1e-30,1e-30);	
			if(bCal_XS){
				if (  Yield_EX->Value>1.0 ){

					//Still Keep the Yield_EX results when xbj > A
					if(Yield_MC->Value>1.0){
						///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						//Note: Live_Time is already corrected when calculating Nf_EX for each run#, Ntg is already in cm^-2			    
						//      Ep is directly caluclate from the center value of Xbj bin, but not the central momentum
						XS_Born_XEMC = gGet_XS_Born_XEMC(Target_Name.Data(),Kin.Data(),Ebeam,Ep,Angle);
						Cross_Section->Value = Yield_EX->Value/Yield_MC->Value * XS_Born_XEMC;
						///////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
						///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						Cross_Section->Stat_Err=Cross_Section->Value 
							* sqrt( pow(Yield_EX->Stat_Err   / Yield_EX->Value, 2)
									+       pow(Yield_MC->Stat_Err / Yield_MC->Value, 2)) ;
						Cross_Section->Sys_Err=Cross_Section->Value 
							* sqrt( pow(Yield_EX->Sys_Err   / Yield_EX->Value, 2)
									+       pow(Yield_MC->Sys_Err / Yield_MC->Value, 2)) ;
						XS_Avg->Value    = 0.0;
						XS_Avg->Stat_Err = 0.0;
						XS_Avg->Sys_Err  = 0.0;	
						int N_XS = 0;
						for(int i=0;i<Chain_Size;i++){
							XS_Chain[i] = Y_EX_Chain[i]/Y_MC_Chain[i] * XS_Born_XEMC; 
							XS_StatChain[i] = XS_Chain[i] 
								* sqrt( pow( Y_EX_StatChain[i] / Y_EX_Chain[i], 2) 
										+       pow( Y_MC_StatChain[i] / Y_MC_Chain[i], 2)) ;
							XS_SysChain[i] = XS_Chain[i] 
								* sqrt( pow( Y_EX_SysChain[i] / Y_EX_Chain[i], 2) 
										+       pow( Y_MC_SysChain[i] / Y_MC_Chain[i], 2)) ;
							//If the Stat_Err is less than 50%, add this run. This will reject those bad runs or short runs.
							if( XS_Chain[i]>1e-33 && (XS_StatChain[i]/XS_Chain[i])<1.0 ){
								N_XS++;
								XS_Avg->Value    += XS_Chain[i];
								XS_Avg->Stat_Err += pow( XS_StatChain[i], 2);
								XS_Avg->Sys_Err  += pow( XS_SysChain[i] , 2);
							}
						}
						if(N_XS>0){
							XS_Avg->Value   /= N_XS;
							XS_Avg->Stat_Err = sqrt(XS_Avg->Stat_Err) / N_XS;
							XS_Avg->Sys_Err  = sqrt(XS_Avg->Sys_Err ) / N_XS;
						}
					}
					else{
						Cross_Section->SetValue(-1.0,1e-30,1e-30);	
						Yield_MC->SetValue(-1.0,1e-30,1e-30);	
					}
				}	// if ( Yield_EX->Value>0)
			}// If(bCal_XS)
			/*}}}*/
			
			/*Output{{{*/
			if(!(bCal_XS&&Cross_Section->Value<-1e-33)){
				outstring=Form("%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e%*e",
						16,Cross_Section->Value,
						16,Cross_Section->Stat_Err,
						16,Cross_Section->Sys_Err,
						16,XS_Born_XEMC,
					   	16,xbj,  //Xbj
						16,Angle, //Deg
						16,Ep,//GeV 
						16,nu,  //Nu
						16,qsq, //Q2
						16,BinCor, //Bin Centering Correction
						16,Livetime->Value,
						16,Ntg->Value,
						16,Ne->Value,
						16,Nf_EX->Value,
						16,Nf_AL->Value,
						16,XS_SAMC_SUM->Value,
						16,Nf_SAMC_Gen_All->Value,
						16,Yield_EX->Value,
						16,Yield_EX->Stat_Err,
						16,Yield_EX->Sys_Err,
						16,Yield_MC->Value,
						16,Yield_MC->Stat_Err,
						16,Yield_MC->Sys_Err,
						16,Yield_AL->Value,
						16,Yield_AL->Stat_Err,
						16,Yield_AL->Sys_Err
							);
				outfile<<outstring<<endl;
			}//if(bCal_XS&&Cross_Section->Value>1e-33){
			/*}}}*/ 

			outlog <<endl
				<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
			outlog << Form("      --> Xbj=%5.3f, Cross_Section=%7.6e / %7.6e / %7.6e",xbj, Cross_Section->Value, Cross_Section->Sys_Err, Cross_Section->Stat_Err)<<endl;
			outlog <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;   
			cout <<endl
				<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
			cout << Form("      --> Xbj=%5.3f, Cross_Section=%7.6e / %7.6e / %7.6e",xbj, Cross_Section->Value, Cross_Section->Sys_Err, Cross_Section->Stat_Err)<<endl;
			cout <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;   
			NoCal++;
			NumOfAccInput++;

			Acc_var_Value+=bin_size;

			T->Fill();	
		}/*}}}Xbj*/
	}//if ( Total_Pion_Rej->Value>0&&Total_E_Eff->Value>0&&Ne->Value>0&&Ntg->Value>0&&Livetime->Value>0&&P0->Value>0 )
	/*}}}*/

	fxs->cd(); T->Write(); fxs->Close(); 

	/*free memory{{{*/
	outfile.close();
	outlog.close();
	delete Total_E_Eff;
	delete Total_Pion_Rej;
	delete Cross_Section;
	delete XS_SAMC_SUM;
	delete Yield_EX;
	delete Yield_AL;
	delete Yield_MC;
	delete Yield_EX_Avg;
	delete Yield_MC_Avg;
	//delete Nf_SAMC;
	delete Nf_EX;
	delete Nf_AL;
	delete Nf_SAMC_Gen_All;
	delete LiveTimeChain;
	delete LiveTimeStatErrChain;
	delete LiveTimeSysErrChain;
	delete NtgChain;
	delete E0Chain;
	delete dEEChain;
	delete AngleChain;
	delete P0Chain;
	delete RunChain;
	delete NeChain;
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
char* gComment=NULL;
int getargs(int argc, char**argv);

/*int main{{{*/
int main(int argc, char**argv){
	int pserr = getargs(argc,argv);
	if(pserr<0){
		outlog<<"****** I can not recognize the parameters!"<<endl;
		cout<<"****** I can not recognize the parameters!"<<endl;
		return -1;
	}

	TString InputFileName = "xs_input/xs_c12_kin3.1.dat"; //outlog <<"Input File Name: "<<endl; cin >> InputFileName; 
	InputFileName = gInput;
	TString LogFileName = "log/xs_c12_kin3.1.log";  // outlog <<"Log File Name: "<<endl; cin >> LogFileName;
	LogFileName = gLog;
	TString CommentName = gComment;

	int err = CrossSection(InputFileName.Data(), LogFileName.Data(),CommentName.Data());

	if(err>0)
		outlog<<"@@@@@@@@@ Code runs successfully!"<<endl;
	cout  <<"@@@@@@@@@ Code runs successfully!"<<endl;

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
				case 'C':
					if(*(++argptr))
						gComment = argptr;
					err=0;
					break;
				case 'h':
					cout <<"================================================"<<endl;
					cout <<"Option: " <<endl;
					cout <<"-I[xs_c12_kin3.1.dat] Input File" <<endl;
					cout <<"-L[c12_kin3.1.log] Log File"<<endl;
					cout <<"-C[Dp0.04] Comment, such as Dp=4%"<<endl;
					cout <<endl<<"================================================"<<endl;
					noStop=false;
					err=-1;
					goto OUT;
					break;
				default:
					cout <<"================================================"<<endl;
					cout <<"Option: " <<endl;
					cout <<"-I[xs_c12_kin3.1.dat] Input File" <<endl;
					cout <<"-L[c12_kin3.1.log] Log File"<<endl;
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
