/*************************************************/
/*           X>2 Dummy Contaimination Study      */
/*         Zhihong Ye, 07/17/2013                */
/*************************************************/
/*Include{{{*/
#include "./SRC/XGT2_Main.h"
#include "./SRC/XGT2_LiveTime.h"
#include "./SRC/XGT2_Beam.h"
#include "./SRC/XGT2_Data.h"
#include "./SRC/XGT2_Cuts.h"
#include "./SRC/XGT2_Target.h"
#include "./SRC/XGT2_Physics.h"
/*}}}*/
 const double CM2M = 1./100.;
 const double M2CM = 100.;
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
		aEp_Min[0]= 1000.0;
		aEp_Max[0]=-1000.0;
	}

}
/*}}}*/

/*int gGet_Yield_Zbin(const TString& aInputFile, const TString& aComment){{{*/
int gGet_Yield_Zbin(const TString& aInputFile, const TString& aComment){
	/////////////////////////////
	// Initialize Input values
	/////////////////////////////
	/*Read InputFile{{{*/
	TString aInputFileName = "./xs_input/"+aInputFile;
	outlog << endl << Form("@@@@@@ Initialize Input values from %s ...",aInputFileName.Data())<<endl;
	cerr   << endl << Form("@@@@@@ Initialize Input values from %s ...",aInputFileName.Data())<<endl;
	Int_t RunNo=-1;
	Int_t aline=0;
	const TString Arm=gGet_InputFile_Var(aInputFileName,aline);aline++;                outlog << "   ARM = " << Arm <<endl;
	const TString Kin=gGet_InputFile_Var(aInputFileName,aline);aline++;                outlog << "   Kin = " << Kin<<endl;
	const TString Target_Name=gGet_InputFile_Var(aInputFileName,aline);aline++;        outlog << "   Target Name = " << Target_Name<<endl;
	const Bool_t IsExtTgt=gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++;     outlog << "   Extended Target? " << IsExtTgt<<endl;
	const Int_t Main_Trigger=gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++;  outlog << "   Main Trigger " << Main_Trigger<<endl;
	const TString S_Tree_Name=gGet_InputFile_Var(aInputFileName,aline);aline++;        outlog << "   Scaler Tree " << S_Tree_Name<<endl;
	const TString T_Tree_Name="T";
	const TString E_Tree_Name="E";
	const TString DECAT_Name="DBB";
	const vector<Int_t> RunNoChain=gGet_RunNoChain(RunNo,gGet_InputFile_Var(aInputFileName,aline));aline++;	
	const Int_t Chain_Size = RunNoChain.size();
	int* RunChain = new int[Chain_Size];
	for(int i=0;i<Chain_Size;i++)
		RunChain[i] = RunNoChain[i];

	gCheckDirectory(Form("%s/results/%s",XGT2_DIR.Data(),gGetScriptname(__FILE__).Data()));
	Beam_Trip_Cut = gGet_InputFile_Var(aInputFileName,aline).Atof();aline++;   outlog << Form("   Beam Tip Cut = %6.2f uA",Beam_Trip_Cut)<<endl;
	/*}}}end of reading input*/

	///////////////////////////////////
	// Number of Incident particles
	///////////////////////////////////
	/*Get Ne{{{*/
	cerr   << endl << "@@@@@@ Calculating Incident Particle ..."<<endl;	
	Double_t runtime;
	Double_t* NeChain = new double[Chain_Size];
	const XGT2_VAR* Ne=gCal_Ne(RunNoChain,Arm,runtime, NeChain);
	cerr   <<Form("     Ne = %4.3e", Ne->Value)<<endl;	
	/*}}}end of Get Ne*/

	//////////////////////////
	// Check PreScalerR
	//////////////////////////
	/*Get PS{{{*/
	cerr   << endl << "@@@@@@ Checking and Storing PreScaler Factors for every run ..." <<endl;	
	Int_t* PSChain = new int[Chain_Size*8];
	gCheck_PS(RunNoChain,PSChain);
	/*}}}end of PS*/

	//////////////////////////
	// LiveTime
	//////////////////////////
	/*LiveTime{{{*/
	cerr    << endl<< "@@@@@@ Checking and Storing Live Time for every run ..." <<endl;
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
	/*}}}end of LiveTime*/

	//////////////////////
	// Check Current
	//////////////////////
	//Check and Print out values of currents from each Arm for each run.
	double *Current = new double[Chain_Size];
	gCheck_Current(RunNoChain, Arm.Data(), Current);

	//////////////////
	// Set Cuts
	//////////////////
	/*Set Cuts{{{*/
	TString ElectronCuts;
	Double_t Target_Length_Cut = 0.0; //Symmetry Cuts, like |VZ|<Target_Length_Cut
	TString AccCut_FileName = Form("./Acceptance/Acc_%s_Table.dat", Arm.Data());
	cerr   << endl << Form("@@@@@@ Obtaining Cuts from %s ...", AccCut_FileName.Data()) <<endl;
	
	ElectronCuts = gGet_XGT2_Cut(Arm,Target_Name,AccCut_FileName,&Target_Length_Cut); 
	cerr <<"--------------------------------------------------------------------------------------"<<endl;
	cerr   <<"     --> "<<ElectronCuts.Data()<<endl;
	cerr   <<"--------------------------------------------------------------------------------------"<<endl;
	
	double Ep_Min=1000.0, Ep_Max=-1000.0;
    const double Ep_Size = 0.02;
	gGet_Ep_Cut(Kin.Data(),Arm.Data(), &Ep_Min,&Ep_Max);
	
	/*}}}end of Cuts*/

	////////////////////////////////////
	// Calculate Yield 
	////////////////////////////////////
	cerr <<endl<< "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl;
	cerr       << "oO0Oo              --- Calculating Yield ---              oO0Oo" <<endl;
	cerr       << "oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo "<<endl<<endl;

	/////////////////////
	// Define variables
	/////////////////////
	/*Define variables{{{*/
	XGT2_VAR* Total_E_Eff=new XGT2_VAR();
	XGT2_VAR* Total_Pion_Rej=new XGT2_VAR();
	//For X>2 Experiment, PID efficiencies are all above 99%, so we ignore the values and put them into sysmetic errors
	Total_E_Eff->SetValue(1,0,0.01); 
	Total_Pion_Rej->SetValue(1,0,0.01);
	XGT2_VAR* Nf_EX=new XGT2_VAR();  //Total Number of real events in each bin after all cuts
    double VZ = -1000.0;
    double Ep = -1000.0;

	ofstream outfile(Form("./output/%s_Yield_%s_%s_ep.dat", Target_Name.Data(), Arm.Data(), Kin.Data()));
	/*}}}*/

	Int_t NoCal=0;
	map<Double_t,Double_t> map_err;

    ////////////////
	// Set Variables 
	////////////////
	/*Set Variables{{{*/
	Double_t Nf_EX_Chain[Chain_Size];	
    Double_t Nf_EX_StatChain[Chain_Size];	
    Double_t Nf_EX_SysChain[Chain_Size];	
    Double_t Y_EX_Chain[Chain_Size];	
    Double_t Y_EX_StatChain[Chain_Size];	
    Double_t Y_EX_SysChain[Chain_Size];	
	/*}}}*/
	
	/*Kinematics Settings{{{*/
	////////////////////////////////////////////
	//Check Beam Energy and Store all P0
	/////////////////////////////////////////////  
	outlog  << endl<< "@@@@@@ Checking Beam Energy for each run ..." <<endl; 
	cerr    << endl<< "@@@@@@ Checking Beam Energy for each run ..." <<endl; 
	Double_t* E0Chain = new double[Chain_Size]; //GeV
	Double_t* dEEChain = new double[Chain_Size];//GeV
	gCheck_E0(RunNoChain,Arm,E_Tree_Name,E0Chain, dEEChain); //Check P0 for each run from the Dipole readings. 

	////////////////////////////////////////////
	//Check Central Momentum and Store all P0
	/////////////////////////////////////////////  
	outlog  << endl<< "@@@@@@ Checking Central Momentum for each run ..." <<endl; 
	cerr    << endl<< "@@@@@@ Checking Central Momentum for each run ..." <<endl; 
	Double_t* P0Chain = new double[Chain_Size];//GeV/c
	gCheck_P0(RunNoChain,Arm,E_Tree_Name,P0Chain); //Check P0 for each run from the Dipole readings. 

	////////////////////////////////////////////
	//Check Scattering Angle and Store all values
	/////////////////////////////////////////////  
	//	Angle = gGet_Angle_Corrected(RunNoChain, Arm, Angle);
	outlog  << endl<< "@@@@@@ Checking Scattering Angle for each run ..." <<endl; 
	cerr    << endl<< "@@@@@@ Checking Scattering Angle for each run ..." <<endl; 
	Double_t* AngleChain = new double[Chain_Size];//Degree
	gCheck_Angle(RunNoChain,Arm,E_Tree_Name,AngleChain); //Check P0 for each run from the Dipole readings. 
	/*}}}*/

	///////////////////////////////////
	// Calculate Yield and Save
	///////////////////////////////////
	double Y_Err = 0.0;
	Target_Length_Cut = 15.*CM2M;
    double offset = 0.0;
	if(Arm=="L"){
		if(Target_Name=="H2"||Target_Name=="He3"||Target_Name=="He4")
		offset = -0.5*CM2M;
	if(Target_Name=="Dummy20cm")
		offset = -0.5*CM2M;
	}

	if(Arm=="R"){
		if(Target_Name=="H2"||Target_Name=="He3"||Target_Name=="He4")
		offset = -0.5*CM2M;
	if(Target_Name=="Dummy20cm")
		offset = -0.5*CM2M;
	}

	double VZ_Max =  Target_Length_Cut-offset; 
	double VZ_Min = -Target_Length_Cut-offset; 
	double Yield_Total, Yield_Err_Total;
	int bin = 200;
	if(Target_Name=="C12"||Target_Name=="Ca40"||Target_Name=="Ca48")
		bin = 1;
	double step = (VZ_Max-VZ_Min)/bin;

	/*Calculate Yield and Save{{{*/
    Ep = Ep_Min + Ep_Size/2.0;
	/*Ep Bin{{{*/
	while (Ep<=Ep_Max){
		VZ = VZ_Min + step/2.0 ;// - offset;
		/*VZ Bin{{{*/
		while ( VZ<=VZ_Max ){
			cerr  <<Form(" +++ Calculating Ep = %f, VZ = %f  ...", Ep, VZ)<<endl;

			/*Get_Nf_Bin{{{*/ 
			Nf_EX = gGet_Nf_EX(RunNoChain,Arm,Ep, Ep_Size,VZ,step,ElectronCuts.Data(),LiveTimeChain,LiveTimeStatErrChain,LiveTimeSysErrChain,
					PSChain,E0Chain, P0Chain, AngleChain,IsExtTgt, Nf_EX_Chain,Nf_EX_StatChain,Nf_EX_SysChain);
			cerr   << "      Total Number of EX Events = " << Nf_EX->Value<<endl;
			/*}}}*/

			cerr  <<Form(" +++ Now VZ = %f  ...", VZ)<<endl;
			/*Yield{{{*/
			Yield_Total = 0.0;
			Yield_Err_Total = 0.0;
			for(int i=0;i<Chain_Size;i++){
				Y_EX_Chain[i]     = 1e-32;
				Y_EX_StatChain[i] = 1e-32;
				Y_EX_SysChain[i]  = 1e-32;
				Y_Err = 1e-32;
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
					Y_Err = sqrt(pow(Y_EX_StatChain[i],2)+pow(Y_EX_SysChain[i],2));

					Yield_Total +=Y_EX_Chain[i];
					Yield_Err_Total = pow(Y_Err,2);
				}
				else{
					Y_EX_Chain[i]     = -1.;
					Y_EX_StatChain[i] = -1.;
					Y_EX_SysChain[i]  = -1.;
					Y_Err = -1.;
				}
				cerr<< " #### VZ="<<VZ<<" RunNo="<<RunChain[i]<<" Y="<<Y_EX_Chain[i]<<endl;
			}
			Yield_Err_Total = sqrt(Yield_Err_Total);
			outfile <<Form("%8.6f %8.6f %12e %12e", Ep, VZ, Yield_Total, Yield_Err_Total )<<endl;
			/*}}}*/

			if((VZ>-0.12&&VZ<-0.08) || (VZ>0.07&&VZ<0.11))
				VZ += step;
			else
				VZ += step*20.0;//Speed up!
			NoCal++;
		}//while(VZ<=VZ_Max)
		/*}}}VZ Bin*/
		Ep += Ep_Size;
	}//while(Ep<=Ep_Max)
	/*}}}Ep Bin*/

	/*free memory{{{*/
	delete Total_E_Eff;
	delete Total_Pion_Rej;
	delete Nf_EX;
	delete LiveTimeChain;
	delete LiveTimeStatErrChain;
	delete LiveTimeSysErrChain;
	delete RunChain;
	delete NeChain;
	delete PSChain;
	/*}}}*/
	return 0;
}
/*}}}*/

/*int main{{{*/
char* gInput;
char* gArm;
int getargs(int argc, char**argv);

int main(int argc, char**argv){
	int pserr = getargs(argc,argv);
	if(pserr<0){
		cerr<<"****** I can not recognize the parameters!"<<endl;
		return -1;
	}

	TString InputFileName = "";  
	InputFileName = gInput;
	TString CommentName = "";

	int err = gGet_Yield_Zbin(InputFileName.Data(), CommentName.Data());

	if(err>0)
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
				case 'h':
					cerr <<"================================================"<<endl;
					cerr <<"Option: " <<endl;
					cerr <<"-I[C12.txt] Input File" <<endl;
					cerr <<endl<<"================================================"<<endl;
					noStop=false;
					err=-1;
					goto OUT;
					break;
				default:
					cerr <<"================================================"<<endl;
					cerr <<"Option: " <<endl;
					cerr <<"-I[C12.txt] Input File" <<endl;
					cerr <<endl<<"================================================"<<endl;
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
