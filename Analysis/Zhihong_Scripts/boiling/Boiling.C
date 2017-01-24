/******************************************/
/*           X>2 Boiling Study            */
/*         Zhihong Ye, 07/17/2013         */
/******************************************/
/*Include{{{*/
#include "./SRC/XGT2_Main.h"
#include "./SRC/XGT2_LiveTime.h"
#include "./SRC/XGT2_Beam.h"
#include "./SRC/XGT2_Data.h"
#include "./SRC/XGT2_Cuts.h"
#include "./SRC/XGT2_Target.h"
#include "./SRC/XGT2_Physics.h"
/*}}}*/

/*int Boiling(const TString& aInputFile){{{*/
int CrossSection(const TString& aInputFileName,const TString& Arm, const TString& aComment){
	/////////////////////////////
	// Initialize Input values
	/////////////////////////////
	/*Read InputFile{{{*/
	cerr   << endl << Form("@@@@@@ Initialize Input values from %s ...",aInputFileName.Data())<<endl;
    vector <double> CurrentChain; 
	vector <int> RunNoChain;
	double current = -1.0; int runno = -1;
	int i=0;
	Double_t* RunChain = new double[20];
    TString TargetName = "C12";
	ifstream infile(aInputFileName.Data());
    infile >> TargetName;
	while(!infile.eof()){
		infile >> current >> runno;
		if(!(runno==RunChain[i-1])){
			CurrentChain.push_back(current);	
			RunNoChain.push_back(runno);	
			RunChain[i]=runno;
		}
		i++;
	}
	infile.close();

	const int Chain_Size = RunNoChain.size();
	Beam_Trip_Cut = 5.0;   cerr << Form("   Beam Tip Cut = %6.2f uA",Beam_Trip_Cut)<<endl;
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
	
	ElectronCuts = gGet_XGT2_Cut(Arm,TargetName,AccCut_FileName,&Target_Length_Cut); 
	cerr <<"--------------------------------------------------------------------------------------"<<endl;
	cerr   <<"     --> "<<ElectronCuts.Data()<<endl;
	cerr   <<"--------------------------------------------------------------------------------------"<<endl;
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

	ofstream outfile(Form("%s_Yield_%s.dat", TargetName.Data(), Arm.Data()));
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

	///////////////////////////////////
	// Calculate Yield and Save
	///////////////////////////////////
	double Y_Err = 0.0;
	Target_Length_Cut = 0.12;
    double offset = 0.0;
	double VZ_Max =  Target_Length_Cut + offset; 
	double VZ_Min = -Target_Length_Cut + offset; 
	int bin = 20;
	if(TargetName=="C12"||TargetName=="Ca40"||TargetName=="Ca48")
		bin = 1;
	double step = (VZ_Max-VZ_Min)/bin;

	/*Calculate Yield and Save{{{*/
	VZ = VZ_Min + step/2.0;
	outfile << Chain_Size;
	for(int i=0;i<Chain_Size;i++){
	   outfile <<"  "<< RunChain[i];
	}
	outfile<<endl;

	while ( VZ<=VZ_Max ){
		cerr  <<Form(" +++ Calculating VZ = %f  ...", VZ)<<endl;

		/*Get_Nf_Bin{{{*/ 
		Nf_EX = gGet_Nf_EX(RunNoChain,Arm,ElectronCuts.Data(),VZ,step,LiveTimeChain,LiveTimeStatErrChain,LiveTimeSysErrChain,
				PSChain,Nf_EX_Chain,Nf_EX_StatChain,Nf_EX_SysChain);
		cerr   << "      Total Number of EX Events = " << Nf_EX->Value<<endl;

		/*}}}*/

		cerr  <<Form(" +++ Now VZ = %f  ...", VZ)<<endl;
		/*Yield{{{*/
		//For Experiment Data
	    outfile << VZ <<endl;
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
			}
			else{
				Y_EX_Chain[i]     = -1.;
				Y_EX_StatChain[i] = -1.;
				Y_EX_SysChain[i]  = -1.;
				Y_Err = -1.;
			}
			cerr<< " #### VZ="<<VZ<<" RunNo="<<RunChain[i]<<" Y="<<Y_EX_Chain[i]<<endl;
			outfile <<"  "<<RunChain[i]<<"  "<< Current[i]<<"  "<<Y_EX_Chain[i]
				    <<"  "<<Y_Err<<endl;
		}
		VZ += step;
		NoCal++;
	}/*}}}VZ*/
	/*}}}*/

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

char* gInput;
char* gArm;
char* gComment=NULL;
int getargs(int argc, char**argv);

/*int main{{{*/
int main(int argc, char**argv){
	int pserr = getargs(argc,argv);
	if(pserr<0){
		cerr<<"****** I can not recognize the parameters!"<<endl;
		return -1;
	}

	TString InputFileName = "C12.txt";  
	InputFileName = gInput;
	TString ArmName = "L"; 
	ArmName = gArm;
	TString CommentName = gComment;

	int err = CrossSection(InputFileName.Data(), ArmName.Data(),CommentName.Data());

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
				case 'A':
					if(*(++argptr))
						gArm = argptr;
					err=0;
					break;
				case 'C':
					if(*(++argptr))
						gComment = argptr;
					err=0;
					break;
				case 'h':
					cerr <<"================================================"<<endl;
					cerr <<"Option: " <<endl;
					cerr <<"-I[C12.txt] Input File" <<endl;
					cerr <<"-A[L] HRS Arm"<<endl;
					cerr <<"-C[ ] Comment, such as Dp=4%"<<endl;
					cerr <<endl<<"================================================"<<endl;
					noStop=false;
					err=-1;
					goto OUT;
					break;
				default:
					cerr <<"================================================"<<endl;
					cerr <<"Option: " <<endl;
					cerr <<"-I[C12.txt] Input File" <<endl;
					cerr <<"-A[L] HRS Arm"<<endl;
					cerr <<"-C[ ] Comment, such as Dp=4%"<<endl;
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
