/******************************************/
/*     X>2 Cross Section Calculation      */
/*         Zhihong Ye, 10/14/2011         */
/*         --Bin by Bin Correction        */
/******************************************/
/*Include{{{*/
#include "/work/halla/e08014/disk1/yez/Ep_Cross_Section/SRC/XGT2_Main.h"
#include "/work/halla/e08014/disk1/yez/Ep_Cross_Section/SRC/XGT2_Beam.h"
/*}}}*/

/*int GetQeI(const TString& aInputFile,const TString& aLogFile){{{*/
int GetQeI(const TString& aInputFile,const TString& aOutputFile,const TString& aLogFile){

    ifstream runlist(aInputFile.Data());	
    ofstream outfile(aOutputFile.Data());	
	Beam_Trip_Cut = 10.0;

	outfile<<Form("%8s %12s %12s %12s %12s %12s %12s",
			"RunNo",
			"I-L",
			"I-L Err",
			"I-R",
			"I-R Err",
			"Q",
			"Q Err"
			)<<endl;

	XGT2_VAR* IL=new XGT2_VAR();
	XGT2_VAR* IR=new XGT2_VAR();
	XGT2_VAR* Ne=new XGT2_VAR();
	Int_t RunNo = 0;    

	int dum=0;
	while(!(runlist.eof())){
		runlist>>RunNo>>dum>>dum>>dum;
	
		/*Get Current{{{*/
		TChain *Tree = gAddTree(RunNo,"T");
		IR=gGet_Current(Tree,"R");
		IL=gGet_Current(Tree,"L");
		/*}}}end of Get Ne*/

		/*Get Ne{{{*/
		TChain *Right_Tree = gAddTree(RunNo,"RIGHT");
		Ne = gGet_BeamCharge(Right_Tree, "R");
		/*}}}end of Get Ne*/

		outfile<<Form("%8d %12e %12e %12e %12e %12e %12e",
			    RunNo,
			    IL->Value,
			    IL->Sys_Err,
			    IR->Value,
			    IR->Sys_Err,
			    Ne->Value,
			    Ne->Sys_Err
				)<<endl;
		cerr<<Form(" === %8d %12e %12e %12e %12e %12e %12e",
				RunNo,
			    IL->Value,
			    IL->Sys_Err,
			    IR->Value,
			    IR->Sys_Err,
			    Ne->Value,
			    Ne->Sys_Err
				)<<endl;


		delete Right_Tree;
		delete Tree;
	}

	runlist.close();
	outfile.close();
	return 0;
}
/*}}}*/

/*int main{{{*/
int main(int argc, char**argv){
	TString InputFileName = "runlist.txt"; //outlog <<"Input File Name: "<<endl; cin >> InputFileName; 
	TString LogFileName = "ntg_qe.log";  // outlog <<"Log File Name: "<<endl; cin >> LogFileName;
	TString OutputFileName = "XGT2_Beam_Charge.dat";  // outlog <<"Log File Name: "<<endl; cin >> LogFileName;

	int err = GetQeI(InputFileName.Data(),OutputFileName.Data(), LogFileName.Data());

	if(err>0)
		outlog<<"@@@@@@@@@ Code runs successfully!"<<endl;
	cerr  <<"@@@@@@@@@ Code runs successfully!"<<endl;

}
/*}}}*/

