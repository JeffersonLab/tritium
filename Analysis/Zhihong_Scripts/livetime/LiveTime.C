/////////////////////////////////////////////////////////
//This code calculate livetime w and w/o beam trip cuts.
//Before use this,"AddCurrent" and "Scaler" are needed
//to run first to add real time currents and triggers.
//--- Zhihong Ye, 08/29/2012
//////////////////////////////////////////////////////////

//Note: I used rootfiles without any cuts, in /work/Rootfiles/e08014_%d.root
//      In the ./Cross_Section/SRC/XGT2_*.h, I used shorter version of rootfiles which can not be used to calculate livetime.
//      Don't make this kind of different versions mistake!  - Zhihong Ye, 09/08/2012
#include "XGT2_Main.h"
#include "XGT2_LiveTime.h"

int main(){
	//Define input run list
	TString RunList = "runlist.txt";
	cerr<<"Tell me the run list (runlist.txt): "; cin >> RunList;
	ifstream infile; infile.open(RunList);
	TString RootFile_Path = "/work/halla/e08014/disk1/Rootfiles/";
	//Define output LiveTime tables
	ofstream left("LeftR_LiveTimeAll.table");
	ofstream right("RightL_LiveTimeAll.table");
	left <<  Form("%6s %6s %12s %12s %12s %12s %12s %12s %12s","Run#","PS", "LT_Avg","LT_Avg_Err","LT","LT_Err","Rate","I_Max","I_Mean")<<endl;
	right << Form("%6s %6s %12s %12s %12s %12s %12s %12s %12s","Run#","PS", "LT_Avg","LT_Avg_Err","LT","LT_Err","Rate","I_Max","I_Mean")<<endl;

	Int_t RunNo,iL, iR, iCom;
	Int_t PS[8], Main_Trigger=0; 
	TString S_Tree_Name, I_Name;
	TString DAQ_Cut = "1";
	Double_t Beam_Trip_CutL = -1.0,Beam_Trip_CutR = -1.0;
	Double_t I_Mean,rms,I_Max,I_Cut;
    Double_t Rate;

	Int_t counter=0; Bool_t bPrint=kFALSE;
	XGT2_VAR* LT = new XGT2_VAR();
	XGT2_VAR* LT_Avg = new XGT2_VAR();

	while(!(infile.eof())){
		infile >> RunNo >> iL >> iR >> iCom;

   //	iL =1; iR=1;
		TString File_Form=Form("e08014_%d.root",RunNo);
		if ( gSystem->FindFile(RootFile_Path.Data(),File_Form) ){
			counter++; bPrint=kFALSE;

		if(!(counter%10)){
			cerr<<Form("--- Count = %d,Working on Run#%d",counter,RunNo)<<endl;
		    bPrint=kTRUE;
		}
			gGet_PS(RunNo, PS);
			// for(int j=0;j<8;j++){
			//   cerr << Form("    PS%d = %d ",j+1,PS[j]);
			// }
			// cerr<<endl;
			TChain *T_Tree = gAddTree(RunNo,"T");
			//if(iR==1){
			if(1){
				Main_Trigger = 1;
				//S_Tree_Name = "RIGHT";
				S_Tree_Name = "LEFT";
				//I_Name = "right_current";
				I_Name = "left_current";
				I_Cut = 0.5; //Cut out current that lower than half of the max current

				TChain *S_Tree = gAddTree(RunNo,S_Tree_Name.Data());
				I_Max = S_Tree->GetMaximum(I_Name.Data());
				if(I_Max>150.0)
					I_Cut = 0.0;//If current values go crazy, don't use beam trip cut
				//Beam_Trip_CutR = I_Cut * I_Max;//I assume the max current is close to the setting value
				Beam_Trip_CutR = 10.0;
				gGet_Leaf_Mean_and_RMS(S_Tree,I_Name.Data(),Form("%s>=0.0&&%s<=200.0",I_Name.Data(),I_Name.Data()),I_Mean,rms); 
				
				gGet_Leaf_Mean_and_RMS(S_Tree,Form("t%dr",Main_Trigger),"",Rate,rms); 

				//Average Live Time for the entire run
				DAQ_Cut = Form("((DBB.evtypebits>>%d)&1)",Main_Trigger);//No Beam Trip Cut
				LT_Avg=gCal_LiveTime_Avg(T_Tree,RunNo,DAQ_Cut.Data(),S_Tree,Form("t%dc",Main_Trigger),PS[Main_Trigger-1]);
                if(LT_Avg->Value>1.0)
					LT_Avg->Value = 1.0;

				//Live Time after remove events with beam trip
				if(I_Max<130.0){
					DAQ_Cut = Form("((DBB.evtypebits>>%d)&1) && %s>=%f",Main_Trigger, I_Name.Data(), Beam_Trip_CutR);
					LT=gCal_LiveTime(T_Tree,RunNo,DAQ_Cut.Data(),S_Tree,Form("rt_t%dc",Main_Trigger),I_Name.Data(),PS[Main_Trigger-1]);
				}
				else{
			       LT->Value = LT_Avg->Value;
			       LT->Stat_Err = LT_Avg->Stat_Err;
				}
	
			    if(LT->Value<0.4||LT->Value>1.0){
				    LT->Value = LT_Avg->Value;
					LT->Stat_Err = LT_Avg->Stat_Err;
				}       
				delete S_Tree;

				if(bPrint)
					cerr << Form("   RIGHT: Run=%4d PS=%2d LT_Avg=%5.4f LT=%5.4f, Rate = %f, I_max = %f, I_mean = %f", RunNo, PS[Main_Trigger-1], LT_Avg->Value, LT->Value, Rate,I_Max,I_Mean)<<endl;
				right << Form("%6d %6d %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f", RunNo, PS[Main_Trigger-1],LT_Avg->Value, LT_Avg->Stat_Err, LT->Value, LT->Stat_Err,Rate,I_Max,I_Mean)<<endl;
			}

			//if(iL==1){
			if(1){
				Main_Trigger = 3;
				//S_Tree_Name = "LEFT";
				S_Tree_Name = "RIGHT";
				//I_Name = "left_current";
				I_Name = "right_current";
				I_Cut = 0.5;

				TChain *S_Tree = gAddTree(RunNo,S_Tree_Name.Data());
			    I_Max =  S_Tree->GetMaximum(I_Name.Data());
				if(I_Max>150.0)
					I_Cut = 0.0;
			    //Beam_Trip_CutL = I_Cut * I_Max;//I assume the max current is close to the setting value
			    Beam_Trip_CutL = 10.0;
				gGet_Leaf_Mean_and_RMS(S_Tree,I_Name.Data(),Form("%s>=0.0&&%s<=500.0",I_Name.Data(),I_Name.Data()),I_Mean,rms); 

				gGet_Leaf_Mean_and_RMS(S_Tree,Form("t%dr",Main_Trigger),"",Rate,rms); 
				
				//Average Live Time for the entire run
				DAQ_Cut = Form("((DBB.evtypebits>>%d)&1)",Main_Trigger);//No Beam Trip Cut
				LT_Avg=gCal_LiveTime_Avg(T_Tree,RunNo,DAQ_Cut.Data(),S_Tree,Form("t%dc",Main_Trigger),PS[Main_Trigger-1]);

				if(LT_Avg->Value>1.0)
					LT_Avg->Value = 1.0;

				//Live Time after remove events with beam trip, have to use RIGHT scaler due to bad LEFT BCM hence bad currents
				if(I_Max<130.0){
				DAQ_Cut = Form("((DBB.evtypebits>>%d)&1) && %s>=%f",Main_Trigger, I_Name.Data(), Beam_Trip_CutL);
				LT=gCal_LiveTime(T_Tree,RunNo,DAQ_Cut.Data(),S_Tree,Form("rt_t%dc",Main_Trigger),I_Name.Data(),PS[Main_Trigger-1]);
				}
				else{
			       LT->Value = LT_Avg->Value;
			       LT->Stat_Err = LT_Avg->Stat_Err;
				}
	
			    if(LT->Value<0.4||LT->Value>1.0){
				    LT->Value = LT_Avg->Value;
					LT->Stat_Err = LT_Avg->Stat_Err;
				}       
				delete S_Tree;

				if(bPrint)
				    cerr << Form("    LEFT: Run=%4d PS=%2d LT_Avg=%5.4f LT=%5.4f, Rate = %f, I_max = %f, I_mean = %f", RunNo, PS[Main_Trigger-1], LT_Avg->Value, LT->Value, Rate,I_Max,I_Mean)<<endl;
				left << Form("%6d %6d %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f", RunNo, PS[Main_Trigger-1],LT_Avg->Value, LT_Avg->Stat_Err, LT->Value, LT->Stat_Err,Rate,I_Max,I_Mean)<<endl;

			}
			delete T_Tree;
		}
	}
	delete LT_Avg;
	delete LT;
	left.close();
	right.close();
}
