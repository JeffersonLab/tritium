//#include "XGT2_Main.h"
#include "XGT2_Beam.h"
#include "XGT2_LiveTime.h"

inline XGT2_VAR* gGet_Yield(TTree *aTree, TTree *aS_Tree, const Int_t aRunNo, const TString& aArm, Int_t *aPS, const TString& aCuts, const Double_t aEff)
{ 
  XGT2_VAR* aYield = new XGT2_VAR();
  Double_t u1,runtime,charge;
 
  XGT2_VAR* aNi=gCal_Ni_Single(aRunNo,aArm,runtime);
  //  aNi->Value/=1e18;
  //  aNi->Stat_Err = sqrt(aNi->Value);
 
  int aMain_Trigger;
  TString aS_Tree_Name;
  if(aArm =="L"){
    aMain_Trigger = 3;
    aS_Tree_Name = "LEFT"; 
  }
  else if(aArm =="R"){
    aMain_Trigger = 1;
    aS_Tree_Name = "RIGHT"; 
  }

  XGT2_VAR* aLiveTime = gCal_LiveTime(aTree,aRunNo,Form("(DBB.evtypebits&(1<<%d))>0&&right_current>5.0",aMain_Trigger),aS_Tree,Form("t%dc",aMain_Trigger),aPS[aMain_Trigger-1]);

  Double_t aNf = 1.0 * aPS[aMain_Trigger-1]*aTree->GetEntries(aCuts.Data())/aLiveTime->Value;
  Double_t yield = aNf/aNi->Value/aEff; 
  outlog<<Form("   Nf = %f, Ni = %e, LiveTime=%f",aNf,aNi->Value,aLiveTime->Value)<<endl;
  cerr<<Form("   Nf = %f, Ni = %e, LiveTime=%f",aNf,aNi->Value,aLiveTime->Value)<<endl;

  delete aNi;
  delete aLiveTime;
  /* delete aTree; */
  /* delete aS_Tree; */

  if(aNf>1e-16){
	  Double_t aNf_Err = 1.0/sqrt(aNf);
	  Double_t aNi_Err = 0.0; //1.0/sqrt(aNi->Value); //Error value when Ni if hugh.
  }
  else{
	  Double_t aNf_Err = 0.0;
	  Double_t aNi_Err = 0.0; //1.0/sqrt(aNi->Value); //Error value when Ni if hugh.
  }

  aYield->Value = yield;
  aYield->Stat_Err = yield*sqrt(aNf_Err*aNf_Err + aNi_Err*aNi_Err);

  return aYield;
}

inline double gGet_Density(const double aDensity, const double aCurrent, const TString& aTarget_Name){
  double Boiling_Factor = 0.0;
  
  //Boiling Factor are calculated at 100uA
  //Values can be updated here.
  if(aTarget_Name == "D2")
    Boiling_Factor = -0.188;
  else if(aTarget_Name == "He3")
    Boiling_Factor = -0.26;
  else if(aTarget_Name == "He4")
    Boiling_Factor = -0.021;
  else
    Boiling_Factor = 0.0;
  
  double aDensity_Real = aDensity * (1.0 - Boiling_Factor * aCurrent /100.00); 
  return aDensity_Real;
}
