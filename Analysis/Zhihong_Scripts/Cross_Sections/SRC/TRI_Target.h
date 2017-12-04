/*inline TRI_VAR* gCal_Ntg(const vector<Int_t>& aRunNoChain, const TString& aArm, const Int_t aTarget_A, const Int_t aTarget_Z){{{*/
inline TRI_VAR* gCal_Ntg(const vector<Int_t>& aRunNoChain,const TString& aArm,const int aTarget_A, const int aTarget_Z, 
        const Double_t aTarget_Thickness,const Double_t aTarget_Thickness_Stat_Err, double* aNtgChain){
  TRI_VAR* aNtg = new TRI_VAR();
  Double_t aNtg_Sum = 0.0, aNtg_Err=0.0;
  for (unsigned int i=0; i<aRunNoChain.size(); i++ ){
      Int_t aRunNo=aRunNoChain[i];
      TChain *aT_Tree = gAddTree(aRunNo,"T"); 
	  Double_t aCurrent = gGet_Current_Avg(aT_Tree, aArm.Data());
	  Double_t aTarget_Thickness_Corrected = 0.0;
	  Double_t aTarget_Thickness_Corrected_Err = 0.0;
      
      aNtgChain[i] = aTarget_Thickness * Na / aTarget_A ;

	  outlog<<Form("      For Run#=%d, Current = %6.2f, rho= %8.6e, Ntg = %4.3e",aRunNo, aCurrent, aTarget_Thickness_Corrected, aNtgChain[i])<<endl;
	  cerr  <<Form("      For Run#=%d, Current = %6.2f, rho= %8.6e, Ntg = %4.3e",aRunNo, aCurrent, aTarget_Thickness_Corrected, aNtgChain[i])<<endl;
      aNtg_Sum += aNtgChain[i];
      aNtg_Err += pow((aTarget_Thickness_Stat_Err/aTarget_Thickness),2);
      delete aT_Tree;
  }
  //Take the average value of all runs
  if(aRunNoChain.size()>1e-16){
	  aNtg->Value    = aNtg_Sum/aRunNoChain.size(); 
	  aNtg->Stat_Err = 0.0; 
	  aNtg->Sys_Err  = aNtg->Value * sqrt(aNtg_Err/aRunNoChain.size()); 
  }
  else
	  aNtg->SetValue(-1.0,0.0,0.0);

  return aNtg;
}
/*}}}*/
