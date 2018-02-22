 //inline Double_t gGetTarget_Boiling(const TString &aArm, const Int_t aTarget_A, const Int_t aTarget_Z, const Double_t aTarget_Density, const Double_t aCurrent)/*{{{*/
 inline Double_t gGetTarget_Boiling(const TString &aArm, const Int_t aTarget_A, const Int_t aTarget_Z, const Double_t aTarget_Density, const Double_t aCurrent){

     Double_t aSlope = 0.0, aOffset = 0.0;
     if(aTarget_A == 1 && aTarget_Z == 1){ //Hydrogen
         aSlope = -0.0051; // +/-0.001
         aOffset = 0.97; // +/-0.015
     }else if(aTarget_A == 2 && aTarget_Z == 1){//Deuteron
         aSlope = -0.004; // +/-0.0005
         aOffset = 1.008; // +/-0.006
     }else if(aTarget_A == 3 && aTarget_Z == 1){//Tritium
         aSlope = -0.0051; // +/-0.0004
         aOffset = 1.008; // +/-0.006
     }else if(aTarget_A == 3 && aTarget_Z == 2){//Helium3
         aSlope = -0.0032; // +/-0.0008
         aOffset = 0.991; // +/-0.01
     }

     Double_t aDensity = aTarget_Density * ( aSlope * aCurrent + aOffset);//aCurrent in uA
     //cout<<Form("--- Target(A=%d, Z=%d): Density(I=%3.1fuA) = %6.3e (%6.3e)",
     //        aTarget_A, aTarget_Z, aCurrent, aDensity, aTarget_Density)<<endl;
     return aDensity;
 }/*}}}*/

/*inline TRI_VAR* gCal_Ntg(const vector<Int_t>& aRunNoChain, const TString& aArm, const Int_t aTarget_A, const Int_t aTarget_Z){{{*/
inline TRI_VAR* gCal_Ntg(const vector<Int_t>& aRunNoChain,const TString& aArm,const int aTarget_A, const int aTarget_Z, 
        const Double_t aTarget_Thickness,const Double_t aTarget_Thickness_Err,const double aTarget_Length, double* aNtgChain){
  TRI_VAR* aNtg = new TRI_VAR();
  Double_t aNtg_Sum = 0.0, aNtg_Err=0.0;
  for (unsigned int i=0; i<aRunNoChain.size(); i++ ){
      Int_t aRunNo=aRunNoChain[i];
	  Double_t aCurrent = gGet_Current_Avg(aRunNo, aArm.Data());
      if(aCurrent>22.5)
          aCurrent=22.5;
	  //Double_t aCurrent = 20.0;//use for boiling effect, has not considered it yet

	  Double_t aTarget_Thickness_Corrected = gGetTarget_Boiling(aArm, aTarget_A, aTarget_Z, aTarget_Thickness, aCurrent) * aTarget_Length/25.00;
	  Double_t aTarget_Thickness_Err_Corrected = aTarget_Thickness_Corrected * aTarget_Thickness_Err/aTarget_Thickness  * aTarget_Length/25.00;
      
      aNtgChain[i] = aTarget_Thickness_Corrected * Na / aTarget_A ;

	  cerr  <<Form("      For Run#=%d, Current = %6.2f, Rho(L=%3.1f cm) = %8.6e, Ntg = %12.4e",aRunNo, aCurrent, aTarget_Length, aTarget_Thickness_Corrected, aNtgChain[i])<<endl;
      aNtg_Sum += aNtgChain[i];
      aNtg_Err += pow((aTarget_Thickness_Err_Corrected/aTarget_Thickness_Corrected),2);
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
