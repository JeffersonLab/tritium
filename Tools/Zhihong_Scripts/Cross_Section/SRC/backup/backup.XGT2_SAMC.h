/*Int_t gGet_Nf_SAMC(TChain* aTree){{{*/
Int_t gGet_Nf_SAMC(TChain* aTree,const Double_t aBin[4],const Double_t aBinCut[4], const Double_t aAccCut[9], const Double_t aCS_Center,const TString &aCuts)
{
 /* Bin Cuts:
      |Dp_tg  -aBin[0]| < aCut[0]; 
      |Nu_tg  -aBin[1]| < aCut[1];
      |Qsq_tg -aBin[2]| < aCut[2];   MeV^2 in SAMC and GeV^2 in Real Data
      |Xbj_tg -aBin[3]| < aCut[3];
    Acc Cuts:
      |x_fp|       < aAccCut[0];
      |y_fp|       < aAccCut[1];
      |th_fp|      < aAccCut[2];
      |ph_fp|      < aAccCut[3];
      |dp_rec|     < aAccCut[4];
      |y_tg_rec|   < aAccCut[5];
      |th_tg_rec|  < aAccCut[6];
      |ph_tg_rec|  < aAccCut[7];
      |reactz_rec| < aAccCut[8];
    Weight:
      w[i] = aCS_Avg / aCS[i];
 */

 Int_t aNf = 0;   
 Int_t aNevt = aTree->GetEntries();

 Double_t x_fp,y_fp,th_fp,ph_fp,dp_rec,y_tg_rec,th_tg_rec,ph_tg_rec,reactz_rec,Angle_rec;
 Double_t E_s,E_p, Theta, Qsq,Xbj,cs_Final;
 Double_t aSum, aWeight;
 Bool_t aIsPassed; 

 aTree->SetBranchAddress("x_fp",         &x_fp);
 aTree->SetBranchAddress("y_fp",         &y_fp);
 aTree->SetBranchAddress("th_fp",        &th_fp);
 aTree->SetBranchAddress("ph_fp",        &ph_fp);
 aTree->SetBranchAddress("dp_rec",       &dp_rec);
 aTree->SetBranchAddress("y_tg_rec",     &y_tg_rec);
 aTree->SetBranchAddress("th_tg_rec",    &th_tg_rec);
 aTree->SetBranchAddress("ph_tg_rec",    &ph_tg_rec);
 aTree->SetBranchAddress("reactz_rec",   &reactz_rec);
 aTree->SetBranchAddress("Angle_rec",    &Angle_rec);
 aTree->SetBranchAddress("E_s",          &E_s);
 aTree->SetBranchAddress("E_p",          &E_p);
 aTree->SetBranchAddress("Theta",          &Theta);
 aTree->SetBranchAddress("Qsq",          &Qsq);
 aTree->SetBranchAddress("Xbj",          &Xbj);
 aTree->SetBranchAddress("cs_Final",     &cs_Final);
 aTree->SetBranchAddress("IsPassed",     &aIsPassed);

 //Define Cuts on each bin
 TString aBinCuts=Form("abs(dp_rec-%f)<=%f&&abs(E_s-E_p-%f)<=%f&&abs(Qsq/1e6-%f)<=%f&&abs(Xbj-%f)<=%f",
                    aBin[0], aBinCut[0],
                    aBin[1], aBinCut[1],
                    aBin[2], aBinCut[2],
                    aBin[3], aBinCut[3]);

 //Get Average Cross Section at the Center Angle
 Double_t aCS_Avg = 0.0, aRMS=0.0;
 TString aCenterCuts = aBinCuts+"&&"+aCuts+"&&"+Form("abs(th_tg_rec - %f)<%f",0.0, 1e-3); //Get CS at Center Scattering Angle for each Bin
 cerr<< "++++  SAMC Cuts for this bin at the Central Angle: "<<endl
      <<"      ------------------------------------------------"<<endl
      << aCenterCuts.Data()<<endl
      <<"      ------------------------------------------------"<<endl<<endl;
 gGet_Leaf_Mean_and_RMS(aTree,"cs_Final",aCenterCuts.Data(),aCS_Avg,aRMS);
 cerr <<Form("       Average Center Angle XS is %f", aCS_Avg)<<endl;

 Int_t aCount = 0;
 aSum =0.0;
 for(int i=0;i<aNevt;i++){
  aTree->GetEntry(i);  
  if(abs(x_fp)       < aAccCut[0] &&
     abs(y_fp)       < aAccCut[1] &&
     abs(th_fp)      < aAccCut[2] &&
     abs(ph_fp)      < aAccCut[3] &&
     abs(dp_rec)     < aAccCut[4] &&
     abs(y_tg_rec)   < aAccCut[5] &&
     abs(th_tg_rec)  < aAccCut[6] &&
     abs(ph_tg_rec)  < aAccCut[7] &&
     abs(reactz_rec) < aAccCut[8] &&
     abs(dp_rec - aBin[0])< aBinCut[0] &&
     abs(E_s-E_p- aBin[1])< aBinCut[1] &&
     abs(Qsq/1e6- aBin[2])< aBinCut[2] &&
     abs(Xbj    - aBin[3])< aBinCut[3] &&
     aIsPassed ){
    //Calculate the value of weight for each event in the cuts:
    aWeight = aCS_Avg / cs_Final;
    aCount++;
    aSum += aWeight;
  }
 }
 cerr <<"      SAMC Events in this bin before correction = " << aCount <<endl;
 cerr <<"      SAMC Events in this bin after Correction  = "<<aSum<<endl<<endl;
 aNf = (Int_t) aSum;
 return aNf ; 
}
/*}}}*/

inline void gGet_XGT2_Acc_Cut(const TString& aArm,const TString& aTarget, const Double_t& aMom, Double_t aAccCuts[9]){
   
  Double_t aX_Fp=0,aY_Fp=0,aTheta_Fp=0,aPhi_Fp=0;      // Focal Plane Cuts
  Double_t aDp_Tg=0,aY_Tg=0,aTheta_Tg=0,aPhi_Tg=0,aRctPt_Z=0;     // Target Plane Cuts
  
  Double_t gMom, gX_Fp, gY_Fp, gTh_Fp, gPh_Fp, gY_Tg, gTh_Tg, gPh_Tg, gDp_Tg,gRctPt_Z;
  TString gKin,gTgt;
  
  if(aArm == "L"){
      //Search the Cut Table
    ifstream InCut("Acc_Left_Table.dat");
    TString comment;
    comment.ReadLine(InCut);
    comment.ReadLine(InCut);
    while(!(InCut.eof())){
      InCut >> gMom >> gKin >> gTgt
	    >> gX_Fp >> gY_Fp >> gTh_Fp >> gPh_Fp
	    >> gDp_Tg >> gY_Tg >> gTh_Tg >> gPh_Tg >> gRctPt_Z;

      if(aTarget.Data() == gTgt && (gMom - aMom) < 0.02 && (gMom - aMom) >-0.02){
	aX_Fp = gX_Fp;   aY_Fp = gY_Fp; aTheta_Fp = gTh_Fp; aPhi_Fp = gPh_Fp; 
	aDp_Tg = gDp_Tg; aY_Tg = gY_Tg; aTheta_Tg = gTh_Tg; aPhi_Tg = gPh_Tg; aRctPt_Z = gRctPt_Z;
	break;
      }
    }
  }
       
  if(aArm == "R"){
    //Search the Cut Table
    ifstream InCut("Acc_Right_Table.dat");
    TString comment;
    comment.ReadLine(InCut);
    comment.ReadLine(InCut);
    while(!(InCut.eof())){
      InCut >> gMom >> gKin >> gTgt
	    >> gX_Fp >> gY_Fp >> gTh_Fp >> gPh_Fp
	    >> gDp_Tg >> gY_Tg >> gTh_Tg >> gPh_Tg >> gRctPt_Z;
     
      if(aTarget.Data() == gTgt.Data() && (gMom - aMom) < 0.002 && (gMom - aMom) >- 0.002){
	aX_Fp = gX_Fp;   aY_Fp = gY_Fp; aTheta_Fp = gTh_Fp; aPhi_Fp = gPh_Fp; 
	aDp_Tg = gDp_Tg; aY_Tg = gY_Tg; aTheta_Tg = gTh_Tg; aPhi_Tg = gPh_Tg; aRctPt_Z = gRctPt_Z;
	break;
      }
    }
  }
  aAccCuts[0] = aX_Fp;  aAccCuts[1] = aY_Fp; aAccCuts[2] = aTheta_Fp; aAccCuts[3] = aPhi_Fp;
  aAccCuts[4] = aDp_Tg; aAccCuts[5] = aY_Tg; aAccCuts[6] = aTheta_Tg; aAccCuts[7] = aPhi_Tg;
  aAccCuts[8] = aRctPt_Z;                         
}

/*inline TString gGet_XGT2_Cut_SAMC(const TString& aArm,const TString& aTarget, const Double_t& aMom){{{*/
inline TString gGet_XGT2_Cut_SAMC(const TString& aArm,const TString& aTarget, const Double_t& aMom)
{
  TString Acc_Cut_Fp, Acc_Cut_Tg, Calo_Cuts, Cer_Cuts, Trigger_Cut;
  Double_t aX_Fp=0,aY_Fp=0,aTheta_Fp=0,aPhi_Fp=0;      // Focal Plane Cuts
  Double_t aDp_Tg=0,aY_Tg=0,aTheta_Tg=0,aPhi_Tg=0,aRctPt_Z=0;     // Target Plane Cuts
  
  Double_t gMom, gX_Fp, gY_Fp, gTh_Fp, gPh_Fp, gY_Tg, gTh_Tg, gPh_Tg, gDp_Tg,gRctPt_Z;
  TString gKin,gTgt;
  
  if(aArm == "L"){
      //Search the Cut Table
    ifstream InCut("Acc_Left_Table.dat");
    TString comment;
    comment.ReadLine(InCut);
    comment.ReadLine(InCut);
    while(!(InCut.eof())){
      InCut >> gMom >> gKin >> gTgt
	    >> gX_Fp >> gY_Fp >> gTh_Fp >> gPh_Fp
	    >> gDp_Tg >> gY_Tg >> gTh_Tg >> gPh_Tg >> gRctPt_Z;
//       cerr << "Mom = "<< gMom <<"  Old Mom = "<<aMom<<endl;
//       cerr << "Tgt = "<< gTgt <<"  Old Tgt = "<<aTarget.Data()<<endl;	

    if(aTarget.Data() == gTgt && (gMom - aMom) < 0.02 && (gMom - aMom) >-0.02){
	aX_Fp = gX_Fp;   aY_Fp = gY_Fp; aTheta_Fp = gTh_Fp; aPhi_Fp = gPh_Fp; 
	aDp_Tg = gDp_Tg; aY_Tg = gY_Tg; aTheta_Tg = gTh_Tg; aPhi_Tg = gPh_Tg; aRctPt_Z = gRctPt_Z;
	cerr <<"   === Applying Accpetance Cuts on SAMC Data: "<<endl
	     <<"   === Focal:   X_Fp = "<<aX_Fp <<", Y_Fp = "<<aY_Fp<<", Theta_Fp = "<<aTheta_Fp<<", Phi_Fp = "<<aPhi_Fp<<endl
	     <<"   === Target: Dp_Tg = "<<aDp_Tg<<", Y_Tg = "<<aY_Tg<<", Theta_Tg = "<<aTheta_Tg<<", Phi_Tg = "<<aPhi_Tg<<", RctPt_Z = "<<aRctPt_Z<<endl;
	break;
      }
    }
  }
       
  if(aArm == "R"){
    //Search the Cut Table
    ifstream InCut("Acc_Right_Table.dat");
    TString comment;
    comment.ReadLine(InCut);
    comment.ReadLine(InCut);
    while(!(InCut.eof())){
      InCut >> gMom >> gKin >> gTgt
	    >> gX_Fp >> gY_Fp >> gTh_Fp >> gPh_Fp
	    >> gDp_Tg >> gY_Tg >> gTh_Tg >> gPh_Tg >> gRctPt_Z;
     
      if(aTarget.Data() == gTgt.Data() && (gMom - aMom) < 0.002 && (gMom - aMom) >- 0.002){
	aX_Fp = gX_Fp;   aY_Fp = gY_Fp; aTheta_Fp = gTh_Fp; aPhi_Fp = gPh_Fp; 
	aDp_Tg = gDp_Tg; aY_Tg = gY_Tg; aTheta_Tg = gTh_Tg; aPhi_Tg = gPh_Tg; aRctPt_Z = gRctPt_Z;
	
	cerr <<"   ^^^ Applying Accpetance Cuts on SAMC Data: "<<endl
	     <<"   ^^^^Focal:   X_Fp = "<<aX_Fp <<" Y_Fp"<<aY_Fp<<" Theta_Fp"<<aTheta_Fp<<"Phi_Fp"<<aPhi_Fp<<endl
	     <<"   ^^^^Target: Dp_Tg = "<<aDp_Tg<<" Y_Tg"<<aY_Tg<<" Theta_Tg"<<aTheta_Tg<<"Phi_Tg"<<aPhi_Tg<<", RctPt_Z = "<<aRctPt_Z<<endl;
	break;
      }
    }
  }

    //Focal Plane Cuts
    Acc_Cut_Fp = Form("abs(x_fp)<=%4.3f && abs(y_fp)<=%4.3f && abs(th_fp)<=%4.3f && abs(ph_fp)<=%4.3f"
		      , aX_Fp, aY_Fp, aTheta_Fp, aPhi_Fp);
    //Target Plane Cuts
    Acc_Cut_Tg = Form("abs(dp_rec)<=%4.3f && abs(y_tg_rec)<=%4.3f && abs(th_tg_rec)<=%4.3f && abs(ph_tg_rec)<=%4.3f && abs(reactz_rec)<%4.3f&&IsPassed",
		       aDp_Tg, aY_Tg, aTheta_Tg, aPhi_Tg, aRctPt_Z);  

  TString aCuts = Acc_Cut_Fp+"&&"+ Acc_Cut_Tg;
  return aCuts;
}
/*}}}*/

inline TChain* gGet_SAMC_Chain(const TString& aSAMC_File)
{
  TChain *aTree = new TChain("SAMC");
  TString aSAMC_Chain = aSAMC_File.Data();
  aTree->Add(aSAMC_Chain.Data());
  for(int i=1;i<10;i++){
    aSAMC_Chain = aSAMC_Chain.ReplaceAll(Form("%d_xem.root",i-1),Form("%d_xem.root",i));
    cerr <<"       Add SAMC Files -> "<<aSAMC_Chain.Data()<<endl;
    aTree->Add(aSAMC_Chain.Data());
  }
  return aTree;
}
