/*
  Weight some variables by cross section.
  -- Zhihong Ye, 08/08/2011
*/

#include "XGT2.h"

//int main()
void refill_full()
{
  TString filename; 
  cerr <<"Which RootFiles? "; cin >> filename;

  gSystem->Exec("cp -vf "+filename+" "+filename+".old");

  TFile *f1 = new TFile(filename,"update");
  assert(f1);

  // hxfp->Delete(); hyfp->Delete(); hthfp->Delete(); hphfp->Delete();
  // hxtg->Delete(); hytg->Delete(); hthtg->Delete(); hphtg->Delete(); hdptg->Delete();

  TTree *t1 =(TTree*)f1->Get("SAMC");  assert(t1);
  
  /***********************************************************/
  /**               Set BranchAddress                       **/
  /***********************************************************/
  t1->SetBranchAddress("Id",&Id); //Event Id
  t1->SetBranchAddress("IsPassed",&IsPassed); //if pass through the magnet, 1=true, 0=false(not pass)
  t1->SetBranchAddress("IsQualified",&IsQualified); //if rec var is in the range of gen, 1=true, 0=false(not
  
  /*only for RadCor{{{*/
  t1->SetBranchAddress("mp_cor",&mp_cor   ); //Multiple Photon Correction
  t1->SetBranchAddress("tail",&tail   ); //radiative tail
  t1->SetBranchAddress("R",&R   ); //Phys.Rev.D 12 1884 A83
 
  t1->SetBranchAddress("cs_qfs",&cs_qfs   );//quasi-free scattering cross section
  t1->SetBranchAddress("cs_del",&cs_del   );//Delta
  t1->SetBranchAddress("cs_x",&cs_x   );  //DIS
  t1->SetBranchAddress("cs_r1500",&cs_r1500   ); //Resonace 1500 MeV
  t1->SetBranchAddress("cs_r1700",&cs_r1700   ); //Resonace 1700 MeV
  t1->SetBranchAddress("cs_2N",&cs_2N   ); //Dip region
  t1->SetBranchAddress("cs_q",&cs_q   );//sum of all the cross section above
 
  t1->SetBranchAddress("cs_p",&cs_p   ); //sigma_p (peak approximation for internal radation)
  t1->SetBranchAddress("cs_b_bt",&cs_b_bt   ); //the part before target in sigma_b
  t1->SetBranchAddress("cs_b_at",&cs_b_at   ); //the part after target in sigma_b
  t1->SetBranchAddress("cs_b",&cs_b   ); //sigma_b (external bremsstrahlung ionization)
  t1->SetBranchAddress("cs_ex",&cs_ex   ); //sigma_ex (extract internal radation)
  t1->SetBranchAddress("cs_qt",&cs_qt   );
  /*}}}*/
  
  /*Member Data from random{{{*/
  t1->SetBranchAddress("beam_x",&beam_x    ); //cm
  t1->SetBranchAddress("beam_y",&beam_y    ); //cm
  t1->SetBranchAddress("reactz_gen",&reactz_gen    ); // interaction z position cm
  t1->SetBranchAddress("E_s",&E_s    );//=incident beam energy (MeV)
  t1->SetBranchAddress("E_p",&E_p    );//=incident beam energy (MeV)
  //gen=generator
  t1->SetBranchAddress("th_tg_gen",&th_tg_gen    ); //theta target
  t1->SetBranchAddress("ph_tg_gen",&ph_tg_gen    ); //phi target
  t1->SetBranchAddress("dp_gen",&dp_gen    ); //dp at target
  /*}}}*/
  
  /*Member Data derived from variables above{{{*/
  t1->SetBranchAddress("x_tg_gen",&x_tg_gen    ); //
  t1->SetBranchAddress("y_tg_gen",&y_tg_gen    ); //
  /*}}}*/
  
  /*Member Data for refinement{{{*/
  t1->SetBranchAddress("dp_ref",&dp_ref    ); //dp for John.LeRose matrix
  t1->SetBranchAddress("x_tg_ref",&x_tg_ref    ); //x tg for John.LeRose matrix
  t1->SetBranchAddress("y_tg_ref",&y_tg_ref    ); //y tg for John.LeRose matrix
  t1->SetBranchAddress("th_tg_ref",&th_tg_ref    ); //th tg for John.LeRose matrix
  t1->SetBranchAddress("ph_tg_ref",&ph_tg_ref    ); //ph tg for John.LeRose matrix
  /*}}}*/
  
  /*Member Data for focal plane {{{*/
  t1->SetBranchAddress("x_fp",&x_fp    ); //
  t1->SetBranchAddress("y_fp",&y_fp    ); //
  t1->SetBranchAddress("th_fp",&th_fp    ); //
  t1->SetBranchAddress("ph_fp",&ph_fp    ); //
  t1->SetBranchAddress("th_fp_no_ort",&th_fp_no_ort    ); //th_fp without TXFIT orthogonalization
 
  t1->SetBranchAddress("q1ex_x",&q1ex_x    );//0:x 1:y
  t1->SetBranchAddress("dent_x",&dent_x    );
  t1->SetBranchAddress("dext_x",&dext_x    );
  t1->SetBranchAddress("q3en_x",&q3en_x    );
  t1->SetBranchAddress("q3ex_x",&q3ex_x    );
 
  t1->SetBranchAddress("q1ex_y",&q1ex_y    );//0:x 1:y
  t1->SetBranchAddress("dent_y",&dent_y    );
  t1->SetBranchAddress("dext_y",&dext_y    );
  t1->SetBranchAddress("q3en_y",&q3en_y    );
  t1->SetBranchAddress("q3ex_y",&q3ex_y    );
  
  t1->SetBranchAddress("IsPassedQ1Ex",&IsPassedQ1Ex    );
  t1->SetBranchAddress("IsPassedDipoleEn",&IsPassedDipoleEn    );
  t1->SetBranchAddress("IsPassedDipoleEx",&IsPassedDipoleEx    );
  t1->SetBranchAddress("IsPassedQ3En",&IsPassedQ3En    );
  t1->SetBranchAddress("IsPassedQ3Ex",&IsPassedQ3Ex    );
  /*}}}*/
  
  /*Member Data for reconstruct target varibles {{{*/
  t1->SetBranchAddress("x_tg_rec",&x_tg_rec    ); //cm
  t1->SetBranchAddress("y_tg_rec",&y_tg_rec    ); //cm
  t1->SetBranchAddress("th_tg_rec",&th_tg_rec    ); //
  t1->SetBranchAddress("ph_tg_rec",&ph_tg_rec    ); //
  t1->SetBranchAddress("dp_rec",&dp_rec    );
  t1->SetBranchAddress("reactz_rec",&reactz_rec    );
  t1->SetBranchAddress("rvalue",&rvalue    );

  t1->SetBranchAddress("cs_M",&cs_M    );//mott cross section
  t1->SetBranchAddress("cs_Final",&cs_Final    );//final cross section
  
  // /*Weight Variables on Focal Plane and Target Plane{{{*/
  // Double_t x_fp_w,y_fp_w,th_fp_w,ph_fp_w;
  // Double_t x_tg_w,y_tg_w,th_tg_w,ph_tg_w,dp_tg_w;
  // t1->SetBranchAddress("x_fp_w",&x_fp_w);
  // t1->SetBranchAddress("y_fp_w",&y_fp_w);
  // t1->SetBranchAddress("th_fp_w",&th_fp_w);
  // t1->SetBranchAddress("ph_fp_w",&ph_fp_w);
  
  // t1->SetBranchAddress("x_tg_w",&x_tg_w);
  // t1->SetBranchAddress("y_tg_w",&y_tg_w);
  // t1->SetBranchAddress("th_tg_w",&th_tg_w);
  // t1->SetBranchAddress("ph_tg_w",&ph_tg_w);
  // t1->SetBranchAddress("dp_tg_w",&dp_tg_w);
  // /*}}}*/
  
/***********************************************************/
/*       End of  Setting BranchAddress                     */
/***********************************************************/

/***********************************************************/
/**                 Add new  Branch                       **/
/***********************************************************/
  // Double_t x_fp_w,y_fp_w,th_fp_w,ph_fp_w;
  // Double_t x_tg_w,y_tg_w,th_tg_w,ph_tg_w,dp_tg_w;
  
  // t1->Branch("x_fp_w",&x_fp_w,"x_fp_w/D");
  // t1->Branch("y_fp_w",&y_fp_w,"y_fp_w/D");
  // t1->Branch("th_fp_w",&th_fp_w,"th_fp_w/D");
  // t1->Branch("ph_fp_w",&ph_fp_w,"ph_fp_w/D");
  
  // t1->Branch("x_tg_w",&x_tg_w,"x_tg_w/D");
  // t1->Branch("y_tg_w",&y_tg_w,"y_tg_w/D");
  // t1->Branch("th_tg_w",&th_tg_w,"tg_fp_w/D");
  // t1->Branch("ph_tg_w",&ph_tg_w,"ph_tg_w/D");
  // t1->Branch("dp_tg_w",&dp_tg_w,"dp_tg_w/D");
  
/***********************************************************/
/**            End of  Adding new  Branch                 **/
/***********************************************************/
  
  //TChain *S = new TChain("SAMC");
  //S->Add("c12_L_xgt2_samc.root");
  //SetTreeb(S);

  Int_t Evt = t1->GetEntries();
  double cs_Average=0;
  int count;

  for(int i=0;i<Evt;i++){
    t1->GetEntry(i);
    if(IsPassed&&cs_Final>0.0){
      cs_Average+=cs_Final*1000.;
      //     cerr << "Final XS = "<<cs_Final<<endl;
      count++;
    }
  }
  
  cs_Average/=count;
  cerr << "Average Cross Section is "<<cs_Average<<endl;
  
  TH1F *hxfp = new TH1F("hxfp","Weighted X_fp",200,-0.6,0.6);
  TH1F *hyfp = new TH1F("hyfp","Weighted Y_fp",200,-0.06,0.06);
  TH1F *hthfp = new TH1F("hthfp","Weighted Theta_fp",200,-0.06,0.06);
  TH1F *hphfp = new TH1F("hphfp","Weighted Phi_fp",200,-0.05,0.05);

  TH1F *hxtg = new TH1F("hxtg","Weighted X_tg",200,-0.01,0.01);
  if(filename.Contains("he3")||filename.Contains("h2")||filename.Contains("he4")){
    TH1F *hytg = new TH1F("hytg","Weighted Y_tg",200,-0.07,0.07);
    TH1F *hrpz = new TH1F("hrpz","Weighted RctPnt_Z",200,-0.12,0.12);
  }
  else{
    TH1F *hytg = new TH1F("hytg","Weighted Y_tg",200,-0.01,0.01);
    TH1F *hrpz = new TH1F("hrpz","Weighted RctPnt_Z",200,-0.01,0.01);
  }
  TH1F *hthtg = new TH1F("hthtg","Weighted Theta_tg",200,-0.06,0.06);
  TH1F *hphtg = new TH1F("hphtg","Weighted Phi_tg",200,-0.05,0.05);
  TH1F *hdptg = new TH1F("hdptg","Weighted Dp_tg",200,-0.06,0.06);

  double weight=1.0;

  for(int j=0;j<Evt;j++){
    t1->GetEntry(j);
    if(IsPassed&&cs_Final>0.0){

      // if(abs(y_tg_rec)<0.025
      // 	 &&abs(th_tg_rec)<0.04
      // 	 &&abs(ph_tg_rec)<0.02
      // 	 &&abs(dp_rec)<0.03){	
      if(1
	 //	 &&abs(y_tg_rec)<0.25
	 //	 &&abs(th_tg_rec)<0.04
	 //	 &&abs(ph_tg_rec)<0.02
      	 &&abs(reactz_rec)<0.07
	 //	 &&abs(dp_rec)<0.035	
	 ){
	weight = cs_Final*1000./cs_Average;
	
	hxfp->Fill(x_fp,weight);
	hyfp->Fill(y_fp,weight);
	hthfp->Fill(th_fp,weight);
	hphfp->Fill(ph_fp,weight);
	
	hxtg->Fill(x_tg_rec,weight);
	hytg->Fill(y_tg_rec,weight);
	hthtg->Fill(th_tg_rec,weight);
	hphtg->Fill(ph_tg_rec,weight);
	hdptg->Fill(dp_rec,weight);
	hrpz ->Fill(reactz_rec,weight);

	// x_fp_w = x_fp/weight;      y_fp_w = y_fp/weight;  
	// th_fp_w = th_fp/weight;    ph_fp_w = ph_fp/weight;  
	
	// x_tg_w = x_tg_rec/weight;      y_tg_w = y_tg_rec/weight;  
	// th_tg_w = th_tg_rec/weight;    ph_tg_w = ph_tg_rec/weight;  
	// dp_tg_w = dp_rec/weight;  
	
	// t1->Fill();
	  }
    }
  }

  // TCanvas *c1 = new TCanvas("c1","c1",800,800);
  // c1->Divide(2,2);
  // c1->cd(1);hxfp->Draw();c1->cd(2);hyfp->Draw();
  // c1->cd(3);hthfp->Draw();c1->cd(4);hphfp->Draw();
  
  // TCanvas *c2 = new TCanvas("c2","c2",800,800);
  // c2->Divide(2,2);
  // c2->cd(1);hdptg->Draw();c2->cd(2);hytg->Draw();
  // c2->cd(3);hthtg->Draw();c2->cd(4);hphtg->Draw();
  
  hxfp->Write(); hyfp->Write(); hthfp->Write(); hphfp->Write();
  hxtg->Write(); hytg->Write(); hthtg->Write(); hphtg->Write(); hdptg->Write();
  hrpz->Write();
  //  t1->Write();
  f1->Close();

}

