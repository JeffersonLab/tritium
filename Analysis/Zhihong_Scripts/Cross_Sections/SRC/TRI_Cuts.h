/////////////////////////////////
//  Cuts 
/////////////////////////////////
//Cuts on both Acceptance and PID
//Needed large modification depending on experiments
/*inline TString gGet_TRI_Cut(const Int_t& aRunNo,const TString& aArm,const TString& aTarget, const TString& aKin, const Double_t& aAngle, const Double_t& aMom,const Int_t& aMain_Trigger,const TString& aAccCut_FileName, const bool aIsExtTgt){{{*/
inline TString gGet_TRI_Cut(const Int_t& aRunNo,const TString& aArm,const TString& aTarget, const TString& aKin, const Double_t& aAngle, const Double_t& aMom,const Int_t& aMain_Trigger,const TString& aAccCut_FileName, const bool aIsExtTgt, double* aTarget_Length_Cut)
{
    TString Acc_Cut_Fp, Acc_Cut_Tg, Calo_Cuts, Cer_Cuts, Trigger_Cut;
    Double_t aX_Fp=0,aY_Fp=0,aTheta_Fp=0,aPhi_Fp=0;      // Focal Plane Cuts
    Double_t aDp_Tg=0,aY_Tg=0,aTheta_Tg=0,aPhi_Tg=0,aRctPt_Z=0;     // Target Plane Cuts
    Double_t aCalo_Cuts=0,aCalo_2nd=0,aCer_Cuts=0; //PID Cuts 
    Double_t aZ_Offset = 0.0, aY_Offset = 0.0;
    Double_t gMom, gX_Fp, gY_Fp, gTh_Fp, gPh_Fp, gY_Tg, gTh_Tg, gPh_Tg, gDp_Tg,gRctPt_Z;
    TString gKin,gTgt;

    if(aRunNo == 1111) // Reserve this for special runs
    {
        outlog << "Not in Use so far" <<endl;
        cerr << "Not in Use so far" <<endl;
    }

    /*If("L"){{{*/
    if(aArm == "L")
    {
        aCalo_Cuts = 0.5;    // E/P Cuts
        aCalo_2nd = 100.0;   // PRL2 Cuts
        aCer_Cuts = 50.0;    // Cer Cuts
        aY_Offset = 0.0;

        Trigger_Cut = Form("%s.edtpl[0]==0&&((%s.evtypebits>>%d)&1)",D_Name.Data(),D_Name.Data(),aMain_Trigger);
        Bool_t bFind=kFALSE; 
        //Search the Cut Table
        ifstream InCut(aAccCut_FileName.Data());
        TString comment;
        comment.ReadLine(InCut);
        comment.ReadLine(InCut);
        while(!(InCut.eof())){
            InCut >> gMom >> gKin >> gTgt
                >> gX_Fp >> gY_Fp >> gTh_Fp >> gPh_Fp
                >> gDp_Tg >> gY_Tg >> gTh_Tg >> gPh_Tg >> gRctPt_Z;
            //			cerr << "Mom = "<< gMom <<"  Old Mom = "<<aMom<<endl;
            //			cerr << "Tgt = "<< gTgt <<"  Old Tgt = "<<aTarget.Data()<<endl;	
            TString bKin ="Kin"+gKin;
            if(aTarget.Data() == gTgt && fabs(gMom - aMom) < 0.02){
                aX_Fp = gX_Fp;   aY_Fp = gY_Fp; aTheta_Fp = gTh_Fp; aPhi_Fp = gPh_Fp; 
                aDp_Tg = gDp_Tg; aY_Tg = gY_Tg; aTheta_Tg = gTh_Tg; aPhi_Tg = gPh_Tg; 
                aRctPt_Z = gRctPt_Z; //The Cut is symmetic but the target center has an offset

                //The optics determins that Y_tg does not change too much on diff angles but ReactionZ does
                //So we use the range of VZ and scale it by the angle, then apply the value on the Ytg cut
                // meanwhile, the ReactionZ cut is disable	-- Zhihong Ye, 02/02/2013
                aY_Tg = fabs(aRctPt_Z * sin(aAngle*TMath::DegToRad())); 
                aY_Offset = -(aZ_Offset) * sin(aAngle*TMath::DegToRad()); 
                *aTarget_Length_Cut = aRctPt_Z;

                outlog <<"   === Applying Accpetance Cuts: "<<endl
                    <<"   === Focal:   X_Fp = "<<aX_Fp <<", Y_Fp = "<<aY_Fp<<", Theta_Fp = "<<aTheta_Fp<<", Phi_Fp = "<<aPhi_Fp<<endl
                    <<"   === Target: Dp_Tg = "<<aDp_Tg<<", Y_Tg = "<<aY_Tg<<", Theta_Tg = "<<aTheta_Tg<<", Phi_Tg = "<<aPhi_Tg<<", RctPt_Z = "<<aRctPt_Z<<endl;
                cerr <<"   === Applying Accpetance Cuts: "<<endl
                    <<"   === Focal:   X_Fp = "<<aX_Fp <<", Y_Fp = "<<aY_Fp<<", Theta_Fp = "<<aTheta_Fp<<", Phi_Fp = "<<aPhi_Fp<<endl
                    <<"   === Target: Dp_Tg = "<<aDp_Tg<<", Y_Tg = "<<aY_Tg<<", Theta_Tg = "<<aTheta_Tg<<", Phi_Tg = "<<aPhi_Tg<<", RctPt_Z = "<<aRctPt_Z<<endl;
                bFind=kTRUE;
                break;
            }
            }
            if(!bFind){	
                outlog <<"***ERROR***, Can not obtain cuts for Left Arm!"<<endl;
                cerr <<"***ERROR***, Can not obtain cuts for Left Arm!"<<endl;
            }

            //Cherenkov Cut
            Cer_Cuts = Form("L.cer.asum_c>=%6.2f",aCer_Cuts);
            //Calorimeter Cuts
            Calo_Cuts = Form("(L.prl1.e+L.prl2.e)/L.tr.p/1000.0>=%f && L.prl2.e>%f",aCalo_Cuts,aCalo_2nd);
            //Calo_Cuts = Form("epL>=%6.5f && L.prl2.e>=%6.2f",aCalo_Cuts,aCalo_2nd);
            //Focal Plane Cuts
            Acc_Cut_Fp = Form("L.tr.n==1 && abs(L.tr.x)<=%6.5f && abs(L.tr.y)<=%6.5f && abs(L.tr.th)<=%6.5f && abs(L.tr.ph)<=%6.5f", aX_Fp, aY_Fp, aTheta_Fp, aPhi_Fp);
            //Target Plane Cuts
            TString aEX_Dp= Form("%s.gold.dp",aArm.Data());
            TString aEX_Th= Form("%s.gold.th",aArm.Data());
            TString aEX_Ph= Form("%s.gold.ph",aArm.Data());
            TString aEX_Y = Form("%s.gold.y",aArm.Data());
            if(aIsExtTgt){
                aEX_Dp= Form("ex%s.dp",aArm.Data());
                aEX_Th= Form("ex%s.th",aArm.Data());
                aEX_Ph= Form("ex%s.ph",aArm.Data());
                aEX_Y = Form("ex%s.y",aArm.Data());
            }
            Acc_Cut_Tg = Form("abs(%s)<=%6.5f && abs(%s-%6.5f)<=%6.5f && abs(%s)<=%6.5f && abs(%s)<=%6.5f", aEX_Dp.Data(),aDp_Tg, aEX_Y.Data(),aY_Offset,aY_Tg, aEX_Th.Data(),aTheta_Tg,aEX_Ph.Data(),aPhi_Tg);
        }
        /*If("L")}}}*/

        /*If("R"){{{*/
        if(aArm == "R")
        {
            aCalo_Cuts = 0.5;    // E/P Cuts
            aCalo_2nd = 200.0;   // PRL2 Cuts
            aCer_Cuts = 50.0;    // Cer Cuts
            aY_Offset = 0.0;

            //I assume we alway care about Trigger 1 and 3, not 6&7, which is danger but OK so far.
            Trigger_Cut = Form("%s.edtpl[0]==0&&((%s.evtypebits>>%d)&1)",D_Name.Data(),D_Name.Data(),aMain_Trigger);
            Bool_t bFind=kFALSE; 
            //Search the Cut Table
            ifstream InCut(aAccCut_FileName.Data());
            TString comment;
            comment.ReadLine(InCut);
            comment.ReadLine(InCut);
            while(!(InCut.eof())){
                InCut >> gMom >> gKin >> gTgt
                    >> gX_Fp >> gY_Fp >> gTh_Fp >> gPh_Fp
                    >> gDp_Tg >> gY_Tg >> gTh_Tg >> gPh_Tg >> gRctPt_Z;
                //			cerr << "Mom = "<< gMom <<"  Old Mom = "<<aMom<<endl;
                //			cerr << "Tgt = "<< gTgt <<"  Old Tgt = "<<aTarget.Data()<<endl;	
                TString bKin ="Kin"+gKin;
                if(aTarget.Data() == gTgt && fabs(gMom - aMom) < 0.02){
                    //if(aTarget.Data() == gTgt.Data() && aKin.Data()==bKin.Data()){
                    aX_Fp = gX_Fp;   aY_Fp = gY_Fp; aTheta_Fp = gTh_Fp; aPhi_Fp = gPh_Fp; 
                    aDp_Tg = gDp_Tg; aY_Tg = gY_Tg; aTheta_Tg = gTh_Tg; aPhi_Tg = gPh_Tg;
                    aRctPt_Z = gRctPt_Z; //The Cut is symmetic but the target center may have an offset

                    //The optics determins that Y_tg does not change too much on diff angles but ReactionZ does
                    //So we use the range of VZ and scale it by the angle, then apply the value on the Ytg cut
                    // meanwhile, the ReactionZ cut is disable	-- Zhihong Ye, 02/02/2013
                    aY_Tg =  fabs(aRctPt_Z * sin((-aAngle)*TMath::DegToRad())); //Angle is - on HRS-R
                    aY_Offset = -(aZ_Offset) * sin((-aAngle)*TMath::DegToRad()); 
                    *aTarget_Length_Cut = aRctPt_Z;

                    outlog <<"   === Applying Accpetance Cuts: "<<endl
                        <<"   === Focal:   X_Fp = "<<aX_Fp <<", Y_Fp = "<<aY_Fp<<", Theta_Fp = "<<aTheta_Fp<<", Phi_Fp = "<<aPhi_Fp<<endl
                        <<"   === Target: Dp_Tg = "<<aDp_Tg<<", Y_Tg = "<<aY_Tg<<", Theta_Tg = "<<aTheta_Tg<<", Phi_Tg = "<<aPhi_Tg<<", RctPt_Z = "<<aRctPt_Z<<endl;
                    cerr <<"   === Applying Accpetance Cuts: "<<endl
                        <<"   === Focal:   X_Fp = "<<aX_Fp <<", Y_Fp = "<<aY_Fp<<", Theta_Fp = "<<aTheta_Fp<<", Phi_Fp = "<<aPhi_Fp<<endl
                        <<"   === Target: Dp_Tg = "<<aDp_Tg<<", Y_Tg = "<<aY_Tg<<", Theta_Tg = "<<aTheta_Tg<<", Phi_Tg = "<<aPhi_Tg<<", RctPt_Z = "<<aRctPt_Z<<endl;
                    bFind=kTRUE;
                    break;
                }
                }
                if(!bFind){	
                    outlog <<"***ERROR***, Can not obtain cuts for Right Arm!"<<endl;
                    cerr <<"***ERROR***, Can not obtain cuts for Right Arm!"<<endl;
                }

                //Cherenkov Cut
                Cer_Cuts = Form("R.cer.asum_c>=%6.2f",aCer_Cuts);
                //Calorimeter Cuts
                //Calo_Cuts = Form("epR>=%6.5f && R.sh.e>=%6.2f",aCalo_Cuts,aCalo_2nd);
                Calo_Cuts = Form("(R.ps.e+R.sh.e)/R.tr.p/1000.0>=%f && R.sh2.e>%f",aCalo_Cuts,aCalo_2nd);
                //Focal Plane Cuts
                Acc_Cut_Fp = Form("R.tr.n==1 && abs(R.tr.x)<=%6.5f && abs(R.tr.y)<=%6.5f && abs(R.tr.th)<=%6.5f && abs(R.tr.ph)<=%6.5f", aX_Fp, aY_Fp, aTheta_Fp, aPhi_Fp);
                //Target Plane Cuts
                TString aEX_Dp= Form("%s.gold.dp",aArm.Data());
                TString aEX_Th= Form("%s.gold.th",aArm.Data());
                TString aEX_Ph= Form("%s.gold.ph",aArm.Data());
                TString aEX_Y = Form("%s.gold.y",aArm.Data());
                if(aIsExtTgt){
                    aEX_Dp= Form("ex%s.dp",aArm.Data());
                    aEX_Th= Form("ex%s.th",aArm.Data());
                    aEX_Ph= Form("ex%s.ph",aArm.Data());
                    aEX_Y = Form("ex%s.y",aArm.Data());
                }
                Acc_Cut_Tg = Form("abs(%s)<=%6.5f && abs(%s-%6.5f)<=%6.5f && abs(%s)<=%6.5f && abs(%s)<=%6.5f", aEX_Dp.Data(),aDp_Tg, aEX_Y.Data(),aY_Offset,aY_Tg, aEX_Th.Data(),aTheta_Tg,aEX_Ph.Data(),aPhi_Tg);
            }
            /*If("R")}}}*/

            //Get the true target length after cuts. We need to treat differently if we cut on different values on up- and down- stream;
            Target_Length = 2.0*aRctPt_Z;

            TString aCuts = Trigger_Cut+"&&"+Cer_Cuts+"&&"+Calo_Cuts+"&&"+Acc_Cut_Fp+"&&"+ Acc_Cut_Tg;
            return aCuts;
}
/*}}}*/

//Only Acceptance Cuts, not PID included, mainly used for MC data 
/*inline void gGet_TRI_Acc_Cut(const TString& aArm,const TString& aTarget, const Double_t& aMom, Double_t aAccCuts[9],const TString& aAccCut_FileName){{{*/
inline void gGet_TRI_Acc_Cut(const TString& aArm,const TString& aTarget, const Double_t& aMom, Double_t aAccCuts[9],const TString& aAccCut_FileName){
    Double_t aX_Fp=0,aY_Fp=0,aTheta_Fp=0,aPhi_Fp=0;      // Focal Plane Cuts
    Double_t aDp_Tg=0,aY_Tg=0,aTheta_Tg=0,aPhi_Tg=0,aRctPt_Z=0;     // Target Plane Cuts

    Double_t gMom, gX_Fp, gY_Fp, gTh_Fp, gPh_Fp, gY_Tg, gTh_Tg, gPh_Tg, gDp_Tg,gRctPt_Z;
    TString gKin,gTgt;
    
    ifstream InCut(aAccCut_FileName.Data());
    TString comment;
    comment.ReadLine(InCut);
    comment.ReadLine(InCut);
    //Search the Cut Table
    while(!(InCut.eof())){
        InCut >> gMom >> gKin >> gTgt
            >> gX_Fp >> gY_Fp >> gTh_Fp >> gPh_Fp
            >> gDp_Tg >> gY_Tg >> gTh_Tg >> gPh_Tg >> gRctPt_Z;

        if(aTarget.Data() == gTgt && fabs(gMom - aMom) < 0.02 ){
            aX_Fp = gX_Fp;   aY_Fp = gY_Fp; aTheta_Fp = gTh_Fp; aPhi_Fp = gPh_Fp; 
            aDp_Tg = gDp_Tg; aY_Tg = gY_Tg; aTheta_Tg = gTh_Tg; aPhi_Tg = gPh_Tg; aRctPt_Z = gRctPt_Z;
            break;
        }
    }

    aAccCuts[0] = aX_Fp;  aAccCuts[1] = aY_Fp; aAccCuts[2] = aTheta_Fp; aAccCuts[3] = aPhi_Fp;
    aAccCuts[4] = aDp_Tg; aAccCuts[5] = aY_Tg; aAccCuts[6] = aTheta_Tg; aAccCuts[7] = aPhi_Tg;
    aAccCuts[8] = aRctPt_Z;                         
}/*}}}*/
