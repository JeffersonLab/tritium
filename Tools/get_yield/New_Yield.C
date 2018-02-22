#include "TRI_Main.h"
#include "TRI_Tools.h"
#include "TRI_Beam.h"
#include "TRI_Target.h"
#include "TRI_Data.h"


int main(int argc, char* argv[]){
//void New_Yield(TString aInputFileName=""){
	const double xB[16] = {0.218, 0.258, 0.298, 0.338, 0.378, 0.418, 0.458, 0.498, 0.538, 0.578, 0.618, 0.658, 0.698,0.738, 0.778, 0.818 };
	//TString aInputFileName="";
    TString aInputFileName = argv[1];
    if(aInputFileName==""){
		cout<<" --- Input File = ";cin >> aInputFileName;
	}
	aInputFileName = "./Runlist/"+aInputFileName;
    
    	/*Read InputFile{{{*/
	//TString aInputFileName = ""; cout<<"--- Input File Name : "; cin >> aInputFileName;
	// int kin = 0; cout<<"--- Kin (1 or 2 or 3 ...) : "; cin >> kin;
	cout   << endl << Form("@@@@@@ Initialize Input values from %s ...",aInputFileName.Data())<<endl;
	Int_t aline=0;
	const TString Target=gGet_InputFile_Var(aInputFileName,aline);aline++; cout << "   Target Name = " << Target<<endl; 
	const int kin=gGet_InputFile_Var(aInputFileName,aline).Atoi();aline++;                cout << Form("   Kin = %d ", kin)<<endl; 
	TString Kin = Form("Kin%d", kin);

	TString Run_String = gGet_InputFile_Var(aInputFileName,aline);
	const vector<Int_t> RunNoChain=gGet_RunNoChain(Run_String);aline++;
	const int NumberOfRun = RunNoChain.size();
	cout << "   Run-List = " ;//<< Run_String<<endl;
    for(int i=0;i<NumberOfRun;i++){
        cout<<RunNoChain[i]<<" ";

    }
	cout<<endl;
	/*}}}*/

	TString Arm="";// = gGet_InputFile_Var(aInputFileName,aline);aline++; 
	if(RunNoChain[0]>50000)
		Arm = "R";
	else
		Arm = "L";

    /*Cuts{{{*/
	TString cut_L = Form("DL.evtypebits>>%d&1",MAIN_TRIGGER_LEFT);
	TString cut_e_L = Form("DL.evtypebits>>%d&1 && L.tr.n==1 && L.cer.asum_c>%4.1f && (L.prl1.e+L.prl2.e)/(L.gold.p*1000.)>%3.2f && (L.prl2.e)>%3.2f && (rpl.z)>%4.3f && (rpl.z)<%4.3f && (L.tr.tg_ph)>%4.3f && (L.tr.tg_ph)<%4.3f && (L.tr.tg_th)>%4.3f && (L.tr.tg_th)<%4.3f && (L.tr.tg_dp)>%4.3f && (L.tr.tg_dp)<%4.3f && evLeftdnew_r>%f",MAIN_TRIGGER_LEFT, GC_Cut, EP_Cut, EC2_Cut, TG_VZ_Min, TG_VZ_Max,TG_Phi_Min, TG_Phi_Max, TG_Theta_Min, TG_Theta_Max, TG_Dp_Min, TG_Dp_Max,BCM_TRIP_CUT);
	TString cut_pi_L = Form("DL.evtypebits>>%d&1 && L.tr.n==1 && L.cer.asum_c<%4.1f && (L.prl1.e+L.prl2.e)/(L.gold.p*1000.)<%3.2f && (rpl.z)>%4.3f && (rpl.z)<%4.3f && (L.tr.tg_ph)>%4.3f && (L.tr.tg_ph)<%4.3f && (L.tr.tg_th)>%4.3f && (L.tr.tg_th)<%4.3f && (L.tr.tg_dp)>%4.3f && (L.tr.tg_dp)<%4.3f && evLeftdnew_r>%f",MAIN_TRIGGER_LEFT, Pi_GC_Cut, Pi_EP_Cut, TG_VZ_Min, TG_VZ_Max,TG_Phi_Min, TG_Phi_Max, TG_Theta_Min, TG_Theta_Max, TG_Dp_Min, TG_Dp_Max,BCM_TRIP_CUT);

	TString cut_R = Form("DR.evtypebits>>%d&1",MAIN_TRIGGER_RIGHT);
	TString cut_e_R = Form("DR.evtypebits>>%d&1 && R.tr.n==1 && R.cer.asum_c>%4.1f && (R.ps.e+R.sh.e)/(R.gold.p*1000.)>%3.2f && (R.sh.e)>%3.2f && (rpr.z)>%4.3f && (rpr.z)<%4.3f && (R.gold.ph)>%4.3f && (R.gold.ph)<%4.3f && (R.gold.th)>%4.3f && (R.gold.th)<%4.3f && (R.gold.dp)>%4.3f && (R.gold.dp)<%4.3f && evRightdnew_r>%f",MAIN_TRIGGER_RIGHT, GC_Cut, 0.95, EC2_Cut, TG_VZ_Min, TG_VZ_Max,TG_Phi_Min, TG_Phi_Max, TG_Theta_Min, TG_Theta_Max, TG_Dp_Min, TG_Dp_Max, BCM_TRIP_CUT);
	TString cut_pi_R = Form("DR.evtypebits>>%d&1 && R.tr.n==1 && R.cer.asum_c<%4.1f && (R.ps.e+R.sh.e)/(R.gold.p*1000.)<%3.2f && (rpr.z)>%4.3f && (rpr.z)<%4.3f && (R.gold.ph)>%4.3f && (R.gold.ph)<%4.3f && (R.gold.th)>%4.3f && (R.gold.th)<%4.3f && (R.gold.dp)>%4.3f && (R.gold.dp)<%4.3f && evRightdnew_r>%f",MAIN_TRIGGER_RIGHT, Pi_GC_Cut, Pi_EP_Cut, TG_VZ_Min, TG_VZ_Max,TG_Phi_Min, TG_Phi_Max, TG_Theta_Min, TG_Theta_Max, TG_Dp_Min, TG_Dp_Max, BCM_TRIP_CUT);

	TString cut = "", cut_e="", cut_pi="";
	int Main_Trigger = 0;
	if(Arm=="L"){
		cut = cut_L;
		cut_e = cut_e_L;
		cut_pi = cut_pi_L;
		Main_Trigger = MAIN_TRIGGER_LEFT;
	}else if(Arm=="R"){
		cut = cut_R;
		cut_e = cut_e_R;
		cut_pi = cut_pi_R;
		Main_Trigger = MAIN_TRIGGER_RIGHT;
	}
/*}}}*/
  
    /*Get Electron and Pion Counts{{{*/
    cout<<"--- Good Electron Cut: "<< cut_e<<endl;
    cout<<"--- Good Pion Cut: "<< cut_pi<<endl;
    
    Double_t* NeleChain = new double[NumberOfRun];
    Double_t* NpiChain = new double[NumberOfRun];
	Double_t* NtotChain = new double[NumberOfRun];
	Double_t* LiveTimeChain = new double[NumberOfRun];
	Double_t* LiveTimeErrChain = new double[NumberOfRun];
    gGetGoodEvents(RunNoChain,Arm, cut_e, cut_pi, cut, NeleChain, NpiChain, NtotChain, LiveTimeChain, LiveTimeErrChain );
   
    int Total_Events= 0.0, Good_Ele=0, Good_Pi=0;
    for(int i=0;i<NumberOfRun;i++){
       Total_Events += NtotChain[i];
       Good_Ele += NeleChain[i];
       Good_Pi  += NpiChain[i];
    }
 
    cout<<endl; 
    cout<<Form("--- For %s at %s on HRS-%s:  ", Target.Data(), Kin.Data(), Arm.Data())<<endl;
    cout<<"--- Total Good Electrons from these runs are = "<<Good_Ele<<endl;
    cout<<endl; 
    cout<<Form("--- For %s at %s on HRS-%s:  ", Target.Data(), Kin.Data(), Arm.Data())<<endl;
    cout<<"--- Total Good Pions from these runs are = "<<Good_Pi<<endl;
    /*}}}*/

	/*Calclate Yield{{{*/
	///////////////////////////////////
	// Number of Incident particles
	///////////////////////////////////
	/*Get Ne{{{*/
	cout   << endl << "@@@@@@ Calculating Incident Particle ..."<<endl;	
	Double_t runtime;
	Double_t* NeChain = new double[NumberOfRun];
	Double_t* CurrentChain = new double[NumberOfRun];
	Double_t* RateChain = new double[NumberOfRun];
	TRI_VAR* Ne=gCal_Ne(RunNoChain,Arm, &runtime, NeChain, CurrentChain, RateChain);
	cout   <<Form("     Ne = %4.3e", Ne->Value)<<endl;	
	/*}}}end of Get Ne*/

	////////////////////////////////
	// Target Luminocity
	////////////////////////////////
	/*Get Ntg{{{*/
	cout   << endl << "@@@@@@ Calculating Target Luminocity ..."<<endl;

	int Target_A=0, Target_Z=0;
	double Target_Thickness=0.0, Target_Thickness_Err=0.0;
	int err = gGet_TargetInfo(Target.Data(), &Target_A, &Target_Z, &Target_Thickness, &Target_Thickness_Err);
	double Target_Length = (TG_VZ_Max - TG_VZ_Min)*100.; //meter to cm
	TRI_VAR* Ntg = new TRI_VAR();
	Double_t* NtgChain = new double[NumberOfRun];
	//Boiling Effect for long cryo-targets will be corrected in this subroutine
	Ntg= gCal_Ntg(RunNoChain,Arm.Data(),Target_A,Target_Z,Target_Thickness,Target_Thickness_Err,Target_Length, NtgChain);
	cout   <<Form("     Avg Ntg = %4.3e; Before Boiling Corrected = %4.3e",Ntg->Value,Target_Thickness/Target_A*Na)<<endl;
	/*}}}end of Get Ntg*/

	const double phase_space = (TG_Theta_Max - TG_Theta_Min) * (TG_Phi_Max - TG_Phi_Min) * (TG_Dp_Max - TG_Dp_Min) * P0; 
	double XS = Good_Ele / Ne->Value / Ntg->Value / phase_space * CM2ToNB;
	//double XS_Err = XS * sqrt( 1/Good_Ele + pow(Ntg->Sys_Err/Ntg->Value, 2) + pow(Ne->Sys_Err/Ne->Value, 2) );
	double XS_Err = XS * sqrt( 1/Good_Ele) ;
	cout<<endl<<"-----------------------------------------------------------"<<endl;
    cout<<"##### Raw Electron Cross Section = "<< XS <<" nb/GeV/sr"<<endl;

    double XS_Pi = Good_Pi / Ne->Value / Ntg->Value / phase_space * CM2ToNB;
	//double XS_Pi_Eerr = XS * sqrt( 1/Good_Pi + pow(Ntg->Sys_Err/Ntg->Value, 2) + pow(Ne->Sys_Err/Ne->Value, 2) );
	double XS_Pi_Err = XS_Pi * sqrt( 1/Good_Pi) ;
	cout<<"##### Raw Pion Cross Section = "<< XS_Pi <<" nb/GeV/sr"<<endl;
	cout<<"-----------------------------------------------------------"<<endl<<endl;
	/*}}}*/
    
    /*Run-By-Run XS{{{*/
	Double_t* XS_ele_Chain = new double[NumberOfRun];
	Double_t* XSErr_ele_Chain = new double[NumberOfRun];
    Double_t* XS_pi_Chain = new double[NumberOfRun];
    Double_t* XSErr_pi_Chain = new double[NumberOfRun];
    
    Double_t* R_ele = new double[NumberOfRun];
    Double_t* RErr_ele = new double[NumberOfRun];
    Double_t* R_pi = new double[NumberOfRun];
    Double_t* RErr_pi = new double[NumberOfRun];
    Double_t* RunN = new double[NumberOfRun];

    for(int i=0;i<NumberOfRun;i++){
        XS_ele_Chain[i]= NeleChain[i]/ NeChain[i] / NtgChain[i] / phase_space  * CM2ToNB;
        XSErr_ele_Chain[i] = XS_ele_Chain[i] * sqrt(1./NeleChain[i] );

        XS_pi_Chain[i] = NpiChain[i] / NeChain[i] / NtgChain[i] / phase_space  * CM2ToNB;
        XSErr_pi_Chain[i] = XS_ele_Chain[i] * sqrt(1./NpiChain[i] );

        if(i==0){
            R_ele[i] = 1.;
            RErr_ele[i] = 0.0;
            R_pi[i] = 1.0;
            RErr_pi[i] = 0.0;
        }else{
            R_ele[i] = XS_ele_Chain[i]/XS_ele_Chain[0];
            RErr_ele[i] =R_ele[i]* sqrt( 1./NeleChain[0] +  1./NeleChain[i] );
            R_pi[i] = XS_pi_Chain[i]/XS_pi_Chain[0];
            RErr_pi[i] = R_pi[i]*sqrt( 1./NpiChain[0] +  1./NpiChain[i] );
        }
        RunN[i] =i+1;
    }
    /*}}}*/
	
    TString Output_TXT = Form("results/%s_%s_%s_%d_T%d.dat", Target.Data(),Kin.Data(),Arm.Data(), RunNoChain[0],Main_Trigger);
    ofstream outf(Output_TXT.Data());

    outf <<Form("%6s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s", "RunNo", "N_tot", "N_ele", "N_pi", "LiveTime", "LT_Err", "Ntg", "Charge","Current", "Rate", "XS_ele","XSErr_ele","XS_pi","XSErr_pi")<<endl;
    for(int i=0;i<NumberOfRun;i++){
        outf <<Form("%6d %10.4e %10.4e %10.4e %10.3f %10.3f %10.4e %10.4e %10.3f %10.2f %10.4e %10.4e %10.4e %10.4e", 
                RunNoChain[i], NtotChain[i], NeleChain[i], NpiChain[i], LiveTimeChain[i], LiveTimeErrChain[i],NtgChain[i],  NeChain[i], CurrentChain[i], RateChain[i], XS_ele_Chain[i], XSErr_ele_Chain[i], XS_pi_Chain[i], XSErr_pi_Chain[i])<<endl;
    }

    outf<<endl;
    outf<<endl;
    outf<<"===================================================================================="<<endl;
    outf<<"XS Unit = nb/GeV/sr "<<endl;
    outf<<"===================================================================================="<<endl;
    outf<<Form( "%14s %14s %14s %14s %14s %14s %14s %14s  %14s","Total-Charge", "Total-Target", "Phase-Space", " Good-Ele", "XS-Ele", "XS-Ele Err", "Good-Pi", "XS-Pi", "XS-Pi Err")<<endl;
    outf<<Form( "%14.4e %14.4e %14.4e %14e %14.4e %14.4e %14e %14.4e %14.4e",Ne->Value, Ntg->Value, phase_space, (double)Good_Ele, XS, XS_Err, (double)Good_Pi, XS_Pi, XS_Pi_Err)<<endl;
    outf<<"===================================================================================="<<endl;
    outf<<endl;
    outf<<endl;
    outf<<"--- Good Electron Cut: "<< cut_e<<endl;
    outf<<"--- Good Pion Cut: "<< cut_pi<<endl;
    
    outf.close();

    /*Making Plots {{{*/
	TChain* T_Tree =(TChain*) gGetTree(RunNoChain, "T");
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TH1F *h1 = new TH1F("h1","",200, -0.2,0.2);
    h1->SetXTitle("ReactZ (m)");
    if(Arm=="L"){
        T_Tree->Draw("rpl.z>>h1",TCut(cut_L));
        T_Tree->SetLineColor(2);
        T_Tree->Draw("rpl.z",TCut(cut_e_L),"same");
    }
    if(Arm=="R"){
        T_Tree->Draw("rpr.z>>h1",TCut(cut_R));
        T_Tree->SetLineColor(2);
        T_Tree->Draw("rpr.z",TCut(cut_e_R),"same");
    }
    Int_t max = h1->GetMaximum();
    h1->GetYaxis()->SetRangeUser(0, max*1.1);

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    t1->SetTextFont(32);
    t1->SetTextSize(0.06);
    t1->SetTextColor(1); 
    t1->DrawLatex(0.15, 0.92, Form("Target = %s at %s (x=%4.3f) on HRS-%s ",Target.Data(), Kin.Data(),xB[kin-1], Arm.Data() ));

    t1->SetTextColor(4); 
    t1->SetTextSize(0.04);
    t1->DrawLatex(0.14, 0.85, Form("Total Good Electrons = %4.1fK (Total %4.1fK T%d Events) ", (Good_Ele/1000.), (Total_Events/1000.), Main_Trigger));

    t1->SetTextColor(4); 
    t1->SetTextSize(0.04);
    t1->DrawLatex(0.16, 0.80, Form("Electron Normalized Yield = %6.3e #pm %2.1f%% nb/GeV/sr", XS, XS_Err/XS*100.));

    t1->SetTextColor(2); 
    t1->SetTextSize(0.04);
    t1->DrawLatex(0.14, 0.70, Form("Total Good Pions = %4.1fK (Total %4.1fK T%d Events) ", (Good_Pi/1000.), (Total_Events/1000.), Main_Trigger));

    t1->SetTextColor(2); 
    t1->SetTextSize(0.04);
    t1->DrawLatex(0.16, 0.65, Form("Pion Normalized Yield = %6.3e #pm %2.1f%% nb/GeV/sr", XS_Pi, XS_Pi_Err/XS_Pi*100.));

    t1->SetTextColor(6); 
    t1->SetTextSize(0.04);
    t1->DrawLatex(0.2, 0.55, Form("Target-Lumi = %6.3e cm^-1",Ntg->Value));
    t1->DrawLatex(0.2, 0.50, Form("Beam-Charge = %6.3e",Ne->Value));

    TString Output_PNG = Form("results/%s_%s_%s_%d_T%d.png", Target.Data(),Kin.Data(),Arm.Data(), RunNoChain[0],Main_Trigger);
    TString Output_PDF = Form("results/%s_%s_%s_%d_T%d.pdf", Target.Data(),Kin.Data(),Arm.Data(), RunNoChain[0],Main_Trigger);
    c1->Print(Output_PNG.Data());
    c1->Print(Output_PDF.Data());
    /*}}}*/

    /*Making Run-By-Run Electron XS Plots{{{*/
    TCanvas *c2 = new TCanvas("c2","c2",800,600);
    TH2F *h2 = new TH2F("h2","",200, 0,NumberOfRun+1, 200, 0.0, 2.0);
    h2->SetXTitle("RunID");
    h2->SetYTitle("XS(RunID!=0)/XS(RunID==0)");
    h2->Draw();

    TGraphErrors *pl_r[NumberOfRun];
    int N_Color = 1, N_Style = 20;
    for(int i=0;i<NumberOfRun;i++){
        pl_r[i] = new TGraphErrors(1,&RunN[i],&R_ele[i],0,&RErr_ele[i]);
        pl_r[i]-> SetMarkerSize(1.5);
        N_Color++;
        if(N_Color>9&&N_Color<37)	N_Color = 37;
        if(N_Color>49)	N_Color = 1;
        if(N_Color==2) N_Color = 28;
        if(R_ele[i]<0.1) N_Color = 2;
        if(N_Color==5)		
            pl_r[i]-> SetMarkerColor(1);
        else
            pl_r[i]-> SetMarkerColor(N_Color);
        N_Style++;
        if(N_Style>30)	N_Style = 2;
        if(N_Style>5&&N_Style<20) N_Style = 20;
        pl_r[i]-> SetMarkerStyle(N_Style);
        pl_r[i]->Draw("P");
    }

    TLatex t2; N_Color = 1;
    t2.DrawLatex(0.2*(1),1.40,Form("%s@ %s on HRS-%s(%d Runs)", Target.Data(),Kin.Data(),Arm.Data(),NumberOfRun));
    for(int i=0;i<NumberOfRun;i++){
        t2.SetTextAngle(90);
        if(NumberOfRun>10)
            t2.SetTextSize(0.02);
        else
            t2.SetTextSize(0.03);
        N_Color++;
        if(N_Color>9&&N_Color<37)	N_Color = 37;
        if(N_Color>49)	N_Color = 1;
        if(N_Color==2) N_Color = 28;
        if(R_ele[i]<0.1) N_Color = 2;
        if(N_Color==5) 
            t2.SetTextColor(1);
        else
            t2.SetTextColor(N_Color);
        if(R_ele[i]<0.5 || R_ele[i]>1.5){
            t2.SetTextColor(2);
            t2.DrawLatex(i+1,0.91,Form("#%d (Bad?)", RunNoChain[i]));
        }
        else if((R_ele[i]+RErr_ele[i])>1.5)
            t2.DrawLatex(i+1,(R_ele[i]-RErr_ele[i])*0.95,Form("#%d", RunNoChain[i]));
        else
            t2.DrawLatex(i+1,(R_ele[i]+RErr_ele[i])*1.01,Form("#%d", RunNoChain[i]));
    }     
    TLine *line1= new TLine(0, 1, NumberOfRun+1, 1);
    line1->SetLineColor(7);
    line1->SetLineWidth(1.5);
    line1->SetLineStyle(9);
    line1->Draw("same");

    TString Output_RUN = Form("results/Runs_Elec_%s_%s_%s_%d_T%d.png", Target.Data(),Kin.Data(),Arm.Data(), RunNoChain[0],Main_Trigger);
    c2->Print(Output_RUN.Data());
    /*}}}*/
    
    /*Making Run-By-Run Pion XS Plots{{{*/
    TCanvas *c3 = new TCanvas("c3","c3",800,600);
    TH2F *h3 = new TH2F("h3","",200, 0,NumberOfRun+1, 200, 0.0, 2.0);
    h3->SetXTitle("RunID");
    h3->SetYTitle("XS_pi(RunID!=0)/XS_pi(RunID==0)");
    h3->Draw();

    TGraphErrors *pl_pi[NumberOfRun];
    N_Color = 1; N_Style = 20;
    for(int i=0;i<NumberOfRun;i++){
        pl_pi[i] = new TGraphErrors(1,&RunN[i],&R_pi[i],0,&RErr_pi[i]);
        pl_pi[i]-> SetMarkerSize(1.5);
        N_Color++;
        if(N_Color>9&&N_Color<37)	N_Color = 37;
        if(N_Color>49)	N_Color = 1;
        if(N_Color==2) N_Color = 28;
        if(R_pi[i]<0.1) N_Color = 2;
        if(N_Color==5)		
            pl_pi[i]-> SetMarkerColor(1);
        else
            pl_pi[i]-> SetMarkerColor(N_Color);
        N_Style++;
        if(N_Style>30)	N_Style = 2;
        if(N_Style>5&&N_Style<20) N_Style = 20;
        pl_pi[i]-> SetMarkerStyle(N_Style);
        pl_pi[i]->Draw("P");
    }

    TLatex t3; N_Color = 1;
    t3.DrawLatex(0.2*(1),1.40,Form("%s@ %s on HRS-%s(%d Runs)", Target.Data(),Kin.Data(),Arm.Data(),NumberOfRun));
    for(int i=0;i<NumberOfRun;i++){
        t3.SetTextAngle(90);
        if(NumberOfRun>10)
            t3.SetTextSize(0.02);
        else
            t3.SetTextSize(0.03);
        N_Color++;
        if(N_Color>9&&N_Color<37)	N_Color = 37;
        if(N_Color>49)	N_Color = 1;
        if(N_Color==2) N_Color = 28;
        if(R_pi[i]<0.1) N_Color = 2;
        if(N_Color==5) 
            t3.SetTextColor(1);
        else
            t3.SetTextColor(N_Color);
        if(R_pi[i]<0.5 || R_pi[i]>1.5){
            t3.SetTextColor(2);
            t3.DrawLatex(i+1,0.91,Form("#%d (Bad?)", RunNoChain[i]));
        }
        else if((R_pi[i]+RErr_pi[i])>1.5)
            t3.DrawLatex(i+1,(R_pi[i]-RErr_pi[i])*0.95,Form("#%d", RunNoChain[i]));
        else
            t3.DrawLatex(i+1,(R_pi[i]+RErr_pi[i])*1.01,Form("#%d", RunNoChain[i]));
    }     
    TLine *line2= new TLine(0, 1, NumberOfRun+1, 1);
    line2->SetLineColor(7);
    line2->SetLineWidth(1.5);
    line2->SetLineStyle(9);
    line2->Draw("same");

    TString Output_RUN2 = Form("results/Runs_Pion_%s_%s_%s_%d_T%d.png", Target.Data(),Kin.Data(),Arm.Data(), RunNoChain[0],Main_Trigger);
    c3->Print(Output_RUN2.Data());
    /*}}}*/
}
