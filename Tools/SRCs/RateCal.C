#include "GetRate.h"
#include "XEMC.h" //XEMC returns XS in nb/(sr*MeV) 

/*GetRate{{{*/
void GetRate(double E0, double P0, double T0, 
	    double* xb_avg, double* W2_avg, 
            double* xs_he3,double* xs_he3_avg,double* rate_he3, 
            double* xs_h3, double* xs_h3_avg,double* rate_h3, 
            double* xs_d2,double* xs_d2_avg,double* rate_d2, double *rate_window){

	bool IsDebug = false;
	//bool IsDebug = true;
	/*
	   double E0 = 11.0; //GeV
	   double P0; //GeV, Scattered Electrons' Central Momentum
	   double T0; //Degree, Scattered Electrons' Central Angle
	   cout<<" --- Beam Energy (E0, GeV) = "; cin >> E0;
	   cout<<" --- Scattered Electron's Central Momentum (P0, GeV) = "; cin >> P0;
	   cout<<" --- Scattered Electron's Central Angle (Theta0, Deg) = "; cin >> T0;
	   */
	int Nevnt = 1000000;
	//int Nevnt = 1000;
	if(IsDebug){cout<<" ----Total Samples = "; cin >> Nevnt;}
	cout<<endl<<endl;

    XEMC* xem = new XEMC();
    xem->Init("input/He3_Input.dat"); 

    Event* evt = new Event();     
    evt->SetSpectrometer("HRS");
    double Phase_Space = evt->GetPhaseSpace() * P0; //sr*GeV //it was deltaTheta * deltaPhi *deltaDp

    /*Define Target Info{{{*/
    double Beam_Current = 0.0;//A

    //Upstream Window:   0.010" x 0.069g/cm-2
    //Downstream Window: 0.011" x 0.075g/cm-2
    const double Window_ArialDensity =  ( 0.069+0.075 );//g/cm-2, Al-2219 Alloy
    double Window_Length = 15.;//cm

    double Target_Density = 0.0;//g*cm-3
    double Target_Thick = 0.0;//cm
    evt->SetVertexRange(-Target_Thick/2.0, Target_Thick/2.0);

    //H2
    Beam_Current = 20*1e-6;//25uA
    Target_Thick = Window_Length;
    Target_Density = 0.124/25.00;//gcm-3
    double Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
    double D2_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/2;//cm-2*s-1
    xem->Process(E0, P0, T0, 2, 1);
    double sigma0_d2 = xem->XS_Born() / MeV2GeV;//form nb/MeV/sr to nb/GeV/sr
    double sigma0_d2_rad = xem->XS_Rad() / MeV2GeV;//form nb/MeV/sr to nb/GeV/sr
    
    //H3
    Beam_Current = 20*1e-6;//25uA
    Target_Thick = Window_Length;
    Target_Density = 0.081/25.00;//g*cm-3
    double H3_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/3;//cm-2*s-1
    xem->Process(E0, P0, T0, 3, 1);
    double sigma0_h3 = xem->XS_Born() / MeV2GeV;//form nb/MeV/sr to nb/GeV/sr
    double sigma0_h3_rad = xem->XS_Rad() / MeV2GeV;//form nb/MeV/sr to nb/GeV/sr
    
    //He3
    Beam_Current = 25*1e-6;//25uA
    Target_Thick = Window_Length;
    Target_Density = 0.093/25.00;//g*cm-3
    double He3_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/3;//cm-2*s-1
    xem->Process(E0, P0, T0, 3, 2);
    double sigma0_he3 = xem->XS_Born() / MeV2GeV;//form nb/MeV/sr to nb/GeV/sr
    double sigma0_he3_rad = xem->XS_Rad() / MeV2GeV;//form nb/MeV/sr to nb/GeV/sr
    
    cout<<Form("H2: XS_RAD = %e,  XS_Born=%e", sigma0_d2_rad, sigma0_d2)<<endl;
    cout<<Form("RadCorr: R(H2) = %f, R(H3) = %f, R(He3) = %f", sigma0_d2_rad/sigma0_d2,sigma0_h3_rad/sigma0_h3,sigma0_he3_rad/sigma0_he3)<<endl;

    xs_d2[0] = sigma0_d2;
    xs_h3[0] = sigma0_h3;
    xs_he3[0] = sigma0_he3;
    /*}}}*/

    const double CosA = cos( T0*Deg2Rad );
    const double SinA = sin( T0*Deg2Rad );
    double x_tg, y_tg, z_tg, dp_tg,theta_tg, phi_tg,Theta, Ep, sigma_d2, sigma_h3, sigma_he3,sigma_al;
    double xb,nu, Q2, W;   
    double D2_XS_Sum = 0.0,H3_XS_Sum = 0.0,He3_XS_Sum = 0.0;
    double XS_Sum_Window = 0.0;

    /*Define ROOT file and tree{{{*/
    TString filename = Form("./results/rate_src_E%d_P%d_T%d.root", (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
    TFile *file = new TFile(filename,"recreate");
    TTree* T = new TTree("T","A New Tree");

    T->Branch("x_tg",&x_tg, "x_tg/D");
    T->Branch("y_tg",&x_tg, "y_tg/D");
    T->Branch("z_tg",&z_tg, "z_tg/D");
    T->Branch("dp_tg",&dp_tg, "dp_tg/D");
    T->Branch("theta_tg",&theta_tg, "theta_tg/D");
    T->Branch("phi_tg",&phi_tg, "phi_tg/D");

    T->Branch("E0",&E0, "E0/D");
    T->Branch("Ep",&Ep, "Ep/D");
    T->Branch("Theta",&Theta, "Theta/D");
    T->Branch("Q2",&Q2, "Q2/D");
    T->Branch("W",&W, "W/D");
    T->Branch("xb",&xb, "xb/D");
    T->Branch("nu",&nu, "nu/D");

    T->Branch("He3_Target_Lumi",&He3_Target_Lumi, "He3_Target_Lumi/D");
    T->Branch("H3_Target_Lumi",&H3_Target_Lumi, "H3_Target_Lumi/D");
    T->Branch("D2_Target_Lumi",&D2_Target_Lumi, "D2_Target_Lumi/D");
    T->Branch("Window_Lumi",&Window_Lumi, "Window_Lumi/D");
    T->Branch("D2_XS_Sum",&D2_XS_Sum, "D2_XS_Sum/D");
    T->Branch("H3_XS_Sum",&H3_XS_Sum, "H3_XS_Sum/D");
    T->Branch("He3_XS_Sum",&He3_XS_Sum, "He3_XS_Sum/D");
    T->Branch("XS_Sum_Window",&XS_Sum_Window, "XS_Sum_Window/D");

    T->Branch("sigma0_d2",&sigma0_d2, "sigma0_d2/D");
    T->Branch("sigma_d2",&sigma_d2, "sigma_d2/D");
    T->Branch("sigma0_h3",&sigma0_h3, "sigma0_h3/D");
    T->Branch("sigma_he3",&sigma_he3, "sigma_he3/D");
    T->Branch("Phase_Space",&Phase_Space, "Phase_Space/D");

    /*}}}*/

    /*Start sampling events{{{*/
    int Count=0;  double xb_sum=0.0; double W2_sum=0.0;
    int err = -1;
    //randomly generate events within the spectrometer phase space and caluclate the XS
    for(int i=0; i<Nevnt; i++){
        x_tg = evt->GetTG_X(); //cm
        y_tg = evt->GetTG_Y(); //cm
        z_tg = evt->GetTG_Z(); //cm

        dp_tg = evt->GetTG_Dp();//GeV
        Ep = P0 * (1.0 + dp_tg);//GeV

        theta_tg = evt->GetTG_Theta(); //rad
        phi_tg = evt->GetTG_Phi(); //rad
        Theta = acos( (CosA-phi_tg*SinA) / sqrt(1.+theta_tg*theta_tg+phi_tg*phi_tg) ) / Deg2Rad;//degree

        nu = E0-Ep;
        Q2 = 4.0*E0*Ep*pow(sin(Theta*Deg2Rad/2.0),2);
        xb = Q2 / (2.0*PROTON_MASS_GeV*nu);
        W = sqrt( PROTON_MASS_GeV*PROTON_MASS_GeV + 2.0*PROTON_MASS_GeV*nu - Q2 );

        if(IsDebug)
            printf("^^^^ Ep=%8.5fGeV, Theta=%8.5fDeg, xb = %4.3f\n",Ep,Theta, xb);

        sigma_d2 = 1e-33;  sigma_h3 = 1e-33; sigma_he3 = 1e-33; sigma_al = 1e-33;
        err = xem->Process(E0, Ep, Theta, 2, 1);
        if(err>=0) sigma_d2 = xem->XS_Born() / MeV2GeV;//form nb/MeV/sr to nb/GeV/sr
        err = xem->Process(E0, Ep, Theta, 3, 1);
        if(err>=0) sigma_h3 = xem->XS_Born() / MeV2GeV;//form nb/MeV/sr to nb/GeV/sr
        err = xem->Process(E0, Ep, Theta, 3, 2);
        if(err>=0) sigma_he3 = xem->XS_Born() / MeV2GeV;//form nb/MeV/sr to nb/GeV/sr
        err = xem->Process(E0, Ep, Theta, 27, 13);
        if(err>=0) sigma_al = xem->XS_Born() / MeV2GeV;//form nb/MeV/sr to nb/GeV/sr
 
        /*Calculate the cross sections{{{*/
        if(IsDebug)
            printf("^^^^ XS_D = %10.4e\n", sigma_d2);

        D2_XS_Sum += sigma_d2 * NBARN_TO_CM2;//cm2/GeV/sr
        H3_XS_Sum += sigma_h3 * NBARN_TO_CM2;//cm2/GeV/sr
        He3_XS_Sum += sigma_he3 * NBARN_TO_CM2;//cm2/GeV/sr

        XS_Sum_Window += sigma_al * NBARN_TO_CM2;//Assuming it is pure Al, Z = 13, N = 14
        xb_sum+=xb * sigma_h3 * NBARN_TO_CM2;;
        W2_sum+=W*W * sigma_h3 * NBARN_TO_CM2;;
        Count++;
        /*}}}*/

        if(IsDebug)
            printf("XS_Sum=%10.4e, XS_Sum_Window=%10.4e\n",H3_XS_Sum, XS_Sum_Window);
        T->Fill();
    }
    /*}}}*/

    /*Calculate rates{{{*/
    double He3_Target_Rate = He3_Target_Lumi * He3_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
    double H3_Target_Rate = H3_Target_Lumi * H3_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
    double D2_Target_Rate = D2_Target_Lumi * D2_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
    double Window_Rate = Window_Lumi * XS_Sum_Window * Phase_Space / Nevnt;  //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz

    xb_avg[0] = xb_sum/H3_XS_Sum;
    W2_avg[0] = W2_sum/H3_XS_Sum;
    xs_d2_avg[0] = D2_XS_Sum/Count / NBARN_TO_CM2;
    xs_h3_avg[0] = H3_XS_Sum/Count / NBARN_TO_CM2;
    xs_he3_avg[0] = He3_XS_Sum/Count / NBARN_TO_CM2;
    /*}}}*/

    rate_d2[0] = D2_Target_Rate;
    rate_h3[0] = H3_Target_Rate;
    rate_he3[0] = He3_Target_Rate;
    rate_window[0] = Window_Rate;

    file->cd(); T->Write(); file->Close();
    delete evt;
    delete xem;
}
/*}}}GetRate*/

/*main{{{*//*{{{*//*}}}*/
int main(){

    const int NBin = 7;
    //const double E0 = 4.4;
    const double Evt_Bin = 25000;
    const double EBeam[NBin] = {4.4, 4.4, 4.4, 4.4, 2.2, 2.2,2.2};
    const double P0[NBin] = {3.710, 3.980, 3.535, 3.790, 1.870, 1.770, 1.665};
    const double T0[NBin] = {17.0, 17.0, 19.0, 19.0, 22.0, 26.2, 30.3};

    double Q2, W2, xbj;
    double rate_he3, rate_h3, rate_d2, rate_w, hour_he3, hour_h3, hour_d2;
    double total_hour,xs_he3, xs_h3, xs_d2,xs_he3_avg, xs_h3_avg, xs_d2_avg, xbj_avg,W2_avg;
    double total_day=0.0;

    FILE *new_table; 
    new_table =fopen(Form("hrs_table_E3.6_lowW2.dat"),"w");

    //only with HRS
    const string spec = "HRS";
    for(int i=0;i<NBin; i++){/*{{{*/
        double E0 = EBeam[i];
        double Ep = P0[i];
        double Theta = T0[i];
        double nu = E0-Ep;
        Q2 = 4.0*E0*Ep*pow(sin(Theta*Deg2Rad/2.0),2);
        xbj = Q2 / (2.0*PROTON_MASS_GeV*nu);
        W2 = PROTON_MASS_GeV*PROTON_MASS_GeV + 2.0*PROTON_MASS_GeV*nu - Q2;

        GetRate(E0, Ep, Theta,&xbj_avg,&W2_avg, &xs_he3, &xs_he3_avg, &rate_he3, &xs_h3,&xs_h3_avg, &rate_h3, &xs_d2,&xs_d2_avg, &rate_d2, &rate_w);

        //Acount for 10% loss of using SOS quards/*{{{*/
        if(spec=="HRS"){
            rate_he3 *=0.9;
            rate_h3 *=0.9;
            rate_d2 *=0.9;
            rate_w *=0.9;
            //80% due to RC, 80% due to efficiecies
            rate_he3 *=0.8*0.8;
            rate_h3 *=0.8*0.8;
            rate_d2 *=0.8*0.8;
            rate_w *=0.8*0.8;
        }
        if(rate_he3>1e-4)
            hour_he3 =  Evt_Bin/(rate_he3*60.*60.);
        else
            hour_he3 = 0;
        if(rate_h3>1e-4)
            hour_h3 =  Evt_Bin/(rate_h3*60.*60.);
        else
            hour_h3 = 0;
        if(rate_d2>1e-4)
            hour_d2 =  Evt_Bin/(rate_d2*60.*60.);
        else
            hour_d2 = 0;
        total_hour = hour_d2+hour_h3+hour_he3;
        //if(hour_d2<1.0&&hour_d2>0.00001)  hour_d2 = 1.0;
        //if(hour_h3<1.0&&hour_h3>0.00001)  hour_h3 = 1.0;
        //if(hour_he3<1.0&&hour_h3>0.00001) hour_he3 = 1.0;
        /*}}}*/
        
        printf("--- P0=%4.3f GeV/c, T0=%3.1f Deg, Xbj = %3.2f, Q2 = %3.2f GeV^2, Rates = %6.2f(H2), %6.2f(H3), %6.2f(He3), %6.2f(Win) Hz\n", P0[i], T0[i],xbj, Q2, rate_d2, rate_h3, rate_he3, rate_w);

        fprintf(new_table,"%4s & %4.2f(%4.2f) & %4.2f(%4.2f) & %4.2f & %4.2f  & %4.2f & %6.4f &  %6.4f & %6.4f & %6.3f(%4.2f) & %6.3f(%4.2f) & %6.3f(%4.2f) & %5.2f  & %6.3f \\\\\\hline \n",
                "HRS", xbj,xbj_avg, W2,W2_avg, Q2, P0[i], T0[i],xs_d2_avg, xs_h3_avg, xs_he3_avg, rate_d2,hour_d2, rate_h3,hour_h3,  rate_he3,hour_he3, total_hour, rate_w);
        total_day +=total_hour;
    }
    /*}}}*/
    total_day /= 24.0;
    fprintf(new_table, " &  &  &  &  &  &  &  &  &  &  &  &  & %5.2f\\\\\\hline \n", total_day);
    fclose(new_table);
}
/*}}}*/
