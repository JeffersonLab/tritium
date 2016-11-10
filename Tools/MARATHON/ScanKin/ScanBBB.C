#include "../GetRate.h"
#include "DIS_Lite.h"
/*GetRate{{{*/
void GetRate(const string spec,  const string Com, double E0, double P0, double T0, 
		double* xb_avg, double* xs_he3,double* xs_he3_avg,double* rate_he3, double* xs_h3,double* xs_h3_avg,double* rate_h3, double* xs_h2,double* xs_h2_avg,double* rate_h2, double *rate_window){
	bool IsDebug = false;
	//bool IsDebug = true;
	/*
	   double E0 = 11.0; //GeV
	   double P0; //GeV, Scattered Electrons' Central Momentum
	   double T0; //Degree, Scattered Electrons' Central Angle
	   string spec = "";
	   cout<<" --- Which Spectrometer (HRS or BB) = "; cin >> spec;
	   cout<<" --- Beam Energy (E0, GeV) = "; cin >> E0;
	   cout<<" --- Scattered Electron's Central Momentum (P0, GeV) = "; cin >> P0;
	   cout<<" --- Scattered Electron's Central Angle (Theta0, Deg) = "; cin >> T0;
	   */
	int Nevnt = 20000000;
	if(IsDebug){cout<<" ----Total Samples = "; cin >> Nevnt;}
	cout<<endl<<endl;

    DIS* dis = new DIS();
    dis->SetKin(E0, P0, T0*Deg2Rad);
 
	Event* evt = new Event();     
	evt->SetSpectrometer(spec);
	if(spec=="BB"){
		if(Com=="narrow")
			evt->SetTGDpRange(-0.045, 0.045);
		else{
			double pmin = (0.6  - P0)/P0;//BB's minimum momentum set to 1GeV
			//double pmin = (1.0  - P0)/P0;//BB's minimum momentum set to 1GeV
			double pmax = (11.0 - P0)/P0;//BB's maximum momentum set to 11GeV
			evt->SetTGDpRange(pmin, pmax);
		}	
	}	
	double Phase_Space = evt->GetPhaseSpace() * P0; //sr*GeV //it was deltaTheta * deltaPhi *deltaDp

	/*Define Target Info{{{*/
	double Beam_Current = 0.0;//A

	//Upstream Window:   0.010" x 0.069g/cm-2
	//Downstream Window: 0.011" x 0.075g/cm-2
	const double Window_ArialDensity =  ( 0.069+0.075 );//g/cm-2, Al-2219 Alloy
	double Window_Length = 25.;//cm
	if(spec=="HRS")
		Window_Length = 15.0;//cm

	double Target_Density = 0.0;//g*cm-3
	double Target_Thick = 0.0;//cm
	evt->SetVertexRange(-Target_Thick/2.0, Target_Thick/2.0);

	Beam_Current = 20*1e-6;//25uA
	Target_Thick = Window_Length;
	Target_Density = 0.124/25.00;//gcm-3

	double Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
	double D2_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/2;//cm-2*s-1
    double sigma0_h2 = dis->Sigma(2,1);

	Beam_Current = 20*1e-6;//25uA
	Target_Thick = Window_Length;
	Target_Density = 0.081/25.00;//g*cm-3

	double H3_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/3;//cm-2*s-1
    double sigma0_h3 = dis->Sigma(3,1);


	Beam_Current = 25*1e-6;//25uA
	Target_Thick = Window_Length;
	Target_Density = 0.093/25.00;//g*cm-3

	double He3_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/3;//cm-2*s-1
    double sigma0_he3 = dis->Sigma(3,2);

	xs_h2[0] = sigma0_h2;
	xs_h3[0] = sigma0_h3;
	xs_he3[0] = sigma0_he3;
	/*}}}*/

	const double CosA = cos( T0*Deg2Rad );
	const double SinA = sin( T0*Deg2Rad );
	double dp_tg,theta_tg, phi_tg,Theta, Ep, sigma_h2, sigma_h3, sigma_he3,sigma_al;
	double xs_p,xs_n,xs_d;
	double xb,nu, Q2, W;   
	double D2_XS_Sum = 0.0,H3_XS_Sum = 0.0,He3_XS_Sum = 0.0;
	double XS_Sum_Window = 0.0;

	/*Define ROOT file and tree{{{*/
    TString filename = Form("./results3/%s_rate_marathon_E%d_P%d_T%d_lowE.root", spec.c_str(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
    if(IsDebug)
        filename = Form("./results3/%s_rate_marathon_E%d_P%d_T%d_debug_lowE.root", spec.c_str(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
    if(spec=="BB" &&Com=="narrow")	
        filename = Form("./results3/%s_rate_marathon_E%d_P%d_T%d_narrow_lowE.root",spec.c_str(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
    TFile *file = new TFile(filename,"recreate");
    TTree* T = new TTree("T","A New Tree");

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
    T->Branch("xs_p",&xs_p, "xs_p/D");
    T->Branch("xs_n",&xs_n, "xs_n/D");
    T->Branch("xs_d",&xs_d, "xs_d/D");

    T->Branch("sigma0_h2",&sigma0_h2, "sigma0_h2/D");
    T->Branch("sigma_h2",&sigma_h2, "sigma_h2/D");
    T->Branch("sigma0_h3",&sigma0_h3, "sigma0_h3/D");
    T->Branch("sigma_he3",&sigma_he3, "sigma_he3/D");
    T->Branch("Phase_Space",&Phase_Space, "Phase_Space/D");

	/*}}}*/

	/*Start sampling events{{{*/
    int Count=0; double xb_sum=0.0;
    //randomly generate events within the spectrometer phase space and caluclate the XS
	for(int i=0; i<Nevnt; i++){
		xs_p = -1e-33; xs_n = -1e-33; xs_d = -1e-33; 

		dp_tg = evt->GetTG_Dp();//GeV
		Ep = P0 * (1.0 + dp_tg);//GeV

		theta_tg = evt->GetTG_Theta(); //rad
		phi_tg = evt->GetTG_Phi(); //rad
		Theta = acos( (CosA-phi_tg*SinA) / sqrt(1.+theta_tg*theta_tg+phi_tg*phi_tg) );//rad

		nu = E0-Ep;
		Q2 = 4.0*E0*Ep*pow(sin(Theta/2.0),2);
		xb = Q2 / (2.0*PROTON_MASS*nu);
		W = sqrt( PROTON_MASS*PROTON_MASS + 2.0*PROTON_MASS*nu - Q2 );

		if(IsDebug)
			printf("^^^^ Ep=%8.5fGeV, Theta=%8.5fDeg, xb = %4.3f\n"
					,Ep,Theta/Deg2Rad, xb);

		sigma_h2 = 0;
		sigma_h3 = 0;
		sigma_he3 = 0;
		sigma_al = 0;
        /*Calculate the cross sections{{{*/
        if(xb<1.0){
            dis->SetKin(E0, Ep, Theta);
            xs_p = dis->Sigma("Proton");
            xs_n = dis->Sigma("Neutron");
            xs_d = dis->Sigma("Deutron");

            if(isnan(xs_p)||xs_p<-1e-20)
                xs_p=1e-66;
            if(isnan(xs_n)||xs_n<-1e-20)
                xs_n=1e-66;
            if(isnan(xs_d)||xs_d<-1e-20)
                xs_d=1e-66;

            if(IsDebug)
                printf("^^^^ XS_P = %10.4e, XS_N =  %10.4e\n", xs_p, xs_n);

            sigma_h2 = xs_d;
            sigma_h3 = xs_d+xs_n;
            sigma_he3 = xs_d+xs_p;
            sigma_al = (xs_p * 13 + xs_n * 14);

            if(IsDebug)
                printf("^^^^ XS_D = %10.4e\n", xs_d);

            //if(W>2.0&&xb<1.00){//Only consider the DIS events with W>2.0 cut, for high-x poins, we may need smaller W cut.
            if((W>2.0&&xb<0.70)||(W>sqrt(3.0)&&xb>0.70)){//Only consider the DIS events with W>2.0 cut, for high-x poins, we may need smaller W cut.
                D2_XS_Sum += sigma_h2 * NBARN_TO_CM2;//cm2/GeV/sr
                H3_XS_Sum += sigma_h3 * NBARN_TO_CM2;//cm2/GeV/sr
                He3_XS_Sum += sigma_he3 * NBARN_TO_CM2;//cm2/GeV/sr
                XS_Sum_Window += sigma_al * NBARN_TO_CM2;//Assuming it is pure Al, Z = 13, N = 14
                xb_sum += xb;
                Count++;
            }
        }
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

    xb_avg[0] = xb_sum/Count;
	xs_h2_avg[0] = D2_XS_Sum/Count / NBARN_TO_CM2;
	xs_h3_avg[0] = H3_XS_Sum/Count / NBARN_TO_CM2;
	xs_he3_avg[0] = He3_XS_Sum/Count / NBARN_TO_CM2;
	/*}}}*/

	rate_h2[0] = D2_Target_Rate;
	rate_h3[0] = H3_Target_Rate;
	rate_he3[0] = He3_Target_Rate;
	rate_window[0] = Window_Rate;

    file->cd(); T->Write(); file->Close();
	delete evt;
    delete dis;
}
/*}}}GetRate*/

/*main{{{*/
int main(){

	//const int N_X = 19;
	const int N_X =1;
	const double X_Bin[N_X] = { 0.83};


	//const int N_T = 4;
	//const int Theta_Bin[N_T] = {20, 21, 22, 23};
	//const int Theta_Bin[N_T] = {24, 25, 26, 27};
	//const int Theta_Bin[N_T] = {28, 29, 30, 31};
	//const int Theta_Bin[N_T] = {32, 33, 34, 35};
	//const int Theta_Bin[N_T] = {42, 47, 52, 57};
	//const int N_T = 28;
	//const int Theta_Bin[N_T] = {30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,46,48,50,52,54,55,56,57,58,59,60,61,62};
	//const int N_T =12;
	//const int Theta_Bin[N_T] = {51, 52, 53, 54, 55, 56, 57, 58,59,60,61,62};

    //const int N_T =4;
    //const int Theta_Bin[N_T]={63,65,67,69};
    const int N_T =10;
    const int Theta_Bin[N_T]={39,41,43,45, 47, 49, 51,53,55,69};


	const double E0 = 11.0;
	double P0, T0, xbj, W2, Q2, nu;
	double xs_he3, xs_h3, rate_he3, rate_h3, xs_he3_avg, xs_h3_avg;
	double xs_h2, rate_h2, xs_h2_avg, rate_w, xbj_avg;
    double rate_he3_bb, rate_h3_bb, rate_h2_bb,rate_w_bb;
    double hour_he3_bb, hour_h3_bb, hour_h2_bb;
	const string spec = "BB";
	double hour_h2,hour_h3 = 0,hour_he3 = 0;
    double total_hour,total_hour_bb;

	//cerr<<"--- Theta0 (degree) = "; cin >> T0;
	for(int j=0;j<N_T;j++){
		T0 = Theta_Bin[j];
		ofstream new_table(Form("hx_%s_table_T%d_lowE.dat",spec.c_str(), (int)T0));
		ofstream new_table1(Form("hx_%s_table_T%d_large_lowE.dat",spec.c_str(), (int)T0));
	
        double SinSQ = pow(sin(T0*Deg2Rad/2.), 2);
		for(int i=0;i<N_X;i++){
			xbj = X_Bin[i];
			P0 = 2.0*PROTON_MASS*E0*xbj/(2.0*PROTON_MASS*xbj + 4.0*E0*SinSQ);
			nu = E0 - P0;
			Q2 = 4.0 * E0 * P0 *SinSQ;
			W2 = PROTON_MASS*PROTON_MASS + 2.0*PROTON_MASS*nu - Q2;

            GetRate(spec.c_str(),"", E0, P0, T0, &xbj_avg, &xs_he3, &xs_he3_avg, &rate_he3, &xs_h3,&xs_h3_avg, &rate_h3, &xs_h2,&xs_h2_avg, &rate_h2, &rate_w);
            if(spec.find("BB",0)!=string::npos)
                GetRate(spec.c_str(),"narrow", E0, P0, T0, &xbj_avg, &xs_he3, &xs_he3_avg, &rate_he3_bb, &xs_h3,&xs_h3_avg, &rate_h3_bb, &xs_h2,&xs_h2_avg, &rate_h2_bb, &rate_w_bb);

			//Consider 90% reduction in acceptance when using SOS quards	

            if(spec.find("HRS",0)!=string::npos){
                rate_h2 *=0.9;			
                rate_h3 *=0.9;			
                rate_he3 *=0.9;			
            }

            //20% RC and 20% Efficiencies
            rate_h2 *=0.8*0.8;			
            rate_h3 *=0.8*0.8;			
            rate_he3 *=0.8*0.8;			

			if(rate_h2>1e-4)
				hour_h2 =  25000./(rate_h2*60.*60.);
			else
				hour_h2 = 0;
			if(rate_he3>1e-4)
				hour_he3 =  25000./(rate_he3*60.*60.);
			else
				hour_he3 = 0;
			if(rate_h3>1e-4)
				hour_h3 =  25000./(rate_h3*60.*60.);
			else
				hour_h3 = 0;
            
            if(rate_he3_bb>1e-4)
                hour_he3_bb =  25000./(rate_he3_bb*60.*60.);
            else
                hour_he3_bb = 0;
            if(rate_h3_bb>1e-4)
                hour_h3_bb =  25000./(rate_h3_bb*60.*60.);
            else
                hour_h3_bb = 0;
            if(rate_h2_bb>1e-4)
                hour_h2_bb =  25000./(rate_h2_bb*60.*60.);
            else
                hour_h2_bb = 0;

			//if(hour_h2<1.0&&hour_h2>0.00001)  hour_h2 = 1.0;
			//if(hour_h3<1.0&&hour_h3>0.00001)  hour_h3 = 1.0;
			//if(hour_he3<1.0&&hour_h3>0.00001) hour_he3 = 1.0;

            total_hour = hour_h3+hour_he3+hour_h2;
            total_hour_bb = hour_h3_bb+hour_he3_bb+hour_h2_bb;

            cout<<Form("&&&& xbj = %4.2f, Ep = %4.2f,  Theta = %4.2f", xbj, P0, T0)<<endl; 
            cout<<Form("                  H2 XS_New = %8.4f,  H2 XS_Avg = %8.4f", xs_h2, xs_h2_avg)<<endl; 
            cout<<Form("                  H3 XS_New = %8.4f,  H3 XS_Avg = %8.4f", xs_h3, xs_h3_avg)<<endl; 
            cout<<Form("                 He3 XS_New = %8.4f, He3 XS_Avg = %8.4f", xs_he3, xs_he3_avg)<<endl; 


            if(spec.find("BB",0)!=string::npos){
                new_table<<Form("%4s & %4.2f(%4.2f) & %4.2f & %4.2f & %4.2f  & %4.2f & %6.4f &  %6.4f & %6.4f & %6.3f(%4.2f) & %6.3f(%4.2f) & %6.3f(%4.2f) & %6.3f  & %5.2f \\\\\\hline ",
                        spec.c_str(),xbj,xbj_avg, W2, Q2, P0, T0, xs_h2_avg, xs_h3_avg, xs_he3_avg, rate_h2_bb,hour_h2_bb, rate_h3_bb,hour_h3_bb,rate_he3_bb,hour_he3_bb,rate_w, total_hour_bb)<<endl;
                new_table1<<Form(" %4s & %4.2f(%4.2f) & %4.2f & %4.2f & %4.2f  & %4.2f & %6.4f &  %6.4f & %6.4f & %6.3f(%4.2f) & %6.3f(%4.2f) & %6.3f(%4.2f) & %6.3f   & %5.2f \\\\\\hline ",
                        spec.c_str(),xbj,xbj_avg, W2, Q2, P0, T0, xs_h2_avg, xs_h3_avg, xs_he3_avg, rate_h2,hour_h2, rate_h3,hour_h3,rate_he3,hour_he3,rate_w, total_hour)<<endl;
            }else
                new_table<<Form(" %4s & %4.2f(%4.2f) & %4.2f & %4.2f & %4.2f  & %4.2f & %6.4f &  %6.4f & %6.4f & %6.3f(%4.2f) & %6.3f(%4.2f) & %6.3f(%4.2f) & %6.3f   & %5.2f \\\\\\hline ",
                        spec.c_str(),xbj,xbj_avg, W2, Q2, P0, T0, xs_h2_avg, xs_h3_avg, xs_he3_avg, rate_h2,hour_h2, rate_h3,hour_h3,rate_he3,hour_he3,rate_w, total_hour)<<endl;
        }
        new_table.close(); 
        new_table1.close(); 
    }
}
/*}}}*/
