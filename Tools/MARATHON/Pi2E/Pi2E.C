#include "../GetRate.h"
#include "../DIS_Lite.h"
/*GetRate{{{*/
void GetRate(const string spec,  const string Com, double E0, double P0, double T0, 
		double* xb_avg, double* xs_he3,double* xs_he3_avg,double* rate_he3, double* xs_h3,double* xs_h3_avg,double* rate_h3, double* xs_d2,double* xs_d2_avg,double* rate_d2, double *rate_window){
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
	int Nevnt = 100000;
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
			double pmin = (1.0  - P0)/P0;//BB's minimum momentum set to 1GeV
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
	//if(spec=="HRS")
		//Window_Length = 15.0;//cm

	double Target_Density = 0.0;//g*cm-3
	double Target_Thick = 0.0;//cm
	evt->SetVertexRange(-Target_Thick/2.0, Target_Thick/2.0);

	Beam_Current = 20*1e-6;//25uA
	Target_Thick = Window_Length;
	Target_Density = 0.124/25.00;//gcm-3

	double Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
	double D2_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/2;//cm-2*s-1
    double sigma0_d2 = dis->Sigma(2,1);

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

	xs_d2[0] = sigma0_d2;
	xs_h3[0] = sigma0_h3;
	xs_he3[0] = sigma0_he3;
	/*}}}*/

	const double CosA = cos( T0*Deg2Rad );
	const double SinA = sin( T0*Deg2Rad );
	double x_tg, y_tg, z_tg, dp_tg,theta_tg, phi_tg,Theta, Ep, sigma_d2, sigma_h3, sigma_he3,sigma_al;
	double xs_p,xs_n,xs_d;
	double xb,nu, Q2, W;   
	double D2_XS_Sum = 0.0,H3_XS_Sum = 0.0,He3_XS_Sum = 0.0;
	double XS_Sum_Window = 0.0;

	/*Start sampling events{{{*/
    int Count=0; double xb_sum=0.0;
    //randomly generate events within the spectrometer phase space and caluclate the XS
	for(int i=0; i<Nevnt; i++){
		xs_p = -1e-33; xs_n = -1e-33; xs_d = -1e-33; 

		x_tg = evt->GetTG_X(); //cm
		y_tg = evt->GetTG_Y(); //cm
		z_tg = evt->GetTG_Z(); //cm

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

		sigma_d2 = 0;
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

            sigma_d2 = xs_d;
            sigma_h3 = xs_d+xs_n;
            sigma_he3 = xs_d+xs_p;
            sigma_al = (xs_p * 13 + xs_n * 14);

            if(IsDebug)
                printf("^^^^ XS_D = %10.4e\n", xs_d);

            //if(W>2.0&&xb<1.00){//Only consider the DIS events with W>2.0 cut, for high-x poins, we may need smaller W cut.
            if((W>2.0&&xb<0.70)||(W>sqrt(3.0)&&xb>0.70)){//Only consider the DIS events with W>2.0 cut, for high-x poins, we may need smaller W cut.
                D2_XS_Sum += sigma_d2 * NBARN_TO_CM2;//cm2/GeV/sr
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
	}
	/*}}}*/

	/*Calculate rates{{{*/
	double He3_Target_Rate = He3_Target_Lumi * He3_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
	double H3_Target_Rate = H3_Target_Lumi * H3_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
	double D2_Target_Rate = D2_Target_Lumi * D2_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
	double Window_Rate = Window_Lumi * XS_Sum_Window * Phase_Space / Nevnt;  //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz

    xb_avg[0] = xb_sum/Count;
	xs_d2_avg[0] = D2_XS_Sum/Count / NBARN_TO_CM2;
	xs_h3_avg[0] = H3_XS_Sum/Count / NBARN_TO_CM2;
	xs_he3_avg[0] = He3_XS_Sum/Count / NBARN_TO_CM2;
	/*}}}*/

	rate_d2[0] = D2_Target_Rate;
	rate_h3[0] = H3_Target_Rate;
	rate_he3[0] = He3_Target_Rate;
	rate_window[0] = Window_Rate;

	delete evt;
    delete dis;
}
/*}}}GetRate*/

/*main{{{*/
int main(){
    ifstream table("marathon_pitoe_part2.dat");
    string com, spec;
	table >> com >> com >> com >> com >> com >> com >> com ;

    const double E0=11.0;//GeV
	double P0, T0, xbj, xb, Q2, W2, xs_e, xs_pim, xs_pi0, r_pi2e;

	double xs_he3, xs_h3, xs_d2, rate_he3, rate_h3, rate_d2;
	double rate_he3_bb, rate_h3_bb, rate_d2_bb, rate_w, rate_w_bb;
    double xs_h3_avg, xs_he3_avg,xs_d2_avg, xb_avg;
    double rate_pim_d2, rate_pim_h3, rate_pim_he3;

	ofstream new_table(Form("pi2e_T47_latex.dat"));
	ofstream new_table1(Form("pi2e_T47.dat"));
    new_table1<<Form("%6s   %6s   %6s   %6s   %6s   %6s   %6s  %6s   %6s   %6s   %6s   %6s",
            spec.c_str(),"xbj","xb_avg", "W2", "Q2", "P0", "Theta", "xs_e", "xs_avg","xs_pim","rate_h3", "rate_pim")<<endl;
     

    TFile *file = new TFile("marathon_pi2e.root","recreate");
    TTree *T = new TTree("T","a new tree");
    T->Branch("P0", &P0 , "P0/D");
    T->Branch("T0", &T0 , "T0/D");
    T->Branch("Q2", &Q2 , "Q2/D");
    T->Branch("W2", &W2 , "W2/D");
    T->Branch("xbj", &xbj , "xbj/D");
    T->Branch("xs_e", &xs_e , "xs_e/D");
    T->Branch("xs_pim", &xs_pim , "xs_pim/D");
    T->Branch("xs_pi0", &xs_pi0 , "xs_pi0/D");
    T->Branch("r_pi2e", &r_pi2e , "r_pi2e/D");

    T->Branch("xs_h3", &xs_h3, "xs_h3/D");
    T->Branch("xs_he3", &xs_he3, "xs_he3/D");
    T->Branch("xs_d2", &xs_d2, "xs_d2/D");

    while(table >> xb >> T0 >> P0 >>xs_e >> xs_pim >> xs_pi0 >> r_pi2e){
	
		double SinSQ = pow(sin(T0*Deg2Rad/2.), 2);
		double nu = E0 - P0;
		Q2 = 4.0 * E0 * P0 *SinSQ;
		W2 = PROTON_MASS*PROTON_MASS + 2.0*PROTON_MASS*nu - Q2;
        xbj = Q2/2.0/nu/PROTON_MASS;
	
		spec ="HRS";

        GetRate(spec.c_str(),"", E0, P0, T0, &xb_avg, &xs_he3, &xs_he3_avg, &rate_he3, &xs_h3,&xs_h3_avg, &rate_h3, &xs_d2,&xs_d2_avg, &rate_d2, &rate_w);

        rate_he3_bb = 0.0;
        rate_h3_bb = 0.0;
        rate_d2_bb = 0.0;
		if(spec=="BB"){
			GetRate(spec.c_str(),"narrow", E0, P0, T0, &xb_avg, &xs_he3,&xs_he3_avg, &rate_he3_bb, &xs_h3,&xs_h3_avg, &rate_h3_bb, &xs_d2,&xs_d2_avg, &rate_d2_bb, &rate_w_bb);
		}

        //Acount for 10% loss of using SOS quards
        if(spec=="HRS"){
            rate_he3 *=0.9;
            rate_h3 *=0.9;
            rate_d2 *=0.9;
        }
        //80% due to RC, 80% due to efficiecies
        rate_he3 *=0.8*0.8;
        rate_h3 *=0.8*0.8;
        rate_d2 *=0.8*0.8;

        rate_pim_d2  = xs_pim/xs_d2_avg * rate_d2;
        rate_pim_h3  = xs_pim/xs_h3_avg * rate_h3;
        rate_pim_he3 = xs_pim/xs_he3_avg* rate_he3;
        //rate_pim_d2 = r_pi2e * rate_d2;
        //rate_pim_h3 = r_pi2e * rate_h3;
        //rate_pim_he3 = r_pi2e * rate_he3;

		cout<<Form("&&&& xbj = %6.2f (%6.2f),  H3 XS_Old = %8.4f,  H3 XS_New = %8.4f,  H3 XS_Avg = %8.4f",xbj,xb, xs_e,xs_h3, xs_h3_avg)<<endl; 

        if(spec=="BB") {
            rate_pim_d2  = xs_pim/xs_d2_avg * rate_d2_bb;
            rate_pim_h3  = xs_pim/xs_h3_avg * rate_h3_bb;
            rate_pim_he3 = xs_pim/xs_he3_avg* rate_he3_bb;
            //rate_pim_d2 = r_pi2e * rate_d2_bb;
            //rate_pim_h3 = r_pi2e * rate_h3_bb;
            //rate_pim_he3 = r_pi2e * rate_he3_bb;

            new_table<<Form("%6s & %6.2f(%6.2f) & %6.2f & %6.2f & %6.2f  & %6.2f & %6.4f &  %6.4f & %6.4f(%6.4f) & %6.4f\\\\\\hline ",
                    spec.c_str(),xbj,xb_avg, W2, Q2, P0, T0, xs_e, xs_h3_avg, rate_h3,rate_h3_bb, rate_pim_h3)<<endl;
            new_table1<<Form("%6s   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f  %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f ",
                    spec.c_str(),xbj,xb_avg, W2, Q2, P0, T0, xs_e, xs_h3_avg,xs_pim,  rate_h3,  rate_h3_bb, rate_pim_h3)<<endl;
                 }
        else{
         new_table<<Form("%6s & %6.2f(%6.2f) & %6.2f & %6.2f & %6.2f  & %6.2f & %6.4f &  %6.4f & %6.4f & %6.4f\\\\\\hline ",
                    spec.c_str(),xbj,xb_avg, W2, Q2, P0, T0, xs_e, xs_h3_avg, rate_h3, rate_pim_h3)<<endl;
         new_table1<<Form("%6s   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f  %6.4f   %6.4f   %6.4f   %6.4f   %6.4f",
                    spec.c_str(),xbj,xb_avg, W2, Q2, P0, T0, xs_e, xs_h3_avg,xs_pim,  rate_h3, rate_pim_h3)<<endl;
         }

       T->Fill();
    }
    new_table.close();
    new_table1.close();
    //
    file->cd(); T->Write(); file->Close();
}
/*}}}*/
