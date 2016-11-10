#include "GetRate.h"
#include "christy_bosted_inelastic.h"
#include "DIS.h"
/*GetRateAll{{{*/
/*void GetRateAll(string spec, double E0, double P0, double T0, string Com, */
		//double* xs_he3,double* xs_he3_avg,double* rate_he3, 
		//double* xs_h3,double* xs_h3_avg,double* rate_h3, 
		//double* xs_d2,double* xs_d2_avg,double* rate_d2, double *rate_window){

    //int Nevnt = 100000;

	//Event* evt = new Event();     
	//evt->SetSpectrometer(spec);
	//if(spec=="BB"){
		//if(Com=="narrow")
			//evt->SetTGDpRange(-0.045, 0.045);
		//else{
			//double pmin = (1.0  - P0)/P0;//BB's minimum momentum set to 1GeV
			//double pmax = (11.0 - P0)/P0;//BB's maximum momentum set to 11GeV
			//evt->SetTGDpRange(pmin, pmax);
		//}	
	//}	
	//double Phase_Space = evt->GetPhaseSpace() * P0; //sr*GeV //it was deltaTheta * deltaPhi *deltaDp

	//[>Define Target Info{{{<]
	//double Beam_Current = 0.0;//A

	////Upstream Window:   0.010" x 0.069g/cm-2
	////Downstream Window: 0.011" x 0.075g/cm-2
	//const double Window_ArialDensity =  ( 0.069+0.075 );//g/cm-2, Al-2219 Alloy
	//double Window_Length = 25.;//cm
	//if(spec=="HRS")
		//Window_Length = 15.0;//cm

	//double Target_Density = 0.0;//g*cm-3
	//double Target_Thick = 0.0;//cm
	//evt->SetVertexRange(-Target_Thick/2.0, Target_Thick/2.0);

	//Beam_Current = 20*1e-6;//25uA
	//Target_Thick = Window_Length;
	//Target_Density = 0.124/25.00;//gcm-3
	//double Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
	//double D2_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/2;//cm-2*s-1
	//double sigma0_d2 = (sigma_p(E0, T0*Deg2Rad, P0)+sigma_n(E0, T0*Deg2Rad, P0)) ;

	//Beam_Current = 20*1e-6;//25uA
	//Target_Thick = Window_Length;
	//Target_Density = 0.081/25.00;//g*cm-3
	//double H3_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/3;//cm-2*s-1
	//double sigma0_h3 = (sigma_p( E0, T0*Deg2Rad, P0)+2.*sigma_n(E0, T0*Deg2Rad, P0)) ;

	//Beam_Current = 25*1e-6;//25uA
	//Target_Thick = Window_Length;
	//Target_Density = 0.093/25.00;//g*cm-3
	//double He3_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/3;//cm-2*s-1
	//double sigma0_he3 = (2.*sigma_p( E0, T0*Deg2Rad, P0)+sigma_n(E0, T0*Deg2Rad, P0)) ;

	//xs_d2[0] = sigma0_d2;
	//xs_h3[0] = sigma0_h3;
	//xs_he3[0] = sigma0_he3;
	//[>}}}<]

	//const double CosA = cos( T0*Deg2Rad );
	//const double SinA = sin( T0*Deg2Rad );
	//double dp_tg,theta_tg, phi_tg,Theta, Ep, sigma_d2, sigma_h3, sigma_he3, sigma_al;
	//double xs_p,xs_n;
	//double xb,nu, Q2, W;   
	//double D2_XS_Sum = 0.0,H3_XS_Sum = 0.0,He3_XS_Sum = 0.0;
	//double XS_Sum_Window = 0.0;
	//double dis_xs_p,dis_xs_n;

	//[>Start sampling events{{{<]
	//int Count=0;
	////randomly generate events within the spectrometer phase space and caluclate the XS
	//for(int i=0; i<Nevnt; i++){
		//xs_p = -1e-33; xs_n = -1e-33;

		//dp_tg = evt->GetTG_Dp();//GeV
		//Ep = P0 * (1.0 + dp_tg);//GeV

		//theta_tg = evt->GetTG_Theta(); //rad
		//phi_tg = evt->GetTG_Phi(); //rad
		//Theta = acos( (CosA-phi_tg*SinA) / sqrt(1.+theta_tg*theta_tg+phi_tg*phi_tg) );//rad

		//nu = E0-Ep;
		//Q2 = 4.0*E0*Ep*pow(sin(Theta/2.0),2);
		//xb = Q2 / (2.0*PROTON_MASS*nu);
		//W = sqrt( PROTON_MASS*PROTON_MASS + 2.0*PROTON_MASS*nu - Q2 );

		//[>Calculate the cross sections{{{<]
		//sigma_h3 = 0;
		//sigma_he3 = 0;
		//xs_p = sigma_p( E0, Theta, Ep); //nb/GeV/sr
		//xs_n = sigma_n( E0, Theta, Ep);//nb/GeV/sr
		//if(isnan(xs_p)||xs_p<-1e-20)
			//xs_p=1e-66;
		//if(isnan(xs_n)||xs_n<-1e-20)
			//xs_n=1e-66;

        //dis_xs_p = dis_sigma_p( E0, Ep, Theta); //nb/GeV/sr
		//dis_xs_n = dis_sigma_n( E0, Ep, Theta);//nb/GeV/sr

        //if(!(i%1000))
            //printf("xb = %3.2f, xs_p= %6.4e/%6.4e, xs_n= %6.4e/%6.4e\n", xb, xs_p, dis_xs_p, xs_n, dis_xs_n);

        //sigma_d2  = dis_xs_p+dis_xs_n;
		//sigma_h3  = dis_xs_p+2.*dis_xs_n;
		//sigma_he3 = 2.*dis_xs_p+dis_xs_n;
        //sigma_al  = (dis_xs_p * 13 + dis_xs_n * 14);
		//[>}}}<]

		////if(W>2.0&&xb<1.00){//Only consider the DIS events with W>2.0 cut, for high-x poins, we may need smaller W cut.
		//if((W>2.0&&xb<0.70)||(W>sqrt(3.0)&&xb>0.70)){//Only consider the DIS events with W>2.0 cut, for high-x poins, we may need smaller W cut.
			//D2_XS_Sum += sigma_d2 * NBARN_TO_CM2;//cm2/GeV/sr
			//H3_XS_Sum += sigma_h3 * NBARN_TO_CM2;//cm2/GeV/sr
			//He3_XS_Sum += sigma_he3 * NBARN_TO_CM2;//cm2/GeV/sr
			//XS_Sum_Window += sigma_al * NBARN_TO_CM2;//Assuming it is pure Al, Z = 13, N = 14
		//}
		//Count++;

	//}
	//[>}}}<]

	//[>Calculate rates{{{<]
	//double He3_Target_Rate = He3_Target_Lumi * He3_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
	//double H3_Target_Rate = H3_Target_Lumi * H3_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
	//double D2_Target_Rate = D2_Target_Lumi * D2_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
	//double Window_Rate = Window_Lumi * XS_Sum_Window * Phase_Space / Nevnt;  //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz

	//xs_d2_avg[0] = D2_XS_Sum/Count / NBARN_TO_CM2;
	//xs_h3_avg[0] = H3_XS_Sum/Count / NBARN_TO_CM2;
	//xs_he3_avg[0] = He3_XS_Sum/Count / NBARN_TO_CM2;
	//[>}}}<]

	//rate_d2[0] = D2_Target_Rate;
	//rate_h3[0] = H3_Target_Rate;
	//rate_he3[0] = He3_Target_Rate;
	//rate_window[0] = Window_Rate;

	//delete evt;
/*}*/
/*}}}GetRate*/

/*GetRate{{{*/
void GetRate(const string spec, const string Model,  const string Com, double E0, double P0, double T0, 
		double* xb_avg, double* W2_avg, double* xs_he3,double* xs_he3_avg,double* rate_he3, double* xs_h3,double* xs_h3_avg,double* rate_h3, double* xs_d2,double* xs_d2_avg,double* rate_d2, double *rate_window){
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
    dis->SetModel(Model.c_str());
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

	/*Define ROOT file and tree{{{*/
    TString filename = Form("./results/%s_rate_marathon_%s_E%d_P%d_T%d.root", Model.c_str(),spec.c_str(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
    if(IsDebug)
        filename = Form("./results/%s_rate_marathon_%s_E%d_P%d_T%d_debug.root", Model.c_str(),spec.c_str(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
    if(spec=="BB" &&Com=="narrow")	
        filename = Form("./results/%s_rate_marathon_%s_E%d_P%d_T%d_narrow.root",Model.c_str(), spec.c_str(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
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
    T->Branch("xs_p",&xs_p, "xs_p/D");
    T->Branch("xs_n",&xs_n, "xs_n/D");
    T->Branch("xs_d",&xs_d, "xs_d/D");

    T->Branch("sigma0_d2",&sigma0_d2, "sigma0_d2/D");
    T->Branch("sigma_d2",&sigma_d2, "sigma_d2/D");
    T->Branch("sigma0_h3",&sigma0_h3, "sigma0_h3/D");
    T->Branch("sigma_he3",&sigma_he3, "sigma_he3/D");
    T->Branch("Phase_Space",&Phase_Space, "Phase_Space/D");

	/*}}}*/

	/*Start sampling events{{{*/
    int Count=0;  double xb_sum=0.0; double W2_sum=0.0;
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
                xb_sum+=xb* sigma_h3 * NBARN_TO_CM2;
                W2_sum+=W*W* sigma_h3 * NBARN_TO_CM2;
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
    
        //xb_avg[0] = xb_sum/Count;
        //W2_avg[0] = W2_sum/Count;
        W2_avg[0] = W2_sum/H3_XS_Sum;
        xb_avg[0] = xb_sum/H3_XS_Sum;
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
    delete dis;
}
/*}}}GetRate*/
/*main{{{*//*{{{*//*}}}*/
int main(){

    const int NBin = 17;
    const double E0 = 11.0;
    const double Evt_Bin = 25000;
    //const double P0[NBin] = {4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   2.07};
    const double P0[NBin] = {3.6,   3.6,   3.6,   3.6,   3.6,   3.6,   3.6,   3.6,   3.6,   3.6,   3.6,   3.6,   3.6,   3.6,   3.6,  3.6,   3.20};
	const double X_Bin[NBin] = { 0.23, 0.27,0.31, 0.35,0.39,0.43,0.47,0.51,0.55, 0.59, 0.63, 0.67, 0.71, 0.75, 0.79, 0.83, 0.87};

    double T0,Q2, W2, xbj;
    double rate_he3, rate_h3, rate_d2, rate_w, hour_he3, hour_h3, hour_d2;
    double total_hour,xs_he3, xs_h3, xs_d2,xs_he3_avg, xs_h3_avg, xs_d2_avg, xbj_avg,W2_avg;
    double total_day=0.0;
    string model = ""; 
    int mod = 0; cout<<"-- Which model (1->CJ, 2->F1F2IN09, 3->F2ALLM or 4->PB) "; cin >> mod;
    if(mod==1) model="CJ";
    else if(mod==2) model="F1F2IN09";
    else if(mod==3) model="F2ALLM";
    else if(mod==4) model="PB";
    else {cerr<<"*** Wrong model#!!!!"<<endl; return -111;}

    FILE *new_table; 
    new_table =fopen(Form("hrs_table_du_%s_E3.6_lowW2.dat",model.c_str()),"w");

    //only with HRS
    const string spec = "HRS";
    for(int i=0;i<NBin; i++){/*{{{*/
        double nu = E0 - P0[NBin-1-i];
        xbj = X_Bin[NBin-1-i];
        Q2 = 2.0 * PROTON_MASS * nu * xbj;
        double SinSQ = Q2 / 4.0 / E0 / P0[NBin-1-i];
        T0 = 2.0 * asin(sqrt(SinSQ)) / Deg2Rad;

        if(abs(xbj-0.87)<0.001)
            T0 = 47.10;

        W2 = PROTON_MASS*PROTON_MASS + 2.0*PROTON_MASS*nu - Q2;

        GetRate(spec.c_str(), model.c_str(),  "",  E0, P0[NBin-1-i], T0,&xbj_avg,&W2_avg, &xs_he3, &xs_he3_avg, &rate_he3, &xs_h3,&xs_h3_avg, &rate_h3, &xs_d2,&xs_d2_avg, &rate_d2, &rate_w);

        printf("--- Xbj = %f, P0=%f, T0=%f Rate_He3 = %f\n", xbj, P0[NBin-1-i], T0, rate_he3);

        //Acount for 10% loss of using SOS quards/*{{{*/
        if(spec=="HRS"){
            rate_he3 *=0.9;
            rate_h3 *=0.9;
            rate_d2 *=0.9;
           //80% due to RC, 80% due to efficiecies
            rate_he3 *=0.8*0.8;
            rate_h3 *=0.8*0.8;
            rate_d2 *=0.8*0.8;
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

        fprintf(new_table,"%4s & %4.2f(%4.2f) & %4.2f(%4.2f) & %4.2f & %4.2f  & %4.2f & %6.4f &  %6.4f & %6.4f & %6.3f(%4.2f) & %6.3f(%4.2f) & %6.3f(%4.2f) & %5.2f \\\\\\hline \n",
                "HRS", xbj,xbj_avg, W2,W2_avg, Q2, P0[NBin-1-i], T0,xs_d2_avg, xs_h3_avg, xs_he3_avg, rate_d2,hour_d2, rate_h3,hour_h3,  rate_he3,hour_he3, total_hour);
        total_day +=total_hour;
    }
    /*}}}*/
    total_day /= 24.0;
    fprintf(new_table, " &  &  &  &  &  &  &  &  &  &  &  & %5.2f\\\\\\hline \n", total_day);
    fclose(new_table);
}
/*}}}*/

/*main{{{*/
/*int main(){*/
        //string spec;
        //const double E0 = 11.0;
        //const double Evt_Bin = 25000.;//25K events per data point
	//double P0, T0, xbj, Q2, W2;

	//double xs_he3, xs_h3, xs_d2, rate_he3, rate_h3, rate_d2;
	//double rate_he3_bb, rate_h3_bb, rate_d2_bb, rate_w, rate_w_bb;
	//double xs_h3_avg, xs_he3_avg,xs_d2_avg;
	////int n_he3_new, n_h3_new, n_d2_new;
	//double hour_d2, hour_d2_bb, hour_he3, hour_he3_bb, hour_h3, hour_h3_bb;

        //cout <<"---- Spectrometer (BB or HRS) "; cin >> spec;
        //cout <<"---- E0 = 11~GeV"<<endl;
        //cout <<"---- Ep (GeV) = "; cin >> P0;
        //cout <<"---- Theta (Degree) = "; cin >> T0;

	//double SinSQ = pow(sin(T0*Deg2Rad/2.), 2);
	//double nu = E0 - P0;
	//Q2 = 4.0 * E0 * P0 *SinSQ;
	//W2 = PROTON_MASS*PROTON_MASS + 2.0*PROTON_MASS*nu - Q2;
    //xbj = Q2 / 2.0 / nu / PROTON_MASS;

	//GetRateAll(spec.c_str(), E0, P0, T0, "", &xs_he3, &xs_he3_avg, &rate_he3, &xs_h3,&xs_h3_avg, &rate_h3, &xs_d2,&xs_d2_avg, &rate_d2, &rate_w);

    //printf("--- %s: Rate_He3 = %f,  Rate_H3 = %f\n",  spec.c_str(), rate_he3, rate_h3);

	//rate_he3_bb = 0.0;
	//rate_h3_bb = 0.0;
	//rate_d2_bb = 0.0;
	//if(spec=="BB"){
		//GetRateAll(spec.c_str(), E0, P0, T0, "narrow", &xs_he3,&xs_he3_avg, &rate_he3_bb, &xs_h3,&xs_h3_avg, &rate_h3_bb, &xs_d2,&xs_d2_avg, &rate_d2_bb, &rate_w_bb);
	//}

	////Acount for 10% loss of using SOS quards
	//if(spec=="HRS"){
		//rate_he3 *=0.9;
		//rate_h3 *=0.9;
		//rate_d2 *=0.9;
	//}

	//if(rate_he3>1e-4)
		//hour_he3 =  Evt_Bin/(rate_he3*60.*60.);
	//else
		//hour_he3 = 0;
	//if(rate_h3>1e-4)
		//hour_h3 =  Evt_Bin/(rate_h3*60.*60.);
	//else
		//hour_h3 = 0;
	//if(rate_d2>1e-4)
		//hour_d2 =  Evt_Bin/(rate_d2*60.*60.);
	//else
		//hour_d2 = 0;
	////if(hour_d2<1.0&&hour_d2>0.00001)  hour_d2 = 1.0;
	////if(hour_h3<1.0&&hour_h3>0.00001)  hour_h3 = 1.0;
	////if(hour_he3<1.0&&hour_h3>0.00001) hour_he3 = 1.0;

	//if(rate_he3_bb>1e-4)
		//hour_he3_bb =  25000./(rate_he3_bb*60.*60.);
	//else
		//hour_he3_bb = 0;
	//if(rate_h3_bb>1e-4)
		//hour_h3_bb =  25000./(rate_h3_bb*60.*60.);
	//else
		//hour_h3_bb = 0;
	//if(rate_d2_bb>1e-4)
		//hour_d2_bb =  25000./(rate_d2_bb*60.*60.);
	//else
		//hour_d2_bb = 0;
	////if(hour_d2_bb<1.0&&hour_d2_bb>0.00001)  hour_d2_bb = 1.0;
	////if(hour_h3_bb<1.0&&hour_h3_bb>0.00001)  hour_h3_bb = 1.0;
	////if(hour_he3_bb<1.0&&hour_h3_bb>0.00001) hour_he3_bb = 1.0;

        //cout<<endl<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
	//if(spec=="BB") {
		//printf("&BigBite: E0 = 11 GeV, Ep = %4.2f GeV/c, Theta = %4.2f Degree\n",P0, T0); 
		//printf("          xbj = %4.2f, W2 = %4.2f GeV2, Q2 = %4.2f GeV2:\n", xbj, W2, Q2); 
		//printf("          (Numbers in the brackets are with a momentum cut+/-4.5percent)\n");
		//printf("--------------------------------------------------------------\n");
		//printf("      He3: XS = %6.4f nb/sr/GeV, Rate = %5.2f(%5.2f) Hz, Time = %5.2f(%4.1f) hour\n",xs_he3, rate_he3,rate_he3_bb, hour_he3,hour_he3_bb); 
		//printf("       H3: XS = %6.4f nb/sr/GeV, Rate = %5.2f(%5.2f) Hz, Time = %5.2f(%4.1f) hour\n",xs_h3, rate_h3,rate_h3_bb, hour_h3,hour_h3_bb); 
		//printf("       D2: XS = %6.4f nb/sr/GeV, Rate = %5.2f(%5.2f) Hz, Time = %5.2f(%4.1f) hour\n",xs_d2, rate_d2,rate_d2_bb, hour_d2,hour_d2_bb); 
		//printf("  Windows: Rate = %5.2f(%4.1f) Hz at 20 uA\n", rate_w, rate_w_bb);
        //}
	//else{
		//printf("&&&  HRS: E0 = 11 GeV, Ep = %4.2f GeV/c, Theta = %4.2f Degree\n", P0, T0 ); 
		//printf("          xbj = %4.2f, W2 = %4.2f GeV2, Q2 = %4.2f GeV2:\n", xbj, W2, Q2); 
		//printf("--------------------------------------------------------------\n");
		//printf("      He3: XS = %6.4f nb/sr/GeV, Rate = %5.2f Hz, Time = %4.1f hour\n",xs_he3, rate_he3,hour_he3); 
		//printf("       H3: XS = %6.4f nb/sr/GeV, Rate = %5.2f Hz, Time = %4.1f hour\n",xs_h3, rate_h3, hour_h3); 
		//printf("       D2: XS = %6.4f nb/sr/GeV, Rate = %5.2f Hz, Time = %4.1f hour\n",xs_d2, rate_d2, hour_d2); 
		//printf("  Windows: Rate = %5.2f Hz at 20 uA\n", rate_w);
	//}
        //cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl<<endl;
/*}*/
/*}}}*/
