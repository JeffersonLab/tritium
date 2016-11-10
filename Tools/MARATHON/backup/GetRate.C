#include "GetRate.h"
#include "christy_bosted_inelastic.h"

/*GetRateAll{{{*/
void GetRateAll(const TString spec, const TString Target, double E0, double P0, double T0, const TString Com, 
		double* xs_he3,double* xs_he3_avg,double* rate_he3, double* xs_h3,double* xs_h3_avg,double* rate_h3, double* xs_d2,double* xs_d2_avg,double* rate_d2, double *rate_window){
	//bool IsDebug = kTRUE;
	bool IsDebug = kFALSE;
	/*
	   double E0 = 11.0; //GeV
	   double P0; //GeV, Scattered Electrons' Central Momentum
	   double T0; //Degree, Scattered Electrons' Central Angle
	   TString spec = "";
	   TString Target = "";
	   cout<<" --- Which Spectrometer (HRS or BB) = "; cin >> spec;
	   cout<<" --- Which Target (H2, D2, H3, He3, or C12) = "; cin >> Target;
	   cout<<" --- Beam Energy (E0, GeV) = "; cin >> E0;
	   cout<<" --- Scattered Electron's Central Momentum (P0, GeV) = "; cin >> P0;
	   cout<<" --- Scattered Electron's Central Angle (Theta0, Deg) = "; cin >> T0;
	   */
	int Nevnt = 100000;
	if(IsDebug){cout<<" ----Total Samples = "; cin >> Nevnt;}
	cout<<endl<<endl;

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

	//else if(Target=="D2"){
	Beam_Current = 20*1e-6;//25uA
	Target_Thick = Window_Length;
	Target_Density = 0.124/25.00;//gcm-3

	double Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
	double D2_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/2;//cm-2*s-1
	double sigma0_d2 = (sigma_d(E0, T0*Deg2Rad, P0)) ;

	//else if(Target=="H3"){
	Beam_Current = 20*1e-6;//25uA
	Target_Thick = Window_Length;
	Target_Density = 0.081/25.00;//g*cm-3

	double H3_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/3;//cm-2*s-1
	double sigma0_h3 = (sigma_d( E0, T0*Deg2Rad, P0)+sigma_n(E0, T0*Deg2Rad, P0)) ;

	//else if(Target=="He3"){
	Beam_Current = 25*1e-6;//25uA
	Target_Thick = Window_Length;
	Target_Density = 0.093/25.00;//g*cm-3

	double He3_Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/3;//cm-2*s-1
	double sigma0_he3 = (sigma_d( E0, T0*Deg2Rad, P0)+sigma_p(E0, T0*Deg2Rad, P0)) ;

	xs_d2[0] = sigma0_d2;
	xs_h3[0] = sigma0_h3;
	xs_he3[0] = sigma0_he3;
	/*}}}*/

	const double CosA = cos( T0*Deg2Rad );
	const double SinA = sin( T0*Deg2Rad );
	double x_tg, y_tg, z_tg, dp_tg,theta_tg, phi_tg,Theta, Ep, sigma_d2, sigma_h3, sigma_he3;
	double xs_p,xs_n,xs_d;
	double xs_p0,xs_n0,xs_d0;
	double xb,nu, Q2, W, eta;   
	double D2_XS_Sum = 0.0,H3_XS_Sum = 0.0,He3_XS_Sum = 0.0;
	double XS_Sum_Window = 0.0;

	/*Define ROOT file and tree{{{*/
	TString filename = Form("./results/rate_marathon_%s_E%d_P%d_T%d.root", spec.Data(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
	if(IsDebug)
		filename = Form("./results/rate_marathon_%s_E%d_P%d_T%d_debug.root", spec.Data(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
	if(spec=="BB" &&Com=="narrow")	
		filename = Form("./results/rate_marathon_%s_E%d_P%d_T%d_narrow.root", spec.Data(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
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
	T->Branch("eta",&eta, "eta/D");

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
	T->Branch("xs_p0",&xs_p0, "xs_p0/D");
	T->Branch("xs_n0",&xs_n0, "xs_n0/D");
	T->Branch("xs_d0",&xs_d0, "xs_d0/D");
	T->Branch("sigma0_d2",&sigma0_d2, "sigma0_d2/D");
	T->Branch("sigma_d2",&sigma_d2, "sigma_d2/D");
	T->Branch("sigma0_h3",&sigma0_h3, "sigma0_h3/D");
	T->Branch("sigma_he3",&sigma_he3, "sigma_he3/D");
	T->Branch("Phase_Space",&Phase_Space, "Phase_Space/D");

	/*}}}*/

	/*Start sampling events{{{*/
	int Count=0;
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
		double tao = Q2/4.0/PROTON_MASS/PROTON_MASS;
		eta = 1./(1.0 + 2.0*(1.0+tao)*pow(tan(Theta/2.0), 2));

		if(IsDebug)
			cout<<Form("^^^^ x=%8.4fcm, y=%8.5fcm, z=%8.5fcm, Ep=%8.5fGeV, Theta=%8.5fDeg"
					,x_tg,y_tg,z_tg,Ep,Theta/Deg2Rad)<<endl;

		/*Calculate the cross sections{{{*/
		sigma_d2 = 0;
		sigma_h3 = 0;
		sigma_he3 = 0;
		xs_p = sigma_p( E0, Theta, Ep); //nb/GeV/sr
		xs_n = sigma_n( E0, Theta, Ep);//nb/GeV/sr
		xs_d = sigma_d( E0, Theta, Ep);//nb/GeV/sr

		xs_p0 = sigma_p( E0, T0, Ep); //nb/GeV/sr
		xs_n0 = sigma_n( E0, T0, Ep);//nb/GeV/sr
		xs_d0 = sigma_d( E0, T0, Ep);//nb/GeV/sr

		if(isnan(xs_p)||xs_p<-1e-20)
			xs_p=1e-66;
		if(isnan(xs_n)||xs_n<-1e-20)
			xs_n=1e-66;
		if(isnan(xs_d)||xs_d<-1e-20)
			xs_d=1e-66;

		if(IsDebug)
			cout<<Form("^^^^ XS_P = %10.4e, XS_N =  %10.4e", xs_p, xs_n)<<endl;

		sigma_d2 = xs_d;
		sigma_h3 = xs_d+xs_n;
		sigma_he3 = xs_d+xs_p;
		if(IsDebug)
			cout<<Form("^^^^ XS_D = %10.4e", xs_d)<<endl;
		/*}}}*/

		//if(W>2.0&&xb<1.00){//Only consider the DIS events with W>2.0 cut, for high-x poins, we may need smaller W cut.
		if((W>2.0&&xb<0.70)||(W>sqrt(3.0)&&xb>0.70)){//Only consider the DIS events with W>2.0 cut, for high-x poins, we may need smaller W cut.
			D2_XS_Sum += sigma_d2 * NBARN_TO_CM2;//cm2/GeV/sr
			H3_XS_Sum += sigma_h3 * NBARN_TO_CM2;//cm2/GeV/sr
			He3_XS_Sum += sigma_he3 * NBARN_TO_CM2;//cm2/GeV/sr
			XS_Sum_Window += (xs_p * 13 + xs_n * 14) * NBARN_TO_CM2;//Assuming it is pure Al, Z = 13, N = 14
		    Count++;
		}

		if(IsDebug)
			cout<<Form("XS_Sum=%10.4e, XS_Sum_Window=%10.4e",H3_XS_Sum, XS_Sum_Window)<<endl;
		T->Fill();
	}
	/*}}}*/

	/*Calculate rates{{{*/
	double He3_Target_Rate = He3_Target_Lumi * He3_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
	double H3_Target_Rate = H3_Target_Lumi * H3_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
	double D2_Target_Rate = D2_Target_Lumi * D2_XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
	double Window_Rate = Window_Lumi * XS_Sum_Window * Phase_Space / Nevnt;  //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz

	xs_d2_avg[0] = D2_XS_Sum/Count;
	xs_h3_avg[0] = H3_XS_Sum/Count;
	xs_he3_avg[0] = He3_XS_Sum/Count;
	
	cout<<Form("&&&& Phase_Space = %6.4f", Phase_Space)<<endl;
	cout<<Form("&&&& Target Rate = %8.2f Hz,  Window Rate = %8.2f Hz, Total Rate = %8.2f Hz", He3_Target_Rate, Window_Rate, He3_Target_Rate+Window_Rate)<<endl; 
	/*}}}*/

	rate_d2[0] = D2_Target_Rate;
	rate_h3[0] = H3_Target_Rate;
	rate_he3[0] = He3_Target_Rate;
	rate_window[0] = Window_Rate;

	file->cd(); T->Write(); file->Close();
	delete evt;
}
/*}}}GetRate*/

/*GetRate{{{*/
void GetRate(const TString spec, const TString Target, double E0, double P0, double T0, const TString Com, double* xs_central,double* rate, double *rate_window){
	//bool IsDebug = kTRUE;
	bool IsDebug = kFALSE;
	/*
	   double E0 = 11.0; //GeV
	   double P0; //GeV, Scattered Electrons' Central Momentum
	   double T0; //Degree, Scattered Electrons' Central Angle
	   TString spec = "";
	   TString Target = "";
	   cout<<" --- Which Spectrometer (HRS or BB) = "; cin >> spec;
	   cout<<" --- Which Target (H2, D2, H3, He3, or C12) = "; cin >> Target;
	   cout<<" --- Beam Energy (E0, GeV) = "; cin >> E0;
	   cout<<" --- Scattered Electron's Central Momentum (P0, GeV) = "; cin >> P0;
	   cout<<" --- Scattered Electron's Central Angle (Theta0, Deg) = "; cin >> T0;
	   */
	int Nevnt = 100000;
	if(IsDebug){cout<<" ----Total Samples = "; cin >> Nevnt;}
	cout<<endl<<endl;

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

	double Window_Lumi = 0.0;
	double Target_Density = 0.0;//g*cm-3
	double Target_Thick = 0.0;//cm
	double Target_Lumi = 0.0;
	double sigma0 = 0.0;

	int A = 0, Z =0;
	if(Target=="H2"){
		A = 1; Z = 1;
		Beam_Current = 20.0*1e-6;//A, 25uA
		Target_Thick = Window_Length;
		Target_Density = 0.062/25.00;//gcm-3

		Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
		Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/A;//cm-2*s-1

		evt->SetVertexRange(-Target_Thick/2.0, Target_Thick/2.0);

		sigma0 = (sigma_p( E0, T0*Deg2Rad, P0)) ;
		cout<<Form("--- Central XS = %10.4f nb/GeV/sr", sigma0 )<<endl;	   
	}
	else if(Target=="D2"){
		A = 2; Z = 1;
		Beam_Current = 20*1e-6;//25uA
		Target_Thick = Window_Length;
		Target_Density = 0.124/25.00;//gcm-3

		Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
		Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/A;//cm-2*s-1

		sigma0 = (sigma_d(E0, T0*Deg2Rad, P0)) ;

		cout<<Form("--- Central XS = %10.4f nb/GeV/sr", sigma0 )<<endl;	   
		evt->SetVertexRange(-Target_Thick/2.0, Target_Thick/2.0);
	}
	else if(Target=="H3"){
		A = 3; Z = 1;
		Beam_Current = 20*1e-6;//25uA
		Target_Thick = Window_Length;
		Target_Density = 0.081/25.00;//g*cm-3
		//Target_Density = 0.00324;//g*cm-3

		Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
		Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/A;//cm-2*s-1

		evt->SetVertexRange(-Target_Thick/2.0, Target_Thick/2.0);

		sigma0 = (sigma_d( E0, T0*Deg2Rad, P0)+sigma_n(E0, T0*Deg2Rad, P0)) ;
		cout<<Form("--- Central XS = %10.4f nb/GeV/sr", sigma0 )<<endl;	   
	}
	else if(Target=="He3"){
		A = 3; Z = 2;
		Beam_Current = 25*1e-6;//25uA
		Target_Thick = Window_Length;
		Target_Density = 0.093/25.00;//g*cm-3

		Window_Lumi= Beam_Current/NQe * Window_ArialDensity*Avogadro/27;//cm-2*s-1
		Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/A;//cm-2*s-1

		evt->SetVertexRange(-Target_Thick/2.0, Target_Thick/2.0);

		sigma0 = (sigma_d( E0, T0*Deg2Rad, P0)+sigma_p(E0, T0*Deg2Rad, P0)) ;
		cout<<Form("--- Central XS = %10.4f nb/GeV/sr", sigma0 )<<endl;	   
	}
	else if(Target=="C12"){
		A = 12; Z = 6;
		Beam_Current = 25*1e-6;//25uA
		Target_Density = 0;//g*cm-3
		Target_Thick = 0;

		Window_Lumi= 0.0;
		Target_Lumi= Beam_Current/NQe * Target_Thick*Target_Density*Avogadro/A;//cm-2*s-1, per nucleus instead of per-nucleon

		evt->SetVertexRange(-Target_Thick/2.0, Target_Thick/2.0);

		sigma0 =(Z*sigma_p( E0, T0*Deg2Rad, P0)+(A-Z)*sigma_n(E0, T0*Deg2Rad, P0));	   
		cout<<Form("--- Central XS = %10.4f nb/GeV/sr", sigma0 )<<endl;	   
	}
	else{
		cerr<<"***ERROR, I don't know this target!!!"<<endl;
	}
	cerr<<Form("--- For %s on %s: E0 = %4.2f, Ep = %4.2f,  Theta=%4.2f", Target.Data(), spec.Data(), E0, P0, T0)<<endl;
	cout<<Form("&&&& Target_Lumi=%10.4e, Window_Lumi=%10.4e",Target_Lumi, Window_Lumi)<<endl;
	cout<<endl<<endl;
	xs_central[0] = sigma0;
	/*}}}*/

	const double CosA = cos( T0*Deg2Rad );
	const double SinA = sin( T0*Deg2Rad );
	double x_tg, y_tg, z_tg, dp_tg,theta_tg, phi_tg,Theta, Ep, sigma;
	double xs_p,xs_n,xs_d;
	double xs_p0,xs_n0,xs_d0;
	double xb,nu, Q2, W, eta;   
	double XS_Sum = 0.0;
	double XS_Sum_Window = 0.0;

	/*Define ROOT file and tree{{{*/
	TString filename = Form("./results/rate_marathon_%s_%s_E%d_P%d_T%d.root", spec.Data(), Target.Data(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
	if(IsDebug)
		filename = Form("./results/rate_marathon_%s_%s_E%d_P%d_T%d_debug.root", spec.Data(), Target.Data(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
	if(spec=="BB" &&Com=="narrow")	
		filename = Form("./results/rate_marathon_%s_%s_E%d_P%d_T%d_narrow.root", spec.Data(), Target.Data(), (int)(E0+0.5), (int)(P0+0.5), (int)(T0+0.5));
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
	T->Branch("eta",&eta, "eta/D");

	T->Branch("Target_Lumi",&Target_Lumi, "Target_Lumi/D");
	T->Branch("Window_Lumi",&Window_Lumi, "Window_Lumi/D");
	T->Branch("XS_Sum",&XS_Sum, "XS_Sum/D");
	T->Branch("XS_Sum_Window",&XS_Sum_Window, "XS_Sum_Window/D");
	T->Branch("xs_p",&xs_p, "xs_p/D");
	T->Branch("xs_n",&xs_n, "xs_n/D");
	T->Branch("xs_d",&xs_d, "xs_d/D");
	T->Branch("xs_p0",&xs_p0, "xs_p0/D");
	T->Branch("xs_n0",&xs_n0, "xs_n0/D");
	T->Branch("xs_d0",&xs_d0, "xs_d0/D");
	T->Branch("sigma0",&sigma0, "sigma0/D");
	T->Branch("sigma",&sigma, "sigma/D");
	T->Branch("Phase_Space",&Phase_Space, "Phase_Space/D");

	/*}}}*/

	/*Start sampling events{{{*/
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
		double tao = Q2/4.0/PROTON_MASS/PROTON_MASS;
		eta = 1./(1.0 + 2.0*(1.0+tao)*pow(tan(Theta/2.0), 2));

		if(IsDebug)
			cout<<Form("^^^^ x=%8.4fcm, y=%8.5fcm, z=%8.5fcm, Ep=%8.5fGeV, Theta=%8.5fDeg"
					,x_tg,y_tg,z_tg,Ep,Theta/Deg2Rad)<<endl;

		/*Calculate the cross sections{{{*/
		sigma = 0;
		xs_p = sigma_p( E0, Theta, Ep); //nb/GeV/sr
		xs_n = sigma_n( E0, Theta, Ep);//nb/GeV/sr
		xs_d = sigma_d( E0, Theta, Ep);//nb/GeV/sr

		xs_p0 = sigma_p( E0, T0, Ep); //nb/GeV/sr
		xs_n0 = sigma_n( E0, T0, Ep);//nb/GeV/sr
		xs_d0 = sigma_d( E0, T0, Ep);//nb/GeV/sr

		if(isnan(xs_p)||xs_p<-1e-20)
			xs_p=1e-66;
		if(isnan(xs_n)||xs_n<-1e-20)
			xs_n=1e-66;
		if(isnan(xs_d)||xs_d<-1e-20)
			xs_d=1e-66;

		if(IsDebug)
			cout<<Form("^^^^ XS_P = %10.4e, XS_N =  %10.4e", xs_p, xs_n)<<endl;

		if(A==1)
			sigma = xs_p;
		else if(A>3)
			sigma = Z*xs_p + (A-Z)*xs_n;
		else{
			if(Target=="D2")
				sigma = xs_d;
			if(Target=="H3")
				sigma = xs_d+xs_n;
			if(Target=="He3")
				sigma = xs_d+xs_p;
			if(IsDebug)
				cout<<Form("^^^^ XS_D = %10.4e", xs_d)<<endl;
		}
		/*}}}*/

		//if(W>2.0&&xb<1.00){//Only consider the DIS events with W>2.0 cut, for high-x poins, we may need smaller W cut.
		if((W>2.0&&xb<0.70)||(W>sqrt(3.0)&&xb>0.70)){//Only consider the DIS events with W>2.0 cut, for high-x poins, we may need smaller W cut.
			XS_Sum += sigma * NBARN_TO_CM2;//cm2/GeV/sr
			XS_Sum_Window += (xs_p * 13 + xs_n * 14) * NBARN_TO_CM2;//Assuming it is pure Al, Z = 13, N = 14
		}

		if(IsDebug)
			cout<<Form("XS_Sum=%10.4e, XS_Sum_Window=%10.4e",XS_Sum, XS_Sum_Window)<<endl;
		T->Fill();
	}
	/*}}}*/

	/*Calculate rates{{{*/
	double Target_Rate = Target_Lumi * XS_Sum * Phase_Space / Nevnt; //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz
	double Window_Rate = Window_Lumi * XS_Sum_Window * Phase_Space / Nevnt;  //cm-2*s-1 * cm2/GeV/sr * sr*GeV = s-1 = Hz

	cout<<Form("&&&& Phase_Space = %6.4f", Phase_Space)<<endl;

	cout<<Form("&&&& Target Rate = %8.2f Hz,  Window Rate = %8.2f Hz, Total Rate = %8.2f Hz", Target_Rate, Window_Rate, Target_Rate+Window_Rate)<<endl; 
	/*}}}*/

	rate[0] = Target_Rate;
	rate_window[0] = Window_Rate;

	file->cd(); T->Write(); file->Close();
	delete evt;
}
/*}}}GetRate*/

/*main{{{*//*{{{*//*}}}*/
int main(){

    const int NBin = 18;
    const double E0 = 11.0;
    const double P0[NBin] = {4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   2.07};
    const double T0[NBin] = {15.10, 16.30, 17.50, 18.60, 19.60, 20.60, 21.60, 22.50, 23.40, 24.20, 25.00, 25.90, 26.70, 27.40, 28.10, 28.90, 47.10};
    double Ep, Theta, xbj;
    int rate_he3, rate_h3, rate_d2, rate_w, hour_he3, hour_h3, hour_d2;

    ofstream new_table("hrs_table_du_new.dat");

    double xs_he3, xs_h3, xs_d2,xs_he3_avg, xs_h3_avg, xs_d2_avg, rate_he3, rate_h3, rate_d2, rate_w;

    //only with HRS
    const TString spec = "HRS";
    for(int i=0;i<NBin; i++){

        GetRateAll(spec.Data(), "He3",E0, P0[i], T0[i], "", &xs_he3, &xs_he3_avg, &rate_he3, &xs_h3,&xs_h3_avg, &rate_h3, &xs_d2,&xs_d2_avg, &rate_d2, &rate_w);

        new_table << Form(" %4.2f & %4.2f & %4.2f & %4.2f &
                %6.4f &  %6.4f & %6.4f &
                %6.3f(%4.2f) & %6.3f(%4.2f) & %6.3f(%4.2f) & %6.3f",
                xbj, W2, Q2, P0[i], T0[i], 
                xs_he3_avg, xs_d2_avg, xs_h3_avg, rate_he3,hour_he3, rate_d2,hour_d2, rate_h3,hour_h3, rate_w) 
            <<endl;
    }

}
/*}}}*/


/*main{{{*//*{{{*//*}}}*/
/*int main(){*/
    //ifstream table("table_du.dat");
    //TString com, spec;
    //table >> com >> com >> com >> com >> com >> com >> com >> com >> com >> com >> com >> com >> com;

    //double E0, P0, T0, xbj;
    //double xs_he3_old, xs_h3_old;
    //double rate_he3_old, rate_h3_old, rate_d2_old;
    //int n_he3, n_h3, n_d2, bt_he3, bt_h3, bt_d2;

    //ofstream new_table("hrs_table_du_w3.dat");
    //ofstream new_table_he3("hrs_table_du_he3_w3.dat");
    //ofstream new_table_h3("hrs_table_du_h3_w3.dat");

    //double xs_he3, xs_h3, xs_d2,xs_he3_avg, xs_h3_avg, xs_d2_avg, rate_he3, rate_h3, rate_d2, rate_w;
    //double rate_he3_bb,rate_h3_bb,rate_d2_bb, rate_w_bb;
    ////int n_he3_new, n_h3_new, n_d2_new;

    //while(!(table.eof())){
        //table >> E0 >> P0 >> T0 >> spec >> xbj >> xs_he3_old >> xs_h3_old >> n_he3 >> n_h3 >> n_d2 >> bt_he3 >> bt_h3 >> bt_d2; 
        //rate_he3_old = n_he3*1.0/60./60.;
        //rate_h3_old = n_h3*1.0/60./60.;
        //rate_d2_old = n_d2*1.0/60./60.;

        ////only with HRS
        //spec = "HRS";

        //GetRateAll(spec.Data(), "He3",E0, P0, T0, "", &xs_he3, &xs_he3_avg, &rate_he3, &xs_h3,&xs_h3_avg, &rate_h3, &xs_d2,&xs_d2_avg, &rate_d2, &rate_w);
        ////GetRate(spec.Data(), "He3",E0, P0, T0, "", &xs_he3, &rate_he3, &rate_w);
        ////GetRate(spec.Data(), "H3" ,E0, P0, T0, "", &xs_h3, &rate_h3, &rate_w);
        ////GetRate(spec.Data(), "D2" ,E0, P0, T0, "", &xs_d2, &rate_d2, &rate_w);

        //if(spec=="BB"){
            //GetRateAll(spec.Data(), "He3",E0, P0, T0, "narrow", &xs_he3,&xs_he3_avg, &rate_he3_bb, &xs_h3,&xs_h3_avg, &rate_h3_bb, &xs_d2,&xs_d2_avg, &rate_d2_bb, &rate_w_bb);
            ////GetRate(spec.Data(), "He3",E0, P0, T0, "narrow", &xs_he3, &rate_he3_bb, &rate_w_bb);
            ////GetRate(spec.Data(), "H3" ,E0, P0, T0, "narrow", &xs_h3, &rate_h3_bb, &rate_w_bb);
            ////GetRate(spec.Data(), "D2" ,E0, P0, T0, "narrow", &xs_d2, &rate_d2_bb, &rate_w_bb);
        //}

        //if(spec=="BB") {
            //new_table << Form("%4s  %4.2f  %4.2f  %6.2f  %6.2f  %6.4f/%6.4f/%6.4f  %6.4f/%6.4f/%6.4f  %5.2f/%5.2f(%4.2f)  %5.2f/%5.2f(%4.2f)  %5.2f/%5.2f(%4.2f)  %6.2f",
                    //spec.Data(),xbj, E0, P0, T0, xs_he3_old, xs_he3, xs_he3_avg, xs_h3_old, xs_h3, xs_h3_avg, rate_he3_old, rate_he3, rate_he3_bb, rate_h3_old, rate_h3, rate_h3_bb, rate_d2_old, rate_d2, rate_d2_bb, rate_w) 
                //<<endl;
            //new_table_he3 << Form("%4s & %4.2f & %4.2f &  %6.2f  & %6.2f  &  %6.4f&%6.4f &%6.4f &  %5.2f & %5.2f(%4.2f) &  %6.2f\\",
                    //spec.Data(),xbj, E0, P0, T0, xs_he3_old, xs_he3, xs_he3_avg, rate_he3_old, rate_he3, rate_he3_bb, rate_w*25./20.) 
                //<<endl;
            //new_table_h3 << Form("%4s & %4.2f & %4.2f &  %6.2f  & %6.2f  &  %6.4f&%6.4f &%6.4f &  %5.2f & %5.2f(%4.2f) &  %6.2f\\",
                    //spec.Data(),xbj, E0, P0, T0, xs_h3_old, xs_h3, xs_h3_avg, rate_h3_old, rate_h3, rate_h3_bb, rate_w) 
                //<<endl;
        //}
        //else{
            //new_table << Form("%4s  %4.2f  %4.2f  %6.2f  %6.2f  %6.4f/%6.4f/%6.4f  %6.4f/%6.4f/%6.4f  %6.2f/%6.2f  %6.2f/%6.2f  %6.2f/%6.2f  %6.2f",
                    //spec.Data(),xbj, E0, P0, T0, xs_he3_old, xs_he3, xs_he3_avg, xs_h3_old, xs_h3, xs_h3_avg, rate_he3_old, rate_he3, rate_h3_old, rate_h3, rate_d2_old, rate_d2, rate_w) 
                //<<endl;
            //new_table_he3 << Form("%4s & %4.2f & %4.2f &  %6.2f  & %6.2f  &  %6.4f&  %6.4f&%6.4f &  %5.2f & %5.2f &  %6.2f\\",
                    //spec.Data(),xbj, E0, P0, T0, xs_he3_old, xs_he3, xs_he3_avg, rate_he3_old, rate_he3,  rate_w*25./20.) 
                //<<endl;
            //new_table_h3 << Form("%4s & %4.2f & %4.2f &  %6.2f  & %6.2f  &  %6.4f&  %6.4f&%6.4f &  %5.2f & %5.2f &  %6.2f\\",
                    //spec.Data(),xbj, E0, P0, T0, xs_h3_old, xs_h3, xs_h3_avg, rate_h3_old, rate_h3, rate_w) 
                //<<endl;
        //}
    //}

/*}*/
/*}}}*/
