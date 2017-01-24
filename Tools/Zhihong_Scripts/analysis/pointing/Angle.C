const Double_t L = 8.458*1e3; //mm

void Angle(){

	Double_t Scatt_Angle = 25.00;//Deg
	cerr<< " --- Set Angle = "; cin >> Scatt_Angle;
	Double_t L_Z_React = 0.00; //mm
	cerr<< " --- L_Z_Off = "; cin >> L_Z_React;
	Double_t R_Z_React = 0.00; //mm
//	cerr<< " --- R_Z_Off = "; cin >> R_Z_React;
    
   	Double_t L_Scatt_Angle = 0.0; //Deg
	Double_t R_Scatt_Angle = 0.0; //Deg
////Multi-C
//	Double_t D_LHRS = L_Z_React - ( 0.2);
//	Double_t D_RHRS = R_Z_React - (-0.3);

	////C12 Foil
	Double_t D_LHRS = L_Z_React - ( 0.20);
	Double_t D_RHRS = R_Z_React - (-0.29);

    L_Scatt_Angle = Scatt_Angle - (D_LHRS/L)*180/3.1415927;
    R_Scatt_Angle = Scatt_Angle - (D_RHRS/L)*180/3.1415927;
   
   	cerr<<Form("For HRS-L, Scatter_Angle: Set = %f, True = %f", Scatt_Angle, L_Scatt_Angle)<<endl;
    cerr<<Form("For HRS-R, Scatter_Angle: Set = %f, True = %f", Scatt_Angle, R_Scatt_Angle)<<endl;
}
