//L-HRS at -21.480, from DVCS 2010-Nov survey
const Double_t L_Dx = -1.25; //mm
const Double_t L_Dy = 1.91; //mm
const Double_t L_Sieve_Z = -876.30;//mm
const Double_t L_Sieve_x = -1.05; //mm
const Double_t L_Sieve_y =  0.20; //mm

//const Double_t L_Sieve_x = -0.67; //mm
//const Double_t L_Sieve_y = -0.07; //mm

//R-HRS at 20.022, from DVCS 2009-Nov survey
const Double_t R_Dx = -0.73; //mm
const Double_t R_Dy = -3.10; //mm
//From A1n, Jin Ge 2009
const Double_t R_Sieve_Z = -892.28;//mm
const Double_t R_Sieve_x = 1.04; //mm
const Double_t R_Sieve_y = 0.05; //mm

//Beam Position for X>2 Experiment
const Double_t Beam_X = -2.668; //mm
const Double_t Beam_Y = 3.022; //mm


void Pointing(){

	Double_t Scatt_Angle = 25.00;//Deg
	cerr<< " --- Set Angle = "; cin >> Scatt_Angle;
	Double_t L_Z_React = 0.00; //mm
	cerr<< " --- L_Z_Off = "; cin >> L_Z_React;
	Double_t R_Z_React = 0.00; //mm
	cerr<< " --- R_Z_Off = "; cin >> R_Z_React;
 	
	Double_t L_Th = 0.0; //rad
	Double_t L_Ph = 0.0; //rad
	Double_t R_Th = 0.0; //rad
	Double_t R_Ph = 0.0; //rad
    Double_t L_Scatt_Angle = 0.0; //Deg
	Double_t R_Scatt_Angle = 0.0; //Deg

	Double_t L_LHRS = 2058.6 - abs(L_Sieve_Z);
	Double_t L_RHRS = 2068.2 - abs(R_Sieve_Z);
	Double_t CosA = cos(Scatt_Angle*3.1415926/180);
	Double_t SinA = sin(Scatt_Angle*3.1415926/180);

	//LHRS
    L_Ph = (L_Sieve_y + L_Dy   - Beam_X*CosA + L_Z_React*SinA)
		 / (L_LHRS - Beam_X*SinA - L_Z_React*CosA);
	L_Th = (L_Sieve_x + L_Dx - Beam_Y)
		 / (L_LHRS - Beam_X*SinA - L_Z_React*CosA);
	L_Scatt_Angle = acos( (CosA - L_Ph*SinA)/sqrt(1.0+L_Ph*L_Ph + L_Th*L_Th)) 
		          * 180/3.1415926;

	SinA *= -1.0;
    R_Ph = (R_Sieve_y + R_Dy   - Beam_X*CosA + R_Z_React*SinA)
		 / (L_RHRS - Beam_X*SinA - R_Z_React*CosA);
	R_Th = (R_Sieve_x + R_Dx - Beam_Y)
		 / (L_RHRS - Beam_X*SinA - R_Z_React*CosA);
	R_Scatt_Angle = acos( (CosA - R_Ph*SinA)/sqrt(1.0+R_Ph*R_Ph + R_Th*R_Th))
	              * 180/3.1415926;

    cerr<<Form("For HRS-L, Scatter_Angle: Set = %f, True = %f", Scatt_Angle, L_Scatt_Angle)<<endl;
    cerr<<Form("For HRS-R, Scatter_Angle: Set = %f, True = %f", Scatt_Angle, R_Scatt_Angle)<<endl;
}
