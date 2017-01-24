#include "XEM_Target.h"
inline void gCal_DIS(double kE0,double kEp,double kTheta, XEM_TGT* kTarget,double kPmax, int kInnt, int kInnp,XEM_SIG* sig)
{
	/*
	//	modified version of BWF's routine bdis
	//	ag,bg,cg,dg are coeff of a 2 gaussian fit to the f(y)
	//	which was converted to n(k)
	*/
	Int_t A,Z,N;
	Double_t f0,a,b,B,alpha,esep,resol;
	A = kTarget->A;
	Z = kTarget->Z;
	N = A - Z;
	esep = kTarget->ESep;
	resol = kTarget->Resol;
	f0 = kTarget->f0;
	B = kTarget->B;
	a = kTarget->a;
	b = kTarget->b;
	alpha = kTarget->alpha;

	f0*=1000;
	B*=1000;
	alpha/=1000;
	a*=1000;
	b*=1000;

	double agd,bgd,cgd,dgd;
	if(A==2){//d2, FIX_HERE,fake values, need to get the real ones
		agd =    0.39167e+01;
		bgd =    0.78468e+02;
		cgd =    0.17924e+00;
		dgd =    0.66511e+01;
	}
	else if(A==3){//he3
		agd =    0.39167e+01;
		bgd =    0.78468e+02;
		cgd =    0.17924e+00;
		dgd =    0.66511e+01;
	}
	else if(A==4){//	!he4
		agd =   0.35160e+01;
		bgd =   0.60049e+02;
		cgd =   0.38390e+00;
		dgd =   0.13421e+02;
	}
	else if(A==12||A==9){ // !carbon
		agd =   0.28757e+01;
		bgd =   0.41922e+02;
		cgd =   0.33801e+00;
		dgd =   0.13824e+02;
	}
	else if(A==27){//     !alum
		agd =   0.25660e+01;
		bgd =   0.36962e+02;
		cgd =   0.41675e+00;
		dgd =   0.13772e+02;
	}
	else if((A==40)||(A==48)){//     !Calcium, FIX_HERE
		agd =   0.28316e+01;
		bgd =    0.44624e+02;
		cgd =    0.37850e+00;
		dgd =    0.12003e+02;
	}
	else if((A==56)||(A==64)){//     !iron
		agd =   0.28316e+01;
		bgd =    0.44624e+02;
		cgd =    0.37850e+00;
		dgd =    0.12003e+02;
	}
	else if(A==197){//!gold, Jerry, Gold!
		agd =   0.24947e+01;
		bgd =   0.30614e+02;
		cgd =   0.42398e+00;
		dgd =   0.12272e+02;
	}
	else{
		cerr << "****** Wrong Target! Check it! --1"<<endl;
		sig->Value=-1000.0;
		return;
	}

	Double_t CS_Theta = cos(kTheta/2.0);
	Double_t SN_Theta = sin(kTheta/2.0);
	Double_t TN_Theta = tan(kTheta/2.0);
	Double_t SN_SQ = SN_Theta*SN_Theta;
	// Double_t CS_SQ = CS_Theta*CS_Theta;
	Double_t TN_SQ = TN_Theta*TN_Theta;
	Double_t Qsq = 4.0*kE0*kEp*SN_SQ;
	Double_t nu = kE0 - kEp;
	//  Double_t xbj = Qsq/2.0/P_MASS/nu;

	Double_t Mass_A_1; // (A-1) mass 
	if(A==2)
		Mass_A_1=0.93827;
	else if(A==3)  
		Mass_A_1 = 1.87609;
	else if (A==4)  
		Mass_A_1 = 2.8094;
	else if (A==9)  
		Mass_A_1 = 7.5402;
	else if (A==12)  
		Mass_A_1 = 10.25553;
	else if (A==27)  
		Mass_A_1 = 24.205;
	else if (A==40)  
		Mass_A_1 = 36.59253; //FIX_HERE
	else if (A==48)  
		Mass_A_1 = 44.09869; //FIX_HERE
	else if (A==56)  
		Mass_A_1 = 51.1743;
	else if (A==64)  
		Mass_A_1 = 57.6886;
	else if (A==197)  
		Mass_A_1 = 182.5394;
	else{
		cerr << "****** Wrong Target! Check it! ---2"<<endl;
		sig->Value=-1000.0;
		return;
	}
	Double_t Mass_A= Mass_A_1 + P_MASS - esep;
	Double_t Qv = sqrt(Qsq+nu*nu);
	Double_t W1P=0.0,W2P=0.0,W1N=0.0,W2N=0.0;
	Double_t du = 2.0/kInnt;
	Double_t dp = kPmax/kInnp;
	Double_t Wp_Max = 0.0;

	Double_t Ei,W1,W2,W1P1,W2P1,W1N1,W2N1,WpMin2,Wp,Nu;
	Double_t p,u,nu_p,w11,w22,weight1,weight2,arg1,arg2,arg21,arg22,radical;
	XEM_VAR3* f1f2 = new XEM_VAR3();

	for( int ip=1;ip<=kInnp;ip++){
		p = (ip - 0.5)*dp;
		W1P1=0.0;W2P1=0.0;
		W1N1=0.0;W2N1=0.0;

		for(int iu=1;iu<=kInnt;iu++){
			u = (iu-0.5)*du-1.0;
			Ei = Mass_A - sqrt(p*p+Mass_A_1*Mass_A_1);
			nu_p = (Ei*nu - p*u*Qv)/P_MASS;
			WpMin2 = pow((P_MASS+0.135),2);
			radical =(Ei*Ei - p*p + 2.*Ei*nu - 2.*p*u*Qv - Qsq);

			if(radical> WpMin2){
				Wp = sqrt(radical);
				w11=0.0;
				w22=0.0;
				Nu = (Wp*Wp-PM_SQ+Qsq)/2.0/P_MASS;

				gCal_F1F2In06(1,1,Qsq,Wp*Wp,f1f2);
				w11= f1f2->First;  w22 = f1f2->Second; 

				if(!(w11==0)&&!(w22==0)){
					w11/=P_MASS; w22/=Nu;
					W1=w11;  W2=w22;

					if (Wp > Wp_Max) 
						Wp_Max = Wp;

					arg1 = W1+(1.0-u*u)*p*p*W2/2.0/PM_SQ;
					arg21 = pow((1.0+p*u*Qsq/P_MASS/nu_p/Qv),2)
						*nu_p*nu_p/nu/nu;
					arg22 = (1.0-u*u)*p*p*Qsq/2.0/Qv/Qv/PM_SQ;
					arg2 = W2*(arg21 + arg22);
					W1P1 += arg1*du;
					W2P1 += arg2*du;

					w11=0.0;    w22=0.0;
					gCal_F1F2In06(1,0,Qsq,Wp*Wp,f1f2);
					w11= f1f2->First;  w22 = f1f2->Second; 

					if(!(w11==0.0)&&!(w22==0.0)){
						w11/=P_MASS;   w22/=Nu;
						W1=w11; W2=w22;

						arg1 = W1 + (1-u*u)*p*p*W2/2.0/PM_SQ;
						arg21 = pow((1.0+p*u*Qsq/P_MASS/nu_p/Qv),2)
							*nu_p*nu_p/nu/nu;
						arg22 = (1.0-u*u)*p*p*Qsq/2.0/Qv/Qv/PM_SQ;
						arg2 = W2*(arg21 + arg22);
						W1N1 += arg1*du;
						W2N1 += arg2*du;
					}// if(!(w11==0.0)&&!(w22==0.0)){
				}//if(!(w11==0)&&!(w22==0)){
				}//if((Ei*Ei - p*p + 2.*Ei*nu - 2.*p*u*Qv - Qsq)> WpMin2){
				} //for(int iu=1;iu<=kInnt;iu++){

				weight1 = exp(-(a*p)*(a*p))
					* (f0-B)*alpha*alpha
					* (a*a*(alpha*alpha+p*p)+1)
					/ pow((alpha*alpha+p*p),2);

				if (A==2){
					if(p<0) 
						weight1=weight1-(b*B*exp(-b*abs(p))/(2*p));
					else
						weight1=weight1+(b*B*exp(-b*abs(p))/(2*p));
				}
				else
					weight1=weight1+(b*b*B*exp(-(b*p)*(b*p)));

				weight1= weight1*p*p*dp/TMath::Pi();
				weight2 = (agd*bgd/TMath::Pi()*exp(-bgd*p*p) + cgd*dgd/TMath::Pi()*exp(-dgd*p*p))*p*p*dp;

				W1P = W1P + W1P1*weight1;
				W2P = W2P + W2P1*weight1;
				W1N = W1N + W1N1*weight1;
				W2N = W2N + W2N1*weight1;
				}//for( int ip=1;ip<=kInnp;ip++){
				delete f1f2;

				Double_t W1a = 2*TMath::Pi()*(Z*W1P + N*W1N);
				Double_t W2a = 2*TMath::Pi()*(Z*W2P + N*W2N);
				Double_t sigm = CS_Theta*Alpha/(2.0*SN_SQ*kE0);
				sigm = pow(HBARC_XEM*sigm*MeVToGeV,2);

				sig->Value = 1e4*sigm*(W2a+2.0*W1a*TN_SQ);
				sig->Factor=0.0; //No in Use in this subroutine
				return;
				}
