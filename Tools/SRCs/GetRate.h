#ifndef GetRate_H
/*Include{{{*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <errno.h>
#include <sstream>
//#include <iterator>
//from root
#include <TMath.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMatrix.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TRandom2.h>
/*}}}*/

/*Constants{{{*/
using namespace::std;
//static const double PI=3.14159265358979323846;
//static const double ALPHA=1/137.036; //fine-structure constant
static const double Deg2Rad = 3.14159265358979323846/180.0;
static const double MeV2GeV = 0.001;
//static const double AP=ALPHA/PI;
//static const double ELECTRON_MASS=0.510998918*MeV2GeV; //GeV
static const double NQe = 1.6*1e-19;
static const double Avogadro = 6.02*1e+23;
//static const double M2=ELECTRON_MASS*ELECTRON_MASS*MeV2GeV; //GeV
//static const double AMU=931.494043*MeV2GeV; //GeV
static const double PROTON_MASS_GeV=938.272029*MeV2GeV; //GeV
//static const double DELTA_MASS=1232.*MeV2GeV; //GeV
//static const double PI_MASS=139.57*MeV2GeV; //GeV
//static const double HBARC=197.326968*MeV2GeV; //hbar*c GeV fm
static const double NBARN_TO_GEV2= 1./3.894 * 1e-5; //nbarn-->GeV^-2
static const double GEV2_TO_NBARN= 3.894 * 1e5; //GeV^-2-->nbarn
static const double NBARN_TO_CM2=1e-33; //nbarn-->CM^2
/*}}}*/

class Event{
	public:
		Event(){
			Init_HRS();
			Init_BB();
		//	fRan = new TRandom2();
		//	fRan->SetSeed(0); 
            seed();
		}
		virtual ~Event(){
		//	delete fRan;
		}

		//double GetTG_Dp(){return TG_Dp_Min + (TG_Dp_Max - TG_Dp_Min) * fRan->Rndm();}//deltaP
		//double GetTG_Theta(){return TG_Theta_Min + (TG_Theta_Max - TG_Theta_Min) * fRan->Rndm();}//rad
		//double GetTG_Phi(){return TG_Phi_Min + (TG_Phi_Max - TG_Phi_Min) * fRan->Rndm();}//rad
		//double GetTG_X(){return TG_X_Min + (TG_X_Max - TG_X_Min) * fRan->Rndm();}//cm
		//double GetTG_Y(){return TG_Y_Min + (TG_Y_Max - TG_Y_Min) * fRan->Rndm();}//cm
		//double GetTG_Z(){return Vertex_Z_Min + (Vertex_Z_Max - Vertex_Z_Min) * fRan->Rndm();}//cm
		double GetTG_Dp(){return TG_Dp_Min + (TG_Dp_Max - TG_Dp_Min) * UnifRand(0.,1.);}//deltaP
		double GetTG_Theta(){return TG_Theta_Min + (TG_Theta_Max - TG_Theta_Min) * UnifRand(0.,1.);}//rad
		double GetTG_Phi(){return TG_Phi_Min + (TG_Phi_Max - TG_Phi_Min) * UnifRand(0.,1.);}//rad
		double GetTG_X(){return TG_X_Min + (TG_X_Max - TG_X_Min) * UnifRand(0.,1.);}//cm
		double GetTG_Y(){return TG_Y_Min + (TG_Y_Max - TG_Y_Min) * UnifRand(0.,1.);}//cm
		double GetTG_Z(){return Vertex_Z_Min + (Vertex_Z_Max - Vertex_Z_Min)  * UnifRand(0.,1.);}//cm


		double GetPhaseSpace(){
			double phase_space = (TG_Theta_Max - TG_Theta_Min) * (TG_Phi_Max - TG_Phi_Min) * (TG_Dp_Max - TG_Dp_Min); 
			return phase_space;	
		}

		void SetSpectrometer(const string spec){
			if(spec.find("LHRS",0)!=string::npos||spec.find("RHRS",0)!=string::npos||spec.find("HRS",0)!=string::npos){
				TG_Dp_Min = HRS_TG_Dp_Min;
				TG_Dp_Max = HRS_TG_Dp_Max;
				TG_Theta_Min = HRS_TG_Theta_Min;
				TG_Theta_Max = HRS_TG_Theta_Max;
				TG_Phi_Min = HRS_TG_Phi_Min;
				TG_Phi_Max = HRS_TG_Phi_Max;
				TG_X_Min = HRS_TG_X_Min;
				TG_X_Max = HRS_TG_X_Max;
				TG_Y_Min = HRS_TG_Y_Min;
				TG_Y_Max = HRS_TG_Y_Max;
			}
			else if(spec.find("BigBite",0)!=string::npos||spec.find("BB",0)!=string::npos){
				TG_Dp_Min = BB_TG_Dp_Min;
				TG_Dp_Max = BB_TG_Dp_Max;
				TG_Theta_Min = BB_TG_Theta_Min;
				TG_Theta_Max = BB_TG_Theta_Max;
				TG_Phi_Min = BB_TG_Phi_Min;
				TG_Phi_Max = BB_TG_Phi_Max;
				TG_X_Min = BB_TG_X_Min;
				TG_X_Max = BB_TG_X_Max;
				TG_Y_Min = BB_TG_Y_Min;
				TG_Y_Max = BB_TG_Y_Max;
			}
			else
				cerr<<"**** Error, I don't know this spectrometer (In GetRate.h, void SetSpectrometer(...))!!!"<<endl;
		}


		bool IsOnTG(double theta, double phi, double dp){
			if(theta<TG_Theta_Min || theta>TG_Theta_Max)
				return 0;
			else if(phi<TG_Phi_Min || phi>TG_Phi_Max)
				return 0;
			else if(dp<TG_Dp_Min || dp>TG_Dp_Max)
				return 0;
			else 
				return 1;   
		}

		bool IsOnFP(double theta, double phi, double x, double y){
			if(theta<FP_Theta_Min || theta>FP_Theta_Max)
				return 0;
			else if(phi<FP_Phi_Min || phi>FP_Phi_Max)
				return 0;
			else if(x<FP_X_Min || x>FP_X_Max)
				return 0;
			else if(y<FP_Y_Min || y>FP_Y_Max)
				return 0;
			else 
				return 1;   
		}

		void SetTGRange(double dp_min, double dp_max, 
				double th_min, double th_max, 
				double ph_min, double ph_max, 
				double x_min, double x_max, 
				double y_min, double y_max){
			TG_Dp_Min = dp_min;
			TG_Dp_Max = dp_max;
			TG_Theta_Min = th_min;
			TG_Theta_Max = th_max;
			TG_Phi_Min = ph_min;
			TG_Phi_Max = ph_max;
			TG_X_Min = x_min;
			TG_X_Max = x_max;
			TG_Y_Min = x_min;
			TG_Y_Max = x_max;
		}

		void SetTGDpRange(double dp_min, double dp_max){
			TG_Dp_Min = dp_min;
			TG_Dp_Max = dp_max;
		}
		void SetTGThetaRange(double th_min, double th_max){
			TG_Theta_Min = th_min;
			TG_Theta_Max = th_max;
		}
		void SetTGPhiRange(double ph_min, double ph_max){
			TG_Phi_Min = ph_min;
			TG_Phi_Max = ph_max;
		}
		void SetVertexRange(double min, double max){
			Vertex_Z_Min = min;
			Vertex_Z_Max = max;
		}

		void SetFPRange(double th_min, double th_max, 
				double ph_min, double ph_max, 
				double x_min, double x_max, 
				double y_min, double y_max){
			FP_Theta_Min = th_min;
			FP_Theta_Max = th_max;
			FP_Phi_Min = ph_min;
			FP_Phi_Max = ph_max;
			FP_X_Min = x_min;
			FP_X_Max = x_max;
			FP_Y_Min = x_min;
			FP_Y_Max = x_max;
		}

	private:
		void Init_HRS(){
			//Define the HRS acceptance range on the target plane,
			//Note that the actual acceptance should be smaller
			HRS_TG_Dp_Min = -0.045; //-%
			HRS_TG_Dp_Max =  0.045; //+%
			HRS_TG_Phi_Min = -0.032;//rad
			HRS_TG_Phi_Max =  0.032;//rad
			HRS_TG_Theta_Min = -0.045;//rad
			HRS_TG_Theta_Max =  0.045;//rad
			HRS_TG_X_Min = -0.0; //
			HRS_TG_X_Max =  0.0; //
			HRS_TG_Y_Min = -0.0; //
			HRS_TG_Y_Max =  0.0; //

			//Define the HRS acceptance range on the focal plane:
			HRS_FP_Phi_Min = -0.09;//rad
			HRS_FP_Phi_Max =  0.09;//rad
			HRS_FP_Theta_Min = -0.50;//rad
			HRS_FP_Theta_Max =  0.50;//rad
			HRS_FP_X_Min = -1.0;//cm
			HRS_FP_X_Max =  1.0;//cm
			HRS_FP_Y_Min = -0.95;//cm
			HRS_FP_Y_Max =  0.95;//cm
		}
		void Init_BB(){
			//Define the  BigBite acceptance range on the target plane,
			//Note that the actual acceptance should be smaller
			BB_TG_Dp_Min =-0.045; //%  
			BB_TG_Dp_Max = 0.045; //% 
			BB_TG_Phi_Min = -0.08;//rad
			BB_TG_Phi_Max =  0.08;//rad
			BB_TG_Theta_Min = -0.17;//rad
			BB_TG_Theta_Max =  0.17;//rad
			BB_TG_X_Min = -0.0; //
			BB_TG_X_Max =  0.0; //
			BB_TG_Y_Min = -0.0; //
			BB_TG_Y_Max =  0.0; //

			//Define the BigBite acceptance range on the focal plane:
			BB_FP_Phi_Min = -0.09;//rad
			BB_FP_Phi_Max =  0.09;//rad
			BB_FP_Theta_Min = -0.50;//rad
			BB_FP_Theta_Max =  0.50;//rad
			BB_FP_X_Min = -1.0;//cm
			BB_FP_X_Max =  1.0;//cm
			BB_FP_Y_Min = -0.95;//cm
			BB_FP_Y_Max =  0.95;//cm
		}
    private:
        double UnifRand(){
            return rand()/double(RAND_MAX);
        }
        double UnifRand(double a, double b){
            return (b-a)*UnifRand()+a;
        }
        void seed(){
            srand((unsigned)time(0)); 
        }


	private:
		//TRandom2* fRan;

		double TG_Dp_Min ; //-6%
		double TG_Dp_Max ; //+6%
		double TG_Phi_Min ;//rad
		double TG_Phi_Max ;//rad
		double TG_Theta_Min ;//rad
		double TG_Theta_Max ;//rad
		double TG_X_Min ; //
		double TG_X_Max ; //
		double TG_Y_Min ; //
		double TG_Y_Max ; //

		double FP_Phi_Min ;//rad
		double FP_Phi_Max ;//rad
		double FP_Theta_Min ;//rad
		double FP_Theta_Max ;//rad
		double FP_X_Min ;//cm
		double FP_X_Max ;//cm
		double FP_Y_Min ;//cm
		double FP_Y_Max ;//cm

		//Define the Vertex_Z Range:
		double Vertex_Z_Min ; //cm
		double Vertex_Z_Max; //cm

	private:
			//Define the HRS acceptance range on the target plane,
			//Note that the actual acceptance should be smaller
			double HRS_TG_Dp_Min ; //-6%
			double HRS_TG_Dp_Max ; //+6%
			double HRS_TG_Phi_Min ;//rad
			double HRS_TG_Phi_Max ;//rad
			double HRS_TG_Theta_Min ;//rad
			double HRS_TG_Theta_Max ;//rad
			double HRS_TG_X_Min ; //
			double HRS_TG_X_Max ; //
			double HRS_TG_Y_Min ; //
			double HRS_TG_Y_Max ; //

			//Define the HRS acceptance range on the focal plane:
			double HRS_FP_Phi_Min ;//rad
			double HRS_FP_Phi_Max ;//rad
			double HRS_FP_Theta_Min ;//rad
			double HRS_FP_Theta_Max ;//rad
			double HRS_FP_X_Min ;//cm
			double HRS_FP_X_Max ;//cm
			double HRS_FP_Y_Min ;//cm
			double HRS_FP_Y_Max ;//cm

			//Define the  BigBite acceptance range on the target plane,
			//Note that the actual acceptance should be smaller
			double BB_TG_Dp_Min ; //-6%
			double BB_TG_Dp_Max ; //+6%
			double BB_TG_Phi_Min;//rad
			double BB_TG_Phi_Max ;//rad
			double BB_TG_Theta_Min ;//rad
			double BB_TG_Theta_Max ;//rad
			double BB_TG_X_Min ; //
			double BB_TG_X_Max ; //
			double BB_TG_Y_Min ; //
			double BB_TG_Y_Max ; //

			//Define the HRS acceptance range on the focal plane:
			double BB_FP_Phi_Min ;//rad
			double BB_FP_Phi_Max ;//rad
			double BB_FP_Theta_Min ;//rad
			double BB_FP_Theta_Max ;//rad
			double BB_FP_X_Min ;//cm
			double BB_FP_X_Max ;//cm
			double BB_FP_Y_Min ;//cm
			double BB_FP_Y_Max ;//cm
};
#endif
