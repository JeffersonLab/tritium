//////////////////////////////////////////////////////
// A Subroutine to call XEMC Cross Section Model,
// or get Cross Section values from a lookup table
//     -- Zhihong Ye, 07/11/2012
//////////////////////////////////////////////////////
//#include "XEMC.h"
//Make sure the bins are the same as ones when generating tables
//Check with the code that generates the tables "./Gen_XS_Table/GenTable.C"
//const int fTheta_Bin_Num = 200;
//const int fEp_Bin_Num = 200;   //For old tables by Kin
//const int fEp_Bin_Num = 1000; //For new tables
//Define Acceptace Range and Bin-Size
//const double Ep_Bin_Size = 0.001;// 1 MeV/c per Bin
//const Double_t Delta_Theta = 0.07*TMath::RadToDeg()*2;     
//const Double_t Size_Theta = Delta_Theta/fTheta_Bin_Num;
const double Ep_Min_21 = 2.730, Ep_Max_21 = 3.240;
const double Ep_Min_23 = 2.680, Ep_Max_23 = 3.220;
const double Ep_Min_25 = 2.350, Ep_Max_25 = 3.180;
const double Ep_Min_28 = 2.680, Ep_Max_28 = 3.020;
class Get_XS
{
	public:
		Get_XS(){
			//do nothing
            fE_Beam = 3.356;//GeV
		}
		virtual ~Get_XS(){
			// delete aXS;
		}

		/* void Load_Table( const TString& kTarget_Name,const TString& kKin){{{*/
		void Load_Table( const TString& kTarget_Name,const TString& kKin){
			//The path or file name might be changed, so double check before use.

			/*Set Kin{{{*/
			Int_t aScat_Angle=0; 
			if(kKin.Contains("3.1")||kKin.Contains("3.2"))  {
				fTheta_Center = 21.00; aScat_Angle=21; fEp_Min = Ep_Min_21; fEp_Max = Ep_Max_21; } 
			else if(kKin.Contains("4.1")||kKin.Contains("4.2"))  {
				fTheta_Center = 23.00; aScat_Angle=23; fEp_Min = Ep_Min_23; fEp_Max = Ep_Max_23; } 
			else if(kKin.Contains("5.0")||kKin.Contains("5.05")||kKin.Contains("5.1")||kKin.Contains("5.2"))  {
				fTheta_Center = 25.00; aScat_Angle=25; fEp_Min = Ep_Min_25; fEp_Max = Ep_Max_25; } 
			else if(kKin.Contains("6.5"))  {
				fTheta_Center = 28.00; aScat_Angle=28; fEp_Min = Ep_Min_28; fEp_Max = Ep_Max_28; } 
			else {
				cerr<< "I don't understand the Kinematics setting!"<<endl;}
			/*}}}*/
			TString Table_Name = Form("XS_%s_%d.table", kTarget_Name.Data(),aScat_Angle); 
			TString Table_File = XEMC_Table_Path+Table_Name;
			cerr<<" --- Reading XS tables from file = "<<Table_File.Data()<<endl;
			ifstream table; table.open(Table_File.Data());
			if(!table.good()){
				outlog<<"*** ERROR, cannot read target table from"<<Table_File.Data()<<endl;
				cerr<<"*** ERROR, cannot read target table from"<<Table_File.Data()<<endl;
			}
			//variables to stores values reading from tables
			TString comment;
			comment.ReadLine(table);
			Double_t aAngle,aEp, aTheta, aXS_Born, aXS_Rad,aTail;
			Int_t aA,aZ, aEp_Bin, aTheta_Bin;
			Int_t aCount=0;
			//Now load the entire table into the 2D Array
            Int_t aEp_Bin_Max = 0, aTheta_Bin_Max = 0;
			while(!table.eof()){	
				table >> aA >> aZ >> aAngle >> aEp >> aEp_Bin >> aTheta >> aTheta_Bin >> aXS_Born >> aXS_Rad >> aTail;

                //===================================================================================================
				//The following section is designed for studying the bin-centering correction, Z. Ye 02/12/2015
				//MAKE-SURE to turn off this section when calculating the normal cross section values
				if(0){//if(0) to turn off
					double aE0 = fE_Beam; //Beam energy fixed at 3.356 GeV
					double aQ2 = aE0 * aEp * pow( sin(aTheta*TMath::Pi()/180./2.) ,2);
					double aXbj = aQ2 / (2.*TARGET_MASS_A*(aE0-aEp));

					aXS_Born *= aXbj; aXS_Rad *= aXbj; //method#1: scale the model by aXbj
			    	aXS_Born *= sqrt(aQ2)/2; aXS_Rad *= sqrt(aQ2)/2; //method#2: scale the model by sqrt(Q2)/2
				   cerr<<"*************************************************************************************"<<endl;
				   cerr<<"*** Attention, I applied a fake correction on the XEMC model. Please remove it!!! ***"<<endl;
				   cerr<<"*************************************************************************************"<<endl;
				   
				   outlog<<"*************************************************************************************"<<endl;
				   outlog<<"*** Attention, I applied a fake correction on the XEMC model. Please remove it!!! ***"<<endl;
				   outlog<<"*************************************************************************************"<<endl;
				}
				//===================================================================================================

				aEp_Table[aTheta_Bin][aEp_Bin] = aEp; 	
				aTheta_Table[aTheta_Bin][aEp_Bin] = aTheta; 	
				aXS_Born_Table[aTheta_Bin][aEp_Bin] = aXS_Born; 	
				aXS_Rad_Table[aTheta_Bin][aEp_Bin] = aXS_Rad; 	
				//aXS_Rad_Table[aTheta_Bin][aEp_Bin] = aXS_Rad-aTail; 	//Test the tail, Z. Ye 01/13/2015

				//cerr<<Form("Ep=%f,Theta=%f, XS_Born[Theta=%d][Ep=%d]=%f",aEp,aTheta,aTheta_Bin,aEp_Bin,aXS_Born)<<endl;
                if(aEp_Bin > aEp_Bin_Max) aEp_Bin_Max = aEp_Bin;
                if(aTheta_Bin > aTheta_Bin_Max) aTheta_Bin_Max = aTheta_Bin;
				aCount++;
			}
			table.close();

            fEp_Min = aEp_Table[0][0];
            fEp_Max = aEp_Table[0][aEp_Bin_Max];
            fEp_Bin_Num = aEp_Bin_Max; 
            fEp_Bin_Size = (fEp_Max - fEp_Min)/fEp_Bin_Num;

            fTheta_Min = aTheta_Table[0][0];
            fTheta_Max = aTheta_Table[aTheta_Bin_Max][0];
            fTheta_Bin_Num = aTheta_Bin_Max; 
            fTheta_Bin_Size = (fTheta_Max - fTheta_Min)/fTheta_Bin_Num;

            cout<<Form("Ep:     Min=%6.4f, Max=%6.5f, Bin=%d, Size=%6.5f", fEp_Min, fEp_Max, fEp_Bin_Num,fEp_Bin_Size)<<endl;
            cout<<Form("Theta:  Min=%6.4f, Max=%6.5f, Bin=%d, Size=%6.5f", fTheta_Min, fTheta_Max, fTheta_Bin_Num,fTheta_Bin_Size)<<endl;
			cerr<<" --- Done!"<<endl;
		}
		/*}}}*/
		
		/*inline double gGet_XS( const double kE0, const double kEp, const double kTheta){{{*/
		inline double gGet_XS( const double kE0, const double kEp, const double kTheta){

			//The following info are from XS_Table generating subroutine,make sure values are the same as ones in that subroutine
			//Decide which element this setting mostly close to
            int iT = (int) ((fabs(kTheta) - fTheta_Min)/fTheta_Bin_Size+0.5);//find the closest point
			int iP = (int) ((kEp - fEp_Min)/fEp_Bin_Size)+1;//find two points that this value falls in between
 

//			cerr<<Form("  For Theta=%f, Ep=%f, Bin_Theta=%d, Bin_Ep=%d", kTheta, kEp, iT,iP)<<endl;
//
			Double_t aDelta_Ep = (aEp_Table[iT][iP-1]-aEp_Table[iT][iP]);
			Double_t aDelta_XS = (aXS_Rad_Table[iT][iP-1]-aXS_Rad_Table[iT][iP]);
			Double_t aXS_Rad = 0.0;

            if(kEp<aEp_Table[iT][iP-1] || kEp>aEp_Table[iT][iP]){
                cerr<<Form("--- ERROR, Finding Ep incorrectly when reading the XS tables (XGT2_XEMC.h, #line162, Ep=%6.4f~[%d, %6.4f,%6.4f]) ", kEp, iP, aEp_Table[iT][iP-1],  aEp_Table[iT][iP])<<endl;
                outlog<<Form("--- ERROR, Finding Ep incorrectly when reading the XS tables (XGT2_XEMC.h, #line162, Ep=%6.4f~[%d, %6.4f,%6.4f]) ", kEp, iP, aEp_Table[iT][iP-1],aEp_Table[iT][iP])<<endl;
                //return -1000;
                exit(-100);
            }
            if(abs(kTheta-aTheta_Table[iT][iP])>fTheta_Bin_Size*0.5){
                cerr<<Form("--- ERROR, Finding Angle incorrectly when reading the XS tables (XGT2_XEMC.h, #line168, T=%6.4f~[%d, %6.4f,%6.4f,%6.4f]) ", kTheta, iT, aTheta_Table[iT-1][iP], aTheta_Table[iT][iP], aTheta_Table[iT+1][iP])<<endl;
                //outlog<<Form("--- ERROR, Finding Angle incorrectly when reading the XS tables (XGT2_XEMC.h, #line168, T=%6.4f~[%d, %6.4f,%6.4f,%6.4f]) ", kTheta, iT, aTheta_Table[iT-1][iP], aTheta_Table[iT][iP], aTheta_Table[iT+1][iP])<<endl;
                return -2000;
                exit(-110);
            }

			if(fabs(aDelta_Ep)>1e-16&&aXS_Rad_Table[iT][iP]>1e-16&&aXS_Rad_Table[iT][iP-1]>1e-16)
				aXS_Rad  = aXS_Rad_Table[iT][iP-1] - (aEp_Table[iT][iP-1]-kEp)/aDelta_Ep*aDelta_XS; 
            else
				aXS_Rad = 0.0;

			return aXS_Rad;
		}
		/*}}}*/

		/*inline double gGet_XS_Born( const double kE0, const double kEp, const double kTheta){{{*/
		inline double gGet_XS_Born( const double kE0, const double kEp, const double kTheta){

			//The following info are from XS_Table generating subroutine,make sure values are the same as ones in that subroutine
			//Decide which element this setting mostly close to
           	int iT = (int) ((fabs(kTheta) - fTheta_Min)/fTheta_Bin_Size+0.5);//find the closest point
			int iP = (int) ((kEp - fEp_Min)/fEp_Bin_Size)+1;//find two points that this value falls in between
 
            Double_t aDelta_Ep = (aEp_Table[iT][iP-1]-aEp_Table[iT][iP]);
            Double_t aDelta_XS = (aXS_Born_Table[iT][iP-1]-aXS_Born_Table[iT][iP]);
			Double_t aXS_Born = 0.0;

		//	outlog<< Form("-----iT=%d, iP=%d, Ep=%f, Th=%f XS_Born=%e, XS_Born=%e ----", 
		//			iT, iP-1,  aEp_Table[iT][iP-1],aTheta_Table[iT][iP-1],aXS_Born_Table[iT][iP-1],aXS_Born_Table[iT][iP])<<endl;

		//	cerr<< Form("-----iT=%d, iP=%d, Ep=%f, Th=%f XS_Born=%e, XS_Born=%e ----", 
		//			iT, iP-1,  aEp_Table[iT][iP-1],aTheta_Table[iT][iP-1],aXS_Born_Table[iT][iP-1],aXS_Born_Table[iT][iP])<<endl;

            if(kEp<aEp_Table[iT][iP-1] || kEp>aEp_Table[iT][iP]){
                cerr<<Form("--- ERROR, Finding Ep incorrectly when reading the XS tables (XGT2_XEMC.h, #line162, Ep=%6.4f~[%d, %6.4f,%6.4f]) ", kEp, iP, aEp_Table[iT][iP-1],  aEp_Table[iT][iP])<<endl;
                outlog<<Form("--- ERROR, Finding Ep incorrectly when reading the XS tables (XGT2_XEMC.h, #line162, Ep=%6.4f~[%d, %6.4f,%6.4f]) ", kEp, iP, aEp_Table[iT][iP-1],aEp_Table[iT][iP])<<endl;
                //return -1000;
                exit(-100);
            }
            if(abs(kTheta-aTheta_Table[iT][iP])>fTheta_Bin_Size*0.5){
                cerr<<Form("--- ERROR, Finding Angle incorrectly when reading the XS tables (XGT2_XEMC.h, #line168, T=%6.4f~[%d, %6.4f,%6.4f,%6.4f]) ", kTheta, iT, aTheta_Table[iT-1][iP], aTheta_Table[iT][iP], aTheta_Table[iT+1][iP])<<endl;
                outlog<<Form("--- ERROR, Finding Angle incorrectly when reading the XS tables (XGT2_XEMC.h, #line168, T=%6.4f~[%d, %6.4f,%6.4f,%6.4f]) ", kTheta, iT, aTheta_Table[iT-1][iP], aTheta_Table[iT][iP], aTheta_Table[iT+1][iP])<<endl;
                //return -2000;
                exit(-110);
            }

			if(fabs(aDelta_Ep)>1e-16&&aXS_Born_Table[iT][iP]>1e-16&&aXS_Born_Table[iT][iP-1]>1e-16)
			     aXS_Born  = aXS_Born_Table[iT][iP-1] - (aEp_Table[iT][iP-1]-kEp)/aDelta_Ep*aDelta_XS; 
            else
				aXS_Born = 0.0;

			return aXS_Born;
		}
		/*}}}*/

		inline double gGet_XS_Center(){
			double aXS = gGet_XS(fE_Beam,fEp_Center,fTheta_Center);
			return aXS;
		}

		inline double gGet_XS_Center(double kEp){
			double aXS = gGet_XS(fE_Beam,kEp,fTheta_Center);
			return aXS;
		}

	private:
		//  XEMC* aXS;
		int aA_Table[1000][1000];
		int aZ_Table[1000][1000]; 
		double aE0_Table[1000][1000];
		double aEp_Table[1000][1000];
		double aEp_1000able[1000][1000];
		double aTheta_Table[1000][1000];
		double aTheta_Bin_Table[1000][1000];
		double aXS_Rad_Table[1000][1000];
		double aXS_Born_Table[1000][1000];
        double fE_Beam;
		double fEp_Center;
        double fEp_Min;
		double fEp_Max;
		double fEp_Bin_Size;
		int fEp_Bin_Num;

		double fTheta_Center;
        double fTheta_Min;
		double fTheta_Max;
		double fTheta_Bin_Size;
        int fTheta_Bin_Num;
};
