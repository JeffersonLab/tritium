//////////////////////////////////////////////////////
// A Subroutine to call XEMC Cross Section Model,
// or get Cross Section values from a lookup table
//     -- Zhihong Ye, 07/11/2012
//////////////////////////////////////////////////////
//#include "XEMC.h"
//Make sure the bins are the same as ones when generating tables
const int Bin_T = 200;
const int Bin_P = 200;   //For old tables by Kin
//const int Bin_P = 1000; //For new tables
const double Ep_Bin_Size = 0.001;// 1 MeV/c per Bin
class Get_XS
{
	public:
		Get_XS(){
			//do nothing
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
				Theta_Center = 21.00; aScat_Angle=21; Ep_Min = 2.730; Ep_Max = 3.240; } 
			else if(kKin.Contains("4.1")||kKin.Contains("4.2"))  {
				Theta_Center = 23.00; aScat_Angle=23; Ep_Min = 2.680; Ep_Max = 3.220; } 
			else if(kKin.Contains("5.0")||kKin.Contains("5.05")||kKin.Contains("5.1")||kKin.Contains("5.2"))  {
				Theta_Center = 25.00; aScat_Angle=25; Ep_Min = 2.350; Ep_Max = 3.180; } 
			else if(kKin.Contains("6.5"))  {
				Theta_Center = 28.00; aScat_Angle=28; Ep_Min = 2.680; Ep_Max = 3.020; } 
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
			while(!table.eof()){	
				table >> aA >> aZ >> aAngle >> aEp >> aEp_Bin >> aTheta >> aTheta_Bin >> aXS_Born >> aXS_Rad >> aTail;

                //===================================================================================================
				//The following section is designed for studying the bin-centering correction, Z. Ye 02/12/2015
				//MAKE-SURE to turn off this section when calculating the normal cross section values
				if(0){//if(0) to turn off
					double aE0 = 3.356; //Beam energy fixed at 3.356 GeV
					double aQ2 = aE0 * aEp * pow( sin(aTheta*TMath::Pi()/180./2.) ,2);
					double aXbj = aQ2 / (2.*TARGET_MASS_A*(aE0-aEp));
					//double aXbj = aQ2 / (2.*PROTON_MASS*(aE0-aEp));

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
				aCount++;
			}
			table.close();

			cerr<<" --- Done!"<<endl;
		}
		/*}}}*/
		
		/* void Load_Table_Angle( const TString& kTarget_Name,const Double_t kAngle){{{*/
		void Load_Table_Angle( const TString& kTarget_Name,const Double_t kAngle){
			//The path or file name might be changed, so double check before use.

			/*Set Kin{{{*/
			Int_t aScat_Angle=0; 
			if(fabs(kAngle-21.)<1.0)  {
				Theta_Center = 21.00; aScat_Angle=21; Ep_Min = 2.730; Ep_Max = 3.240; } 
			else if(fabs(kAngle-23.)<1.0)  {
				Theta_Center = 23.00; aScat_Angle=23; Ep_Min = 2.680; Ep_Max = 3.220; } 
			else if(fabs(kAngle-25.)<1.0)  {
				Theta_Center = 25.00; aScat_Angle=25; Ep_Min = 2.350; Ep_Max = 3.180; } 
			else if(fabs(kAngle-28.)<1.0)  {
				Theta_Center = 28.00; aScat_Angle=28; Ep_Min = 2.680; Ep_Max = 3.020; } 
			else {
				cerr<< "I don't understand the Angle setting!"<<endl;}
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
			while(!table.eof()){	
				table >> aA >> aZ >> aAngle >> aEp >> aEp_Bin >> aTheta >> aTheta_Bin >> aXS_Born >> aXS_Rad >> aTail;

				aEp_Table[aTheta_Bin][aEp_Bin] = aEp; 	
				aTheta_Table[aTheta_Bin][aEp_Bin] = aTheta; 	
				aXS_Born_Table[aTheta_Bin][aEp_Bin] = aXS_Born; 	
				aXS_Rad_Table[aTheta_Bin][aEp_Bin] = aXS_Rad; 	
				//aXS_Rad_Table[aTheta_Bin][aEp_Bin] = aXS_Rad-aTail; 	//Test the tail, Z. Ye 01/13/2015

				//cerr<<Form("Ep=%f,Theta=%f, XS_Born[Theta=%d][Ep=%d]=%f",aEp,aTheta,aTheta_Bin,aEp_Bin,aXS_Born)<<endl;
				aCount++;
			}
			table.close();

			cerr<<" --- Done!"<<endl;
		}
		/*}}}*/

		/*inline double gGet_XS( const double kE0, const double kEp, const double kTheta){{{*/
		inline double gGet_XS( const double kE0, const double kEp, const double kTheta){

			//The following info are from XS_Table generating subroutine,make sure values are the same as ones in that subroutine
			//Define Acceptace Range and Bin-Size
			Double_t Delta_Theta = 0.07*TMath::RadToDeg();     Double_t Size_Theta = Delta_Theta/Bin_T;
			//Decide which element this setting mostly close to
			int iT = (int) ((fabs(kTheta) - Theta_Center + Delta_Theta/2.0)/Size_Theta+0.5);
			int iP = (int) ((kEp - Ep_Min)/Ep_Bin_Size+0.5);
//			cerr<<Form("  For Theta=%f, Ep=%f, Bin_Theta=%d, Bin_Ep=%d", kTheta, kEp, iT,iP)<<endl;
			Double_t aDelta_Ep = (aEp_Table[iT][iP-1]-aEp_Table[iT][iP]);
			Double_t aDelta_XS = (aXS_Rad_Table[iT][iP-1]-aXS_Rad_Table[iT][iP]);
			Double_t aXS_Rad = 0.0;

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
			//Define Acceptace Range and Bin-Size
			Double_t Delta_Theta = 0.07*TMath::RadToDeg();     Double_t Size_Theta = Delta_Theta/Bin_T;
			//Decide which element this setting mostly close to
			int iT = (int) ((fabs(kTheta) - Theta_Center + Delta_Theta/2.0)/Size_Theta+0.5);
			int iP = (int) ((kEp - Ep_Min)/Ep_Bin_Size+0.5);
            Double_t aDelta_Ep = (aEp_Table[iT][iP-1]-aEp_Table[iT][iP]);
			Double_t aDelta_XS = (aXS_Born_Table[iT][iP-1]-aXS_Born_Table[iT][iP]);
			Double_t aXS_Born = 0.0;

		//	outlog<< Form("-----iT=%d, iP=%d, Ep=%f, Th=%f XS_Born=%e, XS_Born=%e ----", 
		//			iT, iP-1,  aEp_Table[iT][iP-1],aTheta_Table[iT][iP-1],aXS_Born_Table[iT][iP-1],aXS_Born_Table[iT][iP])<<endl;

		//	cerr<< Form("-----iT=%d, iP=%d, Ep=%f, Th=%f XS_Born=%e, XS_Born=%e ----", 
		//			iT, iP-1,  aEp_Table[iT][iP-1],aTheta_Table[iT][iP-1],aXS_Born_Table[iT][iP-1],aXS_Born_Table[iT][iP])<<endl;

			if(fabs(aDelta_Ep)>1e-16&&aXS_Born_Table[iT][iP]>1e-16&&aXS_Born_Table[iT][iP-1]>1e-16)
			     aXS_Born  = aXS_Born_Table[iT][iP-1] - (aEp_Table[iT][iP-1]-kEp)/aDelta_Ep*aDelta_XS; 
            else
				aXS_Born = 0.0;

			return aXS_Born;
		}
		/*}}}*/

		inline double gGet_XS_Center(){
			double aXS = gGet_XS(E_Beam,Ep_Center,Theta_Center);
			return aXS;
		}

		inline double gGet_XS_Center(double kEp){
			double aXS = gGet_XS(E_Beam,kEp,Theta_Center);
			return aXS;
		}

	private:
		//  XEMC* aXS;
		int aA_Table[Bin_T][1000];
		int aZ_Table[Bin_T][1000]; 
		double aE0_Table[Bin_T][1000];
		double aEp_Table[Bin_T][1000];
		double aEp_Bin_Table[Bin_T][1000];
		double aTheta_Table[Bin_T][1000];
		double aTheta_Bin_Table[Bin_T][1000];
		double aXS_Rad_Table[Bin_T][1000];
		double aXS_Born_Table[Bin_T][1000];
		double Theta_Center;
		double Ep_Center;
		double Ep_Min;
		double Ep_Max;
		double E_Beam;
};
