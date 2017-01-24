////////////////////////////////////////////////////////////
//    A class to call and calculate cross section from    //
//      XEMC Radiation Cross Section Model                //
//        -- Zhihong Ye, 06/18/2012                       //
////////////////////////////////////////////////////////////
#include "XEMCEvent.h"

class XEMC
{
	public:
		XEMC(){
		}

		/*~XEMC{{{*/
		virtual ~XEMC(){
			Win_Before_Mag_Name.clear();
			Win_Before_Mag_Z.clear();
			Win_Before_Mag_A.clear();
			Win_Before_Mag_L.clear();
			Win_Before_Mag_rho.clear();
			Win_Before_Mag_X0.clear();

			Win_After_Mag_Name.clear();
			Win_After_Mag_Z.clear();
			Win_After_Mag_A.clear();
			Win_After_Mag_L.clear();
			Win_After_Mag_rho.clear();
			Win_After_Mag_X0.clear();

			delete Event;
		}
        /*}}}*/

		XEMC(XEMC const&){};
		XEMC& operator=(XEMC const&){};

		/*int Process(const double kE0, const double kEp, const double kTheta, const int kA, const int kZ){{{*/
		int Process(const double kE0, const double kEp, const double kTheta, const int kA, const int kZ)
		{ 
			if(kE0>12.0||kEp>12.0){
				cerr<<" *** Are you sure? E0 or Ep > 12GeV?!!! *** "<<endl;
				return -1;
			}
			/*{{{ Intialize parameters and calculate Radiated Cross Section*/
			Event->Init(); //Have to initilize every time to clear Win_Before/After_Mag
			Event->E_s=kE0*1000.;   //E0, in MeV
			Event->E_p=kEp*1000.;   //Ep, in MeV
			Event->theta=kTheta; //Theta, in Deg
			Event->Target.A = kA;
			Event->Target.Z = kZ;
			//And re-add materials for Win_Before/After_Mag with new angle setting
			for(unsigned int ii=0;ii<Win_Before_Mag_Name.size();ii++)
				Event->AddOneMaterial(Event->Win_Before_Mag, Win_Before_Mag_X0[ii], Win_Before_Mag_rho[ii],Win_Before_Mag_L[ii],Win_Before_Mag_A[ii],Win_Before_Mag_Z[ii],Win_Before_Mag_Name[ii].c_str());
			for(unsigned int jj=0;jj<Win_After_Mag_Name.size();jj++)
				Event->AddOneMaterial(Event->Win_After_Mag, Win_After_Mag_X0[jj], Win_After_Mag_rho[jj],Win_After_Mag_L[jj],Win_After_Mag_A[jj],Win_After_Mag_Z[jj],Win_After_Mag_Name[jj].c_str());
			Event->reactz_gen = 0 - z0; //For foil targets, assume the reaction points are alwayl at the center
			err=Event->Run();
			if ( err==-1 ){
				cerr << "********** Don't not success to calculate XS!!!"<<endl;
				xs_rad = -1111.00;
				xs_dis = -1112.00;
				xs_qe = -1113.00;
				xs_born = -2222.00;
				tail = -1000.00;
			}
			else{
				xs_born = Event->cs_Born;
				xs_dis = Event->cs_dis;
				xs_qe = Event->cs_qe;
				tail = Event->tail;
			}
			xs_mott = Event->cs_M;

			if ( IsDebug && err==0 )
				Event->Print();

			/*Finished Radiated Cross Section Calculation}}}*/
			return 0;
		}
		/*}}}*/

		/*int Process(const double kE0, const double kEp, const double kTheta, const int kA, const int kZ, const kRctPt){{{*/
		int Process(const double kE0, const double kEp, const double kTheta, const int kA, const int kZ, const double kRctPt)
		{ 
			if(kE0>12.0||kEp>12.0){
				cerr<<" *** Are you sure? E0 or Ep > 12GeV?!!! *** "<<endl;
				return -1;
			}
			/*{{{ Intialize parameters and calculate Radiated Cross Section*/
			Event->Init(); //Have to initilize every time to clear Win_Before/After_Mag
			Event->E_s=kE0*1000.;   //E0, in MeV
			Event->E_p=kEp*1000.;   //Ep, in MeV
			Event->theta=kTheta; //Theta, in Deg
			Event->Target.A = kA;
			Event->Target.Z = kZ;
			//And re-add materials for Win_Before/After_Mag with new angle setting
			for(unsigned int ii=0;ii<Win_Before_Mag_Name.size();ii++)
				Event->AddOneMaterial(Event->Win_Before_Mag, Win_Before_Mag_X0[ii], Win_Before_Mag_rho[ii],Win_Before_Mag_L[ii],Win_Before_Mag_A[ii],Win_Before_Mag_Z[ii],Win_Before_Mag_Name[ii].c_str());
			for(unsigned int jj=0;jj<Win_After_Mag_Name.size();jj++)
				Event->AddOneMaterial(Event->Win_After_Mag, Win_After_Mag_X0[jj], Win_After_Mag_rho[jj],Win_After_Mag_L[jj],Win_After_Mag_A[jj],Win_After_Mag_Z[jj],Win_After_Mag_Name[jj].c_str());
			Event->reactz_gen = kRctPt - z0; //For foil targets, assume the reaction points are alwayl at the center
			err=Event->Run();
			if ( err==-1 ){
				cerr << "********** Don't not success to calculate XS!!!"<<endl;
				xs_rad = -1111.00;
				xs_born = -2222.00;
				xs_dis = -1112.00;
				xs_qe = -1113.00;
				tail = -1000.00;
			}
			else{
				xs_dis = Event->cs_dis;
				xs_qe = Event->cs_qe;
				xs_born = Event->cs_Born;
				xs_rad = Event->cs_Final;
				tail = Event->tail;
			}
			xs_mott = Event->cs_M;

			if ( IsDebug && err==0 )
				Event->Print();

			/*Finished Radiated Cross Section Calculation}}}*/
			return 0;
		}
		/*}}}*/

		/*int Process(const double kE0, const double kEp, const double kTheta){{{*/
		int Process(const double kE0, const double kEp, const double kTheta)
		{
			/*{{{ Intialize kinematic parameters and calculate Radiated Cross Section*/
			if(kE0>12.0||kEp>12.0){
				cerr<<" *** Are you sure? E0 or Ep > 12GeV?!!! *** "<<endl;
				return -1;
			}
			Event->Init(); //Have to initilize every time to clear Win_Before/After_Mag
			Event->E_s=kE0*1000.;   //E0, in MeV
			Event->E_p=kEp*1000.;   //Ep, in MeV
			Event->theta=kTheta; //Theta, in Deg

			//And re-add materials for Win_Before/After_Mag with new angle setting
			for(unsigned int ii=0;ii<Win_Before_Mag_Name.size();ii++)
				Event->AddOneMaterial(Event->Win_Before_Mag, Win_Before_Mag_X0[ii], Win_Before_Mag_rho[ii],Win_Before_Mag_L[ii],Win_Before_Mag_A[ii],Win_Before_Mag_Z[ii],Win_Before_Mag_Name[ii].c_str());
			for(unsigned int jj=0;jj<Win_After_Mag_Name.size();jj++)
				Event->AddOneMaterial(Event->Win_After_Mag, Win_After_Mag_X0[jj], Win_After_Mag_rho[jj],Win_After_Mag_L[jj],Win_After_Mag_A[jj],Win_After_Mag_Z[jj],Win_After_Mag_Name[jj].c_str());

			Event->reactz_gen = 0 - z0; //For foil targets, assume the reaction points are alwayl at the center
			err=Event->Run();
			if ( err==-1 ){
				cerr << "********** Don't not success to calculate XS!!!"<<endl;
				xs_rad = -1111.00;
				xs_dis = -1111.00;
				xs_qe = -1111.00;
				xs_born = -2222.00;
				tail = -1000.00;
			}
			else{
				xs_born = Event->cs_Born;
				xs_dis = Event->cs_dis;
				xs_qe = Event->cs_qe;
				xs_rad = Event->cs_Final;
				tail = Event->tail;
			}
			xs_mott = Event->cs_M;

			if ( IsDebug && err==0 )
				Event->Print();

			/*Finished Radiated Cross Section Calculation}}}*/
			return 0;
		}
		/*}}}*/

		/*int Process(const double kE0, const double kEp, const double kTheta, const double kRctPt){{{*/
		int Process(const double kE0, const double kEp, const double kTheta, const double kRctPt)
		{
			/*{{{ Intialize kinematic parameters and calculate Radiated Cross Section*/
			if(kE0>12.0||kEp>12.0){
				cerr<<" *** Are you sure? E0 or Ep > 12GeV?!!! *** "<<endl;
				return -1;
			}
			Event->Init(); //Have to initilize every time to clear Win_Before/After_Mag
			Event->E_s=kE0*1000.;   //E0, in MeV
			Event->E_p=kEp*1000.;   //Ep, in MeV
			Event->theta=kTheta; //Theta, in Deg
				tail = Event->tail;

			//And re-add materials for Win_Before/After_Mag with new angle setting
			for(unsigned int ii=0;ii<Win_Before_Mag_Name.size();ii++)
				Event->AddOneMaterial(Event->Win_Before_Mag, Win_Before_Mag_X0[ii], Win_Before_Mag_rho[ii],Win_Before_Mag_L[ii],Win_Before_Mag_A[ii],Win_Before_Mag_Z[ii],Win_Before_Mag_Name[ii].c_str());
			for(unsigned int jj=0;jj<Win_After_Mag_Name.size();jj++)
				Event->AddOneMaterial(Event->Win_After_Mag, Win_After_Mag_X0[jj], Win_After_Mag_rho[jj],Win_After_Mag_L[jj],Win_After_Mag_A[jj],Win_After_Mag_Z[jj],Win_After_Mag_Name[jj].c_str());

			Event->reactz_gen = kRctPt - z0; //For foil targets, assume the reaction points are alwayl at the center
			err=Event->Run();
			if ( err==-1 ){
				cerr << "********** Don't not success to calculate XS!!!"<<endl;
				xs_rad = -1111.00;
				xs_born = -2222.00;
				xs_dis = -1111.00;
				xs_qe = -1111.00;
				tail = -1000.00;
			}
			else{
				xs_born = Event->cs_Born;
				xs_dis = Event->cs_dis;
				xs_qe = Event->cs_qe;
				xs_rad = Event->cs_Final;
				tail = Event->tail;
			}
			xs_mott = Event->cs_M;

			if ( IsDebug && err==0 )
				Event->Print();

			/*Finished Radiated Cross Section Calculation}}}*/
			return 0;
		}
		/*}}}*/

		/*Return Variables{{{*/
		double XS_Rad(){
			return xs_rad;
		}

		double XS_Born(){
			return xs_born;
		}
		double XS_QE(){
			return xs_qe;
		}
		double XS_DIS(){
			return xs_dis;
		}
		double XS_Mott(){
			return xs_mott;
		}
		double Elastic_Tail(){
			return tail;
		}
        /*}}}*/

		/*Init(){{{*/
		void Init (const string kTarget_Input){

			cerr <<"---> Initializing Target Parameters from "<<kTarget_Input<<" ..."<<endl;
			int i,j,k;
			vector<string> inputdata;
			/*Read XEMCfile{{{*/
			FILE* XEMCfile;
			XEMCfile=fopen(kTarget_Input.c_str(),"r");
			char buf[CHAR_LEN];
			char data[CHAR_LEN];
			while ( fgets(buf,CHAR_LEN,XEMCfile) )
			{
				i=0;
				while ( buf[i]==' ' )
				{
					i++;
				}
				if ( buf[i]!='#' )
				{
					j=0;
					while ( buf[i]!='#' && buf[i]!='\0' )
					{
						data[j]=buf[i];
						i++; j++;
					}
					data[j]='\0';
					while ( data[--j]==' ' || data[j]=='\t' )
					{
						//remove space or tab at the end of data
						data[j]='\0';
					}
					inputdata.push_back(data);
				}
				//else it's comment, skipped
			}
			fclose(XEMCfile);
			/*}}}*/

			/*Set Global Value{{{*/
			k=0;
			IsDebug=atoi(inputdata[k++].c_str());
			IsEnergyLoss=atoi(inputdata[k++].c_str());
			HRS_L=atof(inputdata[k++].c_str());
			D_x=atof(inputdata[k++].c_str());
			D_y=atof(inputdata[k++].c_str());
			beam_x=atof(inputdata[k++].c_str());
			beam_y=atof(inputdata[k++].c_str());
			//QFS model, not used in XEMC
			IsQFS=atoi(inputdata[k++].c_str());
			IsQFS_Q2dep=atoi(inputdata[k++].c_str());
			Fermi_Moment=atof(inputdata[k++].c_str());
			NIE=atof(inputdata[k++].c_str());
			DEL_SEP=atof(inputdata[k++].c_str());
			Form_Factor_Id=atoi(inputdata[k++].c_str());
			IsInternal_Bre=atoi(inputdata[k++].c_str());
			IsAppro=atoi(inputdata[k++].c_str());
			IsExternal_Bre=atoi(inputdata[k++].c_str());
			IsMulti_Photon_Cor=atoi(inputdata[k++].c_str());
			IsRadCor=atoi(inputdata[k++].c_str());
			DeltaE=atof(inputdata[k++].c_str());
			XEMCFLAG=atoi(inputdata[k++].c_str());
			//Add Target
			Target_Name=inputdata[k++].c_str();
			Target_Z=atoi(inputdata[k++].c_str());
			Target_A=atof(inputdata[k++].c_str());
			Target_T=atof(inputdata[k++].c_str());
			Target_rho=atof(inputdata[k++].c_str());
			z0=atof(inputdata[k++].c_str());
			T_L=atof(inputdata[k++].c_str());
			T_H=atof(inputdata[k++].c_str());
			Target_L=T_L;//in HCS
			//Add Target Cell
			Win_i_Name=inputdata[k++].c_str();
			Win_i_Z=atoi(inputdata[k++].c_str());
			Win_i_A=atof(inputdata[k++].c_str());
			Win_i_T=atof(inputdata[k++].c_str());
			Win_i_rho=atof(inputdata[k++].c_str());
			Win_f_Name=inputdata[k++].c_str();
			Win_f_Z=atoi(inputdata[k++].c_str());
			Win_f_A=atof(inputdata[k++].c_str());
			Win_f_T=atof(inputdata[k++].c_str());
			Win_f_rho=atof(inputdata[k++].c_str());
			Target_Theta=atof(inputdata[k++].c_str());
			IsBump=atoi(inputdata[k++].c_str());

			//Target Chamber and other material before HRS magnetics
			j=k;
			while ( atof(inputdata[j].c_str())>=0 ){
				Win_Before_Mag_Name.push_back(inputdata[j++].c_str());
				Win_Before_Mag_Z.push_back(atoi(inputdata[j++].c_str()));
				Win_Before_Mag_A.push_back(atof(inputdata[j++].c_str()));
				Win_Before_Mag_L.push_back(atof(inputdata[j++].c_str()));
				Win_Before_Mag_rho.push_back(atof(inputdata[j++].c_str()));
				Win_Before_Mag_X0.push_back(atof(inputdata[j++].c_str()));
			}
			j++;
			while ( atof(inputdata[j].c_str())>=0 ){
				Win_After_Mag_Name.push_back(inputdata[j++].c_str());
				Win_After_Mag_Z.push_back(atoi(inputdata[j++].c_str()));
				Win_After_Mag_A.push_back(atof(inputdata[j++].c_str()));
				Win_After_Mag_L.push_back(atof(inputdata[j++].c_str()));
				Win_After_Mag_rho.push_back(atof(inputdata[j++].c_str()));
				Win_After_Mag_X0.push_back(atof(inputdata[j++].c_str()));     
			}

			/*Initialize values in XEMCEvent class{{{*/
			Event = new XEMCEvent();
			//Flags
			Event->IsDebug = IsDebug;
			Event->IsEnergyLoss = IsEnergyLoss;
			Event->HRS_L = HRS_L;
			Event->D_x = D_x;
			Event->D_y = D_y;
			Event->beam_x = beam_x;
			Event->beam_y = beam_y;
			Event->z0 = z0;
			Event->T_L = T_L;
			Event->T_H = T_H;
			//QFS model, not used for XEMC model
			Event->IsQFS = IsQFS;
			Event->IsQFS_Q2dep = IsQFS_Q2dep;
			Event->Fermi_Moment = Fermi_Moment;
			Event->NIE = NIE;
			Event->DEL_SEP = DEL_SEP;
			Event->Form_Factor_Id = Form_Factor_Id;
			Event->IsInternal_Bre = IsInternal_Bre;
			Event->IsAppro = IsAppro;
			Event->IsExternal_Bre = IsExternal_Bre;
			Event->IsMulti_Photon_Cor = IsMulti_Photon_Cor;
			Event->IsRadCor = IsRadCor;
			Event->DeltaE = DeltaE;
			Event->XEMCFLAG = XEMCFLAG;
			//Target
			Event->Target.Name= Target_Name;
			Event->Target.A= Target_A;
			Event->Target.Z= Target_Z;
			Event->Target.T= Target_T;
			Event->Target.rho= Target_rho;
			Event->Target.L=T_L;//in HCS
			Event->IsBump = IsBump;
			//Target Cells
			Event->Win_i.Name= Win_i_Name;
			Event->Win_i.Z= Win_i_Z;
			Event->Win_i.A= Win_i_A;
			Event->Win_i.T= Win_i_T;
			Event->Win_i.rho= Win_i_rho;
			Event->Win_f.Name= Win_f_Name;
			Event->Win_f.Z= Win_f_Z;
			Event->Win_f.A= Win_f_A;
			Event->Win_f.T= Win_f_T;
			Event->Win_f.rho= Win_f_rho;
			Event->T_theta= Target_Theta;
			/*}}}*/

			if ( IsDebug )
				DebugPrint();

			inputdata.clear();
			cerr<<"---> Initialization done!"<<endl;
		}
		/*}}}*/
		/*}}}*/

		/*Clear some variables if needed to re-initialize{{{*/
        void Clear(){
			Win_Before_Mag_Name.clear();
			Win_Before_Mag_Z.clear();
			Win_Before_Mag_A.clear();
			Win_Before_Mag_L.clear();
			Win_Before_Mag_rho.clear();
			Win_Before_Mag_X0.clear();

			Win_After_Mag_Name.clear();
			Win_After_Mag_Z.clear();
			Win_After_Mag_A.clear();
			Win_After_Mag_L.clear();
			Win_After_Mag_rho.clear();
			Win_After_Mag_X0.clear();
		}
        /*}}}*/

	private:
		/*DebugPrint(){{{*/
		void DebugPrint()
		{
			if ( IsEnergyLoss )
				printf("Energy Loss Enabled.\n");
			else
				printf("Energy Loss Disabled.\n");

			printf("%-*s=%*.2f %-*s %-*s\n",15,"HRS_L",            10,HRS_L,                    8,"cm",40,"(Length of HRS)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"D_x",              10,D_x,                      8,"cm",40,"(X offset(in TCS) between TCS and HCS)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"D_y",              10,D_y,                      8,"cm",40,"(Y offset(in TCS) between TCS and HCS)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"beam_x",           10,beam_x,                   8,"cm",40,"(beam x center for generator)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"beam_y",           10,beam_y,                   8,"cm",40,"(beam y center for generator)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"z0",               10,z0,                       8,"cm",40,"(target center for Target Windows)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"T_L",              10,T_L,                      8,"cm",40,"(target length for Target Windowns)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"T_H",              10,T_H,                      8,"cm",40,"(target height for Target Windowns)");
			printf("%-*s=%*.1d %-*s %-*s\n",15,"IsQFS " ,          10,IsQFS,                    8," ",40,"(QFS cross section model, default = false to use XEMC instead)");
			printf("%-*s=%*.1d %-*s %-*s\n",15,"IsQFS_Q2dep" ,     10,IsQFS_Q2dep,              8," ",40,"(QFS cross section model, not used in XEMC)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Fermi_Moment",     10,Fermi_Moment,             8,"MeV",40,"(fermi momentum of target for QFS model, not used in XEMC)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"NIE",              10,NIE,                      8,"MeV",40,"(nucleon separation energy for QFS model, not used in XEMC)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"DEL_SEP",          10,DEL_SEP,                  8,"MeV",40,"(delta separation energy for QFS model, not used in XEMC)");
			printf("%-*s=%*.2d %-*s %-*s\n",15,"Form_Factor",      10,Form_Factor_Id,           8," ",40,"(for QFS model, not used in XEMC)");
			printf("%-*s=%*.1d %-*s %-*s\n",15,"IsInternal_Bre" ,  10,IsInternal_Bre,           8," ",40,"(for QFS model, not used in XEMC)");
			printf("%-*s=%*.1d %-*s %-*s\n",15,"IsAppro",          10,IsAppro,                  8," ",40,"(for QFS model, not used in XEMC)");
			printf("%-*s=%*.1d %-*s %-*s\n",15,"IsExternal_Bre",   10,IsExternal_Bre,           8," ",40,"(for QFS model, not used in XEMC)");
			printf("%-*s=%*.1d %-*s %-*s\n",15,"Multi_Photon_Cor", 10,IsMulti_Photon_Cor,       8," ",40,"(for QFS model, not used in XEMC)");
			printf("%-*s=%*.1d %-*s %-*s\n",15,"IsRadC" ,          10,IsRadCor,                 8," ",40,"(enable raditive correction)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"DeltaE" ,          10,DeltaE,                   8," ",40,"(for QFS model, not used in XEMC)");
			printf("%-*s=%*.1d %-*s %-*s\n",15,"XEMCFLAG" ,        10,XEMCFLAG,                 8," ",40,"(for XEMC model, 1->QE+DIS,2->QE only, 3->DIS only)");

			//    printf("%-*s=%*.*s %-*s %-*s\n",15,"Target Name" ,       10,Target_Name,                   8," ",40,"");
			printf("%-*s=%*.2d %-*s %-*s\n",15,"Target Z" ,          10,Target_Z,               8," ",40,"");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Target A" ,          10,Target_A,               8," ",40,"");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Target T" ,          10,Target_T,               8," ",40,"");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Target rho" ,        10,Target_rho,             8," ",40,"");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Target L" ,          10,Target_L,               8," ",40,"");
			//   printf("%-*s=%*.2s %-*s %-*s\n",15,"Win_i Name" ,       10,Win_i_Name,                   8," ",40,"");
			printf("%-*s=%*.2d %-*s %-*s\n",15,"Win_i Z" ,          10,Win_i_Z,                 8," ",40,"");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Win_i A" ,          10,Win_i_A,                 8," ",40,"");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Win_i T" ,          10,Win_i_T,                 8," ",40,"");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Win_i rho" ,        10,Win_i_rho,               8," ",40,"");
			//   printf("%-*s=%*.2s %-*s %-*s\n",15,"Win_f Name" ,       10,Win_f_Name,               8," ",40,"");
			printf("%-*s=%*.2d %-*s %-*s\n",15,"Win_f Z" ,          10,Win_f_Z,                 8," ",40,"");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Win_f A" ,          10,Win_f_A,                 8," ",40,"");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Win_f T" ,          10,Win_f_T,                 8," ",40,"");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Win_f rho" ,        10,Win_f_rho,               8," ",40,"");

		}
		/*}}}*/

	private:
		/*Some Constants{{{*/		
		//some of them maybe define in other place, so make sure they have the same units.
		static const int CHAR_LEN=1024;
		//default unit
		//MeV cm g rad
		//cross section output unit = nbarn/sr/MeV

		int err;
		XEMCEvent* Event;
		double xs_rad, xs_born, xs_mott,xs_qe,xs_dis;

		bool IsDebug; //output to debugfile
		bool IsEnergyLoss; //enable Energy Loss
		double HRS_L; //Length Of HRS ( Left or Right,here _L=Length) (cm)
		double D_x; //X offset(in TCS) between TCS and HCS
		double D_y; //Y offset(in TCS) between TCS and HCS
		double beam_x; //cm
		double beam_y; //Beam Center cm
		double z0;    //target center for generator(cm)
		double T_L;   //target length for generator(cm)
		double T_H;   //target height for generator(cm)

		bool IsBump; //For special treatment of long targets due to the bump issue  
		bool IsQFS;    //QFS model flag, enable to calculate xs instead of XEMC.
		bool IsQFS_Q2dep; //enable Q2 dependence when calculating Q.E. Peak
		double Fermi_Moment;  //Fermi_Moment,fermi momentum of target (MeV)
		double NIE;   //NIE,nucleon interaction energy (MeV)
		double DEL_SEP;       //DEL_SEP,delta separation energy (MeV)

		bool IsRadCor; //enable Radiative correction
		double tail;
		int XEMCFLAG;
        /*}}}*/

		/*Only for Tail{{{*/
		int Form_Factor_Id; //which form factor used
		bool  IsInternal_Bre; //enable/disable internal bremsstrahlung
		bool  IsExternal_Bre; //enable/disable external bremsstrahlung
		bool  IsMulti_Photon_Cor; //enable/disable multiple photon correction
		bool  IsAppro; //use sigma_p instead of sigma_ex (approximate exact internal bremsstrahlung)
		double DeltaE; //default=5 MeV only for radiative correction

		string Target_Name;
		int Target_Z;
		double Target_A;
		double Target_T;
		double Target_rho;
		double Target_L;
		string Win_i_Name;
		int Win_i_Z;
		double Win_i_A;
		double Win_i_T;
		double Win_i_rho;
		string Win_f_Name;
		int Win_f_Z;
		double Win_f_A;
		double Win_f_T;
		double Win_f_rho;
		double Target_Theta;

		vector<string> Win_Before_Mag_Name;
		vector<int>    Win_Before_Mag_Z;
		vector<double> Win_Before_Mag_A;
		vector<double> Win_Before_Mag_L;
		vector<double> Win_Before_Mag_rho;
		vector<double> Win_Before_Mag_X0;

		vector<string> Win_After_Mag_Name;
		vector<int>    Win_After_Mag_Z;
		vector<double> Win_After_Mag_A;
		vector<double> Win_After_Mag_L;
		vector<double> Win_After_Mag_rho;
		vector<double> Win_After_Mag_X0;
		/*}}}*/

};
