#ifndef SAMC_PHYSICAL_CONSTANT_H
#include "SAMCEvent.h"
#endif

/*void usage(const char* aCommand){{{*/
void usage(const char* aCommand)
{
	printf("--------------------\n");
	printf("%-*s%s -dh filename\n",10,"Usage:",aCommand);
	printf("%-*s -d debug\n",10," ");
	printf("%-*s -h help\n",10," ");
	printf("%-*s filename: absolute file location\n",10," ");
	exit(8);
}
/*}}}*/

/*void getargs(int argc,char** argv){{{*/
void getargs(int argc,char** argv)
{
	const char* cmd=*argv;
	while ( argc-- > 1 )
	{
		const char* opt=*++argv;
		if ( *opt == '-' )
		{
			//option
			while ( *++opt != '\0' )
			{
				switch ( *opt )
				{
					case 'd':
						IsDebug=true;
						break;
					case 'h':
					default:
						usage(cmd);
				}
			}
		}
		else
		{
			if ( File_Name.empty() )
				File_Name=opt;
		}
	}
	if ( File_Name.empty() )
	{
		printf("[Error %s: Line %d] No File Input.\n",__FILE__,__LINE__);
		usage(cmd);
	}
}
/*}}}*/

/*bool IsANumber(string aString) {{{*/
bool IsANumber(string aString) 
{
	istringstream iss(aString);
	ostringstream oss;
	int x;
	iss >> x;
	oss << x;
	if ( aString == oss.str() )
		return true;
	else
		return false;
}
/*}}}*/

/*int main(int argc, char** argv){{{*/
int main(int argc, char** argv)
{
	srand(time(NULL));
	getargs(argc,argv);
	int i,j,k,fail_events=0;

	printf("Start...\n");
	vector<string> inputdata;
	/*Read samcfile{{{*/
	FILE* samcfile;
	samcfile=fopen(File_Name.c_str(),"r");
	char buf[CHAR_LEN];
	char data[CHAR_LEN];
	while ( fgets(buf,CHAR_LEN,samcfile) )
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
	fclose(samcfile);
	/*}}}*/
	/*Set Global Value{{{*/
	printf("Setting Global Value...\n");
	k=0;
	string samc_rootfilename=inputdata[k++];
	string userdefgen_filename=inputdata[k++];
	int Num_Of_Events;
	if ( IsANumber(userdefgen_filename) ) {
		Num_Of_Events=atoi(userdefgen_filename.c_str());
		userdefgen_filename="";
	}
	else {
		Num_Of_Events=atoi(inputdata[k++].c_str());
	}
	
	IsMyModel=atoi(inputdata[k++].c_str());
	IsMultiScat=atoi(inputdata[k++].c_str());
	IsEnergyLoss=atoi(inputdata[k++].c_str());
	Which_Kin=atoi(inputdata[k++].c_str());
	HRS_L=atof(inputdata[k++].c_str());
	Q1_Exit_Eff_L=atof(inputdata[k++].c_str());
	D_Entrance_Eff_L=atof(inputdata[k++].c_str());
	D_Exit_Eff_L=atof(inputdata[k++].c_str());
	Q3_Entrance_Eff_L=atof(inputdata[k++].c_str());
	Q3_Exit_Eff_L=atof(inputdata[k++].c_str());
	FP_Eff_L=atof(inputdata[k++].c_str());
	Q1_Radius=atof(inputdata[k++].c_str());
	D_X_Radius=atof(inputdata[k++].c_str());
	D_Y_L=atof(inputdata[k++].c_str());
	Q3_Entrance_Radius=atof(inputdata[k++].c_str());
	Q3_Exit_Radius=atof(inputdata[k++].c_str());
	VDC_Res_x=atof(inputdata[k++].c_str());
	VDC_Res_y=atof(inputdata[k++].c_str());
	VDC_Res_th=atof(inputdata[k++].c_str());
	VDC_Res_ph=atof(inputdata[k++].c_str());
	D_x=atof(inputdata[k++].c_str());
	D_y=atof(inputdata[k++].c_str());
	delta_dp=atof(inputdata[k++].c_str());
	delta_th=atof(inputdata[k++].c_str());
	delta_ph=atof(inputdata[k++].c_str());
	Which_Beam_Profile=atoi(inputdata[k++].c_str());
	raster_x_size=atof(inputdata[k++].c_str());
	raster_y_size=atof(inputdata[k++].c_str());
	gaus_x_sigma=atof(inputdata[k++].c_str());
	gaus_y_sigma=atof(inputdata[k++].c_str());
	beam_x_center=atof(inputdata[k++].c_str());
	beam_y_center=atof(inputdata[k++].c_str());
	z0=atof(inputdata[k++].c_str());
	T_L=atof(inputdata[k++].c_str());
	T_H=atof(inputdata[k++].c_str());
	Fermi_Moment=atof(inputdata[k++].c_str());
	NIE=atof(inputdata[k++].c_str());
	DEL_SEP=atof(inputdata[k++].c_str());
	IsRadCor=atoi(inputdata[k++].c_str());
	Form_Factor_Id=atoi(inputdata[k++].c_str());
	IsInternal_Bre=atoi(inputdata[k++].c_str());
	IsAppro=atoi(inputdata[k++].c_str());
	IsExternal_Bre=atoi(inputdata[k++].c_str());
	IsMulti_Photon_Cor=atoi(inputdata[k++].c_str());
	DeltaE=atof(inputdata[k++].c_str());
	E0=atof(inputdata[k++].c_str());
	P0=atof(inputdata[k++].c_str());
	if ( IsDebug )
	{
		if ( IsMyModel )
			printf("Use my own one-body model.\n");
		else
			printf("Use old model.\n");
		if ( IsMultiScat )
			printf("Multi-Scattering Enabled.\n");
		else
			printf("Multi-Scattering Disabled.\n");
		if ( IsEnergyLoss )
			printf("Energy Loss Enabled.\n");
		else
			printf("Energy Loss Disabled.\n");
		switch ( Which_Kin )
		{
			case 1:
				printf("Elastic.\n");
				break;
			case 2:
				printf("Quasi-Elastic.\n");
				break;
			case 0:
			default:
				printf("Phase Space.\n");
				break;
		}
		printf("%-*s=%*s %-*s %-*s\n",  15,"samc_rootfilename",10,samc_rootfilename.c_str(),8,"",  40,"(Save Results to Root File Name)");
		if ( !userdefgen_filename.empty() ) {
			printf("%-*s=%*s %-*s %-*s\n",  15,"userdefgen_filename",10,userdefgen_filename.c_str(),8,"",  40,"(User-Defined Generator File Name)");
		}
		printf("%-*s=%*d %-*s %-*s\n",  15,"Num_Of_Events",    10,Num_Of_Events,            8,"",  40,"(Number Of Events)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"HRS_L",            10,HRS_L,                    8,"cm",40,"(Length of HRS)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"Q1_Exit_Eff_L",            10,Q1_Exit_Eff_L,                    8,"cm",40,"(Length of Magnetic)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"D_Entrance_Eff_L",            10,D_Entrance_Eff_L,                    8,"cm",40,"(Length of Magnetic)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"D_Exit_Eff_L",            10,D_Exit_Eff_L,                    8,"cm",40,"(Length of Magnetic)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"Q3_Entrance_Eff_L",            10,Q3_Entrance_Eff_L,                    8,"cm",40,"(Length of Magnetic)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"Q3_Exit_Eff_L",            10,Q3_Exit_Eff_L,                    8,"cm",40,"(Length of Magnetic)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"FP_Eff_L",            10,FP_Eff_L,                    8,"cm",40,"(Length of Magnetic)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"VDC_Res_x",        10,VDC_Res_x,                8,"cm",40,"(Resolution of VDC x)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"VDC_Res_y",        10,VDC_Res_y,                8,"cm",40,"(Resolution of VDC y)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"VDC_Res_th",       10,VDC_Res_th,               8,"mr",40,"(Resolution of VDC th)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"VDC_Res_ph",       10,VDC_Res_ph,               8,"mr",40,"(Resolution of VDC ph)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"delta_dp",         10,delta_dp*100,             8,"%", 40,"(full width of \xce\x94(dp) for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"delta_th",         10,delta_th,                 8,"",40,"(full width of \xce\x94(tan(\xce\xb8tg)) for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"delta_ph",         10,delta_ph,                 8,"",40,"(full width of \xce\x94(tan(\xcf\x86tg)) for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"raster_x_size",    10,raster_x_size,            8,"cm",40,"(raster x full size for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"raster_y_size",    10,raster_y_size,            8,"cm",40,"(raster y full size for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"beam_x_center",    10,beam_x_center,            8,"cm",40,"(beam x center for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"beam_y_center",    10,beam_y_center,            8,"cm",40,"(beam y center for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"z0",               10,z0,                       8,"cm",40,"(target center for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"T_L",              10,T_L,                      8,"cm",40,"(target length for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"T_H",              10,T_H,                      8,"cm",40,"(target height for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"Fermi_Moment",              10,Fermi_Moment,                      8,"MeV",40,"(fermi momentum of target)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"NIE",              10,NIE,                      8,"MeV",40,"(nucleon separation energy)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"DEL_SEP",              10,DEL_SEP,                      8,"MeV",40,"(delta separation energy)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"D_x",              10,D_x,                      8,"cm",40,"(X offset(in TCS) between TCS and HCS)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"D_y",              10,D_y,                      8,"cm",40,"(Y offset(in TCS) between TCS and HCS)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"E0",               10,E0,                       8,"MeV",40,"(Incident Beam Energy for generator)");
		printf("%-*s=%*.2f %-*s %-*s\n",15,"P0",               10,P0,                       8,"MeV",40,"(HRS Setting Momentum for generator)");
	}
	/*}}}*/

	int err;
	double randnum;
	ifstream userdefgen_file;
	string comment;
	double usrgen[6];
	if ( !userdefgen_filename.empty() ) {
		userdefgen_file.open(userdefgen_filename.c_str());
		if ( !userdefgen_file.good() ) {
			printf("[Error %s: Line %d] %s File Not Found.\n",__FILE__,__LINE__,userdefgen_filename.c_str());
			return 10;
		}
		getline(userdefgen_file,comment);//read comment;

	}
	TFile* f=new TFile(samc_rootfilename.c_str(),"recreate");
	TTree* T=new TTree("SAMC","Tree with Acceptance Simulation");
	/*Root output variables{{{*/
	//ro=root output
	int roId; //Event Id
	double rox_tg_gen;
	double roy_tg_gen;
	double roth_tg_gen;
	double roph_tg_gen;
	double roreactz_gen;
	double rodp_gen;
	double rox_tg_ref;
	double roy_tg_ref;
	double roth_tg_ref;
	double roph_tg_ref;
	double rodp_ref;
	double rox_fp;
	double roy_fp;
	double roth_fp;
	double roth_fp_no_ort;
	double roph_fp;
	double rox_tg_rec;
	double roy_tg_rec;
	double roth_tg_rec;
	double roph_tg_rec;
	double roreactz_rec;
	double rodp_rec;
	double roE_s;//=s(3) MeV
	double roE_p;
	double roTheta;
	double rocs_M;
	double rorvalue;
	double roq1ex_x;
	double roq1ex_y;
	double rodent_x;
	double rodent_y;
	double rodext_x;
	double rodext_y;
	double roq3en_x;
	double roq3en_y;
	double roq3ex_x;
	double roq3ex_y;
	int roIsPassedQ1Ex;
	int roIsPassedDipoleEn;
	int roIsPassedDipoleEx;
	int roIsPassedQ3En;
	int roIsPassedQ3Ex;

	double rocs_qfs;//quasi-free scattering cross section
	double rocs_del;//Delta
	double rocs_x;  //DIS
	double rocs_r1500; //Resonace 1500 MeV
	double rocs_r1700; //Resonace 1700 MeV
	double rocs_2N; //Dip region
	double rocs_q;//sum of all the cross section above
	double rocs_p; //sigma_p (peak approximation for internal radation)
	double rocs_b_bt; //the part before target in sigma_b
	double rocs_b_at; //the part after target in sigma_b
	double rocs_b; //sigma_b (external bremsstrahlung ionization)
	double rocs_ex; //sigma_ex (extract internal radation)
	double romp_cor; //Multiple Photon Correction
	double rotail; //radiative tail
	double roR; //Phys.Rev.D 12 1884 A83
	double rocs_qt;
	double rocs_Final;
	int roIsPassed; //if the seed passes through the magnet
	int roIsQualified; //if the reconstruct in the gen range
	double robeam_x;
	double robeam_y;
	double rox_collimator;
	double roy_collimator;


	//int roZ; //Atomic Number of Target
	//double roA; //Atomic Weight of Target(g/mol)
	//double roT; //target thickness(g/cm^2)
	//double roTR; //target thickness(unit in rad_len)
	//int roZ_i; //Atomic Number of init window
	//double roA_i; //Atomic Weight of init window(g/mol)
	//double roT_i; //target thickness of init window(g/cm^2)
	//double roTR_i; //target thickness of init window(g/cm^2)
	//int roZ_f; //Atomic Number of final window
	//double roA_f; //Atomic Weight of final window(g/mol)
	//double roT_f; //target thickness of final window(g/cm^2)
	//double roTR_f; //target thickness of final window(rad_len)
	//double roT_theta; //target angle(deg)
	//double rotheta; //scattering angle(deg)
	//double rotb;//thickness before the scattering in radation length
	//double rota;//thickness after the scattering in radation length
	//double rotheta_rad;//scattering angle in SAMC
	//double rosinsq;//sin(theta/2)^2
	//double roM_T;//Mass of Target MeV
	//double roM_f;//Mass of final hardonic system MeV
	//double roE_p;//=p(3) MeV
	//double rocostheta_s;//cos(theta_s) angle between s and u
	//double rocostheta_p;//cos(theta_p) angle between p and u
	//double roE_p_max; //E' max from E_s MeV
	//double row_s;//A50 MeV
	//double row_p;//A51 MeV
	//double rov_s; //A53
	//double rov_p;
	//double robtb;//b*tb in radiation length
	//double robta;//b*ta in radiation length
	//double robtr;//b*tr (equivalent radiator)
	//double roQ2;// MeV*MeV
	//double roq2;//q2=-Q2 MeV*MeV
	//double rocs_p; //sigma_p (peak approximation for internal radation)
	//double rocs_b_bt; //the part before target in sigma_b
	//double rocs_b_at; //the part after target in sigma_b
	//double rocs_b; //sigma_b (external bremsstrahlung ionization)
	//double rocs_ex; //sigma_ex (extract internal radation)
	//double romp_cor; //Multiple Photon Correction
	//double rotail;
	//double rotail_err; //error for radiative tail
	/*}}}*/
	/*Root output branch{{{*/
	T->Branch("Id",&roId,"roId/I");
	T->Branch("x_tg_gen",&rox_tg_gen,"rox_tg_gen/D");
	T->Branch("y_tg_gen",&roy_tg_gen,"roy_tg_gen/D");
	T->Branch("th_tg_gen",&roth_tg_gen,"roth_tg_gen/D");
	T->Branch("ph_tg_gen",&roph_tg_gen,"roph_tg_gen/D");
	T->Branch("reactz_gen",&roreactz_gen,"roreactz_gen/D");
	T->Branch("dp_gen",&rodp_gen,"rodp_gen/D");
	T->Branch("x_tg_ref",&rox_tg_ref,"rox_tg_ref/D");
	T->Branch("y_tg_ref",&roy_tg_ref,"roy_tg_ref/D");
	T->Branch("th_tg_ref",&roth_tg_ref,"roth_tg_ref/D");
	T->Branch("ph_tg_ref",&roph_tg_ref,"roph_tg_ref/D");
	T->Branch("dp_ref",&rodp_ref,"rodp_ref/D");
	T->Branch("x_fp",&rox_fp,"rox_fp/D");
	T->Branch("y_fp",&roy_fp,"roy_fp/D");
	T->Branch("th_fp",&roth_fp,"roth_fp/D");
	T->Branch("th_fp_no_ort",&roth_fp_no_ort,"roth_fp_no_ort/D");
	T->Branch("ph_fp",&roph_fp,"roph_fp/D");
	T->Branch("x_tg_rec",&rox_tg_rec,"rox_tg_rec/D");
	T->Branch("y_tg_rec",&roy_tg_rec,"roy_tg_rec/D");
	T->Branch("th_tg_rec",&roth_tg_rec,"roth_tg_rec/D");
	T->Branch("ph_tg_rec",&roph_tg_rec,"roph_tg_rec/D");
	T->Branch("reactz_rec",&roreactz_rec,"roreactz_rec/D");
	T->Branch("dp_rec",&rodp_rec,"rodp_rec/D");
	T->Branch("E_s",&roE_s,"roE_s/D");
	T->Branch("E_p",&roE_p,"roE_p/D");
	T->Branch("Theta",&roTheta,"roTheta/D");
	T->Branch("cs_M",&rocs_M,"rocs_M/D");
	T->Branch("rvalue",&rorvalue,"rorvalue/D");
	T->Branch("q1ex_x",&roq1ex_x,"roq1ex_x/D");
	T->Branch("dent_x",&rodent_x,"rodent_x/D");
	T->Branch("dext_x",&rodext_x,"rodext_x/D");
	T->Branch("q3en_x",&roq3en_x,"roq3en_x/D");
	T->Branch("q3ex_x",&roq3ex_x,"roq3ex_x/D");
	T->Branch("q1ex_y",&roq1ex_y,"roq1ex_y/D");
	T->Branch("dent_y",&rodent_y,"rodent_y/D");
	T->Branch("dext_y",&rodext_y,"rodext_y/D");
	T->Branch("q3en_y",&roq3en_y,"roq3en_y/D");
	T->Branch("q3ex_y",&roq3ex_y,"roq3ex_y/D");
	T->Branch("IsPassedQ1Ex",&roIsPassedQ1Ex,"roIsPassedQ1Ex/I");
	T->Branch("IsPassedDipoleEn",&roIsPassedDipoleEn,"roIsPassedDipoleEn/I");
	T->Branch("IsPassedDipoleEx",&roIsPassedDipoleEx,"roIsPassedDipoleEx/I");
	T->Branch("IsPassedQ3En",&roIsPassedQ3En,"roIsPassedQ3En/I");
	T->Branch("IsPassedQ3Ex",&roIsPassedQ3Ex,"roIsPassedQ3Ex/I");

	T->Branch("cs_qfs",&rocs_qfs,"rocs_qfs/D");
	T->Branch("cs_del",&rocs_del,"rocs_del/D");
	T->Branch("cs_x",&rocs_x,"rocs_x/D");
	T->Branch("cs_r1500",&rocs_r1500,"rocs_r1500/D");
	T->Branch("cs_r1700",&rocs_r1700,"rocs_r1700/D");
	T->Branch("cs_2N",&rocs_2N,"rocs_2N/D");
	T->Branch("cs_q",&rocs_q,"rocs_q/D");
	T->Branch("cs_p",&rocs_p,"rocs_p/D");
	T->Branch("cs_b_bt",&rocs_b_bt,"rocs_b_bt/D");
	T->Branch("cs_b_at",&rocs_b_at,"rocs_b_at/D");
	T->Branch("cs_b",&rocs_b,"rocs_b/D");
	T->Branch("cs_ex",&rocs_ex,"rocs_ex/D");
	T->Branch("mp_cor",&romp_cor,"romp_cor/D");
	T->Branch("tail",&rotail,"rotail/D");
	T->Branch("R",&roR,"roR/D");
	T->Branch("cs_qt",&rocs_qt,"rocs_qt/D");
	T->Branch("cs_Final",&rocs_Final,"rocs_Final/D");
	T->Branch("IsPassed",&roIsPassed,"roIsPassed/I");
	T->Branch("IsQualified",&roIsQualified,"roIsQualified/I");
	T->Branch("beam_x",&robeam_x,"robeam_x/D");
	T->Branch("beam_y",&robeam_y,"robeam_y/D");
	T->Branch("x_collimator",&rox_collimator,"rox_collimator/D");
	T->Branch("y_collimator",&roy_collimator,"roy_collimator/D");

	//T->Branch("Z",&roZ,"roZ/I");
	//T->Branch("A",&roA,"roA/D");
	//T->Branch("T",&roT,"roT/D");
	//T->Branch("TR",&roTR,"roTR/D");
	//T->Branch("Z_i",&roZ_i,"roZ_i/I");
	//T->Branch("A_i",&roA_i,"roA_i/D");
	//T->Branch("T_i",&roT_i,"roT_i/D");
	//T->Branch("TR_i",&roTR_i,"roTR_i/D");
	//T->Branch("Z_f",&roZ_f,"roZ_f/I");
	//T->Branch("A_f",&roA_f,"roA_f/D");
	//T->Branch("T_f",&roT_f,"roT_f/D");
	//T->Branch("TR_f",&roTR_f,"roTR_f/D");
	//T->Branch("T_theta",&roT_theta,"roT_theta/D");
	//T->Branch("theta",&rotheta,"rotheta/D");
	//T->Branch("tb",&rotb,"rotb/D");
	//T->Branch("ta",&rota,"rota/D");
	//T->Branch("theta_rad",&rotheta_rad,"rotheta_rad/D");
	//T->Branch("sinsq",&rosinsq,"rosinsq/D");
	//T->Branch("M_T",&roM_T,"roM_T/D");
	//T->Branch("M_f",&roM_f,"roM_f/D");
	//T->Branch("E_p",&roE_p,"roE_p/D");
	//T->Branch("costheta_s",&rocostheta_s,"rocostheta_s/D");
	//T->Branch("costheta_p",&rocostheta_p,"rocostheta_p/D");
	//T->Branch("E_p_max",&roE_p_max,"roE_p_max/D");
	//T->Branch("w_s",&row_s,"row_s/D");
	//T->Branch("w_p",&row_p,"row_p/D");
	//T->Branch("v_s",&rov_s,"rov_s/D");
	//T->Branch("v_p",&rov_p,"rov_p/D");
	//T->Branch("btb",&robtb,"robtb/D");
	//T->Branch("bta",&robta,"robta/D");
	//T->Branch("btr",&robtr,"robtr/D");
	//T->Branch("Q2",&roQ2,"roQ2/D");
	//T->Branch("q2",&roq2,"roq2/D");
	//T->Branch("cs_p",&rocs_p,"rocs_p/D");
	//T->Branch("cs_b_bt",&rocs_b_bt,"rocs_b_bt/D");
	//T->Branch("cs_b_at",&rocs_b_at,"rocs_b_at/D");
	//T->Branch("cs_b",&rocs_b,"rocs_b/D");
	//T->Branch("cs_ex",&rocs_ex,"rocs_ex/D");
	//T->Branch("mp_cor",&romp_cor,"romp_cor/D");
	//T->Branch("tail",&rotail,"rotail/D");
	//T->Branch("tail_err",&rotail,"rotail_err/D");
	/*}}}*/
	/*Process Event by Event{{{*/
	TStopwatch timer;
	timer.Start();
	//only init rfunction once
	left_init_r_function__();
	right_init_r_function__();
	fail_events=0;
	for ( i=0; i<Num_Of_Events; i++ )
	{
		j=k;
		SAMCEvent* Event=new SAMCEvent();
		Event->Id=i;
		Event->E_s=gRandom->Gaus(E0,E0*3e-5);//beam dispersion is 3e-5
		Event->theta=atof(inputdata[j++].c_str());
		Event->Target.Name=inputdata[j++].c_str();
		Event->Target.Z=atoi(inputdata[j++].c_str());
		Event->Target.A=atof(inputdata[j++].c_str());
		Event->Target.T=atof(inputdata[j++].c_str());
		Event->Target.rho=atof(inputdata[j++].c_str());
		Event->Target.L=T_L;//in HCS
		Event->Win_i.Name=inputdata[j++].c_str();
		Event->Win_i.Z=atoi(inputdata[j++].c_str());
		Event->Win_i.A=atof(inputdata[j++].c_str());
		Event->Win_i.T=atof(inputdata[j++].c_str());
		Event->Win_i.rho=atof(inputdata[j++].c_str());
		Event->Win_f.Name=inputdata[j++].c_str();
		Event->Win_f.Z=atoi(inputdata[j++].c_str());
		Event->Win_f.A=atof(inputdata[j++].c_str());
		Event->Win_f.T=atof(inputdata[j++].c_str());
		Event->Win_f.rho=atof(inputdata[j++].c_str());
		Event->T_theta=atof(inputdata[j++].c_str());
		while ( atof(inputdata[j].c_str())>=0 )
		{
			Event->AddOneMaterial(Event->Win_Before_Mag,atof(inputdata[j++].c_str()),atof(inputdata[j++].c_str()),atof(inputdata[j++].c_str()),atof(inputdata[j++].c_str()),atoi(inputdata[j++].c_str()),inputdata[j++].c_str());
		}
		j++;
		while ( atof(inputdata[j].c_str())>=0 )
		{
			Event->AddOneMaterial(Event->Win_After_Mag,atof(inputdata[j++].c_str()),atof(inputdata[j++].c_str()),atof(inputdata[j++].c_str()),atof(inputdata[j++].c_str()),atoi(inputdata[j++].c_str()),inputdata[j++].c_str());
		}

		if ( !userdefgen_filename.empty() ) {
			int buf_i=0;
			for ( buf_i = 0; buf_i < 6; ++buf_i ) {
				userdefgen_file>>usrgen[buf_i];
			}
			Event->beam_x=usrgen[0];
			Event->beam_y=usrgen[1];
			Event->reactz_gen=usrgen[2];
			Event->th_tg_gen=usrgen[3];
			Event->ph_tg_gen=usrgen[4];
			Event->dp_gen=usrgen[5];
		}
		else {
			if ( Which_Beam_Profile==0 ) {
				Event->beam_x=beam_x_center+(gRandom->Rndm()-0.5)*raster_x_size;
				Event->beam_y=beam_y_center+(gRandom->Rndm()-0.5)*raster_y_size;
			}
			else if ( Which_Beam_Profile==1 ) {
				Event->beam_x=gRandom->Gaus(beam_x_center,gaus_x_sigma);
				Event->beam_y=gRandom->Gaus(beam_y_center,gaus_y_sigma);
			}
			else if ( Which_Beam_Profile==2 ) {
				Event->beam_x=beam_x_center+(gRandom->Rndm()-0.5)*raster_x_size;
				Event->beam_y=beam_y_center+(gRandom->Rndm()-0.5)*raster_y_size;
				Event->beam_x=gRandom->Gaus(Event->beam_x,gaus_x_sigma);
				Event->beam_y=gRandom->Gaus(Event->beam_y,gaus_y_sigma);
			}
			Event->reactz_gen=z0+(gRandom->Rndm()-0.5)*T_L*cos(Event->T_theta*DegToRad());
			//cout<<"---DEBUG---"<<endl;
			//cout<<"z0="<<z0<<endl;
			//cout<<"T_L="<<T_L<<endl;
			//cout<<"beam_x="<<Event->beam_x<<endl;
			//cout<<"reactz_gen="<<Event->reactz_gen<<endl;
			//cout<<"T_theta="<<Event->T_theta<<endl;
			//cout<<"-----------"<<endl;
			Event->th_tg_gen=(gRandom->Rndm()-0.5)*delta_th;
			Event->ph_tg_gen=(gRandom->Rndm()-0.5)*delta_ph;
			Event->dp_gen=(gRandom->Rndm()-0.5)*delta_dp;
		}

		err=Event->Process();
		if ( err==-1 )
		{
			exit(err);
		}
		//else if ( err>0 && Event->IsPassed==1 )
		//{
		//	i-=err;//This event is bad, simulate again. To make sure total good number of events == one user inputs.
		//}
		else
		{
			//i-=err;//if err>0 This event is bad, simulate again. To make sure total passed number of events == one user inputs.
			if ( Event->IsPassed==0 && userdefgen_filename.empty() ) {
				++fail_events; //but also record the bad one to calculate acceptance
				--i;
			}
			if ( (i+1)%1000==0 || i==0 || i==(Num_Of_Events-1) ) {
				if ( userdefgen_filename.empty() ) {
					printf("\t%d Good Events Simulated. And %d Bad Events Saved.\n",i+1,fail_events);
				}
				else {
					printf("\t%d User-Defined Events Simulated.\n",i+1);
				}
			}
			/*Save results to TTree{{{*/
			roId=Event->Id;
			rox_tg_gen=Event->x_tg_gen/100;
			roy_tg_gen=Event->y_tg_gen/100;
			roth_tg_gen=Event->th_tg_gen;
			roph_tg_gen=Event->ph_tg_gen;
			roreactz_gen=Event->reactz_gen/100;
			rodp_gen=Event->dp_gen;
			rox_tg_ref=Event->x_tg_ref/100;
			roy_tg_ref=Event->y_tg_ref/100;
			roth_tg_ref=Event->th_tg_ref;
			roph_tg_ref=Event->ph_tg_ref;
			rodp_ref=Event->dp_ref;
			rox_fp=Event->x_fp/100;
			roy_fp=Event->y_fp/100;
			roth_fp=Event->th_fp;
			roth_fp_no_ort=Event->th_fp_no_ort;
			roph_fp=Event->ph_fp;
			rox_tg_rec=Event->x_tg_rec/100;//cm->meter
			roy_tg_rec=Event->y_tg_rec/100;
			roth_tg_rec=Event->th_tg_rec;
			roph_tg_rec=Event->ph_tg_rec;
			roreactz_rec=Event->reactz_rec/100;
			rodp_rec=Event->dp_rec;
			roE_s=Event->E_s/1000;//MeV to GeV
			roE_p=Event->E_p/1000;
			roTheta=Event->Angle*RadToDeg();
			rocs_M=Event->cs_M;
			rorvalue=Event->rvalue;
			roq1ex_x=Event->q1ex[0];
			rodent_x=Event->dent[0];
			rodext_x=Event->dext[0];
			roq3en_x=Event->q3en[0];
			roq3ex_x=Event->q3ex[0];
			roq1ex_y=Event->q1ex[1];
			rodent_y=Event->dent[1];
			rodext_y=Event->dext[1];
			roq3en_y=Event->q3en[1];
			roq3ex_y=Event->q3ex[1];
			roIsPassedQ1Ex=Event->IsPassedQ1Ex;
			roIsPassedDipoleEn=Event->IsPassedDipoleEn;
			roIsPassedDipoleEx=Event->IsPassedDipoleEx;
			roIsPassedQ3En=Event->IsPassedQ3En;
			roIsPassedQ3Ex=Event->IsPassedQ3Ex;

			rocs_qfs=Event->cs_qfs;
			rocs_del=Event->cs_del;
			rocs_x=Event->cs_x;
			rocs_r1500=Event->cs_r1500;
			rocs_r1700=Event->cs_r1700;
			rocs_2N=Event->cs_2N;
			rocs_q=Event->cs_q;
			rocs_p=Event->cs_p;
			rocs_b_bt=Event->cs_b_bt;
			rocs_b_at=Event->cs_b_at;
			rocs_b=Event->cs_b;
			rocs_ex=Event->cs_ex;
			romp_cor=Event->mp_cor;
			rotail=Event->tail;
			roR=Event->R;
			rocs_qt=Event->cs_qt;
			rocs_Final=Event->cs_Final;
			roIsPassed=Event->IsPassed;
			roIsQualified=Event->IsQualified;
			robeam_x=Event->beam_x/100;//cm to m
			robeam_y=Event->beam_y/100;//cm to m
			rox_collimator=Event->x_collimator/100;//cm to m
			roy_collimator=Event->y_collimator/100;//cm to m

			//		roZ=Event->Z;
			//		roA=Event->A;
			//		roT=Event->T;
			//		roTR=Event->TR;
			//		roZ_i=Event->Z_i;
			//		roA_i=Event->A_i;
			//		roT_i=Event->T_i;
			//		roTR_i=Event->TR_i;
			//		roZ_f=Event->Z_f;
			//		roA_f=Event->A_f;
			//		roT_f=Event->T_f;
			//		roTR_f=Event->TR_f;
			//		roT_theta=Event->T_theta;
			//		rotheta=Event->theta;
			//		rotb=Event->tb;
			//		rota=Event->ta;
			//		rotheta_rad=Event->theta_rad;
			//		rosinsq=Event->sinsq;
			//		roM_T=Event->M_T;
			//		roM_f=Event->M_f;
			//		roE_p=Event->E_p;
			//		rocostheta_s=Event->costheta_s;
			//		rocostheta_p=Event->costheta_p;
			//		roE_p_max=Event->E_p_max;
			//		row_s=Event->w_s;
			//		row_p=Event->w_p;
			//		rov_s=Event->v_s;
			//		rov_p=Event->v_p;
			//		robtb=Event->btb;
			//		robta=Event->bta;
			//		robtr=Event->btr;
			//		roQ2=Event->Q2;
			//		roq2=Event->q2;
			//		rocs_p=Event->cs_p;
			//		rocs_b_bt=Event->cs_b_bt;
			//		rocs_b_at=Event->cs_b_at;
			//		rocs_b=Event->cs_b;
			//		rocs_ex=Event->cs_ex;
			//		romp_cor=Event->mp_cor;
			//		rotail=Event->tail;
			//		rotail_err=Event->tail_err;
			T->Fill();
			/*}}}*/
		}

		if ( IsDebug && err==0 )
			Event->Print();

		delete Event;
	}
	if ( !userdefgen_filename.empty() ) {
		userdefgen_file.close();
	}
	/*}}}*/

	T->Write();
	delete T;
	f->Close();
	delete f;
	inputdata.clear();
	printf("Results saved at %s\n",samc_rootfilename.c_str());
	printf("Time at the end of job = %f seconds\n",timer.CpuTime());

	return 0;
}
/*}}}*/
