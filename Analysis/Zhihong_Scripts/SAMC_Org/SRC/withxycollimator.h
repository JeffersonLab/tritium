#ifndef RAD_PHYSICAL_CONSTANT_H
#include "SAMC.h"
//#include "Left_funcs.h"
/*Call fortran codes{{{*/
extern "C"
{
	//must be like the following
	//change l(r)txfit_ in fortran to real function ...
	//otherwise no result is passed.
	//for left arm
	float ltxfit_(float*,int&);
	float ldelta_(float*,int&);
	float ltheta_(float*,int&);
	float lphi_(float*,int&);
	float ly00_(float*,int&);
	double x_e_q1ex__(float*,int&);
	double y_e_q1ex__(float*,int&);
	double x_e_dent__(float*,int&);
	double y_e_dent__(float*,int&);
	double x_e_dext__(float*,int&);
	double y_e_dext__(float*,int&);
	double x_e_q3en__(float*,int&);
	double y_e_q3en__(float*,int&);
	double x_e_q3ex__(float*,int&);
	double y_e_q3ex__(float*,int&);
	double x_e_fp__(float*,int&);
	double y_e_fp__(float*,int&);
	double p_e_fp__(float*,int&);
	double t_e_fp__(float*,int&);
	void left_init_r_function__();
	void left_getindex__(double*,double*,int*);
	double left_rfunction__(float*,float*,float*,float*);

	//for right arm
	float rtxfit_(float*,int&);
	float rdelta_(float*,int&);
	float rtheta_(float*,int&);
	float rphi_(float*,int&);
	float ry00_(float*,int&);
	double x_h_q1ex__(float*,int&);
	double y_h_q1ex__(float*,int&);
	double x_h_dent__(float*,int&);
	double y_h_dent__(float*,int&);
	double x_h_dext__(float*,int&);
	double y_h_dext__(float*,int&);
	double x_h_q3en__(float*,int&);
	double y_h_q3en__(float*,int&);
	double x_h_q3ex__(float*,int&);
	double y_h_q3ex__(float*,int&);
	double x_h_fp__(float*,int&);
	double y_h_fp__(float*,int&);
	double p_h_fp__(float*,int&);
	double t_h_fp__(float*,int&);
	void right_init_r_function__();
	void right_getindex__(double*,double*,int*);
	double right_rfunction__(float*,float*,float*,float*);
}
/*}}}*/
#define MSIZE 5
/*struct Material{{{*/
struct Material
{
	string Name; //Name
	int Z;//Z
	double A;//A
	double M;//Mass
	double T;//Thickness g/cm^2
	double TR;//Thickness in Rad_Len
	double rho;//density g/cm^3
	double bt;//bt
	double X0; //Radiation Length g/cm^2
	double L;//Length:distance along Z axix in TCS
}
;/*}}}*/
#endif

/*class SAMCEvent{{{*/
class SAMCEvent
{
	public:
		/*SAMCEvent(){{{*/
		SAMCEvent()
		{
			Init();
		}
		/*}}}*/
		/*virtual ~SAMCEvent(){{{*/
		virtual ~SAMCEvent()
		{
			Win_Before_Mag.clear();
			Win_After_Mag.clear();
		}
		/*}}}*/

		SAMCEvent(SAMCEvent const&){};
		SAMCEvent& operator=(SAMCEvent const&){};

		/*int Process(){{{*/
		int Process()
		{
			int Num_Event_Add=0;
			Generator();//generate original target variables
			Num_Event_Add=RefineTg();//transport to front of magnetic, then back to get refined target variables
			
			Num_Event_Add+=ToFp(x_tg_ref,y_tg_ref,th_tg_ref,ph_tg_ref,dp_ref);//transfer to focal plane using John.LeRose matrix
			if ( Num_Event_Add>0 ) {
				IsPassed=0;
				IsQualified=0;
			}
			else
				Num_Event_Add=ReconstructTg(x_fp,y_fp,th_fp,ph_fp,x_tg_ref);
			return 0;
		}
		/*}}}*/

		/*void AddOneMaterial(vector<Material>& aWin,const double& aX0,const double& arho,const double& aL,const double& aA,const int& aZ,string aName){{{*/
		void AddOneMaterial(vector<Material>& aWin,const double& aX0,const double& arho,const double& aL,const double& aA,const int& aZ,string aName)
		{
			//To make sure read correctly, so I reverse the order compared with input file.
			Material a;
			a.Name=aName;
			a.Z=aZ;
			a.A=aA;
			a.L=aL;
			a.rho=arho;
			a.T=aL*arho;
			a.X0=aX0;
			if ( fabs(aX0)<1e-10 )
			{
				a.TR=0;
				a.bt=0;
			}
			else
			{
				a.TR=a.T/aX0;
				a.bt=b(aZ)*aL*arho/aX0;
			}
			aWin.push_back(a);
		}
		/*}}}*/


	private:
		/*void Generator(){{{*/
		void Generator()
		{
			//set value for Member Data derived from variables from file
			//File provides E_s,theta,Target.(Z,A,T,rho) Win_i.(Z,A,T,rho)
			//Win_f.(Z,A,T,rho) T_theta
			//Win_Before_Mag(Name,Z,A,L,rho,X0)
			//Win_After_Mag(Name,Z,A,L,rho,X0)
			//beam_x,beam_y,reactz_gen,th_tg_gen,ph_tg_gen,dp_gen
			//z0,HRS_L,VDC_Res(x,y,th,ph),D_(x,y),T_L,P0
			//IsMultiScat,IsEnergyLoss,Which_Kin,FP_Eff_L

			/*Set Material.(Z,A,M,X0,T,TR,bt){{{*/
			SetMaterial(Target);
			SetMaterial(Win_i);
			SetMaterial(Win_f);
			Target.TR=Target.TR;
			Target.bt=b(Target.Z)*Target.TR;
			/*}}}*/

			/*Set Beam Info(HCS and TCS){{{*/
			//know beam_x,beam_y,reactz_gen,E_s,theta,HRS_L
			//Set s,s_TCS,x_tg_gen,y_tg_gen,p_TCS,p_P,p_P_TCS
			theta_rad=theta*DegToRad();//rad
			s(0)=beam_x; //cm
			s(1)=beam_y; //cm
			s(2)=reactz_gen; //cm
			target_center_HCS(0)=0;
			target_center_HCS(1)=0;
			target_center_HCS(2)=z0;
			target_edgepoint_HCS(0)=theta/fabs(theta)*T_H/2*cos(T_theta*DegToRad());//if theta>0,T_H/2, if<0, -T_H/2
			target_edgepoint_HCS(1)=0;
			target_edgepoint_HCS(2)=z0+T_L/2*cos(T_theta*DegToRad());
			TLorentzVector lte;//the position of edge at HCS at T_theta=0
			lte(0)=beam_x;
			lte(1)=beam_y;
			lte(2)=z0+T_L/2*cos(T_theta*DegToRad());//here, lte is edge plane  in HCS
			lte.RotateZ(PI/2);//passive ratation around HCS, so -(-PI/2)
			lte.RotateX(theta_rad);//passive ratation around HCS, so -(-theta_rad)
			lte(0)-=D_x;
			lte(1)-=D_y;//now in TCS

			s_TCS=s;
			s_TCS.RotateZ(PI/2);//passive ratation around HCS, so -(-PI/2)
			s_TCS.RotateX(theta_rad);//passive ratation around HCS, so -(-theta_rad)
			s_TCS(0)-=D_x;
			s_TCS(1)-=D_y;

			target_center_TCS=target_center_HCS;
			target_center_TCS.RotateZ(PI/2);//passive ratation around HCS, so -(-PI/2)
			target_center_TCS.RotateX(theta_rad);//passive ratation around HCS, so -(-theta_rad)
			target_center_TCS(0)-=D_x;
			target_center_TCS(1)-=D_y;

			target_edgepoint_TCS=target_edgepoint_HCS;
			target_edgepoint_TCS.RotateZ(PI/2);//passive ratation around HCS, so -(-PI/2)
			target_edgepoint_TCS.RotateX(theta_rad);//passive ratation around HCS, so -(-theta_rad)
			target_edgepoint_TCS(0)-=D_x;
			target_edgepoint_TCS(1)-=D_y;

			p_inter_point_TCS=s_TCS;
			reactz_TCS=s_TCS.Z();

			s_TCS(0)-=s_TCS.Z()*th_tg_gen;
			s_TCS(1)-=s_TCS.Z()*ph_tg_gen;
			s_TCS(2)-=s_TCS.Z();

			x_tg_gen=s_TCS(0);
			y_tg_gen=s_TCS(1);
			//cout<<"---DEBUG---"<<endl;
			//cout<<"y_tg_gen="<<y_tg_gen<<endl;
			//cout<<"y_tg_gen_cal="<<-D_y-(reactz_gen*(cos(theta_rad)*ph_tg_gen+sin(theta_rad))-beam_x*(cos(theta_rad)-ph_tg_gen*sin(theta_rad)))<<endl;
			//cout<<"reactz_gen="<<reactz_gen<<endl;
			//cout<<"theta_rad="<<theta_rad<<endl;
			//cout<<"ph_tg_gen="<<ph_tg_gen<<endl;
			//cout<<"beam_x="<<beam_x<<endl;
			//cout<<"-----------"<<endl;
			
			//Fix me: I don't think we need this
			//Material halftarget=Target;
			//halftarget.T /= 2;
			//E_s-=Ion_Loss(E_s,Win_i);
			//E_s-=Bremss_Loss(E_s,Win_i.bt+btr);
			//E_s-=Ion_Loss(E_s,halftarget);
			//E_s-=Bremss_Loss(E_s,halftarget.bt+btr);
			s(3)=E_s;//MeV
			s_TCS(3)=E_s;//MeV

			p_TCS.SetVect(s_TCS.Vect());
			TLorentzVector lz(0,0,1,0);//in HCS
			p_P=lz;
			p_P.RotateZ(PI/2);
			p_P.RotateX(theta_rad);
			p_P.RotateY(-atan(th_tg_gen));
			p_P.RotateX(atan(ph_tg_gen));
			Angle=p_P.Angle(lz.Vect());
			Angle_Deg=Angle*RadToDeg();
			sinsq_Angle=sin(Angle/2)*sin(Angle/2);
			sinsq=sin(Angle/2)*sin(Angle/2);//sin(Angle/2)^2
			//p_P_TCS=lz;//Now think it's in TCS
			//p_P_TCS.RotateY(atan(th_tg_gen));//p_P_TCS.X()/p_P_TCS.Z()=th_tg
			//p_P_TCS.RotateX(-atan(ph_tg_gen));//p_P_TCS.Y()/p_P_TCS.Z()=ph_tg
			p_P_TCS.SetX(th_tg_gen);
			p_P_TCS.SetY(ph_tg_gen);
			p_P_TCS.SetZ(1);
			p_P_TCS.SetVect(p_P_TCS.Vect().Unit());
			/*}}}*/

			size_t i;
			size_t imax;
			Material m0;
			/*Set Windows{{{*/
			//know Win_Before_Mag,Win_f,HRS_L,theta
			//Add Target,Win_f,Air if necessary to Win_Before_Mag
			//Win_After_Mag don't need to be changed since it's vacuum.
			//for Ion_Loss Z,A,T,rho
			//for Bremss_Loss bt
			//for MultiScattering TR
			//for Transport L
			//Thickness definition diagram
			//Win_i |<--T-->|

			//Target|<--T-->|

			//Win_f ---------
			//      |<--T
			//      ---------

			//Target diagram,z0,TCS,reactz_gen,T
			//|<----T--->|
			//|  !  !   !|
			//   TCS!   !
			//      z0  !
			//          reactz_gen
			//          <-lL
			vector<Material>::iterator it=Win_Before_Mag.begin();//here only have the material before mag defined in input file.
			double lHL=0; //if lHL<HRS_Lcm or <3.57m(See Hall_A_NIM), i assume it's air
			int FirstInWinBeforeMagBlock=0;
			imax=Win_Before_Mag.size();
			
			for ( i=0; i<(imax-1); i++ )
				lHL+=Win_Before_Mag[i].L;

			//Win_Before_Mag add Win_f
			if ( Win_f.Z!=0 && Win_f.A!=0 && Win_f.rho!=0 && Win_f.TR!=0 )
			{
				m0=Win_f;//Win_f added to Win_Before_Mag
				m0.Name="Win_f";
				Win_Before_Mag.insert(it,m0);
				FirstInWinBeforeMagBlock++;
			}
			//Win_Before_Mag add Target
			it=Win_Before_Mag.begin();

			m0=Target;//target after interaction point added to Win_Before_Mag
			double lL; //assume target is rentangle lL=distance between reactz_gen and edge of target
			lL=lte(2)-p_inter_point_TCS(2);//verticle distance between interaction point and edge plane of target
			lL*=ph_tg_gen;
			lL+=lte(1);//height on edge
			if ( fabs(lL)<fabs(target_edgepoint_TCS(1)) ) {
				//in the target
				lL=(lte(2)-p_inter_point_TCS(2))/cos(atan(ph_tg_gen));
			}
			else {
				//out of target before hiting the edge of target, the edge means the downstream face.
				lL=fabs((lte(2)-p_inter_point_TCS(2))*ph_tg_gen/sin(atan(ph_tg_gen)));
			}
			m0.L=lL;
			Win_Before_Mag.insert(it,m0);
			FirstInWinBeforeMagBlock++;

			//Win_Before_Mag add Air before last material in the input file
			imax=Win_Before_Mag.size();
			if ( (lHL+Win_Before_Mag[imax-1].L)<(HRS_L) )
			{
				//http://pdg.lbl.gov/2009/AtomicNuclearProperties/HTML_PAGES/104.html
				double total=0.000124+0.755267+0.231781+0.012827;//total mass of component of air
				m0.Name="Air";
				m0.Z=int(6*0.000124/total+7*0.755267/total+8*0.231781/total+18*0.012827/total);
				m0.A=m0.Z/0.499;
				m0.L=HRS_L-lHL;
				
				m0.rho=1.205e-03;
				m0.T=m0.L*m0.rho;
				m0.X0=36.66;
				m0.bt=b(m0.Z)*m0.T/m0.X0;
				it=Win_Before_Mag.end();
				it--;
				Win_Before_Mag.insert(it,m0);
			}
			else if ( (lHL+Win_Before_Mag[imax-1].L)>HRS_L )
			{
				printf("[Error %s: Line %d] Total windows length=%f>HRS_L=%f.\n",__FILE__,__LINE__,(lHL+Win_Before_Mag[imax-1].L),HRS_L);
				exit(-3);
			}

			//following is to split the Air to 2 part, so I will check if it hits the collimator at 1.1m
			lHL=0;
			for ( i = 0; i < Win_Before_Mag.size(); ++i ) {
				lHL+=Win_Before_Mag[i].L;
				if ( lHL>110.0 ) {
					break;
				}
			}
			lHL-=Win_Before_Mag[i].L;
			int colli=i;
			Material m=Win_Before_Mag[i];
			m.L=110.-lHL;
			it=Win_Before_Mag.begin();
			for ( i = 0; i < colli; ++i ) {
				it++;
			}
			Win_Before_Mag.insert(it,m);
			Win_Before_Mag[i+1].L-=m.L;
			
			//for ( i = 0; i < Win_Before_Mag.size(); ++i ) {
			//	PrintMaterial(Win_Before_Mag[i]);
			//}
			////Correct First Material.L in Win_Before_Mag Block of inputfile
			//it=Win_Before_Mag.begin();
			//for ( i=0; i<FirstInWinBeforeMagBlock; i++ )
			//	(it+FirstInWinBeforeMagBlock)->L-=(it+i)->L;

			/*}}}*/

			/*Set Target Info(TCS){{{*/
			//know P0,dp_gen,Which_Kin,IsEnergyLoss,E_s,E_p,sinsq,Target
			//Set dp_gen,s_TCS,x_tg_gen,y_tg_gen,p_TCS,p_P,p_P_TCS
			//Set Q2,q2,btr,Win_Before_Mag[0].bt
			//x_tg_gen,y_tg_gen,see Set Beam Info(HCS and TCS)
			TRandom* tr=new TRandom();
			switch ( Which_Kin )
			{
				case 1: //elastic
					E_p=tr->Gaus(0,P0*3e-4);
					E_p+=E_s/(1+2*E_s*sinsq/Target.M);
					break;
				case 2: //quasi-elastic
				case 0: //phase default
				default:
					E_p=P0*(1+dp_gen);
			}
			delete tr;
			Q2=4*E_s*E_p*sinsq;
			q2=-Q2;
			btr=AP*(log(Q2/(ELECTRON_MASS*ELECTRON_MASS))-1);//b*t_r
			//Correct Win_Before_Mag[0](Target) bt
			Win_Before_Mag[0].bt+=btr;

			p_TCS(3)=E_p;//MeV
			p_P(3)=E_p;
			p_P_TCS(3)=E_p;
			dp_gen=(E_p-P0)/P0;


			/*}}}*/

		}
		/*}}}*/

		/*int RefineTg(){{{*/
		int RefineTg()
		{
			//return Number of Event needs to be added. just 1 ^_^
			//refine target variables from Generator because John.LeRose matrix only works for vacuum
			//know Win_Before_Mag,p_TCS,p_P_TCS
			p_TCS_ref=p_inter_point_TCS;
			p_P_TCS_ref=p_P_TCS;
			p_P_TCS_ref(3)=E_p;

			size_t i;
			size_t imax;
			vector<Material> Win_Empty;//tmp use

			Material mixture=GetMixture(Win_Before_Mag);//just need to know the length

			double offset=p_TCS(2)-p_TCS_ref(2);//L=along Z in TCS
			mixture.L=offset;

			if ( IsMyModel )
			{
				/*IsMyModel{{{*/
				TLorentzVector lPos;
				TLorentzVector lMom;
				float matrix[MSIZE];
				int msize=MSIZE;
				double xtest,ytest;

				double ldp=dp_gen;
				if ( IsEnergyLoss )
				{
					//Fix me: Do I need to count the energy loss due to windows after target?
					//I think so. So I add them.
					double lE_p=E_p;
					imax=Win_Before_Mag.size();
					for ( i=0; i<imax; i++ )
					{
						lE_p-=Ion_Loss(lE_p,Win_Before_Mag[i]);
						lE_p-=Bremss_Loss(lE_p,Win_Before_Mag[i].bt);
					}
					ldp=(lE_p-P0)/P0;
				}
				for ( i = 0; i < 2; ++i ) {
					q1ex[i]=0;
					dent[i]=0;
					dext[i]=0;
					q3en[i]=0;
					q3ex[i]=0;
				}
				IsPassedQ1Ex=false;
				IsPassedDipoleEn=false;
				IsPassedDipoleEx=false;
				IsPassedQ3En=false;
				IsPassedQ3Ex=false;
					/*Q1 Exit{{{*/
					lPos=p_TCS_ref;
					lMom=p_P_TCS_ref;
					GetRef_Plane(lPos,lMom,Win_Before_Mag,Win_Empty,Q1_Exit_Eff_L,offset);
					matrix[0]=lPos.X()/100;
					matrix[1]=lMom.X()/lMom.Z();
					matrix[2]=lPos.Y()/100;
					matrix[3]=lMom.Y()/lMom.Z();
					matrix[4]=ldp;
					if ( theta>0 )
					{
						xtest=x_e_q1ex__(matrix,msize)*100; //cm
						ytest=y_e_q1ex__(matrix,msize)*100; //cm
						q1ex[0]=xtest/100;
						q1ex[1]=ytest/100;
						if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q1_Radius*Q1_Radius )
						{
							//(xtets!=xtest)==true to avoid nan(Not a number)
							//printf("Blocked by Q1 Exit.\n");
							return 1;
						}
					}
					else
					{
						xtest=x_h_q1ex__(matrix,msize)*100; //cm
						ytest=y_h_q1ex__(matrix,msize)*100; //cm
						q1ex[0]=xtest/100;
						q1ex[1]=ytest/100;
						//I copy this because I'am afraid the condition may be different for these two arms.
						if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q1_Radius*Q1_Radius )
						{
							//(xtets!=xtest)==true to avoid nan(Not a number)
							//printf("Blocked by Q1 Exit.\n");
							return 1;
						}
					}
					IsPassedQ1Ex=true;
					/*}}}*/
					/*Dipole Entrance{{{*/
					lPos=p_TCS_ref;
					lMom=p_P_TCS_ref;
					GetRef_Plane(lPos,lMom,Win_Before_Mag,Win_Empty,D_Entrance_Eff_L,offset);
					matrix[0]=lPos.X()/100;
					matrix[1]=lMom.X()/lMom.Z();
					matrix[2]=lPos.Y()/100;
					matrix[3]=lMom.Y()/lMom.Z();
					matrix[4]=ldp;
					if ( theta>0 )
					{
						xtest=x_e_dent__(matrix,msize)*100; //cm
						ytest=y_e_dent__(matrix,msize)*100; //cm
						dent[0]=xtest/100;
						dent[1]=ytest/100;
						if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
						{
							//printf("Blocked by Dipole Entrance.\n");
							return 1;
						}
					}
					else
					{
						xtest=x_h_dent__(matrix,msize)*100; //cm
						ytest=y_h_dent__(matrix,msize)*100; //cm
						dent[0]=xtest/100;
						dent[1]=ytest/100;
						if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
						{
							//printf("Blocked by Dipole Entrance.\n");
							return 1;
						}
					}
					IsPassedDipoleEn=true;
					/*}}}*/
					/*Dipole Exit{{{*/
					lPos=p_TCS_ref;
					lMom=p_P_TCS_ref;
					GetRef_Plane(lPos,lMom,Win_Before_Mag,Win_Empty,D_Exit_Eff_L,offset);
					matrix[0]=lPos.X()/100;
					matrix[1]=lMom.X()/lMom.Z();
					matrix[2]=lPos.Y()/100;
					matrix[3]=lMom.Y()/lMom.Z();
					matrix[4]=ldp;
					if ( theta>0 )
					{
						xtest=x_e_dext__(matrix,msize)*100; //cm
						ytest=y_e_dext__(matrix,msize)*100; //cm
						dext[0]=xtest/100;
						dext[1]=ytest/100;
						if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
						{
							//printf("Blocked by Dipole Exit.\n");
							return 1;
						}
					}
					else
					{
						xtest=x_h_dext__(matrix,msize)*100; //cm
						ytest=y_h_dext__(matrix,msize)*100; //cm
						dext[0]=xtest/100;
						dext[1]=ytest/100;
						if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
						{
							//printf("Blocked by Dipole Exit.\n");
							return 1;
						}
					}
					IsPassedDipoleEx=true;
					/*}}}*/
					/*Q3 Entrance{{{*/
					lPos=p_TCS_ref;
					lMom=p_P_TCS_ref;
					GetRef_Plane(lPos,lMom,Win_Before_Mag,Win_Empty,Q3_Entrance_Eff_L,offset);
					matrix[0]=lPos.X()/100;
					matrix[1]=lMom.X()/lMom.Z();
					matrix[2]=lPos.Y()/100;
					matrix[3]=lMom.Y()/lMom.Z();
					matrix[4]=ldp;
					if ( theta>0 )
					{
						xtest=x_e_q3en__(matrix,msize)*100; //cm
						ytest=y_e_q3en__(matrix,msize)*100; //cm
						q3en[0]=xtest/100;
						q3en[1]=ytest/100;
						if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q3_Entrance_Radius*Q3_Entrance_Radius )
						{
							//printf("Blocked by Q3 Entrance.\n");
							return 1;
						}
					}
					else
					{
						xtest=x_h_q3en__(matrix,msize)*100; //cm
						ytest=y_h_q3en__(matrix,msize)*100; //cm
						q3en[0]=xtest/100;
						q3en[1]=ytest/100;
						if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>30*30. )
						{
							//printf("Blocked by Q3 Entrance.\n");
							return 1;
						}
					}
					IsPassedQ3En=true;
					/*}}}*/
					/*Q3 Exit{{{*/
					lPos=p_TCS_ref;
					lMom=p_P_TCS_ref;
					GetRef_Plane(lPos,lMom,Win_Before_Mag,Win_Empty,Q3_Exit_Eff_L,offset);
					matrix[0]=lPos.X()/100;
					matrix[1]=lMom.X()/lMom.Z();
					matrix[2]=lPos.Y()/100;
					matrix[3]=lMom.Y()/lMom.Z();
					matrix[4]=ldp;
					if ( theta>0 )
					{
						xtest=x_e_q3ex__(matrix,msize)*100; //cm
						ytest=y_e_q3ex__(matrix,msize)*100; //cm
						q3ex[0]=xtest/100;//cm to m
						q3ex[1]=ytest/100;
						if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q3_Exit_Radius*Q3_Exit_Radius )
						{
							//printf("Blocked by Q3 Exit.\n");
							return 1;
						}
					}
					else
					{
						xtest=x_h_q3ex__(matrix,msize)*100; //cm
						ytest=y_h_q3ex__(matrix,msize)*100; //cm
						q3ex[0]=xtest/100;
						q3ex[1]=ytest/100;
						if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>30*30. )
						{
							//printf("Blocked by Q3 Exit.\n");
							return 1;
						}
					}
					IsPassedQ3Ex=true;
					/*}}}*/

				GetRef_Plane(p_TCS_ref,p_P_TCS_ref,Win_Before_Mag,Win_After_Mag,FP_Eff_L,offset);
				/*}}}*/
			}
			else
				GetRef_Plane(p_TCS_ref,p_P_TCS_ref,Win_Before_Mag,Win_Empty,0,offset,IsEnergyLoss);

			x_tg_ref=p_TCS_ref(0);
			y_tg_ref=p_TCS_ref(1);
			th_tg_ref=p_P_TCS_ref(0)/p_P_TCS_ref(2);
			ph_tg_ref=p_P_TCS_ref(1)/p_P_TCS_ref(2);
			dp_ref=dp_gen;
			E_p=p_P_TCS_ref(3);
			
			
			if ( E_p<ELECTRON_MASS ) {
				printf("[Warning %s: Line %d] E_p=%g<ELECTRON_MASS=%g\n",__FILE__,__LINE__,E_p,ELECTRON_MASS);
				return 1;
			}
			if ( IsEnergyLoss ) {
				Q2=4*E_s*E_p*sinsq;
				q2=-Q2;
				p_TCS_ref(3)=E_p;//MeV
				p_P_TCS_ref(3)=E_p;
				dp_ref=(E_p-P0)/P0;
			}
			cs_M=sigma_M(E_s,Angle_Deg);
			cs_Final=cs_M;
			
			if ( IsRadCor )
			{
				//all nbarn/sr
				cs_qfs=sigma_qfs(Target.Z,Target.A,E_s,E_s-E_p,Angle_Deg,Fermi_Moment,NIE,true);//Q.E. Peak
				cs_del=sigma_del(Target.A,E_s,E_s-E_p,Angle_Deg,Fermi_Moment,DEL_SEP);//Delta
				cs_x=sigma_x(Target.A,E_s,E_s-E_p,Angle_Deg);//DIS
				cs_r1500=sigma_resonance(1500.,Target.A,E_s,E_s-E_p,Angle_Deg,Fermi_Moment);//Resonace 1500 MeV
				cs_r1700=sigma_resonance(1700.,Target.A,E_s,E_s-E_p,Angle_Deg,Fermi_Moment);//Resonace 1700 MeV
				cs_2N=sigma_2N(Target.Z,Target.A,E_s,E_s-E_p,Angle_Deg,Fermi_Moment);//Dip region
				cs_q=cs_qfs+cs_del+cs_x+cs_r1500+cs_r1700+cs_2N;
				//Printf("-------START--------");
				//double tmp=735.0;
				//Printf("sigma_qfs=%e",sigma_qfs(2,3,4238.6,tmp,15.5,130.0,10.0));
				//Printf("sigma_del=%e",sigma_del(3,4238.6,tmp,15.5,130.0,15.0));
				//Printf("sigma_x=%e",sigma_x(3,4238.6,tmp,15.5));
				//Printf("sigma_resonance_1500=%e",sigma_resonance(1500.,3,4238.6,tmp,15.5,130.0));
				//Printf("sigma_resonance_1700=%e",sigma_resonance(1700.,3,4238.6,tmp,15.5,130.0));
				//Printf("sigma_2N=%e",sigma_2N(2,3,4238.6,tmp,15.5,130.0));
				//Printf("cs_q=%e",cs_q);
				//Printf("-------END--------");
				//Printf("Z=%d,A=%g,E_s=%g,E_p=%g,W=%g,theta=%g,NIE=%g,DEL_SEP=%g,Fermi=%g",Target.Z,Target.A,E_s,E_p,E_s-E_p,Angle_Deg,NIE,DEL_SEP,Fermi_Moment);
				//Printf("sigma_qfs=%e",cs_qfs);
				//Printf("sigma_del=%e",cs_del);
				//Printf("sigma_x=%e",cs_x);
				//Printf("sigma_resonance_1500=%e",cs_r1500);
				//Printf("sigma_resonance_1700=%e",cs_r1700);
				//Printf("sigma_2N=%e",cs_2N);
				//Printf("cs_q=%e",cs_q);

				/*Set Info for RadCor{{{*/
				tail_s.SetXYZT(0,0,sqrt(E_s*E_s-M2),E_s);//assume along z
				tail_p.SetXYZT(sqrt(E_p*E_p-M2)*sin(Angle),0,sqrt(E_p*E_p-M2)*cos(Angle),E_p);//assume along z
				tail_t.SetXYZT(0,0,0,Target.M);
				tail_u=tail_s-tail_p+tail_t;

				costheta_s=cos(tail_s.Angle(tail_u.Vect()));//tail_s.Angle(u.Vect()) return rad
				costheta_p=cos(tail_p.Angle(tail_u.Vect()));

				w_s=E_s-E_p/(1-2*E_p*sinsq_Angle/Target.M); //A50
				w_p=E_s/(1+2*E_s*sinsq_Angle/Target.M)-E_p; //A51
				v_s=w_s/E_s; //A53
				v_p=w_p/(E_p+w_p);

				tb=Win_i.TR+Target.TR/2;
				ta=Win_f.TR+Target.TR/2;
				btb=Win_i.bt+Target.bt/2;//b*tb
				bta=Win_f.bt+Target.bt/2;//b*ta

				E_s_min=E_p/(1-2*E_p*sinsq_Angle/Target.M);
				E_p_max=E_s/(1+2*E_s*sinsq_Angle/Target.M);
				/*}}}*/
				cs_qt=CalcuateQuasiTail(cs_q);
				tail=CalcuateTail();
				cs_Final=cs_qt+tail;//not minus,we simulate the real situation
			}

			return 0;
		}
		/*}}}*/

		/*int ToFp(const double& ax,const double& ay,const double& ath,const double& aph,const double& adp){{{*/
		int ToFp(const double& ax,const double& ay,const double& ath,const double& aph,const double& adp)
		{
			//return Number of Event needs to be added. just 1 ^_^
			//cm->meter, ath,aph,adp are correct
			//must be float, otherwise cannot pass to fortran correctly. float<->dimension
			float matrix[MSIZE]={ax/100.,ath,ay/100,aph,adp};
			int msize=MSIZE;
			double xtest,ytest;
			int i;
			if ( !IsMyModel ) {
				for ( i = 0; i < 2; ++i ) {
					q1ex[i]=0;
					dent[i]=0;
					dext[i]=0;
					q3en[i]=0;
					q3ex[i]=0;
				}
				IsPassedQ1Ex=false;
				IsPassedDipoleEn=false;
				IsPassedDipoleEx=false;
				IsPassedQ3En=false;
				IsPassedQ3Ex=false;
			}
			if ( theta>0 )
			{
				//for left arm
				if ( !IsMyModel )
				{
					/*Q1 Exit{{{*/
					xtest=x_e_q1ex__(matrix,msize)*100; //cm
					ytest=y_e_q1ex__(matrix,msize)*100; //cm
					q1ex[0]=xtest/100;
					q1ex[1]=ytest/100;
					if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q1_Radius*Q1_Radius )
					{
						//(xtets!=xtest)==true to avoid nan(Not a number)
						//printf("Blocked by Q1 Exit.\n");
						return 1;
					}
					IsPassedQ1Ex=true;
					/*}}}*/
					/*Dipole Entrance{{{*/
					// Transport electron to dipole entrance, trapezoid define by (jjl)
					// -40cm<x<40cm (+x is down in HRS frame)
					// y=+-(12.5*(1-(1.25*x/840)) (smallest gap at positive x, x in cm)
					xtest=x_e_dent__(matrix,msize)*100; //cm
					ytest=y_e_dent__(matrix,msize)*100; //cm
					dent[0]=xtest/100;
					dent[1]=ytest/100;
					if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
					{
						//printf("Blocked by Dipole Entrance.\n");
						return 1;
					}
					IsPassedDipoleEn=true;
					/*}}}*/
					/*Dipole Exit{{{*/
					xtest=x_e_dext__(matrix,msize)*100; //cm
					ytest=y_e_dext__(matrix,msize)*100; //cm
					dext[0]=xtest/100;
					dext[1]=ytest/100;
					if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
					{
						//printf("Blocked by Dipole Exit.\n");
						return 1;
					}
					IsPassedDipoleEx=true;
					/*}}}*/
					/*Q3 Entrance{{{*/
					xtest=x_e_q3en__(matrix,msize)*100; //cm
					ytest=y_e_q3en__(matrix,msize)*100; //cm
					q3en[0]=xtest/100;
					q3en[1]=ytest/100;
					if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q3_Entrance_Radius*Q3_Entrance_Radius )
					{
						//printf("Blocked by Q3 Entrance.\n");
						return 1;
					}
					IsPassedQ3En=true;
					/*}}}*/
					/*Q3 Exit{{{*/
					xtest=x_e_q3ex__(matrix,msize)*100; //cm
					ytest=y_e_q3ex__(matrix,msize)*100; //cm
					q3ex[0]=xtest/100;
					q3ex[1]=ytest/100;
					if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q3_Exit_Radius*Q3_Exit_Radius )
					{
						//printf("Blocked by Q3 Exit.\n");
						return 1;
					}
					IsPassedQ3Ex=true;
					/*}}}*/
				}
				x_fp=x_e_fp__(matrix,msize)*100.; //cm
				y_fp=y_e_fp__(matrix,msize)*100.; //cm
				th_fp=t_e_fp__(matrix,msize); //(tantheta)
				ph_fp=p_e_fp__(matrix,msize); //(tanphi)
			}
			else
			{
				//for right arm
				if ( !IsMyModel )
				{
					/*Q1 Exit{{{*/
					xtest=x_h_q1ex__(matrix,msize)*100; //cm
					ytest=y_h_q1ex__(matrix,msize)*100; //cm
					q1ex[0]=xtest/100;
					q1ex[1]=ytest/100;
					if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q1_Radius*Q1_Radius )
					{
						//(xtets!=xtest)==true to avoid nan(Not a number)
						//printf("Blocked by Q1 Exit.\n");
						return 1;
					}
					IsPassedQ1Ex=true;
					/*}}}*/
					/*Dipole Entrance{{{*/
					// Transport electron to dipole entrance, trapezoid define by (jjl)
					// -40cm<x<40cm (+x is down in HRS frame)
					// y=+-(D_Y_L*(1-(1.25*x/840)) (smallest gap at positive x, x in cm)
					xtest=x_h_dent__(matrix,msize)*100; //cm
					ytest=y_h_dent__(matrix,msize)*100; //cm
					dent[0]=xtest/100;
					dent[1]=ytest/100;
					if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
					{
						//printf("Blocked by Dipole Entrance.\n");
						return 1;
					}
					IsPassedDipoleEn=true;
					/*}}}*/
					/*Dipole Exit{{{*/
					xtest=x_h_dext__(matrix,msize)*100; //cm
					ytest=y_h_dext__(matrix,msize)*100; //cm
					dext[0]=xtest/100;
					dext[1]=ytest/100;
					if ( (xtest!=xtest)==true || (ytest!=ytest)==true || fabs(xtest)>D_X_Radius || fabs(ytest)>D_Y_L*(1-1.25*xtest/840) ) //nan
					{
						//printf("Blocked by Dipole Exit.\n");
						return 1;
					}
					IsPassedDipoleEx=true;
					/*}}}*/
					/*Q3 Entrance{{{*/
					xtest=x_h_q3en__(matrix,msize)*100; //cm
					ytest=y_h_q3en__(matrix,msize)*100; //cm
					q3en[0]=xtest/100;
					q3en[1]=ytest/100;
					if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q3_Entrance_Radius*Q3_Entrance_Radius )
					{
						//printf("Blocked by Q3 Entrance.\n");
						return 1;
					}
					IsPassedQ3En=true;
					/*}}}*/
					/*Q3 Exit{{{*/
					xtest=x_h_q3ex__(matrix,msize)*100; //cm
					ytest=y_h_q3ex__(matrix,msize)*100; //cm
					q3ex[0]=xtest/100;
					q3ex[1]=ytest/100;
					if ( (xtest!=xtest)==true || (ytest!=ytest)==true || (xtest*xtest+ytest*ytest)>Q3_Exit_Radius*Q3_Exit_Radius )
					{
						//printf("Blocked by Q3 Exit.\n");
						return 1;
					}
					IsPassedQ3Ex=true;
					/*}}}*/
				}
				x_fp=x_h_fp__(matrix,msize)*100.; //cm
				y_fp=y_h_fp__(matrix,msize)*100.; //cm
				th_fp=t_h_fp__(matrix,msize); //(tantheta)
				ph_fp=p_h_fp__(matrix,msize); //(tanphi)
			}
			if ( IsMultiScat && !IsMyModel )
			{
				//Material mixture=GetMixture(Win_After_Mag);
				//mixture.TR*=sqrt(ph_fp*ph_fp+th_fp*th_fp);
				//mixture.L=0;
				//p_FP.SetXYZT(x_fp,y_fp,0,p_TCS_ref(3));
				//p_P_FP.SetXYZT(th_fp,ph_fp,1,p_TCS_ref(3));
				//Transport(p_FP,p_P_FP,mixture,IsMultiScat);//Just change th_fp,ph_fp
				//th_fp=p_P_FP(0)/p_P_FP(2);
				//ph_fp=p_P_FP(1)/p_P_FP(2);


				//Material mixture=GetMixture(Win_After_Mag);
				p_FP.SetXYZT(x_fp,y_fp,0,p_TCS_ref(3));
				p_P_FP.SetXYZT(th_fp,ph_fp,1,p_TCS_ref(3));
				Material Vacuum;
				Vacuum.L=0;
				for ( i = 0; i < Win_After_Mag.size(); ++i ) {
					Vacuum.L-=Win_After_Mag[i].L;
				}
				Vacuum.Name="Vacuum";
				Vacuum.TR=0;
				Transport(p_FP,p_P_FP,Vacuum,IsMultiScat);//Just change th_fp,ph_fp
				for ( i = 0; i < Win_After_Mag.size(); ++i ) {
					Transport(p_FP,p_P_FP,Win_After_Mag[i],IsMultiScat,IsEnergyLoss);//Just change th_fp,ph_fp
				}
				//Transport(p_FP,p_P_FP,mixture,IsMultiScat);//Just change th_fp,ph_fp
				th_fp=p_P_FP(0)/p_P_FP(2);
				ph_fp=p_P_FP(1)/p_P_FP(2);
				x_fp=p_FP(0);
				y_fp=p_FP(1);
			}
			//Smearing
			TRandom* tr=new TRandom();
			x_fp=tr->Gaus(x_fp,VDC_Res_x);
			y_fp=tr->Gaus(y_fp,VDC_Res_y);
			th_fp=tr->Gaus(th_fp,VDC_Res_th/1000.);
			ph_fp=tr->Gaus(ph_fp,VDC_Res_ph/1000.);
			th_fp_no_ort=th_fp;
			delete tr;
			matrix[0]=x_fp/100.;
			msize=1;
			//orthogonalize theta see John.LeRose Webpage
			//http://hallaweb.jlab.org/news/minutes/tranferfuncs.html
			if ( theta>0 )
			{
				th_fp-=ltxfit_(matrix,msize);
			}
			else
				th_fp-=rtxfit_(matrix,msize);
			return 0;
		}
		/*}}}*/

		/*int ReconstructTg(const double& ax,const double& ay,const double& ath,const double& aph,const double& axtg){{{*/
		int ReconstructTg(const double& ax,const double& ay,const double& ath,const double& aph,const double& axtg)
		{
			//return Number of Event needs to be added. just 1 ^_^
			float matrix[MSIZE]={ax/100.,ath,ay/100,aph,axtg/100};
			int msize=MSIZE;

			x_tg_rec=axtg;
			float rf_y,rf_d,rf_th,rf_ph;
			if ( theta>0 )
			{
				//printf("%f %f %f %f %f\n",matrix[0],matrix[1],matrix[2],matrix[3],matrix[4]);
				y_tg_rec=ly00_(matrix,msize)*100;
				th_tg_rec=ltheta_(matrix,msize);
				ph_tg_rec=lphi_(matrix,msize);
				dp_rec=ldelta_(matrix,msize);
				rf_y=y_tg_rec/100;
				rf_d=dp_rec;
				rf_th=th_tg_rec;
				rf_ph=ph_tg_rec;
				rvalue=left_rfunction__(&rf_y,&rf_d,&rf_th,&rf_ph);
			}
			else
			{
				y_tg_rec=ry00_(matrix,msize)*100;
				th_tg_rec=rtheta_(matrix,msize);
				ph_tg_rec=rphi_(matrix,msize);
				dp_rec=rdelta_(matrix,msize);
				rf_y=y_tg_rec/100;
				rf_d=dp_rec;
				rf_th=th_tg_rec;
				rf_ph=ph_tg_rec;
				rvalue=right_rfunction__(&rf_y,&rf_d,&rf_th,&rf_ph);
			}
			TLorentzVector rec(x_tg_rec,y_tg_rec,0,P0*(1+dp_rec));
			
			rec(0)+=D_x;
			rec(1)+=D_y;
			rec(0)+=reactz_TCS*th_tg_rec;
			rec(1)+=reactz_TCS*ph_tg_rec;
			rec(2)+=reactz_TCS;
			rec.RotateX(-theta_rad);
			rec.RotateZ(-PI/2);
			reactz_rec=rec(2);

			//TVector3 yax(0.,1.,0.);
			//TVector3 beam_org(beam_x,beam_y,0);
			//TVector3 beam_ray(0,0,1);
			//Double_t t;
			//TVector3 v;
			//if( !IntersectPlaneWithRay( p_P.Vect(), yax, rec.Vect(), 
			//			beam_org, beam_ray, t, v ))
			//{
			//	return 1; // Oops, track and beam parallel?
			//}
			//else
			//{
			//	reactz_rec=v(2);
			//}

			if (
					//fabs(reactz_rec/100)<0.1 &&
					//fabs(y_tg_rec/100)<=0.01 &&
					rvalue>-1000&&
					fabs(dp_rec)<delta_dp/2 &&
					fabs(th_tg_rec)<delta_th/2 &&
					fabs(ph_tg_rec)<delta_ph/2
					//fabs(ph_tg_rec)<delta_ph/2
			   )
				IsQualified=1;
			else
				IsQualified=0;
			return 0;
		}
		/*}}}*/

		/*Bool_t IntersectPlaneWithRay( const TVector3& xax,{{{*/
		Bool_t IntersectPlaneWithRay( const TVector3& xax,
				const TVector3& yax,
				const TVector3& org,
				const TVector3& ray_start,
				const TVector3& ray_vect,
				Double_t& length,
				TVector3& intersect )
		{
			// Find intersection point of plane (given by 'xax', 'yax', 'org') with
			// ray (given by 'ray_start', 'ray_vect'). 
			// Returns true if intersection found, else false (ray parallel to plane).
			// Output is in 'length' and 'intersect', where
			//   intersect = ray_start + length*ray_vect
			// 'length' and 'intersect' must be provided by the caller.

			// Calculate explicitly for speed.

			Double_t nom[9], den[9];
			nom[0] = den[0] = xax.X();
			nom[3] = den[3] = xax.Y();
			nom[6] = den[6] = xax.Z();
			nom[1] = den[1] = yax.X();
			nom[4] = den[4] = yax.Y();
			nom[7] = den[7] = yax.Z();
			den[2] = -ray_vect.X();
			den[5] = -ray_vect.Y();
			den[8] = -ray_vect.Z();

			Double_t det1 = den[0]*(den[4]*den[8]-den[7]*den[5])
				-den[3]*(den[1]*den[8]-den[7]*den[2])
				+den[6]*(den[1]*den[5]-den[4]*den[2]);
			if( fabs(det1) < 1e-5 )
				return false;

			nom[2] = ray_start.X()-org.X();
			nom[5] = ray_start.Y()-org.Y();
			nom[8] = ray_start.Z()-org.Z();
			Double_t det2 = nom[0]*(nom[4]*nom[8]-nom[7]*nom[5])
				-nom[3]*(nom[1]*nom[8]-nom[7]*nom[2])
				+nom[6]*(nom[1]*nom[5]-nom[4]*nom[2]);

			length = det2/det1;
			intersect = ray_start + length*ray_vect;
			return true;
		}
		/*}}}*/

		/*double CalcuateTail(){{{*/
		double CalcuateTail()
		{
			double result;

			int flag_internal=1,flag_external=1;
			if ( !IsInternal_Bre )
				flag_internal=0;
			if ( !IsExternal_Bre )
				flag_external=0;

			mp_cor=F_soft();
			double lmulti_photon_cor=mp_cor;
			if ( !IsMulti_Photon_Cor )
				lmulti_photon_cor=1;

			cs_p=sigma_p()*lmulti_photon_cor*MEV2SR_TO_NBARNSR;
			cs_ex=sigma_ex()*lmulti_photon_cor*MEV2SR_TO_NBARNSR;
			cs_b=sigma_b(cs_b_bt,cs_b_at)*lmulti_photon_cor*MEV2SR_TO_NBARNSR;
			cs_b_bt *= lmulti_photon_cor*MEV2SR_TO_NBARNSR;
			cs_b_at *= lmulti_photon_cor*MEV2SR_TO_NBARNSR;

			double lcs_int=cs_ex;
			if ( IsAppro )
				lcs_int=cs_p;

			result=(cs_b*flag_external+lcs_int*flag_internal);
			////Fix me: Do not understand here? See radtail_m.f XERR
			//tail_err=1e-35;
			//tail_err=sqrt(tail_err*tail_err+(0.03*lcs_int*lcs_int*flag_internal+0.1*cs_b*cs_b*flag_external)*lmulti_photon_cor);
			return result;
		}
		/*}}}*/

		/*double CalcuateQuasiTail(const double& aCS_NoRad){{{*/
		double CalcuateQuasiTail(const double& aCS_NoRad)
		{
			//make aCS_NoRad to CS_With_QE_Rad
			R=(Target.M+2*E_s*sinsq_Angle)/(Target.M-2*E_p*sinsq_Angle);
			//Phys.Rev.D 12,1884 A82
			double cs_qt_t0=pow(R*DeltaE/E_s,btb+btr)*pow(DeltaE/E_p,bta+btr)*(1-xi(Target.Z,tb,ta)/(DeltaE*(1-bta-btb-2*btr)))*_F(Q2);

			int Num_Of_Steps=80;

			if ( (int)(E_s-R*DeltaE-E_s_min)>Num_Of_Steps )
			{
				Num_Of_Steps=(int)(E_s-R*DeltaE-E_s_min)/2;
				Num_Of_Steps*=2;
			}
			double cs_qt_t1,cs_qt_t2;
			if ( (E_s-R*DeltaE)>=E_s_min )
			{
				cs_qt_t1=Sim_Integral(1,E_s_min,E_s-R*DeltaE,Num_Of_Steps,&SAMCEvent::sigma_qt_Es_Integrand);
			}
			else
				cs_qt_t1=0;
			//Fix me: do we need 3?
			//add 3 to keep from double subtracting the elastic rad tail, see meziani's code line 402-405
			if ( (int)(E_p_max-E_p-DeltaE)>Num_Of_Steps )
			{
				Num_Of_Steps=(int)(E_p_max-E_p-DeltaE)/2;
				Num_Of_Steps*=2;
			}
			if ( (E_p_max)>=(E_p+DeltaE) )
				cs_qt_t2=Sim_Integral(1,E_p+DeltaE,(E_p_max),Num_Of_Steps,&SAMCEvent::sigma_qt_Ep_Integrand);
			else
				cs_qt_t2=0;
			double result;
			result=aCS_NoRad*cs_qt_t0+cs_qt_t1+cs_qt_t2;
			//Printf("E_s=%g,E_p=%g,theta=%g,result=%g",E_s,E_p,Angle_Deg,result);
			if ( result<0 )
				result=0;
			////Fix me: err calculation
			////copy from meziani's code
			////quasi_cs_err=sqrt(quasi_cs_err*quasi_cs_err+0.05*0.05*(cs_qt_t1*cs_qt_t1+cs_qt_t2*cs_qt_t2));
			return result;
		}
		/*}}}*/
		
		


		//For radiative tail
	private:
		/*double sigma_b(double& ocs_b_bt,double& ocs_b_at){{{*/
		double sigma_b(double& ocs_b_bt,double& ocs_b_at)
		{
			//Phys.Rev.D 12,1884 (A49)
			//External tail ( straggling of target, bremsstrahlung and ionization loss )
			//output cross section b before target(oCS_b_bt)
			//output cross section b after target(oCS_b_at)
			double t0,t1,t2;
			t0=(Target.M+2*(E_s-w_s)*sinsq_Angle)/(Target.M-2*E_p*sinsq_Angle);
			t1=_sigma_el(E_s-w_s)*(btb*phi(v_s)/w_s+xi(Target.Z,tb,ta)/(2*w_s*w_s));
			t2=_sigma_el(E_s)*(bta*phi(v_p)/w_p+xi(Target.Z,tb,ta)/(2*w_p*w_p));
			ocs_b_bt=t0*t1;
			ocs_b_at=t2;
			return ocs_b_bt+ocs_b_at;
		}
		/*}}}*/

		/*double _sigma_el(const double& aE){{{*/
		double _sigma_el(const double& aE)
		{
			double lEPRM=aE/(1+2*aE*sinsq_Angle/Target.M);
			double lQ2=4*aE*lEPRM*sinsq_Angle;
			//Phys.Rev.D 12,1884 (A55)
			return _F(lQ2)*sigma_el(aE);
		}
		/*}}}*/

		/*double sigma_el(const double& aE){{{*/
		double sigma_el(const double& aE)
		{
			//Phys.Rev.D 12,1884 (A3)
			double W_1el,W_2el;
			double lEPRM=aE/(1+2*aE*sinsq_Angle/Target.M);
			double lQ2=4*aE*lEPRM*sinsq_Angle;
			Wel(lQ2,Form_Factor_Id,W_1el,W_2el);
			double t0,t1;
			t0=ALPHA*ALPHA*(1-sinsq_Angle)/(4*aE*aE*sinsq_Angle*sinsq_Angle)*lEPRM/aE;
			t1=W_2el+2*sinsq_Angle*W_1el/(1-sinsq_Angle);
			return t0*t1;
		}
		/*}}}*/

		/*void Wel(const double& aQ2,const int& aChoice,double& oW_1el_Q2,double& oW_2el_Q2) const{{{*/
		void Wel(const double& aQ2,const int& aChoice,double& oW_1el_Q2,double& oW_2el_Q2) const
		{
			//Phys.Rev.D 12,1884, it seems q2=Q2>0
			//Form Factor Function, Get W1_el,W2_el
			//aQ2:MeV^2
			//aChoice=1: Proton Form Factors from  Rev.Mod.Phys.41.205 B3, III2-4
			double lQ2_fm;
			double G_E,G_M;
			double tau;
			double F0_e;//charge form factor for He3 and He4
			double F0_m;//magnetic form factor for He3 and He4

			lQ2_fm=aQ2/(HBARC*HBARC);//translate MeV to fm^-1
			double lQ_fm=sqrt(lQ2_fm);

			if ( aChoice<=4 )
			{
				if ( aChoice==1 || aChoice==4 )
				{
					/*Proton or Deuterium{{{*/
					//Phys.Rev.D 12,1884 (A7-A12)
					//double H={1.0007, 1.01807, 1.05584, 0.836380, 0.6864584, 0.672830};
					//here 1e+6 because of GeV->MeV
					//A7
					//here is proton
					//Fix me: only for comparsion with meziani's code
					//G_E=P(aQ2)/pow((1+aQ2/0.71e+6),2);
					G_E=1./pow((1+aQ2/0.71e+6),2);
					G_M=(1+K_P)*G_E;
					if ( aChoice==4 )
					{
						//A12
						//deuterium
						G_E=P(aQ2)/pow((1+aQ2/0.71e+6),2);
						G_M=(1+K_P)*G_E;
						tau=aQ2/4/(PROTON_MASS*PROTON_MASS); //PROTON_MASS
						double G_Mn=K_N*G_E;
						double F_1n=tau*G_Mn/(1+tau);
						double F_1p=(G_E+tau*G_M)/(1+tau);
						double F_2n=G_Mn/(K_N*(1+tau));
						double F2p=(G_M-G_E)/(K_P*(1+tau));
						double G_P=F_1n+F_1p;
						double G_S=G_P+K_N*F_2n+K_P*F2p;
						double F_d=1.580/lQ_fm*(atan(lQ_fm/0.930)-2*atan(lQ_fm/3.19)+atan(lQ_fm/5.45));

						oW_1el_Q2=F_d*F_d*(2/3.)*tau*G_S*G_S;//W1el(Q2)
						oW_2el_Q2=F_d*F_d*G_P*G_P+oW_1el_Q2;//W2el(Q2)


					}
					/*}}}*/
				}
				else if ( aChoice==2 )
				{
					/*He3{{{*/
					//He3 form factors from Phys.Rev.C.15.1396
					//bet fit values for the parameters of He3 charge distribution a,b,c,d,p,q0 for equ(6)+equ(8)
					double par_e[]={0.675,0.366,0.836,-6.78e-3,0.90,3.98};
					//best fit values for the parameters of He3 magnetic distribution a,b,c equ(6)
					double par_m[]={0.654,0.456,0.821};
					//He3 charge form factors from Phys.Rev.C.15.1396, Equ(6)+(8)
					F0_e=exp(-par_e[0]*par_e[0]*lQ2_fm)-par_e[1]*par_e[1]*lQ2_fm*exp(-par_e[2]*par_e[2]*lQ2_fm);
					F0_e += par_e[3]*exp( -pow((sqrt(lQ2_fm)-par_e[5])/par_e[4],2) );
					//He3 magnetic form factors from Phys.Rev.C.15.1396, Equ(6)
					F0_m=exp(-par_m[0]*par_m[0]*lQ2_fm)-par_m[1]*par_m[1]*lQ2_fm*exp(-par_m[2]*par_m[2]*lQ2_fm);
					/*}}}*/
				}
				else if ( aChoice==3 )
				{
					/*He4{{{*/
					//He4 form factors from Phys.Rev.C.15.1396
					//best fit values for the parameters of He4 charge distribution a,b for equ(10)
					double par_e[]={0.316,0.675};//fm
					//He4 charge form factors from Phys.Rev.C.15.1396, Equ(10)
					F0_e=(1-pow(par_e[0]*par_e[0]*lQ2_fm,6))*exp(-par_e[1]*par_e[1]*lQ2_fm);
					//He4 magnetic form factors
					F0_m=0;
					/*}}}*/
				}
				if ( aChoice==2 || aChoice==3 )
				{
					//For He3 or He4
					//Normalize by He_Z
					int He_Z=2;
					G_E=He_Z*F0_e;
					G_M=He_Z*(1+K_HE)*F0_m;
				}

				if ( aChoice!=4 )
				{
					tau=aQ2/4/(Target.M*Target.M); //not PROTON_MASS
					//printf("tau=%g,Target.M=%g,aQ2=%g\n",tau,Target.M,aQ2);
					//Phys.Rev.D 12,1884 (A6)
					oW_1el_Q2=tau*G_M*G_M;//W1el(Q2)
					oW_2el_Q2=(G_E*G_E+tau*G_M*G_M)/(1+tau);//W2el(Q2)
				}
			}
			else
			{
				/*Models > 4{{{*/
				//Next please see Rev.Mod.Phys.28.214
				//page 218, table I
				//F(x)=F(qa), x=qa, q=lQ_fm
				//a=1.07A^1/3
				//a is the root-mean-sqaure radium of the
				//charge distribution
				double lF;
				//Phys.Rev.D 12,1884
				//equ A18
				double la=1.07*pow(Target.A,1/3.);
				double x=lQ_fm*la;
				if ( aChoice==5 ) //Model Point
				{
					lF=1;
				}
				else if ( aChoice==6 ) //Model Uniform
				{
					lF=5*sqrt(5/3.)*pow(x,-3)*(sin(sqrt(5/3.))*x-sqrt(5/3.)*x*cos(sqrt(5/3.))*x);
				}
				else if ( aChoice==7 ) //Model Gaussian
				{
					lF=exp(-x*x/6.);
				}
				else if ( aChoice==8 ) //Model Exponential
				{
					lF=pow((1+x*x/12.),-2);
				}
				else if ( aChoice==9 ) //Model Shell
				{
					lF=1/x*sin(x);
				}
				else if ( aChoice==10 ) //Model Hollow exponential
				{
					lF=(1-x*x/60.)*pow((1+x*x/20.),-3);
				}
				else if ( aChoice==11 ) //Model ...
				{
					lF=(1-x*x/30.)*pow((1+x*x/30.),-4);
				}
				else if ( aChoice==12 ) //Model Yukawa I
				{
					lF=sqrt(2)/(x*tan(x/sqrt(2)));
				}
				else if ( aChoice==13 ) //Model Yukawa II
				{
					lF=1/(1+x*x/6.);
				}
				else if ( aChoice==14 ) //Model Hollow Gaussian
				{
					lF=(1-x*x/15.)*exp(-x*x/10.);
				}
				else if ( aChoice==15 ) //Model Generalized shell model
				{
					double lalpha=(Target.Z-2)/3.;
					double lk=3*(2+5*lalpha)/(2*(2+3*lalpha));
					lF=(1-lalpha*x*x/(2*lk*lk*(2+3*lalpha)))*exp(-x*x/(4*lk*lk));
				}
				else if ( aChoice==16 ) //Model Modified exponential
				{
					lF=pow((1+x*x/18.),-3);
				}
				else if ( aChoice==17 ) //C,Fe,Pb
				{
					//Phys.Rev.D 12,1884
					//equ A16-A18
					double lb=2.4;//fm
					lF=1./(1+x*x/6.)*exp(-lQ_fm*lQ_fm*lb*lb/6.);
				}
				//Phys.Rev.D 12,1884
				//equ A17
				oW_1el_Q2=0;
				oW_2el_Q2=Target.Z*Target.Z*lF*lF;
				/*}}}*/
			}
		}
		/*}}}*/

		/*double P(const double& aQ2) const{{{*/
		double P(const double& aQ2) const
		{
			//only For proton Form Factor Function
			//deviation of the measured form factors from the dipole expression
			//Phys.Rev.D 12,1884 (A8)
			//q2=Q2>0
			double H[]={1.0007, 1.01807, 1.05584, 0.836380, 0.6864584, 0.672830};
			double P=0;
			double prod;
			double lQ=sqrt(aQ2)/1000.; //GeV
			for ( int i=0; i<6; i++ )
			{
				prod=1;
				for ( int j=0; j<6; j++ )
				{
					if ( i!=j )
						prod *= (lQ-j)/(i-j);
				}
				P += H[i]*prod;
			}
			return P;//just a number
		}
		/*}}}*/

		/*double F_soft() const{{{*/
		double F_soft() const
		{
			//Multiple-photon correction
			//Phys.Rev.D 12,1884 (A58)
			return pow(w_s/E_s,btb+btr)*pow(w_p/(E_p+w_p),btb+btr);
		}
		/*}}}*/

		/*double sigma_p(){{{*/
		double sigma_p()
		{
			//Phys.Rev.D 12,1884 (A56)
			//Internal correction ( peak approximation for sigma_ex )
			double t0,t1,t2;
			t0=(Target.M+2*(E_s-w_s)*sinsq_Angle)/(Target.M-2*E_p*sinsq_Angle);
			t1=_sigma_el(E_s-w_s)*(btr*phi(v_s)/w_s);
			t2=_sigma_el(E_s)*(btr*phi(v_p)/w_p);
			return t0*t1+t2;
		}
		/*}}}*/

		/*double sigma_ex(){{{*/
		double sigma_ex()
		{
			//Phys.Rev.D 12,1884 (A24)
			//exact Internal bremsstrahlung
			double t0,t1;
			t0=ALPHA*ALPHA*ALPHA/(2*PI)*E_p/E_s;
			t1=Sim_Integral(0,-1,1,10000,&SAMCEvent::sigma_ex_Integrand);
			return t0*t1;
		}
		/*}}}*/

		/*double sigma_ex_Integrand(const double& aCOSTHK){{{*/
		double sigma_ex_Integrand(const double& aCOSTHK)
		{
			//Phys.Rev.D 12,1884 (A24)
			//aCOSTHK=cos(thetak)
			//s-k-gc=Pf-t, k*k=0 real photon k=(|w|,w);
			double u_0=tail_u(3);
			double abs_u3=tail_u.P();//3 mean vector
			if ( abs_u3<0 )
				cout<<"DEBUG-->abs_u3="<<abs_u3<<endl;

			//costheta_s costheta_p defined in SetTailMemeber()
			//since they are not calculated everytime aCOSTHK changes
			double s_dot_p=tail_s.Dot(tail_p);
			//double u2=sqrt(u(3)*u(3)-u.P()*u.P());
			double u2=tail_u*tail_u;
			double w=0.5*(u2-Target.M*Target.M)/(u_0-abs_u3*aCOSTHK);
			double lq2=2*M2-2*s_dot_p-2*w*(E_s-E_p)+2*w*abs_u3*aCOSTHK;
			double a=w*(E_p-tail_p.P()*costheta_p*aCOSTHK);
			double aPRM=w*(E_s-tail_s.P()*costheta_s*aCOSTHK);//a'
			double bPRM=-w*tail_p.P()*sqrt(1-costheta_p*costheta_p)*sqrt(1-aCOSTHK*aCOSTHK);
			double v=1/(aPRM-a);
			double x=sqrt(a*a-bPRM*bPRM);
			double y=sqrt(aPRM*aPRM-bPRM*bPRM);

			double t0,t1,t2,t3,t4,t5,t6,l_W_2_q2,l_W_1_q2;
			double l_F_q2,lW_1el,lW_2el;
			//if proton,no l_F_q2
			if ( Form_Factor_Id==1 )
				l_F_q2=1;
			else
			{
				l_F_q2=_F(-lq2);
			}
			Wel(-lq2,Form_Factor_Id,lW_1el,lW_2el);
			l_W_2_q2=l_F_q2*lW_2el;
			l_W_1_q2=l_F_q2*lW_1el;
			t0=2*Target.M*w/(lq2*lq2*(u_0-abs_u3*aCOSTHK));
			t1=-a*M2*(2*E_s*(E_p+w)+lq2/2)/pow(x,3)-aPRM*M2*(2*E_p*(E_s-w)+lq2/2)/pow(y,3);
			t2=-2+2*v*(1/x-1/y)*( M2*(s_dot_p-w*w)+s_dot_p*( 2*E_s*E_p-s_dot_p+w*(E_s-E_p) ) );
			t3=( 2*(E_s*E_p+E_s*w+E_p*E_p)+lq2/2-s_dot_p-M2 )/x;
			t4=-( 2*(E_s*E_p-E_p*w+E_s*E_s)+lq2/2-s_dot_p-M2 )/y;
			t5=(a/pow(x,3)+aPRM/pow(y,3))*M2*(2*M2+lq2)+4+4*v*(1/x-1/y)*s_dot_p*(s_dot_p-2*M2);
			t6=(1/x-1/y)*(2*s_dot_p+2*M2-lq2);

			return t0*( l_W_2_q2*(t1+t2+t3+t4)+l_W_1_q2*(t5+t6) );
		}
		/*}}}*/

		//for quasielastic radiative correction
	private:
		/*double sigma_qt_Es_Integrand(const double& aEs){{{*/
		double sigma_qt_Es_Integrand(const double& aEs)
		{
			double t0,t1;
			double dEs=E_s-aEs;
			t0=_sigma_q(aEs,E_p)*pow(dEs/(E_p*R),bta+btr)*pow(dEs/E_s,btb+btr);
			t1=(btb+btr)/dEs*phi(dEs/E_s)+xi(Target.Z,tb,ta)/(2*dEs*dEs);
			return t0*t1;
		}
		/*}}}*/

		/*double sigma_qt_Ep_Integrand(const double& aEp){{{*/
		double sigma_qt_Ep_Integrand(const double& aEp)
		{
			double t0,t1;
			double dEp=aEp-E_p;
			t0=_sigma_q(E_s,aEp)*pow(dEp/aEp,bta+btr)*pow(dEp*R/E_s,btb+btr);
			t1=(bta+btr)/dEp*phi(dEp/aEp)+xi(Target.Z,tb,ta)/(2*dEp*dEp);
			return t0*t1;
		}
		/*}}}*/

		/*double _sigma_q(const double& aEs,const double& aEp){{{*/
		double _sigma_q(const double& aEs,const double& aEp)
		{
			double lQ2=4*aEs*aEp*sinsq_Angle;
			double lsigma_q=sigma_qfs(Target.Z,Target.A,aEs,aEs-aEp,Angle_Deg,Fermi_Moment,NIE,true);
			lsigma_q+=sigma_del(Target.A,aEs,aEs-aEp,Angle_Deg,Fermi_Moment,DEL_SEP);//Delta
			lsigma_q+=sigma_x(Target.A,aEs,aEs-aEp,Angle_Deg);//DIS
			lsigma_q+=sigma_resonance(1500.,Target.A,aEs,aEs-aEp,Angle_Deg,Fermi_Moment);//Resonace 1500 MeV
			lsigma_q+=sigma_resonance(1700.,Target.A,aEs,aEs-aEp,Angle_Deg,Fermi_Moment);//Resonace 1700 MeV
			lsigma_q+=sigma_2N(Target.Z,Target.A,aEs,aEs-aEp,Angle_Deg,Fermi_Moment);//Dip region
			return _F(lQ2)*lsigma_q;
		}
		/*}}}*/

		//For quasi-elastic radiative correction and elastic radiative tail
	private:
		/*double _F(const double& aQ2){{{*/
		double _F(const double& aQ2)
		{
			//Phys.Rev.D 12,1884 (A44)
			//aQ2=q^2, not Q2
			double t0,t1,t2,t3;
			t0=1+0.5772*(bta+btb);
			t1=2*AP*(-14/9.+13*log(aQ2/ELECTRON_MASS/ELECTRON_MASS)/12.);
			t2=-AP*log(E_s/E_p)*log(E_s/E_p)/2;
			t3=AP*(PI*PI/6.-DiLog(1-sinsq_Angle));
			return t0+t1+t2+t3;
		}
		/*}}}*/

		/*double phi(const double& aMu){{{*/
		double phi(const double& aMu)
		{
			//Phys.Rev.D 12,1884 A54
			return (1-aMu+3/4.*aMu*aMu);
		}
		/*}}}*/

		/*double xi(const int& aZ,const double& aTb,const double& aTa){{{*/
		double xi(const int& aZ,const double& aTb,const double& aTa)
		{
			//Phys.Rev.D 12,1884 A52
			return PI*ELECTRON_MASS/2./ALPHA*(aTb+aTa)/((aZ+eta(aZ))*log(183.*pow(aZ,-1/3.)));
		}
		/*}}}*/

		/*double Sim_Integral(const int& aChoice,const double& aMin,const double& aMax,const int& {{{*/
		double Sim_Integral(const int& aChoice,const double& aMin,const double& aMax,const int& aNum_Of_Step,double (SAMCEvent::*aFunc)(const double&))
		{
			//Simpson Integration
			if ( aChoice==0 )
			{
				/*meziani's code{{{*/
				//copy from meziani's code since it's faster.
				double A=aMin;
				double B=aMax;
				double EPS=(aMax - aMin)/aNum_Of_Step ;
				double DA=B-A;
				double FA=(this->*aFunc)(A);
				double FM=4*(this->*aFunc)(.5*(A+B));
				double FB=(this->*aFunc)(B);
				double AREA=1.0;
				double EST=1.0;
				double LEG[20],F2T[20],FMT[20],F3T[20],F4T[20];
				double FBT[20],DXT[20],X1T[20],X2T[20],ART[20];
				double EPST[20],ES2T[20],ES3T[20],SUM1[20],SUM2[20];
				int L=1;
				double DX,  X1,  X2,  F1,  F2,  F3,  F4,  DX6,  EST1,  EST2,  EST3,  SUM;
				bool IsInLoop1=true;

				//Loop1
				while ( IsInLoop1 )
				{
					IsInLoop1=true;
					DX=DA/3;
					X1=A+DX;
					X2=X1+DX;
					F1=4.*(this->*aFunc)(A+.5*DX);
					F2=(this->*aFunc)(X1);
					F3=(this->*aFunc)(X2);
					F4=4.*(this->*aFunc)(A+2.5*DX);
					DX6=DX/6.;
					EST1=(FA+F1+F2)*DX6;
					EST2=(F2+FM+F3)*DX6;
					EST3=(F3+F4+FB)*DX6;
					AREA=AREA-fabs(EST)+fabs(EST1)+fabs(EST2)+fabs(EST3);
					SUM=EST1+EST2+EST3;

					
					if (
							((fabs(EST-SUM)>EPS*AREA||EST==1.0)&&(L>=20))
							|| (fabs(EST-SUM)<=EPS*AREA&&EST!=1.0)
					   )
					{
						do
						{
							if ( LEG[L]>2 )
							{//7
								SUM1[L]=SUM;
								LEG[L]=2;
								A=X1T[L];
								DA=DXT[L];
								FA=F2T[L];
								FM=FMT[L];
								FB=F3T[L];
								AREA=ART[L];
								EST=ES2T[L];
								IsInLoop1=true;
								break;
							}
							else if ( LEG[L]==2 )
							{//8
								SUM2[L]=SUM;
								LEG[L]=1;
								A=X2T[L];
								DA=DXT[L];
								FA=F3T[L];
								FM=F4T[L];
								FB=FBT[L];
								AREA=ART[L];
								EST=ES3T[L];
								EPS=EPST[L];
								IsInLoop1=true;
								break;
							}
							else
							{//9
								SUM=SUM1[L]+SUM2[L]+SUM;
								L=L-1;
								IsInLoop1=false;
							}
						}while(L>1);
					}
					else
					{
						L=L+1;
						LEG[L]=3;
						F2T[L]=F2;
						FMT[L]=FM;
						F3T[L]=F3;
						F4T[L]=F4;
						FBT[L]=FB;
						DXT[L]=DX;
						X1T[L]=X1;
						X2T[L]=X2;
						ART[L]=AREA;
						EPST[L]=EPS/1.7;
						ES2T[L]=EST2;
						ES3T[L]=EST3;

						DA=DX;
						FM=F1;
						FB=F2;
						EST=EST1;
						EPS=EPST[L];
					}
				}

				return SUM;
				/*}}}*/
			}
			else if ( aChoice==1 )
			{
				double Step = (aMax - aMin)/aNum_Of_Step ;
				Step=(aMax-aMin)/aNum_Of_Step;
				double SI = (this->*aFunc)(aMin);
				SI += (this->*aFunc)(aMax);
				double XX = aMin + Step;
				SI += (this->*aFunc)(XX) *4.;
				double AA=aMin;
				double BB;
				for ( int i=2; i<aNum_Of_Step; i+=2 )
				{
					AA=i*Step;
					AA+=aMin;
					BB=AA+Step;
					SI += 2.*(this->*aFunc)(AA) + 4.*(this->*aFunc)(BB);
				}
				return Step * SI/3.0;
			}
			else
			{
				/*my code{{{*/
				double Step = (aMax - aMin)/aNum_Of_Step ;
				double x = aMin ;
				double xPlus = aMin + 0.5*Step ;
				double mean = ((this->*aFunc)(aMin) + (this->*aFunc)(aMax))*0.5 ;
				double sum = (this->*aFunc)(xPlus) ;
				for(int i=1;i<aNum_Of_Step;i++)
				{
					x     += Step ;
					xPlus += Step ;
					mean  += (this->*aFunc)(x) ;
					sum   += (this->*aFunc)(xPlus) ;
				}
				mean += 2.0*sum ;
				return mean*Step/3.0 ;
				/*}}}*/
			}

		}
		/*}}}*/
		



	private:
		/*double Ion_Loss(const double& aE0,const Material& aMaterial){{{*/
		double Ion_Loss(const double& aE0,const Material& aMaterial)
		{
			//aT: g/cm^2, arho: g/cm^3
			//Particle Booklet Equ(27.9)
			//printf("Z=%d,A=%f,T=%f,rho=%f\n",aMaterial.Z,aMaterial.A,aMaterial.T,aMaterial.rho);
			double lK=0.307075/aMaterial.A;// cm^2/g
			double lbetasq=1-ELECTRON_MASS*ELECTRON_MASS/(aE0*aE0);
			double lxi=lK/2*aMaterial.Z/aMaterial.A*aMaterial.T/lbetasq;//aT: g/cm^2
			double lhbarwsq=28.816*28.816*aMaterial.rho*aMaterial.Z/aMaterial.A*1e-12;//MeV arho is density of absorber
			double j=0.200;
			double Delta_p=lxi*(log(2*ELECTRON_MASS*ELECTRON_MASS*lxi/lhbarwsq)+j);
			double lw=4*lxi;
			double result=0;
			TRandom* tr=new TRandom();
			if ( aMaterial.Z!=0 && aMaterial.A!=0 && aMaterial.T!=0 && aMaterial.rho!=0 )
				result=tr->Landau(Delta_p,lw);
			if ( result>(aE0-ELECTRON_MASS) )
				result=aE0-ELECTRON_MASS;
			if ( result<0 )
				result=0;
			delete tr;
			return result;
		}
		/*}}}*/

		/*double Bremss_Loss(const double& aE0,const double& abt){{{*/
		double Bremss_Loss(const double& aE0,const double& abt)
		{
			//Bremsstrahlung Energy Loss for external and internal(equivalent radiator)
			//Xiaodong Jiang, PhD.thesis Equ 5.15
			//http://filburt.mit.edu/oops/Html/Pub/theses/xjiang.ps
			//*0.999 to avoid lose all energy
			double result=0;
			TRandom* tr=new TRandom();
			if ( abt!=0 )
				result=aE0*pow(tr->Rndm()*0.999,1./abt);
			if ( result>(aE0-ELECTRON_MASS) )
				result=aE0-ELECTRON_MASS;
			if ( result<0 )
				result=0;
			delete tr;
			return result;
		}
		/*}}}*/

		/*double eta(const int& aZ){{{*/
		double eta(const int& aZ)
		{
			//Phys.Rev.D 12,1884 A46
			return log(1440*pow(aZ,-2/3.))/log(183*pow(aZ,-1/3.));
		}
		/*}}}*/

		/*double b(const int& aZ){{{*/
		double b(const int& aZ)
		{
			//Phys.Rev.D 12,1884 A45
			if ( aZ!=0 )
				return 4./3.*( 1+1./9.*(aZ+1)/(aZ+eta(aZ))/log(183*pow(aZ,-1/3.)) );
			else
				return 0;
		}
		/*}}}*/

		/*double Rad_Len(const int& aZ,const double& aA){{{*/
		double Rad_Len(const int& aZ,const double& aA)
		{
			//particle book equation 27.20
			//Lrad=elastic form factor F_el, scattering on the nucleus
			//Lrad_prime=inelastic form factor F_inel, scattering on the shell electrons
			//f(Z)=Coulomb correction
			double Lrad,Lrad_prime,f_Z;
			if ( aZ==1 )
			{
				Lrad=5.31;
				Lrad_prime=6.144;
			}
			else if ( aZ==2 )
			{
				Lrad=4.79;
				Lrad_prime=5.621;
			}
			else if ( aZ==3 )
			{
				Lrad=4.74;
				Lrad_prime=5.805;
			}
			else if ( aZ==4 )
			{
				Lrad=4.71;
				Lrad_prime=5.924;
			}
			else
			{
				Lrad=log(184.15*pow(aZ,-1./3.));
				Lrad_prime=log(1194.*pow(aZ,-2./3.));
			}
			double a=ALPHA*aZ;
			a=a*a;
			f_Z=a*(1./(1+a)+0.20206-0.0369*a+0.0083*a*a-0.002*a*a*a);
			if ( aZ!=0 && aA!=0 )
				return 716.408*aA/(aZ*aZ*(Lrad-f_Z)+aZ*Lrad_prime);
			else
				return 0;
		}
		/*}}}*/

		/*void Transport(TLorentzVector& aPos,TLorentzVector& aMom,const Material& aMaterial,const bool& aIsMultiScatt=false){{{*/
		void Transport(TLorentzVector& aPos,TLorentzVector& aMom,const Material& aMaterial,const bool& aIsMultiScatt=false,const bool& aIsEnergyloss=false)
		{
			//need aMaterial.TR and aMaterial.L
			//aPos and aMom are inputs, also outputs
			//aPos: position aMom: momentum
			double lE=aMom(3);
			double ms_phi,ms_theta;
			if ( aIsMultiScatt && fabs(aMaterial.TR)>1e-15 )
			{
				ms_phi=MultiScattering(lE,aMaterial.TR);//rad
			}
			else
				ms_phi=0;
			TRandom* tr=new TRandom();
			ms_theta=PI/2*tr->Rndm();
			delete tr;
			double lth=aMom.X()/aMom.Z();
			double lph=aMom.Y()/aMom.Z();

			//pass L/2
			aPos.SetX(aPos.X()+aMaterial.L/2*lth);
			aPos.SetY(aPos.Y()+aMaterial.L/2*lph);
			aPos.SetZ(aPos.Z()+aMaterial.L/2);

			//change the angle and pass the rest L/2
			lth+=ms_phi*cos(ms_theta);
			lph+=ms_phi*sin(ms_theta);
			aMom.SetX(lth);
			aMom.SetY(lph);
			aMom.SetZ(1);
			aPos.SetX(aPos.X()+aMaterial.L/2*lth);
			aPos.SetY(aPos.Y()+aMaterial.L/2*lph);
			aPos.SetZ(aPos.Z()+aMaterial.L/2);
			if ( aIsEnergyloss ) {
				lE-=Ion_Loss(lE,aMaterial);
				lE-=Bremss_Loss(lE,aMaterial.bt);
				aMom(3)=lE;
			}
		}
		/*}}}*/

		/*double MultiScattering(const double& aE,const double& aTR){{{*/
		double MultiScattering(const double& aE,const double& aTR)
		{
			//only for electron
			double lPsq=aE*aE-ELECTRON_MASS*ELECTRON_MASS;
			double bcp=lPsq/aE;
			double ltheta0=13.6/bcp*sqrt(aTR)*(1+0.038*log(aTR));
			if ( aTR!=0 )
			{
				//if ( aE>800 )
				//	return gRandom->Gaus(0,ltheta0/2);//rad
				//else if ( aE>550 )
				//	return gRandom->Gaus(0,ltheta0);//rad
				//else
				//	return gRandom->Gaus(0,ltheta0*1.5);//rad
				double result=gRandom->Gaus(0,ltheta0/2);
				return result;//rad
			}
			else
				return 0;
		}
		/*}}}*/

		/*void SetMaterial(Material& aMaterial){{{*/
		void SetMaterial(Material& aMaterial)
		{
			aMaterial.M=aMaterial.A*AMU;//MeV
			aMaterial.X0=Rad_Len(aMaterial.Z,aMaterial.A);
			if ( aMaterial.X0!=0 )
				aMaterial.TR=aMaterial.T/aMaterial.X0;
			else
				aMaterial.TR=0;
			aMaterial.bt=b(aMaterial.Z)*aMaterial.TR;
		}
		/*}}}*/

		/*Material GetMixture(const vector<Material>& aWin){{{*/
		Material GetMixture(const vector<Material>& aWin)
		{
			size_t i;
			size_t imax=aWin.size();
			Material mixture;
			mixture.Name="mixture";
			mixture.L=0;
			mixture.A=0;
			//get mixture TR
			for ( i=0; i<imax; i++ )
			{
				mixture.L+=aWin[i].L;
				mixture.A+=aWin[i].A;
			}
			mixture.TR=0;
			for ( i=0; i<imax; i++ )
			{
				mixture.TR+=aWin[i].A/mixture.A*aWin[i].TR;
			}
			return mixture;
		}
		/*}}}*/

		/*void GetRef_Plane(TLorentzVector& aoPos,TLorentzVector& aoMom,const vector<Material>& aWinBefore,const vector<Material>& aWinAfter,const double& aL,const double& aOffset){{{*/
		void GetRef_Plane(TLorentzVector& aoPos,TLorentzVector& aoMom,const vector<Material>& aWinBefore,const vector<Material>& aWinAfter,const double& aL,const double& aOffset,const bool& aIsEnergyloss=false)
		{
			//Get Refinement aoPos and aoMom for each plane on q1,q2,d,q3,fp
			//ao means input and output
			//aL=Vacuum Length, aOffset=distance between interaction point and z=0 in TCS
			vector<Material> AllWins;
			Material Vacuum;

			AllWins.clear();
			AllWins=aWinBefore;
			Vacuum.Name="Vacuum";
			Vacuum.L=aL;
			Vacuum.A=0;
			Vacuum.TR=0;
			AllWins.push_back(Vacuum);
			size_t i;
			for ( i=0; i<aWinAfter.size(); i++ )
				AllWins.push_back(aWinAfter[i]);
			Vacuum.L=0;
			for ( i = 0; i < AllWins.size(); ++i ) {
				Transport(aoPos,aoMom,AllWins[i],IsMultiScat,aIsEnergyloss);//move to mag
				Vacuum.L-=AllWins[i].L;
				if ( fabs(aoPos(2)-110.)<1 ) {
					x_collimator=aoPos(0);
					y_collimator=aoPos(1);
				}
			}
			Vacuum.L+=aOffset;
			Transport(aoPos,aoMom,Vacuum);//back to ztg=0 in TCS
		}
		/*}}}*/

		/*double sigma_qfs(const int& aZ,const double& aA,const double& aE,const double& aomega,const double& aTheta,const double& aPF,const double& aEPS,const bool& aIsQ2dep=false){{{*/
		double sigma_qfs(const int& aZ,const double& aA,const double& aE,const double& aomega,const double& aTheta,const double& aPF,const double& aEPS,const bool& aIsQ2dep=false)
		{
			//Q.E. Peak
			//aE=MeV,aTheta=deg
			//aPF: fermi momentum
			//aEPS: nucleon interaction energy
			double GAMR=120; //no idea
			double PFR=230; //no idea
			double QMSRQ=4*730.*(730-115.)*pow(sin(37.1*DegToRad()/2),2); //no idea
			double QVSRQ=QMSRQ+115*115; //no idea
			double AP0=840; //no idea
			double AP1=750; //no idea
			double AP;
			int NA=(int)aA;
			if ( NA==1 )
				AP=AP0;
			else if ( NA<4 )
				AP=AP0+(aA-1)*(AP1-AP0)/3;
			else
				AP=AP1;
			double lth=aTheta*DegToRad();
			double lQ2=4*aE*(aE-aomega)*sin(lth/2)*sin(lth/2);
			double lq3mSQ=lQ2+aomega*aomega;

			double sigma_s=sigma_M(aE,aTheta)*Recoil(aE,aTheta,PROTON_MASS);
			double sigma_p=pow(GEp(lQ2,AP),2)+tau(lQ2)*pow(GMp(lQ2,AP),2);
			sigma_p /= (1+tau(lQ2));
			sigma_p += 2*tau(lQ2)*pow(GMp(lQ2,AP),2)*pow(tan(lth/2),2);
			sigma_p *= sigma_s;
			double sigma_n=pow(GEn(lQ2,AP),2)+tau(lQ2)*pow(GMn(lQ2,AP),2);
			sigma_n /= (1+tau(lQ2));
			sigma_n += 2*tau(lQ2)*pow(GMn(lQ2,AP),2)*pow(tan(lth/2),2);
			sigma_n *= sigma_s;

			double lEp=4*aE*aE*sin(lth/2)*sin(lth/2)/(1+aE*(1-cos(lth))/PROTON_MASS);//Q2
			lEp=lEp/(2*PROTON_MASS)+aEPS;//omega
			lEp=aE-lEp;

			double ARG,DEN,GAMQ,result;
			if ( NA==1 )
			{
				ARG=(aE-aomega-lEp)/sqrt(2.);
				DEN=2.51;
			}
			else
			{
				GAMQ=GAMR*aPF/PFR*sqrt(lq3mSQ)/sqrt(QVSRQ);
				ARG=(aE-aomega-lEp)/(1.2*GAMQ/2);
				DEN=2.13*(GAMQ/2);
			}
			int NQ=int(ARG);
			if ( abs(NQ)>10 )
				result=0;
			else
				result=(aZ*sigma_p+(aA-aZ)*sigma_n)*exp(-ARG*ARG)/DEN;
			if ( aIsQ2dep )
			{
				double lpar[2];
				lpar[0]=0.686704; //this is only work for carbon
				lpar[1]=0.0453828; //this is only work for carbon
				if ( Target.Z==26 ) {//for iron
					lpar[0]=1; //just satisfy the data
					lpar[1]=0.0453828; //just satisfy the data
				}
				double lCor=lpar[1]*1e6/lQ2+lpar[0];
				result/=lCor;
			}
			return result;
		}
		/*}}}*/

		/*double sigma_del(const double& aA,const double& aE,const double& aomega,const double& aTheta,const double& aPF,const double& aEPSD){{{*/
		double sigma_del(const double& aA,const double& aE,const double& aomega,const double& aTheta,const double& aPF,const double& aEPSD)
		{
			//Delta
			//aE=MeV,aTheta=deg
			//aPF: fermi momentum
			double AP0=840; //no idea
			double AP1=750; //no idea

			double GAMR=120; //no idea
			double PFR=230; //no idea
			double QMSR=4.*730.*(730.-390.)*pow(sin(37.1*DegToRad()/2.),2);
			double QVSR=QMSR+390.*390;
			double QMSRQ=4*730.*(730-115.)*pow(sin(37.1*DegToRad()/2),2); //no idea
			double QVSRQ=QMSRQ+115*115; //no idea
			double QFDP=1.02e-7;//no idea
			double GAMSPRD=140;
			double GAMDP=110;
			double GSPRDA;
			double AD0=774;
			double AD1=700;
			double AD;
			int NA=(int)aA;
			if ( NA==1 )
			{
				GSPRDA=0;
				AD=AD0;
			}
			else if ( NA<4 )
			{
				GSPRDA=(aA-1)*GAMSPRD/3.;
				AD=AD0+(aA-1)*(AD1-AD0)/3.;
			}
			else
			{
				GSPRDA=GAMSPRD;
				AD=AD1;
			}
			double lth=aTheta*DegToRad();
			double lQ2=4*aE*(aE-aomega)*sin(lth/2)*sin(lth/2);
			double lq3mSQ=lQ2+aomega*aomega;
			double result;
			double GAMQ;
			if ( NA>1 )
				GAMQ=GAMR*aPF*sqrt(lq3mSQ)/PFR/sqrt(QVSRQ);
			else
				GAMQ=0.;
			double GAM=sqrt(GAMDP*GAMDP+GAMQ*GAMQ+GSPRDA*GSPRDA);
			result=QFDP*(GAMDP/GAM);
			double EKAPPA=aomega-lQ2/2./PROTON_MASS;
			double CMTOT2=PROTON_MASS*PROTON_MASS+2*PROTON_MASS*EKAPPA;
			result*=CMTOT2*GAM*GAM;
			double DM=1219;
			result/=(pow((CMTOT2-(DM+aEPSD)*(DM+aEPSD)),2)+CMTOT2*GAM*GAM);
			result*=pow(FD(lQ2,AD),2)/pow(FD(QMSR,AD),2);
			result*=lq3mSQ/QVSR;
			result*=(lQ2/2./lq3mSQ+tan(lth/2.)*tan(lth/2));
			result/=(QMSR/2./QVSR+pow(tan(37.1*DegToRad()/2.),2));
			result*=sigma_M(aE,aTheta)/sigma_M(730,37.1);
			double WTHRESH = 4.*aE*aE*sin(lth/2.)*sin(lth/2.)+PI_MASS*PI_MASS+2.*PI_MASS*PROTON_MASS;
			WTHRESH /= 2.*PROTON_MASS;
			double THRESHD = 1.+aPF/PROTON_MASS+aPF*aPF/2./PROTON_MASS/PROTON_MASS+2.*aE*sin(lth/2.)*sin(lth/2.)/PROTON_MASS;
			WTHRESH = WTHRESH/THRESHD;
			double THRESH;
			double GAMPI=5;
			if ( aomega>WTHRESH )
				THRESH=1.-exp(-(aomega-WTHRESH)/GAMPI);
			else
				THRESH=0.;
			result*=aA*THRESH;

			return result*1e7; //output is nb/sr
		}
		/*}}}*/

		/*double sigma_x(const double& aA,const double& aE,const double& aomega,const double& aTheta){{{*/
		double sigma_x(const double& aA,const double& aE,const double& aomega,const double& aTheta)
		{
			//DIS
			//aE=MeV,aTheta=deg
			//aPF: fermi momentum

			double SIG0=100.e-4;
			double SIG1=54.*1.e-1;
			double GAM0=650.;
			double lth=aTheta*PI/180.;

			double SIGEE=0;
			if ( aomega>1e-5 )
			{
				double QMS=4.*aE*(aE-aomega)*pow(sin(lth/2.),2);
				double ARG0=aomega-QMS/2./PROTON_MASS-PI_MASS-PI_MASS*PI_MASS/2./PROTON_MASS;
				double ARG1=ARG0/GAM0;
				double ARG=ARG1*ARG1/2.;
				double SHAPE;
				if ( ARG1>8 )
					SHAPE=1.+SIG1/SIG0/ARG0;
				else if ( ARG1<1E-5)
					SHAPE=0;
				else if ( ARG1<0.1 )
					SHAPE=SIG1*ARG0/2./GAM0/GAM0/SIG0;
				else
					SHAPE=(1.-exp(-ARG))*(1.+SIG1/SIG0/ARG0);
				double EKAPPA=aomega-QMS/2./PROTON_MASS;
				double SIGGAM=SIG0*SHAPE;
				double QS=QMS+aomega*aomega;
				double EPS=1./(1.+2.*QS*tan(lth/2.)*tan(lth/2)/QMS);
				double FLUX=ALPHA*EKAPPA*(aE-aomega)/2./PI/PI/QMS/aE/(1.-EPS);
				if ( FLUX<1E-20 )
					FLUX=0;
				SIGEE=FLUX*SIGGAM*FPHENOM(QMS)*FPHENOM(QMS);
				double lR=0.56*1.E6/(QMS+PROTON_MASS*PROTON_MASS);
				double FACTOR1=1.+EPS*lR;
				SIGEE*=FACTOR1;
			}

			return aA*SIGEE*1e7;// nbarn/sr
		}
		/*}}}*/

		/*double sigma_2N(const int& aZ,const double& aA,const double& aE,const double& aomega,const double& aTheta,const double& aPF){{{*/
		double sigma_2N(const int& aZ,const double& aA,const double& aE,const double& aomega,const double& aTheta,const double& aPF)
		{
			//Dip region
			double lth=aTheta*DegToRad();
			double DM=1232.;
			double A2=550.;
			double PFR=60.;
			double GAM2N=20.;
			double GAMQFR=40.;
			double GAMREF=300.;
			double GAMR=GAMREF;
			double SIGREF=0.20e-7;
			double QMSR=4.*596.8*(596.8-380.)*pow(sin(60.*DegToRad()/2.),2);
			double QVSR=QMSR+380.*380;
			double SIGKIN=0.5*sigma_M(596.8e0,60.);
			SIGKIN=SIGKIN*(QMSR/2./QVSR+pow(tan(60.*DegToRad()/2.),2));
			SIGKIN=SIGKIN*QVSR*pow(FD(QMSR,A2),2);
			SIGKIN=SIGKIN*GAMR/GAMREF;
			double SIGCON=SIGREF/SIGKIN;
			double QMS=4.*aE*(aE-aomega)*pow(sin(lth/2.),2);
			double QVS=QMS+aomega*aomega;
			double GAMQF=GAMQFR*(aPF/PFR)*(sqrt(QVS)/sqrt(QVSR));
			double EFFMASS=(PROTON_MASS+DM)/2.;
			double SIG=(aZ*(aA-aZ)/aA)*sigma_M(aE,aTheta);
			SIG=SIG*(QMS/2./QVS+pow(tan(lth/2.),2));
			SIG=SIG*QVS*pow(FD(QMS,A2),2);
			double EKAPPA=aomega-QMS/2./PROTON_MASS;
			double CMTOT2=PROTON_MASS*PROTON_MASS+2.*PROTON_MASS*EKAPPA;
			//GAM=sqrt(GAMR*GAMR+GAMQF*GAMQF);
			double GAM=GAMR;
			SIG=SIG*CMTOT2*GAM*GAM;
			SIG=SIG/(pow((CMTOT2-EFFMASS*EFFMASS),2)+CMTOT2*GAM*GAM);
			SIG=SIG*(GAMR/GAM)*SIGCON;
			double result=SIG;
			double WTHRESH=QMS/4./PROTON_MASS;
			double THRESH;
			if ( aomega>WTHRESH )
				THRESH=1.-exp(-(aomega-WTHRESH)/GAM2N);
			else
				THRESH=0.;
			result=result*THRESH;
			return result*1e7; // nbarn/sr
		}
		/*}}}*/

		/*double sigma_resonance(const double& aRes_Energy,const double& aA,const double& aE,const double& aomega,const double& aTheta,const double& aPF){{{*/
		double sigma_resonance(const double& aRes_Energy,const double& aA,const double& aE,const double& aomega,const double& aTheta,const double& aPF)
		{
			//resonace
			double lth=aTheta*DegToRad();
			double PFR=230.;
			double EPSR=0.;
			double AR0;
			double AR1;
			double QFRP;
			double VPHOTONE;
			if ( fabs(aRes_Energy-1500)<1e-3 )
			{
				AR0=1000.;
				AR1=1000.;
				QFRP=1.20e-7;
				VPHOTONE=1240;
			}
			else if ( fabs(aRes_Energy-1700)<1e-3 )
			{
				AR0=1200.;
				AR1=1200.;
				QFRP=0.68e-7;
				VPHOTONE=1520;
			}
			double GAMQFR=120.;
			double GAMSPRD=140.;
			double GAMR=110.;
			double GAMPI=5.;
			double QMSQFR=4.*730.*(730.-115.)*pow(sin(37.1*DegToRad()/2.),2);
			double QVSQFR=QMSQFR+115.*115;
			double QMSRR=4.*10000.*(10000.-VPHOTONE)*pow(sin(6.*DegToRad()/2.),2);
			double QVSRR=QMSRR+VPHOTONE*VPHOTONE;
			double SIGREF=pow(FD(QMSRR,AR0),2)*QVSRR;
			SIGREF=SIGREF*(QMSRR/2./QVSRR+pow(tan(6.*DegToRad()/2.),2));
			SIGREF=SIGREF*sigma_M(10000.e0,6.);
			int NA=int(aA);
			double QFR,GSPRDA,AR;
			if ( NA==1 )
			{
				QFR=QFRP;
				GSPRDA=0.;
				AR=AR0;
			}
			else if ( NA<4 )
			{
				QFR=QFRP;
				GSPRDA=(aA-1.)*GAMSPRD/3.;
				AR=AR0+(aA-1.)*(AR1-AR0)/3.;
			}
			else
			{
				AR=AR1;
				GSPRDA=GAMSPRD;
				QFR=QFRP;
			}
			double QMS=4.*aE*(aE-aomega)*pow(sin(lth/2.),2);
			double QVS=QMS+aomega*aomega;
			double GAMQ;
			if ( NA>1 )
				GAMQ=GAMQFR*aPF*sqrt(QVS)/PFR/sqrt(QVSQFR);
			else
				GAMQ=0.;
			double CMTOT2=PROTON_MASS*PROTON_MASS+2.*PROTON_MASS*aomega-QMS;
			double WTHRESH=4.*aE*aE*pow(sin(lth/2.),2)+PI_MASS*PI_MASS+2.*PI_MASS*PROTON_MASS;
			WTHRESH=WTHRESH/2./PROTON_MASS;
			double THRESHD=1.+aPF/PROTON_MASS+aPF*aPF/2./PROTON_MASS/PROTON_MASS+2.*aE*pow(sin(lth/2.),2)/PROTON_MASS;
			WTHRESH=WTHRESH/THRESHD;
			double THRESH;
			if ( aomega>WTHRESH)
				THRESH=1.-exp(-(aomega-WTHRESH)/GAMPI);
			else
				THRESH=0.;
			double EPR=aE-(aRes_Energy-PROTON_MASS)*(aRes_Energy+PROTON_MASS)/2./PROTON_MASS;
			EPR=EPR/(1.+2.*aE*pow(sin(lth/2.),2)/PROTON_MASS);
			EPR=EPR-EPSR;
			double WR=aE-EPR;
			double GAM=sqrt(GAMR*GAMR+GAMQ*GAMQ+GSPRDA*GSPRDA);
			double result;
			result=QFR*(GAMR/GAM)/SIGREF;
			result*=CMTOT2*GAM*GAM;
			result /= (pow((CMTOT2-(aRes_Energy+EPSR)*(aRes_Energy+EPSR)),2)+CMTOT2*GAM*GAM);
			result*=QVS*pow(FD(QMS,AR),2);
			result*=(QMS/2./QVS+pow(tan(lth/2.),2));
			result*=sigma_M(aE,aTheta);
			result*=aA*THRESH;
			return result*1e7; // nbarn/sr
		}
		/*}}}*/
		


		/*double FPHENOM(const double& aQ2){{{*/
		double FPHENOM(const double& aQ2)
		{
			double A1=.55;
			double A2=20./1.E6;
			double B1=.45;
			double B2=.45/1.E6;
			double C1=0.03;
			double C2=0.2/1.E12;
			double result=A1*exp(-A2*aQ2)+B1*exp(-B2*aQ2);
			result+=C1*exp(-C2*pow((aQ2-4.5E6),2));
			result=sqrt(result);
			return result;
		}
		/*}}}*/

		/*double FD(const double& aQ2,const double aA){{{*/
		double FD(const double& aQ2,const double aA)
		{
			return 1/pow((1+aQ2/(aA*aA)),2);
		}
		/*}}}*/

		/*double sigma_M(const double& aE,const double& aTheta){{{*/
		double sigma_M(const double& aE,const double& aTheta)
		{
			//Mott Cross Section
			//aE=MeV,aTheta=deg
			double ltheta=aTheta/2.*PI/180;
			double mott;
			if ( ltheta!=0 )
				mott=pow(ALPHA*cos(ltheta)/(2*aE*pow(sin(ltheta),2)),2);
			else
				mott=0;
			return mott*MEV2SR_TO_NBARNSR; //nbarn
			//return mott; //nbarn
		}
		/*}}}*/

		/*double Recoil(const double& aE,const double& aTheta,const double& aM){{{*/
		double Recoil(const double& aE,const double& aTheta,const double& aM)
		{
			return 1/(1+aE*(1-cos(aTheta*DegToRad()))/aM);
		}
		/*}}}*/

		/*double GEp(const double& aQ2,const double aPar){{{*/
		double GEp(const double& aQ2,const double aPar)
		{
			return 1/pow((1+aQ2/(aPar*aPar)),2);
		}
		/*}}}*/

		/*double GEn(const double& aQ2,const double aPar){{{*/
		double GEn(const double& aQ2,const double aPar)
		{
			double result = -MU_N;
			result = result * tau(aQ2)/( 1.0+5.6*tau(aQ2) );
			result = result * GEp(aQ2,aPar);
			return result;
		}
		/*}}}*/

		/*double GMp(const double& aQ2,const double aPar){{{*/
		double GMp(const double& aQ2,const double aPar)
		{
			return MU_P * GEp(aQ2,aPar);
		}
		/*}}}*/

		/*double GMn(const double& aQ2,const double aPar){{{*/
		double GMn(const double& aQ2,const double aPar)
		{
			return MU_N * GEp(aQ2,aPar);
		}
		/*}}}*/

		/*double tau(const double& aQ2){{{*/
		double tau(const double& aQ2)
		{
			return aQ2/4/(PROTON_MASS*PROTON_MASS); //PROTON_MASS
		}
		/*}}}*/

	public:
		//Not safe, but it's ok for physicsists
		/*Member Data from file{{{*/
		int Id; //Event Id
		int IsPassed; //if pass through the magnet, 1=true, 0=false(not pass)
		int IsQualified; //if rec var is in the range of gen, 1=true, 0=false(not qualified)
		double theta; //scattering angle(deg)
		Material Target; //target material
		Material Win_i,Win_f;//front[inital] window and back[final] window, stick with target
		vector<Material> Win_Before_Mag;//windows after target and before magnetic
		vector<Material> Win_After_Mag;//windows after target and after magnetic
		double T_theta; //target angle(deg), the central line of target is along beam
		//----
		//|  |
		//|  |--->central line, so the T_theta is defined like phi_tg
		//|  |
		//----
		/*}}}*/

		/*only for RadCor{{{*/
		double cs_qfs;//quasi-free scattering cross section
		double cs_del;//Delta
		double cs_x;  //DIS
		double cs_r1500; //Resonace 1500 MeV
		double cs_r1700; //Resonace 1700 MeV
		double cs_2N; //Dip region
		double cs_q;//sum of all the cross section above


		TLorentzVector tail_s;
		TLorentzVector tail_p;
		TLorentzVector tail_t;
		TLorentzVector tail_u;
		double costheta_s;
		double costheta_p;
		double w_s,w_p;
		double v_s,v_p;
		double tb,ta;
		double btb;//b*tb
		double bta;//b*ta
		double cs_p; //sigma_p (peak approximation for internal radation)
		double cs_b_bt; //the part before target in sigma_b
		double cs_b_at; //the part after target in sigma_b
		double cs_b; //sigma_b (external bremsstrahlung ionization)
		double cs_ex; //sigma_ex (extract internal radation)
		double mp_cor; //Multiple Photon Correction
		double tail; //radiative tail
		double tail_err; //error for radiative tail


		double R; //Phys.Rev.D 12 1884 A83
		double E_s_min; //E min from E_p MeV
		double E_p_max; //Ep max from E_s MeV
		double cs_qt;
		/*}}}*/

		/*Member Data derived from variables above{{{*/
		//TLorentzVector ( vec, t ) = ( k,E ) in ROOT
		//HCS: Hall Coordinate System
		//TCS: Target Coordinate System
		TLorentzVector s; //s(s3,E_s) 4-momentum position of the incident electron in HCS
		TLorentzVector target_center_HCS; // 4-momentum position of the target center in HCS
		TLorentzVector target_center_TCS; // 4-momentum position of the target center in TCS
		TLorentzVector target_edgepoint_HCS; //4-momentum position of the edge point in HCS (if theta>0,T_H/2,if <0,-T_H/2)
		TLorentzVector target_edgepoint_TCS; //4-momentum position of the edge point in TCS
		TLorentzVector s_TCS; //s_TCS(s3,E_s) 4-momentum position of the incident electron in TCS, and ztg=0
		TLorentzVector p_TCS; //p_TCS(p3,E_p) 4-momentum position of the outgoing electron in TCS, and ztg=0
		TLorentzVector p_inter_point_TCS; //p_inter_point_TCS(p3,E_p) 4-momentum position of the outgoing electron in TCS, and ztg=reactz_TCS
		TLorentzVector p_P; //p_P 4-momentum momentum of the outgoing electron in HCS
		TLorentzVector p_P_TCS; //p_P_TCS 4-momentum momentum of the outgoing electron in TCS

		double Angle; //real scattering angel, not from file,rad
		double Angle_Deg; //real scattering angel, not from file,Deg
		double sinsq_Angle;//sin(Angle/2)^2

		double theta_rad;//scattering angle in SAMC
		double sinsq;//sin(theta/2)^2
		double reactz_TCS;//
		double btr;//b*tr (equivalent radiator)

		double Q2;// MeV*MeV
		double q2;//q2=-Q2 MeV*MeV
		/*}}}*/

		/*Member Data from function{{{*/
		/*}}}*/

		/*Member Data from random{{{*/
		double beam_x; //cm
		double beam_y; //cm
		double reactz_gen; // interaction z position cm
		double E_s;//=incident beam energy (MeV)
		double E_p;//=incident beam energy (MeV)
		//gen=generator
		double th_tg_gen;//theta target
		double ph_tg_gen;//phi target
		double dp_gen;//dp at target
		/*}}}*/

		/*Member Data derived from variables above{{{*/
		double x_tg_gen; //
		double y_tg_gen; //
		/*}}}*/

		/*Member Data for refinement{{{*/
		double dp_ref; //dp for John.LeRose matrix
		double x_tg_ref; //x tg for John.LeRose matrix
		double y_tg_ref; //y tg for John.LeRose matrix
		double th_tg_ref; //th tg for John.LeRose matrix
		double ph_tg_ref; //ph tg for John.LeRose matrix
		TLorentzVector p_TCS_ref;//
		TLorentzVector p_P_TCS_ref;//
		double x_collimator;
		double y_collimator;
		/*}}}*/

		/*Member Data for focal plane {{{*/
		TLorentzVector p_FP; //p_FP 4-momentum momentum of the outgoing electron in focal plane
		TLorentzVector p_P_FP; //p_P_FP 4-momentum momentum of the outgoing electron in focal plane
		double x_fp; //
		double y_fp; //
		double th_fp; //
		double ph_fp; //
		double th_fp_no_ort; //th_fp without TXFIT orthogonalization
		double q1ex[2];//0:x 1:y
		double dent[2];
		double dext[2];
		double q3en[2];
		double q3ex[2];
		bool IsPassedQ1Ex;
		bool IsPassedDipoleEn;
		bool IsPassedDipoleEx;
		bool IsPassedQ3En;
		bool IsPassedQ3Ex;
		/*}}}*/

		/*Member Data for reconstruct target varibles {{{*/
		double x_tg_rec; //cm
		double y_tg_rec; //cm
		double th_tg_rec; //
		double ph_tg_rec; //
		double dp_rec;
		double reactz_rec;
		double rvalue;
		double cs_M;//mott cross section
		double cs_Final;//final cross section
		//=cs_M if IsRadCor=0
		//=cs_q-tail-quasi-cor; if IsRadCor=1
		/*}}}*/

		/*void Init(){{{*/
		void Init()
		{
			IsPassed=1;
			IsQualified=1;
			Win_Before_Mag.clear();
			Win_After_Mag.clear();
		};
		/*}}}*/

		/*void Print(){{{*/
		void Print()
		{
			printf("--------------------\n");
			//refer to http://www.fileformat.info/info/unicode/char/search.htm for \x
			printf("%*s=(%*.2f,%*.2f,%*.4f,%*.2f) (4-momentum position of target center [HCS])\n",20,"target_center_HCS(x,y,z,E_s)",10,target_center_HCS(0),10,target_center_HCS(1),10,target_center_HCS(2),10,target_center_HCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.4f,%*.2f) (4-momentum position of target center [TCS])\n",20,"target_center_TCS(x,y,z,E_s)",10,target_center_TCS(0),10,target_center_TCS(1),10,target_center_TCS(2),10,target_center_TCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.4f,%*.2f) (4-momentum position of target edgepoint [HCS])\n",20,"target_edgepoint_HCS(x,y,z,E_s)",10,target_edgepoint_HCS(0),10,target_edgepoint_HCS(1),10,target_edgepoint_HCS(2),10,target_edgepoint_HCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.4f,%*.2f) (4-momentum position of target edgepoint [TCS])\n",20,"target_edgepoint_TCS(x,y,z,E_s)",10,target_edgepoint_TCS(0),10,target_edgepoint_TCS(1),10,target_edgepoint_TCS(2),10,target_edgepoint_TCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of incident electron [HCS])\n",20,"s(x,y,z,E_s)",10,s(0),10,s(1),10,s(2),10,s(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of incident electron [TCS],ztg=0)\n",20,"s_TCS(x,y,z,E_s)",10,s_TCS(0),10,s_TCS(1),10,s_TCS(2),10,s_TCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of outgoing electron [TCS], ztg=0)\n",20,"p_TCS(x,y,z,E_s)",10,p_TCS(0),10,p_TCS(1),10,p_TCS(2),10,p_TCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of outgoing electron [TCS], ztg=reactz_TCS=%.2f)\n",20,"p_inter_point_TCS(x,y,z,E_s)",10,p_inter_point_TCS(0),10,p_inter_point_TCS(1),10,p_inter_point_TCS(2),10,p_inter_point_TCS(3),reactz_TCS);
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [HCS])\n",20,"p_P(x,y,z,E_s)",10,p_P(0),10,p_P(1),10,p_P(2),10,p_P(3));
			printf("%*s=(%*.2e,%*.2e,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [TCS])\n",20,"p_P_TCS(x,y,z,E_s)",10,p_P_TCS(0),10,p_P_TCS(1),10,p_P_TCS(2),10,p_P_TCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of outgoing electron [TCS], ztg=0 after refinement)\n",20,"p_TCS_ref(x,y,z,E_s)",10,p_TCS_ref(0),10,p_TCS_ref(1),10,p_TCS_ref(2),10,p_TCS_ref(3));
			printf("%*s=(%*.2e,%*.2e,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [TCS] after refinement)\n",20,"p_P_TCS_ref(x,y,z,E_s)",10,p_P_TCS_ref(0),10,p_P_TCS_ref(1),10,p_P_TCS_ref(2),10,p_P_TCS_ref(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [FP])\n",20,"p_P(x,y,z,E_s)",10,p_P(0),10,p_P(1),10,p_P(2),10,p_P(3));
			printf("%*s=(%*.2e,%*.2e,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [FP])\n",20,"p_P_FP(x,y,z,E_s)",10,p_P_FP(0),10,p_P_FP(1),10,p_P_FP(2),10,p_P_FP(3));
			printf("%-*s=%*d %-*s %-*s\n",  15,"Id",        10,Id,        8,"",              40,"(Event Id)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"E_s",       10,E_s,       8,"MeV",           40,"(Es Incident Energy)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"E_p",       10,E_p,       8,"MeV",           40,"(Ep Scattering Energy)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"theta",     10,theta,     9,"\xc2\xb0",      40,"(\xce\xb8 scattering angle deg)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"theta_rad", 10,theta_rad, 8,"",              40,"(\xce\xb8 scattering angle rad)");
			printf("%-*s=%*s %-*s %-*s\n",  15,"Target",         10,Target.Name.c_str(),  8,"",              40,"(Target Name)");
			printf("%-*s=%*d %-*s %-*s\n",  15,"Z",         10,Target.Z,  8,"",              40,"(Atomic Number of target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"A",         10,Target.A,  8,"g/mol",         40,"(Mass Number of target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"T",         10,Target.T,         9,"g/cm\xc2\xb2",  40,"(Thickness of target)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"TR",        10,Target.TR,        8,"rad_len",       40,"(Thickness of target in rad_len)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"rho",       10,Target.rho,       9,"g/cm\xc2\xb3",  40,"(\xf0\x9d\x9c\x8c density of target)");
			printf("%-*s=%*s %-*s %-*s\n",  15,"Win_i",         10,Win_i.Name.c_str(),  8,"",              40,"(Win_i Name)");
			printf("%-*s=%*d %-*s %-*s\n",  15,"Z_i",       10,Win_i.Z,       8,"",              40,"(Atomic Number of front[initial] window)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"A_i",       10,Win_i.A,       8,"g/mol",         40,"(Mass Number of front[initial] window)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"T_i",       10,Win_i.T,       9,"g/cm\xc2\xb2",  40,"(Thickness of front[initial] window)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"TR_i",      10,Win_i.TR,      8,"rad_len",       40,"(Thickness of front[initial] window in rad_len)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"rho_i",     10,Win_i.rho,     9,"g/cm\xc2\xb3",  40,"(\xf0\x9d\x9c\x8ci density of front[initial] window)");
			printf("%-*s=%*s %-*s %-*s\n",  15,"Win_f",         10,Win_f.Name.c_str(),  8,"",              40,"(Win_f Name)");
			printf("%-*s=%*d %-*s %-*s\n",  15,"Z_f",       10,Win_f.Z,       8,"",              40,"(Atomic Number of back[final] window)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"A_f",       10,Win_f.A,       8,"g/mol",         40,"(Mass Number of back[final] window)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"T_f",       10,Win_f.T,       9,"g/cm\xc2\xb2",  40,"(Thickness of back[final] window)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"TR_f",      10,Win_f.TR,      8,"rad_len",       40,"(Thickness of back[final] window in rad_len)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"rho_f",     10,Win_f.rho,     9,"g/cm\xc2\xb3",  40,"(\xf0\x9d\x9c\x8ci density of back[final] window)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"T_theta",   10,T_theta,   9,"\xc2\xb0",      40,"(\xce\xb8t target angle)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"sinsq",     10,sinsq,     8,"",              40,"(sin\xc2\xb2(\xce\xb8/2))");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"M_T",       10,Target.M,  8,"MeV",           40,"(Mass of target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"beam_x",    10,beam_x,    8,"cm",            40,"(Beam X)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"beam_y",    10,beam_y,    8,"cm",            40,"(Beam Y)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"reactz_gen",    10,reactz_gen,    8,"cm",            40,"(Z react)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"x_tg_gen",      10,x_tg_gen,      8,"cm",            40,"(x target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"y_tg_gen",      10,y_tg_gen,      8,"cm",            40,"(y target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"th_tg_gen",     10,th_tg_gen,     8,"",              40,"(tan(\xce\xb8tg))");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"ph_tg_gen",     10,ph_tg_gen,     8,"",              40,"(tan(\xcf\x86tg))");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"dp_gen",        10,dp_gen*100,    8,"%",             40,"(dp at target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"th_tg_ref", 10,th_tg_ref, 8,"",              40,"(tan(\xce\xb8tg))[ref_tg]");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"ph_tg_ref", 10,ph_tg_ref, 8,"",              40,"(tan(\xcf\x86tg))[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"dp_ref",    10,dp_ref*100,8,"%",             40,"(dp at target)[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"x_tg_ref",  10,x_tg_ref,  8,"cm",            40,"(x target)[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"y_tg_ref",  10,y_tg_ref,  8,"cm",            40,"(y target)[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"x_fp",  10,x_fp,  8,"cm",            40,"(x target)[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"y_fp",  10,y_fp,  8,"cm",            40,"(y target)[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"th_fp", 10,th_fp, 8,"",              40,"(tan(\xce\xb8tg))[ref_tg]");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"ph_fp", 10,ph_fp, 8,"",              40,"(tan(\xcf\x86tg))[ref_tg]");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"x_tg_rec",      10,x_tg_rec,      8,"cm",            40,"(x target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"y_tg_rec",      10,y_tg_rec,      8,"cm",            40,"(y target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"th_tg_rec",     10,th_tg_rec,     8,"",              40,"(tan(\xce\xb8tg))");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"ph_tg_rec",     10,ph_tg_rec,     8,"",              40,"(tan(\xcf\x86tg))");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"dp_rec",        10,dp_rec*100,    8,"%",             40,"(dp at target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"rvalue",        10,rvalue,    8,"",             40,"(rvalue of target reconstructed variables)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Angle",     10,Angle,     8,"rad",           40,"(real scattering angle[rad])");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Angle",     10,Angle*RadToDeg(),     9,"\xc2\xb0",           40,"(real scattering angle[deg])");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"Q2",        10,Q2,        9,"MeV\xc2\xb2",   40,"(Q\xc2\xb2)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"q2",        10,q2,        9,"MeV\xc2\xb2",   40,"(-Q\xc2\xb2)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"btr",       10,btr,       8,"rad_len",       40,"(b*tr equivalent radiator unit in rad_len)");
			size_t i;
			printf("\t\tWindows List Before Magnetic\n");
			for ( i=0; i<Win_Before_Mag.size(); i++ )
			{
				PrintMaterial(Win_Before_Mag[i]);
			}
			printf("\t\tWindows List After Magnetic\n");
			for ( i=0; i<Win_After_Mag.size(); i++ )
			{
				PrintMaterial(Win_After_Mag[i]);
			}

		};
		/*}}}*/
		/*void PrintMaterial(){{{*/
		void PrintMaterial(const Material& aMaterial)
		{
			printf("%-*s=%*d %-*s %-*s\n",15,Form("%s_Z",aMaterial.Name.c_str()),     10,aMaterial.Z,       8,"",  40,"(Atomic Number)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,Form("%s_A",aMaterial.Name.c_str()),     10,aMaterial.A,       8,"g/mol",  40,"(Atomic Weight)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,Form("%s_L",aMaterial.Name.c_str()),     10,aMaterial.L,       8,"cm",  40,"(length)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,Form("%s_T",aMaterial.Name.c_str()),     10,aMaterial.T,       9,"g/cm\xc2\xb2",  40,"(thickness)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,Form("%s_TR",aMaterial.Name.c_str()),     10,aMaterial.TR,       8,"rad_len",  40,"(thickness in rad_len)");
			printf("%-*s  =%*.2e %-*s %-*s\n",16,Form("%s_\xf0\x9d\x9c\x8c",aMaterial.Name.c_str()),     10,aMaterial.rho,       9,"g/cm\xc2\xb3",  40,"(density)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,Form("%s_X0",aMaterial.Name.c_str()),     10,aMaterial.X0,       9,"g/cm\xc2\xb2",  40,"(Radiation Length)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,Form("%s_bt",aMaterial.Name.c_str()),     10,aMaterial.bt,       8,"",  40,"(bt)");
		}
		;/*}}}*/
		
}
;/*}}}*/
