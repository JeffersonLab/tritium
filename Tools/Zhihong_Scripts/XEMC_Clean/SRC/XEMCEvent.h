////////////////////////////////////////////////////////////
//      XEMC Radiation Cross Section Model                //
//        -- Zhihong Ye, 06/14/2012                       //
//Note:                                                   //
//  The born cross section model is from XEMC_Born.h,that //
//  I converted from old Fortran xem package into C++.    //
//  The radiation part is extract from Huan Yao's C++     //
//  version of SAMC code, where I cleaned up the event    //
//  generation and info related to HRS forward and back   //
//  -ward reconstruction, etc.                            //
////////////////////////////////////////////////////////////
#include "XEMC_Const.h"
#include "XEMC_Born.h"

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
};
/*}}}*/

/*class XEMCEvent{{{*/
class XEMCEvent
{
	public:
		/*XEMCEvent(){{{*/
		XEMCEvent()
		{
		}
		/*}}}*/

		/*virtual ~XEMCEvent(){{{*/
		virtual ~XEMCEvent()
		{
			Win_Before_Mag.clear();
			Win_After_Mag.clear();
		}
		/*}}}*/

		XEMCEvent(XEMCEvent const&){};
		XEMCEvent& operator=(XEMCEvent const&){};

		/*void Init(){{{*/
		void Init()
		{
			Win_Before_Mag.clear();
			Win_After_Mag.clear();
		};
		/*}}}*/

		/*int Run(){{{*/
		int Run()
		{   
			int Num_Event_Add=0;
			SetTarget();//Define Target Materials, like Windowns and Chambers
			Num_Event_Add=CalcXS();//Calculate Cross Sections

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
		/*void SetTarget(){{{*/
		void SetTarget()
		{
			//set value for Member Data derived from variables from file
			//File provides E_s,theta,Target.(Z,A,T,rho) Win_i.(Z,A,T,rho)
			//Win_f.(Z,A,T,rho) T_theta
			//Win_Before_Mag(Name,Z,A,L,rho,X0)
			//Win_After_Mag(Name,Z,A,L,rho,X0)
			//beam_x,beam_y,reactz_gen,th_tg_gen,ph_tg_gen,dp_gen
			//z0,HRS_L,VDC_Res(x,y,th,ph),D_(x,y),T_L,P0
			/*Set Material.(Z,A,M,X0,T,TR,bt){{{*/
			SetMaterial(Target);
			SetMaterial(Win_i);
			SetMaterial(Win_f);
			/*}}}*/

			/*Set Beam Info(HCS and TCS){{{*/
			//know beam_x,beam_y,reactz_gen,E_s,theta,HRS_L
			//Set s,s_TCS,x_tg_gen,y_tg_gen,p_TCS,p_P,p_P_TCS
			theta_rad=theta*DegToRad();//rad
			s(0)=beam_x; //beam offset, cm
			s(1)=beam_y; //beam offset, cm
			reactz_gen+=-(beam_x)*tan(T_theta*DegToRad());
			s(2)=reactz_gen;  //--Z. Ye, 09/11/2012
			//s(2)=0.0;    //Assume reaction always happens at the center of the target, replace reactz_gen; //cm
			target_edgepoint_TRCS(0)=theta/fabs(theta)*T_H/2;//if theta>0,T_H/2, if<0, -T_H/2 in Target Rotation Coordinate System(T_theta=0) not TCS, check Coordinate.svg
			target_edgepoint_TRCS(1)=0;
			target_edgepoint_TRCS(2)=T_L/2;//and z0=0

			TLorentzVector lp_TRCS;//the interaction point in TRCS at T_theta=0
			lp_TRCS=s;
			lp_TRCS(2)-=z0;
			//Printf("s(%g,%g,%g),target_edgepoint_TRCS(%g,%g,%g),lp_TRCS(%g,%g,%g)",s(0),s(1),s(2),target_edgepoint_TRCS(0),target_edgepoint_TRCS(1),target_edgepoint_TRCS(2),lp_TRCS(0),lp_TRCS(1),lp_TRCS(2));
			lp_TRCS.RotateY(-T_theta*DegToRad());

			s_TCS=s;
			s_TCS.RotateZ(PI/2);//passive ratation around HCS, so -(-PI/2)
			s_TCS.RotateX(theta_rad);//passive ratation around HCS, so -(-theta_rad)
			s_TCS(0)-=D_x;
			s_TCS(1)-=D_y;

			p_inter_point_TCS=s_TCS;
			reactz_TCS=s_TCS.Z();

			s_TCS(0)-=s_TCS.Z()*theta_rad; //replace th_tg_gen;
			s_TCS(1)-=s_TCS.Z()*0.0;       //replace ph_tg_gen;
			s_TCS(2)-=s_TCS.Z();

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
			p_P.RotateY(-atan(theta_rad)); //replace (-atan(th_tg_gen));
			p_P.RotateX(atan(0.0));//atan(ph_tg_gen)

			Angle=theta_rad;
			Angle_Deg = theta;
			//     Angle_Deg=Angle*RadToDeg();
			sinsq_Angle=sin(Angle/2)*sin(Angle/2);
			sinsq=sin(Angle/2)*sin(Angle/2);//sin(Angle/2)^2
			//p_P_TCS=lz;//Now think it's in TCS
			//p_P_TCS.RotateY(atan(th_tg_gen));//p_P_TCS.X()/p_P_TCS.Z()=th_tg
			//p_P_TCS.RotateX(-atan(ph_tg_gen));//p_P_TCS.Y()/p_P_TCS.Z()=ph_tg
			p_P_TCS.SetX(theta_rad); //replace th_tg_gen
			p_P_TCS.SetY(0.0);//replace ph_tg_gen
			p_P_TCS.SetZ(1);
			p_P_TCS.SetVect(p_P_TCS.Vect().Unit());

			p_TCS(3)=E_p;//MeV
			p_P(3)=E_p;
			p_P_TCS(3)=E_p;

			/*end of Set Beam}}}*/

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
			size_t i;
			size_t imax;
			Material m0;
			vector<Material>::iterator it=Win_Before_Mag.begin();//here only have the material before mag defined in input file.
			double lHL=0; //if lHL<HRS_Lcm or <3.57m(See Hall_A_NIM), i assume it's air
			unsigned int FirstInWinBeforeMagBlock=0;
			imax=Win_Before_Mag.size();

			for ( i=0; i<imax; i++ ){
				lHL+=Win_Before_Mag[i].L;
				//cerr<< Form("Win#%d, L=%f",i,Win_Before_Mag[i].L)<<endl;
			}
			double lphrad=theta_rad-T_theta*DegToRad();//scattering angle in y axis in TCS or x axis in TRCS
			//it doesn't I add atan(ph_tg_gen) since it's so small for target and windows

			//Win_Before_Mag add Win_f
			if ( Win_f.Z!=0 && Win_f.A!=0 && Win_f.rho!=0 && Win_f.TR!=0 )
			{
				m0=Win_f;//Win_f added to Win_Before_Mag
				m0.Name="Win_f";
				Win_Before_Mag.insert(it,m0);
				FirstInWinBeforeMagBlock++;
			}

			//Add special treatment for E08014 long targets, which have temperature distribution along targets
			//Upstream part from -10cm to -2cm are denser and Downstream from -2cm to +10cm.
			//Basically we treat this system like two targets, each of which has different target lenght and 
			//density, while all others remain the same. So in the caluclation, only m0.L, m0.rho,m0.TR,and m0.bt
			//are need to be updated for each part.	--Z. Ye, 09/11/2012
			/*if(IsBump){{{*/
			if(IsBump){
				//cerr<<"--- Doing special treatment for long targets! "<<endl;
				//FIX_HERE, now we calculate the actual length from the reaction point 
				double lL; //assume target is rentangle lL=distance between reactz_gen and edge of target
				//in TRCS x plane<--> y plane in TCS
				double Target_Boundary = -2.0 - z0; //Currently I take the boundary of two parts at -2cm //FIX_HERE, a potential bug since I am not sure the sign of z0, Z.YE9/19/2012
				double rho_Up = 1.225;
				double rho_Down = 0.85;
				lL=target_edgepoint_TRCS(2)-lp_TRCS(2);//verticle distance between interaction point and edge plane of target
				if ( lL>T_L ) {
					Printf("target_edgepoint_TRCS(%g,%g,%g),lp_TRCS(%g,%g,%g)",target_edgepoint_TRCS(0),target_edgepoint_TRCS(1),target_edgepoint_TRCS(2),lp_TRCS(0),lp_TRCS(1),lp_TRCS(2));
				}
				/*If the reaction point locates at the upstream part,which is densier, we need to treat the target as two independent parts{{{*/
				if(lL>(target_edgepoint_TRCS(2)-Target_Boundary)){
					lL*=tan(lphrad);
					lL+=lp_TRCS(0);//x on edge

					/*in the target and come out from the downstream edge of the target{{{*/
					if ( fabs(lL)<fabs(target_edgepoint_TRCS(0)) ) {//HCS 0=x top view
						//Win_Before_Mag add Target
						it=Win_Before_Mag.begin();

						//Add Upstream Part
						m0=Target;	
						m0.L=T_L - (target_edgepoint_TRCS(2)-Target_Boundary);    //FIX_HERE, from -10cm to -2cm, target is denser but final number need to conform
						m0.rho*=rho_Up; //FIX_HERE, final density are need to conform. 
						m0.T=m0.L*m0.rho;
						m0.bt=b(m0.Z)*m0.T/m0.X0;
						lL=fabs((lp_TRCS(2)-Target_Boundary)/cos(lphrad));//The length from the boundary to the reaction points in the dense area.
						m0.L=lL;
						Win_Before_Mag.insert(it,m0);
						FirstInWinBeforeMagBlock++;

						Target_Up_TR = fabs(lp_TRCS(2)+T_L/2.0+z0)*m0.rho/m0.X0;
						Target_Up_bt =b(m0.Z)*Target_Up_TR;
						Target_Down_TR = m0.L*m0.rho/m0.X0;
						Target_Down_bt = b(m0.Z)*Target_Down_TR;
						it++;
						//Add Downstream Part
						m0=Target;
						m0.L=target_edgepoint_TRCS(2)-Target_Boundary;   //FIX_HERE, need to conform
						m0.rho*=rho_Down;//FIX_HERE, need to conform
						m0.T=m0.L*m0.rho;
						m0.bt=b(m0.Z)*m0.T/m0.X0;
						lL=fabs((target_edgepoint_TRCS(2)-Target_Boundary)/cos(lphrad));//The length from the boundary to the reaction points in the dense area.
						m0.L=lL;
						Win_Before_Mag.insert(it,m0);
						FirstInWinBeforeMagBlock++;

						Target_Down_TR += m0.L*m0.rho/m0.X0;
						Target_Down_bt += b(m0.Z)*Target_Down_TR;
					}/*}}}*/
					/*out of target before hiting the edge of target, the edge means the downstream face{{{*/
					else {
						//Add Upstream Part
						m0=Target;	
						m0.L=T_L - (target_edgepoint_TRCS(2)-Target_Boundary);    //FIX_HERE, from -10cm to -2cm, target is denser but final number need to conform
						m0.rho*=rho_Up; //FIX_HERE, final density are need to conform. 
						m0.T=m0.L*m0.rho;
						m0.bt=b(m0.Z)*m0.T/m0.X0;

						lL=fabs((lp_TRCS(2)-Target_Boundary)/cos(lphrad));//The length from the boundary to the reaction points in the dense area.
						/*If the particle goes out before goes into the downstream part{{{*/
						if(lL*sin(lphrad)>fabs((target_edgepoint_TRCS(0)-lp_TRCS(0))/sin(lphrad))){
							lL=fabs((target_edgepoint_TRCS(0)-lp_TRCS(0))/sin(lphrad));
							m0.L=lL;	
							Win_Before_Mag.insert(it,m0);
							FirstInWinBeforeMagBlock++;

							Target_Up_TR = fabs(lp_TRCS(2)+T_L/2.0+z0)*m0.rho/m0.X0;
							Target_Up_bt =b(m0.Z)*Target_Up_TR;
							Target_Down_TR = m0.L*m0.rho/m0.X0;
							Target_Down_bt = b(m0.Z)*Target_Down_TR;
						}/*}}}*/
						/*If the particle travel through the downstream part, add both{{{*/
						else{	
							m0.L=lL;	
							Win_Before_Mag.insert(it,m0);
							FirstInWinBeforeMagBlock++;

							Target_Up_TR = fabs(lp_TRCS(2)+T_L/2.0+z0)*m0.rho/m0.X0;
							Target_Up_bt =b(m0.Z)*Target_Up_TR;
							Target_Down_TR = m0.L*m0.rho/m0.X0;
							Target_Down_bt = b(m0.Z)*Target_Down_TR;

							it++;
							//Add Downstream Part
							m0=Target;
							m0.L=target_edgepoint_TRCS(2)-Target_Boundary;   //FIX_HERE, need to conform
							m0.rho*=rho_Down;//FIX_HERE, need to conform
							m0.T=m0.L*m0.rho;
							m0.bt=b(m0.Z)*m0.T/m0.X0;
							lL=fabs((target_edgepoint_TRCS(0)-lp_TRCS(0))/sin(lphrad))-fabs((lp_TRCS(2)-Target_Boundary)/cos(lphrad));
							if(lL>=0.0){//If the particle travel through downstream part before It get out, add this part, otherwise only consider the upstream part
								m0.L=lL;
								Win_Before_Mag.insert(it,m0);					
								FirstInWinBeforeMagBlock++;
								Target_Down_TR += m0.L*m0.rho/m0.X0;
								Target_Down_bt += b(m0.Z)*Target_Down_TR;
							}
						}/*}}}*/
					}/*}}}*/
				}/*}}}*/
				/*If the reaction happens at the downstream part, we only need to consider this part only, unless we have angle larger than 90deg, which we don't have{{{*/	
				else{
					//Add Upstream part
					m0=Target;
					m0.rho*=rho_Up;
					Target_Up_TR = fabs(Target_Boundary+T_L/2.0+z0)*m0.rho/m0.X0;
					Target_Up_bt =b(m0.Z)*Target_Up_TR;

					//Add Downstream Part
					m0=Target;
					m0.L=target_edgepoint_TRCS(2)-Target_Boundary;   //FIX_HERE, need to conform
					m0.rho*=rho_Down;//FIX_HERE, need to conform
					m0.T=m0.L*m0.rho;
					m0.bt=b(m0.Z)*m0.T/m0.X0;
					//FIX_HERE, now we calculate the actual length from the reaction point 
					//in TRCS x plane<--> y plane in TCS
					lL=target_edgepoint_TRCS(2)-lp_TRCS(2);//verticle distance between interaction point and edge plane of target
					if ( lL>T_L ) {
						Printf("target_edgepoint_TRCS(%g,%g,%g),lp_TRCS(%g,%g,%g)",target_edgepoint_TRCS(0),target_edgepoint_TRCS(1),target_edgepoint_TRCS(2),lp_TRCS(0),lp_TRCS(1),lp_TRCS(2));
					}
					lL*=tan(lphrad);
					lL+=lp_TRCS(0);//x on edge

					if ( fabs(lL)<fabs(target_edgepoint_TRCS(0)) ) {//HCS 0=x top view
						//in the target
						lL=fabs((target_edgepoint_TRCS(2)-lp_TRCS(2))/cos(lphrad));
					}
					else {
						//out of target before hiting the edge of target, the edge means the downstream face.
						lL=fabs((target_edgepoint_TRCS(0)-lp_TRCS(0))/sin(lphrad));
					}
					m0.L=lL;
					Win_Before_Mag.insert(it,m0);
					FirstInWinBeforeMagBlock++;
					Target_Up_TR += fabs(lp_TRCS(2)-Target_Boundary)*m0.rho/m0.X0;
					Target_Up_bt += b(m0.Z)*Target_Up_TR;
					Target_Down_TR = m0.L*m0.rho/m0.X0;
					Target_Down_bt = b(m0.Z)*Target_Down_TR;
				}/*End of downstream part only}}}*/

			//	cerr<<Form("L=%f, TR_Up=%f, bt_Up=%f, TR_Down=%f, bt_Down=%f", lL,Target_Up_TR, Target_Up_bt, Target_Down_TR, Target_Down_bt)<<endl;
			}/*End of special long target treatment}}}*/
			//No Special Treatment
			else{	
				//Win_Before_Mag add Target
				it=Win_Before_Mag.begin();
				m0=Target;//target after interaction point added to Win_Before_Mag
				double lL; //assume target is rentangle lL=distance between reactz_gen and edge of target
				//in TRCS x plane<--> y plane in TCS
				lL=target_edgepoint_TRCS(2)-lp_TRCS(2);//verticle distance between interaction point and edge plane of target

				if ( lL>T_L ) {
					Printf("target_edgepoint_TRCS(%g,%g,%g),lp_TRCS(%g,%g,%g)",target_edgepoint_TRCS(0),target_edgepoint_TRCS(1),target_edgepoint_TRCS(2),lp_TRCS(0),lp_TRCS(1),lp_TRCS(2));
				}
				lL*=tan(lphrad);
				lL+=lp_TRCS(0);//x on edge

				if ( fabs(lL)<fabs(target_edgepoint_TRCS(0)) ) {//HCS 0=x top view
					//in the target
					lL=fabs((target_edgepoint_TRCS(2)-lp_TRCS(2))/cos(lphrad));
				}
				else {
					//out of target before hiting the edge of target, the edge means the downstream face.
					lL=fabs((target_edgepoint_TRCS(0)-lp_TRCS(0))/sin(lphrad));
				}

				m0.L=lL;
				Win_Before_Mag.insert(it,m0);
				FirstInWinBeforeMagBlock++;
			};

			//Win_Before_Mag add Air before last material in the input file
			imax=Win_Before_Mag.size();
			if ( (lHL)<(HRS_L) )
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
			else if ( (lHL)>HRS_L )
			{
				printf("[Error %s: Line %d] Total windows length=%f>HRS_L=%f.\n",__FILE__,__LINE__,(lHL+Win_Before_Mag[imax-1].L),HRS_L);
				exit(-3);
			}
			//Correct First Material.L in Win_Before_Mag Block of inputfile
			it=Win_Before_Mag.begin();
			(it+FirstInWinBeforeMagBlock)->L+=reactz_TCS;

			for ( i=0; i<FirstInWinBeforeMagBlock; i++ )
				(it+FirstInWinBeforeMagBlock)->L-=(it+i)->L;

			for ( i = 0; i < Win_Before_Mag.size(); ++i ) {
				Win_Before_Mag[i].T=Win_Before_Mag[i].L*Win_Before_Mag[i].rho;
				if ( fabs(Win_Before_Mag[i].X0)<1e-10 ) {
					Win_Before_Mag[i].TR=0;
				}
				else {
					Win_Before_Mag[i].TR=Win_Before_Mag[i].T/Win_Before_Mag[i].X0;
				}
				Win_Before_Mag[i].bt=b(Win_Before_Mag[i].Z)*Win_Before_Mag[i].TR;
			}
			/*end of Set Windows}}}*/

		}
		/*}}}*/

		/*int CalcXS(){{{*/
		int CalcXS()
		{    
			size_t i;
			size_t imax;
			vector<Material> Win_Empty;//tmp use

			Material mixture;

			double offset=p_TCS(2)-p_inter_point_TCS(2);//L=along Z in TCS
			mixture.L=offset;

			if ( E_p<ELECTRON_MASS ) {
				printf("[Warning %s: Line %d] E_p=%g<ELECTRON_MASS=%g\n",__FILE__,__LINE__,E_p,ELECTRON_MASS);
				return 1;
			}

			Q2=4*E_s*E_p*sinsq;
			q2=-Q2;
			btr=AP*(log(Q2/(ELECTRON_MASS*ELECTRON_MASS))-1);//b*t_r

			/*Energy Loss{{{*/
			if ( IsEnergyLoss ) {
				//Fix me: Do I need to count the energy loss due to windows after target?
				//I think so. So I add them.
				imax=Win_Before_Mag.size();
				for ( i=0; i<imax; i++ )
				{
					E_p-=Ion_Loss(E_p,Win_Before_Mag[i]);
					E_p-=Bremss_Loss(E_p,Win_Before_Mag[i].bt);
				}
				imax=Win_After_Mag.size();
				for ( i=0; i<imax; i++ )
				{
					E_p-=Ion_Loss(E_p,Win_After_Mag[i]);
					E_p-=Bremss_Loss(E_p,Win_After_Mag[i].bt);
				}
				if ( E_p<0 ) {
					printf("[Warning %s: Line %d] E_p after Energy Loss=%f<0.\n",__FILE__,__LINE__,E_p);
					return 1;
				}
				Q2=4*E_s*E_p*sinsq;
				q2=-Q2;
			}
			/*}}}*/

			cs_M=sigma_M(E_s,Angle_Deg);

			/*{{{Born Crosss Section Model*/
			//NOTE: I replace the born cross section model with my own XEMC model 
			//   -- Zhihong Ye 01/25/2012
			//The unit of XS is nb/(sr*MeV/c);
			if (IsQFS){
				cs_qfs=sigma_qfs(Target.Z,Target.A,E_s,E_s-E_p,Angle_Deg,Fermi_Moment,NIE,IsQFS_Q2dep);//Q.E. Peak
				cs_del=sigma_del(Target.Z,Target.A,E_s,E_s-E_p,Angle_Deg,Fermi_Moment,DEL_SEP);//Delta
				cs_x=sigma_x(Target.A,E_s,E_s-E_p,Angle_Deg);//DIS
				cs_r1500=sigma_resonance(1500.,Target.A,Target.Z,E_s,E_s-E_p,Angle_Deg,Fermi_Moment);//Resonace 1500 MeV
				cs_r1700=sigma_resonance(1700.,Target.A,Target.Z,E_s,E_s-E_p,Angle_Deg,Fermi_Moment);//Resonace 1700 MeV
				cs_2N=sigma_2N(Target.Z,Target.A,E_s,E_s-E_p,Angle_Deg,Fermi_Moment);//Dip region
				cs_Born=cs_qfs+cs_del+cs_x+cs_r1500+cs_r1700+cs_2N;
			}
			else{ //ZYE
				//Born Cross Section Model from XEM, Last Option : 1->QE+DIS_F1F2IN09, 2->QE Only, 3->DIS_F1F2IN09 only 4->DIS_F1F2IN06_XEM Only
				cs_qe = XEMC_Born(E_s,E_p,Angle_Deg,(int)(Target.A+0.5),(int)(Target.Z+0.5),2); 
				cs_dis = XEMC_Born(E_s,E_p,Angle_Deg,(int)(Target.A+0.5),(int)(Target.Z+0.5),3); 
				cs_Born = cs_qe + cs_dis; 
			}
			cs_Final = cs_Born;
			/*}}}*/

			if ( IsRadCor && cs_Born >= 1e-13 )
			{
				//all nbarn/sr
				cs_q = cs_Final;

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
				if(IsBump){
					tb=Win_i.TR+Target_Up_TR;
					ta=Win_f.TR+Target_Down_TR;
					btb=Win_i.bt+Target_Up_bt;//b*tb
					bta=Win_f.bt+Target_Down_bt;//b*ta
				}
				else{
					tb=Win_i.TR+Target.TR/2;
					ta=Win_f.TR+Target.TR/2;
					btb=Win_i.bt+Target.bt/2;//b*tb
					bta=Win_f.bt+Target.bt/2;//b*ta
				}
				E_s_min=E_p/(1-2*E_p*sinsq_Angle/Target.M);
				E_p_max=E_s/(1+2*E_s*sinsq_Angle/Target.M);
				/*}}}*/
                
				cs_qt=CalcuateQuasiTail(cs_q); //Quasi-Elastic Tail
				tail=CalcuateTail(); //Elastic Tail
				cs_Final=cs_qt+tail;//not minus,we simulate the real situation
				//cerr <<Form("Radiated XS = %e", cs_Final)<<endl;
			}
			return 0;
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

			cerr<<Form("--- Rad-Tail: mp_cor=%e,cs_p=%e,cs_ex=%e,cs_b=%e,cs_b_bt=%e,cs_b_at=%e --- ",mp_cor,cs_p,cs_ex,cs_b,cs_b_bt,cs_b_at)<<endl;

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
				cs_qt_t1=Sim_Integral(1,E_s_min,E_s-R*DeltaE,Num_Of_Steps,&XEMCEvent::sigma_qt_Es_Integrand);
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
				cs_qt_t2=Sim_Integral(1,E_p+DeltaE,(E_p_max),Num_Of_Steps,&XEMCEvent::sigma_qt_Ep_Integrand);
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
		double F_soft() 
		{
			//Multiple-photon correction
			//Phys.Rev.D 12,1884 (A58)
			cerr<<Form(" --- In F_soft: 1=%e,2=%e,3=%e",
					                   w_s/E_s,(btb+btr),w_p/(w_p+E_p))<<endl;
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
			t1=Sim_Integral(0,-1,1,10000,&XEMCEvent::sigma_ex_Integrand);
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
			double lsigma_q = 0.0;

			/*Born Cross Section Model{{{*/
			if(IsQFS){
				lsigma_q=sigma_qfs(Target.Z,Target.A,aEs,aEs-aEp,Angle_Deg,Fermi_Moment,NIE,IsQFS_Q2dep);
				lsigma_q+=sigma_del(Target.Z,Target.A,aEs,aEs-aEp,Angle_Deg,Fermi_Moment,DEL_SEP);//Delta
				lsigma_q+=sigma_x(Target.A,aEs,aEs-aEp,Angle_Deg);//DIS
				lsigma_q+=sigma_resonance(1500.,Target.A,Target.Z,aEs,aEs-aEp,Angle_Deg,Fermi_Moment);//Resonace 1500 MeV
				lsigma_q+=sigma_resonance(1700.,Target.A,Target.Z,aEs,aEs-aEp,Angle_Deg,Fermi_Moment);//Resonace 1700 MeV
				lsigma_q+=sigma_2N(Target.Z,Target.A,aEs,aEs-aEp,Angle_Deg,Fermi_Moment);//Dip region
			}
			else{
				//Born Cross Section Model from XEM, 
				//Last Option : 1->QE+DIS_F1F2IN09, 2->QE Only, 3->DIS_F1F2IN09 only, 4->DIS_F1F2IN06 Only,
				//              5->QE_DIS+F1F2IN09+Coulomb Correction
				lsigma_q = XEMC_Born(aEs,aEp,Angle_Deg,(int)(Target.A+0.5),(int)(Target.Z+0.5),5);
			}
			/*}}}*/
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
		double Sim_Integral(const int& aChoice,const double& aMin,const double& aMax,const int& aNum_Of_Step,double (XEMCEvent::*aFunc)(const double&))
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
							((fabs(EST-SUM)>EPS*AREA||EST==1.0)&&(L>=19))
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
			double lK=0.307075;// cm^2/g for A=1 g/mol
			double lbetasq=1-ELECTRON_MASS*ELECTRON_MASS/(aE0*aE0);
			double lxi=lK/2*aMaterial.Z/aMaterial.A*aMaterial.T/lbetasq;//aT: g/cm^2
			double lhbarwsq=28.816*28.816*aMaterial.rho*aMaterial.Z/aMaterial.A*1e-12;//MeV arho is density of absorber
			double j=0.200;
			double Delta_p=lxi*(log(2*ELECTRON_MASS*lxi/lhbarwsq)+j);
			double lw=4*lxi;
			double result=0;
			if ( aMaterial.Z!=0 && aMaterial.A!=0 && aMaterial.T!=0 && aMaterial.rho!=0 )
				result=gRandom->Landau(Delta_p,lw);
			if ( result>(aE0-ELECTRON_MASS) )
				result=aE0-ELECTRON_MASS;
			if ( result<0 )
				result=0;
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
			if ( abt!=0 )
				result=aE0*pow(gRandom->Rndm()*0.999,1./abt);
			if ( result>(aE0-ELECTRON_MASS) )
				result=aE0-ELECTRON_MASS;
			if ( result<0 )
				result=0;
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

		/*void SetMaterial(Material& aMaterial){{{*/
		void SetMaterial(Material& aMaterial)
		{
			aMaterial.M=aMaterial.A*AMU;//MeV
			if ( aMaterial.L==0 && aMaterial.rho!=0 ) {
				aMaterial.L=aMaterial.T/aMaterial.rho;
			}
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
				//if ( Target.Z==26 ) {//for iron
				//      lpar[0]=1.0; //just satisfy the data
				//      lpar[1]=0.0453828; //just satisfy the data
				//}
				double lCor=lpar[1]*1e6/lQ2+lpar[0];
				result/=lCor;
			}
			return result;
		}
		/*}}}*/

		/*double sigma_del(const double& aA,const double& aE,const double& aomega,const double& aTheta,const double& aPF,const double& aEPSD){{{*/
		double sigma_del(const double& aZ, const double& aA,const double& aE,const double& aomega,const double& aTheta,const double& aPF,const double& aEPSD)
		{
			//Delta
			//aE=MeV,aTheta=deg
			//aPF: fermi momentum
			//double AP0=840; //no idea
			//double AP1=750; //no idea
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
			//  double PFR=60.;
			double GAM2N=20.;
			//  double GAMQFR=40.;
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
			//    double GAMQF=GAMQFR*(aPF/PFR)*(sqrt(QVS)/sqrt(QVSR));
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
		double sigma_resonance(const double& aRes_Energy,const double& aA,const double& aZ,const double& aE,const double& aomega,const double& aTheta,const double& aPF)
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
			//double WR=aE-EPR;
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


		/*double PROD_AND(const double& ax,const double& ay){{{*/
		double PROD_AND(const double& ax,const double& ay)
		{
			double prod;
			prod=TMath::Min(ax,ay);
			//prod=ax+ay-sqrt(ax*ax+ay*ay);
			return prod;
		}
		/*}}}*/


	public:
		//Not safe, but it's ok for physicsists
		/*Member Data from file{{{*/
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


		double cs_M; //Mott Cross Section
		double cs_Final; //Mott Cross Section
		double cs_Born; //Born Cross Section from XEMC
		double cs_qe; //Radiated Cross Section if inlucding rediative effect
		double cs_dis; //Radiated Cross Section if inlucding rediative effect

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

		bool IsDebug; //output to debugfile
		bool IsEnergyLoss; //enable Energy Loss
		double z0;    //target center for generator(cm)
		double T_L;   //target length for generator(cm)
		double T_H;   //target height for generator(cm)
		double Fermi_Moment;  //Fermi_Moment,fermi momentum of target (MeV)
		double NIE;   //NIE,nucleon interaction energy (MeV)
		double DEL_SEP;       //DEL_SEP,delta separation energy (MeV)

		bool IsRadCor; //enable Radiative correction
		bool IsQFS; //enable QFS to calcualte xs instead of XEMC
		bool IsQFS_Q2dep; //enable Q2 dependence when calculating Q.E. Peak
		bool IsBump; 
		double Target_Up_bt;
		double Target_Up_TR;
		double Target_Down_bt;
		double Target_Down_TR;


		/*Only for T ail{{{*/
		int Form_Factor_Id; //which form factor used
		bool  IsInternal_Bre; //enable/disable internal bremsstrahlung
		bool  IsExternal_Bre; //enable/disable external bremsstrahlung
		bool  IsMulti_Photon_Cor; //enable/disable multiple photon correction
		bool  IsAppro; //use sigma_p instead of sigma_ex (approximate exact internal bremsstrahlung)
		double DeltaE; //default=5 MeV only for radiative correction

		int XEMCFLAG;
		double R; //Phys.Rev.D 12 1884 A83
		double E_s_min; //E min from E_p MeV
		double E_p_max; //Ep max from E_s MeV
		double cs_qt;
		double Angle_rec ; //Scattering Angle after reconstruction MeV*MeV
		double Qsq ; //Q2 after reconstruction MeV*MeV
		double Xbj ; //Xbj after reconstruction               
		/*}}}*/

		/*Member Data derived from variables above{{{*/
		//TLorentzVector ( vec, t ) = ( k,E ) in ROOT
		//HCS: Hall Coordinate System
		//TCS: Target Coordinate System
		TLorentzVector s; //s(s3,E_s) 4-momentum position of the incident electron in HCS
		TLorentzVector target_edgepoint_TRCS; //4-momentum position of the edge point in HCS (if theta>0,T_H/2,if <0,-T_H/2)
		TLorentzVector s_TCS; //s_TCS(s3,E_s) 4-momentum position of the incident electron in TCS, and ztg=0
		TLorentzVector p_TCS; //p_TCS(p3,E_p) 4-momentum position of the outgoing electron in TCS, and ztg=0
		TLorentzVector p_inter_point_TCS; //p_inter_point_TCS(p3,E_p) 4-momentum position of the outgoing electron in TCS, and ztg=reactz_TCS
		TLorentzVector p_P; //p_P 4-momentum momentum of the outgoing electron in HCS
		TLorentzVector p_P_TCS; //p_P_TCS 4-momentum momentum of the outgoing electron in TCS

		double Angle; //real scattering angel, not from file,rad
		double Angle_Deg; //real scattering angel, not from file,Deg
		double sinsq_Angle;//sin(Angle/2)^2

		double theta_rad;//scattering angle in XEMC
		double sinsq;//sin(theta/2)^2
		double reactz_TCS;//
		double btr;//b*tr (equivalent radiator)

		double Q2;// MeV*MeV
		double q2;//q2=-Q2 MeV*MeV
		/*}}}*/

		/*Member Data {{{*/
		double beam_x; //cm
		double beam_y; //cm
		double D_x;
		double D_y;
		double HRS_L;
		double E_s;//=incident beam energy (MeV)
		double E_p;//=incident beam energy (MeV)
		double reactz_gen;
		/*void Print(){{{*/
		void Print()
		{
			printf("--------------------\n");
			//refer to http://www.fileformat.info/info/unicode/char/search.htm for \x
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of incident electron [HCS])\n",20,"s(x,y,z,E_s)",10,s(0),10,s(1),10,s(2),10,s(3));
			printf("%*s=(%*.2f,%*.2f,%*.4f,%*.2f) (4-momentum position of target edgepoint [HCS])\n",20,"target_edgepoint_TRCS(x,y,z,E_s)",10,target_edgepoint_TRCS(0),10,target_edgepoint_TRCS(1),10,target_edgepoint_TRCS(2),10,target_edgepoint_TRCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of incident electron [TCS],ztg=0)\n",20,"s_TCS(x,y,z,E_s)",10,s_TCS(0),10,s_TCS(1),10,s_TCS(2),10,s_TCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of outgoing electron [TCS], ztg=0)\n",20,"p_TCS(x,y,z,E_s)",10,p_TCS(0),10,p_TCS(1),10,p_TCS(2),10,p_TCS(3));
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum position of outgoing electron [TCS], ztg=reactz_TCS=%.2f)\n",20,"p_inter_point_TCS(x,y,z,E_s)",10,p_inter_point_TCS(0),10,p_inter_point_TCS(1),10,p_inter_point_TCS(2),10,p_inter_point_TCS(3),reactz_TCS);
			printf("%*s=(%*.2f,%*.2f,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [HCS])\n",20,"p_P(x,y,z,E_s)",10,p_P(0),10,p_P(1),10,p_P(2),10,p_P(3));
			printf("%*s=(%*.2e,%*.2e,%*.2f,%*.2f) (4-momentum momentum of outgoing electron [TCS])\n",20,"p_P_TCS(x,y,z,E_s)",10,p_P_TCS(0),10,p_P_TCS(1),10,p_P_TCS(2),10,p_P_TCS(3));
			printf("%-*s=%*.2f %-*s %-*s\n",15,"E_s",       10,E_s,       8,"MeV",           40,"(Es Incident Energy)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"E_p",       10,E_p,       8,"MeV",           40,"(Ep Scattering Energy)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"theta",     10,theta,     9,"\xc2\xb0",      40,"(\xce\xb8 scattering angle deg)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"theta_rad", 10,theta_rad, 8,"",              40,"(\xce\xb8 scattering angle rad)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"T_theta",   10,T_theta,   9,"\xc2\xb0",      40,"(\xce\xb8t target angle)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"sinsq",     10,sinsq,     8,"",              40,"(sin\xc2\xb2(\xce\xb8/2))");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"M_T",       10,Target.M,  8,"MeV",           40,"(Mass of target)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"beam_x",    10,beam_x,    8,"cm",            40,"(Beam X)");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"beam_y",    10,beam_y,    8,"cm",            40,"(Beam Y)");

			printf("%-*s=%*.2f %-*s %-*s\n",15,"Angle",     10,Angle,     8,"rad",           40,"(real scattering angle[rad])");
			printf("%-*s=%*.2f %-*s %-*s\n",15,"Angle",     10,Angle*RadToDeg(),     9,"\xc2\xb0",           40,"(real scattering angle[deg])");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"Q2",        10,Q2,        9,"MeV\xc2\xb2",   40,"(Q\xc2\xb2)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"q2",        10,q2,        9,"MeV\xc2\xb2",   40,"(-Q\xc2\xb2)");
			printf("%-*s=%*.2e %-*s %-*s\n",15,"btr",       10,btr,       8,"rad_len",       40,"(b*tr equivalent radiator unit in rad_len)");
			PrintMaterial(Target);
			PrintMaterial(Win_i);
			PrintMaterial(Win_f);
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
