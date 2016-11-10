#ifndef DIS_H
#include "christy_bosted_inelastic.h"
static const double GeV2_to_nbarn= 0.3894*1e6; //GeV^-2-->nbarn
        extern "C"/*{{{*/
            {
                /*  ISET             PDF*/
                //
                //  100      | CJ12_min central PDF
                //  101-138  | CJ12_min error PDF sets
                //           |
                //  200      | CJ12_mid central PDF
                //  201-238  | CJ12_mid error PDF sets
                //           |
                //  300      | CJ12_max central PDF
                //  301-338  | CJ12_max error PDF sets
                void setcj_(int*);

                // Flavor index: Iptn = -5:5 for bb, cb, sb, db, ub, g, u, d, s, c, b,
                void cjpdf_(const int*,const double*, const double*);

                //CJpdf_all(x,Q,bb,cb,sb,db,ub,g,u,d,s,c,b) 
                void cjpdf_all_(const double*, const double*, 
                        const double*, const double*, const double*, const double*, const double*, 
                        const double*, const double*, const double*, const double*, const double*, const double*); 
            
                void f1f2in09_(const double*,const double*, const double*, const double*, double*, double*, double*);
            }
            /*}}}*/

class DIS{
    public:
        DIS(){
          fInitCJ = false;//default is not to initialize CJ tables
        }
        virtual ~DIS(){
            //	delete fRan;
        }
    public:
        void Print(){/*{{{*/
            printf("***************************************************************\n");
            printf("**  Run SetCJ(N) before using CJ, N=100->min, 200->mid, 300-max\n");
            printf("**  Change kinematic setting by running SetKin(E0, Ep, Theta*Rad)\n");
            printf("** \n");
            printf("** Three Models to get F2: F1F2IN09, F2ALLM, and CJ\n");
            printf("**  All get F2p and F2n, and F1F2IN09 can also get F2d\n");
            printf("**  Call:  DIS->GetF2(Model, Nucleon),like Proton, Neutron, Deutron\n");
            printf("** \n");
            printf("** Four Models to get DIS XS: PB, F1F2IN09, F2ALLM, and CJ\n");
            printf("**  Call:  DIS->Sigma(Model, Nucleon),like Proton, Neutron, Deutron\n");
            printf("***************************************************************\n");
        }/*}}}*/

        void SetCJ(){/*{{{*/
            /*  ISET             PDF*/
            //
            //  100      | CJ12_min central PDF
            //  101-138  | CJ12_min error PDF sets
            //           |
            //  200      | CJ12_mid central PDF
            //  201-238  | CJ12_mid error PDF sets
            //           |
            //  300      | CJ12_max central PDF
            //  301-338  | CJ12_max error PDF sets
            //int iset = 200;
            //setcj_(&iset);
            int iset =100;
            setcj_(&iset);
            fInitCJ = true;
        }/*}}}*/

        void SetCJ(int iset = 100){/*{{{*/
            /*  ISET             PDF*/
            //
            //  100      | CJ12_min central PDF
            //  101-138  | CJ12_min error PDF sets
            //           |
            //  200      | CJ12_mid central PDF
            //  201-238  | CJ12_mid error PDF sets
            //           |
            //  300      | CJ12_max central PDF
            //  301-338  | CJ12_max error PDF sets
            //int iset = 200;
            //setcj_(&iset);

            setcj_(&iset);
            fInitCJ = true;
        }/*}}}*/

        double GetRnp(const double kXbj, const double kQ2){/*{{{*/
            //arXiv:0805.3116v2, J. Arrington et al
            const double p1 =  0.816;
            const double p2 = -0.661;
            const double p3 =  0.184;
            const double p4 =  5.509;
            const double p5 = -0.034;
            const double p6 =  8.714;
            const double p7 = -0.072;
            const double p8 =  0.450;

            double xi = 2.0*kXbj/(1.0+sqrt(1.0+4.0*pow(0.938272,2)*kXbj*kXbj/kQ2));

            double R = (p1+p2*xi) + p3*exp(-p4*xi) + p5*exp(-p6*(1-xi)) + p7*pow(std::max(0.0, xi-p8), 2);

            return R;
        }/*}}}*/

        void SetModel(const string& kModel){/*{{{*/
            fModel = kModel;
            if(fModel.find("CJ",0)!=string::npos){
                printf("^^^ CJ is used, need to run SetCJ(N).\n");
                printf("^^^ N=100(Min),200(Mid) and 300(Max). Default is 100.\n");
                if(!fInitCJ)
                    SetCJ(100);
            }
        }/*}}}*/

        void SetF1F2Kin(double kXbj, double kQ2){/*{{{*/
            //F1 or F2 only needs xb and Q2;
            fQ2 = kQ2;
            fXbj = kXbj;
            fW2 = pow(0.938272, 2) + fQ2/fXbj - fQ2;
        }/*}}}*/

        void SetKin(double kE0, double kEp, double kTheta){/*{{{*/
            fE0 = kE0;
            fEp = kEp;
            fTheta = kTheta;
            fSinSQ = pow( sin(fTheta * 0.5), 2);
            fCosSQ = pow( cos(fTheta * 0.5), 2);
            fQ2 = 4.0 * fE0 * fEp * fSinSQ;
            fNu = fE0 - fEp;
            fXbj = fQ2 / (2.0 * 0.938272 * fNu);
            fW2 = pow(0.938272, 2) + 2.0*0.938272*fNu - fQ2;
        }/*}}}*/

        void SetKin(string kModel, double kE0, double kEp, double kTheta){/*{{{*/
            fE0 = kE0;
            fEp = kEp;
            fTheta = kTheta;
            fModel = kModel;
            fSinSQ = pow( sin(fTheta * 0.5), 2);
            fCosSQ = pow( cos(fTheta * 0.5), 2);
            fQ2 = 4.0 * fE0 * fEp * fSinSQ;
            fNu = fE0 - fEp;
            fXbj = fQ2 / (2.0 * 0.938272 * fNu);
            fW2 = pow(0.938272, 2) + 2.0*0.938272*fNu - fQ2;

            if(fModel.find("CJ",0)!=string::npos){
                printf("^^^ CJ is used, need to run SetCJ(N).\n");
                printf("^^^ N=100(Min),200(Mid) and 300(Max). Default is 100.\n");
                if(!fInitCJ)
                    SetCJ(100);
            }

        }/*}}}*/

        double GetF2(string kNucleon){/*{{{*/
           //Note: return per-nucleus for A>1
            double kF1 = 0.0;
            double kF2 = 0.0;

            if(fModel.find("F1F2IN09",0)!=std::string::npos){
                if(kNucleon.find("Proton",0)!=std::string::npos)
                    Bosted_F1F209(1, 1, fQ2, fW2, &kF1, &kF2);
                else if(kNucleon.find("Neutron",0)!=std::string::npos)
                    Bosted_F1F209(1, 0, fQ2, fW2, &kF1, &kF2);
                else if(kNucleon.find("Deutron",0)!=std::string::npos)
                    Bosted_F1F209(2, 1, fQ2, fW2, &kF1, &kF2);
                else
                    cout<<"*** ERROR, In GetF2(), I don't know the nucleus"<<endl;
            }
            if(fModel.find("F2ALLM",0)!=std::string::npos){
                double kF2p = GetF2_F2ALLM(fXbj, fQ2, "Proton");
                double kF2n = GetF2_F2ALLM(fXbj, fQ2, "Neutron");
                if(kNucleon.find("Proton",0)!=std::string::npos)
                    kF2 = kF2p;
                else if(kNucleon.find("Neutron",0)!=std::string::npos)
                    kF2 = kF2n;
                else if(kNucleon.find("Deutron",0)!=std::string::npos)
                    kF2 = kF2p+kF2n;
                else
                    cout<<"*** ERROR, In GetF2(), I don't know the nucleus"<<endl;
            }
            if(fModel.find("CJ",0)!=std::string::npos){
                double kF2p = GetF2_CJ(fXbj, fQ2, "Proton");
                double kF2n = GetF2_CJ(fXbj, fQ2, "Neutron");
                if(kNucleon.find("Proton",0)!=std::string::npos)
                    kF2 = kF2p;
                else if(kNucleon.find("Neutron",0)!=std::string::npos)
                    kF2 = kF2n;
                else if(kNucleon.find("Deutron",0)!=std::string::npos)
                    kF2 = kF2p+kF2n;
                else
                    cout<<"*** ERROR, In GetF2(), I don't know the nucleus"<<endl;
            }
            return kF2;
        }/*}}}*/
        
        double GetF2(string kModel, string kNucleon){/*{{{*/
           //Note: return per-nucleus for A>1
            double kF1 = 0.0;
            double kF2 = 0.0;

            if(kModel.find("F1F2IN09",0)!=std::string::npos){
                if(kNucleon.find("Proton",0)!=std::string::npos)
                    Bosted_F1F209(1, 1, fQ2, fW2, &kF1, &kF2);
                else if(kNucleon.find("Neutron",0)!=std::string::npos){
                    double kF1p = 0., kF2p = 0.;
                    Bosted_F1F209(1, 1, fQ2, fW2, &kF1p, &kF2p);
                    double f2nf2p = GetRnp(fXbj,fQ2);//use John's F2n/F2p ratio
                    kF2 = f2nf2p * kF2p; 
                    
                    //Bosted_F1F209(1, 0, fQ2, fW2, &kF1, &kF2);
                }
                else if(kNucleon.find("Deutron",0)!=std::string::npos)
                    Bosted_F1F209(2, 1, fQ2, fW2, &kF1, &kF2);
                else
                    cout<<"*** ERROR, In GetF2(), I don't know the nucleus"<<endl;
            }
            if(kModel.find("F2ALLM",0)!=std::string::npos){
                double kF2p = GetF2_F2ALLM(fXbj, fQ2, "Proton");
                double kF2n = GetF2_F2ALLM(fXbj, fQ2, "Neutron");
                if(kNucleon.find("Proton",0)!=std::string::npos)
                    kF2 = kF2p;
                else if(kNucleon.find("Neutron",0)!=std::string::npos)
                    kF2 = kF2n;
                else if(kNucleon.find("Deutron",0)!=std::string::npos)
                    kF2 = kF2p+kF2n;
                else
                    cout<<"*** ERROR, In GetF2(), I don't know the nucleus"<<endl;
            }
            if(kModel.find("CJ",0)!=std::string::npos){
                double kF2p = GetF2_CJ(fXbj, fQ2, "Proton");
                double kF2n = GetF2_CJ(fXbj, fQ2, "Neutron");
                if(kNucleon.find("Proton",0)!=std::string::npos)
                    kF2 = kF2p;
                else if(kNucleon.find("Neutron",0)!=std::string::npos)
                    kF2 = kF2n;
                else if(kNucleon.find("Deutron",0)!=std::string::npos)
                    kF2 = kF2p+kF2n;
                else
                    cout<<"*** ERROR, In GetF2(), I don't know the nucleus"<<endl;
            }
        return kF2;
        }/*}}}*/

        double GetF1(string kNucleon){/*{{{*/
           //Note: return per-nucleus for A>1
            double kF1 = 0.0;
            double kF2 = 0.0;

            if(fModel.find("F1F2IN09",0)!=std::string::npos){
                if(kNucleon.find("Proton",0)!=std::string::npos)
                    Bosted_F1F209(1, 1, fQ2, fW2, &kF1, &kF2);
                else if(kNucleon.find("Neutron",0)!=std::string::npos){
                    //Bosted_F1F209(1, 0, fQ2, fW2, &kF1, &kF2);
                    double kF1p = 0., kF2p = 0.;
                    Bosted_F1F209(1, 1, fQ2, fW2, &kF1p, &kF2p);
                    double f2nf2p = GetRnp(fXbj,fQ2);//use John's F2n/F2p ratio
                    kF2 = f2nf2p * kF2p; 
                }
                 else if(kNucleon.find("Deutron",0)!=std::string::npos)
                    Bosted_F1F209(2, 1, fQ2, fW2, &kF1, &kF2);
                else
                    cout<<"*** ERROR, In GetF1(), I don't know the nucleus"<<endl;
            }
            else{
                kF2 = GetF2(kNucleon.c_str());
                kF1 = kF2/2.0/fXbj;
            }

            return kF1;
        }/*}}}*/
        
        double GetF1(string kModel, string kNucleon){/*{{{*/
           //Note: return per-nucleus for A>1
            double kF1 = 0.0;
            double kF2 = 0.0;

            if(kModel.find("F1F2IN09",0)!=std::string::npos){
                if(kNucleon.find("Proton",0)!=std::string::npos)
                    Bosted_F1F209(1, 1, fQ2, fW2, &kF1, &kF2);
                else if(kNucleon.find("Neutron",0)!=std::string::npos)
                    Bosted_F1F209(1, 0, fQ2, fW2, &kF1, &kF2);
                else if(kNucleon.find("Deutron",0)!=std::string::npos)
                    Bosted_F1F209(2, 1, fQ2, fW2, &kF1, &kF2);
                else
                    cout<<"*** ERROR, In GetF1(), I don't know the nucleus"<<endl;
            }
            else{
                kF2 = GetF2(kModel.c_str(), kNucleon.c_str());
                kF1 = kF2/2.0/fXbj;
            }

            return kF1;
        }/*}}}*/

        double Sigma(string kNucleon){/*GetDIS_Sigma{{{*/
           //Note: return nb/GeV/sr
            double xs = 0.0;
            if(fModel.find("PB",0)!=std::string::npos){
                if(kNucleon.find("Proton",0)!=std::string::npos)
                    xs = sigma_p(fE0, fTheta, fEp);//Peter Bosted's DIS model
                else if(kNucleon.find("Neutron",0)!=std::string::npos)
                    xs = sigma_n(fE0, fTheta, fEp);//Peter Bosted's DIS model
                else if(kNucleon.find("Deutron",0)!=std::string::npos)
                    xs = sigma_d(fE0, fTheta, fEp);//Peter Bosted's DIS model
                else
                    cout<<"*** ERROR, In DIS->Sigma(), I don't know the nucleus"<<endl;
             }
            else{
                double F2 = GetF2(kNucleon.c_str());
                xs = 4.0 * pow(1/137.036, 2) * fEp*fEp / fQ2/fQ2 
                    * (1./fNu * fCosSQ + 1/0.938272/fXbj * fSinSQ) * F2*GeV2_to_nbarn;
            }
            return xs;
        }
        /*}}}*/
        
        double Sigma(string kModel, string kNucleon){/*GetDIS_Sigma{{{*/
           //Note: return nb/GeV/sr
            double xs = 0.0;
            if(kModel.find("PB",0)!=std::string::npos){
                if(kNucleon.find("Proton",0)!=std::string::npos)
                    xs = sigma_p(fE0, fTheta, fEp);//Peter Bosted's DIS model
                else if(kNucleon.find("Neutron",0)!=std::string::npos)
                    xs = sigma_n(fE0, fTheta, fEp);//Peter Bosted's DIS model
                else if(kNucleon.find("Deutron",0)!=std::string::npos)
                    xs = sigma_d(fE0, fTheta, fEp);//Peter Bosted's DIS model
                else
                    cout<<"*** ERROR, In DIS->Sigma(), I don't know the nucleus"<<endl;
             }
            else{
                double F2 = GetF2(kModel.c_str(), kNucleon.c_str());
                xs = 4.0 * pow(1/137.036, 2) * fEp*fEp / fQ2/fQ2 
                    * (1./fNu * fCosSQ + 1/0.938272/fXbj * fSinSQ) * F2*GeV2_to_nbarn;
            }
            return xs;
        }
        /*}}}*/

        double Sigma( int kA, int kZ){/*{{{*/
           //Note: return nb/GeV/sr
            double xs=0;
            if(kA==3&&kZ==1)//Tritum
                xs = Sigma("Deutron")+Sigma("Neutron");
            else if(kA==3&&kZ==2)//Helium3
                xs = Sigma("Deutron")+Sigma("Proton");
            else
                xs = (kZ*Sigma("Proton")+(kA-kZ)*Sigma("Neutron"));
            return xs;
        }
        /*}}}*/
        
        double Sigma(string kModel, int kA, int kZ){/*{{{*/
           //Note: return nb/GeV/sr
            double xs=0;
            if(kA==3&&kZ==1)//Tritum
                xs = Sigma(kModel.c_str(), "Deutron")+Sigma(kModel.c_str(), "Neutron");
            else if(kA==3&&kZ==2)//Helium3
                xs = Sigma(kModel.c_str(), "Deutron")+Sigma(kModel.c_str(), "Proton");
            else
                xs = (kZ*Sigma(kModel.c_str(), "Proton")+(kA-kZ)*Sigma(kModel.c_str(), "Neutron"));
            return xs;
        }
        /*}}}*/

    private:
        /*CTEQ-JLab{{{*/

        double GetF2_CJ(const double x, const double Q2, const string nucleon){

            //define PDFs:
            double u=0., d=0., s=0., c=0., b=0.,g=0.;
            double ubar=0., dbar=0., sbar=0., cbar=0., bbar=0.;
            double Q = sqrt(Q2);
            cjpdf_all_(&x, &Q, &bbar, &cbar, &sbar, &dbar, &ubar, &g, &u, &d, &s, &c, &b);

            double F2=0;
            if(nucleon.find("Proton",0)!=string::npos)
                F2 = x*( pow(2./3.,2)*(u+ubar) + pow(-1./3.,2)*(d+dbar) +pow(-1./3.,2)*(s+sbar) );
            else if(nucleon.find("Neutron",0)!=string::npos)
                F2 = x*( pow(-1./3.,2)*(u+ubar) + pow(2./3.,2)*(d+dbar) +pow(-1./3.,2)*(s+sbar) );
            else
                cout<<"***ERROR**, No Proton or Neutron in GetF2_CJ!!!"<<endl;
            return F2;
        }
        /*}}}*/

        /*Bosted's F1F209{{{*/
        void Bosted_F1F209(const int kA, const int kZ, const double kQ2, const double kW2, double* kF1, double* kF2)
        {
            /*
            //--------------------------------------------------------------------
            // This is the updated version of F1F2 subroutine. For more info:
            //    http://arxiv.org/abs/1203.2262
            // Original Code locates at:
            //    https:/userweb.jlab.org/~bosted/fits.html
            //
            // Directly call the FORTRAN subroutine instead of converting into C++
            //
            //    --- Zhihong Ye, 12/03/2012
            //--------------------------------------------------------------------
            */

            //=======================================
            // Updated F1F2 in Fortran Subroutines
            //=======================================
            double f1 = -1.1, f2 = -1.2, rc = -1.3;
            double dZ = (double) (kZ); double dA = (double) (kA);
            f1f2in09_(&dZ, &dA, &kQ2, &kW2, &f1, &f2, &rc);
            if(kA==2 && kZ==1 && kW2<-1e-33){
                // DEBUG -- For D2, Can not get F1 and F2 while Wsq < 0, Z. Ye, 12/05/2012
                f1 = 0.0; f2 = 0.0; rc = 0.0;
            }
            kF1[0] = f1;
            kF2[0] = f2;
            return;
        }
        /*}}}*/

        /*F2ALLM97{{{*/
        double GetF2_F2ALLM(const double Xbj, const double Q2, const string Nucleon){
            //      allm97, NMC published measured points fQ2>0.75 GeV2
            //      for values Q<1 use data of E665!
            //      parameterization of F2 , according to
            //      H.Abramowicz and A.Levy, hep-ph/9712415
            //
            //      3*10-6 < x  < 0.85, W2>3GeV2
            //      0.   < Q2 < 5000 GeV2, dof=0.97
            //
            /*Define Constants{{{*/	
            const double S11 = 0.28067;
            const double S12 = 0.22291;
            const double S13 = 2.1979;
            const double A11 = -0.0808;
            const double A12 = -0.44812;
            const double A13 = 1.1709;
            const double B11 = pow(0.60243, 2);
            const double B12 = pow(1.3754, 2);
            const double B13 = 1.8439;
            const double S21 = 0.80107;
            const double S22 = 0.97307;
            const double S23 = 3.4942;
            const double A21 = 0.58400;
            const double A22 = 0.37888;
            const double A23 = 2.6063;
            const double B21 = pow(0.10711, 2);
            const double B22 = pow(1.9386, 2);
            const double B23 = 0.49338;

            const double M02=0.31985;
            const double M12=49.457;
            const double M22   =   0.15052;
            const double LAM2=0.065270;
            const double Q02=0.46017+LAM2;        
            //const double ALFA=112.2;
            const double XMP2=0.8802;//proton mass square
            /*}}}*/
            double W2=Q2*(1./Xbj -1.)+XMP2;
            //double W=sqrt(W2);

            double S, Z, XP, AP, BP, SP, F2P, XR, AR, BR, SR, F2R; 
            if(Q2<=1e-6){//very very low Q2, -->0
                S=0.;
                Z=1.;
                //
                //  POMERON
                //
                XP=1./(1.+(W2-XMP2)/(Q2+M12));
                AP=A11;
                SP=S11;
                F2P=SP*pow(XP, AP);
                //                                              
                //  REGGEON
                //
                XR=1./(1.+(W2-XMP2)/(Q2+M22));
                AR=A21;
                BR=B21;
                SR=S21;
                F2R=SR*pow(XR, AR);
            }
            else{
                S=log(log((Q2+Q02)/LAM2)/log(Q02/LAM2));
                Z=1.-Xbj;
                //
                //  POMERON
                //
                XP=1./(1.+(W2-XMP2)/(Q2+M12));
                AP=A11+(A11-A12)*(1./(1.+pow(S, A13))-1.);
                BP=B11+B12*pow(S, B13);
                SP=S11+(S11-S12)*(1./(1.+pow(S, S13))-1.);
                F2P=SP*pow(XP, AP)*pow(Z, BP);
                //
                //  REGGEON
                //
                XR=1./(1.+(W2-XMP2)/(Q2+M22));
                AR=A21+A22*pow(S, A23);
                BR=B21+B22*pow(S, B23);
                SR=S21+S22*pow(S, S23);
                F2R=SR*pow(XR, AR)*pow(Z, BR);
            }

            double f2allm = 0.0;

            if(Nucleon.find("Proton",0)!=string::npos)
                f2allm=Q2/(Q2+M02)*(F2P+F2R);
            if(Nucleon.find("Neutron",0)!=string::npos){
                //double f2nf2p = 0.976+Xbj*(-1.34+Xbj*(1.319+Xbj*(-2.133+Xbj*1.533)));
                //With Q2-dependence correction
                //double f2nf2p=(0.979+Xbj*(-0.286+Xbj*(-1.554+Xbj*(3.704+Xbj*(-2.443))))
                //       +(0.00297-0.0626*Xbj+0.0968*Xbj*Xbj)*log(Q2))*2.-1.;
                double f2nf2p = GetRnp(Xbj, Q2);//use John's F2n/F2p ratio
                f2allm=Q2/(Q2+M02)*(F2P+F2R)*f2nf2p;
            }
            return f2allm;
        }
        /*}}}*/

    private:
        string fModel;

        double fE0;//GeV
        double fEp;//GeV
        double fTheta;//Rad
        double fSinSQ;
        double fCosSQ;
        double fQ2;
        double fNu;
        double fXbj;
        double fW2;
        bool fInitCJ;
};
#endif
