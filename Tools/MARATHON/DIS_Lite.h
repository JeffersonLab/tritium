#ifndef DIS_H
static const double GeV2_to_nbarn= 0.3894*1e6; //GeV^-2-->nbarn

class DIS{
    public:
        DIS(){
        
        }
        virtual ~DIS(){
            //	delete fRan;
        }
    public:
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

        double GetF2(string kNucleon){/*{{{*/
            //Note: return per-nucleus for A>1
            double kF2 = 0.0;

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

            return kF2;
        }/*}}}*/

        double GetF1(string kNucleon){/*{{{*/
           //Note: return per-nucleus for A>1
            double kF1 = 0.0;
            double kF2 = 0.0;

            kF2 = GetF2(kNucleon.c_str());
            kF1 = kF2/2.0/fXbj;

            return kF1;
        }/*}}}*/
        
        double Sigma(string kNucleon){/*GetDIS_Sigma{{{*/
            //Note: return nb/GeV/sr
            double xs = 0.0;
            double F2 = GetF2(kNucleon.c_str());
            xs = 4.0 * pow(1/137.036, 2) * fEp*fEp / fQ2/fQ2 
                * (1./fNu * fCosSQ + 1/0.938272/fXbj * fSinSQ) * F2*GeV2_to_nbarn;

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
       
        double Xbj(){ return fXbj; }
        double Q2(){ return fQ2; }
        double W2(){ return fW2; }
        double Nu(){ return fNu; }

    private:
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
        double fE0;//GeV
        double fEp;//GeV
        double fTheta;//Rad
        double fSinSQ;
        double fCosSQ;
        double fQ2;
        double fNu;
        double fXbj;
        double fW2;
};
#endif
