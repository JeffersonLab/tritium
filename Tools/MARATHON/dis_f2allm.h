//      allm97, NMC published measured points Q2>0.75 GeV2
//      for values Q<1 use data of E665!
//      parameterization of F2 , according to
//      H.Abramowicz and A.Levy, hep-ph/9712415
//
//      3*10-6 < x  < 0.85, W2>3GeV2
//      0.   < Q2 < 5000 GeV2, dof=0.97
//
static const double ubarn_to_nbarn= 1e6; //ubarn-->nbarn

void f2allm_f2d(double xb, double Q2, double* F2p, double* F2d){/*{{{*/

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
    double W2=Q2*(1./xb -1.)+XMP2;
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
        Z=1.-xb;
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

    double f2allm_f2p=Q2/(Q2+M02)*(F2P+F2R);

    //double f2nf2p = 0.976+xb*(-1.34+xb*(1.319+xb*(-2.133+xb*1.533)));
    //With Q2-dependence correction
    double f2nf2p=(0.979+xb*(-0.286+xb*(-1.554+xb*(3.704+xb*(-2.443))))
            +(0.00297-0.0626*xb+0.0968*xb*xb)*log(Q2))*2.-1.;

    double f2allm_f2d=f2allm_f2p*0.5*(f2nf2p+1.);

    F2p[0] = f2allm_f2p;
    F2d[0] = f2allm_f2d;
}
/*}}}*/


double sigma_f2allm(double E0, double Ep, double Theta, int A, int Z){/*{{{*/

    double SinSQ = pow( sin(Theta * 0.5), 2);
    double CosSQ = pow( sin(Theta * 0.5), 2);
    double Q2 = 4.0 * E0 * Ep * SinSQ;
    double nu = E0 - Ep;
    double xb = Q2 / (2.0 * 0.938272 * nu);

    double F2p = 0.0, F2d = 0.0;
    f2allm_f2d(xb, Q2,  &F2p, &F2d);

    double sigma_p = 4.0 * pow(1/137.036, 2) * Ep*Ep / Q2/Q2
    	* (1./nu * CosSQ + 1/0.938272/xb * SinSQ)
    	* F2p;
    double sigma_d = 4.0 * pow(1/137.036, 2) * Ep*Ep / Q2/Q2
        * (1./nu * CosSQ + 1/0.938272/xb * SinSQ)
    	* F2d;

    if(A==1 && Z==1)
        return sigma_p*ubarn_to_nbarn;
    else if(A==2 && Z==1)
        return sigma_d*ubarn_to_nbarn;
    else
        return (Z*sigma_p+(A-Z)*(sigma_d-sigma_p))*ubarn_to_nbarn;
}
/*}}}*/


double f2allm_sigma_p(double E0, double Ep, double Theta){/*f2allm_sigma_p{{{*/

    double SinSQ = pow( sin(Theta * 0.5), 2);
    double CosSQ = pow( sin(Theta * 0.5), 2);
    double Q2 = 4.0 * E0 * Ep * SinSQ;
    double nu = E0 - Ep;
    double xb = Q2 / (2.0 * 0.938272 * nu);

    double F2p = 0.0, F2d = 0.0;
    f2allm_f2d(xb, Q2,  &F2p, &F2d);

    //double sigma_p = 4.0 * pow(3.1415926, 2) * (1/137.036) / Q2 / (1.-xb)
    //    * (Q2+4.0*pow(0.938272,2)*xb*xb) / Q2
    //    * F2p;
    double sigma_p = 4.0 * pow(1/137.036, 2) * Ep*Ep / Q2/Q2
    	* (1./nu * CosSQ + 1/0.938272/xb * SinSQ)
    	* F2p;

    //printf("** xb=%4.3f, Q2=%5.2f, F2p = %5.3f, F2n = %5.3f\n", xb, Q2, F2p, (2.*F2d-F2p));

    return sigma_p*ubarn_to_nbarn;//ubarn/GeV into nb/GeV
}
/*}}}*/


double f2allm_sigma_d(double E0, double Ep, double Theta){/*f2allm_sigma_d{{{*/

    double SinSQ = pow( sin(Theta* 0.5), 2);
    double CosSQ = pow( sin(Theta* 0.5), 2);
    double Q2 = 4.0 * E0 * Ep * SinSQ;
    double nu = E0 - Ep;
    double xb = Q2 / (2.0 * 0.938272 * nu);

    double F2p = 0.0, F2d = 0.0;
    f2allm_f2d(xb, Q2,  &F2p, &F2d);

    //double sigma_d = 4.0 * pow(3.1415926, 2) * (1/137.036) / Q2 / (1.-xb)
        //* (Q2+4.0*pow(0.938272,2)*xb*xb) / Q2
        //* F2d;
    double sigma_d = 4.0 * pow(1/137.036, 2) * Ep*Ep / Q2/Q2
        * (1./nu * CosSQ + 1/0.938272/xb * SinSQ)
        * F2d;//in ubarn/GeV

    return sigma_d*ubarn_to_nbarn;//ubarn/GeV into nb/GeV
}
/*}}}*/


double f2allm_sigma_n(double E0, double Ep, double Theta){/*f2allm_sigma_n{{{*/

    double SinSQ = pow( sin(Theta * 0.5), 2);
    double CosSQ = pow( sin(Theta * 0.5), 2);
    double Q2 = 4.0 * E0 * Ep * SinSQ;
    double nu = E0 - Ep;
    double xb = Q2 / (2.0 * 0.938272 * nu);

    double F2p = 0.0, F2d = 0.0;
    f2allm_f2d(xb, Q2,  &F2p, &F2d);

    //double sigma_n = 4.0 * pow(3.1415926, 2) * (1/137.036) / Q2 / (1.-xb)
    //            * (Q2+4.0*pow(0.938272,2)*xb*xb) / Q2
    //			* (F2d-F2p);
    double sigma_n = 4.0 * pow(1/137.036, 2) * Ep*Ep / Q2/Q2
        * (1./nu * CosSQ + 1/0.938272/xb * SinSQ)
        * (2.*F2d-F2p);
    
    return sigma_n*ubarn_to_nbarn;//ubarn/GeV into nb/GeV
}/*}}}*/
