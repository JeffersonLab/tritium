	
/////////////////////////////////////////
/// BCM Constant
//////////////////////////////////////////

/*BCM{{{*/
// Note: Only Left U1,D1 are calibrated, but D1 had problem during run# 3835~3888.
//       Right BCM looks like having problems.
const Double_t u_Coeff[2][3]=
  {
    //[2]
    //0=Left 1=Right
    //[3]
    //0=1, 1=3, 2=10 u1,u3,u10
//Alpha: U1    U3    U10         //Updated in 07/20/2012
    {2070.0,    0.0, 0.0},            //u Coeff to calculate Charge for left arm
    {1035.0, 3178.0, 0.0},            //u Coeff to calculate Charge for right arm
  };
//                 Offset        D1     D3    D10        //Updated in 07/20/2012
const Double_t u_Offset[2][3]={{400.0,   0.0, 0.0},    //the offsets of the u1,u3,u10 for left arm
            			       {200.0, 174.0, 0.0},};  //the offsets of the u1,u3,u10 for right arm
const Double_t d_Coeff[2][3]=
  {
    //[2]
    //0=Left 1=Right
    //[3]
    //0=1, 1=3, 2=10 u1,u3,u10
   
//Alpha: D1    D3    D10        //Updated in 07/20/2012
    {1263.0,    0.0, 0.0},      //d Coeff to calculate Charge for left arm
    {1263.0, 3902.0, 0.0},      //d Coeff to calculate Charge for right arm
  };
//                 Offset        D1    D3   D10       //Updated in 07/20/2012
const Double_t d_Offset[2][3]={{20.0,  0.0, 0.0},     //the offsets of the d1,d3,d10 for left arm
			                   {37.0, 22.0, 0.0},};   //the offsets of the d1,d3,d10 for right arm
