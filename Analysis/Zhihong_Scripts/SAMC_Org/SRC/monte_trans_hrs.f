c
c Spectrometer transfer functions based on comm_8_dir.dat for electron and
c hcomm8_dir.dat for hadron. SNAKE has been modified to give the "correct" 
c transverse properties
c                                 4/24/02 -jjl
c
c
c typical call: answer = function(x,5)
c INPUTS: x = 5 or more element array 
c              x(1)=x0  (meters)
c              x(2)=theta0 (really tan(theta0))
c              x(3)=y0   (meters)
c              x(4)=phi0 (really tan(phi0))
c              x(5)=delta (fractional value NOT percent)
c         M=5
c
c OUTPUT: units are the same as inputs
c 
c NOMENCLATURE: function name = prefix + _e_ or _h_ +suffix
c           prefixes:     x means xfinal
c                         t means thetafinal
c                         y means yfinal
c                         p means phifinal
c                         l means pathlength difference from central trajectory
c     
c           suffixes:     fp means target to focus
c                         q1ex means target to Q1 exit
c                         dvac means target to Vacuum box at dipole entrance
c                         dent means target to Dipole entrance
c                         dext means target to dipole exit
c                         q3en means target to Q3 entrance
c                         q3ex means target to Q3 exit
c
c          _e_ is for electron arm
c          _h_ is for hadron arm
c
      function x_e_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/  0.2291023E-02/
      data xmin/
     1 -0.49999E-01,-0.59933E-01,-0.49999E-01,-0.39798E-01,-0.49962E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49871E-01, 0.59944E-01, 0.49999E-01, 0.39798E-01, 0.49911E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23078892E-02, 0.12167830E+00,-0.35640283E-06, 0.22157644E-01,
     +  0.31249220E-02, 0.12923281E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x51 = x5
c
c                  function
c
      x_e_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)*x11            *x51
c
      return
      end
      function t_e_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/ -0.2017684E-03/
      data xmin/
     1 -0.49999E-01,-0.59933E-01,-0.49999E-01,-0.39798E-01,-0.49962E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49871E-01, 0.59944E-01, 0.49999E-01, 0.39798E-01, 0.49911E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23196600E-03,-0.26751672E-02, 0.32225914E-05,-0.25083920E-01,
     + -0.72884180E-06, 0.42211463E-05, 0.28944534E-02, 0.62753790E-06,
     + -0.26580210E-05,-0.19738180E-03,-0.23147034E-03,-0.34753754E-03,
     +  0.11475461E-02,-0.10901972E-03,-0.23569600E-03,-0.17024830E-03,
     + -0.93260860E-04,-0.11148390E-02,-0.40971572E-03,-0.12338403E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_e_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      t_e_q1ex    =t_e_q1ex    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)*x11        *x42    
     8  +coeff( 17)*x12*x21            
      t_e_q1ex    =t_e_q1ex    
     9  +coeff( 18)    *x23*x31*x41    
     1  +coeff( 19)    *x21*x33*x41    
     2  +coeff( 20)    *x23    *x42    
c
      return
      end
      function y_e_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49999E-01,-0.59933E-01,-0.49999E-01,-0.39798E-01,-0.49962E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49871E-01, 0.59944E-01, 0.49999E-01, 0.39798E-01, 0.49911E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.81841350E-01, 0.17506140E+00,-0.16923380E-02,-0.26534970E-02,
     + -0.36940103E-03,-0.44431910E-03, 0.38019860E-03,-0.42201260E-03,
     + -0.32199462E-03, 0.91815942E-04, 0.31874010E-03, 0.12476630E-03,
     + -0.95577874E-04,-0.93289323E-04, 0.10052420E-03,-0.28514240E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_e_q1ex    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31    *x51
     4  +coeff(  4)            *x41*x51
     5  +coeff(  5)    *x22*x31        
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)            *x43    
     8  +coeff(  8)*x11*x21*x31        
      y_e_q1ex    =y_e_q1ex    
     9  +coeff(  9)*x11*x21    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x41*x52
     4  +coeff( 13)*x12    *x31        
     5  +coeff( 14)*x12        *x41    
     6  +coeff( 15)        *x33    *x52
     7  +coeff( 16)*x11*x23    *x41    
c
      return
      end
      function p_e_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49999E-01,-0.59933E-01,-0.49999E-01,-0.39798E-01,-0.49962E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49871E-01, 0.59944E-01, 0.49999E-01, 0.39798E-01, 0.49911E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29974423E-01, 0.88498140E-01,-0.16159070E-02,-0.26054990E-02,
     + -0.46181760E-03,-0.36558634E-03,-0.40714670E-03,-0.38139830E-04,
     + -0.22437180E-03, 0.12091850E-03,-0.11691460E-04, 0.95883780E-04,
     + -0.49089670E-03,-0.37452951E-03,-0.52230473E-03,-0.27665100E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_e_q1ex    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31    *x51
     4  +coeff(  4)            *x41*x51
     5  +coeff(  5)    *x22    *x41    
     6  +coeff(  6)*x11*x21*x31        
     7  +coeff(  7)*x11*x21    *x41    
     8  +coeff(  8)    *x24*x31        
      p_e_q1ex    =p_e_q1ex    
     9  +coeff(  9)    *x22*x31        
     1  +coeff( 10)            *x41*x52
     2  +coeff( 11)*x12        *x41    
     3  +coeff( 12)        *x33    *x52
     4  +coeff( 13)*x11*x23    *x41    
     5  +coeff( 14)*x12*x22*x31        
     6  +coeff( 15)*x12*x22    *x41    
     7  +coeff( 16)*x11*x23*x33        
c
      return
      end
      function l_e_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.2177480E-02/
      data xmin/
     1 -0.49999E-01,-0.59933E-01,-0.49926E-01,-0.39798E-01,-0.49962E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49871E-01, 0.59944E-01, 0.49999E-01, 0.39594E-01, 0.49911E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19429341E-02,-0.34014100E-02,-0.22326310E-02,-0.53287832E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_e_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x22            
     3  +coeff(  3)        *x31*x41    
     4  +coeff(  4)            *x42    
c
      return
      end
      function x_e_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2874856E-02/
      data xmin/
     1 -0.49999E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49871E-01, 0.59944E-01, 0.49999E-01, 0.39017E-01, 0.49911E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25233451E-02,-0.12084104E+00,-0.42483003E-04, 0.53561640E-01,
     +  0.35550600E-02,-0.12656780E-03,-0.72390400E-02, 0.23614281E-05,
     + -0.36916723E-02, 0.63900330E-03, 0.60967891E-03, 0.20578900E-02,
     +  0.25585030E-02,-0.53806020E-02, 0.44937800E-03, 0.58186430E-03,
     +  0.95167820E-03, 0.83940272E-03, 0.10050940E-02, 0.10391260E-02,
     + -0.48442694E-06, 0.36566733E-03, 0.33098330E-02, 0.91023640E-03,
     + -0.19623510E-04, 0.12883171E-02,-0.28023780E-03, 0.30207360E-03,
     +  0.28035312E-02,-0.18483960E-03, 0.37129570E-02, 0.31478223E-02,
     +  0.31431741E-03, 0.11832542E-02, 0.96280004E-04, 0.77530000E-03,
     +  0.42620352E-02, 0.23424890E-02, 0.34652340E-03, 0.33656430E-02,
     +  0.80791683E-02, 0.47177923E-02, 0.79984860E-02, 0.84981741E-02,
     +  0.35976143E-02, 0.58548890E-02, 0.16988450E-02,-0.80807454E-03,
     + -0.60535220E-03,-0.64156510E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_e_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_e_dent    =x_e_dent    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)*x12                
      x_e_dent    =x_e_dent    
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)*x11    *x31*x41    
     2  +coeff( 20)*x11        *x42    
     3  +coeff( 21)            *x44    
     4  +coeff( 22)*x11            *x52
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)*x11*x22*x32        
     7  +coeff( 25)*x11    *x34        
     8  +coeff( 26)*x12*x21    *x42    
      x_e_dent    =x_e_dent    
     9  +coeff( 27)            *x42    
     1  +coeff( 28)*x11    *x32        
     2  +coeff( 29)    *x23*x31*x41    
     3  +coeff( 30)*x12            *x51
     4  +coeff( 31)*x11*x22*x31*x41    
     5  +coeff( 32)*x11*x22    *x42    
     6  +coeff( 33)*x12*x23            
     7  +coeff( 34)*x12*x21*x31*x41    
     8  +coeff( 35)*x11*x21        *x51
      x_e_dent    =x_e_dent    
     9  +coeff( 36)    *x23*x32        
     1  +coeff( 37)    *x21*x31*x43    
     2  +coeff( 38)    *x21    *x44    
     3  +coeff( 39)*x12*x21*x32        
     4  +coeff( 40)    *x23*x33*x41    
     5  +coeff( 41)    *x23*x32*x42    
     6  +coeff( 42)    *x21*x34*x42    
     7  +coeff( 43)    *x23*x31*x43    
     8  +coeff( 44)    *x21*x33*x43    
      x_e_dent    =x_e_dent    
     9  +coeff( 45)    *x23    *x44    
     1  +coeff( 46)    *x21*x32*x44    
     2  +coeff( 47)*x11*x24    *x42    
     3  +coeff( 48)*x12*x21*x31*x41*x51
     4  +coeff( 49)*x12*x21    *x42*x51
     5  +coeff( 50)    *x24            
c
c
c  transform to TRANSPORT x
c
      x_e_dent=-x_e_dent*3.346065 !multiply by -cos(30)/sin(15degrees)

      return
      end
      function t_e_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.3414553E-02/
      data xmin/
     1 -0.49999E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49871E-01, 0.59944E-01, 0.49999E-01, 0.39017E-01, 0.49911E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.31699314E-02, 0.55782970E-01, 0.57176413E-03,-0.32120022E-01,
     +  0.37012402E-02, 0.12714673E-03, 0.12832750E-02, 0.33851390E-02,
     +  0.31713962E-02,-0.35328774E-02,-0.55214011E-03,-0.56918873E-03,
     + -0.16593752E-02, 0.29888840E-02, 0.24384910E-03, 0.56448101E-03,
     + -0.29679110E-03, 0.72874600E-03,-0.30305044E-03,-0.52991113E-03,
     + -0.26682560E-03,-0.62241000E-03, 0.31943460E-03,-0.22649410E-03,
     + -0.16834923E-03, 0.69845010E-03,-0.33952333E-02,-0.15445360E-02,
     +  0.37475652E-03,-0.23353190E-02,-0.51089394E-02,-0.28393601E-02,
     + -0.17549200E-03, 0.36549930E-03, 0.32995911E-03,-0.20893630E-03,
     + -0.24659294E-03,-0.54905060E-03,-0.17805520E-03,-0.15295502E-02,
     + -0.16882151E-03,-0.25671991E-03,-0.11456620E-02,-0.56785362E-03,
     + -0.10977602E-02,-0.30354710E-02,-0.25786220E-02,-0.59414730E-03,
     + -0.15173793E-02,-0.34127733E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_e_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_e_dent    =t_e_dent    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)    *x21        *x52
      t_e_dent    =t_e_dent    
     9  +coeff( 18)*x12                
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)*x11    *x31*x41    
     3  +coeff( 21)*x11        *x42    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)*x11*x21        *x51
     6  +coeff( 24)*x11            *x52
     7  +coeff( 25)*x12*x21            
     8  +coeff( 26)*x11*x23            
      t_e_dent    =t_e_dent    
     9  +coeff( 27)    *x23*x31*x41    
     1  +coeff( 28)    *x21*x33*x41    
     2  +coeff( 29)*x11*x21    *x42    
     3  +coeff( 30)    *x23    *x42    
     4  +coeff( 31)    *x21*x31*x43    
     5  +coeff( 32)    *x21    *x44    
     6  +coeff( 33)*x12            *x51
     7  +coeff( 34)    *x22*x31*x41*x51
     8  +coeff( 35)        *x33*x41*x51
      t_e_dent    =t_e_dent    
     9  +coeff( 36)    *x21*x31*x41*x52
     1  +coeff( 37)    *x21    *x42*x52
     2  +coeff( 38)*x11*x22*x32        
     3  +coeff( 39)*x11    *x34        
     4  +coeff( 40)*x11*x22    *x42    
     5  +coeff( 41)*x12*x21*x32        
     6  +coeff( 42)    *x23*x34        
     7  +coeff( 43)*x12*x21*x31*x41    
     8  +coeff( 44)    *x23*x33*x41    
      t_e_dent    =t_e_dent    
     9  +coeff( 45)*x12*x21    *x42    
     1  +coeff( 46)    *x23*x32*x42    
     2  +coeff( 47)    *x21*x34*x42    
     3  +coeff( 48)*x12*x24            
     4  +coeff( 49)*x11*x22*x33*x41    
     5  +coeff( 50)                *x52
c
      return
      end
      function y_e_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49999E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49871E-01, 0.59944E-01, 0.49999E-01, 0.39017E-01, 0.49911E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35342101E-01, 0.16724030E+00, 0.37616072E-02, 0.47935120E-02,
     +  0.66122040E-02, 0.19396580E-01,-0.17115452E-02,-0.15339180E-02,
     + -0.28046173E-03,-0.21800524E-02,-0.32966890E-02,-0.41445720E-02,
     + -0.33107453E-02,-0.85640954E-03,-0.13779600E-02, 0.16162780E-03,
     +  0.49735410E-03, 0.77961693E-03,-0.34009520E-03,-0.51539680E-03,
     + -0.24144770E-02,-0.32172540E-04, 0.83655140E-05,-0.56629783E-05,
     + -0.57193342E-05, 0.24067014E-04,-0.84576400E-03,-0.56393881E-03,
     + -0.25767040E-02,-0.17083070E-02,-0.25782163E-03,-0.61278270E-03,
     +  0.38545273E-03,-0.14106894E-02,-0.81438160E-03,-0.27806654E-02,
     +  0.46769352E-03,-0.60153182E-03,-0.73768710E-03,-0.38490590E-02,
     + -0.68438181E-03, 0.25428071E-03, 0.18367564E-02,-0.45734184E-04,
     + -0.96687850E-03,-0.10024640E-02,-0.44662570E-03,-0.64379690E-03,
     +  0.68305950E-03, 0.17742160E-02,-0.16139943E-03, 0.19464944E-02,
     +  0.39310284E-03, 0.12085520E-02, 0.10053610E-02, 0.90616390E-03,
     + -0.38299760E-02, 0.16025262E-02,-0.28041270E-02, 0.24713980E-03,
     +  0.46055840E-03,-0.20859192E-02,-0.49049162E-03,-0.15697660E-03,
     + -0.40746101E-03,-0.46272664E-02,-0.87046650E-04,-0.12360782E-03,
     +  0.10361533E-03,-0.82317830E-03, 0.14019710E-02,-0.13942223E-03,
     + -0.54330330E-03,-0.88591300E-03,-0.58707700E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_e_dent    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      y_e_dent    =y_e_dent    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x31*x42    
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)            *x41*x52
     7  +coeff( 16)*x11*x21*x31        
     8  +coeff( 17)*x11*x21    *x41    
      y_e_dent    =y_e_dent    
     9  +coeff( 18)*x11        *x41*x51
     1  +coeff( 19)*x12    *x31        
     2  +coeff( 20)*x12        *x41    
     3  +coeff( 21)    *x22*x32*x41    
     4  +coeff( 22)        *x34*x41    
     5  +coeff( 23)    *x22*x31    *x52
     6  +coeff( 24)        *x33    *x52
     7  +coeff( 25)        *x32*x41*x52
     8  +coeff( 26)*x11    *x33    *x51
      y_e_dent    =y_e_dent    
     9  +coeff( 27)*x12*x22    *x41    
     1  +coeff( 28)    *x22*x34*x41    
     2  +coeff( 29)*x12*x24*x31        
     3  +coeff( 30)        *x32*x41    
     4  +coeff( 31)    *x21*x31    *x51
     5  +coeff( 32)        *x31    *x52
     6  +coeff( 33)*x11    *x31    *x51
     7  +coeff( 34)    *x24    *x41    
     8  +coeff( 35)*x11*x21    *x43*x51
      y_e_dent    =y_e_dent    
     9  +coeff( 36)*x12*x24    *x41    
     1  +coeff( 37)    *x23    *x41    
     2  +coeff( 38)    *x24*x31        
     3  +coeff( 39)*x11*x22    *x41    
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)*x11*x23*x31        
     6  +coeff( 42)*x12*x21    *x41    
     7  +coeff( 43)*x11*x21    *x43    
     8  +coeff( 44)*x12        *x41*x51
      y_e_dent    =y_e_dent    
     9  +coeff( 45)    *x24*x31*x42    
     1  +coeff( 46)    *x22*x33*x42    
     2  +coeff( 47)    *x24    *x43    
     3  +coeff( 48)*x11*x21*x31*x42*x51
     4  +coeff( 49)*x11*x23*x32*x41    
     5  +coeff( 50)*x11*x23*x31*x42    
     6  +coeff( 51)*x11*x21*x33*x42    
     7  +coeff( 52)*x11*x23    *x43    
     8  +coeff( 53)*x11*x21*x31*x44    
      y_e_dent    =y_e_dent    
     9  +coeff( 54)*x12    *x31*x42*x51
     1  +coeff( 55)*x12        *x43*x51
     2  +coeff( 56)    *x24*x34*x41    
     3  +coeff( 57)*x12*x22*x31*x42    
     4  +coeff( 58)    *x24*x33*x42    
     5  +coeff( 59)*x12*x22    *x43    
     6  +coeff( 60)*x11*x23*x32*x43    
     7  +coeff( 61)*x12    *x34*x41*x51
     8  +coeff( 62)*x12*x24*x32*x41    
      y_e_dent    =y_e_dent    
     9  +coeff( 63)*x12*x22*x32*x43    
     1  +coeff( 64)*x11*x22*x31        
     2  +coeff( 65)    *x22*x33        
     3  +coeff( 66)    *x22*x31*x42    
     4  +coeff( 67)*x11        *x41*x52
     5  +coeff( 68)*x11*x21*x33        
     6  +coeff( 69)    *x23*x33        
     7  +coeff( 70)*x11*x23    *x41    
     8  +coeff( 71)*x11*x21*x31*x42    
      y_e_dent    =y_e_dent    
     9  +coeff( 72)    *x24    *x41*x51
     1  +coeff( 73)*x12    *x32*x41    
     2  +coeff( 74)*x12    *x31*x42    
     3  +coeff( 75)*x12        *x43    
c
      return
      end
      function p_e_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49999E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49871E-01, 0.59944E-01, 0.49999E-01, 0.39017E-01, 0.49911E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.59034161E-01,-0.70599030E-01,-0.15304200E-01,-0.33376490E-01,
     +  0.80299163E-02, 0.19418172E-01, 0.87782703E-02, 0.17350020E-01,
     + -0.13830902E-01, 0.18362710E-03,-0.26504332E-02,-0.97379190E-03,
     + -0.17768964E-02, 0.68317750E-02,-0.48337250E-03, 0.17420850E-01,
     +  0.12076070E-02, 0.16719500E-03, 0.33663510E-03,-0.16157210E-02,
     + -0.66434670E-03, 0.11247001E-04,-0.34227620E-02,-0.20680201E-02,
     +  0.57507430E-04, 0.43470750E-03, 0.88458653E-03,-0.76538344E-04,
     +  0.45455510E-04,-0.19892151E-02, 0.14756900E-02,-0.54754554E-02,
     + -0.36184290E-02,-0.48629520E-02,-0.41747340E-02,-0.11017020E-01,
     + -0.58698502E-03,-0.19562890E-01,-0.27165060E-03, 0.21910380E-03,
     + -0.11438112E-01,-0.19797764E-02, 0.25657934E-03,-0.89548982E-03,
     +  0.14315161E-02,-0.12714640E-02, 0.11379063E-02, 0.32490730E-03,
     +  0.10601631E-02,-0.83809363E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_e_dent    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)*x11        *x41    
      p_e_dent    =p_e_dent    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)        *x31    *x52
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)*x11*x21*x31        
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)*x11*x21    *x41    
     8  +coeff( 17)    *x23    *x41    
      p_e_dent    =p_e_dent    
     9  +coeff( 18)    *x21*x32*x41    
     1  +coeff( 19)    *x21*x31*x42    
     2  +coeff( 20)*x11    *x31    *x51
     3  +coeff( 21)    *x22*x31    *x51
     4  +coeff( 22)        *x33    *x51
     5  +coeff( 23)*x11        *x41*x51
     6  +coeff( 24)    *x22    *x41*x51
     7  +coeff( 25)        *x32*x41*x51
     8  +coeff( 26)        *x31*x42*x51
      p_e_dent    =p_e_dent    
     9  +coeff( 27)            *x43*x51
     1  +coeff( 28)    *x21*x31    *x52
     2  +coeff( 29)    *x21    *x41*x52
     3  +coeff( 30)*x12    *x31        
     4  +coeff( 31)*x11*x22*x31        
     5  +coeff( 32)    *x24*x31        
     6  +coeff( 33)    *x22*x33        
     7  +coeff( 34)*x12        *x41    
     8  +coeff( 35)    *x24    *x41    
      p_e_dent    =p_e_dent    
     9  +coeff( 36)    *x22*x32*x41    
     1  +coeff( 37)        *x34*x41    
     2  +coeff( 38)    *x22*x31*x42    
     3  +coeff( 39)        *x33*x42    
     4  +coeff( 40)*x11        *x43    
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)        *x31*x44    
     7  +coeff( 43)    *x23*x31    *x51
     8  +coeff( 44)*x11*x21    *x41*x51
      p_e_dent    =p_e_dent    
     9  +coeff( 45)    *x23    *x41*x51
     1  +coeff( 46)    *x22*x31    *x52
     2  +coeff( 47)        *x33    *x52
     3  +coeff( 48)*x11        *x41*x52
     4  +coeff( 49)        *x32*x41*x52
     5  +coeff( 50)*x12*x21*x31        
c
      return
      end
      function l_e_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 29)
      data ncoeff/ 28/
      data avdat/  0.3626360E-03/
      data xmin/
     1 -0.49999E-01,-0.59666E-01,-0.49926E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49871E-01, 0.59944E-01, 0.49999E-01, 0.38045E-01, 0.49911E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27912494E-03, 0.23325632E+00,-0.10367600E+00, 0.14026960E-01,
     +  0.10400180E-01,-0.17262393E-01, 0.16072744E-01,-0.52257850E-02,
     + -0.98317800E-02, 0.73564784E-04,-0.15044030E-02,-0.67829443E-02,
     + -0.86661680E-02,-0.56482353E-04,-0.96706070E-03,-0.23986980E-02,
     + -0.98699010E-02,-0.17240390E-02, 0.11955470E-02,-0.11848030E-02,
     + -0.50764500E-02, 0.10477270E-02, 0.23649660E-03, 0.92756381E-03,
     + -0.30412630E-02,-0.22963662E-02,-0.25711422E-02,-0.55026950E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_e_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)*x11            *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)*x12                
      l_e_dent    =l_e_dent    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x33*x41    
     2  +coeff( 11)        *x32        
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21*x33*x41    
     6  +coeff( 15)*x11*x24            
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)    *x21*x31*x41    
      l_e_dent    =l_e_dent    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)*x12            *x51
     5  +coeff( 23)            *x41    
     6  +coeff( 24)        *x31*x41*x51
     7  +coeff( 25)*x11*x22            
     8  +coeff( 26)*x11    *x31*x41    
      l_e_dent    =l_e_dent    
     9  +coeff( 27)*x11        *x42    
     1  +coeff( 28)*x11            *x52
c
      return
      end
      function x_e_dext    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 71)
      data ncoeff/ 70/
      data avdat/  0.1110108E-01/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49656E-01, 0.59944E-01, 0.49999E-01, 0.39017E-01, 0.49820E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.99786580E-02, 0.34102730E+00, 0.13299680E+00,-0.19494394E+00,
     +  0.11681910E-01,-0.24411790E-02,-0.69915610E-02, 0.45143330E-01,
     + -0.34517934E-02,-0.87458550E-02, 0.15779644E-03,-0.10184620E-01,
     + -0.10036770E-01, 0.49190083E-02,-0.58382320E-04,-0.15691341E-02,
     +  0.17120610E-02,-0.36907070E-02,-0.55352990E-03, 0.36174384E-04,
     +  0.15912650E-02,-0.13264240E-02, 0.37837901E-03, 0.21401834E-03,
     + -0.12442660E-01,-0.14123584E-01,-0.60592780E-03,-0.92067920E-03,
     + -0.25165320E-02, 0.26429180E-03,-0.27866480E-02,-0.92987600E-03,
     + -0.59723730E-03,-0.10599454E-02,-0.15391281E-01,-0.59127871E-03,
     + -0.13899390E-01,-0.71777671E-03,-0.92284852E-03,-0.81224934E-04,
     + -0.96311332E-03, 0.18224400E-03,-0.11296860E-04, 0.23599683E-02,
     + -0.10613530E-01,-0.78935950E-03,-0.24805394E-01,-0.15506700E-01,
     +  0.46077792E-03,-0.15996800E-02, 0.19473693E-03,-0.13060792E-03,
     + -0.39421110E-02,-0.43243760E-02,-0.20391050E-01,-0.79594952E-02,
     + -0.11929050E-01,-0.16893100E-01,-0.22075381E-02,-0.86813590E-02,
     + -0.81064840E-03,-0.17362682E-02, 0.42708634E-03,-0.11541480E-01,
     +  0.15821232E-03,-0.12015960E-02,-0.27201450E-02, 0.33098060E-03,
     + -0.85397041E-03, 0.16835090E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_e_dext    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_e_dext    =x_e_dext    
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)*x12                
      x_e_dext    =x_e_dext    
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x22*x32        
     2  +coeff( 20)        *x34        
     3  +coeff( 21)*x11*x21        *x51
     4  +coeff( 22)    *x23*x32        
     5  +coeff( 23)    *x21*x34        
     6  +coeff( 24)    *x21*x32    *x52
     7  +coeff( 25)*x11*x22*x31*x41    
     8  +coeff( 26)*x11*x22    *x42    
      x_e_dext    =x_e_dext    
     9  +coeff( 27)*x12*x23            
     1  +coeff( 28)*x12*x21*x32        
     2  +coeff( 29)    *x23*x34        
     3  +coeff( 30)*x11*x22*x34        
     4  +coeff( 31)    *x21*x32        
     5  +coeff( 32)    *x21        *x52
     6  +coeff( 33)*x11    *x32        
     7  +coeff( 34)*x11    *x31*x41    
     8  +coeff( 35)    *x23    *x42    
      x_e_dext    =x_e_dext    
     9  +coeff( 36)*x12            *x51
     1  +coeff( 37)    *x23*x33*x41    
     2  +coeff( 38)        *x32        
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)        *x33*x41    
     5  +coeff( 41)*x11        *x42    
     6  +coeff( 42)*x11            *x52
     7  +coeff( 43)*x12*x21            
     8  +coeff( 44)*x11*x23            
      x_e_dext    =x_e_dext    
     9  +coeff( 45)    *x23*x31*x41    
     1  +coeff( 46)*x11*x21    *x42    
     2  +coeff( 47)    *x21*x31*x43    
     3  +coeff( 48)    *x21    *x44    
     4  +coeff( 49)        *x32*x42*x51
     5  +coeff( 50)*x12*x22            
     6  +coeff( 51)    *x24*x32        
     7  +coeff( 52)    *x24        *x52
     8  +coeff( 53)*x12*x21*x31*x41    
      x_e_dext    =x_e_dext    
     9  +coeff( 54)*x12*x21    *x42    
     1  +coeff( 55)    *x23*x32*x42    
     2  +coeff( 56)    *x21*x34*x42    
     3  +coeff( 57)    *x23*x31*x43    
     4  +coeff( 58)    *x21*x33*x43    
     5  +coeff( 59)    *x23    *x44    
     6  +coeff( 60)    *x21*x32*x44    
     7  +coeff( 61)*x11*x24*x32        
     8  +coeff( 62)    *x24            
      x_e_dext    =x_e_dext    
     9  +coeff( 63)    *x23        *x51
     1  +coeff( 64)    *x21*x32*x42    
     2  +coeff( 65)        *x32    *x53
     3  +coeff( 66)*x11*x24            
     4  +coeff( 67)*x11*x22*x32        
     5  +coeff( 68)*x12        *x42    
     6  +coeff( 69)*x11    *x31*x43    
     7  +coeff( 70)*x12            *x52
c
c
c  transform to TRANSPORT x
c
      x_e_dext=x_e_dext*0.866025
      return
      end
      function t_e_dext    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.7116606E-02/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49656E-01, 0.59944E-01, 0.49999E-01, 0.39017E-01, 0.49820E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.68454640E-02,-0.75662600E-01, 0.22821990E-01, 0.32990200E-01,
     + -0.36332320E-02,-0.49980410E-03, 0.16115524E-02,-0.23635310E-02,
     + -0.40929792E-02,-0.12402804E-02, 0.43888594E-02,-0.99546290E-04,
     +  0.39718850E-03, 0.29430620E-02, 0.49875250E-02,-0.33242523E-02,
     + -0.37521330E-03, 0.37132790E-03, 0.35105790E-03,-0.21101120E-02,
     +  0.44961630E-03,-0.13248561E-02, 0.14617491E-02, 0.38918300E-03,
     +  0.64985650E-03, 0.11145540E-02,-0.28060824E-03, 0.62849890E-03,
     +  0.17622992E-03,-0.14365620E-03, 0.29336880E-03, 0.27817740E-02,
     +  0.21814694E-02, 0.11669032E-03,-0.26509890E-03, 0.22550380E-04,
     +  0.54835260E-03, 0.13229511E-03,-0.47679140E-03, 0.21707754E-03,
     +  0.19940992E-03,-0.37958880E-03, 0.88387064E-03,-0.39739200E-03,
     + -0.65741730E-03,-0.28100342E-03, 0.57288490E-03, 0.59780053E-03,
     +  0.23668760E-03,-0.29964730E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_e_dext    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_e_dext    =t_e_dext    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x22        *x51
      t_e_dext    =t_e_dext    
     9  +coeff( 18)        *x32    *x51
     1  +coeff( 19)        *x31*x41*x51
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x12                
     5  +coeff( 23)*x11*x22            
     6  +coeff( 24)*x11    *x32        
     7  +coeff( 25)*x11    *x31*x41    
     8  +coeff( 26)    *x22    *x42    
      t_e_dext    =t_e_dext    
     9  +coeff( 27)*x11*x21        *x51
     1  +coeff( 28)    *x21    *x42*x51
     2  +coeff( 29)*x11            *x52
     3  +coeff( 30)*x12*x21            
     4  +coeff( 31)*x12            *x51
     5  +coeff( 32)*x11*x22*x31*x41    
     6  +coeff( 33)*x11*x22    *x42    
     7  +coeff( 34)*x11    *x33*x41*x51
     8  +coeff( 35)    *x22*x32        
      t_e_dext    =t_e_dext    
     9  +coeff( 36)    *x23        *x51
     1  +coeff( 37)    *x21*x31*x41*x51
     2  +coeff( 38)    *x22        *x52
     3  +coeff( 39)            *x42*x52
     4  +coeff( 40)*x11*x21*x32        
     5  +coeff( 41)*x11*x21*x31*x41    
     6  +coeff( 42)*x11*x21    *x42    
     7  +coeff( 43)    *x23    *x42    
     8  +coeff( 44)*x11*x22        *x51
      t_e_dext    =t_e_dext    
     9  +coeff( 45)*x11    *x31*x41*x51
     1  +coeff( 46)*x11        *x42*x51
     2  +coeff( 47)    *x22    *x42*x51
     3  +coeff( 48)*x11*x22*x32        
     4  +coeff( 49)*x12*x21        *x51
     5  +coeff( 50)        *x31*x41*x54
c
      return
      end
      function y_e_dext    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49656E-01, 0.59944E-01, 0.49999E-01, 0.39017E-01, 0.49820E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.53764500E-01, 0.13579890E+00,-0.19198070E-01,-0.68191360E-01,
     +  0.23486850E-01, 0.62972970E-01, 0.93214020E-02,-0.10074173E-01,
     +  0.31461420E-01,-0.32107870E-01,-0.55280392E-02,-0.12335190E-01,
     + -0.10352062E-01,-0.24714451E-02,-0.65535310E-02,-0.20833610E-02,
     + -0.45033474E-02, 0.62523020E-02, 0.11765343E-02, 0.41468522E-03,
     +  0.22100731E-01, 0.39300960E-02, 0.58204182E-02, 0.51213754E-02,
     +  0.16983120E-03,-0.20400881E-02,-0.38726230E-02, 0.10596770E-02,
     +  0.14953000E-02, 0.66530670E-03, 0.25302980E-03,-0.18471960E-02,
     + -0.55518560E-05,-0.18037292E-02,-0.83756490E-03,-0.55801840E-02,
     + -0.18023780E-02,-0.55613870E-02, 0.21240350E-03,-0.16010130E-01,
     + -0.16079140E-03,-0.12974630E-01,-0.32611723E-03,-0.48134270E-03,
     + -0.54591011E-04,-0.21364000E-03,-0.14220551E-02,-0.22849403E-02,
     + -0.53227671E-04, 0.42208653E-03, 0.69688220E-02, 0.14049700E-03,
     +  0.71615810E-04, 0.13515230E-02, 0.46105500E-03,-0.36730950E-02,
     + -0.63491730E-02,-0.12648740E-02,-0.42321020E-02,-0.37057930E-02,
     +  0.87875872E-02,-0.14284020E-04, 0.92909430E-02, 0.25466922E-03,
     + -0.21025983E-02, 0.60783330E-02,-0.85534460E-03, 0.25949780E-02,
     + -0.37050020E-03,-0.10911263E-02, 0.18612290E-02,-0.27416991E-02,
     +  0.14094010E-02, 0.65438272E-02, 0.38378760E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_e_dext    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      y_e_dext    =y_e_dext    
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)        *x32*x41    
     3  +coeff( 12)        *x31*x42    
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x21*x31    *x51
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)            *x41*x52
      y_e_dext    =y_e_dext    
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)        *x33    *x51
     8  +coeff( 26)*x11        *x41*x51
      y_e_dext    =y_e_dext    
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)        *x32*x41*x51
     2  +coeff( 29)            *x43*x51
     3  +coeff( 30)    *x21    *x41*x52
     4  +coeff( 31)            *x41*x53
     5  +coeff( 32)*x12    *x31        
     6  +coeff( 33)*x11*x22*x31        
     7  +coeff( 34)    *x24*x31        
     8  +coeff( 35)    *x22*x33        
      y_e_dext    =y_e_dext    
     9  +coeff( 36)*x12        *x41    
     1  +coeff( 37)*x11*x22    *x41    
     2  +coeff( 38)    *x22*x32*x41    
     3  +coeff( 39)*x11    *x31*x42    
     4  +coeff( 40)    *x22*x31*x42    
     5  +coeff( 41)        *x33*x42    
     6  +coeff( 42)    *x22    *x43    
     7  +coeff( 43)*x11*x21*x31    *x51
     8  +coeff( 44)*x11*x21    *x41*x51
      y_e_dext    =y_e_dext    
     9  +coeff( 45)        *x33    *x52
     1  +coeff( 46)*x11        *x41*x52
     2  +coeff( 47)*x11*x23*x31        
     3  +coeff( 48)*x11*x23    *x41    
     4  +coeff( 49)    *x21*x34*x41    
     5  +coeff( 50)    *x21*x33*x42    
     6  +coeff( 51)*x11*x21    *x43    
     7  +coeff( 52)*x11*x22*x31    *x51
     8  +coeff( 53)    *x24*x31    *x51
      y_e_dext    =y_e_dext    
     9  +coeff( 54)*x12        *x41*x51
     1  +coeff( 55)        *x31*x44*x51
     2  +coeff( 56)*x12*x22*x31        
     3  +coeff( 57)*x12*x22    *x41    
     4  +coeff( 58)    *x24*x32*x41    
     5  +coeff( 59)*x12    *x31*x42    
     6  +coeff( 60)*x12        *x43    
     7  +coeff( 61)*x11*x23*x31*x42    
     8  +coeff( 62)*x12*x21    *x43    
      y_e_dext    =y_e_dext    
     9  +coeff( 63)*x11*x23    *x43    
     1  +coeff( 64)*x12    *x31*x42*x51
     2  +coeff( 65)*x12    *x34*x41    
     3  +coeff( 66)*x11*x23*x34*x41    
     4  +coeff( 67)        *x33        
     5  +coeff( 68)    *x21*x32*x41    
     6  +coeff( 69)*x11    *x31    *x51
     7  +coeff( 70)    *x22*x31    *x51
     8  +coeff( 71)        *x31*x42*x51
      y_e_dext    =y_e_dext    
     9  +coeff( 72)    *x24    *x41    
     1  +coeff( 73)*x12*x21    *x41    
     2  +coeff( 74)*x11*x21*x31*x42    
     3  +coeff( 75)*x12    *x31    *x51
c
      return
      end
      function p_e_dext    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49656E-01, 0.59944E-01, 0.49999E-01, 0.39017E-01, 0.49820E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.19728220E-01, 0.72582820E-02,-0.35936710E-02,-0.15191650E-01,
     +  0.43401420E-02, 0.14257683E-01, 0.19975370E-02,-0.20226250E-02,
     +  0.67009330E-02,-0.90548992E-02,-0.21795793E-02,-0.21945731E-02,
     + -0.19519444E-03,-0.36546573E-02, 0.55452660E-04,-0.59271710E-03,
     +  0.61937020E-03, 0.29116243E-03,-0.28679740E-04, 0.65113240E-02,
     +  0.12735462E-02, 0.81536050E-04, 0.14678371E-03,-0.24507832E-03,
     +  0.39825100E-04, 0.21672830E-03,-0.14058450E-02, 0.99433470E-04,
     +  0.41246353E-03, 0.42829141E-03,-0.61753082E-04,-0.27641592E-03,
     + -0.75077921E-04,-0.18786034E-02,-0.11831743E-02, 0.85897001E-04,
     + -0.44765980E-03,-0.10285133E-02,-0.22338321E-02,-0.27482490E-03,
     + -0.27801380E-02,-0.14044650E-02,-0.18061590E-03, 0.22487820E-03,
     +  0.75190230E-03, 0.45483853E-03,-0.24787880E-03, 0.31151230E-03,
     + -0.21619800E-03, 0.37405063E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_e_dext    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      p_e_dext    =p_e_dext    
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x21*x31    *x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)*x11*x21*x31        
      p_e_dext    =p_e_dext    
     9  +coeff( 18)    *x23*x31        
     1  +coeff( 19)    *x21*x33        
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x21*x32*x41    
     5  +coeff( 23)    *x21    *x43    
     6  +coeff( 24)*x11    *x31    *x51
     7  +coeff( 25)        *x33    *x51
     8  +coeff( 26)*x11        *x41*x51
      p_e_dext    =p_e_dext    
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)        *x32*x41*x51
     2  +coeff( 29)        *x31*x42*x51
     3  +coeff( 30)            *x43*x51
     4  +coeff( 31)    *x21*x31    *x52
     5  +coeff( 32)*x11*x22*x31        
     6  +coeff( 33)    *x24*x31        
     7  +coeff( 34)*x12        *x41    
     8  +coeff( 35)*x11*x22    *x41    
      p_e_dext    =p_e_dext    
     9  +coeff( 36)    *x24    *x41    
     1  +coeff( 37)    *x22*x32*x41    
     2  +coeff( 38)        *x34*x41    
     3  +coeff( 39)        *x33*x42    
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)        *x32*x43    
     6  +coeff( 42)        *x31*x44    
     7  +coeff( 43)    *x21*x33    *x51
     8  +coeff( 44)*x11*x21    *x41*x51
      p_e_dext    =p_e_dext    
     9  +coeff( 45)    *x23    *x41*x51
     1  +coeff( 46)    *x21*x32*x41*x51
     2  +coeff( 47)        *x33    *x52
     3  +coeff( 48)    *x22    *x41*x52
     4  +coeff( 49)        *x32*x41*x52
     5  +coeff( 50)    *x21    *x41*x53
c
      return
      end
      function l_e_dext    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/ -0.2903154E-01/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49926E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49656E-01, 0.59944E-01, 0.49999E-01, 0.38045E-01, 0.49820E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27257550E-01,-0.48304301E+00,-0.97353780E-01, 0.25220254E+00,
     + -0.45421522E-01,-0.30781500E-01, 0.27138990E-01,-0.14777870E-01,
     + -0.76393290E-02, 0.29907453E-01, 0.99321771E-02, 0.18703493E-02,
     + -0.53144893E-02,-0.60146674E-02, 0.17131690E-02, 0.62809080E-02,
     +  0.24595394E-01,-0.27443354E-02,-0.16935200E-02, 0.16191592E-02,
     +  0.15276613E-02, 0.41472502E-02, 0.41780890E-02,-0.16757473E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_e_dext    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)*x11            *x51
      l_e_dext    =l_e_dext    
     9  +coeff(  9)*x12                
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)*x11*x22            
     3  +coeff( 12)    *x23*x31*x41    
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)    *x21*x31*x41    
      l_e_dext    =l_e_dext    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)        *x34        
     2  +coeff( 20)    *x23            
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)*x11        *x42    
     6  +coeff( 24)*x11*x21        *x51
c
      return
      end
      function x_e_q3en    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 71)
      data ncoeff/ 70/
      data avdat/  0.1595396E-02/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49656E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49820E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16519014E-02, 0.21749283E+00, 0.13823223E+00,-0.13503350E+00,
     +  0.18138520E-01,-0.23615392E-02,-0.63336040E-02, 0.34997832E-01,
     + -0.58776400E-02,-0.14813502E-01, 0.10949544E-02,-0.72269113E-02,
     + -0.10265183E-01, 0.14052081E-02, 0.26236362E-02, 0.37224720E-04,
     + -0.11855470E-02, 0.30548204E-02,-0.48988931E-02,-0.37633570E-03,
     + -0.55971170E-03, 0.16091260E-04,-0.10396054E-02,-0.96510850E-04,
     + -0.19412644E-03,-0.92452640E-03, 0.56732650E-03,-0.10021084E-01,
     +  0.10742950E-03,-0.72289500E-03,-0.78908511E-03,-0.64139644E-03,
     + -0.27387550E-02, 0.36064183E-03, 0.28322581E-03, 0.78652013E-03,
     + -0.10114030E-01,-0.15280770E-02,-0.31905970E-03,-0.52451720E-02,
     + -0.19297020E-02,-0.22153300E-02,-0.45745040E-02,-0.19725824E-02,
     + -0.33779610E-03,-0.15405882E-02,-0.14725553E-02, 0.51462880E-03,
     + -0.39159730E-03, 0.19596074E-02,-0.15032910E-01,-0.79016843E-02,
     + -0.13252740E-02,-0.42869660E-02, 0.54241920E-03, 0.27101370E-03,
     + -0.75293260E-03, 0.66324850E-03,-0.91537900E-02,-0.14994230E-01,
     + -0.24114912E-01,-0.16033532E-01,-0.22728922E-02,-0.49768190E-02,
     +  0.69678481E-02,-0.10678150E-01,-0.13211560E-01, 0.14446572E-03,
     + -0.38427620E-03, 0.32270312E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_e_q3en    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_e_q3en    =x_e_q3en    
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)            *x42*x51
      x_e_q3en    =x_e_q3en    
     9  +coeff( 18)*x12                
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)*x11    *x32        
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)        *x34        
     5  +coeff( 23)*x11    *x31*x41    
     6  +coeff( 24)*x11        *x42    
     7  +coeff( 25)*x11*x21        *x51
     8  +coeff( 26)    *x23*x32        
      x_e_q3en    =x_e_q3en    
     9  +coeff( 27)    *x21*x34        
     1  +coeff( 28)    *x23    *x42    
     2  +coeff( 29)    *x21*x32    *x52
     3  +coeff( 30)    *x22*x33*x41    
     4  +coeff( 31)    *x23*x34        
     5  +coeff( 32)        *x32        
     6  +coeff( 33)    *x21*x32        
     7  +coeff( 34)                *x53
     8  +coeff( 35)*x11            *x52
      x_e_q3en    =x_e_q3en    
     9  +coeff( 36)*x12*x21            
     1  +coeff( 37)    *x23*x31*x41    
     2  +coeff( 38)*x11*x21    *x42    
     3  +coeff( 39)*x12            *x51
     4  +coeff( 40)*x11*x22    *x42    
     5  +coeff( 41)*x12*x21*x31*x41    
     6  +coeff( 42)*x12*x21    *x42    
     7  +coeff( 43)*x11*x24*x31*x41    
     8  +coeff( 44)*x12*x23*x32        
      x_e_q3en    =x_e_q3en    
     9  +coeff( 45)    *x21        *x52
     1  +coeff( 46)    *x24            
     2  +coeff( 47)    *x22*x31*x41    
     3  +coeff( 48)    *x23        *x51
     4  +coeff( 49)    *x21    *x42*x51
     5  +coeff( 50)*x11*x23            
     6  +coeff( 51)    *x21*x31*x43    
     7  +coeff( 52)    *x21    *x44    
     8  +coeff( 53)*x12*x22            
      x_e_q3en    =x_e_q3en    
     9  +coeff( 54)*x11*x22*x31*x41    
     1  +coeff( 55)*x12        *x42    
     2  +coeff( 56)    *x22*x32*x42    
     3  +coeff( 57)*x11*x21    *x42*x51
     4  +coeff( 58)    *x21*x32*x42*x51
     5  +coeff( 59)    *x23*x32*x42    
     6  +coeff( 60)    *x21*x34*x42    
     7  +coeff( 61)    *x21*x33*x43    
     8  +coeff( 62)    *x21*x32*x44    
      x_e_q3en    =x_e_q3en    
     9  +coeff( 63)*x11*x24*x32        
     1  +coeff( 64)*x11*x24    *x42    
     2  +coeff( 65)    *x23*x34*x42    
     3  +coeff( 66)*x12*x23*x33*x41    
     4  +coeff( 67)*x12*x23*x32*x42    
     5  +coeff( 68)        *x31*x41*x51
     6  +coeff( 69)        *x32*x42    
     7  +coeff( 70)            *x44    
c
      return
      end
      function t_e_q3en    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.6288869E-02/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49656E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49820E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.63027041E-02,-0.75843030E-01, 0.22741720E-01, 0.33414531E-01,
     + -0.49983533E-02, 0.12350831E-03,-0.92987180E-03,-0.52838930E-02,
     + -0.14230100E-02, 0.59034210E-02, 0.68071210E-04, 0.59808100E-03,
     +  0.26538210E-02, 0.25358290E-03,-0.26481420E-02,-0.36160604E-03,
     + -0.14498483E-03, 0.76097314E-03, 0.53654010E-03,-0.18056170E-02,
     +  0.92432054E-03, 0.99415704E-03, 0.87917764E-03,-0.64248740E-04,
     + -0.67105270E-04,-0.52156850E-03,-0.66609622E-03, 0.32759650E-02,
     +  0.41933762E-03, 0.30697302E-02, 0.28215430E-02, 0.23670320E-03,
     +  0.25039681E-03, 0.93339284E-03, 0.72388153E-03, 0.13320210E-03,
     +  0.18564333E-03, 0.27451082E-03,-0.71748970E-03,-0.83323150E-03,
     +  0.20949900E-03, 0.23756930E-02, 0.52129374E-02, 0.29905410E-02,
     +  0.21469660E-03,-0.22615450E-03, 0.20752360E-02, 0.52407550E-03,
     +  0.31158482E-03, 0.50876090E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_e_q3en    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)    *x21        *x51
      t_e_q3en    =t_e_q3en    
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)        *x32    *x51
      t_e_q3en    =t_e_q3en    
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)*x12                
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)*x11        *x42    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)        *x32*x42    
     7  +coeff( 25)            *x44    
     8  +coeff( 26)*x11*x21        *x51
      t_e_q3en    =t_e_q3en    
     9  +coeff( 27)*x11*x21    *x42    
     1  +coeff( 28)    *x21*x32*x42    
     2  +coeff( 29)*x12            *x51
     3  +coeff( 30)*x11*x22*x31*x41    
     4  +coeff( 31)*x11*x22    *x42    
     5  +coeff( 32)        *x32*x44    
     6  +coeff( 33)*x12*x21*x32        
     7  +coeff( 34)*x11*x24*x32        
     8  +coeff( 35)            *x42    
      t_e_q3en    =t_e_q3en    
     9  +coeff( 36)                *x53
     1  +coeff( 37)*x11    *x32        
     2  +coeff( 38)*x11    *x31*x41    
     3  +coeff( 39)    *x21*x31*x41*x51
     4  +coeff( 40)    *x21    *x42*x51
     5  +coeff( 41)*x11*x23            
     6  +coeff( 42)    *x23    *x42    
     7  +coeff( 43)    *x21*x31*x43    
     8  +coeff( 44)    *x21    *x44    
      t_e_q3en    =t_e_q3en    
     9  +coeff( 45)*x11    *x31*x41*x51
     1  +coeff( 46)    *x24*x31*x41    
     2  +coeff( 47)    *x23*x33*x41    
     3  +coeff( 48)*x12*x21    *x42    
     4  +coeff( 49)    *x24*x34        
     5  +coeff( 50)    *x22*x31*x41    
c
      return
      end
      function y_e_q3en    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49656E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49820E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.72916142E-01, 0.14839620E+00,-0.20870130E-01,-0.85577880E-01,
     +  0.29649560E-01, 0.77998824E-01, 0.10028400E-01,-0.11879642E-01,
     + -0.11467423E-02, 0.38998280E-01,-0.39111942E-01,-0.16433780E-01,
     + -0.13690390E-01,-0.11308750E-01,-0.25893112E-02,-0.61044020E-02,
     +  0.70292633E-02, 0.16751580E-02, 0.59016830E-03, 0.27667030E-01,
     +  0.73217451E-02, 0.32120500E-02, 0.79354653E-02, 0.65456310E-02,
     + -0.33713630E-03, 0.31847183E-03,-0.10806410E-02,-0.41034090E-02,
     +  0.12909530E-02, 0.16191974E-02,-0.22831351E-02,-0.13147770E-02,
     + -0.82434960E-02,-0.52428813E-02,-0.64960130E-02,-0.44639222E-03,
     + -0.18339972E-01,-0.95000510E-03,-0.15296204E-01,-0.34948430E-02,
     + -0.16268430E-02, 0.14000390E-03, 0.63091813E-03, 0.54993113E-03,
     + -0.10605084E-02, 0.17487542E-03,-0.27396110E-02,-0.16862082E-03,
     +  0.35127964E-02,-0.35453150E-02,-0.44454063E-04, 0.40405160E-02,
     +  0.25044832E-03, 0.12836200E-02,-0.38951670E-03,-0.11234564E-02,
     + -0.11071254E-02,-0.40418910E-02,-0.95611670E-03,-0.66687972E-02,
     + -0.16256290E-02, 0.20697510E-02,-0.66952990E-03, 0.23153174E-01,
     +  0.18228722E-01, 0.37080240E-02, 0.32681570E-02, 0.76669820E-02,
     + -0.66839074E-02,-0.31987030E-02,-0.15644440E-02, 0.22735884E-02,
     + -0.31997180E-02, 0.14678960E-02, 0.14224200E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_e_q3en    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      y_e_q3en    =y_e_q3en    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x31*x42    
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)*x11*x21*x31        
      y_e_q3en    =y_e_q3en    
     9  +coeff( 18)    *x23*x31        
     1  +coeff( 19)    *x21*x33        
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x21*x32*x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)        *x33    *x51
      y_e_q3en    =y_e_q3en    
     9  +coeff( 27)*x11        *x41*x51
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)        *x32*x41*x51
     3  +coeff( 30)            *x43*x51
     4  +coeff( 31)*x12    *x31        
     5  +coeff( 32)    *x24*x31        
     6  +coeff( 33)*x12        *x41    
     7  +coeff( 34)*x11*x22    *x41    
     8  +coeff( 35)    *x22*x32*x41    
      y_e_q3en    =y_e_q3en    
     9  +coeff( 36)        *x34*x41    
     1  +coeff( 37)    *x22*x31*x42    
     2  +coeff( 38)        *x33*x42    
     3  +coeff( 39)    *x22    *x43    
     4  +coeff( 40)        *x32*x43    
     5  +coeff( 41)        *x31*x44    
     6  +coeff( 42)*x11*x21*x31    *x51
     7  +coeff( 43)    *x23*x31    *x51
     8  +coeff( 44)    *x21*x33    *x51
      y_e_q3en    =y_e_q3en    
     9  +coeff( 45)*x11*x21    *x41*x51
     1  +coeff( 46)    *x21*x31    *x53
     2  +coeff( 47)*x11*x23*x31        
     3  +coeff( 48)*x11*x21*x33        
     4  +coeff( 49)*x12*x21    *x41    
     5  +coeff( 50)*x11*x23    *x41    
     6  +coeff( 51)    *x21*x34*x41    
     7  +coeff( 52)*x11*x21    *x43    
     8  +coeff( 53)    *x24*x31    *x51
      y_e_q3en    =y_e_q3en    
     9  +coeff( 54)*x12        *x41*x51
     1  +coeff( 55)    *x22*x32*x41*x51
     2  +coeff( 56)*x11    *x31*x42*x51
     3  +coeff( 57)        *x33*x42*x51
     4  +coeff( 58)*x12*x22*x31        
     5  +coeff( 59)    *x24*x33        
     6  +coeff( 60)*x12*x22    *x41    
     7  +coeff( 61)    *x22*x34*x41    
     8  +coeff( 62)        *x34*x43    
      y_e_q3en    =y_e_q3en    
     9  +coeff( 63)    *x23*x33    *x51
     1  +coeff( 64)*x11*x23*x31*x42    
     2  +coeff( 65)*x11*x23    *x43    
     3  +coeff( 66)*x12    *x31*x42*x51
     4  +coeff( 67)*x12        *x43*x51
     5  +coeff( 68)*x11*x23*x34*x41    
     6  +coeff( 69)        *x32*x41    
     7  +coeff( 70)    *x21*x31    *x51
     8  +coeff( 71)    *x22*x31    *x51
      y_e_q3en    =y_e_q3en    
     9  +coeff( 72)        *x31*x42*x51
     1  +coeff( 73)    *x24    *x41    
     2  +coeff( 74)    *x23    *x41*x51
     3  +coeff( 75)    *x21*x31*x42*x51
c
      return
      end
      function p_e_q3en    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49656E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49820E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.19616674E-01, 0.13489820E-01,-0.54806590E-02,-0.15902681E-01,
     +  0.47048162E-02, 0.16325220E-01, 0.28256340E-02,-0.27996364E-02,
     + -0.21531720E-03, 0.66884000E-02,-0.10775683E-01,-0.15569243E-02,
     + -0.35613672E-02,-0.31734944E-02,-0.22991130E-02,-0.68307200E-04,
     + -0.58315600E-03, 0.13694011E-02, 0.23639630E-03, 0.23110750E-03,
     +  0.77786590E-02, 0.14303700E-02, 0.13714940E-02, 0.89336550E-03,
     + -0.42584200E-03, 0.67573290E-04,-0.37369310E-03,-0.11277800E-02,
     +  0.17031461E-03, 0.57912763E-03, 0.15458890E-04, 0.30013810E-05,
     + -0.51072851E-03, 0.16113380E-03, 0.12287790E-03,-0.25118393E-02,
     + -0.12325420E-02, 0.16042524E-03,-0.24971013E-03,-0.24073400E-03,
     +  0.50346023E-03, 0.36554920E-03,-0.76887554E-05,-0.75448834E-03,
     +  0.78960641E-03,-0.46750364E-03, 0.58129160E-03, 0.30036410E-03,
     +  0.15139511E-02, 0.15392031E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_e_q3en    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      p_e_q3en    =p_e_q3en    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)            *x41*x52
      p_e_q3en    =p_e_q3en    
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)        *x33    *x51
      p_e_q3en    =p_e_q3en    
     9  +coeff( 27)*x11        *x41*x51
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)        *x32*x41*x51
     3  +coeff( 30)        *x31*x42*x51
     4  +coeff( 31)            *x43*x51
     5  +coeff( 32)        *x31    *x53
     6  +coeff( 33)*x12    *x31        
     7  +coeff( 34)*x11*x22*x31        
     8  +coeff( 35)    *x24*x31        
      p_e_q3en    =p_e_q3en    
     9  +coeff( 36)*x12        *x41    
     1  +coeff( 37)*x11*x22    *x41    
     2  +coeff( 38)*x11*x21*x31    *x51
     3  +coeff( 39)    *x23*x31    *x51
     4  +coeff( 40)    *x21*x33    *x51
     5  +coeff( 41)    *x23    *x41*x51
     6  +coeff( 42)    *x22    *x41*x52
     7  +coeff( 43)*x12*x21*x31        
     8  +coeff( 44)*x11*x23*x31        
      p_e_q3en    =p_e_q3en    
     9  +coeff( 45)*x12*x21    *x41    
     1  +coeff( 46)*x11*x23    *x41    
     2  +coeff( 47)    *x23*x32*x41    
     3  +coeff( 48)    *x21*x34*x41    
     4  +coeff( 49)*x11*x21*x31*x42    
     5  +coeff( 50)*x11*x21    *x43    
c
      return
      end
      function l_e_q3en    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.2068126E-01/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49900E-01,-0.39017E-01,-0.49904E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49656E-01, 0.59522E-01, 0.49999E-01, 0.38045E-01, 0.49820E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.20679340E-01,-0.31027040E+00, 0.15443252E+00,-0.31112690E-01,
     + -0.28305811E-01,-0.20879773E-01, 0.25445750E-01,-0.13035272E-01,
     + -0.73870720E-02,-0.65434740E-02,-0.90784570E-02, 0.13647180E-01,
     +  0.69216573E-02,-0.81310780E-04, 0.73897130E-02,-0.22086561E-02,
     +  0.39317221E-02, 0.11907452E-01,-0.29069520E-02, 0.10529650E-01,
     +  0.16557694E-02, 0.14443000E-02, 0.10161491E-02, 0.31574831E-02,
     +  0.27402260E-02,-0.77327690E-03, 0.86968584E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_e_q3en    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)*x11            *x51
      l_e_q3en    =l_e_q3en    
     9  +coeff(  9)*x12                
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x22*x32        
     6  +coeff( 15)    *x23*x31*x41    
     7  +coeff( 16)        *x32        
     8  +coeff( 17)    *x21*x32        
      l_e_q3en    =l_e_q3en    
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)    *x23    *x42    
     3  +coeff( 21)        *x31*x41*x51
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)*x11    *x31*x41    
     7  +coeff( 25)*x11        *x42    
     8  +coeff( 26)*x11*x21        *x51
      l_e_q3en    =l_e_q3en    
     9  +coeff( 27)*x11            *x52
c
      return
      end
      function x_e_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 71)
      data ncoeff/ 70/
      data avdat/ -0.1220197E-01/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49382E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49548E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49777E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13774750E-01, 0.69178070E-01, 0.30952391E+00,-0.94327840E-01,
     +  0.11229442E-01,-0.78155440E-03,-0.70269880E-02,-0.80083890E-02,
     +  0.32798603E-01,-0.20618270E-01,-0.20648420E-02,-0.14661070E-01,
     + -0.42677730E-02, 0.25802000E-02, 0.17545121E-03,-0.44622190E-03,
     +  0.21291212E-02,-0.16168003E-02,-0.35516852E-02,-0.10998991E-02,
     + -0.17319520E-02, 0.45465080E-02, 0.12600452E-03,-0.15927270E-02,
     + -0.16242980E-02,-0.30383610E-02, 0.73221680E-03, 0.60234391E-02,
     + -0.58413700E-03, 0.11518553E-02,-0.74018610E-02,-0.33079711E-02,
     + -0.65905690E-02,-0.66497894E-02,-0.10220660E-01,-0.13258340E-01,
     + -0.63591571E-02, 0.70022960E-03,-0.64999563E-03,-0.42472891E-02,
     + -0.25858062E-04,-0.24098690E-03, 0.11994871E-02,-0.16607300E-02,
     + -0.37191461E-02, 0.12524100E-02,-0.36635983E-02, 0.55301110E-02,
     +  0.62375911E-03,-0.19990040E-02,-0.85415801E-03, 0.16167020E-02,
     + -0.25451313E-02, 0.28814920E-03, 0.11486682E-02,-0.78727720E-04,
     +  0.53106620E-03, 0.10661330E-02,-0.14029210E-02,-0.40316150E-03,
     +  0.69948873E-03,-0.81440420E-03, 0.28701790E-02,-0.51073631E-03,
     + -0.24056690E-02,-0.49408640E-03, 0.17014733E-02,-0.99533890E-03,
     + -0.14977430E-02, 0.38046650E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_e_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      x_e_q3ex    =x_e_q3ex    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)            *x42*x51
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)                *x53
      x_e_q3ex    =x_e_q3ex    
     9  +coeff( 18)*x12                
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)*x11    *x31*x41    
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)*x11        *x42    
     5  +coeff( 23)        *x32*x42    
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)    *x21*x31*x41*x51
     8  +coeff( 26)    *x21    *x42*x51
      x_e_q3ex    =x_e_q3ex    
     9  +coeff( 27)*x11            *x52
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)    *x21*x34        
     3  +coeff( 30)*x11*x21*x31*x41    
     4  +coeff( 31)    *x23*x31*x41    
     5  +coeff( 32)    *x21*x33*x41    
     6  +coeff( 33)*x11*x21    *x42    
     7  +coeff( 34)    *x23    *x42    
     8  +coeff( 35)    *x21*x32*x42    
      x_e_q3ex    =x_e_q3ex    
     9  +coeff( 36)    *x21*x31*x43    
     1  +coeff( 37)    *x21    *x44    
     2  +coeff( 38)*x12            *x51
     3  +coeff( 39)        *x34    *x51
     4  +coeff( 40)*x12*x22            
     5  +coeff( 41)    *x24*x32        
     6  +coeff( 42)    *x22    *x44    
     7  +coeff( 43)    *x23            
     8  +coeff( 44)    *x21*x32        
      x_e_q3ex    =x_e_q3ex    
     9  +coeff( 45)    *x21*x31*x41    
     1  +coeff( 46)        *x31*x41*x51
     2  +coeff( 47)    *x24            
     3  +coeff( 48)    *x22    *x42    
     4  +coeff( 49)*x12*x21            
     5  +coeff( 50)    *x23*x32        
     6  +coeff( 51)*x12    *x31*x41    
     7  +coeff( 52)*x12        *x42    
     8  +coeff( 53)*x12*x21    *x44    
      x_e_q3ex    =x_e_q3ex    
     9  +coeff( 54)*x11    *x32        
     1  +coeff( 55)    *x21*x32    *x51
     2  +coeff( 56)            *x42*x52
     3  +coeff( 57)    *x21        *x53
     4  +coeff( 58)*x11*x22        *x51
     5  +coeff( 59)    *x24        *x51
     6  +coeff( 60)*x11    *x32    *x51
     7  +coeff( 61)*x11    *x31*x41*x51
     8  +coeff( 62)    *x22*x31*x41*x51
      x_e_q3ex    =x_e_q3ex    
     9  +coeff( 63)    *x22    *x42*x51
     1  +coeff( 64)*x11            *x53
     2  +coeff( 65)*x11*x22    *x42    
     3  +coeff( 66)*x12*x21        *x51
     4  +coeff( 67)*x11*x23        *x51
     5  +coeff( 68)    *x23*x32    *x51
     6  +coeff( 69)*x11*x21    *x42*x51
     7  +coeff( 70)        *x34    *x52
c
      return
      end
      function t_e_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.7701860E-02/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49382E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49548E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49777E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.81624920E-02,-0.27714410E-01, 0.97573430E-01,-0.51057580E-02,
     +  0.30727893E-04,-0.12439821E-03,-0.25775381E-02,-0.17348903E-02,
     +  0.41777333E-02,-0.89875692E-02, 0.30640074E-02, 0.27284133E-03,
     + -0.50389630E-03, 0.14843401E-03,-0.40131030E-02,-0.15779170E-02,
     + -0.40272670E-03, 0.40918382E-03,-0.42735640E-03, 0.47065393E-03,
     + -0.15240182E-02,-0.84575740E-04,-0.21519810E-03, 0.21685630E-02,
     +  0.20976440E-02,-0.22101250E-03,-0.87002932E-03, 0.55685971E-03,
     +  0.18653803E-02, 0.43888210E-03,-0.23799524E-02, 0.33147250E-03,
     +  0.21117033E-03,-0.12708164E-02,-0.56087710E-04, 0.61456710E-03,
     + -0.78197130E-03, 0.64816100E-03,-0.10114333E-02, 0.32941412E-03,
     +  0.28489860E-03, 0.44300330E-03,-0.20771410E-03,-0.32120724E-03,
     + -0.62754670E-03,-0.18258050E-03, 0.12667630E-03,-0.36586984E-03,
     +  0.15180680E-03,-0.35281880E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_e_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_e_q3ex    =t_e_q3ex    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)        *x32    *x51
      t_e_q3ex    =t_e_q3ex    
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)*x12                
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)    *x22*x31*x41    
     6  +coeff( 24)*x11        *x42    
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)*x11*x21        *x51
      t_e_q3ex    =t_e_q3ex    
     9  +coeff( 27)    *x21    *x42*x51
     1  +coeff( 28)*x11            *x52
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)*x11*x21*x31*x41    
     4  +coeff( 31)*x11*x21    *x42    
     5  +coeff( 32)*x12            *x51
     6  +coeff( 33)        *x33*x41*x51
     7  +coeff( 34)*x12*x22            
     8  +coeff( 35)*x11    *x34        
      t_e_q3ex    =t_e_q3ex    
     9  +coeff( 36)*x12        *x42    
     1  +coeff( 37)    *x21*x33*x41*x51
     2  +coeff( 38)        *x31*x41*x51
     3  +coeff( 39)    *x24            
     4  +coeff( 40)*x11    *x32        
     5  +coeff( 41)        *x31*x41*x52
     6  +coeff( 42)*x11*x22        *x51
     7  +coeff( 43)*x11*x22            
     8  +coeff( 44)    *x23        *x51
      t_e_q3ex    =t_e_q3ex    
     9  +coeff( 45)    *x21*x31*x41*x51
     1  +coeff( 46)    *x22        *x52
     2  +coeff( 47)        *x32    *x52
     3  +coeff( 48)            *x42*x52
     4  +coeff( 49)*x12*x21            
     5  +coeff( 50)    *x23*x31*x41    
c
      return
      end
      function y_e_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49382E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49548E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49777E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.78569680E-01, 0.10022683E+00,-0.21710540E-01,-0.76580610E-01,
     +  0.22258980E-01, 0.77015243E-01, 0.10915990E-01,-0.10988423E-01,
     + -0.10554260E-02, 0.33390830E-01,-0.44704213E-01,-0.76790280E-02,
     + -0.15801023E-01,-0.13875892E-01,-0.14062263E-01, 0.38273380E-03,
     + -0.29261080E-02, 0.46472614E-02, 0.11620844E-02, 0.46873441E-03,
     +  0.30276460E-01, 0.76571302E-02, 0.59965660E-02, 0.49887150E-02,
     + -0.13479430E-02, 0.10724790E-02,-0.66065930E-02, 0.15803270E-02,
     +  0.10161731E-03,-0.76738010E-03,-0.17457700E-02, 0.18692961E-03,
     + -0.11494084E-02,-0.92428030E-02,-0.72142160E-02, 0.62064410E-03,
     + -0.56230680E-04,-0.24057290E-03, 0.15704903E-02,-0.13641130E-03,
     +  0.46175280E-02,-0.23639314E-03,-0.59528084E-03,-0.31019572E-03,
     +  0.92174643E-02, 0.94643132E-02,-0.36910432E-02, 0.14366230E-02,
     + -0.23010070E-02,-0.50538270E-02,-0.37750110E-02,-0.54438790E-02,
     + -0.22650720E-02, 0.42763080E-03, 0.10181690E-02,-0.72915432E-03,
     + -0.14524920E-02, 0.28678860E-02,-0.33339890E-03, 0.60389670E-02,
     + -0.23319460E-02, 0.29656490E-02,-0.31837960E-03, 0.64140610E-03,
     +  0.81393340E-03,-0.64004590E-02,-0.45025860E-02, 0.13784783E-02,
     +  0.29125270E-02, 0.22107490E-02, 0.39822430E-03, 0.13817811E-02,
     + -0.59957320E-03, 0.60094130E-03, 0.90991833E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_e_q3ex    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      y_e_q3ex    =y_e_q3ex    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)            *x41*x52
      y_e_q3ex    =y_e_q3ex    
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)*x11        *x41*x51
      y_e_q3ex    =y_e_q3ex    
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)        *x31*x42*x51
     2  +coeff( 29)            *x43*x51
     3  +coeff( 30)            *x41*x53
     4  +coeff( 31)*x12    *x31        
     5  +coeff( 32)*x11*x22*x31        
     6  +coeff( 33)    *x24*x31        
     7  +coeff( 34)*x12        *x41    
     8  +coeff( 35)*x11*x22    *x41    
      y_e_q3ex    =y_e_q3ex    
     9  +coeff( 36)*x11        *x43    
     1  +coeff( 37)    *x23*x31    *x51
     2  +coeff( 38)    *x21*x33    *x51
     3  +coeff( 39)    *x21    *x43*x51
     4  +coeff( 40)*x11*x23*x31        
     5  +coeff( 41)*x12*x21    *x41    
     6  +coeff( 42)*x11*x23    *x41    
     7  +coeff( 43)    *x23*x32*x41    
     8  +coeff( 44)    *x21*x34*x41    
      y_e_q3ex    =y_e_q3ex    
     9  +coeff( 45)*x11*x21*x31*x42    
     1  +coeff( 46)*x11*x21    *x43    
     2  +coeff( 47)    *x24*x31    *x51
     3  +coeff( 48)*x12        *x41*x51
     4  +coeff( 49)*x11*x22    *x41*x51
     5  +coeff( 50)    *x22*x31    *x53
     6  +coeff( 51)*x12*x22*x31        
     7  +coeff( 52)*x12*x22    *x41    
     8  +coeff( 53)    *x24    *x43    
      y_e_q3ex    =y_e_q3ex    
     9  +coeff( 54)*x11*x23*x31    *x51
     1  +coeff( 55)    *x23*x33    *x51
     2  +coeff( 56)    *x22*x31*x42*x52
     3  +coeff( 57)    *x22    *x43*x52
     4  +coeff( 58)*x11*x21*x34*x41    
     5  +coeff( 59)    *x23*x34*x41    
     6  +coeff( 60)    *x24*x31    *x53
     7  +coeff( 61)    *x21*x31    *x51
     8  +coeff( 62)    *x21*x32*x41    
      y_e_q3ex    =y_e_q3ex    
     9  +coeff( 63)    *x21*x31    *x52
     1  +coeff( 64)    *x22*x33        
     2  +coeff( 65)*x11    *x31*x42    
     3  +coeff( 66)    *x22*x31*x42    
     4  +coeff( 67)    *x22    *x43    
     5  +coeff( 68)*x11*x21*x31    *x51
     6  +coeff( 69)    *x23    *x41*x51
     7  +coeff( 70)    *x21*x31*x42*x51
     8  +coeff( 71)*x11    *x31    *x52
      y_e_q3ex    =y_e_q3ex    
     9  +coeff( 72)    *x22    *x41*x52
     1  +coeff( 73)    *x21    *x41*x53
     2  +coeff( 74)        *x34*x41*x51
     3  +coeff( 75)    *x23    *x41*x52
c
      return
      end
      function p_e_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49382E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49548E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49777E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15683910E-01,-0.45329172E-01, 0.47240783E-02, 0.21510554E-01,
     + -0.96127754E-02,-0.16712430E-01,-0.21229733E-02, 0.42605651E-02,
     +  0.41965470E-03,-0.10281521E-01, 0.10750161E-01, 0.21190000E-02,
     +  0.54505230E-02, 0.42349020E-02, 0.58504160E-03, 0.20043990E-02,
     +  0.23136580E-02,-0.20856800E-02,-0.51587942E-03,-0.10541774E-03,
     + -0.54872680E-02,-0.11760670E-02,-0.18356833E-02,-0.10275380E-02,
     + -0.32406661E-03, 0.13329260E-03, 0.17918602E-02,-0.27108610E-04,
     + -0.30272890E-03,-0.69604400E-03, 0.25639744E-03, 0.22238030E-03,
     + -0.44048270E-03,-0.93614053E-03, 0.89365790E-03,-0.11766690E-03,
     +  0.19701523E-02, 0.39204271E-03,-0.29359690E-03, 0.36105320E-03,
     +  0.37405773E-03, 0.43009550E-03, 0.70304691E-03, 0.55237830E-03,
     +  0.44366350E-03, 0.79824580E-03, 0.45050060E-03,-0.94793993E-03,
     + -0.22886721E-02,-0.22699884E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_e_q3ex    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      p_e_q3ex    =p_e_q3ex    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)            *x41*x52
      p_e_q3ex    =p_e_q3ex    
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)    *x22*x31    *x51
      p_e_q3ex    =p_e_q3ex    
     9  +coeff( 27)*x11        *x41*x51
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)        *x32*x41*x51
     3  +coeff( 30)            *x43*x51
     4  +coeff( 31)    *x21*x31    *x52
     5  +coeff( 32)    *x21    *x41*x52
     6  +coeff( 33)        *x31    *x53
     7  +coeff( 34)            *x41*x53
     8  +coeff( 35)*x12    *x31        
      p_e_q3ex    =p_e_q3ex    
     9  +coeff( 36)    *x24*x31        
     1  +coeff( 37)*x12        *x41    
     2  +coeff( 38)*x11*x22    *x41    
     3  +coeff( 39)        *x32*x43    
     4  +coeff( 40)    *x23*x31    *x51
     5  +coeff( 41)    *x21*x33    *x51
     6  +coeff( 42)*x11*x21    *x41*x51
     7  +coeff( 43)        *x31*x42*x52
     8  +coeff( 44)    *x21*x31    *x53
      p_e_q3ex    =p_e_q3ex    
     9  +coeff( 45)        *x31    *x54
     1  +coeff( 46)*x11*x23*x31        
     2  +coeff( 47)*x11*x23    *x41    
     3  +coeff( 48)    *x23*x32*x41    
     4  +coeff( 49)*x11*x21*x31*x42    
     5  +coeff( 50)*x11*x21    *x43    
c
      return
      end
      function l_e_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 29)
      data ncoeff/ 28/
      data avdat/ -0.2786803E-01/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49875E-01,-0.39017E-01,-0.49382E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49548E-01, 0.59522E-01, 0.49999E-01, 0.38045E-01, 0.49777E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27625940E-01,-0.31045360E+00, 0.15375980E+00,-0.31422361E-01,
     + -0.32568650E-01, 0.28328020E-01,-0.14242880E-01,-0.13605430E-01,
     + -0.79342052E-02,-0.65823490E-02,-0.10204131E-01, 0.18741370E-01,
     +  0.75362260E-02,-0.23680560E-02,-0.23673440E-03, 0.68621564E-04,
     +  0.21621072E-02,-0.25449534E-02,-0.57169330E-02, 0.30572053E-02,
     +  0.13676580E-01,-0.23958270E-02, 0.11698753E-02, 0.44583100E-02,
     +  0.17153130E-02, 0.13751140E-02, 0.34827331E-02, 0.14314740E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_e_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)*x11*x21            
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)    *x21        *x51
      l_e_q3ex    =l_e_q3ex    
     9  +coeff(  9)*x12                
     1  +coeff( 10)                *x52
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)*x11*x21        *x51
     6  +coeff( 15)        *x31*x41*x52
     7  +coeff( 16)    *x23*x31*x41    
     8  +coeff( 17)*x12            *x51
      l_e_q3ex    =l_e_q3ex    
     9  +coeff( 18)        *x32        
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)    *x21*x32        
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)                *x53
     6  +coeff( 24)*x11    *x31*x41    
     7  +coeff( 25)        *x31*x41*x51
     8  +coeff( 26)*x11    *x32        
      l_e_q3ex    =l_e_q3ex    
     9  +coeff( 27)*x11        *x42    
     1  +coeff( 28)*x11            *x52
c
      return
      end
      function x_e_fp      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 99)
      data ncoeff/ 98/
      data avdat/ -0.3443713E-01/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49382E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49548E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49777E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.37155032E-01,-0.10795770E-01,-0.10862133E+00, 0.59015540E+00,
     +  0.12794470E-01,-0.13701843E-02,-0.14273820E-01,-0.12867620E-01,
     +  0.31791050E-02, 0.46057004E-01,-0.51412480E-02,-0.93342860E-02,
     + -0.45905023E-01,-0.39116973E-02,-0.26770532E-01,-0.47983870E-02,
     + -0.14957030E-02, 0.10860082E-01, 0.43876940E-02,-0.17392800E-02,
     +  0.10948830E-02,-0.31128560E-02,-0.13197450E-02, 0.23570430E-02,
     +  0.13995260E-02, 0.45618880E-02, 0.11981771E-01, 0.16087020E-01,
     +  0.20314231E-02,-0.13078850E-01,-0.20036192E-02,-0.60705083E-03,
     +  0.45807320E-03,-0.71897820E-03, 0.21150584E-02,-0.26883983E-02,
     +  0.99664570E-03, 0.32924110E-02, 0.14450450E-02,-0.81640250E-02,
     +  0.13990630E-02,-0.41946210E-02,-0.11401880E-01,-0.97624573E-03,
     +  0.32350690E-02, 0.57721263E-02, 0.14408093E-02, 0.40357510E-02,
     +  0.92239370E-02,-0.62572910E-02,-0.48272430E-03, 0.31707854E-03,
     +  0.13551170E-04,-0.39282660E-03,-0.17429814E-02,-0.31722830E-02,
     + -0.35667850E-02,-0.12593740E-02,-0.30908910E-02, 0.10857554E-02,
     + -0.24426590E-02,-0.38973700E-02, 0.36684160E-02,-0.78689280E-02,
     + -0.40648221E-02,-0.11384733E-02,-0.44527670E-03,-0.42579141E-02,
     +  0.95703680E-02,-0.21074640E-02, 0.61807304E-03,-0.15717081E-02,
     +  0.25850570E-02, 0.54039840E-03,-0.62351000E-03,-0.48759960E-03,
     + -0.10140950E-01,-0.10282380E-01,-0.61396700E-02, 0.80258090E-03,
     + -0.53488942E-02,-0.29020560E-02,-0.17483963E-02,-0.29232094E-02,
     + -0.11706253E-02, 0.11320370E-02,-0.65684220E-02, 0.45133850E-02,
     +  0.99152862E-03,-0.11383234E-02, 0.27339190E-02,-0.14313000E-02,
     + -0.22990350E-02,-0.26455831E-02, 0.14904000E-02, 0.37392620E-02,
     + -0.33634493E-02,-0.18704290E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x14 = x13*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_e_fp      =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      x_e_fp      =x_e_fp      
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)*x12                
     3  +coeff( 12)*x11            *x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)*x11    *x31*x41    
      x_e_fp      =x_e_fp      
     9  +coeff( 18)*x11        *x42    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)        *x32    *x51
     3  +coeff( 21)            *x42*x51
     4  +coeff( 22)*x11*x21        *x51
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)*x12            *x51
     7  +coeff( 25)*x11            *x52
     8  +coeff( 26)                *x53
      x_e_fp      =x_e_fp      
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)*x11*x21*x31*x41    
     3  +coeff( 30)*x11*x21    *x42    
     4  +coeff( 31)    *x23*x32        
     5  +coeff( 32)    *x21*x34        
     6  +coeff( 33)        *x33*x41*x51
     7  +coeff( 34)*x12*x24            
     8  +coeff( 35)    *x23            
      x_e_fp      =x_e_fp      
     9  +coeff( 36)    *x21*x32        
     1  +coeff( 37)*x11    *x32        
     2  +coeff( 38)        *x31*x41*x51
     3  +coeff( 39)*x12*x21            
     4  +coeff( 40)    *x24            
     5  +coeff( 41)    *x21*x32    *x51
     6  +coeff( 42)    *x21    *x42*x51
     7  +coeff( 43)*x12*x22            
     8  +coeff( 44)*x12    *x31*x41    
      x_e_fp      =x_e_fp      
     9  +coeff( 45)*x12        *x42    
     1  +coeff( 46)*x11*x22        *x51
     2  +coeff( 47)*x11    *x31*x41*x51
     3  +coeff( 48)*x13*x21            
     4  +coeff( 49)    *x22    *x42*x51
     5  +coeff( 50)*x11*x21    *x42*x51
     6  +coeff( 51)    *x21*x33*x41*x51
     7  +coeff( 52)*x12*x23        *x51
     8  +coeff( 53)    *x24*x33*x41    
      x_e_fp      =x_e_fp      
     9  +coeff( 54)*x13                
     1  +coeff( 55)    *x22*x31*x41    
     2  +coeff( 56)    *x23        *x51
     3  +coeff( 57)    *x21*x31*x41*x51
     4  +coeff( 58)*x11    *x32    *x51
     5  +coeff( 59)    *x22        *x52
     6  +coeff( 60)        *x32    *x52
     7  +coeff( 61)            *x42*x52
     8  +coeff( 62)*x12*x21        *x51
      x_e_fp      =x_e_fp      
     9  +coeff( 63)*x11*x21        *x52
     1  +coeff( 64)    *x23*x31*x41    
     2  +coeff( 65)    *x23    *x42    
     3  +coeff( 66)*x12            *x52
     4  +coeff( 67)                *x54
     5  +coeff( 68)*x11*x22    *x42    
     6  +coeff( 69)*x11*x23        *x51
     7  +coeff( 70)*x11*x21*x32*x42    
     8  +coeff( 71)    *x24    *x42*x51
      x_e_fp      =x_e_fp      
     9  +coeff( 72)*x11*x25        *x51
     1  +coeff( 73)*x12*x24        *x51
     2  +coeff( 74)    *x22*x32        
     3  +coeff( 75)*x11*x21*x32        
     4  +coeff( 76)*x14                
     5  +coeff( 77)    *x21*x32*x42    
     6  +coeff( 78)    *x21*x31*x43    
     7  +coeff( 79)    *x21    *x44    
     8  +coeff( 80)*x13            *x51
      x_e_fp      =x_e_fp      
     9  +coeff( 81)    *x24        *x51
     1  +coeff( 82)*x12*x21*x31*x41    
     2  +coeff( 83)*x12*x21    *x42    
     3  +coeff( 84)    *x23        *x52
     4  +coeff( 85)    *x21*x32    *x52
     5  +coeff( 86)    *x21    *x42*x52
     6  +coeff( 87)*x12*x22        *x51
     7  +coeff( 88)*x11*x22        *x52
     8  +coeff( 89)*x11    *x32    *x52
      x_e_fp      =x_e_fp      
     9  +coeff( 90)            *x42*x53
     1  +coeff( 91)    *x22*x32*x42    
     2  +coeff( 92)*x12*x21        *x52
     3  +coeff( 93)    *x21    *x44*x51
     4  +coeff( 94)*x12*x22    *x42    
     5  +coeff( 95)*x12    *x32*x42    
     6  +coeff( 96)    *x22    *x42*x52
     7  +coeff( 97)*x11*x21    *x42*x52
     8  +coeff( 98)    *x21*x35*x41    
c
      return
      end
      function t_e_fp      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.7701758E-02/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49382E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49548E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49777E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.81675900E-02,-0.27736354E-01, 0.97573094E-01,-0.51001930E-02,
     +  0.15689012E-04,-0.12602712E-03,-0.25627131E-02,-0.17744784E-02,
     +  0.41822350E-02,-0.89948440E-02, 0.30736250E-02, 0.29798030E-03,
     + -0.49082940E-03,-0.35448911E-05,-0.39611002E-02,-0.15759330E-02,
     + -0.40961140E-03, 0.42070680E-03,-0.41418450E-03, 0.47379353E-03,
     + -0.15297280E-02,-0.98346594E-04,-0.27851123E-03, 0.21638721E-02,
     +  0.21449564E-02,-0.17100160E-03,-0.88800681E-03, 0.54135170E-03,
     +  0.18661210E-02, 0.45042580E-03,-0.24743811E-02, 0.32730930E-03,
     +  0.20511262E-03,-0.12664140E-02,-0.28046960E-04, 0.66318723E-03,
     + -0.81379380E-03, 0.63841380E-03,-0.99996151E-03, 0.30493061E-03,
     +  0.28953500E-03, 0.45558000E-03,-0.10580624E-02,-0.20358310E-03,
     + -0.33107610E-03,-0.59000484E-03,-0.18790720E-03, 0.13143791E-03,
     + -0.28493290E-03, 0.15165310E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_e_fp      =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_e_fp      =t_e_fp      
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)        *x32    *x51
      t_e_fp      =t_e_fp      
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)*x12                
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)    *x22*x31*x41    
     6  +coeff( 24)*x11        *x42    
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)*x11*x21        *x51
      t_e_fp      =t_e_fp      
     9  +coeff( 27)    *x21    *x42*x51
     1  +coeff( 28)*x11            *x52
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)*x11*x21*x31*x41    
     4  +coeff( 31)*x11*x21    *x42    
     5  +coeff( 32)*x12            *x51
     6  +coeff( 33)        *x33*x41*x51
     7  +coeff( 34)*x12*x22            
     8  +coeff( 35)*x11    *x34        
      t_e_fp      =t_e_fp      
     9  +coeff( 36)*x12        *x42    
     1  +coeff( 37)    *x21*x33*x41*x51
     2  +coeff( 38)        *x31*x41*x51
     3  +coeff( 39)    *x24            
     4  +coeff( 40)*x11    *x32        
     5  +coeff( 41)        *x31*x41*x52
     6  +coeff( 42)*x11*x22        *x51
     7  +coeff( 43)*x11*x21    *x42*x53
     8  +coeff( 44)*x11*x22            
      t_e_fp      =t_e_fp      
     9  +coeff( 45)    *x23        *x51
     1  +coeff( 46)    *x21*x31*x41*x51
     2  +coeff( 47)    *x22        *x52
     3  +coeff( 48)        *x32    *x52
     4  +coeff( 49)            *x42*x52
     5  +coeff( 50)*x12*x21            
c
      return
      end
      function y_e_fp      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49382E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49548E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49777E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.32913632E-01,-0.29685720E-01,-0.81914350E-02,-0.13443930E-01,
     + -0.41986380E-02, 0.29732933E-01, 0.46267810E-02,-0.11502980E-02,
     +  0.61857572E-04, 0.36575960E-02,-0.16593800E-01,-0.17056693E-02,
     + -0.15918974E-02,-0.31741410E-02,-0.11255603E-01, 0.50170314E-02,
     +  0.50080921E-02,-0.30877540E-03,-0.42647020E-03, 0.47030120E-03,
     +  0.15827363E-01, 0.31421160E-02, 0.45247450E-03, 0.10516723E-02,
     + -0.14805350E-02,-0.30792530E-02, 0.42882250E-02,-0.72004362E-02,
     + -0.14548230E-02,-0.61118390E-03,-0.33805693E-02, 0.10225764E-03,
     +  0.29503300E-03,-0.25052353E-02,-0.34494120E-02, 0.15580150E-03,
     +  0.99453400E-03, 0.41377870E-03,-0.85272770E-04,-0.44857603E-02,
     + -0.46826801E-02, 0.11747243E-02, 0.18950340E-04, 0.12908454E-02,
     +  0.42989840E-03, 0.22046091E-02,-0.97590270E-03,-0.21433460E-03,
     +  0.22797930E-02, 0.29127071E-02,-0.12159710E-02, 0.16683212E-02,
     +  0.67305760E-03, 0.96093220E-03, 0.12833820E-02, 0.23796040E-02,
     + -0.30213890E-02, 0.23250480E-02,-0.22437204E-02,-0.26154870E-03,
     + -0.34936063E-03, 0.30383453E-02,-0.16253163E-02,-0.76994742E-03,
     +  0.23006324E-03,-0.26845410E-02, 0.15645333E-02,-0.48591571E-03,
     + -0.14327161E-02, 0.58747570E-03,-0.22989662E-02,-0.63732500E-03,
     +  0.37344420E-02,-0.23595730E-03, 0.14310773E-03, 0.58180210E-03,
     +  0.26523560E-02, 0.10913573E-02,-0.97333040E-03, 0.26622070E-02,
     +  0.68718061E-03, 0.56164310E-03, 0.16311430E-02,-0.11635500E-02,
     + -0.15842470E-02,-0.26282270E-02, 0.79617014E-03, 0.17206370E-03,
     +  0.42587830E-03,-0.44570590E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_e_fp      =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      y_e_fp      =y_e_fp      
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)            *x41*x52
      y_e_fp      =y_e_fp      
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x32*x41    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)    *x22*x31    *x51
      y_e_fp      =y_e_fp      
     9  +coeff( 27)*x11        *x41*x51
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)        *x32*x41*x51
     3  +coeff( 30)        *x31*x42*x51
     4  +coeff( 31)            *x43*x51
     5  +coeff( 32)    *x21*x31    *x52
     6  +coeff( 33)    *x21    *x41*x52
     7  +coeff( 34)        *x31    *x53
     8  +coeff( 35)            *x41*x53
      y_e_fp      =y_e_fp      
     9  +coeff( 36)*x12    *x31        
     1  +coeff( 37)*x11*x22*x31        
     2  +coeff( 38)    *x24*x31        
     3  +coeff( 39)*x11    *x33        
     4  +coeff( 40)*x12        *x41    
     5  +coeff( 41)*x11*x22    *x41    
     6  +coeff( 42)    *x24    *x41    
     7  +coeff( 43)        *x34*x41    
     8  +coeff( 44)*x11    *x31*x42    
      y_e_fp      =y_e_fp      
     9  +coeff( 45)*x11        *x43    
     1  +coeff( 46)*x11*x21*x31    *x51
     2  +coeff( 47)    *x23*x31    *x51
     3  +coeff( 48)    *x21*x33    *x51
     4  +coeff( 49)*x11*x21    *x41*x51
     5  +coeff( 50)    *x23    *x41*x51
     6  +coeff( 51)    *x21*x32*x41*x51
     7  +coeff( 52)    *x21    *x43*x51
     8  +coeff( 53)*x11    *x31    *x52
      y_e_fp      =y_e_fp      
     9  +coeff( 54)    *x22*x31    *x52
     1  +coeff( 55)    *x22    *x41*x52
     2  +coeff( 56)        *x31*x42*x52
     3  +coeff( 57)    *x21    *x41*x53
     4  +coeff( 58)        *x31    *x54
     5  +coeff( 59)            *x41*x54
     6  +coeff( 60)*x12*x21*x31        
     7  +coeff( 61)*x11*x21*x33        
     8  +coeff( 62)*x12*x21    *x41    
      y_e_fp      =y_e_fp      
     9  +coeff( 63)*x11*x23    *x41    
     1  +coeff( 64)*x12    *x31    *x51
     2  +coeff( 65)*x12        *x41*x51
     3  +coeff( 66)*x11*x22    *x41*x51
     4  +coeff( 67)*x11    *x32*x41*x51
     5  +coeff( 68)    *x22*x32*x41*x51
     6  +coeff( 69)        *x33*x42*x51
     7  +coeff( 70)*x11*x21    *x41*x52
     8  +coeff( 71)*x11    *x31    *x53
      y_e_fp      =y_e_fp      
     9  +coeff( 72)        *x33    *x53
     1  +coeff( 73)*x11        *x41*x53
     2  +coeff( 74)*x12*x22*x31        
     3  +coeff( 75)*x11    *x34*x41    
     4  +coeff( 76)*x12*x21    *x41*x51
     5  +coeff( 77)    *x21*x33*x42*x51
     6  +coeff( 78)*x12    *x31    *x52
     7  +coeff( 79)*x12        *x41*x52
     8  +coeff( 80)    *x23*x31    *x53
      y_e_fp      =y_e_fp      
     9  +coeff( 81)*x11    *x31    *x54
     1  +coeff( 82)        *x33    *x54
     2  +coeff( 83)*x11*x23    *x43    
     3  +coeff( 84)*x12*x22*x31    *x51
     4  +coeff( 85)*x12*x22    *x41*x51
     5  +coeff( 86)*x11*x23*x31    *x52
     6  +coeff( 87)    *x24*x31    *x53
     7  +coeff( 88)    *x22*x33        
     8  +coeff( 89)*x11    *x32*x41    
      y_e_fp      =y_e_fp      
     9  +coeff( 90)    *x22*x32*x41    
c
      return
      end
      function p_e_fp      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49999E-01,-0.39017E-01,-0.49382E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49548E-01, 0.59522E-01, 0.49999E-01, 0.39017E-01, 0.49777E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15683950E-01,-0.45328780E-01, 0.47251110E-02, 0.21510580E-01,
     + -0.96115730E-02,-0.16709933E-01,-0.21226932E-02, 0.42602160E-02,
     +  0.41997690E-03,-0.10281050E-01, 0.10751363E-01, 0.21192760E-02,
     +  0.54506133E-02, 0.42353550E-02, 0.58259110E-03, 0.20028632E-02,
     +  0.23060804E-02,-0.20858694E-02,-0.51687530E-03,-0.10642402E-03,
     + -0.54871410E-02,-0.11772260E-02,-0.18411360E-02,-0.10291303E-02,
     + -0.32885810E-03, 0.13445952E-03, 0.17969170E-02,-0.27114443E-04,
     + -0.30402504E-03,-0.69690660E-03, 0.25982144E-03, 0.22648890E-03,
     + -0.44754600E-03,-0.94445050E-03, 0.89366873E-03,-0.11772910E-03,
     +  0.19699824E-02, 0.39296233E-03,-0.29693701E-03, 0.36224551E-03,
     +  0.37290391E-03, 0.42719670E-03, 0.70249290E-03, 0.55551272E-03,
     +  0.45923080E-03, 0.79812703E-03, 0.45003223E-03,-0.95158640E-03,
     + -0.22874441E-02,-0.22714461E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_e_fp      =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      p_e_fp      =p_e_fp      
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)            *x41*x52
      p_e_fp      =p_e_fp      
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)    *x22*x31    *x51
      p_e_fp      =p_e_fp      
     9  +coeff( 27)*x11        *x41*x51
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)        *x32*x41*x51
     3  +coeff( 30)            *x43*x51
     4  +coeff( 31)    *x21*x31    *x52
     5  +coeff( 32)    *x21    *x41*x52
     6  +coeff( 33)        *x31    *x53
     7  +coeff( 34)            *x41*x53
     8  +coeff( 35)*x12    *x31        
      p_e_fp      =p_e_fp      
     9  +coeff( 36)    *x24*x31        
     1  +coeff( 37)*x12        *x41    
     2  +coeff( 38)*x11*x22    *x41    
     3  +coeff( 39)        *x32*x43    
     4  +coeff( 40)    *x23*x31    *x51
     5  +coeff( 41)    *x21*x33    *x51
     6  +coeff( 42)*x11*x21    *x41*x51
     7  +coeff( 43)        *x31*x42*x52
     8  +coeff( 44)    *x21*x31    *x53
      p_e_fp      =p_e_fp      
     9  +coeff( 45)        *x31    *x54
     1  +coeff( 46)*x11*x23*x31        
     2  +coeff( 47)*x11*x23    *x41    
     3  +coeff( 48)    *x23*x32*x41    
     4  +coeff( 49)*x11*x21*x31*x42    
     5  +coeff( 50)*x11*x21    *x43    
c
      return
      end
      function l_e_fp      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.3295795E-01/
      data xmin/
     1 -0.49884E-01,-0.59666E-01,-0.49875E-01,-0.39017E-01,-0.49382E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49548E-01, 0.59522E-01, 0.49999E-01, 0.38045E-01, 0.49777E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.32771060E-01,-0.31099164E+00, 0.15420630E+00,-0.32572094E-01,
     + -0.33649533E-01,-0.20670500E-01, 0.27861652E-01,-0.13091820E-01,
     + -0.80019701E-02,-0.12782110E-01,-0.56333480E-02, 0.21613173E-01,
     + -0.30136280E-02,-0.40285503E-02, 0.13067961E-01, 0.47722733E-02,
     +  0.76484300E-02,-0.46282312E-02, 0.33041180E-02, 0.31001810E-02,
     +  0.31465990E-02, 0.47435070E-02,-0.11720681E-02, 0.14726460E-02,
     +  0.12310660E-02, 0.21317450E-02, 0.30913040E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_e_fp      =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)*x11            *x51
      l_e_fp      =l_e_fp      
     9  +coeff(  9)*x12                
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)        *x32        
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)                *x53
     8  +coeff( 17)*x11*x22            
      l_e_fp      =l_e_fp      
     9  +coeff( 18)*x11*x21        *x51
     1  +coeff( 19)*x12            *x51
     2  +coeff( 20)    *x21*x32        
     3  +coeff( 21)        *x31*x41*x51
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)        *x32    *x51
     7  +coeff( 25)*x11    *x32        
     8  +coeff( 26)*x11        *x42    
      l_e_fp      =l_e_fp      
     9  +coeff( 27)    *x21    *x42*x51
c
      return
      end
      function x_h_q1ex      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.3375489E-03/
      data xmin/
     1 -0.49925E-01,-0.59938E-01,-0.49982E-01,-0.39996E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49930E-01, 0.59958E-01, 0.49982E-01, 0.39996E-01, 0.49981E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35849900E-03, 0.12244262E+00, 0.15007774E-05, 0.22364171E-01,
     +  0.31374450E-02, 0.12903060E-02,-0.65926543E-03,-0.15062403E-03,
     + -0.21326651E-03,-0.40822600E-04,-0.10107400E-03, 0.38427213E-04,
     + -0.17633401E-03,-0.19893560E-03,-0.62794460E-03,-0.25165473E-03,
     + -0.62113700E-04,-0.81925441E-04,-0.43756240E-03,-0.23815364E-06,
     + -0.50492550E-03,-0.19327270E-04, 0.17702460E-05,-0.21925533E-03,
     + -0.54374330E-03,-0.12292171E-03,-0.41188520E-04,-0.40384511E-04,
     + -0.41558740E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_h_q1ex      =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)*x11            *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x21        *x52
      x_h_q1ex      =x_h_q1ex      
     9  +coeff(  9)*x11*x22            
     1  +coeff( 10)    *x21*x33*x41    
     2  +coeff( 11)    *x25*x32        
     3  +coeff( 12)    *x21*x35*x41    
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)*x11        *x42    
     8  +coeff( 17)*x11            *x52
      x_h_q1ex      =x_h_q1ex      
     9  +coeff( 18)*x12*x21            
     1  +coeff( 19)*x11*x22*x31*x41    
     2  +coeff( 20)*x11    *x33*x41    
     3  +coeff( 21)*x11*x22    *x42    
     4  +coeff( 22)*x13    *x32        
     5  +coeff( 23)*x13    *x31*x41    
     6  +coeff( 24)*x11    *x31*x41    
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)*x11*x22*x32        
      x_h_q1ex      =x_h_q1ex      
     9  +coeff( 27)    *x25*x31*x41    
     1  +coeff( 28)*x11    *x32        
     2  +coeff( 29)    *x23*x31*x41    
c
      return
      end
      function t_h_q1ex      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1013849E-03/
      data xmin/
     1 -0.49925E-01,-0.59938E-01,-0.49982E-01,-0.39996E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49930E-01, 0.59958E-01, 0.49982E-01, 0.39996E-01, 0.49981E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10004510E-03,-0.27876730E-02, 0.14566250E-05,-0.24970682E-01,
     + -0.16891700E-05,-0.70258101E-06,-0.13733160E-05,-0.14873832E-05,
     +  0.29021590E-02,-0.29398990E-07,-0.14902230E-05,-0.18963133E-03,
     + -0.34513391E-03, 0.11444013E-02,-0.13421362E-03,-0.14102200E-03,
     + -0.14221430E-03,-0.52341394E-04,-0.92076000E-04,-0.19629693E-03,
     + -0.74997770E-03, 0.39205704E-04,-0.10215490E-02,-0.11979093E-02,
     + -0.12678012E-03,-0.11548054E-02,-0.43242150E-03,-0.62793132E-03,
     + -0.33738990E-03,-0.86620962E-03,-0.25112820E-04,-0.12803700E-03,
     + -0.30933500E-03,-0.45673190E-04,-0.11021982E-03,-0.34900770E-03,
     +  0.76455120E-05,-0.95542210E-03,-0.56038383E-03,-0.24665250E-03,
     +  0.12653240E-03, 0.12842193E-03,-0.26372753E-03,-0.18481800E-02,
     + -0.95761573E-03,-0.16171470E-02,-0.16091020E-02,-0.60661500E-03,
     + -0.10225111E-02,-0.24892800E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_h_q1ex      =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_h_q1ex      =t_h_q1ex      
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)*x11        *x42    
      t_h_q1ex      =t_h_q1ex      
     9  +coeff( 18)*x11            *x52
     1  +coeff( 19)*x12*x21            
     2  +coeff( 20)    *x23*x32        
     3  +coeff( 21)    *x23*x31*x41    
     4  +coeff( 22)    *x21*x33*x41    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)*x11*x22*x31*x41    
     7  +coeff( 25)*x11    *x33*x41    
     8  +coeff( 26)*x11*x22    *x42    
      t_h_q1ex      =t_h_q1ex      
     9  +coeff( 27)*x11    *x32*x42    
     1  +coeff( 28)*x11    *x31*x43    
     2  +coeff( 29)*x11        *x44    
     3  +coeff( 30)    *x23*x33*x41    
     4  +coeff( 31)*x11*x24*x32        
     5  +coeff( 32)    *x21*x32        
     6  +coeff( 33)    *x21*x31*x41    
     7  +coeff( 34)*x11    *x32        
     8  +coeff( 35)*x11    *x31*x41    
      t_h_q1ex      =t_h_q1ex      
     9  +coeff( 36)*x12*x21    *x42    
     1  +coeff( 37)    *x21    *x42*x51
     2  +coeff( 38)    *x21*x31*x43    
     3  +coeff( 39)    *x21    *x44    
     4  +coeff( 40)*x11*x22*x32        
     5  +coeff( 41)    *x23*x31*x41*x51
     6  +coeff( 42)    *x23    *x42*x51
     7  +coeff( 43)*x12*x21*x31*x41    
     8  +coeff( 44)    *x23*x32*x42    
      t_h_q1ex      =t_h_q1ex      
     9  +coeff( 45)    *x21*x34*x42    
     1  +coeff( 46)    *x23*x31*x43    
     2  +coeff( 47)    *x21*x33*x43    
     3  +coeff( 48)    *x23    *x44    
     4  +coeff( 49)    *x21*x32*x44    
     5  +coeff( 50)*x11*x24    *x42    
c
      return
      end
      function y_h_q1ex      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49925E-01,-0.59938E-01,-0.49982E-01,-0.39996E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49930E-01, 0.59958E-01, 0.49982E-01, 0.39996E-01, 0.49981E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.81667990E-01, 0.17635662E+00,-0.16872921E-02,-0.26787910E-02,
     + -0.71292473E-04,-0.39347110E-04, 0.37725013E-03,-0.16700114E-03,
     + -0.13148132E-03, 0.28443850E-03, 0.13932040E-03,-0.88420483E-05,
     +  0.23325583E-04, 0.30903561E-06, 0.89804052E-04,-0.55435881E-03,
     + -0.79695892E-03,-0.40256010E-03,-0.57762870E-03,-0.41430720E-05,
     +  0.58805642E-04, 0.28063252E-04,-0.26085250E-03,-0.37284163E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_h_q1ex      =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31    *x51
     4  +coeff(  4)            *x41*x51
     5  +coeff(  5)    *x22*x31        
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)            *x43    
     8  +coeff(  8)*x11*x21*x31        
      y_h_q1ex      =y_h_q1ex      
     9  +coeff(  9)*x11*x21    *x41    
     1  +coeff( 10)        *x31*x42    
     2  +coeff( 11)            *x41*x52
     3  +coeff( 12)*x12    *x31        
     4  +coeff( 13)*x12        *x41    
     5  +coeff( 14)    *x22*x31    *x52
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)*x11*x23*x31        
     8  +coeff( 17)*x11*x23    *x41    
      y_h_q1ex      =y_h_q1ex      
     9  +coeff( 18)*x12*x22*x31        
     1  +coeff( 19)*x12*x22    *x41    
     2  +coeff( 20)        *x33        
     3  +coeff( 21)        *x32*x41    
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)    *x24*x31        
     6  +coeff( 24)    *x24    *x41    
c
      return
      end
      function p_h_q1ex      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 49)
      data ncoeff/ 48/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49925E-01,-0.59938E-01,-0.49982E-01,-0.39996E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49930E-01, 0.59958E-01, 0.49982E-01, 0.39996E-01, 0.49981E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29884570E-01, 0.89014761E-01,-0.16226782E-02,-0.26435700E-02,
     + -0.19364574E-03,-0.78411170E-04,-0.12668280E-03,-0.19235490E-03,
     + -0.10036650E-03, 0.13932560E-03,-0.17473510E-04,-0.52379960E-05,
     + -0.88935110E-03,-0.40997040E-03,-0.61950442E-03, 0.16740550E-04,
     + -0.11327050E-03, 0.88133510E-04,-0.33659310E-03,-0.58752542E-03,
     + -0.79981720E-04, 0.57252582E-04, 0.30412284E-04, 0.39193230E-04,
     +  0.61300401E-04,-0.23371774E-04,-0.37409530E-04,-0.15823790E-04,
     +  0.18078353E-04,-0.45513610E-03,-0.69376590E-03,-0.39820700E-03,
     +  0.78391400E-04, 0.84996320E-04,-0.66514380E-06, 0.43463690E-05,
     + -0.28101670E-03,-0.40066360E-03,-0.18537843E-03, 0.10897012E-04,
     + -0.11210360E-03,-0.11546600E-04,-0.23433492E-04,-0.89770040E-04,
     + -0.24001250E-04, 0.74796820E-04, 0.18843044E-04,-0.84056744E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_h_q1ex      =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31    *x51
     4  +coeff(  4)            *x41*x51
     5  +coeff(  5)    *x22    *x41    
     6  +coeff(  6)*x11*x21*x31        
     7  +coeff(  7)*x11*x21    *x41    
     8  +coeff(  8)    *x24*x31        
      p_h_q1ex      =p_h_q1ex      
     9  +coeff(  9)    *x22*x31        
     1  +coeff( 10)            *x41*x52
     2  +coeff( 11)*x12        *x41    
     3  +coeff( 12)    *x22*x31    *x52
     4  +coeff( 13)*x11*x23    *x41    
     5  +coeff( 14)*x12*x22*x31        
     6  +coeff( 15)*x12*x22    *x41    
     7  +coeff( 16)*x11*x23*x33        
     8  +coeff( 17)        *x31*x42    
      p_h_q1ex      =p_h_q1ex      
     9  +coeff( 18)        *x31    *x52
     1  +coeff( 19)    *x24    *x41    
     2  +coeff( 20)*x11*x23*x31        
     3  +coeff( 21)        *x32*x41    
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)*x11*x21    *x41*x51
     6  +coeff( 24)    *x24*x31    *x51
     7  +coeff( 25)*x11*x23*x31    *x51
     8  +coeff( 26)        *x33        
      p_h_q1ex      =p_h_q1ex      
     9  +coeff( 27)            *x43    
     1  +coeff( 28)*x12    *x31        
     2  +coeff( 29)*x11*x21*x31    *x51
     3  +coeff( 30)*x11*x21*x32*x41    
     4  +coeff( 31)*x11*x21*x31*x42    
     5  +coeff( 32)*x11*x21    *x43    
     6  +coeff( 33)*x12*x22*x31    *x51
     7  +coeff( 34)*x12*x22    *x41*x51
     8  +coeff( 35)    *x21    *x41    
      p_h_q1ex      =p_h_q1ex      
     9  +coeff( 36)*x11*x22*x31        
     1  +coeff( 37)    *x22*x32*x41    
     2  +coeff( 38)    *x22*x31*x42    
     3  +coeff( 39)    *x22    *x43    
     4  +coeff( 40)    *x21*x33    *x51
     5  +coeff( 41)*x11*x21*x33        
     6  +coeff( 42)    *x23*x33        
     7  +coeff( 43)    *x23*x32*x41    
     8  +coeff( 44)    *x24*x33        
      p_h_q1ex      =p_h_q1ex      
     9  +coeff( 45)*x12    *x31*x42    
     1  +coeff( 46)*x11*x23    *x41*x51
     2  +coeff( 47)    *x21*x34*x41*x51
     3  +coeff( 48)*x11*x23*x31*x42    
c
      return
      end
      function l_h_q1ex      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/ -0.2214059E-02/
      data xmin/
     1 -0.49925E-01,-0.59938E-01,-0.49982E-01,-0.39950E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49930E-01, 0.59958E-01, 0.49967E-01, 0.39996E-01, 0.49981E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22165873E-02,-0.34503700E-02,-0.25650572E-02,-0.55814771E-02,
     + -0.51993100E-06, 0.20245570E-03,-0.40181690E-03,-0.29725960E-03,
     +  0.69063720E-07,-0.21261230E-04, 0.21917883E-03, 0.57522140E-04,
     + -0.68454570E-06, 0.16941414E-05, 0.78157250E-06, 0.42951710E-04,
     +  0.20618923E-03, 0.27167070E-04,-0.52664780E-05,-0.36626530E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_h_q1ex      =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x22            
     3  +coeff(  3)        *x31*x41    
     4  +coeff(  4)            *x42    
     5  +coeff(  5)*x12    *x32        
     6  +coeff(  6)*x11*x21            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)*x12                
      l_h_q1ex      =l_h_q1ex      
     9  +coeff(  9)                *x51
     1  +coeff( 10)    *x22        *x51
     2  +coeff( 11)            *x42*x51
     3  +coeff( 12)*x11*x21        *x51
     4  +coeff( 13)    *x22*x31*x41*x51
     5  +coeff( 14)        *x33*x41*x51
     6  +coeff( 15)        *x34*x42*x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      l_h_q1ex      =l_h_q1ex      
     9  +coeff( 18)*x12            *x51
     1  +coeff( 19)            *x41    
     2  +coeff( 20)                *x52
c
      return
      end
      function x_h_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 99)
      data ncoeff/ 98/
      data avdat/ -0.1535021E-02/
      data xmin/
     1 -0.49916E-01,-0.59181E-01,-0.49982E-01,-0.38578E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49982E-01, 0.38578E-01, 0.49981E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.71629610E-03,-0.12025470E+00,-0.65138570E-04, 0.52803613E-01,
     +  0.35387484E-02,-0.72120530E-04,-0.22658590E-03,-0.69988890E-02,
     +  0.81707330E-05,-0.35969461E-02,-0.52806534E-02, 0.55605070E-03,
     +  0.14491813E-02, 0.17013032E-02, 0.42726870E-03, 0.91785733E-03,
     +  0.58462780E-03, 0.85782600E-03, 0.50587934E-03, 0.48707930E-03,
     +  0.36417280E-03, 0.76238793E-04,-0.17766140E-03, 0.44009080E-03,
     +  0.17207160E-03, 0.34883054E-04, 0.34301590E-03, 0.25561271E-03,
     +  0.37209670E-03, 0.14418410E-03, 0.13123790E-02, 0.30981460E-02,
     +  0.22610030E-02, 0.21547640E-02, 0.85680664E-03, 0.11347250E-03,
     + -0.60969201E-03,-0.59370684E-03,-0.98023660E-04, 0.91469344E-04,
     +  0.14130212E-02, 0.73982780E-02, 0.44211000E-02, 0.18245580E-04,
     + -0.35277451E-04, 0.75117260E-03, 0.17516000E-02, 0.53244310E-02,
     +  0.82594493E-03, 0.72684930E-02, 0.17887650E-02, 0.37047574E-02,
     +  0.24282434E-02, 0.16103800E-02, 0.15692810E-02, 0.20947344E-02,
     +  0.29782264E-03, 0.67238940E-04,-0.58694300E-04,-0.21776861E-04,
     + -0.15544540E-03,-0.46403420E-04,-0.11180051E-03,-0.12809144E-03,
     +  0.52029230E-04, 0.50373280E-02, 0.18721644E-02, 0.32340750E-03,
     +  0.84178271E-03, 0.24512280E-03, 0.25591190E-02,-0.62623294E-04,
     +  0.33630210E-02, 0.63579453E-03, 0.58125610E-03, 0.86432700E-03,
     +  0.14870920E-03, 0.25814950E-02,-0.97775682E-05,-0.31750610E-04,
     + -0.23635354E-03,-0.24027430E-03,-0.13250412E-03,-0.94744270E-04,
     +  0.81627530E-04, 0.83806931E-04, 0.15077203E-02, 0.61471900E-04,
     +  0.29993991E-04,-0.41774881E-04, 0.41720623E-03, 0.10783160E-02,
     +  0.10103050E-02,-0.58429272E-04,-0.15328750E-03, 0.15702560E-02,
     +  0.42300960E-03, 0.69075444E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_h_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_h_dent    =x_h_dent    
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)*x11            *x51
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)*x12                
     8  +coeff( 17)    *x21        *x52
      x_h_dent    =x_h_dent    
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)*x11    *x31*x41    
     2  +coeff( 20)*x11        *x42    
     3  +coeff( 21)*x11            *x52
     4  +coeff( 22)*x12*x21            
     5  +coeff( 23)*x12            *x51
     6  +coeff( 24)    *x23*x32        
     7  +coeff( 25)    *x21*x34        
     8  +coeff( 26)*x11    *x34        
      x_h_dent    =x_h_dent    
     9  +coeff( 27)    *x25*x32        
     1  +coeff( 28)    *x23*x34        
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)*x11    *x32        
     4  +coeff( 31)    *x23    *x42    
     5  +coeff( 32)*x11*x22*x31*x41    
     6  +coeff( 33)*x11*x22    *x42    
     7  +coeff( 34)    *x25*x31*x41    
     8  +coeff( 35)*x11*x24*x32        
      x_h_dent    =x_h_dent    
     9  +coeff( 36)*x11*x21        *x51
     1  +coeff( 37)    *x21*x31*x41*x51
     2  +coeff( 38)    *x21    *x42*x51
     3  +coeff( 39)*x11*x22        *x51
     4  +coeff( 40)*x13                
     5  +coeff( 41)    *x23*x31*x41    
     6  +coeff( 42)    *x21*x31*x43    
     7  +coeff( 43)    *x21    *x44    
     8  +coeff( 44)    *x25        *x51
      x_h_dent    =x_h_dent    
     9  +coeff( 45)    *x23*x32    *x51
     1  +coeff( 46)*x12*x21    *x42    
     2  +coeff( 47)    *x23*x33*x41    
     3  +coeff( 48)    *x23*x32*x42    
     4  +coeff( 49)    *x21*x34*x42    
     5  +coeff( 50)    *x23*x31*x43    
     6  +coeff( 51)    *x21*x33*x43    
     7  +coeff( 52)    *x23    *x44    
     8  +coeff( 53)    *x21*x32*x44    
      x_h_dent    =x_h_dent    
     9  +coeff( 54)    *x21*x31*x45    
     1  +coeff( 55)*x12*x23*x31*x41    
     2  +coeff( 56)*x12*x23    *x42    
     3  +coeff( 57)*x13*x22*x32        
     4  +coeff( 58)    *x25    *x44    
     5  +coeff( 59)            *x42*x51
     6  +coeff( 60)    *x24            
     7  +coeff( 61)    *x21*x32    *x51
     8  +coeff( 62)    *x21        *x53
      x_h_dent    =x_h_dent    
     9  +coeff( 63)*x11    *x31*x41*x51
     1  +coeff( 64)*x11        *x42*x51
     2  +coeff( 65)    *x25            
     3  +coeff( 66)    *x21*x32*x42    
     4  +coeff( 67)*x11    *x31*x43    
     5  +coeff( 68)*x12*x21*x32        
     6  +coeff( 69)*x12*x21*x31*x41    
     7  +coeff( 70)    *x21*x35*x41    
     8  +coeff( 71)    *x25    *x42    
      x_h_dent    =x_h_dent    
     9  +coeff( 72)*x11    *x32    *x53
     1  +coeff( 73)*x11*x24    *x42    
     2  +coeff( 74)*x11    *x34*x42    
     3  +coeff( 75)*x11    *x33*x43    
     4  +coeff( 76)*x11*x22    *x44    
     5  +coeff( 77)*x11    *x32*x44    
     6  +coeff( 78)*x11*x24*x31*x43    
     7  +coeff( 79)        *x32        
     8  +coeff( 80)        *x31*x41*x51
      x_h_dent    =x_h_dent    
     9  +coeff( 81)    *x22*x31*x41    
     1  +coeff( 82)    *x22    *x42    
     2  +coeff( 83)    *x23        *x51
     3  +coeff( 84)*x11*x23            
     4  +coeff( 85)*x11*x21*x31*x41    
     5  +coeff( 86)*x11*x21    *x42    
     6  +coeff( 87)    *x21*x33*x41    
     7  +coeff( 88)*x12*x22            
     8  +coeff( 89)*x12        *x42    
      x_h_dent    =x_h_dent    
     9  +coeff( 90)*x12*x21        *x51
     1  +coeff( 91)*x11*x22*x32        
     2  +coeff( 92)*x11    *x32*x42    
     3  +coeff( 93)*x11        *x44    
     4  +coeff( 94)    *x24*x32        
     5  +coeff( 95)*x11*x22*x31*x41*x51
     6  +coeff( 96)*x11*x24*x31*x41    
     7  +coeff( 97)*x11    *x35*x41    
     8  +coeff( 98)*x11*x22*x32*x42    
c
c
c
c  transform to TRANSPORT x
c
      x_h_dent=-x_h_dent*3.346065 !multiply by -cos(30)/sin(15degrees)

      return
      end
      function t_h_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.2898941E-02/
      data xmin/
     1 -0.49916E-01,-0.59181E-01,-0.49982E-01,-0.38578E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49982E-01, 0.38578E-01, 0.49981E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24971072E-02, 0.55540412E-01, 0.57472864E-03,-0.31754672E-01,
     +  0.42544110E-02, 0.13580210E-03, 0.14423782E-02, 0.34904391E-02,
     +  0.31496793E-02,-0.33593740E-04,-0.35757510E-02,-0.57309800E-03,
     + -0.72933750E-03,-0.31520740E-02,-0.35621451E-02, 0.29423970E-02,
     +  0.21364680E-03, 0.63483130E-03,-0.32891391E-03, 0.76313290E-03,
     + -0.18012740E-03,-0.34251040E-03,-0.11230130E-03,-0.13157070E-02,
     +  0.34142981E-03,-0.21231750E-03,-0.18781810E-03, 0.10215720E-02,
     + -0.18785230E-03, 0.13397860E-03, 0.43849560E-04,-0.67037623E-03,
     +  0.14233410E-05,-0.24699040E-02,-0.19067483E-03,-0.66799414E-03,
     +  0.60462100E-04, 0.33203430E-03,-0.67167780E-03,-0.10301604E-03,
     + -0.84923201E-03, 0.33529612E-03, 0.30509500E-03, 0.47144250E-03,
     + -0.57186260E-03, 0.10521050E-03,-0.53864152E-03,-0.20163240E-03,
     + -0.18214020E-02,-0.31696290E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_h_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_h_dent    =t_h_dent    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x22        *x51
      t_h_dent    =t_h_dent    
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)*x12                
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)*x11        *x42    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)*x11*x21        *x51
     8  +coeff( 26)*x11            *x52
      t_h_dent    =t_h_dent    
     9  +coeff( 27)*x12*x21            
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)*x12            *x51
     3  +coeff( 30)    *x22*x31*x41*x51
     4  +coeff( 31)        *x33*x41*x51
     5  +coeff( 32)*x12*x22            
     6  +coeff( 33)*x11    *x34        
     7  +coeff( 34)*x11*x22    *x42    
     8  +coeff( 35)*x11*x23    *x42    
      t_h_dent    =t_h_dent    
     9  +coeff( 36)*x11*x24*x31*x41    
     1  +coeff( 37)        *x32    *x51
     2  +coeff( 38)        *x31*x41*x51
     3  +coeff( 39)    *x24            
     4  +coeff( 40)*x11    *x32        
     5  +coeff( 41)    *x22*x31*x41    
     6  +coeff( 42)    *x21*x31*x41*x51
     7  +coeff( 43)    *x21    *x42*x51
     8  +coeff( 44)*x11*x21    *x42    
      t_h_dent    =t_h_dent    
     9  +coeff( 45)    *x23    *x42    
     1  +coeff( 46)*x11*x22        *x51
     2  +coeff( 47)*x11*x22*x32        
     3  +coeff( 48)    *x24*x32        
     4  +coeff( 49)*x11*x22*x31*x41    
     5  +coeff( 50)*x12*x23            
c
      return
      end
      function y_h_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49916E-01,-0.59181E-01,-0.49982E-01,-0.38578E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49982E-01, 0.38578E-01, 0.49981E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.36056954E-01, 0.16785410E+00, 0.36235090E-02, 0.45011010E-02,
     +  0.18828533E-01,-0.16534873E-02,-0.18899680E-02,-0.31576170E-03,
     + -0.20804640E-02,-0.33639660E-02,-0.41465461E-02,-0.31985461E-02,
     + -0.27105960E-03,-0.13457190E-02, 0.15237920E-03, 0.12797672E-04,
     +  0.92349730E-04, 0.37827133E-03,-0.29552410E-03, 0.64744264E-04,
     +  0.73304000E-03,-0.23977670E-04, 0.31390282E-03, 0.63810760E-03,
     +  0.40176731E-03, 0.23505250E-04, 0.44442290E-04,-0.36995910E-03,
     + -0.65421982E-03,-0.21505240E-02, 0.63809260E-04, 0.36733243E-04,
     +  0.53751720E-04, 0.14283140E-03, 0.23496973E-04, 0.99980560E-04,
     +  0.82317560E-04, 0.79150340E-05,-0.81800730E-03, 0.12686604E-03,
     +  0.64674520E-03, 0.61792280E-03, 0.21093220E-03, 0.70311652E-04,
     +  0.24263510E-03, 0.43051200E-04, 0.29930853E-03, 0.48388243E-03,
     +  0.21783380E-05,-0.41889121E-04, 0.14387460E-03, 0.16577790E-03,
     + -0.86729550E-03,-0.16902401E-03,-0.39435290E-04, 0.57332250E-05,
     +  0.11544480E-04,-0.16459260E-03,-0.12111694E-02, 0.31366530E-04,
     +  0.19911360E-04,-0.73710434E-04,-0.91396940E-03,-0.97569794E-03,
     + -0.77486690E-04,-0.28425020E-02, 0.28339820E-03, 0.16893320E-02,
     +  0.64087202E-02,-0.20052744E-02,-0.89633200E-03,-0.61981380E-03,
     + -0.76346990E-03,-0.12133560E-02,-0.38692220E-02, 0.18126860E-03,
     +  0.52967382E-03,-0.25016530E-03,-0.55087320E-02,-0.42472230E-02,
     + -0.14603622E-03, 0.24066910E-03, 0.34402690E-03, 0.18256500E-02,
     + -0.63241994E-03,-0.75250460E-03,-0.66119461E-03, 0.35816081E-02,
     +  0.23006321E-02, 0.14070290E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_h_dent    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)*x11    *x31        
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)        *x33        
      y_h_dent    =y_h_dent    
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x21*x31    *x51
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)    *x21*x33        
     8  +coeff( 17)*x11*x21    *x41    
      y_h_dent    =y_h_dent    
     9  +coeff( 18)*x11    *x31    *x51
     1  +coeff( 19)    *x22*x31    *x51
     2  +coeff( 20)        *x33    *x51
     3  +coeff( 21)*x11        *x41*x51
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)        *x32*x41*x51
     6  +coeff( 24)        *x31*x42*x51
     7  +coeff( 25)            *x43*x51
     8  +coeff( 26)        *x31    *x53
      y_h_dent    =y_h_dent    
     9  +coeff( 27)            *x41*x53
     1  +coeff( 28)*x12    *x31        
     2  +coeff( 29)*x12        *x41    
     3  +coeff( 30)    *x22*x32*x41    
     4  +coeff( 31)        *x34*x41    
     5  +coeff( 32)        *x32*x43    
     6  +coeff( 33)    *x23*x31    *x51
     7  +coeff( 34)    *x23    *x41*x51
     8  +coeff( 35)    *x21*x32*x41*x51
      y_h_dent    =y_h_dent    
     9  +coeff( 36)    *x21*x31*x42*x51
     1  +coeff( 37)    *x21    *x43*x51
     2  +coeff( 38)    *x22*x31    *x52
     3  +coeff( 39)*x11*x23*x31        
     4  +coeff( 40)*x12    *x31    *x51
     5  +coeff( 41)    *x24*x31    *x51
     6  +coeff( 42)    *x22*x33    *x51
     7  +coeff( 43)*x12        *x41*x51
     8  +coeff( 44)        *x34*x41*x51
      y_h_dent    =y_h_dent    
     9  +coeff( 45)    *x22*x31*x42*x51
     1  +coeff( 46)        *x33*x42*x51
     2  +coeff( 47)    *x22    *x43*x51
     3  +coeff( 48)    *x22*x31    *x53
     4  +coeff( 49)        *x33    *x53
     5  +coeff( 50)        *x32*x41*x53
     6  +coeff( 51)        *x31*x42*x53
     7  +coeff( 52)            *x43*x53
     8  +coeff( 53)*x12*x22    *x41    
      y_h_dent    =y_h_dent    
     9  +coeff( 54)    *x22*x34*x41    
     1  +coeff( 55)        *x34*x43    
     2  +coeff( 56)*x12    *x31    *x52
     3  +coeff( 57)        *x33    *x54
     4  +coeff( 58)*x12    *x33    *x51
     5  +coeff( 59)    *x24*x33    *x51
     6  +coeff( 60)    *x22*x34*x41*x51
     7  +coeff( 61)    *x22*x33*x42*x51
     8  +coeff( 62)*x12    *x31    *x53
      y_h_dent    =y_h_dent    
     9  +coeff( 63)    *x24*x31    *x53
     1  +coeff( 64)    *x22*x33    *x53
     2  +coeff( 65)        *x33*x42*x53
     3  +coeff( 66)*x12*x24*x31        
     4  +coeff( 67)*x12    *x33    *x53
     5  +coeff( 68)    *x24*x33    *x53
     6  +coeff( 69)        *x31    *x51
     7  +coeff( 70)        *x32*x41    
     8  +coeff( 71)    *x21    *x41*x51
      y_h_dent    =y_h_dent    
     9  +coeff( 72)        *x31    *x52
     1  +coeff( 73)*x11*x22    *x41    
     2  +coeff( 74)    *x24    *x41    
     3  +coeff( 75)*x12*x24    *x41    
     4  +coeff( 76)*x11*x21*x31        
     5  +coeff( 77)    *x23    *x41    
     6  +coeff( 78)*x11*x22*x31        
     7  +coeff( 79)    *x22*x31*x42    
     8  +coeff( 80)    *x22    *x43    
      y_h_dent    =y_h_dent    
     9  +coeff( 81)*x11*x21    *x41*x51
     1  +coeff( 82)*x12*x21    *x41    
     2  +coeff( 83)*x11*x21*x32*x41    
     3  +coeff( 84)*x11*x21    *x43    
     4  +coeff( 85)    *x24*x33        
     5  +coeff( 86)*x12    *x31*x42    
     6  +coeff( 87)*x12        *x43    
     7  +coeff( 88)*x11*x23*x31*x42    
     8  +coeff( 89)*x11*x23    *x43    
      y_h_dent    =y_h_dent    
     9  +coeff( 90)*x11*x21*x31*x44    
c
      return
      end
      function p_h_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49916E-01,-0.59181E-01,-0.49982E-01,-0.38578E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49982E-01, 0.38578E-01, 0.49981E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.57035330E-01,-0.66649760E-01,-0.15605124E-01,-0.33370193E-01,
     +  0.79647090E-02, 0.18842274E-01, 0.86366360E-02,-0.72480850E-02,
     + -0.16314914E-03, 0.16924310E-01,-0.19928230E-01,-0.17442510E-02,
     + -0.48058440E-02,-0.41857920E-02, 0.22856484E-04, 0.59966840E-03,
     + -0.63869600E-03,-0.13773163E-02, 0.66865090E-02,-0.19101364E-03,
     +  0.28429910E-03, 0.17307520E-01, 0.14068870E-02, 0.14886234E-02,
     +  0.27483403E-02, 0.18652243E-02,-0.14780182E-02,-0.68956951E-03,
     + -0.16195170E-03,-0.33577650E-02,-0.15930490E-02,-0.18122250E-02,
     +  0.12215053E-02,-0.44475040E-02, 0.15350863E-02,-0.19898840E-02,
     + -0.14988950E-02, 0.33134693E-03,-0.66879970E-03,-0.26880903E-02,
     +  0.14459890E-02, 0.10785383E-02,-0.80949040E-03,-0.60654660E-03,
     + -0.55008672E-03, 0.12630110E-02, 0.36996061E-03, 0.53073873E-03,
     + -0.15329880E-03,-0.14889900E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_h_dent    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      p_h_dent    =p_h_dent    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)        *x31    *x52
      p_h_dent    =p_h_dent    
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)*x11*x21*x31        
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)    *x21*x33        
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x21*x32*x41    
     7  +coeff( 25)    *x21*x31*x42    
     8  +coeff( 26)    *x21    *x43    
      p_h_dent    =p_h_dent    
     9  +coeff( 27)*x11    *x31    *x51
     1  +coeff( 28)    *x22*x31    *x51
     2  +coeff( 29)        *x33    *x51
     3  +coeff( 30)*x11        *x41*x51
     4  +coeff( 31)    *x22    *x41*x51
     5  +coeff( 32)*x12    *x31        
     6  +coeff( 33)*x11*x22*x31        
     7  +coeff( 34)*x12        *x41    
     8  +coeff( 35)    *x24    *x41    
      p_h_dent    =p_h_dent    
     9  +coeff( 36)    *x22    *x43    
     1  +coeff( 37)*x11*x21    *x41*x51
     2  +coeff( 38)*x11        *x41*x52
     3  +coeff( 39)*x12*x21*x31        
     4  +coeff( 40)*x11*x23    *x41    
     5  +coeff( 41)*x11*x21    *x43    
     6  +coeff( 42)*x12        *x41*x51
     7  +coeff( 43)*x12*x22*x31        
     8  +coeff( 44)*x11*x23*x31    *x51
      p_h_dent    =p_h_dent    
     9  +coeff( 45)*x11*x21*x33    *x51
     1  +coeff( 46)*x11*x21*x33*x42    
     2  +coeff( 47)*x12*x22*x31    *x51
     3  +coeff( 48)*x12    *x33    *x51
     4  +coeff( 49)    *x22*x32*x41    
     5  +coeff( 50)    *x22*x31*x42    
c
      return
      end
      function l_h_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.2246371E-02/
      data xmin/
     1 -0.49916E-01,-0.59181E-01,-0.49982E-01,-0.37177E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38578E-01, 0.49981E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.37948470E-02, 0.23173892E+00,-0.10206040E+00,-0.17124252E-01,
     +  0.14008970E-01, 0.15702890E-01, 0.10251931E-01,-0.50747920E-02,
     + -0.96750260E-02, 0.24996130E-04,-0.14276160E-02,-0.66068270E-02,
     + -0.68480242E-02,-0.67017910E-03,-0.35562871E-02,-0.74062390E-02,
     + -0.16536230E-03,-0.88131222E-02, 0.10884250E-02,-0.20167560E-02,
     + -0.15303213E-02,-0.40561590E-02,-0.33076631E-03,-0.64736610E-02,
     + -0.20652404E-02,-0.24095120E-02,-0.79453680E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_h_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)*x11*x21            
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)*x12                
      l_h_dent    =l_h_dent    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x33*x41    
     2  +coeff( 11)        *x32        
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x23*x31*x41    
     8  +coeff( 17)    *x21*x33*x41    
      l_h_dent    =l_h_dent    
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)*x12            *x51
     2  +coeff( 20)    *x24        *x51
     3  +coeff( 21)    *x21*x32    *x52
     4  +coeff( 22)    *x23*x34        
     5  +coeff( 23)            *x41    
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)*x11    *x31*x41    
     8  +coeff( 26)*x11        *x42    
      l_h_dent    =l_h_dent    
     9  +coeff( 27)*x11            *x52
c
      return
      end
      function x_h_dext    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 98)
      data ncoeff/ 97/
      data avdat/  0.9306630E-02/
      data xmin/
     1 -0.49916E-01,-0.58959E-01,-0.49982E-01,-0.38531E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49982E-01, 0.38531E-01, 0.49826E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.65509951E-02, 0.33843564E+00, 0.13330030E+00,-0.19341273E+00,
     +  0.11760033E-01,-0.68127710E-03,-0.71706240E-02, 0.44610530E-01,
     + -0.33686220E-02,-0.91793043E-02, 0.49357502E-02,-0.13458090E-02,
     + -0.60881082E-02,-0.81093800E-02, 0.65664310E-04,-0.17224123E-02,
     +  0.17639430E-02,-0.88922830E-03,-0.41402860E-02, 0.79062140E-03,
     +  0.10392813E-02, 0.14350280E-02,-0.12002340E-02,-0.13102302E-03,
     + -0.25189230E-03,-0.70312333E-03,-0.13841351E-01,-0.26033523E-04,
     + -0.71405484E-04,-0.11403022E-02,-0.24609852E-02,-0.18406561E-02,
     + -0.15960760E-04,-0.60296780E-03, 0.27188532E-02,-0.13216090E-02,
     + -0.37336200E-02,-0.14688734E-01,-0.18700620E-02,-0.10161433E-01,
     + -0.10688210E-01,-0.12883740E-02,-0.25070520E-03, 0.40832670E-04,
     +  0.17898553E-03, 0.55863772E-03, 0.80347014E-03, 0.79814570E-03,
     + -0.23060790E-03,-0.81397520E-02,-0.25635004E-01,-0.35789653E-01,
     + -0.18598620E-01,-0.23829701E-02, 0.56152550E-03,-0.50462163E-04,
     + -0.36264434E-02,-0.29698144E-02,-0.70642461E-02,-0.18890830E-03,
     +  0.10599090E-02, 0.13056040E-03, 0.25252460E-03,-0.18390400E-03,
     + -0.41455123E-03,-0.28816190E-03,-0.16006453E-03,-0.92927791E-03,
     + -0.24169800E-03,-0.29987632E-03,-0.49774880E-02,-0.67174260E-02,
     + -0.35247621E-02, 0.21602324E-05,-0.59878043E-02, 0.25901120E-03,
     +  0.38279700E-04, 0.30085640E-03, 0.11645530E-03,-0.36528780E-03,
     + -0.38135160E-03,-0.15801420E-03, 0.73823850E-03, 0.24482160E-03,
     + -0.25510940E-03,-0.17707874E-02, 0.37687180E-03,-0.43854292E-03,
     + -0.93689432E-03,-0.43322882E-03,-0.38333600E-03,-0.32136933E-02,
     +  0.55109942E-03, 0.45981420E-03,-0.45082940E-02,-0.35784791E-02,
     + -0.10936801E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_h_dext    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_h_dext    =x_h_dext    
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)*x11            *x51
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)*x12                
      x_h_dext    =x_h_dext    
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)*x11    *x31*x41    
     3  +coeff( 21)*x11        *x42    
     4  +coeff( 22)*x11*x21        *x51
     5  +coeff( 23)    *x22*x31*x41    
     6  +coeff( 24)        *x33*x41    
     7  +coeff( 25)        *x31*x43    
     8  +coeff( 26)    *x25            
      x_h_dext    =x_h_dext    
     9  +coeff( 27)    *x23    *x42    
     1  +coeff( 28)        *x35*x41    
     2  +coeff( 29)*x13    *x32        
     3  +coeff( 30)    *x25*x31*x41    
     4  +coeff( 31)        *x31*x41    
     5  +coeff( 32)    *x24            
     6  +coeff( 33)*x12*x21            
     7  +coeff( 34)*x12            *x51
     8  +coeff( 35)*x11*x23            
      x_h_dext    =x_h_dext    
     9  +coeff( 36)*x11*x21    *x42    
     1  +coeff( 37)    *x23*x32        
     2  +coeff( 38)    *x23*x31*x41    
     3  +coeff( 39)*x12*x22            
     4  +coeff( 40)*x11*x22*x31*x41    
     5  +coeff( 41)*x11*x22    *x42    
     6  +coeff( 42)*x11*x24*x32        
     7  +coeff( 43)        *x31*x41*x51
     8  +coeff( 44)*x11    *x32        
      x_h_dext    =x_h_dext    
     9  +coeff( 45)*x11            *x52
     1  +coeff( 46)    *x23        *x51
     2  +coeff( 47)    *x21*x31*x41*x51
     3  +coeff( 48)    *x21    *x42*x51
     4  +coeff( 49)            *x42*x52
     5  +coeff( 50)    *x21*x33*x41    
     6  +coeff( 51)    *x21*x32*x42    
     7  +coeff( 52)    *x21*x31*x43    
     8  +coeff( 53)    *x21    *x44    
      x_h_dext    =x_h_dext    
     9  +coeff( 54)*x11*x22*x32        
     1  +coeff( 55)*x13*x21            
     2  +coeff( 56)    *x22*x34        
     3  +coeff( 57)*x12*x21*x31*x41    
     4  +coeff( 58)*x12*x21    *x42    
     5  +coeff( 59)*x11*x24    *x42    
     6  +coeff( 60)*x12*x23*x32        
     7  +coeff( 61)    *x23            
     8  +coeff( 62)                *x53
      x_h_dext    =x_h_dext    
     9  +coeff( 63)    *x21*x32    *x51
     1  +coeff( 64)    *x22        *x52
     2  +coeff( 65)*x11*x21*x31*x41    
     3  +coeff( 66)*x11        *x42*x51
     4  +coeff( 67)*x13                
     5  +coeff( 68)    *x21*x34        
     6  +coeff( 69)    *x24        *x51
     7  +coeff( 70)*x11*x24            
     8  +coeff( 71)*x11    *x32*x42    
      x_h_dext    =x_h_dext    
     9  +coeff( 72)*x11    *x31*x43    
     1  +coeff( 73)*x11        *x44    
     2  +coeff( 74)        *x32    *x53
     3  +coeff( 75)*x11*x24*x31*x41    
     4  +coeff( 76)*x11*x22*x33*x41    
     5  +coeff( 77)*x11    *x35*x41    
     6  +coeff( 78)*x12        *x42*x52
     7  +coeff( 79)        *x32    *x51
     8  +coeff( 80)    *x22*x32        
      x_h_dext    =x_h_dext    
     9  +coeff( 81)        *x32*x42    
     1  +coeff( 82)        *x31*x41*x52
     2  +coeff( 83)    *x22    *x42*x51
     3  +coeff( 84)*x12        *x42    
     4  +coeff( 85)*x11    *x34        
     5  +coeff( 86)*x11    *x33*x41    
     6  +coeff( 87)*x11*x23        *x51
     7  +coeff( 88)*x12*x23            
     8  +coeff( 89)*x12*x21*x32        
      x_h_dext    =x_h_dext    
     9  +coeff( 90)*x13    *x31*x41    
     1  +coeff( 91)*x13        *x42    
     2  +coeff( 92)    *x25    *x42    
     3  +coeff( 93)    *x24*x31*x41*x51
     4  +coeff( 94)*x12*x22*x31*x41    
     5  +coeff( 95)*x12*x23    *x42    
     6  +coeff( 96)*x12*x25*x31*x41    
     7  +coeff( 97)*x13*x24*x32        
c
c
c
c  transform to TRANSPORT x
c
      x_h_dext=x_h_dext*0.866025

      return
      end
      function t_h_dext    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.6216832E-02/
      data xmin/
     1 -0.49916E-01,-0.58959E-01,-0.49982E-01,-0.38531E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49982E-01, 0.38531E-01, 0.49826E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.55481381E-02,-0.75245261E-01, 0.22721853E-01, 0.32675940E-01,
     + -0.35593710E-02,-0.46012131E-03, 0.16156964E-02,-0.26590062E-02,
     + -0.40323124E-02,-0.12267964E-02, 0.43484060E-02, 0.41600770E-03,
     +  0.29218893E-02, 0.53515210E-02,-0.32217593E-02,-0.29680712E-03,
     +  0.36674780E-03, 0.27162884E-03,-0.20359710E-02, 0.45333130E-03,
     + -0.13420510E-02, 0.13435070E-02, 0.33513930E-03,-0.17462710E-03,
     +  0.53759251E-03, 0.30378360E-03, 0.13470554E-02,-0.28976790E-03,
     + -0.12154070E-04, 0.18845072E-03, 0.30899743E-03,-0.42196712E-03,
     + -0.57354423E-03,-0.40518560E-03, 0.29123220E-02, 0.26330391E-02,
     +  0.10436640E-03,-0.26693860E-03, 0.41942360E-03,-0.40633682E-03,
     + -0.11988254E-03, 0.20725620E-03, 0.13602550E-03,-0.44014630E-03,
     + -0.15294323E-03, 0.62806080E-03, 0.22084390E-03, 0.19981210E-03,
     +  0.50560030E-03,-0.30015211E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_h_dext    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_h_dext    =t_h_dext    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)        *x32    *x51
      t_h_dext    =t_h_dext    
     9  +coeff( 18)        *x31*x41*x51
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)*x12                
     4  +coeff( 22)*x11*x22            
     5  +coeff( 23)*x11    *x32        
     6  +coeff( 24)    *x22*x32        
     7  +coeff( 25)*x11    *x31*x41    
     8  +coeff( 26)    *x22*x31*x41    
      t_h_dext    =t_h_dext    
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)*x11*x21        *x51
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)*x11            *x52
     4  +coeff( 31)*x12            *x51
     5  +coeff( 32)*x11*x22        *x51
     6  +coeff( 33)*x11    *x31*x41*x51
     7  +coeff( 34)*x11        *x42*x51
     8  +coeff( 35)*x11*x22*x31*x41    
      t_h_dext    =t_h_dext    
     9  +coeff( 36)*x11*x22    *x42    
     1  +coeff( 37)                *x53
     2  +coeff( 38)*x11        *x42    
     3  +coeff( 39)    *x21    *x42*x51
     4  +coeff( 40)            *x42*x52
     5  +coeff( 41)*x12*x21            
     6  +coeff( 42)*x11*x21*x32        
     7  +coeff( 43)*x11*x21*x31*x41    
     8  +coeff( 44)*x11*x21    *x42    
      t_h_dext    =t_h_dext    
     9  +coeff( 45)*x11    *x32    *x51
     1  +coeff( 46)*x11*x22*x32        
     2  +coeff( 47)*x12        *x42    
     3  +coeff( 48)*x12*x21        *x51
     4  +coeff( 49)    *x23*x31*x41*x51
     5  +coeff( 50)        *x31*x41*x54
c
      return
      end
      function y_h_dext    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49916E-01,-0.58959E-01,-0.49982E-01,-0.38531E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49982E-01, 0.38531E-01, 0.49826E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.50602160E-01, 0.14259460E+00,-0.19309570E-01,-0.67861070E-01,
     +  0.23020700E-01, 0.61587952E-01, 0.92712044E-02,-0.10923792E-01,
     + -0.86397920E-03, 0.31046760E-01,-0.32080110E-01,-0.55115194E-02,
     + -0.12290904E-01,-0.10068313E-01,-0.24390330E-02,-0.64786160E-02,
     + -0.20461811E-02,-0.43895510E-02, 0.59891150E-02, 0.96624053E-03,
     +  0.28479730E-03, 0.20380392E-01, 0.38850070E-02, 0.49907670E-02,
     +  0.46247170E-02, 0.67321980E-04,-0.21309391E-02,-0.37479230E-02,
     +  0.23133833E-02, 0.15228940E-02,-0.20574660E-02,-0.60559350E-02,
     + -0.15164044E-02,-0.10643273E-01,-0.24650574E-03,-0.52107190E-03,
     + -0.12000083E-02, 0.12912870E-02,-0.91066374E-03, 0.68999370E-03,
     +  0.41072684E-03, 0.83717713E-02, 0.23167203E-02, 0.42932620E-03,
     + -0.13522792E-03, 0.36317302E-04, 0.12017140E-02,-0.62238381E-04,
     + -0.31412404E-03, 0.17716480E-03,-0.23442781E-02, 0.82263723E-02,
     +  0.71081672E-02, 0.34165594E-02, 0.34245674E-02, 0.86206141E-02,
     + -0.45324690E-02,-0.37954543E-03,-0.98243660E-03, 0.10811181E-02,
     +  0.73461230E-03,-0.22511123E-02,-0.97766570E-02, 0.29357280E-02,
     +  0.68341814E-02, 0.17972020E-02, 0.48677740E-03, 0.21626980E-02,
     +  0.39435690E-02,-0.27764220E-02,-0.26515000E-02,-0.21660900E-02,
     + -0.59620390E-03,-0.35663654E-02,-0.70409760E-02,-0.77030621E-02,
     + -0.11113794E-02,-0.47033450E-02,-0.13506013E-01,-0.92746740E-02,
     + -0.23812013E-02, 0.16478110E-02, 0.14628200E-03, 0.28099540E-03,
     + -0.81511580E-03,-0.51080950E-03,-0.64381130E-02, 0.28851000E-03,
     + -0.14178170E-03, 0.21029103E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_h_dext    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      y_h_dext    =y_h_dext    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)        *x31    *x52
      y_h_dext    =y_h_dext    
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)*x11*x21*x31        
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)    *x21*x33        
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x21*x31*x42    
     7  +coeff( 25)    *x21    *x43    
     8  +coeff( 26)        *x33    *x51
      y_h_dext    =y_h_dext    
     9  +coeff( 27)*x11        *x41*x51
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)        *x31*x42*x51
     3  +coeff( 30)            *x43*x51
     4  +coeff( 31)*x12    *x31        
     5  +coeff( 32)*x12        *x41    
     6  +coeff( 33)*x11*x22    *x41    
     7  +coeff( 34)    *x22    *x43    
     8  +coeff( 35)*x11*x21*x31    *x51
      y_h_dext    =y_h_dext    
     9  +coeff( 36)*x11*x21    *x41*x51
     1  +coeff( 37)*x11*x23*x31        
     2  +coeff( 38)*x12*x21    *x41    
     3  +coeff( 39)*x11*x23    *x41    
     4  +coeff( 40)    *x23*x32*x41    
     5  +coeff( 41)    *x21*x34*x41    
     6  +coeff( 42)*x11*x21    *x43    
     7  +coeff( 43)    *x21*x32*x43    
     8  +coeff( 44)*x12    *x31    *x51
      y_h_dext    =y_h_dext    
     9  +coeff( 45)    *x24*x31    *x51
     1  +coeff( 46)*x11    *x33    *x51
     2  +coeff( 47)*x12        *x41*x51
     3  +coeff( 48)        *x34*x41*x51
     4  +coeff( 49)*x12*x22*x31        
     5  +coeff( 50)*x12*x22    *x41    
     6  +coeff( 51)    *x22*x31*x44    
     7  +coeff( 52)*x11*x23*x31*x42    
     8  +coeff( 53)*x11*x23    *x43    
      y_h_dext    =y_h_dext    
     9  +coeff( 54)*x11*x21*x31*x44    
     1  +coeff( 55)    *x21*x33*x44    
     2  +coeff( 56)*x11*x23*x32*x43    
     3  +coeff( 57)*x11*x23*x34*x43    
     4  +coeff( 58)*x11    *x31    *x51
     5  +coeff( 59)    *x22*x31    *x51
     6  +coeff( 60)        *x32*x41*x51
     7  +coeff( 61)    *x21    *x41*x52
     8  +coeff( 62)    *x24    *x41    
      y_h_dext    =y_h_dext    
     9  +coeff( 63)    *x22*x31*x42    
     1  +coeff( 64)*x11*x21*x32*x41    
     2  +coeff( 65)*x11*x21*x31*x42    
     3  +coeff( 66)    *x23*x31*x42    
     4  +coeff( 67)    *x23*x31    *x52
     5  +coeff( 68)    *x24*x32*x41    
     6  +coeff( 69)    *x22*x34*x41    
     7  +coeff( 70)*x12    *x31*x42    
     8  +coeff( 71)    *x24*x31*x42    
      y_h_dext    =y_h_dext    
     9  +coeff( 72)*x12        *x43    
     1  +coeff( 73)    *x24    *x43    
     2  +coeff( 74)*x12*x24*x31        
     3  +coeff( 75)*x12*x24    *x41    
     4  +coeff( 76)*x12*x22*x32*x41    
     5  +coeff( 77)*x12    *x34*x41    
     6  +coeff( 78)    *x24*x34*x41    
     7  +coeff( 79)*x12*x22*x31*x42    
     8  +coeff( 80)*x12*x22    *x43    
      y_h_dext    =y_h_dext    
     9  +coeff( 81)*x12*x24*x33        
     1  +coeff( 82)    *x21*x32*x41    
     2  +coeff( 83)        *x31    *x53
     3  +coeff( 84)            *x41*x53
     4  +coeff( 85)    *x24*x31        
     5  +coeff( 86)    *x22*x33        
     6  +coeff( 87)    *x22*x32*x41    
     7  +coeff( 88)*x11    *x31*x42    
     8  +coeff( 89)*x11        *x41*x52
      y_h_dext    =y_h_dext    
     9  +coeff( 90)*x12*x21*x31        
c
      return
      end
      function p_h_dext    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49916E-01,-0.58959E-01,-0.49982E-01,-0.38531E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49982E-01, 0.38531E-01, 0.49826E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.18960190E-01, 0.91240390E-02,-0.38409530E-02,-0.15517221E-01,
     +  0.44016381E-02, 0.14154990E-01, 0.19985410E-02,-0.18697390E-02,
     + -0.19497550E-03, 0.66359303E-02,-0.89951930E-02,-0.13051284E-02,
     + -0.31395950E-02,-0.23363700E-02,-0.32824501E-03,-0.33513973E-02,
     + -0.87877080E-04,-0.59719590E-03, 0.91248740E-03, 0.39616713E-03,
     +  0.41717210E-04, 0.69428351E-02, 0.16841380E-02, 0.14068700E-02,
     +  0.60504320E-03,-0.22675300E-03,-0.29366012E-03,-0.26322513E-04,
     +  0.22568450E-03,-0.13656911E-02, 0.21174810E-03, 0.36852710E-03,
     + -0.32235190E-03,-0.15265940E-03,-0.18141680E-03,-0.19532610E-02,
     + -0.15825614E-02,-0.11857440E-02,-0.82992450E-03, 0.78477023E-04,
     +  0.12707780E-03, 0.43840892E-03, 0.16141180E-03, 0.15281300E-03,
     + -0.89017110E-03, 0.93947830E-03,-0.10936330E-02, 0.17347280E-03,
     +  0.40316913E-03, 0.35242780E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_h_dext    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      p_h_dext    =p_h_dext    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)        *x31    *x52
      p_h_dext    =p_h_dext    
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)*x11*x21*x31        
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)    *x21*x33        
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x21*x31*x42    
     7  +coeff( 25)    *x21    *x43    
     8  +coeff( 26)*x11    *x31    *x51
      p_h_dext    =p_h_dext    
     9  +coeff( 27)    *x22*x31    *x51
     1  +coeff( 28)        *x33    *x51
     2  +coeff( 29)*x11        *x41*x51
     3  +coeff( 30)    *x22    *x41*x51
     4  +coeff( 31)        *x31*x42*x51
     5  +coeff( 32)            *x43*x51
     6  +coeff( 33)    *x21    *x41*x52
     7  +coeff( 34)            *x41*x53
     8  +coeff( 35)*x11*x22*x31        
      p_h_dext    =p_h_dext    
     9  +coeff( 36)*x12        *x41    
     1  +coeff( 37)*x11*x22    *x41    
     2  +coeff( 38)    *x22*x31*x42    
     3  +coeff( 39)    *x22    *x43    
     4  +coeff( 40)    *x21*x33    *x51
     5  +coeff( 41)*x11*x21    *x41*x51
     6  +coeff( 42)    *x23    *x41*x51
     7  +coeff( 43)*x11        *x41*x52
     8  +coeff( 44)    *x22    *x41*x52
      p_h_dext    =p_h_dext    
     9  +coeff( 45)*x11*x23*x31        
     1  +coeff( 46)*x12*x21    *x41    
     2  +coeff( 47)*x11*x23    *x41    
     3  +coeff( 48)    *x23*x32*x41    
     4  +coeff( 49)    *x21*x34*x41    
     5  +coeff( 50)*x11*x21*x31*x42    
c
      return
      end
      function l_h_dext    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/ -0.2567878E-01/
      data xmin/
     1 -0.49916E-01,-0.58959E-01,-0.49982E-01,-0.37177E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49826E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21673450E-01,-0.47936591E+00,-0.97724340E-01, 0.25024491E+00,
     + -0.45174933E-01,-0.30488220E-01, 0.26749360E-01,-0.14659230E-01,
     + -0.74208200E-02, 0.27898020E-01, 0.93553730E-02, 0.14100160E-02,
     + -0.51825220E-02, 0.16935980E-02, 0.58924853E-02, 0.24027714E-01,
     + -0.28118474E-02, 0.47213863E-02,-0.41700880E-03, 0.70988390E-02,
     + -0.16353040E-02,-0.48241741E-02, 0.20944613E-02, 0.29364982E-02,
     + -0.17072630E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_h_dext    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)*x11            *x51
      l_h_dext    =l_h_dext    
     9  +coeff(  9)*x12                
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)*x11*x22            
     3  +coeff( 12)    *x23*x31*x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x22        *x51
      l_h_dext    =l_h_dext    
     9  +coeff( 18)*x11    *x31*x41    
     1  +coeff( 19)    *x22*x31*x41    
     2  +coeff( 20)*x11*x22    *x42    
     3  +coeff( 21)        *x32        
     4  +coeff( 22)        *x31*x41    
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)*x11        *x42    
     7  +coeff( 25)*x11*x21        *x51
c
      return
      end
      function x_h_q3en    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 99)
      data ncoeff/ 98/
      data avdat/  0.1545243E-02/
      data xmin/
     1 -0.49916E-01,-0.58628E-01,-0.49967E-01,-0.38531E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49826E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.77329730E-03, 0.21630164E+00, 0.13871500E+00,-0.13417391E+00,
     +  0.18286220E-01,-0.65355771E-03,-0.64108250E-02, 0.34822683E-01,
     + -0.58757100E-02,-0.15288072E-01, 0.14458730E-02,-0.12527800E-02,
     + -0.41857040E-02,-0.72001153E-02, 0.28312101E-02, 0.50542290E-04,
     + -0.15452080E-02, 0.30631350E-02,-0.31271530E-03,-0.42366003E-02,
     +  0.13704060E-02,-0.19493120E-02,-0.16946010E-03,-0.29157810E-03,
     + -0.90850890E-02,-0.17438424E-01, 0.19607710E-04,-0.58228490E-02,
     + -0.10971324E-02,-0.77365800E-02, 0.45200133E-04, 0.14995654E-03,
     + -0.49659230E-05,-0.23602300E-02, 0.38652100E-03,-0.38155730E-03,
     + -0.17574740E-02, 0.23796290E-03,-0.25252540E-03, 0.24281230E-02,
     + -0.19537531E-02,-0.10753870E-01,-0.24287880E-01,-0.12654160E-01,
     + -0.18565881E-02,-0.22522770E-02, 0.17036274E-03,-0.38788552E-04,
     +  0.16088390E-02, 0.35058160E-03,-0.52110570E-03, 0.50882860E-03,
     +  0.83464244E-03,-0.34807130E-03,-0.23871834E-02,-0.56810840E-02,
     +  0.12454251E-02, 0.53988554E-03,-0.55283020E-03,-0.81479240E-03,
     +  0.60194000E-03,-0.25197370E-02,-0.17424570E-02,-0.11878252E-02,
     +  0.91649992E-04, 0.17085333E-03, 0.22636790E-03, 0.27248600E-03,
     + -0.36953170E-03,-0.60302300E-03, 0.13393240E-03,-0.22083550E-03,
     + -0.33212904E-02,-0.44209514E-02,-0.23393730E-02, 0.64589430E-03,
     + -0.23565990E-02, 0.69624000E-03,-0.49677550E-02,-0.34095013E-02,
     + -0.96234300E-03,-0.42267510E-02,-0.35640510E-03,-0.18494242E-03,
     + -0.27679893E-03,-0.10697121E-03,-0.18181423E-03,-0.42278863E-03,
     +  0.65333490E-04,-0.13015144E-03, 0.27800531E-03, 0.75959462E-04,
     + -0.33249620E-03,-0.38549370E-03,-0.15633754E-03,-0.21882240E-03,
     +  0.80657552E-03,-0.16002500E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_h_q3en    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      x_h_q3en    =x_h_q3en    
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)*x11            *x51
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      x_h_q3en    =x_h_q3en    
     9  +coeff( 18)*x12                
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)*x11        *x42    
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)        *x33*x41    
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)    *x21*x32*x42    
      x_h_q3en    =x_h_q3en    
     9  +coeff( 27)        *x31*x43*x51
     1  +coeff( 28)*x11*x22*x31*x41    
     2  +coeff( 29)*x11    *x33*x41    
     3  +coeff( 30)*x11*x22    *x42    
     4  +coeff( 31)        *x35*x41    
     5  +coeff( 32)*x11*x21*x32*x42    
     6  +coeff( 33)*x13    *x32        
     7  +coeff( 34)        *x31*x41    
     8  +coeff( 35)                *x53
      x_h_q3en    =x_h_q3en    
     9  +coeff( 36)*x11*x21        *x51
     1  +coeff( 37)    *x24            
     2  +coeff( 38)*x11            *x52
     3  +coeff( 39)*x12            *x51
     4  +coeff( 40)*x11*x23            
     5  +coeff( 41)*x11*x21    *x42    
     6  +coeff( 42)    *x23*x31*x41    
     7  +coeff( 43)    *x21*x31*x43    
     8  +coeff( 44)    *x21    *x44    
      x_h_q3en    =x_h_q3en    
     9  +coeff( 45)*x12*x22            
     1  +coeff( 46)*x11*x22*x32        
     2  +coeff( 47)    *x21*x35*x41    
     3  +coeff( 48)    *x25*x34        
     4  +coeff( 49)    *x23            
     5  +coeff( 50)*x11    *x31*x41    
     6  +coeff( 51)    *x22*x32        
     7  +coeff( 52)    *x23        *x51
     8  +coeff( 53)*x12*x21            
      x_h_q3en    =x_h_q3en    
     9  +coeff( 54)*x13                
     1  +coeff( 55)    *x23*x32        
     2  +coeff( 56)    *x21*x33*x41    
     3  +coeff( 57)    *x22    *x42*x51
     4  +coeff( 58)*x12        *x42    
     5  +coeff( 59)*x11*x24            
     6  +coeff( 60)*x11*x21    *x42*x51
     7  +coeff( 61)*x13*x21            
     8  +coeff( 62)*x12*x21*x31*x41    
      x_h_q3en    =x_h_q3en    
     9  +coeff( 63)*x12*x21    *x42    
     1  +coeff( 64)*x12*x23*x32        
     2  +coeff( 65)        *x32    *x51
     3  +coeff( 66)*x11    *x32        
     4  +coeff( 67)    *x21*x32    *x51
     5  +coeff( 68)    *x21*x31*x41*x51
     6  +coeff( 69)    *x25            
     7  +coeff( 70)    *x21*x34        
     8  +coeff( 71)*x11*x21        *x52
      x_h_q3en    =x_h_q3en    
     9  +coeff( 72)*x12*x21        *x51
     1  +coeff( 73)*x11    *x32*x42    
     2  +coeff( 74)*x11    *x31*x43    
     3  +coeff( 75)*x11        *x44    
     4  +coeff( 76)*x11*x23        *x51
     5  +coeff( 77)    *x25    *x42    
     6  +coeff( 78)    *x24*x31*x41*x51
     7  +coeff( 79)*x11*x24*x31*x41    
     8  +coeff( 80)*x11*x24    *x42    
      x_h_q3en    =x_h_q3en    
     9  +coeff( 81)*x11*x23*x31*x41*x51
     1  +coeff( 82)*x12*x23    *x42    
     2  +coeff( 83)        *x32*x42    
     3  +coeff( 84)        *x31*x41*x52
     4  +coeff( 85)            *x42*x52
     5  +coeff( 86)*x11    *x32    *x51
     6  +coeff( 87)*x11        *x42*x51
     7  +coeff( 88)    *x24        *x51
     8  +coeff( 89)*x12    *x31*x41    
      x_h_q3en    =x_h_q3en    
     9  +coeff( 90)*x11    *x34        
     1  +coeff( 91)*x11    *x31*x41*x52
     2  +coeff( 92)*x13            *x51
     3  +coeff( 93)*x12*x21*x32        
     4  +coeff( 94)*x12*x22        *x51
     5  +coeff( 95)*x13    *x31*x41    
     6  +coeff( 96)    *x21*x31*x41*x54
     7  +coeff( 97)*x11*x22    *x42*x52
     8  +coeff( 98)*x12*x23*x31*x41    
c
      return
      end
      function t_h_q3en    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.5499774E-02/
      data xmin/
     1 -0.49916E-01,-0.58628E-01,-0.49967E-01,-0.38531E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49826E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.45958170E-02,-0.75638250E-01, 0.22695930E-01, 0.33088760E-01,
     + -0.48401914E-02, 0.14259594E-03,-0.79394143E-03, 0.87711310E-03,
     + -0.51435540E-02,-0.14174223E-02, 0.57899510E-02, 0.17929410E-03,
     +  0.75097180E-03, 0.28025873E-02, 0.16121480E-02,-0.26298850E-02,
     + -0.35001503E-03, 0.76226510E-03, 0.56232354E-03,-0.17191210E-02,
     +  0.10170110E-02, 0.10618010E-02,-0.53446730E-03,-0.93313492E-03,
     +  0.43135620E-03, 0.18926540E-02, 0.15249264E-03, 0.22483570E-02,
     +  0.78318660E-04, 0.17763690E-03, 0.47250190E-02, 0.15365864E-03,
     + -0.14175492E-03, 0.13326133E-03, 0.19970170E-03, 0.32369880E-03,
     +  0.70541421E-03,-0.18538880E-03,-0.89714792E-03, 0.14642410E-03,
     +  0.37185533E-03, 0.12494022E-03, 0.12160890E-02,-0.57595554E-03,
     +  0.82102441E-03, 0.57828180E-03, 0.20995593E-03,-0.23970040E-03,
     +  0.41444031E-02, 0.72250050E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_h_q3en    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_h_q3en    =t_h_q3en    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x22        *x51
      t_h_q3en    =t_h_q3en    
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)*x12                
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)*x11        *x42    
     5  +coeff( 23)*x11*x21        *x51
     6  +coeff( 24)    *x21    *x42*x51
     7  +coeff( 25)*x12            *x51
     8  +coeff( 26)*x11*x22*x31*x41    
      t_h_q3en    =t_h_q3en    
     9  +coeff( 27)*x11    *x33*x41    
     1  +coeff( 28)*x11*x22    *x42    
     2  +coeff( 29)    *x21*x33*x41*x51
     3  +coeff( 30)*x12*x21*x32        
     4  +coeff( 31)    *x23    *x44    
     5  +coeff( 32)*x11*x22*x34        
     6  +coeff( 33)        *x32    *x51
     7  +coeff( 34)                *x53
     8  +coeff( 35)*x11    *x32        
      t_h_q3en    =t_h_q3en    
     9  +coeff( 36)*x11    *x31*x41    
     1  +coeff( 37)    *x22    *x42    
     2  +coeff( 38)    *x23        *x51
     3  +coeff( 39)    *x21*x31*x41*x51
     4  +coeff( 40)*x12*x21            
     5  +coeff( 41)*x11*x23            
     6  +coeff( 42)*x11*x21*x31*x41    
     7  +coeff( 43)    *x23*x31*x41    
     8  +coeff( 44)*x11*x21    *x42    
      t_h_q3en    =t_h_q3en    
     9  +coeff( 45)    *x23    *x42    
     1  +coeff( 46)    *x21*x31*x43    
     2  +coeff( 47)*x11    *x31*x41*x51
     3  +coeff( 48)*x12*x22            
     4  +coeff( 49)    *x23*x32*x42    
     5  +coeff( 50)    *x23*x31*x43    
c
      return
      end
      function y_h_q3en    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49916E-01,-0.58628E-01,-0.49967E-01,-0.38531E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49826E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.69441124E-01, 0.15574613E+00,-0.21445790E-01,-0.86071953E-01,
     +  0.29134720E-01, 0.75872264E-01, 0.10116431E-01,-0.10453602E-01,
     + -0.10129090E-02, 0.38963861E-01,-0.37689924E-01,-0.59417830E-02,
     + -0.14034762E-01,-0.11806081E-01,-0.11126561E-01,-0.25257414E-02,
     + -0.58982870E-02, 0.59704690E-02, 0.19283990E-02, 0.35974330E-03,
     +  0.24963650E-01, 0.72971894E-02, 0.77002262E-02, 0.55632740E-02,
     + -0.43591670E-03,-0.14168580E-02, 0.17879210E-03,-0.10871050E-02,
     + -0.36907230E-02, 0.28407101E-02, 0.21834683E-02, 0.59430140E-03,
     + -0.17074394E-02,-0.66472981E-02,-0.60869110E-02,-0.18683081E-01,
     +  0.44839413E-04,-0.14822310E-01,-0.11909910E-02,-0.73235231E-03,
     + -0.31804721E-03,-0.13912690E-02, 0.24714670E-04,-0.19766682E-03,
     + -0.21552830E-03,-0.11182590E-02, 0.34651241E-02, 0.63698800E-03,
     +  0.32086140E-02, 0.47902480E-03,-0.45431820E-03, 0.11627894E-01,
     +  0.11079360E-01, 0.15944360E-02, 0.15570040E-02,-0.20913030E-03,
     + -0.50010820E-02,-0.52018640E-02,-0.94476220E-03, 0.86238130E-03,
     +  0.36660220E-03,-0.22559960E-02, 0.27327020E-02, 0.14854740E-02,
     + -0.24889840E-02,-0.12276862E-02,-0.33713611E-02,-0.73008383E-02,
     +  0.25013112E-03, 0.26925310E-03, 0.45178091E-03,-0.57215914E-02,
     + -0.41018320E-02,-0.45882920E-02,-0.20982760E-02,-0.19152960E-03,
     +  0.20025444E-02, 0.25066160E-03,-0.59026360E-03,-0.10331650E-02,
     +  0.13223890E-02, 0.10555380E-03, 0.44068850E-03,-0.99222770E-03,
     +  0.24392781E-03,-0.43366840E-03,-0.23166004E-02, 0.16895231E-02,
     + -0.41576864E-03,-0.39771283E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_h_q3en    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      y_h_q3en    =y_h_q3en    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)            *x41*x52
      y_h_q3en    =y_h_q3en    
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)    *x22*x31    *x51
      y_h_q3en    =y_h_q3en    
     9  +coeff( 27)        *x33    *x51
     1  +coeff( 28)*x11        *x41*x51
     2  +coeff( 29)    *x22    *x41*x51
     3  +coeff( 30)        *x31*x42*x51
     4  +coeff( 31)            *x43*x51
     5  +coeff( 32)    *x21    *x41*x52
     6  +coeff( 33)*x12    *x31        
     7  +coeff( 34)*x12        *x41    
     8  +coeff( 35)*x11*x22    *x41    
      y_h_q3en    =y_h_q3en    
     9  +coeff( 36)    *x22*x31*x42    
     1  +coeff( 37)*x11        *x43    
     2  +coeff( 38)    *x22    *x43    
     3  +coeff( 39)        *x31*x44    
     4  +coeff( 40)    *x23*x31    *x51
     5  +coeff( 41)    *x21*x33    *x51
     6  +coeff( 42)*x11*x21    *x41*x51
     7  +coeff( 43)    *x21*x31*x42*x51
     8  +coeff( 44)*x11        *x41*x52
      y_h_q3en    =y_h_q3en    
     9  +coeff( 45)    *x21*x31    *x53
     1  +coeff( 46)*x11*x23*x31        
     2  +coeff( 47)*x12*x21    *x41    
     3  +coeff( 48)*x11*x23    *x41    
     4  +coeff( 49)*x11*x21*x32*x41    
     5  +coeff( 50)    *x23*x32*x41    
     6  +coeff( 51)    *x21*x34*x41    
     7  +coeff( 52)*x11*x21*x31*x42    
     8  +coeff( 53)*x11*x21    *x43    
      y_h_q3en    =y_h_q3en    
     9  +coeff( 54)    *x21*x32*x43    
     1  +coeff( 55)*x12        *x41*x51
     2  +coeff( 56)        *x34*x41*x51
     3  +coeff( 57)*x12*x22*x31        
     4  +coeff( 58)*x12*x22    *x41    
     5  +coeff( 59)    *x24*x32*x41    
     6  +coeff( 60)    *x23*x33    *x51
     7  +coeff( 61)    *x23*x31    *x53
     8  +coeff( 62)    *x21*x31    *x51
      y_h_q3en    =y_h_q3en    
     9  +coeff( 63)    *x21*x32*x41    
     1  +coeff( 64)        *x32*x41*x51
     2  +coeff( 65)    *x24*x31        
     3  +coeff( 66)    *x22*x33        
     4  +coeff( 67)    *x24    *x41    
     5  +coeff( 68)    *x22*x32*x41    
     6  +coeff( 69)*x11*x21*x31    *x51
     7  +coeff( 70)*x12    *x31    *x51
     8  +coeff( 71)*x11*x22*x31    *x51
      y_h_q3en    =y_h_q3en    
     9  +coeff( 72)*x12    *x31*x42    
     1  +coeff( 73)*x12        *x43    
     2  +coeff( 74)*x12*x24    *x41    
     3  +coeff( 75)*x12*x22*x32*x41    
     4  +coeff( 76)*x12    *x34*x41    
     5  +coeff( 77)*x12*x22*x34*x41    
     6  +coeff( 78)            *x41*x53
     7  +coeff( 79)*x11*x22*x31        
     8  +coeff( 80)        *x32*x43    
      y_h_q3en    =y_h_q3en    
     9  +coeff( 81)    *x23    *x41*x51
     1  +coeff( 82)*x11    *x31    *x52
     2  +coeff( 83)*x12*x21*x31        
     3  +coeff( 84)*x11*x22    *x41*x51
     4  +coeff( 85)*x11    *x32*x41*x51
     5  +coeff( 86)*x11    *x31*x42*x51
     6  +coeff( 87)*x12    *x32*x41    
     7  +coeff( 88)*x11*x22    *x43    
     8  +coeff( 89)*x12*x21*x31    *x51
      y_h_q3en    =y_h_q3en    
     9  +coeff( 90)*x12        *x41*x52
c
      return
      end
      function p_h_q3en    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49916E-01,-0.58628E-01,-0.49967E-01,-0.38531E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49826E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.18889920E-01, 0.15121280E-01,-0.54687900E-02,-0.15849500E-01,
     +  0.46873474E-02, 0.16209310E-01, 0.28210720E-02,-0.27242521E-02,
     + -0.19715450E-03, 0.66788000E-02,-0.10472200E-01,-0.14575920E-02,
     + -0.34933670E-02,-0.28676760E-02,-0.42121460E-03,-0.23965780E-02,
     + -0.67393630E-04,-0.56079903E-03, 0.16231250E-02, 0.31699620E-03,
     +  0.11468532E-03, 0.79985493E-02, 0.13957672E-02, 0.10673890E-02,
     +  0.50950660E-03,-0.37304930E-03, 0.62636150E-05,-0.33016441E-03,
     + -0.98392380E-03, 0.51532380E-03, 0.14685624E-03,-0.12710530E-03,
     + -0.57542731E-03, 0.14252450E-03,-0.24774684E-02,-0.15237030E-02,
     + -0.60778763E-03, 0.17635370E-03,-0.36724103E-03, 0.53371820E-03,
     + -0.11664393E-04, 0.27549650E-03, 0.27311750E-03, 0.23329333E-03,
     + -0.79141850E-05,-0.73091051E-03,-0.17507610E-03, 0.91353000E-03,
     + -0.73992120E-03, 0.84025674E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_h_q3en    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      p_h_q3en    =p_h_q3en    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)        *x31    *x52
      p_h_q3en    =p_h_q3en    
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)*x11*x21*x31        
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)    *x21*x33        
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x21*x31*x42    
     7  +coeff( 25)    *x21    *x43    
     8  +coeff( 26)*x11    *x31    *x51
      p_h_q3en    =p_h_q3en    
     9  +coeff( 27)        *x33    *x51
     1  +coeff( 28)*x11        *x41*x51
     2  +coeff( 29)    *x22    *x41*x51
     3  +coeff( 30)        *x31*x42*x51
     4  +coeff( 31)    *x21    *x41*x52
     5  +coeff( 32)            *x41*x53
     6  +coeff( 33)*x12    *x31        
     7  +coeff( 34)*x11*x22*x31        
     8  +coeff( 35)*x12        *x41    
      p_h_q3en    =p_h_q3en    
     9  +coeff( 36)*x11*x22    *x41    
     1  +coeff( 37)    *x22    *x43    
     2  +coeff( 38)*x11*x21*x31    *x51
     3  +coeff( 39)*x11*x21    *x41*x51
     4  +coeff( 40)    *x23    *x41*x51
     5  +coeff( 41)*x11        *x41*x52
     6  +coeff( 42)    *x22    *x41*x52
     7  +coeff( 43)        *x31*x42*x52
     8  +coeff( 44)    *x21*x31    *x53
      p_h_q3en    =p_h_q3en    
     9  +coeff( 45)*x12*x21*x31        
     1  +coeff( 46)*x11*x23*x31        
     2  +coeff( 47)*x11*x21*x33        
     3  +coeff( 48)*x12*x21    *x41    
     4  +coeff( 49)*x11*x23    *x41    
     5  +coeff( 50)    *x23*x32*x41    
c
      return
      end
      function l_h_q3en    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 29)
      data ncoeff/ 28/
      data avdat/ -0.1907076E-01/
      data xmin/
     1 -0.49916E-01,-0.58628E-01,-0.49917E-01,-0.37177E-01,-0.49961E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49826E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15514072E-01,-0.30898400E+00, 0.15272383E+00,-0.31423691E-01,
     + -0.20870191E-01, 0.25153854E-01,-0.12929610E-01, 0.86179380E-03,
     + -0.28693664E-01,-0.71570300E-02,-0.80879100E-02, 0.17082344E-01,
     + -0.32651112E-02, 0.72952923E-02, 0.27650460E-03,-0.52839960E-03,
     +  0.99919950E-03,-0.20584962E-02,-0.61424644E-02, 0.31451084E-02,
     +  0.13729922E-01, 0.39322050E-02, 0.34898680E-02, 0.12647962E-02,
     +  0.17908140E-02, 0.13308650E-02, 0.10856330E-02, 0.10297382E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_h_q3en    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)*x11*x21            
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)    *x24            
      l_h_q3en    =l_h_q3en    
     9  +coeff(  9)    *x22            
     1  +coeff( 10)*x12                
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)        *x33*x41    
     7  +coeff( 16)        *x32*x42    
     8  +coeff( 17)    *x23*x31*x41    
      l_h_q3en    =l_h_q3en    
     9  +coeff( 18)        *x32        
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)    *x21*x32        
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)*x11        *x42    
     6  +coeff( 24)*x11            *x52
     7  +coeff( 25)        *x31*x41*x51
     8  +coeff( 26)            *x42*x51
      l_h_q3en    =l_h_q3en    
     9  +coeff( 27)    *x21        *x52
     1  +coeff( 28)*x11    *x32        
c
      return
      end
      function x_h_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 99)
      data ncoeff/ 98/
      data avdat/ -0.7321482E-02/
      data xmin/
     1 -0.49916E-01,-0.58195E-01,-0.49967E-01,-0.38531E-01,-0.49743E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49756E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.85125593E-02, 0.68434163E-01, 0.31092770E+00,-0.94056554E-01,
     +  0.11734900E-01,-0.87931650E-03,-0.72611410E-02,-0.81220390E-02,
     +  0.33308480E-01,-0.20670920E-01,-0.35338094E-02,-0.43853670E-02,
     +  0.10755421E-02,-0.15083140E-01, 0.30371972E-02, 0.12657180E-02,
     +  0.76557470E-04,-0.14904210E-02, 0.24118702E-04,-0.27431342E-02,
     + -0.10004043E-02, 0.50260510E-02, 0.21329990E-02,-0.19646360E-02,
     + -0.16624020E-02, 0.70194522E-02, 0.63076760E-03,-0.16438710E-02,
     + -0.24316061E-02, 0.99582540E-03, 0.70932712E-02, 0.76309870E-03,
     + -0.86443973E-02,-0.11634913E-02,-0.67415480E-02,-0.33011562E-02,
     + -0.29497502E-02,-0.10285804E-01,-0.13009724E-01,-0.64565570E-02,
     + -0.52371430E-02, 0.20100780E-02,-0.13938114E-02,-0.19119360E-02,
     + -0.27929721E-02,-0.51761954E-03,-0.38457680E-02, 0.61458611E-03,
     + -0.24599420E-03, 0.44297250E-02, 0.15430660E-02,-0.11495313E-02,
     +  0.42247481E-03, 0.72377902E-03, 0.71023820E-04, 0.24363050E-02,
     +  0.48219350E-03,-0.24503240E-02,-0.19756704E-03, 0.43940641E-02,
     + -0.36211710E-03,-0.33455200E-03,-0.13994810E-02,-0.61251500E-02,
     + -0.29966381E-02,-0.26027884E-03,-0.24723920E-02,-0.20179950E-02,
     + -0.22093213E-03,-0.37895271E-02, 0.28508420E-02,-0.15402040E-02,
     + -0.79270750E-03,-0.50895903E-02,-0.76460830E-03,-0.30015134E-02,
     +  0.34114983E-03,-0.10134723E-02, 0.28219510E-03,-0.31957310E-03,
     +  0.11437492E-02,-0.89473481E-03,-0.16379980E-02,-0.41833793E-03,
     + -0.15020744E-03, 0.34763600E-02,-0.27403670E-02,-0.95796530E-03,
     +  0.13305180E-02, 0.17644452E-02,-0.55500020E-03,-0.82447200E-03,
     + -0.27176330E-03,-0.32278930E-03, 0.33347294E-03,-0.22114190E-03,
     +  0.19784591E-03, 0.10015310E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_h_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      x_h_q3ex    =x_h_q3ex    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)*x11            *x51
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      x_h_q3ex    =x_h_q3ex    
     9  +coeff( 18)*x12                
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)*x11    *x31*x41    
     4  +coeff( 22)*x11        *x42    
     5  +coeff( 23)                *x53
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)    *x22*x31*x41    
     8  +coeff( 26)    *x22    *x42    
      x_h_q3ex    =x_h_q3ex    
     9  +coeff( 27)*x11            *x52
     1  +coeff( 28)    *x21*x31*x41*x51
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)*x12            *x51
     4  +coeff( 31)*x11*x23            
     5  +coeff( 32)*x11*x21*x31*x41    
     6  +coeff( 33)*x11*x21    *x42    
     7  +coeff( 34)    *x23*x32        
     8  +coeff( 35)    *x23*x31*x41    
      x_h_q3ex    =x_h_q3ex    
     9  +coeff( 36)    *x21*x33*x41    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x21*x32*x42    
     3  +coeff( 39)    *x21*x31*x43    
     4  +coeff( 40)    *x21    *x44    
     5  +coeff( 41)*x12*x22            
     6  +coeff( 42)*x12        *x42    
     7  +coeff( 43)    *x24    *x42    
     8  +coeff( 44)    *x21*x32        
      x_h_q3ex    =x_h_q3ex    
     9  +coeff( 45)    *x21*x31*x41    
     1  +coeff( 46)        *x32    *x51
     2  +coeff( 47)    *x24            
     3  +coeff( 48)*x12*x21            
     4  +coeff( 49)*x13                
     5  +coeff( 50)*x11*x23        *x51
     6  +coeff( 51)*x13*x21            
     7  +coeff( 52)*x11    *x32*x42*x51
     8  +coeff( 53)*x11    *x32        
      x_h_q3ex    =x_h_q3ex    
     9  +coeff( 54)    *x21*x32    *x51
     1  +coeff( 55)    *x21        *x53
     2  +coeff( 56)*x11*x22        *x51
     3  +coeff( 57)*x11    *x31*x41*x51
     4  +coeff( 58)    *x24        *x51
     5  +coeff( 59)    *x22*x31*x41*x51
     6  +coeff( 60)    *x22    *x42*x51
     7  +coeff( 61)*x12    *x31*x41    
     8  +coeff( 62)*x11            *x53
      x_h_q3ex    =x_h_q3ex    
     9  +coeff( 63)*x12*x21        *x51
     1  +coeff( 64)*x11*x22    *x42    
     2  +coeff( 65)*x11*x21    *x42*x51
     3  +coeff( 66)    *x25        *x51
     4  +coeff( 67)*x12*x21*x31*x41    
     5  +coeff( 68)*x12*x22        *x51
     6  +coeff( 69)*x11    *x34    *x51
     7  +coeff( 70)*x12*x22    *x42    
     8  +coeff( 71)    *x23    *x42*x52
      x_h_q3ex    =x_h_q3ex    
     9  +coeff( 72)*x11*x24*x32        
     1  +coeff( 73)*x11*x22*x34        
     2  +coeff( 74)*x11*x24*x31*x41    
     3  +coeff( 75)*x12*x25            
     4  +coeff( 76)*x12*x23    *x42    
     5  +coeff( 77)*x11*x22        *x54
     6  +coeff( 78)    *x23        *x51
     7  +coeff( 79)        *x32    *x52
     8  +coeff( 80)*x11*x21*x32        
      x_h_q3ex    =x_h_q3ex    
     9  +coeff( 81)*x11*x21        *x52
     1  +coeff( 82)        *x32*x42*x51
     2  +coeff( 83)    *x23        *x52
     3  +coeff( 84)    *x21*x32    *x52
     4  +coeff( 85)    *x24        *x52
     5  +coeff( 86)*x11*x23    *x42    
     6  +coeff( 87)    *x25    *x42    
     7  +coeff( 88)*x11*x21    *x42*x52
     8  +coeff( 89)*x13*x21    *x42    
      x_h_q3ex    =x_h_q3ex    
     9  +coeff( 90)*x11*x24        *x52
     1  +coeff( 91)*x11*x23        *x53
     2  +coeff( 92)    *x22        *x52
     3  +coeff( 93)            *x42*x52
     4  +coeff( 94)*x11    *x32    *x51
     5  +coeff( 95)*x11        *x42*x51
     6  +coeff( 96)    *x21*x34        
     7  +coeff( 97)        *x34    *x51
     8  +coeff( 98)*x12    *x32        
c
      return
      end
      function t_h_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.5844384E-02/
      data xmin/
     1 -0.49916E-01,-0.58195E-01,-0.49967E-01,-0.38531E-01,-0.49743E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49756E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.54963952E-02,-0.27682254E-01, 0.97763434E-01,-0.50869930E-02,
     +  0.21789073E-04,-0.18883910E-03,-0.25550744E-02,-0.17508114E-02,
     +  0.41388063E-02,-0.90707880E-02, 0.29365450E-02, 0.28634641E-03,
     + -0.42606680E-03, 0.76937161E-04,-0.39965950E-02,-0.15629410E-02,
     + -0.34280680E-03, 0.92364114E-03, 0.21363160E-03,-0.31765043E-03,
     +  0.54733870E-03,-0.14843130E-02,-0.23120630E-03,-0.13836820E-03,
     + -0.24534290E-04, 0.21115874E-02, 0.24522260E-02, 0.49547130E-04,
     + -0.77218360E-03, 0.48066931E-03, 0.13156800E-03, 0.18360110E-02,
     +  0.32994104E-03,-0.25274210E-02, 0.35574834E-03,-0.12976574E-02,
     +  0.14069363E-04, 0.70910702E-03,-0.93831564E-03,-0.31838050E-03,
     + -0.80899880E-03,-0.38829501E-03, 0.17607193E-03, 0.40015563E-03,
     +  0.40938850E-03, 0.87551143E-03, 0.24421524E-03, 0.14840050E-03,
     + -0.23383460E-03, 0.13768620E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_h_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_h_q3ex    =t_h_q3ex    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)        *x32    *x51
      t_h_q3ex    =t_h_q3ex    
     9  +coeff( 18)        *x31*x41*x51
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)                *x53
     4  +coeff( 22)*x12                
     5  +coeff( 23)*x11*x22            
     6  +coeff( 24)*x11    *x31*x41    
     7  +coeff( 25)    *x22*x31*x41    
     8  +coeff( 26)*x11        *x42    
      t_h_q3ex    =t_h_q3ex    
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)    *x21*x32    *x51
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)*x11            *x52
     4  +coeff( 31)        *x31*x41*x52
     5  +coeff( 32)*x11*x23            
     6  +coeff( 33)*x11*x21*x31*x41    
     7  +coeff( 34)*x11*x21    *x42    
     8  +coeff( 35)*x12            *x51
      t_h_q3ex    =t_h_q3ex    
     9  +coeff( 36)*x12*x22            
     1  +coeff( 37)*x11    *x34        
     2  +coeff( 38)*x12        *x42    
     3  +coeff( 39)    *x24            
     4  +coeff( 40)*x11*x21        *x51
     5  +coeff( 41)    *x21*x31*x41*x51
     6  +coeff( 42)            *x42*x52
     7  +coeff( 43)*x12*x21            
     8  +coeff( 44)*x11*x22        *x51
      t_h_q3ex    =t_h_q3ex    
     9  +coeff( 45)*x11    *x31*x41*x51
     1  +coeff( 46)    *x22    *x42*x51
     2  +coeff( 47)*x11    *x32        
     3  +coeff( 48)    *x22*x32        
     4  +coeff( 49)    *x23        *x51
     5  +coeff( 50)        *x32    *x52
c
      return
      end
      function y_h_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49916E-01,-0.58195E-01,-0.49967E-01,-0.38531E-01,-0.49743E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49756E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.75681210E-01, 0.10657150E+00,-0.22240722E-01,-0.77030550E-01,
     +  0.22430712E-01, 0.75910170E-01, 0.11034391E-01,-0.10709131E-01,
     + -0.85496070E-03, 0.33496331E-01,-0.39938550E-01,-0.56499140E-02,
     + -0.12406274E-01,-0.11045943E-01,-0.21724950E-02,-0.13665831E-01,
     +  0.33233631E-03,-0.30131554E-02, 0.50723054E-02, 0.13440820E-02,
     +  0.56121504E-03, 0.29378034E-01, 0.74553360E-02, 0.47593950E-02,
     + -0.12677330E-02,-0.27212100E-02,-0.46217610E-04, 0.88208550E-03,
     + -0.59218620E-02, 0.53458620E-03, 0.18945274E-02,-0.54950494E-03,
     + -0.15412670E-02,-0.18020684E-03,-0.81687590E-02,-0.77346754E-02,
     + -0.13857380E-01, 0.12011213E-02,-0.92394422E-03, 0.44515100E-02,
     + -0.56048482E-03, 0.88020954E-02, 0.33490060E-02, 0.12367663E-02,
     +  0.13868804E-02,-0.17419140E-02, 0.74343563E-03,-0.76939640E-03,
     + -0.32224750E-03,-0.15456860E-02,-0.25849684E-02,-0.49354580E-02,
     +  0.45309440E-02,-0.32787274E-02,-0.21461750E-02, 0.12526890E-02,
     +  0.44156910E-03, 0.36502280E-02, 0.21306774E-02,-0.18715903E-02,
     + -0.21365291E-02, 0.16697732E-02, 0.51073110E-02, 0.44088700E-03,
     + -0.61198930E-03,-0.64094113E-02,-0.13107290E-01, 0.20124660E-02,
     +  0.13842991E-02,-0.33315242E-03, 0.29415860E-02, 0.78758860E-02,
     + -0.10816360E-02,-0.21272480E-02, 0.34478350E-02,-0.39066094E-02,
     + -0.60984710E-02, 0.62294481E-02, 0.61708150E-03,-0.58081970E-04,
     + -0.19363470E-02, 0.43334370E-03,-0.43461920E-03,-0.74869780E-03,
     +  0.23726891E-03, 0.32948583E-03,-0.23127480E-03, 0.15451860E-03,
     +  0.14801084E-02, 0.19137264E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_h_q3ex    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      y_h_q3ex    =y_h_q3ex    
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)        *x31    *x52
      y_h_q3ex    =y_h_q3ex    
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)*x11*x21*x31        
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)    *x21*x33        
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)    *x22*x31    *x51
      y_h_q3ex    =y_h_q3ex    
     9  +coeff( 27)        *x33    *x51
     1  +coeff( 28)*x11        *x41*x51
     2  +coeff( 29)    *x22    *x41*x51
     3  +coeff( 30)        *x32*x41*x51
     4  +coeff( 31)        *x31*x42*x51
     5  +coeff( 32)        *x31    *x53
     6  +coeff( 33)*x12    *x31        
     7  +coeff( 34)*x11*x22*x31        
     8  +coeff( 35)*x12        *x41    
      y_h_q3ex    =y_h_q3ex    
     9  +coeff( 36)*x11*x22    *x41    
     1  +coeff( 37)    *x22    *x43    
     2  +coeff( 38)*x11*x21*x31    *x51
     3  +coeff( 39)*x11*x23*x31        
     4  +coeff( 40)*x12*x21    *x41    
     5  +coeff( 41)*x11*x23    *x41    
     6  +coeff( 42)*x11*x21    *x43    
     7  +coeff( 43)    *x21*x32*x43    
     8  +coeff( 44)    *x21*x31*x44    
      y_h_q3ex    =y_h_q3ex    
     9  +coeff( 45)*x12        *x41*x51
     1  +coeff( 46)*x11*x22    *x41*x51
     2  +coeff( 47)*x11    *x32*x41*x51
     3  +coeff( 48)*x11*x21    *x41*x52
     4  +coeff( 49)    *x22    *x41*x53
     5  +coeff( 50)*x12*x22*x31        
     6  +coeff( 51)*x12*x22    *x41    
     7  +coeff( 52)*x12    *x31*x42    
     8  +coeff( 53)*x11*x22*x31*x42    
      y_h_q3ex    =y_h_q3ex    
     9  +coeff( 54)*x12        *x43    
     1  +coeff( 55)    *x22*x31*x44    
     2  +coeff( 56)    *x23    *x43*x51
     3  +coeff( 57)    *x21*x31*x44*x51
     4  +coeff( 58)*x11*x23*x31*x42    
     5  +coeff( 59)*x11*x23    *x43    
     6  +coeff( 60)    *x21*x34*x43    
     7  +coeff( 61)*x11*x22    *x43*x51
     8  +coeff( 62)    *x21*x32*x41    
      y_h_q3ex    =y_h_q3ex    
     9  +coeff( 63)    *x21*x31*x42    
     1  +coeff( 64)            *x43*x51
     2  +coeff( 65)            *x41*x53
     3  +coeff( 66)    *x22*x32*x41    
     4  +coeff( 67)    *x22*x31*x42    
     5  +coeff( 68)    *x23    *x41*x51
     6  +coeff( 69)    *x22    *x41*x52
     7  +coeff( 70)*x11*x21*x33        
     8  +coeff( 71)*x11*x21*x32*x41    
      y_h_q3ex    =y_h_q3ex    
     9  +coeff( 72)*x11*x21*x31*x42    
     1  +coeff( 73)    *x24*x33        
     2  +coeff( 74)*x12    *x32*x41    
     3  +coeff( 75)    *x24    *x43    
     4  +coeff( 76)*x12*x24*x31        
     5  +coeff( 77)*x12*x24    *x41    
     6  +coeff( 78)*x11*x24    *x43    
     7  +coeff( 79)    *x21    *x41*x52
     8  +coeff( 80)*x11    *x33        
      y_h_q3ex    =y_h_q3ex    
     9  +coeff( 81)    *x24    *x41    
     1  +coeff( 82)*x11    *x32*x41    
     2  +coeff( 83)*x11*x21    *x41*x51
     3  +coeff( 84)    *x21*x32*x41*x51
     4  +coeff( 85)*x11    *x31    *x52
     5  +coeff( 86)    *x22*x31    *x52
     6  +coeff( 87)*x11        *x41*x52
     7  +coeff( 88)*x12*x21*x31        
     8  +coeff( 89)    *x23*x32*x41    
      y_h_q3ex    =y_h_q3ex    
     9  +coeff( 90)    *x23*x31*x42    
c
      return
      end
      function p_h_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49916E-01,-0.58195E-01,-0.49967E-01,-0.38531E-01,-0.49743E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49756E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14936370E-01,-0.46864733E-01, 0.47950350E-02, 0.21509730E-01,
     + -0.92439160E-02,-0.16271980E-01,-0.21930410E-02, 0.37552001E-02,
     + -0.10050082E-01, 0.10082950E-01, 0.20532694E-02, 0.47332630E-02,
     +  0.33862213E-02, 0.63777441E-03, 0.77115593E-03, 0.17383930E-02,
     +  0.20994194E-02,-0.23388052E-02,-0.45452473E-03,-0.20862860E-04,
     + -0.59344852E-02,-0.10176131E-02,-0.18157960E-02,-0.91152620E-03,
     + -0.32465794E-03,-0.36283694E-04, 0.17671894E-02,-0.62300880E-03,
     + -0.77135802E-03, 0.29464901E-03, 0.32284262E-03,-0.51993253E-03,
     + -0.63222902E-03, 0.95491902E-03, 0.30479930E-03, 0.19569951E-02,
     +  0.39472174E-03, 0.18872321E-02, 0.39103240E-03,-0.43351613E-03,
     +  0.20023183E-02, 0.61856410E-03, 0.46814130E-03,-0.38682282E-03,
     + -0.50739630E-03, 0.41080290E-03, 0.10954382E-02,-0.32291901E-03,
     +  0.12491972E-02,-0.71833360E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_h_q3ex    =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      p_h_q3ex    =p_h_q3ex    
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)        *x32*x41    
     3  +coeff( 12)        *x31*x42    
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x21*x31    *x51
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)            *x41*x52
      p_h_q3ex    =p_h_q3ex    
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)        *x33    *x51
      p_h_q3ex    =p_h_q3ex    
     9  +coeff( 27)*x11        *x41*x51
     1  +coeff( 28)        *x31*x42*x51
     2  +coeff( 29)            *x43*x51
     3  +coeff( 30)    *x21*x31    *x52
     4  +coeff( 31)    *x21    *x41*x52
     5  +coeff( 32)        *x31    *x53
     6  +coeff( 33)            *x41*x53
     7  +coeff( 34)*x12    *x31        
     8  +coeff( 35)    *x22*x33        
      p_h_q3ex    =p_h_q3ex    
     9  +coeff( 36)*x12        *x41    
     1  +coeff( 37)*x11*x22    *x41    
     2  +coeff( 38)    *x22*x31*x42    
     3  +coeff( 39)        *x33*x42    
     4  +coeff( 40)*x11        *x43    
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)*x11*x21    *x41*x51
     7  +coeff( 43)        *x33    *x52
     8  +coeff( 44)*x11        *x41*x52
      p_h_q3ex    =p_h_q3ex    
     9  +coeff( 45)    *x21    *x41*x53
     1  +coeff( 46)        *x31    *x54
     2  +coeff( 47)*x11*x23*x31        
     3  +coeff( 48)*x12*x21    *x41    
     4  +coeff( 49)*x11*x23    *x41    
     5  +coeff( 50)    *x23*x32*x41    
c
      return
      end
      function l_h_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.2672948E-01/
      data xmin/
     1 -0.49916E-01,-0.58195E-01,-0.49917E-01,-0.37177E-01,-0.49743E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49756E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22107854E-01,-0.30847504E+00, 0.15290800E+00,-0.31179944E-01,
     + -0.31978590E-01, 0.27905374E-01,-0.14168570E-01,-0.13572700E-01,
     + -0.77936481E-02,-0.66827070E-02,-0.94821294E-02, 0.17930630E-01,
     +  0.77609630E-02,-0.21823553E-02, 0.21199763E-02,-0.10485412E-02,
     +  0.41734660E-02, 0.35401273E-02,-0.24117173E-02,-0.48804520E-02,
     +  0.30531974E-02, 0.11623742E-01,-0.25939930E-02, 0.21840800E-02,
     +  0.39946171E-02, 0.33901224E-02, 0.80294843E-03, 0.11819710E-02,
     +  0.12330300E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_h_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)*x11*x21            
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)    *x21        *x51
      l_h_q3ex    =l_h_q3ex    
     9  +coeff(  9)*x12                
     1  +coeff( 10)                *x52
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)*x11*x21        *x51
     6  +coeff( 15)    *x23*x31*x41    
     7  +coeff( 16)*x12    *x31*x41    
     8  +coeff( 17)*x12*x21*x31*x41    
      l_h_q3ex    =l_h_q3ex    
     9  +coeff( 18)*x12            *x53
     1  +coeff( 19)        *x32        
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)    *x21*x32        
     4  +coeff( 22)    *x21*x31*x41    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)        *x31*x41*x51
     7  +coeff( 25)*x11    *x31*x41    
     8  +coeff( 26)*x11        *x42    
      l_h_q3ex    =l_h_q3ex    
     9  +coeff( 27)        *x32    *x51
     1  +coeff( 28)*x11    *x32        
     2  +coeff( 29)*x11            *x52
c
      return
      end
      function x_h_fp      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 99)
      data ncoeff/ 98/
      data avdat/ -0.2421117E-01/
      data xmin/
     1 -0.49916E-01,-0.58195E-01,-0.49967E-01,-0.38531E-01,-0.49743E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49756E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.24286630E-01,-0.11277183E-01,-0.10859660E+00, 0.59249620E+00,
     +  0.12496204E-01,-0.14831861E-02,-0.14665860E-01,-0.13090062E-01,
     +  0.34681590E-02, 0.45957520E-01,-0.54860440E-02,-0.86097093E-02,
     + -0.46639030E-01, 0.31433020E-02,-0.59460750E-02,-0.31191650E-01,
     + -0.58705173E-02,-0.16037671E-02, 0.11341150E-01, 0.45085350E-02,
     + -0.15722800E-02, 0.38508982E-02, 0.25951260E-03,-0.32607380E-02,
     + -0.95008080E-03, 0.21960210E-02, 0.17570650E-02, 0.44938200E-02,
     +  0.15387231E-01, 0.14245940E-01, 0.15404970E-02,-0.19113622E-01,
     + -0.10600310E-01, 0.41892640E-02,-0.26794180E-03, 0.76001830E-04,
     + -0.36377783E-02,-0.40231761E-02, 0.13905220E-02, 0.23079680E-02,
     + -0.71062800E-02,-0.32591684E-02,-0.46009410E-02,-0.20842393E-02,
     + -0.22494360E-02, 0.32760151E-02, 0.10416992E-01, 0.90548820E-02,
     + -0.51135080E-03, 0.70694403E-03,-0.18311900E-02, 0.13294960E-02,
     +  0.60759242E-02, 0.15932000E-02, 0.97269820E-03,-0.36732791E-02,
     +  0.31879700E-02,-0.15310190E-02,-0.59692293E-03,-0.54608400E-02,
     + -0.70034070E-02,-0.65153581E-02, 0.44548910E-02,-0.62765060E-03,
     + -0.47643654E-03,-0.17954690E-02,-0.75981393E-03, 0.61074470E-03,
     +  0.37425310E-02, 0.72558142E-03,-0.54057384E-03,-0.30779060E-02,
     + -0.23942730E-02,-0.49346610E-02,-0.72191434E-03,-0.71411422E-03,
     +  0.38361761E-02,-0.16903660E-02,-0.67832092E-02,-0.49816244E-02,
     + -0.36847890E-02,-0.41711190E-02,-0.11455430E-02, 0.22485114E-02,
     +  0.66712530E-03,-0.10404200E-02,-0.14746483E-02, 0.85169180E-02,
     + -0.70671700E-02, 0.31207110E-02, 0.36115380E-02, 0.14189422E-02,
     +  0.23267643E-02,-0.20222320E-02,-0.38290620E-02,-0.68946190E-02,
     + -0.54307410E-02, 0.48593930E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x14 = x13*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_h_fp      =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      x_h_fp      =x_h_fp      
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)*x12                
     3  +coeff( 12)*x11            *x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x22            
      x_h_fp      =x_h_fp      
     9  +coeff( 18)*x11    *x31*x41    
     1  +coeff( 19)*x11        *x42    
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)        *x32    *x51
     4  +coeff( 22)        *x31*x41*x51
     5  +coeff( 23)            *x42*x51
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)    *x21        *x52
     8  +coeff( 26)*x12            *x51
      x_h_fp      =x_h_fp      
     9  +coeff( 27)*x11            *x52
     1  +coeff( 28)                *x53
     2  +coeff( 29)    *x22    *x42    
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)*x11*x21*x31*x41    
     5  +coeff( 32)*x11*x21    *x42    
     6  +coeff( 33)*x12*x22            
     7  +coeff( 34)*x12        *x42    
     8  +coeff( 35)    *x23*x32        
      x_h_fp      =x_h_fp      
     9  +coeff( 36)    *x21*x34        
     1  +coeff( 37)    *x24    *x42    
     2  +coeff( 38)    *x21*x32        
     3  +coeff( 39)*x11    *x32        
     4  +coeff( 40)*x12*x21            
     5  +coeff( 41)    *x24            
     6  +coeff( 42)    *x21*x31*x41*x51
     7  +coeff( 43)    *x21    *x42*x51
     8  +coeff( 44)*x11    *x32    *x51
      x_h_fp      =x_h_fp      
     9  +coeff( 45)            *x42*x52
     1  +coeff( 46)*x13*x21            
     2  +coeff( 47)    *x22    *x42*x51
     3  +coeff( 48)*x11*x23        *x51
     4  +coeff( 49)*x13                
     5  +coeff( 50)    *x22*x32        
     6  +coeff( 51)    *x22*x31*x41    
     7  +coeff( 52)    *x21*x32    *x51
     8  +coeff( 53)*x11*x22        *x51
      x_h_fp      =x_h_fp      
     9  +coeff( 54)*x11    *x31*x41*x51
     1  +coeff( 55)        *x32    *x52
     2  +coeff( 56)*x12*x21        *x51
     3  +coeff( 57)*x11*x21        *x52
     4  +coeff( 58)    *x25            
     5  +coeff( 59)*x11            *x53
     6  +coeff( 60)    *x24        *x51
     7  +coeff( 61)*x11*x21    *x42*x51
     8  +coeff( 62)*x12*x22        *x51
      x_h_fp      =x_h_fp      
     9  +coeff( 63)*x11*x22        *x52
     1  +coeff( 64)*x11*x23*x32        
     2  +coeff( 65)    *x25        *x51
     3  +coeff( 66)*x12    *x33*x41    
     4  +coeff( 67)    *x24        *x52
     5  +coeff( 68)    *x25        *x52
     6  +coeff( 69)    *x23    *x42*x52
     7  +coeff( 70)        *x33*x41    
     8  +coeff( 71)*x11*x21*x32        
      x_h_fp      =x_h_fp      
     9  +coeff( 72)    *x23        *x51
     1  +coeff( 73)    *x22        *x52
     2  +coeff( 74)    *x23*x31*x41    
     3  +coeff( 75)    *x23    *x42    
     4  +coeff( 76)*x12            *x52
     5  +coeff( 77)*x11*x24            
     6  +coeff( 78)*x11*x22*x32        
     7  +coeff( 79)*x11*x22    *x42    
     8  +coeff( 80)*x12*x23            
      x_h_fp      =x_h_fp      
     9  +coeff( 81)*x12*x21*x31*x41    
     1  +coeff( 82)    *x23        *x52
     2  +coeff( 83)    *x21*x32    *x52
     3  +coeff( 84)*x13*x22            
     4  +coeff( 85)*x12        *x42*x51
     5  +coeff( 86)            *x42*x53
     6  +coeff( 87)*x12*x21        *x52
     7  +coeff( 88)*x11*x23    *x42    
     8  +coeff( 89)*x12*x22    *x42    
      x_h_fp      =x_h_fp      
     9  +coeff( 90)    *x22    *x42*x52
     1  +coeff( 91)*x13*x21    *x42    
     2  +coeff( 92)*x13    *x32    *x51
     3  +coeff( 93)*x13*x23        *x51
     4  +coeff( 94)    *x25*x31*x41*x51
     5  +coeff( 95)*x13*x21    *x42*x52
     6  +coeff( 96)*x14*x21*x32*x42    
     7  +coeff( 97)*x14*x21*x31*x43    
     8  +coeff( 98)        *x32*x42    
c
      return
      end
      function t_h_fp      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/ -0.5844322E-02/
      data xmin/
     1 -0.49916E-01,-0.58195E-01,-0.49967E-01,-0.38531E-01,-0.49743E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49756E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.54843751E-02,-0.27646191E-01, 0.97765021E-01,-0.51363823E-02,
     +  0.61653430E-04,-0.19645714E-03,-0.24542822E-02,-0.17126270E-02,
     +  0.42407730E-02,-0.90651450E-02, 0.29384554E-02, 0.18060430E-03,
     + -0.35691890E-03,-0.29223934E-04,-0.43115060E-02,-0.14355490E-02,
     + -0.27786981E-03, 0.95334474E-03, 0.41842952E-03,-0.31061280E-03,
     +  0.57739671E-03,-0.14564020E-02,-0.97997610E-04,-0.16124702E-03,
     + -0.91851510E-04, 0.21728580E-02, 0.25497530E-02, 0.16333610E-03,
     + -0.89549330E-03, 0.55646970E-03, 0.15125292E-03, 0.20847920E-02,
     +  0.30517260E-03,-0.28069843E-02, 0.29385284E-03,-0.14666090E-02,
     + -0.13779630E-04, 0.57385710E-03,-0.99419630E-03,-0.35386310E-03,
     + -0.82771420E-03,-0.39031500E-03, 0.14334842E-03, 0.77096372E-03,
     +  0.31635013E-03, 0.16754570E-02, 0.32423634E-03, 0.18343710E-03,
     + -0.49363240E-03, 0.23073650E-03, 0.54446310E-03,-0.47648110E-03,
     + -0.22594150E-03,-0.23816962E-03,-0.18600502E-03,-0.54669050E-03,
     + -0.56294713E-04,-0.26814520E-03, 0.45920930E-03,-0.45518612E-03,
     +  0.71013240E-03,-0.23603964E-03,-0.11879680E-02,-0.29523460E-03,
     + -0.37242731E-03,-0.18197672E-03, 0.25505572E-03,-0.86351741E-04,
     + -0.21774700E-03, 0.11809463E-03,-0.12353701E-03, 0.45827444E-03,
     +  0.51171040E-03, 0.12548650E-03,-0.11003960E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_h_fp      =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_h_fp      =t_h_fp      
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)        *x32    *x51
      t_h_fp      =t_h_fp      
     9  +coeff( 18)        *x31*x41*x51
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)                *x53
     4  +coeff( 22)*x12                
     5  +coeff( 23)*x11*x22            
     6  +coeff( 24)*x11    *x31*x41    
     7  +coeff( 25)    *x22*x31*x41    
     8  +coeff( 26)*x11        *x42    
      t_h_fp      =t_h_fp      
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)    *x21*x32    *x51
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)*x11            *x52
     4  +coeff( 31)        *x31*x41*x52
     5  +coeff( 32)*x11*x23            
     6  +coeff( 33)*x11*x21*x31*x41    
     7  +coeff( 34)*x11*x21    *x42    
     8  +coeff( 35)*x12            *x51
      t_h_fp      =t_h_fp      
     9  +coeff( 36)*x12*x22            
     1  +coeff( 37)*x11    *x34        
     2  +coeff( 38)*x12        *x42    
     3  +coeff( 39)    *x24            
     4  +coeff( 40)*x11*x21        *x51
     5  +coeff( 41)    *x21*x31*x41*x51
     6  +coeff( 42)            *x42*x52
     7  +coeff( 43)*x12*x21            
     8  +coeff( 44)*x11*x22        *x51
      t_h_fp      =t_h_fp      
     9  +coeff( 45)*x11    *x31*x41*x51
     1  +coeff( 46)    *x22    *x42*x51
     2  +coeff( 47)*x11    *x32        
     3  +coeff( 48)    *x22*x32        
     4  +coeff( 49)    *x23        *x51
     5  +coeff( 50)        *x32    *x52
     6  +coeff( 51)    *x23    *x42    
     7  +coeff( 52)    *x24        *x51
     8  +coeff( 53)*x11    *x32    *x51
      t_h_fp      =t_h_fp      
     9  +coeff( 54)    *x21*x32    *x52
     1  +coeff( 55)*x11            *x53
     2  +coeff( 56)            *x42*x53
     3  +coeff( 57)*x11*x22*x32        
     4  +coeff( 58)*x12    *x31*x41    
     5  +coeff( 59)*x11*x22*x31*x41    
     6  +coeff( 60)*x12*x21        *x51
     7  +coeff( 61)*x11*x23        *x51
     8  +coeff( 62)*x11*x21*x31*x41*x51
      t_h_fp      =t_h_fp      
     9  +coeff( 63)*x11*x21    *x42*x51
     1  +coeff( 64)    *x24        *x52
     2  +coeff( 65)*x11    *x31*x41*x52
     3  +coeff( 66)*x11*x23*x32        
     4  +coeff( 67)*x12*x21        *x53
     5  +coeff( 68)*x11*x21*x32        
     6  +coeff( 69)        *x32*x42*x51
     7  +coeff( 70)*x11*x21        *x52
     8  +coeff( 71)    *x23        *x52
      t_h_fp      =t_h_fp      
     9  +coeff( 72)    *x21*x31*x41*x52
     1  +coeff( 73)    *x21    *x42*x52
     2  +coeff( 74)    *x22        *x53
     3  +coeff( 75)        *x32    *x53
c
      return
      end
      function y_h_fp      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 98)
      data ncoeff/ 97/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49916E-01,-0.58195E-01,-0.49967E-01,-0.38531E-01,-0.49743E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49756E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.32664854E-01,-0.27110470E-01,-0.80249812E-02,-0.13638180E-01,
     + -0.36577780E-02, 0.30838350E-01, 0.45124143E-02,-0.97184232E-03,
     +  0.10133010E-03, 0.37141940E-02,-0.16654760E-01,-0.17257540E-02,
     + -0.13546671E-02,-0.31549200E-02,-0.10392862E-01, 0.45335381E-02,
     +  0.30435093E-02,-0.21609204E-03,-0.32155810E-03, 0.44277694E-03,
     +  0.15457490E-01, 0.30886460E-02, 0.98068201E-04,-0.68175910E-03,
     +  0.91300610E-03,-0.12765450E-02,-0.30583320E-02, 0.53497264E-02,
     + -0.69239140E-02,-0.11144150E-02, 0.22052592E-04,-0.34343281E-02,
     + -0.14479930E-03, 0.29090282E-03,-0.25475313E-02,-0.34920503E-02,
     +  0.40173104E-05, 0.95370210E-03,-0.57150620E-04,-0.47486040E-02,
     + -0.48595820E-02, 0.64424780E-03, 0.64555840E-03, 0.11216262E-02,
     +  0.20805874E-02,-0.14803220E-02, 0.16992632E-02,-0.12927451E-02,
     +  0.25328111E-02, 0.80010312E-03, 0.35177790E-03, 0.32413094E-02,
     +  0.27949463E-02,-0.29035610E-02, 0.29208260E-02,-0.44625400E-03,
     + -0.53526390E-04, 0.30632710E-02, 0.79099690E-03,-0.11023111E-02,
     +  0.27634320E-03,-0.20515700E-02, 0.92544310E-03,-0.11389500E-02,
     + -0.20015270E-02,-0.19284032E-02,-0.14276450E-02,-0.17919500E-02,
     + -0.52070990E-04, 0.41026830E-02,-0.23301200E-02,-0.62318984E-03,
     +  0.17075300E-02, 0.58927682E-04,-0.81121020E-03,-0.80683652E-03,
     +  0.29966023E-02,-0.24646550E-02, 0.13643281E-02, 0.12631430E-02,
     + -0.65830163E-03, 0.20348210E-02, 0.19196280E-02, 0.16181800E-02,
     +  0.56382400E-03,-0.15025583E-02,-0.21865332E-02,-0.26810483E-03,
     + -0.13724993E-02,-0.60403250E-03,-0.19757540E-02, 0.12929380E-02,
     + -0.42735203E-02, 0.15066593E-02,-0.75198010E-02, 0.27450644E-02,
     +  0.29428450E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_h_fp      =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      y_h_fp      =y_h_fp      
     9  +coeff(  9)        *x33        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)            *x41*x52
      y_h_fp      =y_h_fp      
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x32*x41    
     6  +coeff( 24)    *x21*x31*x42    
     7  +coeff( 25)    *x21    *x43    
     8  +coeff( 26)*x11    *x31    *x51
      y_h_fp      =y_h_fp      
     9  +coeff( 27)    *x22*x31    *x51
     1  +coeff( 28)*x11        *x41*x51
     2  +coeff( 29)    *x22    *x41*x51
     3  +coeff( 30)        *x32*x41*x51
     4  +coeff( 31)        *x31*x42*x51
     5  +coeff( 32)            *x43*x51
     6  +coeff( 33)    *x21*x31    *x52
     7  +coeff( 34)    *x21    *x41*x52
     8  +coeff( 35)        *x31    *x53
      y_h_fp      =y_h_fp      
     9  +coeff( 36)            *x41*x53
     1  +coeff( 37)*x12    *x31        
     2  +coeff( 38)*x11*x22*x31        
     3  +coeff( 39)*x11    *x33        
     4  +coeff( 40)*x12        *x41    
     5  +coeff( 41)*x11*x22    *x41    
     6  +coeff( 42)    *x24    *x41    
     7  +coeff( 43)*x11    *x32*x41    
     8  +coeff( 44)*x11    *x31*x42    
      y_h_fp      =y_h_fp      
     9  +coeff( 45)*x11*x21*x31    *x51
     1  +coeff( 46)    *x23*x31    *x51
     2  +coeff( 47)*x11*x21    *x41*x51
     3  +coeff( 48)    *x21*x32*x41*x51
     4  +coeff( 49)    *x21*x31*x42*x51
     5  +coeff( 50)*x11    *x31    *x52
     6  +coeff( 51)*x11        *x41*x52
     7  +coeff( 52)    *x22    *x41*x52
     8  +coeff( 53)        *x31*x42*x52
      y_h_fp      =y_h_fp      
     9  +coeff( 54)    *x21    *x41*x53
     1  +coeff( 55)        *x31    *x54
     2  +coeff( 56)*x12*x21*x31        
     3  +coeff( 57)*x11*x21*x33        
     4  +coeff( 58)*x12*x21    *x41    
     5  +coeff( 59)*x11*x21    *x43    
     6  +coeff( 60)*x12    *x31    *x51
     7  +coeff( 61)    *x22*x33    *x51
     8  +coeff( 62)*x11*x22    *x41*x51
      y_h_fp      =y_h_fp      
     9  +coeff( 63)*x11    *x32*x41*x51
     1  +coeff( 64)    *x22*x32*x41*x51
     2  +coeff( 65)*x11    *x31*x42*x51
     3  +coeff( 66)    *x22*x31*x42*x51
     4  +coeff( 67)        *x33*x42*x51
     5  +coeff( 68)*x11        *x43*x51
     6  +coeff( 69)*x11*x21    *x41*x52
     7  +coeff( 70)    *x21*x31*x42*x52
     8  +coeff( 71)*x11    *x31    *x53
      y_h_fp      =y_h_fp      
     9  +coeff( 72)        *x33    *x53
     1  +coeff( 73)*x11        *x41*x53
     2  +coeff( 74)*x12*x21*x31    *x51
     3  +coeff( 75)    *x23*x32*x41*x51
     4  +coeff( 76)    *x23*x31*x42*x51
     5  +coeff( 77)    *x23    *x43*x51
     6  +coeff( 78)    *x21*x31*x44*x51
     7  +coeff( 79)*x12    *x31    *x52
     8  +coeff( 80)    *x22*x33    *x52
      y_h_fp      =y_h_fp      
     9  +coeff( 81)*x11    *x32*x41*x52
     1  +coeff( 82)    *x22*x31*x42*x52
     2  +coeff( 83)    *x23*x31    *x53
     3  +coeff( 84)    *x23    *x41*x53
     4  +coeff( 85)*x11    *x31    *x54
     5  +coeff( 86)    *x22    *x41*x54
     6  +coeff( 87)            *x43*x54
     7  +coeff( 88)    *x24*x33    *x51
     8  +coeff( 89)*x11*x23*x31    *x52
      y_h_fp      =y_h_fp      
     9  +coeff( 90)*x11*x21*x33    *x52
     1  +coeff( 91)*x11*x23    *x41*x52
     2  +coeff( 92)    *x23*x32*x41*x52
     3  +coeff( 93)*x11*x21*x31*x42*x52
     4  +coeff( 94)    *x21*x32*x43*x52
     5  +coeff( 95)    *x21*x31*x44*x52
     6  +coeff( 96)    *x22*x31*x42*x53
     7  +coeff( 97)*x11        *x43*x53
c
      return
      end
      function p_h_fp      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 81)
      data ncoeff/ 80/
      data avdat/  0.0000000E+00/
      data xmin/
     1 -0.49916E-01,-0.58195E-01,-0.49967E-01,-0.38531E-01,-0.49743E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49756E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14913730E-01,-0.46598341E-01, 0.48219110E-02, 0.22023600E-01,
     + -0.90491883E-02,-0.15737414E-01,-0.22045160E-02, 0.28688720E-02,
     + -0.10218570E-01, 0.80243991E-02, 0.16763001E-02, 0.41528040E-02,
     +  0.31400264E-02, 0.68979820E-03, 0.88253541E-03, 0.17601693E-02,
     +  0.23403551E-02,-0.17674210E-02,-0.42190020E-03,-0.53926930E-04,
     + -0.46801210E-02,-0.15127932E-02,-0.22280880E-02,-0.15084081E-02,
     + -0.17637011E-03,-0.16993590E-03, 0.13376441E-02,-0.10028590E-02,
     + -0.11599902E-02, 0.19823660E-03,-0.62130530E-04,-0.66171912E-03,
     + -0.90016383E-03, 0.70608490E-03, 0.28582022E-03, 0.15560820E-02,
     +  0.72096654E-03, 0.49933830E-02, 0.82143851E-04,-0.28301772E-03,
     +  0.39931721E-02, 0.71987794E-03, 0.15204402E-03, 0.71109680E-04,
     + -0.86504450E-03, 0.62545240E-03, 0.89502012E-04,-0.40322100E-03,
     + -0.33327680E-03,-0.75740812E-04, 0.99999230E-04,-0.35557453E-02,
     + -0.30864262E-02,-0.19182800E-03, 0.30848332E-04,-0.52606540E-03,
     + -0.35745100E-03, 0.28057070E-03,-0.40272710E-03, 0.27669351E-04,
     + -0.33508880E-03, 0.86389100E-03, 0.12938490E-02,-0.71639122E-04,
     +  0.21983201E-02, 0.19149984E-03, 0.62488892E-03, 0.22901363E-04,
     +  0.55461863E-04, 0.29514724E-03,-0.81318280E-03,-0.57721850E-03,
     +  0.79645500E-03, 0.14032900E-02, 0.18697611E-02, 0.22159700E-03,
     +  0.13120240E-03, 0.63010820E-03,-0.43396240E-03,-0.10622890E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_h_fp      =avdat
     1  +coeff(  1)        *x31        
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x22*x31        
      p_h_fp      =p_h_fp      
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)        *x32*x41    
     3  +coeff( 12)        *x31*x42    
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x21*x31    *x51
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)            *x41*x52
      p_h_fp      =p_h_fp      
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21*x33        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)*x11    *x31    *x51
     8  +coeff( 26)        *x33    *x51
      p_h_fp      =p_h_fp      
     9  +coeff( 27)*x11        *x41*x51
     1  +coeff( 28)        *x31*x42*x51
     2  +coeff( 29)            *x43*x51
     3  +coeff( 30)    *x21*x31    *x52
     4  +coeff( 31)    *x21    *x41*x52
     5  +coeff( 32)        *x31    *x53
     6  +coeff( 33)            *x41*x53
     7  +coeff( 34)*x12    *x31        
     8  +coeff( 35)    *x22*x33        
      p_h_fp      =p_h_fp      
     9  +coeff( 36)*x12        *x41    
     1  +coeff( 37)*x11*x22    *x41    
     2  +coeff( 38)    *x22*x31*x42    
     3  +coeff( 39)        *x33*x42    
     4  +coeff( 40)*x11        *x43    
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)*x11*x21    *x41*x51
     7  +coeff( 43)        *x33    *x52
     8  +coeff( 44)*x11        *x41*x52
      p_h_fp      =p_h_fp      
     9  +coeff( 45)    *x21    *x41*x53
     1  +coeff( 46)        *x31    *x54
     2  +coeff( 47)*x11*x23*x31        
     3  +coeff( 48)*x12*x21    *x41    
     4  +coeff( 49)*x11*x23    *x41    
     5  +coeff( 50)    *x23*x32*x41    
     6  +coeff( 51)    *x21*x34*x41    
     7  +coeff( 52)*x11*x21*x31*x42    
     8  +coeff( 53)*x11*x21    *x43    
      p_h_fp      =p_h_fp      
     9  +coeff( 54)*x12    *x31    *x51
     1  +coeff( 55)    *x22*x33    *x51
     2  +coeff( 56)*x12        *x41*x51
     3  +coeff( 57)    *x24    *x41*x51
     4  +coeff( 58)*x11    *x32*x41*x51
     5  +coeff( 59)    *x22*x32*x41*x51
     6  +coeff( 60)        *x34*x41*x51
     7  +coeff( 61)*x11    *x31    *x53
     8  +coeff( 62)*x11        *x41*x53
      p_h_fp      =p_h_fp      
     9  +coeff( 63)*x12*x22*x31        
     1  +coeff( 64)*x12    *x33        
     2  +coeff( 65)*x12*x22    *x41    
     3  +coeff( 66)    *x24*x32*x41    
     4  +coeff( 67)*x12    *x31*x42    
     5  +coeff( 68)*x11    *x33*x42    
     6  +coeff( 69)*x11*x21*x34*x41    
     7  +coeff( 70)        *x33        
     8  +coeff( 71)    *x21*x32*x41    
      p_h_fp      =p_h_fp      
     9  +coeff( 72)        *x32*x41*x51
     1  +coeff( 73)    *x24*x31        
     2  +coeff( 74)    *x24    *x41    
     3  +coeff( 75)    *x22*x32*x41    
     4  +coeff( 76)*x11    *x31*x42    
     5  +coeff( 77)*x11    *x31    *x52
     6  +coeff( 78)        *x31*x42*x52
     7  +coeff( 79)            *x41*x54
     8  +coeff( 80)*x11*x21*x32*x41    
c
      return
      end
      function l_h_fp      (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10) !,xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3203921E-01/
      data xmin/
     1 -0.49916E-01,-0.58195E-01,-0.49917E-01,-0.37177E-01,-0.49743E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49854E-01, 0.59958E-01, 0.49967E-01, 0.38531E-01, 0.49756E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27463080E-01,-0.30836250E+00, 0.15310410E+00,-0.31949450E-01,
     + -0.33311330E-01,-0.20824270E-01, 0.27692070E-01,-0.12713600E-01,
     + -0.78404510E-02,-0.12467630E-01,-0.49907410E-02, 0.20204463E-01,
     + -0.26692180E-02,-0.34709521E-02, 0.12782970E-01, 0.39222310E-02,
     +  0.77088270E-02,-0.31097140E-02, 0.24455094E-02, 0.29524910E-02,
     +  0.29942710E-02, 0.29395252E-02, 0.71436060E-03,-0.44706181E-03,
     +  0.88388730E-03,-0.19145610E-02, 0.12304940E-02,-0.43720854E-04,
     + -0.11005070E-02, 0.33721160E-03,-0.16628884E-02,-0.41280113E-03,
     + -0.13056130E-02, 0.16923742E-02,-0.23426151E-04,-0.16700950E-02,
     + -0.20155801E-02,-0.11102090E-02, 0.94766374E-02, 0.10786190E-01,
     +  0.23618114E-02,-0.81451100E-04, 0.72094211E-04, 0.86190550E-04,
     +  0.80853851E-03,-0.87509641E-03, 0.17573331E-02,-0.82882400E-03,
     + -0.23616370E-03, 0.17779940E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_h_fp      =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)*x11            *x51
      l_h_fp      =l_h_fp      
     9  +coeff(  9)*x12                
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)        *x32        
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)                *x53
     8  +coeff( 17)*x11*x22            
      l_h_fp      =l_h_fp      
     9  +coeff( 18)*x11*x21        *x51
     1  +coeff( 19)*x12            *x51
     2  +coeff( 20)    *x21*x32        
     3  +coeff( 21)        *x31*x41*x51
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)        *x34    *x51
     6  +coeff( 24)            *x41    
     7  +coeff( 25)    *x21    *x41    
     8  +coeff( 26)    *x22        *x51
      l_h_fp      =l_h_fp      
     9  +coeff( 27)*x11    *x32        
     1  +coeff( 28)*x11        *x42    
     2  +coeff( 29)*x12*x21            
     3  +coeff( 30)    *x21*x31        
     4  +coeff( 31)    *x21        *x52
     5  +coeff( 32)    *x22*x32        
     6  +coeff( 33)    *x23        *x51
     7  +coeff( 34)*x11            *x52
     8  +coeff( 35)*x11*x21    *x42    
      l_h_fp      =l_h_fp      
     9  +coeff( 36)*x11    *x31*x41*x51
     1  +coeff( 37)*x11        *x42*x51
     2  +coeff( 38)*x11*x24            
     3  +coeff( 39)*x11*x22*x31*x41    
     4  +coeff( 40)*x11*x22    *x42    
     5  +coeff( 41)*x11*x24*x32        
     6  +coeff( 42)        *x31        
     7  +coeff( 43)        *x31    *x51
     8  +coeff( 44)*x11    *x31        
      l_h_fp      =l_h_fp      
     9  +coeff( 45)        *x32    *x51
     1  +coeff( 46)            *x42*x51
     2  +coeff( 47)    *x22    *x42    
     3  +coeff( 48)        *x32*x42    
     4  +coeff( 49)    *x21*x31*x41*x51
     5  +coeff( 50)    *x21    *x42*x51
c
      return
      end
