      subroutine left_init_r_function
c      
c initializes arrays etc in the common block of the r-function package
c contains some diagnostic write statements you can comment them out if you want
c You must call left_init_r_function once before calling left_rfunction(which you can then
c do over and over again without reinitializing).
c                                           -JJL 6/6/01
      parameter (narray=1400)
      common /lines/rline(12,narray,6),tmin(narray),tmax(narray)
     +              ,phimin(narray),phimax(narray),dy,dd,npts
      integer ndum
      character*24 daet

c
c  for each of the 12 lines defining the theta, phi acceptance of the form:
c   0=A*theta+B*phi+C
c    rline(i,n,1)= yn (y value of nth set) i refers to the line number
c    rline(i,n,2)= deltan (delta value for nth set)
c    rline(i,n,3)=A (for the nth set)
c    rline(i,n,4)=B
c    rline(i,n,5)=C
c    rline(i,n,6)=R^2 or 2 for vertical lines or 3 for lines through the origin
c
c  read in the files with the line parameters
c
      open(11,file='./New-Left/no-collimator/line1.dat',status='old')
      open(12,file='./New-Left/no-collimator/line2.dat',status='old')
      open(13,file='./New-Left/no-collimator/line3.dat',status='old')
      open(14,file='./New-Left/no-collimator/line4.dat',status='old')
      open(15,file='./New-Left/no-collimator/line5.dat',status='old')
      open(16,file='./New-Left/no-collimator/line6.dat',status='old')
      open(17,file='./New-Left/no-collimator/line7.dat',status='old')
      open(18,file='./New-Left/no-collimator/line8.dat',status='old')
      open(19,file='./New-Left/no-collimator/line9.dat',status='old')
      open(20,file='./New-Left/no-collimator/line10.dat',status='old')
      open(21,file='./New-Left/no-collimator/line11.dat',status='old')
      open(22,file='./New-Left/no-collimator/line12.dat',status='old')
      open(23,file='./New-Left/no-collimator/dy-dd.dat',status='old')
      open(30,file='./New-Left/no-collimator/min-max.dat',status='old')
      read(23,*)dy,dd   !y and delta grid spacing
      close(23)
      !write(6,*)"initializing r-function arrays"
      open(31,file='./New-Left/no-collimator/information.dat'
     &,status='old')




      read(31,'(a24)')daet
      !write(6,*)"Files generated on:"
      !write(6,'(a24)')daet
      read(31,'(a24)')daet
      !write(6,'(a24)')daet
      read(31,'(a24)')daet
      !write(6,'(a24)')daet
      close(31)
      do 10 i=1,narray
      
      do 11 n=1,12
      nin=n+10
      read(nin,*,end=12)ndum,(rline(n,i,j),j=1,6)
   11 continue
      read(30,*)ydum,ddum,tmin(i),tmax(i),phimin(i),phimax(i)
   10 continue
   12 npts=i
      if(npts.ge.narray)then
      write(6,*)"ERROR in init-r-function: Need to increase narray"
      endif
      return
      end
      
      subroutine left_getindex(y,d,index)
      
c for given y and delta (y,d) find index for line parameters
      parameter (narray=1400)                  
      common /lines/rline(12,narray,6),tmin(narray),tmax(narray)
     +              ,phimin(narray),phimax(narray),dy,dd,npts
      real y,d,ddummy,testa
      integer index
      ddummy=100000.
      index=0
c      write(6,*)"npts=",npts
      !print*,y,d
      do i=1,npts
       a=((y-rline(1,i,1))**2)+((d-rline(1,i,2))**2)
       if(a.lt.ddummy)then
        index=i
        ddummy=a
       endif
      enddo
c      write(6,*)index
c test to see that a is not too big
      ddummy=sqrt(ddummy)
      testa=sqrt(dd**2+dy**2)/2
c      write(6,*)dd,dy,testa,ddummy
      if(ddummy.gt.testa) index=0   !outside the grid
      return
      end
      
      real function left_rfunction(y,d,theta,phi)
c
c  sample left_rfunction call
c   y is in meters (y0)
c   d (delta) is a fraction
c   theta and phi are radians (theta_0 and phi_0)
c      
c  sample left_rfunction call:  answer=left_rfunction(y,d,theta,phi)
c      
c  answer is the distance in phi/theta space of the given  
c  theta and phi point from the edge of the acceptance for the given y and d
c  negative numbers are outside the acceptance (-1000. is the default outside value)
c      
      
      parameter (narray=1400)                  
      common /lines/rline(12,narray,6),tmin(narray),tmax(narray)
     +              ,phimin(narray),phimax(narray),dy,dd,npts
      real y,d,theta,phi,dline(16)
      integer index
      left_rfunction=100000.
      !print*,"in F,y=",y," d=",d," theta=",theta," phi=",phi
      call left_getindex(y,d,index)
      !write(6,*)"index=",index,y,d 
      if(index.gt.0)then
      do 10 i=1,12   !loop through the 12 lines
c
c set up the signs  (see notes: 5/25-30/01)
c
      if(rline(i,index,4).eq.0.) go to 10  !ignore vertical lines, the box takes care of it      
      if(rline(i,index,6).eq.3.) go to 10  !ignore lines through the origin
      em=-rline(i,index,3)/rline(i,index,4) !conventional slope of the line
      isign=+1
      if((tmax(index).gt.0.).and.(tmin(index).lt.0.))then  !origin is inside the acceptance
       go to 7 !all signs positive
       else
       left_rfunction=-1000.
       return
      endif
c
c  dline is the distance from the point to the line
c       
  7   dline(i)=isign*(
     +         (rline(i,index,3)*theta)
     +        +(rline(i,index,4)*phi)
     +        + rline(i,index,5))
     +        /sqrt(rline(i,index,3)**2+rline(i,index,4)**2)
c
c  then find the minimum
c      
      if(dline(i).lt.left_rfunction)left_rfunction=dline(i)
  10  continue
c
c make the simple box 
c     
      dline(13)=tmax(index)-theta
      if(dline(13).lt.left_rfunction)left_rfunction=dline(13)
      dline(14)=theta-tmin(index)
      if(dline(14).lt.left_rfunction)left_rfunction=dline(14)
      dline(15)=phimax(index)-phi
      if(dline(15).lt.left_rfunction)left_rfunction=dline(15)
      dline(16)=phi-phimin(index)
      if(dline(16).lt.left_rfunction)left_rfunction=dline(16)

      else
      left_rfunction=-1000.  !default no acceptance value
      endif
      return
      end
      
      
