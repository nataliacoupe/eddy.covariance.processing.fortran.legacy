!===================================================================================================
!                                      lba code:  lfluxvar.for          
!===================================================================================================
!     sets of parameter declarations and variables used in eddy flux calculation in lfluxes.for
!---------------------------------------------------------------------------------------------------
!     eddy input vars:     
      real*8 po,pres, to,temp
!     storage variables for one aggregation interval:
      real*8 co2(maxhigh), h2o(maxhigh), u(maxhigh) 
      real*8 h2o_v(maxhigh), co2_v(maxhigh) ! voltages to keep for EB*.FLX output nrc 170125
      real*8 v(maxhigh), w(maxhigh), t(maxhigh)
      real*8 stat(maxhigh), cal(maxhigh)
      real*8 ht(maxhigh), hw(maxhigh)       ! storage for water-smeared t,w.
      real*8 sw(maxhigh)                    ! temp storage for co2-smeared w.
      integer n                             ! number of data points actually input (max n =maxhigh).
!     slow data in 201 file:
!     re-use v() for voltages, store(,) for storing volts
!     flux calculation and output variables 
      logical calbit, unexpcal   ! unexpcal added 30-may-02
      integer calib, endzer
      integer idateMM, idate     ! nrc 181010
      real*8 store(nwdim,10)     ! changed from 7 (5-nov-01)
!     1-4: u,v,w,t; 
!     5-6: co2,h2o; 
!     7  : stat;
!     8  : jday;
!     9  : secs;
!     10 : cal.
      real*8 flux(3,3),avgs(3),vars(3)
!     flux(1,*) = temp flux (smeared)
!     flux(2,*) = water flux h' sm(w')
!     flux(3,*) = co2 flux
      real*8 f9avg(3),f9std(3),f9max(3),f9lag(3),f9bus(4)

!     double-precision for key flux-calc variables:
      real*8 avg,var,aint,slope
      real*8 uint,uslope,vint,vslope,wint,wslope,tint,tslope
      real*8 hint,hslope,cint,cslope
      real*8 xint,yint         		     ! voltages to keep for EB*.FLX output nrc 170125
      real*8 st(3,3),angles(3),cov(3),uvw(3),rot(3,3)  ! , work(13)
      parameter(nlag=6)                      ! number of + and - lags to do correlation at

!     declaration of bad sonic fraction variables
!     -------------------------------------------
!     badsonic--fraction of half hour interval with bad sonic data
!     numbadsonic--lines of -99.999's in .200 file
      real*8 fracnotused  ! badsonic, 
      integer numbadsonic, numfiltsonic, numbadli, numeitherbad

!...................................................................................................
!    added numbadsonic & numnotused to common block flx4byt
!    09-may-2002 it was removed, because now it's passed as argument to subroutines
!...................................................................................................

!     variables in common for flux calculation subroutines:
      common /flx8byt/ rot
      common /flx4byt/ u,v,w,t,n

!     declaration of mirror file parameters:
!     --------------------------------------
      real*8 mirfreq,mirlen
      integer mirnum
!     where:
!     mirfreq = frequency of mirror data
!     length = averaging interval for mirror data in seconds (1/2 hour is 1800 seconds)
!     mirnum = number of points to be aggregated      
      parameter(mirfreq=0.5, mirlen=1800.)
      parameter(mirnum=mirfreq*mirlen)

!     variable declarations:
!     ----------------------
      real*8 mirt(mirnum)
      real*8 mirdew(mirnum)
      real*8 mirsumt,mirsumdew,miragdew,miragt
!     where:
!     mirt--mirror temperature readings		2008->temp
!     mirdew--mirror dew point readings		2008->RH
!     miragt,miragdew--aggregated temp and dew point for period
!     mirsumt, mirsumdew--summation variables for averaging
      real*8 statcm, statcal(0:3)
      
!     common declarations:
!     --------------------
!     mircon--array of conversion factors for chilled mirror data
!     conversion from volts to degrees. see documention in lmirror.for file under the 
!     getmirfacs() routine.
      real*8 mircon(2,2,2)
      common /mirconv/ mircon

!     status bit variables
!     --------------------
!     bitarray--array of status bits 15,14,13,12 for a given status word
!     bitarray(1)= 1 or 0 {bit15}
!     bitarray(2)= 1 or 0 {bit14}
!     bitarray(3)= 1 or 0 {bit13}
!     bitarray(4)= 1 or 0 {bit12}
!
!     numbitarray--array of number of flags (i.e. 1's) of each bit in a given half hour.
!     numbitarray(1) = 0 to 14,400 {bit15}
!     numbitarray(2) = 0 to 14,400 {bit14}
!     numbitarray(3) = 0 to 14,400 {bit13}
!     numbitarray(4) = 0 to 14,400 {bit12}
      
      integer diagarray(4)
      integer numbitarray(4)
      common /statusblock/ diagarray       

! {note: all variables must be declared before any assignment statements are used}

