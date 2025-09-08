!===================================================================================================
!                                     lba program: lco2vars6                                       !
!===================================================================================================
!      used to declare variables in common for lcalpf module for calibrating the profile system.
!      note:  parameters defined in lparams7 'include' module.
!
!      15-dec-2004:  change all real*4 to real*8  : especially needed for times
!      01-mar-2009:  includes dv2x2 one of the cal tanks is missing
!===================================================================================================

      real*8 tzero(nzdim), vzco2(nzdim),vzh2o(nzdim)
      real*8 pcellz(nzdim), tcellz(nzdim)
      real*8 flsamz(nzdim), flcalz(nzdim), prconz(nzdim) ! new 10-oct-2001
      real*8 vzair(nzdim), vzairh(nzdim), pzair(nzdim), tzair(nzdim) !%%% new 3-nov-2001
      integer zerflag(nzdim)  ! licor read-light status for zeros
      ! =0 if zero ok, =1 if licor ready light off
      ! =100 if this zero segment is too short
!     calibrations:
      real*8 tcalave(ncdim), tcal(ncdim,7)  ! average time, specific times 1-7 as in vcal
      real*8 vcal(ncdim,11) 
!     vcal(n,*):  avg voltage of each segment (n) of each cal-sequence 
!     n: 1=ls;    2=ms;        3=hs;        4=zero2;      5=surv std;    6=dilution;   7=air  (co2)
!        8=zero2; 9= (dummy); 10=dilution; 11=air  (h2o)
      real*8 pcellc(ncdim,4), tcellc(ncdim,4)                 ! cals (was p/tcal for lfluxes)
      real*8 flcal(ncdim,4),  flsam(ncdim,4), prcon(ncdim)    ! flow for cal and sample periods
!     p/tcellc(*,n):  n:  1=main cal, 2=hi-flow zero, 3=dilution, 4=air
!     2nd dim changed to 4 on 4-oct-01 in order to allow for separate data on air-segment
!     previous assumption that dilution p,t applies to air sample as well 
!     was not necessarily good!  (see fig ed_cals.010915_detail.ps in python)

      integer calflag(ncdim)  ! licor ready-light status for cals
      !=0 if zero ok, 
      ! =xx1 if licor ready light off
      ! =x1x if cal low, high-span signals indistinguishable 
      ! =1xx if this cal-sequence has a segment that is too short to use.

!     surveillance standards:
      real*8 tss(nsdim), vss(nsdim), vssh(nsdim), pcells(nsdim), 
     &       tcells(nsdim)
      !  tss=time, vss=co2 voltage, vssh=h2o voltage , p, t
      real*8 pcons(nsdim), flsams(nsdim), flcals(nsdim)
      integer liflag(nsdim)  ! licor ready-light status

!     fitting co2 data 
      integer licoroff                              ! num times licor off in interval
      real*8 coef1(ncdim),coef2(ncdim),coef3(ncdim) ! calibration coefs 1,2,3
      real*8 v0(3),dv3x3(3,3),coef(3),coefss(3)     ! zero, volatages, coefs for single fit, 
						    ! coefss as to force values of surface standard cal tank (nrc 180120)
!nrc................................................................................................
      real*8 cref2x2(2),dv2x2(2,2),coef2x2(2)       ! just in case one of the cal tanks is missing
!nrc................................................................................................
!  variables for eddy/profile co2 calibrations
      integer*4 iecal, izero, iss, nofecal, nofzero, nofss, ssflag
!
      common /co28byt/ tzero, tss, tcalave, tcal  
      common /co24byt/ vcal, vzco2, vzh2o, vss, vssh ! vssh added 23-oct-01 (missing)
      common /co2ints/ iecal, izero, iss, ssflag     ! rl

!  indices to vcal array for calibrations:        
      parameter (ihs=1,ims=2,ils=3)   !   hi, mid, low span
      !   skip normal zero (which would come next)
      parameter(iz2=4)   !   zero#2 (hi-flow)
      parameter(isurv=5) !   surveillance standard (or second zero2)
      parameter(idil=6)  !   dilution
      parameter(iair=7)  !   post-cal air sample
      parameter(iz2h=8,isurvh=9,idilh=10)
      !   zero2, surveil. std, dilution (water signal)
      parameter(iairh=11)!   post-cal air (water signal) 

!  indices to v201 data array (201 data) (indices common to eddy and profile):
! 
      parameter(ipcon=1)          ! pressure at p-controller
      parameter(ifsam=2, ifcal=3) ! sample and cal flow
      parameter(ipcell=4)         ! cell pressure
      parameter(irl = 5)          ! ready-light
      parameter(itcell=6)         ! cell temperature
      
! -----------
!     now, define vars used only locally (in lcaled):
!     input arrays for eddy data:
      integer*4 is200(2), is201ed(4)  ! number of status numbers
!     sid cal , sid cal mux cm
!     note:  in this program (as opposed to lsplit) individual status 
!     bits are not stored separtely, but are grouped as above
      real*8 v200(7), v201ed(10) ! number of voltage fields on 200, 201 lines
!     200: co2 h2o u v w t cstat 
!     201: pcon flsam flcal pcell rl tcell tdetect tpump tdew tamb

!     working array:
      real*8 work(nwdim,2)  ! 1=co2, 2=h2o
      ! nwdim default=9600 (num fast points per cal seq)
      real*8 pcell, tcell     ! pressure, temp during cal
      real*8 pcon,fsam,fcal
      real*8 vrl              ! voltage of rl (typically ~1.9 if on)
      real*8 co2air, co2dil, h2oair, h2odil   ! co2 (ppm) and h2o (mmol/mol) 
      !  air and dilution samples (passed to getconc)

!     subroutines of lcaled:
      integer*4 nfirst, ifirst, istore  ! ssflag
      real*8 vco2, zco2, zh2o, z2co2, z2h2o
      common /ed4byt/ work

!     input arrays for profile data:  
      real*8 v201PF(12)       ! 201-input line voltages (see l1split for details)
      integer*4 is201PF(24)   ! 201-input line status bits
!     scalar vars for profile data
      real*8 xpco2,xph2o      ! co2 (ppm), h2o (mmol/mol)
      real*8 xco2ss, xh2oss   ! co2 from surveillance standard (forced) (nrc 180120)
      real*8 xpco2_v,xph2o_v  !keep track of volatages for output -calc cal and h2o correction factors
      integer*4 hhl, mml      ! low-span hours and minutes
      logical newstat         ! check for new status in input data
!     sizes for working arrays (max data points per sampling segment):
      parameter (npfdim=2*npwdim) 
      !** increase the array size to account for the longer cal/zero
      !** segments in the re-installed instruments. 
      real*8 work1(npfdim),work2(npfdim),work3(npfdim),work4(npfdim)
      real*8 work5(npfdim),work6(npfdim),work7(npfdim)
!     for aggregating (1-4): co2, h2o, p, t (new:  fsam, fcal, pcon)
!     npwdim= (freqlo*60*thubseq/10)*3 = 180, set in 'lparams.for' module)
!        
!     output arrays for co2 samples (level, tprfly, vpco2, vph2o, pcprf, tcprf)
!     are still defined within lcalpf because of their large size
!     --  arrays sized for 2-week data period:  30/hr * 24 * 14 days = 10080 (default) 
    
