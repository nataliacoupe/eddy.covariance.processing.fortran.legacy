!===========================================================================================
!                                  lba code: lparams.for
!===========================================================================================
!     as per 2008
      parameter (freqhi=8, freqlo=0.5, thubseq=20, tcalint=720, tedzint=120)  !nrc 200829
!      parameter (freqhi=8, freqlo=0.5, thubseq=20, tcalint=360, tedzint=120) ! 2002-2006 nrc 180302
!                                                                               no need as tcalint is not used
!     where:
!     freqhi           = number observations per second high frequency
!     freqhi           = number observations per second low frequency (freqlo=0.5 every 2 seconds)
!     thubseq          = 20 min (* 60 s/min = 1200 s)
!     tcalint(min)     = cal-interval time 
!     tedzint(min)     = zero-interval time for eddy  (cal is same as thubseq)

      parameter (nsec=7510, maxhigh=nsec*freqhi, maxlow=nsec*freqlo)
!     nsec    = max seconds of data to store (3600 for 1 hour, 1800 for 1/2 hour, etc).
!     maxhigh = max high-frequency for 1-hr of flux aggregation = 28880.
!     maxlow  = max low-frequency 1805

!     parameters for co2 calibration (sized for 2-weeks data period):
!     ---------------------------------------------------------------
      parameter (ncdim=14000, nezdim=5600)		!nrc increased nezdim from 1400 to 5600 for longer/more zeros

!     where:
!     nezdim = (1440/tedzint)*14
!     ncdim  = (1440/tcalint)*14

      parameter (npzdim=(1440/thubseq)*15)!1440=60*24 sec/day
!     ncdim  = size for cals:             1440/tcalint (def=4/day) * 14 days =  56 default
!     nezdim = size for eddy zeros:       (1440 min/day)/(min/zero)[=12] *14 = 168 default
!     npzdim = size for profile zeros:    (1440/day)/(min/hubseq)*14 d = 1008 default
!-------------------------------------------------------------------------------------------
!     After reinstallation:
!     ncdim  = size for cals:             1440/tcalint (def=2/day) * 14 days =  28 default
!-------------------------------------------------------------------------------------------

      parameter(vminrl = 0.5)     ! minimum voltage for rl for licor to be on.

!     parameters for data-storage arrays: 

!      parameter(tgap201 = 1.5)   ! minimum time (in min) missing in 201 file to be flagged 
      parameter(tgap201 = 3.)     ! minimum time (in min) missing in 201 file to be flagged nrc 180515
                                  ! as a gap.

!     eddy:
      parameter (nsdim=10, nwdim=6*freqhi*thubseq*60)       ! default: 19200
!     where:
!     nsdim = surveillance stdandard tests: up to 3/wk (only 1 anticipated tough).
!     nwdim = working array (maximum fast data points per calibration sequence)
!     default:  1200 sec/hubseq * 8 hz = 9600
!     (note: for eddy, mostly only 8 of 10 segments of hub-seq stored: zero,hs,ms,ls,zero, 
!     zero2, dilution, surveilance standard)

!     profile:
      parameter (npwdim=(freqlo*60*thubseq/10)*3)           ! default = 180.
!     sizes for working arrays (maximum data points per sampling segment).
!     hub-sequence divided to 10 segments (def=2 min each) at freqlo (def=0.5).
!     normally, it will sample 1 segment and, at end of calibration sequence, it will 
!     sample 2, but uses 3, to allow extra points (default nwdim=3*.5*120=180 points).

      parameter (npdim=60/(thubseq/10)* 24 * 30) 
!     default=10800  (changed from 15 to 30 days - 120823 ktw)
!     co2 sample:  arrays sized for 2-week data period: 30/hr*24*14 days = 10080 (default) 

!     parameters characterizing datafile structure:

!     number of variables on input file lines (non-time, non-status bit):
      parameter(n1PF=12, n12PF=32, n13PF=42)       
!     PF: 201, 201+202, 201+202+203 be: n1PF=12,n12PF=32, n13PF=42

      parameter(n1GD=7, n12GD=31, n13GD=39)        
!     GD: same as above. datafile parameters are symbolic integer constants (it allows them
!     to be used in 'case' statement of function xsum).

!     number of status variables on input files: 
      data nspf/24/                              ! number of status bits (profile)
      data nsgd/18/                              ! number of status bits (ground)

!     number of mux'd input variables:
      data nmuxpf/0/                             ! no mux'd variables anymore.
      data nmuxgd/0/

!     output-file parameters:
      data noutpf/51/        ! total output variables (hkp.clm, incl status, not incl date).
      data nhkppf/21/        ! number of output that are housekeeping (including status).
      data noutgd/50/        ! ...and similarly for ground.
      data nhkpgd/25/

!     input-file parameters:
      include 'lpaths.for'   ! path definitions (depends, in part, on windows or unix os)
                             ! defines root pathnames rootsplit, rootprocin, rootprocout and
                             ! data directories, e.g.:
                             ! data dirsplit/'split'/      <- input from *.20? files
                             ! data dirprocin/'process'/   <- input from cal/zero files *.zer, *.cal
                             ! data dirprocout/'process'/  <- output for processed data
!     system variables (input from 'lsysdat.inp' file)

!     general variables / variables used in common

      real*8 cref(3)                              ! reference-gas concentrations(for calibrations): 
                                                  ! high span, middle span, and low span.
      real*8 fcoef(3)                             ! fake calibration coefficients -only when no cal tanks available 3xPF 3xEB nrc
      real*8 cvapor                               ! multiplier for h2o correction nrc 170127 
      real*8 csurv                                ! surveillance standard concentration.
      real*8 flsamed, flsampf, flcaled, flcalpf   ! standard flow rates.
      real*8 zflsam, zflcal, zpcon, zpcell        ! zero voltages for flow/pressure-meters
      real*8 gflsam,gflcal,gpcon, gpcell, gtcell       
      real*8 clicor(5),hlicor(3),t0co2,t0h2o
      real*8 tlagco2,tau,tauh                     ! in seconds
      real*8 gu,gv,gw,gt,dirofu
!     where:
!     g* = gain of u, v, w and t (sonic data)
!     dirofu = direction of u, degrees from north clockwise

!     eddy variables
      real*8 tdelayed,tavecal,tflstep,tovlap   
!     minimum time waiting lag (after solenoid set) before measurement time to average together for 
!     calibration.
!     total time calibration-gas flowing during a single calibration step. times correspond to 
!     number of points (below):

      integer*4 ndelayed,navecal,nflstep,novlap           ! fast data (8 hz - nominal)
!     minimum number of points waiting lag (after solenoid set) before measurement.
!     number of points to average together for calibration, total number of points cal-gas 
!     flowing during a single calibration step.

      integer*4 nlaglo, navecalo,nflsteplo                 ! slow data (1/2 hz nominal)
!     variables corresponding to above, but for slow data.

      integer*4 kco2, kh2o
 
      character*5 filtcode

      integer filtset(5)          ! user-specified filter (incl. 5 flags)
      integer filtvec(maxhigh)    ! sonic filter vector
      integer badlicorrl(maxlow)  ! licor-bad indicator (slow data from 201 file)
      integer licorrlvec(maxhigh) ! licor filter vector (stretched for fast data)

!     initialization to default to filter on second bit

      data filtcode/'11110'/          ! default filter, 09-may-2002:  changed from '01000'
      data filtset/5*0/               ! actual filter (numeric), will be set by "filtcode"
      data filtvec/maxhigh*0/

!     profile variables 
      real*8 tavepf,tdelaypf      ! time(s): average for measurement,lag before measure,tot segment.
      integer*4 navepf, ndelaypf  ! number of points corresponding to tavepf, tdelaypf.
      integer*4 lagco2,lagh2o     ! number of data points to lag between licor, sonic.

!     data variables (mostly from ground system) aggregated by laggr
      real*8 cc(60),offset(60)    ! conversion coefficients for converting voltages to physical 
                                  ! values (see lconv.inp)
      real*8 sxo                  ! sd of output summary variable
      integer*4 nvarin, nmux
      integer*4 hh,mm             ! timing variables 
      real*8 ss                   ! used to convert fraction of day to hour:min:sec
      integer*4 is201(24), stbit(24,maxlow), kcal(0:3)
      integer*4 n1, n2, n3, ns, nskip  
      real*8 xo(0:60)             ! used to write to output files (note that zero-element is not an 
                                  ! output)
!     data file / command line housekeeping variables:

      character*6 platform                            ! watcom or gnufor

!     files names:
      character*8  lname                              ! name of program module
      character*4  lyear                              ! year of data 
      character*6  yymmdd, ftypes, tempstr
      character*2  dtype
      character*12  ifilen(8)                         ! up to 8 input files of 8 chars each
      character*120 fnamein, fname, fnameout          ! name of input, output files 
!      character*12 ifilen(6)                         ! up to 6 input files of 8 chars each
!      character*96 fnamein, fname, fnameout          ! name of input, output files 
      integer*4     fnum,enum                         ! file number, eddy number

!     command line:
      character*32  arg
      character*100 cmdline
      integer*4 fgetcmd,cmdlen 

!     timing:
      real*8 time0,time1,time2,time3,dtime,period     ! double-precision for timing.
      real*8 tstart,tstart1,tnext,tmid,tnext2,tnext3  ! tnext2,3 added 20-jul-04 for 
                                                      ! year-crossing test.
      real*8 year, days00, yrdays                     ! days since 1/1/00, days in the current year.

!     lsplit log reporting, ascii:
      logical verbose, ascii

!     common declarations:
      common /sys8byt/ tdelayed, tavecal, tflstep, tovlap, tavepf, 
     &                 tdelaypf, time0, days00, year
      common /sys4byt/ cref, csurv, flsamed, flsampf, flcaled, flcalpf, 
     &                 zflsam, zflcal, zpcon, zpcell, gflsam, gflcal,
     &                 gpcon, gpcell, gtcell, clicor, hlicor, t0co2, 
     &                 t0h2o, tlagco2, tau, tauh, gu, gv, gw, gt, sxo, 
     &                 dirofu, cc, offset, ss, fcoef, cvapor
      common /sysints/ ndelayed,navecal,nflstep, novlap, nlaglo, 
     &                 navecalo,nflsteplo, kco2, kh2o, nlagco2,nlagh2o,
     &                 navepf, ndelaypf, nvarin, nmux, hh, mm

