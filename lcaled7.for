!===================================================================================================
!                                     lba code: lcaled6
!===================================================================================================      
!  Created by Scott Saleska       (started Jul 2000)
! 
!   LcalEd5      : so named to match version number in Lflux5 (01-July-04)
!    16-Nov-04   : added routine to read cal tank numbers from file
!                  Lcaltank.inp, and select correct date.
!                  Also, changed headers for output files.                
!   
! 
!   LcalEd4      : Fix boundary-crossing problem for year (19-Jul-04)
!                  (NB still need to fix surveillance std)
!                  kind of inelegant:  relies on cals and zeros at new year's 
!                  midnight being mis-labeled as
!                  day-of-year 366.00xx (367.00xx in leap year).
!  
!   LcalEd3
!    13-Dec-01   : format line 517 (at ~812): increased precision of cal coefs (to 6)
!    05-Nov-01   : corrected header-line of *.ZER file to 9 (not 6)
!    23-Oct-01(c): added vssH to common block (was missing, so not calc water for Surv.Std)
!      *.ZER:  added prcon (pressure cal? control), flsam (flow during sample),flcal (cal) output 
!      *.CAL:  - corrected order of interpolated CO2-zero output to low,mid,high
!              - corrected value of air,dilution voltages (to make NOT zero-adjusted)
!              - corrected: store co2(v), not h2o, in co2 array [work(,1)] for air sample,
!                and store h2o(v) in h2o array [work(,2)] instead of nothing
!                ->line ~400:  'work(istore,2) = v200(ih2o)'  [was mistakenly 'work(istore,1)']
!                  (co2air values calculated before 23-Oct were based on water voltage!)
!    11-Oct-01(b):  few slight changes to run on python
!     7-Oct-01:  completed re-doing of Licor ready-light check
!       - new format for *.CAL, *.SS output files 
!          (includes lots more diagnostic info about cals, including flow info,
!       - *.SS:  corrected order of output columns:  P,T and co2,h2o were reversed
!     24-Sept-01:  add common (across L*.for modules) front-end 'Linitial.FOR'
c
!   LcalEd2
!     09-Sep-01: change name of gain variables for Licor temp & pressure
!               (gTecell -> gTcell;  gPecell --> gPcell)
!               (and similarly in LcalPF for profile Licor)
!     early Aug-01: add calFlag routines, and output calFlag too
!   Upgrades (incorporated in version stored as LcalEd1.for):
!     13-Jul-01: detect if Licor is on from ready light, avoid calc. cal coefs if
!                signal is too small
!     12-Jul-01: finish up date to allow Lparams file to set input/output paths
!                fix 10x output of diluted, undiluted sample; in *.CAL file (line ~640)
!     13-Jun-01: replace "{T/N}lag{ed/pf}" with "{T/N}delay{ed/pf}"
!     18-May-01: 'path' parameter added for finding input files (set in Lparams include file)
!                modified to output zero/cal times correctly
!                (line 323:  time0 = time1, b/c getzero/cal subroutines use time0)
!  
!    version 1.0 compiled 03-Apr-01 
!        (default directories set for field processing)
! 
!    preliminary versions:
!     version 0.3 compiled 01-Mar-01 (moved key subroutines to include file)
!     version 0.2 compiled 26-Feb-01 (accomodate change in LSysdat.inp, fix
!                                     to 'getconc' subroutine)
!     version 0.1 compiled 03-Nov-00
!     (adapted from part of Boreas code B2SPLIT by S.C. Wofsy and S.-M. Fan)
! 
!     ---------------------------------------------------------------------------------------
!     calibrates eddy data
!     invoke: 'lcaled7 yymmdd -files 12'
!  
!     reads through high-frequency eddy data to extract the calibration and zeroing data, estimates 
!     calibration coefficients, and write to a calibration file for later application in deriving 
!     co2 concentrations for samples.
!
!     input data files:
!       .\lsysdat.inp = input parameters describing lba system for calibrations
!       .\yymmdd\split\e*yymmdd.200 = scaled data for id=200 on large disk (c:)
!       .\yymmdd\split\e*yymmdd.201 = scaled data for id=201 on large disk (c:)
!     ouput data files:  (* is a or b for eddy1, eddy2)
!       .\yymmdd\process\e*yymmdd.zer = eddy zero data (mainly zero of co2 and h2o)
!       .\yymmdd\process\e*yymmdd.cal = eddy cal files (mainly cal-coefs for calculating co2)
!       .\yymmdd\process\e*yymmdd.ss  = eddy surveillance std data (if present)
!
!     note: input data (*.200) consists primarily of continuous high-frequency data (8 hz), 
!           interrupted by:
!     (1) a 2-minute zero segment (every 2 hours--contrast every 20-min in profile)
!     (2) a (nominal) 20-min cal sequence (every 6 hrs), broken into the following 2-min segments:
!                                                                 cal bits ('iscal')
!          0. zero (coincides with end of profile sequence)        01 
!          1. high-span for 'tflstep' (120) seconds                01
!          2. mid-span for 'tflstep'  (120) seconds                01
!          3. low-span for 'tflstep'  (120) seconds                01
!          4. zero                                                 01
!          5. zero #2 (zero at hi flow = 8000 sccm)                11
!          6. zero #2 (or surv std if time -- 1 per week)          11 (01)
!          7. dilution (4000 sccm zero/ 4000 sccm sample air)      10
!          8. sample (to finish out sequence of 10 segments)       00
!          9. sample (to finish out sequence of 10 segments)       00
!         10. sample (to finish out sequence of 10 segments)       00
!  
!     (3) the zero and calibration segments of (1) and (2) are divided into:
!         a.  lag interval (tdelayed seconds for transition), and 
!         b.  a sample interval (tavecal seconds of stable data)
!             nrc 01-Nov-2017 increased the transition as low values of H2O are 
!             embeded on the sample (lsyst.inp)
!
!     overall program logic contains 5 steps:
!
!     loop through input files (EA/ebyymmdd.200):
!     [step 1]: sets up for input file (open, get right data from lsysdata.inp)
!     ---------------------------------------------------------------------------------------
! 
!     [step 2]: averages high-frequency data to get one point for each zero/calibration segment:
!     ---------------------------------------------------------------------------------------
!     loop through all of e*yymmdd.200 file; for each cal/zero sequence:
!     (a) stores all co2 (n=1), h2o (n=2) voltages in a temporary array:  work(n, istore) from 
!         point where iscal=01 until iscal=00 again.
!     (b) averages up each segment of the work array (not including first tdelayed seconds), 
!         put them into these arrays using 'getzero' and 'getcal' subroutines:
!         - zero data (normal zeros only)
!           vzco2(izero), izero = 1 to nofzero; 
!           tzero(izero)  time of zero (in decimal days)
!         - cal data:
!           vecal(n, iecal), where iecal = 1 to nofecal
!           n= 1: low-span    2: mid-span,        3: high-span   
!              4: zero#2      5: surv std/zer#2   6: dilution    7: air sample
!              8-11:  as above, but water signal
!           tecal(n, iecal) = times for each of above
!           tecalave(iecal)  = average of time for hs, ms, ls
!         - surveillance standard data (once per week):
!           vss(iss), where iss = 1 to nofss
!           tss(iss) = time of surveilance standard (ss)
!  
!     [step 3]: retrieves pressure and temperature data for each cal sequence from *.201 file
!     ---------------------------------------------------------------------------------------
!     (a) zeros:  loop through tzero(izero) times, get corresponding p and t from .201, 
!                 puts average torr, degc into: pcellz(izero), tcellz(izero)
!     (b) cals:   loop through tecal(iecal) times, get corresponding p and t, puts average torr,
!                 degc into: pcell(n,iecal), tcell(n, iecal)
!     (c) surveillance standards: loop through tss(iss) times, get corresponding p and t, puts 
!                                 average torr, degc into: pcells(iss), tcells(iss)
!     (d) writes the zero output file (to E*.zer)
!
!     [step 4]: fit the calibration coefficients
!     ---------------------------------------------------------------------------------------
!     loop through each calibration (1 to nofecal):
!     (a) gets the zeros for this cal (by interpolation)
!     (b) fits the cal coefficients -- store in: coef1(iecal), coef2(), coef3()
!     (c) calculates h2o, co2 concentrations for dilution and post-cal sample (use cal-seq zero, 
!         uninterpolated)
!     (d) writes the calibration output file (to E*.cal)
!
!     step 5: calculate surveillance standard concentration
!     -----------------------------------------------------
!     loop through each surv. std. measurement (1 to nofss):
!     (a) interpolates zeros
!     (b) interpolates cal coefficients
!     (c) calculates concentration
!     (d) writes the results (to E*.ss)
!============================================================================================
      program lcaled 
      include 'lparams.for' 

! parameters for array-sizes + sys vars
!    ! as of March 2001, Lparams.for included:
!     parameters for CO2 Calibration:  --  sized for 2-week data period :
!        parameter (ncdim=(1440/tcalint)*14, nezdim=(1440/tedzint)*14 )
!        parameter (npzdim=(1440/thubseq)*14 )
        !  ncdim=size for cals:  1440/tcalint (def=4/day) * 14 days      = 0056 default
        !  nezdim=size for eddy zeros: (1440 min/day)/(min/zero)[=12]*14 = 0168 default
        !  npzdim=size for profile zeros: (1440/day)/(min/hubseq)*14 d   = 1008 default
!    EDDY:
        ! cal-intervals (at most, one hub-sequence worth of fast data):
!      parameter (nsdim=10, nwdim=freqhi*thubseq*60)   ! default:  9600
        !  nsdim=surveil. std. tests: up to 3/wk (though only 1 anticipated): 
        !  nwdim=working array (max fast data points per cal sequence):
        !    default:  1200 sec/hubseq * 8 hz = 9600
        !    (note:  for eddy at most only 8 of 10 segments of hub-seq stored:
        !        zero, hs, ms, ls, zero, zero2, dilution, surv. std)
!    ! Lparams.for also defines, and puts into common-block 'sysvars':
!     SYSTEM VARIABLES (input from 'LSYSDAT.inp' file)
!    ! Lparams.for also defines DATA FILE/COMMAND LINE HOUSEKEEPING VARS
!      (used in retrieving command-line info and file-name info

      parameter(nzdim=nezdim)
      include 'lco2vars7.for'

!     variables for doing CO2 conc, including:
!     zeros:
!        real vzco2(nzdim),  vzh2o(nzdim), tzero(nzdim), 
!        real pcellz(nzdim), tcellz(nzdim)
!     cals:
!        real tcalave(ncdim), tcal(ncdim,7)  ! average time, specific times 1-7 as in vecal
!        real vcal(ncdim,11), pecell(3,ncdim), tecell(3,ncdim)
          !  vcal(n,*):  avg voltage of each segment (n) of each cal-sequence (*)
          !  n: 1=hs; 2=ms; 3=ls; 4=zero2; 5=zero2; 6=dilution; 7=air (CO2)
          !                       8=zero2; 9=dummy;  10=dilution; 11=air                     (H2O)
!        real pcellc(ncdim,3), tcellc(ncdim,3)  ! cals (was p/tcal for lfluxes)
          !  p/tcellc(*,n):  n:  1=main cal, 2=hi-flow zero, 3=dilution, 4=air
!     fitting CO2 data 
!        real*8 coef1(ncdim),coef2(ncdim),coef3(ncdim) ! calibration coefs 1,2,3
!        real*8 v0(3),dV3x3(3,3),coef(3)               ! zero, volatages, coefs for single fit
!     eddy cals:
        ! integer*4 iecal,izero,iss,nofecal,nofzero,nofss, ssflag
        ! common /CO2vars/ iecal,izero,iss,time0,vcal,tcal,tcalave, 
!                          vzero,tzero,vss,tss,ssflag
!     ---------------------------------------------------------------------------------------
!     now, define vars used only locally (in lcaled):
!     Input arrays for eddy data:
!        integer*4 is200(2), is201(4)  ! number of status numbers
          !        sid cal mux cm
          !  note: in this program (as opposed to lsplit) individual status 
          !         bits are not stored separtely, but are grouped as above
!        real*4 v200(7), v201(10) ! number of voltage fields on 200, 201 lines
!               200: co2 h2o u v w t cstat sid cal mux cm
!               201: pcon flsam flcal pcell RL tcell tdetect tpump tdew tamb
! 
!     working array:
!        real*4  work(nwdim,2)  ! 1=co2, 2=h2o
!        real*8  pcell, tcell
!        logical LicorOFF
!        common  work

      logical done, debug 
      real*8  yrdays1, yrdays0, tdoy   
      integer doy, flagtank, segcal, j, k
      real*8  gmtsecs
      real*8  a1(nzdim),a2(nzdim),a3(nzdim),a4(nzdim),a5(nzdim)
      real*8  a6(nzdim),a7(nzdim),a8(nzdim),a9(nzdim)  			!nrc 170727
      integer idate                                  			!flag based on the date of the collection
      common /ewg/ gmtsecs, doy
!     set default initial parameters ----------------------------------------------------------------
      lname = 'lcaled'
      ftypes= '12'
      include 'linitial7.for'    ! common initialization code. sets directory structure (root paths,
                                 ! dirsplit, dirprocin, dirprocout), retrieves command line
                                 ! file names for processing [ to ifilen(1..nifiles) ] 
                                 ! displays initial program message.

!     indices to input voltage arrays:
      ico2 = 1                   ! indices to v200 for location of needed data
      ih2o = 2
!...................................................................................................
! nrc 100827 to deal with those periods when cal has double in time (e.g. 100820)
! variables were fixed, now can be changed for when cal is longer: 3 min cals: 1440 points
! exept for dil zero cal is 2 min long: 960 points
      if8 = 0.                    
      if4 = 960.                !dil cal time (calflag=10)
      flagtank = 0.	        !flag indicates no calibration tanks if 1
      segcal   = 3.		!original cal sequence was 898 seconds ZR-HS-MS-LS-ZR 180s (3min) per step
      				!increase to 5min if cal sequence is 40min instead of 15min or 2 if 10min total
! nrc ..............................................................................................

!---------------------------------------------------------------------------------------------------
!     eddy co2 calibration:
!     yymmddxx.dat = raw data from cr10 in binary form
!     xxyymmdd.200 = scaled data for id=200
!     xxyymmdd.201 = scaled data for id=201
!---------------------------------------------------------------------------------------------------


!===================================================================================================
!                               main program loop starts here
!===================================================================================================
      if (nifiles.eq.1) close(10)   ! will not be needing more from param file.
      do infile=1, nifiles          ! loops over each input file.

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!     first: opens the input file (E*.200) and output files:
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
      dtype= ifilen(infile)(1:2) ! dtype is data type string
      enum=1

!     parameter files:

      if(dtype.eq.'EB') then 
       enum = 2
      endif  
      print *,'  Reading parameter file lsysdat.inp ...'
      call readsysd(enum)            ! reads data from lsysdata.inp for eddy system 2.
      call readconv(enum)            ! reads data from lconv.inp for conversion coefficients.

!     raw data input (*.200,201 data):
      fnamein = ifilen(infile)
      call getfullpath( rootsplit,yymmdd, dirsplit, fnamein, lfnin) 

!     calibration data (processed *.cal, *.zer files) 
      fname = ifilen(infile)
      call getfullpath( rootprocin, yymmdd, dirprocout, fname, lfn)
      print *, '  Opening input files for calibration ',fnamein(1:lfnin)//'.200"..'
      open(20,file=fnamein(1:lfnin)//'.200',status='old') 
      read(20,*)                          ! header lines
      read(20,*)
!     opens output files:
      open(12, file=fname(1:lfn)//'.zer',status='unknown')
      open(13, file=fname(1:lfn)//'.cal',status='unknown')
      open(14, file=fname(1:lfn)//'.ss', status='unknown')
      open(99, file=fname(1:lfn)//'.erc',status='unknown')  ! error log file

!     puts headers in the output files:
!
!---------------------------------------------------------------------------------------------------
!     zero file
!---------------------------------------------------------------------------------------------------
      write(12,120) enum,yymmdd,platform,lname,iyear,imonth,iday,ihrs,
     &              imin
120   format('2  10  LBA Eddy', i1,', zero on: ', a6,' (',a6,' ',a7,
     &          ' run: ',i4,'-',i2.2,'-',i2.2,', ',i2.2,':',i2.2,')' )
      write(12,121) 
121   format('jdstart.gmt   doy   zco2.v   zh2o.v   pz.torr.   tz.c   pcon.t   flsam.sccm   fcal   zeroflag')   

!---------------------------------------------------------------------------------------------------
!     calibration file
!---------------------------------------------------------------------------------------------------
      write(13,130) enum,yymmdd,platform,lname,iyear,imonth,iday,ihrs,
     &              imin
130   format('6  44  LBA Eddy',i1,', cal on: ',a6, ' (',a6,' ',a7,
     &           ' run: ',i4,'-',i2.2,'-',i2.2,', ',i2.2,':',i2.2,')' )
      write(13,131) 
131   format('jdstart.gmt   doy   zlow.v   zmid.v   zhigh.v   dvlo   dvmid   dvhigh   calflag   ')  
      write(13,132)  
132   format('zco2.v   zh2o.v   z2co2   z2h2o   vdil   vair   vdilh   vairh   ')   
      write(13,133) 
133   format('pcal.torr   pz2.t   pdil.t   pair.t   pcon.t   tcal.c   tz2.c   tdil.c   tair.c')  
      write(13,134)  
134   format('fsamcal   fsamz2   fsamdil   fsamair   fcalcal   fcalz2   fcaldil   fcalair   ')   
      write(13,135)   
135   format('coef1   coef2   coef3   co2d.dil   co2air.umol.mol   h2od.dil   h2oair.mmol.mol   co2.hs   co2.ms   co2.ls') 


!---------------------------------------------------------------------------------------------------
!     surveillance standard
!---------------------------------------------------------------------------------------------------
      write(14,140) enum,yymmdd,platform,lname,iyear,imonth,iday,ihrs,
     &              imin
140   format('3  20  LBA Eddy',i1,', surveil. on: ',a6,' (',a6,' ',
     &         a7,' run: ',i4,'-',i2.2,'-',i2.2,', ',i2.2,':',i2.2,')' )
      write(14,141)   
141   format('jdstart.gmt   doy   zco2.v   zh2o.v   dvco2.v   dvh2o.v   coef1   coef2   coef3  ')   
      write(14,142)  
142   format('pss.t   pcon.t   tss.c   flsam.sccm   flcal   co2.ppm   h2o.mmol.mol   liflag   co2.ss   coef1.ss  co2.ppm.ss')   

      write(99,55) platform, lname, iyear, imonth, iday,ihrs,imin,isec
55    format(a6,' ',a7,':  on date (yy-mm-dd) ', i4,'-',i2.2,'-',i2.2,
     &             ' (',i2.2,':',i2.2,':',i2.2,'), process:' )
      write(99,*) '   ', fnamein(1:lfnin)//'.200'
110   format('  First GMT time in file is ',f8.4,' (',i5,' sec =',i3.2,
     &            ':',i2.2,':',f5.2,')')
      debug = .false. 
      call read20n(20,200,dtype,time0,v200,1,7,is200,1,2,done,yrdays1,
     &             debug) 
      if (done) goto 300 
      call writetimestamp(99,'  First GMT time in *.200 file is ',32,
     &                    time0)
      write(99,111) 
111   format('  Reading through file to find segments with calstatus', 
     &      ' bit set: ',/, 
     &      '     - store all data from zero and cal segments',/, 
     &      '     - then write zero file (EByymmdd.zer)',/, 
     &      '     - then fit cal-coefficients and write cal ',
     &      'file (EByymmdd.cal)')  
      write(99,*)
      yrdays0 = yrdays1                        ! saves yrdays at beginning of interval.

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!     second -- loop through 200 file to get & store zero and cal segments:
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
      filend = 0
      izero  = 0
      iecal  = 0
      iss    = 0
      ssflag = 0

!     file-loop:
!     ----------
      do while((iecal.le.ncdim).and.(izero.le.nzdim))! goes on until end-of-file causes transfer out 
                                                     ! or (unlikely) until array-size limits 
                                                     ! (nedim and nzdim) are reached.
                                                     ! at end, goes to label 290 if in middle of 
                                                     ! calibration sequence. otherwise, goes to 300.
      iscal=0
      do while(iscal.eq.0)                           ! loop until calibration sequence is found (cal=01)

200   format(i4,1x,f9.3,1x,2(f7.4,1x),4(f7.3,1x),f6.0,1x,i3,1x,i2) 
      call read20n(20,200,dtype,time1,v200,1,7,is200,1,2,done,yrdays1,
     &             debug) 
      if (done) goto 300 
      iscal = is200(2)                               ! cal= 00 : normal sampling

!     errors are reported in read20n subroutine.

205   continue
      enddo                    !  until (iscal.ge.01)
                               !      01 : normal cal or zero
                               !      11 : high-flow zero (zero #2)
                               !      10 : dilution zero (1/2 sample, 1/2 zero)

      call gethhmmss(time1)
      time0  = time1           !  getcal and getzero subs need time0
      istore = 0
      ssflag = 0
      write(99,210) time1, int(gmtsecs), hh,mm,ss
210   format('  Start: ',f10.4,' (',i5,' sec = ',i2.2,':',i2.2,':',
     &                  f5.2 ,') of cal/zero sequence.')
      do while(iscal.ne.0)     !  loop:  store all data from this calibraton/zero sequence.
       istore = istore+1
       work(istore,1) = v200(ico2)
       work(istore,2) = v200(ih2o)
       call read20n(20,200,dtype,time2,v200,1,7,is200,1,2,done,yrdays1,
     &             debug) 
       if (done) goto 290 
       oldcal = iscal
       iscal  = is200(2) 
       if ((oldcal.eq.11) .and. (iscal.eq.01)) then
        ssflag = 1            !   transition:  11--> 01 means surveillance standard
       endif
      enddo                   !   until (iscal.eq.00). ends when returns to normal sampling.

!...................................................................................................
!     Keep track the lenght of the calibration nrc 170110
!...................................................................................................
      read (yymmdd,'(i10)') idate
      if  ((idate.ge.080101).and.
     &     (istore.lt.1430) .or.
     &     ((istore.gt.1480).and.(istore.lt.11000)).or.
     &     (istore.gt.12000))  then
       print *,'Check cal lenght, see istore ', istore, ' should be 7200(Z) or 11520(full cal)'
      endif

!---------------------------------------------------------------------------------------------------
! 210614 Achtung !!! fow controlers changed
!        print*, 'start output to check'
!        print*, 'zflsam:', zflsam, ' gflsam:', gflsam, ' fsam:',fsam, ' zflsam:', zflsam, ' navecalo:',navecalo
!	 print*, 'zpcell:',zpcell, ' gpcell:', gpcell, ' pcell:', pcell
!        print*, 'zflcal:', zflcal,'  gflcal', gflcal, ' fsam:',fcal, ' zflsam:', zflcal, ' gpcon:',gpcon

!---------------------------------------------------------------------------------------------------

!...................................................................................................
!     note from nrc: starts 100824 as calibration time double after power outrages. see 100817 and 
!     original full cal: flag 01: Z-HS-MS-LS-Z flag 11: Z-Zhigh flag 10: Dil
!     100820 skips the first Z (only 1or2 secs) flag 01: HS-MS-LS-Z flag 11: Z-Zhigh flag 10: Dil: 
!            length h-m-l-z-zhigh: 19200 data points + skips single zero calibration
!            nsegs=7 for the "NEW" sequence
!...................................................................................................
!     !Total cal time ~40 minutes= 40*8*60 = 19200 (24+6+6+4)*60*8min  HS-MS-LS-Z + Z-Zhigh + DIL
!     instead of (15+6+2)*60*8min = 11520    Z-HS-MS-LS-Z(cal 01) + ZRhigh(cal 11)  + DIL(cal 10) 
   
      if ((idate.ge.080101).and.(istore.gt.18720).and.(istore.lt.19220).and.(if8.eq.0)) then
       segcal  = 6.   		      ! 4 step sequence HS-MS-LS-Z-Zhigh generally 3 min, now 6 min
       nflstep = segcal*60.*freqhi  ! originally 11520=(15+3+6)*60*8
       if8     = 1.
       if4     = 60*freqhi*4.       ! 960 = 60*8*2  as dilution takes now 4 min
       tflstep = 360.               ! 6 min * 60
       nflstep = tflstep*freqhi     ! number of points for each calibration step (hi-frequency)
       nflsteplo = tflstep*freqlo   ! number of points for each calibration step (low-freq)
       tflstep   = tflstep/86400.d0 ! eddy
       tovlap    = (tflstep + tdelayed)
       tovlap    = (tflstep + tdelayed)
       novlap    = int(tovlap*86400.*freqhi)   
       print*,'  Warning long cal: sequence is taking: ',  istore/(8*60), ' min instead of 23 min.'
      endif
!     !Year 2001 to 2006 original calibration was
!     !Total cal time ~16 minutes= 16*8*60 = 7680 (2+2+2+2+6+2)*60*8min  Z-HS-MS-LS-Z-Z + Z-Zhigh + DIL
!     !Zero cal (abscent since 2008) every 2 hours aditional 10times/day ~2 minutes = 2*8*60 = 980  
!      print*, 'CHECK tavecal: ', tavecal,' tdelayed: ',tdelayed
      if (idate.le.080101) then
         segcal  = 2.   	        ! 7 step sequence Z-HS-MS-LS-Z-Z-Zhigh 2 min each
         nflstep = segcal*60.*freqhi   ! originally 11520=(15+3+6)*60*8
         if8     = 1.
         if4     = 60*freqhi*2.        ! 960 = 60*8*2 as dilution takes now 2 min
         tflstep = 120.                ! 2 min * 60, now 180
         tavecal = 52.		        ! 2 min - 60 sec delay at low frequency  60sec/2sec. 26 points	
         tdelayed  = 60.               ! time delay... remove after calibration. 120 s as per 2008 onwards
         nflstep   = tflstep*freqhi    ! number of points for each calibration step (hi-frequency)
         nflsteplo = tflstep*freqlo    ! number of points for each calibration step (low-freq)
      	 navecalo  = tavecal*freqlo    ! number of points to average for slow (201) eddy cal data -"real" time to average for each cal 
         tflstep   = tflstep/86400.d0  ! eddy
         tavecal   = tavecal/86400.d0
         tdelayed  = tdelayed/86400.d0  
         tovlap    = (tflstep + tdelayed)
         print*,'  Initial setup 2000-2006: cal sequence is taking: ',  istore/(8*60), ' min'
       endif
!       print*, 'CHECK tavecal: ', tavecal,' tdelayed: ',tdelayed
!...................................................................................................
! To "return" to standard calibration time... see if8 going from 1 to 0 
! Zero cal 1440 Total cal 11520   
!...................................................................................................
      if ((idate.ge.080101).and.(((istore.gt.1430).and.(istore.lt.1480)).or.
     &     ((istore.gt.10000).and.(istore.lt.11540))).and.(if8.eq.1)) then        
       segcal  = 3.                  ! 5 step sequence ZR-HS-MS-LS-ZR now 3 min
       nflstep = segcal*60.*freqhi   ! originally 1440=60*24min(ZR-HS-MS-LS-ZR(cal 01)-ZRhigh(cal 11)-dilution(cal 10)
       if8 = 0.
       if4 = 60*freqhi*2.            ! 960 = 60*8*2 dilution takes 2 min
       tflstep = 180.                ! 3 min * 60
       nflstep = tflstep*freqhi      ! number of points for each calibration step (hi-frequency)
       nflsteplo = tflstep*freqlo         ! number of points for each calibration step (low-freq)
       tflstep   = tflstep/86400.d0       ! eddy
       tovlap    = (tflstep + tdelayed)
       tdelayed  = 120
       navecal   = tavecal*freqhi    ! number of points to average for a calibra
       print *,'  Warning: calibration sequence is back to 20 minutes.'
      endif

!...................................................................................................
!................................................................................................... 
      nsegs = nint( float(istore)/nflstep) 
      print*, 'number of segments (nearest int) in this sequence:', nsegs 
!     number of segments (nearest int) in this sequence:
!     nsegs=1: routine, isolated zero
!     nsegs=8: standard cal-sequence (zero,hs,ms,ls,zer,zer2,ss,dil)
!     nsegs=7: cal-seq, but no leading zero (should not happen usually)                       
!     nsegs=7: after re-installation calibration sequence is: hs ->ms ->ls ->zer ->ss ->zer2 ->dil
!     time2 = dble(julday) + gmt/86400.d0 --> time at end of calibration sequence:

      call gethhmmss(time2)
      write(99,230) time2,int(gmtsecs),hh,mm,ss,istore,nsegs, nflstep
230   format('  End: ',f10.4,' (',i5,' sec = ',i2.2,':',i2.2,':',
     &        f5.2,').',i5,' points stored (=',i2,' segments of ',
     &        i5,' points)')

!-----note from ewg (ago/2008):---------------------------------------------------------------------
!     For re-installation of the eddy system, the calibration sequence has changed: the 
!     calibrations are made each 3 min (i.e. 1440 points) long except for dilution zero, 
!     at the end of the sequence which is 2 min long (i.e., 960 points). 
!     Adjust here to get the correct number of segments.
!     Also need to account for the single zeros, which are 3 minutes long.
!---------------------------------------------------------------------------------------------------
      if (nsegs.eq.1) then 
       ansegs = float(istore)/nflstep 
      else  
       ansegs = float(istore-if4)/nflstep + 1.0		!nrc modified, before it was hardwired to 960 instead of if4
      endif 

!...................................................................................................
!     note from nrc: 100827 before fixed coefficient to 960 instead of variable "if4"
!     nflstep fixed at lsubs to nflstep=tflstep*freqhi now derived acording to the length to the cal sequence
!...................................................................................................
      if (abs(nsegs-ansegs).gt.0.1) then
       print*, 'QA Non-integral number of segments: ',abs(nsegs-ansegs),' skipping to next sequence'
       print*, 'number of cal segments nsegs: ',nsegs, ' and number of obs/time per segment ansegs: ', ansegs
       write(99,233) float(istore)/nflstep
233    format('  Non-integral number of segments (',f6.2, 
     &       ') in this sequence.',' skipping to next sequence...')
      else

!     Routine, isolated zero, izero as ID counter :
       if (nsegs.eq.1) then
        izero = izero+1
        call getzero(1)

!     Extracts mean zero voltages from 'work*' arrays, starting with beginning of array (element 1) 
!     and skipping over the first 'ndelayed' points. then, puts results in: vzero(1,izero) (co2) and
!     vzero(2,izero) (h2o).
!     Also, puts time1+tdeayed into  tzero(izero)  (delay accounts for transition):

        write(99,*) '  Zero(',izero,') time = ',tzero(izero)
       else

!     gets one more segment of normal air sample after end of calibration sequence in order to 
!     compare to dilution.

!...................................................................................................
!     note: flcal was fixed in 23-oct-2001: flcal was mistakenly writed as "work(istore,1)"
!...................................................................................................
       do i=1,nflstep
        istore = istore+1
        work(istore,1) = v200(ico2)
        work(istore,2) = v200(ih2o)
        call read20n(20,200,dtype,time2,v200,1,7,is200,1,2,done,yrdays1,
     &             debug) 
        if (done) goto 290 
       enddo
       write(99,234) nflstep,istore,ssflag
234    format('  ',i4,' more points stored for air sample = ',
     &       i6,' total (ssflag = ',i2,')')

!     -----------------------------------------------------------------------
!     routine: zero, hs, ms, ls, zero, zero2, ss, dilution (before tree fall)
!     routine: zero, hs, ms, ls, zero, zero2, zero2, dilution (after reinstallation)
!     nrc 170730 there is no ss between cals at 2002 I am not sure about this comment
!     0:zero, 1:hs, 2:ms, 3:ls, 4:zero, 5:zero2, 6:zero2, 7:dilution
!     -----------------------------------------------------------------------
      if (nsegs.eq.8) then  
!     first, gets the initial zero:

!     then, starts with element 1, extracts mean zero voltages from "work" arrays. 
!     results in: tzero(izero), vzero(1,izero), and vzero(2,izero)
       izero = izero+1
       call getzero( 1 )
       write(99,*) '  Zero(',izero,') time = ',tzero(izero)

!     then, calls the routine calibration sequence:
       iecal = iecal + 1
       call getcal( 1+nflstep )  	   !nflstep is the length of the first zero cal*8*60, if 3 min 1440 if 2 960
!     now, takes care of high, mid, low, and zero, plus post-calibration sequence (zero2, dilution,
!     and ss). 
!     extracts means, starting with work(#, 1+nflstep), and puts results in vecal(n, iecal) and 
!     tcal(n, iecal).
!     n=1-3 (hs, ms, ls);      
!     ??not true nrc 170730?? n=4-7 (zero2,ss,dilution,air);       n=8-11 (4-7, h2o)
!     n=4-7 (zero,zero2,zero2,dilution);
!     if surveilance standard ssflag=1, sets tss(iss) and vss(iss)
!     puts results of zero in: tzero(izero), vzero(1-2,izero) (and increments izero):
       write(99,*) '   Calibration (',iecal,'): HS time = ',tcal(iecal,1)
       write(99,*) '                          : MS time = ',tcal(iecal,2)
       write(99,*) '                          : LS time = ',tcal(iecal,3)
       write(99,*) '   Zero(',izero,') time............ = ',tzero(izero)
!...............................................................................................
!     routine-sequency, no leading zero:
      elseif (nsegs.eq.7) then 
       iecal = iecal + 1
       call getcal(1)
!......nrc added 170120.........................................................................
       print*, 'cal with no leading zero'
       write(99,*) '   Calibration (',iecal,'): HS time = ',tcal(iecal,1)

       write(99,*) '                          : MS time = ',tcal(iecal,2)
       write(99,*) '                          : LS time = ',tcal(iecal,3)
       write(99,*) '   Zero(',izero,') time............ = ',tzero(izero)
!...............................................................................................
!     takes care of high, mid, low, zero, and post-calibration sequence (zero2, dilution, ss) just 
!     as above.
      else 
       write(99,235) nsegs
235    format('  Incorrect number of segments (',i3,
     &       ') in this sequence.',' skipping to next sequence.')
      endif
      endif
      endif
      enddo  

!---------------------------------------------------------------------------------------------------
!     end and file loop. after this cal-seq stored, goes look for next one.
!---------------------------------------------------------------------------------------------------
      write(99,236) izero, nzdim, iecal,ncdim 
236   format('  Too many cals or zeros. End processing of *.200 ',
     &       'file after: ',/, 
     &       '     ->',i3,' zeros, max allowed =',i3,'(nzdim).',/,
     &       '     ->',i3,' cals,  max allowed =',i3,'(ncdim).')   
      goto 300
      
!---------------------------------------------------------------------------------------------------
!     third --  gets p and t for zero and calibration segments from *.201 file
!---------------------------------------------------------------------------------------------------
290   continue              ! jumps here if end-of-file is in middle of calibration sequence.
      write(99,*) '  End-of-file in middle of this cal-sequence'
300   continue              ! jumps here if end-of-file in is normal sample mode.
      close(20)
      nofzero = izero       ! saves counters.
      nofecal = iecal
      nofss   = iss
      write(99,302) nofzero, nofecal, nofss
      write(*,302)  nofzero, nofecal, nofss
302   format('  End of 200-file: ',i4,' zeros, ',i3,' cals, ',
     &            i3,' surveillance standard.') 
      write(99,*)
      open(21,file=fnamein(1:lfnin)//'.201',status='old')
      write(99,*) '  Opened input file ', fnamein(1:lfnin)//'.201', 
     &            '        to check licor status, get p and t.' 
      read(21,*)        !  saved for 1  header lines
      read(21,*)        !  2
      read(21,*)        !  3
      read(21,*)        !  4
!    
201   format(i4,1x,f9.3,/,8f8.3,/,f9.4,1x,f9.4, 1x,i3,1x,i2,1x,i3,1x,i1)
      call read20n(21,201,dtype,time1,v201ed,1,10,is201,1,4,done,
     &             yrdays1,debug) 
      if (done) goto 350 
      call gethhmmss(time1)
      call writetimestamp(99,'  first gmt time in *.201 file is ',
     &                    32,time1)   

!...................................................................................................
!     3a) loop through zero segments, get p and t, and check licor for these.
!...................................................................................................
      do izero=1,nofzero
       time0 = tzero(izero)         ! time of zero (after lagging from transition)
       write(99,207)izero,time0
207    format('  Search zero ',i3,' (start time =',f9.4,')...')
       do while(time1.lt.time0)
        call read20n(21,201,dtype,time1,v201ed,1,10,is201,1,4,done,
     &              yrdays1,debug) 
        if (done) goto 350 
       enddo                        ! until (time1 .ge. time0). found same time in *.201 file
       licoroff = 0
       pcell = 0
       tcell = 0                    ! initialize (?)
       pcon  = 0
       fsam  = 0
       fcal  = 0
       do j=1,navecalo              ! number of points at low-frequency (navecal was hi-frequency).
        call read20n(21,201,dtype,time2,v201ed,1,10,is201,1,4,done,
     &              yrdays1,debug) 
        if (done) goto 350 
        pcell = pcell + v201ed(ipcell)
        tcell = tcell + v201ed(itcell)
        pcon  = pcon  + v201ed(ipcon)  ! pressure at the controller (downstream)
        fsam  = fsam  + v201ed(ifsam)  ! sample flow
        fcal  = fcal  + v201ed(ifcal)  ! calibration flow
        vrl   = abs(v201ed(irl))       ! check licor ready-light
        if(vrl.le.vminrl) licoroff = licoroff + 1
       enddo
       pcellz(izero)  = gpcell*pcell/float(navecalo)   ! average over zero segment.
       tcellz(izero)  = gtcell*tcell/float(navecalo)
       prconz(izero)  = gpcon *(pcon/float(navecalo)-zpcon)
       flsamz(izero)  = gflsam*(fsam/float(navecalo)-zflsam)
!               print*, ' gflsam:', gflsam, ' fsam:',fsam, ' navecalo:',navecalo,' zflsam:', zflsam
       flcalz(izero)  = gflcal*(fcal/float(navecalo)-zflcal)
       zerflag(izero) = 0                             ! sets zero flag to ok.
       if(licoroff.gt.2) then
        call gethhmmss(time1) 
        write(99,330) 'zero', izero, time1, int(gmtsecs), hh,mm,ss
330     format('  Warning: LiCor was off during ',a4,' # ', i3,
     &        ' at time ', f9.4,' (',i5,' sec = ',i3.2,':',i2.2,
     &        ':',f5.2,')')
        zerflag(izero)=1                             ! zero flag set:  licor off
       endif
      enddo
      write(*,340) nofzero,  ' zeros ...  '
      write(99,340) nofzero, ' zeros ...  '
340   format('  Finished checking *.201 file for ',i3, a12)

	print*, 'after zero'
        print*, 'pcellz:', pcellz(izero), ' tcellz:', tcellz(izero), ' prconz:',prconz(izero), ' flsamz:', flsamz(izero), ' flcalz:',flcalz(izero), '  flsamz:', flsamz(izero), '   zflsam:', zflsam, '   gflcal:', gflcal
        print*, 'pcell:', pcell, ' tcell:', tcell, ' pcon:',pcon, ' fsam:', fsam, ' fcal:',fcal, ' vrl:',vrl

!...................................................................................................
!     3b) loop through calibration segments, get p and t for these.
!...................................................................................................
      rewind 21
      read(21,*)  !  1  header lines
      read(21,*)  !  2
      read(21,*)  !  3
      read(21,*)  !  4

      call read20n(21,201,dtype,time1,v201ed,1,10,is201,1,4,done,
     &             yrdays1,debug) 
      if (done) goto 350 

      do iecal=1,nofecal

!     iecal loop, part 1:   gets p and t, and averages for main calibration sequence.
      pcell = 0                           ! initialize for this calibration.
      tcell = 0
      pcon  = 0
      fsam  = 0
      fcal  = 0
      calflag(iecal) = 0                     ! sets calflag to ok.
      do i=1,3                               ! goes over main cal-sequence (hs,ms,ls).
       licoroff = 0                          ! licor checked separately for each span.
       time0 = tcal(iecal,i)                 ! time for this segment, after tdelayed.
       write(99,206)iecal,i,time0
206    format('  Search calibration ',i3,', segment ',i2,
     &       ' (start-time =',f9.4,')...')
       dowhile (time1.lt.time0)              ! loops to find start time 'time0'
       call read20n(21,201,dtype,time1,v201ed,1,10,is201,1,4,done,
     &             yrdays1,debug) 
       if (done) goto 350 

!     year-crossing check:
      enddo            ! until (time1 .ge. time0). founds same time in *.201 file.
      do j=1,navecalo  ! number of points at low-frequency (navecal is hi-frequency)
      call read20n(21,201,dtype,time1,v201ed,1,10,is201,1,4,done,
     &             yrdays1,debug) 
      if (done) goto 350 
!      print*, ' Check the PRESURE pcell during sample v201ed(ipcell): ', v201ed(ipcell)
       pcell = pcell + v201ed(ipcell)        ! pressure of licor cell
       tcell = tcell + v201ed(itcell)        ! temperature of licor cell
       pcon  = pcon + v201ed(ipcon)          ! pressure at the controller (downstream)
       fsam  = fsam + v201ed(ifsam)          ! sample flow
       fcal  = fcal + v201ed(ifcal)          ! calibration flow
       vrl   = abs(v201ed(irl))              ! check licor ready-light
       if(vrl.le.vminrl) licoroff = licoroff + 1

!     if ready-light voltage is below threshold, licor is off.
       enddo
       if(licoroff.gt.2) then 
        call gethhmmss(time1) 
        write(99,330) 'cal ', iecal, time1, int(gmtsecs), hh,mm,ss
        calflag(iecal) = calflag(iecal) + 1000 ! up to 3000 if out for all three spans.
       endif
      enddo                               ! i=3,1,-1 (main cal sequence)
!      print*, ' Check the PRESURE pcell during sample pcell: ', pcell

!     This needs to be modified when cal is longer                         nrc 160715
!     originally pcellc(iecal,1)=gpcell*(pcell/float(3*navecalo)-zpcell) ! averages over 3-segs of main
!      print*, ' Check the PRESURE pcellc during cal gpcell: ', gpcell,' pcell: ', pcell,' segcal: ', segcal
!      print*, ' navecalo', navecalo, ' zpcell', zpcell
!      print*, ' Check the PRESURE prcon during cal gpcon: ', gpcon,' pcon: ', pcon
!      print*,'segcal: ', segcal
!      pcellc(iecal,1) = gpcell*(pcell/float(segcal*navecalo)-zpcell)     ! averages over 3-min(?) of main
!      tcellc(iecal,1) = gtcell*tcell/float(segcal*navecalo)              ! calibration sequence.
!      prcon(iecal)    = gpcon*(pcon/float(segcal*navecalo)-zpcon)
!      flsam(iecal,1)  = gflsam*(fsam/float(segcal*navecalo)-zflsam)
!      flcal(iecal,1)  = gflcal*(fcal/float(segcal*navecalo)-zflcal)
      pcellc(iecal,1) = gpcell*(pcell/float(3*navecalo)-zpcell)     ! averages over 3-min(?) of main
!      print*, 'gpcell:', gpcell ,' pcell:', pcell, 'navecalo:', navecalo
      tcellc(iecal,1) = gtcell*tcell/float(3*navecalo)              ! calibration sequence.
      prcon(iecal)    = gpcon*(pcon/float(3*navecalo)-zpcon)
      flsam(iecal,1)  = gflsam*(fsam/float(3*navecalo)-zflsam)
      flcal(iecal,1)  = gflcal*(fcal/float(3*navecalo)-zflcal)

!     gains (g*) and zero voltages (z*) from lsysdat.inp file
!     iecal loop, part 2: get p,t average for:  zero2, dilution and air segments
      do i=2,4                 ! do 2=zero2, 3=dilution, 4=air (put values in pecell(,i), tecell(,i)
       it = 2*i                ! time index: 4=zero2, 6=dilution, 7=air sample
                               ! tcal(iecal,n) n=1:HS n=2:MS n=3:LS n=5:HZ 
       if(i.eq.4) it=7         !             7=air sample
       time0 = tcal(iecal,it)  ! time zero2 (tcal(,4)) or dil (tcal(,6)) or air (tcal(,7))
       do while(time1.lt.time0)
        call read20n(21,201,dtype,time1,v201ed,1,10,is201,1,4,done,
     &             yrdays1,debug) 
       enddo                   !  until (time1 .ge. time0) = found same time in *.201 file
       pcell = 0
       tcell = 0 
       pcon  = 0   
       fsam  = 0
       fcal  = 0               ! initialize
       do j = 1,navecalo  
        call read20n(21,201,dtype,time1,v201ed,1,10,is201,1,4,done,
     &             yrdays1,debug) 
        if (done) goto 350 
!        print*, 'CHECK the pressure of air v201ed(ipcell) :', v201ed(ipcell), ' i as in 2zero2 3dilution 4air is: ', i 
!        print*, 'CHECK time1: ', time1, '  time0: ', time0
        pcell = pcell + v201ed(ipcell)
        tcell = tcell + v201ed(itcell)
        fsam  = fsam + v201ed(ifsam)
        fcal  = fcal + v201ed(ifcal) 
       enddo
       pcellc(iecal,i) = gpcell*(pcell/float(navecalo)-zpcell) ! averages over points of this segment
!      print*, 'CHECK the pressure of air pcellc :', pcellc(iecal,i), 'i as in 2zero2 3dilution 4air :', i, ' navecalo: ', navecalo 
       tcellc(iecal,i) = gtcell*tcell/float(navecalo)
       flsam(iecal,i)  = gflsam*(fsam/float(navecalo)-zflsam)
       flcal(iecal,i)  = gflcal*(fcal/float(navecalo)-zflcal)  

!        print*, 'after cal'
!        print*, 'zflsam:', zflsam, ' gflsam:', gflsam, ' fsam:',fsam, ' zflsam:', zflsam, ' navecalo:',navecalo
!	 print*, 'zpcell:',zpcell, ' gpcell:', gpcell, ' pcell:', pcell, '   flsam(iecal,i):', flsam(iecal,i)
!        print*, 'zflcal:', zflcal,'  gflcal', gflcal, ' fcal:',fcal, ' zflcal:', zflcal, ' gpcon:',gpcon, '   tcell:', tcell
       
!      puts zero2, dilution, and air pressures in array elements 2,3, and 4.
       enddo                                ! i = 2,4 (zero2, dilution, air loop)
      enddo                                 ! iecal=1,nofecal

      write(99,340) nofecal, ' cals.    '
      write(*,340)  nofecal, ' cals.    '

!---------------------------------------------------------------------------------------------------
!     3c) loop through surveillance segments, get p,t for these
!---------------------------------------------------------------------------------------------------
      rewind 21
      read(21,*)  !  1  header lines
      read(21,*)  !  2
      read(21,*)  !  3
      read(21,*)  !  4
      call read20n(21,201,dtype,time1,v201ed,1,10,is201,1,4,done,
     &             yrdays1,debug) 
      if (done) goto 350 
      do iss=1,nofss
       pcell = 0
       tcell = 0                      ! initialize
       pcon  = 0   
       fsam  = 0
       fcal  = 0                      ! initialize
       licoroff = 0
       liflag(iss) = 0                ! sets calflag to ok
       time0 = tss(iss)               ! time of zero (after lag)
       dowhile( time1.lt.time0)      ! loop
       call read20n(21,201,dtype,time1,v201ed,1,10,is201,1,4,done,
     &             yrdays1,debug) 
       if (done) goto 350 
      enddo                         ! until (time1 .ge. time0). found same time in *.201 file.
      do j=1,navecalo                ! number of points at low-frequency.
       call read20n(21,201,dtype,time1,v201ed,1,10,is201,1,4,done,
     &             yrdays1,debug) 
       if (done) goto 350 
       pcell = pcell+v201ed(ipcell)
       tcell = pcell+v201ed(itcell)
       pcon  = pcon+v201ed(ipcon)
       fsam  = fsam+v201ed(ifsam)
       fcal  = fcal+v201ed(ifcal) 
       vrl   = abs(v201ed(irl))          ! checks licor ready-light
       if(vrl.le.vminrl) licoroff = licoroff + 1
      enddo
      pcells(iss) = gpcell*(pcell/float(navecalo)-zpcell)
      tcells(iss) = gtcell*tcell/float(navecalo)
      pcons(iss)  = gpcon*(pcon/float(navecalo)-zpcon)
      flsams(iss) = gflsam*(fsam/float(navecalo)-zflsam)
      flcals(iss) = gflcal*(fcal/float(navecalo)-zflcal)         

      if(licoroff.gt.2) then 
       call gethhmmss(time1) 
       write(99,330) 'surv', iss, time1, int(gmtsecs),hh,mm,ss
       liflag(iss) = 1  
      endif
      enddo 
      write(99,340) nofss, ' Surv. Std.'
      write(*,340) nofss,  ' Surv. Std.'

!     all is done. now extracts ps and ts:
      write(99,*) '  All finished with P and T from 201 file.' 
      write(*,*)  '  All finished with P and T from 201 file.' 
      close(21)
      goto 400
      
350   continue               ! comes here if finds unexpected end-of-file, i.e. while reading *.201
                             ! line to get time corresponding to *.200 line already extracted).
      write(99,355) ifilen(infile), time1, time0
355   format('  Error: encountered end of ',a8,'.201 at jday: ',f8.4,
     &        /,'          while looking for time ', f8.4)
      goto 900
!---------------------------------------------------------------------------------------------------
!     3d) writes the zero output file (E*.zer)
!---------------------------------------------------------------------------------------------------
400   continue               ! jumps here if it is done reading *.201 file.

!     writes the zero-data
!...................................................................................................
!     Before 2008 PF calibrated for zero 16 times between 0:00 and 12:00 and EA only 6 times
!     nrc 170726 "artificially" removed zeros when the fsam (sample flow) is greater than 100ssc (negative when cal occurs)
!      print*,'flsamz(izero): ', flsamz(izero),'  flcalz(izero): ', flcalz(izero)
!      print*, 'gflsam:',gflsam,'  fsam:',fsam,' zflsam:',zflsam
      izero = 1
      j = 0
      k = nofzero
      if ((idate.ge.080101).and.(idate.lt.240401))  then
       do while (izero.le.nofzero)
        ! nrc 180820 asume the sample flow is low < 100 scc and the cal flow is relatively high
        print*, 'sample flow is low < 100 scc and the cal flow is relatively high'
        print*, 'izero:', izero, ' nofzero:', nofzero   
        print*, 'flsamz(izero):', flsamz(izero), ' flcalz(izero):', flcalz(izero)   
        if ((flsamz(izero).le.1000).and.(flcalz(izero).gt.2000)) then		
         j = j+1
         a1(j) = tzero(izero)
         a2(j) = vzco2(izero)
         a3(j) = vzh2o(izero)
         a4(j) = pcellz(izero)
         a5(j) = tcellz(izero) 
         a6(j) = prconz(izero)
         a7(j) = flsamz(izero)
         a8(j) = flcalz(izero)
         a9(j) = zerflag(izero)
         izero = izero+1
         k     = k-1
        endif
        izero = izero+1
       enddo
       if (k.lt.nofzero) then
        print*, 'Warning: when sample flow was too high during cal, values were removed'
        j = 1
        do while (j.le.nofzero)
         tzero(j)   = a1(j)
         vzco2(j)   = a2(j)
         vzh2o(j)   = a3(j)
         pcellz(j)  = a4(j)
         tcellz(j)  = a5(j) 
         prconz(j)  = a6(j)
         flsamz(j)  = a7(j)
         flcalz(j)  = a8(j)
         zerflag(j) = a9(j)
         j = j+1
        enddo
        nofzero = k
       endif
      endif

      !............................................................................................
      ! Broken solenoid............................................................................
      !if (idate.gt.190101) then
      ! vzco2(izero) = vzco2(izero)-0.403+0.6423
      ! vzh2o(izero) = vzh2o(izero)+0.933       
      !endif
      ! ...........................................................................................
!---------------------------------------------------------------------------------------------------
! 210614 Achtung !!!
!---------------------------------------------------------------------------------------------------
! check zflsam --when we changed the flow controlers
        print*, 'vzco2(izero): ',  vzco2(izero), '     vzh2o(izero): ',  vzh2o(izero)
!---------------------------------------------------------------------------------------------------

!     writes the zero-data
      do izero=1,nofzero
       tdoy = tzero(izero) - days00 - yrdays1
      !nrc added set boundaries of what is possible given the write format constrains
       if ((vzco2(izero).lt.-99.).or.(vzco2(izero).gt.999.))       vzco2(izero)=-99.999	!f8.4
       if ((vzh2o(izero).lt.-99.).or.(vzh2o(izero).gt.999.))       vzh2o(izero)=-99.999	!f8.4
       if ((pcellz(izero).lt.-999.).or.(pcellz(izero).gt.9999.))   pcellz(izero)=-999.9  !f6.1
       if ((tcellz(izero).lt.-999.).or.(tcellz(izero).gt.9999.))   tcellz(izero)=-999.99 !f7.2
       if ((prconz(izero).lt.-999.).or.(prconz(izero).gt.9999.))   prconz(izero)=-999.99 !f7.2
       if ((flsamz(izero).lt.-9999.).or.(flsamz(izero).gt.99999.)) flsamz(izero)=-9999.9 !f7.1
      !..................................................................................
      !no more crazy dates
       if ((tdoy.ge.1.).and.(tdoy.lt.367.).and.(tzero(izero).ge.1.).and.(tzero(izero).lt.100000.)) then		                
        write(12,412) tzero(izero), tdoy, vzco2(izero), vzh2o(izero),
     &              pcellz(izero), tcellz(izero), prconz(izero), 
     &              flsamz(izero), flcalz(izero), zerflag(izero)
412     format(f10.4,1x,f9.4,2(1x,f8.4),1x,f6.1,1x,f7.2,1x, 
     &       f7.2,2(1x,f7.1),1x,i3)
       endif                                                             !endif no more crazy dates
      enddo
      close(12)
      write(99,*) '  Zero file writen       ', fname(1:lfn)//'.zer' 
      write(*, *) '  Zero file writen       ', fname(1:lfn)//'.zer' 

!---------------------------------------------------------------------------------------------------
!     fourth --  does calculations for calibration, and writes output file (E*.cal).
!                (interpolates zeros to calibration voltages, fits calibibration curve)
!---------------------------------------------------------------------------------------------------
500   continue       ! continues here when done with writing zero file.

!     loops through each calibration:
      jlo=1
      write(99,*)'  Fitting calibration coefficients...'
      write(*, *)'  Fitting calibration coefficients...'
      do iecal=1,nofecal

!...................................................................................................
!     4a. gets the zeros for this calibration [ v0(i), i=1,2,3 --> high,mid,low ]
!...................................................................................................
!      print*, 'CHECK at EDDY first time1: ', time1
      time1 = tcal(iecal,ihs)                 ! time for the hs
  
      call readcal(time1)                     ! readcal -- gets the appropriate cal tank values for 
                                              !            this segment.
! 
!      print*, 'CHECK at EDDY second izero: ',izero , ' nofzero: ', nofzero, ' time1: ', time1, ' jlo: ',jlo

      ! Broken solenoid............................................................................
!      if (idate.gt.190101) then
!       vzco2(izero) = vzco2(izero)-0.403+0.6423
!       vzh2o(izero) = vzh2o(izero)+0.933       
!       vzco2(izero+1) = vzco2(izero+1)-0.403+0.6423
!       vzh2o(izero+1) = vzh2o(izero+1)+0.933       
!      endif
      ! ...........................................................................................

      ! tzero times when there is zero calibration
      izero = 0                               ! nrc 170207 added initialize
      if (nofzero.gt.1)  izero=ihuntf(tzero,nofzero,time1,jlo)  
                                              ! izero:  index to last tzero before time1.
      if (izero.eq.0) then                    ! before the range of recorded zeros
       izero=1
       do i=1,3
        v0(i)=vzco2(izero)                    ! uses co2 zero 1.
       enddo
      if( zerflag(izero).gt.0) calflag(iecal)=calflag(iecal) + 3

!     setting calflag=3 means:  unable to extract nearby zero (reminder:  zerflag = 0 if zero ok.
!                               otherwise zerflag = 1)
      elseif (izero.eq.nofzero) then          ! after the range of recorded zeros.
       do i=1,3
        v0(i) = vzco2(nofzero)                ! use last co2 zero. 
       enddo
      if(zerflag(izero).gt.0) calflag(iecal) = calflag(iecal) + 3
       elseif (izero .gt. 0) then             ! within the range.
        jlo    = izero
        deltav = vzco2(izero+1)-vzco2(izero)
        deltat = tzero(izero+1)-tzero(izero)
        dvdt   = deltav/deltat
       if(zerflag(izero).eq.0) then 

!     if previous zero is good, check the next one.
        if(zerflag(izero+1).eq.0) then  

!     if both zeros are ok...
         do i=1,3                             ! hs,ms,ls
          delt  = tcal(iecal,i)-tzero(izero)  ! use interpolated value
          print*, 'time this zero cal tcal(iecal,i):', tcal(iecal,i),'  time past zero cal tzero(izero):', tzero(izero)
          print*, 'vzco2(izero):', vzco2(izero),'  dvdt:', dvdt, '   delt:',delt
          v0(i) = vzco2(izero)+(dvdt*delt)
         enddo                                ! (these should all be very close)
         print*, 'Both zeros are OK at calibration'  
        else   
!     previous zero ok, but next is bad --> don't interpolate.
         do i=1,3                             ! hs,ms,ls
          v0(i) = vzco2(izero)
         enddo                                ! (these should all be very close)
         calflag(iecal) = calflag(iecal) + 2  ! 2=next zero bad
         print*, 'Eddy previous zero is ok and the next is bad'  
        endif
       elseif(zerflag(izero+1).eq.0) then 
!     if previous zero bad, but next one ok, use next
        do i=1,3                              ! hs,ms,ls
         v0(i) = vzco2(izero+1)
        enddo                                 ! (these should all be very close)
        calflag(iecal) = calflag(iecal) + 1   ! 1=prev zero bad
        print*, 'Eddy previous zero bad, but next one ok, use next' 
!     co2 zero:
       else
!     both zeros are bad
        do i=1,3                              ! hs,ms,ls
         v0(i) = 0                            ! will not be used
        enddo                                 ! (these should all be very close)
        calflag(iecal) = calflag(iecal) + 3   ! 3=both zeros bad
        print*, 'Eddy both zeros are bad' 
       endif
      endif

!---------------------------------------------------------------------------------------------------
! 210614 Achtung !!!
!---------------------------------------------------------------------------------------------------
!  needed to hardwire zflsam --when we changed the flow controlers
!        print*, 'vzco2(izero)',  vzco2(izero)
!        if ((idate.ge.210401).and.(idate.le.220101)) vzco2(izero) = vzco2(izero)+0.16    
!        print*, 'vzco2(izero)',  vzco2(izero)
!---------------------------------------------------------------------------------------------------

!...................................................................................................
!     4b. gets the calibration curve coefficients for this calibration
!...................................................................................................
      do i=1,3                              ! high,mid,low
      print*,' eddy check voltages during zero v0(i): ', v0(i), '  for i:', i,'  dv3x3(i,1): ', vcal(iecal,i)
      ! IRGA's voltage correlation coefficients were reset to zero after lightning changes required for lparameters
      ! Hardwired IRG3-1046 CO2 cal C2 T=32.82, C2 K=18009, C2 A=1.36268E-1, C2 B=4.64832E-6, 
      ! C2 C=8.20996E-9, C2 D = -1.04493E-12, C2 E=6.3323E-17
      ! H2O cal values H2 T = 33.22, H2 K=15621, H2 A= 6.00539E-3, H2 B = 2.50887E-6, H2 C = 3.11171E-11
      ! if ((idate.ge.210730).and.(idate.lt.211231)) 
      !		vcal(iecal,i) = vcal(iecal,i)-0.15 	

       dv3x3(i,1) = vcal(iecal,i)-v0(i)     ! subtract off drift-corrected zeros
       dv3x3(i,2) = dv3x3(i,1)*dv3x3(i,1)   ! voltage squared
       dv3x3(i,3) = dv3x3(i,2)*dv3x3(i,1)   ! voltage cubed
       coef(i)    = 0                       ! initialize coefficient
      enddo

505   format('  Calibration # ',i2,':  cal-matrix',3(/,10x,'|',f7.1,
     &              '|   |',3(1x,f7.3),'|  |coefn|'))
      dvhilo = abs( dv3x3(1,1) - dv3x3(3,1) )  ! signal diff between high and low
      print*,'CHECK difference between high and low voltages, dvhilo: ', dvhilo
      print*,'typically ~1.3 for eddy licors, and ~0.4 for profile'
      
!     dvhilo is typically ~1.3 for eddy licors, and ~0.4 for profile

      if(dvhilo.lt.0.1) calflag(iecal) = calflag(iecal)+10

!nrc added.........................................................................................
      if((cref(1).gt.990.).or.(cref(2).gt.990.).or.(cref(3).gt.990.)
     &           .and.(dvhilo.lt.0.1)) then
       calflag(iecal) = calflag(iecal) - 10
       print*,'changes to calflag based on differences between voltage during hi and low (dvhilo)'
       print *,'  No cal tanks -borrowed coefficients'
      endif
!nrc................................................................................................

      print*,'Enter the cal coeff derivation, where cal tanks concentration are: ', cref, '  calflag:',calflag(iecal)
      if(calflag(iecal).lt.3) then 
!     only invert matrix if calflag was not raised, or if it was raised only due to bad zeros.
!     calflag = 1 if previouse zero bad, 2 if next, 3 if no nearby zero accesible
!nrc................................................................................................
!nrc  cref concentration in ppm, if cref==999 cal gas is off fill the matrix with other 
!     available tanks. If any tank concentration is >990
!nrc................................................................................................
      if((cref(1).gt.990.).or.(cref(2).gt.990.).or.
     &           (cref(3).gt.990.)) then
!...................................................................................................
       if((cref(1).lt.990.).and.(cref(2).gt.990.)
     &                    .and.(cref(3).gt.990.)) then
        coef(1) = cref(1)/dv3x3(1,1)
        coef(2) = 0.0
        coef(3) = 0.0
        print*,'  Exclusively HS calibration.'
       endif
!...................................................................................................
       if((cref(1).gt.990.).and.(cref(2).lt.990.)
     &                    .and.(cref(3).gt.990.)) then
        coef(1) = cref(2)/dv3x3(2,1)
        coef(2) = 0.0
        coef(3) = 0.0
        print*,'  Exclusively MS calibration.'
       endif
!...................................................................................................
       if((cref(1).gt.990.).and.(cref(2).gt.990.)
     &                    .and.(cref(3).lt.990.)) then
        coef(1) = cref(3)/dv3x3(3,1)
        coef(2) = 0.0
        coef(3) = 0.0
        print *,'  Exclusively LS calibration.'
       endif
!...................................................................................................
       if((cref(1).lt.990.).and.(cref(2).gt.990.)
     &                    .and.(cref(3).lt.990.)) then
        dv2x2(1,1) = dv3x3(1,1)                  ! subtracts off drift-corrected zeros
        dv2x2(2,1) = dv3x3(3,1)
        do i=1,2
         dv2x2(i,2) = dv2x2(i,1)*dv2x2(i,1)      ! voltage squared
         coef2x2(i) = 0.0
        enddo                              
        cref2x2(1) = cref(1)
        cref2x2(2) = cref(3)
        print *,'cref(2)',cref(2),'cref(3)',cref(3)
	call solvetwo(dv2x2,coef2x2,cref2x2)
        coef(1) = coef2x2(1)
        coef(2) = coef2x2(2)
        coef(3) = 0.0
        print *,'coef(1)',coef(1),'coef(2)',coef(2),'coef(3)',coef(3)
        print *,'  Eddy Exclusively HS and LS calibration.'
       endif
!...................................................................................................
       if((cref(1).gt.990.).and.(cref(2).lt.990.)
     &                    .and.(cref(3).lt.990.)) then
        dv2x2(1,1) = dv3x3(2,1)                  ! subtracts off drift-corrected zeros
        dv2x2(2,1) = dv3x3(3,1)
        do i=1,2
         dv2x2(i,2) = dv2x2(i,1)*dv2x2(i,1)      ! voltage squared
         coef2x2(i) = 0.0
        enddo
        cref2x2(1)=cref(2)
        cref2x2(2)=cref(3)
	call solvetwo(dv2x2, coef2x2, cref2x2)
        coef(1) = coef2x2(1)
        coef(2) = coef2x2(2)
        coef(3) = 0.0
        print *,'  Exclusively MS and LS calibration.'
       endif
!...................................................................................................
       if((cref(1).lt.990.).and.(cref(2).lt.990.)
     &                    .and.(cref(3).gt.990.)) then
        dv2x2(1,1) = dv3x3(1,1)                   ! subtract off drift-corrected zeros.
        dv2x2(2,1) = dv3x3(2,1)
	do i=1,2
	 dv2x2(i,2) = dv2x2(i,1)*dv2x2(i,1)       ! voltage squared.
     	 coef2x2(i) = 0
	enddo
        cref2x2(1) = cref(1)
        cref2x2(2) = cref(2)
	call solvetwo(dv2x2, coef2x2, cref2x2)
        coef(1) = coef2x2(1)
        coef(2) = coef2x2(2)
        coef(3) = 0.0
        print *,'  Exclusively HS and MS calibration.'
       endif
!...................................................................................................
       if((cref(1).gt.990.).and.(cref(2).gt.990.)
     &                    .and.(cref(3).gt.990.)) then
        print*, 'no CO2 cal tanks available reading values from ',dtype, 'calcoeff'
        call readco2cal(time0,dtype) 
        coef(1) = fcoef(1)
        coef(2) = fcoef(2)
        coef(3) = fcoef(3)
!        print*,'coef(1) at EB cal: ', coef(1)
        print *,'  No cal tanks -borrowed coefficients: ', coef
        flagtank = 1.
      endif
      else				     ! all tanks are in
!===========================================================================================
!        print *,'  CHECK calibration : dv3x3', dv3x3
       call solve(dv3x3, coef, cref)
!===========================================================================================
!     problem:  dv(3x3) * coef(1:3) = cref(1:3)
!     solution: coef = dv^-1 * cref
      endif				     ! endif when at least one cal tank is missing
!nrc end.............................................................................................
      endif				     ! end calflag(iecal).lt.3
      coef1(iecal) = coef(1)                 ! save coefficients for later use. 
      coef2(iecal) = coef(2)
      coef3(iecal) = coef(3)
      call gethhmmss(time1)
      write(99,510) iecal, time1, int(mod(time1,1d0)*86400),
     &              hh,mm,ss, coef, calflag(iecal)
510   format('  Calibration # ',i2,' hs at ', f9.4,' (',i5,' sec=',i3.2,
     &       ':',i2.2,':',f5.2, '), coefs: ',
     &       3(1x,1pe11.4), ' (calflag =',i3,')')

!...................................................................................................
!     4c. calculate the concentrations for dilution and air sample.
!...................................................................................................
!     izero+1 = zero after 3-span set = zero just before dilution
      co2dil = vcal(iecal,idil)-vzco2(izero+1)  
      h2odil = vcal(iecal,idilh)-vzh2o(izero+1)
      
      call getconc(h2odil,co2dil,pcellc(iecal,3),pcellc(iecal,1), 
     &                 tcellc(iecal,3),tcellc(iecal,1),coef) 

!     p,t index=1 <-> cal; index=3 <-> dilution (pdil should be higher b/c high-flow)

      co2air = vcal(iecal,iair) -vzco2(izero+1)
      h2oair = vcal(iecal,iairh)-vzh2o(izero+1)

!        print*,'coef at EB getconc: ', coef
!        print*,'co2air: ', co2air
!        print*,'h2oair: ', h2oair
!        print*,'tratio: ', (273.d0+tcellc(iecal,4))/(273.d0+t0h2o)
!        print*,'pcellc4: ', pcellc(iecal,4)
!        print*,'pcellc1: ', pcellc(iecal,1)
!        print*,'tcellc4: ', tcellc(iecal,4)
!        print*,'tcellc1: ', tcellc(iecal,1)
!        pcellc(iecal,1) = pcellc(iecal,4)
      print*,'before conc, h2oair: ', h2oair,' co2air: ', co2air,' pcellc(4):', pcellc(iecal,4), '  pcellc(1):', pcellc(iecal,1)
      print*,'tcellc:', tcellc(iecal,4),'     tcellc:', tcellc(iecal,1)
      call getconc(h2oair,co2air,pcellc(iecal,4),pcellc(iecal,1), 
     &                 tcellc(iecal,4),tcellc(iecal,1),coef)
      print*,'co2air after getconc: ', co2air,' vcal(iecal,iair):', vcal(iecal,iair) ,' -vzco2(izero+1)', -vzco2(izero+1)

!...................................................................................................
!     4d. write results.
!...................................................................................................
!    note:  wrong order of v0 (was 1,2,3 = hs,ms,ls) changed 23-oct-01
!...................................................................................................
      tdoy = tcalave(iecal) - days00 - yrdays1
      !nrc added set boundaries of what is possible given the write format constrains
      if ((dv3x3(3,1).lt.-9.9999).or.(dv3x3(3,1).gt.99.9999))   dv3x3(3,1)=-9.9999	!f7.4
      if ((dv3x3(2,1).lt.-9.9999).or.(dv3x3(2,1).gt.99.9999))   dv3x3(2,1)=-9.9999	!f7.4
      if ((dv3x3(1,1).lt.-9.9999).or.(dv3x3(1,1).gt.99.9999))   dv3x3(1,1)=-9.9999	!f7.4
      if ((co2dil.lt.-999.).or.(co2dil.gt.9999.))   co2dil=-999.99	!0PF7.2
      if ((co2air.lt.-999.).or.(co2air.gt.9999.))   co2air=-999.99	!0PF7.2
      if ((h2odil.lt.-99.).or.(h2odil.gt.999.))     h2odil=-99.999	!f7.3
      if ((h2oair.lt.-99.).or.(h2oair.gt.999.))     h2oair=-99.999	!f7.3
      do i=1,4 
       if ((tcellc(iecal,i).lt.-99.99).or.(tcellc(iecal,i).gt.99.99))  tcellc(iecal,i)=-99.99
      enddo
      !no more crazy dates
      if ((tdoy.ge.1.).and.(tdoy.lt.367.).and.(tcalave(iecal).ge.1.).and.(tcalave(iecal).lt.100000.)) then	
      !..................................................................................
       write(13,513) tcalave(iecal),tdoy,v0(3),v0(2),v0(1),
     &              dv3x3(3,1),dv3x3(2,1),dv3x3(1,1), calflag(iecal)  ! low,mid,high dv
       write(13,514) vzco2(izero+1),vzh2o(izero+1),vcal(iecal,iz2),
     &              vcal(iecal,iz2h),vcal(iecal,idil),vcal(iecal,iair),
     &              vcal(iecal,idilh),vcal(iecal,iairh)
       write(13,515) (pcellc(iecal,i),i=1,4),prcon(iecal), 
     &              (tcellc(iecal,i),i=1,4)
!     added temperature during cal that gave cal coefs.  Before only up to i=3 tcellc(iecal,i),i=1,3)
       write(13,516) (flsam(iecal,i),i=1,4), (flcal(iecal,i),i=1,4)
       write(13,517) coef,co2dil,co2air,h2odil,h2oair,                ! coef = array of 3 cal-tank concentrations.
     &              cref(1),cref(2),cref(3) 
513    format(f10.4,1x,f8.4,1x, 6(1x,f7.4), 1x, i4 )
514    format(1x,f7.4,1x,7(1x,f7.4) )
515    format(1x,f7.2,1x,4(1x,f7.2), 4(1x,f6.2))
516    format(1x, 8(1x,f7.1))
517    format(1x, 3(1x,1pe13.6),2(1x,0PF7.2), 2(1x,f7.3), 3(1x,f7.2))
      endif		                                              !endif no more crazy dates

!...................................................................................................
!     note:
!     12-13-01: changed from 1pe11.4 to 1pe13.6
!     initial coef format changed. it was: e11.5
!     note on possible formatting confusion for 216-format line:
!     the format-editor 1pe11.4 was used to output cal-coefficients
!     ->  1p, when used with e, the base is multiplied by 10, the exponent reduced by 1;
!         [ the default format for e-editing is 0.xyzeab, this produces x.yze(ab-1) ]
!         but the 1p-setting remains for the whole line, unless changed.  
!     ->  1p, when used with f, multiplies the output by 10
!         thus, without the 0p in front of f7.2, the other values would be multiplied by 10
!         (see p editing, watcom lang. ref. manual, p. 274)
!...................................................................................................
      enddo             ! for "iecal=1, nofecal"
      close(13)
      write(99,*) '  Calibration file writen on ', fname(1:lfn)//'.cal' 
      write(* ,*) '  Calibration file writen on ', fname(1:lfn)//'.cal' 

!---------------------------------------------------------------------------------------------------
!     fifth --  performs calculations for surveillance standard, write output file (E*.ss)
!               (interpolates zeros, interpolates calibration coeffs, calculates conccentrations)
!---------------------------------------------------------------------------------------------------
600   continue                      ! jumps here if done writing *.cal file.
      jlo=1
      klo=1

!     loops through each surveillance standard (probably only 1).
      do iss=1, nofss
      time0 = tss(iss)

!...................................................................................................
!     5a. zero interpolation:
!...................................................................................................
      izero = 0                             ! nrc 170207 added
      if (nofzero.gt.1) izero=ihuntf(tzero,nofzero,time0,jlo)

       if (izero.eq.0) then                 ! before the range of recorded zeros.
        v0h2o = vzh2o(1)
        v0co2 = vzco2(1)
        izero = 1
        jlo   = izero
        if(zerflag(izero).gt.0)  liflag(iss)=liflag(iss) + 10    ! 10 if zero bad, 0 otherwise
         print*, 'At the eddy zero --zeros were bad'  
!     setting lvlflag = x1x means:  unable to extract nearby zero.
!     (reminder:  zerflag =0 if zero ok, 1 otherwise)
!     lvlflag = 100 already if licor is off during the level measurement.

       elseif (izero .eq. nofzero) then     ! after the range of recorded zeros
        v0h2o  = vzh2o(nofzero)
        v0co2  = vzco2(nofzero)
        if(zerflag(izero).gt.0)  liflag(iss)=liflag(iss)+10       ! 10 if zero bad, 0 otherwise
       elseif ((izero.gt.0).and.(izero.lt.nofzero)) then         ! within range
        jlo    = izero
        deltat = tzero(izero+1)-tzero(izero)
       if(zerflag(izero).eq.0) then 

!     if previous zero good, check the next one
       if(zerflag(izero+1).eq.0) then 
        delt = time0-tzero(izero) 
       else                                  ! if the next zero is bad, doesn't interpolate.
        delt = 0                             ! setting delt=0 prevents interpolation below.
        liflag(iss) = liflag(iss) + 2        ! 2 = next zero bad
         print*, 'Initial zeros is OK following not, skiping interpolation'  
       endif
      print*, 'CHECK voltage zero, vzh2o(izero): ', vzh2o(izero)  
!     water zero
        deltav = vzh2o(izero+1)-vzh2o(izero)
        dvdt   = deltav/deltat
        v0h2o  = vzh2o(izero)+dvdt*delt

!     co2 zero
        deltav = vzco2(izero+1)-vzco2(izero)
        dvdt   = deltav/deltat
        v0co2  = vzco2(izero)+dvdt*delt
       elseif(zerflag(izero+1).eq.0) then 

!     if previous zero was bad, but the next one is good, use next.
        liflag(iss) = liflag(iss) + 1         ! 1=prev zero bad

!     water zero
        v0h2o = vzh2o(izero+1)

!     co2 zero
        v0co2 = vzco2(izero+1)
         print*, 'Eddy zero interpolation -previous zero was bad, but the next one is good'  
      else

!     when both are bad:
        liflag(iss) = liflag(iss) + 10         ! 10=both zeros bad
        v0h2o = 0                              ! will not be used
        v0co2 = 0
        print*, 'Eddy zero interpolation -both zeros are bad'  
       endif
      endif
      vss(iss)  = vss(iss) -v0co2            ! co2 zero is subtracted
      vssh(iss) = vssh(iss)-v0h2o            ! h2o

!...................................................................................................
!     5b. calibration coefficient interpolation:
!...................................................................................................
      iecal=0                                ! nrc 170207 added iecal will change is nofcal>1 in next line
      if (nofecal.gt.1) iecal = ihuntf(tcalave,nofecal,time0,klo)  ! index to last calibration before time0
!       print*, 'CHECK cal coefficient interpolation iecal: ', iecal, ' calflag(iecal): ', calflag(iecal)
  
      if (iecal.eq.0) then
       klo   = 1
       iecal = 1
      endif
      coef(1) = coef1(iecal)                 ! don't need to interpolate, b/c ss is part of this cal
      coef(2) = coef2(iecal)
      coef(3) = coef3(iecal)

      if(calflag(iecal).gt.0) liflag(iss) = liflag(iss) + 20

      call gethhmmss(time0) 
      write(99,607) iss, time0, int(mod(time0,1d0)*86400), hh,mm,ss
607   format('  Surveillance standard # ',i2,' at ', f9.4, ' ('
     &       ,i5,' sec = ',i3.2,':',i2.2,':',f5.2,')') 
      call gethhmmss(dble(tcalave(iecal)))
      write(99,608) iecal, tcalave(iecal), hh,mm,ss
608   format('  --> using coefficients from cal #', i3,' at t =', 
     &       f9.4, ' (',i3.2,':',i2.2,':',f5.2,')')

!...................................................................................................
!     5c. derives concentrations from calibration curve:
!     corrects for pressure and temperature deviations from calibration conditions:
!...................................................................................................
      pratio = pcellc(iecal,1)/pcells(iss)
      tratio = (tcells(iss)+273.d0)/(tcellc(iecal,1)+273.d0)
      xpco2  = vss(iss)                     ! zero was subtracted above
!     izero+1 = zero after 3-span set = zero just before dilution
      xph2o  = vssh(iss)                     ! zero was already subtracted 
!===========================================================================================
!     170907 nrc get the coef from the SS and the zero, such as only one tank running
!     the system missed cal tanks in 2015 only the SS was available
      xh2oss = xph2o*1.
      xco2ss = xpco2*1.
      call getconc(xph2o,xpco2,pcells(iss),pcellc(iecal,1),tcells(iss),tcellc(iecal,1),coef)

!     co2 from surface standard (forced) (nrc 180120)
        fdilut = 1.d0/(1.d0-(xph2o/1000.d0))
!        coef(1) = cref(1)/dv3x3(1,1)
        coefss(1) = csurv*.1/(xco2ss*pratio*tratio*fdilut)
        coefss(2) = 0.0
        coefss(3) = 0.0

!     co2 from surface standard (forced) (nrc 180120)
      call getconc(xh2oss,xco2ss,pcells(iss),pcellc(iecal,1),tcells(iss),tcellc(iecal,1),coefss) 
      print*, 'xpco2ss final concentration at survey standard: ', xpco2ss
!===========================================================================================
!...................................................................................................
!     5d. write the results ss:
!...................................................................................................
       tdoy = time0 - days00 - yrdays1
       if ((xpco2.lt.-999.99).or.(xpco2.gt.999.99))	xpco2=999.99	!f7.2 nrc 180111
       if ((xph2o.lt.-999.99).or.(xph2o.gt.999.99))	xph2o=999.99	!f7.2
       if ((csurv.lt.0.0).or.(csurv.gt.999.99))	        xph2o=999.99	!f7.2
       if ((coefss(1).lt.-99999.99).or.(coefss(1).gt.99999.99))	coefss(1)=999.99	!f7.2
       if ((xco2ss.lt.-999.99).or.(xco2ss.gt.999.99))	xco2ss=999.99	!f7.2

       if ((tdoy.ge.1.).and.(tdoy.lt.367.).and.(time0.ge.1.).and.(time0.lt.100000.)) then
        write(14,624) time0,tdoy,v0co2,v0h2o,vss(iss),vssh(iss),coef
        write(14,625) pcells(iss),pcons(iss),tcells(iss),flsams(iss),flcals(iss),xpco2,xph2o,liflag(iss),csurv,coefss(1),xco2ss
       endif
!...................................................................................................
!     note:
!     10-4-01:  re-ordered p,t and co2, h2o (was wrong-way round)
!...................................................................................................
624    format(f10.4,1x,f9.4,2(1x,f7.4),2(1x,f7.4),3(1x,1pe11.4))
625    format( 1x,3(1x,f7.2),2(1x,f7.1),2(1x,f7.2),1x,i4,1x,f7.2,1x,1pe11.4,1x,f7.2)
       enddo
       close(14)
       write(99,*) '  Surveillance standard file writen on  ',
     &             fname(1:lfn)//'.ss' 
       write(* ,*) '  Surveillance standard file writen on  ', 
     &             fname(1:lfn)//'.ss' 
       enddo                 ! from "do infile=1,nifiles "
       goto 990               ! exit this program.  

!     error processing:
900    continue

!     error processing code goes here:
       write(99,*) '  Exiting program LCALED due to error.'
       print *,    '  Exiting program LCALED due to error.'
990    continue
       write(99,*) '  End of program LCALED.'
       print *,    '  End of program LCALED.'
       close(99)
       end 
!===================================================================================================
!     end of program
!===================================================================================================


!===================================================================================================
!                                    subroutines and functions                                     !
!===================================================================================================
      include 'lsubs.for' 

!===================================================================================================
!                                       subroutine getcal                                          !
!===================================================================================================
!      given:
!      work(w,n) [common] = array of input voltages (w=1: co2, w=2: h2o)
!      nfirst = first n to start with for taking averages
!      time0 [common] = time (decimal julian days) of work(w,n=0) 
!      iecal [common] = index to vecal, tecal
!      ssflag = ssflag=1 if get surveil. std.  
!      return (all in common):
!      vcal(w,iecal) = summary mean voltages for cals 
!                   (1-3: high-low; 4-7: zero2, dil, ss, air; 8-11: h2o)
!      tcal(w, iecal) = times for each segment (as above)
!      tcalave(iecal) = mean time of cal-levels 
!                       and also return (if ssflag set)
!      vss(w,iss) = voltage (w=1-3: co2, p, t) during surv. std
!      tss(iss) = time of surv. std.
!--------------------------------------------------------------------------------------------------      
      subroutine getcal( nfirst )
      include 'lparams.for'
      parameter(nzdim=nezdim)
      include 'lco2vars7.for'
!     gets calibration (3 standard spans) :
      write(99,'(a14,i3)') '    --> cal # ',iecal
      do i=1,3                                             ! span-gas cals:  high, mid, low
       vco2 = 0
       ifirst = (nfirst+ndelayed) + (i-1)*nflstep
!      points from beg of 'work' to start of this calibration (including lag).
       do istore = ifirst, (ifirst+navecal-1)
        vco2 = vco2 + work(istore,1)
       enddo
       vcal(iecal,i) = vco2/float(navecal)
       tcal(iecal,i) = time0 + (ifirst/freqhi)/86400. 
!      time0 offset by how far into array (accounts for lag)
      enddo
      tcalave(iecal) = (tcal(iecal,1)+tcal(iecal,2)+
     &                 tcal(iecal,3) )/3                   ! mean cal time
!     gets zero:
      izero = izero+1
      call getzero( nfirst+3*nflstep ) 
!     stores post-main cal sequence:  
!     zero2 (at high-flow), zero2/ss, dilution, air sample (both co2 and h2o)

!...................................................................................................
!     note:
!     modified in 08-11-2008:  the calibrations are now 3 minutes long, except for the diluted zero,
!     which is only 2 minutes long.  the factor, if6, accounts for this.
!     note that an additional 60 sec of sample data is read to allow for the extra time.    
!...................................................................................................
      if6 = 0 
      do i=iz2,iair                    ! 4=zero2, 5=zero2/surv, 6=dilution, 7=air
       z2co2=0                         ! use for accumul: zero #2, surv, dilution for co2
       z2h2o=0                         ! zero #2, surv, dilution for h2o
       if (i.eq.6) if6 = -60  
        ifirst = (nfirst+ndelayed+if6) + i*nflstep     ! n.b.:  not (i-1) b/c zero is there
!       skips over low, mid, high, and zero to start of this calibration.
        do istore = ifirst, (ifirst+navecal-1)
         z2co2 = z2co2 + work(istore,1)
!        gets h2o
         z2h2o = z2h2o + work(istore,2)
        enddo
        vcal(iecal,i)   = z2co2/float(navecal)         ! i   = 4, 5,  6,  7  
        vcal(iecal,i+4) = z2h2o/float(navecal)         ! i+4 = 8, 9, 10, 11  (water)
        tcal(iecal,i)   = time0 + (ifirst/freqhi)/86400. 
      enddo
      if (ssflag.eq.1) then
       iss = iss + 1
       vss(iss)  = vcal(iecal,isurv)                   ! mean reading (isurv=5).
       vssh(iss) = vcal(iecal,isurv+4)                 ! water.
       tss(iss)  = tcal(iecal,isurv)
       ssflag = 0
      endif
      return
      end
!===================================================================================================
!                                   end of subroutine getcal                                       !
!===================================================================================================



!===================================================================================================
!                                      subroutine getzero                                          !
!===================================================================================================
!     given the array of input voltages "work", gets average of values ending with element "ilast".
!---------------------------------------------------------------------------------------------------
      subroutine getzero( nfirst )
      include 'lparams.for'
      parameter(nzdim=nezdim)
      include 'lco2vars7.for'
      zco2 = 0
      zh2o = 0
      ifirst = nfirst+ndelayed               ! lag (eddy) for solenoid switching/tube travel effect.
      do istore = ifirst, (ifirst+navecal-1) ! navecal= number of points in one calibration.
       zco2 = zco2 + work(istore,1)
       zh2o = zh2o + work(istore,2)
      enddo
      if (navecal.gt.0) then		     ! added by nrc 180111 to avoid inf
       vzco2(izero) = zco2/navecal
       vzh2o(izero) = zh2o/navecal
       tzero(izero) = time0 + (ifirst/freqhi)/86400.0
      endif
      return
      end
!===================================================================================================
!                                   end of subroutine getzero                                      !
!===================================================================================================

