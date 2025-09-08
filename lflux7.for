!===================================================================================================
!                                   lba code:  lfluxes7.for                                        !
!===================================================================================================
!     created by scott saleska (started january 2001)
!     (adapted from boreas code fluxes4 by s.-m. fan & doug sutton)
!---------------------------------------------------------------------------------------------------
!     calculates km67 tower fluxes
!     invoke:  lfluxes yymmdd -files ff -min mm 
! 
!     input data files:
!     .\yymmdd\split\e*yymmdd.200 = scaled eddy system data for id=200 (hi-freq data)
!     .\yymmdd\process\e*yymmdd.zer, and e*yymmdd.cal = zero and calibration files
!     output files:
!     .\yymmdd\process\e*yymmdd.flx = derived flux data for this eddy system
!
!     this code calculates eddy correlation fluxes:
!     for each eddy system, do:
!     1. loop:  read dt (default=1/2 hr) of data
!     2. apply calibration based on calibration files:
!          e*yymmdd.zer, and e*yymmdd.cal
!     3. calculate stress tensor for wind data:
!           a. subtract linear trend from raw data
!           b. assemble matrix from covariances
!     4. perform axis rotation on wind data
!     5. calculate heat (bouyancy) flux, <w't'>:
!           a. detrend t
!           b. calculate covariances (incl. lags to +/-4)
!     6. calculate flux correction factor to adjust for smearing of 
!         the co2, h2o signals:
!         flux correction = <t'w'>/< sm(t') sm(w') >,
!         where sm(x_t) = a* sm(x_t-1) + (1-a)x_t is a recursive 
!         numerical exponential smoothing function to simulate 
!         smearing caused by the tubing
!     7. estimate h2o, co2 fluxes 
!     8. write output to file for this eddy system
!         end loop
!     note: for each aggregation interval, cal=0 if no cal/zero in interval.
!           otherwise:              
!     calib:  1 = ends with zero, 2 = zero continuation, 3 = zero tail, 
!             4 = starts with new zero (shouldn't usually happen)
!            10 = full-cal sequence (covariance calculated, but user should ignore it
!     above+100 = licor was off at last zero (and may still be off now)
!---------------------------------------------------------------------------------------------------
      program lfluxes
      include 'lparams.for' 
      parameter(nzdim=nezdim)                   ! number of zeros for eddy zeroing.
      include 'lco2vars7.for'                   ! variables for doing co2 conccentration.
      real*8 dum4(4), dum5(5), dum6(6), dum8(8) ! dummy arrays.
      character*1 bo                            ! debug: indicates whether bounds checking.
      include 'lfluxvar.for'                    ! variables for flux calculation related to 
                                                ! filtering and chilled mirror.
      logical done, newyrflag, debug
      real*8  tdoy, yrdays0

!     water-vapor scale factor:  (added in 16-jul-2004)
      parameter(nwvsf = 50000)                  ! enough for 5+ years of hourly data.
      real*8 tscalewv(nwvsf)                    ! num water-vapor scale factors.
      real*8 scalewv(nwvsf)                     ! corresponding vector of hourly scale factors.
      real*8 tfind
      real*8 wvsf				

!     co2/h2o lag factors:  (added in 12-dec-2004)
      parameter(nlagtot = 50000)              ! enough for 5+ years of hourly data.
      real*8  tlagfac(nlagtot)                ! time of lag factors.
      integer co2co2(nlagtot)                 ! corresponding vector of co2 lags (number of
                                              ! high-frequency samples).
      integer co2h2o(nlagtot)                 ! corresponding vector of h2o lags .

!     sets  initial parameters, get command-line, etc.:
      lname = 'lflux6c '
      ftypes= '1'
      debug = .false.
!===================================================================================================
      call getfilterset(filtcode, filtset)  
!===================================================================================================
      include 'linitial7.for'                   ! common initialization code 
                                                ! sets directory structure (root paths, dirsplit, 
                                                ! dirprocin, dirprocout);
                                                ! retrieves command line, file names for processing
                                                ! [ to ifilen(1..nifiles) ]; 
                                                ! displays initial program message.

      if (navg.eq.0) navg = freqhi*60.*period   ! number of points to aggregate
      print *, '  Using filter code '//filtcode

!===================================================================================================
!                               main program loop starts here
!===================================================================================================

      do infile=1, nifiles                      ! loop over each input file
!---------------------------------------------------------------------------------------------------
!     first --  opens the input file (E*.200) and output files
!---------------------------------------------------------------------------------------------------
      read (yymmdd,'(i10)') idate		! nrc 180302 set calibration coefficient for humidity sensor
      dtype= ifilen(infile)(1:2)                ! dtype is data type string
      enum=1
!     parameter files:
      if (dtype.eq.'EB') then 
      enum=2
      end if  
      print *,'  Reading parameter file lsysdat.inp ...'
!===================================================================================================
      call readsysd(enum)   ! read data from lsysdata.inp for eddy system enum
!===================================================================================================
      call readconv(enum)   ! read data from lconv.inp for conversion coefs
!===================================================================================================
!     gets this loop's filename:
!     raw data input (*.200,201 data)
      fnamein = ifilen(infile)
      call getfullpath(rootsplit,yymmdd, dirsplit, fnamein, lfnin) 

!     calibration data input (processed *.cal, *.zer files)
      fname = ifilen(infile)
!===================================================================================================
      call getfullpath( rootprocin, yymmdd, dirprocin, fname, lfn)
!===================================================================================================
!     data output (*.flx, *.erf files)
      fnameout = ifilen(infile)
!===================================================================================================
      call getfullpath(rootprocout,yymmdd,dirprocout,fnameout,lfnout)
!===================================================================================================
      print *, '  Opening input files "',fnamein(1:lfnin)//'" and '
      print *, '                      "',fname(1:lfn)//'".'
      print *, '  Sending output to   "',fnameout(1:lfnout)

!     set-up error/status file:
      open(99,file=fnameout(1:lfnout)//'.erf',status='unknown')             ! error file
      write(99,55) platform, lname, iyear, imonth, iday,ihrs,imin,isec
55    format(a6,' ',a7,':  on date (yy-mm-dd) ', i4,'-',i2.2,'-',i2.2,
     &       ' (',i2.2,':',i2.2,':',i2.2,'), process:' )
      write(99,*) '  '//fnamein(1:lfnin)//'.200      using filter code '
     &                //filtcode
      print 8
      write(99,8)
8     format( '  Retrieving housekeeping (& chilled mirror) data from:')
      print *,    '  ', fname(1:lfn)//'.201'
      write(99,*) '  ', fname(1:lfn)//'.201'
      print 9, period, navg, dirprocout
      print *
      write(99,9) period, navg, dirprocout
9     format('  Aggregation period: ', f5.2, ' minutes (',
     &       i6, ' datapoints )',/, '  output files in directory: ', a8)
      open(11,file=fname(1:lfn)//'.zer',status='old')  
      read(11,*)                                                            ! header info
      read(11,*)
      print *,   '  Loading zero-file (*.zer) data...'
      write(99,*)'  Loading zero-file data from '//fname(1:lfn)//'.zer'

!     load zero data:
      i = 1
      nofzero=0
      nbad=0
      do while(i.le.nzdim)
110    continue
       read(11,111,end=112) tzero(i), tdoy, vzco2(i), vzh2o(i), 
     &                     pcellz(i), tcellz(i), prconz(i), 
     &                     flsamz(i), flcalz(i), zerflag(i)
111    format( f10.4,1x,f9.4,2(1x,f8.4),1x,f6.1,1x,f7.2,1x,
     &        f7.2,2(1x,f7.1),1x,i3)

       if (zerflag(i). gt. 0) then 
        nbad = nbad+1
       end if
       nofzero=i
       i = i+1
      end do
112   continue                                                  ! reached end of zero file.
      close(11)
      print *,    '  Read in', nofzero,' zeros (including ',nbad,
     &            ' bad zeros).'
      write(99,*) '  Read in', nofzero,' zeros (including ',nbad,
     &            ' bad zeros).'

!     last co2 time, zero:
!     calibration file
      open(12,file=fname(1:lfn)//'.cal',status='old') 
113   continue  
      read(12,'(i2,i3)') nhead, ncol
      do i=2,nhead
      read(12,*)                                                ! reads number of lines in header.
      end do
      print *,   '  Loading cal-file (*.cal) data.'
      write(99,*)'  Loading cal-file data from '//fname(1:lfn)//'.cal'
      i = 1
      nofcal=0
      nbad=0
      do while(i.le.ncdim)                       ! should be more than needed for 7*4=28 cals/week.

120   continue
      read(12,513,end=132) tcalave(i),tdoy, dum6, calflag(i)     
      read(12,514,end=132) vzco2(i),vzh2o(i)
      read(12,515,end=132) (pcellc(i,j),j=1,4),prcon(i), 
     &                     (tcellc(i,j),j=1,3)

!     p and t for j= 1: cal, 2: zero2, 3: dilution, 4: air-sample 
!     include j=4 (air sample) for p but didn't bother for temp (temperature doesn't change).

      read(12,516,end=132) (flsam(i,j),j=1,4), (flcal(i,j),j=1,4)
      read(12,517,end=132) coef1(i), coef2(i), coef3(i)

513   format(f10.4,1x,f8.4,1x, 6(1x,f7.4), 1x, i4 )
514   format(1x,f7.4,1x,7(1x,f7.4) )
515   format(1x,f7.2,1x,4(1x,f7.2), 3(1x,f6.2))
516   format(1x, 8(1x,f7.1))
517   format(1x, 3(1x,1pe13.6),2(1x,0PF7.2), 2(1x,f7.3), 3(1x,f7.2))

      if (calflag(i) .ge. 3) then 
      nbad = nbad+1
      end if
      nofcal=i
      i = i + 1
      end do
132   continue                                                 ! end of calibration file.
      close(12)
      print *,   '  Read in',nofcal,' calibrations (including ',nbad,
     &           ' bad calibrations).'
      write(99,*)'  Read in',nofcal,' calibrations (including ',nbad,
     &           ' bad calibrations).'

!     200 data file:
      open(20,file=fnamein(1:lfnin)//'.200',status='old')
      read(20, *)                                              ! header line 1   
      read(20, *)                                              ! header line 2
      write(99,*) '  Opened input file:  '//fnamein(1:lfnin)//'.200.'

!     hkp and chilled mirror data is in the .201
      open(88,file=fnamein(1:lfnin)//'.201',status='old')
      call skipmirhead(88)                                     ! skips header of eddy 201 file.
      write(99,*) '  Opened input file:  '//fnamein(1:lfnin)//'.201.'

!     allowing for pseudofilter 1 to open:
      if (filtset(5) .eq. 1) then
      call open200file(55,yymmdd,enum)     
      call skip200header(55)
      end if

!     set-up output files        
!     main output files
      open(19, file=fnameout(1:lfnout)//'.flx',status='unknown')
      write(99,*) '  Sending output to:  '//fnameout(1:lfnout)//
     &            '.flx, .erf'
      write(99,*)
      
!     headers for output file:   E*.flx
      write(19,140) enum, yymmdd, period, filtcode, platform, lname, 
     &              iyear,imonth,iday, ihrs,imin,isec
!     where 14 is number of rows and 76 is the number of variables
140   format('14 76  LBA Eddy',i1,': ',a6, ' (', f6.3, ' min), filt=', 
     &        a5, ' (',a6,' ',a8,' on: ',i4, '-', i2.2, '-', i2.2, ', ',
     &        i2.2,':',i2.2,':',i2.2,')' )
      write(19,*) 'jdstart.gmt  doy   cal    wspd.m.s    wdson.deg    hratio      fratio    '
      write(19,*) 'suu   suv    suw   svv    sww '
      write(19,*) 'tu    hu     cu    tv     hv          cv  '
      write(19,*) 'tw    hw     cw    phi    theta       beta '
      write(19,*) 'twlag_bar    hwlag_bar    cwlag_bar   twlag_max    hwlag_max   cwlag_max    '
      write(19,*) 'twlag_sd hwlag_sd    cwlag_sd   maxlagtw       maxlaghw        maxlagcw    '
      write(19,*) 't.c      h2o.mmol.m  co2.ppm    st    sh2o     sco2    '
      write(19,*) 'u_int    v_int    w_int     t_int     h_int    c_int   '
      write(19,*) 'u_slp    v_slp    w_slp     t_slp     h_slp    c_slp   '
      write(19,*) 'tamb.c   rh.tdew  plicor.t  tlicor.c    '
      write(19,*) 'h2o.v    co2.v    pres.tor  temp.c    '           
      write(19,*) 'u        v        w         u_std	 v_std	  w_std    '           
      write(19,*) 'n.99     ndsos    nnolock   nsiglow   nsighi   nlirl    notused.percent    n'           
      
!---------------------------------------------------------------------------------------------------
!     second - loop through input files to get & store data segments:
!---------------------------------------------------------------------------------------------------
!
!     a. initialize to read this file
      time0 = 0
!===================================================================================================
      call read20n(20, 200, dtype, time0, v200, 1, 7, is200, 1, 2, done, 
     &             yrdays, debug)
!===================================================================================================
!     input args :    fileptr (20), line-num (200), dtype (EA/EB)
!     output args:    time0, v200(1..7), is200(1,2), done (flag), debug reset to false
!     (time0 = decimal days since 1/1/2000)
!     also checks year-boundary crossings:  
!     yrdays set to number of days of year if year-boundary encountered
      yrdays0 = yrdays  
      if (done) goto 290                                  ! if at end of file.
      iscal= is200(2)                                     ! cal = 00 means normal sampling.
!===================================================================================================
      call gethhmmss(time0)
!===================================================================================================
      write(99,202) time0, int(gmt), hh,mm,ss
202   format( '  First GMT time in 200 file is ' , f8.4 , ' (' , i5 , 
     &        ' sec =' , i3.2 , ':' , i2.2 , ':' , f5.2 , ')' )
      call startime(time0, tstart1, period)

      call gethhmmss(tstart1)
      write(99,203) tstart1, hh,mm,ss
      write(*,203) tstart1, hh,mm,ss
203   format( '  Start-time is              ',f10.4,' (           ',
     &       i3.2 , ':' , i2.2 , ':' , f5.2 , ')' )

!...................................................................................................
!     note:  when ready to load in more extensive zero/cal(to prevent problems at ends)
!     can move zero/cal loading routine here (already done in lflux4r code)
!...................................................................................................

      write(99,204) (tzero(i), i=1,nofzero)
      write(99,205) (tcalave(i), i=1,nofcal)
204   format('  Zero times: ', 11f10.3)        
205   format('  Calibration times: ', 3f10.3)        
      write(99,*)

!     gets water-cal scale factor
      lrt = lentrim0(rootinp, len(rootinp))
      open(12, file=rootinp(1:lrt)//slash//dtype//'h2ocal.inp',
     &     status='old')                                       ! wv scale factors created in splus.

!     re-use file 12 (was used for co2 cal file)
      read(12,*)        
      print *,   '  Loading water-vapor scale factors.'

      write(99,*)'  Loading cal-file data from '//dtype//'h2ocal.inp.'
      i = 1
      tscalewv(i) = 0 
      nscale=0
      nbad=0

!     moves to desired start time in wv.cal file (quarter day before start of eddy data), acquires 
!     initial wv.cal(i=1):
      do while(tscalewv(i).lt.(time0-.25))
      read(12,*, end=133) tscalewv(i), scalewv(i)
      end do
      nscale=i                   ! i.e., 1
      i = 2                      ! gets set to read cal 2  
      do while(i.le.nwvsf)       ! should be more than needed for 7*4=28 cals/week

!     use list-directed input for universal cal file:
      read(12,*,end=133) tscalewv(i), scalewv(i)
      nscale=i
      i = i + 1
      end do       ! ending "do while (i.le.ncdim)".
133   continue     ! end of calibration file.
      close(12)
      print *,   '  Read in',nscale,' wv scale factors.'
      write(99,*)'  Read in',nscale,' wv scale factors.'
!     finished loading wv scale-factors 

!---------------------------------------------------------------------------------------------------
!     gets lag-correlation factors:
!---------------------------------------------------------------------------------------------------
      lrt = lentrim0(rootinp, len(rootinp))
      open(12, file=rootinp(1:lrt)//slash//dtype//'lags.inp',
     &     status='old')                                  ! lag factors created from llagcor.
!     re-use file 12 
      read(12,*)                                          ! header line
      print *,   '  Loading lag factors.'
      write(99,*)'  Loading lag data from '//dtype//'lags.inp.' 
      i = 1
      tlagfac(i) = 0 
      nlagfac = 0
      nbad    = 0
!     moves to desired start time in lag file (one day before start of eddy data), and acquires 
!     initial lag.
      do while(tlagfac(i).lt.(time0-1.0))
      read(12,*, end=134) tlagfac(i),iyear,idoy,co2co2(i),co2h2o(i)
!      print *,   '  i:',i,' tlagfac(i):', tlagfac(i),' iyear:',iyear,
!     &     ' co2co2(i): ', co2co2(i),' co2h2o(i): ', co2h2o(i)         !nrc 170517
      end do
      nlagfac = i                 ! i.e., 1
      i = 2                       ! gets set to read cal 2.
      write(99,*)'  Starting to read lags for day ',idoy,' of year ',
     &            iyear
      do while(i.le.nlagtot)      ! should be more than needed for 7*4=28 cals/week.

!     uses list-directed input for universal calibration file:
      read(12,*,end=134) tlagfac(i),iyear,idoy,co2co2(i),co2h2o(i)
      nlagfac=i
      if (co2co2(i).le.-999) co2co2(i) = nlagco2  ! default value from lsysdat.inp
      if (co2h2o(i).le.-999) co2h2o(i) = nlagh2o  ! ditto
      i = i + 1
      end do           ! ending "do while (i.le.ncdim)"
134   continue         ! end of calibration file.
      close(12)
      write(99,*)'  Finished reading lags with day ',idoy,' of year ',
     &            iyear
      print *,   '  Read in',nlagfac,' lags for co2 and h2o.'
      write(99,*)'  Read in',nlagfac,' lags for co2 and h2o.'
!     finished lag-facs 

!     move to desired start time in each file, acquire initial values:
!     first, save tovlap of data just *before* start time
!     (in case first full segment has zero at end)
      write(99,*)'  Storing an initial set of data points to ',
     &           'depository: '
      if (time0.lt.(tstart1-tovlap)) then     ! acquire full segment
      do while(time0.lt.(tstart1-tovlap)) 
!===================================================================================================
      call read20n(20,200,dtype,time0,v200,1,7,is200,1,2,done,
     &             yrdays,debug)
!===================================================================================================
!     input arguments :    fileptr (20), line-num (200), dtype (EA/EB)
!     output arguments:    time0, v200(1..7), is200(1,2), done (flag), debug reset to false
      if (done) goto 291         ! end-of-file here, then we are done before we started.
      end do
      else
      write(99,*)'  (Note: Less than 02:50 min of data to store before',
     &           ' start time!)'
      end if            ! if not enough data for full tovlap segment, just start with where we are.
      n=0
      time1 = time0     ! (to record start time of initial interval).
      do while(time0.lt.tstart1)
      n=n+1             ! should accumulate to 3 min worth this initial segment.
      co2(n)  = v200(1)
      h2o(n)  = v200(2)
      u(n)    = v200(3)
      v(n)    = v200(4)
      w(n)    = v200(5)
      t(n)    = v200(6)
      stat(n) = v200(7)
!===================================================================================================
      call read20n(20, 200, dtype, time0, v200, 1, 7, is200, 1, 2, done,
     &             yrdays, debug)
!===================================================================================================
!     input arguments :    fileptr (20), linenum (200), dtype (EA/EB)
!     output arguments:    time0, v200(1..7), is200(1,2), done (flag), debug resets to false.
      if (done) goto 300
      end do        
!     stores n values at end of store(1..novlap,.) array
!     (fills values from bottom-up, in case n <> novlap)
!===================================================================================================
      call gethhmmss(time1)
!===================================================================================================
      write(99,207) time1, hh,mm,ss, n, n/freqhi
207   format('  Starting at ', f10.4 , '(',i2.2,':',i2.2,':',f5.2,
     &       '), storing ', i5,' points (',f7.3,' secs)')
      do i = novla p, (novlap-n+1), -1

!     if-statement fix added 19-may-02:
      if (i.le.0) then                ! check to avoid accessing store zero-element (if n>novlap).
      write(99,*) '  Attempted to store too many initial points before',
     &            ' start !'
      write(99,*) '  (of n=',n,' loaded, can store only last novlap = ',
     &            novlap,')'
      write(99,*)                    ! when this happens, its usually just one point.
      goto 135    

!     skips trying to go back further in time, we have a full novlap.
      end if
      j = i - (novlap-n) 
      store(i,1) = u(j)
      store(i,2) = v(j)
      store(i,3) = w(j)
      store(i,4) = t(j)
      store(i,5) = co2(j)
      store(i,6) = h2o(j)
      store(i,7) = stat(j)
      end do
135   continue
      tnext    = tstart1
      numint   = 0          ! number of aggregation intervals processed.
      calib    = 0          ! whether this interval starts with a cal/zero.
      endzer   = 0          ! whether this interval ends with a zero.
      nseg     = 0 
      calbit   = .false.
      unexpcal = .false.    ! unexpected calibration.
      jlo = 1               ! initialize search variables for licor zeros 
      klo = 1               ! and for calibrations.
      mlo = 1               ! water-vapor scale factor.
      nlo = 1               ! lag factor.

!     b. file-loop to go through this file.
!===================================================================================================
210   continue             ! loop until this input file is exhausted;
                           ! goes on until end-of-file causes transfer out ;
                           ! at end, goto label 301,302,303 (depending on which file ends first)

!     data aggregating, by time:
!     tstart= starting time for this interval, tnext= start for next interval.
!     loop to read data until time of data-line >= tnext
!     zero/cal handling procedure:
!     zero (2-mins):  grab 2-min from next 1/2 hour by adding tskip (=2min) to tnext

!     (1) initialize data arrays
      tstart=tnext                                    ! start time for this interval.
      tnext = tnext + period/1440.d0                  ! minutes to decimal julian days.
      tmid  = tstart + (period/2)/1440.d0             ! mid-point of this period.    
!===================================================================================================
      call gethhmmss(tstart)
!===================================================================================================
      write(99,220) numint+1, tstart,hh,mm,ss
220   format('  Aggregation inteval: ',i4 / '  nominal start-time: ',
     &        f10.4,' (',i2.2,':',i2.2,':',f5.2,')')
!===================================================================================================
      call gethhmmss(tnext)
!===================================================================================================
      write(99,221) tnext,hh,mm,ss
221   format('  End-time :  ',f10.4,' (',i2.2,':',i2.2,':',f5.2,')')

!     (2) loop to read and store one aggregation interval:
10    continue
      n=0
      time1 = time0           ! start time of this segment (could be tstart+tovlap).
!===================================================================================================
      call gethhmmss(time0)
!===================================================================================================
      write(99,230) time0, hh,mm,ss, calib, nseg
230   format('  Starting to read file at:   ',f10.4,
     &       ' (',i2.2,':',i2.2,':',f5.2,'); ',
     &       ' (calib = ',i4,'; nseg = ', i1,')' )

      do while( time0.lt.tstart+tovlap*(nseg+1) )
      n = n+1                 ! should accumulate to 3 min to worth this initial segment.
      co2(n)  = v200(1)
      h2o(n)  = v200(2)
      u(n)    = v200(3)
      v(n)    = v200(4)
      w(n)    = v200(5)
      t(n)    = v200(6)
      stat(n) = v200(7)
      cal(n)  = is200(2)
      if ( (cal(n).ne.0) .and. (.not.calbit)) then
      calbit=.true.           ! note if calibration is in initial segment .

      time1 = time0
!===================================================================================================
      call gethhmmss(time1)
!===================================================================================================
      write(99,235) time1,hh,mm,ss
235   format('  Calibration/zero sequence at start: ',f10.4,
     &       ' (',i2.2,':',i2.2,':',f5.2 ,').')
      end if
!===================================================================================================
      call read20n(20, 200, dtype, time0, v200, 1, 7, is200, 1, 2, done,
     &             yrdays, debug)
!===================================================================================================
!     input args :    fileptr (20), line-num (200), dtype (EA/EB)
!     output args:    time0, v200(1..7), is200(1,2), done (flag), debug reset to false
      if (done) goto 300
      end do
      call gethhmmss(time0)

237   format('  Stopped reading file at:  ',f10.4,
     &       ' (',i2.2,':',i2.2,':',f5.2,'), cum. points= ',i5,'.') 

!     (b) check if this initial segment contains a calibration
      if (calbit) then

!     calibration bit was set in data, so skip over this segment 
      calbit = .false.
      if (calib.eq.0) then  ! no calibration set yet (first time through)
      calib = 4             ! calib = 4 means that this half-hour starts with new zero.
!...................................................................................................
!     note: calib:  1 = ends with zero, 2 = zero continuation, 3 = zero tail, 
!            4 = starts with new zero (shouldn't usually happen)
!            100 = full-cal sequence (changed from 10 on 9-oct-01)
!...................................................................................................
      if (endzer.eq.1) then 
      endzer = 0
      calib  = 2  ! calib=2:  this half-hr continues a zero
      end if
      nseg = 1
      goto 10     ! since first seg was for cal/zero, go read the next 3-min segment
      else        ! if calib already set, then cal bit set in 2nd 3 min, so must be full cal
      calib = 100 ! calib=100:  this half-hr starts with a full cal
      nseg = 0    ! don't need to read past end of this 1/2 hour any more
      end if
      elseif (endzer.eq.1) then ! if calbit not set ...
!     if last segment ended with zero that doesn't carry into this,
!     then force a delay to account for re-equilibration time in tube after zero.
      endzer = 0  ! reset endzer because it has served its purpose
      calib  = 3  ! calib:  3 = tail on zero that started at end of last segment
      nseg   = 1  ! number of segments to carry-over
      goto 10     ! go back to read next segment
      end if

!     (c) read remainder of this segment
15    continue
!===================================================================================================
      call gethhmmss(time0)
!===================================================================================================
      write(99,240) time0,hh,mm,ss, n+1
240   format('  Continuing to read file at: ',f10.4,
     &       ' (',i2.2,':',i2.2,':',f5.2,'), data-point:  ',i5,'.')
      do while(time0.lt.(tnext+tovlap*nseg) )
       n = n+1                   ! should accumulate to navg
       co2(n) = v200(1)
       h2o(n) = v200(2)
       u(n) = v200(3)
       v(n) = v200(4)
       w(n) = v200(5)
       t(n) = v200(6)
       stat(n) = v200(7)
       cal(n)  = is200(2)
      if ((cal(n).ne.0) .and. (calib.ne.100) .and. (.not.calbit) 
     &                 .and. (.not.unexpcal)) then
!     if cal in remaining segment (presumably at end)
!     unless calib=100 (if cal is near beg then is full cal-sequence, so calib should =100)
!     unless calbit not set (first time noticed, since after first time, calbit=true)
      time1 = time0 ! float(julday) + gmt/86400.d0 ! get time at front of cal-segment
!===================================================================================================
      call gethhmmss(time1)
!===================================================================================================
      if ((time1 + 1.5*tovlap) .lt. tnext) then   ! if this cal bit is not at end i.e. if it is more
                                                  ! than 1.5 overlap intervals before the end, then:
      unexpcal=.true.               
      write(99,241) time1, hh,mm,ss
241   format('  Unexpected calibration found at: ',
     &        f10.4,'(',i2.2,':',i2.2,':',f5.2 ,').')
      write(99,*) '  (setting calib=100, nseg=0 for this interval)'
      nseg   = 0
      calib  = 100                
      else                                ! if it is at the end, it is an anticipated zero-sequence.
      calbit = .true.   
      write(99,242) time1, hh,mm,ss
242   format('  Zero-sequence at end: ',f8.4,
     &       ' (',i2.2,':',i2.2,':',f5.2 ,').')
      end if
      end if
!===================================================================================================
      call read20n(20, 200, dtype, time0, v200, 1, 7, is200, 1, 2, done,
     &             yrdays, debug)
!===================================================================================================
!     input arguments :    fileptr (20), linenum(200), dtype (EA/EB)
!     output arguments:    time0, v200(1..7), is200(1,2), done (flag), debug reset to false.
      if (done) goto 300
      end do               ! add tovlap segment if was zero at beginning.
!                          (nseg=1:  add one additional segment).
      call gethhmmss(time0)
      write(99,237) time0, hh,mm,ss,n   ! when stopped reading file at end of interval
      if (unexpcal) then
      unexpcal = .false.
      end if

!    (d) tests if zero at end (if cal/zero at end, it is 'zero', not 'calibration')
      if (.not.calbit) then             ! no zero at the end of this segment.

!     saves last novlap points
      time2 = time0-1/(freqhi*86400.d0) ! true end: 2nd-last line just read.
      j = n - novlap                    ! no zero at the end, store last tovlap points.
      if (n.ge.novlap) then 
      do i = 1, novlap
       j = j + 1
       store(i,1) = u(j)
       store(i,2) = v(j)
       store(i,3) = w(j)
       store(i,4) = t(j)
       store(i,5) = co2(j)
       store(i,6) = h2o(j)
       store(i,7) = stat(j)
      end do
      write(99,244) n-novlap+1, n, tovlap*86400.
244   format('  Storing points ',i5,' - ',i5,' (',f8.4,
     &       ' seconds) in data depository...')
      else
      write(99,*) '  Warning: insufficient points to store in data',
     &            ' depository.'
      end if
      else            ! there is a zero starting at end of this segment.

!     move data in arrays down, restore from store to first part of array.
      calbit = .false.
      endzer = 1
      calib  = 1                       ! calib=1 --> zero at end of segment (1= end, 10=full cal)
      if (n.ge.novlap) then
      time1  = tstart - tovlap          ! true start is now earlier by tovlap
      time2  = tnext - tovlap           ! true end is earlier by same amount.
      write(99,245) novlap, tovlap*86400.
245   format('  Moving points in data array down by ',i5,
     &       ' (',f7.3,' secs)...')
      do i   = 0, n-novlap-1             ! moves existing data down over calibration part.
       j = n-novlap-i
       u(n-i) = u(j)
       v(n-i) = v(j)
       w(n-i) = w(j)
       t(n-i) = t(j)
       co2(n-i)  = co2(j)
       h2o(n-i)  = h2o(j)
       stat(n-i) = stat(j)
      end do

!     restores last segment of data from depository 
      write(99,247) novlap, tovlap*86400.
247   format('  Restoring ',i5,' points (',f8.4,
     &       ' secs) from depository to first part of array...')
      do i = 1, novlap                  ! (from end of last 1/2 hr to beginning of this 1/2 hr).
       u(i) = store(i,1)
       v(i) = store(i,2)
       w(i) = store(i,3)
       t(i) = store(i,4)
       co2(i)  = store(i,5)
       h2o(i)  = store(i,6)
       stat(i) = store(i,7)
      end do
      end if
      end if
!
      if (calib.gt.0 .or. endzer.gt.0) then
      write(99,255) calib, endzer
!===================================================================================================
      call gethhmmss(time1)
!===================================================================================================
      write(99,256) time1, hh,mm,ss
!===================================================================================================
      call gethhmmss(time2)
!===================================================================================================
      write(99,257) time2, hh,mm,ss
255   format('  Note:  cal/zero in this interval (calib='
     &        , i4,'; endzer=',i2,' )' )
256   format('  True start-time :   ',f10.4,
     &       ' (',i2.2,':',i2.2,':',f5.2,')')
257   format('  true end-time   :   ',f10.4,
     &       ' (',i2.2,':',i2.2,':',f5.2,')')
      end if
      if (n.le.navg/2 .or. n .ge. navg*2) then 
      write(99,259) n, navg
259   format('  Warning: bad number of values read:  n= ',i5,
     &       ' (expected navg= ',i5,').')
      write(99,*)'  Skipping this interval...'
      write(99,*)
      goto 160                              ! skips to end of loop for processing this interval.
      end if
      write(99,270) n, 200
270   format('  Aggregating: ',i6,' values from ',i3, ' file...')
      write(99,*)

!---------------------------------------------------------------------------------------------------
!     third -- calculates the flux
!---------------------------------------------------------------------------------------------------

!     step 0:  check sonic status bits (and licor ready-light)
!     status bit additions
!     determines the number of points within a half hour where the diagnosis status bits 
!     (bits 12-15) are flagged.

!     initialization
      do k=1,4
      diagarray(k) = 0     ! status bits 15..12
      numbitarray(k) = 0   ! counter for each status bit
      end do

!     loop for status variable to determine the diagarray for each status word
      do i=1,n   ! for each fast data point in the interval:
!===================================================================================================
      call statusbitarray(int(stat(i)))  ! unpack sonic status word to diagarray
!===================================================================================================
!     diagarray(1...4) <--> status bits 15...12
      filtvec(i) = 0                        ! default = keep this datapoint (inserted 8-may-02).
      do k=1,4                              ! checks each status bit of this point.
      if (diagarray(k).eq.1) then
      numbitarray(k) = numbitarray(k)+1     ! increments counts for each bit.
      end if

!     within the inside loop, set the filtering vector flags:               
      if (filtset(k).eq.1 .and. filtset(5).ne.1) then
!     filtset was set by 'getfilterset' from command line and it is "1" for each bit that the user 
!     may want to filter on.
      if (diagarray(k).eq.1) then
      filtvec(i) = 1                        ! filter out this data-point.

!...................................................................................................
!     note:  if we were filtering on more than one flag, subsequent flags not set would negate
!     the filtering of a previous flag that was set.
!...................................................................................................

      end if
      end if     
      end do
      end do
      if (filtset(5).gt.0) then              
!===================================================================================================
      call pseudofilter(filtvec, filtset(5),yymmdd,55,n,tstart)
!===================================================================================================
      end if      

!---------------------------------------------------------------------------------------------------
!     Chilled mirror data and licor status light:
!     Part of the code to calculate chilled mirror variables has been moved forward because 
!     licor status light filter is an additional filter added after the status bit filter 
!     was written, and because the code already written for it serves more than one 
!     function, I have 'patched' in the licor filter, rather than incorporating it into the 
!     filtering in the section just above.  (srs)
!     pPatch uses licoroff counter and licorstatus variables.
!
!     relocated code--in order to define the filtering before calling the detrend subroutine.
!
!     step 0b:  get chilled mirror and housekeeping inputs from *.201 file
!                 -added starting july 2001 (chilled mirror) and 
!                  sept '01 (housekeeping data)
!                 -relocated sept '01 {dmm}
!
!     chilled mirror additions
!     ------------------------
!     chilled mirror segment does four operations
!     a) matches file pointer in .201 time to time of the loop
!     b) reads in 1 period worth of points (e.g. 900 points,for 30 min at 1/2 hz)
!     c) converts to temperature units
!     d) averages
!
!     conversion (operation "c") requires use of the common block conversion array calculated in 
!     subroutine getmirfacs() called outside the main program loop. critical to the workings of
!     these operations is the flow of the cursor through the .201 file. 
!     it is admittedly hard to follow and not entirely robust--thus the user should be wary of 
!     reading errors within this segment, as they can cause the whole loop to terminate.
!
!     for simplicity, code initially lets eddy 1 and eddy 2 both use conversion factor for eddy 1. 
!     if the conversion factors become different for some reason this code will need to be modified. 
!     variables have been declared in lfluxvars.for. subroutines in lfluxsub.for together with 
!     filter routines.
!
!     nrc 2017 April 2017
!     chilled mirror was "replaced" by an HC2-S3L humidity and temperature sensor
!     nice and simple 
!---------------------------------------------------------------------------------------------------

!     (a) Line up dates

      n1 = n1PF-2                ! number of fields in eddy 201 file is 10. 2 less than in profile.
      ns = 9

!     Initializes:
      kavg = freqlo*period*60.   ! expected number of values to averge
      k = 1                      ! k = index to position within interval (1/2 hr -> 1..kavg points)
      do i = 1,n1                ! n1 = number of input variables
       xo(i) = 0.0d0             ! xo(i) = ith output
      end do
      do i = 1,n
       licorrlvec(i) = 0
      end do
      licoroff = 0         
      do i = 0,3                 ! calibration status:
       statcal(i) = 0.0          ! 0 = frequency of 00 (normal sample).
      end do                     ! 1 = frequency of 01 (normal zero/cal).
                                 ! 2 = frequency of 10 (zero-dilution).
                                 ! 3 = frequency of 11 (hi-flow zero).
      statcm=0.0                 ! frequency of chilled mirror status flag.

!     (b) Reads in one period worth of data:

!     Matches time:
!===================================================================================================
      call read20n(88, 201, dtype, time1, v201ed, 1, 10, is201, 1, ns,
     &             done, yrdays, debug)
!===================================================================================================
      do while(time1.lt.tstart)
!===================================================================================================
      call read20n(88, 201, dtype, time1, v201ed, 1, 10, is201, 1, ns,
     &             done, yrdays, debug)
!===================================================================================================
      if (done) goto 899
!     input arguments :    fileptr (88), linesum (201), dtype (EA/EB)
!     output arguments:    time1, v201(1..10), is201(1..ns), done (flag), debug reset to false
      end do
!     v201ed used instead of v, starting 7-nov-01 (here and below)
!===================================================================================================
      call gethhmmss(time1)
!===================================================================================================
      write(99,650) time1, hh,mm,ss
650   format('  Reading 201 file at:   ',f10.4,
     &       ' (',i2.2,':',i2.2,':',f5.2,')')
      write(99,*) '  Initial v201ed (1-10): ', v201ed
      write(99,*) '  Initial is201 (1-',ns,'): ', is201

!     gets rest of interval:
      if (time1.lt.tnext) then 
       do while (time1.lt.tnext)
!     stores previous record's data:
        do i=1,n1                         ! converts voltages 1..n1 for point k. 
!     sum of each over k (after converting to correct units):
!.....cc and offset are read from lsubs here modified for humidity sensor nrc 180302......... 
         if ((idate.ge.080101).and.(i.ge.9)) then 
           cc(i) = 1
           offset(i) = 0 
         end if
!............................................................................................       
         xo(i) = xo(i)+(v201ed(i)*cc(i))+offset(i) 
        end do

!     extracts relevant status bits:

!     calibration-status
        ical = 2*is201(4) + is201(5)       ! 00 = 0 (data sample)
                                           ! 01 = 1 (normal cal, incl zero)
                                           ! 10 = 2 (dilution)
                                           ! 11 = 3 (hi-flow zero)
        statcal(ical) = statcal(ical) + 1  ! number of times in each calibration status

!     chilled-mirror status
        statcm = statcm + is201(9)
!     licor ready light
!     variables
!     ---------
!     vrl = voltage licor ready light; the voltage of the ready light
!     vminrl     = the minimum voltage for a positive ready light
!     badlicorrl = true/false licor ready light; an array of ones and zeroes indicating
!     note: prior 5/9/02, this was licorrlvec (now is stretched version for fast data) whether the 
!           licor ready light was on (0) or off (1).
!     licoroff   = number of times licor ready light is off

       vrl = -1*v201ed(irl)*cc(irl) + offset(irl)                
       if (vrl.lt.vminrl) then
        badlicorrl(k) = 1              ! ready light off, we filter (k=this data point)
        licoroff = licoroff + 1
       else
        badlicorrl(k) = 0              ! ready light on, we do not filter
       end if              

!     reads next-line of data in this interval
!===================================================================================================
       call read20n(88, 201, dtype, time1, v201ed, 1, 10, is201, 1, ns, 
     &             done, yrdays, debug)
!===================================================================================================
       if (done) goto 899

!     input arguments :    fileptr (88), line-num (201), dtype (EA/EB)
!     output arguments:    time1, v201(1..10), is201(1..ns), done (flag), debug reset to false
!                          time1=jd1+gmt1/dble(86400.)
       k = k+1              ! number of samples in this interval.
       end do           
       k = k-1
      end if
!===================================================================================================
      call gethhmmss(time1)
!===================================================================================================
      write(99,652) time1, hh,mm,ss, k
652   format('  Stopped reading 201 file at:  ',f10.4,
     &       ' (',i2.2,':',i2.2,':',f5.2,'), after ',i5,' records.')
      goto 901
      
899   continue

!     inserts error code for premature end-of file here:
      write(99,*) '  Warning: Premature end-of-file for 201 file.' 
      write(*,*)  '  Warning: Premature end-of-file for 201 file.' 

901   continue

      if (k.le.kavg/2 .or. k.ge.kavg*2) then 
      write(99,259) k, kavg
655   format('  Warning: bad number of values read:  k= ',i5,
     &       ' (expected kavg= ',i5,').')
      if (k.le.kavg/20) then
      write(99,*)'  Skipping this interval...'
      write(99,*)          ! doesn/t skip unless we absolutely have to (k<45 when it should be 900).
      goto 160             ! skips to end of loop for processing this interval.
      end if
      end if
      write(99,270) k, 201                 ! aggregating k values.
      write(99,*)

!     now take averages:
      do i =1,n1                           ! for each var of input line.
       xo(i) = xo(i)/k
      end do
      do i=1,3                             ! for zero/cal status bits.
       statcal(i) = statcal(i)/k
      end do
      statcm   = statcm/k                  ! chilled mirror status
      miragdew = xo(9)                     ! mirsumdew/mirnum (subroutine lfluxvar) 2008-> RH
      miragt = xo(10)                      ! mirsumt/mirnum   (subroutine lfluxvar) 2008-> tair deg C    
      pres   = xo(ipcell)                  ! licor cell pressure in this half-hour
      temp   = xo(itcell)                  ! licor cell temperature in this half-hour

!     End of housekeeping/ chilled mirror code within loop / HC2-S23L Humidity and temperature sensor 

!     Modify the master filtering vector:  
!     Discard this data-point if either right sonic flag set or if licor ready-light off.

!     Implementation on 9-may-02 (srs):
!     Because the previous scheme was wrong, not taking into account that licor status 
!     is slow data (1/2 hz) but filtvec is matched to fast data (8hz).
!     Also, we want to keep track of bad licor and bad sonic separately (if licor is out, 
!     we could still have good sonic data).

      nrep = int(n/k+0.5) 
!     where:    
!     nrep = number PF replicates: expand to fill fast-day vector

      nrem = n - k*nrep   
!     where:    
!     nrem = number of remaining points
!     n    = number of fast points, k = number of slow points

!     nominally 16 fast points per slow (for 8 hz and 1/2 hz)
      do i = 1,k
      if (badlicorrl(i).eq.1) then
      ibase = (i-1)*nrep               ! baseline fast point corresponding to slow-point i
      do j = ibase+1,ibase+nrep
      if (j.le.n) licorrlvec(j) = 1    ! check if n/k not integer (nrem= -)
      end do                           ! end of the fast segment corresponding to slow-point i
      if (i.eq.k .and. nrem.gt.0) then ! check if m/k not integer (nrem= +)
      do j=1,nrem
      licorrlvec(ibase+nrep+j) = 1
      end do
      end if                           ! if remaining points at end.
      end if
      end do                           ! goes on to the next slow point i.

!     code to calculate fraction of bad sonic. we determine bad sonic readings by determining number
!     of -99.999's in temperature on the .201 file. choice of temperature was arbitrary. 
!     u, v, w would also have worked.

!...................................................................................................
!     note: one must run this before the detrending subroutine, as the detrending routine modifies 
!           the values.
!           as an addition, we will find the number of points that have been filtered.
!...................................................................................................

      numbadsonic  = 0               ! hard sonic errors
      numfiltsonic = 0               ! filtered sonic errors + hard sonic errors.
      numbadli = 0                   ! licor ready light off                (added in 9-may-2002)
      numeitherbad = 0               ! intersection of sonic and licor bad  (added in 9-may-2002)
      do i=1,n
      if (t(i).eq.-99.999) then
      numbadsonic = numbadsonic+1
      end if
      if (filtvec(i).eq.1 .or. t(i).eq.-99.999 ) then
      numfiltsonic = numfiltsonic+1  ! numfiltered=numfiltered+1

!...................................................................................................
!     note:  hard sonic errors (-99.999) also separately removed in detrend subroutine.
!...................................................................................................

      end if
      if (licorrlvec(i).eq.1) then
      numbadli = numbadli+1                        ! added in 9-may-2002
      end if
      if (t(i).eq.-99.999 .or. filtvec(i).eq.1 
     &                    .or. licorrlvec(i).eq.1) then
      numeitherbad = numeitherbad+1                ! added in 9-may-2002
      end if
      end do
      fracnotused = numeitherbad/dble(n)           ! numnotused/dble(n) (change 9-may-02)


!---------------------------------------------------------------------------------------------------
!     step 1: calculates the h2o and co2 mixing ratios;
!             obtains information from the correct calibration period;
!             interpolates zero and calibration information for the time of the aggregation period:
!             calculates the h2o concentration and then the co2 concentration.
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!     1(a). zero interpolation:
!     skips zero interpolation for refit.  use zeros from calibrations only.
!---------------------------------------------------------------------------------------------------
      !added "if", as some collection data missing zero cal nrc 170207
      izero=0
      if (nofzero.gt.0) izero = ihuntf(tzero,nofzero,tmid,jlo) ! tzero(jlo)= last time <= tmid
      if (izero.eq.0) then                                     ! before the range of recorded zeros
      v0co2 = vzco2(1)
      v0h2o = vzh2o(1)
      izero = 1
      jlo   = 1
      if (zerflag(izero).gt.0) calib = calib + 200             ! 200 if zero bad, 0 otherwise

!     setting calib =2xx means:  unable to extract nearby zero
!     (reminder:  zerflag =0 if zero ok, 1 otherwise)
      elseif (izero .eq. nofzero) then                         ! after the range of recorded zeros
      v0co2 = vzco2(nofzero)
      v0h2o = vzh2o(nofzero)
      if (zerflag(nofzero).gt.0) calib = calib + 200
      elseif (izero .gt. 0 .and. izero .lt. nofzero) then      ! within range
      jlo    = izero
      deltat = tzero(izero+1)-tzero(izero)
      if (zerflag(izero).eq.0) then 

!     if previous zero good, check the next one:
      if (zerflag(izero+1).eq.0) then 
      delt = tmid-tzero(izero)        ! 7-nov-01:  interpolate to tmid (not time0)
      else      ! next zero bad, don't interpolate
      delt  = 0 ! setting delt=0 prevents interpolation below
      calib = calib + 20              ! x2x = next zero bad
      end if

!     water zero
      deltav = vzh2o(izero+1)-vzh2o(izero)
      dvdt   = deltav/deltat
      v0h2o  = vzh2o(izero)+dvdt*delt

!     co2 zero
      deltav = vzco2(izero+1)-vzco2(izero)
      dvdt   = deltav/deltat
      v0co2  = vzco2(izero)+dvdt*delt
      elseif (zerflag(izero+1).eq.0) then 

!     if previous zero bad, but next one good, use next
      calib = calib + 10            ! x1x=prev zero bad

!     water zero:
      v0h2o = vzh2o(izero+1)

!     co2 zero:
      v0co2 = vzco2(izero+1)
      else

!     both are bad:
      calib = calib + 200           ! 2xx=both zeros bad
      v0h2o = 0                     ! will not be used.  (?)
      v0co2 = 0
      end if
      end if 
      
      write(99,610) izero,nofzero,v0h2o,v0co2,calib,jlo
610   format('  After zero #', i3,' (of ',i3,
     &       '): interpolated h2o,co2 zeros=',2f8.4,
     &       ', calib=', i3,' jlo=',i3)
      if (calib.ge.10) then 
      write(99,*)'  Warning: LiCor was off for at least one nearby',
     &           ' calibration or zero'
      end if

!---------------------------------------------------------------------------------------------------
!     1(b). calibration coefficient interpolation:
!     note: logic of cal use is:  both nearby cals good, use.  if one is bad, use the other
!           if both are bad (or there is only one nearby and it is bad), then don't do.
!---------------------------------------------------------------------------------------------------
      !added "if", as some collection data missing zero cal nrc 170207
215   continue
      iecal=0
      if (nofzero.gt.0) iecal=ihuntf(tcalave,nofcal,tmid,klo)    ! index to last cal before tmid
      if (iecal .eq. 0) then
       coef(1) = coef1(1)
       coef(2) = coef2(1)
       coef(3) = coef3(1)
       po = pcellc(1,1)                          ! p/tcal(*,1) = calibratio p/t
       to = tcellc(1,1)
!       po = 496.94				 ! QA
!       to = 43.56
       iecal = 1
       klo   = iecal
      if (calflag(iecal).ge.3) calib = calib + 400
      elseif (iecal .eq. nofcal) then
       coef(1) = coef1(nofcal)
       coef(2) = coef2(nofcal)
       coef(3) = coef3(nofcal)
       po      = pcellc(nofcal,1)                ! p/tcal(*,1) = calibration p/t
       to      = tcellc(nofcal,1)
!       po = 496.94				 ! QA
!       to = 43.56
      if (calflag(iecal).ge.3) calib = calib + 400
      elseif (iecal .gt. 0 .and. iecal .lt. nofcal) then
       klo     = iecal

!---------------------------------------------------------------------------------------------------
!     possible values of calflag:
!     calflag =  0 if everything ok with licor during cal
!             = xxxn: n=1: previous zero bad (next ok); 
!                     n=2: next zero bad (prev ok); 
!                     n=3: both bad
!             = xxnx: n=1 if low, high spans indistinguishable; n=2 if time mismatch
!             = xnxx: n= # of cal-seq segments that were not long enough to average
!             = nxxx: n= # of span levels (1-3) during which licor ready-light off
!              (nb:  licor status during hi-flow zero, dilution not checked)
!---------------------------------------------------------------------------------------------------

      if (calflag(iecal).lt.3) then
!     if previous calibration is good, check next one :
       deltat = tcalave(iecal+1)-tcalave(iecal)  
      if (calflag(iecal+1).lt.3) then     ! if next calibratin good: interpolate.
       delt  = tmid-tcalave(iecal)
      else                         ! next calibration bad: don't interpolate
       delt  = 0                   ! setting delt=0 prevents interpolation below
       calib = calib + 60          !             =x6x:  next calibration bad, previous good
      end if

!---------------------------------------------------------------------------------------------------
!     when skipping zero interp:
!              dcdt = (vzco2(iecal+1)-vzco2(iecal))/deltat
!              v0co2= vzco2(iecal)+ dcdt*delt   ! 
!              dcdt = (vzh2o(iecal+1)-vzh2o(iecal))/deltat
!              v0h2o= vzh2o(iecal)+ dcdt*delt   ! 
!---------------------------------------------------------------------------------------------------

!     calibration coefficients:
       dcdt    = (coef1(iecal+1)-coef1(iecal))/deltat ! coef 1
       coef(1) = coef1(iecal) + (dcdt*delt)
       dcdt    = (coef2(iecal+1)-coef2(iecal))/deltat ! coef 2
       coef(2) = coef2(iecal) + (dcdt*delt)
       dcdt    = (coef3(iecal+1)-coef3(iecal))/deltat ! coef 3
       coef(3) = coef3(iecal) + (dcdt*delt)

!     p and t conditions:
       dcdt = (pcellc(iecal+1,1)-pcellc(iecal,1))/deltat ! po [= pcal(*,1)]
       po   = pcellc(iecal,1) + (dcdt*delt)
       dcdt = (tcellc(iecal+1,1)-tcellc(iecal,1))/deltat ! to [= tcal(*,1)]
       to   = tcellc(iecal,1)  + (dcdt*delt)
!       po = 496.94			   ! QA
!       to = 43.56
      elseif (calflag(iecal+1).lt.3) then 
       coef(1) = coef1(iecal+1)
       coef(2) = coef2(iecal+1)
       coef(3) = coef3(iecal+1)
       po      = pcellc(iecal+1,1) 
       to      = tcellc(iecal+1,1) 
       calib   = calib + 30                ! =x3x:  previous cal bad, next good
!       po = 496.94			   ! QA
!       to = 43.56
      else 
       coef(1) = 0
       coef(2) = 0
       coef(3) = 0
       po      = 0
       to      = 0
       calib   = calib + 400              ! = 4xx: both cals bad
      end if                              ! if calflag(iecal) < 3 (good), 
                                          ! calflag(iecal+1) < 3 (next good), or neither good
                                          ! if iecal = 0, nofcal, or something in between.
      end if      
      write(99,615) iecal, coef, calib
615   format('  After Eddy calibration #',i3,
     &       ': interpolated coefficients =', 3(1x,1pe11.4),
     &       ', calib =',i4)
      write(99,620) po, pres, to, temp
620   format(' : interp po, current pres, interp to, current temp =', 
     &        4f6.1)

!---------------------------------------------------------------------------------------------------
!     in 16-july-04:
!     finds water cal scale factor
!     uses:
!     tscalewv(1..nscale) -> vector of times (long time period, absolute decimal days since 
!                                             beginning 1/1/01)
!     scalewv(1..nscale) corresponding vector of hourly scale factors
!     tfind (real)
!     sets:
!     wvsf = water-vapor scale factor for this interval
!---------------------------------------------------------------------------------------------------        
!..................................................................................nrc 170127
      call readh2ocal(time0,dtype)
!...........................................................................................
       wvsf = dble(cvapor)
      write(99,621),'time and water vapour calibration factor', time0, wvsf
!       read(yymmdd, '(i2,i2,i2)') iyy,imm,idd  ! name of this file 
!       tfind  = tmid                           ! days since 1/1/01
!       iscale = 0                              ! nrc added following if (ihuntf) 170207 
!      if (nscale.gt.0)        iscale=ihuntf(tscalewv,nscale,tfind,mlo) ! tscalewv(jlo)= last time <= tmid
!      if (iscale.eq.0) then                    ! before the range of scale factors
!       wvsf   = scalewv(1)
!       iscale = 1
!       mlo    = 1
!      elseif (iscale.eq.nscale) then           ! after the range of scale factors
!       wvsf = scalewv(nscale)
!      elseif ((iscale.gt.0).and.(iscale.lt.nscale)) then   ! within range
!       mlo  = iscale
!       wvsf = scalewv(iscale)
!      end if
!      write(99,621) tfind, wvsf
621   format(' : water vapor scale factor for jd ',f10.6,': ',f10.6,1x,f6.3,1x,f6.3)
      read(yymmdd, '(i2,i2,i2)') iyy,imm,idd   ! name of this file 
      tfind = tmid                             ! days since 1/1/01
      ilag  = 0                                ! nrc added following if (ihuntf) 170207 
      if (nlagfac.gt.0) ilag = ihuntf(tlagfac,nlagfac,tfind,nlo) ! tscalewv(jlo)= last time <= tmid
      if (ilag.eq.0) then                      ! before the range of scale factors
       ilag   = 1
       lagco2 = co2co2(ilag)
       lagh2o = co2h2o(ilag)
       nlo    = 1
      elseif (ilag.eq.nlagfac) then           ! after the range of scale factors
       lagco2 = co2co2(ilag)
       lagh2o = co2h2o(ilag)
      elseif ((ilag.gt.0).and.(ilag.lt.nlagfac)) then   ! within range
       nlo    = ilag
       lagco2 = co2co2(ilag)
       lagh2o = co2h2o(ilag)
       nlo    = ilag
      end if
      write(99,622) tfind, lagco2, lagh2o
622   format(' : on JD ',f8.2,', lags = ',i3,'(CO2), ',i3,'(H2O)')

!---------------------------------------------------------------------------------------------------
!     1(c). derive concentrations from cal-curve for each point in this interval
!---------------------------------------------------------------------------------------------------
      do i = 1,n
       co2(i) = co2(i)-v0co2 ! *pres/po   ! %% zero-adjusted voltages (at pair=pres)
       h2o(i) = h2o(i)-v0h2o ! *pres/po
       if ((i.le.2).or.(i.ge.(n-2))) then
        write(99,625) i, co2(i), h2o(i)
       end if
       nmux = i
625    format('       ',i5,':  CO2, H2O (input V): ', 2(1x,f6.3) )
       h2o_v(i) = h2o(i)     ! keep voltage for later output nrc
       co2_v(i) = co2(i)
!       xint = h2o_v(i)
!       yint = co2_v(i)

!===================================================================================================
       call getconc2(h2o(i),co2(i),pres,po,temp,to,coef,wvsf)  ! 16-jul-04
!===================================================================================================
!     converts zero-adjusted voltages to mole fractions:
!     co2=ppm, h2o=mmol/mol
       if ((i.le.2).or.(i.ge.(n-2))) then
        write(99,627) co2(i), h2o(i)
       end if
627    format('            CO2, H2O (output concentration): ',2(1x,f6.2))
      end do
      write(99,629) h2o(n), co2(n)
629   format('                   : last H2O, CO2 concentration:', 2f6.1)

!---------------------------------------------------------------------------------------------------
!     step 2: calculates the stress tensor (wind covariance matrix).  
!            a. detrend the data for u, v, w.
!            b. record average speed during the half-hour into
!            uvw(1) for u, uvw(2) for v, uvw(3) for w.  
!            c. record the slopes and intercepts of the linear detrending line into ?int and ?slope.  
!            d. calculate the stress tensor.
!---------------------------------------------------------------------------------------------------
!     detrends u :
      write(99,*) '  Flux calculations:'
      write(99,*) '  Detrend wind.'

!===================================================================================================
      call detrend(n,u,avg,var,aint,slope,filtvec,numfiltsonic) 
!===================================================================================================
!     regress u against (1:n) (exclude bad points id'd by filtvec), return: 
!     u = residuals (u detrended) (filtered points set =0)
!     avg = average raw u
!     var = variance of detrended u (residual variance)
!     aint= intercept of detrend line
!     slope= slope of detrend line
      uvw(1)  = avg*gu     ! <u> (gain included to convert to m/sec)
      st(1,1) = var*gu*gu  ! st= wind covariance matrix, 1,1: sigma_uu
      uint    = aint
      uslope  = slope

!     detrends v:
!===================================================================================================
      call detrend(n,v,avg,var,aint,slope,filtvec,numfiltsonic) 
!===================================================================================================
      uvw(2)  = avg*gv     ! <v>
      st(2,2) = var*gv*gv  ! sigma_vv
      vint    = aint
      vslope  = slope

!     detrends w:
!===================================================================================================
      call detrend(n,w,avg,var,aint,slope,filtvec,numfiltsonic) 
!===================================================================================================
      uvw(3)  = avg*gw           ! <w>
      st(3,3) = var*gw*gw        ! sigma_ww
      wint    = aint
      wslope  = slope
      st(1,2) = 0.d0             ! off-diagonal (covariance) terms of wind covariance matrix
      st(1,3) = 0.d0
      st(2,1) = 0.d0
      st(2,3) = 0.d0
      st(3,1) = 0.d0
      st(3,2) = 0.d0
      do i = 1,n                 ! accumulate covariances (except bad points set=0, don't add)
       st(1,2) = st(1,2)+dble(u(i))*dble(v(i))
       st(1,3) = st(1,3)+dble(u(i))*dble(w(i))  ! momentum flux
       st(2,3) = st(2,3)+dble(v(i))*dble(w(i))
      end do

!     normalized by number of points - number of bad sonic points
      st(1,2) = st(1,2)*gu*gv/dble(n-numfiltsonic) ! normalize
      st(1,3) = st(1,3)*gu*gw/dble(n-numfiltsonic)
      st(2,3) = st(2,3)*gv*gw/dble(n-numfiltsonic)
      st(2,1) = st(1,2)          ! set symmetric terms equal
      st(3,1) = st(1,3)
      st(3,2) = st(2,3)

!     step 3: rotate the wind data so that 1) the direction is correct and 
!                                          2) that the mean vertical wind is 0.0 m/s.
!---------------------------------------------------------------------------------------------------
      write(99,*) '  Rotating wind.'
!===================================================================================================
      call rotate(angles,st,uvw,rot)
!===================================================================================================
!     input:  uvw(3), st(3x3) = wind mean vector and covariance matrix
!     output: st(3x3)= rotated covariance matrix
!             angles(3) = phi (about z), theta (about y), beta (about x) 
!                     (all in degrees)
!             rot(3x3) = rotation matrix (in common, accessed by covar)

!---------------------------------------------------------------------------------------------------
!     step 4: calculate the heat flux before any smearing is done.  this
!             will give the true heat flux and average temperature data.
!             vertical heat flux will be stored in flux(1,3).
!---------------------------------------------------------------------------------------------------
      do i = 1,3
       avgs(i)  = 0.
       vars(i)  = 0.
       f9avg(i) = 0.
       f9max(i) = 0.
       f9std(i) = 0.
       do j = 1,3
        flux(i,j) = 0.0
       end do
      end do
! 
!     calculates sensible heat fluxes:
!===================================================================================================
      call detrend(n,t,avg,var,aint,slope,filtvec,numfiltsonic)
!===================================================================================================

!     regress t against (1:n), and return:   t    = residuals (t detrended)
!                                            avg  = average raw t
!                                            var  = variance of detrended t (residual variance)
!                                            aint = intercept of detrend line
!                                            slope= slope of detrend line
      avgs(1) = avg*gt
      vars(1) = var*gt*gt
      tint    = aint
      tslope  = slope
      write(99,*) '  Calculating sensible heat flux.'
!===================================================================================================
      call covar(cov, f9bus, 0, gt, numfiltsonic)
!===================================================================================================
      ! uses u,v,w,t in common to calculate <u't'>, <v't'>, <w't'>
      ! and store result in cov(3) (note: uses rotated winds)
      !   f9bus(4) is lag correlation summary across 9 lags (+/-4)
      !   3rd arg (=0 for heat flux) is the lag between t and wind
      !   4th arg (gt) is the gain to convert t to true temperature
      !   5th arg (numfiltsonic):  points to throw out when normalizing covar

      flux(1,1) = cov(1)  ! <u't'> (rotated)
      flux(1,2) = cov(2)  ! <v't'>
      flux(1,3) = cov(3)  ! <w't'>
      fluxt     = cov(3)  ! store vertical flux term again (won't be smeared)
      f9avg(1)  = f9bus(1)
      f9lag(1)  = f9bus(2)
      f9max(1)  = f9bus(3)
      f9std(1)  = f9bus(4)

!---------------------------------------------------------------------------------------------------
!     step 5: data smearing (filtering) to get co2 and h2o correction factor:
!---------------------------------------------------------------------------------------------------     
!     <t'w'> / <smear(t')*smear(w')>
!     where smear(x) [=sm(x)] is recursive exponential filter: 
!     smear(x)_t = alpha*smear(x)_t-1 + (1-alpha)*x_t
!     note:  raw flux = <w'c'>
!     corrected flux = <sm(w')*c'>* (<w't'>/<sm(w')*sm(t')>)
      do j = 1,n
       ht(j) = t(j)         ! replicates for water flux
       hw(j) = w(j)
      end do
      write(99,*) '  Filter wind and temperature for correction factor.'
!===================================================================================================
      call smear(u,n,tau)   ! returns u smeared for co2
!===================================================================================================
      call smear(v,n,tau)   ! smear v for co2
!===================================================================================================
      call smear(w,n,tau)   ! smear w for co2
!===================================================================================================
      call smear(t,n,tau)   ! smear t for co2
!===================================================================================================
      call smear(hw,n,tauh) ! smear w for water
!===================================================================================================
      call smear(ht,n,tauh) ! smear t for water
!===================================================================================================
      do j = 1,n
       sw(j) = w(j)         ! store smeared w for co2
      end do
!===================================================================================================
      call covar(cov, f9bus, 0, gt, numfiltsonic)  ! smeared heat flux (for co2)
!===================================================================================================
      ratio = 1.
      if (fluxt .ne. 0.) ratio = cov(3)/fluxt 
!     ratio for co2:   <sm(t')*sm(w')>/<t'w'>
      if (abs(ratio) .gt. 10.) ratio = -999.
!     legacy from boreas:  boreas code stores heat flux as the smeared heat flux.  
!     what we really want for heat flux is un-smeared version (though this could
!     be recovered by dividing by the flux correction factor).
!        flux(1,1)=cov(1)  ! <sm(u')*sm(t')>  (smeared for co2)
!        flux(1,2)=cov(2)  ! <sm(v')*sm(t')>
!        flux(1,3)=cov(3)  ! <sm(w')*sm(t')>
!     smear for water
      do i=1,n
       t(i) = ht(i)        ! t smeared for h2o
       w(i) = hw(i)        ! w smeared for h2o
      end do
!===================================================================================================
      call covar(cov,f9bus,0,gt,numfiltsonic) ! smeared heat flux (for h2o)
!===================================================================================================
      hratio = 1.
      if (fluxt .ne. 0.) hratio = cov(3)/fluxt  
      ! ratio for h2o:  <sm(t')*sm(w')>/<t'w'>
      if (abs(ratio) .gt. 10.) hratio = -999.

!---------------------------------------------------------------------------------------------------
!     step 6: calculate the covariances for h2o and co2 using the smeared data.
!---------------------------------------------------------------------------------------------------
!...................................................................................nrc 170201 start
      call detrend(n,h2o_v,avg,var,aint,slope,licorrlvec,numbadli)       ! voltages h2o and co2o out
      xint    = aint      
      call detrend(n,co2_v,avg,var,aint,slope,licorrlvec,numbadli)
      yint    = aint
!............................................................................................nrc end

630   format('  ',i5,': H2O concentration: ', (1x,1pe10.4) )
      write(99,*) '  Detrend H2O for water flux.'
!===================================================================================================
      call detrend(n,h2o,avg,var,aint,slope,licorrlvec,numbadli)  ! filtvec,numfiltsonic)
!===================================================================================================
!     licorrlvec indicates licor ready light off             ! (changed 9-may-02)     
      write(99,*)'  (mean H2O =',avg
      avgs(2) = avg
      vars(2) = var
      hint    = aint
      hslope  = slope
      do i=1,n
       t(i) = h2o(i)                                         ! h' to calculate water fluxes
!     note:  w already smeared for water.
      end do
      write(99,*) '  Calculating water flux, using lagh2o =',lagh2o
!===================================================================================================
      call covar(cov,f9bus,lagh2o,1.d0,numeitherbad)
!===================================================================================================
!     numeitherbad counts intersection of bad sonic with bad licor
      flux(2,1) = cov(1)  ! <h' * sm(u')>  (note:  u,v still smeared for co2)
      flux(2,2) = cov(2)  ! <h' * sm(v')>  
      flux(2,3) = cov(3)  !  = <h' * sm(w')> 
      f9avg(2)  = f9bus(1)
      f9lag(2)  = f9bus(2)
      f9max(2)  = f9bus(3)
      f9std(2)  = f9bus(4)
      do j=1,n
       w(j) = sw(j)       ! reset w to w smeared for co2
      end do
632   format('  ',i5,': CO2 conc: ', (1x,1pe10.4) )
      write(99,*)'  Detrend CO2 for flux.'
!===================================================================================================
      call detrend(n,co2,avg,var,aint,slope,licorrlvec,numbadli)  !filtvec,numnotused)
!===================================================================================================
      avgs(3) = avg
      vars(3) = var
      cint    = aint
      cslope  = slope
      write(99,*)' (mean CO2 =',avg
      do i = 1,n
       t(i) = co2(i)
      end do
      write(99,*) '  Calculating CO2 flux, using lagco2 =',lagco2
!===================================================================================================
      call covar(cov,f9bus,lagco2,1.d0,numeitherbad)
!===================================================================================================
      flux(3,1) = cov(1)   ! <c' * sm(u')>  (u,v still smeared for co2)
      flux(3,2) = cov(2)   ! <c' * sm(v')>
      flux(3,3) = cov(3)   ! <c' * sm(w')>
      f9avg(3)  = f9bus(1)
      f9lag(3)  = f9bus(2)
      f9max(3)  = f9bus(3)
      f9std(3)  = f9bus(4)

!---------------------------------------------------------------------------------------------------
!     step 7: write the output:
!---------------------------------------------------------------------------------------------------
50    speed = sqrt((uvw(1)*uvw(1))+(uvw(2)*uvw(2))+(uvw(3)*uvw(3)))
      write(99,*)'  Write the output.'
      if (uvw(1) .gt. 0.0) then
       direc = dirofu-angles(1)
      else
       direc = dirofu+angles(1)-180.0
      end if
      if (direc .ge. 360.) direc=direc-360.
      if (direc .lt. 0.)   direc=direc+360.
      tdoy = tstart-days00-yrdays0 
      !..................................................................................
      !nrc added set boundaries of what is possible given the write format constrains
      if ((speed.lt.-99.999).or.(speed.gt.999.999).or.(speed.ne.speed*1.00))    speed = -999.99	!f8.3
      if ((direc.lt.-99.999).or.(direc.gt.999.999).or.(direc.ne.direc*1.00))    direc = -999.99	!f8.3
      if ((hratio.lt.-9999.999).or.(hratio.gt.99999.999).or.(hratio.ne.hratio*1.00))  hratio = -999.99	!1pe11.4
      if ((ratio.lt.-9999.999).or.(ratio.gt.99999.999).or.(ratio.ne.ratio*1.00))      ratio = -999.99	!1pe11.4

      if ((st(2,2).lt.-99.999).or.(st(2,2).gt.999.999).or.(st(2,2).ne.st(2,2)*1.00)) st(2,2) = -999.99	!f8.3
      if ((st(3,3).lt.-99.999).or.(st(3,3).gt.999.999).or.(st(3,3).ne.st(3,3)*1.00)) st(3,3) = -999.99	!f8.3

      do i=1,3
        if ((st(1,i).lt.-999.999).or.(st(1,i).gt.999.999).or.(st(1,i).ne.st(1,i)*1.00)) st(1,i) = 999.999	!f8.3
        if ((flux(1,i).lt.-9999.999).or.(flux(1,i).gt.9999.999).or.(flux(1,i).ne.flux(1,i)*1.00)) flux(1,i) = -999.99!f8.3
        if ((flux(2,i).lt.-9999.999).or.(flux(2,i).gt.9999.999).or.(flux(2,i).ne.flux(2,i)*1.00)) flux(2,i) = -999.99!f8.3
        if ((flux(3,i).lt.-9999.999).or.(flux(3,i).gt.9999.999).or.(flux(3,i).ne.flux(3,i)*1.00)) flux(3,i) = -999.99!f8.3
        if ((angles(i).lt.-9999.999).or.(angles(i).gt.9999.999).or.(angles(i).ne.angles(i)*1.00)) angles(i) = -999.99!f8.3
        if ((f9avg(i).lt.-9999.999).or.(f9avg(i).gt.9999.999).or.(f9avg(i).ne.f9avg(i)*1.00)) f9avg(i) = -999.99	!f8.3
        if ((f9std(i).lt.-9999.999).or.(f9std(i).gt.9999.999).or.(f9std(i).ne.f9std(i)*1.00)) f9std(i) = -999.99	!f8.3
        if ((f9max(i).lt.-9999.999).or.(f9max(i).gt.9999.999).or.(f9max(i).ne.f9max(i)*1.00)) f9max(i) = -999.99	!f8.3
        if ((f9lag(i).lt.-9999.999).or.(f9lag(i).gt.9999.999).or.(f9lag(i).ne.f9lag(i)*1.00)) f9lag(i) = -999.99	!f8.3
        if ((uvw(i).lt.-999.999).or.(uvw(i).gt.999.999).or.(uvw(i).ne.uvw(i)*1.00))             uvw(i) =  99.99999	!f9.6
        if ((st(i,i).lt.-99.999).or.(st(i,i).gt.99.999).or.(st(i,i).ne.st(i,i)*1.00))       st(i,i) = 99.999	!f8.3
        if ((avgs(i).lt.-9999.999).or.(avgs(i).gt.9999.999).or.(avgs(i).ne.avgs(i)*1.00))  avgs(i) = -999.99	!f11.4 the h20 co2
        if ((vars(i).lt.-9999.999).or.(vars(i).gt.9999.999).or.(vars(i).ne.vars(i)*1.00))  vars(i) = -999.99	!f11.4
        if (numbitarray(i).ne.numbitarray(i)*1.00)  numbitarray(i) = 99
      enddo

      if ((tslope.lt.-9.999) .or.(tslope.gt.99.999).or.(tslope.ne.tslope*1.00))     tslope = -999.99	!f9.6
      if ((wslope.lt.-9999.999).or.(wslope.gt.9999.999).or.(wslope.ne.wslope*1.00))   wslope = -999.99	!f8.3
      if ((cslope.lt.-9999.999).or.(cslope.gt.9999.999).or.(cslope.ne.cslope*1.00))   cslope = -999.99	!f8.3
      if ((uslope.lt.-9999.999).or.(uslope.gt.9999.999).or.(uslope.ne.uslope*1.00))   uslope = -999.99	!f8.3
      if ((hslope.lt.-9999.999).or.(hslope.gt.9999.999).or.(hslope.ne.hslope*1.00))   hslope = -999.99	!f8.3
      if ((vslope.lt.-9999.999).or.(vslope.gt.9999.999).or.(vslope.ne.vslope*1.00))   vslope = -999.99	!f8.3
   
      if ((uint.lt.-9999.999).or.(uint.gt.9999.999).or.(uint.ne.uint*1.00))   uint = -999.999	        !f8.3
      if ((vint.lt.-9999.999).or.(vint.gt.9999.999).or.(vint.ne.vint*1.00))   vint = -999.999	!f8.3
      if ((wint.lt.-9999.999).or.(wint.gt.9999.999).or.(wint.ne.wint*1.00))   wint = -999.999	!f8.3
      if ((tint.lt.-9999.999).or.(tint.gt.9999.999).or.(tint.ne.tint*1.00))   tint = -999.999	        !f8.3
      if ((hint.lt.-9999.999).or.(hint.gt.9999.999).or.(hint.ne.hint*1.00))   hint = -999.999	        !f8.3
      if ((cint.lt.-9999.999).or.(cint.gt.9999.999).or.(cint.ne.cint*1.00))   cint = -999.999	        !f8.3
      
      if ((yint.lt.-99.999) .or.(yint.gt.99.999).or.(yint.ne.yint*1.00))       yint = -9.99	        !f9.6
      if ((xint.lt.-99.999) .or.(xint.gt.99.999).or.(xint.ne.xint*1.00))       xint = -9.99	        !f9.6

      if ((temp.lt.-99.999)    .or.(temp.gt.999.999).or.(temp.ne.temp*1.00))             temp  = -99.99	!f7.3
      if ((miragdew.lt.-99.999).or.(miragdew.gt.99.999).or.(miragdew.ne.miragdew*1.00))  miragdew = -99.99	!f7.3
      if ((miragt.lt.-99.999).or.(miragt.gt.99.999).or.(miragt.ne.miragt*1.00))          miragt   = -99.99	!f7.3

      if (numbadsonic.ne.numbadsonic*1.00)        numbadsonic = 99
      if (numbitarray(4).ne.numbitarray(4)*1.00)  numbitarray(4) = 99
      if (numbadli.ne.numbadli*1.00)              numbadli = 99
      if (fracnotused.ne.fracnotused*1.00)        fracnotused = 99

      !no more crazy dates
      if ((tdoy.ge.1.).and.(tdoy.lt.367.).and.(tstart.ge.1.).and.(tstart.lt.100000.).or.(tstart.eq.tstart*1.00)) then
      !..................................................................................
       write(19,70)tstart,tdoy,calib,speed,direc,hratio,ratio
       write(19,80)st(1,1),st(1,2),st(1,3),st(2,2),st(3,3) 
       write(19,80)flux,angles  ! flux(3x3)=9 elements; angles=3 elements: on 2 lines
       write(19,80)f9avg,f9max  ! 3, 3
       write(19,80)f9std,f9lag  ! 3, 3
       write(19,80)avgs,vars    ! 3, 3
       write(19,80)uint,vint,wint,tint,hint,cint
       write(19,80)uslope,vslope,wslope,tslope,hslope,cslope
!     nrc added to obtain H2Ocal fn Tair.vs.Tson and CO2 cal without tanks
70     format(f10.4,1x,f8.4,i4,2(1x,f8.3),1x,1pe10.3,1x,1pe10.3)
80     format(6(1x,1pe11.4))
!     2002-2006 chilled mirror file output and bad sonic file output
!     2008-2017 humidity and temperature output and bad sonic file output
       if (idate.ge.080101) then
        write(19,86)miragdew,miragt,pres,temp
       else
        write(19,86)miragt,miragdew,pres,temp
       end if
       write(19,87)xint,yint,po,to	              		! h2o.v co2.v p.sample t.sample nrc added 170126
!       write(19,87)h2o_v,co2_v,po,to	              		! h2o.v co2.v p.sample t.sample nrc added 170126
       write(19,88)uvw(1),uvw(2),uvw(3),st(1,1),st(2,2),st(3,3)	! u v w u_std v_std w_std  nrc added 170510 wind speed mean and std
       write(19,89) numbadsonic,numbitarray(1),numbitarray(2),
     &             numbitarray(3),numbitarray(4), numbadli, 
     &             fracnotused*100.d0, n  ! licoroff
!     ^ counts out of 14400 (fast) not 900 (slow)
86     format(1x,f7.3, 5x,f7.3, 4x, f9.3, 4x, f7.3) 
87     format(1x,f9.6, 1x,f9.6, 1x, f9.3, 4x, f7.3)   ! nrc added 170126
88     format(6(1x,f9.6))   			      ! nrc added 170510
89     format(1x, 6(1x,i5), 1x, f9.4, 1x, i5)
!     end of flux calculation (steps 1-7)
      endif		                              ! endif no more crazy dates
!---------------------------------------------------------------------------------------------------
!     fourth -- continue loop, then close up files all used here.
!---------------------------------------------------------------------------------------------------
160   continue
      numint= numint+1                       	! number of aggregation intervals processed.
      tstart= tnext                          	! start time for next interval.
      tnext = tstart + period/1440.d0        	! next time for next interval.

!     minutes to decimal julian days:
      tmid = tstart + (period/2)/1440.d0     	! mid-point of next period.   
      write(99,350) calib, endzer, nseg 
350   format('  End the loop.  (calib = ',i4,'; endzer = ',
     &         i2,'; nseg = ',i1,')')
      write(99,*)
       
!===================================================================================================
      call gethhmmss(tstart)
!===================================================================================================
      write(99,220) numint+1, tstart,hh,mm,ss   
      call gethhmmss(tnext)
      write(99,221) tnext,hh,mm,ss
!     err found 29-may-02 :  above causes problem when calib set to, e.g. 3x or 6x (nearby cal bad)
!                            because fails to restore data from depository if zero at beginning of 
!                            this (if nearby calibration is bad, and zero at beg, calib=62, and does
!                            not pass if-statement above)
      if (calib.ge.1 .and. nseg.eq.1 .and. endzer.eq.0) then  
!     zero at beginning, when:        calib >=1    calib = x1, x2, or x3 if zero    and 
!                                     nseg=1       set to go one couple min segment past end
!                                     (note:  for full cal-seq, nseg is set = 0)
!                                     endzer = 0     no zero at end (probably redundant)
!                                     then:  restore last segment of data from depository 
      write(99,247) novlap, tovlap*86400.  ! notify:  restoring to top of array from end of this 1/2
      do i = 1, novlap                     ! hour (which went over to compensate for calibration at 
       u(i) = store(i,1)                   ! beginning) to the beginning of the next 1/2 hr (so it
       v(i) = store(i,2)                   ! will start on 1/2 hour boundary). results in several 
       w(i) = store(i,3)                   ! minutes of overlap between 1/2 hours.
       t(i) = store(i,4)
       co2(i)  = store(i,5)
       h2o(i)  = store(i,6)
       stat(i) = store(i,7)
      end do
      n = novlap    ! set n to be the number of points just recovered from depository
      nseg=0        ! reset num of extra 'tovlap' segments to read in
      calib=0       ! reset cal flag
      goto 15       ! file loop:  read next aggregation interval
!     continue by reading remainder of file (after first few minutes)
      end if
      calib=0
      goto 10                     ! file loop:  read next aggregation interval
!     continues by reading first few minutes again (repeats until the end-of-file is reached with
!                                                   read statment)

!     exit from one of the file-read loops

!    (1) exit due to premature end-of-file error:
290   continue
      write(99,*)
      write(99,*)'  Error:  end-of-file on first read for: ',
     &    fnamein(1:lfnin)//'.200'
      print *,   '  Error:  end-of-file on first read for: ',
     &    fnamein(1:lfnin)//'.200'
      goto 400
      
291   continue
      write(99,*)
      write(99,*)'  Error:  end-of-file before start time in file: ',
     &    fnamein(1:lfnin)//'.200' 
      print *,   '  Error:  end-of-file before start time in file: ',
     &    fnamein(1:lfnin)//'.200' 
      goto 400
 
!     (2) exit during later read (probably normal)
300   continue
      write(99,*)
      write(99,*)'  Discard interval:  reached end of: ',
     &    fnamein(1:lfnin)//'.200' 
      time0 = julday+gmt/86400.d0
!===================================================================================================
      call gethhmmss(time0)
!===================================================================================================
      write(99,237) time0, hh,mm,ss, n
      print *
      print *, '  Reached end-of-file in file: ',
     &    fnamein(1:lfnin)//'.200' 
      write(*,237)  time0, hh,mm,ss, n

400   continue
      tnext = tnext - period/1440.d0       ! minutes to decimal julian days
      write(99,*)
      write(99,*)'  Finished processing "',fnamein(1:lfnin), 
     &           '" input file.'
      print *,   '  Finished processing "',fnamein(1:lfnin), 
     &           '" input file.'
      write( *,410) tnext-tstart1, numint, period
      write(99,410) tnext-tstart1, numint, period
410   format('  (',f7.3,' days of data aggregated into ',i3,1x,f6.2,
     &         '-min intervals)')
      write(*,*)                           ! blank line
      write(99,*)
      close(19)                            ! close output file (*.flx)
      close(16)                            ! close *.hkp
      close(20)                            ! close input files
      close(88)                            ! close chilled mirror file (.201)
      close(99)                            ! close error/log file
      if (filtset(5) .eq. 1) then          ! pseudofilter 1 file
      close(55)
      end if        
      end do                ! ending 'do infile=1, nifiles  (loop over each set of input files)
990   continue
      print *, '  End of program LFLUX'
      print *
      end  

!===================================================================================================
!     end of program
!===================================================================================================


!===================================================================================================
!                                    subroutines and functions                                     !
!===================================================================================================
      include 'lsubs.for' 
      include 'lfluxsub.for'                             ! (contains code from lmirror and lfilter)

!===================================================================================================
!                                      subroutine detrend                                          !
!=================================================================================================== 
!     input:  x (raw data), ndata (number of data)
!     added filter-control vector for filtering, a 1 indicates to filter on that point ntossed
!     number of points filtered or bad ("tossed" out).
!     return: x (detrended data ( = 0 if filter set), xmean, var, a = intercept, b = slope.
!     remove bad data points before detrending (bad if value =-99.999)
!---------------------------------------------------------------------------------------------------
      subroutine detrend(ndata,x,xmean,var,a,b,filter,ntossed)
      real*8 x(ndata)
      real*8 xmean,var,a,b
      real*8 sx,ss,si,st2,sioss,vara
      real*8 aa,ba ,y
      integer filter(ndata)
      integer ntossed  
      sx=0.d0
!     modified ss to subtract off the number of points not used      
      ss    = dble(ndata)-dble(ntossed)
      si    = ss*(ss+1.)/2.d0
      st2   = 0.d0
      ba    = 0.d0
      vara  = 0.d0
      sioss = si/ss
      do 14 i=1,ndata
       if (x(i) .ne. -99.999 .and. filter(i).ne.1) then
        sx  = sx+dble(x(i))
        y   = dble(i)-sioss
        st2 = st2+y*y
        ba  = ba+y*dble(x(i))
       end if     
14    continue
      xmean = sx/ss
      ba    = ba/st2
      aa    = (sx-si*ba)/ss
      do 16 i = 1,ndata
      if(x(i) .ne. -99.999 .and. filter(i).ne.1) then
       x(i) = x(i)-aa-ba*dble(i)
      else
       x(i) = 0.0
      end if
      vara =vara+dble(x(i))*dble(x(i))            
16    continue
      var  = vara/ss
      a    = aa
      b    = ba
      return
      end
!===================================================================================================
!                                     end of subroutine detrend                                    !
!===================================================================================================



!===================================================================================================
!                                         subroutine rotate                                        !
!===================================================================================================
!     this program does axis rotation adopted from moore and fitz
!     originally from mcmillan
!     the rotation is done in two steps:
!     first,  rotate an angle of "phi" about z and then an angle of "theta" about new y, 
!             so that the new x (u) axis represents the streamline;
!     second, rotate an angle of "beta" about the new x to make y-component surface 
!             st zero:  [v'w'] = 0.
!
!     enter:  mean (u,v,w), stress tensor, scalar fluxes
!     return: corrected stress tensor and rotations
!---------------------------------------------------------------------------------------------------
      subroutine rotate(angles, st, uvw, rot)
      implicit real*8 (a-h,o-z)
      real*8 st(3,3), angles(3), cov(3), uvw(3), rot(3,3)
      real*8 bet(3,3)
      real*8 tmp(3,3)
      rad2deg=57.29578

!     next define rotation angles: phi, cta (theta)
      usq=uvw(1)**2
      vsq=uvw(2)**2
      wsq=uvw(3)**2
      wind1=sqrt(usq+vsq)
      wind2=sqrt(usq+vsq+wsq)
      cphi=uvw(1)/wind1
      sphi=uvw(2)/wind1
      ccta=wind1/wind2
      scta=uvw(3)/wind2

!     next define the first rotation matrix
      rot(1,1) = cphi*ccta
      rot(1,2) = sphi*ccta
      rot(1,3) = scta
      rot(2,1) = -sphi
      rot(2,2) = cphi
      rot(2,3) = 0.0 
      rot(3,1) = -cphi*scta
      rot(3,2) = -sphi*scta 
      rot(3,3) = ccta

!     next calculate stress tensor after first rotation
      do 2 i=1,3
      do 2 j=1,3
      tmp(i,j) = 0.0
      do 2 k=1,3
      do 2 l=1,3
      tmp(i,j) = tmp(i,j)+rot(i,k)*rot(j,l)*st(k,l)
2     continue

!     next define rotation angle: beta
      diff = tmp(2,2)-tmp(3,3)
      if (diff.eq. 0.0) then
      beta = 1.5708
      else
      tan2b = 2.0*tmp(2,3)/(tmp(2,2)-tmp(3,3))
      beta  = 0.5*atan(tan2b)
      endif 
      cbeta = cos(beta)
      sbeta = sin(beta)

!     next define the second rotation matrix
      bet(1,1) = 1.0
      bet(1,2) = 0.0
      bet(1,3) = 0.0
      bet(2,1) = 0.0
      bet(2,2) = cbeta
      bet(2,3) = sbeta
      bet(3,1) = 0.0
      bet(3,2) = -sbeta
      bet(3,3) = cbeta

!     next calculate the stress tensor after second rotation
      do 4 i=1,3
      do 4 j=1,3
      st(i,j)=0.0
      do 4 k=1,3
      do 4 l=1,3
      st(i,j)=st(i,j)+bet(i,k)*bet(j,l)*tmp(k,l)
  4   continue
      angles(1)=asin(sphi)*rad2deg
      angles(2)=asin(scta)*rad2deg
      angles(3)=beta*rad2deg

!     returned angles are in degrees
!     next correct scaler fluxes  
      return
      end
!===================================================================================================
!                                     end of subroutine rotate                                     !
!===================================================================================================



!===================================================================================================
!                                        subroutine covar                                          !
!===================================================================================================
!     subroutine for calculating covariance using detrended data.
!
!     uses u,v,w,t in common to calculate <u't'>, <v't'>, <w't'>  (t can be re-set beforehand to be 
!                                                                  any scalar).
!     ntossed inserted 5/9/02 to replace numnotused
!     (to allow for separate filtering on sonic and licor data)
!     and store result in cov(3) (note: uses rotated winds)
!     f9bus(nlag*2+1) is lag correlation summary across +/- nlag lags
!     lag (=0 for heat flux) is the lag between t and wind
!     gain = (default=1) the gain to convert 't' to correct units
!---------------------------------------------------------------------------------------------------
      subroutine covar(cov,f9bus,lag,gain,ntossed)
      real*8 wcavg, wclag, wcmax, wcstd, uc,vc,wc
      include 'lparams.for'

!...................................................................................................
!     note:  lparams includes in common:  gu,gv,gw,gt
!            (all read in from lsysdat.inp file, typically =1)
!...................................................................................................

      include 'lfluxvar.for'
      real*8 work(2*nlag+1)
      real*8 gain
      cov(1)= 0.d0
      cov(2)= 0.d0
      cov(3)= 0.d0
      wcavg = 0.d0
      wclag = 0.d0
      wcmax = 0.d0
      wcstd = 0.d0
      nlags = (nlag*2+1)
      do ilag =-nlag,nlag
      uc=0.d0
      vc=0.d0
      wc=0.d0
      lags  = ilag+lag                    ! offset this lag (lags) by ilag from center lag.
      nmlag = n-abs(lags)                 ! reducse number of samples by lag.
!     accumulates covariances (sum of product part):
      if (lags .lt. 0) then   
      do i=1,nmlag
      j   = i-lags
      uc  = uc+u(j)*t(i)
      vc  = vc+v(j)*t(i)
      wc  = wc+w(j)*t(i)
      end do 
      else
      do i=1,nmlag
      j  = i+lags
      uc = uc+u(i)*t(j)
      vc = vc+v(i)*t(j)
      wc = wc+w(i)*t(j)                   ! vertical flux
      end do 
      end if

!     normalize:
!     have now removed the number bad sonic points fromthe normalization. this could introduce a 
!     slight error given that when lags are figured, one covariance point is lost for each lag. if 
!     that point is a bad point, we'll have excluded it twice. max points of error are 6 out 14,400
!     so for ease of code, we've accepted a 0.1% potential error.

      uc=uc/dble(nmlag-ntossed)*gu*gain  ! 5/9/02 change: numnotused --> ntossed
      vc=vc/dble(nmlag-ntossed)*gv*gain
      wc=wc/dble(nmlag-ntossed)*gw*gain

!     matrix multiply: r-matrix by covariance vector
      if (ilag .eq. 0) then
      cov(1)=rot(1,1)*uc+rot(1,2)*vc+rot(1,3)*wc              ! rotated <u't'>
      cov(2)=rot(2,1)*uc+rot(2,2)*vc+rot(2,3)*wc              ! rotated <v't'>
      cov(3)=rot(3,1)*uc+rot(3,2)*vc+rot(3,3)*wc              ! rotated <w't'>
      end if
      wc=rot(3,1)*uc+rot(3,2)*vc+rot(3,3)*wc                  ! rotated w't' for all lags
      wcavg=wcavg+wc
      wcstd=wcstd+wc*wc
      work(ilag+nlag+1)=wc                                    ! store this <w't'>
      end do                                                  ! do ilag = -nlag, +nlag
      wcavg=wcavg/nlags
      delta=wcstd/nlags-wcavg*wcavg
      if (delta .le. 0.) wcstd=0.
      if (delta .gt. 0.) wcstd=sqrt(delta)
      wcmax=wcavg

!     finds maximum covariance magnitude:
      if (wcavg .ge. 0.) then

!     chooses most positive for max:
      do i=1,nlags
       if (wcmax .lt. work(i)) then
        wcmax=work(i)    

!     saves the lag at which maximum occurs:
        wclag=float(i-nlag-1)
       end if
205   end do
      else

!     chooses most negative for maximum:
      do i=1,nlags
       if (wcmax .gt. work(i)) then
        wcmax=work(i)
        wclag=float(i-nlag-1)
       end if
210   end do
      end if
      f9bus(1)=wcavg
      f9bus(2)=wclag
      f9bus(3)=wcmax
      f9bus(4)=wcstd
      return
      end
!===================================================================================================
!                                       end of subroutine covar                                    !
!===================================================================================================



!===================================================================================================
!                                         subroutine smear                                         !
!===================================================================================================
      subroutine smear(x,n,tau)
      real*8 x(n), tau
      real*8 alpha, alpha1
      alpha = exp(-2.5d-1/dble(tau))
      alpha1 = 1. - alpha
      x(1) = x(1)
      do i = 2, n
      x(i) = x(i-1)*alpha + x(i)*alpha1
      end do
      return
      end
!===================================================================================================
!                                      end of subroutine smear                                     !
!===================================================================================================



!===================================================================================================
!                                        subroutine read200                                        !
!===================================================================================================
!     dtime added for year-crossing check (srs, 20-jul-04)
!---------------------------------------------------------------------------------------------------
      subroutine read200( julday, gmt, v200, is200, done)
      include 'lparams.for'
      parameter(nzdim=nezdim)                   ! number of zeros for eddy zeroing.
      include 'lco2vars7.for'                   ! includes v200, is200.
      integer julday
      real*8 gmt 
      logical done
10    continue
      read(20,200, end=90, err=30) julday, gmt, v200, is200
200   format(i4,1x,f9.3,1x,2(f7.4,1x),4(f7.3,1x),f6.0,1x,i3,1x,i2)
      dtime = julday+gmt/86400.d0 - time0       ! time difference since previous line read at
                                                ! year-cross, dtime= 1.000-365.999 = -364.999 
                                                ! (or -365.999 in leap year).
                                                ! dtime in common block defined in lparams in
      time0 = julday+gmt/86400.d0               ! common block 'sys8byt'.
      return
!     error:
30    call gethhmmss(time0)
      write(99,105) julday, gmt,  time0, hh,mm,ss
      write(*, 105) julday, gmt,  time0, hh,mm,ss
105   format('  Warning: error in reading E*.200:  jday: ', i3, ', ',
     &       f9.3,' s' ,/, '  [ decimal jday = ', f8.4, ' (',i2.2,':',
     &       i2.2,':',f5.2,') ]' )
      goto 10                                   ! continues until good value is found.
!     end-of-file:
90    done = .true.
      end
!===================================================================================================
!                                    end of subroutine read200                                     !
!===================================================================================================


