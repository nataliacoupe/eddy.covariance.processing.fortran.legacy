!===================================================================================================
!                                     lba code: lcalpf6
!===================================================================================================   
!	created by scott saleska (started in july, 2000)
!	version 1.0 compiled 03-apr-01, saved as lcalpf.
!     adapted from xtlicor boreas code written by song-miao fan.
!
!	calibrates profile data. invoke 'lcalpf yymmdd'.
!	extracts zero and cal data, estimates calibration coefficients, and calculates concentrations 
!     for profile (level and time) data.
!
!	input data files:
!	.\yymmdd\split\pfyymmdd.201 	= scaled data for id=201
!	.\lsysdat.inp 			= general input parameters for lba processing
!
!	ouput data files:
! 	(1) calibration data
!	.\yymmdd\process\pfyymmdd.zer = zeros of profile instrument at each zero-time
!	.\yymmdd\process\pfyymmdd.cal = cal data (coefs, etc.) from each cal-time
!	(2) output data
!	.\yymmdd\process\pfyymmdd.co2 = concentration data calc from interpolated zeros and coefs
!	.\yymmdd\process\pfyymmdd.erc = error/log file for calibration
!	.\yymmdd\process\pfyymmdd.fin   (later) 
!
!	note:  input data structure is:
!	(1) a 20-min profile/hub sequence, broken up into
!	(2) 10 2-minute segments sampled at rate of freqlo (0.5 hz), broken into:
!	(3) a lag interval (tdelaypf secs for transition), and a sample interval (tavepf secs of good, 
!         stable data).
!
!     overall program steps:
!
!	-------------------------------------------------------------------------------
!	step 1: averages 0.5 hz raw data down to get one point for each sample interval
!	-------------------------------------------------------------------------------
!	(a-c) read all data from pfyymmdd.201, put averages of voltages from the sample interval
!           (i.e. not incl first delaypf secs) into these arrays:
!	- cal data:
!		vcal(ipcal,n), ipcal = 1 to nofpcal, n= order of cal-sequence, i.e.:
!		n=1: high-span,         	
!               n=2: mid-span,
!               n=3: low-span,   
!		n=4: dilution,
!               n=5: dilution (h2o)   
!		n=6: level 1 (air)
!               n=7: level 1 (h2o)
!		tcal(ipcal,n) = time of each of above
!		tcellc(), pcellc() = t, p of licor cell during calibration (from mid-span)
!     - zero data:
!		vzco2(izero), vzh2o(izero); izero = 1 to nofzero
!	- sample data
!		vpco2(iprfl), vph2o(iprfl); iprfl= 1 to nofprfl
!		tcprf(), pcprf() = t,p during sample
!	(d) writes zero-data to pfyymmdd.zer , including ready-light flag, zerflag:
!		(zerflag = 0 if all ok with licor during this zero, 1 if ready-light off for two 
!            sample points or more)
!
!	--------------------------------------------------------------------
!	step 2: estimates calibration coefficients from each set of cal data
!	--------------------------------------------------------------------
!
!	loop thru all cal sequence and for each, do:
!	(a) gets zeros corresponding to cal measurements  
!		- interpolates zero voltages --> zlow, zmid, zhigh
! 		- gets zero-adjusted span voltages dvlo, dvmid, dvhigh
!	(b) calculates cal coeffs, coef(i), i=1,2,3, put into coefx(ipcal), x=1,2,3 assuming:  
!           cref_xx = coef(1)*dv_xx + coef(2)*dv_xx^2 + coef(3)*dv_xx^3, wwhere xx is lo, mid, hi, 
!           and cref is the concentration of span gas, in ppm)
!	(c) calculates h2o, co2 concentrations for cal dilution, level1 samples
!        (diagnostics to check cal:  dilution = 1/2 level1 )
!	(d) writes cal-data to file pfyymmdd.cal (now including calflag:
!        (calflag =  0 if everything ok with licor
!         = xxxn: n=1: previous zero bad (next ok); 
!                 n=2: next zero bad (prev ok); 
!                 n=3: both bad
!         = xxnx: n=1 if low, high spans indistinguishable; 
!                 n=2 if time mismatch
!         = xnxx: n= # of cal-seq segments not long enough to average
!         = nxxx: n= # of span levels (1-3) during which licor ready-light if off allows for typical
!                    calflag
!        filter:  if calflag<3, is usable with one zero 
!
!	----------------------------------------------------------
!	step 3: calculates concentrations for each sample interval
!	----------------------------------------------------------
!
!    (a) interpolates zeros and coefficientss to each sampling measurement time
!    (b) calculate h2o (xph2o), co2 (xpco2) concentrations:
!   	- calculates h2o concentration from licor factory coefs
!    	- calculates co2 concentration using pressure, temp ratios relative to our cal, and h2o 
!       corrections due to pressure broadening and dilution) 
!    (c) writes concentration data to file pfyymmdd.co2 (now including liflag):
!        liflag =  0 if everything ok with licor
!               = xxn: n=1: previous zero bad (next ok); 
!                      n=2: next zero bad (last ok); 
!                      n=3: previous cal bad (next ok); 
!                      n=6: next cal bad (last ok)
!                     ( possible n:  1,2,3,4,5,6,7,8)
!               = xmx: m=1: both (or sole possible) zeros bad; 
!                      m=2: both (or sole possible) cals bad; 
!                      m=3 zeros and cals all bad
!               = 1xx: licor ready-light off for some portion of this measurement
!
!	typical filter to use with liflag:  if liflag<10, usable with at least one cal,zero
!---------------------------------------------------------------------------------------------------
!	note:  cal/zero sequence is as follows:
!	zero happens every 20 minutes each profile cycle (contrast w. every two hours in eddy system)
!	cal sequence happens every 6 hrs, w. following steps:
!                                                           cal bits ('iscal')
!		 0. zero (as end of previous profile sequence)           01
!		 1. high-span for 'tflstep' (120) seconds                01
!		 2. mid-span  for 'tflstep' (120) seconds                01
!		 3. low-span  for 'tflstep' (120) seconds                01
!		 4. zero                                                 01
!		 5. zero #2 (same for profile system)                    11
!		 6. zero #2 (or surv std if time -- 1 per week)          11 (01)
!		 7. dilution (500 sccm zero/ 500 sccm level 1)           10
!	 	 8. level 1 (to compare to dilution)                     00
!		 9. column integral  (to finish out sequence)            00
!		10. column integral                                      00
!       
!       nrc:    0. zero skiped during 2015 use 4. or 5.
!---------------------------------------------------------------------------------------------------
!	note: normal profile sequence is as follows (total 20 min):
!		 1. level 1 (top of tower)  (for 'tflstep' seconds)
!		 2. level 2
!		 3. level 3        	current level is in 9 status bits (16-24)
!		 4. level 4    
!		 5. level 5            	and
!		 6. level 6       
!		 7. level 7        	is converted to level number 1-9 in 'islvl'
!		 8. level 8
!		 9. column integral
!		10. zero
!---------------------------------------------------------------------------------------------------
!
!	input 'yymmddpf.201' data fields:
!	data:  13 variables
!	header has 4 lines:
!	 lba profile 201 data, 010126:  slow channels (1/2 hz), scaled to volts
!	julday, gmt(s),
!	pcon  flsam   flcal   pcell   rl      tcell   tdetect   tpump
!	co2   h2o     spare   spare   is201PF(24)
!
!	24 is201PF bits are:
!                 cal calgas hub levels
!	sid cal mux con zlmhs 1234 5678 col
!	011  00 000  00 00000 1000 0000 0
!
!	relevant status bits are converted in this program to:
!	islvl = level currently being sampled (0=none, 1-8 levels, 9=column)
!	iscal = cal status:  0= normal sample, 1= normal cal, 2=dilution, 3=zero #2
!	(for profile, zero #2 is same as zero; for eddy zero#2 is at high-flow)
!	ispan = cal gas status: -1= none, 0=zero air, 1=low-span, 2=mid-span, 3=hi-span
!-------------------------------------------------------------------------------------------
      program lcalpf 
      include 'lparams.for' 
!-------------------------------------------------------------------------------------------
!     lparams7 brings:
!	parameters for array-sizes/sys vars, including:
!    	parameters for co2 calibration: sized for 2-week data period :
!    	parameter (ncdim=(1440/tcalint)*14, nezdim=(1440/tedzint)*14 )
!    	parameter (npzdim=(1440/thubseq)*14 )
!	ncdim =size for cals:  1440/tcalint (def=4/day) * 14 days      	=   56 default
!	nezdim=size for eddy zeros: (1440 min/day)/(min/zero)[=12]*14 	=  168 default
!	npzdim=size for profile zeros: (1440/day)/(min/hubseq)*14 d   	= 1008 default
!       number of seconds per day 60*24=1440
!	profile: 
!     parameter (npwdim=3*freqlo*60*thubseq/10)  ! default=180
!    	sizes for working arrays (max data points per sampling segment)
!    	hub-sequence divided to 10 segments (def=2 min each) at freqlo (def=0.5)
!     normally, will sample 1 segment, but at end of cal-seq, will sample 2,
!     but use 3 to allow extra (default nwdim=3*.5*120=180 points)
!     parameter (npdim=60/(thubseq/10)* 24 * 14) ! default=10080
!     co2 sample:  --  arrays sized for 2-week data period
!                      30/hr * 24 * 14 days = 10080 (default) 
!	lparams7.for also defines, and puts into common-block 'sysvars':
!   	system variables (input from 'lsysdat.inp' file)
!	lparams7.for also defines data file/command line housekeeping vars
!    	(used in retrieving command-line info and file-name info
!-------------------------------------------------------------------------------------------

!	sets num zeros parameter for use in profile system:
      parameter(nzdim=npzdim)
      include 'lco2vars7.for' 
!-------------------------------------------------------------------------------------------
!     lco2vars6 bring:
!	variables for doing co2 concentration, including:
!	zeros:
!		real tzero(nzdim), vzco2(nzdim),vzh2o(nzdim)
!		real pcellz(nzdim), tcellz(nzdim)
!	cals:
!	real tcalave(ncdim), tcal(ncdim,7)  ! average time, specific times 1-7 as in vecal
!	real vcal(ncdim,11)  , pecell(3,ncdim), tecell(3,ncdim)
!   	vcal(n,*): average voltage of each segment (n) of each caluibration-sequence (*)
!			n: 	1  = ls;         2 = ms; 	3 = hs; 
!				4  = zero2; 	 5 = ss std; 	6 = dilution; 
!				7  = air (co2)	 8 = zero2;	9 = (dummy);
!				10 = dilution;	11 = air (h2o)
!        	real pcellc(ncdim,3), tcellc(ncdim,3)  ! cals (was p/tcal for lfluxes)
!		p/tcellc(*,n):  n:  1=main cal, 2=hi-flow zero, 3=dilution
!
!     -----------
!     local input arrays for profile data:  
!     (parts commented out are now included in lco2vars):
!     real*4 v201PF(12)       ! 201-input line voltages (see l1split for details)
!     integer*4 is201PF(24)   ! 201-input line status bits
!     logical ioldcal         ! check for new status in input data
!     sizes for working arrays (max data points per sampling segment):
!     real*4 work1(npwdim),work2(npwdim),work3(npwdim),work4(npwdim)
!     for aggregating co2, h2o, p, t
!     (npwdim= 2 120-sec/interval * 0.5 hz = 120, set by 'lparams.for' module)
!        
!     output arrays for co2 samples 
!     --  arrays sized for 2-week data period
!        30/hr * 24 * 14 days = 10080 (default) 
!-------------------------------------------------------------------------------------------
      real*8 pres, po, temp, to  
      integer*4 level(npdim)
      real*8  tprfl(npdim),vpco2(npdim),vph2o(npdim)  ! time, voltages
      real*8  pcprf(npdim),tcprf(npdim),fsampf(npdim),fcalpf(npdim)  ! pressure, temp
      integer lvlflag(npdim)                         ! licor ready-light status for profile levels
      logical done, debug 
      real*8  yrdays1, yrdays0, tdoy
      integer idate                                  !flag to carry on when profile collected every 7 or 8 sec
      real*8  wvsf      			     !water vapour caliration nrc 170127

!-------------------------------------------------------------------------------------------
!   step 0: set-up.  gets date for file to process, gets system parameters
!                    opens input & output files, writes headers
!                    gets date ('yymmdd') from command line for input data file
!-------------------------------------------------------------------------------------------
!     sets default initial parameters:
!
      lname  = 'lcalpf'
      ftypes = '3'
      include  'linitial7.for' !  common initialization code 

!---------------------------------------------------------------------------------------------------
!     linitial7 brings:
!     sets directory structure (root path, dirsplit, dirprocin, dirprocout), retrieves command line, 
!     file names for processing [ to ifilen(1..nifiles) ] and displays initial program message.
!---------------------------------------------------------------------------------------------------

      ico2 = 9         ! indices to v201PF for location of needed data
      ih2o = 10        ! (see module lsplit for details)
!     reads data from system parameter file
!===================================================================================================
      call readsysd(3) ! 3= read eddy system data (see 'lsubs.for' for code)
      tpfstep   = tflstep                	! nrc added 180820 as never defined before and needed to define calflag
!===================================================================================================
!     nrc from August 2015 to May 2016 the profile collected every 8 sec instead of 1 s
!     we modified the tavepf acordingly
!     if the SDM-INT is broken, the profile goes every 8 sec. God knows why!
      read (yymmdd,'(i10)') idate
      if  (((idate.ge.150101).and.(idate.lt.160601)).or.
     &   ((idate.ge.140820).and.(idate.le.140830)).or.
     &   ((idate.ge.130228).and.(idate.le.130401)).or.
     &   ((idate.ge.100817).and.(idate.le.100820)).or.
     &   ((idate.ge.171009).and.(idate.le.171030))) then
!       freqlo   = freqlo/4			! not 0.5 hz (2 sec) but 8 sec -can be changed as it is a parameter
       tavepf   = tavepf/4 	     		! all variables that depend on freqlo need to be recalculated
       tdelaypf = tavepf/4 	     	
       navepf   = navepf/4 
       ndelaypf = ndelaypf/4
       ndelayed = tdelayed*freqhi 		! number of points to wait before averaging
       print*, 'PF sampling every 8 sec at: ',idate
      end if

      if  (idate.lt.080101) then
       tpfstep  = 120.                   	! 2-minute cals 
       tdelaypf = 60
       tavepf   = 52.
       nflsteplo = tflstep*freqlo        	! number of points for each calibration step (low-freq)
       navecalo  = tavecal*freqlo        	! number of points to average for slow (201) eddy cal data 
       navepf    = tavepf*freqlo         	! number of points for one profile level average
       ndelaypf  = tdelaypf*freqlo
       navepf    = tavepf*freqlo         	! num for one profile level avg
       tdelayed  = tdelayed/86400.d0  
       tavepf    = tavepf/86400.d0
       tpfstep   = tpfstep/86400.d0
       print*,'  Initial setup 2000-2006: cal sequence was taking 2 min'
      endif

!............................................................................................

!     gets the full-path filenames:
!     raw data input (*.201 data)
      fnamein = ifilen(1)

!===========================================================================================
      call getfullpath( rootsplit,yymmdd, dirsplit, fnamein, lfnin) 
!===========================================================================================

!     calibration data (processed *.cal, *.zer files)
      fname = ifilen(1)

!===========================================================================================
      call getfullpath( rootprocin, yymmdd, dirprocout, fname, lfn)
!===========================================================================================

!...........................................................................................
!     note:
!     use of 'dirprocout' (will be '..in' from lfluxes)
!...........................................................................................

      write(* ,*) '  Opening input file ', fnamein(1:lfnin)//'.201'
      open(21,file=fnamein(1:lfnin)//'.201',status='old')  ! action='read')

!     opens output files
      open(12, file=fname(1:lfn)//'.zer',status='unknown')
      open(13, file=fname(1:lfn)//'.cal',status='unknown')
      open(14, file=fname(1:lfn)//'.co2',status='unknown')
      open(99, file=fname(1:lfn)//'.erc',status='unknown')  ! error file
!     note:
!     status unknown means any existing files will be over-written.
 
      write(99,55) platform, lname, iyear, imonth, iday,ihrs,imin,isec
55    format(a6,' ',a7,':  on date (yy-mm-dd) ', i4,'-',i2.2,'-',
     &             i2.2,' (',i2.2,':',i2.2,':',i2.2,'), process:' )
      write(99,*) '  ', fnamein(1:lfnin)//'.201'

!     puts headers in the output files

!     zero file 
      write(12,121) yymmdd,platform,lname,iyear,imonth,iday, ihrs,imin
121   format('2  10  LBA Profile zero, from: ', a6,' (',a6,' ',a7,
     &          ' on: ',i4,'-',i2.2,'-',i2.2,', ',i2.2,':',i2.2,')' )
      write(12,122) 
122   format('jdstart.gmt   doy   zco2.v   zh2o.v   pz.torr   ',
     &       'tzer.c   pcon.t   fsam   fcal   zerflag') 

!     calibration file 
      write(13,123) yymmdd,platform,lname,iyear,imonth,iday, ihrs,imin
123   format('5  35  LBA Profile cal, from: ', a6,' (',a6,' ',a7,
     &            ' on: ',i4,'-',i2.2,'-',i2.2,', ',i2.2,':',i2.2,')' )
      write(13,124) 
124   format('jdstart.gmt   doy   zlow.v   zmid.v   zhigh   ',
     &       'dvlow   dvmid   dvhigh   dvdil   dvair  ') 
      write(13,126) 
126   format('dvdilh   dvairh   pcal.t   pdil   pair   pcon   ',
     &               'tcal.c   tair   calflag   ') 
      write(13,127) 
127   format('fsamcal   fsamdil   fsamair   fcalcal   fcaldil   ',
     &       'fcalair   ') 
      write(13,128) 
128   format('coef1   coef2   coef3   co2d.ppm   co2air   ',
     &       'h2od.ppt   h2o.mmol.mol   co2.hs   co2.ms   co2.ls   ') 

!     co2 file 
      write(14,129) yymmdd,platform,lname,iyear,imonth,iday, ihrs,imin
129   format('4  23  LBA Profile CO2, from: ', a6,' (',a6,' ',a7,
     &          ' on: ',i4,'-',i2.2,'-',i2.2,', ',i2.2,':',i2.2,')' )
      write(14,131)  
131   format('jdstart.gmt   doy   level   z.v   dv   coef1   ',
     &       'coef2   coef3   zh2o   ')
      write(14,132)  
132   format('dvh2o   co2.ppm   h2o.ppt   pratio   tratio    ',
     &       'pair.t   fsam.sccm   flcal   liflag')  
      write(14,133)  
133   format('h2o.v    co2.v    po.tor    temp.C   to.C')  

!---------------------------------------------------------------------------------------------------
!     step 1:  extract all raw data from 201 file, average to get 1 value for each sample interval.
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!     step 1a. initialize for loop: read/store first data line from *.201 file
!---------------------------------------------------------------------------------------------------
!     sets counter indices so that arrays will load element 1 first time around
      iprfl = 1     !  counter for # profile samples
      ipcal = 1     !  counter for # profile calibrations
      izero = 1     !  counter for # zeros
      iseg  = 0     !  counter for segments read 
      licoroff = 0
      read(21,*)    !  skip over header lines
      read(21,*)    !  2
      read(21,*)    !  3
      read(21,*)    !  4
      dtype = "PF"
      debug = .false. 
      nsv   = 12
      nsi   = 24

!===================================================================================================
      call read20n(21,201,dtype,time0,v201PF,1,nsv,is201PF,1,nsi, 
     &             done,yrdays1,debug) 
!===================================================================================================
      tnext = time0

!===================================================================================================
      call readcal(time0) 
!===================================================================================================

201     format(i4,1x,f8.1,/,8f8.3,/,f9.4,3f8.4,1x,2(3i1,1x,2i1,1x),  
     &         5i1,1x,2(4i1,1x),i1) ! profile format 
!     must be identical to line 345 (?) of program lsplit (the format for creating the 201 file in 
!     the first place)
      iread = 1                            ! index for inputing data
!
!     store first line of data (raw voltages):
      work1(iread) = v201PF(ico2)          ! co2_profile       (iread=1 for these)
      work2(iread) = v201PF(ih2o)          ! h2o_profile
      work3(iread) = v201PF(itcell)        ! tcell_profile
      work4(iread) = v201PF(ipcell)        ! pcell_profile
      work5(iread) = v201PF(ifsam)
      work6(iread) = v201PF(ifcal)
      zerflag(izero) = 0                   ! izero=1            for telling when licor rl is off
      calflag(ipcal) = 0                   ! ipcal=1            counters
      lvlflag(iprfl) = 0                   ! iprfl=1

!===================================================================================================
      call gethhmmss(time0)
!===================================================================================================

      write(99,30) time0, int(gmt), hh,mm,ss
30    format('  First GMT time in file is ',f8.4,' (',i5,' secs) = ',
     &                  i3.2,':',i2.2,':',f5.2,' (hh:mm:ss.s)')
      write(99,31) fname(1:lfn),fname(1:lfn)
31    format('  Reading sequence of 2-min segments from file: ',/,
     &          '     - look for and store data from zeros and cals',/,
     &          '     - then write zero file: ',a,'.zer',/,
     &          '     - then write cal  file: ',a,'.cal',/)

!..........................................................................................nrc start
! after 170109 CC is wrong, needs to be the same as CAL PFxxxxx.201
!..........................................................................................
      if ((idate.ge.170109).and.(idate.lt.171109)) then	      
       if (is201PF(12).eq.0.) is201PF(12)=5.
       if ((is201PF(4).eq.0.).and.(is201PF(5).eq.1.).and.(is201PF(9).eq.0.).and.(is201PF(10).eq.1.)) then 
        is201PF(11)=1.
        do i=12,24
         is201PF(i) = 0.
        end do
       end if
       if (is201PF(12).eq.1.) is201PF(12)=0.
       if (is201PF(12).eq.5.) is201PF(12)=1.
       is201PF(9)  = 1*is201PF(4)
       is201PF(10) = 1*is201PF(5)
      end if
!............................................................................................nrc end

!    status info:
      islvl=max(is201PF(16)*1,is201PF(17)*2,is201PF(18)*3,is201PF(19)*4,
     &          is201PF(20)*5,is201PF(21)*6,is201PF(22)*7,is201PF(23)*8, 
     &          is201PF(24)*9)                          ! =1-8 levels, 9 if col, 0 if none

      ispan=-1+max(is201PF(11)*1,is201PF(12)*2,is201PF(13)*3,
     &             is201PF(14)*4,is201PF(15)*5)         ! =0 (zero), =1 (ls), =2(ms) =3(hs), =4(ss)
                                                        ! =-1 if no span gas is open
      
      iscal=(is201PF(9)*2 + is201PF(10)*1 )             ! = 0 if sampling, = 1 normal calibration
                                                        ! = 2 if dilution, = 3 if hi-flow zero (only
                                                        !                        relevant for eddy).
      
      ioldcal = iscal
      newstat = .false.

!---------------------------------------------------------------------------------------------------
!     step 1b: main loop (line 130) to read and average all segments in .201 file
!---------------------------------------------------------------------------------------------------

130   continue                           ! goes on to read and store average of next segment.

!     1b.i  loop to read and store one measurement segment in .201 file.
!           reads/stores in workn arrays until there is a change in status, indicating the end of 
!           this 2-min measurement segment.

!     continues reading data-file:
      dowhile(.not. newstat)             ! loop 
      time1 = time0  
!===================================================================================================
      call read20n(21,201,dtype,time0,v201PF,1,nsv,is201PF,1,nsi, 
     &                     done,yrdays1,debug) 
!===================================================================================================
      if (done) goto 135 
      dtime  = time0 - time1                                  ! time difference since previous line. 
                                                              ! dtime:   2.31481481E-05
!...................................................................................................
!      if ((idate.ge.170108).and.(idate.lt.170120)) then	      ! after 170109 CC is wrong, needs to be the same as CAL
!       if (is201PF(12).eq.0.) is201PF(12)=5.
!       if ((is201PF(4).eq.0.).and.(is201PF(5).eq.1.).and.(is201PF(9).eq.0.).and.(is201PF(10).eq.1.)) then 
!        is201PF(11)=1.
!        do i=12, 24
!         is201PF(i) = 0.
!        end do
!       end if
!       if (is201PF(12).eq.1.) is201PF(12) = 0.
!       if (is201PF(12).eq.5.) is201PF(12) = 1.
!       is201PF(9)  = 1*is201PF(4)
!       is201PF(10) = 1*is201PF(5)
!      end if
!...................................................................................................

      islvl0 = max(is201PF(16)*1,is201PF(17)*2,is201PF(18)*3,
     &             is201PF(19)*4,is201PF(20)*5,is201PF(21)*6,
     &             is201PF(22)*7,is201PF(23)*8,is201PF(24)*9) ! =1-8 levels, 9 if col, 0 if none
      ispan0 = -1+max(is201PF(11)*1,is201PF(12)*2,is201PF(13)*3,
     &                is201PF(14)*4,is201PF(15)*5)            ! =0 (zero), =1 (ls), =2(ms) =3(hs), =4(ss)      
!    where:
!    ispan = 0 (zero), 
!          = 1 (LS), 
!          = 2 (MS), 
!          = 3 (HS), 
!          = 4 (SS)
!          =-1 if no span gases open.
      iscal0 = (is201PF(4)*2 + is201PF(5)*1 )             

!     where:
!     iscal0 = 0 if sample, 
!            = 1 normal calibration,
!            = 2 dilution, 
!            = 3 hi-flow zero (only relevant for eddy)	      ! useful if the initial cal zero is missed 2015
      newstat=(islvl.ne.islvl0 .or. ispan.ne.ispan0 .or. iscal.ne.iscal0 
     &              .or. 1440*dtime.ge.tgap201)

!     dtime > tgap201 (minimum time to be considered a gap, in min)
      if(.not. newstat) then                       ! stores these voltages if status is unchanged.
       iread = iread+1
       work1(iread) = v201PF(ico2)
       work2(iread) = v201PF(ih2o)
       work3(iread) = v201PF(itcell)
       work4(iread) = v201PF(ipcell)
       work5(iread) = v201PF(ifsam)
       work6(iread) = v201PF(ifcal)
       work7(iread) = v201PF(ipcon)
       vrl = abs(v201PF(irl))                      ! check licor ready-light.
       if(vrl.lt.vminrl) licoroff = licoroff + 1
       if (iread.ge.npfdim) newstat = .true.       ! make sure array-size not exceeded.
      endif
      end do                                       ! until ioldcal.

!     1b.ii: averages last part of workn arrays down to one point each, allowing for 'tdelaypf' 
!            seconds of flushing period before averaging. stores averages in data arrays according
!            to status.
!
      iseg   = iseg+1   ! number of the segment just completed.
      tstart = tnext    ! start time of segment just completed.
      tnext  = time0    ! time of last point read = first point of next segment just starting.

!===================================================================================================
      call gethhmmss(tstart)
!===================================================================================================

!     check number of points in each segment.
!     sometimes there are too many points:
      if (iread .ge. npfdim) then
      write(99,50) iseg,tstart,int(mod(tstart,1d0)*86400),hh,mm,ss,
     &                      iscal, iread
50    format('  Error in segment #',i4,':  time = ',f8.4,' (',i5,
     &           ' sec) = ',i3.2,':',i2.2,':',f5.2,'; calstatus=',i1,
     &           '):',/,'      -->  ',i4,' values in segment ', 
     &           '  -- "workn" array size-limit reached !') 
      endif

!...................................................................................................
!     sometimes there is an insufficient number of points:  
      if (iread .lt. navepf+ndelaypf+1) then
       print*,'Insufficient number of points for PF calibration:', iread, ', should be: ',navepf+ndelaypf+1
       print*,'navepf: ', navepf, ', ndelaypf: ',ndelaypf
       write(99,57) iseg,tstart,int(mod(tstart,1d0)*86400),hh,mm,ss,
     &                      iscal, iread, navepf+ndelaypf+1
57     format('  Error in segment #',i4,':  time =',f8.4,' (',i5,
     &           ' sec) = ',i3.2,':',i2.2,':',f5.2,'; calstatus=',i1,
     &           '):',/, '      -->  ',i3,' values -- fewer data than ',
     &           i2,' required !') 
!   ir
       if(iscal.eq.1) then                         ! if this is a calibration or zero, can't use it.
        if((islvl.eq.0).and.(ispan.eq.0)) then     ! if zero (or zero #2)
         write(99,*) '      -->  this is a zero; disallow this zero'
        elseif(ispan.ge.1) then
         write(99,*) '      -->  this is a cal; disallow this cal'
         calflag(ipcal) = calflag(ipcal) + 100

!     100s place of calflag = n that insufficient points collected needed for calibration because 
!     there are several segments in cal-sequence, so even though this segment may be ignored, it 
!     contaminates the whole sequence.

        endif
       endif
      else                        ! proceeds if number of points in this segment is sufficient for 

       iend = (iread-1)           ! measurement, but makes end one less to be 1 away from transition.

       istart = iend-navepf+1     ! backup over navepf number needed for average.
       nave   = iend-istart+1     ! note: total points=iend-istart+1=navepf
       avg1 = 0                       
       avg2 = 0
       avg3 = 0
       avg4 = 0
       avg5 = 0
       avg6 = 0
       avg7 = 0
       do i =  istart, iend       ! calculates the average of the navepf.
        avg1 = avg1+work1(i)
        avg2 = avg2+work2(i)
        avg3 = avg3+work3(i)
        avg4 = avg4+work4(i)
        avg5 = avg5+work5(i)
        avg6 = avg6+work6(i)
        avg7 = avg7+work7(i)
       end do 
       avg1 = avg1/float(nave)    ! points ('tavepf' time)
       avg2 = avg2/float(nave)
       avg3 = avg3/float(nave)
       avg4 = avg4/float(nave)
       avg5 = avg5/float(nave)
       avg6 = avg6/float(nave)
       avg7 = avg7/float(nave) 
!     stores the averaged data in the right place:
!     zero:  no level sampled, and span gas = 0 
!     n.b.:  iscal can be zero or nonzero for profile zero

      if ((islvl.eq.0).and.(ispan.eq.0)) then      ! if zero (or zero #2).
       tzero(izero)  = tnext-tavepf/2.             ! mid-point of averaging period.

       vzco2(izero)  = avg1
       vzh2o(izero)  = avg2                        ! zero#2 is same for profile system.
       tcellz(izero) = avg3*gtcell
       pcellz(izero) = (avg4-zpcell)*gpcell
       flsamz(izero) = (avg5-zflsam)*gflsam
       flcalz(izero) = (avg6-zflcal)*gflcal
       prconz(izero) = (avg7-zpcon )*gpcon

!===================================================================================================
       call gethhmmss(tstart)
!===================================================================================================

       write(99,110) izero, iseg, iread, tstart,
     &              int(mod(tstart,1d0)*86400),hh,mm,ss
110    format('  Zero # izero=', i4,' (segment ',i4, ', ',i4,
     &       ' points.  time(start)= ', 
     &        f9.4, ' (',i5,' sec) = ',i3.2,':',i2.2,':',f5.2,')')
       if(licoroff.gt.2) then
        write(99,115) 'zero   '
115     format('  Warning:  LiCor ready-light is off during this ',a7)
        zerflag(izero) = zerflag(izero)+1               ! zero flag set: licor is off.
!     n.b.:  zeros stand alone, zerflag should be zero before this.
       endif
      endif
!     cal:  iscal=1, and span set at some level
      if (iscal.eq.1 .and. ispan.eq.3) then              ! if high-span
       tcal(ipcal,ihs) = tnext-tavepf/2.                 ! time = middle of average period.
       vcal(ipcal,ihs) = avg1                            ! stores average
       if(licoroff.gt.2) then
        write(99,115) '  HS calibration - sets calflag to "licor off"'
        calflag(ipcal) = calflag(ipcal) + 1000           ! sets calflag to "licor off".
       endif
      endif 
      if (iscal.eq.1 .and. ispan.eq.2) then              ! if middle-span
       tcal(ipcal,ims) = tnext-tavepf/2.              
       vcal(ipcal,ims) = avg1                        
       tcellc(ipcal,1) = avg3*gtcell                     ! p and t during mid-span are
       pcellc(ipcal,1) = (avg4-zpcell)*gpcell            ! used for this calibration.
       flsam(ipcal,1)  = (avg5-zflsam)*gflsam
       flcal(ipcal,1)  = (avg6-zflcal)*gflcal
       prcon(ipcal)    = (avg7-zpcon )*gpcon
       if(licoroff.gt.2) then
        write(99,115) '  MS calibration - sets calflag to "licor off"'
        calflag(ipcal) = calflag(ipcal) + 1000         ! sets calflag to "licor off".
!     n.b.: increments (by 1000's) for each span that licor is off.
       endif
      endif 
      if (iscal.eq.1 .and. ispan.eq.1) then            ! if low-span (last in sequence)
       tcal(ipcal,ils) = tnext-tavepf/2.              
       vcal(ipcal,ils) = avg1
       tcalave(ipcal)  = 
     &       (tcal(ipcal,ils)+tcal(ipcal,ims)+tcal(ipcal,ihs))/3.

!===================================================================================================
       call gethhmmss(tstart)
!===================================================================================================

       write(99,120) ipcal, iseg, iread, tstart,
     &              int(mod(tstart,1d0)*86400),hh,mm,ss
120    format('   -> cal # ipcal=', i4,' (low-span segment= ',i4, 
     &             ', ',i4,' points. time(startls)= ', f9.4, ' (',
     &               i5,' sec) = ',i3.2,':',i2.2,':',f5.2,')')
      if(licoroff.gt.2) then
       write(99,115) '  LS calibration - sets calflag to "licor off"'
       calflag(ipcal) = calflag(ipcal) + 1000  ! sets calflag to "licor off"
!     n.b.:  calflag could have been set to 100 during high or mid-spans if
!     one of those segments was of insufficient length
       endif
      endif 
      if (iscal.eq.2 .and. ispan.eq.0 .and. islvl.ge.1) then      ! if dilution.
       tcal(ipcal,idil)  = tnext-tavepf/2.              
       vcal(ipcal,idil)  = avg1                        
       vcal(ipcal,idilh) = avg2
       tcellc(ipcal,3)   = avg3*gtcell                            ! p and t during dilution are
       pcellc(ipcal,3)   = (avg4-zpcell)*gpcell                   ! used for this calibration.
       flsam(ipcal,3)    = (avg5-zflsam)*gflsam
       flcal(ipcal,3)    = (avg6-zflcal)*gflcal
      endif 
      if (iscal.eq.0 .and. islvl.gt.0 .and. ioldcal.eq.2 ) then     ! if cal-air

!     if this one is normal, and previous calibration was 2 (binary 10 = dilution)
!     then store in vcal array (nb: it will *also* store as profile data, below.)
       tcal(ipcal,iair)  = tnext-tavepf/2.              
       vcal(ipcal,iair)  = avg1                        
       vcal(ipcal,iairh) = avg2
       tcellc(ipcal,4)   = avg3*gtcell                      ! p and t during air sample are
       pcellc(ipcal,4)   = (avg4-zpcell)*gpcell             ! used for this calibration.
       flsam(ipcal,4)    = (avg5-zflsam)*gflsam
       flcal(ipcal,4)    = (avg6-zflcal)*gflcal
      endif 
      if (iscal.eq.0 .and. islvl.gt.0) then                ! if profile level
       tprfl(iprfl)  = tnext-tavepf/2.                     ! (nb: level > 0 prevents confusing
       level(iprfl)  = islvl                               ! with zeros at end of each 
       vpco2(iprfl)  = avg1                                ! hub-sequence)
       vph2o(iprfl)  = avg2
       tcprf(iprfl)  = avg3*gtcell
       pcprf(iprfl)  = (avg4-zpcell)*gpcell
       fsampf(iprfl) = (avg5-zflsam)*gflsam
       fcalpf(iprfl) = (avg6-zflcal)*gflcal
       if(licoroff.gt.2) then
        write(99,115) '  Level '//char(islvl+48)    ! level number converted to char via ascii
        lvlflag(iprfl) = lvlflag(iprfl) + 100       ! sets lvlflag to "licor off".
!     n.b.: levels stand alone, lvlflag should be zero before this.
       endif
      endif
      if (iscal.eq.1 .and. ispan.eq.4) then         ! if surveillance standard.
       tprfl(iprfl)  = tnext-tavepf/2.      
       level(iprfl)  = 0                            ! treats like profile level 0
       vpco2(iprfl)  = avg1
       vph2o(iprfl)  = avg2
       tcprf(iprfl)  = avg3*gtcell
       pcprf(iprfl)  = (avg4-zpcell)*gpcell
       fsampf(iprfl) = (avg5-zflsam)*gflsam
       fcalpf(iprfl) = (avg6-zflcal)*gflcal
       call gethhmmss(tstart) 
       write(99,125) iprfl, iseg,iread,tstart,int(mod(tstart,1d0)*86400),
     &              hh,mm,ss
125    format('  Surveilance standard (level=0) in iprfl=', i4,
     &       ' (segment ',i4, ', ',i4,' points. time(start)= ', 
     &       f9.4,' (',i5,' sec) = ',i3.2,':',i2.2,':',f5.2,')')
       if(licoroff.gt.2) then
        write(99,115) '  Survstd'
        lvlflag(iprfl) = lvlflag(iprfl) + 100       ! sets lvlflag to "licor off"
!     n.b.: surveillance standards stand alone, lvlflag should be zero before this.
       endif
      endif
!     counters are of single-segment samples, so increments are based on status of current segment:
      if (iscal.eq.0 .and. islvl.gt.0 .and. ioldcal.eq.2 ) then      ! end of calibration-sequence.
!     ii iscal=0 and islvl>0 = normal air sample, but last iscal (ioldcal)=2 
       ipcal = ipcal+1   
       calflag(ipcal) = 0                           ! sets calflag to ok for next calibration.
      endif
      if (ispan.eq.0 .and. islvl.eq.0) then         ! zero segment
       izero = izero+1                              ! span-gas = zero and no level
       zerflag(izero) = 0
      endif
      if (iscal.eq.0 .and. islvl.gt.0) then         ! profile level segment; i.e., if 
                                                    ! normal sample and at some level.
       iprfl = iprfl+1   
       lvlflag(iprfl) = 0
      endif
      if (ispan.eq.4) then                          ! surveillance standard segment.
       iprfl = iprfl+1   
       lvlflag(iprfl) = 0
!     if this is surveillance standard, also increment profile counter:
      endif
      if( ipcal.gt.ncdim .or. izero.gt.nzdim .or. iprfl.gt.npdim ) then
       write(* ,*)'  Reached array-size limit for cals, zeros, or data.'
       write(* ,*)'  Ending data analysis before end of file ',
     &            fnamein(1:lfnin)//'.201'
       write(99,*)'  Reached array-size limit for cals, zeros, or data.'
       write(99,*)'  Ending data analysis before end of file ',
     &            fnamein(1:lfnin)//'.201'
       goto 140
      endif
      endif      ! end else: if there are sufficient data points read into, takes average.

!     1b.iii: reinitializes work arrays (1 to 4) for the new intake period.
!             goes to 130 to do the next iteration of loop.

!     check time to next segment:
      if(1440*dtime .gt. tgap201) then

!===================================================================================================
      call gethhmmss(dtime)
!===================================================================================================

      write(99,40) int(86400*dtime), hh,mm,ss
40         format('  -------------',/,'   warning: ',i6,' sec = ',i3.2,
     &            ':',i2.2,':',f5.2,
     &             ' (hh:mm:ss) time-gap in *.201 file:')

!===================================================================================================
      call gethhmmss(tstart)
!===================================================================================================

      write(99,42) iseg,tstart,int(mod(tstart,1d0)*86400), hh,mm,ss, 
     &             iread,iscal
42    format('  Last segment #',i4,':  start-time = ',f8.4,' (',i5,
     &       ' sec) = ',i3.2,':',i2.2,':',f5.2,
     &       ' (', i4,' points, calstatus=',i1,')')

!===================================================================================================
      call gethhmmss(tnext)
!===================================================================================================

      write(99,44) iseg+1, tnext, int(mod(tnext,1d0)*86400), hh,mm,ss
44    format('  Next segment #',i4,':  start-time = ',f8.4,' (',i5,
     &       ' sec) = ',i3.2,':',i2.2,':',f5.2,
     &             /,'   -------------')
      if(1440*dtime .ge. 6) then 
!      if long-gap (long enough for full-cal sequence), wipes out any partially done cal-sequence.
       calflag(ipcal) = 0
       tcalave(ipcal) = 0
      endif
      endif
!     reinitializes for next segment
      iread=1
      work1(iread) = v201PF(ico2)
      work2(iread) = v201PF(ih2o)
      work3(iread) = v201PF(itcell)
      work4(iread) = v201PF(ipcell)
      work5(iread) = v201PF(ifsam)
      work6(iread) = v201PF(ifcal)
      work7(iread) = v201PF(ipcon)
      licoroff = 0
      ioldcal = iscal     ! preserves old cal status (for testing for cal-seq air).
      iscal   = iscal0    ! resets current status to be the new status.
      islvl   = islvl0
      ispan   = ispan0
      newstat = .false.
      goto 130            ! goes back and read a new data-segment.

!-------------------------------------------------------------------------------------------
!     step 1c: done reading 201 file. tie-up loose ends.
!---------------------------------------------------------------------------------------------------
135   continue            ! jumps here from read near line 130 if file at end.
      write(*,*)  '  Reached end of input file ', 
     &           fnamein(1:lfnin)//'.201'
      write(99,*) '  Reached end of input file ', 
     &           fnamein(1:lfnin)//'.201'
140   continue

      close(21)           ! closes the input file.
!
      nofpcal = ipcal-1   ! sets size of arrays containing averaged data.
      nofzero = izero-1
      nofprfl = iprfl-1
      write(*,142) iseg, nofzero, nofpcal, nofprfl
      write(99,142) iseg, nofzero, nofpcal, nofprfl
142   format('  Processed ',i5,' segments: ',i4,' zeros, ',i3,' cals, ',
     &       i5,' profile samples') 

!===========================================================================================
      call gethhmmss(tnext)
!===========================================================================================

      write(*,145) tnext, hh, mm, ss
      write(99,145) tnext, hh, mm, ss
145   format('  End time: ',f12.4,' (of last finished segment: ',i2.2,
     &       ':',i2.2,':',f5.2 ,').')

!-------------------------------------------------------------------------------------------
!     step 1d: writes zero-data to the output file.
!-------------------------------------------------------------------------------------------
      do izero = 1,nofzero
       tdoy = tzero(izero) - days00 - yrdays1  
       !nrc added set boundaries of what is possible given the write format constrains
       if ((vzco2(izero).lt.-9.).or.(vzco2(izero).gt.99.))   	vzco2(izero)  = -9.9999	!f7.4
       if ((vzh2o(izero).lt.-9.).or.(vzh2o(izero).gt.99.))   	vzh2o(izero)  = -9.9999	!f7.4
       if ((pcellz(izero).lt.-999.).or.(pcellz(izero).gt.999.))  pcellz(izero) = -999.99	!f7.2
       if ((tcellz(izero).lt.-999.).or.(tcellz(izero).gt.999.))  tcellz(izero) = -999.99	!f7.2
       !discard unrealistic dates and times
       if ((tdoy.ge.1.).and.(tdoy.lt.367.).and.(tzero(izero).gt.0.).and.(tzero(izero).lt.1000000.)) then		
      !..................................................................................
        write(12,212) tzero(izero), tdoy, vzco2(izero), vzh2o(izero), 
     &              pcellz(izero), tcellz(izero), prconz(izero), 
     &              flsamz(izero), flcalz(izero), zerflag(izero)
212     format(f10.4,1x,f8.4,2(1x,f7.4),1x,3(1x,f7.2),2(1x,f7.1),1x,i3)
       end if 		                				!no more crazy dates
      end do
      close(12)
      write(99,*) '  Zero file writen on ', fname(1:lfn)//'.zer' 
      write(*, *) '  Zero file writen on ', fname(1:lfn)//'.zer' 

!-------------------------------------------------------------------------------------------
!     step 2: estimates calibration coefficients from fitting calibration data to reference
!             gas concentrations.  write results to file pfyymmdd.cal.
!             includes updating calflag to indicate licor status:
!             calflag =  0 if everything ok with licor
!             = xxxn: n=1: previous zero bad; n=2: next zero bad; n=3: both bad
!             = xxnx: n=1 if low & high spans indistinguishable; n=2 if times don't match
!             = xnxx: n= # of cal-seq segments that were not long enough to average
!             = nxxx: n= # of span levels (1-3) during which licor ready-light off
!             these last two would have been set in previous section.
!-------------------------------------------------------------------------------------------
      jlo=1
      write(99,*)'  Fitting calibration coefficients...'
      write(*, *)'  Fitting calibration coefficients...'
      do ipcal=1,nofpcal
      time0 = dble( tcalave(ipcal) )

!-------------------------------------------------------------------------------------------
!      step 2a:  get the zeros for this calibration
!-------------------------------------------------------------------------------------------
      izero = 0                            ! nrc 170207 added as well as next if condition initialize

      if (nofzero.gt.1)       izero = ihuntf(tzero,nofzero,time0,jlo) ! izero = index to last tzero before time0
      if (izero.eq.0) then                 ! before the range of recorded zeros.
       izero = 1
       v0hs  = vzco2(izero)
       v0ms  = v0hs
       v0ls  = v0hs

       if (zerflag(izero).gt.0) calflag(ipcal) = calflag(ipcal) + 3
        print*,'Profile calibration unable to extract nearby zero'   
!      setting calflag=3 means: unable to extract nearby zero
!      (reminder: zerflag =0 if zero ok, 1 otherwise)

       elseif (izero .eq. nofzero) then    ! after the range of recorded zeros.
        v0hs = vzco2(nofzero)
        v0ms = v0hs
        v0ls = v0hs
       if (zerflag(izero).gt.0) calflag(ipcal) = calflag(ipcal) + 3
        print*,'Profile calibration unable to extract nearby zero within the range'
!       calflag=3 : unable to extract nearby zero within the range.

      else
       jlo   = izero
       delth = tcal(ipcal,ihs)-tzero(izero)  ! hs interp
       deltm = tcal(ipcal,ims)-tzero(izero)  ! ms interp
       deltl = tcal(ipcal,ils)-tzero(izero)  ! ls interp
!     now get slope: 
       deltat = tzero(izero+1)-tzero(izero)
       deltav = vzco2(izero+1)-vzco2(izero)
       dvdt   = deltav/deltat
       if (zerflag(izero).eq.0) then 
!       if previous zero good, check the next one

        if (zerflag(izero+1).ne.0) then  
!       previous zero is ok and the next is bad, don't interpolate.
         delt  = 0                               ! setting delt=0 prevents interpolation below.
         delth = delt
         deltm = delt
         deltl = delt
         calflag(ipcal) = calflag(ipcal) + 2     ! 2 = next zero is bad
        print*,'Profile calibration previous zero is ok and the next is bad'
        endif
!       co2 zero.
        v0hs = vzco2(izero)+dvdt*delth
        v0ms = vzco2(izero)+dvdt*deltm
        v0ls = vzco2(izero)+dvdt*deltl

       elseif(zerflag(izero+1).eq.0) then 
!     if previous zero is bad, but next one ok, use next co2 zero
        v0hs = vzco2(izero+1)
        v0ms = vzco2(izero+1)
        v0ls = vzco2(izero+1)
        calflag(ipcal) = calflag(ipcal) + 1     ! 1 = previous zero was bad
       print*,'Profile calibration previous zero is bad'

       else
!     if both zeros are bad... interpolates anyway, to see what they are.
        v0hs = vzco2(izero)+dvdt*delth 
        v0ms = vzco2(izero)+dvdt*deltm
        v0ls = vzco2(izero)+dvdt*deltl
        calflag(ipcal) = calflag(ipcal) + 3     ! 3 = both zeros are bad
       print*,'Profile calibration both zeros are bad'
      endif
!     if the previous zero (izero) is bad, calflag = 1
!     if the next zero (izero+1) is bad, calflag = 2 (if both are bad, will get 3)
!     note: calflag already >= 1000 if off during this calibration)
      end if
!     if the previous zero (# izero) is bad, set calflag (1-> previous zero bad)

!-------------------------------------------------------------------------------------------
!     step 2b: gets the calibration curve coefficients for this calibration
!              (dv3x3 array is defined in lco2vars module)
!-------------------------------------------------------------------------------------------
!      print*,' v0hs:', v0hs, ' v0ms: ', v0ms,' v0ls: ', v0ls
      dv3x3(1,1) = vcal(ipcal,ihs)-v0hs         ! subtracts off drift-corrected zeros.
      dv3x3(2,1) = vcal(ipcal,ims)-v0ms
      dv3x3(3,1) = vcal(ipcal,ils)-v0ls
      do i=1,3                                   ! form matrix of voltages.
       dv3x3(i,2) = dv3x3(i,1)*dv3x3(i,1)        ! voltage squared.
       dv3x3(i,3) = dv3x3(i,2)*dv3x3(i,1)        ! voltage cubed.
       coef(i)    = 0                            ! initialize coefficient
      end do
      dvhilo = abs( dv3x3(1,1) - dv3x3(3,1) )    ! signal difference between hi and low
!     dvhilo is typically = ~1.3 for eddy licors, and ~0.4 for profile
      if(dvhilo.lt.0.1) then 
       calflag(ipcal) = calflag(ipcal)+10        ! we want calflag=xx1x if low and high spans are 
                                                 ! indistinguishable  .
       print*,'Profile calibration low and high spans are indistinguishable'
      endif
!nrc........................................................................................
!     when a missing tank is present and low and high spans are indistinguishable
      if(((cref(1).gt.990.).or.(cref(2).gt.990.).or.(cref(3).gt.990.))
     &            .and.(dvhilo.lt.0.1)) then
       dvhilo = 0.2			         ! forced to be ~ val reorded 160706
       calflag(ipcal) = calflag(ipcal)-10
       iread = navepf+ndelaypf+1
       print*,'Profile calibration no available tanks'
      endif
!nrc........................................................................................

      if(calflag(ipcal).lt.3) then  
!     only inverts matrix if calflag was not raised, or if raised only due to bad zeros.
!     calflag = 1 if previouse zero bad, 2 if next, 3 if no nearby zero accesible.
!...........................................................................................
!nrc  start.................................................................................
!nrc  cref concentration in ppm, if cref==999 cal gas is off fill the matrix with other 
!     available tanks. If any tank concentration is >990
      if((cref(1).gt.990.).or.(cref(2).gt.990.)
     &                    .or.(cref(3).gt.990.)) then
!MS and LS concentrations are greater than 990
       if((cref(1).lt.990.).and.(cref(2).gt.990.)
     &                    .and.(cref(3).gt.990.)) then
        coef(1) = cref(1)/dv3x3(1,1)
        coef(2) = 0.0
        coef(3) = 0.0
        print *,'  Exclusively HS calibration.'
       endif
!...........................................................................................
       if((cref(1).gt.990.).and.(cref(2).lt.990.)
     &                    .and.(cref(3).gt.990.)) then
        coef(1) = cref(2)/dv3x3(2,1)
        coef(2) = 0.0
        coef(3) = 0.0
        print *,'  Exclusively MS calibration.'
       endif
!...........................................................................................
       if((cref(1).gt.990.).and.(cref(2).gt.990.)
     &                    .and.(cref(3).lt.990.)) then
        coef(1) = cref(3)/dv3x3(3,1)
        coef(2) = 0.0
        coef(3) = 0.0
        print *,'  Exclusively LS calibration.'
       endif
!...........................................................................................
       if((cref(1).lt.990.).and.(cref(2).gt.990.)
     &                    .and.(cref(3).lt.990.)) then
        dv2x2(1,1) = dv3x3(1,1)              ! subtract off drift-corrected zeros.
        dv2x2(2,1) = dv3x3(3,1)
        do i=1,2
         dv2x2(i,2) = dv2x2(i,1)*dv2x2(i,1)   ! voltage squared.
         coef2x2(i) = 0
        end do					
        cref2x2(1) = cref(1)
        cref2x2(2) = cref(3)
        call solvetwo(dv2x2, coef2x2, cref2x2)
        coef(1) = coef2x2(1)
        coef(2) = coef2x2(2)
        coef(3) = 0.0
!        print *,'coef(1)',coef(1),'coef(2)',coef(2),'coef(3)',coef(3)
        print *,'  PF Exclusively HS and LS calibration.'
       endif
!...........................................................................................
       if((cref(1).gt.990.).and.(cref(2).lt.990.)
     &                    .and.(cref(3).lt.990.)) then
        dv2x2(1,1) = dv3x3(2,1)           ! subtract off drift-corrected zeros.
        dv2x2(2,1) = dv3x3(3,1)
       do i=1,2
        dv2x2(i,2) = dv2x2(i,1)*dv2x2(i,1)  ! voltage squared.
        coef2x2(i) = 0
       end do
        cref2x2(1) = cref(2)
        cref2x2(2) = cref(3)
        call solvetwo(dv2x2, coef2x2, cref2x2)
        coef(1) = coef2x2(1)
        coef(2) = coef2x2(2)
        coef(3) = 0.0
        print *,'  PF Exclusively MS and LS calibration.'
       endif
!...........................................................................................
       if((cref(1).lt.990.).and.(cref(2).lt.990.)
     &                    .and.(cref(3).gt.990.)) then
        dv2x2(1,1) = dv3x3(1,1)           ! subtract off drift-corrected zeros.
        dv2x2(2,1) = dv3x3(2,1)
        do i=1,2
         dv2x2(i,2) = dv2x2(i,1)*dv2x2(i,1)         ! voltage squared.
         coef2x2(i) = 0
        end do
        cref2x2(1) = cref(1)
        cref2x2(2) = cref(2)
        call solvetwo(dv2x2, coef2x2, cref2x2)
        coef(1) = coef2x2(1)
        coef(2) = coef2x2(2)
        coef(3) = 0.0
        print *,'  PF Exclusively HS and MS calibration.'
       endif
!...........................................................................................
       if((cref(1).gt.990.).and.(cref(2).gt.990.)
     &                    .and.(cref(3).gt.990.)) then
        call readco2cal(time0,'PF') 
        coef(1) = fcoef(1)
        coef(2) = fcoef(2)
        coef(3) = fcoef(3)
        print *,'  No cal tanks -borrowed coefficients: ', coef
       end if
      else				! all tanks available
!===========================================================================================
      call solve(dv3x3, coef, cref) 
!===========================================================================================

!     problem:  dv(3x3) * coef(1:3) = cref(1:3)
!     solution: coef = dv^-1 * cref
!     note: order is 1=hs, 2=ms, 3=ls


      endif				! end if when at least one cal tank is missing
!nrc end.............................................................................................
      endif				! end calflag(iecal).lt.3
      coef1(ipcal) = coef(1)		! save coefficients for later use
      coef2(ipcal) = coef(2)
      coef3(ipcal) = coef(3)
!===========================================================================================
      call gethhmmss(time0) 
!===========================================================================================

       write(99,510) ipcal, time0, int(mod(time0,1d0)*86400), hh,mm,ss, 
     &               coef, calflag(ipcal)
510   format('  Calibration # ',i2,' at t =', f9.4, ' (',i5,' sec) = ',
     &        i3.2,':',i2.2,':',f5.2,
     &       ', coefs: ',3(1x,1pe11.4), ' (calflag =',i4,')')  ! 0pi?

!     checks is timebase for this calibration is consistent:
      deltat=abs(tcal(ipcal,ils)-tcal(ipcal,ihs))
!       tpfstep = tpfstep*1.5  ! change due to 3-minute cals oct 2006. nrc 180516 commented 
                                          ! as tpfstep recalculated at lcaed  
      if (deltat .gt. (3.*tpfstep)) then  ! deltat should ~= 2*tpfstep.
      calflag(ipcal) = calflag(ipcal)+20  ! we want calflag = xx2x if times do not match.
       print*,'Profile calibration if times do not match'
       print*,' check parameters deltat: ', deltat, '     tpfstep: ', tpfstep

!===========================================================================================
      call gethhmmss( dble(tcal(ipcal,ils)) )
!===========================================================================================

      hhl = hh
      mml = mm
      ssl = ss

!===========================================================================================
      call gethhmmss( dble(tcal(ipcal,ihs)) )
!===========================================================================================

      write(99,220) ipcal, tcal(ipcal,ils), int(mod(tcal(ipcal,ils),
     &              1.0)*86400), hhl,mml,ssl,tcal(ipcal,ihs), 
     &              int(mod(tcal(ipcal,ihs),1.0)*86400), hh,mm,ss
220   format('  Warning for calibration # ',i3, ':  time (low- span) =',
     &             f9.4, ' (',i5,' sec) = ',i3.2,':',i2.2,':',f5.2,';',
     &             /,'                            time (high-span) =', 
     &             f9.4, ' (',i5,' sec) = ',i3.2,':',i2.2,':',f5.2,')')
      write(99,221) 
221   format(12x,'  High-Span and Low-Span times do not match.') 
      endif

!-------------------------------------------------------------------------------------------
!     step 2c:  calculate concentrations of dilution & air measurements
!-------------------------------------------------------------------------------------------
      vcal(ipcal,idil) = vcal(ipcal,idil)-vzco2(izero+1)            
!     izero+1 = zero after 3-span set = zero just before dilution
      co2dil = vcal(ipcal,idil)
      vcal(ipcal,idilh) = vcal(ipcal,idilh)-vzh2o(izero+1)
      h2odil = vcal(ipcal,idilh)
      po   = pcellc(ipcal,1)     ! assumes this applies yo measurements as well as calibrations 
      to   = tcellc(ipcal,1)     ! (i.e., no p and t corrections applied).
      pres = pcellc(ipcal,3)
      temp = tcellc(ipcal,3)

!===========================================================================================
      call getconc(h2odil,co2dil,pres,po,temp,to,coef)
!===========================================================================================

      vcal(ipcal,iair) = vcal(ipcal,iair)-vzco2(izero+1)
      co2air = vcal(ipcal,iair)
      vcal(ipcal,iairh) = vcal(ipcal,iairh)-vzh2o(izero+1)
      h2oair = vcal(ipcal,iairh)
      pres   = pcellc(ipcal,4)      
      temp   = tcellc(ipcal,4)

!===========================================================================================
      call getconc(h2oair,co2air,pres,po,temp,to,coef)
!===========================================================================================

!-------------------------------------------------------------------------------------------
!     step 2d:  write results to the calibration file (pfyymmdd.cal)
!-------------------------------------------------------------------------------------------
      tdoy = tcalave(ipcal) - days00 - yrdays1  
      !nrc added set boundaries of what is possible given the write format constrains
      if ((dv3x3(3,1).lt.-9.9999).or.(dv3x3(3,1).gt.9.9999))   dv3x3(3,1)=9.9999	!f7.4
      if ((dv3x3(2,1).lt.-9.9999).or.(dv3x3(2,1).gt.9.9999))   dv3x3(2,1)=9.9999	!f7.4
      if ((dv3x3(1,1).lt.-9.9999).or.(dv3x3(1,1).gt.9.9999))   dv3x3(1,1)=9.9999	!f7.4
      if ((vcal(ipcal,idil).lt.-9.9999).or.(vcal(ipcal,idil).gt.9.9999))   vcal(ipcal,idil)=9.9999	!f7.4
      if ((vcal(ipcal,iair).lt.-9.9999).or.(vcal(ipcal,iair).gt.9.9999))   vcal(ipcal,iair)=9.9999	!f7.4
      if ((co2dil.lt.-999.).or.(co2dil.gt.9999.))   co2dil = -999.99	!0PF7.2
      if ((co2air.lt.-999.).or.(co2air.gt.9999.))   co2air = -999.99	!0PF7.2
      if ((h2odil.lt.-99.).or.(h2odil.gt.999.))     h2odil = -99.999	!f7.3
      if ((h2oair.lt.-99.).or.(h2oair.gt.999.))     h2oair = -99.999	!f7.3
      if ((vcal(ipcal,idilh).lt.-9.).or.(vcal(ipcal,idilh).gt.99.))     vcal(ipcal,idilh) = -9.9999	!f7.4
      if ((vcal(ipcal,iairh).lt.-9.).or.(vcal(ipcal,iairh).gt.99.))     vcal(ipcal,iairh) = -9.9999	!f7.4
      if ((tdoy.ge.1.).and.(tdoy.lt.367.)) then		                !no more crazy dates
      !..................................................................................
      write(13,213) tcalave(ipcal),tdoy,v0ls,v0ms,v0hs, dv3x3(3,1),
     &                    dv3x3(2,1),dv3x3(1,1), 
     &                    vcal(ipcal,idil),vcal(ipcal,iair)       	! low,mid,high,dil,air dv; 
      write(13,214) vcal(ipcal,idilh),vcal(ipcal,iairh),pcellc(ipcal,1),
     &                   pcellc(ipcal,3),pcellc(ipcal,4),prcon(ipcal),
     &                   tcellc(ipcal,1),tcellc(ipcal,4),calflag(ipcal)
      write(13,215) flsam(ipcal,1),flsam(ipcal,3),flsam(ipcal,4), 
     &                   flcal(ipcal,1),flcal(ipcal,3),flcal(ipcal,4)   ! 1=cal, 3=dil, 4=air
      write(13,216) coef,co2dil,co2air,h2odil,h2oair,                   ! coef = array of 3
     &                  cref(1),cref(2),cref(3)                         ! cal tank conc.   
!
213   format( f10.4,1x,f8.4,8(1x,f7.4) )
214   format( 2(1x,f7.4),1x,6(1x,f7.2),4x,i4 )
215   format( 6(1x,f7.1) )
216   format( 1x,3(1x,1pe11.4),2(1x,0PF7.2),2(1x,f7.3),3(1x,f7.2) )
      endif		                                                !endif no more crazy dates

!...........................................................................................
!     note on possible formatting confusion for 216-format line:
!     the format-editor 1pe11.4 was used to output cal-coefficients
!      ->  1p, when used with e, the base is multiplied by 10, the exponent reduced by 1;
!          [ the default format for e-editing is 0.xyzeab, this produces x.yze(ab-1) ]
!          but the 1p-setting remains for the whole line, unless changed.  
!      ->  1p, when used with f, multiplies the output by 10.
!          thus, without the 0p in front of f7.2, the other values would be multiplied by 10
!     (see p editing, watcom lang. ref. manual, p. 274)
!...........................................................................................

      end do             ! ipcal=1,nofpcal
      close(13)
      write(99,*) '  Calibration file writen on ', fname(1:lfn)//'.cal' 
      write(* ,*) '  Calibration file writen on ', fname(1:lfn)//'.cal' 

!-------------------------------------------------------------------------------------------
!     step 3:  calculates concentrations for the profile data.
!              writes the results to pfyymmdd.co2
!              includes updating liflag to indicate licor status:
!              liflag  =  0 if everything is ok with licor
!                      = xxxn: n=1: just previous zero bad; n=2: just next zero bad; 
!                              n=3: just previous cal bad;  n=6: just next cal bad
!                              ( possible n:  1,2,3,4,5,6,7,8)
!                      = xxmx: m=1: both (or sole possible) zeros bad; 
!                              m=2: both (or sole possible) cals bad; m=3 both zers and cals
!                      = x1xx: licor ready-light off for some portion of this measurement
!                 this allows filter:  if liflag<10, usable with at least one cal,zero
!-------------------------------------------------------------------------------------------
       jlo = 1
       klo = 1
!     loop through all profile measurements:
      do iprfl = 1,nofprfl
       time0   = tprfl(iprfl)

!---------------------------------------------------------------------------------------------------
!       step 3a.  interpolate zero data to the time of this profile measurement
!-------------------------------------------------------------------------------------------
      izero=0                                 ! nrc 170207 added
      if (nofzero>1)       izero = ihuntf(tzero,nofzero,time0,jlo)
      if (izero .eq. 0) then                  ! before the range of recorded zeros.
       v0h2o = vzh2o(1)
       v0co2 = vzco2(1)
       izero = 1
       jlo   = izero
      if (zerflag(izero).gt.0) lvlflag(iprfl) = lvlflag(iprfl) + 10  ! 10 if zero bad, 0 otherwise.
!     setting lvlflag=x1x means it was unable to extract nearby zero
!     (reminder:  zerflag =0 if zero ok, 1 otherwise)
!                 lvlflag = 100 already if licor went off during the level measurements.
      elseif (izero .eq. nofzero) then        ! after the range of recorded zeros
       v0h2o = vzh2o(nofzero)
       v0co2 = vzco2(nofzero)
      if (zerflag(izero).gt.0) lvlflag(iprfl) = lvlflag(iprfl) + 10  ! 10 if zero bad, 0 otherwise.
      elseif ((izero.gt.0) .and. (izero.lt.nofzero)) then            ! within range.
       jlo    = izero
       deltat = tzero(izero+1)-tzero(izero)
      if (zerflag(izero).eq.0) then 
!     if previous zero was good, check the next one.
      if (zerflag(izero+1).eq.0) then 
       delt = time0-tzero(izero)
      else                                    ! if next zero is bad, don't interpolate yet
       delt = 0                               ! setting delt=0 prevents interpolation below
       lvlflag(iprfl) = lvlflag(iprfl) + 2    ! 2=next zero bad
      endif
!     water zero
       deltav = vzh2o(izero+1)-vzh2o(izero)
       dvdt   = deltav/deltat
       v0h2o  = vzh2o(izero)+(dvdt*delt)
!     co2 zero
       deltav = vzco2(izero+1)-vzco2(izero)
       dvdt   = deltav/deltat
       v0co2  = vzco2(izero)+(dvdt*delt)
      elseif(zerflag(izero+1).eq.0) then 
!     if previous zero bad, but next one good, use next
       lvlflag(iprfl) = lvlflag(iprfl) + 1  ! 1=prev zero bad
!     water zero
       v0h2o = vzh2o(izero+1)
!     co2 zero
       v0co2 = vzco2(izero+1)
      else
!     if both are bad:
       lvlflag(iprfl) = lvlflag(iprfl) + 10       ! 10 = both zeros are bad
       v0h2o = vzh2o(izero)                       ! will not be used, but see what it is anyway
       v0co2 = vzco2(izero)
      endif
      end if
      if (vph2o(iprfl).ne.v0h2o) vph2o(iprfl) = vph2o(iprfl)-v0h2o          ! h2o zero is subtracted
      if (vpco2(iprfl).ne.v0co2) vpco2(iprfl) = vpco2(iprfl)-v0co2          ! co2 zero is subtracted

!---------------------------------------------------------------------------------------------------
!     step 3b:  interpolates calibration coefficients to this measurement
!---------------------------------------------------------------------------------------------------
      ipcal=0                                ! nrc 170207 added
      if (nofpcal>1)       ipcal = ihuntf(tcalave,nofpcal,time0,klo)  ! index to last tcal before time0
      if (ipcal .eq. 0) then
       coef(1) = coef1(1)
       coef(2) = coef2(1)
       coef(3) = coef3(1)
       ipcal   = 1
       klo     = ipcal
      if (calflag(ipcal).ge.3) lvlflag(iprfl) = lvlflag(iprfl) + 20
!     lvlflag xxmx, m = 2 if calibrations are not retrievable.
      else if (ipcal .eq. nofpcal) then
       coef(1) = coef1(nofpcal)
       coef(2) = coef2(nofpcal)
       coef(3) = coef3(nofpcal)
      if (calflag(ipcal).ge.3) lvlflag(iprfl) = lvlflag(iprfl) + 20
!     lvlflag xxmx, m = 2 if calibrations are not retrievable.
      else                    ! in the range (0 > ipcal > nofpcal)
       klo=ipcal
       deltat = tcalave(ipcal+1)-tcalave(ipcal)
      if (calflag(ipcal).lt.3) then
!     if previous calibration was good, check the next one.
      if (calflag(ipcal+1).lt.3) then
       delt = time0-tcalave(ipcal)
      else                    ! previous calibration was good, but next is bad, so don't interpolate.
       delt=0                 ! setting delt=0 prevents interpolation below
       lvlflag(iprfl) = lvlflag(iprfl) + 6          ! 6 = just next calibration is bad
      endif
       dcdt    = (coef1(ipcal+1)-coef1(ipcal))/deltat
       coef(1) = coef1(ipcal)+(dcdt*delt)
       dcdt    = (coef2(ipcal+1)-coef2(ipcal))/deltat
       coef(2) = coef2(ipcal)+(dcdt*delt)
       dcdt    = (coef3(ipcal+1)-coef3(ipcal))/deltat
       coef(3) = coef3(ipcal)+(dcdt*delt)
      elseif (calflag(ipcal+1).lt.3) then
!     previous calibration was bad, bit next is good, so use only next.
       coef(1) = coef1(ipcal+1)
       coef(2) = coef2(ipcal+1)
       coef(3) = coef3(ipcal+1)
       lvlflag(iprfl) = lvlflag(iprfl) + 3     ! 3 = just previous calibration was bad
      else 
       lvlflag(iprfl) = lvlflag(iprfl) + 20    ! 20 = both calibrations are bad
       coef(1) = coef1(ipcal)
       coef(2) = coef2(ipcal)
       coef(3) = coef3(ipcal)
       endif
      end if

!---------------------------------------------------------------------------------------------------
!nrc start: get coef from forced when no tanks available
      if((cref(1).gt.990.).and.(cref(2).gt.990.)
     &                    .and.(cref(3).gt.990.)) then
!==================================================================================================
       call readco2cal(time0,'PF') 
!==================================================================================================
       coef(1) = fcoef(1)
       coef(2) = fcoef(2)
       coef(3) = fcoef(3)
      end if
!nrc end
!-------------------------------------------------------------------------------------------
!     step 3c.  calculate concentrations from interpolated cal-curve
!-------------------------------------------------------------------------------------------
!     wvsf    = 1.d0 			!commented nrc 170125
!     gets water-cal scale factor
!..................................................................................nrc 170127
      call readh2ocal(time0,'PF')
!...........................................................................................
      wvsf = dble(cvapor)
!      print*, 'time0: ', time0, '  wvsf: ', wvsf

!     sets h2o and co2 voltages, call subroutine getconc2 to get the concentrations.
      pratio  = pcellc(ipcal,1)/pcprf(iprfl)  ! po/p
      tratio  = (273.0+tcprf(iprfl))/(273.0+tcellc(ipcal,1)) ! t/to
      xph2o   = vph2o(iprfl)
      xpco2   = vpco2(iprfl)
      xph2o_v = xph2o			!added to keep record of voltages nrc 170125
      xpco2_v = xpco2

!===========================================================================================
      call getconc2(xph2o,xpco2,pcprf(iprfl),pcellc(ipcal,1),
     &                          tcprf(iprfl),tcellc(ipcal,1),coef,wvsf)
!===========================================================================================
!     (?) question:  interpolate P0 and T0 in future? (probably not needed)

!-------------------------------------------------------------------------------------------
!     step 3d.  write the results
!-------------------------------------------------------------------------------------------
      tdoy = time0 - days00 - yrdays1  
      !nrc added set boundaries of what is possible given the write format constrains
      if ((v0co2.lt.-999.).or.(v0co2.gt.9999.))               v0co2 = -999.9999		!f9.4
      if ((vpco2(iprfl).lt.-999.).or.(vpco2(iprfl).gt.9999.)) vpco2(iprfl) = -999.9999	!f9.4
      if ((v0h2o.lt.-.9999).or.(v0h2o.gt.9.9999))             v0h2o = 9.9999		!f7.4
      if ((vph2o(iprfl).lt.-9.).or.(vph2o(iprfl).gt.99.))     vph2o(iprfl) = -9.9999	!f7.4
      if ((xpco2.lt.-999.).or.(xpco2.gt.9999.))               xpco2 = -999.99		!f7.2
      if ((xph2o.lt.-99.).or.(xph2o.gt.999.))                 xph2o = -99.999		!f7.3
      if ((pratio.lt.-99.).or.(pratio.gt.999.).or.(pcprf(iprfl).eq.0.)) pratio = -99.9999	!f8.4
      if ((tratio.lt.-99.).or.(tratio.gt.999.))               tratio = -99.9999		!f8.4
      if ((xph2o_v.lt.-9.).or.(xph2o_v.gt.99.)) xph2o_v=-9.9999	!f8.4
      if ((xpco2_v.lt.-9.).or.(xpco2_v.gt.99.)) xpco2_v=-9.9999	!f8.4
      if ((pcellc(ipcal,1).lt.-99.).or.(pcellc(ipcal,1).gt.999.)) pratio = -99.9999	!f8.4
      if ((tcprf(iprfl).lt.-9.).or.(tcprf(iprfl).gt.99.))         tcprf(iprfl) = -9.9999	!f8.4
      !no more crazy dates
      if ((tdoy.ge.1.).and.(tdoy.lt.367.).and.(time0.ge.1.).and.(time0.lt.100000.)) then
      !..................................................................................
      write(14,224) time0,tdoy,level(iprfl),v0co2,vpco2(iprfl),
     &              coef,v0h2o
      write(14,225) vph2o(iprfl),xpco2,xph2o,pratio,tratio, 
     &              pcprf(iprfl), fsampf(iprfl), fcalpf(iprfl), 
     &              lvlflag(iprfl)
      write(14,226) xph2o_v,xpco2_v,pcellc(ipcal,1),tcprf(iprfl),tcellc(ipcal,1)
224   format(f10.4,1x,f8.4,1x,i2,3x,2(1x,f9.4),1x,3(1x,1pe11.4),1x,f7.4)
225   format(1x,f7.4,1x,f7.2,1x,f7.3,1x,2(f8.4),1x,f6.1,1x,f8.1,3x,f8.1,
     &       1x,i4 )
226   format(1x,f10.7,1x,f10.7,1x,f9.4,1x,f9.6,1x,f9.6)
      endif		                                                !endif no more crazy dates
      end do
      close(14)
      write(99,*) '  CO2 data file writen on: ', fname(1:lfn)//'.co2' 
      write(* ,*) '  CO2 data file writen on: ', fname(1:lfn)//'.co2' 
      write(99,*)
      write(99,*) '  End of program LCALPF.'
      print *,    '  End of program LCALPF.'
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

