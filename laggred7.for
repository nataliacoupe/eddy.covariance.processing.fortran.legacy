!=======================================================================
!                        lba code:  laggred7.for          
!=======================================================================
!     created by     scott saleska  (started sept 2001)
!     lba aggregate eddy data
!     (start with just housekeeping)
!     laggr 001212 -files 34 -min 30 (20 min is default)
!     (20 min default corresponds to calibration-sequence length) 
!
!     Program to aggregate low-frequency data 0.5 Hz data to specified 
!     intervals.
!     Ideal: do 2-minute intervals during calibration/zero times, and do
!     standard aggregation elsewhere.
!     ------------------------------------------------------------------
      program laggred

!     Parameters for array-sizes + system variables:
      include 'lparams.for'   

!     Number of zeros for eddy zeroing:
      parameter(nzdim=nezdim)

!     Variables for doing CO2 concentrations:
      include 'lco2vars7.for'

!     Filtering and housekeeping variables:
      include 'lfluxvar.for'

      logical done, debug  
      real*8 yrdays1, tdoy  

!     Sets initial parameters, gets command-line, etc.:
      lname = 'laggred'
      ftypes='12'
      fnum=1
      include 'linitial7.for'  
!     'LINITIAL" is a common initialization code. It sets the directory 
!     structure (root paths, dirsplit, dirprocin, dirprocout), then
!     retrieves command line, file names for processing 
!     [ to ifilen(1..nifiles) ] and displays initial program message.
      if (navg.eq.0) navg = freqlo*60.*period

!     Form input file name(s):
      if(fnum.eq.1) then 
      dtype='EA'
      else
      dtype='EB'
      endif
      enum=fnum
      print *,'  Reading parameter file lsysdat.inp ...'
      call readsysd(enum)  ! read generic vars, plus eddy data

!     ------------------------------------------------------------------
      call readconv(enum)  ! read conversion data, from eddy part 
!     raw data input (*.201 data)
      fnamein = ifilen(1)
      call getfullpath( rootsplit,yymmdd, dirsplit, fnamein, lfnin) 
!     data output (*.hkp, *.ere)
      fnameout = ifilen(1)
      call getfullpath( rootprocout,yymmdd,dirprocout,fnameout,lfnout )
      print *, '  -> opening input files "'//fnamein(1:lfnin)//'"...'
      print *, '  -> sending output to   "'//fnameout(1:lfnout)//'"...'

!===================================================================================================
!                                  main program starts here    
!===================================================================================================
      read (yymmdd,'(i10)') idate		! nrc 180302 set calibration coefficient for humidity sensor
!     ------------------------------------------------------------------
!     first --  opens the input file (E*.201) and output files
!     ------------------------------------------------------------------

!     opens input file (hkp and chilled mirror data is in the .201)
      open(21,file=fnamein(1:lfnin)//'.201',status='old')
      call skipmirhead(21)    ! skip header of eddy 201 file

!     opens output files
      open(16,file=fnameout(1:lfnout)//'.hkp',status='unknown')         ! housekeeping
      open(99,file=fnameout(1:lfnout)//'.ere',status='unknown')         ! errorfile

!     watcom-specific date functions
      write(99,55) iyear, imonth, iday,ihrs,imin,isec
      write(99,*) 
55    format('  LAGGRED version 6:  processing date (yy-mm-dd): ', i4,
     &       '-',i2.2,'-',i2.2,' (',i2.2,':',i2.2,':',i3.2,')' )
      write(99,*) '  Opened input file ', fnamein(1:lfnin)//'.201'
      write(*,*)  '  Processing sample file ',fnamein(1:lfnin)//'.201'
      write(16,141) enum,yymmdd,period,iyear,imonth,iday, ihrs,imin

141   format('3 16 LBA Eddy',i1,' housekeeping from: ',a6,' (',f5.1,  
     &       ' min) (on ',i4,'-',i3.2,'-',i2.2,', ',i2.2,':',i2.2,')')
      write(16,142) 

142   format('jdstart.gmt   doy   pcon   flsam    flcal   pcell   rl   tcell   tdetect   tpump   ')  
      write(16,*) 'tdew.rh  tamb   c01     c10      c11   statcm'

!     ------------------------------------------------------------------
!     second:  initialize values to read the file
!     ------------------------------------------------------------------

      n1 = n1PF-2  ! number of fields in eddy 201 file is 10 (2 less than profile)
      ns=9
      debug= .false. 
      call read20n(21, 201, dtype, time0, v, 1, n1, is201, 1, ns,
     &             done, yrdays1, debug)
      if (done) goto 399
211   format(i4, 1x, f9.3, /, 8f8.3, /, f9.4, 1x, f9.4, 1x, 3i1, 1x, 
     &       2i1, 1x, 3i1, 1x, i1) ! eddy
!     same format as output in lsplit
      call gethhmmss(time0)
      write(99,202) time0, hh,mm,ss
202   format('  First gmtime in 201 file is ',f10.4,' (',i3.2,':',
     &        i2.2,':',f5.2,')')

      call startime(time0, tstart1, period)  ! get start-time tstart1
                                             ! default: next integral 'period'
      call gethhmmss(tstart1)
      write(99,203) tstart1, hh,mm,ss
203   format('  Start-time is   ',f10.4,' (',i3.2,':',i2.2,':',f5.2,')')

!     move to desired start time in file, acquire initial values:
      do  while(time0.lt.tstart1)
      call read20n(21,201,dtype,time0,v,1,n1,is201,1,ns,
     &             done,yrdays1,debug)
      if (done) goto 399
      end do
      tnext = tstart1
      numint=0

!     ------------------------------------------------------------------
!     third:  main loop for reading this file
!     ------------------------------------------------------------------

210   continue ! loops until this input file is exhausted;
               ! goes on until end-of-file causes transfer out;
               ! at end, goes to label 899
               ! data aggregating, by time:
               ! tstart= starting time for this interval, 
               ! tnext= start for next interval
               ! loops to read data until time of data-line >= tnext
!     ------------------------------------------------------------------
!     (1) initialize for this loop:
!     ------------------------------------------------------------------

      tstart = tnext                         ! start time for this interval
      tnext  = tnext + period/1440.d0        ! minutes to decimal julian days 1440=60*24 sec/day
      tmid   = tstart + (period/2)/1440.d0   ! mid-point of this period 1440=60*24 sec/day     
      call gethhmmss(tstart)
      write(99,220) numint+1, tstart,hh,mm,ss
220   format('  Aggregation interval: ',i4 / '          start-time: ',
     &        f10.4,' (',i2.2,':',i2.2,':',f5.2,')')
      call gethhmmss(tnext)
      write(99,221) tnext,hh,mm,ss
221   format('  Expected end-time :  ',f10.4,' (',i2.2,':',i2.2,':'
     &       ,f5.2,')')
      k=1
      do  i=1,n1               ! n1 = number of input variables
       xo(i) = 0.0d0           ! xo(i) = ith output
      end do
      do  i=0,3                ! cal status
       statcal(i)=0.0          ! 0 = frequency of 00 (normal sample)
      end do                   ! 1 = freq of 01 (normal zero/cal)
                               ! 2 = freq of 10 (zero-dilution)
                               ! 3 = freq of 11 (hi-flow zero)
      statcm=0.0               ! freq of chilled mirror status flag

!     ---------------------------------------------------------------------
!     (2) read in one period worth of data:
!     ---------------------------------------------------------------------
      if(time0.lt.tnext) then 
      do while (time0.lt.tnext) 

!     stores previous record's data: (i) voltages
      do  i=1,n1                                 ! converts voltages 1..n1 for point k
       xi = v(i)                                 ! ith input variable
       if(i.eq.5) then                           ! ready-light
        if(abs(v(i)) .ge. vminrl) then 
         xi=1.0d0                                ! > min voltage, ready light on
        else xi=0.0d0                            ! otherwise, licor off
        endif
       endif
       if (i.ge.7 .and. i.le.8) then             ! thermistor (tdetect or tpump)
        xi = tempf( -v(i) )                      ! thermistor+pull-up sensors
       endif
       if ((idate.ge.080101).and.(i.ge.9)) then  ! after 2017 we are using and humidity and temperature sensor
         cc(i) = 1                               ! not the chilling mirror
         offset(i) = 0
       endif
!     sum of each over k (after converting to correct units):
       xo(i) = xo(i) + xi*cc(i) + offset(i)
      end do
!          cal-status
!          ical = 2*is201(4) + is201(5)    ! 00 = 0 (data sample)
                                           ! 01 = 1 (normal calibration, including zero)
                                           ! 10 = 2 (dilution)
                                           ! 11 = 3 (hi-flow zero)
      icalin = 10*is201(4) + is201(5) 
      ical   = icalin
      if(icalin.gt.1) then                 ! icalin = 10 or 11
       ical  = 2
       if(icalin.eq.11) ical = 3
      endif
      statcal(ical) = statcal(ical) + 1    ! number of times in each calibration status.
!     chilled-mirror status
      statcm = statcm + icm
!     read next-line of data in this interval:
      call read20n(21,201,dtype,time0,v,1,n1,is201,1,ns,
     &             done,yrdays1,debug)
      if (done) goto 399      
       k = k+1                             ! number of samples in this interval.
      end do                    
       k = k-1
      endif
      call gethhmmss(time0)
      write(99,652) time0, hh,mm,ss, k
652   format('  Stopped reading 201 file at:  ',f10.4,' (',
     &        i2.2,':',i2.2,':',f5.2,'), cum. records= ',i5,'.')
      goto 901

!     ------------------------------------------------------------------
!     (3) aggregates this period's data:
!     ------------------------------------------------------------------
901   continue
      write(99,270) k, 201            
270   format('  Aggregating: ',i6,' values from ',i3, 'file...')
!     now take averages:
      do  i = 1, n1                        ! for each var of input line
       xo(i) = xo(i)/k
      end do
      do  i=1,3                            ! for zero/cal status bits
       statcal(i) = statcal(i)/k
      end do
      statcm   = statcm/k                  ! chilled mirror status
      miragdew = xo(9)                     ! mirsumdew/mirnum
      miragt   = xo(10)                    ! mirsumt/mirnum       

!     ------------------------------------------------------------------
!     (4) write housekeeping output:
!     ------------------------------------------------------------------
      tdoy = tstart - days00 - yrdays1
      !..................................................................................
      !nrc added set boundaries of what is possible given the write format constrains
      do i =1, 10
       if ((xo(i).lt.-999.999).or.(xo(i).gt.999.999))     xo(i) = 999.999	!f7.3
      enddo
      if ((tdoy.ge.1.).and.(tdoy.lt.367.).and.(tstart.ge.1.).and.(tstart.lt.100000.)) then
       write(16,161) tstart, tdoy, (xo(i), i=1,10), (statcal(i), i=1,3), 
     &               statcm

161    format(f10.4,1x,f8.4,1x,f8.2,3(1x,f8.2),f6.2,3(1x,f7.2),/,   
     &          1x,2(1x,f7.3),1x,3(f5.2),1x,f5.2 ) 
      endif
!     ------------------------------------------------------------------
!     (5) continue looping:
!     ------------------------------------------------------------------
!     Number of aggregation intervals processed:
      numint = numint+1
      write(99,*)'  End the loop.'
      write(99,*)

!     Continues reading the *.201 file:
      goto 210  

!     ------------------------------------------------------------------
!     (6) Close-up shop and exit when end-of-file reached:
!     ------------------------------------------------------------------
399   continue
      write(99,*)'  Reached end-of-file in file: ',
     &    fnamein(1:lfnin)//'.201' 
      print *,   '  Reached end-of-file in file: ',
     &    fnamein(1:lfnin)//'.201' 

400   continue

!     Converts minutes to decimal julian days (starting 2000, not 1970!)
      tnext = tnext - period/1440.d0 
      write(99,*)'  Finished processing "',
     &            fnamein(1:lfnin), '" input file.'
      print *,'  Finished processing "',
     &            fnamein(1:lfnin), '" input file.'
      write( *,410) tnext-tstart1, numint, period
      write(99,410) tnext-tstart1, numint, period
410   format('    (',f7.3,' days of data aggregated into ',i3,1x,f6.2,
     &          '-min intervals)')

!     Writes a blank line:
      write(*,*)
      write(99,*)

!     Closes housekeeping input file (*.201):
      close(21)

!     Closes housekeeping output file (*.hkp):
      close(16)

!     Closes error/log files:
      close(99)

990   continue
      print *, '  End of program LAGGRED.'
      print *

!     ------------------------------------------------------------------
!     End of program:
      end
!     ------------------------------------------------------------------


!     ------------------------------------------------------------------
!                                subroutines and functions
!     ------------------------------------------------------------------

      include 'lsubs.for'
      include 'lfluxsub.for'  ! (contains code from lmirror and lfilter)
