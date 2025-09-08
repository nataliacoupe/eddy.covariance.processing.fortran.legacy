!===================================================================================================
!                                       lba code:  laggr7.for                
!===================================================================================================
!     created by scott saleska  (started august 2000)
!
!     version 1.0 compiled on 2-april-01
!     (adapted from part of boreas code b3split by s.c. wofsy and s.-m. fan)
!
!     lba aggregate slow data
!     laggr 001212 -files 34 -min 30 (30 min is default)
!
!     program to aggregate low-freq data 0.5 hz data to 
!     specified intervals (10 min default)
!
!       input data files:
!     .\yymmdd\split\pfyymmdd.201 = scaled profile system data for id=201
!     .\yymmdd\split\pfyymmdd.202 = scaled profile system data for id=202
!     .\yymmdd\split\pfyymmdd.203 = scaled profile system data for id=203
!       and/or:
!     .\yymmdd\split\gdyymmdd.201 = scaled ground system date
!     .\yymmdd\split\gdyymmdd.202 
!     .\yymmdd\split\gdyymmdd.203 
!
!       ouput data files:
!     .\yymmdd\process\pfyymmdd.hkp = profile system housekeeping stuff       
!     .\yymmdd\process\pfyymmdd.clm = profile system climate data
!       and/or:
!     .\yymmdd\process\gdyymmdd.hkp = ground system housekeeping stuff       
!     .\yymmdd\process\gdyymmdd.clm = ground system climate data
!===================================================================================================
      program laggr
      include 'lparams.for'             ! parameters for array-sizes + sys vars
!     input variables:
      real*8 v(10+22+10)                ! sensor voltages
      real*8 store(42,maxlow), x(maxlow)! stored voltages  ! stmux(3,8,maxlow/8) 
!     output variables:
      integer iomap(60), sdmap(60)      ! index arrays to xo for mapping input index (to store, 
                                        ! stmux, or stbit ) to output index (to xo array)
      integer nout, nhkp
      integer f201, f202, f203          ! labels for input formats (201,202, 203)
      integer f16, f17, f18             ! labels for output format (hkp,.clm, sun)
      logical newyrflag, done, debug
      real*8 tdoy, yrdays1, yrdays2, yrdays3
      real*8 tempf, xmed, xout, xsum 
!     default initial parameters:
      lname = 'laggr6'
      ftypes= '34'
      include 'linitial7.for'           ! common initialization code.
                                        ! sets directory structure (root paths, dirsplit, dirprocin,
                                        ! dirprocout), retrieves command line, file names for 
                                        ! processing [ to ifilen(1..nifiles) ] and displays initial 
                                        ! program message.
      if (navg.eq.0) navg = freqlo*60.*period
      call readsysd(3)                  ! 3 = read profile system data (see 'lsubs7.for' for code)
      call readconv(3)                  ! conversion coeffs from lconv.inp, starting with profile

!===================================================================================================
!                                 main program loop starts here   
!===================================================================================================

      do infile=1, nifiles              ! loop over each input file
!     ----------------------------------------------------------------------------------------------
!     first >> setup: open the input and output files
!     ---------------------------------------------------------------------------------------------
      dtype= ifilen(infile)(1:2)        ! dtype is data type string
      include 'lfileagr.for'  
!     lfileagr opens and sets input and output files,depending on whether profile(PF) or ground(GD)
!     * open input files:          * opens output files:
!       file pointr 21 = *.201,      99 = *.era
!       file pointr 22 = *.202,      16 = *.hkp
!       file pointr 23 = *.203       17 = *.clm
!                                    18 = *.sun (if profile)
!     * sets file-dependent input and output parameters:  
!       n1,n2,n3 = number of sensor inputs on 201,202,203 lines
!       ns       = number of status bits on 201 line
!       nmux     = num of mux'd variables for this file
!       nvarin   = total num of non-mux'd vars
!       nskip    = number of status bits to skip in output
!       nout,nhkp = # of total output vars, number that are housekeeping
!       f16,f17,f18 = format identifiers for output to files 16,17,18 
!     * sets input-output map, iomap(1..60):
!       format of input-output map:
!       i=indx to iomap            'iomap(i)'  (indx to xo)
!       (order of input)           (order of output)
!       ----------------           ------------------
!       1 to nvarin;               output order of non-mux'd measurement
!       nvarin+1 to nvarin+nmux*8  where to output this mux'd measurement
!       nvarin+nmux*8+1  to        where to put this status output
!       nvarin+nmux*8+nos          (all status go to hkp file)
!       (nos=num of status)

!     ----------------------------------------------------------------------------------------------
!     second >> loop through input files to get & store data segments
!     ----------------------------------------------------------------------------------------------

!     ----------------------------------------------------------------------------------------------
!     a. initialize to read this file:
!     ----------------------------------------------------------------------------------------------
      write(99,*) '  Start reading input files...'
      debug=.false.
      yrdays1 = 0.         ! days in this year to add to time* var in event of new year
      yrdays2 = 0.
      yrdays3 = 0.
      call read20n(21,201,dtype,time1,v,1,n1,is201,1,ns,done, 
     &             yrdays1, debug)

!     input arguments :    fileptr (21), line-num (201,202,203), dtype (PF/GD)
!     output arguments:    time1, v(1..n1), is201(1..ns), done (flag), debug reset to false, and
!     also checks year-boundary crossings:  
!     yrdays1 set to num days in yr if year-boundary encountered

      if(done) goto 290
      call read20n(22,202,dtype,time2,v,n1+1,n1+n2,is201,1,1,done,
     &             yrdays2, debug)
      call read20n(23,203,dtype,time3,v,n1+n2+1,n1+n2+n3,is201,1,1,done,
     &             yrdays3, debug)

!     read first line of each file (error exit to 290 for file 201) initializes time1, time2, time3:
      call startime(time1, tstart, period)

!     unit 6 is terminal, 38 is string length:
      call writetimestamp(6,  '  First time in *.201 file is      ',38,
     &                    time1) 
      call writetimestamp(6,  '  Start-time in all files will be  ',38,
     &                    tstart)
      call writetimestamp(99, '  First time in *.201 file is      ',38,
     &                    time1)
      call writetimestamp(99, '  Start-time in all files will be  ',38,
     &                    tstart)

!     move to desired start time in each file, acquire initial values:
      debug=.true.
      if (time1.lt.tstart) then
      dowhile (time1.lt.tstart)
      call read20n(21, 201, dtype, time1, v, 1, n1, is201, 1, ns, done, 
     &             yrdays1,debug)
      if(done) goto 291
      enddo 
      endif   
      dowhile (time2.lt.tstart)
      call read20n(22, 202, dtype, time2, v, n1+1, n1+n2, is201, 1, 1, 
     &             done, yrdays2,debug)
      if(done) goto 292
      enddo 
      dowhile (time3.lt.tstart)
      call read20n(23, 203, dtype, time3, v, n1+n2+1, n1+n2+n3, is201, 
     &             1, 1, done, yrdays3, debug)
      if(done) goto 293
      enddo

!...................................................................................................
!     note:  
!     now v(i), is201(i) have initial values across all files
!...................................................................................................

      tnext = tstart
      numint = 0             ! number of aggregation intervals processed

!     ----------------------------------------------------------------------------------------------
!     b. file-loop to go through this file
!     ----------------------------------------------------------------------------------------------
200   continue               ! loop until this input file is exhausted; 
                             ! goes on until end-of-file causes transfer out;
                             ! at end, goto label 301,302,303 (depending on which file ends first).

!     ----------------------------------------------------------------------------------------------
!     (b.1) initialize data arrays:
!     ----------------------------------------------------------------------------------------------
      tnext = tnext + period/dble(1440.)  ! minutes to decimal julian days.
      yrdays0 = yrdays1                   ! saves yrdays in effect at beginning of interval.
      do i=0,3
      kcal(i)=0
      enddo
      ios=0
      ncaldif=0
      x1_4=0
      x5_8=0
      write(99,*) '  Getting data for next interval.'

!     ----------------------------------------------------------------------------------------------
!     (b.2) loop to read and store one aggregation interval:
!     (each file can be on different time-base; k1,k2,k3 will contain the number of samples in the 
!      interval for each file)
!     ----------------------------------------------------------------------------------------------
!     ---------
!     201 file:
!     ---------
      k1=0
      if(time1.lt.tnext) then               ! loop until all datapoints in this interval collected.
      dowhile (time1.lt.tnext) 
      k1=k1+1
      do i=1,n1                             ! stores voltages 1..n1 in for point k1.
      store(i,k1) = v(i)
      enddo
      do i=1,ns
      stbit(i,k1)=is201(i)                  ! stores status bits 1..ns.
      enddo
      call read20n(21, 201, dtype, time1, v, 1, n1, is201, 1, ns, done,
     &             yrdays1, debug)
      if(done) goto 301
      enddo                                 ! leaves time1, v(i), is201(i) with initial.
      endif                                 ! values for next aggregation interval.
!     ---------
!     202 file
!     ---------
      k2=0
      if (time2.lt.tnext) then              ! loop until all datapoints in this interval collected.
      dowhile (time2.lt.tnext)
      k2=k2+1
      do i=1,n2                             ! stores voltages in for point k2.
      store(n1+i,k2) = v(n1+i)
      enddo
      call read20n(22, 202, dtype, time2, v, n1+1, n1+n2, is201, 1, 1,
     &             done, yrdays2,debug)
      if(done) goto 302
      enddo 
      endif
!     ---------
!     203 file
!     ---------
      k3=0
      if (time3.lt.tnext) then              ! loop until all datapoints in this interval collected. 
      dowhile (time3.lt.tnext)
      k3=k3+1
      do i=1,n3                             ! stores voltages in for point k3
      store(n1+n2+i,k3) = v(n1+n2+i)
      enddo
      call read20n(23, 203, dtype, time3, v, n1+n2+1, n1+n2+n3, is201,
     &             1, 1, done, yrdays3, debug)
      if(done) goto 303
      enddo 
      endif

!...................................................................................................
!     note:  2(a-c) ends with:
!     store(1..n1 , 1..k1) = the k1 201-file voltages in this interval
!     store(n1+1..n1+n2, 1..k2) = the k2 202-file voltages
!     store(n1+n2..n1+n2+n3, 1..k3) = the k3 203-file voltages
!     stbit(1..ns, 1..k1) = the k1 bits for ns statuses in 201-file 
!   
!     v(1..n1+n2+n3) = the initial voltages from the next interval
!     is201(1..n2) = the initial statuses from the next interval
!     tnext = the time at the end of the next interval
!...................................................................................................

!     ----------------------------------------------------------------------------------------------
!    (b.3) find the summary statistic (e.g. mean) for this aggregation interval:
!     ----------------------------------------------------------------------------------------------
!     summarize all non-mux'd variables:
      call gethhmmss(time1)
      write(99,'(a21,i6,a13,f10.4,a2,i2.2,a1,i2.2,a1,f5.2,a2)') 
     &      ' aggregation interval: ', numint+1,'; next time: ',time1,
     &      ' (',hh,':',mm,':',ss,') '
      write(99,105) k1,k2,k3
105   format('  Aggregating ',i3,' values from *201; ',i3,' from *202; '
     &       ,i3,' from *203...')
      do i=1,nvarin               ! all non-mux'd vars (num vars input = nvarin=n1+n2+n3).
      if(i.le.n1) then            ! 201-file variables.
      navg=k1
      elseif (i.le.n1+n2) then    ! 202-file variables.
      navg=k2
      else                        ! 203-file variables.
      navg=k3
      endif
      do k=1,navg
      x(k)=store(i,k)             ! retrieves this variable to x.
      end do
      xo(iomap(i))=xsum(x,navg,i,dtype)
!...................................................................................................
!     1. maps ith series of store(i,.) to x,;
!     2. then, summarizes the series x
!     3. then, maps the summary to xo using iomap(i) 
!     (n.b:  if input variable 'i' is not slated for output, iomap(i)=0, and function 'xsum' 
!     converts voltage to correct units (based on index i and data type dtype (PF or GD), and 
!     returns appropriate summary statistics.
!     if navg=0, xsum returns -999.99)
!...................................................................................................
      xo(sdmap(i))= sxo   ! sxo common variable calculated by xsum
                          ! if sdmap(i) is not 0, then save standard deviation statistics for output 
      enddo  
!     ----------------------------------------------------------------------------------------------
!     (b.4) summarizes status bits from 201 line:
!     ----------------------------------------------------------------------------------------------
      write(99,*) '  Summarize status bits'
!     (b.4a): part that need special treatment because several status bits need to be combined to 
!             produce one output var (first nskip status bits).
!     (?)

!     (b.4a.i): system designator (bits 1-3 -- both 'PF' and 'GD')
      ios = 1                  ! io offset for status output (increment by 1 for each status output)
      i = nvarin+nmux*8 + ios  ! offset relative to other variables 
					 ! as of 9/2001:  i=nvarin+1=43, iomap(i)=9
      xo(iomap(i)) = ismux(stbit(3,1),stbit(2,1),stbit(1,1))-1
!     system designator (1st 3 status bits, k=1 index doesn't change).
!     (n.b.:  profile=011 (sys 3); ground=100 (sys 4). use function ismux here but note: ismux makes 
!      000 = 1)

!     (b.4a.ii): calibration status (bits 4-5) or cal-control (9-10) (if 'PF')
      if(dtype.eq.'PF') then
      do k=1,k1
      ical = 2*stbit(9,k) + stbit(10,k)  	  ! 00 = 0 (data sample)
                                              ! 01 = 1 (normal cal, incl zero)
                                              ! 10 = 2 (dilution)
                                              ! 11 = 3 (hi-flow zero)
      kcal(ical) = kcal(ical) + 1  		      ! num of times in each cal status
						                      !  is cal control (bits 9-10) diff from cal status?
      icalcon = 2*stbit(9,k) + stbit(10,k)        ! calibration-control
      if(ical.ne.icalcon) ncaldif = ncaldif + 1   ! counts number of times the calibration status is
                                                  ! different from cal-control.
      end do
!     output:  cal-statuses (outputs 10-12)
      do ical = 1,3                               ! summarize different cal statuses (1,2,3).
      ios = ios + 1                               ! increment io offset for status output.
      i = nvarin+nmux*8 + ios                     ! ios= iomap offset relative to previous vars.
      if(k1.eq.0) then 
      xo(iomap(i)) = 0.
      else
      xo(iomap(i)) = kcal(ical)/float(k1)         ! fraction of time spent in this cal state.
      endif
      enddo
!     output:  is cal-control different from cal? (profile only)   (outputs 13)
      ios = ios + 1                               ! increment io offset for status output
      i = nvarin+nmux*8 + ios                     ! offset relative to previous variables.
      if(k1.eq.0) then 
      xo(iomap(i)) = 0.
      else
      xo(iomap(i)) = ncaldif/float(k1)            ! fraction of time cal-control differs from cal.
      endif
      endif
!     (b.4b):  status bits that are an output onto themselves (PF and GD)
      write(99,*) '  Last part of status bits.'
      do is=(nskip+1),ns                       ! nskip+1 = 11 = beginning of individual status bits
      do k=1, k1                               ! ns = nspf = 24 for profile
      x(k) = float(stbit(is, k))               ! retrieves this interval's bits (as real).
      enddo
      xtemp = xsum(x,k1,nvarin+nmux+(is-nskip), dtype)  ! summarize interval (default: fraction of
                                                        ! of interval bit is set).

      if (k1.eq.0) xtemp=0.                             ! override result if no data in interval.
!     status index ranges from 1 to (ns-nskip) beyond number of voltage variables.
!     gets temporary x-summary: fraction of time bit 'is' is set.
       if ((dtype.eq.'PF').and.(is.ge.16)) then  ! if 'is'>=16 ---> status indicates level
!     special treatment here: aggregates several status bits for common output. each bit's sum, 
!     xtemp= 0-9 = number of times this level was sampled, rounded to nearest whole number.
!     aggregates levels 1-4 into 4-digit real; likewise for levels 5-8.
      if (is.le.19) then 
      x1_4 = x1_4 + int(xtemp+.5) * 10**(19-is)          ! each level is a power of 10.
      else if(is.le.23) then
      x5_8 = x5_8 + int(xtemp+.5) * 10**(23-is)
      endif
      else                              ! default: generates output corresponding to this status.
      ios = ios + 1
      i = nvarin+nmux*8 + ios
      xo(iomap(i)) = xtemp              ! fraction of interval that this bit is set.
      endif
      enddo
      if (dtype.eq.'PF') then           ! at end of status bit gives a summary for PF.
      ios = ios+1    
      xo(iomap(nvarin+nmux*8+ios)) = x1_4   ! levels 1-4.
      ios = ios+1
      xo(iomap(nvarin+nmux*8+ios)) = x5_8   ! levels 5-8.
      ios = ios+1
      xo(iomap(nvarin+nmux*8+ios)) = xtemp  ! column sum.
      endif
!     ----------------------------------------------------------------------------------------------
!     (b.5) write this interval's data to output files:
!     ----------------------------------------------------------------------------------------------
      write(99,*) '  Write to file'
      tbeg = tnext - (period)/1440.   ! initial time for this period
      tdoy = tbeg - days00 - yrdays0  ! doy of initial time (yrdays0 was read at start of interval).
      if ((tdoy.ge.1.).and.(tdoy.lt.367.)) then		      ! nrc 161212 to avoid non-realistic dates
       write(16,f16) tbeg, tdoy, (xo(i), i=1,nhkp)            ! housekeeping
!.....nrc to avoid improper format..................................................................
       do i=(nhkp+1),nout
        if ((xo(i).lt.-999.99).or.(xo(i).gt.999.99)) then
         xo(i)=-999.99
        endif
       enddo
!....................................................................................................
       write(17,f17) tbeg, tdoy, (xo(i), i=nhkp+1,nout)       ! climate data
       if (dtype .eq. 'PF') then                              ! writes sunshine sensor data  
        write(18,f18) tbeg, tdoy, (xo(i), i=nout+1,nout+4)    ! 2 means, then 2 sd's
       endif  
      end if
!     profile system output formats: ---------------------------------------------------------------
161   format(f11.5,1x,f9.5, f8.2,3(1x,f8.2),1x,f7.2,3(1x,f7.2),/,f4.0,
     &       1x,3(f5.2),1x,f5.2,2x,2(1x,f5.2),3(f5.2),1x,2(f6.0),f3.0)  ! *.hkp
171   format(f11.5,1x,f9.5, 7(1x,f8.3),/,8(1x,f8.3),/,8(1x,f8.3),/,
     &       7(1x,f8.3) )                                               ! *.clm
181   format(f11.5,1x,f9.5, 1x,4(1x,f10.3) )                            ! *.sun

!     ground system output formats: ----------------------------------------------------------------
162   format(f11.5,1x,f9.5, 1x,f7.2,4(1x,f8.2),/,7(1x,f8.2),/,
     &       1x,f3.0, 5(f5.2),1x,3(1x,f5.2),1x,3(1x,f4.2),1x,f5.2)      ! *.hkp
172   format(f11.5,1x,f9.5, 1x,f8.3,6(1x,f8.3),/,8(1x,f8.3),/,
     &       8(1x,f8.3),/,2(1x,f8.3))                                   ! *.clm
      numint= numint+1                ! number of aggregation intervals processed
      goto 200  

!===================================================================================================
!                                 main program loop ends here   
!                                (not the end of the program!)
!===================================================================================================      
!     file loop: 
!     reads this set of input files and repeats until end-of-file is reached with read statement.
!     ----------------------------------------------------------------------------------------------
!     third >> close up all lfiles used here and exit from one of the file-read loops.
!     ----------------------------------------------------------------------------------------------

!     (1) exit due to premature end-of-file error:
290   continue
      write(99,*)'  Error:  end-of-file on first read for: ',
     &    fnamein(1:lfnin)//'.201'
      print *, '  Error:  end-of-file on first read for: ',
     &    fnamein(1:lfnin)//'.201'
      goto 400
      
291   continue
      write(99,*)'  Error:  end-of-file before start time in file: ',
     &    fnamein(1:lfnin)//'.201' 
      print *, '  Error:  end-of-file before start time in file: ',
     &    fnamein(1:lfnin)//'.201' 
      goto 400
292   continue
      write(99,*)'  Error:  end-of-file before start time in file: ',
     &    fnamein(1:lfnin)//'.202' 
      print *, '  Error:  end-of-file before start time in file: ',
     &    fnamein(1:lfnin)//'.202' 
      goto 400
293     continue
      write(99,*)'  Error:  end-of-file before start time in file: ',
     &    fnamein(1:lfnin)//'.203' 
      print *, '  Error:  end-of-file before start time in file: ',
     &    fnamein(1:lfnin)//'.203' 
      goto 400

!     (2) exit during later read (probably normal):
301   continue
      write(99,*)'  Reached end-of-file in file: ',
     &    fnamein(1:lfnin)//'.201' 
      print *, '  Reached end-of-file in file: ',
     &    fnamein(1:lfnin)//'.201' 
      goto 400
302     continue
      write(99,*)'  Reached end-of-file in file: ',
     &    fnamein(1:lfnin)//'.202' 
      print *, '  Reached end-of-file in file: ',
     &    fnamein(1:lfnin)//'.202' 
      goto 400
303     continue
      write(99,*)'  Reached end-of-file in file: ',
     &    fnamein(1:lfnin)//'.203' 
      print *, '  Reached end-of-file in file: ',
     &    fnamein(1:lfnin)//'.203' 
400   continue
      tnext = tnext - period/1440.                           ! minutes to decimal julian days 1440=60*24 sec/day
      write(99,*)'  Finished processing "',fnamein(1:lfnin), 
     &           '" input file(s).'
      print *,   '  Finished processing "',fnamein(1:lfnin), 
     &           '" input file(s).'
      write( *,410) tnext-tstart, numint, period
      write(99,410) tnext-tstart, numint, period
410   format('  (',f7.3,' days of data aggregated into ',i3,1x,f6.2,
     &          '-min intervals)')
      close(16)                                             ! closes output files (*.hkp)
      close(17)                                             ! (*.clm)
      if(dtype.eq.'PF') then 
      close(18)                                             ! closes profile output file (*.sun).
      endif
      close(21)                                             ! closes input files.
      close(22)      
      close(23)
      close(99)      
      enddo           ! does infile=1, nifiles (loop over each set of input files.
990   continue
      print *
      print *, '  End of program LAGGR.'
      print *
      end  

!===================================================================================================
!                                         end of program
!===================================================================================================      


!===================================================================================================
!                                   subroutines and functions:      
!===================================================================================================

      include 'lsubs.for' 

!===================================================================================================
!                                           function xsum                                          !
!===================================================================================================
!     gets summary statistics for array of 'navg' input voltages in 'x';
!     does conversion based on indx, dtype ('PF' or 'GD') and based on knowing order of input file 
!     variables. 
!     if navg=0, then returns -999.99.
!     note:  indx ranges are
!     1 to nvarin:  non-mux inputs  
!                  (num on 201,202,203 input lines)
!     nvarin+1 to nvarin+nmux:             mux'd variables 1 to nmux 
!     nvarin+nmux+1 to nvarin+nmux+nstat:  status variables 1 to nstat

      function xsum( x, navg, indx, dtype )
      include 'lparams.for'                     ! parameters for array-sizes + sys vars
      real*8 tempf, xmed, xout, xsum, x(navg) 
      real*8 x1, xn, sx1                        ! needs double precision to calculate variance.
      real*8 rad

!     number of (non-time, non-status bit) variables on input file lines:
!     (need constant expressions for these in select case statement in function 'xsum')

      if(navg.eq.0) then
      write(99,*) '  Note: zero values to summarize in this interval'
      xout=-999.99
      sxo =-999.99
      else
      x1 =0d0
      sx1=0d0
!     default:  just calculates mean and standard deviation after applying linear conversion
!     (if these are status bits, then cc=1, offset=0 and x(i)= 0 or 1, so returns fraction of 
!     aggregation interval that status bit=1))
      do i=1,navg                               ! accumulates total of converted signals.
      xn = x(i) 
      x1 = x1 + xn
      sx1 = sx1 + xn*xn
      end do
      x1 = x1/dble(navg)                        ! double precision mean.
      xout = x1 *cc(indx) + offset(indx)        ! converts to physical units.
      sxo = sdev(sx1,x1,navg) *cc(indx)  
10    format('  xo(',i2 ,' ) = ',f9.3,', sdev = ',f11.8,', navg = ',i3,
     &       ' (sign: ',f3.0,')')

!     specific exceptions --------------------------------------------------------------------------

!     --------------
!     profile system
!     --------------
      if(dtype.eq.'PF') then 
      select case(indx)                         ! indx is variable number in input file.
!     exception 1:  ready-light
!     -------------------------
      case(5)                                   ! 5 = ready light
      xout=0
      do i=1,navg                               ! gets mean of 1's and 0's
      if(abs(x(i)) .ge. vminrl) xout=xout+1     ! if > min voltage, ready light "on".
      end do
      xout = xout/float(navg)                   ! fraction of time ready light was "on".

!     exception 2:  thermistors
!     -------------------------
      case(7:8,(n1PF+12):(n1PF+20))             ! 201: 7 (tdetect), 8 (tpump); 
                                                ! 202: 12 (tair_cal), 13-20 (tair1-tair8)
      xout = tempf( -xmed(x,navg) )             ! uses median for thermistor+pull-up sensors;
                                                ! calculates with polynomial in tempf.
                                                ! note: rxmeverse signal
      xout = xout*cc(indx) + offset(indx)       ! provides linear calibration correction (if any).
 
!     exception 3: adjustment to wind-direction
!     -----------------------------------------
      case(n1PF+1)                              ! wind direction (line 202, first position)
      x1 =0d0                                   ! (re-calculates avoiding windir changes > 300 deg)
      sx1=0d0      
      do i=1,navg  
      wd = x(i) *cc(indx)+offset(indx)
      if (i.gt.1) then 
      delwd = wd - oldwd
      if(delwd.gt.+300) wd = wd-360.
      if(delwd.lt.-300) wd = wd+360.
      endif
      oldwd = wd
      x1 = x1 + wd
      sx1 = sx1 + wd*wd
      end do
      x1 = x1/dble(navg)
      xout = x1
      sxo = sdev(sx1,x1,navg)                   ! function to get sdev from sum-of-square sx1.

!     exception 4:  net radiometer calculation
!     ----------------------------------------
      case((n1PF+3):(n1PF+4))                   ! 202: 3 (netrad) & 4 (netrad_cal)
      x1=0d0
      sx1=0d0
      do i=1,navg
!     commented the following lines as when signal es negative it should be negative radiation....
!      if(-x(i).ge.0.) then       ! signal = -x
!      rad = -x(i) *cc(indx)      ! signal positive
!      else 
!      rad = -x(i) *offset(indx)  ! signal negative
!      endif
      rad = -x(i) *cc(indx)      ! signal positive
!     conversion factor is different if signal voltage is positive or negative.
!...................................................................................................
!     to add:  rad = rad * wind-adjust factor          (?)                           
!...................................................................................................
      x1 = x1 + rad
      sx1 = sx1 + rad*rad
      end do
      x1 = x1/dble(navg)
      xout = x1 
      sxo = sdev(sx1,x1,navg)

!     exception 5:  rain should be summed, not averaged
!     -------------------------------------------------
      case(n12PF+9 )                            ! 203:  9 = rain tips
      xout = xout*float(navg)                   ! retrieves sum from average.

!     exception 6:  a few of the status bits (when indx > n13PF)
!     ----------------------------------------------------------
!     note:  1-5 (solenoid level) does by default:  fraction of time bit=1
      case( (n13PF+6):(n13PF+14) )            ! 6-14:  levels 1-8, column
      xout = xout*float(navg)/nflsteplo       ! restores sum from average, then divide by nflstep
                                              ! to get number of segments at a single level.
      end select
!     -------------
!     ground system
!     -------------
      elseif (dtype.eq.'GD') then               ! dtype = 'GD'

!     only exception for GD is thermistors:
!     -------------------------------------
      if(indx.eq.5 .or. (indx.ge.7 .and. indx.le.23 )) then  ! changed range 5-23 from 8-23
                                                             ! exclude 6 (ttrap_co)
      xout = tempf( - xmed(x,navg)  )                        ! (5-7 sensors in shack)
      xout = xout*cc(indx) + offset(indx)                    ! linear calibration corr. (if any)
      endif

!     all the rest are simple averages, or status-bit counts.   
      endif                                           ! end for "if(dtype.eq.'PF') else ('GD')"
      endif                                           ! enf for "if(navg.eq.0)  else ..."
      xsum = xout
      return
      end                   
!===================================================================================================
!                                      end of function xsum                                        !
!===================================================================================================


!===================================================================================================
!                                          function sdev                                           !
!===================================================================================================
!     calculates standard-deviation from sum-of-squares of series and from the mean of the series.
!     checks for:
!     -  series of length 1, 
!     -  [sum-of-squares - n mean^2] <0, that can happen rarely due to round-off error, when there 
!        is no variation in series.
!---------------------------------------------------------------------------------------------------

      function sdev(ss,xbar,n)
      real*8 ss, ssdev                 ! sum-of-squares = sum( xi^2) over all i (double precision).
      real*8 xbar                      ! double-precision mean.
      real*8 sxo,sdev
      integer n
      ssdev = (ss - n*xbar*xbar)       ! unnormalized sum-of-squared deviations
      if(n.eq.1 .or. ssdev.lt.0.) then ! rarely, sx1<0 when it should be 0, due to round-off error
                                       ! even when using double precision.
      sxo=0.
      elseif(n.gt.1) then
      sxo = sqrt( ssdev/dble(n-1) ) 
      else 
      sxo = -999.99
      endif
      sdev = sxo
      return 
      end
!===================================================================================================
!                                      end of function sdev                                        !
!===================================================================================================



!===================================================================================================
!                                         function ismux                                           !
!===================================================================================================
      function ismux(i1,i2,i3)

!     ismux =  1    2    3    4    5    6    7    8           
!        i1 =  0    1    0    1    0    1    0    1
!        i2 =  0    0    1    1    0    0    1    1
!        i3 =  0    0    0    0    1    1    1    1

      ismux=0
      if (i3 .eq. 0) then
      if (i1 .eq. 0 .and. i2 .eq. 0) ismux=1
      if (i1 .eq. 1 .and. i2 .eq. 0) ismux=2
      if (i1 .eq. 0 .and. i2 .eq. 1) ismux=3
      if (i1 .eq. 1 .and. i2 .eq. 1) ismux=4
      else
      if (i1 .eq. 0 .and. i2 .eq. 0) ismux=5
      if (i1 .eq. 1 .and. i2 .eq. 0) ismux=6
      if (i1 .eq. 0 .and. i2 .eq. 1) ismux=7
      if (i1 .eq. 1 .and. i2 .eq. 1) ismux=8
      end if
      return
      end
!===================================================================================================
!                                     end of function ismux                                        !
!===================================================================================================

