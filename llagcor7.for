!===================================================================================================
!                                      lba code:  llagcor7.for
!===================================================================================================
!     modified from boreas code lagcor.for by scott r. saleska in september 2004.
!     lagcor.for made by s.-m. fan in 27-jul-1993.
!
!     calculates correlation functions between vertical wind (w) and other variables (t, co2, h2o).
!
!     history:
!     12-nov-2004 - llagcor modified by ewg to retrieve bad sonic flags and licor ready light flags.
!===================================================================================================

      program llagcor
      include 'lparams.for'                     ! parameters for array-sizes + sys variables.
      parameter(nzdim=nezdim)
      include 'lco2vars7.for'                   ! variables for doing co2 concentration, including:
      include 'lfluxvar.for'                    ! filtering and housekeeping variables
      parameter (n1024=16, m=1024*n1024)        ! 34.13 minutes w. 8 hz data.
      real*4 cor1(801),cor2(801),cor3(801)
      real*4 cor4(801),cor5(801)
      real*8 wvar, tvar, sco2, sh2o 
      logical done
      real*4 sday(30), stime(7) 
      data stime/0.4375, 0.5208, 0.6042, 0.6875, 0.7708, 0.8542, 0.9375/  
      ! time(gmt) 1030    1230    1430    1630    1830    2030    2230
      ! time(lt)  0630    0830    1030    1230    1430    1630    1830

!     sets initial parameters, gets command-line etc:
      lname = 'llagcor'
      ftypes='1'                               ! eddy 1.
      call getfilterset(filtcode, filtset)     ! uses filtcode (default set in lparams) to establish
                                               ! default filtset.
      include 'linitial7.for'                  ! common initialization code .
                                               ! sets directory structure (root paths, dirsplit, 
                                               ! dirprocin, dirprocout), retrieves command line, 
                                               ! file names for processing [ to ifilen(1..nifiles) ]
                                               ! and displays initial program message.
      if (navg.eq.0) navg = freqhi*60.*period  ! number of points to aggregate.
      print *, '  using filter code '//filtcode

!===================================================================================================
!                                      main code starts here
!===================================================================================================
      infile=1
      dtype= ifilen(infile)(1:2)             ! dtype is data type string.
      enum=1
!     parameter files: 
      if(dtype.eq.'EB') then
      enum=2
      endif  
      print *,'  reading parameter file lsysdat.inp ... for enum ', enum
      call readsysd(enum)                   ! reads data from lsysdata.inp for eddy system enum.
      print *,'  pass 1'
      call readconv(enum)                   ! reads data from lconv.inp for conversion coefficients.
      print *,'  pass 2'

!     gets filename: 
!     raw data input (*.200, 201 data)
      fnamein = ifilen(infile)
      print *,'  getfillpath'
      call getfullpath( rootsplit,yymmdd, dirsplit, fnamein, lfnin) 

!     data output (*.cor and *.erl files) 
      fnameout = ifilen(infile)
      call getfullpath( rootprocout,yymmdd,dirprocout,fnameout,lfnout)
      print *, '  opening input files "',fnamein(1:lfnin)//'.200" '
      print *, '  sending output to   "',fnameout(1:lfnout)//'.cor"'
      print *, '                      "',fnameout(1:lfnout)//'.erl"'

!     opens files.
!     sets-up error/status file.
      open(99,file=fnameout(1:lfnout)//'.erl',status='unknown')  

!     error file for lags:
      write(99,55) platform, lname, iyear, imonth, iday,ihrs,imin,isec
55    format(a6,' ',a7,':  on date (yy-mm-dd) ', i4,'-',i2.2,'-',i2.2,
     &       ' (',i2.2,':',i2.2,':',i2.2,'), process:' )
      write(99,*) '  '//fnamein(1:lfnin)//'.200',' using filter code '
     &                //filtcode

!     200 data file:
      open(18,file=fnamein(1:lfnin)//'.200',status='old')                ! removed action=read
      read(18, *)                                                        ! header line 1   
      read(18, *)                                                        ! header line 2
      write(99,*) '  opened input file:  '//fnamein(1:lfnin)//'.200.'

!     201 data file:
      open(88,file=fnamein(1:lfnin)//'.201',status='old')                ! removed action=read
      read(88, *)                                                        ! header line 1   
      read(88, *)                                                        ! header line 2   
      read(88, *)                                                        ! header line 3
      read(88, *)                                                        ! header line 4
      write(99,*) '  opened input file:  '//fnamein(1:lfnin)//'.201.'

!     main output files:
      open(19, file=fnameout(1:lfnout)//'.cor',status='unknown')
      write(99,*) '  sending output to: '//fnameout(1:lfnout)//'.cor'
      write(99,*)

!     writes header to output file:
      write(19,140) enum,yymmdd,filtcode,platform,lname,iyear,imonth,
     &              iday,ihrs,imin,isec
140   format('2  18  LBA Eddy',i1,': ',a6, ' , filt=', a5, ' (',a6,' ',
     &       a8,' on: ',i4,
     &          '-',i2.2,'-',i2.2,', ',i2.2,':',i2.2,':',i2.2,')' )
      write(19,141)  
141   format('gmtstart    nseg    lag    wt    wc    wh    ct    ',
     &       'ht    n1    n2    n3    n4    nbads    nbadli    ',  
     &       'npt    tflx    co2flx    h2oflx')  

!     moves to start of fast data:
      call read200(julday, gmt, v200, is200, done)
!                     1     1     7     2
      if(done) goto 290                                 ! if at end of file.
      iscal = is200(2)                                  ! cal = 00 --> normal sampling.
      call gethhmmss(time0)
      write(99,202) time0, int(gmt), hh,mm,ss
202   format('  First gmt time in 200 file is ',f8.4,
     &       ' (',i5,' sec =',i3,':',i2,':',f5.2,')')

!-------------------------------------------------------------------------------------------
!     sets up the times for which to calculate the lag correlations
!-------------------------------------------------------------------------------------------
      npts = 0 
      do i=1, 4
       do j=1, 6 
        npts = npts + 1
        sday(npts) = julday + (i-1) + stime(j) 
       enddo
      enddo 
      istart = 1 
      do j=1, npts-1
       if (time0 .gt. sday(j) .and. time0 .lt. sday(j+1)) istart = j 
      enddo 
      print *, '  Starting point: ',istart,sday(istart)  
      
!     calls startime(time0, tstart1, period)
!     tstart1 determines desired start-time from function:
      nperiod=0
      do j = istart, npts
      tstart1 = sday(j)  
      call gethhmmss(tstart1)
      write(99,203) tstart1, hh,mm,ss
      write(*,203)  tstart1, hh,mm,ss
203   format('  Start-time is             ',f8.4,' (           ',i3.2,
     &       ':',i2.2,':',f5.2,')')
      write(99,*)'  moving to start time.'
      do while(time0.lt.tstart1) 
       call read200(julday, gmt, v200, is200, done)
       if(done) goto 291              ! end-of-file here, then we are done before we started.
      enddo

!     collects data:
      ncal=0                                            ! init num cal points.
20    continue
      n=0                                               ! init num data points.
      npostc=0                                          ! init num post-cal points.
      calbit    = .false.                               ! calibration bit indicator.
      nsonflags = 0                                     ! number of sonic flags.

!     checks sonic status bits : initialization
      do k=1,4
       diagarray(k)=0                                   ! status bits 15..12
       numbitarray(k)=0                                 ! counters for each status bit
      enddo
      numbadsonic=0                                     ! hard sonic errors
    
!     skip any cals, plus allow for delay after calibration
      do while(is200(2).ne.0)                     ! keeps reading while this is a calibration.
       calbit = .true.
       call read200(julday, gmt, v200, is200, done)
       if(done) goto 300                          ! end-of-file here when this is not a calibration.
      enddo
      if(calbit) then
       time0 = float(julday) + gmt/86400.d0       ! gets time at point right after cal/zer ended.
       call gethhmmss(time0)
       write(99,218) time0,hh,mm,ss
218    format('  Cal/Zero interruption finished at: ',f8.4,
     &       ' (',i2,':',i2,':',f5.2,')')
       write(99,*) '  (skipping extra 50 sec more to equilibrate flow)'
       do while(npostc.lt.ndelayed)               ! ndelayed = number of points for 50 seconds (set 
                                                  !            in readsysd subroutine).
        npostc = npostc+1
        if(is200(2).ne.0) goto 20                 ! if calibration encountered again, goes back 
                                                  ! (should not happen unless manual calibration is
                                                  ! forced, for instance).
        call read200(julday, gmt, v200, is200, done)
        if(done) goto 300                         ! end-of-file here
       enddo
       calbit=.false.
      endif

      tstart = float(julday) + gmt/86400.d0       ! gets time at beginning of this segment           <<<<<<<<<<<<
      tnext  = tstart+(m/8)/86400.d0
      call gethhmmss(tstart)
      write(99,220) nperiod+1, tstart, hh, mm, ss
220   format('  Aggregation segment: ',i4 / '  nominal start-time: ',
     &        f8.4,' (',i2,':',i2,':',f5.2,')')
      call gethhmmss(tnext)
      write(99,221) tnext,hh,mm,ss
221   format('  End-time :  ',f8.4,' (',i2.2,':',i2.2,':',f5.2,')')

      do while(n.lt.m)
      n=n+1                                 ! accumulates n points.
      time0  = float(julday) + gmt/86400.d0 ! gets time at beginning of this segment.
      co2(n) = -v200(1)                     ! reverse inversion of mux board on licor voltage. 
      h2o(n) = -v200(2)                     ! same (note: only does if using raw voltage, not if
                                            ! using cal-func).
      u(n)   = v200(3)
      v(n)   = v200(4)
      w(n)   = v200(5)
      t(n)   = v200(6)
      stat(n)= v200(7)
      cal(n) = is200(2)
      if(cal(n).ne.0) then
      call gethhmmss(time0)
      write(99,222) time0,hh,mm,ss
222   format('  Calibration bits interrupt segment at:     ',f8.4,
     &       ' (',i2.2,':',i2.2,':',f5.2,')')
      goto 20                ! if this is a calibration point, go back and flush until no more cals.
      endif

!     checks sonic status bits
!     determines status variable for each fast data point in the interval
      call statusbitarray(int(stat(n)))     ! unpacks sonic status word to diagarray, where
                                            ! diagarray(1...4) <--> status bits 15...12.
      do k=1,4                              ! checks each status bit of this point.
      if (diagarray(k).eq.1) then
      numbitarray(k)=numbitarray(k)+1       ! increment counts for each bit.
      endif
      enddo
                                       
!     determines bad sonic readings by -99.999 in temperature. the choice of using temperature is 
!     arbitrary; u, v or w would also work.
      if (t(n).eq.-99.999) then                                
       numbadsonic=numbadsonic+1                              
      endif
      call read200(julday, gmt, v200, is200, done)
      if (done) goto 300                     ! end-of-file here.
      enddo

!     licor status light
!     gets housekeeping inputs from *.201 file
!       a) matches file pointer in .201 time to time of the loop;
!       b) reads in 1 period worth of points (e.g. 1024 points,for 34.13 min at 1/2 hz)
!     variables declared in lfluxvars7.for and subroutines in lfluxsub.for.
!     as of 09-sept-2001: sets-up to allow access to all variables in *.201 file.

!     line up dates:
      call matchmirtm(tstart,88)            ! cursor positioned at dew point reading for tstart 
                                            ! after this line.
      n1 = n1pf-2                           ! number of fields in eddy 201 file is 10 -- 2 less than
                                            ! in profile.
      ns=9
!     initializes:
      kavg = freqlo * m/8.                  ! where m is # of points in 34.13 min @ 8hz 
      k=1                                   ! k = position index within interval 
      do i=1,n
       licorrlvec(i) = 0
      enddo
      licoroff = 0   

      call gethhmmss(tstart)
      write(99,650) tstart, hh,mm,ss
650   format('  starting to read 201 file at:  ',f8.4,
     &       ' (',i2.2,':',i2.2,':',f5.2,')')

!     reads in one period worth of data: 
!     reads first line:
      read(88,211,end=899)jd1,gmt1,v201ed,(is201(i),i=1,ns)  
211   format(i4 , 1x , f9.3 , / , 8f8.3 , / , f9.4 , 1x , f9.4 , 1x , 
     &       3i1 , 1x , 2i1 , 1x , 3i1 , 1x , i1) 
      time1=jd1+gmt1/dble(86400.)
!     then gets rest of interval:
      if(time1.lt.tnext) then
      do while (time1.lt.tnext) 

!     licor ready light filtering
!
!     variables:
!     vrl...............voltage licor ready light; the voltage of the ready light.
!     vminrl............the minimum voltage for a positive ready light.
!     badlicorrl........true/false licor ready light; array of ones and zeroes indicating whether 
!                       the ready light was on (0) or off (1).
!     licoroff..........number of times licor ready light is off (slow data).

      vrl = -1*v201ed(irl)*cc(irl) + offset(irl)                
      if (vrl.lt.vminrl) then
      badlicorrl(k) = 1                     ! ready light off, filter this point.
      licoroff = licoroff + 1
      else
      badlicorrl(k) = 0                     ! ready light on, do not filter.
      endif              

!     reads next-line of data in this interval:
      read(88,211,end=899)jd1,gmt1,v201ed,(is201(i),i=1,ns)
      time1=jd1+gmt1/dble(86400.)
      k=k+1                                 ! number of samples in this interval.        
      enddo  
      k=k-1
      endif
      call gethhmmss(time1)
      write(99,652) time1, hh,mm,ss, k
652   format('  Stopped reading 201 file at:  ',f8.4,
     &       ' (',i2.2,':',i2.2,':',f5.2,'), after ',i5,' records.')
      goto 901

!     inserts error code for premature end-of file here.
899   continue    
      write(99,*) '  Warning: Premature end-of-file for 201 file.' 
      write(*,*)  '  Warning: Premature end-of-file for 201 file.' 
901   continue

      if(k.le.kavg/2 .or. k.ge.kavg*2) then 
      write(99,655) k, kavg
655   format('  Warning: bad number of values read:  k= ',i5,
     &       ' (expected kavg= ',i5,').')
      write(99,*)
      endif

!     end of housekeeping code ---------------------------------------------------------------------
      
!     finds the number of points with licor ready light (rl) off
      numbadli= licoroff*16              ! licor ready light off               
      write(99,270) n, 200
270   format('  Calculating correlations with ',i6,' values from *.'
     &        ,i3, ' file...')
      write(99,*)

!     Detrends and calculates correlations:
      nperiod=nperiod+1
      call detrend(m,w,avg,wvar)  
      call detrend(m,t,avg,tvar)
      call detrend(m,co2,avg,sco2)
      call detrend(m,h2o,avg,sh2o)
      call calcor(t,  w,m,cor1,801)       
      call calcor(co2,w,m,cor2,801)       
      call calcor(h2o,w,m,cor3,801)       
      call calcor(co2,t,m,cor4,801)       
      call calcor(h2o,t,m,cor5,801)       

!     output results:
      do 80 i=1,801
      lag = i-401
      tw = sqrt(tvar*wvar)*cor1(i)  
      co2w = sqrt(sco2*wvar)*cor2(i)  
      h2ow = sqrt(sh2o*wvar)*cor3(i)  
      write(19,70) tstart,nperiod,lag,cor1(i),cor2(i),cor3(i),cor4(i),
     &             cor5(i),(numbitarray(k),k=1,4),numbadsonic,numbadli,
     &             m,tw,co2w,h2ow    
70    format(f8.4,1x,i4,1x,i4,5(1x,f8.4),7i7,1x,f12.4,2(1x,f10.4))
80    continue
      enddo                                               ! skips to the next selected time period   

!===== exit from one of the file-read loops ========================================================

!===== (1) exit due to premature end-of-file error

290   continue
      write(99,*)
      write(99,*)'  error:  end-of-file on first read for: ',
     &              fnamein(1:lfnin)//'.200'
      print *,   '  error:  end-of-file on first read for: ',
     &              fnamein(1:lfnin)//'.200'
      goto 400
      
291   continue
      write(99,*)
      write(99,*)'  error:  end-of-file before start time in file: ',
     &              fnamein(1:lfnin)//'.200' 
      print *,   '  error:  end-of-file before start time in file: ',
     &              fnamein(1:lfnin)//'.200' 
      goto 400

!===== (2) exit during later read (probably normal)

300   continue
      write(99,*)
      write(99,*)'  discard interval: reached end of: ',
     &              fnamein(1:lfnin)//'.200' 
      time0 = julday+gmt/86400.d0
      call gethhmmss(time0)
400   continue
      write(99,*)
      write(99,*)'  finished processing "',
     &              fnamein(1:lfnin), '" input file.'
      print *,   '  finished processing "',
     &              fnamein(1:lfnin), '" input file.'
      write(*,*)                                    ! blank line
      write(99,*)
      close(19)                                     ! close output file
      close(18)                                     ! close input files
      close(99)                                     ! close error/log file
990   continue
      print *, '  llagcor '//yymmdd//': program finished.'
      print *
      end                                          ! end of the program

!===================================================================================================
!                                    subroutines and functions
!===================================================================================================
      include 'lsubs.for'  
      include 'lfluxsub.for'
!===================================================================================================
!                               subroutine detrend(ndata,x,xmean,var)                              !
!===================================================================================================
!     input:  x (raw data), ndata (number of data)                                                 !
!     return: x (detrended data), xmean, var                                                       !
!---------------------------------------------------------------------------------------------------
      subroutine detrend (ndata,x,xmean,var)
      real*8 xmean, var
      real*8 x(ndata)
      sx=0.
      ss=float(ndata)
      si=ss*(ss+1.)/2.
      st2=0.
      b=0.
      var=0.
      sioss=si/ss
      do 14 i=1,ndata
      sx=sx+x(i)
      y=float(i)-sioss
      st2=st2+y*y
      b=b+y*x(i)
14    continue
      xmean=sx/ss
      b=b/st2
      a=(sx-si*b)/ss
      do 16 i=1,ndata
      x(i)=x(i)-a-b*float(i)
      var=var+x(i)*x(i)
16    continue
      var=var/ss
      return
      end
!===================================================================================================
!                           end of subroutine detrend(ndata,x,xmean,var)                           !
!===================================================================================================



!===================================================================================================
!                                         subroutine calcor
!===================================================================================================
      subroutine calcor(x,y,m,corr,mcor)
      parameter (nmax=65536)
      real*8 x(m),y(m)
      real*4 corr(mcor)
      real*4 a1(nmax),a2(nmax)
      m2=2*m
      mm=(mcor+1)/2
      kk=m-mm+1
      rm=1./float(m)
      call twofft(x,y,a1,a2,m)     
!     sets amplitude at low frequencies to zero (a high-pass filter):
      a1(1)=0.0 
      a1(2)=0.0 
      a1(3)=0.0 
      a1(4)=0.0 
      a1(5)=0.0 
      a1(6)=0.0 
      a1(7)=0.0 
      a1(8)=0.0 
      a1(9)=0.0 
      a1(10)=0.0 
      a1(11)=0.0 
      a1(12)=0.0 
      a1(m2-9)=0.0 
      a1(m2-8)=0.0 
      a1(m2-7)=0.0 
      a1(m2-6)=0.0 
      a1(m2-5)=0.0 
      a1(m2-4)=0.0 
      a1(m2-3)=0.0 
      a1(m2-2)=0.0 
      a1(m2-1)=0.0  
      a1(m2-0)=0.0  
      a2(1)=0.0 
      a2(2)=0.0 
      a2(3)=0.0 
      a2(4)=0.0 
      a2(5)=0.0 
      a2(6)=0.0 
      a2(7)=0.0 
      a2(8)=0.0 
      a2(9)=0.0 
      a2(10)=0.0 
      a2(11)=0.0 
      a2(12)=0.0 
      a2(m2-9)=0.0  
      a2(m2-8)=0.0  
      a2(m2-7)=0.0  
      a2(m2-6)=0.0  
      a2(m2-5)=0.0  
      a2(m2-4)=0.0  
      a2(m2-3)=0.0  
      a2(m2-2)=0.0  
      a2(m2-1)=0.0  
      a2(m2-0)=0.0  
!     end of high-pass filter.
      call realft(a1,m,-1) 
      call realft(a2,m,-1) 
      j=0     
      do 110 i=1,m2,2  
      j=j+1
      x(j)=rm*a1(i)  
      y(j)=rm*a2(i)  
110   continue
      call correl(x,x,m,a1)       
      v0=a1(1) 
      call correl(y,y,m,a1)       
      v1=sqrt(v0*a1(1))
      if (v1 .eq. 0.) goto 130
      call correl(x,y,m,a1)       
      corr(mm)=a1(1)/v1
      do 120 i=1,mm-1
      corr(i)=a1(i+kk)/v1
      corr(i+mm)=a1(i+1)/v1
120   continue
      return
130   print *, '  error: v1=0.0 in calcor subroutine.'
      return
      end
!===================================================================================================
!                                    end of subroutine calcor                                      !
!===================================================================================================



!===================================================================================================
!                                        subroutine correl                                         !
!===================================================================================================
      subroutine correl(data1,data2,n,ans)
      parameter(nmax=65536)
      real*8 data1(n), data2(n)
      complex fft(nmax),ans(n)
      call twofft(data1,data2,fft,ans,n)
      no2=float(n)/2.0
      do 11 i=1,n/2+1
      ans(i)=fft(i)*conjg(ans(i))/no2
11    continue
      ans(1)=cmplx(real(ans(1)),real(ans(n/2+1)))
      call realft(ans,n/2,-1)
      return
      end
!===================================================================================================
!                                    end of subroutine correl                                      !
!===================================================================================================



!===================================================================================================
!                                        subroutine realft                                         !
!===================================================================================================
      subroutine realft(data,n,isign)
      real*8 wr,wi,wpr,wpi,wtemp,theta
      dimension data(*)
      theta=6.28318530717959d0/2.0d0/dble(n)
      c1=0.5
      if (isign.eq.1) then
      c2=-0.5
      call four1(data,n,+1)
      else
      c2=0.5
      theta=-theta
      endif
      wpr=-2.0d0*dsin(0.5d0*theta)**2
      wpi=dsin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=2*n+3
      do 11 i=2,n/2+1
      i1=2*i-1
      i2=i1+1
      i3=n2p3-i2
      i4=i3+1
      wrs=sngl(wr)
      wis=sngl(wi)
      h1r=c1*(data(i1)+data(i3))
      h1i=c1*(data(i2)-data(i4))
      h2r=-c2*(data(i2)+data(i4))
      h2i=c2*(data(i1)-data(i3))
      data(i1)=h1r+wrs*h2r-wis*h2i
      data(i2)=h1i+wrs*h2i+wis*h2r
      data(i3)=h1r-wrs*h2r+wis*h2i
      data(i4)=-h1i+wrs*h2i+wis*h2r
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
      h1r=data(1)
      data(1)=h1r+data(2)
      data(2)=h1r-data(2)
      else
      h1r=data(1)
      data(1)=c1*(h1r+data(2))
      data(2)=c1*(h1r-data(2))
      call four1(data,n,-1)
      endif
      return
      end
!===================================================================================================
!                                     end of subroutine realft                                     !
!===================================================================================================



!===================================================================================================
!                                         subroutine twofft                                        !
!===================================================================================================
      subroutine twofft(data1,data2,fft1,fft2,n)
      real*8 data1(n), data2(n)
      complex fft1(n),fft2(n),h1,h2,c1,c2
      c1=cmplx(0.5,0.0)
      c2=cmplx(0.0,-0.5)
      do 11 j=1,n
      fft1(j)=cmplx(data1(j),data2(j))
11    continue
      call four1(fft1,n,1)
      fft2(1)=cmplx(aimag(fft1(1)),0.0)
      fft1(1)=cmplx(real(fft1(1)),0.0)
      n2=n+2
      do 12 j=2,n/2+1
      h1=c1*(fft1(j)+conjg(fft1(n2-j)))
      h2=c2*(fft1(j)-conjg(fft1(n2-j)))
      fft1(j)=h1
      fft1(n2-j)=conjg(h1)
      fft2(j)=h2
      fft2(n2-j)=conjg(h2)
12    continue
      return
      end
!===================================================================================================
!                                      end of subroutine twofft                                    !
!===================================================================================================



!===================================================================================================
!                                          subroutine four1                                        !
!===================================================================================================
      subroutine four1(data,nn,isign)
      real*8 wr,wi,wpr,wpi,wtemp,theta
      dimension data(*)
      n=2*nn
      j=1
      do 11 i=1,n,2
      if(j.gt.i)then
      tempr=data(j)
      tempi=data(j+1)
      data(j)=data(i)
      data(j+1)=data(i+1)
      data(i)=tempr
      data(i+1)=tempi
      endif
      m=n/2
1     if ((m.ge.2).and.(j.gt.m)) then
      j=j-m
      m=m/2
      go to 1
      endif
      j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
      istep=2*mmax
      theta=6.28318530717959d0/(isign*mmax)
      wpr=-2.d0*dsin(0.5d0*theta)**2
      wpi=dsin(theta)
      wr=1.d0
      wi=0.d0
      do 13 m=1,mmax,2
      do 12 i=m,n,istep
      j=i+mmax
      tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
      tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
      data(j)=data(i)-tempr
      data(j+1)=data(i+1)-tempi
      data(i)=data(i)+tempr
      data(i+1)=data(i+1)+tempi
12    continue
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
13    continue
      mmax=istep
      go to 2
      endif
      return
      end
!===================================================================================================
!                                     end of subroutine four1                                      !
!===================================================================================================


!===================================================================================================
!                                        subroutine read200                                        !
!===================================================================================================
      subroutine read200( julday, gmt, v200, is200, done)
      include 'lparams.for'
      parameter(nzdim=nezdim)                               ! number of zeros for eddy 'zeroing'.
      include 'lco2vars7.for'                               ! includes v200, is200.
      integer julday
      real*4 gmt
      logical done
10    continue
      read(18,200, end=90, err=30) julday, gmt, v200, is200
!     (?) isn't 20 in standard lba software set?  (? so why it is 18 here?)
200   format(i4,1x,f9.3,1x,2(f7.4,1x),4(f7.3,1x),f6.0,1x,i3,1x,i2)
      time0 = julday+gmt/86400.d0                          ! in common block 'sys8byt'.
      return
!     error:
30    call gethhmmss(time0)
      write(99,105) julday, gmt,  time0, hh,mm,ss
      write(*, 105) julday, gmt,  time0, hh,mm,ss
105   format('  warning:  error in reading e*.200:  jday: ', i3, ', ',
     &        f9.3,' s' ,/, '    [ decimal jday = ', f8.4, ' (',
     &        i2.2,':',i2.2,':',f5.2,') ]' )
      goto 10                                              ! continues until a good value is found.
90    done = .true.
      end
!===================================================================================================
!                                       end of subroutine read200                                  !
!===================================================================================================
