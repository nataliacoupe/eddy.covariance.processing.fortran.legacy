!===================================================================================================
!                                        lsplit6.eddy
!===================================================================================================
!     (adapted from b1split boreas code written by song-miao fan)
!     lba split raw data
!     invoke:  'lsplit -date yymmdd -files 1234'
!     takes cr10 binary data, converts them into ascii form,and then separates into readings
!     from individual sensors, applies scale factors to obtain original sensor voltages from
!     cr10-logged voltage-divided mv.
! 
!     input files:
!     .\yymmdd\yymmddxx.dat
!     output files:
!     .\yymmdd\split\xxyymmdd.200, 201, 202, 203
!     error/log file to:
!     .\yymmdd\process\xxyymmdd.err
!
!     6 possible kinds of cr10 input data from files named:  yymmddxx.dat
!     where xx depends on site_id of system from which data comes:
!          xx:  EA       EB       PF         GD        sn        nt
!     site_id:  1=eddy1, 2=eddy2, 3=profile, 4=ground, 5=sonic2, 6=natal
!
!     input cr10 data fields (by site, and data-line id #):
!
!     site = 1 or 2 (EA or EB):
!     id=100 (initialization loop header) - 6 fields:
!     100, yr, day, hhmm, ss, site_id
!     id=200 (8 hz) --  1+7 = 8 fields:
!        200, co2, h2o, u, v, w, c, status {, tfw}
!
!     id=201 (0.5 hz) - 1+3+2+10 = 16 fields: 
!        201, day,      hhmm,       ss,        status1,    status2, 
!             pcon,     flsam,      flcal,     pcell,      rl,        tcell,
!             tdetect,  tpump,      tamb       rh.tdew
!
!     id=202 (0.5 hz) - 1+3+6 = 10 fields:
!        202, day,      hhmm,       ss,         
!             spare,    spare,      spare,     spare
!
!     id=203 (0.5 hz) - 1+3+10 = 14 fields:
!        203, day,      hhmm,       ss,         
!             cmstatus, spare,      spare,     spare,      spare,     spare,
!             spare,    spare,      spare,     spare
!
!     201 line was set to this (cmstatus at end, now in 203):
!     id=201 (0.5 hz) - 1+3+2+11 = 17 fields: 
!        201, day,      hhmm,       ss,        status1,    status2, 
!             pcon,     flsam,      flcal,     pcell,      rl,        tcell,
!             tdetect,  tpump,      tamb,      rh.tdew,    cmstatus
!
!     notes:
!     201: status1 is a spare
!     202: 4 spares are mux inputs where, e.g., sonic analog outputs could go
!     203: spares:  7 are digital status signals, last 2 are pulse inputsc id=200 (4 hz) -- 
!                   1+7 fields:
!     in output, put tamb, tdew.rh in 201 line, put cmstatus in status bits
!
!     site = 3 (PF = profile):
!     id=100 (initialization loop header) - 1+5 fields:
!     100, yr, day, hhmm, ss, site_id
!     id=201 (0.5 hz) - 1+5+10=16 fields
!        201, day,       hhmm,      ss,       status1,     status2,
!             pcon,      flsam,     flcal,    pcell,       rl,         tcell,
!             tcinlet,   tpump,     co2,      h2o,         spare,      spare
!     id=202 (0.5 hz) - 1+3+20 = 24 fields
!        202, day,       hhmm,      ss,
!             wdi,       spare,       
!             netrad,    netrad_cal,
!             par1_up,   par1_down, par2_up,  par_cal,     spare,      spare,
!             spare,     spare,     tair_1,   tair_2,      tair_3,     tair_4,
!             tair_5,    tair_6,    tair_7,   tair_8
!     id=203 (0.5 hz) - 1+3+10=14 fields
!        203, day,       hhmm,      ss,
!             ws1,       ws2,       ws3,      ws4,         spare,      spare,
!             spare,     spare,     rain,     spare
!
!     site = 4 (GD = ground):
!     id=100 (initialization loop header) - 6 fields:
!        100, yr, day, hhmm, ss, site_id
!     id=201 (0.5 hz) - 1+5+16=22 fields
!        201, day,       hhmm,       ss,         status1,     status2,    
!             p_cohs,    p_col,      p_ss,       p_hs,        p_ms,     p_ls,
!             p_zag,     tpump_co,   co,         tcell_co,    pcon_co,  pamb,
!             tcatalyst, par(1-8),   tsoil(1-8), tsoil(9-16)
!    
!     site_id = 5 (sn = sonic):
!     id=100 (initialization loop header) - 6 fields:
!        100, yr, day, hhmm, ss, site_id
!     id=200 (4 hz) --  1+5 fields:
!        200,  u, v, w, t, status
!
!     all cr10 units:
!     id=300 (low power output) - 4 fields: 
!        300, day, hhmm, ss, ??
!---------------------------------------------------------------------------------------------------

      program lsplit
      include 'lparams.for'                 ! to be compatible with 'linitial7.for'
      character*1 cr10(4096), ctmp
      character*96 basename, fnamerr        ! name of output file
      integer nifiles, infile                 
      integer id, hhmm, status1, status2, icsat 
      integer*2 icr10(2048), itmp
      integer isol1(12), isol2(12)          ! 12-bit status words 
      integer is200(5)                      ! 24 bits is max needed for is201 
      integer f200, f201, f202, f203        ! format line labels 
      logical lskip, output                 ! whether to skip, ascii format, output data
      real startnum                         ! start.numid: start= number of label 100 to start with
      integer nid, nidstart, nidend         ! number of labels, num to start, num to end with
      integer jrec                          ! index to block number (block = 2048 2-byte words)
      integer norec                         ! number of records (= # discrete cr10 lines)

!...................................................................................................
!     note:  3 kinds of errors:  
!     1. bad-data error (this word uninterpretable according to cr10 binary encoding scheme)
!     2. invalid id# error:  expected line id# (100,200,201,etc) but found something else
!     3. time-mismatch error:  this time comes before previous time 
!...................................................................................................

      integer noerr   ! number of bad-data errors within this block
      real*4 toterr   ! total number of bad-data errors in this file
      integer nblkerr ! number of blocks with >=1024 bad-data errors
      integer niderr, totiderr  ! number of invalid id errors
      integer ntierr, tottierr ! number of time mismatch errors
      integer eof     ! check for end-of-file for direct access file 
      
!      parameter (nd200=8, nd201=13, nd202=24, nd203=10)
      parameter (nd200=8, nd201=15, nd202=24, nd203=10)
      real*4 acr10(2058), rcr10(2048)      ! input ascii array for record block
      real*4 a200(nd200), a201(nd201), a202(nd202), a203(nd203)
                                           ! input ascii arrays (not incl. ID/date/time/status data)
      real*4 v200(nd200), v201(nd201), v202(nd202), v203(nd203)
                                           ! output voltage arrays (not incl. ID/date/time/status data) 
      real*4 sf200(nd200),sf201(nd201),sf202(nd202),sf203(nd203)
                                           ! scale-factor arrays  ( v20n = a20n / sf20n ) 
      real*4 sf200SN(nd200)  

      equivalence (cr10(1),icr10(1))       ! arrays co-located in memory

      data sf200/2*500., 6*1./ 
      data sf200SN/8*1./  !no gain for the for csat3
      real*8 ttsn
      real*4 sf201E(10)                    ! 201=0.5 Hz (eddy & profile combined)   **modified from 12 to 10 by nrc 180302**
                                           ! before     Pcon Flsam Flcal Pcell RL Tcell Tdetect Tpump CO2  H2O   spare
                                           ! 2016 after pcon flsam flcal pcell RL tcell tdetect tpump tamb rh.tdew
      data sf201E/6*500., 2*1000., 2*1./   ! **modified by nrc 180302**
      real*4 sf203E(1)      
      data sf203E/1000./                   ! 203 data
      real*4 sf201PF(12)                   !  **new for PF nrc 180302**
      data sf201PF/6*500., 2*1000., 2*500., 2*1000./ ! profile includes two extra compare to eddy       **new for PF nrc 180302**
      real*4 sf202PF(20)                   ! 202 profile line: wdi spare netrad(x2) par(x4) spare(x3) tsoil(x9)
      data sf202PF/2*500., 2*1000., 4*1000., 3*1000., 9*1000./
      real*4 sf203PF(10)                   ! 203 line: ws(x4), spare(x4), rain, spare 
      data sf203PF/8*1000., 2*1./ 
      real*4 sf201GD(7)                    ! 201 line ground: CO pcon fl stat tcell ttrap tpump
      data sf201GD/4*500., 3*1000./        ! fl, stat still tbd
      real*4 sf202GD(24)                   ! 202 line ground: tsoil(x16)  par(x8)
      data sf202GD/16*1000.,8*1000./
      real*4 sf203GD(8)                    ! 203 data  
      data sf203GD/ 7*500., 250./
      real*4 sf201SN(3)       
      data sf201SN/3*1./   
      data is200/5*0/   
      data is201/24*0/   

!nrc added to sort EB200 file...........................................................
!      integer idcount
      integer   indEB, indEB1, count200EBmax
      parameter (count200EBmax=20*24*60*60*8)			!assume 20 days data
      real    var2(count200EBmax)
      real    var3(count200EBmax), var4(count200EBmax), var5(count200EBmax)
      real    var6(count200EBmax), var7(count200EBmax), var8(count200EBmax)
      real    var13(count200EBmax),var14(count200EBmax),var15(count200EBmax), var16(count200EBmax)
      integer ina1(count200EBmax), ina12(count200EBmax)
      integer ina9(count200EBmax), ina10(count200EBmax), ina11(count200EBmax)
      real    ina1r(count200EBmax), ina12r(count200EBmax)
      real    ina9r(count200EBmax), ina10r(count200EBmax), ina11r(count200EBmax)
      real    keyEB(count200EBmax), keyEBk(count200EBmax)
      integer jdmax, inewyr, irecord

!.......................................................................................

!     input data are in milivolts, sf*1000 --> output in volts:
!     get user input (from command line)
      jrstart = 1           ! first record reading binary data, default
      jrend   = 100000      ! last  record reading binary data, default
      startid = 0.          ! option:  id# to keyEB off of:  output only a certain number
                            ! of the specified id#.  e.g. startid=100 (default=0:  
                            ! output all id's)

      ascii = .false.

!     set initial parameters, get command-line, etc.  

      lname = 'lsplit'
      ftypes= '1234'           ! default file types to process: 
                               ! 1=eddy1, 2=eddy2, 3=profile, 4=ground, 5=sonic
      include 'linitial7.for'  ! common initialization code 

      call infile_spl(yymmdd, ftypes, ifilen, nifiles)  

      print *
      print 9, nifiles 
9     format('  LSPLIT: processing ', i1, ' input files:' ) 
      do i=1,nifiles
      print *, '  ', ifilen(i)
      end do 

!     ==============================================================================================
!                                main program loop starts here
!     ==============================================================================================
      do infile=1, nifiles    ! loop over each input file
!
!     generic initializations
      startid= startidin         ! id# to start with (0 if none specified)
      nid=0                      ! init number of startid's found in file
      nidstart=1
      nidend=999
      if(startid.eq.0) then      ! don't allow output if we have to find startid first
      output= .true. 
      else 
      output= .false.
      if(startnum.ne.0) then 
      nidstart = int(startnum)   ! number of startid's to go by before reading
      nidend = int(10*(startnum-float(nidstart)))     ! read up to 9 labels
      print *, '  # before decimal = ',nidstart
      print *, '  # after  decimal = ',nidend
      nidend = nidstart + nidend 
      endif
      endif

!     get this loop's filename:
!     raw binary data input (*.dat files)
      dtype= ifilen(infile)(7:8)              ! dtype is data type string

      if (dtype .eq. 'EA') idtype = 1 
      if (dtype .eq. 'EB') idtype = 2 
      if (dtype .eq. 'PF') idtype = 3 
      if (dtype .eq. 'GD') idtype = 4 
      if (dtype .eq. 'SN') idtype = 5 
      if (dtype .eq. 'NT') idtype = 6 

      read (yymmdd,'(i10)') idate			    ! nrc 170207 add to hard wire problems with loggers 2012-2013
      idateMM = 100*(idate-(10000*int(idate/10000)))/100  ! nrc 181010 add to hard wire problems with loggers 2018
      print*, 'Month of the year, idateMM: ',imonth

      basename = yymmdd//dtype//'.dat'
      lrt = lentrim0(rootsplit, len(rootsplit))
      fname = rootsplit(1:lrt)//slash//yymmdd//slash//basename
      lfn = lentrim0(fname, len(fname))

      fnameout = dtype//yymmdd
      call getfullpath(rootsplit,yymmdd,dirsplit,fnameout,lfnout)
      
      fnamerr = dtype//yymmdd
      call getfullpath(rootprocout,yymmdd,dirprocout,fnamerr,lfnerr)

      print *, '  Opening input file "',fname(1:lfn)//'"' 
      print *, '  for output to "',fnameout(1:lfnout)//'.2*" files.'

      if(.not. ascii) then                                    ! default:  binary format
      open(11,file=fname(1:lfn),status='old',
     &       form='unformatted',access='direct',recl=4096) 
      else   
      open(11,file=fname(1:lfn),status='old') 
      endif

!......................................................................................
!    set-up output files:
!......................................................................................
      open(99,file=fnamerr(1:lfnerr)//'.err',status='unknown')  

      write(99,55) platform, lname, iyear, imonth, iday,ihrs,imin,isec
55    format(a6,' ',a7,':  on date (yy-mm-dd) ', i4,'-',i2.2,
     &       '-',i2.2,' (',i2.2,':',i2.2,':',i2.2,'), process:' )
      write(99,*) '   ', fname(1:lfn)  
      if(verbose) write(99,*) '  Report mode = verbose'

      if(ascii) write(99,*) '  (to be read as ascii file)'

      if (dtype .ne. 'SN') then  
      open(21,file=fnameout(1:lfnout)//'.201',status='unknown') 
      endif   
      if(startid.ne.0) then
      write(99,'(a25,i4,a10,f5.0)')'  Start with occurence # ',
     &      nidstart,' of label=',startid
      write(*, '(a25,i4,a10,f5.0)')'  Start with occurence # ',
     &      nidstart,' of label=',startid
      if(nidend.lt.999) then
      write(99,'(a23,i4)')'  End after occurence #',nidend
      write(*, '(a23,i4)')'  End after occurence #',nidend
      endif
      endif

      select case(idtype)
!.......................................................................................
      case(1:2)  
      enum=1
      if(dtype.eq.'EB') enum=2
      idsys = enum  ! set system id
      nf200 = 7     ! number of input fields (not incl line id#)
      nf201 = 15    ! includes 5:  day, hhmm, ss, status1, status2
      nf202 = 7     ! includes 3:  day, hhmm, ss  (remainder all spares)
      nf203 = 13    ! includes 3:        "        (remainder: cmstatus + spares)
      ns200 = 5     ! number of status bits to output on 200 line
      ns201 = 9     !                                 on 201 line
      do i=1,nf201-5       ! scale factor array (don't include day/time/status)
       sf201(i) = sf201E(i)
      enddo
       sf203(1) = sf203E(1)! only 1 var (cmstatus) on 203 line
      do i=1,nf203-3       ! debug test: use profile scale-factors
       sf203(i) = sf203PF(i)
      enddo
      assign 201 to f200 ! format label numbers for writing data to files
      assign 211 to f201
      open(20,file=fnameout(1:lfnout)//'.200',status='unknown')
!     200 file:
      write(20,01) enum,yymmdd
01    format('LBA Eddy',i1, ', ',a6,' fast channels (8 hz), 200 data')
      write(20,*)'julday  gmt.s  co2.v  h2o.v  u.m.s  v.m.s',
     &           '  w.m.s  t.c     stat sid cal'
!     format for data output (used in output section to match headers here):
201   format(i4,1x,f9.3,1x,2(f7.4,1x),3(f7.3,1x),f7.3,1x,i6,
     &       1x,3i1,1x,2i1) ! 200 line 
!     201 file:
      write(21,11) enum,yymmdd
11    format('LBA Eddy',i1, ', ',a6,' slow data (1/2 hz),',
     &       ' scaled to Volts')
      write(21,*)'julday  gmt.s '
      write(21,*)
     & '   pcon   flsam   flcal   pcell   rl   tcell   tdetect   tpump'
      write(21,*)
     & '   tamb   rh.tdew   sid   cal   mux   cm'
211   format(i4,1x,f9.3,/,8f8.3,/,f9.4,1x,f9.4,1x,3i1,1x,
     &       2i1,1x,3i1,1x,i1)

!.......................................................................................
!     case 2:  profile data:
      case(3)
      idsys = 3    
      nf200 = 7              ! not used for profile
      nf201 = 17             ! not incl id#, including 5:  day, hhmm, ss, status1,2
      nf202 = 23             ! not incl id#, including 3:  day, hhmm, ss
      nf203 = 13             ! not incl id#, including 3:  day, hhmm, ss
      ns201 = 24             ! number status bits to output for 201 file
      do i=1,nf201-5         ! scale factor array (not incl 5 date/time/status vars)
      sf201(i) = sf201PF(i)  ! changed from sf201E to sf201PF 180302
      enddo
      do i=1,nf202-3         ! scale factor array (not incl 3 date/time vars)
      sf202(i) = sf202PF(i)
      enddo
      do i=1,nf203-3         ! scale factor array (not incl 3 date/time vars)
      sf203(i) = sf203PF(i)
      enddo
      assign 212 to f201  ! format label number:  212 -> 201.2 = 201 label #2
      assign 222 to f202  ! 202 label #2 (after eddy)
      assign 232 to f203  ! 203 label #2
      open(22,file=fnameout(1:lfnout)//'.202',status='unknown') 
      open(23,file=fnameout(1:lfnout)//'.203',status='unknown')
!     201 file:
      write(21,12) yymmdd
12    format('LBA Profile 201 data, ',a6, 
     &       ': slow channels (1/2 hz), scaled to V')
      write(21,*) 'julday   gmt.s   '
      write(21,*) 'pcon   flsam   flcal   pcell   rl   tcell   ',
     &            'tdetect  tpump   '
      write(21,*) 'co2   h2o   spare   spare   sid   cal   mux   cc',
     &            '   csol   1234   5678   c' 
212   format(i4,1x,f9.3,/,8f8.3,/,4(f9.4,1x),2(3i1,1x,2i1,1x),      
     &          5i1,1x,2(4i1,1x),i1)                     ! profile 201
!     202 file:
      write(22,22) yymmdd
22    format('LBA Profile 202 data, ',a6, 
     & ': slow channels (1/2 hz), scaled to V')
      write(22,*) 'julday   gmt.s   '
      write(22,*) 'wdi   spare   netrad   nrad_cal   par1_up    ',
     &            'par1_dwn   '
      write(22,*) 'par2_up   par_cal  spare    spare    spare    ',
     &            'tair_cal   '
      write(22,*) 'tair_1   tair_2   tair_3   tair_4   tair_5   ',
     &            'tair_6   tair_7   tair_8'
222   format(i4,1x,f9.3,/,6(f8.3,1x)/,6(f8.3,1x),/,8(f8.3,1x) ) 
!     203 file:
      write(23,32) yymmdd
32    format('LBA Profile 203 data, ',a6,
     & ': slow channels (1/2 hz), scaled to V')
      write(23,*)'julday  gmt.s   '
      write(23,*)'ws1    ws2    ws3    ws4    spare    spare    ',
     &           'spare    spare    rain    spare'  
232   format(i4,1x,f9.3,/,10f7.3) ! profile 203

!.......................................................................................
!     case 3:  ground data
      case(4)  
      idsys = 4
      nf201 = 12
      nf202 = 27
      nf203 = 11
      nf300 = 4              ! power loss header
      ns201 = 18
      do i=1,nf201-5         ! scale factor array (not incl 5 date/time/status vars)
      sf201(i) = sf201GD(i)
      enddo
      do i=1,nf202-3         ! scale factor array
      sf202(i) = sf202GD(i)
      enddo
      do i=1,nf203-3         ! scale factor array
      sf203(i) = sf203GD(i)
      enddo
      assign 213 to f201     ! format label  (213-> 201.3, 201 label #3)
      assign 223 to f202     ! 222 -> 202 label #3
      assign 233 to f203     ! 232 -> 203 label #3
      open(22,file=fnameout(1:lfnout)//'.202',status='unknown') 
      open(23,file=fnameout(1:lfnout)//'.203',status='unknown')
!     write headers to output files:
!     201 file:
      write(21,13) yymmdd
13    format('LBA Ground 201 data, ',a6,
     & ': slow channels (1/2 hz), scaled to V')
      write(21,*) 'julday  gmt.s   '
      write(21,*) 'co   pcon_co   fl_co   stat_co   tcell_co   ttrap  ',
     &            'tpmp_d   sid   mux   gp21   cdpz   lhzi'
!     format for outputting ground 201 data corresponding to these headers:     
213   format(i4,1x,f9.3,/,7f8.3,1x,2(3i1,1x),3(4i1,1x))       ! ground 201
!     202 file:
      write(22,23) yymmdd
23    format('LBA Ground 202 data, ',a6,': 5-min avg, scaled to V')
      write(22,*)'julday  gmt.s    '
      write(22,*)'ts1    ts2    ts3    ts4    ts5    ts6    ts7    ts8 '
      write(22,*)'ts9    ts10   ts11   ts12   ts13   ts14   ts15   ts16'
      write(22,*)'par1   par2   par3   par4   par5   par6   par7   par8'
!     format for outputting ground 202 data corresponding to these headers:     
223   format(i4,1x,f8.0,3(/,8f8.3))         ! ground 202 (5-min avgs)
!     203 file:
      write(23,33) yymmdd
33    format('LBA Ground 203 data, ',a6,': 5-min avg, scaled to volts')
      write(23,*) 'julday   gmt.s    p_ss    p_hs   p_ms   p_ls   ',
     &            'p_cohs   p_cols   p_zag   pamb'

!     format for outputting ground 203 data corresponding to these headers:     
233   format(i4,1x,f8.0,8f8.3)             ! ground 203 (5-min avgs)

!.......................................................................................
!     case 4:  backup sonic data
      case(5)  ! sonic data : starting sept '2006  
      idsys=5
      nf200=5  
      ns200=5
      do i=1,nf200 
      sf200(i) = sf200SN(i)  
      enddo
      nf201=3
      assign 214 to f200 !               
      open(20,file=fnameout(1:lfnout)//'.200',status='unknown') 
!     200 file:
      write(20,14) yymmdd
14    format('LBA Sonic, ',a6,' fast channels (8 hz), 200 data')
      write(20,*)'julday   gmt.s   u.m.s   v.m.s   w.m.s   ',
     &           't.c   stat'
!     format for data output (used in output section to match headers here):
214   format(i4,1x,f9.3,1x,3(f7.3,1x),f7.3,1x,i6) ! 200 line for sonic
!     format for outputting sonic 200 data corresponding to headers:     
202   format(i4,1x,f8.2,1x,5(f8.4,1x),3i1,1x,i1)  ! 200 line for backup sonic

      case(6)
      idsys=6
      nf201=15     ! fill-in more when we know what natal system is like
      end select

!     ==============================================================================================
!     reads yymmdd.dat & convert cr10 binary data into ascii form
!     ==============================================================================================

      julday=-999
!     relevant for eddy data file:
      nfast=-1         ! number of fast data lines read (-1 if no 201 lines read yet)
      cmlast=10.0      ! last chilled-mirror status (10 v if none read yet)
!                      ! (from 203 input line but output as status bit on 201 line)
      noerr=0          ! number of 2-byte words with errors in current input block (0-2048)
      noblkerr=0       ! number of blocks with significant num of errors (ie noerr>=1024)
      toterr=0.0      ! total number of words with errors in whole file
      niderr=0         ! number invalid id# errors in this output block
      totiderr=0       ! total number of invalid id errors in this file
      ntierr=0         ! number of time-mismatch errors in this output block
      tottierr=0       ! total time-mismatch errors in this file
      norec=0          ! number of records (cr10-data lines) read
      n201=0           ! number of records (201-lines only) read
      nbadinrow=0      ! number of bad values or records in a row
      ngoodinrow=0     ! number of good 201-line records in a row
      tsave=0.
      k=0
      itmp=0
      n12=0  
      jrec=jrstart-1
      id=0
      indEB=1	       !counts the number of records on EB***.200
      indEB1=1	       !counts the number of records on EB***.201

!     get next 4096-byte record

100   jrec=jrec+1
      noerr=0          ! re-zero number of errors in current record
      if (jrec.gt.jrend) goto 980         ! normal exit of reading this file
      read (11,rec=jrec,iostat=eof) cr10  ! cr10: 4096 1-byte chars
      if (eof.lt.0) goto 990  
      do i=1,4095,2                       ! flip bytes within each 2-byte word
      ctmp=cr10(i+1)
      cr10(i+1)=cr10(i)
      cr10(i)=ctmp
      enddo

!     decodes binary cr10 record (in icr10) into ascii (in acr10) ==================================
!     j is index to icr10 array (the data input block)
!     k is index to acr10 ( the corresponding decoded output block)
!     in principle k more or less = j, except:
!     * 4-byte value increments j twice, k once
!     * when bad-data error, throw it out (increments j, not k)
!     * offset between 2048-array and data-lines

      j=0  
110   j=j+1
      if (j.gt.2048) then                                ! at end of this record
       if(verbose) then 
        write(99,111) jrec, icr10(1), icr10(2048), noerr ! , norec, julday, tlast
111     format('  Input block ',i5,': icr10(1),(2048)= ',i6,', ',i6,
     &       '; bad-data errors= ',i4)
       elseif(noerr.gt.0) then
        write(99,112) jrec, noerr
112     format('  Input block ',i5,': ',i5,
     &       ' of 2048 words have bad-data error')
       endif
       goto 100                         ! go to get next record
      endif
      if(ascii) then                   ! if ascii file, no decoding necessary
       k=k+1
       acr10(k) = rcr10(j)
       if (k .eq. 2048) goto 160
       goto 110
      endif
      if ( iand(icr10(j),  7168) .ne.   7168 ) goto 120  ! a 2-byte word of data 
      if ( iand(icr10(j),  -512) .eq.  -1024 ) goto 130  ! start of array
      if ( iand(icr10(j), 15360) .eq.   7168 ) goto 140  ! first 2 bytes of 4-byte value
      if ( iand(icr10(j), -1024) .eq.  15360 ) goto 150  ! last 2 bytes of 4-byte value

      noerr = noerr+1
      toterr = toterr+1.0

      if ((nblkerr.eq.0).or.(noerr.eq.1).or.(noerr.eq.1024)) then
       write(99,114)  noerr,jrec, j, icr10(j) 
114    format('  Bad data error #',i5,' in block jrec=',i5,
     &       ':  icr10(j =',i5,') = ',i6)
      endif
      if (noerr.eq.1024) nblkerr=nblkerr+1        ! number of records with 1024 errors or more.
      goto 110
!     2-bytes data:
120   k=k+1
      if ( iand(icr10(j), -32768) .ne. 0) then
       sf=-1.
      else
       sf=1.
      end if
      if (iand(icr10(j),16384).ne.0) sf=sf*0.01
      if (iand(icr10(j),8192).ne.0)  sf=sf*0.1
      acr10(k) = float(iand(icr10(j),8191))*sf
      if (k.ge.2048) goto 160
      goto 110
!     start of array:
130   k=k+1
      acr10(k) = float(iand(icr10(j),511))		! variable file id 20X
      if (k.ge.2048) goto 160
      goto 110
!     first 2 bytes of 4-byte value:
140   continue       ! the k=k+1
      itmp = icr10(j)
      n12 = n12+1    ! n12 is number of sequential first halves without a second half.  
      if (n12.eq.1) k=k+1 
      goto 110
!     third-4th bytes of a 4-byte value:
150   if (n12.eq.0) then 
      write(99,151) 
151   format(2x,'  Warning: found second half of 4-byte with no ',
     &          'first half')  
      goto 110       ! debug:  only if first value (added 9-25-01)
      endif  
      if (n12.gt.1) then 
      write(99,152) n12 
152   format(2x,'  Warning: found second half of 4-byte after ',i4,
     &          ' first halves')  
      endif  
      n12 = 0  
      if (iand(itmp,16384).ne.0) then
       sf = -1.
      else
       sf = 1.
      end if
      if (iand(itmp,512).ne.0) sf = sf*0.0001 
      if (iand(itmp,256).ne.0) sf = sf*0.01 
      if (iand(itmp,-32768).ne.0) sf = sf*0.1 
      acr10(k) = sf*float( 256*iand(icr10(j),256) + 
     &              256*iand(itmp, 255) + iand(icr10(j), 255) )
      itmp=0
      if (k.ge.2048) goto 160
      goto 110
!---------------------------------------------------------------------------------------------------
!     end of decode binary icr10 into ascii acr10
!---------------------------------------------------------------------------------------------------
!     next step:  have a full decoded output block (acr10), so break into data lines:
!---------------------------------------------------------------------------------------------------
160   jw1=1                           ! start at beginning of ascii array

165   continue
      if (jw1.ge.2049) then         ! at end, go get next ***changed from .eq. *** (1/23/02)
      k = 0
      if(verbose.or.noerr.gt.0) then 
      write(99,298) jrec, j, norec, julday, tlast
298   format('  Output block end (in block ',i5,', j=',i4,') at rec# ',
     &       i7,', day ',i3,', ',f9.3,' s')
      endif
      niderr=0                        ! re-zero number of id errors in record
      ntierr=0                        ! re-zero number of time-mismatch errors in record
      goto 110
      end if
      if (startid.ne.0) then          ! check if we are looking for specific labels
      if (acr10(jw1).eq.startid) then
      nid= nid+1
      endif
      if(nid.ge.nidstart .and. (.not.output) ) then
      output = .true.
      write(99,*)'  Starting data output.'
      write(*,*) '  Starting data output.'
      elseif(nid.gt.nidend) then
      write(99,*)'  Reached end of ID-block. Exiting this file...'
      write(*,*) '  Reached end of ID-block. Exiting this file...'
      goto 980                        ! exit reading this file
      endif
      endif

      if (acr10(jw1).eq.200.) then    ! data-line id = 200.  (eddy1,2 or backup sonic)
      jw2 = jw1+nf200                 ! check end of this data-line to see
      if (jw2.gt.2048) goto 170       ! if acr10 ends in middle of this line
      id = 200

      do i=1,nf200                    ! loop thru data fields on this line
      jw=jw1+i                        ! no date/time fields (compare to 201 etc lines)
      a200(i) = acr10(jw)
      end do
      jw1=jw2+1                           ! set jw1 to beg of next line
      if(output) goto 180             	  ! goto write 200 data
      elseif (acr10(jw1).eq.201.) then    ! data-line id = 201 (eddy1,2 or profile or ground)
      jw2=jw1+nf201
      if (jw2.gt.2048) goto 170
      id = 201
      julday = int(acr10(jw1+1))          ! first 5 vars are date/time/status
      hhmm   = int(acr10(jw1+2))
      ssec   = acr10(jw1+3)
      if (dtype.ne.'SN') then 
      status1 = int(acr10(jw1+4))
      status2 = int(acr10(jw1+5))

!.....nrc 170207.................................................................................
      if ((idate.ge.120321).and.(idate.lt.160601)) then	 ! after 170109 problems with the loggers in the years 2012-2013
       if (status1.eq.6999) status1 = 0
       if (status2.eq.6999) status2 = 1024
      end if
!................................................................................................

      do i=1,nf201-5                  ! n fields left to store less by 5 (date/time/status)
      jw = jw1+i+5                    ! jw index starts with 1 after 5 date/time/status vars
      a201(i) = acr10(jw) 
      end do
      endif  
      jw1 = jw2+1                          ! set jw1 to beg of next line
      if(output) goto 185                  ! goto write 201 data
      elseif (acr10(jw1).eq.202.) then   ! eddy, profile, ground with 202 line
      jw2 = jw1+nf202
      if (jw2.gt.2048) goto 170
      id     = 202
      julday = int( acr10(jw1+1) )
      hhmm   = int( acr10(jw1+2) )
      ssec   = acr10(jw1+3)
      do i=1,nf202-3
       jw      = jw1+i+3
       a202(i) = acr10(jw)
      end do
      jw1 = jw2+1                            ! set jw1 to beg of next line
      if(output) goto 190                  ! goto write 202 data
      elseif (acr10(jw1).eq.203.) then   ! eddy, profile, ground with 203 line
      jw2 = jw1+nf203
      if (jw2.gt.2048) goto 170
      id     = 203
      julday = int(acr10(jw1+1))
      hhmm   = int(acr10(jw1+2))
      ssec   = acr10(jw1+3)
      do i=1,nf203-3
      jw      = jw1+i+3
      a203(i) = acr10(jw)
      end do
      jw1 = jw2+1                         ! set jw1 to beg of next line
      if(output) goto 195                 ! goto write 203 data
      elseif (acr10(jw1) .eq. 100.) then  ! initialization loop header
      jw2 = jw1+5
      if (jw2 .gt. 2048) goto 170
      id     = 100
      iyear  = int(acr10(jw1+1))
      julday = int(acr10(jw1+2))
      hhmm   = int(acr10(jw1+3))
      hh     = int(hhmm/100)
      mm     = hhmm-hh*100
      ssec   = acr10(jw1+4)
      idsite = int(acr10(jw1+5))
      write (*,301) nid,iyear, julday, hh, mm, ssec,jrec, jw1
      write (99,301) nid,iyear, julday, hh, mm, ssec,jrec, jw1
301   format('  ID=100 (initialization header # ',i3,') at year=',i5,
     &    ' day:',i4,i3,':',i2,':',f5.2,/,9x,'(rec ',i5,',indx ',i4,')')
      jw1 = jw2+1
      goto 165                             ! get next data line
      elseif (acr10(jw1) .eq. 300.) then   ! low-power output
      jw2 = jw1+4
      if (jw2 .gt. 2048) goto 170
      julday  = int( acr10(jw1+1) )
      hhmm    = int( acr10(jw1+2) )
      ssec    = acr10(jw1+3)
      hh      = int(hhmm/100)
      mm      = hhmm-hh*100
      tpwroff = acr10(jw1+4)*2             ! time that power was off (min)
      write (*,303) julday, hh, mm, ssec, tpwroff
      write (99,303)julday, hh, mm, ssec, hh*3600+mm*60+ssec, tpwroff
303   format('  ID=300 (low-power header), at day=',i4,i3,':',i2,':',
     &       f5.2,' (,',f8.1,') after ',f6.0,' min. no power)')
      jw1 = jw2+1
      goto 165                             ! get next data line
      else
      write(99,309) acr10(jw1), jrec, jw1  ! unrecognized data line
309   format('  Invalid ID# (', e11.4, ' not = 100,200..203,300)',
     1       ' in/near input block # ', i6, ' (jw1= ', i4,')')
      niderr     = niderr+1
      totiderr   = totiderr+1
      nbadinrow  = nbadinrow+1             ! also incremented if bad rec (e.g. time mismatch) error 
      ngoodinrow = 0
      jw1 = jw1+1
      endif
      goto 165                       ! get next data line
!     goto 170 (from above) if jw2>2048 (at end of 2048-integer record), but have part of data-line.
170   k = 0      
      if (jw1 .gt. 2048) goto 110    ! if <=, we're still in the middle of a data-line
      do js=jw1,2048           
      k = k+1                        ! save remaining elements of acr10, move
      acr10(k) = acr10(js)           ! to beginning, and set k = # remaining
      end do
      goto 110                       ! go read another record from file


!---------------------------------------------------------------------------------------------------
!      write ascii data into output files:
!---------------------------------------------------------------------------------------------------
!     4 hz data :: id 200 
!     (eddy systems and backup sonic)
!---------------------------------------------------------------------------------------------------
180   norec=norec+1
      if(nfast.lt.0) goto 165               ! abort if no time markers read yet
      if (lskip) goto 165                   ! or if last time-marker no good
      nfast=nfast+1
      call vecdiv(a200,sf200,v200,nf200)    ! return: v200 = a200 / sf200
      tt=tlast+float(nfast-1)/freqhi        ! tt (secs) interpolated from last 20n line

!...................................................................................................
!     note:
!     the eddy system was ~8.07 min behind... the midnight crossings are fixed later (nrc).
!...................................................................................................
      if ((dtype.eq.'EB').and.(idate.ge.080729).and.(idate.le.081219)) then
       tt=tt+615     !+522
      endif 
      if ((idate.ge.081220).and.(idate.le.090414)) then
       tt=tt+610     !+522
      endif 
!...................................................................................................
!     note:
!     the eddy system was 2:54 min behind... the midnight crossings are fixed later (nrc).
!...................................................................................................
      if ((idate.ge.090415).and.(idate.le.090502)) then
       tt = tt-15
      endif 
!...................................................................................................
!     note:
!     the eddy system was 2:54 min behind... the midnight crossings are fixed later (nrc).
!...................................................................................................
      if ((idate.ge.090503).and.(idate.le.090522)) then
       tt = tt-120
      endif 
      if ((idate.ge.090523).and.(idate.le.091101)) then
       tt = tt+140
      endif 
      if ((idate.ge.091101).and.(idate.le.091201)) then
       tt = tt+100
      endif 
      if ((idate.ge.091201).and.(idate.le.091231)) then
       tt = tt+0
      endif 

!...................................................................................................
182   format('  200-line time-error:  200-line time (secs)   =',f9.3,/,
     &       '                        interpolation from 201 =',f9.3,
     &       ' (diff =',f7.3,')')
!...................................................................................................
      jd=julday
      if(tt.ge.86400.) then             ! fix midnight crossings
       jd = jd+1
       tt = tt-86400.
       if ((jd.eq.367).and.(((int(year/4))-(year/4)).eq.0)) then
         print*, 'leap year'
         jd = 1    	   ! ewg: fix end-of-year crossing 
         year = year+1    ! nrc: fix end-of-year crossing 
       end if
       if ((jd.eq.366).and.(((int(year/4))-(year/4)).ne.0)) then
         jd = 1     	   ! ewg: fix end-of-year crossing 
         year = year+1    ! nrc: fix end-of-year crossing 
       end if
      endif
      nf = nf200
      if((dtype.eq.'EA').or.(dtype.eq.'EB').or.(dtype.eq.'SN')) then  
       icsat=int(a200(nf200))
       nf = nf200-1                      ! don't count csat status as real
      endif

!...................................................................................................
      if ((is200(3).eq.1).and.(idate.ge.080101)) then
        print*, ' Eddy system incorrectly set to EA -negative H2O and CO2 voltages during:',idate
        idsys = 2	       ! nrc EB is now system 1
        enum  = 2	
        dtype ='EB'
        is200(1) = 0;  is200(2) = 1;  is200(3) = 0	        
      end if
!...................................................................................................

!...................................................................................................

!...................................................................................................
!     note:
!     the eddy system co2 voltage zero was -1. forced the timeseries to be closer to 0
!     licor 210401 to 220401 nrc
!...................................................................................................
!      if ((idate.ge.210401).and.(idate.le.210610)) then
!       v200(1) = 2*(v200(1)+1)       !co2
!       v200(2) = 2*(v200(2)+0.5)     !h2o
!      endif 
!      if ((idate.ge.210610).and.(idate.le.210729)) then
!       v200(1) = 1.5*(v200(1)+1.)         !co2
!       v200(2) = 1.5*(v200(2)+1.)         !h2o
!      endif 
!      if ((idate.ge.210729).and.(idate.le.310729)) then
!       v200(1) = 1.5*(v200(1)+1.)         !co2
!       v200(2) = 2.0*(v200(2)+1.1)         !h2o
!      endif 
! IRGA's voltage correlation coefficients were reset to zero after lightning
! Hardwired IRG3-1046  
! C2 T=32.82, C2 K=18009, C2 A=1.36268E-1, C2 B=4.64832E-6, C2 C=8.20996E-9, C2 D = -1.04493E-12, C2 E=6.3323E-17
! H2O cal values H2 T = 33.22, H2 K=15621, H2 A= 6.00539E-3, H2 B = 2.50887E-6, H2 C = 3.11171E-11  
! F(V) = a 1V + a2V2 + a3V3 + a4V4 + a5V5  for H2O
! µmol/mol = A(mV) + B(mV) 2 + C(mV)3 + D(mV) 4 + E(mV)5.
!      if ((idate.ge.210730).and.(idate.lt.211231).and.(v200(1).lt.0)) then 
!	v200(1) = (clicor(1)*v200(1))+(clicor(2)*v200(1)*v200(1))+(clicor(3)*v200(1)*v200(1)*v200(1))+(clicor(4)*v200(1)*v200(1)*v200(1)*v200(1))+!!(clicor(5)*v200(1)*v200(1)*v200(1)*v200(1)*v200(1)*v200(1))
!      	v200(2) = (hlicor(1)*v200(2))+(hlicor(2)*v200(2)*v200(2))+(hlicor(3)*v200(2)*v200(2)*v200(2))
!      endif
!      if ((dtype.eq.'EA').or.(dtype.eq.'EB')) then 
!       write(20,f200)jd,tt,(v200(i),i=1,nf),icsat,(is200(i),i=1,ns200)
!      endif 
!...................................................................................................

      if ((dtype.eq.'EA').or.(dtype.eq.'EB')) then 
       indEB=indEB+1                     ! nrc file length counter added 

      !needed to hardwire positive voltages 111219 - 120104
        if ((idate.ge.111219).and.(idate.le.120104)) v200(1) = -1.*v200(1)        

      !set boundaries of what is possible for u.m.s  v.m.s  w.m.s  t.c given the write format constrains
        if ((v200(1).lt.-9.99).or.(v200(1).gt.99.99))   v200(1)=-9.99
        if ((v200(2).lt.-9.99).or.(v200(2).gt.99.99))   v200(2)=-9.99
        if ((v200(3).lt.-99.99).or.(v200(3).gt.999.99)) v200(3)=-99.99
        if ((v200(4).lt.-99.99).or.(v200(4).gt.999.99)) v200(4)=-99.99
        if ((v200(5).lt.-99.99).or.(v200(5).gt.999.99)) v200(5)=-99.99
        if ((v200(6).lt.-99.99).or.(v200(6).gt.999.99)) v200(6)=-99.99
       
!This has been forced, should be properly fixed
      if (icsat.ge.4032)  icsat=icsat-4032
       write(20,f200)jd,tt,(v200(i),i=1,nf),icsat,(is200(i),i=1,ns200)
      endif 

      if (dtype.eq.'SN') then 
       write(20,f200)jd,tt,(v200(i),i=1,nf),icsat
      endif 
!     f200  check file-specific initializations for definition of format labels f200
      goto 165
!---------------------------------------------------------------------------------------------------
!     0.5 hz data :: id 201
!     (eddy, profile, ground, and sonic)
!---------------------------------------------------------------------------------------------------
185   norec=norec+1
      n201 = n201+1  				! num 201 lines
      if(n201.le.1) lastday=julday  		! initialize
      lskip=.false.
      nfast=0
      hh=int(hhmm/100)
      mm=hhmm-hh*100
      tlast=float(hh*3600 + mm*60) + ssec  	! time in secs only
      if (dtype .eq. 'SN') goto 165  
!     error checking moved further below:
!     2 status words (12-bit) in every 201 line:
!     index# 12    11    10     9      8     7     6     5     4     3     2     1
!     2048  1024   512   256    128    64    32    16     8     4     2     1
!     eddy:
!     #1 spare spare spare spare  spare spare spare spare spare spare spare spare
!     #2   sd3   sd2   sd1  cal2   cal1  mux2  mux1  mux0 spare spare spare spare 
!     profile:
!     #1    ms    hs   ss  lvl_9  lvl_8 lvl_7 lvl_6 lvl_5 lvl_4 lvl_3 lvl_2 lvl_1
!     #2   sd3   sd2   sd1  cal2   cal1  mux2  mux1  mux0 calc2 calc1 zero  ls
!     ground:
!     #1 pwrco pwrdp pwrpc pwrzag spare spare spare spare co_ls co_hs co_ze co_in
!     #2   sd3   sd2   sd1 power  spare  mux2  mux1  mux0 spare pwrpf pwre2 pwre1
!
!     status will be output as follows
!     is200 (taken from last 201 line)
!     e*  5 bits:  1-3=sysdesig; 4-5=cal
!     sn  1 bit:   1= sonic status
!     is201:
!     E* 9 bits:   1-3=sysdesig; 4-5=cal; 6-8=mux; 9=cm hygr. status
!     PF 24 bits:  1-3=sysdesig; 4-5=cal; 6-8=mux; 9-10=calcon; 11-15=cal sols
!                  16-24= level sols 
!     GD 17 bits:  1-3=sysdesig; 4-6=mux; 7-13=pwr stat; 14-17=co cal sols
!     SN  4 bits:  1-3=sysdesig; 4=sonic status
!
!     extract status bits and put into integer arrays:
      call istat(status1,isol1)    		! status bits from 201 line
      call istat(status2,isol2)

      nf=nf201-5   	                        ! leave out day, hhmm, ss, status1, status2
!....................................................................................
!     eddy 1 or 2
      if(dtype(1:1).eq.'E') then 
       indEB1 = indEB1+1                 ! nrc file length counter added  
       do i=1,5
        is200(i)   = isol2(13-i)         ! sysid 3,2,1;  cal2,cal1
        is201(i)   = is200(i)            ! 1st 5 bits same for 200, 201
        is201(i+5) = isol2(8-i)          ! mux2, mux1, mux0, + 2 spare bits
        ! note:  last 2 (9,10) are spares; 9 replaced below
       enddo
       is201(9) = 1 ! chilled-mirror status:  assume high=ok
       if(cmlast.lt.2.5)     is201(9) = 0   ! low if mirror dirty
       if(cmlast.ge.9.9)     is201(9) = 2   ! no 203 lines read yet
       if (idate.ge.080101)  is201(9) = 1   ! chilled-mirror status:  set ok after 2008
       ! note:  cmlast taken from last eddy 203 line
      endif
!....................................................................................
!     profile
      if(dtype.eq.'PF') then 
       do i=1,12
        is201(i) = isol2(13-i)        ! sysid 3,2,1; cal2,1; mux2,1,0; calc2,1; zero, ls
       enddo
       do i=13,15
        is201(i) = isol1(12+(13-i))   ! ms,hs,ss
       enddo
       do i=16,24
        is201(i) = isol1(i-15)        ! level 1-9, column
       enddo
      end if
!....................................................................................
!     ground 
      if(dtype.eq.'GD') then 
       do i=1,3 
        is201(i)    = isol2(13-i)   ! sysid 3,2,1; 
        is201(i+3)  = isol2(8-i)    ! mux 2,1,0
        is201(i+7)  = isol2(4-i)    ! pwr: PF, e2, e1
        is201(i+10) = isol1(13-i)   ! pwr: co, dumppump, pc
        is201(i+14) = isol1(5-i)    ! co sol:  ls, hs, zag
       enddo
       is201(7)  = isol2(9)         ! ground (main) power
       is201(14) = isol1(9)         ! pwr: zag
       is201(18) = isol1(1)         ! solenoid: co inlet
      endif
!....................................................................................
!     backup sonic
      if(dtype.eq.'SN') then 
       do i=1,3
        is201(i) = isol2(13-i)
       enddo
       is201(4)  = isol1(2)  ! sonic status
       is200(1)  = isol1(2)
      endif
!....................................................................................
      id201 = is201(1)*4 + is201(2)*2 + is201(3)               ! get system id stored in this line
!      if ((dtype.eq.'EB').and.(id201.ne.1)) is201(1)=0        ! needs to check why this is happening CHECK
!     1=eddy1, 2=eddy2, 3=profile, 4=ground, etc
!       print*, ' check id201: ', id201
!       print*, ' check is201(1): ', is201(1)
!       print*, ' check is201(2): ', is201(2)
!       print*, ' check is201(3): ', is201(3)
!...............................................................................
      if ((id201.eq.1).and.(idate.ge.080101)) then
        print*, ' Eddy system incorrectly set to EA -negative H2O and CO2 voltages during:',idate
        id201 = 2	       ! nrc EB is now system 1
        dtype ='EB'
        is201(1) = 0;  is201(2) = 1;  is201(3) = 0	        
      end if
!...............................................................................
!end nrc 171205.................................................................
!     ready to output, first check if this record ok:
      recerr = checkrec(julday, tlast, lastday, tprev, -0.5, jdmax)   ! check: record time ok?
!     elapsed time must be greater than 0
      if(id201.ne.idsys)  recerr=1.     ! if sys id on 201 line not right, skip also
       if (recerr.ne.0) then
       write(99,350) 201, norec, jrec, jw1-nf201-1, julday, tlast, 
     &              lastday, tprev
350     format(' ',i3,': mismatch, rec #',i7,' (block #', i5,
     &       ', jw1=',i5,'): jday ',i4,
     &       ', ',f9.3,' sec (prev= ',i4,1x,f9.3,')')
        write(99,*)'   (num 201 = ',n201,'; sysid =',id201,')'
        ntierr   = ntierr+1
        tottierr = tottierr+1
!     now, leave only for id error:      
        lskip = .true.
        goto 165
       endif
       ngoodinrow = ngoodinrow+1
!     if multiple bad lines in a row, wait for more than one good line before outputting
      if(nbadinrow.gt.1) then
      if(ngoodinrow.ge.1) nbadinrow = 0 ! will output next time around
      write(99,352) 201, norec, jrec, jw1-nf201-1, julday, tlast
352   format('  ',i3,': OK, but no output, rec #',i7,' (block #', i5,
     &       ', jw1=',i5,'): jday ',i4,
     &        ', ',f9.3,' sec')
       write(99,*)'  (want 2 good 201 lines in a row, now have ',
     &           ngoodinrow,')'
       lskip = .true.
       goto 165
      endif
      nbadinrow = 0
      lastday   = julday                ! keep record of this day
      tprev     = tlast                 ! keep record of this time (secs)
      call vecdiv(a201,sf201,v201,nf)   ! return: v201 = a201 / sf201
      jd = julday

!...................................................................................................
!nrc 171205.........................................................................................
!...................................................................................................
!     note:
!     The eddy system was ~8.07 min behind the profile during 2008. Therefore, the calibration comes first at the eddy
!     and then to the profile... the midnight crossings are fixed later (nrc).
!...................................................................................................
      tt = tlast
!...................................................................................................
!     note:
!     the eddy system was ~8.07 min behind... the midnight crossings are fixed later (nrc).
!...................................................................................................
      if ((dtype.eq.'EB').and.(idate.ge.080729).and.(idate.le.081219)) then
       tt=tt+615     !+522
      endif 
      if ((idate.ge.081220).and.(idate.le.090414)) then
       tt=tt+610     !+522
      endif 
!...................................................................................................
!     note:
!     the eddy system was 2:54 min behind... the midnight crossings are fixed later (nrc).
!...................................................................................................
      if ((idate.ge.090415).and.(idate.le.090502)) then
       tt = tt-15
      endif 
!...................................................................................................
!     note:
!     the eddy system was 2:54 min behind... the midnight crossings are fixed later (nrc).
!...................................................................................................
      if ((idate.ge.090503).and.(idate.le.090522)) then
       tt = tt-120
      endif 
      if ((idate.ge.090523).and.(idate.le.091101)) then
       tt = tt+140
      endif 
      if ((idate.ge.091101).and.(idate.le.091201)) then
       tt = tt+100
      endif 
      if ((idate.ge.091201).and.(idate.le.091231)) then
       tt = tt+0
      endif 
      
!...................................................................................................
      if(tt.ge.86400.) then             ! fix midnight crossings
       jd = jd+1
       tt = tt-86400.
       endif                            ! nrc 210521
       if ((jd.eq.367).and.(((int(year/4))-(year/4)).eq.0)) then
         jd = 1    	   	! ewg: fix end-of-year crossing 
         year = year+1    	! nrc: fix end-of-year crossing 
       end if
       if ((jd.eq.366).and.(((int(year/4))-(year/4)).ne.0)) then
         jd = 1    	   	! ewg: fix end-of-year crossing 
         year = year+1    	! nrc: fix end-of-year crossing 
       end if
!      endif
!...............................................................................
      if ((id201.eq.1).and.(idate.ge.080101)) then
        print*, ' Eddy system incorrectly set to EA -negative H2O and CO2 voltages during:',idate
        id201 = 2	       ! nrc EB is now system 1
        dtype ='EB'
        is201(1) = 0;  is201(2) = 1;  is201(3) = 0	        
      end if
!...............................................................................
!end nrc 171205.................................................................

 !     !...................................................................................................
 !     ! In 210401 we changed flow controler --changed treshold flux  				!nrc 210805
 !     !............................................................................................
 !     if ((idate.ge.210729).and.(idate.le.310729)) then
 !      v201(3) = v201(3)-0.03         !flux      1:pcon   2:flsam   3:flcal   4: pcell  
 !     endif 
      !............................................................................................


!      write(21,f201)jd,tlast,(v201(i),i=1,nf),(is201(i),i=1,ns201)
       if ((tt.le.86400).and.(tt.ge.0).and.(jd.le.366).and.(jd.ge.1)) then
        write(21,f201)jd,tt,(v201(i),i=1,nf),(is201(i),i=1,ns201)
       endif

!     f201  check file-specific initializations for output formats defined by f201
      if(verbose.and.(ntierr.gt.0 .or. niderr.gt.0) ) then
      write(99,354) 201, norec, jrec, jw1-nf201-1, julday, tlast
354   format('  ',i3,' CR10 record #',i7,' (block #', i5,', jw1=',i5,
     &       '): jday ',i3,
     &            ', ',f9.3,' sec' )
      endif                
      goto 165

!---------------------------------------------------------------------------------------------------
!     0.5 hz data :: id 202
!     (eddy and profile and ground)
!---------------------------------------------------------------------------------------------------
190   norec = norec+1
      lskip = .false.
      nfast = 0
      hh    = int(hhmm/100)
      mm    = hhmm-hh*100
      tlast = float(hh*3600 + mm*60) + ssec  ! time in secs only
!     do not write repeated data (due, e.g., to errors in cr10 data dump):
      recerr = checkrec(julday, tlast, lastday, tprev,-0.5,jdmax)  ! check: record ok?
!     (elapsed time (from 201) must be > -0.5 sec)
      if(id201.ne.idsys)  recerr=1.   ! if sys id on last 201 line was not right, skip also
      if(recerr.ne.0.) then
       write(99,350) 202, norec, jrec, jw1-nf202-1, julday, tlast, 
     &              lastday, tprev
       write(99,*)'   (last 201 # = ',n201,'; sysid =',id201,')'
       lskip = .true.
       goto 165
      endif
! if several bad 201 lines in a row, don't output until multiple good 201 lines found
      if(nbadinrow.gt.1) then   
       write(99,352) 202, norec, jrec, jw1-nf201-1, julday, tlast
       write(99,*)'  (Want 2 good 201 lines in a row, now have ',
     &           ngoodinrow,')'
       lskip = .true.
       goto 165
      endif
      if(verbose .and. (ntierr.gt.0 .or. niderr.gt.0) ) then
       write(99,354) 202, norec, jrec, jw1-nf201-1, julday, tlast
      endif
      if(dtype(1:1).eq.'E') then
       goto 165  ! eddy 202 line is all spares, skip output
      endif
      nf = nf202-3  ! nb:  3 less -- don't do day, hhmm, sec
      call vecdiv(a202,sf202,v202,nf)   !  return: v202 = a202 / sf202
      write(22,f202)julday,tlast,(v202(i),i=1,nf)  ! changed file unit # from 202
!     f202  check file-specific initializations for output formats defined by f202
      goto 165

!---------------------------------------------------------------------------------------------------
!     0.5 hz data :: id 203
!     (eddy and profile and ground systems)
!---------------------------------------------------------------------------------------------------
195   norec = norec+1
      lskip = .false.
      nfast = 0
      hh    = int(hhmm/100)
      mm    = hhmm-hh*100
      tlast = float(hh*3600 + mm*60) + ssec  ! time in secs only
!     do not write repeated data (due, e.g., to errors in cr10 data dump):
      recerr = checkrec(julday, tlast, lastday, tprev,-0.5,jdmax)  ! check: record ok?
      ! (elapsed time (from 201) must be > -0.5 sec)
      if(id201.ne.idsys)  recerr=1.   ! if sys id on last 201 line was not right, skip also
      if( recerr.ne.0.) then
      write(99,350) 203, norec, jrec, jw1-nf202-1, julday, tlast, 
     &              lastday, tprev
      write(99,*)'   (last 201 # = ',n201,'; sysid =',id201,')'
      lskip = .true.
      goto 165
      endif
      ! if several bad 201 lines in a row, don't output until multiple good 201 lines found
      if(nbadinrow.gt.1) then   
      write(99,352) 203, norec, jrec, jw1-nf201-1, julday, tlast
      write(99,*)'  (Want 2 good 201 lines in a row, now have ',
     &           ngoodinrow,')'
      lskip = .true.
      goto 165
      endif

!     don't let non-201 lines affect future timestamp criteria:
      if(verbose.and.(ntierr.gt.0.or.niderr.gt.0)) then
      write(99,354) 203, norec, jrec, jw1-nf201-1, julday, tlast
      endif
      if(dtype(1:1).eq.'E') then 
       cmlast = a203(1)/sf203(1)  ! save chilled mir status 
      ! for output on next 201 line
       goto 165  ! remainder of eddy 203 line is all spares, skip output
      endif
      nf=nf203-3  ! nb:  3 less -- don't do day, hhmm, sec
      call vecdiv(a203,sf203,v203,nf)   !  return: v203 = a203 / sf203
      write(23,f203)julday,tlast,(v203(i),i=1,nf)
!     f203  check file-specific initializations for output formats defined by f203

      goto 165

!---------------------------------------------------------------------------------------------------
!    error exit (defunct)
!---------------------------------------------------------------------------------------------------
900   write(99,*) ' '
      write(99,901) ierrlim, ifilen(infile)
901   format     ('  Error exit (num errors >=',i5,') from ', a12)
      write(99,*) '  Abort processing: '
      write(99,*) '  at record (4096 bytes/rec) # jrec =',jrec
      write(99,*) '  # lines (cr10 record) read: norec = ',norec
      write(99,*) '  # errors reading, decoding: noerr = ',noerr
!     terminal output:
      write(*,901) ierrlim, ifilen(infile)
      write(*,*)  '  Abort processing at error #       = ', noerr
      goto 992

!---------------------------------------------------------------------------------------------------
!     normal exit
!     user-defined exit (record #jrend reached)
!---------------------------------------------------------------------------------------------------

980   write(99,*) ' '
      write(99,*) '  Warning:  exit from ', ifilen(infile)
      write(99,*) '  because record limit exceeded at record ',
     &            '(4096 bytes/rec) # jrec =',jrec
      write(99,*) '  Blocks (of 2048 2-byte words) read: jrec = ',jrec
      write(99,*) '  CR10 data-lines read:              norec = ',norec
      write(99,*) '  Number of 201 cr10 lines            n201 = ',n201
      write(99,*) '  Total # of 2-byte words read = jrec*2048 = ',
     &             jrec*2048
      write(99,*) '  Error summary (2-byte words): '
      write(99,*) '  Number of blocks with >=1024 bad-data errors = ',
     &             nblkerr,' (',nblkerr/jrec*100, '%)'
      write(99,*) '  Total number of bad-data errors          = ',
     &             toterr,' (',toterr/(jrec*2048)*100,'%)'
      write(99,*) '  Error summary (within records): '
      write(99,*) '  Number of invalid-id errors              = ',
     &             totiderr
      write(99,*) '  Number of 201 time-base/mismatch errors  = ',
     &             tottierr
      goto 992

!     exit due to end-of-file:
990   write(99,*) ' '
      write(99,*) '  Normal exit at end of ',ifilen(infile)
      write(99,*) '  (reached end-of-file)'
      write(99,*) '  Blocks (of 2048 2-byte words) read: jrec     = ',
     &               jrec
      write(99,*) '  CR10 data-lines (records) read:    norec     = ',
     &               norec
      write(99,*) '  Number of 201-line records          n201     = ',
     &               n201
      write(99,*) '  Total # of 2-byte words read = jrec*2048     = ',
     &               jrec*2048
      write(99,*) '  Error summary (2-byte words): '
      write(99,*) '  Number of blocks with >=1024 bad-data errors = ',
     &               nblkerr,' (',nblkerr/jrec*100, '%)'
      write(99,*) '  Total number of bad-data errors              = ',
     &               toterr,' (',toterr/(jrec*2048)*100,'%)'
      write(99,*) '  Error summary (within records): '
      write(99,*) '  Number of invalid-id errors                  = ',
     &               totiderr
      write(99,*) '  Number of 201 time-base/mismatch errors      = ',
     &               tottierr

992   close(99)     

!     end of main program loop:
995   end do ! do infile=1, nifiles (loop to next input file)

!---------------------------------------------------------------------------------------------------
!     nrc:160902 	opens the input file (E*.200) and sorts it
!     The ugliest program I have ever written -nrc 160920
!---------------------------------------------------------------------------------------------------
      if ((idate.ge.150101).and.(idate.le.160530)) then	      ! after 170109 CC is wrong, needs to be the same as CAL
!       print *,'  Reading E*.200 and sorting its time vector dtype: ',dtype
       dtype  = 'EB'
       enum   = 2
       infile = 1
       fnameout = dtype//yymmdd
       
       call getfullpath( rootsplit,yymmdd, dirsplit, fnameout, lfnin) 
       print *,' file name:', fnameout(1:lfnout)//'.200'
       open(44,file=fnameout(1:lfnout)//'.200',status='old')  
       read(44,*)                                                            ! header info
       read(44,*)
       i = 1
       do while(i.le.indEB)
510     continue
        read(44,511,end=512) ina1(i), var2(i), var3(i), var4(i), 
     &                     var5(i), var6(i), var7(i), var8(i),  
     &                     ina9(i), ina10(i), ina11(i)
511     format(i4,1x,f9.3,1x,2(f7.4,1x),3(f7.3,1x),f7.3,1x,
     &             i6,1x,i3.3,1x,i2.2) ! 200 line 
!        print*,'ina9(i): ',ina9(i), ' ina10(i): ', ina10(i), ' ina11(i): ',ina11(i)
        i=i+1
       end do
512    continue   
       close(44, status='delete')
      
      !creates an index for it to be sorted and converts to reals integer 
       jdmax = 365
       if(((int(year/4))-(year/4)).eq.0) jdmax=366
       print*,'Maximum numbers of days in  year: ',year, ' jdmax: ', jdmax

       inewyr=0
       i=1
       do while (i.le.indEB)
        if (ina1(i).eq.1)  then				! flags start of the year
!        if ((ina1(i).eq.1).and.(idateMM.gt.11))  then	! flags start of the year
!         ina1(i)=ina1(i)+jdmax			        ! adds the last JD e.g. 365 or 366
         inewyr=1
        end if     

        if (ina1(i).gt.jdmax) then		         ! flags start of the year
         ina1(i) = ina1(i)-jdmax		 	 ! substracts the last JD e.g. 365 or 366
        end if

        ina1r(i)  = dble(ina1(i))      !float
        ina9r(i)  = dble(ina9(i))
        ina10r(i) = dble(ina10(i))
        ina11r(i) = dble(ina11(i))
        keyEB(i)  = dble(var2(i))+dble(100000.000*(ina1r(i)-ina1r(1)))  ! counter id to sort the file
        if (ina1(i).lt.0.)  keyEB(i)=990000.000+dble(var2(i))		! kick the strange values to the end
        if ((inewyr.eq.1.).and.(ina1r(i).lt.10)) keyEB(i)=keyEB(i)+(100000.00*370.)
        if (keyEB(i).eq.0.) keyEB(i)=990000.000+dble(var2(i))		! kick the strange values to the end
        keyEBk(i) = keyEB(i)
        i=i+1
       end do
!       print *,'  keyEBk: ',keyEBk
!       pause
      
       i=1
       irecord=0
       do while (i.le.n-1)
        if (keyEBk(i).gt.keyEBk(i+1)) irecord=1+irecord
        i=i+1
       end do
       print*, 'Times the gmt time vector goes back and forth: ', irecord

       call sort2(indEB,keyEB,var2) 
       call back2(indEB,keyEB,keyEBk)
       call sort2(indEB,keyEB,var3) 
       call back2(indEB,keyEB,keyEBk)
       call sort2(indEB,keyEB,var4) 
       call back2(indEB,keyEB,keyEBk)
       call sort2(indEB,keyEB,var5) 
       call back2(indEB,keyEB,keyEBk)
       call sort2(indEB,keyEB,var6) 
       call back2(indEB,keyEB,keyEBk)
       call sort2(indEB,keyEB,var7) 
       call back2(indEB,keyEB,keyEBk)
       call sort2(indEB,keyEB,var8) 
       call back2(indEB,keyEB,keyEBk)
       call sort2(indEB,keyEB,ina1r) 
       call back2(indEB,keyEB,keyEBk)
       call sort2(indEB,keyEB,ina9r) 
       call back2(indEB,keyEB,keyEBk)
       call sort2(indEB,keyEB,ina10r) 
       call back2(indEB,keyEB,keyEBk)
       call sort2(indEB,keyEB,ina11r) 
       call back2(indEB,keyEB,keyEBk)

       print*,'indEB: ',indEB
       i=1
       do while (i.le.indEB)
        ina1(i)  = int(ina1r(i))
        ina9(i)  = int(ina9r(i))
        ina10(i) = int(ina10r(i))
        ina11(i) = int(ina11r(i))
        i=i+1
!        print*,'ina9(i): ',ina9(i), ' ina10(i): ', ina10(i), ' ina11(i): ',ina11(i)
       end do

       i=1
       if (inewyr.eq.1) then
        do while (i.le.indEB)
         if (ina1(i).gt.jdmax) then			       !flags start of the year
          ina1(i)=ina1(i)-jdmax				       !substracts the last JD e.g. 365 or 366
         endif
         i=i+1
        enddo
       endif

!!     200 file:
       open(44,file=fnameout(1:lfnout)//'.200', status='new')  
       write(44,613) enum, yymmdd
613    format('LBA Eddy',i1, ', ',a6,' fast channels (8 hz), 200 data')
       write(44,*)'julday  gmt.s  co2.v  h2o.v  u.m.s  v.m.s',
     &            '  w.m.s  t.c     stat sid cal'
       
       i=1
       do 612 while(i.le.indEB)
        if (ina1(i).eq.0) goto 615		! skips jd==0
        write(44,614)ina1(i),var2(i),var3(i),var4(i),var5(i),var6(i),var7(i),var8(i),
     &             int(ina9(i)),int(ina10(i)),int(ina11(i))
614     format(i4,1x,f9.3,1x,2(f7.4,1x),3(f7.3,1x),f7.3,1x,
     &             i6,1x,i3.3,1x,i2.2)		! 200 line 
615     i=i+1
612     continue   
        close(44)
!       end do
       print *, '  File EB*.200 is sorted'
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!     nrc:160902 	Opens the input file (E*.201) and sorts it
!---------------------------------------------------------------------------------------------------
       print *,'  Reading E*.201 and sorting its time vector'
       dtype='EB'
       enum=2
       infile=1
       fnameout = dtype//yymmdd

       call getfullpath( rootsplit,yymmdd, dirsplit, fnameout, lfnin) 
       print *,' file name:', fnameout(1:lfnout)//'.201'
       open(45,file=fnameout(1:lfnout)//'.201',status='old')  
       read(45,*)                                                            ! 4 lines header info
       read(45,*)
       read(45,*)
       read(45,*)
       i=1
       do while(i.le.indEB1)
710     continue
        read(45,711,end=712) ina1(i),var2(i),  
     &                     var3(i),var4(i),var5(i),var6(i),var7(i),var8(i),var13(i),var14(i),
     &                     var15(i), var16(i),ina9(i), ina10(i), ina11(i), ina12(i)
711     format(i4,1x,f9.3,/,8f8.3,/,f9.4,1x,f9.4,1x,i3,1x,
     &       i2,1x,i3,1x,i1) ! 201 line 
        i=i+1
       enddo
712    continue   
       close(45, status='delete')
      
!      jdmax=valmax(indEB1,ina1)
      !creates an index for it to be sorted and converts to reals integer 
       jdmax=365
       if(((int(year/4))-(year/4)).eq.0) jdmax=366
       inewyr=0
       i=1
       do while (i.le.indEB1)
        if (ina1(i).eq.1) then				! flags start of the year
!        if ((ina1(i).eq.1).and.(idateMM.gt.11)) then	! flags start of the year
!         ina1(i)=ina1(i)+jdmax				! adds the last JD e.g. 365 or 366
         inewyr=1
        end if    

        if (ina1(i).gt.jdmax) then			! flags start of the year
         ina1(i) = ina1(i)-jdmax			! substracts the last JD e.g. 365 or 366
         print*,'check last day of the year -- jdmax: ',jdmax
        end if
 
        ina1r(i)  = dble(ina1(i))			! float
        ina9r(i)  = dble(ina9(i))
        ina10r(i) = dble(ina10(i))
        ina11r(i) = dble(ina11(i))
        ina12r(i) = dble(ina12(i))
        keyEB(i)  = dble(var2(i))+dble(100000.000*(ina1r(i)-ina1r(1)))  ! counter id to sort the file
        if (ina1(i).le.0.) keyEB(i)=990000.000+dble(var2(i))		! kick the strange values to the end
        if ((inewyr.eq.1.).and.(ina1r(i).lt.10)) keyEB(i)=keyEB(i)+(100000.00*370.)
        if (keyEB(i).eq.0.) keyEB(i)=990000.000+dble(var2(i))		! kick the strange values to the end
        keyEBk(i) = keyEB(i)
        i=i+1
       enddo
      
       print*,'indEB1: ', indEB1
       irecord=0
       i=1
       do while (i.le.indEB1-1)     !(i.le.n-1)
        if (keyEBk(i).gt.keyEBk(i+1)) irecord=1+irecord
        i=i+1
       enddo
       print*, 'Times the gmt time vector goes back and forth: ', irecord

       call sort2(indEB1,keyEB,var2) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,var3) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,var4) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,var5) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,var6) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,var7) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,var8) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,var13) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,var14) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,var15) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,var16) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,ina1r)
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,ina9r)
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,ina10r) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,ina11r) 
       call back2(indEB1,keyEB,keyEBk)
       call sort2(indEB1,keyEB,ina12r) 
       call back2(indEB1,keyEB,keyEBk)

       i=1
       do while (i.le.indEB1)
        ina1(i)=int(ina1r(i))
        ina9(i)=int(ina9r(i))
        ina10(i)=int(ina10r(i))
        ina11(i)=int(ina11r(i))
        ina12(i)=int(ina12r(i))
        i=i+1
       enddo

       i=1
       if (inewyr.eq.1) then
        do while (i.le.indEB1)
         if (ina1(i).gt.jdmax) then			!flags start of the year
          ina1(i)=ina1(i)-jdmax				!substracts the last JD e.g. 365 or 366
         end if
         i=i+1
        enddo
       endif

!!     201 file:
       open(45,file=fnameout(1:lfnout)//'.201', status='new')  
       write(45,713) enum, yymmdd
713    format('LBA Eddy',i1, ', ',a6,' slow data (1/2 hz),',
     &       ' scaled to Volts')
       write(45,*)'julday  gmt.s '
       write(45,*)
     & '   pcon   flsam   flcal   pcell   rl   tcell   tdetect   tpump'
       write(45,*)
     & '   tamb   rh.tdew   sid   cal   mux   cm'
       
       i=1
       do 812 while(i.le.indEB1)
        if (ina1(i).eq.0) goto 813		       	! skips jd==0
        write(45,814)ina1(i),var2(i),
     &             var3(i),var4(i),var5(i),var6(i),var7(i),var8(i),var13(i),var14(i),
     &             var15(i),var16(i),ina9(i),ina10(i),ina11(i),ina12(i)
814     format(i4,1x,f9.3,/,8f8.3,/,f9.4,1x,f9.4,1x,i3.3,1x,
     &       i2.2,1x,i3.3,1x,i1)
813     i=i+1
812     continue
!       end do
       close(45)
       
       print *, '  File EB*.201 is sorted'
      end if    ! date brackets to sort file
!---------------------------------------------------------------------------------------------------

      print *, '  End of program LSPLIT EDDY'
      end       
!     end of program


!===================================================================================================
!                                  subroutines and functions
!===================================================================================================

!===================================================================================================
!                                      subroutine vecdiv
! enable the approximation for vectorized division (nrc)
!===================================================================================================
      subroutine vecdiv(a,b,ab,n)
      real*4 a(n),b(n),ab(n)
      do i=1,n
      ab(i)=a(i)/b(i)
      if( abs(ab(i)) .gt. 999.) ab(i)=-99.999
      end do
      return
      end
      function idotprod(a,b,n)
      integer a(n),b(n)
      idotprod=0.
      do i=1,n
      idotprod=idotprod+a(i)*b(i)
      end do
      return
      end
!===================================================================================================
!                                   end of subroutine vecdiv
!===================================================================================================



!===================================================================================================
!                                     subroutine istat
!===================================================================================================
      subroutine istat(is1,isol)
!     extracts status bits from integer is1, puts into integer array isol
      integer is1,isol(12)
      do i=1,12
       isol(i)=0
       if (btest(is1,i-1)) isol(i)=1
      enddo
      return
      end
!===================================================================================================
!                                    end of subroutine istat
!===================================================================================================



!===================================================================================================
!                                   subroutine infile_spl
!===================================================================================================
      subroutine infile_spl(ymd,types,names,num)
      character*6 ymd, types
      character*12 names(6)
      character*12 fchars
      integer fnum
      fchars='EAEBPFGDSNNT'
      num= lentrim0(types,len(types) )
!     need to add error checking (make sure no duplicates, etc)
      do i=1,num
      fnum = ichar(types(i:i))-48
      names(i)=ymd//fchars( ((fnum-1)*2+1):((fnum-1)*2+2) )//'.dat'
      enddo
      return
      end
!===================================================================================================
!                                   end of subroutine infile_spl
!===================================================================================================



!===================================================================================================
!                                        function checkrec
!===================================================================================================

      function checkrec(jday, tsec, prevday, tsecprev, tallow, jdmax)
      logical newyear  
      integer prevday, jdmax
      real*4 checkrec, error, tsec, tsecprev, tallow
      real*8 tnow,tprev, tdiff
      error=0.
      tnow = jday + tsec/86400.                    ! times in decimal days
      tprev = prevday + tsecprev/86400.
      newyear = .false.  
      if (jday .eq. 1 .and. prevday.ge.jdmax) newyear = .true.
      tdiff = 86400.*(tnow-tprev)                  ! elapsed time (in seconds)
!     criteria for error:
      if((jday.lt.0).or.(jday.gt.jdmax+1)) err=1.  ! day out of range
      if((tsec.lt.0.).or.(tsec.gt.86410.)) err=1.  ! time-of-day out of range
      checkrec = error
      end

!===================================================================================================
!                                   end of function checkrec
!===================================================================================================


!===================================================================================================
!                                     subroutine sort2                                      !
!===================================================================================================
      subroutine sort2(n,arr,brr)
      integer n,m,nstack
      real arr(n),brr(n)
      parameter (m=7,nstack=50)
      !Sorts an array arr(1:n) into ascending order using quicksort, while making the corresponding
      !rearrangement of the array brr(1:n).
      integer i,ir,j,jstack,k,l,istack(nstack)
      real a,b,temp
      jstack=0
      l=1
      ir=n

1     if(ir-l.lt.m)then 			!Insertion sort when subarray small enough.
      do j=l+1,ir	!12
      a=arr(j)
      b=brr(j)
      do i=j-1,l,-1	!11 
      if(arr(i).le.a)goto 2
      arr(i+1)=arr(i)
      brr(i+1)=brr(i)
      enddo 		!11
      i=l-1
2     arr(i+1)=a
      brr(i+1)=b
      enddo 		!12
      if(jstack.eq.0)return
      ir=istack(jstack) 			!Pop stack and begin a new round of partitioning.
      l=istack(jstack-1)
      jstack=jstack-2
      else
      k=(l+ir)/2 
!     Choose median of left, center and right elements as partitioning
!     element a. Also rearrange so that a(l) ≤ a(l+1) ≤ a(ir).
      temp=arr(k)
      arr(k)=arr(l+1)
      arr(l+1)=temp
      temp=brr(k)
      brr(k)=brr(l+1)
      brr(l+1)=temp
      if(arr(l).gt.arr(ir))then
      temp=arr(l)
      arr(l)=arr(ir)
      arr(ir)=temp
      temp=brr(l)
      brr(l)=brr(ir)
      brr(ir)=temp
      endif
      if(arr(l+1).gt.arr(ir))then
      temp=arr(l+1)
      arr(l+1)=arr(ir)
      arr(ir)=temp
      temp=brr(l+1)
      brr(l+1)=brr(ir)
      brr(ir)=temp
      endif
      if(arr(l).gt.arr(l+1))then
      temp=arr(l)
      arr(l)=arr(l+1)
      arr(l+1)=temp
      temp=brr(l)
      brr(l)=brr(l+1)
      brr(l+1)=temp
      endif
      i=l+1 				!Initialize pointers for partitioning.
      j=ir
      a=arr(l+1) 			!Partitioning element.
      b=brr(l+1)
3     continue 				!Beginning of innermost loop.
      i=i+1 				!Scan up to find element > a.
      if(arr(i).lt.a)goto 3
4     continue
      j=j-1 				!Scan down to find element < a.
      if(arr(j).gt.a)goto 4
      if(j.lt.i)goto 5 			!Pointers crossed. Exit with partitioning complete.
      temp=arr(i) 			!Exchange elements of both arrays.
      arr(i)=arr(j)
      arr(j)=temp
      temp=brr(i)
      brr(i)=brr(j)
      brr(j)=temp
      goto 3 				!End of innermost loop.
5     arr(l+1)=arr(j) 			!Insert partitioning element in both arrays.
      arr(j)=a
      brr(l+1)=brr(j)
      brr(j)=b
      jstack=jstack+2
!     Push pointers to larger subarray on stack, process smaller subarray immediately.
      if(jstack.gt.nstack)pause 'nstack too small in sort2'
      if(ir-i+1.ge.j-l)then
      istack(jstack)=ir
      istack(jstack-1)=i
      ir=j-1
      else
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
      endif
      endif
      goto 1

      end

!===================================================================================================
!                                   end of subroutine sort2
!===================================================================================================


!===================================================================================================
!                                   subroutine back2
!===================================================================================================
      subroutine back2(n,arr,crr)
      integer n, k
      real crr(n),arr(n)
      k=1
      do while (k.le.n)
      arr(k)=crr(k)
      k=k+1
      end do
      end
!===================================================================================================
!                                   end of subroutine back2
!===================================================================================================

!===================================================================================================
!                                   function valmax
!===================================================================================================
      function valmax(n,arr)
      integer n, valmax, k
      integer arr(n)
      valmax=arr(1)
      k=1
      do while (k.le.n)
      if (valmax.gt.arr(k)) valmax=arr(k)
      k=k+1
      end do
      end
!===================================================================================================
!                                   end of subroutine back2
!===================================================================================================


      include 'lsubs.for'


