!===========================================================================================
!                                  LBA code:  lsplit7.for
!===========================================================================================
!
!     Takes CR10-X binary data, converts them into ascii form, and then separates into 
!     readings from individual sensors, applies scale factors to obtain original sensor 
!     voltages from CR10-logged voltage-divided mV.
!     Invoke:  'lsplit -date yymmdd -files 1234'
!
!     Originally created by Scott Saleska       (started jul 2000)
!     Adapted from b1split Boreas code written by Song-Miao Fan
!
!     Input files:            .\yymmdd\yymmddxx.dat
!     Output files:           .\yymmdd\split\xxyymmdd.200, 201, 202, 203
!     Error/log file to:      .\yymmdd\process\xxyymmdd.err
!
!     There are six possible kinds of CR10 input data from files: yymmddxx.dat, 
!     where xx depends on site_id of system from which data comes:
!     xx:            EA       EB         PF        GD        sn       nt
!     site_id:  1=eddy1  2=eddy2  3=profile  4=ground  5=sonic2  6=natal
!
!     ======================================================================================
!     Input CR10-X data fields (by site, and data-line id #):
!     ======================================================================================
!
!     Site = 1 or 2 (Eddy System A (EA) or Eddy System B (EB)):
!     ---------------------------------------------------------
!
!     id=100 (initialization loop header) - 6 fields:
!     100, yr, day, hhmm, ss, site_id
!
!     id=200 (8 Hz) --  1+7 = 8 fields:
!     200, co2, h2o, u, v, w, c, status {, tfw}
!
!     id=201 (0.5 Hz) - 1+3+2+10 = 16 fields: 
!     201, day, hhmm, ss, status1, status2, 
!     pcon, flsam, flcal, pcell, rl, tcell,
!     tdetect, tpump, tdew.rh, tamb
!
!     id=202 (0.5 Hz) - 1+3+6 = 10 fields:
!     202, day, hhmm, ss,         
!     spare, spare, spare, spare
!
!     id=203 (0.5 hz) - 1+3+10 = 14 fields:
!     203, day, hhmm, ss,         
!     cmstatus, spare, spare, spare, spare, spare,
!     spare, spare, spare, spare
!     
!     Notes:
!     a) 201: status1 is actually a spare
!     b) 202: 4 spares are mux inputs where, e.g., sonic analog outputs could go
!     c) 203: spares:  7 are digital status signals, last 2 are pulse inputs
!        id=200 (4 hz) --  1+7 fields
!     d) in output, put tamb, tdew.rh in 201 line, put cmstatus in status bits
!     e) line 201 was set to this (cmstatus at end, now in 203):
!        id=201 (0.5 hz) - 1+3+2+11 = 17 fields: 
!        201, day, hhmm, ss, status1, status2, 
!        pcon, flsam, flcal, pcell, rl, tcell,
!        tdetect, tpump, tdew.rh, tamb, cmstatus
!
!     Site = 3 (PF = Profile System):
!     -------------------------------
!     id =      100 (initialization loop header) - 1+5 = 6 fields
!               100, yr, day, hhmm, ss, site_id
!
!     id =      201 (0.5 hz) - 1+5+12 = 18 fields
!                201, day, hhmm, ss, status1, status2,
!               pcon, flsam, flcal, pcell, rl, tcell,
!               tcinlet, tpump, co2, h2o, spare, spare
!
!     id =      202 (0.5 hz) - 1+3+20 fields
!                202, day, hhmm, ss, wdi, spare, netrad, netrad_cal
!               par1_up, par1_down, par2_up, par_cal, spare, spare
!               spare, spare, tair_1, tair_2, tair_3, tair_4,
!               tair_5, tair_6, tair_7, tair_8
!
!     id =      203 (0.5 hz) - 1+3+10=14 fields
!                203, day, hhmm, ss,
!               ws1, ws2, ws3, ws4, spare, spare
!               spare, spare, rain, spare
!
!     site = 4 (GD = ground):
!     -----------------------
!     id =      100 (initialization loop header) - 6 fields:
!               100, yr, day, hhmm, ss, site_id
!     id =      201 (0.5 hz) - 1+5+16=22 fields
!               201, day, hhmm, ss, status1, status2,    
!               p_cohs, p_col, p_ss, p_hs, p_ms, p_ls, p_zag, tpump_co, co, 
!               tcell_co, pcon_co, pamb, tcatalyst, par(1-8), 
!               tsoil(1-8) tsoil(9-16)
!
!     site = 5 (sn = sonic):
!     ----------------------
!     id =      100 (initialization loop header) - 6 fields:
!               100, yr, day, hhmm, ss, site_id
!     id =      200 (4 hz) --  1+5 fields:
!               200,  u, v, w, t, status
!
!     all cr10 units:
!     ---------------
!     id =      300 (low power output) - 4 fields: 
!               300, day, hhmm, ss, ??
!-------------------------------------------------------------------------------------------

      program lsplit

      include 'lparams.for'          ! To be compatible with 'linitial7.for'

!     initial declarations:=================================================================
      character*1 cr10(4096),ctmp
      character*96 basename, fnamerr ! Name of output file
      integer nifiles, infile 
      integer id,hhmm,status1,status2, icsat 
      integer*2 icr10(2048),itmp
      integer isol1(12),isol2(12)    ! 12-bit status words
      integer is200(5)               ! 24 bits is max needed for is201 (?)
      integer f200, f201, f202, f203 ! format line labels 
      logical lskip, output          ! whether to skip, ascii format, output data
      real    startnum               ! start.numid: start= number of label 100 to start with
      integer nid,nidstart,nidend    ! number of labels, num to start, num to end with
      integer jrec                   ! index to block number (block = 2048 2-byte words)
      integer norec                  ! number of records (= # discrete cr10 lines)
                                     ! (# of records within each blockrecords cross blocks)

!...........................................................................................
!     *** note
!      that there are 3 possible errors:  
!      1) bad-data error:
!        this word is uninterpretable according to cr10 binary encoding scheme.
!      2) invalid id# error:
!        expected line id# (100,200,201,etc) but found something else
!      3) time-mismatch error:
!        this time comes before previous time 
!...........................................................................................

      integer noerr                   ! number of bad-data errors within this block
      real*4  toterr                  ! total number of bad-data errors in this file
      integer nblkerr                 ! number of blocks with >=1024 bad-data errors
      integer niderr, totiderr        ! number of invalid id errors
      integer ntierr, tottierr        ! number of time mismatch errors
      integer eof                     ! check for end-of-file for direct access file 
        
!...........................................................................................
!     *** Note
!     Actual data array needs (not incl. id/date/time/status data):
!                                                              maximum
!      *.200:  e*  7 fields                                     7
!      *.201:  e* 10 fields; PF 12 fields; GD 13 fields        13
!      *.202:  e*  4 fields; PF 20 fields; GD 24 fields        24
!      *.203:  e* 10 fields; PF 10 fields, GD  8 fields        10
!...........................................................................................

!     Parameters to define max dimensions for input/output arrays:
      parameter (nd200=8, nd201=13, nd202=24, nd203=10)

!     Input ascii array for record block. 
!     Note that the size of "acr10" was increased by 10 to allow slack for errors:
      real*4  acr10(2058), rcr10(2048)                     
      
!     Input ascii arrays (not incl. id/date/time/status data):
      real*4  a200(nd200), a201(nd201), a202(nd202), a203(nd203)

!     Output voltage arrays (not incl. id/date/time/status data):
      real*4  v200(nd200), v201(nd201), v202(nd202), v203(nd203)

!     Scale-factor arrays (v20* = a20*/sf20*):
      real*4 sf200(nd200),sf201(nd201),sf202(nd202),sf203(nd203)
      real*4 sf200SN(nd200)        
                              
!     Arrays co-located in memory:
      equivalence ( cr10(1) , icr10(1) )

!...........................................................................................
!     *** note
!     About scale factors sf200, sf201, sf202, sf203:
!     The scale factors divide the input signal, giving signal output by sensor (in v).
!     Account for gain factor in data system, and cr10x recording in mV. Then, a signal 
!     with gain=1 implies in an scale factor = 1000 to convert to v).
!...........................................................................................

!...........................................................................................
!     *** note
!      remember to change scale factors if:
!     (1) output amplification/reduction factors are changed on board,
!          or
!     (2) data sequence in cr10 data fields is changed.
!...........................................................................................

!     scale factors for eddy system: =======================================================

!     gain x 1000 to convert to volts (not valid for csat3):
      data sf200/2*500., 6*1./

!     no gain for the csat3:
      data sf200SN/8*1./

      real*8 ttsn
      real*4 sf201E(12)     
!          pcon flsam flcal pcell RL tcell tdetect tpump tdew.rh tamb              ! Eddy
!          pcon flsam flcal pcell RL tcell tdetect tpump CO2  H2O  spare spare  ! ProFile
      data sf201E/6*500., 2*1000., 2*500., 2*1000./ 
!                                     ^ tdew.rh/tamb (ED) and CO2/H2O (PF) all have gain 1/2
      real*4 sf203E(1)                                ! 203 data  (only cm status right now)
      data sf203E/1000./

!...........................................................................................
!     *** note                                    
!     201 = 0.5 hz (eddy & profile combined)
!     pcon flsam flcal pcell rl tcell tdetect tpump tdew.rh tamb             --> Eddy
!     pcon flsam flcal pcell rl tcell tdetect tpump co2  h2o  spare spare    --> PF has 2 extra
!     tdew.rh, tamb (ed) and co2, h2o (PF) all have gain 0.5.
!...........................................................................................

!     scale factor to profile: =============================================================

      real*4 sf202PF(20)                  ! 202 line

!     profile:      wdi spare  netrad(x2)  par(x4)  spare(x3)  tsoil(x9)
      data sf202PF/   2*500.,    2*1000.,  4*1000.,  3*1000.,   9*1000./

      real*4 sf203PF(10)                  ! 203 line
                                          ! profile:  ws(x4), spare(x4), rain, spare       
      data sf203PF/8*1000., 2*1./         ! don't scale rain-tips by 1000

!...........................................................................................
!     *** note
!     netrad was left unscaled, because gain of 20 gives good range (-.5 to 2.0 v)
!     31-may-2002: error! netrad scale factor should be 10k ( 1000 mv/v * 10(gain) )
!     (now that it is here, leave it, adjust with conversion factor in lconv.inp file)
!     with sf=1000, will have values of 0.25 to 1 volt
!     par is left un-scaled, b/c gain = 200 gives range 0 to 1.933 v.
!     in 31-may-2002: par is not right either, but it is working ok with existing 
!     conversion factors (should be 100k b/c gain was 100; i think did this deliberately 
!     so the conversion factors would not get too big. ("i" = srs?)
!     tsoil is done by 30k thermistor, 30.1k pull-up (1 for now).
!...........................................................................................


!     scale factors to ground:==============================================================
      real*4 sf201GD(7)                ! 201 line
                                       ! ground: co pcon fl stat tcell ttrap tpump
      data sf201GD/4*500., 3*1000./    ! fl, stat still tbd
      real*4 sf202GD(24)               ! 202 line
                                       ! ground:  tsoil(x16)  par(x8)
      data sf202GD/16*1000.,8*1000./
      real*4 sf203GD(8)                ! 203 data
                                       ! ground:  p_ss p_hs p_ms p_ls pcohs pcols p_zag pamb  
      data sf203GD/ 7*500., 250./

!     scale factors to sonic:===============================================================
      real*4 sf201SN(3)                ! sonic
      data sf201SN/3*1./               ! time stamp only  
      data is200/5*0/                  ! max # status bits
      data is201/24*0/                 ! max # 201 status bits

!...........................................................................................
!     *** note
!     input data are in milivolts, sf*1000 makes output in volts.
!     get user input (from command line)
!...........................................................................................

      jrstart = 1          ! first record reading binary data, default
      jrend   = 100000     ! last  record reading binary data, default
      startid = 0.         ! option: id# to key off of output only a certain # of the 
                           ! specified id#. 
                           ! e.g. startid=100 (default=0:  output all id's)
      ascii = .false.

!     set initial parameters, get command-line, etc.: ======================================

      lname = 'lsplit'
      ftypes= '1234'          ! default file types to process: 
                              ! 1=eddy1, 2=eddy2, 3=profile, 4=ground, 5=sonic
      include 'linitial7.for' ! common initialization code 
                              ! sets directory structure 
                              ! (root paths,dirsplit,dirprocin,dirprocout)
                              ! retrieves command line, file names for processing 
                              ! [ to ifilen(1..nifiles) ] displays initial program message.

      call infile_spl(yymmdd, ftypes, ifilen, nifiles) ! gets 'nfiles' input file names 
                                                       ! (ifilen) from: yymmdd and ftypes.
                                                       ! filenames are:  
                                                       ! 1 = yymmddea; 2 = yymmddeb; 
                                                       ! 3 = yymmddpf; 4 = yymmddgd;
                                                       ! 5 = yymmddsn; 6 = yymmddnt.
      print *
      print 9, nifiles 
9     format('  LSPLIT: processing ', i1, ' input files:') 
      do i=1,nifiles
      print *, '  ', ifilen(i)
      end do 

!     opens input and output files for id=200,201,202, and respectively:
!     set up for reading raw data from drive
!     write the initial split files to local drive (run from d:)
!     write the final outputs to c:, the smaller hard disk  200 mb
!
!     yymmddxx.dat = raw data from cr10 in binary form
!     xxyymmdd.200=scaled data for id=200
!     xxyymmdd.201=scaled data for id=201
!     xxyymmdd.202=scaled data for id=202
!     xxyymmdd.203=scaled data for id=203
!
!     ======================================================================================
!                              main program loop starts here
!     ======================================================================================

      do infile=1, nifiles        ! loops over each input file

!     generic initializations===============================================================
      startid = startidin         ! id# to start with (0 if none specified)
      nid=0                       ! init number of startid's found in file
      nidstart=1
      nidend=999
      if(startid.eq.0) then       ! does not allow output if we have to find startid first
      output= .true. 
      else 
      output= .false.
      if(startnum.ne.0) then 
      nidstart = int(startnum)                     ! number of start id's to go by 
                                                   ! before reading.
      nidend = int(10*(startnum-float(nidstart)))  ! reads up to 9 labels.
      print *, ' # before decimal = ',nidstart
      print *, ' # after decimal  = ',nidend
      nidend = nidstart + nidend 
      endif
      endif

!     gets this loop's filename: raw binary data input (*.dat files)========================

      dtype= ifilen(infile)(7:8)                         ! dtype is data type string
      if (dtype .eq. 'EA') idtype = 1 
      if (dtype .eq. 'EB') idtype = 2 
      if (dtype .eq. 'PF') idtype = 3 
      if (dtype .eq. 'GD') idtype = 4 
      if (dtype .eq. 'SN') idtype = 5 
      if (dtype .eq. 'NT') idtype = 6 
      basename = yymmdd//dtype//'.dat'
      lrt = lentrim0(rootsplit, len(rootsplit))
      fname = rootsplit(1:lrt)//slash//yymmdd//slash//basename
      lfn = lentrim0(fname, len(fname))
      read (yymmdd,'(i10)') idate			! nrc 170207 add to hard wire problems with loggers 2012-2013
 
!     data output: =========================================================================
      fnameout = dtype//yymmdd
      call getfullpath(rootsplit, yymmdd, dirsplit, fnameout, lfnout)
        
!     error/log file: ======================================================================
      fnamerr = dtype//yymmdd
      call getfullpath(rootprocout, yymmdd, dirprocout, fnamerr, lfnerr)

!     set-up input files: ==================================================================
      print *,'  Opening input file "',fname(1:lfn)//'"' ,
     &        '  for output to      "',fnameout(1:lfnout)//'.2*" files'
!      print *,'  for output to      "',fnameout(1:lfnout)//'.2*" files'

      if(.not. ascii) then    ! the default is binary format.
      open(11,file=fname(1:lfn), status='old',
     &        form='unformatted',access='direct',recl=4096)
      else   
!     option to do an ascii text file (doesn't really work):            (???)
      open(11,file=fname(1:lfn),status='old')
      endif

!     set-up output files: =================================================================
      open(99,file=fnamerr(1:lfnerr)//'.err',status='unknown')        ! error/log file

!     inserted 1-26-02:.....................................................................
      write(99,55) platform, lname, iyear, imonth, iday,ihrs,imin,isec
55    format(a6,' ',a7,': on date (yy-mm-dd) ',i4,'-',i2.2,'-',i2.2,
     &       ' (',i2.2,':',i2.2,':',i2.2,'), process:' )
      write(99,*) '   ', fname(1:lfn)  
      if(verbose) write(99,*) '   report mode = verbose'
! *** end of insert..........................................................................

      if(ascii) write(99,*) '   (to be read as ascii file)'
      if (dtype .ne. 'sn') then  
      open(21,file=fnameout(1:lfnout)//'.201',status='unknown')       ! first output file
      endif   
!     note that all files have a 201 line - except sn sonic files
      if(startid.ne.0) then
      write(99,'(a25,i4,a10,f5.0)')' start with occurence # ',nidstart,
     &                             ' of label=',startid
      write(*, '(a25,i4,a10,f5.0)')' start with occurence # ',nidstart,
     &                             ' of label=',startid
      if(nidend.lt.999) then
      write(99,'(a23,i4)')' end after occurence #',nidend
      write(*, '(a23,i4)')' end after occurence #',nidend
      endif
      endif

!     input file-specific initializations: =================================================
      select case(idtype)

!     ======================================================================================
!      case 1 and 2: eddy data
!     ======================================================================================

      case(1:2)                                
       enum=1.                    ! if 'EA', enum=1
       if(dtype.eq.'EB') enum=2.    
       idsys = enum               ! set system id
       nf200 = 7                  ! number of input fields   (not incl line id#)
                                  ! nf200=7 if no date/time  (default)
       nf201 = 15                 ! includes 5: day,hhmm,ss, status1, status2
       nf202 = 7                  ! includes 3: day,hhmm,ss  (remainder all spares)
       nf203 = 13                 ! includes 3:       "      (remainder: cmstatus + spares)
       ns200 = 5                  ! number of status bits to output on 200 line
       ns201 = 9                  !                                 on 201 line
       do i=1,nf201-5              ! scale factor array (don't include day/time/status).
       sf201(i) = sf201E(i)
       print*, 'scale factor sf201E(i): ', sf201E(i)
      enddo
      sf203(1) = sf203E(1)        ! only 1 var (cmstatus) on 203 line
      do i=1,nf203-3              ! debug test: use profile scale-factors
       sf203(i) = sf203PF(i)
      enddo
      assign 201 to f200          ! format label numbers for writing data to files
      assign 211 to f201
      open(20,file=fnameout(1:lfnout)//'.200',status='unknown')
 
!     write headers for eddy output files: =================================================
!     200 file:
      write(20,01) enum,yymmdd
01    format('  LBA Eddy',i1, ', ',a6,' fast channels (8 Hz), 200 data')
      write(20,*)'julday   gmt.s   co2.v   h2o.v   u.m.s   v.m.s',
     &           '   w.m.s   t.c   stat   sid   cal'

!     format for data output (used in output section to match headers here):
201   format(i4,1x,f9.3,1x,2(f7.4,1x),3(f7.3,1x),f7.3,1x,i6,1x,3i1,1x,
     &       2i1) ! 200 line 
!     200 line output format #1 --> 201 (format #2 is backup sonic output)

!     201 file:
      write(21,11) enum,yymmdd
11    format('LBA Eddy',i1, ', ',a6,' slow data (1/2 Hz), scaled to V')
      write(21,*)'julday   gmt.s   '
      write(21,*)'pcon   flsam   flcal   pcell   rl   tcell   ',
     &           'tdetect   tpump'
      write(21,*)'tdew.rh   tamb   sid   cal   mux   cm'

!     format for outputting data corresponding to these headers:     
211   format(i4,1x,f9.3,/,8f8.3,/,f9.4,1x,f9.4,
     &       1x,3i1,1x,2i1,1x,3i1,1x,i1)
!      modified 030806: inserted space between tdew.rh and tamb so don't run together if error

!     ======================================================================================
!      case 3: profile data
!     ======================================================================================
      case(3)
       idsys=3    
       nf200=7                     ! not used for profile
       nf201=17                    ! not incl id#, including 5:  day, hhmm, ss, status1,2
       nf202=23                    ! not incl id#, including 3:  day, hhmm, ss
       nf203=13                    ! not incl id#, including 3:        "
       ns201=24                    ! number status bits to output for 201 file
       do i=1,nf201-5              ! scale factor array (not incl 5 date/time/status vars)
        sf201(i) = sf201E(i)
       enddo
       do i=1,nf202-3              ! scale factor array (not incl 3 date/time vars)
        sf202(i) = sf202PF(i)
       enddo
       do i=1,nf203-3              ! scale factor array (not incl 3 date/time vars)
        sf203(i) = sf203PF(i)
       enddo
       assign 212 to f201          ! format label number:  212 -> 201.2 = 201 label #2
       assign 222 to f202          ! 202 label #2 (after eddy)
       assign 232 to f203          ! 203 label #2
       open(22,file=fnameout(1:lfnout)//'.202',status='unknown') 
       open(23,file=fnameout(1:lfnout)//'.203',status='unknown')

!     write headers to output files: =======================================================
!     201 file:
       write(21,12) yymmdd
12     format('LBA Profile 201 data, ',a6,': slow channels (1/2 Hz),',
     &       ' scaled to V')
       write(21,*)'julday   gmt.s   '
       write(21,*)'pcon   flsam   flcal   pcell   rl   tcell   ',
     &           'tdetect   tpump'
       write(21,*)'co2   h2o   spare   spare   sid   cal   mux   cc   ',
     &           'csol   1234   5678   c'

!     format for outputting data corresponding to these headers:     
212    format(i4,1x,f9.3,/,8f8.3,/,4(f9.4,1x),2(3i1,1x,2i1,1x),5i1,1x,
     &       2(4i1,1x),i1)     ! profile 201

!     202 file:
       write(22,22) yymmdd
22     format('LBA Profile 202 data, ',a6, ': slow channels (1/2 Hz),',
     &       'scaled to V')
       write(22,*)'julday   gmt.s   '
       write(22,*)'wdi   spare   netrad   nrad_cal   par1_up   ',
     &           'par1_dwn   '
       write(22,*)'par2_up   par_cal   totrad.umol   difrad.umol   spare   ',
     &           'tair_cal   '
       write(22,*)'tair_1   tair_2   tair_3   tair_4   tair_5 ',
     &           'tair_6   tair_7   tair_8   '
!     format for outputting data corresponding to these headers:     
222    format(i4,1x,f9.3,/,6(f8.3,1x)/,6(f8.3,1x),/,8(f8.3,1x) )               ! profile 202
!     203 file:
       write(23,32) yymmdd
32     format('LBA Profile 203 data, ',a6,
     &       ': slow channels (1/2 Hz), scaled to V')
       write(23,*)'julday   gmt.s   '
       write(23,*)'ws1   ws2   ws3   ws4   spare   spare   spare   ',
     &           'spare   rain   spare'

!     format for outputting 203 data corresponding to these headers:     
232    format( i4,1x,f9.3,/,10f7.3 )                                           ! profile 203

!     ======================================================================================
!      case 4: ground data
!     ======================================================================================
      case(4)  
       idsys=4
       nf201=12
       nf202=27
       nf203=11
       nf300=4                       ! power loss header
       ns201=18
       do i=1,nf201-5                ! scale factor array (not incl 5 date/time/status vars)
        sf201(i)=sf201GD(i)
       enddo
       do i=1,nf202-3                ! scale factor array
        sf202(i)=sf202GD(i)
       enddo
       do i=1,nf203-3                ! scale factor array
        sf203(i)=sf203GD(i)
       enddo
       assign 213 to f201            ! format label  (213-> 201.3, 201 label #3)
       assign 223 to f202            ! 222 -> 202label #3
       assign 233 to f203            ! 232 -> 203 label #3
       open(22,file=fnameout(1:lfnout)//'.202',status='unknown') 
       open(23,file=fnameout(1:lfnout)//'.203',status='unknown')

!     write headers to output files: =======================================================
!     201 file:
       write(21,13) yymmdd
13     format('LBA Ground 201 data, ',a6, ': slow channels (1/2 Hz),',
     &       ' scaled to V')
       write(21,*)'julday   gmt.s   '
       write(21,*)'co   pcon_co   fl_co   stat_co   tcell_co   ttrap   ',
     &           'tpmp_d   sid   mux   gp21   cdpz   lhzi'
!     format for outputting ground 201 data corresponding to these headers:     
213    format(i4,1x,f9.3,/,7f8.3,1x,2(3i1,1x),3(4i1,1x))       ! ground 201
!     202 file:
       write(22,23) yymmdd
23     format('LBA Ground 202 data, ',a6,': 5-min avg, scaled to Volts')
       write(22,*) 'julday   gmt.s   '
       write(22,*) 'ts1   ts2   ts3   ts4   ts5   ts6   ts7   ts8   '
       write(22,*) 'ts9   ts10   ts11   ts12   ts13   ts14   ts15   ts16'
       write(22,*) '   par1   par2   par3   par4   par5   par6',
     &            '   par7   par8   '
!     format for outputting ground 202 data corresponding to these headers:     
223    format( i4,1x,f8.0,3(/,8f8.3) )                         ! ground 202 (5-min avgs)
!     203 file:
       write(23,33) yymmdd
33     format('LBA Ground 203 data, ',a6,': 5-min avg, scaled to V')
       write(23,*) 'julday gmt.s   p_ss   p_hs   p_ms   p_ls   ',
     &            'p_cohs   p_cols   p_zag   pamb'
!     format for outputting ground 203 data corresponding to these headers:     
233    format(i4,1x,f8.0, 8f8.3)                               ! ground 203 (5-min avgs)

!     ======================================================================================
!      case 5: backup sonic data
!     ======================================================================================
      case(5)                                           ! sonic data : starting sept '2006  
       idsys=5
       nf200=5  
       ns200=5
       do i=1,nf200 
        sf200(i) = sf200SN(i)  
       enddo
       nf201=3
       assign 214 to f200 !               
       open(20,file=fnameout(1:lfnout)//'.200',status='unknown') 
 
!     add: write headers to output files
!     200 file:
       write(20,14) yymmdd
14     format('LBA Sonic, ',a6,' fast channels (8 hz), 200 data')
       write(20,*)'julday   gmt.s   u.m.s   v.m.s   w.m.s   ',
     &           't.c   stat'
!     format for data output (used in output section to match headers here):
214    format(i4,1x,f9.3,1x,3(f7.3,1x),f7.3,1x,i6) ! 200 line for sonic
!     format for outputting sonic 200 data corresponding to headers:     
202    format(i4,1x,f8.2,1x,5(f8.4,1x),3i1,1x,1x)  ! 200 line for backup sonic

!     ======================================================================================
!      case 6: natal data
!     ======================================================================================
      case(6)
       idsys=6
       nf201=15     ! fill-in more when we know what natal system is like

!     add: write headers to output files
      end select

!     ======================================================================================
!     reads yymmdd.dat & convert cr10 binary data into ascii form
!     ======================================================================================
      julday     = -999     ! relevant for eddy data file:
      nfast      = -1       ! number of fast data lines read (-1 if no 201 lines read yet)
      cmlast     = 10.0     ! last chilled-mirror status (10 v if none read yet) from 203 
                            ! input line but output as status bit on 201 line)
      noerr      = 0        ! number of 2-byte words with errors in current input 
                            ! block (0-2048)
      noblkerr   = 0        ! number of blocks with significant number of errors 
                            ! (i.e., noerr >= 1024)
      toterr     = 0.0      ! total number of words with errors in whole file
      niderr     = 0        ! number invalid id# errors in this output block
      totiderr   = 0        ! total number of invalid id errors in this file
      ntierr     = 0        ! number of time-mismatch errors in this output block
      tottierr   = 0        ! total time-mismatch errors in this file
      norec      = 0        ! number of records (cr10-data lines) read
      n201       = 0        ! number of records (201-lines only) read
      nbadinrow  = 0        ! number of bad values or records in a row
      ngoodinrow = 0        ! number of good 201-line records in a row
      tsave      = 0.
      k          = 0
      itmp       = 0
      n12        = 0  
      jrec       = jrstart-1

!     gets the next 4096-byte record =======================================================
100   jrec = jrec+1
      noerr = 0                          ! re-zero number of errors in current record
      if (jrec .gt. jrend) goto 980      ! normal exit of reading this file
      read (11,rec=jrec,iostat=eof) cr10 ! cr10: 4096 1-byte chars
      if (eof .lt. 0) goto 990           ! cr10 co-located with icr10 (2048 2-byte integers)
      do i=1,4095,2                      ! flips bytes within each 2-byte word
       ctmp      = cr10(i+1)
       cr10(i+1) = cr10(i)
       cr10(i)   = ctmp
      enddo

!     decodes binary cr10 record (in icr10) into ascii (in acr10) ==========================
!     j is index to icr10 array (the data input block)
!     k is index to acr10 ( the corresponding decoded output block)
!     in principle k more or less = j, except:
!          * 4-byte value increments j twice, k once
!          * when bad-data error, throw it out (increments j, not k)
!          * offset between 2048-array and data-lines
      j = 0
110   j = j+1
      if (j .gt. 2048) then            ! at end of this record
       if(verbose) then 
        write(99,111) jrec, icr10(1), icr10(2048), noerr
111     format('  Input block ',i5,': icr10(1),(2048)= ',i6,
     &       ', ',i6,'; bad-data errors= ',i4)
       elseif(noerr.gt.0) then
        write(99,112) jrec, noerr
112     format('  Input block ',i5,': ',i5,' of 2048 words have bad data',
     &       ' error')
       endif
       goto 100                        ! go to get next record
      endif
      if(ascii) then                   ! if ascii file, no decoding necessary
       k=k+1
       acr10(k) = rcr10(j)
       if (k .eq. 2048) goto 160
       goto 110
      endif
      if ( iand(icr10(j),  7168) .ne.   7168 ) goto 120     ! a 2-byte word of data 
      if ( iand(icr10(j),  -512) .eq.  -1024 ) goto 130     ! start of array
      if ( iand(icr10(j), 15360) .eq.   7168 ) goto 140     ! first 2 bytes of 4-byte value
      if ( iand(icr10(j), -1024) .eq.  15360 ) goto 150     ! last 2 bytes of 4-byte value

!     error in this 2-byte word (i.e. none of the above) ===================================
!     added in 23-jan-2002:
      noerr = noerr+1
      toterr = toterr+1.0
      if ((nblkerr.eq.0).or.(noerr.eq.1).or.(noerr.eq.1024)) then

!     if nblkerr (first block with errors) or error # 1, 1024 then output error message
       write(99,114)  noerr,jrec, j, icr10(j) 
114    format('  Bad data error #',i5,' in block jrec=',i5,
     &  ':  icr10(j =',i5,') = ',i6)
      endif
      if (noerr .eq. 1024) nblkerr = nblkerr + 1      ! number of records with 1024+ errors
      goto 110

!     2-bytes data:
120   k=k+1
      if ( iand(icr10(j), -32768) .ne. 0) then
       sf = -1.
      else
       sf = 1.
      end if
      if ( iand(icr10(j),  16384) .ne. 0) sf = sf * 0.01
      if ( iand(icr10(j),   8192) .ne. 0) sf = sf * 0.1
      acr10(k) = float( iand(icr10(j), 8191) ) * sf
      if (k.ge.2048) goto 160
      goto 110

!     start of array:
130   k=k+1
      acr10(k) = float( iand(icr10(j), 511) )
      if (k.ge.2048) goto 160
      goto 110

!     first 2 bytes of 4-byte value:
140   continue                      ! k=k+1
      itmp = icr10(j)
      n12  = n12 + 1                 ! n12 is number of sequential first halves 
                                    ! without a second half  

      if (n12 .eq. 1) k = k + 1 
      goto 110

!     third-4th bytes of a 4-byte value:
150   if (n12.eq.0) then 
      write(99,151) 
151   format(2x,'  Warning: found 2nd half of 4-byte with no 1st half')  
      goto 110                      ! debug:  only if first value
      endif  
      if (n12.gt.1) then 
      write(99,152) n12 
152   format(2x,'  Warning: found 2nd half of 4-byte after ',i4,
     &           ' first halves')  
      endif  
      n12 = 0  
      if ( iand(itmp,  16384) .ne. 0 ) then
      sf = -1.
      else
      sf = 1.
      end if
      if ( iand(itmp,    512) .ne. 0) sf = sf * 0.0001 
      if ( iand(itmp,    256) .ne. 0) sf = sf * 0.01 
      if ( iand(itmp, -32768) .ne. 0) sf = sf * 0.1 
      acr10(k) = sf*float( 256*iand(icr10(j), 256) + 
     &                     256*iand(itmp, 255) + iand(icr10(j), 255) )
      itmp=0
      if (k .ge. 2048) goto 160
      goto 110
!     end of binary decoding ===============================================================

!     next step:  have a full decoded output block (acr10), so break into data lines:

160   jw1=1                         ! starts at beginning of ascii array
165   continue
      if (jw1 .ge. 2049) then       ! at end, go get next ***changed from .eq. *** (1/23/02)
       k = 0
       if(verbose .or. noerr.gt.0) then 
        write(99,298) jrec, j, norec, julday, tlast
298     format('  Output block end (in block ',i5,', j=',i4,') at record',
     &         ' number ',i7,', day ',i3,', ',f9.3,' s')
       endif
       niderr = 0                       ! re-zero number of id errors in record
       ntierr = 0                       ! re-zero number of time-mismatch errors in record
       goto 110
      end if
      if (startid.ne.0) then            ! checks if we are looking for specific labels
       if (acr10(jw1).eq.startid) then
       nid = nid+1
       endif
       if(nid.ge.nidstart .and. (.not.output) ) then
        output = .true.
        write(99,*)'  Starting data output.'
        write(*,*) '  Starting data output.'
       elseif(nid.gt.nidend) then
        write(99,*)'  Reached end of ID-block.  Exiting this file.'
        write(*,*) '  Reached end of ID-block.  Exiting this file.'
        goto 980                         ! exit reading this file
       endif
      endif

      if (acr10(jw1) .eq. 200.) then     ! data-line id = 200.  (eddy1,2 or backup sonic)
       jw2 = jw1 + nf200                  ! check end of this data-line to see if acr10 ends
       if (jw2 .gt. 2048) goto 170        ! in middle of this line.
       id = 200
       do i=1,nf200                       ! loop thru data fields on this line.
        jw = jw1+i                         ! no date/time fields (compare to 201 etc lines).
        a200(i) = acr10(jw)
       end do
       jw1=jw2+1                          ! set jw1 to beg of next line
       if(output) goto 180                ! goes to write 200 data
       elseif (acr10(jw1) .eq. 201.) then ! data-line id = 201 (eddy1,2 or profile or ground)
        jw2 = jw1+nf201
       if (jw2.gt.2048) goto 170
        id = 201
        julday = int(acr10(jw1+1))         ! first 5 variables are date/time/status
        hhmm   = int(acr10(jw1+2))
        ssec   = acr10(jw1+3)
       if (dtype.ne.'SN') then 
        status1 = int(acr10(jw1+4))
        status2 = int(acr10(jw1+5))
        
!.....nrc 170207.................................................................................
      if ((idate.ge.120321).and.(idate.lt.151231)) then	 ! after 170109 problems with the loggers in the years 2012-2013
       if (status1.eq.3842) status1 = 256
       if (status2.eq.2311) status2 = 1536
      end if
!................................................................................................

      do i=1,nf201-5                     ! n fields left to store less by 5 (date/time/status)
      jw      = jw1+i+5                  ! jw index starts with 1 after 5 date/time/status vars
      a201(i) = acr10(jw) 
      end do
      endif  
      jw1     = jw2+1                         ! sets jw1 to beg of next line
      if(output) goto 185                     ! goes to write 201 data
      elseif (acr10(jw1).eq.202.) then        ! eddy, profile, ground with 202 line
      jw2    = jw1+nf202
      if (jw2 .gt. 2048) goto 170
      id     = 202
      julday = int( acr10(jw1+1) )
      hhmm   = int( acr10(jw1+2) )
      ssec   = acr10(jw1+3)
      do i=1,nf202-3
      jw      = jw1+i+3
      a202(i) = acr10(jw)
      end do
      jw1 = jw2+1                             ! sets jw1 to beg of next line
      if(output) goto 190                     ! goes to write 202 data
      elseif (acr10(jw1).eq.203.) then        ! eddy, profile, ground with 203 line
       jw2 = jw1+nf203
      if (jw2 .gt. 2048) goto 170
      id = 203
      julday = int( acr10(jw1+1) )
      hhmm   = int( acr10(jw1+2) )
      ssec   = acr10(jw1+3)
      do i=1,nf203-3
       jw = jw1+i+3
       a203(i) = acr10(jw)
      end do
      jw1 = jw2+1                             ! sets jw1 to beg of next line
      if(output) goto 195                     ! goes to write 203 data
      elseif (acr10(jw1) .eq. 100.) then      ! initialization loop header
      jw2=jw1+5
      if (jw2 .gt. 2048) goto 170
      id = 100
      iyear  = int( acr10(jw1+1) )
      julday = int( acr10(jw1+2) )
      hhmm   = int( acr10(jw1+3) )
      hh     = int(hhmm/100)
      mm     = hhmm-hh*100
      ssec   = acr10(jw1+4)
      idsite = int(acr10(jw1+5))
      write (*,301) nid,iyear, julday, hh, mm, ssec,jrec, jw1
      write (99,301) nid,iyear, julday, hh, mm, ssec,jrec, jw1
301   format('  ID=100 (initialization header # ',i3,') at year=',i5,
     &' day:',i4,i3,':',i2,':',f5.2,/,9x,'(rec ',i5,',indx ',i4,')')
      jw1 = jw2+1
      goto 165                                ! gets next data line.
      elseif (acr10(jw1) .eq. 300.) then      ! low-power output.
      jw2 = jw1+4
      if (jw2 .gt. 2048) goto 170
      julday = int( acr10(jw1+1) )
      hhmm   = int( acr10(jw1+2) )
      ssec   = acr10(jw1+3)
      hh=int(hhmm/100)
      mm      = hhmm-hh*100
      tpwroff = acr10(jw1+4)*2                        ! time that power was off (min)
      write(*,303) julday, hh, mm, ssec, tpwroff
      

      write (99,303)julday, hh, mm, ssec, hh*3600+mm*60+ssec, tpwroff
303   format('  ID=300 (low-power header), at day=',i4,i3,':',i2,':',
     &        f5.2,' (,',f8.1,') after ',f6.0,' min. no power)')
      jw1=jw2+1
      goto 165                                         ! gets next data line
      else
      write(99,309) acr10(jw1), jrec, jw1              ! unrecognized data line
309   format('  Invalid ID number (', e11.4, ' not = 100,200..203,300)',
     &  ' in/near input block # ', i6, ' (jw1= ', i4,')')
      niderr=niderr+1
      totiderr = totiderr+1
      nbadinrow = nbadinrow+1      ! also incremented if bad record error 
                                   ! occurs (e.g. time mismatch)
      ngoodinrow = 0
      jw1 = jw1+1
      endif
      goto 165    
!     gets next data line, goes to 170 (from above) if jw2>2048 (at end of 2048-integer 
!     record), but have part of data-line.
170   k=0      
      if (jw1.gt.2048) goto 110   ! if <= we're still in the middle of a 
                                    ! data-line save remaining.
      do js=jw1,2048                ! elements of acr10, move to beginning, 
                                    ! and set k = # remaining.
      k = k+1
      acr10(k) = acr10(js)
      end do
      goto 110                      ! go read another record from file


!     write ascii data into output files:

!     4 hz data :: id 200 
!     eddy systems and backup sonic:

180    norec = norec+1
      if (nfast.lt.0) goto 165             ! abort if no time markers read yet
      if (lskip) goto 165                  ! or if last time-marker no good
       nfast = nfast+1
       call vecdiv(a200,sf200,v200,nf200)   ! return: v200 = a200 / sf200
       tt    = tlast+float(nfast-1)/freqhi  ! tt (secs) interpolated from last 20n line

!     debug:  gets time directly from input line
182   format('   200-line time-error:  200-line time (secs)   =',f9.3,/,
     &       '   interpolation from 201 =',f9.3,' (diff =',f7.3,')')
      jd=julday
!      if(tt.ge.86400.) then                ! added in 13-jun-2001: fixes midnight crossings.
!       jd = jd+1
!       tt = tt-86400.
!       if (jd .eq. 366) jd = 1             ! fixes end-of-year crossing.
!      endif
      if(tt.ge.86400.) then                 ! fix midnight crossings
       jd = jd+1
       tt = tt-86400.
      endif                                 ! nrc 210521
       if ((jd.eq.367).and.(((int(year/4))-(year/4)).eq.0)) then
         print*, 'leap year'
         jd = 1    	   ! ewg: fix end-of-year crossing 
         year = year+1     ! nrc: fix end-of-year crossing 
       end if
       if ((jd.eq.366).and.(((int(year/4))-(year/4)).ne.0)) then
         jd = 1     	   ! ewg: fix end-of-year crossing 
         year = year+1     ! nrc: fix end-of-year crossing 
       end if
      !endif
            
      nf = nf200
      if((dtype.eq.'EA').or.(dtype.eq.'EB').or.(dtype.eq.'SN')) then  
      icsat=int(a200(nf200))
      nf = nf200-1                         ! doesn't count csat status as real.
      endif

      if ((dtype.eq.'EA').or.(dtype.eq.'EB').and.(v200(1).lt.0)) then 
       write(20,f200)jd,tt,(v200(i),i=1,nf),icsat,(is200(i),i=1,ns200)
      endif 
      if ((dtype.eq.'EA').and.(v200(1).gt.0)) then 
       write(20,f200)jd,tt,-v200(1),(v200(i),i=2,nf),icsat,(is200(i),i=1,ns200)
      endif 
      if (dtype.eq.'SN') then 
       write(20,f200)jd,tt,(v200(i),i=1,nf),icsat
      endif 
      goto 165

!     0.5 hz data :: id 201 ......................................................
!     eddy, profile, ground and sonic 
185   norec=norec+1
      n201 = n201+1                         ! number 201 lines
      if(n201.le.1) lastday=julday          ! initializes
      lskip=.false.
      nfast=0
      hh=int(hhmm/100)
      mm=hhmm-(hh*100)
      tlast=float(hh*3600 + mm*60) + ssec   ! time in secs only
      if (dtype .eq. 'SN') goto 165  

!     2 status words (12-bit) in every 201 line:
!     index# 12    11    10     9      8     7     6     5     4     3     2     1
!          2048  1024   512   256    128    64    32    16     8     4     2     1
!     eddy:
!      #1 spare spare spare spare  spare spare spare spare spare spare spare spare
!      #2   sd3   sd2   sd1  cal2   cal1  mux2  mux1  mux0 spare spare spare spare 
!     profile:
!      #1    ms    hs   ss  lvl_9  lvl_8 lvl_7 lvl_6 lvl_5 lvl_4 lvl_3 lvl_2 lvl_1
!      #2   sd3   sd2   sd1  cal2   cal1  mux2  mux1  mux0 calc2 calc1 zero  ls
!     ground:
!      #1 pwrco pwrdp pwrpc pwrzag spare spare spare spare co_ls co_hs co_ze co_in
!      #2   sd3   sd2   sd1 power  spare  mux2  mux1  mux0 spare pwrpf pwre2 pwre1
!
!     status will be output as follows
!     is200 (taken from last 201 line)
!     E*  5 bits:  1-3=sysdesig; 4-5=cal
!     SN  1 bit:   1= sonic status
!     is201:
!     E* 9 bits:   1-3=sysdesig; 4-5=cal; 6-8=mux; 9=cm hygr. status
!     PF 24 bits:  1-3=sysdesig; 4-5=cal; 6-8=mux; 9-10=calcon; 11-15=cal sols
!                  16-24= level sols 
!     GD 17 bits:  1-3=sysdesig; 4-6=mux; 7-13=pwr stat; 14-17=co cal sols
!     SN  4 bits:  1-3=sysdesig; 4=sonic status

!     Extracts status bits and put into integer arrays: 
      call istat(status1,isol1)           ! status bits from 201 line
      call istat(status2,isol2)
      nf=nf201-5                          ! leave out day, hhmm, ss, status1, status2

!     Eddy 1 or 2
      if(dtype(1:1).eq.'E') then  
      do i=1,5
       is200(i)=isol2(13-i)                ! sysid 3,2,1;  cal2,cal1
       is201(i)=is200(i)                   ! 1st 5 bits same for 200, 201
       is201(i+5)= isol2(8-i)              ! mux2, mux1, mux0, + 2 spare bits

!...........................................................................................
!     note:  last 2 (9,10) are spares; 9 replaced below
!...........................................................................................

      enddo
      is201(9) = 1                        ! chilled-mirror status:  assume high=ok
      if(cmlast .lt. 2.5)  is201(9)=0     ! low if mirror dirty
      if(cmlast .ge. 9.9)  is201(9)=2     ! no 203 lines read yet

!...........................................................................................
!     note:  cmlast taken from last eddy 203 line
!...........................................................................................

      endif

!     Profile
      if(dtype.eq.'PF') then 
       do i=1,12
        is201(i)=isol2(13-i)                ! sysid 3,2,1; cal2,1; mux2,1,0; calc2,1; zero, ls
       enddo
       do i=13,15
        is201(i)=isol1(12+(13-i))           ! ms,hs,ss
       enddo
       do i=16,24
        is201(i)=isol1(i-15)                ! level 1-9, column
       enddo
      end if

!     Ground
      if(dtype.eq.'GD') then 
       do i=1,3 
        is201(i)=isol2(13-i)                ! sysid 3,2,1; 
        is201(i+3)=isol2(8-i)               ! mux 2,1,0
        is201(i+7)=isol2(4-i)               ! pwr: PF, e2, e1
        is201(i+10)=isol1(13-i)             ! pwr: co, dumppump, pc
        is201(i+14)=isol1(5-i)              ! co sol:  ls, hs, zag
       enddo
       is201(7)=isol2(9)                   ! ground (main) power
       is201(14)=isol1(9)                  ! pwr: zag
       is201(18)=isol1(1)                  ! solenoid: co inlet
      endif

!     Backup sonic
      if(dtype.eq.'SN') then 
       do i=1,3
        is201(i)=isol2(13-i)
       enddo
       is201(4) = isol1(2)                 ! sonic status
       is200(1) = isol1(2)
      endif
      id201 = is201(1)*4 + is201(2)*2 + is201(3)     ! get system id stored in this line
                                                     ! 1=eddy1, 2=eddy2, 3=profile, 
                                                     ! 4=ground, etc

!     Ready to output, first check if this record is OK:
!     Checks if record time is ok -> Elapsed time must be greater than 0:
      recerr = checkrec(julday, tlast, lastday, tprev, 0.0)  
      if (recerr.ne.0) then
       write(99,350) 201, norec, jrec, jw1-nf201-1, julday, tlast, 
     &                   lastday,tprev
350    format(' ',i3,': mismatch, rec #',i7,' (block #', i5,', jw1=',i5,
     & '): jday ',i4,', ',f9.3,' sec (prev= ',i4,1x,f9.3,')')
       write(99,*)'   (num 201 = ',n201,'; sysid =',id201,')'
       ntierr=ntierr+1
       tottierr=tottierr+1

!     Now, leave only for id error:
       lskip = .true.
       goto 165
      endif

!     included by nrc:......................................................................
      if(id201.ne.idsys) then         ! if sys id on 201 line not right, skip also
       lskip = .true.
       goto 165
      endif
!     end of inclusion......................................................................

      ngoodinrow = ngoodinrow+1       ! if there are multiple bad lines in a row, waits for 
                                      ! more than one good line before outputting.
       if(nbadinrow.gt.1) then
       if(ngoodinrow.ge.1) nbadinrow=0     ! it will output next time around
       write(99,352) 201, norec, jrec, jw1-nf201-1, julday, tlast
352    format(' ',i3,': OK, but no output, record #',i7,' (block #', i5,
     &       ', jw1=',i5,'): jday ',i4,', ',f9.3,' sec')
       write(99,*)   '  (wanted 2 good 201 lines in a row, and have  ',
     &       ngoodinrow,')'
       lskip = .true.
       goto 165
      endif
      nbadinrow=0
      lastday = julday                   ! keeps record of this day
      tprev = tlast                      ! keeps record of this time (secs)
      call vecdiv(a201,sf201,v201,nf)    ! returns: v201 = a201 / sf201
      jd = julday

      !needed to hardwire positive voltages 111219 - 120104
        if ((idate.ge.111219).and.(idate.le.120104)) v200(1) = -1.*v200(1)        

      !set boundaries of what is possible for u.m.s  v.m.s  w.m.s  t.c given the write format constrains
      if ((tlast.ge.0).and.(tlast.le.86400).and.(jd.le.366).and.(jd.ge.1)) then
       write(21,f201)jd,tlast,(v201(i),i=1,nf),(is201(i),i=1,ns201)
      endif
      if(verbose .and. (ntierr.gt.0 .or. niderr.gt.0) ) then
      write(99,354) 201, norec, jrec, jw1-nf201-1, julday, tlast
354   format(' ',i3,' CR10 record #',i7,' (block #', i5,', jw1=',i5,
     &       '): jday ',i3,', ',f9.3,' sec' )
      endif                
      goto 165

!     0.5 hz data :: id 202 ............................................
!     Eddy, Profile, and Ground
190   norec=norec+1
      lskip=.false.
      nfast=0
      hh=int(hhmm/100)
      mm=hhmm-hh*100
      tlast=float(hh*3600 + mm*60) + ssec              ! time in seconds only

!     do not write repeated data (due, e.g., to errors in cr10 data dump):
!     checks if record is ok. elapsed time (from 201) must be greater than -0.5 seconds).
!     if sys id on last 201 line was not right,skip also.
      recerr = checkrec(julday, tlast, lastday, tprev,-0.5)        
      if(id201.ne.idsys)  recerr=1.   
      if( recerr.ne.0.) then
       write(99,350) 202, norec, jrec, jw1-nf202-1, julday, tlast, lastday,tprev
       write(99,*)'   (last 201 # = ',n201,'; sysid =',id201,')'
       lskip = .true.
       goto 165
      endif

!     if several bad 201 lines in a row, don't output until multiple good 201 lines found
      if(nbadinrow.gt.1) then   
       write(99,352) 202, norec, jrec, jw1-nf201-1, julday, tlast
       write(99,*)'   (want 2 good 201 lines in a row, now have ',
     &            ngoodinrow,')'
       lskip = .true.
       goto 165
      endif

!     don't let non-201 lines affect future timestamp criteria:
      if(verbose .and. (ntierr.gt.0 .or. niderr.gt.0) ) then
       write(99,354) 202, norec, jrec, jw1-nf201-1, julday, tlast
      endif
      if(dtype(1:1).eq.'E') then
       goto 165                                   ! eddy 202 line is all spares, skip output
      endif
      nf=nf202-3                                  ! nb:  3 less -- don't do day, hhmm, sec
      call vecdiv(a202,sf202,v202,nf)             ! return: v202 = a202 / sf202
      if ((tlast.le.86400).and.(tlast.ge.0).and.(julday.le.366).and.(julday.ge.1)) then
       write(22,f202)julday,tlast,(v202(i),i=1,nf) ! changed file unit # from 202
      endif
      goto 165

!     0.5 hz data :: id 203
!     eddy, profile, and ground
195   norec=norec+1
      lskip=.false.
      nfast=0
      hh=int(hhmm/100)
      mm=hhmm-hh*100
      tlast=float(hh*3600 + mm*60) + ssec                          ! time in secs only

!     do not write repeated data (due, e.g., to errors in cr10 data dump):
!     checks if record is ok. elapsed time (from 201) must be greater than -0.5 seconds).
!     if sys id on last 201 line was not right,skip also.
      recerr = checkrec(julday, tlast, lastday, tprev,-0.5)  
      if(id201.ne.idsys)  recerr=1.                                 
      if(recerr.ne.0.) then
      write(99,350) 203, norec, jrec, jw1-nf202-1, julday, tlast, lastday,tprev
      write(99,*)'   (last 201 # = ',n201,'; sysid =',id201,')'
      lskip = .true.
      goto 165
      endif

!     if several bad 201 lines in a row, don't output until multiple good 201 lines found.
      if(nbadinrow.gt.1) then   
       write(99,352) 203, norec, jrec, jw1-nf201-1, julday, tlast
       write(99,*)'   (want 2 good 201 lines in a row, now have ',
     &            ngoodinrow,')'
       lskip = .true.
       goto 165
      endif

!     doesn't let non-201 lines affect future timestamp criteria:
      if(verbose .and. (ntierr.gt.0 .or. niderr.gt.0) ) then
       write(99,354) 203, norec, jrec, jw1-nf201-1, julday, tlast
      endif
      if(dtype(1:1).eq.'E') then 
       cmlast = a203(1)/sf203(1)     ! save chilled mirror status. 
                                     ! for output on next 201 line
        print*, 'cmlast split:  ', cmlast

      goto 165                       ! remainder of eddy 203 line is all spares, skip output.
      endif
      nf=nf203-3                            ! nb:  3 less -- don't do day, hhmm, sec
      call vecdiv(a203,sf203,v203,nf)       ! returns v203 = a203 / sf203
      if ((tlast.le.86400).and.(tlast.ge.0).and.(julday.le.366).and.(julday.ge.1)) then
       write(23,f203)julday,tlast,(v203(i),i=1,nf)
      endif
      goto 165

!     error exit (defunct)
900   write(99,*) ' '
      write(99,901) ierrlim, ifilen(infile)
901   format('  Error exit (number of errors >=',i5,') from ', a12)
      write(99,*) '  Abort processing: '
      write(99,*) '   at record (4096 bytes/rec) # jrec =',jrec
      write(99,*) '   # lines (cr10 record) read: norec = ',norec
      write(99,*) '   # errors reading, decoding: noerr = ',noerr        ! terminal output
      write(*,901) ierrlim, ifilen(infile)
      write(*,*)  '  Abort processing at error # ', noerr
      goto 992

!     normal exit
!     user-defined exit (record #jrend reached)
980   write(99,*) ' '
      write(99,*) '  Warning:  exit from ',ifilen(infile)
      write(99,*) '  because record limit exceeded at record ',
     &            '(4096 bytes/rec) # jrec =',jrec
      write(99,*) '  Blocks (of 2048 2-byte words) read: jrec     = ',
     &             jrec
      write(99,*) '  CR10 data-lines read:              norec     = ',
     &             norec
      write(99,*) '  Number of 201 cr10 lines            n201     = ',
     &             n201
      write(99,*) '  Total # of 2-byte words read = jrec*2048     = ',
     &             jrec*2048
      write(99,*) '  Error summary (2-byte words): '
      write(99,*) '  Number of blocks with >=1024 bad-data errors =',
     &             nblkerr, ' (',nblkerr/jrec*100, '%)'
      write(99,*) '  Total number of bad-data errors              =',
     &             toterr,' (',toterr/(jrec*2048)*100,'%)'
      write(99,*) '  Error summary (within records): '
      write(99,*) '  Number of invalid-id errors                  = ',
     &             totiderr
      write(99,*) '  Number of 201 time-base/mismatch errors      = ',
     &             tottierr
      goto 992

!     exit due to end-of-file
990   write(99,*) ' '
      write(99,*) '  Normal exit at end of ',ifilen(infile)
      write(99,*) '  (reached end-of-file)'
      write(99,*) '  Blocks (of 2048 2-byte words) read: jrec  = ',jrec
      write(99,*) '  CR10 data-lines (records) read: norec     = ',norec
      write(99,*) '  Number of 201-line records          n201  = ',n201
      write(99,*) '  Total # of 2-byte words read = jrec*2048  =',
     &             jrec*2048
      write(99,*) '  Error summary (2-byte words):'
      write(99,*) '  Number blocks with >=1024 bad-data errors = ',
     &             nblkerr,'(',nblkerr/jrec*100, '%)'
      write(99,*) '  Total number of bad-data errors           = ',
     &             toterr,' (',toterr/(jrec*2048)*100,'%)'
      write(99,*) '  Error summary (within records):'
      write(99,*) '  Number of invalid-id errors               = ',
     &             totiderr
      write(99,*) '  Number of 201 time-base/mismatch errors   = ',
     &             tottierr

992   close(99)     

!     end of main program loop
995   end do                     ! end of 'do infile=1, nifiles' (loop to next input file).
      print *, '   End of program LSPLIT'
      end       


!===========================================================================================
!                              subroutines and functions
!===========================================================================================

!===========================================================================================
!                               subroutine vecdiv                                          !
!-------------------------------------------------------------------------------------------
      subroutine vecdiv(a,b,ab,n)
      real*4 a(n),b(n),ab(n)
      do i=1,n
       ab(i)=a(i)/b(i)
       if( abs(ab(i)) .gt. 999.) ab(i)=-99.999
      enddo
      return
      end
      function idotprod(a,b,n)
      integer a(n),b(n)
      idotprod=0.
      do i=1,n
       idotprod=idotprod+(a(i)*b(i))
      enddo
      return
      end
!-------------------------------------------------------------------------------------------
!                            end of subroutine vecdiv                                      !
!===========================================================================================



!===========================================================================================
!                                subroutine istat                                          !
!-------------------------------------------------------------------------------------------
!     Extracts status bits from integer is1, puts into integer array isol
!-------------------------------------------------------------------------------------------
      subroutine istat(is1,isol)
      integer is1,isol(12)
      do i=1,12
       isol(i)=0
       if (btest(is1,i-1)) isol(i)=1
      enddo
      return
      end
!-------------------------------------------------------------------------------------------
!                              end of subroutine istat                                     !
!===========================================================================================



!===========================================================================================
!                               subroutine infiles                                         !
!-------------------------------------------------------------------------------------------
      subroutine infile_spl(ymd,types,names,num)
      character*6 ymd, types
      character*12 names(6)
      character*12 fchars
      integer fnum
      fchars='EAEBPFGDSNNT'
      num= lentrim0(types,len(types) )
      do i=1,num
       fnum = ichar(types(i:i))-48
       names(i)=ymd//fchars( ((fnum-1)*2+1):((fnum-1)*2+2) )//'.dat'
      enddo
      return
      end
!-------------------------------------------------------------------------------------------
!                            end of subroutine infiles                                     !
!===========================================================================================



!===========================================================================================
!                                fuction checkrec                                          !
!-------------------------------------------------------------------------------------------
!     checks if the time just read (jday, tsec) is valid, and consistent with previous time 
!     stamp (prevday, tsecprev), i.e., elapsed time must be > tallow (sec). returns value=0 
!     if ok, nonzero otherwise.
!     nrc changed the var name tprev for tbr to avoid confusion
!     nrc increased the precision of the real variables 4 to 6 and 8 to 10
!-------------------------------------------------------------------------------------------
      function checkrec(jday,tsec,prevday,tsecprev,tallow)
      logical newyear  
      integer jday, prevday                               !nrc added jday
      real*4 checkrec, error, tsec, tsecprev, tallow
      real*8 tnow, tbfr, tdiff
      tsec = dble(tsec)                                   !added by nrc 091009
      tsecprev = dble(tsecprev)                           !added by nrc 091009
      error = 0.
      tnow = tsec + (jday*86400.)                         ! times in decimal days
      tbfr = tsecprev + (prevday*86400.)                  ! (now in seconds) nrc
      newyear = .false.  
      if (jday .eq. 1 .and. prevday .ge. 366) newyear = .true.
      if (jday .eq. 1 .and. prevday .ge. 366) newyear = .true.
      tdiff = tnow-tbfr                                   ! elapsed time (in seconds)

!     criteria for error: ==================================================================
      if(jday.lt.0 .or. jday.gt.366)     error=1.           ! day out of range
      if(tsec.lt.0. .or. tsec.gt.86410.) error=1.           ! time-of-day out of range
      if(tdiff.le.tallow .and. .not.newyear) error=1.       ! elapsed must be > tallow
      if(tdiff.gt. (7.*(86400.)) )       error=1.           ! elapsed time too big (> 7 days)
      checkrec = error
      end
!-------------------------------------------------------------------------------------------
!                            end of fuction checkrec                                       !
!===========================================================================================
!
!-------------------------------------------------------------------------------------------
      include 'lsubs.for' 
!-------------------------------------------------------------------------------------------
!
