!===================================================================================================
!                                         lfileagr.for          
!===================================================================================================
!  'include' file to open input and output files for laggr6.for:
!       - get all file names, including appropriate directory paths
!       - open input files (*.201, *.202, *.203)
!       - open output files (*.hkp, *.clm, *.era) and output initial headers, etc
!---------------------------------------------------------------------------------------------------
!     gets filename:
!     raw data input (*.201,202,203 data)
      fnamein = ifilen(infile)
!---------------------------------------------------------------------------------------------------
      call getfullpath( rootsplit,yymmdd, dirsplit, fnamein, lfnin) 
!---------------------------------------------------------------------------------------------------
!     output data (processed *.hkp, *.clm files) 
      fname = ifilen(infile)
!---------------------------------------------------------------------------------------------------
      call getfullpath( rootprocin, yymmdd, dirprocout, fname, lfn)
!---------------------------------------------------------------------------------------------------

!...................................................................................................
!     note: use of 'dirprocout' directory.  (?)
!...................................................................................................

      write(* ,*) '  Opening input files ', fnamein(1:lfnin)

!     input files (common to PF, GD)

!     201 file
      open(21,file=fnamein(1:lfnin)//'.201',status='old')  
      read(21, *)  ! header line 1   GD has 3-line header
      read(21, *)  ! header line 2   ( 1 more in PF )
      read(21, *)  ! header line 3

!     202 file
      open(22,file=fnamein(1:lfnin)//'.202',status='old') 
      read(22, *)  ! header line 1   both PF & GD have 5-line header
      read(22, *)  ! header line 2
      read(22, *)  ! header line 3
      read(22, *)  ! header line 4
      read(22, *)  ! header line 5

!     203 file
      open(23,file=fnamein(1:lfnin)//'.203',status='old') 
      read(23, *)  ! header line 1   GD has 2-line header 
      read(23, *)  ! header line 2   (1 more for PF)

!    output files: 
      open(99, file=fname(1:lfn)//'.era',status='unknown')  ! error file.
      open(16, file=fname(1:lfn)//'.hkp',status='unknown')  ! housekeeping data.
      open(17, file=fname(1:lfn)//'.clm',status='unknown')  ! climate data.

      write(99,55) platform, lname, iyear, imonth, iday,ihrs,imin,isec
55    format(a6,' ',a7,':  on date (yy-mm-dd) ', i4,'-',i2.2,'-',i2.2,
     &       ' (',i2.2,':',i2.2,':',i2.2,'), process:' )
      write(99,*) '  Opened input files ', fnamein(1:lfnin)//
     &            '.201, 202, 203'

!     file-specific set-up:

!     profile system set-up:
      if(dtype.eq.'PF') then   

!     opens additional output if this is profile data for sunshine sensor:
      open(18, file=fname(1:lfn)//'.sun',status='unknown')

!     reads a few more header lines...
      read(21, *)  ! header line 4 in 201 file
      read(23, *)  ! header line 3 in 203 file

!     headers for profile output files: 

!     1) PF*.hkp:
      write(16,61) yymmdd,period,platform, lname, iyear,imonth,iday, 
     &             ihrs,imin, isec
61    format(' 3 23 LBA Profile housekeep: ',a6,' (', f5.2, ' min) (',
     &       a6,' ',a8,' on: ', i4,'-',i2.2,'-',i2.2,', ',i2.2,':',
     &       i2.2,':',i2.2,')' )
      write(16,*) 'jdstart.gmt    doy    pcon    flsam    flcal    pcell    rl    tcell    tdetect    tpump'
      write(16,*) 'sid    c01    c10    c11    ccon    zag    ls    ms    hs    ss    lv1-4    lv5-8    col'

!     2) PF*.clm:
      write(17,71) yymmdd, period, platform, lname, iyear, imonth, iday,
     &             ihrs,imin, isec
71    format('5  32  LBA Profile climate: ',a6, ' (', f5.2, ' min) (',
     &        a6,' ',a8,' on: ', i4,'-',i2.2,'-',i2.2,', ',i2.2,':',
     &        i2.2,':',i2.2,')' )
      write(17,*) 'jdstart.gmt    doy    netrad    snetrad    par1up    spar1up    par1dn    spar1dn    par2up'
      write(17,*) 'spar2up    tair1.c    tair2    tair3    tair4    tair5    tair6    tair7'
      write(17,*) 'tair8    wd    swd    ws1.m.s    sws1    ws2     sws2    ws3'
      write(17,*) 'sws3    ws4    sws4    rain    netrad_cal    par_cal    tair_cal'       

!     3)PF*.sun
      write(18,81) yymmdd, period, platform, lname, iyear, imonth, iday,
     &             ihrs,imin
81    format('2  6  LBA Profile sunshine:  ',a6, ' (', f5.2,
     &       ' min) (',a6,' ',a8,' on: ', i4,'-',i2.2,'-',i2.2,', ',
     &       i2.2,':',i2.2,':',i2.2,')' )
      write(18,*) 'jdstart.gmt    doy    totrad.umol    difrad.umol    stot    sdiff' 
       
       
!     other parameters for profile
!     input-related parameters:
      n1=n1PF                       ! number of sensor inputs on 201 line (not including date).
      n2=n12PF-n1PF                 ! number of sensor inputs on 202 line.
      n3=n13PF-n12PF                ! number on 203 line (lines 1-3 minus 1-2).
      ns=nspf                       ! number status bits (end of 201 line).
      nmux=nmuxpf                   ! no mux'd variables for profile system.
      nvarin=n13PF- nmux            ! total num non-mux'd variables.

!     output-related parameters:
      nskip=10                      ! number of status to skip for status output loop .
                                    ! (sid(3), cal(2), mux(3), & cal-con(2).
      nout=noutpf                   ! number of output variables (split between *hkp and .clm).
      nhkp=nhkppf                   ! number of output variabless that are housekeeping.
      assign 161 to f16
      assign 171 to f17
      assign 181 to f18             ! for sunshine sensor data.

!     input-output map:
!     format of input-output map:
!
!     i=indx to iomap            'iomap(i)'  (indx to xo)
!     (order of input)           (order of output)
!     ----------------           ------------------
!     1 to nvarin;               output order of non-mux'd measurement
!     nvarin+1 to nvarin+nmux*8  where to output this mux'd measurement
!     nvarin+nmux*8+1  to        where to put this status output
!     nvarin+nmux*8+nos          (all status go to hkp file)
!     nos = number of status

      do i=1,60     ! zero output map 
      iomap(i)=0    ! (i.e. default is to not output input i)
      end do

!     reference:  profile system vars, 1-nvarin:          
!     id=201     1  pcon                 15  netrad
!                2  flsam                16  netrad_cal
!                3  flcal                17  par1_up
!                4  pcell                18  par1_down
!                5  rl                   19  par2_up
!                6  tcell                20  par_cal
!                7  tdetect           21-23  spare      (21 & 22 are now sunshine sensor)
!                8  tpump                24  tair_cal
!                9  co2               25-32  tair1 to tair8
!               10  h2o        id=203 33-36  ws1 to ws4
!            11-12  spare             37-40  spare

!     id=202    13  wdir                 41  rain
!               14  spare      varin:    42  spare
!                             (status    43  sid    composite of is201(1-3)
!                              extras    44  cal=01 composite of is201(4-5)
!                                        45  cal=10 composite of is201(4-5)
!                                        46  cal=11 composite of is201(4-5)
!                                        47  cal-cntrl (composite of 4-5,9-10)
!                                     48-52  span gas      from is201(11-15)
!                                     53-62  level sample  from is201(16-24)

!     fill-in according to output order (i.e. right-hand side of below):

!     first, do hkp output
      do i=1,8                        ! 1st 8 output vars (to hkp) are sensor ouputs
      iomap(i)=i                      ! inputs 1-8 are outputs 1-8
      end do                           !  pcon,flsam,flcal,pcell,rl,tcell,tcinlet,tpump
      ns0=nvarin+nmux*8
      do i=1,13                       ! next 13 (also to hkp) are all the status vars
      iomap(ns0+i) = 8+i              ! inputs ns0+i=42+(1-13)=43-55 -> outputs 9-21
      end do

!     now, do.clm output: (offsetting from ones already done 8+13=21=nhkp)
      iomap(15) = nhkp+1              ! input 15=netrad is output nhkp+1
      sdmap(15) = nhkp+2              ! sd of input 15 is output nhkp+2
      ioffs=nhkp+2
      npar=3
      do i=1,npar                     ! input = 17-19 (par1_up,dwn, par2)
      iomap(16+i) = ioffs+(i-1)*2+1   ! mean: last output + (1,3,5)
      sdmap(16+i) = ioffs+(i-1)*2+2   ! sdev: last output + (2,4,6)
      end do
      ioffs= ioffs+npar*2 
      nair=8
      do i=1,nair                     ! input = 8 air temps
      iomap(24+i) = ioffs + i         ! output= last + 1-8
      end do
      ioffs = ioffs + nair
      iomap(13) = ioffs + 1           ! input 13=wind dir is output ioffs+1
      sdmap(13) = ioffs + 2           ! sdev of wind dir
      ioffs = ioffs +2                ! accumulated output offset
      nws=4
      do i=1,nws                      ! input 33-37 = 4 windspeeds (ws1-ws4) + 4 sdev's
      iomap(32+i) = (ioffs)+(i-1)*2+1 ! mean: last+(1,3,5,7)
      sdmap(32+i) = (ioffs)+(i-1)*2+2 ! sdev: last+(2,4,6,8)
      end do
      ioffs = ioffs + nws*2           ! accumulated output offset
      iomap(41) = ioffs + 1           ! 41 = rain

!     calibration outputs (need to add 3 to nout if outputting these)         
      iomap(16) = ioffs + 2           ! 16 = netrad cal input
      iomap(20) = ioffs + 3           ! 20 = par_cal
      iomap(24) = ioffs + 4           ! 24 = tair_cal

!     sunshine sensor outputs
      iomap(21) = nout+1
      iomap(22) = nout+2    
      sdmap(21) = nout+3
      sdmap(22) = nout+4     
       
!     ground system set-up:
      elseif(dtype.eq.'GD') then  

!     gets conversion coefficients from lconv.inp file
!---------------------------------------------------------------------------------------------------
      call readconv(4)
!---------------------------------------------------------------------------------------------------
!     input-related parameters:
      n1=n1GD                         ! number of sensor inputs on 201 line.
      ns=nsgd                         ! number of status bits on 201 line.
      n2=n12GD-n1GD                   ! number of fields on 202 line.
      n3=n13GD-n12GD                  ! number of fields on 203 line.
      nmux=nmuxgd                     ! last 3 inputs are mux'd vars (8 sensors each).
      nvarin=n13GD - nmux             ! total number non-mux'd variables.

!     output-related parameters:
      nskip=6                         ! number status to skip for status output (sid, mux).
      nout=noutgd                     ! number of output variables (including status).
      nhkp=nhkpgd                     ! number of output vars that are housekeeping (incl status).
      assign 162 to f16               ! output format label number for GD hkp.
      assign 172 to f17 

!   headers for ground-system output files: 

!   GD*.hkp
      write(16,62) yymmdd,period,platform, lname, iyear,imonth,iday, 
     &             ihrs,imin,isec
62    format('4 27  LBA Ground housekeep: ',a6, ' (', f5.2,' min)  (',
     &        a6,' ',a8,'on: ',
     &        i4,'-',i2.2,'-',i2.2,', ',i2.2,':',i2.2,':',i2.2,')' )
      write(16,*) 'jdstart.gmt    doy    pconco    flco    tcellco    ttrap    tpump_d    '
      write(16,*) 'pss    phs    pms    pls    phsco    plsco    pzag    '
      write(16,*) 'sid    pwr    pwrp    pwr2    pwr1    pwrco    pwrpmp    pwrpc    pwrzg    lsco    hsco    zagco    samco'


      write(17,72) yymmdd,period,platform, lname, iyear,imonth,iday, 
     &             ihrs,imin,isec
72    format('5 27  LBA Ground climate: ',a6, ' (', f5.2,' min)  (',a6,
     &       ' ',a8,'on: ',i4,'-',i2.2,'-',i2.2,', ',i2.2,':',i2.2,':',
     &       i2.2,')' )
      write(17,*) 'jdstart.gmt  doy     pamb    pars1    pars2    pars3    pars4    pars5    pars6'
      write(17,*) 'pars7    pars8     tsa1.c  tsa2     tsa3    tsa4     tsa5     tsa6'
      write(17,*) 'tsa7     tsa8     tsb1     tsb2     tsb3    tsb4     tsb5     tsb6'
      write(17,*) 'tsb7     tsb8'

!     set-up input-output map for ground system:
      do i=1,60                       ! zero output map 
      iomap(i)=0                      ! (i.e. default is to not output input i)
      end do
!
!     reference:  ground system input-var index numbers (1-nvarin):
!     id=201     1  co             id=203  32  p_ss  (cylinder pressures)
!                2  pcon_co                33  p_hs
!                3  fl_co                  34  p_ms
!                4  stat_co                35  p_ls
!                5  tcell_co               36  p_cohs
!                6  ttrap_co               37  p_cols
!                7  tpump_dump             38  p_zag
!     id=202  8-23  ts1 - ts16             39  p_ambient
!            24-31  par1 - par8   (status  40  sid (combined)
!                                          41-48  power status
!                                          49-52  co solenoid status
!
!     nb:  in order of output (number on right-hand-side below) 
!                             (ref headers above, for order of output)
!     number as index to iomap (on left-hand side below) is input order
!                             (ref input-var index numbers above)

!     first, do hkp output:
      do i=1,2                        ! 1st 4 output variables (to hkp) are sensor ouputs.
      iomap(i+1)=i                    ! inputs i+1=2-3 are outputs 1-2 (skip 1st input: co).
      iomap(i+4)=i+2                  ! inputs i+4=5-6 are outputs 3-4 (skip 4th:  stat_co).
      end do                          ! pcon_co, fl_co, tcell_co, ttrap_co.
      iomap(7) = 5                    ! input 7 (tpump_dmp) is output 5.
      ioffs = 5
      do i=1,7                        ! next 7 (also to hkp) are cylinder pressures.
      iomap(31+i) = ioffs+i           ! inputs 31+i=32-38 are outputs 6-12.
      end do
      ioffs = ioffs+7
      ns0 = nvarin+nmux*8             ! zero for status bits (last of all other inputs).
      do i=1,13                       ! 13 status outputs.
      iomap(ns0+i) = ioffs+i          ! inputs 39+i=40-52 are outputs 12+i=13-25.
      end do

!     then, do.clm output: (offsetting from ones already done 5+7+13=25 =nhkp)
      ioffs =nhkp
      iomap(39) = ioffs+1             ! 39 = ambient pressure is.clm output 1.
      ioffs=ioffs+1
      do i=1,8                        ! next 8 outptus are the 8 surface par sensors.
      iomap(23+i) = ioffs+i           ! inputs 24-31 are next 8 outputs.
      end do
      ioffs= ioffs+8
      do i=1,16                       ! next 16 outputs are the 16 soil temperature probes.
      iomap(7+i) = ioffs+i
      end do
      endif                           ! profile/ground set-up options.

