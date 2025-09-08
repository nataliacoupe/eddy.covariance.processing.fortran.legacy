!===================================================================================================
!                                lba program: linitial7.for
!===================================================================================================
!     by scott saleska
!
!     common code for initializing variables in lba program modules. started 22-sept-01
!
!     26-sep-06:  modified for the 'new' tower installation sept '06.
!     no eddy systems, and we are now using 'sn' for sonic data.
!     15-oct-04 : changed root path to allow for yearly directory structure modified (ewg).
!     10-may-02 : increased length of arg and trimmed spaces before use allows change in arg length 
!                 definition (lparams7.for).
!     09-may-02 : slight change to user message.
!     23-jan-02 : add ability to modify dirprocin via command line.
!---------------------------------------------------------------------------------------------------

!     default initial parameters
!     intialize directory structure (defaults set in lpaths7.for)
      rootinp = root                 ! default root paths for directories
      rootsplit = root
      rootprocin = root
      rootprocout = root
      dirsplit = split               ! default data-directory names
      dirprocin = process
      dirprocout = process
      verbose = .false. 
      data period/30.0/              ! averaging interval (minutes)
      data cc/60*1./, offset/60*0./  ! conversion coefficients, offsets
!     get command line info
      errflag=0
      pflag=0
      navg=0
!---------------------------------------------------------------------------------------------------
      call readcmdline(cmdline,cmdlen)
!---------------------------------------------------------------------------------------------------
!     defined in include file lcmdline.for, which is platform dependent.
!     lcmdline on watcom platform will call fgetcmd
!     on gnu fortran will use: 
!     function iargc() to get # of args
!     and subroutine getarg(i,arg) to retrieve the ith arg
!     readcmdline:  returns command line in cmdline, length in cmdlen

      narg=iargf(cmdline)  ! get num args in cmdline

!     if command line not right, display info, stop:
      if (narg .lt. 1 .or. cmdlen.lt.6) then
      print *, '  Error: incorrect cmdline input.'
      goto 5
      end if
!     retrieve contents of command line 
!     (note:  defaults ftypes set in main code)
      errflag=0
      pflag=0
!---------------------------------------------------------------------------------------------------
      call mygetarg(cmdline,1,arg)
!---------------------------------------------------------------------------------------------------

!.....inserted 1-26-02 to allow for possibility of '-date' command-line.............................
      if(arg(1:1).ne.'-'.and.lentrim0(arg,len(arg)).eq.6) then
      yymmdd=arg(1:6)
      ifirst=2
      end if
!.....end of insert.................................................................................

      if (narg .ne. 1) then    ! if =1, is just yymmdd, ok
!     option:  insert check that narg-1 is even
      do ith=2,narg,2

!.....inserted 1-26-02:  (for compatibility with lsplit)............................................

!---------------------------------------------------------------------------------------------------
      call mygetarg(cmdline,ith,arg)
!---------------------------------------------------------------------------------------------------
      larg = lentrim0(arg,len(arg))          ! ewg/srs 
      if (arg(1:larg) .eq. '-date') then
!---------------------------------------------------------------------------------------------------
      call mygetarg(cmdline,ith+1,arg)
!---------------------------------------------------------------------------------------------------
      yymmdd=arg(1:6)
      end if
      if (arg(1:larg) .eq. '-start') then    ! if specified, block num to start with
!---------------------------------------------------------------------------------------------------
      call mygetarg(cmdline,ith+1,arg)
!---------------------------------------------------------------------------------------------------
      read(arg,*)jrstart
      end if
      if (arg(1:larg) .eq. '-end') then      ! if specified, block num to end with
!---------------------------------------------------------------------------------------------------
      call mygetarg(cmdline,ith+1,arg)
!---------------------------------------------------------------------------------------------------
      read(arg,*)jrend
      end if
      if (arg(1:larg) .eq. '-label') then    ! don't start until this label # found
!---------------------------------------------------------------------------------------------------
      call mygetarg(cmdline,ith+1,arg)
!---------------------------------------------------------------------------------------------------
      read(arg,*)startidin                   ! note:  -label arg can be paired with -numlab arg
      end if
      if (arg(1:larg) .eq. '-numlab') then   ! number of labels to take 
!---------------------------------------------------------------------------------------------------
      call mygetarg(cmdline,ith+1,arg)
!---------------------------------------------------------------------------------------------------
      read(arg,*) startnum                   ! to be split into start and num parts
!     startnum:  start = number of labels to find before outputting data.
!     numid = number of additional labels before stopping output.
      end if
      if (arg(1:larg) .eq. '-ascii') then    ! input files are ascii, not binary
      ascii = .true.                         ! (n.b. this does not quite work)
      end if
      if (arg(1:larg) .eq. '-verbose') then
      verbose = .true.                       ! more reporting in log files
      end if
!.....end of insert.................................................................................

      if (arg(1:larg) .eq.'-files'.or. arg(1:larg).eq.'-file') then
!---------------------------------------------------------------------------------------------------
      call mygetarg(cmdline,ith+1,arg)
!---------------------------------------------------------------------------------------------------
      tempstr=arg(1:5)//'    '                ! '123456' = all file types (currently, max 5 files possible)
      do i=1,lentrim0(tempstr, len(tempstr))  
      fnum = ichar(tempstr(i:i))-48           ! converts to integer
      if ((fnum .lt. 1).or.(fnum.gt.5)) then
      errflag=1
      tempstr(i:i) = ' '  
      end if
      end do
      if (lentrim0(tempstr, len(tempstr)).ge.1) ftypes=tempstr
!     if acceptable file types, override default
!---------------------------------------------------------------------------------------------------
      elseif ((pflag.eq.0) .and. 
     &       ((arg(1:larg).eq.'-minutes').or.(arg(1:larg).eq.'-min'))) 
     &              then
      call mygetarg(cmdline,ith+1,arg)
      read(arg,*)period
      pflag=1                          ! period set as given in minutes
!---------------------------------------------------------------------------------------------------
      elseif ((pflag.eq.0) .and. 
     &       ((arg(1:larg).eq.'-seconds').or.(arg(1:larg).eq.'-sec'))) 
     &              then
      call mygetarg(cmdline,ith+1,arg)
      read(arg,*)period
      period = period/dble(60.)        ! if period is in seconds convert to minutes
!---------------------------------------------------------------------------------------------------
      elseif ((pflag.eq.0) .and. (arg(1:larg).eq.'-raw')) then
      call mygetarg(cmdline,ith+1,arg)
      read(arg,*)period
      navg = period
      period = navg/dble(freqhi*60.)   ! if period is in hertz convert to minutes  
!     getting filtering mask code
!---------------------------------------------------------------------------------------------------
      elseif (arg(1:larg).eq. '-filter') then
      call mygetarg(cmdline,ith+1,arg)
      filtcode=arg(1:5)
      call getfilterset(filtcode,filtset)
!---------------------------------------------------------------------------------------------------
      elseif (arg(1:larg).eq. '-dirin') then
      call mygetarg(cmdline,ith+1,arg)
      dirprocin = arg
!---------------------------------------------------------------------------------------------------
      elseif (arg(1:larg).eq. '-dirout') then
      call mygetarg(cmdline,ith+1,arg)
      dirprocout = arg
!---------------------------------------------------------------------------------------------------
      elseif (arg(1:larg).eq. '-split') then
      call mygetarg(cmdline,ith+1,arg)
      dirsplit = arg
!---------------------------------------------------------------------------------------------------
      elseif (arg(1:larg).eq. '-root') then
      call mygetarg(cmdline,ith+1,arg)
      print *, ' arg(-root) = "'//arg//'"'
      rootinp = arg
      rootsplit = arg
      rootprocin = arg
      rootprocout = arg
!---------------------------------------------------------------------------------------------------
      elseif (arg(1:larg).eq. '-rtinp') then
      call mygetarg(cmdline,ith+1,arg)
      rootinp = arg
!---------------------------------------------------------------------------------------------------
      elseif (arg(1:larg).eq. '-rtsplit') then
      call mygetarg(cmdline,ith+1,arg)
      rootsplit = arg
!---------------------------------------------------------------------------------------------------
      elseif (arg(1:larg).eq. '-rtproci') then
      call mygetarg(cmdline,ith+1,arg)
      rootprocin = arg
!---------------------------------------------------------------------------------------------------
      elseif (arg(1:larg).eq. '-rtproco') then
      call mygetarg(cmdline,ith+1,arg)
      rootprocout = arg
!---------------------------------------------------------------------------------------------------
      else errflag=1
      end if          
      end do
      end if
!---------------------------------------------------------------------------------------------------
      if (period .gt. 120) then
      errflag=1
      end if
      if (errflag.ne.0) then          
      print *, '  Unrecognized or unallowed arguments in cmdline'
      goto 5                               ! prints error message and stop.
      end if   

!     get year information
      lyear = '20'//yymmdd(1:2) 
      rootsplit = root//slash//lyear
      rootprocin = root//slash//lyear
      rootprocout = root//slash//lyear
      print *, ' rootsplit: ', rootsplit  
      print *, ' rootprocin: ', rootprocin  
      print *, ' rootprocout: ', rootprocout  
      year = realchar(lyear,4)             ! year in common block.
      leapdays = int(((year-1)-2000)/4)+1  ! number of leapdays in preceding years.
      days00 = (year-2000)*365 + leapdays  ! number of days in preceding years (back to 1/1/2000)
      yrdays = 0.

!     form input file name(s)
!---------------------------------------------------------------------------------------------------
      call infiles(yymmdd, ftypes, ifilen, nifiles)  
!---------------------------------------------------------------------------------------------------
!     get 'nfiles' input file names (ifilen) from: yymmdd & ftypes
!     filenames are:  1= EAyymmdd; 2= EB, 3=PF, 4=GD, 5=SN  
!---------------------------------------------------------------------------------------------------
      call getdatetime(iyear,imonth, iday, ihrs,imin,isec, platform)
!---------------------------------------------------------------------------------------------------
!     calls platform-specific data functions (see lcmdline7.for)
      print *
      write(*,3) platform, lname, iyear, imonth, iday,ihrs,imin,isec 
3     format(a6,' ',a8,':  processing on: ', i4,'-',i2.2,'-',i2.2,' (',
     &       i2.2,':',i2.2,':',i2.2,'):')
      write(*,4) nifiles 
4     format('  ', i1, ' set(s) of input files:' ) 
      do i=1,nifiles
      print *, '   ', ifilen(i)
      end do 
      goto 6                               ! returns to main code.
5     continue 
      print *, '  Syntax: '//lname//' yymmdd'
      print *, '  or:    '//lname//' yymmdd <options>'
      print *, '  where options are:'
      if(lname.eq.'lsplit  ') then
      print *, '   -files ff  [ff is file: 1=eddy1, 2=eddy2, ',
     &         '3=profile, 4=ground, 1234=all (default)]'
      print *, '   -start 1 -end 10000 [ start= which 4096-byte block',
     &         ' to start with  '
      print *, '                         end = how many blocks to take',
     &         ' before ending ]'
      print *, '   -label xx  -numlab start.num (debug options)'
      print *, '    where:  xx = label to start outputting with ',
     &         '(e.g. label 100)'
      print *, '            start.num = decimal number, with start = ',
     &         'which label to start'
      print *, '                        and num = how many such labels',
     &         ' to take before ending'
      elseif(lname.eq.'lcaled  ') then
      print *, '   -files ff [ff is file:  1=eddy1, 2=eddy2,  12=both ',
     &         '(default)]'
      elseif(lname(1:5).eq.'laggr') then
      if(lname.eq.'laggr   ') then
      print *, '    -files ff   [ff:  3=profile, 4=ground, 34=both ',
     &         '(default)]'
      end if
      if(lname.eq.'laggred ') then
      print *, '    -files ff   [ff:  1=eddy1 (default), 2=eddy2 ]'
      end if
      print *, '    -min mm   [mm: aggregation period, in minutes ',
     &         '(default=30, max is 60)]'
      print *, '    (aggregation period can also be given as:'
      print *, '    -sec ss (where ss is seconds to aggregate)'
      print *, '    -raw nn (num of raw datapoints to aggreg.)'
      elseif(lname.eq.'lfluxes ') then
      print *, '    -files ff  [ff: 1=eddy1, 2=eddy2,  12=both ',
     &         '(default)]'
      print *, '    -min mm  [mm: aggregation period, in minutes ',
     &         '(default=30, max is 60)]'
      print *, '    (aggregation period can also be given as:'
      print *, '    -sec ss (where ss is seconds to aggregate)'
      print *, '    -raw nn (num of raw datapoints to aggreg.)'
      print *, '    -filter fffff (a filter mask code of zeroes and ',
     &         'ones)'
      print *, '     where: 1xxxx = filter if bit 15 set ',
     &         '(speed-of-sound error)'
      print *, '            x1xxx = filter if bit 14 set ',
     &         '(poor signal lock)'
      print *, '            xx1xx = filter if bit 13 set ',
     &         '(signal amplitude too high)'
      print *, '            xxx1x = filter if bit 12 set ',
     &         '(signal amplitude too low)'
      print *, '            xxxx1 = filter this sonic based on bit 14 ',
     &         'of other sonic '
      end if

!     options for all program modules:
      print *, '    -dirout <dir> (dir is output data directory.',
     &         '           [default="process"] )'
      print *, '    -dirin  <dir> (dir is input processed-data ',
     &         'directory.  [default="process"] )'
      print *, '    -split  <dir> (dir is input raw-data directory.',
     &         '        [default="split"] )'
      print *, '    -root   <path> (path is root path, use for all ',
     &         'directories)'
      print *, '         (default ='//root//')'
      print *, '    -rtsplit <path> (root path for split files)'
      print *, '    -rtproci <path> (root path for input process files',
     &         ' -- i.e., calibrated data)'
      print *, '    -rtproco <path> (root path for output process ',
     &         'files)'
      stop
6       continue  
!     continue on with main code:
