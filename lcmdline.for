!===================================================================================================
!                                         lba: lcmdline.for
!===================================================================================================
!     the purpose of this file is portability. we want to define two functions to facilitate getting
!     the command line into an appropriate variable within fortran programs.  this routine is needed
!     because of the differences between the watcom and the gnu fortran platforms. this program will
!     be an included module within the lsubs.for module which itself is included in lfluxes.for.
!
!     each function is defined twice. first it is defined for watcom, and then for gnu fortran. 
!     the definition for the platform on which this file is running will have the appropriate code 
!     available for reading. the other code is intended to be commented out. 
!            
!     important: when porting the code, simply change the commented lines.      
!     ********************************************************************
!
!     note: 
!     originally made to define subroutine 'readcmdline', which calls platform-specific functions to
!     access the command line, should now be used for any platform-specific function.
!
!     general strategy:
!     use this file to define a generalized function that, depending on the platform (linux or 
!     windows), calls the specific function to achieve that function. this file is thus the only one
!     that needs to be different when moving from one platform to another.  
!
!     written by daniel m. matross and scott r. saleska
!
!     history:  
!     09-oct-01:  modified by srs to also include subroutine 'getdatetime': calls platform-specific 
!                 functions to get date and time.
!     23-oct-01:  remove debug line to print out command-line args in readcmdline
!     10-may-02:  changed readcmdline to allow arguments longer (up to 16) than 8 characters.
!     09-sep-04:  changed readcmdline to allow arguments to 32 chars(also required changes length of
!                 directory path names in lpaths.for and to subroutines mygetarg and getfullpath in 
!                 lsubs.for).
!---------------------------------------------------------------------------------------------------


!===================================================================================================
!                                     subroutine: readcmdline
!===================================================================================================
!     usage:  call readcmdline(cmdline, cmdlen)
!     subroutine that returns a character string that contains the command line. other functions in 
!     'lsubs' parse the command line based on spaces, so all spaces must be included. 
!     in watcom, this is simply a wrapper function, because the fgetcmd function is built in. 
!     however, in gnu fortran we must build the string command line from each of the elements.
!---------------------------------------------------------------------------------------------------

!-----1) watcom version-----------------------------------------------------------------------------
!     subroutine readcmdline(cmdline,cmdlen)
!     character*80 cmdline
!     integer*4 cmdlen,fgetcmd
!     cmdlen=fgetcmd(cmdline)
!     end
      
!-----2) gnu version -------------------------------------------------------------------------------
      subroutine readcmdline(cmdline,cmdlen)
      
      character*100 cmdline
      character*33 arg
      integer*4 cmdlen, lencl, lenarg

!     gets first argument of command line:
      if (iargc() .ge. 1) then   
      call getarg(1,arg)
      cmdline = arg                 ! initializes command line with first argument
      lencl = lentrim0(cmdline, len(cmdline)) 
!     initializes length of cmd line with length of first argument.
      else 
      cmdline = char(0)             ! else, there is no command line.
      cmdlen = 0
      return
      end if
!     gets remaining arguments: 
      do i=2,iargc()
      call getarg(i,arg)
      lenarg = lentrim0(arg,len(arg))
      cmdline=cmdline(1:lencl)//' '//arg(1:lenarg)
!     accumulates:  arg1 + ' ' + arg2 + ' '...
      lencl = lencl + 1 + lenarg               ! keeps track of lenght of cmdline, including spaces.
      end do
      cmdlen=lencl
      end
!===================================================================================================
!                                end of subroutine readcmdline
!===================================================================================================      



!===================================================================================================
!                                   subroutine: getdatetime
!===================================================================================================      
!     gets date and time (duh!)
!---------------------------------------------------------------------------------------------------

!------1) watcom version----------------------------------------------------------------------------
!
!        subroutine getdatetime(iyear,imonth, iday, ihrs,imin,isec,platform)
!        integer*4 iyear,imonth, iday, ihrs,imin, isec, ihsec
!        character*6 platform
!        platform = 'watcom'
!        call getdat(iyear, imonth, iday)
!        call gettim(ihrs, imin, isec, ihsec)  
!        end

!------2) gnu version -------------------------------------------------------------------------------

      subroutine getdatetime(iyear,imonth,iday,ihrs,imin,isec,platform)
      integer*4 iyear,imonth, iday, ihrs,imin, isec, ihsec
      character*6 platform
      character*8 date
      character*9 time
      character*5 zone
      integer*4 ivalues(8)
      platform = 'gnufor'
      call date_and_time(date, time, zone, ivalues)
      iyear  = 1000*(ichar(date(1:1))-48) + 100*(ichar(date(2:2))-48)
     &        + 10*(ichar(date(3:3))-48)  + (ichar(date(4:4))-48) 
      imonth = 10*(ichar(date(5:5))-48)   + ichar(date(6:6))-48
      iday   = 10*(ichar(date(7:7))-48)   + ichar(date(8:8))-48 
      ihrs   = 10*(ichar(time(1:1))-48)   + ichar(time(2:2))-48
      imin   = 10*(ichar(time(3:3))-48)   + ichar(time(4:4))-48
      isec   = 10*(ichar(time(5:5))-48)   + ichar(time(6:6))-48
      end

!===================================================================================================
!                                  end of subroutine getdatetime
!===================================================================================================    
