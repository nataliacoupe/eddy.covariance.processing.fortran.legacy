!===================================================================================================
!					       Program lba: lfluxsub.for
!===================================================================================================
!
!     Implemented in 05-Sep-2001.
!     Mombines lmirror and lfilter files written by D. Matross (these are only used by lflux module,
!     don't need to have separate files).
!     From file: lmirror.for
!     Written by: Daniel M. Matross
!     Date: 7/9/2001
!
!     Purpose is to introduce the chilled mirror data into lfluxes. This file is intended as a 
!     module to include all the necessary code to read the chilled mirror data.
!     Modified to include the status bit functions as well, since they are specific to the lfluxes 
!     code.
!---------------------------------------------------------------------------------------------------

!===================================================================================================
!                                      Subroutine skipmirhead
!===================================================================================================
!     Usage: skipmirhead (88)
!     Short psubroutine to read through the header lines of the .201 files that contain the chilled
!     mirror data.
!---------------------------------------------------------------------------------------------------
      subroutine skipmirhead(infp)
      integer i,infp
      do i=1,4
      read(infp, *)
      end do       
      end
!===================================================================================================
!                                   End of subroutine skipmirhead
!===================================================================================================

!===================================================================================================
!                                      Subroutine getmirfacs
!===================================================================================================
!     Usage: getconvfacs()
! 
!     In order to work, a common array of dimensions 2x2x2 must be called in the program unit where 
!     this procedure is called. Procedure opens the file where "lba mircon" factors are stored and 
!     retrieves the ones relavant to the chilled mirror data. Factors are stored in a block of 
!     common memory.
!
!     Factors are entered into a 2x2x2 common mircon array of form (a,b,c).
!     Locations correspond to:	a --> ed1 {1} or ed2 {2}
!      					b --> dew point {1} or ambient temp {2}
!      					c --> multiplier {1} or offset {2}
!---------------------------------------------------------------------------------------------------
      subroutine getmirfacs()
      include 'lparams.for'
      real*8 mircon(2,2,2)
      character*9 name
      common /mirconv/ mircon
      name='lconv.inp'
      lrt = lentrim0(rootinp, len(rootinp))
      open(81,file=rootinp(1:lrt)//slash//name,status='old')
      read(81,*)
      read(81,*)
      read(81,*)
      read(81,*)
      read(81,801) mircon(1,1,1), mircon(1,1,2)
801   format(//f8.1,f8.1)
      read(81,802) mircon(1,2,1), mircon(1,2,2)
802   format(f8.1,f8.1)
      read(81,803) mircon(2,1,1), mircon(2,1,2)
803   format(///f8.1,f8.1)
      read(81,804) mircon(2,2,1), mircon(2,2,2)
804   format(f8.1,f8.1)
      close(81)                
      end
!===================================================================================================
!                                   End of subroutine getmirfacs
!===================================================================================================



!===================================================================================================
!                                        Subroutine rdmirtm
!===================================================================================================
!     Subroutine to read a time within the .201 file, for chilled mirror data.
!	If end of the file is reached, returns -999.99.
!---------------------------------------------------------------------------------------------------
      subroutine rdmirtm(mirtime,infp)
      integer infp
      real*8 mirtime
      real*8 tod,day
      read(infp,811,end=812) day,tod
811   format(f4.0,f9.3)
      mirtime=day+tod/86400.0
      return
812   mirtime= -999.99       
      end
!===================================================================================================
!                                    end of subroutine rdmirtm
!===================================================================================================



!===================================================================================================
!                                      subroutine matchmirtm
!===================================================================================================
!      Usage: matchmirtm(time, infp)
!      Subroutine to match the time in the mirror file with the designated start time.
!      03-sept-2001:  modified to allow reading other vars in c.201 file
!---------------------------------------------------------------------------------------------------
      subroutine matchmirtm(tstart,infp)
      integer infp
      real*8 mirtime
      real*8 tstart 
      call rdmirtm(mirtime,infp)  
      do while (mirtime.lt.tstart)               
      read (infp,*)
      read (infp,*)
      call rdmirtm(mirtime,infp)
      end do
      backspace infp
      end
!===================================================================================================
!                                    end of subroutine matchmirtm
!===================================================================================================



!===================================================================================================
!                                        function voltstodeg
!===================================================================================================
!     Usage: temp = voltstodeg(volts, convmultiplier, conv)
!     Subroutine to convert volts into degrees using linear interpolation.
!     One must input the multiplier and offset (slope and y-intercept) of the conversion factor.
!     Conversion parameters are retrieved with getmirfacs(). For use with chilled mirror data.
!---------------------------------------------------------------------------------------------------
      function voltstodeg(volts, mult, offset)
      real*8 volts,mult,offset
      real*8 voltstodeg
      voltstodeg=mult*volts+offset
      end
!===================================================================================================
!                                     end of function voltstodeg
!===================================================================================================



!===================================================================================================
!                                     subroutine statusbitarray
!===================================================================================================
!     a quick hexadecimal decoding is performed to disregard bits 0=11,leaving only the last 4 bits,
!     which are the input into the diagnstic word.
!---------------------------------------------------------------------------------------------------
      subroutine statusbitarray(status)
      integer status,diag,bit15,bit14,bit13,bit12
      integer diagarray(4)
      common /statusblock/ diagarray
      diag=status/4096
      bit15=diag/8
      bit14=(mod(diag,8))/4
      bit13=(mod(mod(diag,8),4))/2
      bit12= mod(mod(mod(diag,8),4),2)    
      diagarray(1)=bit15
      diagarray(2)=bit14
      diagarray(3)=bit13
      diagarray(4)=bit12
	end
!===================================================================================================
!                                  end of subroutine statusbitarray
!===================================================================================================



!===================================================================================================
!                                     subroutine pseudofilter
!===================================================================================================
!     usage: call pseudofilter(filtvec, filtset(5),enum,n,tstart)
!     this is a function that runs a pseudofilter.
!	it is framed in such a way as to allow multiple types of pseudofilters.	for now, it can only 
!     call one.
!---------------------------------------------------------------------------------------------------
      subroutine pseudofilter(filtvec, type, yymmdd, infp, ndata,tstart)
      integer type, eddynum,infp
      real*8 tstart
      integer filtvec(ndata)
      character*6 yymmdd
      if (type.eq.1) then
      call filtoppnolock(filtvec,infp,yymmdd,ndata,tstart)
      endif
      end
!===================================================================================================
!                                  end of subroutine pseudofilter
!===================================================================================================



!===================================================================================================
!                                     subroutine filtoppnolock
!===================================================================================================
!     usage: call filtoppnolock(vec,1,14400,time)
!     sets a pseudofilter to filter on the no signal lock flag from the opposite eddy system. does 
!     so by opening opposite eddy .200 file and retrieving the status word, then converting it and 
!     returning the control vector as that which has bit14 (the no signal lock bit) flagged.
!---------------------------------------------------------------------------------------------------
!     .200 has already been opened in main program
      subroutine filtoppnolock(filtvec,infp,yymmdd,ndata,tstart)    
      integer filtvec(ndata),filtdiagarray(4)
      integer i,status,iflag,infp
      real*8 tstart
      character*6 yymmdd
      call match200time(infp, tstart) 
      do i=1,ndata
      read(infp,555,end=550) status     
555   format(t64,i6)
      call getstatflag(iflag,status,2)
      filtvec(i)=iflag
      end do                       
550   continue
      end
!===================================================================================================
!                                    end of subroutine filtoppnolock
!===================================================================================================



!===================================================================================================
!                                    subroutine skip200header(infp)
!===================================================================================================
!     subroutine to skip header lines of a 200 file.
!---------------------------------------------------------------------------------------------------
      subroutine skip200header(infp)
      integer infp
      read(infp,*)
      read(infp,*)
      end
!===================================================================================================
!                                end of subroutine skip200header(infp)
!===================================================================================================



!===================================================================================================
!                                     subroutine match200time
!===================================================================================================
!     matches time in the .200 file given by the pointer to tstart given
!---------------------------------------------------------------------------------------------------
	subroutine match200time(infp,tstart)
	real*8 time, tstart
	integer infp
	call read200time(time,infp)
	backspace infp
	do while (time.lt.tstart)               
	call read200time(time,infp)
	end do
	end
!===================================================================================================
!                                  end of subroutine match200time
!===================================================================================================



!===================================================================================================
!                                     function read200time
!===================================================================================================
!     usage: time=read200time(infp)
!     with the cursor properly positioned, this function reads the time and day from a .200 and 
!     returns it as a double with the decimal being the fraction of a day.
!---------------------------------------------------------------------------------------------------
      subroutine read200time(time,infp)
	integer infp
      double precision day, tod, time
      read(infp,552,end=553) day,tod
552   format(f4.0,f9.3)
	time=day+tod/86400.0
      return
553	time= -999.99       
	end
!===================================================================================================
!                                   end of function read200time
!===================================================================================================



!===================================================================================================
!                                     subroutine open200file
!===================================================================================================
!     usage: call open200file(501, yymmdd, eddynum)
!     a subroutine to open the .200 file with the given infile pointer. takes care of names.
!---------------------------------------------------------------------------------------------------
      subroutine open200file(infp,yymmdd,eddynum)
      include 'lpaths.for'
      integer infp, eddynum, lenfn
      character*48 filename
      character*2  edaorb
      character*6  yymmdd
      if (eddynum.eq.1) then
      edaorb = 'EA'
      else 
      edaorb  = 'EB'
      endif
      filename = edaorb//yymmdd//'.200'
      call getfullpath( rootsplit,yymmdd, dirsplit, filename, lenfn) 
      open(infp,file=filename(1:lenfn),status='old')
      return
      end
!===================================================================================================
!                                    end of subroutine open200file
!===================================================================================================



!===================================================================================================
!                                        function getstatflag
!===================================================================================================
!     usage: nolock = getstatflag(4094,2)
!     function to retrieve the integer (0 or 1) value of the status bit that has been decoded. 
!     since we're only interested in the status words, we only need the last four bits 
!     (bits 15,14,13,12) from the status word, thus we use the system:
!     1: bit 15 -- difference in speed of sound
!     2: bit 14 -- no signal lock
!     3: bit 13 -- sonic amplitude too low
!     4: bit 12 -- sonic amplitude too high
!     function relies on integer division.
!     note: this is logically a function, but my watcom compiler won't let work for me as a function
!     {error cp-01}. i'm convinced something is wrong with the compiler. however, i redid it as a 
!     subroutine. a patch, but i'm rather frustrated with fortran in these matters.
!---------------------------------------------------------------------------------------------------
      subroutine getstatflag(iflag,status,whichflag)
      integer status,whichflag,diag,bit15,bit14,bit13,bit12
      integer flagarray(4)
!     we do a quick hex decode to disregard bits 0-11, leaving us only with the last 4 bits, which 
!     are input into the diagonistic word.
      diag=status/4096
      bit15=diag/8
      bit14=(mod(diag,8))/4
      bit13=(mod(mod(diag,8),4))/2
      bit12=mod(mod(mod(diag,8),4),2)
      flagarray(1)=bit15
      flagarray(2)=bit14
      flagarray(3)=bit13
      flagarray(4)=bit12
      if (whichflag.lt.1 .or. whichflag.gt.4) then
      iflag = 0
      else 
      iflag = flagarray(whichflag)        
      endif
      end   
!===================================================================================================
!                                   end of function getstatflag
!=================================================================================================== 
