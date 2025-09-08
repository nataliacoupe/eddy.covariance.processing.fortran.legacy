!===================================================================================================
!                                    program lba: lpaths6.for
!===================================================================================================
!	path definitions (depends, in part, if os is windows or unix).
!	defines root pathnames rootsplit, rootprocin, rootprocout and data directories, for example:
!	data dirsplit/'split'/      <-- input from *.20? files
!	data dirprocin/'process'/   <-- input from cal/zero files *.zer, *.cal
!	data dirprocout/'process'/  <-- output for processed data
!---------------------------------------------------------------------------------------------------

!     parameters characterizing datafile structure:

!---- 1. root directories:--------------------------------------------------------------------------
!	global default root directory and data-directory names.
!	root directory for all is the default root (but each can be modified by command line).

	character*1 slash
	parameter(slash='/')
	character*8 split, process
	parameter(split='split', process='process')
	character*43 root
	parameter(root='/media/ncoupe/seagate/Data/tapajos/km67.H2O/')
	character*58 rootinp, rootsplit, rootprocin, rootprocout
!	character*29 root
!	parameter(root='/home/ncoupe/tapajos/km67.H2O/')                  
!	character*44 rootinp, rootsplit, rootprocin, rootprocout
	
!---- 2. data directories (platform-independent):---------------------------------------------------

!	character*44 dirsplit, dirprocin, dirprocout
	character*58 dirsplit, dirprocin, dirprocout

! 	names of directories:
!	raw data input (200,201 files); cal data (*.cal, zer), and output files, e.g.:
!	split input file (if name of file is 'infile'):
!	fnamein = rootsplit//slash//yymmdd//slash//dirsplit(1:ldir)//slash//infile(1:8)//'.200'  

	common /dirnames/ rootinp, rootsplit, rootprocin, rootprocout, dirsplit, dirprocin, dirprocout
