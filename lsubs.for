!===================================================================================================
!                                           lsubs.for
!
!      core subroutines and functions used in lba program modules.
!      lcalpf, lcaled, lfluxes, laggred
!===================================================================================================

      include 'lcmdline.for'          ! subroutine def for getting command line (platform-dependent)


!===================================================================================================
!                                      subroutine getfullpath                                      !
!===================================================================================================
!     puts together the pieces to make a full path for file 'fname'.
!     returns  fname = root + yymmdd + datadir + fname
!              lfn   = trimmed lenght of filename
!---------------------------------------------------------------------------------------------------
      subroutine getfullpath( rootpath, yymmdd, datadir, fname, lfn) 
      include 'lpaths.for'
!      character*34 rootpath, datadir   ! lengths consistent with defs of corresp. vars 
      character*54 rootpath, datadir   ! lengths consistent with defs of corresp. vars 
      character*6 yymmdd               ! in lparams, lpaths
      character*96 fname
      character*8 basename
      integer ldir,lrt,lfn
      basename = fname(1:8)
      ldir = lentrim0(datadir, len(datadir))  ! default = split
      lrt = lentrim0(rootpath, len(rootpath))
      fname = rootpath(1:lrt)//slash//yymmdd//slash//datadir(1:ldir)
     &        //slash//basename
      lfn = lentrim0(fname, len(fname))    ! length of fname
      return
      end
!===================================================================================================
!                                  end of subroutine getfullpath                                   !
!===================================================================================================


!===================================================================================================
!                                      subroutine readsysd
!===================================================================================================
!     reads system data parameters from text file lsysdat.inp
!	if nread =      1:  read eddy 1, then exit
!                       2:  read eddy 2 data,  exit
!                       3:  read profile data, and exit
!---------------------------------------------------------------------------------------------------
      subroutine readsysd(nread)
      integer nread
      include 'lparams.for'
      lrt = lentrim0(rootinp, len(rootinp))
      print *, '    file "'//rootinp(1:lrt)//slash//
     &         'lsysdat.inp": read generic parameters...'
      open(10, file=rootinp(1:lrt)//slash//'lsysdat.inp',status='old')
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*) cref(1), cref(2), cref(3), csurv
      read(10,*) !~475   ~400    ~325    ~376  -->  cal-gas concentrations
      read(10,*)
      read(10,*) tflstep, tavecal, tdelayed, tavepf, tdelaypf
      read(10,*)                        ! 120    ed: 60   20?  PF: 60    50
      read(10,*)                        ! tflstep = time at each cal-gas/profile level  (all in secs)
      read(10,*) flsamed, flsampf, flcaled, flcalpf  ! use flsam, flcal here
      read(10,*)                        ! 8000           1000
      nflstep  = tflstep*freqhi         ! number of points for each calibration step (hi-frequency)
      navecal  = tavecal*freqhi         ! number of points to average for a calibration
      ndelayed = tdelayed*freqhi        ! number of points to wait before averaging
      if ((ndelayed+navecal) .ge. nflstep) then
      print *,' Averaging time too long. '
      stop 
      endif
!     low-frequency parameters:
      nflsteplo = tflstep*freqlo         ! number of points for each calibration step (low-freq)
      navecalo  = tavecal*freqlo         ! number of points to average for slow (201) eddy cal data 
      navepf    = tavepf*freqlo          ! number of points for one profile level average
      ndelaypf  = tdelaypf*freqlo
!     converts times to fraction of day (high-precision variables):
!     sampling at 24*60*60
      tflstep   = tflstep/86400.d0       ! eddy
      tavecal   = tavecal/86400.d0  
      tdelayed  = tdelayed/86400.d0  
      tavepf    = tavepf/86400.d0        ! profile 
!     overlaps to allow between half-hour segments because of zeros:
      tovlap    = (tflstep + tdelayed)
      novlap    = int(tovlap*86400.*freqhi)    
      if(nread.eq.1) then   
      goto 10                        ! goes directly to read eddy 1
      elseif(nread.eq.2) then 
      do i=1,19                      ! skips over eddy system 1
      read(10,*)
      enddo
      goto 20
      elseif(nread.eq.3) then
      do i=1,19*2                    ! skips over eddy system 1 & 2
      read(10,*)
      enddo
      goto 30
      else
      print *,'Bad parameter "nread" for subroutine readsysd.'
      goto 40
      endif 
!     eddy system 1:  
10    continue
      print *, '    file "'//rootinp(1:lrt)//slash//
     &         'lsysdat.inp": read eddy 1 parameters...'
      read(10,*)
      read(10,*) 
      read(10,*) zflsam, zflcal, zpcon, zpcell  ! only zpcell used
      read(10,*) 
!     where:
!     zflsam = sample flowmeter zero (v)
!     zflcal = calibration flowmeter zero (v)
!     zpcon  = pressure controller zero (v)
!     zpcell = cell pressure zero (v) this changed in 20210401 as flow controler was changed
      read(10,*)
      read(10,*) gflsam,gflcal,gpcon,gpcell,gtcell 
!     gains for flow, pressure and temperature        
      gpcell = gpcell * 760/14.7         ! converts from psi to torr
      read(10,*) 
      read(10,*)
      read(10,*) clicor,t0co2, kco2      ! licor factory coefs, co2
      read(10,*)                         ! x5
      read(10,*)
      read(10,*) hlicor,t0h2o, kh2o      ! licor factory coefs, h2o
      read(10,*)                         ! x3
      read(10,*)  
      read(10,*) tlagco2, tau, tlagh2o, tauh  ! lag between licor and sonic, smearing coefs.

      read(10,*) 
      read(10,*) 
      read(10,*) gu, gv, gw, gt, dirofu  ! g* = gain of u,v,w,t sonic data
      read(10,*)                         ! dirofu = u direction, degrees from north clockwise
      nlagco2 = int(tlagco2*freqhi)      ! number of samples in lag
      nlagh2o = int(tlagh2o*freqhi)
      goto 40      
!     eddy system 2:
20    continue
      print *, '    file "'//rootinp(1:lrt)//slash//
     &         'lsysdat.inp": read eddy 2 parameters...'
      read(10,*)
      read(10,*)
      read(10,*) zflsam, zflcal, zpcon, zpcell  ! not used
      read(10,*)
      read(10,*)
!     gains for flow, pressure and temperature  
      read(10,*) gflsam,gflcal,gpcon,gpcell,gtcell      
      gpcell = gpcell * 760/14.7          ! convert from psi to torr
      read(10,*) 
      read(10,*)
      read(10,*) clicor,t0co2, kco2       ! licor factory coefs, co2
      read(10,*)                          ! x5
      read(10,*)
      read(10,*) hlicor,t0h2o, kh2o       ! licor factory coefs, h2o
      read(10,*)                          ! x3
      read(10,*)  
      read(10,*) tlagco2, tau, tlagh2o, tauh  ! lag between licor and sonic, smearing coefs
      read(10,*) 
      read(10,*) 
      read(10,*) gu, gv, gw, gt, dirofu   ! g* = gain of u,v,w,t sonic data
      read(10,*)                          ! dirofu = u direction, degrees from north clockwise
      nlagco2 = int(tlagco2*freqhi)       ! number of samples in lag
      nlagh2o = int(tlagh2o*freqhi)
      goto 40
!     profile system:
30    continue
      print *, '    file "'//rootinp(1:lrt)//slash//
     &         'lsysdat.inp": read profile parameters...'
      read(10,*)
      read(10,*)
      read(10,*) zflsam, zflcal, zpcon, zpcell  
      read(10,*)
      read(10,*)
!     gains for flow, pressure and temperature
      read(10,*) gflsam,gflcal,gpcon,gpcell,gtcell 
      gpcell = gpcell * 760/14.7                   ! convert from psi to torr
      read(10,*) 
      read(10,*)
      read(10,*) clicor,t0co2, kco2                ! licor factory coefs, co2 (x5)
      read(10,*)                                   
      read(10,*)
      read(10,*) hlicor,t0h2o, kh2o                ! licor factory coefs, h2o (x3)
                                                   
40    continue
      close(10)
      return
      end
!===================================================================================================
!                                    end of subroutine readsysd                                    !
!===================================================================================================



!===================================================================================================
!                                        subroutine readconv                                       !
!===================================================================================================
!     reads conversion data from text file lconv.inp
!     if nread = 1: reads eddy 1 convertion factors, then exit.
!     if nread = 2: reads eddy 2 data, then exit.
!     if nread = 3: reads profile data, then exit.
!     if nread = 4: reads ground system.					
!---------------------------------------------------------------------------------------------------
      subroutine readconv(nread)
      integer nread
      include 'lparams.for' ! among other things, defines variables in common-block 'sysvars', which
                            ! are read in by this subroutine. includes the conversion coefficient
                            ! arrays 'cc' and 'offset' as well as the system file parameters
                            !   (e.g. gpcon, zpcon, etc) input by readsysd from lsysdat.inp

      lrt = lentrim0(rootinp, len(rootinp))
      open(10, file=rootinp(1:lrt)//slash//'lconv.inp',status='old')

      do i=1,4
      read(10,*)  ! get header
      enddo
!     201 line (for eddies or profile):
!     get initial set from values read-in from lsysdat.inp:
!     (will be eddy 1, 2, or profile, depending on system)
      cc(1) = gpcon                   ! pressure controller
      offset(1) = -zpcon*cc(1)        ! (zero to subtract)
      cc(2) = gflsam                  ! sample-flow
      offset(2) = -zflsam*cc(2) 
      cc(3) = gflcal                  ! cal-flow
      offset(3) = -zflcal*cc(3) 
      cc(4) = gpcell                  ! licor cell pressure
      offset(4) = -zpcell*cc(4)
      cc(5) = 1                       ! ready-light
      cc(6) = gtcell                  ! licor cell temperature
!     get remaining 201 for eddy 1:
      do i=6+1, 10  ! finish 201 line
      read(10,*) cc(i),offset(i)
      end do
!     if system 2 or greater, override eddy 1 with eddy 2 
      if(nread.ge.2) then
      read(10,*)   ! skip header for eddy 2
      do i=6+1, 10  ! overwrite with values for eddy 2
      read(10,*) cc(i),offset(i)
      end do
      endif
!     if system 3 or greater, override the eddy 2 values in favor of profile:
      if(nread.ge.3) then 
      read(10,*)      ! skip header for profile
      do i=6+1, n1PF  ! finish profile 201 line
      read(10,*) cc(i),offset(i)
      end do
!     profile 202 line:
      read(10,*)  ! skip 202 header
      do i=n1PF+1, n12PF
      read(10,*) cc(i),offset(i)
      end do
!     profile 203 line:
      read(10,*)  ! skip 203 header
      do i=n12PF+1, n13PF
      read(10,*) cc(i),offset(i)
      end do
      endif
!     if system =4 (ground), override the profile values in favor of ground:
      if(nread.ge.4) then
!     ground 201 line:
      read(10,*)  ! skip 201 header (nb. start over with i=1)
      do i=1, n1GD  
      read(10,*) cc(i),offset(i)
      end do
!     ground 202 line:
      read(10,*)  ! 202 header
      do i=n1GD+1, n12GD
      read(10,*) cc(i),offset(i)
      end do
!     ground 203 line:
      read(10,*)  ! 203 header
      do i=n12GD+1, n13GD
      read(10,*) cc(i),offset(i)
      end do
      do i=n13GD+1, n13GD+10  ! last few to do status bits
      cc(i)=1.
      offset(i)=0.
      enddo
      endif
      if(nread.ge.5) then
      print *,'Bad parameter "nread" for subroutine readconv.'
      endif
!     file all read
      close(10)
      return
      end
!===================================================================================================
!                                    end of subroutine readconv                                    !
!===================================================================================================



!===================================================================================================
!                                       subroutine startime                                        !
!===================================================================================================
!    gets desired start-time based on the initial time. reads from the first line of the input file.
!---------------------------------------------------------------------------------------------------
      subroutine startime(firstime, nextime, nextmin)
      real*8 firstime                      ! initial time in decimal julian days
      real*8 rnumsegs,fracseg              ! function value:  time to start with
      real*8 eps, nextmin, segspday	
!     where:
!     nextmin = next minute to start at.   60  --> starts on the hour
!                                          30  --> starts on the 1/2 hour
!                                          1/2 --> starts on the 30 seconds
      real*8 nextime
!     parameters:
      eps = 2.d0                                             ! epsilon = error allowed (in seconds).
      eps= (eps/60.)/nextmin                                 ! error allowed (fraction of 'nextmin' 
                                                             ! segment).
      segspday = 1440./nextmin                               ! number of 'nextmin' segments per day 1440=60*24 sec/day
      rnumsegs = (firstime - int(firstime))*dble(segspday)   ! number of (real) 'nextmin'-minute 
                                                             ! segments so far in this day.
      fracseg= rnumsegs-int(rnumsegs)                        ! fraction of this segment we are into
      if(fracseg.gt.eps) then                                ! sets 'startime' to start of next 
                                                             ! segment.
      nextime= firstime + (1-fracseg)/dble(segspday)
      else 
      nextime = firstime                                     ! keeps startime at initial time. 
      endif
      return
      end
!===================================================================================================
!                                    end of subroutine startime                                    !
!===================================================================================================



!===================================================================================================
!                                           funtion iargf                                          !
!===================================================================================================
      function iargf(cmdlin)
      character*100 cmdlin
      character*1 oldarg,newarg
      iargf = 0
      oldarg = ' '
      do 20 i = 1, 100
      newarg = cmdlin(i:i)
      if (newarg .ne. ' ') newarg = 'y'
      if (newarg .eq. oldarg) goto 20
      if (newarg .eq. ' ') then
      oldarg = ' '
      else
      oldarg = 'y'
      iargf = iargf + 1
      end if
20    continue
      return
      end
!===================================================================================================
!                                       end of funtion iargf                                       !
!===================================================================================================



!===================================================================================================
!                                       subroutine mygetarg                                        !
!===================================================================================================
      subroutine mygetarg(cmdlin,itharg,argi)
      character*100 cmdlin
      character*32 argi
      character*1 oldarg,newarg
      argi = '                ' 
      iarg = 0
      oldarg = ' '
      do 20 i = 1, 100
      newarg = cmdlin(i:i)
      if (newarg .ne. ' ') newarg = 'y'
      if (newarg .eq. oldarg) goto 20
      if (newarg .eq. ' ') then
      oldarg = ' '
      else
      oldarg = 'y'
      iarg = iarg + 1
      j = i
      if (iarg .eq. itharg) goto 40
      end if
20     continue
      goto 80
40     argi(1:1) = cmdlin(j:j)
      do 60 i = j+1, j+len(argi)-1
      newarg = cmdlin(i:i)
      if (newarg .eq. ' ') goto 80
      argi((i-j+1):(i-j+1)) = cmdlin(i:i)
60    continue
80    return
      end

!===================================================================================================
!                                   end of subroutine mygetarg                                     !
!===================================================================================================



!===================================================================================================
!                                        function lentrim0                                         !
!===================================================================================================
!      note: same as lentrim in watcom, or len_trim in gnu fortran.
!---------------------------------------------------------------------------------------------------
      function lentrim0(xx,nn)
      character*80 xx
      ii=1
      do i=1,nn
      if(xx(i:i).eq.' ' )then
      ii=i-1
      go to 3
      endif
      enddo
      ii=nn
3     lentrim0=ii
      return
      end
!===================================================================================================
!                                    end of function lentrim0                                      !
!===================================================================================================



!===================================================================================================
!                                        subroutine infiles                                        !
!===================================================================================================
!     accepts ymd ('yymmdd') and types (some subset of '123456') and returns 'names(1:num)': a set 
!     of corresponding file names, where num = number of file names returned:
!     types = 1 --> yymmddea; 
!             2 --> yymmddeb; 
!             3 --> yymmddpf;
!             4 --> yymmddgd;
!             5 --> yymmddsn;
!             6 --> ntyymmdd.
!---------------------------------------------------------------------------------------------------
      subroutine infiles(ymd,types,names,num)
      character*6 ymd, types
      character*12 names(6)
      character*12 fchars
      integer fnum, pos
      fchars='EAEBPFGDSNNT'
!     where the substrings are:
!     EA = eddy1, EB=eddy2, PF=profile, GD=ground, sn=backup sonic, nt = natal system.
      num= lentrim0(types, len(types))
!     where:
!     types = some subset of '123456'
!	need to add error checking (make sure no duplicates, etc) 
      do i=1,num
      fnum = ichar(types(i:i))-48           ! converts numbers as characters to integer
      pos = (fnum-1)*2 + 1                  ! position within fchars
      names(i)=fchars( pos:(pos+1) )//ymd   ! take 2 chars of fchars
      enddo
      return
      end
!===================================================================================================
!                                   end of subroutine infiles                                      !
!===================================================================================================


! note: there are two subroutines to calculate concentrations from voltages: getconc (profile dilution  
! & air measurements and eddy calibration) and getconc2 (eddy).

!===================================================================================================
!                                        subroutine getconc                                        !
!===========================================================================================
!     calculates co2 and h2o concentrations from voltages and calibration coefficients.
!     input args: h2o = voltage of h2o signal
!                 co2 = voltage of co2 signal
!                 pres = licor cell pressure (torr) of measurement
!                 po   = pressure during cal that gave cal coeffs
!                 temp = licor cell temperature (c) of measurement
!                 to   = temperature during cal that gave cal coefs
!                 coefs(1:3) = interpolated cal coefficients for licor co2 polynomial
!     input in common:        hlicor = factory cal coefs for h2o
!                             t0h2o  = factory cal temp for h2o   
!     output as args:         h2o    = mole fraction of water vapor (in mmol/mol)
!                             co2    = co2 concentration (in ppm)
!-------------------------------------------------------------------------------------------
      subroutine getconc(h2o, co2, pres, po, temp, to, coefs)
          ! calculate co2, h2o concentrations from voltages & cal coeffs
          ! use for dilution and other cal values not ambient air concentrations
          !   Input arguments:
          !     h2o = voltage of h2o signal
          !     co2 = voltage of co2 signal
          !     pres = Licor cell pressure (torr) of measurement
          !     po   = pressure during cal that gave cal coeffs
          !     temp = Licor cell temperature (C) of measurement
          !     to   = temperature during cal that gave cal coefs
          !     coefs(1:3) = interpolated calibration coefficients 
          !              for Licor CO2 polynomial
          !   Input in common:
          !     hLicor = factory cal coefs for h2o
          !     t0h2o  = factory cal temp for h2o
          !     
          !   Output as arguments:
          !     h2o = mole fraction water vapor (in mmol/mol)
          !     co2 = CO2 conc (in ppm)
      include 'lparams.for'             ! among other things, defines vars in common block 
                                        ! 'sysvars', including:
                                        ! real*4 clicor(5),hlicor(3),t0co2,t0h2o
                                        ! (used in calculating co2, water concentrations

!     dummy arguments:
      real*8 h2o, co2, pres, po, temp, to
      real*8 coefs(3)
!     internal arguments:
      real*8 vx, vsq, vcb               ! 1st,2nd,3rd order of voltage
      real*8 pratio, tratio, wratio, wv, chi, fdilut
!      print*,'h2o: ',h20, '   co2: ', co2, '  pres: ', pres
!      print*,'po: ',po, '  temp: ', temp, ' to: ', to, '  t0h2o: ', t0h2o
!      print*,'coefs(1): ', coefs(1),' coefs(2)', coefs(2),' coefs(3)', coefs(3)
!      print*,'hlicor(1): ', hlicor(1),' hlicor(2)', hlicor(2),' hlicor(3)', hlicor(3)
!-------------------------------------------------------------------------------------------
!     water vapour (wv):  
! 	wv = fw[ vh2o * (po/p)^0.9 ] * (t/to)    (wv in mmol/mol).
! 	where fw[x] = a*x + b*x^2 + c*x^3        (according to licor-6262 manual, sect 3.4).
!     applies licor factory cal coefficients for water vapor and, hence, normalize to 
!     factory calibration conditions: 
!	po = 760 torr, to = t0h2o.
!-------------------------------------------------------------------------------------------
      pratio=1
      if( pres.gt.0.1 ) pratio= 760.d0/pres       ! po/p
      pratio= pratio**9.d-1                       ! (po/p)^0.9 [ water only, see p. 3-6 ]
      tratio=(273.d0+temp)/(273.d0+t0h2o)         ! t/to for water
      vx=h2o*pratio
      vsq=vx*vx
      vcb=vx*vsq
      wv=((hlicor(1)*vx) + (hlicor(2)*vsq) + (hlicor(3)*vcb))*tratio  ! mmol/mol
      wratio=18.d0/(29.d0*(1.d0-1.d-3*wv) + 1.8d-2*wv)  ! molec wt(h2o)/molec wt(moist air)
      if(nmux.eq.1) then
       write(99,*)'               water (wv) =' , wv
      endif
      h2o=wv                                      ! *wratio ! gh2o/kg air  ! leave as molar fraction
!-------------------------------------------------------------------------------------------
!     co2:  
!     [co2] = chi * fdilut * f[ vco2/chi * (po/p) ] * (t/to)
!     where f[x] = a*x + b*x^2 + c*x^3 is calibration polynomial for dry air.
!     applies field-calibration coefficients for co2 and, hence, normalize to field-cal 
!     conditions given by po and to arguments.
!     chi    = water pressure-broadening correction
!     fdilut = water-dilution correction   (see licor-6262 manual, sect. 3.5 and 3.6)
!-------------------------------------------------------------------------------------------
      chi    = 1.d0 + (5.d-1*wv/1000.d0)          ! pressure-broadening correction
      fdilut = 1.d0/(1.d0-wv/1000.d0)             ! correction for water dilution
      pratio = 1 
      if(pres.gt.0.1) pratio = po/pres
       tratio = (temp+273.d0)/(to+273.d0)         ! (nb: '+273' correction added 26-feb-01)

!        print*,'tratio: ', tratio
!        print*,'co2: ', co2
!        print*,'vx: ', vx
!        print*,'fdilut: ', fdilut
!        print*,'pratio: ', pratio
!        print*,'chi: ', chi
!        print*,'co2: ', co2
!        print*,'wv: ', wv

       vx  = co2*pratio/chi
       vsq = vx*vx
       vcb = vx*vsq
       co2 = (coefs(1)*vx)+(coefs(2)*vsq)+(coefs(3)*vcb) ! applies interpolated coefficents
       co2 = chi*fdilut*co2*tratio                ! applies remaining water corrections
       return
      end
!===========================================================================================
!                                   end of subroutine getconc  
!===========================================================================================



!===========================================================================================
!                                      subroutine getconc2       
!===========================================================================================
!     calculates co2 and h2o concentrations from voltages & cal coefficients *includes wvscale*
!     input args: h2o = voltage of h2o signal
!                 co2 = voltage of co2 signal
!                 pres = licor cell pressure (torr) of measurement
!                 po   = pressure during cal that gave cal coeffs
!                 temp = licor cell temperature (c) of measurement
!                 to   = temperature during cal that gave cal coefs
!                 coefs(1:3) = interpolated calibration coefficients for licor co2 polynomial
!                 wvscale = water vapor scale factor     		(added 040717 updated nrc 170125)
!	input in common:
!			hlicor = factory cal coefs for h2o
!			t0h2o  = factory cal temp for h2o
!	output as args:
!			h2o = mole fraction water vapor (in mmol/mol)
!			co2 = co2 conc (in ppm)
!     note in 09-dec-04:
!     make h2o output mixing ratio w.r.t. dry air (mmol h2o/mol dry air)  (it was per mole wet air)
!---------------------------------------------------------------------------------------------------
      subroutine getconc2(h2o, co2, pres, po, temp, to, coefs, wvscale)
      include 'lparams.for'       ! among other things, defines variables in common block 'sysvars', 
                                  ! including:
                                  ! real*4 clicor(5), hlicor(3), t0co2, t0h2o (used to calculate co2
                                  ! and water concentrations).
!     dummy arguments:
      real*8 h2o, co2, pres, po, temp, to
      real*8 coefs(3)
      real*8 wvscale
!     internal arguments:
      real*8 vx, vsq, vcb         ! 1st,2nd,3rd order of voltage
      real*8 pratio, tratio, wratio, wv, chi, fdilut
!---------------------------------------------------------------------------------------------------
!     water vapour (wv):  
! 	wv = fw[ vh2o * (po/p)^0.9 ] * (t/to)       (wv in mmol/mol).
! 	where fw[x] = a*x + b*x^2 + c*x^3           (according to licor-6262 manual, sect 3.4).
! 	applies licor factory calibration coefficients for water vapor and, hence, normalize to 
!       factory calibration conditions: 
!	po = 760 torr, to = t0h2o.
!---------------------------------------------------------------------------------------------------
      pratio=1
      if( pres.gt.0.1 ) pratio= 760.d0/pres                     ! po/p
      pratio = pratio**9.d-1                                    ! (po/p)^0.9 [water only, see p. 3-6]
      tratio = (273.d0+temp)/(273.d0+t0h2o)                     ! t/to for water
      vx  = h2o*pratio
      vsq = vx*vx
      vcb = vx*vsq
      wv  = ((hlicor(1)*vx) + (hlicor(2)*vsq) + (hlicor(3)*vcb))*tratio ! mmol/mol
!...................................................................................................
!     added in 17-jul-2004 (srs) modified 25-jan-17 (nrc)
!       print*, 'wvscale: ', wvscale, '  wv: ', wv
      wv = wv*wvscale
!       print*, 'wvscale: ', wvscale, '  wv: ', wv
      fdilut = 1.d0/(1.d0 - wv/1000.d0)                         ! correction for water dilution
!...................................................................................................
      wratio = 18.d0/(29.d0*(1.d0-1.d-3*wv) + 1.8d-2*wv)        ! molec wt(h2o) / molec wt(moist air)
      if(nmux.eq.1) then
       write(99,*)'               water (wv) =' , wv
      endif
      h2o = wv*fdilut                                           ! *wratio ! gh2o/kg air  ! leave as molar fraction
!---------------------------------------------------------------------------------------------------
!     co2:  
!     [co2] = chi * fdilut * f[ vco2/chi * (po/p) ] * (t/to)
!     where f[x] = a*x + b*x^2 + c*x^3 is calibration polynomial for dry air.
!     applies field-calibration coefficients for co2 and, hence, normalize to field-calibration 
!     conditions given by po and to arguments.
!     chi    = water pressure-broadening correction
!     fdilut = water-dilution correction   (see licor-6262 manual, sect. 3.5 and 3.6)
!---------------------------------------------------------------------------------------------------
      chi    = 1.d0 + (5.d-1*wv/1000.d0)                  ! pressure-broadening correction
      pratio = 1 
      if(pres.gt.0.1) pratio = po/pres
       tratio = (temp+273.d0)/(to+273.d0)                 ! (nb: '+273' correction added 26-feb-01)
       vx  = co2 * pratio/chi
       vsq = vx*vx
       vcb = vx*vsq
       co2 = (coefs(1)*vx)+(coefs(2)*vsq)+(coefs(3)*vcb)  ! apply interpolated coefficents
       co2 = chi*fdilut*co2*tratio                        ! apply remaining water corrections
       return
      end
!===================================================================================================
!                                  end of subroutine getconc2                                      !
!===================================================================================================



!===================================================================================================
!                                        function ihuntf                                           !
!===================================================================================================
!     extracted from "numerical recipes", edition of 1992.
!     given monotonic array xx(1:n) and a value x, returns value = jlo such that x is between 
!     xx(jlo) and xx(jlo+1).  
!     0,n = exceed range. jlo1 on input is the initial guess;  used for searches with rolling
!     best guess.
!---------------------------------------------------------------------------------------------------
      function ihuntf(xx,n,x,jlo1)
      integer jlo,n,inc,jhi,jm,ihuntf,jlo1
      real*8 x,xx(n)
      logical ascnd
      jlo=jlo1
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then 
!     go to bisection straight off
      jlo=0
      jhi=n+1
      goto 3
      end if
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then 
!     hunt up, increment inc
1     jhi=jlo+inc
      if(jhi.gt.n)then
      jhi=n+1
      else if(x.ge.xx(jhi).eqv.ascnd)then 
!     not done hunting
      jlo=jhi
      inc=inc+inc 
!     double increment
      goto 1
      end if
      else                           
!     hunt down
      jhi=jlo
2     jlo=jhi-inc
      if(jlo.lt.1)then
      jlo=0
      else if(x.lt.xx(jlo).eqv.ascnd)then
      jhi=jlo
      inc=inc+inc
      goto 2
      end if
      end if                          
!     end of hunt for bracketting indexes
3     if(jhi-jlo.eq.1)goto 30
      jm=(jhi+jlo)/2         
!     begin bisection
      if(x.gt.xx(jm).eqv.ascnd)then
      jlo=jm
      else
      jhi=jm
      end if
      goto 3
30      ihuntf=jlo
      return
      end
!===================================================================================================
!                                     end of function ihuntf                                       !
!===================================================================================================



!===================================================================================================
!                                       subroutine gaussj                                          !
!===================================================================================================
!     extracted from "numerical recipes", 2nd edition (p. 30)
!     (do-loop labels removed, replaced with 'end do')
!     solve a linear system: a x = b,  
!     by getting x = a^-1 b 
!     using the gauss-jordan elimination method.
!     inputs:
!     a(1:n, 1:n) = matrix in a np:np physical array.
!     b(i:n, 1:m) = m different right-hand n-vectors in a mp:mp array.
!     outputs:
!     a = the inverse of the input matrix.
!     b = the set of m solution vectors.
!---------------------------------------------------------------------------------------------------
      subroutine gaussj(a,n,np,b,m,mp)
      integer m,mp,n,np,nmax
      real*8 a(np,np),b(np,mp)                                          ! make double precision
      parameter (nmax=50)
      integer i,icol,irow,j,k,l,ll,indxc(nmax),indxr(nmax),ipiv(nmax)
      real*8 big,dum,pivinv
      do j=1,n
       ipiv(j)=0
      enddo
      do i=1,n
       big=0.
      do j=1,n
       if(ipiv(j).ne.1)then
        do k=1,n
         if (ipiv(k).eq.0) then
          if (abs(a(j,k)).ge.big)then
           big=abs(a(j,k))
           irow=j
           icol=k
          endif
         else if (ipiv(k).gt.1) then
          write(99,*)'    Error:  calibration matrix is singular.'
          pause '    Error:  calibration matrix is singular.'
         endif
        enddo
       endif
      enddo
      ipiv(icol)=ipiv(icol)+1
      if (irow.ne.icol) then
       do l=1,n
        dum=a(irow,l)
        a(irow,l)=a(icol,l)
        a(icol,l)=dum
       enddo
       do l=1,m
        dum=b(irow,l)
        b(irow,l)=b(icol,l)
        b(icol,l)=dum
       enddo
      endif
      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol).eq.0.) pause '  Singular matrix in gaussj.'
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.
      do l=1,n
       a(icol,l)=a(icol,l)*pivinv
      enddo
      do l=1,m
       b(icol,l)=b(icol,l)*pivinv
      enddo
      do ll=1,n
       if(ll.ne.icol)then
        dum=a(ll,icol)
        a(ll,icol)=0.
        do l=1,n
         a(ll,l)=a(ll,l)-a(icol,l)*dum
        enddo
        do l=1,m
         b(ll,l)=b(ll,l)-b(icol,l)*dum
        enddo
       endif
      enddo
      enddo
      do l=n,1,-1
       if(indxr(l).ne.indxc(l))then
        do k=1,n
         dum=a(k,indxr(l))
         a(k,indxr(l))=a(k,indxc(l))
         a(k,indxc(l))=dum
        enddo
       endif
      enddo
      return
      end
!===================================================================================================
!                                    end of subroutine gaussj                                      !
!===================================================================================================


!---------------------------------------------------------------------------------------------------
! note: there are two subroutines to calculate the solution of a linear system of kind a*x=b: solve
! and solvetwo.
! i am checking if they should be only one and, if 'yes', i will change their call inside the whole
! code. (ktw 120302).
! solve is when you have 3 calibration tanks available and solvetwo is when you have only two (nrc160707)
!===================================================================================================
!                                        subroutine solve                                          !
!===================================================================================================
!     given a linear system:  a*x = b, gets the solution: x = a^-1 b, and returns x.
!---------------------------------------------------------------------------------------------------
      subroutine solve(a, x, b) 
      real*8 a(3,3), x(3)
      real*8 acall(3,3), bcall(3)
      real*8 b(3)
      n  = 3
      np = 3
      m  = 1
      mp = 1
      do i = 1,3                                                 ! don't want inputs to be changed
      acall(1,i) = a(1,i) 
      acall(2,i) = a(2,i) 
      acall(3,i) = a(3,i) 
      bcall(i)   = b(i)
      enddo
!     get this from recipies
      call gaussj(acall,n,np,bcall,m,mp)
      do i=1,3
       x(i)=bcall(i)                                              ! returns x
      enddo
      return
      end
!===================================================================================================
!                                     end of subroutine solve                                      !
!===================================================================================================



!===================================================================================================
!                                       subroutine solvetwo                                        !
!===================================================================================================
!     given a linear system:  a*x = b, gets the solution: x = a^-1 b, and returns x.
!---------------------------------------------------------------------------------------------------
      subroutine solvetwo(a, x, b)
      real*8 a(2,2), x(2)
      real*8 acall(2,2), bcall(2)
      real*8 b(2)
      n  = 2
      np = 2
      m  = 1
      mp = 1
      do i = 1,2                               ! don't want inputs to be changed
       acall(1,i)= a(1,i)
       acall(2,i)= a(2,i)
       bcall(i)= b(i)
      enddo
!     get this from recipies
      call gaussj(acall,n,np,bcall,m,mp)
      do i=1,2
       x(i)=bcall(i)                            ! returns x
      enddo
      return
      end
!===================================================================================================
!                                   end of subroutine solvetwo                                     !
!===================================================================================================



!===================================================================================================
!                                      subroutine gethhmmss                                        !
!===================================================================================================
      subroutine gethhmmss(tdoy)
      real*8 tdoy  ! decimal day-of year
      real*8 gmt   ! fraction of this day
      include 'lparams.for'                      ! among other things, defines vars in common block
                                                 ! 'sysvars', including hh,mm,ss.
      gmt = mod(tdoy,1d0)*86400.d0               ! seconds of start-time
      hh = int(gmt/3600.d0)
      mm = mod(gmt/3600.d0, 1.d0)*60 
      ss = mod(gmt, 60.d0)                       !  gmt-hh*3600-mm*60
      if(ss.ge.59.99) then
      ss = ss-60
      mm = mm+1
      endif
      if(mm.ge.60) then
      mm=mm-60
      hh=hh+1
      endif
      return
      end
!===================================================================================================
!                                  end of subroutine gethhmmss                                     !
!===================================================================================================



!===================================================================================================
!                                   subroutine writetimestamp                                      !
!===================================================================================================
      subroutine writetimestamp(filepoint, msgstring, lenmsg, decjd)
      include 'lparams.for'
      character*40 msgstring
      integer filepoint
      real*8 decjd
      leapdays = int(((year-1)-2000)/4)+1             ! number of leapdays in preceding years (year
                                                      ! in common block).
      daysprev = (year-2000)*365 + leapdays           ! number of days in preceding years (back to
                                                      ! 1/1/2000)
      intdoy = int(decjd - daysprev)   
      call gethhmmss(decjd)
      write(filepoint,10) msgstring(1:lenmsg),decjd,year,intdoy,hh,mm,ss
10    format(a40,'julday ',f8.3,' (',f5.0,' doy ',i3,', ',i2.2,':',
     &       i2.2,':',f5.2,' gmt)')
!     e.g.:  julday 1706.344 (2004 doy 245, 21:00:00.62 gmt)
      end
!===================================================================================================
!                                end of subroutine writetimestamp                                  !
!===================================================================================================



!===================================================================================================
!                                         function tempf                                           !
!===================================================================================================
      function tempf(x)
      real*8 tempf,x,xx
      real*8 y,y2,y3,z
!     thermistor volt --> deg c
      if (x .le. 0.0 .or. x .ge. 2.5) then
       tempf=-999.99
      return
      end if
      xx=30100.*x/(2.5-x)
      y=log(xx)
      y2=y*y
      y3=y2*y
      z=9.3548e-4+2.2071e-4*y+6.8076e-8*y2+1.2407e-7*y3
      tempf=1./z-273.15
      return
      end
!===================================================================================================
!                                     end of function tempf                                        !
!===================================================================================================




!===================================================================================================
!                                          function xmed                                           !
!===================================================================================================
      function xmed(x,n)
      real*8 x(n), xmed 
      call sort(n,x)
      n2=n/2
      if(2*n2.eq.n)then
      xmed=0.5*(x(n2)+x(n2+1))
      else
      xmed=x(n2+1)
      endif
      return
      end
!===================================================================================================
!                                       end of function xmed                                       !
!===================================================================================================



!===================================================================================================
!                                         subroutine sort                                          !
!===================================================================================================
!     nrc 161001 to orgnize scrambeled eddy high frequency data collected 2015 to 2016
!     sorts array 'ra' into ascending numerical order using headsort algorithm.  
!     n  = the length of input array ra.
!     ra is replaced on output by its sorted rearrangement (see hpsort "numerical recipies", p. 329).
!---------------------------------------------------------------------------------------------------
      subroutine sort(n,ra)
      integer i,ir,j,l
      real*8 rra
      real*8 ra(n)
      l=n/2+1
      ir=n
      if(n.lt.2)return             ! returns single element unchanged (important!)
10    continue
      if(l.gt.1)then
      l=l-1
      rra=ra(l)
      else
      rra=ra(ir)
      ra(ir)=ra(1)
      ir=ir-1
      if(ir.eq.1)then
      ra(1)=rra
      return
      endif
      endif
      i=l
      j=l+l
20    dowhile(j.le.ir)              ! loop
      if(j.lt.ir)then
      if(ra(j).lt.ra(j+1))j=j+1
      endif
      if(rra.lt.ra(j))then
      ra(i)=ra(j)
      i=j
      j=j+j
      else
      j=ir+1
      endif
      enddo                         ! until (j.gt.ir)
      ra(i)=rra
      goto 10
      return
      end
!===================================================================================================
!                                     end of subroutine sort                                      !
!===================================================================================================



!===================================================================================================
!                                     subroutine getfilterset                                      !
!===================================================================================================
!     procedure: getfilterset   (moved from lfluxsub on 24-sept-01)
!     usage: call getfilterset(filtcode,filtset)
!     a procedure that has a dual purpose:
!     first, it converts the string code input into an integer array. 
!     second, it performs some error checking in order to convert any user misteps into a 'safe 
!     mode' of no-filtering. it also assures the supremacy of the pseudo-filter over any actual 
!     filtering.
!---------------------------------------------------------------------------------------------------
      subroutine getfilterset(filtcode, filtset)
      character filtcode(5)
      integer filtset(5)
      if (filtcode(5) .gt. '0' .and. filtcode(5).lt.'9') then
      print *,'  Pseudofilter enabled.'
      do k=1,4
      filtset(k) = 0
      filtcode(k)= '0'
      enddo
      if (filtcode(5) .eq. '1') then
      print *, ' Using opposite eddy no lock status',
     &         ' bit as pseudofilter.'
      filtset(5)=1
      filtcode(5)='1'
      endif    
      else
      do k=1,4
      if(filtcode(k) .eq. '1') then
      filtset(k) =1
      filtcode(k) ='1'
      else
      filtset(k) =0
      filtcode(k) ='0'
      endif
      enddo
      filtset(5)=0
      filtcode(5)='0'          
      endif
      end
!===================================================================================================
!                                  end of subroutine getfilterset                                  !
!===================================================================================================



!===================================================================================================
!                                      subroutine realchar                                         !
!===================================================================================================
!     converts character string 'ch' (of length lench) to real (lench is necessary because length of
!     'ch' can vary).
!     warning: it will return a value even if 'ch' is not representing number!
!---------------------------------------------------------------------------------------------------
      function realchar(ch, lench)
      character*15 ch
      real*4 realchar, n, pwr
      n = 0
      ndigits = lentrim0(ch, lench)  
      do i=1, ndigits
      pwr = 10**(ndigits-i)
      n = n + (ichar(ch(i:i))-48)*pwr
!     ichar returns ascii value, 48 is ascii value of '0' ( and 49 is ascii of '1', etc)
      enddo
      realchar = n
      end
!===================================================================================================
!                                    end of subroutine realchar                                    !
!===================================================================================================



!===================================================================================================
!                                      subroutine chareal                                          !
!===================================================================================================
!     converts real to character string
!---------------------------------------------------------------------------------------------------
      subroutine chareal(n,nd,ch)
      real n
      character*15 ch 
      ndigits = nd
      if(nd.eq.0) ndigits = int(log10(n)+1)     ! number of places we will need
      do i=ndigits,1,-1                         ! up to 15 places
      j = ndigits-i+1                           ! reverse order index to character string
                                                ! (most significant decimal goes to position one)
      idigit = int( (n - int(n/10**(i))*10**i)/10**(i-1) )
      ch(j:j) = char(idigit+48)
      enddo
      end 
!===================================================================================================
!                                    end of subroutine chareal                                     !
!===================================================================================================



!===================================================================================================
!                                      subroutine chareal8                                         !
!===================================================================================================
!     converts real to character string
!---------------------------------------------------------------------------------------------------
      subroutine chareal8(n,nd,ch)
      real*8 n
      character*15 ch 
      ndigits = nd
      if(nd.eq.0) ndigits = int(log10(n)+1)  ! number of places we will need
      do i=ndigits,1,-1                      ! up to 15 places
      j = ndigits-i+1                        ! reverse order index to character string
                                             ! (most significant decimal goes to position one)
      idigit = int( (n - int(n/10**(i))*10**i)/10**(i-1) )
      ch(j:j) = char(idigit+48)
      enddo
      end 
!===================================================================================================
!                                    end of subroutine chareal8                                    !
!===================================================================================================


!===================================================================================================
!                                      subroutine read20n                                          !
!===================================================================================================
      subroutine read20n(filepoint,line20n,dtype,decjd,v20n,nv1,nv2,
     &                             is201,ns1,ns2,done, yrdays, debug)
      include 'lparams.for'            ! among other things, defines vars in common block 'sysvars', 
                                       ! including hh,mm,ss.
!     arguments:
      logical done, debug
      integer filepoint, line20n, nv1, nv2
      real*8 v20n(10+22+10)
      real*8 decjd                     ! decimal julian day (days since 1-jan-00)
!     internal vars:
      integer doy,nerr,f20n
      real*8 gmtsecs    
      common /ewg/ gmtsecs,doy
!  
      if(dtype.eq.'PF') then
       if(line20n.eq.201) then
        assign 211 to f20n               ! format label for inputting 201 line
       elseif(line20n.eq.202) then
        assign 221 to f20n
       elseif(line20n.eq.203) then
        assign 231 to f20n
       else 
        goto 30
       endif
      elseif(dtype.eq.'GD') then
       if(line20n.eq.201) then
        assign 212 to f20n
       elseif(line20n.eq.202) then
        assign 222 to f20n
       elseif(line20n.eq.203) then 
        assign 232 to f20n
       else 
        goto 30
       endif
      elseif(dtype.eq.'EA' .or. dtype.eq.'EB') then
       if(line20n.eq.200) then
        assign 200 to f20n
       elseif(line20n.eq.201) then
        assign 201 to f20n
       else 
        goto 30
       endif
      else
       goto 30
      endif
      done = .false.
      nerr=0
10    continue
      if(line20n.eq.201) then
       read(filepoint,f20n,end=90,err=40) doy,gmtsecs,
     &    (v20n(i),i=nv1,nv2),(is201(i),i=ns1,ns2)
!     201 line always has status bits, others do not
      elseif (line20n .eq. 200) then 
       read(filepoint,f20n,end=90,err=40) doy,gmtsecs,
     &    (v20n(i),i=nv1,nv2),(is201(i),i=ns1,ns2)
!     200 line for eddy 1/2 has status bits
      else
       read(filepoint,f20n,end=90,err=40) doy,gmtsecs,(v20n(i),i=nv1,nv2)
      endif
!     year boundary crossing check 
      dtime= days00 + yrdays + doy+gmtsecs/dble(86400.) - decjd  
!     dtime is the difference between this time and last time;
!     days00 is number of days in preceding years, 
!     yrdays=0 or 365 if year-boundary crossed.
      if(dtime.lt.-300) then
      yrdays = int((-dtime)+0.1)                 ! number of days in this year depends 
                                                 ! on height of cliff we fell off of
                                                 ! (initialized to zero in linitial)
      write(99,15) line20n,yrdays 
15    format('  New year boundary detected in *.',i3,
     &       ' file (yrdays = ',f8.2,')' )  
      endif
      decjd = days00 + yrdays + doy+gmtsecs/dble(86400.)  ! decimal julian days since 1-jan-00
!     days00=days in preceding years set in linitial7.for,
!     yrdays= 0 unless we crossed a year boundary 
!     (in which case, it is = 365 or 366)
      if(debug) then
       write(99,100) dtype, line20n, doy, gmt    
       call writetimestamp(99, '        --> decimal ', 20, decjd)
       debug = .false.
      endif
100   format('  Subroutine read20n: on file ',a2,'*.',i3,
     &       ':  doy: ', i3, ', ', f9.3,' sec')
      return
!     errors:
30    continue
      write(99,103) dtype, line20n 
      write(*, 103) dtype, line20n
      call writetimestamp(99, '        --> called  ', 20, decjd)
      call writetimestamp(6 , '        --> called  ', 20, decjd)  
103   format('  Warning:Bad specification of dtype or line number in: ',
     &       a2,'*.',i3)
      return
      
40    continue
      nerr = nerr+1
      if (mod(nerr,10) .eq. 1) then 
       write(99,105) dtype, line20n, doy, gmt, nerr 
       write(*, 105) dtype, line20n, doy, gmt, nerr 
       call writetimestamp(99, '        --> decimal ', 20, decjd)
       call writetimestamp(6 , '        --> decimal ', 20, decjd)  
      endif   
!     unit 6 (same as *) = terminal
105   format('   Warning: Error in subroutine read20n: file ',a2,
     &       '*.',i3,':  doy: ', i3, ', ',f10.3,' sec','  nerr = ',i7)
      goto 10           ! continue until good value found 
!     end-of-file:
90    done = .true.
      return
!     formats for input-files (same as output formats in l1split):
200   format(i4,1x,f9.3,1x,2(f7.4,1x),4(f7.3,1x),f6.0,1x,i3,1x,i2)       ! eddy 200 
201   format(i4,1x,f9.3,/,8f8.3,/,2(f9.4,1x),1x,3i1,1x,2i1,1x,3i1,1x,i1) ! eddy  
211   format(i4,1x,f8.1,/,8f8.3,/,4(f9.4,1x),2(3i1,1x,2i1,1x),      
     &       5i1,1x,2(4i1,1x),i1)                                        ! profile 201
212   format(i4,1x,f8.1,/,7f8.3,1x,2(3i1,1x),3(4i1,1x))                  ! ground  201
221   format(i4,1x,f8.1,/,6(f8.3,1x)/,6(f8.3,1x),/,8(f8.3,1x) )          ! profile 202
222   format(i4,1x,f8.0,3(/,8f8.3))                                      ! (5-min avgs) ground  202 
231   format(i4,1x,f8.1,/,10f7.3)                                        ! profile 203
232   format(i4,1x,f8.0, 8f8.3)                                          ! (5-min avgs) ground  203 
      end
!===================================================================================================
!                                   end of subroutine read20n                                      !
!===================================================================================================



!===================================================================================================
!                                      subroutine readcal                                          !
!===================================================================================================
!     reads calibration tank parameters from text file lcaltank.inp, finds the correct caltank value
!     based on the time of the calibration sequence.
!---------------------------------------------------------------------------------------------------
      subroutine readcal(time1)
      include 'lparams.for'      ! among other things, defines variables in common-block 'sysvars', 
                                 ! cref and csurv are read here.  (as of november, 2004)  
                                 ! the remaining parameters wer read in readsysd:  
                                 ! common /sysvars/ cref,csurv,flsam,flsampf,flcal,flcalpf,
                                 !     x  zflsam,zflcal,zpcon,zpcell,gflsam,gflcal,gpcon,gpecell,
                                 !     x  gtecell, gppcell,gtpcell,clicor,hlicor,t0co2,toh2o,
                                 !     x  tlag,lagco2,lagh2o,tau, tauh, gu,gv,gw,gt,dirofu,               
                                 !     x  tdelayed,tavecal,tflstep, ndelayed,navecal,nflstep, 
                                 !     x  nlaglo,navecalo,nflsteplo,
                                 !     x  tavepf,tdelaypf,navepf,ndelaypf
                                 ! and also directory structure:
                                 ! common /dirnames/ rootinp, rootsplit, rootprocin, rootprocout, 
                                 !     x             dirsplit, dirprocin, dirprocout
      character dumline*80, pnds*1    
      real*8 sdate1, sdate2, jday1, jday2 
      lrt = lentrim0(rootinp, len(rootinp))
!      print *, '    File "'//rootinp(1:lrt)//slash//'lcaltank.inp":',
!     &         ' read cal tank concentrations ...'
      open(10, file=rootinp(1:lrt)//slash//'lcaltank.inp',status='old')
      read(10,*)
      read(10,*)
      read(10,*) sdate1, jday1, cref11, cref21, cref31, csurv1
!                               475     400     325     376   --> cal-gas ppm  
20    read(10,900,end=30) pnds,dumline
900   format(a1,a80)  
 
      if (pnds .eq. "#") then 
      goto 20  
      else 
      read(dumline,*)  sdate2, jday2, cref12, cref22, cref32, csurv2 
      endif   
      if (time1 .gt. jday1 .and. time1 .lt. jday2) then 
      goto 30  
      else 
      sdate1 = sdate2    
      jday1  = jday2    
      cref11 = cref12  
      cref21 = cref22   
      cref31 = cref32   
      csurv1 = csurv2 
      goto 20   
      endif  
30    cref(1) = cref11 
      cref(2) = cref21 
      cref(3) = cref31 
      csurv = csurv1 
40     continue
!       write(6, 41) time1, jday1, cref(1), cref(2), cref(3), csurv  
!41      format(' Time: ', f10.4,1X,' Calibration tank concentrations: ',f10.4,2X,4(f7.2,1X)) 
       close(10)
       return
       end
!===================================================================================================
!                                    end of subroutine readcal                                     !
!===================================================================================================   



!===================================================================================================
!                                      subroutine readco2cal                                       !
!===================================================================================================
!     reads calibration coefficients from text file lcalcoeff.inp, finds the correct calcoefficient
!     value based on the time of the calibration sequence. HS:coeff1, MS:coeff2, LS: coeff3
!---------------------------------------------------------------------------------------------------
      subroutine readco2cal(time1,types)
      include 'lparams.for'      ! fcoef are read here.  (as of june, 2016)  
                                 ! requires directory structure:
                                 ! common /dirnames/ rootinp, rootsplit, rootprocin, rootprocout, 
                                 !     x             dirsplit, dirprocin, dirprocout
      character*2 types
      character dumline*180, pnds*1    
      real*8 jday1, jday2 
      lrt = lentrim0(rootinp, len(rootinp))
!      print*, ' File:', rootinp(1:lrt)//slash//types//'lcalcoeff.inp'
!     print *, ' read cal values to be used if no CO2 tank concentrations are available...'
      open(10, file=rootinp(1:lrt)//slash//types//'lcalcoeff.inp',status='old')
      read(10,*)
      read(10,*)
!      read(10,*) jday1, fcoef31, fcoef21, fcoef11
      read(10,*) jday1, fcoef11, fcoef21, fcoef31
!                               coef1    coef2    coef3		>-ini cal factors  
20    read(10,900,end=30) pnds, dumline
900   format(a1,a100)
      if (pnds .eq. "#") then 
       goto 20  
      else
!       read(dumline,*)  jday2, fcoef32, fcoef22, fcoef12
       read(dumline,*)  jday2, fcoef12, fcoef22, fcoef32
!       print*,'jday2, fcoef12, fcoef22, fcoef32: ',jday2, ' ', fcoef12, ' ', fcoef22, ' ', fcoef32
      endif
      if ((time1.gt.jday1).and.(time1.lt.jday2)) then 
       goto 30  
      else 
       jday1 = jday2    
       fcoef11 = fcoef12  
       fcoef21 = fcoef22   
       fcoef31 = fcoef32   
       goto 20   
      endif  
30    fcoef(1) = fcoef11 
      fcoef(2) = fcoef21 
      fcoef(3) = fcoef31 
40    continue
!      write(6, 41) time1, jday1, fcoef(1), fcoef(2), fcoef(3)  
!41    format(' Time:',F10.4,1X,'Forced cal coefficients: ',F10.4,2X,3(F7.2,1X)) 
      close(10)
      return
      end
!===================================================================================================
!                                 end of subroutine readco2cal                                     !
!===================================================================================================   


!===================================================================================================
!                                      subroutine readh2ocal                                          !
!===================================================================================================
!     reads calibration tank parameters from text file lcaltank.inp, finds the correct caltank value

!     based on the time of the calibration sequence.
!---------------------------------------------------------------------------------------------------
      subroutine readh2ocal(time1,types)
      include 'lparams.for'      ! among other things records directory structure:
                                 ! common /dirnames/ rootinp, rootsplit, rootprocin, rootprocout, 
                                 !     x             dirsplit, dirprocin, dirprocout
      character*2 types
      character   dumline*180    
      real*8      jday1,jday2
      integer*4   j
      real*8      cvapor11,cvapor12

      lrt = lentrim0(rootinp, len(rootinp))
!     print *, '    File:', rootinp(1:lrt)//slash//'h20cal'//types//'.inp'
!     print *, ' read H2O scaling factor'
      open(10, file=rootinp(1:lrt)//slash//types//'h2ocal.inp',status='old')
      read(10,*)
      read(10,*) jday1, cvapor11
      j=1
      do while(j.le.int(time1*2))              ! loop
       read(10,900) dumline
900    format(a60)
       read(dumline,*) jday2, cvapor12
       if ((time1.ge.jday1) .and. (time1.le.jday2)) then 
        j      = j+int(time1*2)  
        cvapor = cvapor11
       else 
        jday1    = jday2   
        cvapor11 = cvapor12
       end if  
       j=j+1
      end do
        cvapor = cvapor11
      close(10)
      end
!===================================================================================================
!                                    end of subroutine readh2ocal                                     !
!===================================================================================================   

