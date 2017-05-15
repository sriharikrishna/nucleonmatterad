c *id* header **********************************************************
c subroutine for time, date, and system information
c calf denotes alpha-specific functions
c cray denotes cray-specific functions
c cibm denotes ibm-specific functions
c crs6 denotes rs/6000-specific functions
c csun denotes sun-specific functions
c cvax denotes vax-specific functions
c ----------------------------------------------------------------------
      subroutine header(sysdat,timdat)
      implicit real*8 (a-h,o-z)
      character*20 timdat
      character*50 macdat,sysdat
      character*9 thedat,thetim
cibm  integer idtx(14)
      character*26 rsdate
csun  integer idtx(3)
      data thedat/'today'/,thetim/'now'/
cray  data macdat/'NCSA  Cray Y-MP/464  UNICOS 6.0  CFT77 5.0'/
cray  data macdat/'NERSC  Cray-2S/8-128  UNICOS 6.0  CFT77 5.0'/
cray  data macdat/'NERSC  Cray Y-MP/C916  UNICOS 7.C  CFT77 5.0'/
crs6  data macdat/'Theory  IBM RS/6000-397  AIX 4.3.2  XL Fortran 5.1'/
crs6  data macdat/'Theory  IBM RS/6000-370  AIX 3.2.5  XL Fortran 3.2'/
crs6  data macdat/'Theory2 IBM RS/6000-320  AIX 3.2.5  XL Fortran 3.2'/
csp1  data macdat/'ANL-MCS  IBM SP-1  AIX 3.2.5  XL Fortran 3.2'/
cqsp  data macdat/'ANL-MCS  Quad SP   AIX 4.2    XL Fortran 4.1'/
chib  data macdat/'ANL-MCS Chiba City  LINUX  Portland Group F90'/
calf  data macdat/'CFNUL DEC AXP  OSF/1 V1.3  DEC Fortran X3.2-379'/
ctl8  data macdat/'Theoryl8 Pentium-4  Linux 7.1  IFC 5.0'/
cjaz  data macdat/'JAZZ Pentium-4  Redhat Enterprise Linux  Intel-8.1'/
cmbp  data macdat/'MacBook Pro Core 2 Duo  MacOSX 10.4  g95 0.9'/
cfus  data macdat/'FUSION Nehalem qc  Redhat Linux 5.4  Intel-11.1'/
cfd   data macdat/'Full Disclosure  SiCortex  pathf95'/
      data macdat/'Theoryl18 i7-920  SciLinux 6.1  IFC 13.0'/
      sysdat=macdat
calf  call time(thetim)
calf  call date(thedat)
cray  write(thetim,1) clock()
cray  write(thedat,1) date()
cray1 format(a8)
cibm  call datimx(idtx)
cibm  write(thetim,1) idtx(5),idtx(4),idtx(3)
cibm1 format(i2,':',i2,':',i2)
cibm  write(thedat,2) idtx(7),idtx(6),idtx(14)
cibm2 format(i2,'/',i2,'/',i2)
crs6  call fdate_(rsdate)
crs6  write(thetim,1) rsdate(12:19)
crs61 format(a8)
crs6  write(thedat,2) rsdate(9:10),rsdate(5:7),rsdate(23:24)
crs62 format(a2,'-',a3,'-',a2)
csun  call itime(idtx)
csun  write(thetim,1) idtx
csun1 format(i2,':',i2,':',i2)
csun  call idate(idtx)
csun  idtx(3)=idtx(3)-1900
csun  write(thedat,2) idtx
csun2 format(i2,'/',i2,'/',i2)
cvax  call time(thetim)
cvax  call date(thedat)
      call fdate_(rsdate)
      write(thetim,1) rsdate(12:19)
    1 format(a8)
      write(thedat,2) rsdate(9:10),rsdate(5:7),rsdate(23:24)
    2 format(a2,'-',a3,'-',a2)
      timdat=thetim//'  '//thedat
      return
      end
c *id* timer ***********************************************************
c function for elapsed cpu time
c calf denotes alpha-specific functions
c cray denotes cray-specific functions
c cibm denotes ibm 370-specific functions
c crs6 denotes ibm rs/6000-specific functions
c csun denotes sun-specific functions
c cvax denotes vax-specific functions
c ----------------------------------------------------------------------
      function timer(arg)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*4 sec(2),x(2)
      real*8 dtime
      if (arg.eq.0.) then
        timer=1.e-6
calf    call etime(sec)
cray    no call
cibm    call cputime(time,ircode)
crs6    no call
csun    call etime(sec)
cvax    call lib$init_timer
      else
        timer=dtime(x)-arg
calf    call etime(sec)
calf    timer=sec(1)-arg
cray    timer=second(0)-arg
cibm    call cputime(time,ircode)
cibm    timer=1.e-6*time-arg
crs6    isec=mclock()
crs6    timer=1.e-2*isec-arg
csun    call etime(sec)
csun    timer=sec(1)-arg
cvax    call lib$stat_timer(2,itime)
cvax    timer=.01*float(itime)
      end if
      return
      end
c *id* fdate_ **********************************************************
      subroutine fdate_ ( rsdate )
      character*26 rsdate
      character*8 date
      character*10 time
      character*5 zone
      integer values(8)
      character*3 months(12)
      data months / 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
     &   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /
!
      call date_and_time ( date, time, zone, values )
!
      rsdate = "1234MMM8dd1hh:mm:ss012yy56"
      rsdate(5:7) = months(values(2))
      rsdate(9:10) = date(7:8)
      rsdate(12:19) = time(1:2) // ":" // time(3:4) // ":"
     &   // time(5:6)
      rsdate(23:24) = date(3:4)
      return
      end
