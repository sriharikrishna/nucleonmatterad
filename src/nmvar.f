      module nmvar
      use params
      !implicit real*8 (a-h,o-z)
      !implicit integer*4 (i-n)
      
      integer*4, parameter:: nu=4/nm
      integer*4,parameter :: nlog=0
      integer*4,parameter :: nin=5
      integer*4,parameter :: nout=6

      !common minim 
      real*8     :: econ
      integer*4  :: ncon,ntype
      
      !common consts
      real*8     :: rho,acn,ast,atn,als,cn,cne
      real*8     :: dt,dr,evx,h2m,h2mcs,pi,s
      real*8  :: kf

      !common rslate 
      real*8, dimension(lgrid) :: r,ri,rs,sl,sls,slp,slps
      real*8, dimension(lgrid) :: sldp,sltp,rllp,rlssx,rsdsl

      !common correl 
      real*8     :: f(lgrid,8),fp(lgrid,8),fds(lgrid,8),v(lgrid,14)

c      !common amatrx 
c      !real*8     :: aa,ab,ad,ae,af,ak,al,as,at,ax
c      !real*8      :: aa(8),ab(8),ad(8,8),ae(6,2),
c     !&      af(8),ak(8,8,8),al(6,6,6),
c     !&      as(6),at(8,8),ax(6,6,6)

      real*8, dimension(8) :: aa = (/1.,3.,3.,9.,6.,18.,1.,3./)
      real*8, dimension(8) :: ab = (/1.,3.,3.,9.,1.,3.,1.,3./)
      real*8, dimension(8,8) :: ad = reshape([0.,0.,0.,0.,0.,0.,0.,0.,0.,
     &,12.,0.,12.,0.,12.,0.,12.,0.,0.,12.,12.,12.,12.,6.,6.,0.,12.,12.
     &,8.,12.,8.,6.,10.,0.,0.,12.,12.,12.,12.,6.,6.,0.,12.,12.,
     &8.,12.,8.,6.,10.
     &,0.,0.,6.,6.,6.,6.,3.,3.,0.,12.,6.,
     &10.,6.,10.,3.,11.],[8,8])
      real*8, dimension(6,2) :: ae = reshape([0.,0.
     &,6.,6.,6.,6.,0.,0.,6.,-2.,6.,-2.],[6,2])
      real*8, dimension(8) :: af = (/1.,1.,1.,1.,0.,0.,0.,0./)
      real*8, dimension(8,8,8) :: ak = reshape([1.,0.,0.,0.,0.,0.,0.,0.,
     & 0.,3.,0.,0.,0.,0.,0.,0.,
     & 0.,0.,3.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,9.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,6.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,18.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,1.,0.,
     & 0.,0.,0.,0.,0.,0.,0.,3.,
     & 0.,1.,0.,0.,0.,0.,0.,0.,
     & 1.,-2.,0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,3.,0.,0.,0.,0.,
     & 0.,0.,3.,-6.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,6.,0.,0.,
     & 0.,0.,0.,0.,6.,-12.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,0.,1.,
     & 0.,0.,0.,0.,0.,0.,1.,-2.,  
     & 0.,0.,1.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,3.,0.,0.,0.,0.,
     & 1.,0.,-2.,0.,0.,0.,0.,0.,
     & 0.,3.,0.,-6.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,2.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,6.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,.333,0.,
     & 0.,0.,0.,0.,0.,0.,0.,1.,
     & 0.,0.,0.,1.,0.,0.,0.,0.,
     & 0.,0.,1.,-2.,0.,0.,0.,0.,
     & 0.,1.,0.,-2.,0.,0.,0.,0.,
     & 1.,-2.,-2.,4.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,2.,0.,0.,
     & 0.,0.,0.,0., 2.,-4.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,0.,.333,
     & 0.,0.,0.,0.,0.,0.,.333,-.667,
     & 0.,0.,0.,0.,1.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,3.,0.,0.,
     & 0.,0.,0.,0.,1.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,3.,0.,0.,
     & 1.,0.,1.,0.,-2.,0.,1.,0.,
     & 0.,3.,0.,3.,0.,-6.,0.,3.,
     & 0.,0.,0.,0.,1.,0.,1.,0.,
     & 0.,0.,0.,0.,0.,3.,0.,3., 
     & 0.,0.,0.,0.,0.,1.,0.,0.,
     & 0.,0.,0.,0.,1.,-2.,0.,0.,
     & 0.,0.,0.,0.,0.,1.,0.,0.,
     & 0.,0.,0.,0.,1.,-2.,0.,0.,
     & 0.,1.,0.,1.,0.,-2.,0.,1.,
     & 1.,-2.,1.,-2.,-2.,4.,1.,-2.,
     & 0.,0.,0.,0.,0.,1.,0.,1.,
     & 0.,0.,0.,0.,1.,-2.,1.,-2.,
     & 0.,0.,0.,0.,0.,0.,1.,0.,
     & 0.,0.,0.,0.,0.,0.,0.,3.,
     & 0.,0.,0.,0.,0.,0.,1.,0.,
     & 0.,0.,0.,0.,0.,0.,0.,3.,
     & 0.,0.,0.,0.,0.,0.,1.,0.,
     & 0.,0.,0.,0.,0.,0.,0.,3.,
     & 1.,0.,1.,0.,1.,0.,1.,0.,
     & 0.,3.,0.,3.,0.,3.,0.,3.,
     & 0.,0.,0.,0.,0.,0.,0.,1.,
     & 0.,0.,0.,0.,0.,0.,1.,-2.,
     & 0.,0.,0.,0.,0.,0.,0.,1.,
     & 0.,0.,0.,0.,0.,0.,1.,-2.,
     & 0.,0.,0.,0.,0.,0.,0.,1.,
     & 0.,0.,0.,0.,0.,0.,1.,-2.,
     & 0.,1.,0.,1.,0.,1.,0.,1.,
     & 1.,-2.,1.,-2.,1.,-2.,1.,-2.],[8,8,8])
      real*8, dimension(6,6,6) :: al = reshape([1.,0.,0.,0.,0.,0.,
     & 0.,3.,0.,0.,0.,0.,
     & 0.,0.,3.,0.,0.,0.,
     & 0.,0.,0.,9.,0.,0.,
     & 0.,0.,0.,0.,6.,0.,
     & 0.,0.,0.,0.,0.,18.,
     & 0.,3.,0.,0.,0.,0.,
     & 3.,6.,0.,0.,0.,0.,
     & 0.,0.,0.,9.,0.,0.,
     & 0.,0.,0.,0.,9.,18.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,18.,0.,0.,0.,0.,
     & 0.,0.,18.,36.,0.,0.,
     & 3.,0.,0.,0.,0.,0.,
     & 0.,9.,0.,0.,3.,0.,
     & 6.,0.,0.,0.,0.,9.,
     & 0.,18.,0.,0.,0.,0.,
     & 0.,0.,-6.,0.,0.,0.,
     & 0.,0.,0.,-18.,0.,0.,
     & 0.,9.,0.,0.,0.,0.,
     & 9.,18.,0.,0.,0.,9.,
     & 0.,18.,0.,0.,9.,18.,
     & 18.,36.,0.,0.,0.,0.,
     & 0.,0.,0.,-18.,0.,0.,
     & 0.,0.,-18.,-36.,0.,0.,
     & 0.,0.,6.,0.,0.,0.,
     & 0.,0.,0.,18.,0.,0.,
     & 0.,0.,-6.,0.,0.,0.,
     & 0.,0.,0.,-18.,6.,0.,
     & -6.,0.,12.,0.,0.,18.,
     & 0.,-18.,0.,36.,0.,0.,
     & 0.,0.,0.,18.,0.,0.,
     & 0.,0.,18.,36.,0.,0.,
     & 0.,0.,0.,-18.,0.,0.,
     & 0.,0.,0.,0.,-18.,-36.,
     & 0.,18.,0.,-18.,0.,36.,
     & 18.,36.,-18.,-36.,36.,72.],[6,6,6])
      real*8, dimension(6) :: as = (/2.25,2.25,1.25,1.25,1.,1./)
      real*8, dimension(8,8) :: at = reshape([1.,0.,0.,0.,0.,0.,0.,0.,
     & 0.,1.,0.,0.,0.,0.,0.,0.,
     & 0.,0.,1.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,1.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,1.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,1.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,1.,0.,
     & 0.,0.,0.,0.,0.,0.,0.,1.],[8,8])
      real*8, dimension(6,6,6) :: ax = reshape([
     & 1.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,1.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,1.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,1.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,1.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,1.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,1.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,1.,0.,1.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,1.,
     & 0.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,1.,0.,1.],[6,6,6])
      
      !common gchain 
      real*8, dimension(lgrid,6) :: gca,gcb,gdd,gde,gee
      real*8, dimension(lgrid)   :: gl,gx,gy,gz
      real*8, dimension(lgrid,14):: gnn

      !common mocfun 
      real*8, dimension(lgrid,6) :: gfdd,gfde,gfed,gfcc,ghdd
      real*8, dimension(lgrid,6) :: ghde,ghed,ghcc
      real*8, dimension(lgrid,6) :: grdc,grdd,grde,gred,gree
      real*8, dimension(lgrid,6) :: grfc,grfd,grfe,grmd,grme

      !common tbpots 
      real*8     :: v3cc(lgrid,6,2),v3dd(lgrid,6,2),v3de(lgrid,6,2)
      real*8     :: v3ee(lgrid,6,2)


      !common sorfun 
      real*8     :: bj(8,6),bk(4,3),bq(6,2),vc(6,3,3)
      real*8     :: bcc(lgrid,3),bde(lgrid,3)

      !common tbcnst 
      real*8     :: u,uf,up,tnia,tnic,tniu,tnix
      real*8     :: cut,cut0,w3v0,w3v1,w3va,w3vc

      !common hotted 
      real*8     :: temp,chmpot,entrpy
      integer*4 :: mstar,ksav,kqav

      !common eblock 
      real*8     :: ev6,evb,evq,ek6,ekb,ef6,ej6,ejb,ep6
      real*8, dimension(14,10) :: wx
      real*8, dimension(8,10) ::wcx,wcdx,wcmx,wcrx
      real*8, dimension(6,4,2) ::w3x

      !common hotfun 
      integer*4, parameter :: ngrid=(20*lgrid+1)
      real*8, dimension(ngrid) :: rx, slx,slpx,sldpx,sltpx
      
      !common echain
      integer, parameter :: legrid=lgrid*(lgrid**2+1)/2 
      real*8, dimension(lgrid,6) :: eca,ecb,edd,ede,eee
      real*8, dimension(legrid) :: sccd,scce,sddd,sdde,sdee,seee

      !common tbfunc 
      real*8, dimension(lgrid) :: tpi,ypi,tpi2
      real*8, dimension(lgrid) :: xt0,xt1,xt2,xt3

      !common pionic 
      real*8     :: eav,fsof,plm,qmin,qmax
        
      !common parhol 
      real*8     :: xph,yph

      !common angle 
      real*8, dimension(legrid) :: xtheta,ytheta,ztheta,stheta
      integer*4  :: index(lgrid,lgrid,lgrid)

      end module nmvar
