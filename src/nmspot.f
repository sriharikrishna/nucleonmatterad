c ********************************
c nmspot
c single-particle potential driver
c ********************************
      program nmspot
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (nlog=0,nin=5,nout=6)
      real*8 kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,h2m,h2mcs,pi,s
      common /consts/ kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,
     &       h2m,h2mcs,pi,s
      real*8 aa(8),ab(8),ad(8,8),ae(6,2),af(8),ak(8,8,8),al(6,6,6),
     &       as(6),at(8,8),ax(6,6,6)
      common /amatrx/ aa,ab,ad,ae,af,ak,al,as,at,ax
      real*8 u,uf,up,tnia,tnic,tniu,tnix,cut,cut0,w3v0,w3v1,w3va,w3vc
      common /tbcnst/ u,uf,up,
     &       tnia,tnic,tniu,tnix,cut,cut0,w3v0,w3v1,w3va,w3vc
      real*8 eav,fsof,plm,qmin,qmax
      common /pionic/ eav,fsof,plm,qmin,qmax
      real*8 temp,mstar,chmpot,entrpy,ksav,kqav
      common /hotted/ temp,mstar,chmpot,entrpy,ksav,kqav
      real*8 xph,yph
      common /parhol/ xph,yph
      dimension e(5,20),gint(6),x(5),y(100)
      character*8 mname(2)
      character*16 pname(10),tname(0:5)
      character*20 timdat
      character*50 sysdat
      data mname/'Nuclear ','Neutron '/
      data pname/'Malfliet-Tjon V','Reid v8','Argonne v8','Argonne v12'
     &,'Argonne v14','Argonne v18','Urbana v14','Urbana v14 + TNI'
     &,'SSC(C)','Paris'/
      data tname/'                ',' + Urbana Vijk  '
     &          ,' + Tucson Vijk  ',' + Brazil Vijk  '
     &          ,' + DD TNR       ',' + DD TNR & TNA '/
c
      timeit=timer(0.)
      call header(sysdat,timdat)
      write(nlog,6) sysdat,timdat
      write(nout,6) sysdat,timdat
    6 format(/80('*')//1x,a50,8x,a20)
c
      pi=acos(-1.)
      do 1 i=1,8
      do 1 j=1,8
    1 ad(i,j)=-ad(i,j)/9.
      do 2 i=1,6
      do 2 j=1,2
    2 ae(i,j)=ae(i,j)/9.
c   ------------------
c   read in parameters
c   ------------------
      read(nin,1001) nm,np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
 1001 format(5x,i2)
      read(nin,1002) kf,rho,dor,acn,ast,atn,als,bst,btn,bls,cn,cne
 1002 format(5x,f10.4)
      read(nin,1002) temp,mstar,tnia,tnic,tniu,tnix,cut,cut0
      read(nin,1001) npi,npf
      read(nin,1002) eav,fsof,plm,qmin,qmax
      write(nlog,1010) mname(nm)
      write(nout,1010) mname(nm)
 1010 format(/4x,a8,'Matter')
      if (temp.gt.0.) then
        write(nlog,1011) temp
        write(nout,1011) temp
 1011   format(/4x,'at t =',f6.2,' mev')
      end if
      write(nlog,1012) pname(np),tname(nt)
      write(nout,1012) pname(np),tname(nt)
 1012 format(/4x,2a16)
      s=float(4/nm)
      if (kf.eq.0.) kf=(6.*pi**2*rho/s)**(1./3.)
      if (rho.eq.0.) rho=s*kf**3/6./pi**2
c  ------------------------
c  ground-state calculation
c  ------------------------
      nix=npi
      niy=npf
      ekf=eav
      xmin=fsof
      xmax=plm
      ymin=qmin
      ymax=qmax
      npi=0
      ns=0
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk,
     &            dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      e0=efree
      no=0
      ns=1
      write(nout,1111)
 1111 format('1')
      write(nlog,1101)
      write(nout,1101)
 1101 format(/4x,'compute single particle potenial function e(x,y)')
c  ----------------------------------------------------
c  cycle through yph=qmin,qmax,qinc=(qmax-qmin)/(npi-1)
c  ----------------------------------------------------
      if (nix.eq.1) then
        xinc=0
      else
        xinc=(xmax-xmin)/(nix-1)
      end if
      if (niy.eq.1) then
        yinc=0
      else
        yinc=(ymax-ymin)/(niy-1)
      end if
      do 20 ix=1,nix
        x(ix)=xmin+(ix-1)*xinc
        do 10 iy=1,niy
          y(iy)=ymin+(iy-1)*yinc
          if (y(iy).eq.0.) y(iy)=.01
          xph=x(ix)
          yph=y(iy)
          if (yph.lt.kf) xph=-x(ix)
          call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk,
     &                dor,bst,btn,bls,npi,npf, gint, endiff, efree)
          e(ix,iy)=efree
          eph=(e(ix,iy)-e0)/xph+ekf
          ue=eph-.5*h2m*yph**2
          write(nlog,1102) e(ix,iy),x(ix),y(iy),eph,ue
          write(nout,1102) e(ix,iy),x(ix),y(iy),eph,ue
 1102     format(/2x,'e(x,y) = ',f8.3,' x = ',f8.3,' y = ',f8.3
     &          ,' eph = ',f8.3,' ue = ',f8.3)
   10   continue
   20 continue
      timeit=timer(timeit)
      write(nlog,999) timeit
      write(nout,999) timeit
  999 format(/5x,'job time =',f8.3,' seconds')
      stop
      end
c *********************
c block data subprogram
c *********************
      block data
c
      real*8 aa(8),ab(8),ad(8,8),ae(6,2),af(8),ak(8,8,8),al(6,6,6),
     &       as(6),at(8,8),ax(6,6,6)
      common /amatrx/ aa,ab,ad,ae,af,ak,al,as,at,ax
c
      data aa/1.,3.,3.,9.,6.,18.,1.,3./
      data ab/1.,3.,3.,9.,1.,3.,1.,3./
      data ad/9*0.,12.,0.,12.,0.,12.,0.,12.,2*0.,4*12.,2*6.,0.,2*12.
     &,8.,12.,8.,6.,10.,2*0.,4*12.,2*6.,0.,2*12.,8.,12.,8.,6.,10.
     &,2*0.,4*6.,2*3.,0.,12.,6.,10.,6.,10.,3.,11./
      data ae/2*0.,4*6.,2*0.,6.,-2.,6.,-2./
      data af/4*1.,4*0./
      data ak/1.,8*0.,3.,8*0.,3.,8*0.,9.,8*0.,6.,8*0.,18.,8*0.,1.,8*0.
     &,3.,  0.,1.,6*0.,1.,-2.,9*0.,3.,6*0.,3.,-6.,9*0.,6.,6*0.,6.,-12.
     &,9*0.,1.,6*0.,1.,-2.,  2*0.,1.,8*0.,3.,4*0.,1.,0.,-2.,6*0.,3.
     &,0.,-6.,8*0.,2.,8*0.,6.,8*0.,.333,8*0.,1.
     &,3*0.,1.,6*0.,1.,-2.,5*0.,1.,0.,-2.,4*0.,1.,2*-2.,4.,9*0.,2.,6*0.
     &,2.,-4.,9*0.,.333,6*0.,.333,-.667,  4*0.,1.,8*0.,3.,6*0.,1.,8*0.
     &,3.,2*0.,1.,0.,1.,0.,-2.,0.,1.,2*0.,3.,0.,3.,0.,-6.,0.,3.,4*0.
     &,1.,0.,1.,6*0.,3.,0.,3.,  5*0.,1.,6*0.,1.,-2.,7*0.,1.,6*0.,1.,-2.
     &,3*0.,1.,0.,1.,0.,-2.,0.,2*1.,-2.,1.,2*-2.,4.,1.,-2.,5*0.,1.,0.,1.
     &,4*0.,1.,-2.,1.,-2.,  6*0.,1.,8*0.,3.,6*0.,1.,8*0.,3.,6*0.,1.,8*0.
     &,3.,1.,0.,1.,0.,1.,0.,1.,2*0.,3.,0.,3.,0.,3.,0.,3., 7*0.,1.,6*0.
     &,1.,-2.,7*0.,1.,6*0.,1.,-2.,7*0.,1.,6*0.,1.,-2.,0.,1.,0.,1.,0.,1.
     &,0.,2*1.,-2.,1.,-2.,1.,-2.,1.,-2./
      data al/1.,6*0.,3.,6*0.,3.,6*0.,9.,6*0.,6.,6*0.,18.
     &,0.,3.,4*0.,3.,6.,7*0.,9.,4*0.,9.,18.,7*0.,18.,4*0.,18.,36.
     &,2*0.,3.,6*0.,9.,2*0.,3.,0.,6.,4*0.,9.,0.,18.,6*0.,-6.,6*0.,-18.
     &,3*0.,9.,4*0.,9.,18.,3*0.,9.,0.,18.,2*0.,9.,2*18.,36.,7*0.,-18.
     &,4*0.,-18.,-36.,4*0.,6.,6*0.,18.,4*0.,-6.,6*0.,-18.,6.,0.,-6.,0.
     &,12.,2*0.,18.,0.,-18.,0.,36.,5*0.,18.,4*0.,18.,36.,5*0.,-18.,4*0.
     &,-18.,-36.,0.,18.,0.,-18.,0.,36.,18.,36.,-18.,-36.,36.,72./
      data as/2*2.25,2*1.25,2*1./
      data at/1.,8*0.,1.,8*0.,1.,8*0.,1.,8*0.,1.,8*0.,1.,8*0.,1.,8*0.,1.
     &/
      data ax/1.,42*0.,1.,42*0.,1.,13*0.,1.,28*0.,1.,13*0.,1.
     &,16*0.,1.,9*0.,1.,0.,1.,30*0.,1.,9*0.,1.,0.,1./
      end
