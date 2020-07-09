c *id* nucmat **********************************************************
c subroutine for driving nuclear/neutron matter code
c ----------------------------------------------------------------------
#if defined (BFGS) && defined (ALLOW_TAPENADE)
      subroutine nucmat(x,n,f,flocald,lprt)
#else
      subroutine nucmat(x,n,f,lprt)
#endif
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      INCLUDE 'DIFFSIZES.inc'
      parameter (nlog=0,nin=5,nout=6)
      logical lprt
      real*8 x(n)
      common /minim/ econ,ncon,ntype
c
      real*8 kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,h2m,h2mcs,pi,s
      common /consts/ kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,
     &       h2m,h2mcs,pi,s
#ifdef ALLOW_TAPENADE
#if defined (DO_ALL)
      COMMON /consts_dv/ astd, atnd, alsd, dtd, drd, evxd
      real*8 xperturb(30,7)
#else
      common /consts_d/ kfd, rhod, acnd, astd, atnd, alsd, cnd, cned,
     & dtd, drd, evxd, h2md, h2mcsd, pid, sd
#endif
#endif
#if defined (DO_ALL) && defined (ALLOW_TAPENADE)
      real*8 kfd(nbdirsmax), rhod(nbdirsmax), acnd(nbdirsmax),
     & astd(nbdirsmax), atnd(nbdirsmax), alsd(nbdirsmax),
     & cnd(nbdirsmax), cned(nbdirsmax),
     & dtd(nbdirsmax), drd(nbdirsmax), evxd(nbdirsmax),
     & h2md(nbdirsmax), h2mcsd(nbdirsmax), pid(nbdirsmax),
     & sd(nbdirsmax)
      real*8 bstd(nbdirsmax), btnd(nbdirsmax),
     & blsd(nbdirsmax), dord(nbdirsmax)
#endif

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
      save nmlocal,np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk,
     &npi,npf,bst,btn,bls,nosave,npisav
c
      real*8 gint(6)
#if defined (DO_ALL) && defined (ALLOW_TAPENADE)
      real*8 gintd(nbdirsmax,6),endiffd(nbdirsmax),
     &efreed(nbdirsmax), flocald(nbdirsmax)
#else
      real*8 gintd
#endif
      character*20 pname(30),tname(0:5),ptnnam
      character*20 timdat
      character*32 mname(4)
      character*50 sysdat
      data mname/'Nuclear matter','Neutron matter'
     &          ,'Not implemented at this time'
     &          ,'Spin-polarized neutron matter'/
      data pname/'Malfliet-Tjon V','Reid v8','Urbana v14'
     &,'Argonne v8','Argonne v9','Argonne v14'
     &,'Argonne v18-csbl','Argonne v18-csbs','Argonne v18'
     &,'Argonne v8p','Argonne v6p','Argonne v4p','Argonne v2p'
     &,'Argonne v1p','Argonne vxp','Argonne v9p (1D)','Argonne v9p (Da)'
     &,'undefined','undefined','undefined'
     &,'Av18 v1.9','Av8p v1.9','Av18 v1.7','Av8p v1.7'
     &,'undefined','Argonne v18p','Argonne v18pq'
     &,'SSC(C) v14','SSC(C) v8p mod','Paris'/
      data tname/'                ',' + Urbana Vijk  '
     &          ,' + Tucson Vijk  ',' + Brazil Vijk  '
     &          ,' + DD TNR       ',' + DD TNR & TNA '/
c
      dor=x(1)
      if (n.ge.2) then
        ast=x(2)
        atn=ast
        als=ast
      end if
      if (n.eq.4) then
        bst=x(3)
        btn=x(4)
#ifdef CASE_SNM
        bls=bst
#else
        bls=0.
#endif
      end if
#ifdef ALLOW_TAPENADE
      astd=0.0
      atnd=0.0
      alsd=0.0
      bstd=0.0
      btnd=0.0
      blsd=0.0
      dord=0.0
      gintd =0.0
      flocald = 0.0
      gint = 0.0
#ifdef DO_ALL
      ndirs=1
      dord(ndirs)=1.0
      if (n.ge.2) then
        ndirs=ndirs+1
        astd(ndirs)=1.0
        atnd=astd
        alsd=astd
      end if
      if (n.eq.4) then
        ndirs=ndirs+1
        bstd(ndirs)=1.0
        ndirs=ndirs+1
        btnd(ndirs)=1.0
c        bls=bst
        blsd=bstd
      end if
#endif
#endif
#ifndef ALLOW_TAPENADE
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
        g2=g2+(gint(l)+1.)**2
      end do
      f=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
#else
#ifndef DO_ALL
      call NMMAINAD_D(np, nv, nt, ni, nie, no, ns, lf, lc, ls, lt
     &               , ll, lg, le, l3, lk, dor, dord, bst, bstd,
     &               btn, btnd, bls, blsd, npi, npf, gint, gintd
     &               , endiff, endiffd, efree, efreed, flocal,
     &               flocald, nmlocal)
#else
      call NMMAINAD_DV(np, nv, nt, ni, nie, no, ns, lf, lc, ls, lt
     &                  , ll, lg, le, l3, lk, dor, dord, bst, bstd
     &                  , btn, btnd, bls, blsd, npi, npf, gint,
     &                   endiff, efree, flocal, flocald, nmlocal,
     &                     nbdirsmax)
#endif
      f=flocal
#endif
       no=0
#if defined (ALLOW_TAPENADE)
#ifdef DO_ALL
      do i=1,nbdirsmax
        write(nlog,*) "flocald%d", flocald(i)
        write(nout,*) "flocald%d", flocald(i)
      end do
#else
      write(nlog,*) "flocald%d", flocald
      write(nout,*) "flocald%d", flocald
#endif
#endif
      return
c ************************
c entry for initialization
c ************************
      entry nminit(x,n)
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
      read(nin,1001) nmlocal,np,nv,nt,ni,nie,no,ns,lf,
     & lc,ls,lt,ll,lg,le,l3,lk
 1001 format(5x,i3)
      read(nin,1002) kf,rho,dor,acn,ast,atn,als,bst,btn,bls,cn,cne
 1002 format(5x,f10.4)
      read(nin,1002) temp,mstar,tnia,tnic,tniu,tnix,cut,cut0
      read(nin,1001) npi,npf
      read(nin,1002) eav,fsof,plm,qmin,qmax
      write(nlog,1010) mname(nmlocal)
      write(nout,1010) mname(nmlocal)
 1010 format(/4x,a32)
      if (temp.gt.0.) then
        write(nlog,1011) temp
        write(nout,1011) temp
 1011   format(/4x,'at T =',f6.2,' MeV')
      end if
      ptnnam=' '
      if (np.le.100) then
        ptnnam(1:20)=pname(np)
      else if (np.gt.100 .and. np.le.200) then
        ptnnam(1:13)='Norfolk v19 #'
        write(ptnnam(14:16),'(i3)') np
      end if
      write(nlog,1012) ptnnam,tname(nt)
      write(nout,1012) ptnnam,tname(nt)
 1012 format(/4x,2a20)
      s=float(4/nmlocal)
      do i=1, 30
        read(nin,*) (xperturb(i,j),j=1,7)
      end do
      read(nin,*) nperturb
      read(nin,*) delta
      write(nlog,*) "delta", delta
      dor=dor*(1+delta*xperturb(nperturb,1))
      ast=ast*(1+delta*xperturb(nperturb,2))
      als=als*(1+delta*xperturb(nperturb,3))
      atn=atn*(1+delta*xperturb(nperturb,4))
      bst=bst*(1+delta*xperturb(nperturb,5))
      bls=bls*(1+delta*xperturb(nperturb,6))
      btn=btn*(1+delta*xperturb(nperturb,7))
      x(1)=dor
      if (n.ge.2) x(2)=ast
      if (n.eq.4) then
        x(3)=bst
        x(4)=btn
      end if
      nosave=no
      npisav=npi
      no=1
      npi=0
      return
c *******************
c entry for final run
c *******************
      entry nmfin(x,n)
      ni=ni+5
      if (nie.gt.0) nie=nie+1
      no=nosave
      npi=npisav
c     lc=2*lc
c     ls=2*ls
c     lt=2*lt
c     ll=2*ll
c     lg=2*lg
c     le=2*le
c     l3=2*l3
      dor=x(1)
      if (n.ge.2) then
        ast=x(2)
        atn=ast
        als=ast
      end if
      if (n.eq.4) then
        bst=x(3)
        btn=x(4)
#ifdef CASE_SNM
        bls=bst
#else
        bls=0.
#endif
      end if
#ifdef ALLOW_TAPENADE
      astd=0.0
      atnd=0.0
      alsd=0.0
      bstd=0.0
      btnd=0.0
      blsd=0.0
      dord=0.0
      gintd =0.0
      flocald = 0.0
      gint = 0.0
#ifdef DO_ALL
      ndirs=1
      dord(ndirs)=1.0
      if (n.ge.2) then
        ndirs=ndirs+1
        astd(ndirs)=1.0
        atnd=astd
        alsd=astd
      end if
      if (n.eq.4) then
        ndirs=ndirs+1
        bstd(ndirs)=1.0
        ndirs=ndirs+1
        btnd(ndirs)=1.0
c        bls=bst
        blsd=bstd
      end if
#endif
#endif
#ifndef ALLOW_TAPENADE
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do 995 l=1,2,nmlocal
  995 g2=g2+(gint(l)+1.)**2
      final=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      fplus=final+abs(endiff)
#else
#ifndef DO_ALL
      call NMMAINAD_D(np, nv, nt, ni, nie, no, ns, lf, lc, ls, lt
     &               , ll, lg, le, l3, lk, dor, dord, bst, bstd,
     &               btn, btnd, bls, blsd, npi, npf, gint, gintd
     &               , endiff, endiffd, efree, efreed, flocal,
     &               flocald, nmlocal)
#else
      call NMMAINAD_DV(np, nv, nt, ni, nie, no, ns, lf, lc, ls, lt
     &                  , ll, lg, le, l3, lk, dor, dord, bst, bstd
     &                  , btn, btnd, bls, blsd, npi, npf, gint,
     &                   endiff, efree, flocal, flocald, nmlocal,
     &                     nbdirsmax)
#endif
      final = flocal
      fplus=flocal+abs(endiff)
#endif

      write(nlog,1095) final,fplus
      write(nout,1095) final,fplus
 1095 format(/2x,'final   fplus',/2f8.3)
      if (no.le.1) go to 999
      call nmout(le,lg,lt,l3,nie,no,nt,nv)
  999 return
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