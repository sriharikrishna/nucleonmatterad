c **********************************************************************
c nm
c program for nuclear/neutron matter
c **********************************************************************
      program nmprog
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
#ifdef DO_FULLX
      include "DIFFSIZES.inc"
#endif
      include "nclude/params.f"
      external nucmat
      parameter (nlog=0,nin=5,nout=6)
      parameter (maxdim=15)
      character*50 sysdat
      character*20 timdat
      logical lprt
      real*8 x(maxdim),scale(maxdim)
      common /minim/ econ,ncon,ntype
      character*4 etype(3)
      data etype/' jf ',' av ',' pb '/
      timeit=timer(0.)
      call header(sysdat,timdat)
      write(nlog,6) sysdat,timdat
      write(nout,6) sysdat,timdat
    6 format(/80('*')//1x,a50,8x,a20)
      read(nin,*) maxcl1,tol1,maxcl2,tol2,alpha,beta,gamma,ndim
      read(nin,*) (scale(i),i=1,ndim)
#if !defined(ONLY_NUCMAT) && !defined(BFGS) && defined(DO_FULLX)
      do i=ndim+1, nbdirsmax
        scale(i) = scale(2)
      end do
#endif
      read(nin,*) econ,ncon,ntype
      call nminit(x,ndim)
      if (maxcl1.ge.1) then
        write(nlog,7) etype(ntype+2),ndim
        write(nout,7) etype(ntype+2),ndim
    7   format(/4x,a4,'energy minimization in',i2,' dimensions')
        if (econ.ne.0.) then
          write(nlog,8) econ,ncon
          write(nout,8) econ,ncon
    8     format(5x,'with constraint',f6.0,
     &              '*sqrt(sum((1+gint(k))**2))**',i1)
        end if
        write(nlog,9) tol1,alpha,beta,gamma,tol2,(scale(i),i=1,ndim)
        write(nout,9) tol1,alpha,beta,gamma,tol2,(scale(i),i=1,ndim)
    9   format(5x,'simplex tolerance is',f6.2,' with coefficients of',
     &        /5x,'reflection ',f6.2,' contraction ',f6.2,
     &          ' expansion ',f6.2,
     &        /5x,'quadratic tolerance is',f6.2,
     &        /5x,'scale factors are',9f6.2)
#ifndef ONLY_NUCMAT
#ifndef BFGS
#ifndef DO_FULLX
        call minimi(maxcl1,tol1,maxcl2,tol2,alpha,beta,gamma,scale,x,
     &              fbest,nucmat,ndim)
#else
        call minimi(maxcl1,tol1,maxcl2,tol2,alpha,beta,gamma,scale,x,
     &              fbest,nucmat,nbdirsmax)
#endif
#else
        mbfgs=6
#ifndef DO_FULLX
        call sdrive(ndim,mbfgs,x(1:ndim),nucmat)
#else
        call sdrive(nbdirsmax,ndim,mbfgs,x(1:nbdirsmax),nucmat)
#endif
#endif
#else
#ifdef ALLOW_TAPENADE
#ifndef DO_FULLX
      call nucmat(x(1:ndim),ndim)
#else
      call nucmat(x(1:nbdirsmax),nbdirsmax)
#endif
#else
      f=0.0
#ifndef DO_FULLX
      call nucmat(x(1:ndim),ndim,f,.TRUE.)
#else
      call nucmat(x(1:nbdirsmax),nbdirsmax,0.0,.TRUE.)
#endif
#endif
#endif
      end if
!#if !defined (ALLOW_TAPENADE) && !defined (ONLY_NUCMAT)
      call nmfin(x,ndim)
!#endif
      timeit=timer(timeit)
      write(nlog,999) timeit
      write(nout,999) timeit
  999 format(/5x,'job time =',f8.3,' seconds')
      stop
      end
#ifndef ALLOW_TAPENADE
c *id* nucmat **********************************************************
c subroutine for driving nuclear/neutron matter code
c ----------------------------------------------------------------------
      subroutine nucmat(x,ndim,f,lprt)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (nlog=0,nin=5,nout=6,nres=7)
      logical lprt
      real*8 x(ndim)
      common /minim/ econ,ncon,ntype
c
      real*8 kf,rho,acn,ast,atn,als,al2,als2,bst,btn,bls,
     &       cn,cne,dt,dr,evx,h2m,h2mcs,pi,s
      common /consts/ kf,rho,acn,ast,atn,als,al2,als2,bst,btn,bls,
     &       cn,cne,dt,dr,evx,h2m,h2mcs,pi,s

      real*8 aa(8),ab(8),ad(8,8),ae(6,2),af(8),ak(8,8,8),al(6,6,6),
     &       as(6),at(8,8),ax(6,6,6)
      common /amatrx/ aa,ab,ad,ae,af,ak,al,as,at,ax
      real*8 u,uf,up,tnia,tnic,tnis,tniu,tnix,cut,cut0,
     &       w3va,w3vc,w3vs,w3vu,w3vx
      common /tbcnst/ u,uf,up,tnia,tnic,tnis,tniu,tnix,cut,cut0,
     &       w3va,w3vc,w3vs,w3vu,w3vx
      real*8 eav,fsof,plm,qmin,qmax
      common /pionic/ eav,fsof,plm,qmin,qmax
      real*8 temp,mstar,chmpot,entrpy,ksav,kqav
      character*24 nffile
      common /files/ nffile
      common /hotted/ temp,mstar,chmpot,entrpy,ksav,kqav
      save nmlocal,np,nv,nt,ni,nie,nio,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
c     &,npi,npf,bst,btn,bls,nosave,npisav
     &,npi,npf,nosave,npisav
c
      real*8 gint(6)
      character*20 pname(30),tname(0:5),ptnnam
      character*20 timdat
      character*32 mname(4)
      character*50 sysdat
      character*160 fname
      character*7 outpre
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
      if (ndim.eq.2) then
        ast=x(2)
        atn=ast
        als=ast
      else if (ndim.eq.3) then
        ast=x(2)
        atn=x(3)
        als=ast
      else if (ndim.eq.4) then
        ast=x(2)
        atn=ast
        als=ast
        bst=x(3)
        btn=x(4)
        if (bls.ne.0.) bls=bst
      else if (ndim.eq.7) then
        ast=x(2)
        atn=x(3)
        als=x(4)
        bst=x(5)
        btn=x(6)
        bls=x(7)
      else if (ndim.eq.9) then
        ast=x(2)
        atn=x(3)
        als=x(4)
        al2=x(5)
        als2=x(6)
        bst=x(7)
        btn=x(8)
        bls=x(9)
      end if
      call nmmain(np,nv,nt,ni,nie,nio,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
        g2=g2+(gint(l)+1.)**2
      end do
      cons=econ*sqrt(g2)**ncon
      f=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
c $$$$$$$$$$$$$$$$$$$$$$$
c NEW stability condition
c $$$$$$$$$$$$$$$$$$$$$$$
      f=f+abs(endiff)/2
      no=0
#if defined CUSTOM_INPUTS && defined ONLY_NUCMAT
      write(nres,*) f
#ifndef DO_FULLX
     & ,(x(i),i=1,ndim)
#else
     &, (x(i),i=1,nbdirsmax)
#endif
#endif
      return
c ************************
c entry for initialization
c ************************
      entry nminit(x,ndim)
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
      read(nin,1001) nmlocal,np,nv,nt,ni,nie,nio,no,ns
 1001 format(5x,i3)
      read(nin,1001) lf,lc,ls,lt,ll,lg,le,l3,lk
      read(nin,1002) kf,rho,dor
 1002 format(5x,f10.4)
      read(nin,1002) acn,ast,atn,als,al2,als2
      read(nin,1002) bst,btn,bls,cn,cne
c      read(nin,1002) temp,mstar,tnia,tnic,tnis
c      read(nin,1002) tniu,tnix,cut,cut0
      read(nin,1002) temp,mstar,tnia,tnic,tnis,tniu,tnix,cut,cut0
      read(nin,1001) npi,npf
      read(nin,1002) eav,fsof,plm,qmin,qmax
      if (ns.eq.3) read(nin,1003) nffile
 1003 format(5x,a24)
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
        write(nlog,1012) ptnnam,tname(nt)
        write(nout,1012) ptnnam,tname(nt)
 1012   format(/4x,2a20)
      else if (np.gt.100 .and. np.le.200) then
        ptnnam(1:12)='Norfolk vij '
        if (nt.gt.100) then
          ptnnam(13:18)='+ Vijk'
          write(nlog,1014) ptnnam,np,nt
          write(nout,1014) ptnnam,np,nt
 1014     format(/4x,a20,'#',i3,' + #',i3)
        else
          write(nlog,1016) ptnnam,np
          write(nout,1016) ptnnam,np
 1016     format(/4x,a12,'#',i3)
        end if
      end if
      s=float(4/nmlocal)
            x(1)=dor
      if (ndim.eq.2) then
        x(2)=ast
      else if (ndim.eq.3) then
        x(2)=ast
        x(3)=atn
      else if (ndim.eq.4) then
        x(2)=ast
        x(3)=bst
        x(4)=btn
        if (bls.ne.0.) bls=bst
      else if (ndim.eq.7) then
        x(2)=ast
        x(3)=atn
        x(4)=als
        x(5)=bst
        x(6)=btn
        x(7)=bls
      else if (ndim.eq.9) then
        x(2)=ast
        x(3)=atn
        x(4)=als
        x(5)=al2
        x(6)=als2
        x(7)=bst
        x(8)=btn
        x(9)=bls
      else
        write(nlog,666)
        write(nout,666)
  666   format(/6x,'stopping because of non-standard ndim')
        stop 666
      end if
      nosave=no
      npisav=npi
      no=1
      npi=0
#ifdef CUSTOM_INPUTS
#if defined (CASE_SNM)
      outpre = "out_snm"
#else
      outpre = "out_pnm"
#endif
#ifndef DO_FULLX
      read(nin,*) (x(i),i=1,ndim)
#ifdef FD
      read(nin,*) nind, xptert
      if(nind.ne.0) x(nind)=x(nind)+x(nind)*xptert
#endif
      if (ndim.eq.2) then
      write(fname,"(A7,2(A1,F5.3),A1,F5.3,3(A1,I2),A4)")
     &outpre, ("_",abs(x(i)),i=1,ndim),"_",rho,"_",lc,
     & "_",ls,"_",lt,".txt"
      else if (ndim.eq.4) then
      write(fname,"(A7,4(A1,F5.3),A1,F5.3,3(A1,I2),A4)")
     &outpre, ("_",abs(x(i)),i=1,ndim),"_",rho,"_",lc,
     & "_",ls,"_",lt,".txt"
      else if (ndim.eq.7) then
      write(fname,"(A7,7(A1,F5.3),A1,F5.3,3(A1,I2),A4)")
     &outpre, ("_",abs(x(i)),i=1,ndim),"_",rho,"_",lc,
     & "_",ls,"_",lt,".txt"
      else if (ndim.eq.9) then
      write(fname,"(A7,9(A1,F5.3),A1,F5.3,3(A1,I2),A4)")
     &outpre, ("_",abs(x(i)),i=1,ndim),"_",rho,"_",lc,
     & "_",ls,"_",lt,".txt"
      end if
#else
      nbdirsmax=7
      read(nin,*) (x(i),i=1,nbdirsmax)
      write(fname,"(A7,7(A1,F19.17),A1,F19.17,3(A1,I2),A4)")
     &outpre, ("_",abs(x(i)),i=1,nbdirsmax),"_",rho,"_",lc,
     & "_",ls,"_",lt,".txt"
#endif
      open(unit=nres,file=fname,action="WRITE")
#endif
      return
c *******************
c entry for final run
c *******************
      entry nmfin(x,ndim)
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
      if (ndim.eq.2) then
        ast=x(2)
        atn=ast
        als=ast
      else if (ndim.eq.3) then
        ast=x(2)
        atn=x(3)
        als=ast
      else if (ndim.eq.4) then
        ast=x(2)
        atn=ast
        als=ast
        bst=x(3)
        btn=x(4)
        if (bls.ne.0.) bls=bst
      else if (ndim.eq.7) then
        ast=x(2)
        atn=x(3)
        als=x(4)
        bst=x(5)
        btn=x(6)
        bls=x(7)
      else if (ndim.eq.9) then
        ast=x(2)
        atn=x(3)
        als=x(4)
        al2=x(5)
        als2=x(6)
        bst=x(7)
        btn=x(8)
        bls=x(9)
      end if
      call nmmain(np,nv,nt,ni,nie,nio,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,npi,npf, gint, endiff, efree)
c     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do 995 l=1,2,nmlocal
  995 g2=g2+(gint(l)+1.)**2
      final=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      fplus=final+abs(endiff)/2
      write(nlog,*) final,fplus
      write(nout,*) final,fplus
 1095 format(/2x,'econs   eplus',/2f17.9)
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
#endif
