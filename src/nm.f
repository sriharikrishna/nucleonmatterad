c **********************************************************************
c nmprog
c program for nuclear/neutron matter
c **********************************************************************
      program nmprog
      use nmvar
      !implicit real*8 (a-h,o-z)
      !implicit integer*4 (i-n)
      external nucmat
      parameter (maxdim=15)
      character*50 sysdat
      character*20 timdat
      logical lprt
      real*8 x(maxdim),scale(maxdim)
      !common /minim/ econ,ncon,ntype
      character*4 etype(3)
      data etype/' jf ',' av ',' pb '/
      timeit=timer(0.)
      call header(sysdat,timdat)
      write(nlog,6) sysdat,timdat
      write(nout,6) sysdat,timdat
    6 format(/80('*')//1x,a50,8x,a20)
      read(nin,*) maxcl1,tol1,maxcl2,tol2,alpha,beta,gamma,n
      read(nin,*) (scale(i),i=1,n)
      read(nin,*) econ,ncon,ntype
      call nminit(x,n)
      if (maxcl1.ge.1) then
        write(nlog,7) etype(ntype+2),n
        write(nout,7) etype(ntype+2),n
    7   format(/4x,a4,'energy minimization in',i2,' dimensions')
        if (econ.ne.0.) then
          write(nlog,8) econ,ncon
          write(nout,8) econ,ncon
    8     format(5x,'with constraint',f6.0,
     &              '*sqrt(sum((1+gint(k))**2))**',i1)
        end if
        write(nlog,9) tol1,alpha,beta,gamma,tol2,(scale(i),i=1,n)
        write(nout,9) tol1,alpha,beta,gamma,tol2,(scale(i),i=1,n)
    9   format(5x,'simplex tolerance is',f7.3,' with coefficients of',
     &        /5x,'reflection ',f7.3,' contraction ',f7.3,
     &          ' expansion ',f7.3,
     &        /5x,'quadratic tolerance is',f7.3,
     &        /5x,'scale factors are',5f7.3)
        call minimi(maxcl1,tol1,maxcl2,tol2,alpha,beta,gamma,scale,x,
     &              fbest,nucmat,n)
      end if
      call nmfin(x,n)
      timeit=timer(timeit)
      write(nlog,999) timeit
      write(nout,999) timeit
  999 format(/5x,'job time =',f8.3,' seconds')
      stop
      end
c *id* nucmat **********************************************************
c subroutine for driving nuclear/neutron matter code
c ----------------------------------------------------------------------
      subroutine nucmat(x,n,flocal,lprt)
      use nmvar
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      logical lprt
      real*8 x(n)
      !common /minim/ econ,ncon,ntype
c
      !real*8 kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,h2m,h2mcs,pi,s
      !common /consts/ kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,
      !&       h2m,h2mcs,pi,s
      !real*8 aa(8),ab(8),ad(8,8),ae(6,2),af(8),ak(8,8,8),al(6,6,6),
      !&      as(6),at(8,8),ax(6,6,6)
      !common /amatrx/ aa,ab,ad,ae,af,ak,al,as,at,ax
      !real*8 u,uf,up,tnia,tnic,tniu,tnix,cut,cut0,w3v0,w3v1,w3va,w3vc
      !common /tbcnst/ u,uf,up,
      !&       tnia,tnic,tniu,tnix,cut,cut0,w3v0,w3v1,w3va,w3vc
      !real*8 eav,fsof,plm,qmin,qmax
      !common /pionic/ eav,fsof,plm,qmin,qmax
      !real*8 temp,mstar,chmpot,entrpy,ksav,kqav
      !common /hotted/ temp,mstar,chmpot,entrpy,ksav,kqav
      save nmlocal,np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt
     &,ll,lg,le,l3,lk,npi,npf
     &,bst,btn,bls,nosave,npisav
c
      real*8 gint(6)
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
        bls=bst
      end if
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do 5 l=1,2,nmlocal
    5 g2=g2+(gint(l)+1.)**2
      flocal=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      no=0
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
     &lc,ls,lt,ll,lg,le,l3,lk
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
        bls=bst
      end if
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do 995 l=1,2,nmlocal
  995 g2=g2+(gint(l)+1.)**2
      final=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      fplus=final+abs(endiff)
      write(nlog,1095) final,fplus
      write(nout,1095) final,fplus
 1095 format(/2x,'final   fplus',/2f8.3)
      if (no.le.1) go to 999
      call nmout(le,lg,lt,l3,nie,no,nt,nv)
  999 return
      end
