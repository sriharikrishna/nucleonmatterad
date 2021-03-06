c **********************************************************************
c nmprog
c program for nuclear/neutron matter
c **********************************************************************
      program nmprog
      use nmvar
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      external nucmat
      parameter (maxdim=15)
      external nmvarinit
      character*50 sysdat
      character*20 timdat
      logical lprt
      real*8 x(maxdim),scale(maxdim)
      !common /minim/ econ,ncon,ntype
      character*4 etype(3)
      data etype/' jf ',' av ',' pb '/
      call processargs(argval)
      timeit=timer(0.)
      call nmvarinit()
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
c *id* processargs *****************************************************
c subroutine for processing command line arguments for OpenAD Case
c ----------------------------------------------------------------------
      subroutine processargs(argval)
      implicit none
      character(len=32)  :: argval
      integer argpos, totarg
      totarg = iargc()  
      if (totarg .ne. 1) then
        stop "ERROR: This program takes only one option"
      end if
      DO argpos = 1, totarg
        CALL getarg(argpos, argval)
      END DO
      end subroutine processargs
c *id* nucmat **********************************************************
c subroutine for driving nuclear/neutron matter code
c ----------------------------------------------------------------------
      subroutine nucmat(x,n,flocal,lprt)
      use nmvar
#ifndef ALLOW_OPENAD
      use nmvarcopypassive
#endif
#ifdef ALLOW_OPENAD
      use OAD_tape
      use OAD_rev
      use OAD_active
      use OAD_cp
#endif
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      logical lprt
      real*8 x(n)
      real*8 h
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
#ifdef ALLOW_OPENAD
      type(active) :: bst,btn,bls,dor
      type(active) :: gint(6),final,flocaload,endiff
#else
      !real*8:: bst,btn,bls,dor
      !real*8 :: final,flocal,endiff

#endif
      save nmlocal,np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt
     &,ll,lg,le,l3,lk,npi,npf
     &,bst,btn,bls,nosave,npisav
c
#ifndef ALLOW_OPENAD
      real*8 gint(6)
      real*8 dor_d, bst_d, btn_d, bls_d, ast_d, atn_d, als_d
      real*8 tmp_dor, tmp_bst, tmp_btn, tmp_bls
      real*8 tmp_ast, tmp_atn, tmp_als
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
      
#ifndef ALLOW_OPENAD
      dor=x(1)
#else
      dor%v=x(1)
      flocaload%v = flocal 
#endif
      if (n.ge.2) then
#ifndef ALLOW_OPENAD
        ast=x(2)
        atn=ast
        als=ast
#else
        ast%v=x(2)
        atn%v=ast%v
        als%v=ast%v
#endif
      end if
      if (n.eq.4) then
#ifndef ALLOW_OPENAD
        bst=x(3)
        btn=x(4)
        bls=bst
#else
        bst%v=x(3)
        btn%v=x(4)
        bls%v=bst%v
#endif
      end if
#ifndef ALLOW_OPENAD
      if (argval .eq. "p" .or. argval .eq. "") then
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do 5 l=1,2,nmlocal
    5 g2=g2+(gint(l)+1.)**2
      flocal=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      else if (argval .eq. "f") then
      write(*,*) "dor%v", dor, "bls%v", bls
      h = 0.0000001
      tmp_dor=dor    
      tmp_bst=bst      
      tmp_btn=btn      
      tmp_bls=bls      
      tmp_ast=ast      
      tmp_atn=atn      
      tmp_als=als      
      
      call var_transfer_store()
      
!! dor
      call var_transfer_restore()
      dor = dor + dor * h 
      !dor = dor + h*10000
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
      g2=g2+(gint(l)+1.)**2
      end do
      dor_d=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      dor=tmp_dor    
      bst=tmp_bst      
      btn=tmp_btn      
      bls=tmp_bls      
      ast=tmp_ast      
      atn=tmp_atn      
      als=tmp_als      
!!
!! bst
      call var_transfer_restore()
      bst = bst + bst * h
      !bst = bst +  h
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
      g2=g2+(gint(l)+1.)**2
      end do
      bst_d=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      dor=tmp_dor    
      bst=tmp_bst      
      btn=tmp_btn      
      bls=tmp_bls      
      ast=tmp_ast      
      atn=tmp_atn      
      als=tmp_als      
!!
!! btn
      call var_transfer_restore()
      btn = btn + btn * h
      !btn = btn + h
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
      g2=g2+(gint(l)+1.)**2
      end do
      btn_d=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      dor=tmp_dor    
      bst=tmp_bst      
      btn=tmp_btn      
      bls=tmp_bls      
      ast=tmp_ast      
      atn=tmp_atn      
      als=tmp_als      
!!
!! bls
      call var_transfer_restore()
      if(bls.ne.0.0) then
        bls = bls + bls * h
      else
        bls = bls+ h
      end if 
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
      g2=g2+(gint(l)+1.)**2
      end do
      bls_d=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      dor=tmp_dor    
      bst=tmp_bst      
      btn=tmp_btn      
      bls=tmp_bls      
      ast=tmp_ast      
      atn=tmp_atn      
      als=tmp_als      
!!
!! ast
      call var_transfer_restore()
      ast = ast + ast * h
      !ast = ast + h
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
      g2=g2+(gint(l)+1.)**2
      end do
      ast_d=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      dor=tmp_dor    
      bst=tmp_bst      
      btn=tmp_btn      
      bls=tmp_bls      
      ast=tmp_ast      
      atn=tmp_atn      
      als=tmp_als      

!!
!! atn
      call var_transfer_restore()
      atn = atn + atn * h
      !atn = atn + h
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
      g2=g2+(gint(l)+1.)**2
      end do
      atn_d=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      dor=tmp_dor    
      bst=tmp_bst      
      btn=tmp_btn      
      bls=tmp_bls      
      ast=tmp_ast      
      atn=tmp_atn      
      als=tmp_als      

!!
!! als
      call var_transfer_restore()
      als = als + als * h
      !als = als + h
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
      g2=g2+(gint(l)+1.)**2
      end do
      als_d=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      bst=tmp_bst      
      btn=tmp_btn      
      bls=tmp_bls      
      ast=tmp_ast      
      atn=tmp_atn      
      als=tmp_als      
!!
      call var_transfer_restore()
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
      g2=g2+(gint(l)+1.)**2
      end do
      flocal=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      dor_d=(dor_d-flocal)/(tmp_dor*h)      
      bst_d=(bst_d-flocal)/(tmp_bst*h)      
      btn_d=(btn_d-flocal)/(tmp_btn*h)     
      if(tmp_bls.ne.0.0) then  
        bls_d=(bls_d-flocal)/(tmp_bls*h)      
      else
        bls_d=(bls_d-flocal)/(h)      
      end if
      ast_d=(ast_d-flocal)/(tmp_ast*h)      
      atn_d=(atn_d-flocal)/(tmp_atn*h)      
      als_d=(als_d-flocal)/(tmp_als*h)      
!!
      write(*,*) "dor%d", dor_d
      write(*,*) "bst%d", bst_d
      write(*,*) "btn%d", btn_d
      write(*,*) "bls%d", bls_d
      write(*,*) "ast%d", ast_d
      write(*,*) "atn%d", atn_d
      write(*,*) "als%d", als_d
      else
        stop ("ERROR : Argument must be 'p' or 'f'")
      end if
#else

      if (argval .eq. "p") then
      our_rev_mode%plain=.TRUE.
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%tape=.FALSE.
      our_rev_mode%adjoint=.FALSE.
      our_rev_mode%topsplit=.FALSE.
      call nmmainad(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree,flocaload
     &           ,nmlocal)
      flocal=flocaload%v
#ifdef ALLOW_OPENAD_FORWARD
      else if (argval .eq. "f") then
      call flush(6)
      dor%d(1) = 1.0
      bst%d(2) = 1.0
      btn%d(3) = 1.0
      bls%d(4) = 1.0
      ast%d(5) = 1.0
      atn%d(6) = 1.0
      als%d(7) = 1.0
      flocaload%d = 0.0
      call nmmainad(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree,flocaload
     &           ,nmlocal)
      flocal=flocaload%v
      write(*,*) "dor%d", flocaload%d(1)
      write(*,*) "bst%d", flocaload%d(2)
      write(*,*) "btn%d", flocaload%d(3)
      write(*,*) "bls%d", flocaload%d(4)
      write(*,*) "ast%d", flocaload%d(5)
      write(*,*) "atn%d", flocaload%d(6)
      write(*,*) "als%d", flocaload%d(7)
      call flush(6)
#else
      else if (argval .eq. "a") then
      !our_rev_mode%plain=.TRUE.
      !our_rev_mode%arg_store=.TRUE.
      our_rev_mode%plain=.FALSE.
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%tape=.TRUE.
      our_rev_mode%adjoint=.FALSE.
      our_rev_mode%topsplit=.FALSE.
      call flush(6)
      call cp_init()
      call oad_tape_init()
      flocaload%v=flocal
      call nmmainad(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree,flocaload
     &           ,nmlocal)
      flocal=flocaload%v
      call flush(6)
      dor%d = 0.0
      bst%d = 0.0
      btn%d = 0.0
      bls%d = 0.0
      ast%d = 0.0
      atn%d = 0.0
      als%d = 0.0
      flocaload%d = 1.0
      our_rev_mode%plain=.FALSE.
      our_rev_mode%arg_store=.FALSE.
      !our_rev_mode%arg_restore=.TRUE.
      !our_rev_mode%tape=.TRUE.
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%tape=.FALSE.
      !our_rev_mode%adjoint=.FALSE.
      our_rev_mode%adjoint=.TRUE.
      our_rev_mode%topsplit=.FALSE.
      call nmmainad(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree,flocaload
     &           ,nmlocal)
      write(*,*) "dor%d", dor%d
      write(*,*) "bst%d", bst%d
      write(*,*) "btn%d", btn%d
      write(*,*) "bls%d", bls%d
      write(*,*) "ast%d", ast%d
      write(*,*) "atn%d", atn%d
      write(*,*) "als%d", als%d
      call flush(6)
      call oad_tape_delete()
#endif
      else
        stop ("ERROR : Argument must be 'a' or 'p'")
      end if
#endif
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
#ifndef ALLOW_OPENAD
      read(nin,1002) kf,rho,dor,acn,ast,atn,als,bst,btn,bls,cn,cne
#else
      read(nin,1002) kf,rho,dor%v,acn,ast%v,atn%v,als%v,
     &bst%v,btn%v,bls%v,cn,cne
#endif
 1002 format(5x,f10.4)
#ifndef ALLOW_OPENAD
      read(nin,1002) temp,mstar,tnia,tnic,tniu,tnix,cut,cut0
#else
      read(nin,1002) temp%v,mstar,tnia,tnic,tniu,tnix,cut,cut0
#endif
      read(nin,1001) npi,npf
      read(nin,1002) eav,fsof,plm,qmin,qmax
      write(nlog,1010) mname(nmlocal)
      write(nout,1010) mname(nmlocal)
 1010 format(/4x,a32)
#ifndef ALLOW_OPENAD
      if (temp.gt.0.) then
#else
      if (temp%v.gt.0.) then
#endif
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
#ifndef ALLOW_OPENAD
      x(1)=dor
      if (n.ge.2) x(2)=ast
#else
      x(1)=dor%v
      if (n.ge.2) x(2)=ast%v
#endif
      if (n.eq.4) then
#ifndef ALLOW_OPENAD
        x(3)=bst
        x(4)=btn
#else
        x(3)=bst%v
        x(4)=btn%v
#endif
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
#ifndef ALLOW_OPENAD
      dor=x(1)
#else
      dor%v=x(1)
#endif
      if (n.ge.2) then
#ifndef ALLOW_OPENAD
        ast=x(2)
        atn=ast
        als=ast
#else
        ast%v=x(2)
        atn%v=ast%v
        als%v=ast%v
#endif
      end if
      if (n.eq.4) then
#ifndef ALLOW_OPENAD
        bst=x(3)
        btn=x(4)
        bls=bst
#else
        bst%v=x(3)
        btn%v=x(4)
        bls%v=bst%v
#endif
      end if
#ifndef ALLOW_OPENAD
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do 995 l=1,2,nmlocal
  995 g2=g2+(gint(l)+1.)**2
      final=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      fplus=final+abs(endiff)
#else
      our_rev_mode%plain=.TRUE.
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%tape=.FALSE.
      our_rev_mode%adjoint=.FALSE.
      our_rev_mode%topsplit=.FALSE.
      call nmmainad(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree
     &           ,final,nmlocal)
      fplus=final%v+abs(endiff%v)
#endif
#ifndef ALLOW_OPENAD
      write(nlog,1095) final,fplus
      write(nout,1095) final,fplus
#else
      write(nlog,1095) final%v,fplus
      write(nout,1095) final%v,fplus
#endif
 1095 format(/2x,'final   fplus',/2f8.3)
      if (no.le.1) go to 999
      call nmout(le,lg,lt,l3,nie,no,nt,nv)
  999 return
      end

      subroutine nmvarinit
      use params
      use nmvar

      !common consts
      rho=0.
      kf=0.

      aa = (/1.,3.,3.,9.,6.,18.,1.,3./)
      ab = (/1.,3.,3.,9.,1.,3.,1.,3./)
      ad = reshape(
     &          [0., 0., 0., 0., 0., 0.,0., 0.,
     &           0.,12., 0.,12., 0.,12.,0.,12.,
     &           0., 0.,12.,12.,12.,12.,6., 6.,
     &           0.,12.,12., 8.,12., 8.,6.,10.,
     &           0., 0.,12.,12.,12.,12.,6., 6.,
     &           0.,12.,12., 8.,12., 8.,6.,10.,
     &           0., 0., 6., 6., 6., 6.,3., 3.,
     &           0.,12., 6.,10., 6.,10.,3.,11.],[8,8])
      ae = reshape([0.,0.
     &,6.,6.,6.,6.,0.,0.,6.,-2.,6.,-2.],[6,2])
      af = (/1.,1.,1.,1.,0.,0.,0.,0./)
      ak = reshape(
     & [1.,0.,0.,0.,0.,0.,0.,0.,
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
      al = reshape([1.,0.,0.,0.,0.,0.,
     & 0.,3.,0.,0.,0.,0.,
     & 0.,0.,3.,0.,0.,0.,
     & 0.,0.,0.,9.,0.,0.,
     & 0.,0.,0.,0.,6.,0.,
     & 0.,0.,0.,0.,0.,18.,
     & 0.,3.,0.,0.,0.,0.,
     & 3.,6.,0.,0.,0.,0.,
     & 0.,0.,0.,9.,0.,0.,
     & 0.,0.,9.,18.,0.,0.,
     & 0.,0.,0.,0.,0.,18.,
     & 0.,0.,0.,0.,18.,36.,
     & 0.,0.,3.,0.,0.,0.,
     & 0.,0.,0.,9.,0.,0.,
     & 3.,0.,6.,0.,0.,0.,
     & 0.,9.,0.,18.,0.,0.,
     & 0.,0.,0.,0.,-6.,0.,
     & 0.,0.,0.,0.,0.,-18.,
     & 0.,0.,0.,9.,0.,0.,
     & 0.,0.,9.,18.,0.,0.,
     & 0.,9.,0.,18.,0.,0.,
     & 9.,18.,18.,36.,0.,0.,
     & 0.,0.,0.,0.,0.,-18.,
     & 0.,0.,0.,0.,-18.,-36.,
     & 0.,0.,0.,0.,6.,0.,
     & 0.,0.,0.,0.,0.,18.,
     & 0.,0.,0.,0.,-6.,0.,
     & 0.,0.,0.,0.,0.,-18.,
     & 6.,0.,-6.,0.,12.,0.,
     & 0.,18.,0.,-18.,0.,36.,
     & 0.,0.,0.,0.,0.,18.,
     & 0.,0.,0.,0.,18.,36.,
     & 0.,0.,0.,0.,0.,-18.,
     & 0.,0.,0.,0.,-18.,-36.,
     & 0.,18.,0.,-18.,0.,36.,
     & 18.,36.,-18.,-36.,36.,72.],[6,6,6])
      as = (/2.25,2.25,1.25,1.25,1.,1./)
      at = reshape([1.,0.,0.,0.,0.,0.,0.,0.,
     & 0.,1.,0.,0.,0.,0.,0.,0.,
     & 0.,0.,1.,0.,0.,0.,0.,0.,
     & 0.,0.,0.,1.,0.,0.,0.,0.,
     & 0.,0.,0.,0.,1.,0.,0.,0.,
     & 0.,0.,0.,0.,0.,1.,0.,0.,
     & 0.,0.,0.,0.,0.,0.,1.,0.,
     & 0.,0.,0.,0.,0.,0.,0.,1.],[8,8])
      ax = reshape([
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

#ifndef ALLOW_OPENAD
       entrpy = 0.
#else
       entrpy%v = 0.
#endif
      end subroutine nmvarinit
