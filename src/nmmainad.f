c *id* nmmain **********************************************************
c main subroutine for calculating properties of nuclear or
c neutron matter
c **********************************************************************
      subroutine nmmainad(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,
     &                    lk,dor,bst,btn,bls,npi,npf, gint, endiff, 
     &                    efree,flocal,nmlocal)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter (nlog=0,nin=5,nout=6)
      common /minim/ econ,ncon,ntype
      real*8 kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,h2m,h2mcs,pi,s
      common /consts/ kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,
     &       h2m,h2mcs,pi,s
      real*8 :: dor,bst,btn,bls,endiff,efree,gint(6),flocal
      integer*4 :: np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
      integer*4 :: npi,npf,l1,l2, nmlocal
      real*8 :: g2
      integer*4 :: l
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
        g2=g2+(gint(l)+1.)**2
      end do 
      flocal=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
      end
