c *id* nmmain **********************************************************
c main subroutine for calculating properties of nuclear or
c neutron matter
c **********************************************************************
      subroutine nmmainad(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,
     &                    lk,dor,bst,btn,bls,npi,npf, gint, endiff, 
     &                    efree,flocal,nmlocal)
      use nmvar
      use nmtbimod
      use nmhncmod
      implicit none
c ----------------------------------------------------------------------
      real*8 :: dor,bst,btn,bls,endiff,efree,gint(6),flocal
      integer*4 :: np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
      integer*4 :: npi,npf,l1,l2, nmlocal
      real*8 :: g2
      integer*4 :: l
!$openad INDEPENDENT(dor)
!$openad INDEPENDENT(bst)
!$openad INDEPENDENT(btn)
!$openad INDEPENDENT(bls)
!$openad INDEPENDENT(ast)
!$openad INDEPENDENT(atn)
!$openad INDEPENDENT(als)
      call nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &           ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      g2=0.
      do l=1,2,nmlocal
        g2=g2+(gint(l)+1.)**2
      end do 
      flocal=efree+ntype*endiff/2+econ*sqrt(g2)**ncon
! DEPENDENT(endiff)
! DEPENDENT(efree)
! DEPENDENT(gint)
!$openad DEPENDENT(flocal)
      end
