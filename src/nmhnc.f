      module nmhncmod
      use nmvar
      use nmsubmod
      implicit none


      private
      real*8 afe(6),afi(6,6,6),afj(6,6,6),afk(6,6,6),ahi(6,6,6)
     &,ahj(6,6,6),ahk(6,6,6),vco(6,3,3),vce(6,3),xcc(lgrid)
     &,xca(lgrid,6),xdd(lgrid,6),xde(lgrid,6),xee(lgrid,6),xgcc(lgrid)
     &,xgcb(lgrid,6),xgdd(lgrid,6),xgde(lgrid,6),xgee(lgrid,6)
     &,xgca(lgrid,6)
     &,qcc(6),qdd(6),qde(6),qee(6),pcc(6),pdd(6),pde(6),pee(6)
     &,zca(6),zcb(6),zdd(6),zde(6),zee(6),v3st(6)
      integer lit(6),mit(6),nit(6)

      real*8 :: eijk
      public :: nmhnc
      
      contains

c *id* nmhnc ***********************************************************
c subroutine for fhnc/soc equations
c **********************************************************************
      subroutine nmhnc(lg,le,l3,ni,nie,no,nt,nv)
      implicit none
      integer*4, parameter :: n3s=5-nm
      integer*4, parameter :: n3t=7-nm
      integer*4, parameter :: nphi=5
      integer*4, parameter :: lx=512
      integer*4 lg,le,l3,ni,nie,no,nt,nv
c
      real*8 :: xa(4),ya(4),yb(4),yc(4),jx(0:lx),rsx(0:lx)
     &,sddx(0:lx),sdex(0:lx),seex(0:lx),cphi(nphi)
     &,sccr(lgrid),sddr(lgrid),sder(lgrid),seer(lgrid)
      real*8 :: coe,dphi,q3,q3cne,q4,dx,dxi,co,q1,q2,q2cn,fc2
      real*8 :: c3,s3,s4,yddd,ydde,dya,dyb,dyc,ydee,yeee,yccd,ycce,cg
      real*8 :: s3s4,rdddde,rdeeee,reeeee,rddddd,rdeede,reeede,rcccc,dln
      real*8 :: x1,x2,x3,x4,x5,x6,axi,axj,axk,rde,rdd,ree,rcut0,xlm,xlm2
      real*8 :: y1,y2,y3,y4,y5,y6,z1,z2,z3,z4,z5,z6,xmu,pa,pb,pc,rcut
      integer*4 :: lmax,kji,ka,kb,ljk,ijk,kl,it
      real*8 :: a,b,c4,c3c4,r34s,xp,xm,sddp,sdep,seep,vpd,vpe,vpp,vfd
      real*8 :: vfe,vfp,ffl,hl,fl2,hlh,ei,eij,x,y,z,xijk,xi,xj,xk,rxmu
      real*8 :: rmu,rlm,ermu,erlm,pac,pap
      integer*4 :: i,j,k,l,m,n,il,ia,jn,ma,mb,mli,np,ikj,jik,lj,dlk,dkn

c -------------------
c  data
c -------------------
      lit = [0,2,2,4,4,4]
      mit = [0,2,4,2,4,4]
      nit = [0,2,4,4,2,4]
c -------------------
c statement functions
c -------------------
c ---------------------
c setup sines & cosines
c ---------------------
      lmax=max(lg,le,l3)
      kji=0
      do 15 i=1,lmax
        do 10 j=1,lmax
          ka=iabs(i-j)+1
          kb=min(i+j-1,lgrid)
          do 5 k=ka,kb
            kji=kji+1
            index(k,j,i)=kji
            xtheta(kji)=(rs(i)+rs(j)-rs(k))/(2*r(i)*r(j))
            ytheta(kji)=(rs(i)+rs(k)-rs(j))/(2*r(i)*r(k))
            ztheta(kji)=(rs(j)+rs(k)-rs(i))/(2*r(j)*r(k))
            stheta(kji)=sqrt(1-xtheta(kji)**2)
    5     continue
   10   continue
   15 continue
c ---------------------------------------
c set up fhnc/4 integration parameters
c interpolation grid and index for r34**2
c ---------------------------------------
      if (nie.ge.1) then
        coe=1-cne
        dphi=pi/nphi
        q3=2*pi*rho*dr*dr
        q3cne=q3*cne
        q4=2*rho*dr*dr*dphi
        do 20 n=1,nphi
          cphi(n)=cos((float(n)-.5)*dphi)
   20   continue
        dx=((le+lg)*dr)**2/lx
        dxi=1/dx
        do 25 i=0,lx
          rsx(i)=i*dx
          call locate(rs,lmax,rsx(i),j)
          if (j.le.1) j=2
          if (j.ge.lmax-1) j=lmax-2
          jx(i)=j
   25   continue
      end if
c -------------
c zero matrices
c -------------
      do 30 i=1,lgrid
        gx(i)=1
        gy(i)=0
        gz(i)=0
        gl(i)=sl(i)
   30 continue
      do 35 i=1,legrid
        sccd(i)=0
        scce(i)=0
        sddd(i)=0
        sdde(i)=0
        sdee(i)=0
        seee(i)=0
   35 continue
      do 40 il=1,lgrid*6
        gca(il,1)=0
        gcb(il,1)=0
        gdd(il,1)=0
        gde(il,1)=0
        gee(il,1)=0
        eca(il,1)=0
        ecb(il,1)=0
        edd(il,1)=0
        ede(il,1)=0
        eee(il,1)=0
   40 continue
      v3cc(:,:,:)=0
      v3dd(:,:,:)=0
      v3de(:,:,:)=0
      v3ee(:,:,:)=0
      do 50 il=1,lgrid*3
        bcc(il,1)=0
        bde(il,1)=0
   50 continue
      do 55 l=1,6
        afe(l)=acex(l,1,l)
   55 continue
      do 60 ljk=1,6*3*3
        vc(ljk,1,1)=1
   60 continue
c -------------------------------------
c main iteration loop for hnc equations
c -------------------------------------
      co=1-cn
      q1=4*pi*rho*dr
      q2=2*pi*rho*dr*dr
      q2cn=q2*cn
c =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>
      do 490 ia=1,ni
c ---------------
c input functions
c ---------------
      do 70 i=1,lmax
        fc2=f(i,1)**2
        xgdd(i,1)=fc2*gx(i)-1
        xgde(i,1)=fc2*gy(i)*gx(i)
        xgee(i,1)=fc2*(gy(i)**2+gz(i)-gl(i)**2/nu)*gx(i)
        xgcb(i,1)=-fc2*gx(i)*gl(i)/nu
        xgcc(i)=xgcb(i,1)
        xgca(i,1)=xgcc(i)
        xca(i,1)=0
        xcc(i)=xgcc(i)+gl(i)/nu
        sddr(i)=xgdd(i,1)*r(i)
        sder(i)=xgde(i,1)*r(i)
        seer(i)=xgee(i,1)*r(i)
        sccr(i)=xgcb(i,1)*r(i)
   70 continue
c --------------------------------------------
c skip elementary diagrams in first iterations
c --------------------------------------------
      if (ia.le.(ni-nie)) go to 200
c ----------------------------------------
c construct interpolated table of sxx(r34)
c ----------------------------------------
      do 90 i=0,lx
        do 80 n=1,4
          jn=jx(i)-2+n
          xa(n)=rs(jn)
          ya(n)=xgdd(jn,1)
          yb(n)=xgde(jn,1)
          yc(n)=xgee(jn,1)
   80   continue
        call polint(xa,ya,4,rsx(i),sddx(i),dya)
        call polint(xa,yb,4,rsx(i),sdex(i),dyb)
        call polint(xa,yc,4,rsx(i),seex(i),dyc)
   90 continue
c ----------------------------
c construct 3-point superbonds
c ----------------------------
!$OMP PARALLEL SHARED(sddd,sdde,sdee,seee,sccd,scce)
!$OMP& PRIVATE(i,j,k,l,m,np)
!$OMP DO SCHEDULE(DYNAMIC)
      do 160 i=1,le                              ! r12
        do 150 j=1,le                            ! r13
          ka=iabs(i-j)+1
          kb=min(i+j-1,le)
          do 140 k=ka,kb                         ! r23
            kji=index(k,j,i)
            c3=xtheta(kji)                       ! cos(theta3)
            s3=stheta(kji)                       ! sin(theta3)
            yddd=0.
            ydde=0.
            ydee=0.
            yeee=0.
            yccd=0.
            ycce=0.
            do 130 l=1,lg                        ! r14
              a=rs(j)+rs(l)
              b=2*r(j)*r(l)
              ma=iabs(i-l)+1
              mb=min(i+l-1,lg)
              do 120 m=ma,mb                     ! r24
                mli=index(m,l,i)
                c4=xtheta(mli)                   ! cos(theta4)
                s4=stheta(mli)                   ! sin(theta4)
                c3c4=c3*c4
                s3s4=s3*s4
                rdddde=sddr(l)*sddr(m)
                rdeeee=sddr(l)*sder(m)
                reeeee=sder(l)*sder(m)
                rddddd=rdddde+rdeeee+sder(l)*sddr(m)
                rdeede=rdeeee+reeeee+sddr(l)*seer(m)
                reeede=reeeee+(sder(l)*seer(m)+seer(l)*sder(m))
                rcccc=sccr(l)*sccr(m)
                do 110 np=1,nphi                 ! phi4
                  cg=c3c4+s3s4*cphi(np)          ! cos(gamma)
                  r34s=a-b*cg                    ! r34**2
                  n=int(r34s*dxi)
                  xp=(r34s-rsx(n))*dxi
                  xm=(rsx(n+1)-r34s)*dxi
                  sddp=xm*sddx(n)+xp*sddx(n+1)
                  sdep=xm*sdex(n)+xp*sdex(n+1)
                  seep=xm*seex(n)+xp*seex(n+1)
                  yddd=yddd+rddddd*sddp+rdddde*sdep
                  ydde=ydde+rddddd*sdep+rdddde*seep
                  ydee=ydee+rdeede*sdep+rdeeee*seep
                  yeee=yeee+reeede*sdep+reeeee*seep
                  yccd=yccd+rcccc*sddp
                  ycce=ycce+rcccc*sdep
  110           continue
  120         continue
  130       continue
            ijk=index(i,j,k)
            sddd(ijk)=q4*ri(i)*yddd
            sdde(ijk)=q4*ri(i)*ydde
            sdee(ijk)=q4*ri(i)*ydee
            seee(ijk)=q4*ri(i)*yeee
            sccd(ijk)=q4*ri(i)*yccd
            scce(ijk)=q4*ri(i)*ycce
  140     continue
  150   continue
  160 continue
!$OMP END DO NOWAIT
!$OMP END PARALLEL
c -----------------------
c set soc input functions
c -----------------------
  200 do 220 l=1+nm,6,nm
        vpd=vc(l,2,1)
        vpe=vc(l,2,2)
        vpp=vc(l,2,3)
        vfd=vc(l,3,1)
        vfe=vc(l,3,2)
        vfp=vc(l,3,3)
        do 210 i=1,lmax
          fc2=f(i,1)**2
          ffl=2*f(i,1)*f(i,l)
          hl=ffl+fc2*(gdd(i,l)+edd(i,l))
          xgca(i,l)=(-hl*vpe*gl(i)/nu+fc2*(gca(i,l)+eca(i,l)))*gx(i)
          xgcb(i,l)=(-hl*vpe*gl(i)/nu+fc2*(gcb(i,l)+ecb(i,l)))*gx(i)
          xgdd(i,l)=hl*gx(i)
          xgde(i,l)=(hl*gy(i)+fc2*(gde(i,l)+ede(i,l)))*gx(i)
          xgee(i,l)=(hl*(gy(i)**2+gz(i))+fc2*(gee(i,l)+eee(i,l)
     &     -af(l)*gl(i)**2/nu+2*(gde(i,l)+ede(i,l))*gy(i)))*gx(i)
          xca(i,l)=xgca(i,l)-gca(i,l)
          xdd(i,l)=xgdd(i,l)-gdd(i,l)
          xde(i,l)=xgde(i,l)-gde(i,l)
          xee(i,l)=xgee(i,l)-gee(i,l)
c -----------------------
          fl2=f(i,l)**2
          hlh=ffl+.5*fc2*gdd(i,l)
          xgdd(i,1)=xgdd(i,1)
     &     +aa(l)*(fl2*vfd**2+hlh*gdd(i,l)*vpd**2)*gx(i)
          xgde(i,1)=xgde(i,1)
     &     +aa(l)*((fl2*vfd*vfp+hlh*gdd(i,l)*vpd*vpp)*gy(i)
     &     +hl*gde(i,l)*vpd*vpe)*gx(i)
          xgee(i,1)=xgee(i,1)
     &     +aa(l)*((fl2*vfp**2+hlh*gdd(i,l)*vpp**2)*(gy(i)**2
     &     +gz(i))+hl*((gee(i,l)-af(l)*gl(i)**2/nu)*vpe**2
     &     +2*gde(i,l)*gy(i)*vpe*vpp)+fc2*(gde(i,l)**2*vpe**2
     &     +2*af(l)*(gca(i,l)+gcb(i,l))*gl(i)*vpe))*gx(i)
     &     -afe(l)*fl2*gl(i)**2*gx(i)*vfe**2/nu
          xgcb(i,1)=xgcb(i,1)+aa(l)*af(l)*vpe*fc2*gx(i)*gca(i,l)
          xgca(i,1)=xgca(i,1)+aa(l)*af(l)*vpe*(fc2*gx(i)-1)*gcb(i,l)
          xca(i,1)=xca(i,1)+aa(l)*af(l)*vpe*(fc2*gx(i)-1)*gcb(i,l)
  210   continue
  220 continue
      do 230 i=1,lmax
        xdd(i,1)=xgdd(i,1)-gdd(i,1)
        xde(i,1)=xgde(i,1)-gde(i,1)
        xee(i,1)=xgee(i,1)-gee(i,1)
  230 continue
c --------------------------------
c three-body fhnc/soc integrations
c --------------------------------
      do 240 il=1,lgrid*6
        gdd(il,1)=gdd(il,1)*co
        gde(il,1)=gde(il,1)*co
        gee(il,1)=gee(il,1)*co
        gca(il,1)=gca(il,1)*co
        gcb(il,1)=gcb(il,1)*co
  240 continue
!$OMP PARALLEL SHARED(gdd,gde,gee,gca,gcb)
!$OMP& PRIVATE(i,j,k,l)
!$OMP DO SCHEDULE(DYNAMIC)
      do 290 i=1,lg
        ei=q2cn*ri(i)
        do 280 j=1,lg
          eij=ei*r(j)
          ka=iabs(i-j)+1
          kb=min(i+j-1,lg)
          do 270 l=1,2,nm
            do 260 k=ka,kb
              eijk=eij*r(k)
              kji=index(k,j,i)
              x=xtheta(kji)
              y=ytheta(kji)
              z=ztheta(kji)
              xijk=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
              xi=3*z**2-1
              xj=1.5*y**2-.5
              xk=1.5*x**2-.5
              gdd(i,l)=gdd(i,l)+sdd(j,k,l,l)
              gde(i,l)=gde(i,l)+sde(j,k,l,l)
              gee(i,l)=gee(i,l)+see(j,k,l,l)
              gca(i,l)=gca(i,l)+eijk*xca(j,l)*xgcc(k)
              gcb(i,l)=gcb(i,l)+eijk*xcc(j)*xgcb(k,l)
              gdd(i,l+2)=gdd(i,l+2)+sdd(j,k,l+2,l+2)+xi*sdd(j,k,l+4,l+4)
              gde(i,l+2)=gde(i,l+2)+sde(j,k,l+2,l+2)+xi*sde(j,k,l+4,l+4)
              gee(i,l+2)=gee(i,l+2)+see(j,k,l+2,l+2)+xi*see(j,k,l+4,l+4)
              gca(i,l+2)=gca(i,l+2)+eijk*xca(j,l+2)*xgcc(k)
              gcb(i,l+2)=gcb(i,l+2)+eijk*xcc(j)*xgcb(k,l+2)
              gdd(i,l+4)=gdd(i,l+4)+xijk*sdd(j,k,l+4,l+4)
     &                  +xk*sdd(j,k,l+4,l+2)+xj*sdd(j,k,l+2,l+4)
              gde(i,l+4)=gde(i,l+4)+xijk*sde(j,k,l+4,l+4)
     &                  +xk*sde(j,k,l+4,l+2)+xj*sde(j,k,l+2,l+4)
              gee(i,l+4)=gee(i,l+4)+xijk*see(j,k,l+4,l+4)
     &                  +xk*see(j,k,l+4,l+2)+xj*see(j,k,l+2,l+4)
              gca(i,l+4)=gca(i,l+4)+eijk*xk*xca(j,l+4)*xgcc(k)
              gcb(i,l+4)=gcb(i,l+4)+eijk*xj*xcc(j)*xgcb(k,l+4)
  260       continue
  270     continue
  280   continue
  290 continue
!$OMP END DO NOWAIT
!$OMP END PARALLEL
c -------------------------------
c elementary diagram integrations
c -------------------------------
      do 340 il=1,lgrid*6
        edd(il,1)=edd(il,1)*coe
        ede(il,1)=ede(il,1)*coe
        eee(il,1)=eee(il,1)*coe
        eca(il,1)=eca(il,1)*coe
        ecb(il,1)=ecb(il,1)*coe
  340 continue
!$OMP PARALLEL SHARED(edd,ede,eee,eca,ecb)
!$OMP& PRIVATE(i,j,k,l)
!$OMP DO SCHEDULE(DYNAMIC)
      do 390 i=1,le
        ei=q3cne*ri(i)
        do 380 j=1,le
          eij=ei*r(j)
          ka=iabs(i-j)+1
          kb=min(i+j-1,le)
c         do 370 l=1,2,nm
          do 370 l=1,1,nm
            do 360 k=ka,kb
              eijk=eij*r(k)
              ijk=index(i,j,k)
              ikj=index(i,k,j)
              jik=index(j,i,k)
              kji=index(k,j,i)
              x=xtheta(kji)
              y=ytheta(kji)
              z=ztheta(kji)
              xijk=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
              xi=3*z**2-1
              xj=1.5*y**2-.5
              xk=1.5*x**2-.5
              edd(i,l)=edd(i,l)+eijk*
     &            ( sddd(ijk)*(xgdd(j,l)*(xgdd(k,l)*vc(l,2,1)
     &                                   +xgde(k,l)*vc(l,2,2))
     &                        +xgde(j,l)* xgdd(k,l)*vc(l,2,2))
     &             +sdde(ijk)* xgdd(j,l)* xgdd(k,l)*vc(l,2,3) )
              ede(i,l)=ede(i,l)+eijk*
     &            ( sddd(ijk)*(xgdd(j,l)*(xgde(k,l)*vc(l,2,1)
     &                                   +xgee(k,l)*vc(l,2,2))
     &                        +xgde(j,l)* xgde(k,l)*vc(l,2,2))
     &             +sdde(jik)*(xgdd(j,l)*(xgdd(k,l)*vc(l,2,1)
     &                                   +xgde(k,l)*vc(l,2,2))
     &                        +xgde(j,l)* xgdd(k,l)*vc(l,2,2))
     &             +sdde(ijk)* xgdd(j,l)* xgde(k,l)*vc(l,2,3)
     &             +sdee(ijk)* xgdd(j,l)* xgdd(k,l)*vc(l,2,3)
     &        -2*nu*sccd(kji)* xgdd(j,l)* xgcc(k)  *vc(l,2,2) )
              eee(i,l)=eee(i,l)+eijk*
     &            ( sddd(ijk)*(xgde(j,l)*(xgde(k,l)*vc(l,2,1)
     &                                   +xgee(k,l)*vc(l,2,2))
     &                        +xgee(j,l)* xgde(k,l)*vc(l,2,2))
     &           +2*sdde(jik)*(xgde(j,l)*(xgdd(k,l)*vc(l,2,1)
     &                                   +xgde(k,l)*vc(l,2,2))
     &                        +xgee(j,l)* xgdd(k,l)*vc(l,2,2))
     &             +sdde(ijk)* xgde(j,l)* xgde(k,l)*vc(l,2,3)
     &           +2*sdee(ijk)* xgde(j,l)* xgdd(k,l)*vc(l,2,3)
     &             +sdee(kji)*(xgdd(j,l)*(xgdd(k,l)*vc(l,2,1)
     &                                   +xgde(k,l)*vc(l,2,2))
     &                        +xgde(j,l)* xgdd(k,l)*vc(l,2,2))
     &             +seee(ijk)* xgdd(j,l)* xgdd(k,l)*vc(l,2,3)
     &       -2*nu*(sccd(kji)* xgde(j,l)* xgcc(k)  *vc(l,2,2)
     &             +scce(kji)* xgdd(j,l)* xgcc(k)  *vc(l,2,2)
     &             +sccd(jik)* xgcc(j)*   xgde(k,l)*vc(l,2,2)
     &             +scce(jik)* xgcc(j)*   xgdd(k,l)*vc(l,2,2)
     &             +sccd(ijk)* xgcc(j)*   xgcc(k)) )
              eca(i,l)=eca(i,l)+eijk*
     &            ( sddd(ijk)* xgca(j,l)* xgcc(k)
     &             +sccd(jik)* xgdd(j,l)* xgcc(k) )
              ecb(i,l)=ecb(i,l)+eijk*
     &            ( sddd(ijk)* xgcc(j)*   xgcb(k,l)
     &             +sccd(kji)* xgcc(j)*   xgdd(k,l) )
c             edd(i,l+2)=edd(i,l+2)+eijk*
c    &            ( sddd(ijk)*(xgdd(j,l+2)*(xgdd(k,l+2)*vc(l+2,2,1)
c    &                                     +xgde(k,l+2)*vc(l+2,2,2))
c    &                        +xgde(j,l+2)* xgdd(k,l+2)*vc(l+2,2,2))
c    &             +sdde(ijk)* xgdd(j,l+2)* xgdd(k,l+2)*vc(l+2,2,3) )
c    &                             +eijk*xi*
c    &            ( sddd(ijk)*(xgdd(j,l+4)*(xgdd(k,l+4)*vc(l+4,2,1)
c    &                                     +xgde(k,l+4)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgdd(k,l+4)*vc(l+4,2,2))
c    &             +sdde(ijk)* xgdd(j,l+4)* xgdd(k,l+4)*vc(l+4,2,3) )
c             ede(i,l+2)=ede(i,l+2)+eijk*
c    &            ( sddd(ijk)*(xgdd(j,l+2)*(xgde(k,l+2)*vc(l+2,2,1)
c    &                                     +xgee(k,l+2)*vc(l+2,2,2))
c    &                        +xgde(j,l+2)* xgde(k,l+2)*vc(l+2,2,2))
c    &             +sdde(jik)*(xgdd(j,l+2)*(xgdd(k,l+2)*vc(l+2,2,1)
c    &                                     +xgde(k,l+2)*vc(l+2,2,2))
c    &                        +xgde(j,l+2)* xgdd(k,l+2)*vc(l+2,2,2))
c    &             +sdde(ijk)* xgdd(j,l+2)* xgde(k,l+2)*vc(l+2,2,3)
c    &             +sdee(ijk)* xgdd(j,l+2)* xgdd(k,l+2)*vc(l+2,2,3)
c    &        -2*nu*sccd(kji)* xgdd(j,l+2)* xgcc(k)    *vc(l+2,2,2) )
c             ede(i,l+2)=ede(i,l+2)+eijk*xi*
c    &            ( sddd(ijk)*(xgdd(j,l+4)*(xgde(k,l+4)*vc(l+4,2,1)
c    &                                     +xgee(k,l+4)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgde(k,l+4)*vc(l+4,2,2))
c    &             +sdde(jik)*(xgdd(j,l+4)*(xgdd(k,l+4)*vc(l+4,2,1)
c    &                                     +xgde(k,l+4)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgdd(k,l+4)*vc(l+4,2,2))
c    &             +sdde(ijk)* xgdd(j,l+4)* xgde(k,l+4)*vc(l+4,2,3)
c    &             +sdee(ijk)* xgdd(j,l+4)* xgdd(k,l+4)*vc(l+4,2,3) )
c             eee(i,l+2)=eee(i,l+2)+eijk*
c    &            ( sddd(ijk)*(xgde(j,l+2)*(xgde(k,l+2)*vc(l+2,2,1)
c    &                                     +xgee(k,l+2)*vc(l+2,2,2))
c    &                        +xgee(j,l+2)* xgde(k,l+2)*vc(l+2,2,2))
c    &           +2*sdde(jik)*(xgde(j,l+2)*(xgdd(k,l+2)*vc(l+2,2,1)
c    &                                     +xgde(k,l+2)*vc(l+2,2,2))
c    &                        +xgee(j,l+2)* xgdd(k,l+2)*vc(l+2,2,2))
c    &             +sdde(ijk)* xgde(j,l+2)* xgde(k,l+2)*vc(l+2,2,3)
c    &           +2*sdee(ijk)* xgde(j,l+2)* xgdd(k,l+2)*vc(l+2,2,3)
c    &             +sdee(kji)*(xgdd(j,l+2)*(xgdd(k,l+2)*vc(l+2,2,1)
c    &                                     +xgde(k,l+2)*vc(l+2,2,2))
c    &                        +xgde(j,l+2)* xgdd(k,l+2)*vc(l+2,2,2))
c    &             +seee(ijk)* xgdd(j,l+2)* xgdd(k,l+2)*vc(l+2,2,3)
c    &       -2*nu*(sccd(kji)* xgde(j,l+2)* xgcc(k)    *vc(l+2,2,2)
c    &             +scce(kji)* xgdd(j,l+2)* xgcc(k)    *vc(l+2,2,2)
c    &             +sccd(jik)* xgcc(j)*     xgde(k,l+2)*vc(l+2,2,2)
c    &             +scce(jik)* xgcc(j)*     xgdd(k,l+2)*vc(l+2,2,2)
c    &             +sccd(ijk)* xgcc(j)*     xgcc(k)) )
c             eee(i,l+2)=eee(i,l+2)+eijk*xi*
c    &            ( sddd(ijk)*(xgde(j,l+4)*(xgde(k,l+4)*vc(l+4,2,1)
c    &                                     +xgee(k,l+4)*vc(l+4,2,2))
c    &                        +xgee(j,l+4)* xgde(k,l+4)*vc(l+4,2,2))
c    &           +2*sdde(jik)*(xgde(j,l+4)*(xgdd(k,l+4)*vc(l+4,2,1)
c    &                                     +xgde(k,l+4)*vc(l+4,2,2))
c    &                        +xgee(j,l+4)* xgdd(k,l+4)*vc(l+4,2,2))
c    &             +sdde(ijk)* xgde(j,l+4)* xgde(k,l+4)*vc(l+4,2,3)
c    &           +2*sdee(ijk)* xgde(j,l+4)* xgdd(k,l+4)*vc(l+4,2,3)
c    &             +sdee(kji)*(xgdd(j,l+4)*(xgdd(k,l+4)*vc(l+4,2,1)
c    &                                     +xgde(k,l+4)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgdd(k,l+4)*vc(l+4,2,2))
c    &             +seee(ijk)* xgdd(j,l+4)* xgdd(k,l+4)*vc(l+4,2,3) )
c             eca(i,l+2)=eca(i,l+2)+eijk*
c    &            ( sddd(ijk)* xgca(j,l+2)* xgcc(k)
c    &             +sccd(jik)* xgdd(j,l+2)* xgcc(k) )
c             ecb(i,l+2)=ecb(i,l+2)+eijk*
c    &            ( sddd(ijk)* xgcc(j)*     xgcb(k,l+2)
c    &             +sccd(kji)* xgcc(j)*     xgdd(k,l+2) )
c             edd(i,l+4)=edd(i,l+4)+eijk*xijk*
c    &            ( sddd(ijk)*(xgdd(j,l+4)*(xgdd(k,l+4)*vc(l+4,2,1)
c    &                                     +xgde(k,l+4)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgdd(k,l+4)*vc(l+4,2,2))
c    &             +sdde(ijk)* xgdd(j,l+4)* xgdd(k,l+4)*vc(l+4,2,3) )
c    &                             +eijk*xk*
c    &            ( sddd(ijk)*(xgdd(j,l+4)*(xgdd(k,l+2)*vc(l+4,2,1)
c    &                                     +xgde(k,l+2)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgdd(k,l+2)*vc(l+4,2,2))
c    &             +sdde(ijk)* xgdd(j,l+4)* xgdd(k,l+2)*vc(l+4,2,3) )
c    &                             +eijk*xj*
c    &            ( sddd(ijk)*(xgdd(j,l+2)*(xgdd(k,l+4)*vc(l+2,2,1)
c    &                                     +xgde(k,l+4)*vc(l+2,2,2))
c    &                        +xgde(j,l+2)* xgdd(k,l+4)*vc(l+2,2,2))
c    &             +sdde(ijk)* xgdd(j,l+2)* xgdd(k,l+4)*vc(l+2,2,3) )
c             ede(i,l+4)=ede(i,l+4)+eijk*xijk*
c    &            ( sddd(ijk)*(xgdd(j,l+4)*(xgde(k,l+4)*vc(l+4,2,1)
c    &                                     +xgee(k,l+4)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgde(k,l+4)*vc(l+4,2,2))
c    &             +sdde(jik)*(xgdd(j,l+4)*(xgdd(k,l+4)*vc(l+4,2,1)
c    &                                     +xgde(k,l+4)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgdd(k,l+4)*vc(l+4,2,2))
c    &             +sdde(ijk)* xgdd(j,l+4)* xgde(k,l+4)*vc(l+4,2,3)
c    &             +sdee(ijk)* xgdd(j,l+4)* xgdd(k,l+4)*vc(l+4,2,3) )
c             ede(i,l+4)=ede(i,l+4)+eijk*xk*
c    &            ( sddd(ijk)*(xgdd(j,l+4)*(xgde(k,l+2)*vc(l+4,2,1)
c    &                                     +xgee(k,l+2)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgde(k,l+2)*vc(l+4,2,2))
c    &             +sdde(jik)*(xgdd(j,l+4)*(xgdd(k,l+2)*vc(l+4,2,1)
c    &                                     +xgde(k,l+2)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgdd(k,l+2)*vc(l+4,2,2))
c    &             +sdde(ijk)* xgdd(j,l+4)* xgde(k,l+2)*vc(l+4,2,3)
c    &             +sdee(ijk)* xgdd(j,l+4)* xgdd(k,l+2)*vc(l+4,2,3)
c    &        -2*nu*sccd(kji)* xgdd(j,l+4)* xgcc(k)    *vc(l+4,2,2) )
c             ede(i,l+4)=ede(i,l+4)+eijk*xj*
c    &            ( sddd(ijk)*(xgdd(j,l+2)*(xgde(k,l+4)*vc(l+2,2,1)
c    &                                     +xgee(k,l+4)*vc(l+2,2,2))
c    &                        +xgde(j,l+2)* xgde(k,l+4)*vc(l+2,2,2))
c    &             +sdde(jik)*(xgdd(j,l+2)*(xgdd(k,l+4)*vc(l+2,2,1)
c    &                                     +xgde(k,l+4)*vc(l+2,2,2))
c    &                        +xgde(j,l+2)* xgdd(k,l+4)*vc(l+2,2,2))
c    &             +sdde(ijk)* xgdd(j,l+2)* xgde(k,l+4)*vc(l+2,2,3)
c    &             +sdee(ijk)* xgdd(j,l+2)* xgdd(k,l+4)*vc(l+2,2,3) )
c             eee(i,l+4)=eee(i,l+4)+eijk*xijk*
c    &            ( sddd(ijk)*(xgde(j,l+4)*(xgde(k,l+4)*vc(l+4,2,1)
c    &                                     +xgee(k,l+4)*vc(l+4,2,2))
c    &                        +xgee(j,l+4)* xgde(k,l+4)*vc(l+4,2,2))
c    &           +2*sdde(jik)*(xgde(j,l+4)*(xgdd(k,l+4)*vc(l+4,2,1)
c    &                                     +xgde(k,l+4)*vc(l+4,2,2))
c    &                        +xgee(j,l+4)* xgdd(k,l+4)*vc(l+4,2,2))
c    &             +sdde(ijk)* xgde(j,l+4)* xgde(k,l+4)*vc(l+4,2,3)
c    &           +2*sdee(ijk)* xgde(j,l+4)* xgdd(k,l+4)*vc(l+4,2,3)
c    &             +sdee(kji)*(xgdd(j,l+4)*(xgdd(k,l+4)*vc(l+4,2,1)
c    &                                     +xgde(k,l+4)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgdd(k,l+4)*vc(l+4,2,2))
c    &             +seee(ijk)* xgdd(j,l+4)* xgdd(k,l+4)*vc(l+4,2,3) )
c             eee(i,l+4)=eee(i,l+4)+eijk*xk*
c    &            ( sddd(ijk)*(xgde(j,l+4)*(xgde(k,l+2)*vc(l+4,2,1)
c    &                                     +xgee(k,l+2)*vc(l+4,2,2))
c    &                        +xgee(j,l+4)* xgde(k,l+2)*vc(l+4,2,2))
c    &           +2*sdde(jik)*(xgde(j,l+4)*(xgdd(k,l+2)*vc(l+4,2,1)
c    &                                     +xgde(k,l+2)*vc(l+4,2,2))
c    &                        +xgee(j,l+4)* xgdd(k,l+2)*vc(l+4,2,2))
c    &             +sdde(ijk)* xgde(j,l+4)* xgde(k,l+2)*vc(l+4,2,3)
c    &           +2*sdee(ijk)* xgde(j,l+4)* xgdd(k,l+2)*vc(l+4,2,3)
c    &             +sdee(kji)*(xgdd(j,l+4)*(xgdd(k,l+2)*vc(l+4,2,1)
c    &                                     +xgde(k,l+2)*vc(l+4,2,2))
c    &                        +xgde(j,l+4)* xgdd(k,l+2)*vc(l+4,2,2))
c    &             +seee(ijk)* xgdd(j,l+4)* xgdd(k,l+2)*vc(l+4,2,3)
c    &       -2*nu*(sccd(kji)* xgde(j,l+4)* xgcc(k)    *vc(l+2,2,2)
c    &             +scce(kji)* xgdd(j,l+4)* xgcc(k)    *vc(l+2,2,2)) )
c             eee(i,l+4)=eee(i,l+4)+eijk*xj*
c    &            ( sddd(ijk)*(xgde(j,l+2)*(xgde(k,l+4)*vc(l+2,2,1)
c    &                                     +xgee(k,l+4)*vc(l+2,2,2))
c    &                        +xgee(j,l+2)* xgde(k,l+4)*vc(l+2,2,2))
c    &           +2*sdde(jik)*(xgde(j,l+2)*(xgdd(k,l+4)*vc(l+2,2,1)
c    &                                     +xgde(k,l+4)*vc(l+2,2,2))
c    &                        +xgee(j,l+2)* xgdd(k,l+4)*vc(l+2,2,2))
c    &             +sdde(ijk)* xgde(j,l+2)* xgde(k,l+4)*vc(l+2,2,3)
c    &           +2*sdee(ijk)* xgde(j,l+2)* xgdd(k,l+4)*vc(l+2,2,3)
c    &             +sdee(kji)*(xgdd(j,l+2)*(xgdd(k,l+4)*vc(l+2,2,1)
c    &                                     +xgde(k,l+4)*vc(l+2,2,2))
c    &                        +xgde(j,l+2)* xgdd(k,l+4)*vc(l+2,2,2))
c    &             +seee(ijk)* xgdd(j,l+2)* xgdd(k,l+4)*vc(l+2,2,3)
c    &       -2*nu*(sccd(jik)* xgcc(j)*     xgde(k,l+4)*vc(l+2,2,2)
c    &             +scce(jik)* xgcc(j)*     xgdd(k,l+4)*vc(l+2,2,2)) )
c             eca(i,l+4)=eca(i,l+4)+eijk*xk*
c    &            ( sddd(ijk)* xgca(j,l+4)* xgcc(k)
c    &             +sccd(jik)* xgdd(j,l+4)* xgcc(k) )
c             ecb(i,l+4)=ecb(i,l+4)+eijk*xj*
c    &            ( sddd(ijk)* xgcc(j)*     xgcb(k,l+4)
c    &             +sccd(kji)* xgcc(j)*     xgdd(k,l+4) )
  360       continue
  370     continue
  380   continue
  390 continue
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      do 395 i=1,lg
        eca(i,1)=.5*eca(i,1)
        ecb(i,1)=.5*ecb(i,1)
        edd(i,1)=.5*edd(i,1)
        ede(i,1)=.5*ede(i,1)
        eee(i,1)=.5*eee(i,1)
        gl(i)=sl(i)-nu*(gca(i,1)+gcb(i,1)+eca(i,1)+ecb(i,1))
        gx(i)=exp(gdd(i,1)+edd(i,1))
        gy(i)=gde(i,1)+ede(i,1)
        gz(i)=gee(i,1)+eee(i,1)
  395 continue
c --------------------------------
c calculate vertex factors
c bj(l,1)=J^l(d,d,f^2) - RMP(6.29)
c bj(l,2)=J^l(e,d,f^2) - RMP(6.30)
c bj(l,3)=J^l(d,d,P)   - RMP(6.31)
c bj(l,4)=J^l(e,d,P)   - RMP(6.32)
c bj(l,5)=J^l(d,e,P)   - RMP(6.33)
c bj(l,6)=K^l(d,e,f^2) - W80(3.14)
c --------------------------------
      do 400 lj=1,8*6
        bj(lj,1)=0
  400 continue
      do 420 l=1,6,nm
        if (l.eq.1) go to 420
        vpd=vc(l,2,1)
        vpe=vc(l,2,2)
        vpp=vc(l,2,3)
        vfd=vc(l,3,1)
        vfe=vc(l,3,2)
        vfp=vc(l,3,3)
        do 410 i=1,lg
          x=q1*rs(i)*gx(i)
          fl2=f(i,l)**2
          ffl=2*f(i,1)*f(i,l)
          bj(l,1)=bj(l,1)+aa(l)*fl2*((1+gy(i))*vfd
     &     +(gy(i)+gy(i)**2+gz(i))*vfp)*x
          bj(l,2)=bj(l,2)+aa(l)*fl2*(vfd+gy(i)*vfp)*x
          bj(l,3)=bj(l,3)+.5*aa(l)*ffl*((gdd(i,l)+edd(i,l))
     &     *((1+gy(i))*vpd+(gy(i)+gy(i)**2+gz(i))*vpp)
     &     +(gde(i,l)+ede(i,l))*(1+gy(i))*vpe)*x
          bj(l,4)=bj(l,4)+.5*aa(l)*ffl*((gdd(i,l)+edd(i,l))
     &     *(vpd+gy(i)*vpp)+(gde(i,l)+ede(i,l))*vpe)*x
          bj(l,5)=bj(l,5)-aa(l)*af(l)*(ffl+f(i,1)**2
     &     *(gdd(i,l)+edd(i,l)))*gl(i)**2*vpe*x/nu
          bj(l,6)=bj(l,6)-afe(l)*fl2*gl(i)**2*vfe*x/nu
  410   continue
  420 continue
c --------------------------------
c --------------------------------
      do 470 k=1+nm,6,nm
        do 430 l=1,3
          vce(k,l)=0
          vco(k,1,l)=vc(k,1,l)
          vco(k,2,l)=vc(k,2,l)
          vco(k,3,l)=vc(k,3,l)
          vc(k,1,l)=0
          vc(k,2,l)=0
          vc(k,3,l)=0
  430   continue
        do 450 l=1+nm,6,nm
          dlk=ad(l,k)
          vc(k,1,1)=vc(k,1,1)+dlk*(bj(l,1)/2+bj(l,3)/3)
          vce(k,1)=vce(k,1)+dlk*(bj(l,5)/4+bj(l,6)/2)
          vc(k,1,2)=vc(k,1,2)+dlk*(bj(l,2)+bj(l,4)/2)
          vc(k,1,3)=vc(k,1,3)+dlk*(bj(l,2)/2+bj(l,4)/3)
          vc(k,2,1)=vc(k,2,1)+dlk*(5*bj(l,1)/12+bj(l,3)/3)
          vce(k,2)=vce(k,2)+dlk*(bj(l,5)/3+5*bj(l,6)/12)
          vc(k,2,3)=vc(k,2,3)+dlk*(5*bj(l,2)/12+bj(l,4)/3)
          vc(k,3,1)=vc(k,3,1)+dlk*(bj(l,1)/2+5*bj(l,3)/12)
          vce(k,3)=vce(k,3)+dlk*(bj(l,5)/2+bj(l,6)/2)
          vc(k,3,2)=vc(k,3,2)+dlk*(bj(l,2)/2+5*bj(l,4)/12)
          vc(k,3,3)=vc(k,3,3)+dlk*(bj(l,2)/2+5*bj(l,4)/12)
          do 440 n=1,4,nm
            dkn=ad(k,n)
            dln=ad(l,n)
            vce(k,2)=vce(k,2)+dkn*(ak(l,l,n)*aa(n)/afe(l))*bj(l,6)/8
            vce(k,3)=vce(k,3)+dkn*(ak(l,l,n)*aa(n)/afe(l))*bj(l,6)/4
            vc(k,3,2)=vc(k,3,2)+dln*(ak(k,k,n)*aa(n)/afe(k))*(bj(l,2)/4
     &       +bj(l,4)/8)
  440     continue
  450   continue
        vc(k,2,2)=vc(k,1,3)
        do 460 l=1,3
        do 460 m=1,3
          vc(k,l,m)=((1+at(m,1)*vce(k,l))*exp(vc(k,l,m)))*cn
     &             +vco(k,l,m)*co
  460   continue
  470 continue
  490 continue
c <=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=
c -------------------------------------
c add elementary diagrams to SOC chains
c for remainder of program
c -------------------------------------
      do 500 l=1+nm,6,nm
      do 500 i=1,lg
        gca(i,l)=gca(i,l)+eca(i,l)
        gcb(i,l)=gcb(i,l)+ecb(i,l)
        gdd(i,l)=gdd(i,l)+edd(i,l)
        gde(i,l)=gde(i,l)+ede(i,l)
        gee(i,l)=gee(i,l)+eee(i,l)
  500 continue
c -------------
c moc integrals
c -------------
      do 505 k=1,lgrid
        xcc(k)=0
  505 continue
      do 510 kl=1,lgrid*6
        gfcc(kl,1)=0
        gfdd(kl,1)=0
        gfde(kl,1)=0
        gfed(kl,1)=0
        ghcc(kl,1)=0
        ghdd(kl,1)=0
        ghde(kl,1)=0
        ghed(kl,1)=0
        xca(kl,1)=0
        xdd(kl,1)=0
        xee(kl,1)=0
        xgcb(kl,1)=0
        xgdd(kl,1)=0
        xgde(kl,1)=0
        xgee(kl,1)=0
  510 continue
      do 530 l=1+nm,6,nm
        do 512 k=1,lg
          xee(k,l)=xee(k,l)+aa(l)*2*f(k,1)*f(k,l)*(-gl(k)**2/nu)*gx(k)
          xgee(k,l)=xgee(k,l)+aa(l)*2*f(k,1)*f(k,l)*(-gl(k)**2/nu)*gx(k)
  512   continue
        do 526 i=1+nm,6,nm
          do 524 j=1+nm,6,nm
            x1=al(i,j,l)
            x2=ak(i,j,l)*aa(l)
            x3=x2*af(l)
            x4=0
            x5=0
            x6=0
            do 516 n=1,4,nm
              do 514 m=1,6,nm
                x4=x4+ak(i,n,m)*al(m,j,l)+ak(j,n,m)*al(m,i,l)
                x5=x5+.5*(ak(i,n,m)*al(m,j,l)+ak(j,n,m)*al(m,i,l)
     &                   +ak(i,j,m)*(al(m,n,l)+ak(m,n,l)*aa(l)))
  514         continue
              x6=x6+aa(l)*ak(i,j,n)*aa(n)
  516       continue
            do 518 k=1,lg
              xca(k,l)=xca(k,l)-x4*f(k,i)*f(k,j)*gl(k)*gx(k)/nu
              xcc(k)=xcc(k)-x3*f(k,i)*f(k,j)*gl(k)*gx(k)/nu
     &              *vc(i,3,2)*vc(j,3,2)
              xdd(k,l)=xdd(k,l)+x1*f(k,i)*f(k,j)*gx(k)
              xee(k,l)=xee(k,l)+x1*2*f(k,1)*f(k,i)*(-af(j)*gl(k)**2/nu)
     &                *gx(k)
              xgcb(k,l)=xgcb(k,l)-x6*f(k,i)*f(k,j)*gl(k)*gx(k)/nu
              xgdd(k,l)=xgdd(k,l)+x2*f(k,i)*f(k,j)*gx(k)
              xgde(k,l)=xgde(k,l)-x5*f(k,i)*f(k,j)*gl(k)*gx(k)/nu
              xgee(k,l)=xgee(k,l)
     &                 +x2*2*f(k,1)*f(k,i)*(-af(j)*gl(k)**2/nu)*gx(k)
  518       continue
            afi(l,i,j)=0
            afj(l,i,j)=0
            afk(l,i,j)=0
            ahi(l,i,j)=0
            ahj(l,i,j)=0
            ahk(l,i,j)=0
            do 522 n=1,4,nm
            do 522 np=1,4,nm
            do 522 k=1,6,nm
              do 520 m=1,6,nm
                afi(l,i,j)=afi(l,i,j)+(ak(n,l,k)*aa(k)+al(n,l,k))
     &           *(ak(np,i,m)*ax(m,j,k)+ak(np,j,m)*ax(i,m,k))/12
                afj(l,i,j)=afj(l,i,j)+(ak(n,i,k)*aa(k)+al(n,i,k))
     &           *(ak(np,j,m)*ax(m,l,k)+ak(np,l,m)*ax(j,m,k))/12
                afk(l,i,j)=afk(l,i,j)+(ak(n,j,k)*aa(k)+al(n,j,k))
     &           *(ak(np,l,m)*ax(m,i,k)+ak(np,i,m)*ax(l,m,k))/12
                ahi(l,i,j)=ahi(l,i,j)+(2*al(n,l,k)+ak(n,l,k)*aa(k))
     &           *(ak(np,i,m)*ax(m,j,k)+ak(np,j,m)*ax(i,m,k))/16
                ahj(l,i,j)=ahj(l,i,j)+((2*ak(n,i,k)*aa(k)+al(n,i,k))
     &           *ak(np,j,m)*ax(m,l,k)+(ak(n,i,k)*aa(k)+al(n,i,k))
     &           *ak(np,l,m)*ax(j,m,k))/16
                ahk(l,i,j)=ahk(l,i,j)+((2*ak(n,j,k)*aa(k)+al(n,j,k))
     &           *ak(np,i,m)*ax(l,m,k)+(ak(n,j,k)*aa(k)+al(n,j,k))
     &           *ak(np,l,m)*ax(m,i,k))/16
  520         continue
  522       continue
  524     continue
  526   continue
  530 continue
      do 550 i=1,lg
        ei=q2*ri(i)
        do 548 j=1,lg
          eij=ei*r(j)
          do 534 l=1+nm,6,nm
            pcc(l)=-af(l)*gl(j)/nu
            qcc(l)=f(j,1)**2*gx(j)*pcc(l)
            qdd(l)=2*f(j,1)*f(j,l)*gx(j)
            qde(l)=-qdd(l)*gl(j)/nu
            qee(l)=qcc(l)*gl(j)
  534     continue
          ka=iabs(i-j)+1
          kb=min(i+j-1,lg)
          do 540 il=1,2,nm
            do 538 it=3-il,6
              l=il+lit(it)
              m=il+mit(it)
              n=il+nit(it)
              do 536 k=ka,kb
                eijk=eij*r(k)
                kji=index(k,j,i)
                x=xtheta(kji)
                y=ytheta(kji)
                z=ztheta(kji)
                xi=1.5*z**2-.5
                xj=1.5*y**2-.5
                xk=1.5*x**2-.5
                xijk=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
                axi=at(1,it)+at(2,it)+2*at(3,it)*xi+at(4,it)*xj
     &             +at(5,it)*xk+at(6,it)*xijk
                axj=at(1,it)+at(2,it)+at(3,it)*xi+2*at(4,it)*xj
     &             +at(5,it)*xk+at(6,it)*xijk
                axk=at(1,it)+at(2,it)+at(3,it)*xi+at(4,it)*xj
     &           +2*at(5,it)*xk+at(6,it)*xijk
                gfcc(i,l)=gfcc(i,l)+eijk*axk
     &           *(qcc(m)*xgde(k,n)-pcc(m)*xgcb(k,n)*af(l))*vc(l,3,3)
                gfdd(i,l)=gfdd(i,l)+eijk*axk
     &           *qdd(m)*(5*xdd(k,n)+7*xgdd(k,n))*vc(l,3,1)/12
                gfed(i,l)=gfed(i,l)+eijk*axk
     &           *qdd(m)*(xee(k,n)+2*xgee(k,n))*vc(l,2,3)/3
                gfde(i,l)=gfde(i,l)+eijk*axk
     &           *qee(m)*.5*(xdd(k,n)+xgdd(k,n))*vc(l,3,3)
                ghcc(i,l)=ghcc(i,l)+eijk*axk
     &           *(qcc(m)*xca(k,n)-pcc(m)*xgcb(k,n)*af(l))*vc(l,3,3)
                ghdd(i,l)=ghdd(i,l)+eijk*axk
     &           *qdd(m)*.5*(xdd(k,n)+xgdd(k,n))*vc(l,3,1)
                ghed(i,l)=ghed(i,l)+eijk*axk
     &           *qdd(m)*.25*(xee(k,n)+3*xgee(k,n))*vc(l,2,3)
                ghde(i,l)=ghde(i,l)+eijk*axk
     &           *qee(m)*xdd(k,n)*vc(l,3,3)
                rde=-2*f(k,1)*f(k,n)*gl(k)*gx(k)/nu
                gfcc(i,l)=gfcc(i,l)+eijk*(axi*afi(l,m,n)+axj*afj(l,m,n)
     &           +axk*afk(l,m,n))*qde(m)*rde*vc(l,2,3)
                ghcc(i,l)=ghcc(i,l)+eijk*(axi*ahi(l,m,n)+axj*ahj(l,m,n)
     &           +axk*ahk(l,m,n))*qde(m)*rde*vc(l,2,3)
  536         continue
  538       continue
  540     continue
          do 546 k=ka,kb
            eijk=eij*r(k)
            ghcc(i,1)=ghcc(i,1)+eijk*(2*qcc(3)-pcc(3))*xcc(k)
            do 544 n=1+nm,4,nm
              rde=-2*f(k,1)*f(k,n)*gl(k)*gx(k)/nu
              do 542 m=1+nm,4,nm
                ghcc(i,1)=ghcc(i,1)+eijk*aa(m)*aa(n)*(1+.5*ad(m,n))
     &           *qde(m)*rde*vc(m,2,2)*vc(n,2,2)
  542         continue
  544       continue
  546     continue
  548   continue
  550 continue
c -------------------
c sor-chain integrals
c -------------------
      do 555 kl=1,lgrid*6
        grdc(kl,1)=0
        grdd(kl,1)=0
        grde(kl,1)=0
        gred(kl,1)=0
        gree(kl,1)=0
        grfc(kl,1)=0
        grfd(kl,1)=0
        grfe(kl,1)=0
        grmd(kl,1)=0
        grme(kl,1)=0
  555 continue
      do 590 i=1,lg
        ei=q2*ri(i)
        do 580 j=1,lg
          eij=ei*r(j)
          do 560 l=1+nm,6,nm
            pdd(l)=aa(l)*(2*f(j,1)*f(j,l)+.5*f(j,1)**2*gdd(j,l))
     &       *gdd(j,l)*gx(j)
            pde(l)=aa(l)*(2*f(j,1)*f(j,l)+f(j,1)**2*gdd(j,l))*gde(j,l)
     &       *gx(j)
            pee(l)=aa(l)*((2*f(j,1)*f(j,l)+f(j,1)**2*gdd(j,l))*(gee(j,l)
     &       -af(l)*gl(j)**2/nu)+f(j,1)**2*(gde(j,l)**2+af(l)*(gca(j,l)
     &       +gcb(j,l))*gl(j)/vc(l,2,2)))*gx(j)
            qde(l)=pdd(l)*gy(j)
            qee(l)=pdd(l)*gz(j)+qde(l)*gy(j)
            qcc(l)=pde(l)*gy(j)
            zca(l)=-f(j,l)**2*gx(j)*gl(j)**2/nu
            zdd(l)=aa(l)*f(j,l)**2*gx(j)
            zde(l)=zdd(l)*gy(j)
            zee(l)=zdd(l)*gz(j)+zde(l)*gy(j)
  560     continue
          ka=iabs(i-j)+1
          kb=min(i+j-1,lg)
          do 570 l=1+nm,6,nm
            do 565 k=ka,kb
              eijk=eij*r(k)
              rdd=f(k,1)**2*gx(k)
              rde=rdd*gy(k)
              ree=rdd*(gz(k)+gy(k)**2)
              rdd=rdd-1
              grdd(i,l)=grdd(i,l)+eijk*(pdd(l)*(rdd*vc(l,2,1)+rde
     &         *vc(l,2,3))+(pde(l)*vc(l,2,2)+qde(l)*vc(l,2,3))*rdd)
              grde(i,l)=grde(i,l)+eijk*(pdd(l)*(rde*vc(l,2,1)+ree
     &         *vc(l,2,3))+(pde(l)*vc(l,2,2)+2*qde(l)*vc(l,2,3))*rde
     &         +qde(l)*rdd*vc(l,2,1)+(qee(l)*vc(l,2,3)+qcc(l)*vc(l,2,2))
     &         *rdd)
              grdc(i,l)=grdc(i,l)+eijk*(qde(l)*(rde*vc(l,2,1)
     &         +ree*vc(l,2,3))+(qcc(l)*vc(l,2,2)
     &         +qee(l)*vc(l,2,3))*rde)
              gred(i,l)=gred(i,l)+eijk*(pde(l)*(rdd*vc(l,2,1)
     &         +rde*vc(l,2,3))+(pee(l)*vc(l,2,2)+qcc(l)*vc(l,2,3))*rdd)
              gree(i,l)=gree(i,l)+eijk*(pde(l)*(rde*vc(l,2,1)
     &         +ree*vc(l,2,3))+(pee(l)*vc(l,2,2)+qcc(l)*vc(l,2,3))*rde)
              grfd(i,l)=grfd(i,l)+eijk*(zdd(l)*(rdd*vc(l,3,1)
     &         +rde*vc(l,3,3))+zde(l)*rdd*vc(l,3,3))
              grfe(i,l)=grfe(i,l)+eijk*(zdd(l)*(rde*vc(l,3,1)
     &         +ree*vc(l,3,3))+zde(l)*(rdd*vc(l,3,1)+2*rde*vc(l,3,3))
     &         +zee(l)*rdd*vc(l,3,3))
              grfc(i,l)=grfc(i,l)+eijk*(zde(l)*(rde*vc(l,3,1)
     &         +ree*vc(l,3,3))+zee(l)*rde*vc(l,3,3))
              grmd(i,l)=grmd(i,l)+eijk*zca(l)*rdd*vc(l,3,2)
              grme(i,l)=grme(i,l)+eijk*zca(l)*rde*vc(l,3,2)
  565       continue
  570     continue
  580   continue
  590 continue
c -------------
c l.s integrals
c -------------
      if (nv.le.6) go to 690
      do 605 l=7,8,nm
        do 600 i=1,lg
          bj(l,1)=bj(l,1)+f(i,l)**2*rs(i)**2*gx(i)
  600   continue
        bj(l,1)=ksav*q1*aa(l)*bj(l,1)/6
  605 continue
      do 615 l=1+nm,4,nm
        bk(l,2)=0
        bk(l,3)=0
        do 610 i=1,lg
          x=q1*rs(i)*gx(i)/6
          y=aa(l)*(2*f(i,1)*f(i,l)+f(i,1)**2*gdd(i,l))
          bk(l,2)=bk(l,2)-2*x*y*slps(i)/nu
          bk(l,3)=bk(l,3)+2*x*y*sl(i)*(sldp(i)+2*slp(i)*ri(i))/nu
  610   continue
        bk(l,1)=bk(l,2)-bk(l,3)
  615 continue
      do 630 i=1,lg
        ei=q2*ri(i)
        do 625 j=1,lg
          eij=ei*r(j)
          ka=iabs(i-j)+1
          kb=min(i+j-1,lg)
          y1=f(j,1)**2*gx(j)-1
          y2=y1*(sldp(j)+slp(j)*ri(j))
          y3=y1*(sldp(j)-slp(j)*ri(j))
          y4=y1*slp(j)
          y5=y1*sl(j)
          y6=y1*(slps(j)+sl(j)*sldp(j)+2*sl(j)*slp(j)*ri(j))
c ----------------------------------------------------------
c Arya Akmal version:
c         z1=f(j,1)**2*gx(j)
c         z2=gl(j)
c         z3=0.25*(sldp(j)-slp(j)*ri(j))
c         z4=0.25*(sldp(j)+slp(j)*ri(j))
c         z5=slp(j)*ri(i)
c ----------------------------------------------------------
          do 620 k=ka,kb
            eijk=eij*r(k)
            kji=index(k,j,i)
            x=xtheta(kji)
            y=ytheta(kji)
            z=ztheta(kji)
            z1=f(k,1)**2*gx(k)
            z2=sldp(k)+slp(k)*ri(k)-(sldp(k)-slp(k)*ri(k))*y**2
            bde(i,1)=bde(i,1)+eijk*.5*y1*z1*z2*sl(k)
            bde(i,2)=bde(i,2)-eijk*(.5*y1*z1*slps(k)*(1-y**2)+y6/3)
            bde(i,3)=bde(i,3)-eijk*y1*z1*sl(k)*slp(k)*y
            bcc(i,1)=bcc(i,1)+eijk*.5*(y2-y3*x**2)*z1*sl(k)
            bcc(i,2)=bcc(i,2)+eijk*.5*(y4*(z+x*y)*(1+z1)*slp(k)-y5*z2)
            bcc(i,3)=bcc(i,3)-eijk*.25*y4*z1*sl(k)*x
c ----------------------------------------------------------
c Arya Akmal version:
c           y1=f(k,1)**2*gx(k)-1
c           y2=gl(k)
c           bde(i,1)=bde(i,1)-eijk*2.*y1*z1*z2*(x**2*z3-z4+x*z5)
c           bcc(i,1)=bcc(i,1)-eijk*2.*y1*z1*y2*(x**2*z3-z4+x*z5)
c ----------------------------------------------------------
  620     continue
  625   continue
        bde(i,1)=rs(i)*bde(i,1)/nu
        bde(i,2)=rs(i)*bde(i,2)/nu
        bde(i,3)=r(i)*bde(i,3)/nu
        bcc(i,1)=sl(i)*rs(i)*bcc(i,1)/nu**2
        bcc(i,2)=sl(i)*rs(i)*bcc(i,2)/nu**2
        bcc(i,3)=sl(i)*r(i)*bcc(i,3)/nu**2
c ----------------------------------------------------------
c Arya Akmal version [defines only bde(1) & bcc(1)]:
c       bde(i,1)=rs(i)*bde(i,1)/nu
c       bcc(i,1)=gl(i)*rs(i)*bcc(i,1)/nu**2
c ----------------------------------------------------------
  630 continue
      if (nv.le.8) go to 690
c --------------
c l**2 integrals
c --------------
      do 670 i=1,6,nm
        bq(i,1)=0
        bq(i,2)=0
        if (i.eq.1) then
          do 655 j=1,lg
            x=rs(j)*gx(j)
            bq(1,1)=bq(1,1)+f(j,1)*fds(j,1)*x
            bq(1,2)=bq(1,2)+(f(j,1)*fds(j,1)*sls(j)
     &       +2*f(j,1)*fp(j,1)*sl(j)*slp(j)+f(j,1)**2*(sl(j)*sldp(j)
     &       +slps(j)+2*sl(j)*slp(j)*ri(j)))*x
  655     continue
        else if (i.ge.5) then
          do 660 j=1,lg
            x=rs(j)*gx(j)
            bq(i,1)=bq(i,1)+f(j,i)*fds(j,i)*vc(i,3,1)*x
  660     continue
        else
          do 665 j=1,lg
            x=rs(j)*gx(j)
            bq(i,1)=bq(i,1)+f(j,i)*fds(j,i)*vc(i,3,1)*x
            bq(i,2)=bq(i,2)+((f(j,i)*fds(j,1)+f(j,1)*fds(j,i))*sls(j)
     &       +2*(f(j,i)*fp(j,1)+f(j,1)*fp(j,i))*sl(j)*slp(j)+2*f(j,1)
     &       *f(j,i)*(sl(j)*sldp(j)+slps(j)+2*sl(j)*slp(j)*ri(j)))
     &       *vc(i,2,2)*x
  665     continue
        end if
        bq(i,1)=-(2./3.)*q1*aa(i)*bq(i,1)
        bq(i,2)=(2./3.)*q1*aa(i)*bq(i,2)/nu
  670 continue
c print ================================================================
  690 if (no.eq.0) go to 700
      write(nlog,1000) ni,cn,co
      write(nout,1000) ni,cn,co
 1000 format(/4x,'# chain iterations (new/old) = ',i3,' (',f3.2
     &      ,'/',f3.2,')')
      if (nie.ge.1) then
        write(nlog,1002) nie,cne,coe,nphi
        write(nout,1002) nie,cne,coe,nphi
 1002   format(/4x,'# elementary diagram iterations (new/old) = ',i2
     &        ,' (',f3.2,'/',f3.2,')',' nphi = ',i1)
      end if
      write(nout,1006)
 1006 format(/4x,'j(i:ddf2,edf2,ddp,edp,dep,def2)')
      write(nout,1008) bj
 1008 format(8f8.3)
      if (nv.gt.8) then
        write(nout,1010) bq
 1010   format(/4x,'q(i:ddf2,dep)',/,(6f8.3))
      end if
      if (no.ge.1) write(nout,1020) vc
 1020 format(/4x,'m(i,xt(i,p,f2),xy(d,e,ep))',9(/,6f8.3))
c ======================================================================
c ------------------------------
c set three-body force functions
c ------------------------------
  700 if (nt.eq.0) go to 900
      if (nt.le.1.or.nt.ge.4) then
        xmu=.7
        do 710 j=1,lmax
          rxmu=xmu*r(j)
          rcut=1-exp(-cut*rs(j))
          ypi(j)=exp(-rxmu)*rcut/rxmu
          tpi(j)=(1+3/rxmu+3/rxmu**2)*ypi(j)*rcut
          rcut0=1-exp(-cut0*rs(j))
          tpi2(j)=(tpi(j)*(rcut0/rcut)**2)**2
          xt1(j)=(ypi(j)-tpi(j))*r(j)/3
  710   continue
      else if (nt.eq.2.or.nt.eq.3) then
        pa=0.1102*(3-nt)+0.1023*(nt-2)
        pb=-0.2517*(3-nt)-0.2230*(nt-2)
        pc=0.09754*(3-nt)
        xmu=139.6/197.3
        xlm=plm
        xlm2=xlm**2
        do 720 j=1,lmax
          rmu=xmu*r(j)
          rlm=xlm*rmu
          ermu=exp(-rmu)/rmu
          erlm=exp(-rlm)/rmu
          xt0(j)=.5*(xlm2-1)**2*erlm*rmu
          xt1(j)=(1+1/rmu)*ermu-xlm*(1+1/rlm)*erlm-.5*(xlm2-1)*erlm*rmu
          xt2(j)=(1+3/rmu+3/rmu**2)*ermu-xlm2*(1+3/rlm+3/rlm**2)*erlm
     &          -.5*xlm*(xlm2-1)*(1+1/rlm)*erlm*rmu
          xt3(j)=ermu-xlm2*erlm-.5*xlm*(xlm2-1)*(1-2/rlm)*erlm*rmu
  720   continue
      end if
c -------------------------------
c calculate effective two-body v3
c -------------------------------
      if (nt.ge.4) go to 900
      do 730 i=1,l3
        xgdd(i,1)=f(i,1)**2*gx(i)
        xgde(i,1)=xgdd(i,1)*gy(i)
        xgee(i,1)=xgdd(i,1)*(gy(i)**2+gz(i)-gl(i)**2/nu)
        xgcc(i)=-xgdd(i,1)*gl(i)/nu
  730 continue
      do 790 i=1,l3
        ei=q2/r(i)
        do 780 j=1,l3
          eij=ei*r(j)
          ka=iabs(i-j)+1
          kb=min(i+j-1,l3)
c ---------
c urbana v3
c ---------
          if (nt.eq.1) then
cdir$ ivdep
          do 750 k=ka,kb
            eijk=eij*r(k)
            kji=index(k,j,i)
            x=xtheta(kji)
            y=ytheta(kji)
            z=ztheta(kji)
            xi=3*z**2-1
            xj=1.5*y**2-.5
            xk=1.5*x**2-.5
            xijk=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
            v3st(n3s)=4*tnia*(ypi(j)*ypi(k)+tpi(j)*tpi(k)*xi)
            v3st(n3t)=4*tnia*(ypi(j)*tpi(k)*xj+tpi(j)*ypi(k)*xk
     &                       +tpi(j)*tpi(k)*xijk)
            do 740 l=n3s,n3t,2
              v3dd(i,l,1)=v3dd(i,l,1)+v3st(l)*udd(j,k,1,1)
              v3de(i,l,1)=v3de(i,l,1)+v3st(l)*ude(j,k,1,1)
              v3ee(i,l,1)=v3ee(i,l,1)+v3st(l)*uee(j,k,1,1)
              v3cc(i,l,1)=v3cc(i,l,1)+v3st(l)*eijk*xgcc(j)*xgcc(k)
  740       continue
c ------------
c s-wave parts
c ------------
            if (tnix.eq.0.) go to 750
            pap=tnix*xt1(j)*xt1(k)/z
            v3st(n3s)=pap*(xi+1)
            v3st(n3t)=pap*(xijk+xj+xk)
            do 745 l=n3s,n3t,2
              v3dd(i,l,2)=v3dd(i,l,2)+v3st(l)*udd(j,k,1,1)
              v3de(i,l,2)=v3de(i,l,2)+v3st(l)*ude(j,k,1,1)
              v3ee(i,l,2)=v3ee(i,l,2)+v3st(l)*uee(j,k,1,1)
              v3cc(i,l,2)=v3cc(i,l,2)+v3st(l)*eijk*xgcc(j)*xgcc(k)
  745       continue
  750     continue
c -------------------
c tucson-melbourne v3
c -------------------
          else if (nt.eq.2.or.nt.eq.3) then
cdir$ ivdep
          do 770 k=ka,kb
            eijk=eij*r(k)
            kji=index(k,j,i)
            x=xtheta(kji)
            y=ytheta(kji)
            z=ztheta(kji)
            xi=3*z**2-1
            xj=1.5*y**2-.5
            xk=1.5*x**2-.5
            xijk=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
            v3st(n3s)=pb*(xt3(j)*xt3(k)+xt2(j)*xt2(k)*xi)
            v3st(n3t)=pb*(xt3(j)*xt2(k)*xj+xt2(j)*xt3(k)*xk
     &                   +xt2(j)*xt2(k)*xijk)
            do 760 l=n3s,n3t,2
              v3dd(i,l,1)=v3dd(i,l,1)+v3st(l)*udd(j,k,1,1)
              v3de(i,l,1)=v3de(i,l,1)+v3st(l)*ude(j,k,1,1)
              v3ee(i,l,1)=v3ee(i,l,1)+v3st(l)*uee(j,k,1,1)
              v3cc(i,l,1)=v3cc(i,l,1)+v3st(l)*eijk*xgcc(j)*xgcc(k)
  760       continue
c ------------
c s-wave parts
c ------------
            if (tnix.eq.0.) go to 770
            pac=((pa-2*pc)*xt1(j)*xt1(k)
     &         +pc*(xt0(j)*xt1(k)+xt1(j)*xt0(k)))/z
            v3st(n3s)=pac*(xi+1)
            v3st(n3t)=pac*(xijk+xj+xk)
            do 765 l=n3s,n3t,2
              v3dd(i,l,2)=v3dd(i,l,2)+v3st(l)*udd(j,k,1,1)
              v3de(i,l,2)=v3de(i,l,2)+v3st(l)*ude(j,k,1,1)
              v3ee(i,l,2)=v3ee(i,l,2)+v3st(l)*uee(j,k,1,1)
              v3cc(i,l,2)=v3cc(i,l,2)+v3st(l)*eijk*xgcc(j)*xgcc(k)
  765       continue
  770     continue
          end if
  780   continue
  790 continue
c --------------------------
c fhnc/4 contributions to v3
c --------------------------
      if (nie.eq.0) go to 900
      do 890 i=1,le
        ei=q2/r(i)
        do 880 j=1,le
          eij=ei*r(j)
          ka=iabs(i-j)+1
          kb=min(i+j-1,le)
c ---------
c urbana v3
c ---------
          if (nt.eq.1) then
cdir$ ivdep
          do 850 k=ka,kb                         ! r23
            eijk=eij*r(k)
            ijk=index(i,j,k)
            jik=index(j,i,k)
            kji=index(k,j,i)
            x=xtheta(kji)
            y=ytheta(kji)
            z=ztheta(kji)
            xi=3*z**2-1
            xj=1.5*y**2-.5
            xk=1.5*x**2-.5
            xijk=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
            v3st(n3s)=4*tnia*(ypi(j)*ypi(k)+tpi(j)*tpi(k)*xi)
            v3st(n3t)=4*tnia*(ypi(j)*tpi(k)*xj+tpi(j)*ypi(k)*xk
     &                       +tpi(j)*tpi(k)*xijk)
            do 840 l=n3s,n3t,2
              v3dd(i,l,1)=v3dd(i,l,1)+v3st(l)*eijk*
     &              ( sddd(ijk)*xgdd(j,1)*(xgdd(k,1)+2*xgde(k,1))
     &               +sdde(ijk)*xgdd(j,1)*xgdd(k,1) )
              v3de(i,l,1)=v3de(i,l,1)+v3st(l)*eijk*
     &              ( sddd(ijk)*(xgdd(j,1)*(xgde(k,1)+xgee(k,1))
     &                          +xgde(j,1)*xgde(k,1))
     &               +sdde(jik)*(xgdd(j,1)*(xgdd(k,1)+xgde(k,1))
     &                          +xgde(j,1)*xgdd(k,1))
     &               +sdde(ijk)*xgdd(j,1)*xgde(k,1)
     &               +sdee(ijk)*xgdd(j,1)*xgdd(k,1)
     &          -2*nu*sccd(kji)*xgdd(j,1)*xgcc(k) )
              v3ee(i,l,1)=v3ee(i,l,1)+v3st(l)*eijk*
     &              ( sddd(ijk)*xgde(j,1)*(xgde(k,1)+2*xgee(k,1))
     &             +2*sdde(jik)*(xgde(j,1)*(xgdd(k,1)+xgde(k,1))
     &                          +xgee(j,1)*xgdd(k,1))
     &               +sdde(ijk)*xgde(j,1)*xgde(k,1)
     &             +2*sdee(ijk)*xgde(j,1)*xgdd(k,1)
     &               +sdee(kji)*xgdd(j,1)*(xgdd(k,1)+2*xgde(k,1))
     &               +seee(ijk)*xgdd(j,1)*xgdd(k,1)
     &       -2*nu*(2*sccd(kji)*xgde(j,1)*xgcc(k)
     &             +2*scce(kji)*xgdd(j,1)*xgcc(k)
     &               +sccd(ijk)*xgcc(j)*xgcc(k)) )
              v3cc(i,l,1)=v3cc(i,l,1)+v3st(l)*eijk*
     &            ( 2*sddd(ijk)*xgcc(j)*xgcc(k)
     &             +2*sccd(kji)*xgcc(j)*xgdd(k,1) )
  840       continue
  850     continue
c -------------------
c tucson-melbourne v3
c -------------------
          else if (nt.eq.2.or.nt.eq.3) then
cdir$ ivdep
          do 870 k=ka,kb
            eijk=eij*r(k)
            ijk=index(i,j,k)
            jik=index(j,i,k)
            kji=index(k,j,i)
            x=xtheta(kji)
            y=ytheta(kji)
            z=ztheta(kji)
            xi=3*z**2-1
            xj=1.5*y**2-.5
            xk=1.5*x**2-.5
            xijk=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
            pac=((pa-2*pc)*xt1(j)*xt1(k)
     &         +pc*(xt0(j)*xt1(k)+xt1(j)*xt0(k)))/z
            v3st(n3s)=pb*(xt3(j)*xt3(k)+xt2(j)*xt2(k)*xi)+pac*(xi+1)
            v3st(n3t)=pb*(xt3(j)*xt2(k)*xj+xt2(j)*xt3(k)*xk
     &                   +xt2(j)*xt2(k)*xijk)+pac*(xijk+xj+xk)
            do 860 l=n3s,n3t,2
              v3dd(i,l,1)=v3dd(i,l,1)+v3st(l)*eijk*
     &              ( sddd(ijk)*xgdd(j,1)*(xgdd(k,1)+2*xgde(k,1))
     &               +sdde(ijk)*xgdd(j,1)*xgdd(k,1) )
              v3de(i,l,1)=v3de(i,l,1)+v3st(l)*eijk*
     &              ( sddd(ijk)*(xgdd(j,1)*(xgde(k,1)+xgee(k,1))
     &                            +xgde(j,1)*xgde(k,1))
     &               +sdde(jik)*(xgdd(j,1)*(xgdd(k,1)+xgde(k,1))
     &                            +xgde(j,1)*xgdd(k,1))
     &               +sdde(ijk)*xgdd(j,1)*xgde(k,1)
     &               +sdee(ijk)*xgdd(j,1)*xgdd(k,1)
     &          -2*nu*sccd(kji)*xgdd(j,1)*xgcc(k) )
              v3ee(i,l,1)=v3ee(i,l,1)+v3st(l)*eijk*
     &              ( sddd(ijk)*xgde(j,1)*(xgde(k,1)+2*xgee(k,1))
     &             +2*sdde(jik)*(xgde(j,1)*(xgdd(k,1)+xgde(k,1))
     &                            +xgee(j,1)*xgdd(k,1))
     &               +sdde(ijk)*xgde(j,1)*xgde(k,1)
     &             +2*sdee(ijk)*xgde(j,1)*xgdd(k,1)
     &               +sdee(kji)*xgdd(j,1)*(xgdd(k,1)+2*xgde(k,1))
     &               +seee(ijk)*xgdd(j,1)*xgdd(k,1)
     &       -2*nu*(2*sccd(kji)*xgde(j,1)*xgcc(k)
     &             +2*scce(kji)*xgdd(j,1)*xgcc(k)
     &               +sccd(ijk)*xgcc(j)*xgcc(k)) )
              v3cc(i,l,1)=v3cc(i,l,1)+v3st(l)*eijk*
     &              ( 2*sddd(ijk)*xgcc(j)*xgcc(k)
     &               +2*sccd(kji)*xgcc(j)*xgdd(k,1) )
  860       continue
  870     continue
          end if
  880   continue
  890 continue
  900 continue
      return
      end subroutine nmhnc

      function sdd(j,k,l,m)
        use nmvar
        implicit none
        real*8 :: sdd 
        integer*4 :: j,k,l,m
        sdd = eijk*(xdd(j,l)*xgdd(k,m)*vc(l,2,1)
     &            +(xdd(j,l)*xgde(k,m)+xde(j,l)*xgdd(k,m))*vc(l,2,2))
      end function

      function sde(j,k,l,m)
        use nmvar
        implicit none
        real*8 :: sde 
        integer*4 :: j,k,l,m
        sde = eijk*(xdd(j,l)*xgde(k,m)*vc(l,2,1)
     &            +(xdd(j,l)*xgee(k,m)+xde(j,l)*xgde(k,m))*vc(l,2,2))
      end function


      function see(j,k,l,m)
        use nmvar
        implicit none
        real*8 :: see 
        integer*4 :: j,k,l,m
        see = eijk*(xde(j,l)*xgde(k,m)*vc(l,2,1)
     &            +(xde(j,l)*xgee(k,m)+xee(j,l)*xgde(k,m))*vc(l,2,2))
      end function


      function udd(j,k,l,m)
        use nmvar
        implicit none
        real*8 udd 
        integer*4 :: j,k,l,m
        udd = eijk*(xgdd(j,l)*xgdd(k,m)*vc(l,2,1)
     &            +(xgdd(j,l)*xgde(k,m)+xgde(j,l)*xgdd(k,m))*vc(l,2,2))
      end function


      function ude(j,k,l,m)
        use nmvar
        implicit none
        real*8 ude 
        integer*4 :: j,k,l,m
        ude = eijk*(xgdd(j,l)*xgde(k,m)*vc(l,2,1)
     &            +(xgdd(j,l)*xgee(k,m)+xgde(j,l)*xgde(k,m))*vc(l,2,2))
      end function

      function uee(j,k,l,m)
        use nmvar
        implicit none
        real*8 uee 
        integer*4 :: j,k,l,m 
        uee = eijk*(xgde(j,l)*xgde(k,m)*vc(l,2,1)
     &            +(xgde(j,l)*xgee(k,m)+xgee(j,l)*xgde(k,m))*vc(l,2,2))
      end function

      end module nmhncmod

