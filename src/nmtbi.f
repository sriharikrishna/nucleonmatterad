
      module nmtbimod
      use nmvar
      implicit none


      private:: rrr,rry,ryr,rrz,zzz,zzr,zrz,rbybr,rbrtr,rbrty,rbytr,
     &ybrtr,rtrtr,rtrty,rtytr,zbzbz,zbztz
      
      integer*4, parameter :: nu=4/nm
      integer*4, parameter :: n3s=5-nm
      integer*4, parameter :: n3t=7-nm
      integer*4, parameter :: nlog=0
      integer*4, parameter :: nin=5
      integer*4, parameter :: nout=6
      real*8 afe(6),w3vm(10),w3vx(6,2),v3(2,2)
     &,rcc(lgrid),rdd(lgrid),rde(lgrid),ree(lgrid),sdd(lgrid)
     &,ycc(lgrid),ydd(lgrid),yde(lgrid),yee(lgrid),ybcc(lgrid)
     &,ybdd(lgrid),ybde(lgrid),ybee(lgrid),ytcc(lgrid),ytee(lgrid)
     &,zdd(lgrid,6),zde(lgrid,6),zee(lgrid,6),zcc(lgrid,6)
     &,zpdd(lgrid,6),zpde(lgrid,6),zpee(lgrid,6),zpcc(lgrid,6)
     &,zbdd(lgrid,8),zbde(lgrid,8),zbee(lgrid,8),zbcc(lgrid,8)
     &,ztee(lgrid,6)

      public :: nmtbi
      
      contains

c *id* nmtbi ***********************************************************
c nmtbi
c subprogram for three-body integrations
c ----------------------------------------------------------------------
      subroutine nmtbi(lt,lg,le,l3,nie,no,nt)
      implicit none
      integer*4 :: lt,lg,le,l3,nie,no,nt
      integer*4 :: i,j,k,l,kb,ltd,fc2,fc2p,fl2,ffl,fl2p,fflp,ijk,ikj
      integer*4 :: jik,kji,ka
      real*8 :: vid,vie,vip,vpd,vpe,vpp,vfd,vfe,vfp,pa,pap,pb,pc,pd,qv
      real*8 :: v0,ac,qvi,w3vma,w3vmb,x,y,z,qttt,qtts,qtst,qstt,acx,u1
      real*8 :: u2,u3,u4,uf1,uf2,uf3,uf4,uf5,up1,up2,up3,qu,fpi,fpj
      real*8 :: qui,u4a,u4b,u4c,u4d,qpts,qtps,qpps,qspt,qpst,qtpt,qptt
      real*8 :: qppt
      real*8 :: acex !entry point
c

c -------------------------
c statement functions
c error in rbrtr fixed 7/07
c -------------------------
      do 5 l=1,6
    5 afe(l)=acex(l,1,l)
      do 10 i=1,10
   10 w3vm(i)=0
c ----------------------
c set up r,y,z functions
c ----------------------
      ltd=min(2*max(lt,lg,l3),lgrid)
      do 50 i=1,ltd
c ---------------------------
c rxx, rbxx, rtxx definitions
c ---------------------------
        rdd(i)=f(i,1)**2*gx(i)
        rde(i)=rdd(i)*gy(i)
        ree(i)=rdd(i)*(gy(i)**2+gz(i)-gl(i)**2/nu)
        rcc(i)=-rdd(i)*gl(i)/nu
        sdd(i)=rdd(i)-1
        ydd(i)=0
        yde(i)=0
        yee(i)=0
        ycc(i)=0
        ybdd(i)=0
        ybde(i)=0
        ybee(i)=0
        ybcc(i)=0
        ytee(i)=0
        ytcc(i)=0
        zdd(i,1)=0
        zde(i,1)=0
        zee(i,1)=0
        zcc(i,1)=0
        zbdd(i,1)=0
        zbde(i,1)=0
        zbee(i,1)=0
        zbcc(i,1)=0
        ztee(i,1)=0
        fc2=f(i,1)**2
        fc2p=2*f(i,1)*fp(i,1)
        do 40 l=1+nm,6,nm
          vid=vc(l,1,1)
          vie=vc(l,1,2)
          vip=vc(l,1,3)
          vpd=vc(l,2,1)
          vpe=vc(l,2,2)
          vpp=vc(l,2,3)
          vfd=vc(l,3,1)
          vfe=vc(l,3,2)
          vfp=vc(l,3,3)
          fl2=f(i,l)**2
          ffl=2*f(i,1)*f(i,l)
          fl2p=2*f(i,l)*fp(i,l)
          fflp=2*(f(i,1)*fp(i,l)+fp(i,1)*f(i,l))
c ---------------------------
c yxx, ybxx, ytxx definitions
c ---------------------------
          ydd(i)=ydd(i)+aa(l)*( fl2*vfd**2
     &                   +ffl*gdd(i,l)*vpd**2
     &                   +.5*fc2*gdd(i,l)**2*vpd**2 )*gx(i)
          yde(i)=yde(i)+aa(l)*( (fl2*vfd*vfp
     &                    +ffl*gdd(i,l)*vpd*vpp
     &                    +.5*fc2*gdd(i,l)**2*vpd*vpp)*gy(i)
     &                   +(ffl*gde(i,l)*vpe*vpd
     &                    +fc2*gdd(i,l)*gde(i,l)*vpe*vpd) )*gx(i)
          yee(i)=yee(i)+aa(l)*( (fl2*vfp**2
     &                    +ffl*gdd(i,l)*vpp**2
     &                    +.5*fc2*gdd(i,l)**2*vpp**2)*(gy(i)**2+gz(i))
     &                 +2*(ffl*gde(i,l)*vpe*vpp
     &                    +fc2*gdd(i,l)*gde(i,l)*vpe*vpp)*gy(i)
     &                 +ffl*(gee(i,l)-af(l)*gl(i)**2/nu)*vpe**2
     &                 +fc2*gdd(i,l)*(gee(i,l)-af(l)*gl(i)**2/nu)*vpe**2
     &                 +fc2*gde(i,l)**2*vpe**2
     &               +2*fc2*af(l)*(gca(i,l)+gcb(i,l))*gl(i)*vpe )*gx(i)
     &           -afe(l)*fl2*gl(i)**2*vfe**2*gx(i)/nu
          ycc(i)=ycc(i)-aa(l)*af(l)*( ffl*gl(i)*vpe**2/nu
     &                         +fc2*gdd(i,l)*gl(i)*vpe**2/nu
     &                         -fc2*(gca(i,l)+gcb(i,l))*vpe )*gx(i)
     &           -afe(l)*fl2*gl(i)*vfe**2*gx(i)/nu
          ybdd(i)=ybdd(i)+aa(l)*( fl2p*vfd**2
     &                     +fflp*gdd(i,l)*vpd**2
     &                     +.5*fc2p*gdd(i,l)**2*vpd**2 )*gx(i)
          ybde(i)=ybde(i)+aa(l)*( (fl2p*vfd*vfp
     &                      +fflp*gdd(i,l)*vpd*vpp
     &                      +.5*fc2p*gdd(i,l)**2*vpd*vpp)*gy(i)
     &                     +(fflp*gde(i,l)*vpe*vpd
     &                      +fc2p*gdd(i,l)*gde(i,l)*vpe*vpd) )*gx(i)
          ybee(i)=ybee(i)+aa(l)*( (fl2p*vfp**2
     &                     +fflp*gdd(i,l)*vpp**2
     &                     +.5*fc2p*gdd(i,l)**2*vpp**2)*(gy(i)**2+gz(i))
     &                 +2*(fflp*gde(i,l)*vpe*vpp
     &                    +fc2p*gdd(i,l)*gde(i,l)*vpe*vpp)*gy(i)
     &                +fflp*(gee(i,l)-af(l)*gl(i)**2/nu)*vpe**2
     &                +fc2p*gdd(i,l)*(gee(i,l)-af(l)*gl(i)**2/nu)*vpe**2
     &                +fc2p*gde(i,l)**2*vpe**2
     &              +2*fc2p*af(l)*(gca(i,l)+gcb(i,l))*gl(i)*vpe )*gx(i)
     &             -afe(l)*fl2p*gl(i)**2*vfe**2*gx(i)/nu
          ybcc(i)=ybcc(i)-aa(l)*af(l)*( fflp*gl(i)*vpe**2/nu
     &                           +fc2p*gdd(i,l)*gl(i)*vpe**2/nu
     &                           -fc2p*(gca(i,l)+gcb(i,l))*vpe )*gx(i)
     &             -afe(l)*fl2p*gl(i)*vfe**2*gx(i)/nu
          ytee(i)=ytee(i)-aa(l)*af(l)*( ffl*gl(i)*vpe**2/nu
     &                      +fc2*gdd(i,l)*gl(i)*vpe**2/nu )*slp(i)*gx(i)
     &         +aa(l)*af(l)*(fc2*gx(i)-1)*(gca(i,l)+gcb(i,l))*vpe*slp(i)
     &         -afe(l)*fl2*gl(i)*slp(i)*vfe**2*gx(i)/nu
          ytcc(i)=ytcc(i)+aa(l)*af(l)*( ffl*vpe**2
     &                           +fc2*gdd(i,l)*vpe**2 )*slp(i)*gx(i)
     &             +afe(l)*fl2*slp(i)*vfe**2*gx(i)
c ---------------------------
c zxx, zbxx, ztxx definitions
c ---------------------------
          zdd(i,l)=(ffl+fc2*gdd(i,l))*gx(i)*vid**2
          zde(i,l)=(ffl+fc2*gdd(i,l))*gx(i)*gy(i)*vid*vip
     &       +fc2*gde(i,l)*gx(i)*vid*vie
          zee(i,l)=(ffl+fc2*gdd(i,l))*gx(i)*(gy(i)**2+gz(i))*vip**2
     &     +2*fc2*gde(i,l)*gx(i)*gy(i)*vip*vie
     &       +fc2*(gee(i,l)-af(l)*gl(i)**2/nu)*gx(i)*vie**2
          zcc(i,l)=rcc(i)*af(l)*vie**2
          zpdd(i,l)=(ffl+fc2*gdd(i,l))*gx(i)*vpd**2
          zpde(i,l)=(ffl+fc2*gdd(i,l))*gx(i)*gy(i)*vpd*vpp
     &       +fc2*gde(i,l)*gx(i)*vpd*vpe
          zpee(i,l)=(ffl+fc2*gdd(i,l))*gx(i)*(gy(i)**2+gz(i))*vpp**2
     &     +2*fc2*gde(i,l)*gx(i)*gy(i)*vpp*vpe
     &       +fc2*(gee(i,l)-af(l)*gl(i)**2/nu)*gx(i)*vpe**2
          zpcc(i,l)=rcc(i)*af(l)*vpe**2
          zbdd(i,l)=(fflp+fc2p*gdd(i,l))*gx(i)*vpd**2
          zbde(i,l)=(fflp+fc2p*gdd(i,l))*gx(i)*gy(i)*vpd*vpp
     &         +fc2p*gde(i,l)*gx(i)*vpd*vpe
          zbee(i,l)=(fflp+fc2p*gdd(i,l))*gx(i)*(gy(i)**2+gz(i))*vpp**2
     &       +2*fc2p*gde(i,l)*gx(i)*gy(i)*vpp*vpe
     &         +fc2p*(gee(i,l)-af(l)*gl(i)**2/nu)*gx(i)*vpe**2
          ztee(i,l)=-af(l)*rdd(i)*gl(i)*slp(i)*vpe**2/nu
          if (l.ge.5) then
            zbdd(i,l+2)=ffl*ri(i)*gx(i)*vpd**2
            zbde(i,l+2)=ffl*ri(i)*gx(i)*gy(i)*vpd*vpp
            zbee(i,l+2)=ffl*ri(i)*gx(i)*(gy(i)**2+gz(i))*vpp**2
          end if
   40   continue
   50 continue
c ------------------------
c calculate w3
c nt = 1 Urbana type
c      2 Tucson
c      3 Brazil
c      4 density-dependent
c s = 4/nm
c ------------------------
      if (nt.eq.0.or.nt.ge.4) then
        go to 500
      else if (nt.eq.1) then
        pap=4*tnix*.5*(s-1)
        pb=4*tnia*.5*(s-1)
        pd=-16*tnic*.375*(s-2)
      else if (nt.eq.2) then
        pa=0.1102
        pb=-0.2517*.5*(s-1)
        pc=0.09754
        pd=0.2937*.375*(s-2)
      else if (nt.eq.3) then
        pa=0.1023
        pb=-0.2230*.5*(s-1)
        pc=0
        pd=0.2987*.375*(s-2)
      end if
      qv=4*(rho*pi)**2*dr**3
c ------------------------------------------
c integration loops: i=r(mo),j=r(mn),k=r(no)
c ------------------------------------------
      do 300 i=1,l3
        do 290 j=1,l3
          ka=iabs(i-j)+1
          kb=min(i+j-1,ltd)
          if (nt.eq.1) then
            v0=tniu*tpi2(i)*tpi2(j)
            ac=pap*xt1(i)*xt1(j)
            v3(1,1)=ypi(i)*ypi(j)
            v3(1,2)=ypi(i)*tpi(j)
            v3(2,1)=tpi(i)*ypi(j)
            v3(2,2)=tpi(i)*tpi(j)
          else if (nt.eq.2.or.nt.eq.3) then
            v0=0
            ac=(pa-2*pc)*xt1(i)*xt1(j)+pc*(xt0(i)*xt1(j)+xt1(i)*xt0(j))
            v3(1,1)=xt3(i)*xt3(j)
            v3(1,2)=xt3(i)*xt2(j)
            v3(2,1)=xt2(i)*xt3(j)
            v3(2,2)=xt2(i)*xt2(j)
          end if
          if (tniu.eq.0.) go to 240
c ----------------------
c central v3 integration
c w3vm(1) := diagram 3.1
c w3vm(2) := diagram 3.2
c w3vm(3) := diagram 3.3
c ----------------------
          do 210 k=ka,kb
            ijk=index(i,j,k)
            ikj=index(i,k,j)
            jik=index(j,i,k)
            kji=index(k,j,i)
            qvi=qv*r(i)*r(j)*r(k)
            w3vm(1)=w3vm(1)+qvi*v0*rrr(i,j,k,ijk,ikj,jik,kji)
            w3vm(2)=w3vm(2)+qvi*v0*(rry(i,j,k,ijk,ikj,jik,kji)
     &                           +2*ryr(i,j,k,ijk,ikj,jik,kji))
  210     continue
          do 230 l=1,2,nm
            w3vma=0
            w3vmb=0
            do 220 k=ka,kb
              ijk=index(i,j,k)
              ikj=index(i,k,j)
              jik=index(j,i,k)
              kji=index(k,j,i)
              x=xtheta(kji)
              y=ytheta(kji)
              z=ztheta(kji)
              qvi=qv*r(i)*r(j)*r(k)
              qttt=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
              qtts=3*x**2-1
              qtst=1.5*y**2-.5
              qstt=1.5*z**2-.5
              w3vma=w3vma+qvi*v0*
     &         (aa(l)*zzz(i,j,k,l,l,l,ijk,ikj,jik,kji)/vc(l,1,1)**3
     &         +aa(l+2)*(zzz(i,j,k,l+2,l+2,l+2,ijk,ikj,jik,kji)
     &                  +zzz(i,j,k,l+4,l+4,l+2,ijk,ikj,jik,kji)*qtts)
     &                  /vc(l+2,1,1)**3)
              w3vmb=w3vmb+qvi*v0*
     &         (aa(l+4)*(zzz(i,j,k,l+2,l+4,l+4,ijk,ikj,jik,kji)*qstt
     &                  +zzz(i,j,k,l+4,l+2,l+4,ijk,ikj,jik,kji)*qtst
     &                  +zzz(i,j,k,l+4,l+4,l+4,ijk,ikj,jik,kji)*qttt)
     &                  /vc(l+4,1,1)**3)
  220       continue
            w3vm(3)=w3vm(3)+w3vma+w3vmb
  230     continue
c ----------------------------
c 2pi-exchange v3 integrations
c ----------------------------
c ----------------------------
c anticommutator zzr zzz terms
c   w3vm(4) := diagram 2.2
c   w3vm(5)  ~ diagram 2.4
c   w3vm(9) := diagram 2.2-s
c   w3vm(10) ~ diagram 2.4-s
c ----------------------------
  240     do 250 k=ka,kb
            ijk=index(i,j,k)
            ikj=index(i,k,j)
            jik=index(j,i,k)
            kji=index(k,j,i)
            x=xtheta(kji)
            y=ytheta(kji)
            z=ztheta(kji)
            qvi=qv*r(i)*r(j)*r(k)
            qttt=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
            qtts=3*x**2-1
            qtst=1.5*y**2-.5
            qstt=1.5*z**2-.5
            w3vm(4)=w3vm(4)+qvi*pb*
     &       ( (  6*v3(1,1)+6*qtts*v3(2,2) )
     &             *zzr(i,j,k,n3s,n3s,ijk,ikj,jik,kji)
     &        +( 12*v3(2,1)+6*qtts*(v3(1,2)+v3(2,2)) )
     &             *zzr(i,j,k,n3t,n3s,ijk,ikj,jik,kji)
     &        +( 12*v3(1,2)+6*qtts*(v3(2,1)+v3(2,2)) )
     &             *zzr(i,j,k,n3s,n3t,ijk,ikj,jik,kji)
     &        +( 24*v3(2,2)+6*qtts*(v3(1,1)+v3(1,2)+v3(2,1)+v3(2,2)) )
     &             *zzr(i,j,k,n3t,n3t,ijk,ikj,jik,kji))/vc(n3s,1,1)
            w3vm(9)=w3vm(9)+qvi*pb*18*ac*x
     &          *(  zzr(i,j,k,n3s,n3s,ijk,ikj,jik,kji)
     &           +2*zzr(i,j,k,n3t,n3s,ijk,ikj,jik,kji)
     &           +2*zzr(i,j,k,n3s,n3t,ijk,ikj,jik,kji)
     &           +4*zzr(i,j,k,n3t,n3t,ijk,ikj,jik,kji))/vc(n3s,1,1)
            w3vm(5)=w3vm(5)+qvi*pb*
     &       ((24*(v3(1,1)+acx)+24*qtts*(v3(2,2)+acx))
     &             *zdd(i,n3s)*zdd(j,n3s)
     &       +(-24*(v3(2,1)+acx)-12*qtts*(v3(1,2)+v3(2,2)+2*acx))
     &             *zdd(i,n3t)*zdd(j,n3s)
     &       +(-24*(v3(1,2)+acx)-12*qtts*(v3(2,1)+v3(2,2)+2*acx))
     &             *zdd(i,n3s)*zdd(j,n3t)
     &       +(72*(v3(2,2)+acx)-12*qtts*(v3(2,1)+v3(1,2)-2*v3(1,1)))
     &             *zdd(i,n3t)*zdd(j,n3t))*zdd(k,n3s)/vc(n3s,1,1)**3
            w3vm(5)=w3vm(5)+qvi*pb*
     &       ((-24*qtst*(v3(2,1)+acx)-24*qstt*(v3(1,2)+acx)
     &         -24*qttt*(v3(2,2)+acx))*zdd(i,n3s)*zdd(j,n3s)
     &       +((12*qtts+24*(qstt+qtst-1))*(v3(1,2)+acx)
     &        +(-12*qtts+24*(qttt-qtst+1))*(v3(2,2)+acx)
     &        -24*qtst*(v3(1,1)-2*v3(2,1)-acx))*zdd(i,n3t)*zdd(j,n3s)
     &       +((12*qtts+24*(qstt+qtst-1))*(v3(2,1)+acx)
     &        +(-12*qtts+24*(qttt-qstt+1))*(v3(2,2)+acx)
     &        -24*qstt*(v3(1,1)-2*v3(1,2)-acx))*zdd(i,n3s)*zdd(j,n3t)
     &       +((24*qtts-48*(qttt+1))*(v3(2,2)+acx)
     &        +(-12*qtts+24*(qttt-qtst+1))*(v3(2,1)+acx)
     &        +(-12*qtts+24*(qttt-qstt+1))*(v3(1,2)+acx)
     &        -24*qttt*(v3(1,1)+acx))*zdd(i,n3t)*zdd(j,n3t))
     &             *zdd(k,n3t)/vc(n3t,1,1)**3
  250     continue
          if (nm.eq.2) go to 290
c ------------------------
c commutator zzr zrz terms
c ------------------------
          do 260 k=ka,kb
            ijk=index(i,j,k)
            ikj=index(i,k,j)
            jik=index(j,i,k)
            kji=index(k,j,i)
            x=xtheta(kji)
            y=ytheta(kji)
            z=ztheta(kji)
            qvi=qv*r(i)*r(j)*r(k)
            qttt=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
            qtts=3*x**2-1
            qtst=1.5*y**2-.5
            qstt=1.5*z**2-.5
            w3vm(6)=w3vm(6)+qvi*pd*
     &       ((-12*v3(1,1)+6*qtts*v3(2,2))
     &             *zzr(i,j,k,n3s,n3s,ijk,ikj,jik,kji)
     &       +(-24*v3(2,1)+6*qtts*(v3(1,2)+v3(2,2)))
     &             *zzr(i,j,k,n3t,n3s,ijk,ikj,jik,kji)
     &       +(-24*v3(1,2)+6*qtts*(v3(2,1)+v3(2,2)))
     &             *zzr(i,j,k,n3s,n3t,ijk,ikj,jik,kji)
     &       +(6*qtts*(v3(1,1)+v3(1,2)+v3(2,1)+v3(2,2))-48*v3(2,2))
     &             *zzr(i,j,k,n3t,n3t,ijk,ikj,jik,kji))/vc(n3s,1,1)
c --------------------------------------
c following expression corrected 5/30/99
c --------------------------------------
            w3vm(7)=w3vm(7)+2*qvi*pd*
     &       ((-12*v3(1,1)+6*qtts*v3(2,2))
     &             *zrz(i,j,k,n3s,n3s,ijk,ikj,jik,kji)
     &       +(12*v3(2,1)+6*qtts*(v3(1,2)-2*v3(2,2)))
     &             *zrz(i,j,k,n3t,n3s,ijk,ikj,jik,kji)
     &       +(12*qtst*v3(2,1)-24*qstt*v3(1,2)+12*qttt*v3(2,2))
     &             *zrz(i,j,k,n3s,n3t,ijk,ikj,jik,kji)
     &       +(12*qtst*(v3(1,1)-2*v3(2,1))+12*qttt*(v3(1,2)-2*v3(2,2))
     &        +24*qstt*v3(2,2))
     &             *zrz(i,j,k,n3t,n3t,ijk,ikj,jik,kji))/vc(n3s,1,1)
  260     continue
c ------------------------------
c direct three-correlation terms
c ------------------------------
          do 270 k=ka,kb
            kji=index(k,j,i)
            x=xtheta(kji)
            y=ytheta(kji)
            z=ztheta(kji)
            acx=ac/x
            qvi=qv*r(i)*r(j)*r(k)
            qttt=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
            qtts=3*x**2-1
            qtst=1.5*y**2-.5
            qstt=1.5*z**2-.5
            w3vm(8)=w3vm(8)+qvi*pd*
     &       ((-12*v3(1,1)+6*qtts*v3(2,2))
     &             *zdd(i,n3s)*zdd(j,n3s)
     &       +(12*v3(2,1)+6*qtts*(v3(2,2)-2*v3(1,2)))
     &             *zdd(i,n3t)*zdd(j,n3s)
     &       +(12*v3(1,2)+6*qtts*(v3(2,2)-2*v3(2,1)))
     &             *zdd(i,n3s)*zdd(j,n3t)
     &       +(-36*v3(2,2)+6*qtts*(v3(1,1)+v3(2,1)+v3(1,2)))
     &             *zdd(i,n3t)*zdd(j,n3t))*zdd(k,n3s)/vc(n3s,1,1)**3
            w3vm(8)=w3vm(8)+qvi*pd*
     &       ((12*qtst*v3(2,1)+12*qstt*v3(1,2)
     &        -6*(qtts+2*(qtst+qstt-1))*v3(2,2))*zdd(i,n3s)*zdd(j,n3s)
     &       +(12*qtst*(v3(1,1)-2*v3(2,1))+12*qttt*v3(1,2)
     &        +6*(qtts-2*(qttt-qtst+1))*v3(2,2))*zdd(i,n3t)*zdd(j,n3s)
     &       +(12*qstt*(v3(1,1)-2*v3(1,2))+12*qttt*v3(2,1)
     &        +6*(qtts-2*(qttt-qstt+1))*v3(2,2))*zdd(i,n3s)*zdd(j,n3t)
     &       +(-6*(qtts+2*(qtst+qstt-1))*v3(1,1)
     &        +6*(qtts-2*(qttt-qtst+1))*v3(2,1)
     &        +6*(qtts-2*(qttt-qstt+1))*v3(1,2)
     &        -12*(qtts-2*(qttt+1))*v3(2,2))*zdd(i,n3t)*zdd(j,n3t))
     &             *zdd(k,n3t)/vc(n3t,1,1)**3
  270     continue
  290   continue
  300 continue
      w3v0=w3vm(1)+w3vm(2)+w3vm(3)
c w3vx not defined any where so this is commented out 2/13/16
c     w3vx(1,1)=0
c     fac=tnix/tniu
c     w3v1=-.25*(fac*w3v0+w3vx(2,1)+w3vx(3,1)+w3vx(4,1)
c    &               +w3vx(2,2)+w3vx(3,2)+w3vx(4,2))
c -----------------------------------------------------------
      w3va=w3va+w3vm(4)+w3vm(5)
      w3vc=w3vm(6)+w3vm(7)+w3vm(8)
c print ================================================================
      if (no.eq.0) go to 500
      write(nlog,973)
      write(nout,973)
  973 format(/4x,'w3v0(g3):',15x,'w3va(g3):',7x,'w3vc(g3):'/4x,'rrr'
     &,5x,'yrr',5x,'zzz',5x,'zzr',5x,'zzz',5x,'zzr',5x,'zrz',5x,'zzz')
      write(nlog,918) w3vm
      write(nout,918) w3vm
  918 format(8f8.3)
      write(nlog,972)
      write(nout,972)
  972 format(/4x,'w3vx:rrz,zzr'
     &       /4x,'c',7x,'t',7x,'s',7x,'st',6x,'tn',6x,'tnt')
      write(nlog,916) w3vx
      write(nout,916) w3vx
  916 format(6f8.3)
      write(nlog,974) w3v0,w3v1,w3va,w3vc
      write(nout,974) w3v0,w3v1,w3va,w3vc
  974 format(/4x,'w3v0',4x,'w3v1',4x,'w3va',4x,'w3vc'/4f8.3)
c ======================================================================
c --------------------
c calculate u, uf & up
c --------------------
  500 u1=0
      u2=0
      u3=0
      u4=0
      uf1=0
      uf2=0
      uf3=0
      uf4=0
      uf5=0
      up1=0
      up2=0
      up3=0
      qu=-h2m*(rho*pi)**2*dr**3
c ------------------------------------------
c integration loops: i=r(mo),j=r(mn),k=r(no)
c ------------------------------------------
      do 550 i=1,lg
        fpi=2*fp(i,1)/f(i,1)
        do 540 j=1,lg
          fpj=2*fp(j,1)/f(j,1)
          ka=iabs(i-j)+1
          kb=min(i+j-1,ltd)
c ---------------------------------
c rrr, rry, ryr, & yrr integrations
c ---------------------------------
          do 510 k=ka,kb
            ijk=index(i,j,k)
            ikj=index(i,k,j)
            jik=index(j,i,k)
            kji=index(k,j,i)
            qui=qu*r(i)*r(j)*r(k)*xtheta(kji)
            u1=u1+qui*fpi*fpj*rrr(i,j,k,ijk,ikj,jik,kji)
            u2=u2+qui*fpi*fpj*rry(i,j,k,ijk,ikj,jik,kji)
            u3=u3+2*qui*fpi*rbybr(i,j,k,ijk,ikj,jik,kji)
            uf1=uf1+4*qui*fpi*slp(j)*rbrtr(i,j,k,ijk,ikj,jik,kji)
            uf2=uf2+4*qui*fpi*slp(j)*rbrty(i,j,k,ijk,ikj,jik,kji)
            uf3=uf3+4*qui*slp(j)*ybrtr(i,j,k,ijk,ikj,jik,kji)
            uf4=uf4+4*qui*fpi*rbytr(i,j,k,ijk,ikj,jik,kji)
            up1=up1+2*qui*slp(i)*slp(j)*rtrtr(i,j,k,ijk,kji)/nu
            up2=up2+2*qui*slp(i)*slp(j)*rtrty(i,j,k,ijk,kji)/nu
            up3=up3+4*qui*slp(i)*rtytr(i,j,k,ijk,kji)/nu
  510     continue
c ----------------
c zzz integrations
c ----------------
          do 530 l=1,2,nm
            u4a=0
            u4b=0
            u4c=0
            u4d=0
            do 520 k=ka,kb
              ijk=index(i,j,k)
              ikj=index(i,k,j)
              jik=index(j,i,k)
              kji=index(k,j,i)
              x=xtheta(kji)
              y=ytheta(kji)
              z=ztheta(kji)
              qui=qu*r(i)*r(j)*r(k)*x
              qtst=1.5*y**2-.5
              qstt=1.5*z**2-.5
              qtts=3*x**2-1
              qttt=-4.5*x*y*z-1.5*(x**2+y**2+z**2)+1
              qpts=6*(1-x**2)
              qtps=6*(1-x**2)
              qpps=12*x**2
              qspt=-3*(y/x+z)*z
              qpst=-3*(z/x+y)*y
              qtpt=9*x*y*z+4.5*y**2+3*(x**2+z**2)-1.5*y*z/x-3
              qptt=9*x*y*z+4.5*z**2+3*(x**2+y**2)-1.5*y*z/x-3
              qppt=-18*x*y*z-6*x**2-9*(y**2+z**2)-4.5*y*z/x+4.5
              u4a=u4a+qui*
     &         (aa(l)*zbzbz(i,j,k,l,l,l,ijk,ikj,jik,kji)/vc(l,2,1)**3
     &         +aa(l+2)*(zbzbz(i,j,k,l+2,l+2,l+2,ijk,ikj,jik,kji)
     &                  +zbzbz(i,j,k,l+4,l+4,l+2,ijk,ikj,jik,kji)*qtts)
     &                  /vc(l+2,2,1)**3)
              u4b=u4b+qui*
     &         (aa(l+4)*(zbzbz(i,j,k,l+2,l+4,l+4,ijk,ikj,jik,kji)*qstt
     &                  +zbzbz(i,j,k,l+4,l+2,l+4,ijk,ikj,jik,kji)*qtst
     &                  +zbzbz(i,j,k,l+6,l+2,l+4,ijk,ikj,jik,kji)*qpst
     &                  +zbzbz(i,j,k,l+4,l+4,l+4,ijk,ikj,jik,kji)*qttt)
     &                  /vc(l+4,2,1)**3)
              u4c=u4c+qui*
     &         (aa(l+2)*(zbzbz(i,j,k,l+4,l+6,l+2,ijk,ikj,jik,kji)*qtps
     &                  +zbzbz(i,j,k,l+6,l+4,l+2,ijk,ikj,jik,kji)*qpts
     &                  +zbzbz(i,j,k,l+6,l+6,l+2,ijk,ikj,jik,kji)*qpps)
     &                  /vc(l+2,2,1)**3)
              u4d=u4d+qui*
     &         (aa(l+4)*(zbzbz(i,j,k,l+2,l+6,l+4,ijk,ikj,jik,kji)*qspt
     &                  +zbzbz(i,j,k,l+6,l+4,l+4,ijk,ikj,jik,kji)*qptt
     &                  +zbzbz(i,j,k,l+4,l+6,l+4,ijk,ikj,jik,kji)*qtpt
     &                  +zbzbz(i,j,k,l+6,l+6,l+4,ijk,ikj,jik,kji)*qppt)
     &                  /vc(l+4,2,1)**3)
              uf5=uf5+4*qui*
     &         (aa(l)*zbztz(i,j,k,l,l,l,ijk,jik)/vc(l,2,1)**3
     &         +aa(l+2)*zbztz(i,j,k,l+2,l+2,l+2,ijk,jik)/vc(l+2,2,1)**3
     &         +aa(l+4)*(zbztz(i,j,k,l+4,l+2,l+4,ijk,jik)*qtst
     &                  +zbztz(i,j,k,l+6,l+2,l+4,ijk,jik)*qpst)
     &                  /vc(l+4,2,1)**3)
  520       continue
            u4=u4+u4a+u4b+u4c+u4d
  530     continue
  540   continue
  550 continue
c print ================================================================
  600 if (no.eq.0) go to 999
      write(nlog,970) u1,u2,u3,u4
      write(nout,970) u1,u2,u3,u4
  970 format(/4x,'u1',6x,'u2',6x,'u3',6x,'u4'/4f8.3)
      write(nlog,971) uf1,uf2,uf3,uf4,uf5
      write(nout,971) uf1,uf2,uf3,uf4,uf5
  971 format(/4x,'uf1',5x,'uf2',5x,'uf3',5x,'uf4',5x,'uf5'/5f8.3)
      write(nlog,975) up1,up2,up3
      write(nout,975) up1,up2,up3
  975 format(/4x,'up1',5x,'up2',5x,'up3'/3f8.3)
c ======================================================================
  999 u=u1+u2+u3+u4
      uf=uf1+uf2+uf3+uf4+uf5
      up=up1+up2+up3
      return
      end subroutine nmtbi

      function rrr(i,j,k,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: rrr
        integer*4 :: i,j,k,ijk,ikj,jik,kji
        rrr =
     &  rdd(i)*rdd(j)*rdd(k)*(1+sddd(ijk)+sdde(ijk)+sdde(jik)
     &                         +sdde(kji)+sdee(ijk)+sdee(ikj)
     &                         +sdee(kji)+seee(ijk))
     & +rde(i)*rdd(j)*rdd(k)*(2+2*sddd(ijk)+sdde(kji)+sdde(jik)
     &                         +2*sdde(ijk)+sdee(ijk)+sdee(ikj))
     & +rdd(i)*rde(j)*rdd(k)*(2+2*sddd(ijk)+sdde(kji)+sdde(ijk)
     &                         +2*sdde(jik)+sdee(ijk)+sdee(kji))
     & +rdd(i)*rdd(j)*rde(k)*(2+2*sddd(ijk)+sdde(jik)+sdde(ijk)
     &                         +2*sdde(kji)+sdee(ikj)+sdee(kji))
     & +ree(i)*rdd(j)*rdd(k)*(1+sddd(ijk)+sdde(ijk))
     & +rdd(i)*ree(j)*rdd(k)*(1+sddd(ijk)+sdde(jik))
     & +rdd(i)*rdd(j)*ree(k)*(1+sddd(ijk)+sdde(kji))
     & +(rde(i)*rde(j)*rdd(k)+rde(i)*rdd(j)*rde(k)+rdd(i)*rde(j)*rde(k))
     &  *(3+3*sddd(ijk)+sdde(ijk)+sdde(jik)+sdde(kji))
     & +(rdd(i)*(rde(j)*ree(k)+ree(j)*rde(k))
     &  +rdd(j)*(rde(k)*ree(i)+ree(k)*rde(i))
     &  +rdd(k)*(rde(i)*ree(j)+ree(i)*rde(j))
     &  +2*rde(i)*rde(j)*rde(k))*(1+sddd(ijk))
     &  -2*nu*(rcc(i)*rcc(j)*rcc(k)*(1+sddd(ijk))
     &        +rcc(i)*rdd(j)*rdd(k)*(sccd(ijk)+scce(ijk))
     &        +rdd(i)*rcc(j)*rdd(k)*(sccd(jik)+scce(jik))
     &        +rdd(i)*rdd(j)*rcc(k)*(sccd(kji)+scce(kji))
     &        +rcc(i)*(rde(j)*rdd(k)+rdd(j)*rde(k))*sccd(ijk)
     &        +rcc(j)*(rde(k)*rdd(i)+rdd(k)*rde(i))*sccd(jik)
     &        +rcc(k)*(rde(i)*rdd(j)+rdd(i)*rde(j))*sccd(kji)
     &        +rdd(i)*rcc(j)*rcc(k)*sccd(ijk)
     &        +rcc(i)*rdd(j)*rcc(k)*sccd(jik)
     &        +rcc(i)*rcc(j)*rdd(k)*sccd(kji))
      end function rrr

      function rry(i,j,k,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: rry
        integer*4 :: i,j,k,ijk,ikj,jik,kji
        rry=
     &  rdd(i)*rdd(j)*ydd(k)*(1+sddd(ijk)+sdde(ijk)+sdde(jik)
     &                         +sdde(kji)+sdee(ijk)+sdee(ikj)
     &                         +sdee(kji)+seee(ijk))
     & +rde(i)*rdd(j)*ydd(k)*(2+2*sddd(ijk)+sdde(kji)+sdde(jik)
     &                         +2*sdde(ijk)+sdee(ijk)+sdee(ikj))
     & +rdd(i)*rde(j)*ydd(k)*(2+2*sddd(ijk)+sdde(kji)+sdde(ijk)
     &                         +2*sdde(jik)+sdee(ijk)+sdee(kji))
     & +rdd(i)*rdd(j)*yde(k)*(2+2*sddd(ijk)+sdde(jik)+sdde(ijk)
     &                         +2*sdde(kji)+sdee(ikj)+sdee(kji))
     & +ree(i)*rdd(j)*ydd(k)*(1+sddd(ijk)+sdde(ijk))
     & +rdd(i)*ree(j)*ydd(k)*(1+sddd(ijk)+sdde(jik))
     & +rdd(i)*rdd(j)*yee(k)*(1+sddd(ijk)+sdde(kji))
     & +(rde(i)*rde(j)*ydd(k)+rde(i)*rdd(j)*yde(k)+rdd(i)*rde(j)*yde(k))
     &  *(3+3*sddd(ijk)+sdde(ijk)+sdde(jik)+sdde(kji))
     & +(rdd(i)*(rde(j)*yee(k)+ree(j)*yde(k))
     &  +rdd(j)*(yde(k)*ree(i)+yee(k)*rde(i))
     &  +ydd(k)*(rde(i)*ree(j)+ree(i)*rde(j))
     &  +2*rde(i)*rde(j)*yde(k))*(1+sddd(ijk))
     &  -2*nu*(rcc(i)*rcc(j)*ycc(k)*(1+sddd(ijk))
     &        +rcc(i)*rdd(j)*ydd(k)*(sccd(ijk)+scce(ijk))
     &        +rdd(i)*rcc(j)*ydd(k)*(sccd(jik)+scce(jik))
     &        +rdd(i)*rdd(j)*ycc(k)*(sccd(kji)+scce(kji))
     &        +rcc(i)*(rde(j)*ydd(k)+rdd(j)*yde(k))*sccd(ijk)
     &        +rcc(j)*(yde(k)*rdd(i)+ydd(k)*rde(i))*sccd(jik)
     &        +ycc(k)*(rde(i)*rdd(j)+rdd(i)*rde(j))*sccd(kji)
     &        +rdd(i)*rcc(j)*ycc(k)*sccd(ijk)
     &        +rcc(i)*rdd(j)*ycc(k)*sccd(jik)
     &        +rcc(i)*rcc(j)*ydd(k)*sccd(kji))
      end function rry

      function ryr(i,j,k,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: ryr
        integer*4 :: i,j,k,ijk,ikj,jik,kji
        ryr=
     &  rdd(i)*ydd(j)*rdd(k)*(1+sddd(ijk)+sdde(ijk)+sdde(jik)
     &                         +sdde(kji)+sdee(ijk)+sdee(ikj)
     &                         +sdee(kji)+seee(ijk))
     & +rde(i)*ydd(j)*rdd(k)*(2+2*sddd(ijk)+sdde(kji)+sdde(jik)
     &                         +2*sdde(ijk)+sdee(ijk)+sdee(ikj))
     & +rdd(i)*yde(j)*rdd(k)*(2+2*sddd(ijk)+sdde(kji)+sdde(ijk)
     &                         +2*sdde(jik)+sdee(ijk)+sdee(kji))
     & +rdd(i)*ydd(j)*rde(k)*(2+2*sddd(ijk)+sdde(jik)+sdde(ijk)
     &                         +2*sdde(kji)+sdee(ikj)+sdee(kji))
     & +ree(i)*ydd(j)*rdd(k)*(1+sddd(ijk)+sdde(ijk))
     & +rdd(i)*yee(j)*rdd(k)*(1+sddd(ijk)+sdde(jik))
     & +rdd(i)*ydd(j)*ree(k)*(1+sddd(ijk)+sdde(kji))
     & +(rde(i)*yde(j)*rdd(k)+rde(i)*ydd(j)*rde(k)+rdd(i)*yde(j)*rde(k))
     &  *(3+3*sddd(ijk)+sdde(ijk)+sdde(jik)+sdde(kji))
     & +(rdd(i)*(yde(j)*ree(k)+yee(j)*rde(k))
     &  +ydd(j)*(rde(k)*ree(i)+ree(k)*rde(i))
     &  +rdd(k)*(rde(i)*yee(j)+ree(i)*yde(j))
     &  +2*rde(i)*yde(j)*rde(k))*(1+sddd(ijk))
     &  -2*nu*(rcc(i)*ycc(j)*rcc(k)*(1+sddd(ijk))
     &        +rcc(i)*ydd(j)*rdd(k)*(sccd(ijk)+scce(ijk))
     &        +rdd(i)*ycc(j)*rdd(k)*(sccd(jik)+scce(jik))
     &        +rdd(i)*ydd(j)*rcc(k)*(sccd(kji)+scce(kji))
     &        +rcc(i)*(yde(j)*rdd(k)+ydd(j)*rde(k))*sccd(ijk)
     &        +ycc(j)*(rde(k)*rdd(i)+rdd(k)*rde(i))*sccd(jik)
     &        +rcc(k)*(rde(i)*ydd(j)+rdd(i)*yde(j))*sccd(kji)
     &        +rdd(i)*ycc(j)*rcc(k)*sccd(ijk)
     &        +rcc(i)*ydd(j)*rcc(k)*sccd(jik)
     &        +rcc(i)*ycc(j)*rdd(k)*sccd(kji))
      end function ryr

      function rrz(i,j,k,l,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: rrz
        integer*4 :: i,j,k,l,ijk,ikj,jik,kji
        rrz=
     &  rdd(i)*rdd(j)*zdd(k,l)*(1+sddd(ijk)+sdde(ijk)+sdde(jik)
     &                           +sdde(kji)+sdee(ijk)+sdee(ikj)
     &                           +sdee(kji)+seee(ijk))
     & +rde(i)*rdd(j)*zdd(k,l)*(2+2*sddd(ijk)+sdde(kji)+sdde(jik)
     &                           +2*sdde(ijk)+sdee(ijk)+sdee(ikj))
     & +rdd(i)*rde(j)*zdd(k,l)*(2+2*sddd(ijk)+sdde(kji)+sdde(ijk)
     &                           +2*sdde(jik)+sdee(ijk)+sdee(kji))
     & +rdd(i)*rdd(j)*zde(k,l)*(2+2*sddd(ijk)+sdde(jik)+sdde(ijk)
     &                           +2*sdde(kji)+sdee(ikj)+sdee(kji))
     & +ree(i)*rdd(j)*zdd(k,l)*(1+sddd(ijk)+sdde(ijk))
     & +rdd(i)*ree(j)*zdd(k,l)*(1+sddd(ijk)+sdde(jik))
     & +rdd(i)*rdd(j)*zee(k,l)*(1+sddd(ijk)+sdde(kji))
     & +(rde(i)*rde(j)*zdd(k,l)+rde(i)*rdd(j)*zde(k,l)
     &  +rdd(i)*rde(j)*zde(k,l))
     &  *(3+3*sddd(ijk)+sdde(ijk)+sdde(jik)+sdde(kji))
     & +(rdd(i)*(rde(j)*zee(k,l)+ree(j)*zde(k,l))
     &  +rdd(j)*(zde(k,l)*ree(i)+zee(k,l)*rde(i))
     &  +zdd(k,l)*(rde(i)*ree(j)+ree(i)*rde(j))
     &  +2*rde(i)*rde(j)*zde(k,l))*(1+sddd(ijk))
     &  -2*nu*(rcc(i)*rcc(j)*rcc(k)*af(l)*(1+sddd(ijk))
     &        +rcc(i)*rdd(j)*zdd(k,l)*(sccd(ijk)+scce(ijk))
     &        +rdd(i)*rcc(j)*zdd(k,l)*(sccd(jik)+scce(jik))
     &        +rdd(i)*rdd(j)*rcc(k)*af(l)*(sccd(kji)+scce(kji))
     &        +rcc(i)*(rde(j)*zdd(k,l)+rdd(j)*zde(k,l))*sccd(ijk)
     &        +rcc(j)*(zde(k,l)*rdd(i)+zdd(k,l)*rde(i))*sccd(jik)
     &        +rcc(k)*af(l)*(rde(i)*rdd(j)+rdd(i)*rde(j))*sccd(kji)
     &        +rdd(i)*rcc(j)*rcc(k)*af(l)*sccd(ijk)
     &        +rcc(i)*rdd(j)*rcc(k)*af(l)*sccd(jik)
     &        +rcc(i)*rcc(j)*zdd(k,l)*sccd(kji))
      end function rrz

      function zzz(i,j,k,l,m,n,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: zzz
        integer*4 :: i,j,k,l,m,n,ijk,ikj,jik,kji
        zzz=
     & zdd(i,l)*zdd(j,m)*zdd(k,n)*(1+sddd(ijk)+sdde(ijk)+sdde(jik)
     &          +sdde(kji)+sdee(ijk)+sdee(ikj)+sdee(kji)+seee(ijk))
     &+zde(i,l)*zdd(j,m)*zdd(k,n)*(2+2*sddd(ijk)+sdde(kji)
     &          +sdde(jik)+2*sdde(ijk)+sdee(ijk)+sdee(ikj))
     &+zdd(i,l)*zde(j,m)*zdd(k,n)*(2+2*sddd(ijk)+sdde(kji)
     &          +sdde(ijk)+2*sdde(jik)+sdee(ijk)+sdee(kji))
     &+zdd(i,l)*zdd(j,m)*zde(k,n)*(2+2*sddd(ijk)+sdde(jik)
     &          +sdde(ijk)+2*sdde(kji)+sdee(ikj)+sdee(kji))
     & +zee(i,l)*zdd(j,m)*zdd(k,n)*(1+sddd(ijk)+sdde(ijk))
     & +zdd(i,l)*zee(j,m)*zdd(k,n)*(1+sddd(ijk)+sdde(jik))
     & +zdd(i,l)*zdd(j,m)*zee(k,n)*(1+sddd(ijk)+sdde(kji))
     & +(zde(i,l)*zde(j,m)*zdd(k,n)+zde(i,l)*zdd(j,m)*zde(k,n)
     &  +zdd(i,l)*zde(j,m)*zde(k,n))
     &  *(3+3*sddd(ijk)+sdde(ijk)+sdde(jik)+sdde(kji))
     & +(zdd(i,l)*(zde(j,m)*zee(k,n)+zee(j,m)*zde(k,n))
     &  +zdd(j,m)*(zde(k,n)*zee(i,l)+zee(k,n)*zde(i,l))
     &  +zdd(k,n)*(zde(i,l)*zee(j,m)+zee(i,l)*zde(j,m))
     &  +2*zde(i,l)*zde(j,m)*zde(k,n))*(1+sddd(ijk))
      end function zzz

      function zzr(i,j,k,l,m,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: zzr
        integer*4 :: i,j,k,l,m,ijk,ikj,jik,kji
        zzr=
     & zdd(i,l)*zdd(j,m)*rdd(k)*(1+sddd(ijk)+sdde(ijk)+sdde(jik)
     &        +sdde(kji)+sdee(ijk)+sdee(ikj)+sdee(kji)+seee(ijk))
     &+zde(i,l)*zdd(j,m)*rdd(k)*(2+2*sddd(ijk)+sdde(kji)
     &        +sdde(jik)+2*sdde(ijk)+sdee(ijk)+sdee(ikj))
     &+zdd(i,l)*zde(j,m)*rdd(k)*(2+2*sddd(ijk)+sdde(kji)
     &        +sdde(ijk)+2*sdde(jik)+sdee(ijk)+sdee(kji))
     &+zdd(i,l)*zdd(j,m)*rde(k)*(2+2*sddd(ijk)+sdde(jik)
     &        +sdde(ijk)+2*sdde(kji)+sdee(ikj)+sdee(kji))
     & +zee(i,l)*zdd(j,m)*rdd(k)*(1+sddd(ijk)+sdde(ijk))
     & +zdd(i,l)*zee(j,m)*rdd(k)*(1+sddd(ijk)+sdde(jik))
     & +zdd(i,l)*zdd(j,m)*ree(k)*(1+sddd(ijk)+sdde(kji))
     & +(zde(i,l)*zde(j,m)*rdd(k)+zde(i,l)*zdd(j,m)*rde(k)
     &  +zdd(i,l)*zde(j,m)*rde(k))
     &  *(3+3*sddd(ijk)+sdde(ijk)+sdde(jik)+sdde(kji))
     & +(zdd(i,l)*(zde(j,m)*ree(k)+zee(j,m)*rde(k))
     &  +zdd(j,m)*(rde(k)*zee(i,l)+ree(k)*zde(i,l))
     &  +rdd(k)*(zde(i,l)*zee(j,m)+zee(i,l)*zde(j,m))
     &  +2*zde(i,l)*zde(j,m)*rde(k)-2*nu*zcc(i,l)*zcc(j,m)*rcc(k))
     & *(1+sddd(ijk))
      end function zzr

      function zrz(i,j,k,l,n,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: zrz
        integer*4 :: i,j,k,l,n,ijk,ikj,jik,kji
        zrz=
     & zdd(i,l)*rdd(j)*zdd(k,n)*(1+sddd(ijk)+sdde(ijk)+sdde(jik)
     &        +sdde(kji)+sdee(ijk)+sdee(ikj)+sdee(kji)+seee(ijk))
     &+zde(i,l)*rdd(j)*zdd(k,n)*(2+2*sddd(ijk)+sdde(kji)
     &        +sdde(jik)+2*sdde(ijk)+sdee(ijk)+sdee(ikj))
     &+zdd(i,l)*rde(j)*zdd(k,n)*(2+2*sddd(ijk)+sdde(kji)
     &        +sdde(ijk)+2*sdde(jik)+sdee(ijk)+sdee(kji))
     &+zdd(i,l)*rdd(j)*zde(k,n)*(2+2*sddd(ijk)+sdde(jik)
     &        +sdde(ijk)+2*sdde(kji)+sdee(ikj)+sdee(kji))
     & +zee(i,l)*rdd(j)*zdd(k,n)*(1+sddd(ijk)+sdde(ijk))
     & +zdd(i,l)*ree(j)*zdd(k,n)*(1+sddd(ijk)+sdde(jik))
     & +zdd(i,l)*rdd(j)*zee(k,n)*(1+sddd(ijk)+sdde(kji))
     & +(zde(i,l)*rde(j)*zdd(k,n)+zde(i,l)*rdd(j)*zde(k,n)
     &  +zdd(i,l)*rde(j)*zde(k,n))
     &  *(3+3*sddd(ijk)+sdde(ijk)+sdde(jik)+sdde(kji))
     & +(zdd(i,l)*(rde(j)*zee(k,n)+ree(j)*zde(k,n))
     &  +rdd(j)*(zde(k,n)*zee(i,l)+zee(k,n)*zde(i,l))
     &  +zdd(k,n)*(zde(i,l)*ree(j)+zee(i,l)*rde(j))
     &  +2*zde(i,l)*rde(j)*zde(k,n))*(1+sddd(ijk))
      end function zrz

      function rbybr(i,j,k,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: rbybr
        integer*4 :: i,j,k,ijk,ikj,jik,kji
        rbybr =
     &  rdd(i)*ybdd(j)*rdd(k)*(1+sddd(ijk)+sdde(ijk)+sdde(jik)
     &  +sdde(kji)+sdee(ijk)+sdee(ikj)+sdee(kji)+seee(ijk))
     & +rde(i)*ybdd(j)*rdd(k)*(2+2*sddd(ijk)+sdde(kji)+sdde(jik)
     &                          +2*sdde(ijk)+sdee(ijk)+sdee(ikj))
     & +rdd(i)*ybde(j)*rdd(k)*(2+2*sddd(ijk)+sdde(kji)+sdde(ijk)
     &                          +2*sdde(jik)+sdee(ijk)+sdee(kji))
     & +rdd(i)*ybdd(j)*rde(k)*(2+2*sddd(ijk)+sdde(jik)+sdde(ijk)
     &                          +2*sdde(kji)+sdee(ikj)+sdee(kji))
     & +ree(i)*ybdd(j)*rdd(k)*(1+sddd(ijk)+sdde(ijk))
     & +rdd(i)*ybee(j)*rdd(k)*(1+sddd(ijk)+sdde(jik))
     & +rdd(i)*ybdd(j)*ree(k)*(1+sddd(ijk)+sdde(kji))
     & +(rde(i)*ybde(j)*rdd(k)+rde(i)*ybdd(j)*rde(k)
     &  +rdd(i)*ybde(j)*rde(k))
     &  *(3+3*sddd(ijk)+sdde(ijk)+sdde(jik)+sdde(kji))
     & +(rdd(i)*(ybde(j)*ree(k)+ybee(j)*rde(k))
     &  +ybdd(j)*(rde(k)*ree(i)+ree(k)*rde(i))
     &  +rdd(k)*(rde(i)*ybee(j)+ree(i)*ybde(j))
     &  +2*rde(i)*ybde(j)*rde(k))*(1+sddd(ijk))
     &  -2*nu*(rcc(i)*ybcc(j)*rcc(k)*(1+sddd(ijk))
     &        +rcc(i)*ybdd(j)*rdd(k)*(sccd(ijk)+scce(ijk))
     &        +rdd(i)*ybcc(j)*rdd(k)*(sccd(jik)+scce(jik))
     &        +rdd(i)*ybdd(j)*rcc(k)*(sccd(kji)+scce(kji))
     &        +rcc(i)*(ybde(j)*rdd(k)+ybdd(j)*rde(k))*sccd(ijk)
     &        +ybcc(j)*(rde(k)*rdd(i)+rdd(k)*rde(i))*sccd(jik)
     &        +rcc(k)*(rde(i)*ybdd(j)+rdd(i)*ybde(j))*sccd(kji)
     &        +rdd(i)*ybcc(j)*rcc(k)*sccd(ijk)
     &        +rcc(i)*ybdd(j)*rcc(k)*sccd(jik)
     &        +rcc(i)*ybcc(j)*rdd(k)*sccd(kji))
      end function rbybr

      function rbrtr(i,j,k,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: rbrtr
        integer*4 :: i,j,k,ijk,ikj,jik,kji
        rbrtr =
     &  (rdd(i)+rde(i))*rcc(j)*sdd(k)
     & +(rdd(i)*rcc(j)*rde(k)+rcc(i)*sdd(j)*rcc(k))*(1+sddd(ijk))
     &  +rdd(i)*rcc(j)*rdd(k)*(sddd(ijk)+sdde(jik))
     & +(rde(i)*rcc(j)*rdd(k)+rcc(i)*rcc(k))*sddd(ijk)
     &  +rcc(i)*(sdd(j)*rdd(k)+sdd(k))*sccd(kji)
     &  +rdd(i)*(sdd(j)*rdd(k)+sdd(k))*(sccd(jik)+scce(jik))
     & +(rdd(i)*rdd(j)*rde(k)+rde(i)*(sdd(j)*rdd(k)+sdd(k)))*sccd(jik)
     &  +rdd(i)*rdd(j)*rcc(k)*sccd(ijk)
      end function rbrtr

      function rbrty(i,j,k,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: rbrty
        integer*4 :: i,j,k,ijk,ikj,jik,kji
        rbrty =
     &   rdd(i)*rcc(j)*ydd(k)*(1+sddd(ijk)+sdde(jik))
     & +(rdd(i)*rcc(j)*yde(k)+rde(i)*rcc(j)*ydd(k)
     &  +rcc(i)*sdd(j)*ycc(k))*(1+sddd(ijk))
     &  +rcc(i)*ycc(k)*sddd(ijk)
     &  +rcc(i)*rdd(j)*ydd(k)*sccd(kji)
     &  +rdd(i)*rdd(j)*ydd(k)*(sccd(jik)+scce(jik))
     & +(rdd(i)*rdd(j)*yde(k)+rde(i)*rdd(j)*ydd(k))*sccd(jik)
     &  +rdd(i)*rdd(j)*ycc(k)*sccd(ijk)
      end function rbrty

      function rbytr(i,j,k,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: rbytr
        integer*4 :: i,j,k,ijk,ikj,jik,kji
        rbytr =
     &   rdd(i)*ytee(j)*rdd(k)*(1+sddd(ijk)+sdde(jik))
     & +(rdd(i)*ytee(j)*rde(k)+rde(i)*ytee(j)*rdd(k)
     &  +rcc(i)*ytcc(j)*rcc(k))*(1+sddd(ijk))
     &  +rcc(i)*ytcc(j)*rdd(k)*sccd(kji)
     &  +rdd(i)*ytcc(j)*rdd(k)*(sccd(jik)+scce(jik))
     & +(rdd(i)*ytcc(j)*rde(k)+rde(i)*ytcc(j)*rdd(k))*sccd(jik)
     &  +rdd(i)*ytcc(j)*rcc(k)*sccd(ijk)
      end function rbytr

      function ybrtr(i,j,k,ijk,ikj,jik,kji)
        use nmvar
        implicit none
        real*8 :: ybrtr
        integer*4 :: i,j,k,ijk,ikj,jik,kji
        ybrtr =
     &   ybdd(i)*rcc(j)*rdd(k)*(1+sddd(ijk)+sdde(jik))
     & +(ybdd(i)*rcc(j)*rde(k)+ybde(i)*rcc(j)*rdd(k)
     &  +ybcc(i)*sdd(j)*rcc(k))*(1+sddd(ijk))
     &  +ybcc(i)*rcc(k)*sddd(ijk)
     &  +ybcc(i)*(sdd(j)*rdd(k)+sdd(k))*sccd(kji)
     &  +ybdd(i)*(sdd(j)*rdd(k)+sdd(k))*(sccd(jik)+scce(jik))
     & +(ybdd(i)*rdd(j)*rde(k)+ybde(i)*(sdd(j)*rdd(k)+sdd(k)))*sccd(jik)
     &  +ybdd(i)*rdd(j)*rcc(k)*sccd(ijk)
      end function ybrtr

      function rtrtr(i,j,k,ijk,kji)
        use nmvar
        implicit none
        real*8 :: rtrtr
        integer*4 :: i,j,k,ijk,kji
        rtrtr =sdd(i)*sdd(j)*rcc(k)
     & +rdd(i)*rdd(j)*rcc(k)*sddd(ijk)
     & +rdd(i)*rdd(j)*rdd(k)*sccd(kji)
      end function rtrtr

      function rtrty(i,j,k,ijk,kji)
        use nmvar
        implicit none
        real*8 :: rtrty
        integer*4 :: i,j,k,ijk,kji
        rtrty = sdd(i)*sdd(j)*ycc(k)
     & +rdd(i)*rdd(j)*ycc(k)*sddd(ijk)
     & +rdd(i)*rdd(j)*ydd(k)*sccd(kji)
      end function rtrty

      function rtytr(i,j,k,ijk,kji)
        use nmvar
        implicit none
        real*8 :: rtytr
        integer*4 :: i,j,k,ijk,kji
        rtytr = sdd(i)*ytcc(j)*rcc(k)
     & +rdd(i)*ytcc(j)*rcc(k)*sddd(ijk)
     & +rdd(i)*ytcc(j)*rdd(k)*sccd(kji)
      end function rtytr

      function zbzbz(i,j,k,l,m,n,ijk,ikj,jik,kji)
        implicit none
        real*8 :: zbzbz
        integer*4 :: i,j,k,l,m,n,ijk,ikj,jik,kji
        zbzbz=
     & zbdd(i,l)*zbdd(j,m)*zpdd(k,n)*(1+sddd(ijk)+sdde(ijk)+sdde(jik)
     &             +sdde(kji)+sdee(ijk)+sdee(ikj)+sdee(kji)+seee(ijk))
     &+zbde(i,l)*zbdd(j,m)*zpdd(k,n)*(2+2*sddd(ijk)+sdde(kji)
     &             +sdde(jik)+2*sdde(ijk)+sdee(ijk)+sdee(ikj))
     &+zbdd(i,l)*zbde(j,m)*zpdd(k,n)*(2+2*sddd(ijk)+sdde(kji)
     &             +sdde(ijk)+2*sdde(jik)+sdee(ijk)+sdee(kji))
     &+zbdd(i,l)*zbdd(j,m)*zpde(k,n)*(2+2*sddd(ijk)+sdde(jik)
     &             +sdde(ijk)+2*sdde(kji)+sdee(ikj)+sdee(kji))
     & +zbee(i,l)*zbdd(j,m)*zpdd(k,n)*(1+sddd(ijk)+sdde(ijk))
     & +zbdd(i,l)*zbee(j,m)*zpdd(k,n)*(1+sddd(ijk)+sdde(jik))
     & +zbdd(i,l)*zbdd(j,m)*zpee(k,n)*(1+sddd(ijk)+sdde(kji))
     & +(zbde(i,l)*zbde(j,m)*zpdd(k,n)+zbde(i,l)*zbdd(j,m)*zpde(k,n)
     &  +zbdd(i,l)*zbde(j,m)*zpde(k,n))
     &  *(3+3*sddd(ijk)+sdde(ijk)+sdde(jik)+sdde(kji))
     & +(zbdd(i,l)*(zbde(j,m)*zpee(k,n)+zbee(j,m)*zpde(k,n))
     &  +zbdd(j,m)*(zpde(k,n)*zbee(i,l)+zpee(k,n)*zbde(i,l))
     &  +zpdd(k,n)*(zbde(i,l)*zbee(j,m)+zbee(i,l)*zbde(j,m))
     &  +2*zbde(i,l)*zbde(j,m)*zpde(k,n))*(1+sddd(ijk))
      end function zbzbz

      function zbztz(i,j,k,l,m,n,ijk,jik)
        implicit none
        real*8 :: zbztz
        integer*4 :: i,j,k,l,m,n,ijk,jik
        zbztz=
     & +zbdd(i,l)*ztee(j,m)*zpdd(k,n)*(1+sddd(ijk)+sdde(jik))
     & +(zbdd(i,l)*ztee(j,m)*zpde(k,n)+zbde(i,l)*ztee(j,m)*zpdd(k,n))
     & *(1+sddd(ijk))
      end function zbztz

      end module nmtbimod
