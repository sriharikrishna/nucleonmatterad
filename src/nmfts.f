c *id* nmfts ***********************************************************
c subroutine for finding correlation functions
c in spin-isospin channels for v14 problem
c short-range behavior fixed 2/16/2016
c **********************************************************************
      subroutine nmfts(lc,ls,lt,ll,lf,no,np,nt,nv)
      use nmvar
      implicit none
      !implicit real*8 (a-h,o-z)
      !implicit integer*4 (i-n)
c
      !include "params.f"
      
      real*8 rr,vv(22),vp(12),vw(10),rv(6)
      real*8 rix(ngrid),rsx(ngrid),slsx(ngrid),slpsx(ngrid),
     & chi(ngrid),phir(10,ngrid),pm(8,ngrid),psi(8,ngrid),vx(14,ngrid),
     & rlm(8,ngrid),rlx(2,ngrid),ets(2,2),c(3,8),ca(3,8),blm(8)
      real*8 blmn,blmo,c0,chi0,chim,cm,dldifo,e1,e2,fac,fbn,fcm,ftp
      real*8 gam,gchi0,gchim,gp,gpsi0,gpsim,h,h2cp,del,dldif
      real*8 h0,hm,hp,cp,de,fbm,h2,omg,phip,pqpcs,pqqpqs,psi0
      real*8 psim,rdlsx,rllpx,rllpy,rsy,rt2,rt5,ry,sldpy
      real*8 pqpc,slpy,spsy,sly,small,trpidr,x,xfb,xft,xmn,y
      real*8 yb,ybb,yc,yq,yt,z,slsy,xx
      integer*4 i,iloop,j,jm,jp,k,kqq,l,lcx,lfh,lfk,kbb,kp
      integer*4 kq,lkf,llx,loop,ltx,m,n
      integer*4 lc,ls,lt,ll,lf,no,np,nt,nv
      real*8 havetobreaka,havetobreakb,havetobreakc,havetobreakd
      !data ets/4*0./
      ets(1,1) = 0.
      ets(1,2) = 0.
      ets(2,1) = 0.
      ets(2,2) = 0.
      rt2=sqrt(2.)
      rt5=sqrt(5.)
      lfh=lf/2
      ltx=lf*lt
      h=dt/float(ltx)
      if (np.le.100) then
        gam=1. ; del=1. ; xmn=1. ; chi=1. ; omg=1. ; ftp=1.
        if (nt.ge.4) then
          gam=tniu
          del=rho
        end if
        call setpot(np,xmn,gam,del,chi,omg,ftp,h2m,h2mcs)
      else if (np.gt.100 .and. np.le.200) then
        call setpot_chiral(np,h2m,h2mcs)
      end if
      if (nm.eq.2) h2m=h2m-h2mcs
      h2=h*h/12/h2m
c   --------------------
c   set up r,sl,phi,etc.
c   --------------------
      do i=1,ltx+1
        rx(i)=h*float(i)
        rix(i)=1/rx(i)
        rsx(i)=rx(i)*rx(i)
        x=kf*rx(i)
        xx=x*x
        if (temp.eq.0.) then
          if (x.lt.0.125) then
            slx(i)=1-xx/10+xx**2/280
            slpx(i)=(-xx/5+xx**2/70-xx**3/2520)*rix(i)
            sldpx(i)=kf**2*(-1./5+3*xx/70-xx**2/504)
            sltpx(i)=kf**2*(3*xx/35-xx**2/126+xx**3/3960)*rix(i)
          else
            y=sin(x)
            z=cos(x)
            slx(i)=3*(y-x*z)/x**3
            slpx(i)=3*kf*(-slx(i)+y/x)/x
            sldpx(i)=3*(-slpx(i)+(slx(i)+z-2*y/x)*rix(i))*rix(i)
            sltpx(i)=3*(-sldpx(i)*rx(i)+2*slpx(i)
     &       -2*(slx(i)+z-2*y/x)*rix(i)-kf*(y+2*z/x-2*y/xx))/rsx(i)
          end if
        else
          call nmhot(lkf,ltx+1,1)
        end if
        slsx(i)=slx(i)*slx(i)
        slpsx(i)=slpx(i)**2
        rllpx=rx(i)*slx(i)*slpx(i)
        rdlsx=rllpx+rsx(i)*(slpsx(i)+slx(i)*sldpx(i))
        phir(1,i)=sqrt(1.+slsx(i))*rx(i)
        phir(2,i)=sqrt(1.-slsx(i))*rx(i)
        phir(3,i)=phir(2,i)
        phir(4,i)=phir(1,i)
        phir(5,i)=phir(2,i)
        phir(6,i)=phir(1,i)
        phir(7,i)=sqrt((1./5.)*xx-rllpx)*rx(i)
        phir(8,i)=sqrt((1./5.)*xx+rllpx)*rx(i)
        phir(9,i)=sqrt((12./175.)*xx**2+(2./5.)*xx-rdlsx)*rx(i)
        phir(10,i)=sqrt((12./175.)*xx**2+(2./5.)*xx+rdlsx)*rx(i)
        pm(1,i)=h2m*(slsx(i)*slpsx(i)/(1.+slsx(i))-slpsx(i)
     &   -slx(i)*sldpx(i)-2.*slx(i)*slpx(i)*rix(i))/(1.+slsx(i))
        pm(3,i)=h2m*(slsx(i)*slpsx(i)/(1.-slsx(i))+slpsx(i)
     &   +slx(i)*sldpx(i)+2.*slx(i)*slpx(i)*rix(i))/(1.-slsx(i))
        pm(5,i)=pm(3,i)
        phip=((1./5.)*xx-.5*rdlsx)/phir(7,i)
        pm(7,i)=h2m*((phip**2*rsx(i)-(1./5.)*xx+rsx(i)*(slpsx(i)
     &   +slx(i)*sldpx(i)+(1.5*slpx(i)*sldpx(i)+.5*slx(i)*sltpx(i))
     &   *rx(i)))/phir(7,i)**2-2.*phip/phir(7,i))
        psi(1,i)=phir(1,i)
        psi(3,i)=phir(3,i)
        psi(5,i)=0.
        psi(7,i)=0.
        rr=rx(i)
        if (np.le.100) then
          call pot(0,rr,vv,vp,vw)
        else if (np.gt.100 .and. np.le.200) then
          call pot_chiral(rr,vv,vp)
        end if
        if (nm.eq.1) then
          vx(1,i)=acn*vv(1)+ast*(vv(2)-3*(vv(3)+vv(4)))
          vx(2,i)=acn*vv(1)-3*ast*(vv(2)+vv(3)-3.*vv(4))
          vx(3,i)=acn*vv(1)+ast*(vv(2)+vv(3)+vv(4))
          vx(4,i)=acn*vv(1)+ast*(vv(3)-3*(vv(2)+vv(4)))
          vx(5,i)=atn*(vv(5)+vv(6))
          vx(6,i)=atn*(vv(5)-3*vv(6))
          vx(7,i)=als*(vv(7)+vv(8))
          vx(8,i)=als*(vv(7)-3*vv(8))
          vx(9,i)=acn*vv(9)+ast*(vv(10)-3*(vv(11)+vv(12)))
          vx(10,i)=acn*vv(9)-3*ast*(vv(10)+vv(11)-3.*vv(12))
          vx(11,i)=acn*vv(9)+ast*(vv(10)+vv(11)+vv(12))
          vx(12,i)=acn*vv(9)+ast*(vv(11)-3*(vv(10)+vv(12)))
          vx(13,i)=acn*(vv(13)+vv(14))
          vx(14,i)=acn*(vv(13)-3*vv(14))
        else if (nm.eq.2) then
          vx(1,i)=acn*(vv(1)+vv(2)+2*(vv(15)-vv(19)))
     &         -3*ast*(vv(3)+vv(4)+2*vv(16))
          vx(3,i)=acn*(vv(1)+vv(2)+2*(vv(15)-vv(19)))
     &           +ast*(vv(3)+vv(4)+2*vv(16))
          vx(5,i)=atn*(vv(5)+vv(6)+2*vv(17))
          vx(7,i)=als*(vv(7)+vv(8)+2*vv(18))
          vx(9,i)=acn*(vv(9)+vv(10))-3*ast*(vv(11)+vv(12))
          vx(11,i)=acn*(vv(9)+vv(10))+ast*(vv(11)+vv(12))
          vx(13,i)=acn*(vv(13)+vv(14))
        end if
        if (nm.eq.1) then
          pm(2,i)=pm(3,i)
          pm(4,i)=pm(1,i)
          pm(6,i)=pm(1,i)
          phip=((1./5.)*xx+.5*rdlsx)/phir(8,i)
          pm(8,i)=h2m*((phip**2*rsx(i)-(1./5.)*xx-rsx(i)*(slpsx(i)
     &     +slx(i)*sldpx(i)+(1.5*slpx(i)*sldpx(i)+.5*slx(i)*sltpx(i))
     &     *rx(i)))/phir(8,i)**2-2.*phip/phir(8,i))
          psi(2,i)=phir(2,i)
          psi(4,i)=phir(4,i)
          psi(6,i)=0.
          psi(8,i)=0.
          rlm(2,i)=vx(2,i)+vx(10,i)*(phir(7,i)/phir(2,i))**2
          rlx(2,i)=vx(4,i)+(vx(12,i)+(2./3.)*vx(14,i))
     &     *(phir(8,i)/phir(4,i))**2
        end if
        rlm(1,i)=vx(1,i)+vx(9,i)*(phir(8,i)/phir(1,i))**2
        rlx(1,i)=vx(3,i)+(vx(11,i)+(2./3.)*vx(13,i))
     &   *(phir(7,i)/phir(3,i))**2
      end do
      rr=0.
      if (np.le.100) then
        call pot(1,rr,vv,vp,vw)
        if (nm.eq.1) then
          rv(1)=acn*vv(1)+ast*(vv(2)-3*(vv(3)+vv(4)))
          rv(4)=acn*vv(1)+ast*(vv(3)-3*(vv(2)+vv(4)))
          rv(6)=atn*(vv(5)-3*vv(6))
        else if (nm.eq.2) then
          rv(1)=acn*(vv(1)+vv(2)+2*(vv(15)-vv(19)))
     &       -3*ast*(vv(3)+vv(4)+2*vv(16))
        end if
      else
        rv(:)=0.
      end if
c   ---------------------------
c   single-channel psi equation
c   ---------------------------
      do k=1,2,nm
        kp=9-k
        kq=k+8
        blm(k)=0.
        small=1.e-10
        lcx=lf*lc
        if (k.eq.2) lcx=lf*ls
        havetobreaka = 0.0
        do loop=1,20
          do i=1,lcx
            rlm(k,i)=blm(k)
          end do
          psim=0.
          fcm=3*psi(k,1)/phir(k,1)-3*psi(k,2)/phir(k,2)
     &         +psi(k,3)/phir(k,3)
          if (k.eq.1) gpsim=h2*rt2*rv(k)*fcm
          if (k.eq.2) gpsim=h2*h2m*2*kf*fcm/rt5
          psi0=psi(k,1)
          gpsi0=h2*(vx(k,1)-pm(k,1)-blm(k)
     &             +vx(kq,1)*(phir(kp,1)/phir(k,1))**2)*psi0
          do j=2,lcx+1
            gp=h2*(vx(k,j)-pm(k,j)-blm(k)
     &            +vx(kq,j)*(phir(kp,j)/phir(k,j))**2)
            psi(k,j)=(2*psi0+10*gpsi0-psim+gpsim)/(1-gp)
            psim=psi0
            gpsim=gpsi0
            psi0=psi(k,j)
            gpsi0=gp*psi0
          end do
          dldif=phir(k,lcx)*(psi(k,lcx+1)-psi(k,lcx-1))
     &         -psi(k,lcx)*(phir(k,lcx+1)-phir(k,lcx-1))
          if (loop.eq.1) then
            dldifo=dldif
            blmo=blm(k)
            blm(k)=(-1.)**k
          else
            if (abs(dldifo-dldif).le.small) then
              havetobreaka = 1.0
            else
              if(havetobreaka.eq.0.0) then
                blmn=(dldifo*blm(k)-dldif*blmo)/(dldifo-dldif)
                dldifo=dldif
                blmo=blm(k)
                blm(k)=blmn
              end if
            end if
          end if
        end do
        fac=phir(k,lcx)/psi(k,lcx)
        do j=1,lcx
          psi(k,j)=fac*psi(k,j)
        end do
        psi(k,lcx+1)=phir(k,lcx+1)
      end do
c -----------------------------
c coupled-channel psi equations
c -----------------------------
      do l=3,4,nm
        m=l+2
        n=l+4
        k=l-2
        kqq=l+6
        kq=l+8
        kbb=l+10
        blm(l)=0.
        blm(m)=0.
        blm(n)=0.
        lcx=lf*ls
        if (k.eq.2) lcx=lf*lc
        llx=lcx
        do i=1,lcx
          rlm(l,i)=0.
        end do
        do i=lcx+1,ltx+1
          rlm(l,i)=rlx(l-2,i)
        end do
        do i=1,ltx
          rlm(m,i)=0.
        end do
        pqpcs=(phir(n,ltx+1)/phir(l,ltx+1))**2
        rlm(m,ltx+1)=vx(m,ltx+1)-(1./12.)*pqpcs*vx(kbb,ltx+1)
        do i=1,llx
          rlm(n,i)=0.
        end do
        do i=llx+1,ltx+1
          rlm(n,i)=vx(n,i)-.5*vx(kbb,i)
        end do
        do iloop=1,10
c ---------------
c central channel
c ---------------
          havetobreakb=0
          do loop=1,20
            do i=1,lcx
              rlm(l,i)=blm(l)
            end do
            psim=0.
            chim=0.
            if (iloop.eq.1.and.loop.eq.1) then
              psi0=8*h
              psi(l,1)=psi0
              fcm=0
            else
              psi0=psi(l,1)
              fcm=3*psi(l,1)/phir(l,1)-3*psi(l,2)/phir(l,2)
     &             +psi(l,3)/phir(l,3)
            end if
            if (l.eq.3) gpsim=h2*h2m*2*kf*fcm/rt5
            if (l.eq.4) gpsim=h2*rt2*rv(l)*fcm
            gchim=gpsim
            hm=0.
            chi0=psi0
            chi(1)=chi0
            yc=vx(l,1)-rlm(l,1)
            yt=vx(m,1)-rlm(m,1)
            yb=vx(n,1)-rlm(n,1)
            yq=vx(kq,1)
            ybb=vx(kbb,1)
            pqpc=phir(n,1)/phir(l,1)
            pqpcs=pqpc*pqpc
            gpsi0=h2*(yc+pqpcs*(yq+(2./3.)*ybb)-pm(l,1))*psi0
            gchi0=h2*(yc+pqpcs*(yq+(2./3.)*ybb)-pm(l,1))*chi0
            h0=h2*((8.*yt-(2./3.)*pqpcs*ybb)*psi(m,1)
     &            +(2./3.)*pqpc*(yb-.5*ybb)*psi(n,1))
            do j=2,lcx+1
              yc=vx(l,j)-rlm(l,j)
              yt=vx(m,j)-rlm(m,j)
              yb=vx(n,j)-rlm(n,j)
              yq=vx(kq,j)
              ybb=vx(kbb,j)
              pqpc=phir(n,j)/phir(l,j)
              pqpcs=pqpc*pqpc
              gp=h2*(yc+pqpcs*(yq+(2./3.)*ybb)-pm(l,j))
              hp=h2*((8.*yt-(2./3.)*pqpcs*ybb)*psi(m,j)
     &              +(2./3.)*pqpc*(yb-.5*ybb)*psi(n,j))
              psi(l,j)=(2*psi0+10*gpsi0-psim+gpsim+hp+10*h0+hm)/(1-gp)
              psim=psi0
              gpsim=gpsi0
              hm=h0
              psi0=psi(l,j)
              gpsi0=gp*psi0
              h0=hp
              chi(j)=(2*chi0+10*gchi0-chim+gchim)/(1-gp)
              chim=chi0
              gchim=gchi0
              chi0=chi(j)
              gchi0=gp*chi0
            end do
            fac=(phir(l,lcx)-psi(l,lcx))/chi(lcx)
            do j=1,lcx+1
              psi(l,j)=psi(l,j)+fac*chi(j)
            end do
            dldif=phir(l,lcx)*(psi(l,lcx+1)-psi(l,lcx-1))
     &           -psi(l,lcx)*(phir(l,lcx+1)-phir(l,lcx-1))
            if (loop.eq.1) then
              dldifo=dldif
              blmo=blm(l)
              blm(l)=blmo+(-1.)**(l-1)
            else
              if (abs(dldifo-dldif).le.small) then
                havetobreakb = 1.0
              else
                if(havetobreakb.eq.0.0) then
                  blmn=(dldifo*blm(l)-dldif*blmo)/(dldifo-dldif)
                  dldifo=dldif
                  blmo=blm(l)
                  blm(l)=blmn
                end if
              end if
            end if
          end do
          psi(l,lcx+1)=phir(l,lcx+1)
c --------------
c tensor channel
c --------------
          if (.not.(nv.le.4)) then
          havetobreakc = 0.0
          do loop=1,20
            do i=lcx+1,ltx+1
              xft=psi(m,i)/phir(m,i)
              xfb=psi(n,i)/phir(n,i)
              rlm(l,i)=rlx(l-2,i)+8.*(vx(m,i)-rlm(m,i))*xft
     &         +(2./3.)*((vx(n,i)-rlm(n,i)-.5*vx(kbb,i))*xfb
     &         -vx(kbb,i)*xft)*(phir(n,i)/phir(l,i))**2
            end do
            do i=1,ltx
              rlm(m,i)=blm(m)
            end do
            do i=llx+1,ltx+1
              xft=psi(m,i)/phir(m,i)
              rlm(n,i)=vx(n,i)-.5*vx(kbb,i)*(1.-4.*xft)/(1.-xft)
            end do
            psim=0.
            chim=0.
            gpsim=0.
            gchim=0.
            if (iloop.eq.1.and.loop.eq.1) then
              psi0=h
              psi(m,1)=psi0
            else
              psi0=psi(m,1)
            end if
            fcm=3*psi(l,1)/phir(l,1)-3*psi(l,2)/phir(l,2)
     &           +psi(l,3)/phir(l,3)
            if (l.eq.3) hm=0.
            if (l.eq.4) hm=h2*rv(m)*fcm*rt2
            chi0=psi0
            chi(1)=chi0
            yc=vx(l,1)-rlm(l,1)
            yt=vx(m,1)-rlm(m,1)
            yb=vx(n,1)-rlm(n,1)
            yq=vx(kq,1)
            ybb=vx(kbb,1)
            pqpc=phir(n,1)/phir(l,1)
            pqpcs=pqpc*pqpc
            gpsi0=h2*(h2m*6./rsx(1)+yc-2.*yt-3.*yb+6.*yq+9.*ybb
     &       +pqpcs*(yq+(5./6.)*ybb)-pm(m,1))*psi0
            gchi0=h2*(h2m*6./rsx(1)+yc-2.*yt-3.*yb+6.*yq+9.*ybb
     &       +pqpcs*(yq+(5./6.)*ybb)-pm(m,1))*chi0
            h0=h2*((yt-(1./12.)*pqpcs*ybb)*psi(l,1)
     &          -(1./12.)*pqpc*(yb-2.*ybb)*psi(n,1))
            do j=2,ltx+1
              yc=vx(l,j)-rlm(l,j)
              yt=vx(m,j)-rlm(m,j)
              yb=vx(n,j)-rlm(n,j)
              yq=vx(kq,j)
              ybb=vx(kbb,j)
              pqpc=phir(n,j)/phir(l,j)
              pqpcs=pqpc*pqpc
              gp=h2*(h2m*6./rsx(j)+yc-2.*yt-3.*yb+6.*yq+9.*ybb
     &              +pqpcs*(yq+(5./6.)*ybb)-pm(m,j))
              hp=h2*((yt-(1./12.)*pqpcs*ybb)*psi(l,j)
     &            -(1./12.)*pqpc*(yb-2.*ybb)*psi(n,j))
              psi(m,j)=(2*psi0+10*gpsi0-psim+gpsim+hp+10*h0+hm)/(1-gp)
              psim=psi0
              gpsim=gpsi0
              hm=h0
              psi0=psi(m,j)
              gpsi0=gp*psi0
              h0=hp
              chi(j)=(2*chi0+10*gchi0-chim+gchim)/(1-gp)
              chim=chi0
              gchim=gchi0
              chi0=chi(j)
              gchi0=gp*chi0
            end do
            fac=-psi(m,ltx)/chi(ltx)
            do j=1,ltx+1
              psi(m,j)=psi(m,j)+fac*chi(j)
            end do
            dldif=phir(m,ltx)*(psi(m,ltx+1)-psi(m,ltx-1))
     &           -psi(m,ltx)*(phir(m,ltx+1)-phir(m,ltx-1))
            if (loop.eq.1) then
              dldifo=dldif
              blmo=blm(m)
              blm(m)=blmo+(-1.)**(l-1)
            else
              if ((abs(dldifo-dldif).le.small)) then
                havetobreakc = 1.0
              else
                if(havetobreakc.eq.0.0) then
                  blmn=(dldifo*blm(m)-dldif*blmo)/(dldifo-dldif)
                  dldifo=dldif
                  blmo=blm(m)
                  blm(m)=blmn
                end if
              end if
            end if
          end do
          psi(m,ltx+1)=0.
c ------------------
c spin-orbit channel
c ------------------
          if (.not.(nv.le.6)) then
          havetobreakd = 0.0
          do loop=1,20
            do i=lcx+1,ltx+1
              xft=psi(m,i)/phir(m,i)
              xfb=psi(n,i)/phir(n,i)
              rlm(l,i)=rlx(l-2,i)+8*(vx(m,i)-rlm(m,i))*xft
     &         +(2./3.)*((vx(n,i)-rlm(n,i)-.5*vx(kbb,i))*xfb
     &         -vx(kbb,i)*xft)*(phir(n,i)/phir(l,i))**2
            end do
            do i=1,llx
              rlm(n,i)=blm(n)
            end do
            psim=0
            chim=0
            if (iloop.eq.1.and.loop.eq.1) then
              psi0=8*h
              psi(n,1)=psi0
              fbm=0
            else
              psi0=psi(n,1)
              fbm=3*psi(n,1)/phir(n,1)-3*psi(n,2)/phir(n,2)
     &             +psi(n,3)/phir(n,3)
            end if
            if (l.eq.3) gpsim=h2*h2m*2*kf*fbm*rt2/rt5
            if (l.eq.4) gpsim=0
            gchim=gpsim
            hm=0
            chi0=psi0
            chi(1)=chi0
            yc=vx(l,1)-rlm(l,1)
            yt=vx(m,1)-rlm(m,1)
            yb=vx(n,1)-rlm(n,1)
            yq=vx(kq,1)
            ybb=vx(kbb,1)
            pqpc=phir(n,1)/phir(l,1)
            pqpcs=pqpc*pqpc
            pqqpqs=(phir(kqq,1)/phir(n,1))**2
            gpsi0=h2*(yc-yt-.5*yb+pqqpqs*(yq+ybb)-pm(n,1))*psi0
            gchi0=h2*(yc-yt-.5*yb+pqqpqs*(yq+ybb)-pm(n,1))*chi0
            h0=h2*pqpc*((yb-.5*ybb)*psi(l,1)-(yb-2*ybb)*psi(m,1))
            do j=2,llx+1
              yc=vx(l,j)-rlm(l,j)
              yt=vx(m,j)-rlm(m,j)
              yb=vx(n,j)-rlm(n,j)
              yq=vx(kq,j)
              ybb=vx(kbb,j)
              pqpc=phir(n,j)/phir(l,j)
              pqpcs=pqpc*pqpc
              pqqpqs=(phir(kqq,j)/phir(n,j))**2
              gp=h2*(yc-yt-.5*yb+pqqpqs*(yq+ybb)-pm(n,j))
              hp=h2*pqpc*((yb-.5*ybb)*psi(l,j)-(yb-2*ybb)*psi(m,j))
              psi(n,j)=(2*psi0+10*gpsi0-psim+gpsim+hp+10*h0+hm)/(1-gp)
              psim=psi0
              gpsim=gpsi0
              hm=h0
              psi0=psi(n,j)
              gpsi0=gp*psi0
              h0=hp
              chi(j)=(2*chi0+10*gchi0-chim+gchim)/(1-gp)
              chim=chi0
              gchim=gchi0
              chi0=chi(j)
              gchi0=gp*chi0
            end do
            fac=-psi(n,llx)/chi(llx)
            do j=1,llx+1
              psi(n,j)=psi(n,j)+fac*chi(j)
            end do
            dldif=phir(n,llx)*(psi(n,llx+1)-psi(n,llx-1))
     &           -psi(n,llx)*(phir(n,llx+1)-phir(n,llx-1))
            if (loop.eq.1) then
              dldifo=dldif
              blmo=blm(n)
              blm(n)=blmo+(-1)**k
            else
              if (abs(dldifo-dldif).gt.small) then
                if(havetobreakd.eq.0.0) then
                  blmn=(dldifo*blm(n)-dldif*blmo)/(dldifo-dldif)
                  dldifo=dldif
                  blmo=blm(n)
                  blm(n)=blmn
                end if
              else
                 havetobreakd = 1.0
              end if
            end if
          end do
          psi(n,llx+1)=0.
          end if
          end if
        end do
      end do
c print ================================================================
      if (.not.(no.eq.0)) then
      write(nlog,950) ltx,lt
      write(nout,950) ltx,lt
  950 format(/4x,'wave equations solved in t,s projected channels'
     & /4x,'fine grid has ',i3,' points in dt, main grid has ',i2
     & ,' points in dt')
      write(nlog,951) blm
      write(nout,951) blm
  951 format (/4x,'lambda'/4x,'1s',6x,'1p',4x,'3p(c)',3x,'3s(c)',3x
     & ,'3p(t)',3x,'3s(t)',3x,'3p(b)',3x,'3s(b)'/8f8.3)
      if (no.ge.3) then
        write(nout,952)
  952   format(/4x,'correlations on fine grid:'/4x,'r',7x,'1s',7x,'1p'
     &         ,7x,'3p(c)',4x,'3s(c)',4x,'3p(t)',4x,'3s(t)'
     &         ,4x,'3p(b)',4x,'3s(b)')
        write(nout,953) (rx(j),(psi(k,j)/phir(k,j),k=1,8),j=1,ltx)
  953   format(f8.5,8f9.5)
      end if
      if (no.ge.4) then
        open(unit=16,file='fnmout',status='unknown',form='unformatted')
        write(16) (rx(j),(psi(k,j)/phir(k,j),k=1,8),j=1,ltx)
      end if
      end if
c ======================================================================
c -------------------
c projections for r<d
c -------------------
      do i=1,lt
        j=lf*i-lfh
        jp=j+1
        jm=j-1
        rr=rx(j)
        do k=1,8,nm
          blm(k)=rlm(k,j)
          cp=psi(k,jp)/phir(k,jp)
          c0=psi(k,j)/phir(k,j)
          cm=psi(k,jm)/phir(k,jm)
          c(1,k)=c0
          c(2,k)=.5*(cp-cm)/h
          c(3,k)=(cp-2*c0+cm)/h**2+2*c(2,k)*rix(j)
        end do
        do k=1,2,nm
          l=k+2
          m=k+4
          n=k+6
          pqpc=phir(n,j)/phir(l,j)
          ets(k,1)=ets(k,1)+blm(k)*psi(k,j)**2
          ets(k,2)=ets(k,2)+psi(l,j)*(blm(l)*psi(l,j)+16*blm(m)*psi(m,j)
     &     +(4./3.)*pqpc*blm(n)*psi(n,j))+psi(m,j)*(8*(blm(l)-2*blm(m)
     &     -3*blm(n))*psi(m,j)-(4./3.)*pqpc*blm(n)*psi(n,j))
     &     +(2./3.)*psi(n,j)**2*(blm(l)-blm(m)-.5*blm(n))
        end do
        do k=1,3
          ca(k,1)=.25*(3*c(k,3)+c(k,1))
          ca(k,3)=.25*(c(k,3)-c(k,1))
          ca(k,5)=c(k,5)
          ca(k,7)=c(k,7)
          if (.not.(nm.eq.2)) then
          ca(k,2)=.25*ca(k,1)-.0625*(3*c(k,4)+c(k,2))
          ca(k,4)=.25*ca(k,3)-.0625*(c(k,4)-c(k,2))
          ca(k,6)=.25*(ca(k,5)-c(k,6))
          ca(k,8)=.25*(ca(k,7)-c(k,8))
          ca(k,1)=.75*ca(k,1)+.0625*(3*c(k,4)+c(k,2))
          ca(k,3)=.75*ca(k,3)+.0625*(c(k,4)-c(k,2))
          ca(k,5)=.25*(3*ca(k,5)+c(k,6))
          ca(k,7)=.25*(3*ca(k,7)+c(k,8))
          end if
        end do
        do k=1,8,nm
          f(i,k)=ca(1,k)
          fp(i,k)=ca(2,k)
          if (k.le.4) then
            fds(i,k)=ca(3,k)
          else if (k.ge.7) then
            fds(i,k)=ca(3,k)+2*ca(2,k)*rix(j)
          else
            fds(i,k)=ca(3,k)-6*ca(1,k)/rsx(j)
          end if
        end do
        if (np.le.100) then
c -------------------------------------
c special case to use av8' correlations
c but evaluate av18 potentials
c -------------------------------------
c         if (np.eq.10 .and. nv.eq.14) then
c           call setpot(9,xmn,gam,del,chi,omg,ftp,h2m,h2mcs)
c         end if
c ----------------------------------------
          call pot(0,rr,vv,vp,vw)
        else if (np.gt.100 .and. np.le.200) then
          call pot_chiral(rr,vv,vp)
        end if
        do l=1,nv,nm
          v(i,l)=vv(l)
        end do
        r(i)=rx(j)
        ri(i)=rix(j)
        rs(i)=rsx(j)
        x=kf*r(i)
        xx=x*x
        y=sin(x)
        z=cos(x)
        sl(i)=3*(y-x*z)/x**3
        sls(i)=sl(i)**2
        slp(i)=3*kf*(-sl(i)+y/x)/x
        slps(i)=slp(i)**2
        sldp(i)=3*(-slp(i)+(sl(i)+z-2*y/x)*ri(i))*ri(i)
        sltp(i)=3*(-sldp(i)*r(i)+2*slp(i)
     &         -2*(sl(i)+z-2*y/x)*ri(i)-kf*(y+2*z/x-2*y/xx))/rs(i)
        if (.not.(nm.eq.1)) then
        do l=1,nv,2
          v(i,l)=v(i,l)+vv(l+1)
          v(i,l+1)=0
        end do
        v(i,1)=v(i,1)+2*(vv(15)-vv(19))
        v(i,3)=v(i,3)+2*vv(16)
        v(i,5)=v(i,5)+2*vv(17)
        v(i,7)=v(i,7)+2*vv(18)
        end if
      end do
c -------------------
c projections for r>d
c -------------------
      do i=lt+1,lgrid
        ry=h*float(lf*i-lfh)
        rsy=ry*ry
        x=kf*ry
        xx=x*x
        y=sin(x)
        z=cos(x)
        sly=3*(y-x*z)/x**3
        slsy=sly*sly
        slpy=3*kf*(-sly/x+y/xx)
        rllpy=ry*sly*slpy
        rr=ry
        if (np.le.100) then
          call pot(0,rr,vv,vp,vw)
        else if (np.gt.100 .and. np.le.200) then
          call pot_chiral(rr,vv,vp)
        end if
        ets(1,1)=ets(1,1)+rsy*((vv(1)+vv(2)-3*(vv(3)+vv(4)))*(1+slsy)
     &   +(vv(9)+vv(10)-3*(vv(11)+vv(12)))*(.2*xx+rllpy))
        ets(1,2)=ets(1,2)+rsy*((vv(1)+vv(2)+vv(3)+vv(4))*(1-slsy)+(vv(9)
     &   +vv(10)+vv(11)+vv(12)+2*(vv(13)+vv(14))/3)*(.2*xx-rllpy))
        if (nm.eq.1) then
        ets(2,1)=ets(2,1)+rsy*((vv(1)-3*(vv(2)+vv(3)-3*vv(4)))*(1-slsy)
     &   +(vv(9)-3*(vv(10)+vv(11)-3*vv(12)))*(.2*xx-rllpy))
        ets(2,2)=ets(2,2)+rsy*((vv(1)+vv(3)-3*(vv(2)+vv(4)))*(1+slsy)
     &   +(vv(9)+vv(11)-3*(vv(10)+vv(12))+2*(vv(13)-3*vv(14))/3)
     &   *(.2*xx+rllpy))
        end if
        do l=1,nv,nm
          v(i,l)=vv(l)
        end do
        if (nm.eq.2) then
          do l=1,nv,2
            v(i,l)=v(i,l)+vv(l+1)
            v(i,l+1)=0
          end do
          v(i,1)=v(i,1)+2*(vv(15)-vv(19))
          v(i,3)=v(i,3)+2*vv(16)
          v(i,5)=v(i,5)+2*vv(17)
          v(i,7)=v(i,7)+2*vv(18)
        end if
        f(i,1)=1.
        r(i)=ry
        ri(i)=1/ry
        rs(i)=rsy
        sl(i)=sly
        sls(i)=slsy
        slp(i)=slpy
        slps(i)=slpy**2
        sldpy=3*(-slpy/ry+(sly+z-2*y/x)/rsy)
        sldp(i)=sldpy
        sltp(i)=3*(-sldpy*ry+2*slpy-2*(sly+z-2*y/x)/ry
     &   -kf*(y+2*z/x-2*y/xx))/rsy
      end do
      evx=.25_8*ets(1,1)+.75_8*ets(1,2)
      if (nm.eq.1) evx=.75_8*evx+.0625_8*ets(2,1)+.1875_8*ets(2,2)
      e1=.3*h2m*kf**2
      trpidr=2*rho*pi*lf*h
      ets(1,1)=ets(1,1)*trpidr*.25
      ets(1,2)=ets(1,2)*trpidr*.75
      e2=ets(1,1)+ets(1,2)
      if (nm.eq.1) then
        ets(1,1)=.75*ets(1,1)
        ets(1,2)=.75*ets(1,2)
        ets(2,1)=ets(2,1)*trpidr*.0625
        ets(2,2)=ets(2,2)*trpidr*.1875
      e2=ets(1,1)+ets(1,2)+ets(2,1)+ets(2,2)
      end if
      evx=e2-trpidr*evx
c write ================================================================
      if (no.ge.1) then
        write(nout,954) e1,e2,ets,evx
  954   format (/4x,'c1',4x,'c2(ts)',2x,'c2(10)',2x,'c2(00)',2x,'c2(11)'
     &          ,2x,'c2(01)',4x,'cvx'/7f8.3)
      end if
c ======================================================================
      return
      end
