c *id* nmfkts **********************************************************
c subroutine for finding momentum-dependent correlation
c functions in partial-wave channels for v14 problem
c **********************************************************************
      subroutine nmfkts(nm,lc,ls,lt,ll,lf,no,np,nv)
      use parameters
      use wavefunc
      use consts
      use perturb
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
c
c     include 'param.f'
      parameter(kgrid=10,kxgrid=20,kxxgrid=40)
c     real*8 kf,rho,acn,ast,atn,als,cn,dt,dr,evx,h2m,h2mcs,pi,s
c     common /consts/ kf,rho,acn,ast,atn,als,cn,dt,dr,evx,h2m,h2mcs,pi,s
c     real*8 r(lgrid),rs(lgrid),sl(lgrid),sls(lgrid),slp(lgrid),
c    &       slps(lgrid),sldp(lgrid),sltp(lgrid)
c     common /rslate/ r,rs,sl,sls,slp,slps,sldp,sltp
c     real*8 f(lgrid,8),fp(lgrid,8),fds(lgrid,8),v(lgrid,14)
c     common /correl/ f,fp,fds,v
      real*8 rr,vv(18),ww(3),rv(6)
      real*8 rx(ngrid),rsx(ngrid), psis(4,ngrid,kgrid),
     &psit(4,ngrid,kgrid), vx(14,ngrid), k(kgrid), jlx(8,ngrid,kgrid),
     &jlpx(4,ngrid,kgrid),clm(2,4,kgrid),ks(kgrid),
     &psitm(3,ngrid,kgrid,2),psitp(3,ngrid,kgrid,2),almm(3,kgrid,2),
     &almp(3,kgrid,2),chi(ngrid)
      real*8 psi(8,ngrid),phir(10,ngrid),blm(8),uk(8,ngrid,kgrid),
     &utk(3,2,ngrid,kgrid),wtk(3,2,ngrid,kgrid),
     &utkdp(3,2,lgrid,kgrid),wtkdp(3,2,lgrid,kgrid),
     &unm(15,ngrid,kgrid),utnm(5,2,ngrid,kgrid),wtnm(5,2,ngrid,kgrid),
     &ukdp(8,lgrid,kgrid),unmdp(15,lgrid,kgrid),
     &utnmdp(5,2,lgrid,kgrid),wtnmdp(5,2,lgrid,kgrid)
      real*8 vtsj(15,lgrid),vtsjm(5,lgrid),vtsj0(5,lgrid),
     &vtsjp(5,lgrid)
      real*8 c2k(15,kgrid),c2nm(15,kgrid),
     &c2ck(5,2,kgrid),c2cnm(5,2,kgrid)
      common/correlx/ psi,phir,blm
      real*8 w(kgrid),wp(kgrid),kx(kxgrid),kxx(kxxgrid),
     &wxx(kxxgrid),kint(0:kxxgrid),c2(15),c2c(5,2)
      real*8 xvs(ngrid),xvd(ngrid),xvg(ngrid),xvp(ngrid),xvf(ngrid),
     &xvt(ngrid),xu(ngrid),xw(ngrid),xh(ngrid)
      rt2=sqrt(2.)
      rt5=sqrt(5.)
      lfh=lf/2
      ltx=lf*lt
      ltxp=ltx+1
      ltxm=ltx-1
      ltp=lt+1
      ltd=2*lt
      lcx=lf*lc
      lcxp=lcx+1
      lcxm=lcx-1
      llx=lf*ll
      llxp=llx+1
      llxm=llx-1
      lmax=int(10*float(lt)/dt)
      if (lmax.lt.lgrid) lmax=lgrid
      h=dt/float(ltx)
      call setpot(np,0.,rho,h2m,h2mcs)
      if (nm.eq.2) h2m=h2m-h2mcs
      h2=h*h/12/h2m
      dk=kf/float(kgrid)
c   --------------------
c   set up r,j,etc.
c   --------------------
      do i=1,ngrid
        rx(i)=h*float(i)
        rsx(i)=rx(i)*rx(i)
      end do

      do 50 i=1,ngrid
       do 55 j=1,kgrid
	k(j)=(j-0.5)*dk
	ks(j)=k(j)*k(j)
        x=k(j)*rx(i)
        xx=x*x
	xxx=x*xx
	xxxx=xx*xx
	xxxxx=x*xxxx
	xxxxxx=x*xxxxx
          if (x.lt.0.125) then
	    jlx(1,i,j)= 1.-xx/6.
	    jlx(2,i,j)= x*(1.-xx/10.)/3.
	    jlx(3,i,j)= xx*(1.-xx/14.)/15.
	    jlx(4,i,j)= xxx*(1.-xx/18.)/105.
	    jlx(5,i,j)= xxxx*(1.-xx/22.)/945.
	    jlx(6,i,j)= xxxxx*(1.-xx/26.)/10395.
	    jlx(7,i,j)= xxxxxx*(1.-xx/30.)/135135.
	    jlx(8,i,j)= xxxxxx*x*(1.-xx/34.)/2027025.
            jlpx(1,i,j)= k(j)*x*(-1.+xx/10.)/3.
      	    jlpx(2,i,j)= k(j)*(1./3.-xx/10.+xx*xx/168.)
      	    jlpx(3,i,j)= k(j)*x*(1.-xx/7.+xx*xx/170.)*2./15.
      	    jlpx(4,i,j)= k(j)*(3.*xx-5.*xx*xx/18.+7.*xx**3/792.)/105.
          else
            y=sin(x)
            z=cos(x)
	    jlx(1,i,j)= y/x
	    jlx(2,i,j)= y/xx - z/x
	    jlx(3,i,j)= (3./xx -1.)*y/x -3.*z/xx
	    jlx(4,i,j)= (15./xx - 6.)*y/xx -(15./xx -1.)*z/x
	    jlx(5,i,j)= (105./xxxx -45./xx +1.)*y/x-(105./xx -10.)*z/xx
	    jlx(6,i,j)= (945./xxxx -420./xx +15.)*y/xx
     &                 -(945./xxxx -105./xx +1. )*z/x
	    jlx(7,i,j)= (10395./xxxxxx -4725./xxxx +210./xx -1.)*y/x
     &                 -(10395./xxxx   -1260./xx   +21.)*z/xx
	    jlx(8,i,j)= (135135./xxxxxx -62370./xxxx +3150./xx-28.)*y/xx
     &                 -(135135./xxxxxx -17325./xxxx +378./xx -1.)*z/x
	    jlpx(1,i,j)= -k(j)*jlx(2,i,j)
	    jlpx(2,i,j)=  k(j)*(-2.*y/xx +2.*z/x +y)/x
	    jlpx(3,i,j)=  k(j)*((-9./xx +4.)*y/xx +(9./xx -1.)*z/x)
	    jlpx(4,i,j)=  k(j)*((-60./xx/xx +27./xx -1.)*y/x
     &                         +(60./xx -7.)*z/xx)
          end if

55     continue

        rr=rx(i)
        call pot(0,rr,vv,ww)
        if (nm.eq.1) then
          vx(1,i)=acn*vv(1)+ast*(vv(2)-3*(vv(3)+vv(4)))
          vx(2,i)=acn*vv(1)-3*ast*(vv(2)+vv(3)-3.*vv(4))
          vx(3,i)=acn*vv(1)+ast*(vv(2)+vv(3)+vv(4))
          vx(4,i)=acn*vv(1)+ast*(vv(3)-3*(vv(2)+vv(4)))
          vx(5,i)=atn*(vv(5)+vv(6))
          vx(6,i)=atn*(vv(5)-3*vv(6))
          vx(7,i)=als*(vv(7)+vv(8))
          vx(8,i)=als*(vv(7)-3*vv(8))
          vx(9,i)=1.0*vv(9)+ast*(vv(10)-3*(vv(11)+vv(12)))
          vx(10,i)=1.0*vv(9)-3*ast*(vv(10)+vv(11)-3.*vv(12))
          vx(11,i)=1.0*vv(9)+ast*(vv(10)+vv(11)+vv(12))
          vx(12,i)=1.0*vv(9)+ast*(vv(11)-3*(vv(10)+vv(12)))
          vx(13,i)=1.0*(vv(13)+vv(14))
          vx(14,i)=1.0*(vv(13)-3*vv(14))
        else if (nm.eq.2) then
          vx(1,i)=acn*(vv(1)+vv(2)+2*(vv(15)-vv(18)))
     &         -3*ast*(vv(3)+vv(4)+2*vv(16))
          vx(3,i)=acn*(vv(1)+vv(2)+2*(vv(15)-vv(18)))
     &           +ast*(vv(3)+vv(4)+2*vv(16))
          vx(5,i)=atn*(vv(5)+vv(6)+2*vv(17))
          vx(7,i)=als*(vv(7)+vv(8))
          vx(9,i)=1.0*(vv(9)+vv(10))-3*ast*(vv(11)+vv(12))
          vx(11,i)=1.0*(vv(9)+vv(10))+ast*(vv(11)+vv(12))
          vx(13,i)=1.0*(vv(13)+vv(14))
	  vx(2,i)=0.
	  vx(4,i)=0.
	  vx(6,i)=0.
	  vx(8,i)=0.
	  vx(10,i)=0.
	  vx(12,i)=0.
	  vx(14,i)=0.
        end if
   50 continue

c     rr=0.
c     call pot(1,rr,vv,ww)
c     if (nm.eq.1) then
c       rv(1)=acn*vv(1)+ast*(vv(2)-3*(vv(3)+vv(4)))
c       rv(4)=acn*vv(1)+ast*(vv(3)-3*(vv(2)+vv(4)))
c       rv(6)=atn*(vv(5)-3*vv(6))
c     else if (nm.eq.2) then
c       rv(1)=acn*(vv(1)+vv(2)+2*(vv(15)-vv(18)))
c    &     -3*ast*(vv(3)+vv(4)+2*vv(16))
c     end if



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  solve for psi=r*f*j in uncoupled channels
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       do ik=1,kgrid
       do ir=1,ngrid
       do ip=1,4
         psis(ip,ir,ik)=jlx(ip,ir,ik)*rx(ir)
        if(ip.ne.1)then
         psit(ip,ir,ik)=jlx(ip,ir,ik)*rx(ir)
        else
         psit(ip,ir,ik)=jlx(2,ir,ik)*rx(ir)
        end if
       end do
	if(nm.eq.2) then
	 psis(2,ir,ik)=0.
	 psis(4,ir,ik)=0.
	 psit(3,ir,ik)=0.
        end if
       end do
       end do


C jp: j+1. jp:odd: tau=1 (ivs=1), jp:even: tau=0 (ivs=2) (for s=0)
C          jp:odd: tau=0 (ivt=2), jp:even: tau=1 (ivt=1) (for s=1)

      do 500 jp=1,4

      if(mod(jp,2).eq.1)then
      ivs=1
      ivt=2
      else
      ivs=2
      ivt=1
      end if

      do 400 ik=1,kgrid

      clm(1,jp,ik)=1.0
      clm(2,jp,ik)=1.0
      small=1.e-7
      tiny=1.e-10

      if(nm.eq.2.and.mod(jp,2).eq.0) goto 125
C     singlet s loop

      algd=jlpx(jp,lcx,ik)/jlx(jp,lcx,ik) +1./rx(lcx)
      a=(jlx(jp,lcxp,ik)-jlx(jp,lcxm,ik))/(2.*h)
      b=jlx(jp,lcx,ik)-a*rx(lcx)
      l0=-b/(a*h)
      l0m=l0-1
      l0p=l0+1
      rdelt=-b/a -h*float(l0)

      do 110 loop1=1,20

      psim=0.
      gpsim=0.
      psi0=h/10.
      gpsi0=h2*((h2m*(jp-1)*jp)/rx(1)**2 -h2m*k(ik)**2 +vx(ivs,1)
     &   +jp*(jp-1)*vx(ivs+8,1)-clm(1,jp,ik))*psi0

C     Integrate out PSIS

      psis(jp,1,ik)=psi0

      if(abs(algd).lt.10) then
       do 100 ir=2,lcxp

       gp=h2*((h2m*(jp-1)*jp)/rx(ir)**2 -h2m*k(ik)**2 +vx(ivs,ir)
     &    +jp*(jp-1)* vx(ivs+8,ir)-clm(1,jp,ik))

       psis(jp,ir,ik)= (10.*gpsi0 +2.*psi0 +gpsim -psim)/(1.-gp)

       psim=psi0
       gpsim=gpsi0
       psi0=psis(jp,ir,ik)
       gpsi0=gp*psi0

100    end do
      else
       do 101 ir=2,l0p

       gp=h2*((h2m*(jp-1)*jp)/rx(ir)**2 -h2m*k(ik)**2 +vx(ivs,ir)
     &    +jp*(jp-1)* vx(ivs+8,ir)-clm(1,jp,ik))

       psis(jp,ir,ik)= (10.*gpsi0 +2.*psi0 +gpsim -psim)/(1.-gp)

       psim=psi0
       gpsim=gpsi0
       psi0=psis(jp,ir,ik)
       gpsi0=gp*psi0

101    end do
      end if

C     Match the logarithmic derivative to bessel functions
C     if at a node, match the value

      if(abs(algd).lt.10.0) then
         slope= (psis(jp,lcxp,ik)-psis(jp,lcxm,ik))/psis(jp,lcx,ik)
     &  -2.*h*(jlpx(jp,lcx,ik)/jlx(jp,lcx,ik)+ 1./rx(lcx))
      else
         slope= psis(jp,l0,ik)+(psis(jp,l0p,ik)-psis(jp,l0,ik))*rdelt/h
      end if

      if(loop1.eq.1) then
	slopeo=slope
	clmo=clm(1,jp,ik)
	clm(1,jp,ik)=0.9*clmo

	if(ik.gt.1) then
	 clmo=clm(1,jp,ik-1)
	 clm(1,jp,ik)=0.9*clmo
        end if

      else
c       write(7,3700) jp,ik,loop1,clm(1,jp,ik),abs(slopeo-slope)
	if(abs(slopeo-slope).le.small) goto 120
	clmn=(slopeo*clm(1,jp,ik)-slope*clmo)/(slopeo-slope)
	slopeo=slope
	clmo=clm(1,jp,ik)
	clm(1,jp,ik)=clmn
      end if

110   continue

120   if(abs(algd).lt.10.)then
         fac=rx(lcx)*jlx(jp,lcx,ik)/psis(jp,lcx,ik)
         do ir=1,lcx
	  psis(jp,ir,ik)= fac*psis(jp,ir,ik)
         end do
         psis(jp,lcxp,ik)=rx(lcxp)*jlx(jp,lcxp,ik)
       else
	 fac=(rx(l0p)*jlx(jp,l0p,ik)-rx(l0)*jlx(jp,l0,ik))/
     &       (psis(jp,l0p,ik)-psis(jp,l0,ik))
	 do ir=1,l0p
	  psis(jp,ir,ik)=fac*psis(jp,ir,ik)
         end do
       end if


125   if(jp.eq.1) then

C     3P0 loop

      algd=jlpx(jp+1,ltx,ik)/jlx(jp+1,ltx,ik) + 1./rx(ltx)
      a=(jlx(jp+1,ltxp,ik)-jlx(jp+1,ltxm,ik))/(2.*h)
      b=jlx(jp+1,ltx,ik)-a*rx(ltx)
      l0=-b/(a*h)
      l0m=l0-1
      l0p=l0+1
      rdelt=-b/a -h*float(l0)


      do 210 loop2=1,20

      psim=0.
      gpsim=0.
      psi0=h/10.
      gpsi0=h2*(h2m*2./rx(1)**2 -h2m*k(ik)**2 +vx(3,1)
     &+2* vx(11,1) -4.*vx(5,1) -2.*vx(7,1)
     &+4.*vx(13,1) -clm(2,1,ik))*psi0

C     Integrate out PSIT

      psit(jp,1,ik)=psi0

      if(abs(algd).lt.10.0)then
       do 200 ir=2,ltxp

       gp=h2*(h2m*2./rx(ir)**2 -h2m*k(ik)**2 +vx(3,ir)
     & +2.* vx(11,ir) -4.*vx(5,ir) -2.*vx(7,ir)
     & +4.*vx(13,ir) -clm(2,1,ik))

       psit(1,ir,ik)= (10.*gpsi0 +2.*psi0 +gpsim -psim)/(1.-gp)

       psim=psi0
       gpsim=gpsi0
       psi0=psit(jp,ir,ik)
       gpsi0=gp*psi0

200    end do
      else
       do 201 ir=2,l0p

       gp=h2*(h2m*2./rx(ir)**2 -h2m*k(ik)**2 +vx(3,ir)
     & +2.* vx(11,ir) -4.*vx(5,ir) -2.*vx(7,ir)
     & +4.*vx(13,ir) -clm(2,1,ik))

       psit(1,ir,ik)= (10.*gpsi0 +2.*psi0 +gpsim -psim)/(1.-gp)

       psim=psi0
       gpsim=gpsi0
       psi0=psit(jp,ir,ik)
       gpsi0=gp*psi0

201    end do
      end if



C     Match the logarithmic derivative to bessel functions
C     If at a node, match the value

      if(abs(algd).lt.10.0) then
         slope= (psit(jp,ltxp,ik)-psit(jp,ltxm,ik))/psit(jp,ltx,ik)
     &  -2.*h*(jlpx(jp+1,ltx,ik)/jlx(jp+1,ltx,ik)+ 1./rx(ltx))
      else
	 slope= psit(jp,l0,ik)+(psit(jp,l0p,ik)-psit(jp,l0,ik))*rdelt/h
      end if

      if(loop2.eq.1) then
	slopeo=slope
	clmo=clm(2,jp,ik)
	clm(2,jp,ik)=0.9*clmo

        if(ik.gt.1) then
         clmo=clm(2,jp,ik-1)
         clm(2,jp,ik)=0.9*clmo
        end if

      else
c       write(8,3700) jp,ik,loop2,clm(2,jp,ik),abs(slopeo-slope)
	if(abs(slopeo-slope).le.small) goto 220
	clmn=(slopeo*clm(2,jp,ik)-slope*clmo)/(slopeo-slope)
	slopeo=slope
	clmo=clm(2,jp,ik)
	clm(2,jp,ik)=clmn
      end if

210   continue

220    if(abs(algd).lt.10.)then
         fac=rx(ltx)*jlx(jp+1,ltx,ik)/psit(jp,ltx,ik)
         do ir=1,ltx
           psit(jp,ir,ik)=fac*psit(jp,ir,ik)
         end do
         psit(jp,ltxp,ik)=rx(ltxp)*jlx(jp+1,ltxp,ik)
      else
	 fac=(rx(l0p)*jlx(jp+1,l0p,ik)-rx(l0)*jlx(jp+1,l0,ik))/
     &       (psit(jp,l0p,ik)-psit(jp,l0,ik))
         do ir=1,l0p
           psit(jp,ir,ik)=fac*psit(jp,ir,ik)
         end do
      end if


      goto 400
      end if

      if(nm.eq.2.and.jp.eq.3) goto 400
C     uncoupled triplet s loop

      algd=jlpx(jp,ltx,ik)/jlx(jp,ltx,ik) + 1./rx(ltx)
      a=(jlx(jp,ltxp,ik)-jlx(jp,ltxm,ik))/(2.*h)
      b=jlx(jp,ltx,ik)-a*rx(ltx)
      l0=-b/(a*h)
      l0m=l0-1
      l0p=l0+1
      rdelt=-b/a -h*float(l0)

      do 310 loop3=1,20

      psim=0.
      gpsim=0.
      psi0=h/10.
      gpsi0=h2*((h2m*(jp-1)*jp)/rx(1)**2 -h2m*k(ik)**2 +vx(ivt+2,1)
     &+jp*(jp-1)* vx(ivt+10,1) +2.*vx(ivt+4,1) -vx(ivt+6,1)
     &+vx(ivt+12,1) -clm(2,jp,ik))*psi0

C     Integrate out PSIT

      psit(jp,1,ik)=psi0

      if(abs(algd).lt.10) then
       do 300 ir=2,ltxp

       gp=h2*((h2m*(jp-1)*jp)/rx(ir)**2 -h2m*k(ik)**2 +vx(ivt+2,ir)
     & +jp*(jp-1)* vx(ivt+10,ir) +2.*vx(ivt+4,ir) -vx(ivt+6,ir)
     & +vx(ivt+12,ir) -clm(2,jp,ik))

       psit(jp,ir,ik)= (10.*gpsi0 +2.*psi0 +gpsim -psim)/(1.-gp)

       psim=psi0
       gpsim=gpsi0
       psi0=psit(jp,ir,ik)
       gpsi0=gp*psi0

300    end do
      else
       do 301 ir=2,l0p

       gp=h2*((h2m*(jp-1)*jp)/rx(ir)**2 -h2m*k(ik)**2 +vx(ivt+2,ir)
     & +jp*(jp-1)* vx(ivt+10,ir) +2.*vx(ivt+4,ir) -vx(ivt+6,ir)
     & +vx(ivt+12,ir) -clm(2,jp,ik))

       psit(jp,ir,ik)= (10.*gpsi0 +2.*psi0 +gpsim -psim)/(1.-gp)

       psim=psi0
       gpsim=gpsi0
       psi0=psit(jp,ir,ik)
       gpsi0=gp*psi0

301    end do
      end if

C     Match the logarithmic derivative to bessel functions
C     if at a node, match the value


      if(abs(algd).lt.10.0) then
         slope= (psit(jp,ltxp,ik)-psit(jp,ltxm,ik))/psit(jp,ltx,ik)
     &  -2.*h*(jlpx(jp,ltx,ik)/jlx(jp,ltx,ik)+ 1./rx(ltx))
      else
         slope= psit(jp,l0,ik)+(psit(jp,l0p,ik)-psit(jp,l0,ik))*rdelt/h
      end if

      if(loop3.eq.1) then
	slopeo=slope
	clmo=clm(2,jp,ik)
	clm(2,jp,ik)=0.9*clmo

        if(ik.gt.1) then
         clmo=clm(2,jp,ik-1)
         clm(2,jp,ik)=0.9*clmo
        end if

      else
c       write(9,3700) jp,ik,loop3,clm(2,jp,ik),abs(slopeo-slope)
	if(abs(slopeo-slope).le.small) goto 320
	clmn=(slopeo*clm(2,jp,ik)-slope*clmo)/(slopeo-slope)
	slopeo=slope
	clmo=clm(2,jp,ik)
	clm(2,jp,ik)=clmn
      end if

310   continue

320   if(abs(algd).lt.10.)then
         fac=rx(ltx)*jlx(jp,ltx,ik)/psit(jp,ltx,ik)
	 do ir=1,ltx
	  psit(jp,ir,ik)=fac*psit(jp,ir,ik)
         end do
	 psit(jp,ltxp,ik)=rx(ltxp)*jlx(jp,ltxp,ik)
      else
         fac=(rx(l0p)*jlx(jp,l0p,ik)-rx(l0)*jlx(jp,l0,ik))/
     &       (psit(jp,l0p,ik)-psit(jp,l0,ik))
	 do ir=1,l0p
	  psit(jp,ir,ik)=fac*psit(jp,ir,ik)
         end do
      end if


400   continue
500   continue


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c solve for correlations in 3S1-3D1,3P2-3F2 and 3D3-3G3 coupled channels
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      psitm(1:3,1:ngrid,1:kgrid,1:2)=0.
      psitp(1:3,1:ngrid,1:kgrid,1:2)=0.

      xu(1:ngrid)=0.
      xw(1:ngrid)=0.
      xh(1:ngrid)=0.

      if(nm.eq.2) goto 5201
      goto 5201
      do 5001 i=1,ngrid
      xvs(i)=vx(4,i)
      xvd(i)=xvs(i)-2.*vx(6,i)-3.*vx(8,i)+6.*vx(12,i)
     &      +9.*vx(14,i)+h2m*6./rsx(i)
      xvt(i)=vx(6,i)
5001  continue
      h2r=h2*sqrt(8.)

      do 5002 ik=1,10
      xek=h2m*k(ik)*k(ik)
c     xek=0.
      if(ik.eq.1)then
      xld=-10.
      xlt=-0.1
      end if

      xu(1:ngrid)=0.
      do i=1,40
      xw(i)=rsx(i)*rx(i)
      end do

      do 5200 id=1,10

      do i=1,40
      xh(i)=xw(i)
      end do

      do 5103 iw=1,20
      if(iw.eq.20) then
      write(nlog,*) "iw over"," id=",id," 3D1"
      end if
      xwm=xw(39)
      xw0=xw(40)
      gm=h2*(xvd(39)-xek-xld)*xwm
      gp=h2*(xvd(40)-xek-xld)
      g0=gp*xw0
      qm=h2r*(xvt(39)-xlt)*xu(39)
      q0=h2r*(xvt(40)-xlt)*xu(40)
      xhm=xh(39)
      xh0=xh(40)
      ghm=h2*(xvd(39)-xek-xld)*xhm
      gh0=gp*xh0


      do 5104  ir=41,ltxp
      gp=h2*(xvd(ir)-xek-xld)
      qp=h2r*(xvt(ir)-xlt)*xu(ir)
      xw(ir)=(10.*g0+gm+2.*xw0-xwm+10.*q0+qm+qp)/(1.-gp)
      xwm=xw0
      xw0=xw(ir)
      gm=g0
      g0=gp*xw0
      xh(ir)=(10.*gh0+ghm+2.*xh0-xhm)/(1.-gp)
      xhm=xh0
      xh0=xh(ir)
      ghm=gh0
      gh0=gp*xh0
      qm=q0
      q0=qp
5104  continue

c     print*,ik,id,iw,xld
      fac=(rx(ltx)*jlx(3,ltx,ik)-xw(ltx))/xh(ltx)
      diff=xw(ltxp)-xw(ltxm)+fac*(xh(ltxp)-xh(ltxm))
     &    -rx(ltxp)*jlx(3,ltxp,ik)+rx(ltxm)*jlx(3,ltxm,ik)

      if(iw.eq.1)then
      odiff=diff
      oxld=xld
      xld=xld*.9
      goto 5105
      else
      if(abs(odiff-diff).le.small)goto 5106
      xldn=(odiff*xld-diff*oxld)/(odiff-diff)
      oxld=xld
      odiff=diff
      xld=xldn
      end if
5105  continue
5103  continue
5106  continue

      do 5109 ir=1,ltxp
5109  xw(ir)=xw(ir)+fac*xh(ir)

      xu(1)=xw(1)
      xh(1)=xu(1)
      do 5003 iu=1,80
      if(iu.eq.20) then
      write(nlog,*)"iu over"," id=",id," 3D1"
      end if
      xum=0
      xu0=xu(1)
      gm=0.
      gp=h2*(xvs(1)-xek-xld)
      g0=gp*xu0
      qm=0.
      q0=h2r*(xvt(1)-xlt)*xw(1)
      xhm=0.
      xh0=xh(1)
c     xh0=h/10.
      ghm=0.
      gh0=gp*xh0
c     xh(1)=xh0

      do 5004 ir=2,lcxp
      gp=h2*(xvs(ir)-xek-xld)
      qp=h2r*(xvt(ir)-xlt)*xw(ir)
      xu(ir)=(10.*g0+gm+2.*xu0-xum+10.*q0+qm+qp)/(1.-gp)
      xum=xu0
      xu0=xu(ir)
      gm=g0
      g0=gp*xu0
      xh(ir)=(10.*gh0+ghm+2.*xh0-xhm)/(1.-gp)
      xhm=xh0
      xh0=xh(ir)
      ghm=gh0
      gh0=gp*xh0
      qm=q0
      q0=qp
5004  continue

c     print*,ik,id,iu,xlt
      fac=-xu(lcx)/xh(lcx)
      diff=xu(lcxp)-xu(lcxm)+fac*(xh(lcxp)-xh(lcxm))
c     print*,ik,id,iu,xlt,diff
      if(iu.eq.1) then
      odiff=diff
      oxlt=xlt
      xlt=xlt*.9
      goto 5005
      else
      if(abs(odiff-diff).le.tiny) goto 5006
      xltn=(odiff*xlt-diff*oxlt)/(odiff-diff)
      oxlt=xlt
      odiff=diff
      xlt=xltn
      end if
5005  continue
5003  continue
5006  continue


      do 5009 ir=1,lcxp
      xu(ir)=xu(ir)+fac*xh(ir)
5009  continue

5200  continue

      do ir=1,ltxp
      psitm(1,ir,ik,2)=xu(ir)
      psitp(1,ir,ik,2)=xw(ir)
      end do

5002  continue
5201  continue


      do 6001 i=1,ngrid
      xvp(i)=vx(3,i)-(2./5.)*vx(5,i)+vx(7,i)
     &      +2.*vx(11,i)+vx(13,i)+h2m*2./rsx(i)
      xvf(i)=vx(3,i)-(8./5.)*vx(5,i)-4.*vx(7,i)
     &      +12.*vx(11,i)+16.*vx(13,i)+h2m*12./rsx(i)
      xvt(i)=vx(5,i)
6001  continue

      h2r=h2*sqrt(6.)*6./5.
      xlp=-10.
      xlt=-0.1
      xw(1:ngrid)=0.
      do i=1,20
      xu(i)=rsx(i)
      end do

      do 6002 ik=1,kgrid
      xek=h2m*k(ik)*k(ik)
c     xek=0.
      pnode=4.49340945791/k(ik)
      lnode=pnode/h
      lmatch=ltx
      if(lnode.lt.lmatch) then
      lmatch=lnode
      end if
      lmatchp=lmatch+1
      lmatchm=lmatch-1

      do 6200 id=1,10

      do 6003 iu=1,20
       do i=1,20
       xh(i)=xu(i)
       end do
      if(iu.eq.20) then
      write(nlog,*)"iu over"," id=",id
      end if
      xum=xu(19)
      xu0=xu(20)
      gm=h2*(xvp(19)-xek-xlp)*xum
      gp=h2*(xvp(20)-xek-xlp)
      g0=gp*xu0
      qm=h2r*(xvt(19)-xlt)*xw(19)
      q0=h2r*(xvt(20)-xlt)*xw(20)
      xhm=xh(19)
      xh0=xh(20)
      ghm=h2*(xvp(19)-xek-xlp)*xhm
      gh0=gp*xh0

      do 6004 ir=21,ltxp
      gp=h2*(xvp(ir)-xek-xlp)
      qp=h2r*(xvt(ir)-xlt)*xw(ir)
      xu(ir)=(10.*g0+gm+2.*xu0-xum+10.*q0+qm+qp)/(1.-gp)
      xum=xu0
      xu0=xu(ir)
      gm=g0
      g0=gp*xu0
      xh(ir)=(10.*gh0+ghm+2.*xh0-xhm)/(1.-gp)
      xhm=xh0
      xh0=xh(ir)
      ghm=gh0
      gh0=gp*xh0
      qm=q0
      q0=qp
6004  continue

c     print*,ik,id,iu,xlp
      fac=(rx(ltx)*jlx(2,ltx,ik)-xu(ltx))/xh(ltx)
      diff=xu(ltxp)-xu(ltxm)+fac*(xh(ltxp)-xh(ltxm))
     &    -rx(ltxp)*jlx(2,ltxp,ik)+rx(ltxm)*jlx(2,ltxm,ik)
      if(iu.eq.1) then
      odiff=diff
      oxlp=xlp
      xlp=xlp*.9
      goto 6005
      else
      if(abs(odiff-diff).le.small) goto 6006
      xlpn=(odiff*xlp-diff*oxlp)/(odiff-diff)
      oxlp=xlp
      odiff=diff
      xlp=xlpn
      end if
6005  continue
6003  continue
6006  continue


      do 6009 ir=1,ltxp
      xu(ir)=xu(ir)+fac*xh(ir)
6009  continue

      if(ik.eq.1.and.id.eq.1)then
      do i=1,20
      xw(i)=rsx(i)*rsx(i)
      end do
      end if

      do i=1,20
      xh(i)=xw(i)
      end do

      do 6103 iw=1,20
      if(iw.eq.20) then
      write(nlog,*) "iw over"," id=",id
      end if
      xwm=xw(19)
      xw0=xw(20)
      gm=h2*(xvf(19)-xek-xlp)*xwm
      gp=h2*(xvf(20)-xek-xlp)
      g0=gp*xw0
      qm=h2r*(xvt(19)-xlt)*xu(19)
      q0=h2r*(xvt(20)-xlt)*xu(20)
      xhm=xh(19)
      xh0=xh(20)
      ghm=h2*(xvf(19)-xek-xlp)*xhm
      gh0=gp*xh0


      do 6104  ir=21,lmatchp
      gp=h2*(xvf(ir)-xek-xlp)
      qp=h2r*(xvt(ir)-xlt)*xu(ir)
      xw(ir)=(10.*g0+gm+2.*xw0-xwm+10.*q0+qm+qp)/(1.-gp)
      xwm=xw0
      xw0=xw(ir)
      gm=g0
      g0=gp*xw0
      xh(ir)=(10.*gh0+ghm+2.*xh0-xhm)/(1.-gp)
      xhm=xh0
      xh0=xh(ir)
      ghm=gh0
      gh0=gp*xh0
      qm=q0
      q0=qp
6104  continue

c     print*,ik,id,iw,xlt
      fac=-xw(lmatch)/xh(lmatch)
      diff=xw(lmatchp)-xw(lmatchm)+fac*(xh(lmatchp)-xh(lmatchm))
      if(iw.eq.1)then
      odiff=diff
      oxlt=xlt
      xlt=xlt*.9
      goto 6105
      else
      if(abs(odiff-diff).le.small)goto 6106
      xltn=(odiff*xlt-diff*oxlt)/(odiff-diff)
      oxlt=xlt
      odiff=diff
      xlt=xltn
      end if
6105  continue
6103  continue
6106  continue

      do 6109 ir=1,lmatchp
6109  xw(ir)=xw(ir)+fac*xh(ir)
      if(lmatch.ne.ltx)then
      xw(lmatchp+1:ltxp)=0.
      end if

6200  continue

      do ir=1,ltxp
      psitm(2,ir,ik,1)=xu(ir)
      psitp(2,ir,ik,1)=xw(ir)
      end do

c     print*,ik,xlp,xlt
      if(ik.eq.2)then
       xlpm=xlp
       xltm=xlt
      elseif(ik.ge.3)then
       xlpmm=xlpm
       xlpm=xlp
       xltmm=xltm
       xltm=xlt
       xlp=2.*xlpm-xlpmm
       xlt=2.*xltm-xltmm
      end if
6002  continue

      if(nm.eq.2) goto 7201
      do 7001 i=1,ngrid
      xvd(i)=vx(4,i)-(4./7.)*vx(6,i)+2.*vx(8,i)
     &      +6.*vx(12,i)+4.*vx(14,i)+h2m*6./rsx(i)
      xvg(i)=vx(4,i)-(10./7.)*vx(6,i)-5.*vx(8,i)
     &      +20.*vx(12,i)+25.*vx(14,i)+h2m*20./rsx(i)
      xvt(i)=vx(6,i)
7001  continue

      h2r=h2*sqrt(12.)*6./7.
      xld=-10.
      xlt=-0.1
      xw(1:ngrid)=0.
      do i=1,40
      xu(i)=rx(i)*rx(i)*rx(i)
      end do

      do 7002 ik=1,kgrid
      xek=h2m*k(ik)*k(ik)
c     xek=0.
      dnode=5.76345919689/k(ik)
      lnode=dnode/h
      lmatch=ltx
      if(lnode.lt.lmatch) then
      lmatch=lnode
      end if
      lmatchp=lmatch+1
      lmatchm=lmatch-1

      do 7200 id=1,10

      do 7003 iu=1,20
       do i=1,40
       xh(i)=xu(i)
       end do
      if(iu.eq.20) then
      write(nlog,*)"iu over"," id=",id
      end if
      xum=xu(39)
      xu0=xu(40)
      gm=h2*(xvd(39)-xek-xld)*xum
      gp=h2*(xvd(40)-xek-xld)
      g0=gp*xu0
      qm=h2r*(xvt(39)-xlt)*xw(39)
      q0=h2r*(xvt(40)-xlt)*xw(40)
      xhm=xh(39)
      xh0=xh(40)
      ghm=h2*(xvd(39)-xek-xld)*xhm
      gh0=gp*xh0

      do 7004 ir=41,ltxp
      gp=h2*(xvd(ir)-xek-xld)
      qp=h2r*(xvt(ir)-xlt)*xw(ir)
      xu(ir)=(10.*g0+gm+2.*xu0-xum+10.*q0+qm+qp)/(1.-gp)
      xum=xu0
      xu0=xu(ir)
      gm=g0
      g0=gp*xu0
      xh(ir)=(10.*gh0+ghm+2.*xh0-xhm)/(1.-gp)
      xhm=xh0
      xh0=xh(ir)
      ghm=gh0
      gh0=gp*xh0
      qm=q0
      q0=qp
7004  continue

c     print*,ik,id,iu,xld
      fac=(rx(ltx)*jlx(3,ltx,ik)-xu(ltx))/xh(ltx)
      diff=xu(ltxp)-xu(ltxm)+fac*(xh(ltxp)-xh(ltxm))
     &    -rx(ltxp)*jlx(3,ltxp,ik)+rx(ltxm)*jlx(3,ltxm,ik)
      if(iu.eq.1) then
      odiff=diff
      oxld=xld
      xld=xld*.9
      goto 7005
      else
      if(abs(odiff-diff).le.small) goto 7006
      xldn=(odiff*xld-diff*oxld)/(odiff-diff)
      oxld=xld
      odiff=diff
      xld=xldn
      end if
7005  continue
7003  continue
7006  continue


      do 7009 ir=1,ltxp
      xu(ir)=xu(ir)+fac*xh(ir)
7009  continue


      if(ik.eq.1.and.id.eq.1)then
      do i=1,40
      xw(i)=rsx(i)*rsx(i)*rx(i)
      end do
      end if

      do i=1,40
      xh(i)=xw(i)
      end do

      do 7103 iw=1,20
      if(iw.eq.20) then
      write(nlog,*) "iw over"," id=",id
      end if
      xwm=xw(39)
      xw0=xw(40)
      gm=h2*(xvg(39)-xek-xld)*xwm
      gp=h2*(xvg(40)-xek-xld)
      g0=gp*xw0
      qm=h2r*(xvt(39)-xlt)*xu(39)
      q0=h2r*(xvt(40)-xlt)*xu(40)
      xhm=xh(39)
      xh0=xh(40)
      ghm=h2*(xvg(39)-xek-xld)*xhm
      gh0=gp*xh0


      do 7104  ir=41,lmatchp
      gp=h2*(xvg(ir)-xek-xld)
      qp=h2r*(xvt(ir)-xlt)*xu(ir)
      xw(ir)=(10.*g0+gm+2.*xw0-xwm+10.*q0+qm+qp)/(1.-gp)
      xwm=xw0
      xw0=xw(ir)
      gm=g0
      g0=gp*xw0
      xh(ir)=(10.*gh0+ghm+2.*xh0-xhm)/(1.-gp)
      xhm=xh0
      xh0=xh(ir)
      ghm=gh0
      gh0=gp*xh0
      qm=q0
      q0=qp
7104  continue

c     print*,ik,id,iw,xlt
      fac=-xw(lmatch)/xh(lmatch)
      diff=xw(lmatchp)-xw(lmatchm)+fac*(xh(lmatchp)-xh(lmatchm))
      if(iw.eq.1)then
      odiff=diff
      oxlt=xlt
      xlt=xlt*.9
      goto 7105
      else
      if(abs(odiff-diff).le.small)goto 7106
      xltn=(odiff*xlt-diff*oxlt)/(odiff-diff)
      oxlt=xlt
      odiff=diff
      xlt=xltn
      end if
7105  continue
7103  continue
7106  continue

      do 7109 ir=1,lmatchp
7109  xw(ir)=xw(ir)+fac*xh(ir)
      if(lmatch.ne.ltx)then
      xw(lmatchp+1:ltxp)=0.
      end if

7200  continue

      do ir=1,ltxp
      psitm(3,ir,ik,1)=xu(ir)
      psitp(3,ir,ik,1)=xw(ir)
      end do

c     print*,ik,xld,xlt
      if(ik.eq.2)then
       xldm=xld
       xltm=xlt
      elseif(ik.ge.3)then
       xldmm=xldm
       xldm=xld
       xltmm=xltm
       xltm=xlt
       xld=2.*xldm-xldmm
       xlt=2.*xltm-xltmm
      end if
7002  continue
7201  continue



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 2-body cluster contributions for k-averaged and
c unaveraged correlations
c uk(1-8): 1S0,1P1,1D2,1F3,3P0,3P1,3D2,3F3 unaveraged
c unm(1-15): 1S0,1P1,1D2,1F3,3P0,3P1,3D2,3F3,1G4,1H5,1I6,1J7,
c            3G4,3H5,3I6 k-averaged
c utk(1-3;1-2): 3S1(u-,u+),3P2(u-,u+),3D3(u-,u+) unaveraged
c wtk(1-3;1-2): 3D1(w-,w+),3F2(w-,w+),3G3(w-,w+) unaveraged
c utnm(1-5;1,2): 3S1(u-,u+),3P2(u-,u+),3D3(u-,u+),
c                3F4(u-,u+),3G5(u-,u+)            k-averaged
c wtnm(1-5;1,2): 3D1(w-,w+),3F2(w-,w+),3G3(w-,w+),
c                3H4(w-,w+),3I5(w-,w+)            k-averaged
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c off-diagonal matrix element factors:
      a1=2.*rt2
      a2=6.*sqrt(6.)/5.
      a3=12.*sqrt(3.)/7.
      a4=4.*rt5/3.
      a5=6.*sqrt(30.)/11.


      do ik=1,kgrid
      do ir=1,ngrid

       if(ir.le.ltx) then
	if(nm.eq.1) then
        fc10=psi(1,ir)/phir(1,ir)
        fc00=psi(2,ir)/phir(2,ir)
        fc11=psi(3,ir)/phir(3,ir)
        fc01=psi(4,ir)/phir(4,ir)
        ft11=psi(5,ir)/phir(5,ir)
        ft01=psi(6,ir)/phir(6,ir)
        fb11=psi(7,ir)/phir(7,ir)
        fb01=psi(8,ir)/phir(8,ir)
	elseif(nm.eq.2) then
	fc10=psi(1,ir)/phir(1,ir)
	fc00=0.
	fc11=psi(3,ir)/phir(3,ir)
	fc01=0.
	ft11=psi(5,ir)/phir(5,ir)
	ft01=0.
	fb11=psi(7,ir)/phir(7,ir)
	fb01=0.
	end if
       end if

        do ip=1,4
         uk(ip,ir,ik)=psis(ip,ir,ik)
         uk(ip+4,ir,ik)=psit(ip,ir,ik)
        end do

	do ip=1,3
	do ix=1,2
	 utk(ip,ix,ir,ik)=psitm(ip,ir,ik,ix)
	 wtk(ip,ix,ir,ik)=psitp(ip,ir,ik,ix)
        end do
	end do

        if(ir.le.lcx) then
         unm(1,ir,ik)=rx(ir)*jlx(1,ir,ik)*fc10
         unm(2,ir,ik)=rx(ir)*jlx(2,ir,ik)*fc00
         unm(3,ir,ik)=rx(ir)*jlx(3,ir,ik)*fc10
         unm(4,ir,ik)=rx(ir)*jlx(4,ir,ik)*fc00
         unm(9,ir,ik)=rx(ir)*jlx(5,ir,ik)*fc10
         unm(10,ir,ik)=rx(ir)*jlx(6,ir,ik)*fc00
	 unm(11,ir,ik)=rx(ir)*jlx(7,ir,ik)*fc10
	 unm(12,ir,ik)=rx(ir)*jlx(8,ir,ik)*fc00
         utnm(1,1,ir,ik)=rx(ir)*jlx(1,ir,ik)*fc01
        else
         unm(1,ir,ik)=rx(ir)*jlx(1,ir,ik)
         unm(2,ir,ik)=rx(ir)*jlx(2,ir,ik)
         unm(3,ir,ik)=rx(ir)*jlx(3,ir,ik)
         unm(4,ir,ik)=rx(ir)*jlx(4,ir,ik)
         unm(9,ir,ik)=rx(ir)*jlx(5,ir,ik)
         unm(10,ir,ik)=rx(ir)*jlx(6,ir,ik)
	 unm(11,ir,ik)=rx(ir)*jlx(7,ir,ik)
	 unm(12,ir,ik)=rx(ir)*jlx(8,ir,ik)
         utnm(1,1,ir,ik)=rx(ir)*jlx(1,ir,ik)
	 if(nm.eq.2) then
	 unm(2,ir,ik)=0.
	 unm(4,ir,ik)=0.
	 unm(10,ir,ik)=0.
	 unm(12,ir,ik)=0.
	 utnm(1,1,ir,ik)=0.
	 end if
        end if

        if(ir.le.ltx) then
         unm(5,ir,ik)=rx(ir)*jlx(2,ir,ik)*(fc11-4.*ft11-2.*fb11)
         unm(6,ir,ik)=rx(ir)*jlx(2,ir,ik)*(fc11+2.*ft11 -  fb11)
         unm(7,ir,ik)=rx(ir)*jlx(3,ir,ik)*(fc01+2.*ft01 -  fb01)
         unm(8,ir,ik)=rx(ir)*jlx(4,ir,ik)*(fc11+2.*ft11 -  fb11)
	 unm(13,ir,ik)=rx(ir)*jlx(5,ir,ik)*(fc01+2.*ft01 -  fb01)
	 unm(14,ir,ik)=rx(ir)*jlx(6,ir,ik)*(fc11+2.*ft11 -  fb11)
	 unm(15,ir,ik)=rx(ir)*jlx(7,ir,ik)*(fc01+2.*ft01 -  fb01)
         utnm(1,2,ir,ik)=a1*rx(ir)*jlx(3,ir,ik)*ft01
         utnm(2,2,ir,ik)=a2*rx(ir)*jlx(4,ir,ik)*ft11
         utnm(3,2,ir,ik)=a3*rx(ir)*jlx(5,ir,ik)*ft01
	 utnm(4,2,ir,ik)=a4*rx(ir)*jlx(6,ir,ik)*ft11
	 utnm(5,2,ir,ik)=a5*rx(ir)*jlx(7,ir,ik)*ft01
       utnm(2,1,ir,ik)=rx(ir)*jlx(2,ir,ik)*(fc11-(2./5.)  *ft11+   fb11)
       utnm(3,1,ir,ik)=rx(ir)*jlx(3,ir,ik)*(fc01-(4./7.)  *ft01+2.*fb01)
       utnm(4,1,ir,ik)=rx(ir)*jlx(4,ir,ik)*(fc11-(2./3.)  *ft11+3.*fb11)
       utnm(5,1,ir,ik)=rx(ir)*jlx(5,ir,ik)*(fc01-(8./11.) *ft01+4.*fb01)
       wtnm(1,2,ir,ik)=rx(ir)*jlx(3,ir,ik)*(fc01- 2.      *ft01-3.*fb01)
       wtnm(2,2,ir,ik)=rx(ir)*jlx(4,ir,ik)*(fc11-(8./5.)  *ft11-4.*fb11)
       wtnm(3,2,ir,ik)=rx(ir)*jlx(5,ir,ik)*(fc01-(10./7.) *ft01-5.*fb01)
       wtnm(4,2,ir,ik)=rx(ir)*jlx(6,ir,ik)*(fc11-(4./3.)  *ft11-6.*fb11)
       wtnm(5,2,ir,ik)=rx(ir)*jlx(7,ir,ik)*(fc01-(14./11.)*ft01-7.*fb01)
         wtnm(1,1,ir,ik)=a1*rx(ir)*jlx(1,ir,ik)*ft01
         wtnm(2,1,ir,ik)=a2*rx(ir)*jlx(2,ir,ik)*ft11
         wtnm(3,1,ir,ik)=a3*rx(ir)*jlx(3,ir,ik)*ft01
	 wtnm(4,1,ir,ik)=a4*rx(ir)*jlx(4,ir,ik)*ft11
	 wtnm(5,1,ir,ik)=a5*rx(ir)*jlx(5,ir,ik)*ft01
        else
         unm(5,ir,ik)=rx(ir)*jlx(2,ir,ik)
         unm(6,ir,ik)=rx(ir)*jlx(2,ir,ik)
         unm(7,ir,ik)=rx(ir)*jlx(3,ir,ik)
         unm(8,ir,ik)=rx(ir)*jlx(4,ir,ik)
	 unm(13,ir,ik)=rx(ir)*jlx(5,ir,ik)
	 unm(14,ir,ik)=rx(ir)*jlx(6,ir,ik)
	 unm(15,ir,ik)=rx(ir)*jlx(7,ir,ik)
         utnm(1:5,2,ir,ik)=0.
         utnm(2,1,ir,ik)=rx(ir)*jlx(2,ir,ik)
         utnm(3,1,ir,ik)=rx(ir)*jlx(3,ir,ik)
	 utnm(4,1,ir,ik)=rx(ir)*jlx(4,ir,ik)
	 utnm(5,1,ir,ik)=rx(ir)*jlx(5,ir,ik)
         wtnm(1,2,ir,ik)=rx(ir)*jlx(3,ir,ik)
         wtnm(2,2,ir,ik)=rx(ir)*jlx(4,ir,ik)
         wtnm(3,2,ir,ik)=rx(ir)*jlx(5,ir,ik)
	 wtnm(4,2,ir,ik)=rx(ir)*jlx(6,ir,ik)
	 wtnm(5,2,ir,ik)=rx(ir)*jlx(7,ir,ik)
         wtnm(1:5,1,ir,ik)=0.
	 if(nm.eq.2)then
	 unm(7,ir,ik)=0.
	 unm(13,ir,ik)=0.
	 unm(15,ir,ik)=0.
	 utnm(1,2,ir,ik)=0.
	 utnm(3,2,ir,ik)=0.
	 utnm(5,2,ir,ik)=0.
	 utnm(3,1,ir,ik)=0.
	 utnm(5,1,ir,ik)=0.
	 wtnm(1,2,ir,ik)=0.
	 wtnm(3,2,ir,ik)=0.
	 wtnm(5,2,ir,ik)=0.
	 wtnm(1,1,ir,ik)=0.
	 wtnm(3,1,ir,ik)=0.
	 wtnm(5,1,ir,ik)=0.
	 end if
        end if

      end do
      end do

      c2k(1:15,1:kgrid)=0.
      c2nm(1:15,1:kgrid)=0.
      c2ck(1:5,1:2,1:kgrid)=0.
      c2cnm(1:5,1:2,1:kgrid)=0.

      do irc=1,lgrid

c      if(irc.gt.lt) then
        if(nm.eq.1)then
         vc10=v(irc,1)+v(irc,2)-3.*(v(irc,3)+v(irc,4))
         vc00=v(irc,1)-3.*(v(irc,2)+v(irc,3)-3.*v(irc,4))
         vc11=v(irc,1)+v(irc,2)+v(irc,3)+v(irc,4)
         vc01=v(irc,1)+v(irc,3)-3.*(v(irc,2)+v(irc,4))
         vt11=v(irc,5)+v(irc,6)
         vt01=v(irc,5)-3.*v(irc,6)
         vb11=v(irc,7)+v(irc,8)
         vb01=v(irc,7)-3.*v(irc,8)
         vq10=v(irc,9)+v(irc,10)-3.*(v(irc,11)+v(irc,12))
         vq00=v(irc,9)-3.*(v(irc,10)+v(irc,11)-3.*v(irc,12))
         vq11=v(irc,9)+v(irc,10)+v(irc,11)+v(irc,12)
         vq01=v(irc,9)+v(irc,11)-3.*(v(irc,10)+v(irc,12))
         vbb11=v(irc,13)+v(irc,14)
         vbb01=v(irc,13)-3.*v(irc,14)
        elseif(nm.eq.2)then
	 ir=lf*irc-lfh
	 rr=h*float(ir)
	 call pot(0,rr,vv,ww)
          vc10=vv(1)+vv(2)+2*(vv(15)-vv(18))
     &         -3*(vv(3)+vv(4)+2*vv(16))
          vc00=0.
          vc11=vv(1)+vv(2)+2*(vv(15)-vv(18))
     &           +vv(3)+vv(4)+2*vv(16)
          vc01=0.
          vt11=vv(5)+vv(6)+2*vv(17)
          vt01=0.
          vb11=vv(7)+vv(8)
          vb01=0.
          vq10=vv(9)+vv(10)-3*(vv(11)+vv(12))
	  vq00=0.
          vq11=vv(9)+vv(10)+vv(11)+vv(12)
	  vq01=0.
          vbb11=vv(13)+vv(14)
	  vbb01=0.
        end if

	
c      else
c	if(nm.eq.1)then
c        vc10=v(irc,1)+ast*(v(irc,2)-3.*(v(irc,3)+v(irc,4)))
c        vc00=v(irc,1)-3.*ast*(v(irc,2)+v(irc,3)-3.*v(irc,4))
c        vc11=v(irc,1)+ast*(v(irc,2)+v(irc,3)+v(irc,4))
c        vc01=v(irc,1)+ast*(v(irc,3)-3.*(v(irc,2)+v(irc,4)))
c        vt11=atn*(v(irc,5)+v(irc,6))
c        vt01=atn*(v(irc,5)-3.*v(irc,6))
c        vb11=als*(v(irc,7)+v(irc,8))
c        vb01=als*(v(irc,7)-3.*v(irc,8))
c        vq10=v(irc,9)+ast*(v(irc,10)-3.*(v(irc,11)+v(irc,12)))
c        vq00=v(irc,9)-3.*ast*(v(irc,10)+v(irc,11)-3.*v(irc,12))
c        vq11=v(irc,9)+ast*(v(irc,10)+v(irc,11)+v(irc,12))
c        vq01=v(irc,9)+ast*(v(irc,11)-3.*(v(irc,10)+v(irc,12)))
c        vbb11=v(irc,13)+v(irc,14)
c        vbb01=v(irc,13)-3.*v(irc,14)
c       elseif(nm.eq.2)then
c	 ir=lf*irc-lfh
c	 rr=h*float(ir)
c	 call pot(0,rr,vv,ww)
c         vc10=acn*(vv(1)+vv(2)+2*(vv(15)-vv(18)))
c    &         -3*ast*(vv(3)+vv(4)+2*vv(16))
c         vc11=acn*(vv(1)+vv(2)+2*(vv(15)-vv(18)))
c    &           +ast*(vv(3)+vv(4)+2*vv(16))
c         vt11=atn*(vv(5)+vv(6)+2*vv(17))
c         vb11=als*(vv(7)+vv(8))
c         vq10=1.0*(vv(9)+vv(10))-3*ast*(vv(11)+vv(12))
c         vq11=1.0*(vv(9)+vv(10))+ast*(vv(11)+vv(12))
c         vbb11=1.0*(vv(13)+vv(14))
c         vc00=0.
c         vc01=0.
c         vt01=0.
c         vb01=0.
c         vq00=0.
c         vq01=0.
c         vbb01=0.
c       end if
c      end if

         vtsj(1,irc) =vc10
         vtsj(2,irc) =vc00+ 2.*vq00
         vtsj(3,irc) =vc10+ 6.*vq10
         vtsj(4,irc) =vc00+12.*vq00
         vtsj(5,irc) =vc11+ 2.*vq11-4.*vt11-2.*vb11+4.*vbb11
         vtsj(6,irc) =vc11+ 2.*vq11+2.*vt11-   vb11+   vbb11
         vtsj(7,irc) =vc01+ 6.*vq01+2.*vt01-   vb01+   vbb01
         vtsj(8,irc) =vc11+12.*vq11+2.*vt11-   vb11+   vbb11
	 vtsj(9,irc) =vc10+20.*vq10
	 vtsj(10,irc)=vc00+30.*vq00
	 vtsj(11,irc)=vc10+42.*vq10
	 vtsj(12,irc)=vc00+56.*vq00
	 vtsj(13,irc)=vc01+20.*vq01+2.*vt01-   vb01+   vbb01
	 vtsj(14,irc)=vc11+30.*vq11+2.*vt11-   vb11+   vbb11
	 vtsj(15,irc)=vc01+42.*vq01+2.*vt01-   vb01+   vbb01

         vtsjm(1,irc)=vc01
         vtsjp(1,irc)=vc01+ 6.*vq01 -2.       *vt01 -3.*vb01 + 9.*vbb01
         vtsj0(1,irc)=vt01
         vtsjm(2,irc)=vc11+ 2.*vq11 -(2./5.)  *vt11 +   vb11 +    vbb11
         vtsjp(2,irc)=vc11+12.*vq11 -(8./5.)  *vt11 -4.*vb11 +16.*vbb11
         vtsj0(2,irc)=vt11
         vtsjm(3,irc)=vc01+ 6.*vq01 -(4./7.)  *vt01 +2.*vb01 + 4.*vbb01
         vtsjp(3,irc)=vc01+20.*vq01 -(10./7.) *vt01 -5.*vb01 +25.*vbb01
         vtsj0(3,irc)=vt01
	 vtsjm(4,irc)=vc11+12.*vq11 -(2./3.)  *vt11 +3.*vb11 + 9.*vbb11
	 vtsjp(4,irc)=vc11+30.*vq11 -(4./3.)  *vt11 -6.*vb11 +36.*vbb11
	 vtsj0(4,irc)=vt11
	 vtsjm(5,irc)=vc01+20.*vq01 -(8./11.) *vt01 +4.*vb01 +16.*vbb01
	 vtsjp(5,irc)=vc01+42.*vq01 -(14./11.)*vt01 -7.*vb01 +49.*vbb01
	 vtsj0(5,irc)=vt01


         ir=lf*irc-lfh

         do ik=1,kgrid
         do ip=1,8

          ukdp(ip,irc,ik)=(uk(ip,ir+1,ik)-2.*uk(ip,ir,ik)
     &                    +uk(ip,ir-1,ik))/(h*h)
	  unmdp(ip,irc,ik)=(unm(ip,ir+1,ik)-2.*unm(ip,ir,ik)
     &                    +unm(ip,ir-1,ik))/(h*h)
	   if(ip.le.4)then
	   il=ip-1
	   ij=il
	   elseif(ip.eq.5)then
	   il=ip-4
	   ij=il-1
	   elseif(ip.gt.5)then
	   il=ip-5
	   ij=il
	   end if

	   if(ip.eq.2.or.ip.eq.4.or.ip.eq.7)then
	   it=0
	   else
	   it=1
	   end if


           if(nm.eq.1)then
            xjt=(2.*ij+1.)*(2.*it+1.)*dr
           elseif(nm.eq.2)then
	    xjt=(2.*ij+1.)*dr
           end if

	   c2k(ip,ik)=c2k(ip,ik) + xjt*
     & uk(ip,ir,ik)*(
     & -h2m*(ukdp(ip,irc,ik)+(ks(ik)-(il+1.)*il/rs(irc))
     & *uk(ip,ir,ik))
     & +vtsj(ip,irc)*uk(ip,ir,ik) )

	   c2nm(ip,ik)=c2nm(ip,ik) + xjt*
     & unm(ip,ir,ik)*(
     & -h2m*(unmdp(ip,irc,ik)+(ks(ik)-(il+1.)*il/rs(irc))
     & *unm(ip,ir,ik))
     & +vtsj(ip,irc)*unm(ip,ir,ik) )

         end do

	 do ip=9,15
          unmdp(ip,irc,ik)=(unm(ip,ir+1,ik)-2.*unm(ip,ir,ik)
     &                    +unm(ip,ir-1,ik))/(h*h)
	  if(ip.le.12)then
	    il=ip-5
          else
	    il=ip-9
          end if

	  ij=il
	  if(ip.eq.9.or.ip.eq.11.or.ip.eq.14)then
	    it=1
          elseif(ip.eq.10.or.ip.eq.12.or.ip.eq.13.or.ip.eq.15)then
	    it=0
          endif

	  vp=vtsj(ip,irc)

           if(nm.eq.1)then
            xjt=(2.*ij+1.)*(2.*it+1.)*dr
           elseif(nm.eq.2)then
            xjt=(2.*ij+1.)*dr
           end if

           c2nm(ip,ik)=c2nm(ip,ik) + xjt*
     & unm(ip,ir,ik)*(
     & -h2m*(unmdp(ip,irc,ik)+(ks(ik)-(il+1.)*il/rs(irc))
     & *unm(ip,ir,ik))
     & +vp*unm(ip,ir,ik) )
         end do

         do ij=1,5
	 do jx=1,2
          utnmdp(ij,jx,irc,ik)= (utnm(ij,jx,ir+1,ik)
     &         -2.*utnm(ij,jx,ir,ik)+utnm(ij,jx,ir-1,ik))/(h*h)
	  wtnmdp(ij,jx,irc,ik)= (wtnm(ij,jx,ir+1,ik)
     &         -2.*wtnm(ij,jx,ir,ik)+wtnm(ij,jx,ir-1,ik))/(h*h)

	  if(ij.le.3)then
          utkdp(ij,jx,irc,ik)= (utk(ij,jx,ir+1,ik)
     &         -2.*utk(ij,jx,ir,ik)+utk(ij,jx,ir-1,ik))/(h*h)
          wtkdp(ij,jx,irc,ik)= (wtk(ij,jx,ir+1,ik)
     &         -2.*wtk(ij,jx,ir,ik)+wtk(ij,jx,ir-1,ik))/(h*h)
          end if

          ijm=ij-1
	  ijp=ij+1
	  if(ij.eq.1.or.ij.eq.3.or.ij.eq.5)then
	   it=0
          elseif(ij.eq.2.or.ij.eq.4)then
	   it=1
          end if

          if(nm.eq.1)then
           xjt=(2.*ij+1.)*(2.*it+1.)*dr
          elseif(nm.eq.2)then
           xjt=(2.*ij+1.)*dr
          end if

          c2cnm(ij,jx,ik)=c2cnm(ij,jx,ik)+ xjt*
     &    (
     &     utnm(ij,jx,ir,ik)*(
     &          -h2m*(
     &                utnmdp(ij,jx,irc,ik)
     &               -((ijm*(ijm+1.)/rs(irc))-ks(ik))*utnm(ij,jx,ir,ik)
     &               )
     &          +vtsjm(ij,irc)*utnm(ij,jx,ir,ik)
     &                       )
     &    +wtnm(ij,jx,ir,ik)*(
     &          -h2m*(
     &                wtnmdp(ij,jx,irc,ik)
     &               -((ijp*(ijp+1.)/rs(irc))-ks(ik))*wtnm(ij,jx,ir,ik)
     &               )
     &          +vtsjp(ij,irc)*wtnm(ij,jx,ir,ik)
     &                       )
     &    +utnm(ij,jx,ir,ik)*vtsj0(ij,irc)*wtnm(ij,jx,ir,ik)*
     &                      (12.*sqrt(ij*(ij+1.))/(2.*ij+1.))
     &    )

          if(ij.le.3)then
          c2ck(ij,jx,ik)=c2ck(ij,jx,ik)+ xjt*
     &    (
     &     utk(ij,jx,ir,ik)*(
     &          -h2m*(
     &                utkdp(ij,jx,irc,ik)
     &               -((ijm*(ijm+1.)/rs(irc))-ks(ik))*utk(ij,jx,ir,ik)
     &               )
     &          +vtsjm(ij,irc)*utk(ij,jx,ir,ik)
     &                       )
     &    +wtk(ij,jx,ir,ik)*(
     &          -h2m*(
     &                wtkdp(ij,jx,irc,ik)
     &               -((ijp*(ijp+1.)/rs(irc))-ks(ik))*wtk(ij,jx,ir,ik)
     &               )
     &          +vtsjp(ij,irc)*wtk(ij,jx,ir,ik)
     &                       )
     &    +utk(ij,jx,ir,ik)*vtsj0(ij,irc)*wtk(ij,jx,ir,ik)*
     &                      (12.*sqrt(ij*(ij+1.))/(2.*ij+1.))
     &    )
          end if
	 end do
	 end do

         end do
      end do

      dkx=kf/float(kxgrid)
      dkxx=kf/float(kxxgrid)

      c2(1:15)=0.
      do ip=1,15

      do i=1,kgrid
       w(i)=c2nm(ip,i)-c2k(ip,i)
      end do
      do i=1,kxgrid
       kx(i)=dkx*(i-0.5)
      end do
      call spline(k,w,kgrid,wp)
      do i=1,kxxgrid
       kxx(i)=dkxx*(i-0.5)
       call splint(k,w,wp,kgrid,kxx(i),wxx(i))
      end do

      kint(0)=0.
      do i=1,kxxgrid
       kint(i)=kint(i-1)+kxx(i)*wxx(i)*dkxx
      end do

      do m=1,kxgrid
      do n=1,kxgrid
      c2(ip)=c2(ip)+(2.*dkx*dkx/(pi*pi*pi*rho))*kx(m)*kx(n)*
     & (kint(n+m-1)-kint(abs(n-m)))
      end do
      end do

      end do

      c2c(1:5,1:2)=0.
      do ij=1,5
      do jx=1,2

      do i=1,kgrid
       w(i)=c2cnm(ij,jx,i)-c2ck(ij,jx,i)
      end do
      call spline(k,w,kgrid,wp)
      do i=1,kxxgrid
       call splint(k,w,wp,kgrid,kxx(i),wxx(i))
      end do

      kint(0)=0.
      do i=1,kxxgrid
       kint(i)=kint(i-1)+kxx(i)*wxx(i)*dkxx
      end do

      do m=1,kxgrid
      do n=1,kxgrid
      c2c(ij,jx)=c2c(ij,jx)+(2.*dkx*dkx/(pi*pi*pi*rho))*kx(m)*kx(n)*
     & (kint(n+m-1)-kint(abs(n-m)))
      end do
      end do

      end do
      end do

      c210=c2(1)+c2(3)+c2(9)+c2(11)
      c200=c2(2)+c2(4)+c2(10)+c2(12)
      c211=c2(5)+c2(6)+c2(8)+c2(14)+c2c(2,1)+c2c(4,1)
     &+c2c(2,2)+c2c(4,2)
      c201=c2(7)+c2(13)+c2(15)+c2c(1,1)+c2c(3,1)+c2c(5,1)
     &+c2c(1,2)+c2c(3,2)+c2c(5,2)
      c2ts=c210+c200+c211+c201

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c print out results
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i=1,ngrid
c     if(i.le.lcx) then
c     write(10,3600) rx(i),(psis(1,i,kg)/rx(i),kg=1,10)
c     write(14,3600) rx(i),(psit(1,i,kg)/rx(i),kg=1,10)
c     write(15,3600) rx(i),(psit(2,i,kg)/rx(i),kg=1,10)
c     write(16,3600) rx(i),(psit(3,i,kg)/rx(i),kg=1,10)
c     write(17,3600) rx(i),(psit(4,i,kg)/rx(i),kg=1,10)
c     end if
c     write(50,3600) rx(i),(jlx(1,i,kg),kg=1,10)
c     write(51,3600) rx(i),(jlx(2,i,kg),kg=1,10)
c     write(52,3600) rx(i),(jlx(3,i,kg),kg=1,10)
c     write(53,3600) rx(i),(jlx(4,i,kg),kg=1,10)
c     end do

      if(i.le.lcx) then
c     write(10,3600) rx(i),(psis(1,i,kg)/rx(i),kg=1,10)
c     write(11,3600) rx(i),(psis(2,i,kg)/rx(i),kg=1,10)
c     write(12,3600) rx(i),(psis(3,i,kg)/rx(i),kg=1,10)
c     write(13,3600) rx(i),(psis(4,i,kg)/rx(i),kg=1,10)
      end if
      if(i.le.ltx) then
c     write(14,3600) rx(i),(psit(1,i,kg)/rx(i),kg=1,10)
c     write(15,3600) rx(i),(psit(2,i,kg)/rx(i),kg=1,10)
c     write(16,3600) rx(i),(psit(3,i,kg)/rx(i),kg=1,10)
c     write(17,3600) rx(i),(psit(4,i,kg)/rx(i),kg=1,10)
c     write(18,3600) rx(i),(psitm(1,i,kg,1)/rx(i),kg=1,10)
c     write(19,3600) rx(i),(psitp(1,i,kg,1)/rx(i),kg=1,10)
c     write(20,3600) rx(i),(psitm(1,i,kg,2)/rx(i),kg=1,10)
c     write(21,3600) rx(i),(psitp(1,i,kg,2)/rx(i),kg=1,10)
c     write(22,3600) rx(i),(psitm(2,i,kg,1)/rx(i),kg=1,10)
c     write(23,3600) rx(i),(psitp(2,i,kg,1)/rx(i),kg=1,10)
c     write(24,3600) rx(i),(psitm(2,i,kg,2)/rx(i),kg=1,10)
c     write(25,3600) rx(i),(psitp(2,i,kg,2)/rx(i),kg=1,10)
c     write(26,3600) rx(i),(psitm(3,i,kg,1)/rx(i),kg=1,10)
c     write(27,3600) rx(i),(psitp(3,i,kg,1)/rx(i),kg=1,10)
c     write(28,3600) rx(i),(utnm(1,1,i,kg)/rx(i),kg=1,10)
c     write(29,3600) rx(i),(wtnm(1,1,i,kg)/rx(i),kg=1,10)
c     write(30,3600) rx(i),(utnm(1,2,i,kg)/rx(i),kg=1,10)
c     write(31,3600) rx(i),(wtnm(1,2,i,kg)/rx(i),kg=1,10)
c     write(32,3600) rx(i),(utnm(2,1,i,kg)/rx(i),kg=1,10)
c     write(33,3600) rx(i),(wtnm(2,1,i,kg)/rx(i),kg=1,10)
c     write(34,3600) rx(i),(utnm(2,2,i,kg)/rx(i),kg=1,10)
c     write(35,3600) rx(i),(wtnm(2,2,i,kg)/rx(i),kg=1,10)
c     write(36,3600) rx(i),(utnm(3,1,i,kg)/rx(i),kg=1,10)
c     write(37,3600) rx(i),(wtnm(3,1,i,kg)/rx(i),kg=1,10)
c     write(40,3600) rx(i),(unm(1,i,kg)/rx(i),kg=1,10)
c     write(41,3600) rx(i),(unm(2,i,kg)/rx(i),kg=1,10)
c     write(42,3600) rx(i),(unm(3,i,kg)/rx(i),kg=1,10)
c     write(43,3600) rx(i),(unm(4,i,kg)/rx(i),kg=1,10)
c     write(44,3600) rx(i),(unm(5,i,kg)/rx(i),kg=1,10)
c     write(45,3600) rx(i),(unm(6,i,kg)/rx(i),kg=1,10)
c     write(46,3600) rx(i),(unm(7,i,kg)/rx(i),kg=1,10)
c     write(47,3600) rx(i),(unm(8,i,kg)/rx(i),kg=1,10)
      end if

      end do

c     do ik=1,10
c     write(32,4100) (c2k(kg,ik),kg=1,4)
c     end do
c     do ik=1,10
c     write(32,4100) (c2k(kg,ik),kg=5,8)
c     end do
c     do ik=1,10
c     write(32,4100) (c2nm(kg,ik),kg=1,4)
c     end do
c     do ik=1,10
c     write(32,4100) (c2nm(kg,ik),kg=5,8)
c     end do
c     do ik=1,10
c     write(32,4100) (c2cnm(kg,ik),kg=1,3)
c     end do

      write(nout,4300) 0.,c2ts,c210,c200,c211,c201

      write(nout,4201) (c2(i),i=1,4)
      write(nout,4202) (c2(i),i=5,8)
      xjsum=0.d0
      do i=1,8
        xjsum=xjsum+c2(i)
      end do
c.....write(nout,4203) (c2(i),i=9,12)
c.....write(nout,4204) (c2(i),i=13,15)
      write(nout,4205) (c2c(i,1),i=1,5)
        xjsum=xjsum+c2c(2,1)
      write(nout,4206)xjsum
c.....write(nout,4205) (c2c(i,2),i=1,5)

3600   format(11f8.4)
3605   format(9f10.4)
3700   format(3i3,2f14.7)
4000   format(10f8.4)
4100   format(4f14.7)
4200   format(f14.7,1x,f14.7,1x,f14.7,1x,f14.7)
4201   format(/6x,'1S0',11x,'1P1',11x,'1D2',11x,'1F3'/5f14.7)
4202   format(/6x,'3P0',11x,'3P1',11x,'3D2',11x,'3F3'/5f14.7)
4203   format(/6x,'1G4',11x,'1H5',11x,'1I6',11x,'1J7'/5f14.7)
4204   format(/6x,'3G4',11x,'3H5',11x,'3I6',11x,'3J7'/5f14.7)
4205   format(/6x,'3S1-3D1',7x,'3P2-3F2',7x,'3D3-3G3',7x,'3F4-3H4',
     &         7x,'3G5-3I5'/5f14.7)
4206   format(/6x,'de2b',/f14.7)
4300   format(7f8.3)

      end
