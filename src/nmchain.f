c *id* nmchain *********************************************************
c subroutine for calculating chain contributions
c ----------------------------------------------------------------------
      subroutine nmchain(nv,nt,no,lg,l3)
      use nmvar
      implicit none
      !implicit real*8 (a-h,o-z)
      !implicit integer*4 (i-n)
      !include "params.f"
      integer*4 :: nv,nt,no,lg,l3
      integer*4 :: li,kj,j,lj,l,k,i,m,n,mp,ir,ka,kb,jk
      real*8    :: wvc,wvcs,wfcdd,wfccc,wfcsdd,wfcscc,wpcdd,wpccc
      real*8    :: wpcsdd,wpcscc,x,vd,ve,vp,qx
      real*8    :: z,ze,zf,zfe,zk,zke,zp,zpe,zj,zje
      real*8    :: x0,x1,x2,x3,x4,x5,x6,x7,x8,xcc,xdd,xde,xee
      real*8    :: y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11
      real*8    :: ydd,yca,ycb,yde,yee,qe,qd,g1,g2,g3,g4,g5
      real*8    :: wkc,wkcs,wjc,wjcs,wfc,wfcs,wpc,wpcs,wvcb,wkcb,wjcb
      real*8    :: wvcd,wvcds,wfcddd,wfcdcc,wfcdsd,wfcdsc,wpcddd,wpcdcc
      real*8    :: wpcdsd,wpcdsc,wkcd,wkcds,wjcd,wjcds,wfcd,wfcds
      real*8    :: wpcd,wpcds,wvcm,wvcms,wfcmcc,wfcmdd,wfcmsc,wfcmsd
      real*8    :: wpcmcc,wpcmdd,wpcmsc,wpcmsd,wkcm,wkcms,wjcm,wjcms
      real*8    :: wfcm,wfcms,wpcm,wpcms,wvcr,wfcr,wpcr,wkcr,wjcr
      real*8    :: xx,yy,zz,z1,z2,z3,z4,z5
      real*8    :: di1dk1,ditdkt,dli,dlj,dlk,d3
      real*8    :: acex !entry point
      real*8    :: ac   !function
c ----------------------------------------------------------------------
      real*8 wckx(10),wcjx(10),wcdkx(10),wcdjx(10),wcmjx(10),wcmkx(10)
     &,wcrjx(10),wcrkx(10)
      li=max(lg,l3)
      wvc=0
      wvcs=0
      wfcdd=0
      wfccc=0
      wfcsdd=0
      wfcscc=0
      wpcdd=0
      wpccc=0
      wpcsdd=0
      wpcscc=0
      w3v0=0
      w3v1=0
      w3va=0
      w3vc=0
      do 50 kj=1,8*10
        wcx(kj,1)=0
        wcdx(kj,1)=0
        wcmx(kj,1)=0
        wcrx(kj,1)=0
   50 continue
      do 52 j=1,10
        wckx(j)=0
        wcjx(j)=0
        wcdkx(j)=0
        wcdjx(j)=0
        wcmkx(j)=0
        wcmjx(j)=0
        wcrkx(j)=0
        wcrjx(j)=0
   52 continue
      x=2*rho*pi*dr
      do 550 j=1,6,nm
        lj=1
        if (j.eq.1) lj=2
        do 540 l=1+nm,6,nm
          vd=vc(l,lj,1)
          ve=vc(l,lj,2)
          vp=vc(l,lj,3)
          do 530 i=1,6,nm
          do 530 k=1,i,nm
            qx=2*x
            if (i.eq.k) qx=x
            x1=0
            x2=0
            x3=0
            x4=0
            x5=0
            x6=0
            x7=0
            if (i.ge.5.and.k.ge.5.and.l.lt.5) x7=6
            if (i.ge.5.and.k.ge.5.and.l.ge.5) x7=3
            do 510 m=1,6,nm
              y1=ak(i,j,m)*ak(k,l,m)*aa(m)
              y2=ak(i,j,m)*al(k,l,m)
              y3=ak(j,k,m)*al(i,l,m)
              y4=ak(i,k,m)*al(j,l,m)
              x1=x1+(11*y1+5*(y2+y3)+3*y4)/24
              x2=x2+.5*y1+.25*(y2+y3)
              x3=x3+y1
              do 500 n=1,4,nm
                if (l.le.4) x6=x6+ak(j,k,m)*ak(n,i,m)*aa(m)*aa(l)
                do 500 mp=1,6,nm
                  y5=ak(j,k,m)*ak(n,i,mp)*al(mp,l,m)
                  y6=ak(i,j,m)*ak(m,k,mp)*al(n,l,mp)
                  y7=ak(k,n,m)*ak(m,i,mp)*al(j,l,mp)
                  y8=ak(i,j,m)*ak(k,n,mp)*al(mp,l,m)
                  y9=ak(n,l,m)*ak(i,j,mp)*ak(m,mp,k)*aa(k)
                  y10=ak(j,k,m)*ak(m,n,mp)*al(i,l,mp)
                  y11=ak(n,i,m)*ak(m,j,mp)*al(k,l,mp)
                  x4=x4+.125*(y5+y6+y7+y8)+y9/3+(y10+y11)/12
                  x5=x5+.25*(y6+y9)+.125*(y5+y8+y10+y11)
  500           continue
  510       continue
c ------------------------
c exchange potential terms
c ------------------------
            if (x4.eq.0..and.x5.eq.0..and.x6.eq.0.) go to 521
            do 520 ir=1,li
              xdd=gx(ir)
              xcc=gl(ir)
              ydd=gdd(ir,l)
              yca=gca(ir,l)
              ycb=yca-gcb(ir,l)
              z=f(ir,i)*v(ir,j)*f(ir,k)*xdd*qx*rs(ir)
              wcx(j,4)=wcx(j,4)-z*ydd*x4*xcc**2/nu
              wcx(j,9)=wcx(j,9)-z*ydd*x4*xcc**2*ve**2/nu
              wcx(j,5)=wcx(j,5)+z*(4*yca*x5-2*ycb*x6)*xcc
              wcx(j,10)=wcx(j,10)+z*(4*yca*x5-2*ycb*x6)*xcc*ve**2
c ----------------------------
c effective 2-body part of tni
c ----------------------------
              z3=f(ir,i)*f(ir,k)*xdd*qx*rs(ir)
              w3x(j,4,1)=w3x(j,4,1)
     &         -z3*(ydd*x4*(v3dd(ir,j,1)*xcc**2/s-2*v3cc(ir,j,1)*xcc)
     &         -(4*yca*x5-2*ycb*x6)*(v3dd(ir,j,1)*xcc-s*v3cc(ir,j,1)))
     &         *ve**2
c -----------
c s-wave part
c -----------
              w3x(j,4,2)=w3x(j,4,2)
     &         -z3*(ydd*x4*(v3dd(ir,j,2)*xcc**2/s-2*v3cc(ir,j,2)*xcc)
     &         -(4*yca*x5-2*ycb*x6)*(v3dd(ir,j,2)*xcc-s*v3cc(ir,j,2)))
     &         *ve**2
c ----------------------------
              gnn(ir,j)=gnn(ir,j)+f(ir,i)*f(ir,k)*xdd*(-ydd*x4*xcc**2
     &         *ve**2/nu+(4*yca*x5-2*ycb*x6)*xcc*ve**2)*qx/x
  520       continue
c ----------------------
c direct potential terms
c ----------------------
  521       if (x1.eq.0..and.x2.eq.0..and.x3.eq.0.) go to 523
            do 522 ir=1,li
              xdd=gx(ir)
              xde=gy(ir)
              xee=gz(ir)
              ydd=gdd(ir,l)
              yde=gde(ir,l)
              yee=gee(ir,l)
              z=f(ir,i)*v(ir,j)*f(ir,k)*xdd*qx*rs(ir)
              wcx(j,1)=wcx(j,1)+z*ydd*x1*((1+xde)**2+xee)
              wcx(j,6)=wcx(j,6)+z*ydd*x1*(vd**2+2*xde*vd*vp
     &         +(xde**2+xee)*vp**2)
              wcx(j,2)=wcx(j,2)+2*z*yde*x2*(1+xde)
              wcx(j,7)=wcx(j,7)+2*z*yde*x2*(vd+xde*vp)*ve
              wcx(j,3)=wcx(j,3)+z*yee*x3
              wcx(j,8)=wcx(j,8)+z*yee*x3*ve**2
c ----------------------------
c effective 2-body part of tni
c ----------------------------
              z3=f(ir,i)*f(ir,k)*xdd*qx*rs(ir)
              w3x(j,3,1)=w3x(j,3,1)
     &         +z3*(ydd*x1*(v3dd(ir,j,1)*((1+xde)**2+xee)
     &         *vd**2+2*v3de(ir,j,1)*(1+xde)*vd*vp+v3ee(ir,j,1)*vp**2)
     &         +2*yde*x2*(v3dd(ir,j,1)*(1+xde)*vd*ve+v3de(ir,j,1)*vp*ve)
     &         +yee*x3*v3dd(ir,j,1)*ve**2)
c -----------
c s-wave part
c -----------
              w3x(j,3,2)=w3x(j,3,2)
     &         +z3*(ydd*x1*(v3dd(ir,j,2)*((1+xde)**2+xee)
     &         *vd**2+2*v3de(ir,j,2)*(1+xde)*vd*vp+v3ee(ir,j,2)*vp**2)
     &         +2*yde*x2*(v3dd(ir,j,2)*(1+xde)*vd*ve+v3de(ir,j,2)*vp*ve)
     &         +yee*x3*v3dd(ir,j,2)*ve**2)
c ----------------------------
              gnn(ir,j)=gnn(ir,j)+f(ir,i)*f(ir,k)*xdd*
     &         (ydd*x1*(vd**2+2*xde*vd*vp+(xde**2+xee)*vp**2)
     &         +2*yde*x2*(vd+xde*vp)*ve
     &         +yee*x3*ve**2)*qx/x
  522       continue
c -------------------------
c exchange kinetic energies
c -------------------------
  523 if (j.gt.1) go to 530
      if (x4.eq.0..and.x5.eq.0..and.x6.eq.0.) go to 525
      do 524 ir=1,li
      xdd=gx(ir)
      xcc=gl(ir)
      ydd=gdd(ir,l)
      yca=gca(ir,l)
      ycb=yca-gcb(ir,l)
      zf=h2m*(f(ir,i)*fp(ir,k)+f(ir,k)*fp(ir,i))*xdd*slp(ir)*rs(ir)*qx
      wfcdd=wfcdd+zf*ydd*x4*xcc/nu
      wfcsdd=wfcsdd+zf*ydd*x4*xcc*(ve**2-1)/nu
      wfccc=wfccc-zf*(2*yca*x5-ycb*x6)
      wfcscc=wfcscc-zf*(2*yca*x5-ycb*x6)*(ve**2-1)
      zk=-.5*h2m*(f(ir,i)*fds(ir,k)+f(ir,k)*fds(ir,i))*xdd*qx*rs(ir)
      wckx(4)=wckx(4)-zk*ydd*x4*xcc**2/nu
      wckx(9)=wckx(9)-zk*ydd*x4*xcc**2*ve**2/nu
      wckx(5)=wckx(5)+zk*(4*yca*x5-2*ycb*x6)*xcc
      wckx(10)=wckx(10)+zk*(4*yca*x5-2*ycb*x6)*xcc*ve**2
      zp=.5*h2m*(f(ir,i)*f(ir,k)*xdd-at(i,1)*at(k,1))*qx
      wpcdd=wpcdd-zp*ydd*x4*rlssx(ir)/nu
      wpcsdd=wpcsdd-zp*ydd*x4*(ve**2-1)*rlssx(ir)/nu
      wpccc=wpccc+zp*(2*yca*x5-ycb*x6)*rsdsl(ir)
      wpcscc=wpcscc+zp*(2*yca*x5-ycb*x6)*(ve**2-1)*rsdsl(ir)
      zj=h2m*(fp(ir,i)*fp(ir,k)*rs(ir)+x7*f(ir,i)*f(ir,k))*qx*xdd
      wcjx(4)=wcjx(4)-zj*ydd*x4*xcc**2/nu
      wcjx(9)=wcjx(9)-zj*ydd*x4*xcc**2*ve**2/nu
      wcjx(5)=wcjx(5)+zj*(4*yca*x5-2*ycb*x6)*xcc
      wcjx(10)=wcjx(10)+zj*(4*yca*x5-2*ycb*x6)*xcc*ve**2
  524 continue
c -----------------------
c direct kinetic energies
c -----------------------
  525 if (x1.eq.0..and.x2.eq.0..and.x3.eq.0.) go to 530
      do 526 ir=1,li
      xdd=gx(ir)
      xde=gy(ir)
      xee=gz(ir)
      ydd=gdd(ir,l)
      yde=gde(ir,l)
      yee=gee(ir,l)
      zk=-.5*h2m*(f(ir,i)*fds(ir,k)+f(ir,k)*fds(ir,i))*xdd*qx*rs(ir)
      wckx(1)=wckx(1)+zk*ydd*x1*((1+xde)**2+xee)
      wckx(6)=wckx(6)+zk*ydd*x1*(vd**2+2*xde*vd*vp
     & +(xde**2+xee)*vp**2)
      wckx(2)=wckx(2)+2*zk*yde*x2*(1+xde)
      wckx(7)=wckx(7)+2*zk*yde*x2*(vd+xde*vp)*ve
      wckx(3)=wckx(3)+zk*yee*x3
      wckx(8)=wckx(8)+zk*yee*x3*ve**2
      zj=h2m*(fp(ir,i)*fp(ir,k)*rs(ir)+x7*f(ir,i)*f(ir,k))*qx*xdd
      wcjx(1)=wcjx(1)+zj*ydd*x1*((1+xde)**2+xee)
      wcjx(6)=wcjx(6)+zj*ydd*x1*(vd**2+2*xde*vd*vp
     & +(xde**2+xee)*vp**2)
      wcjx(2)=wcjx(2)+2*zj*yde*x2*(1+xde)
      wcjx(7)=wcjx(7)+2*zj*yde*x2*(vd+xde*vp)*ve
      wcjx(3)=wcjx(3)+zj*yee*x3
      wcjx(8)=wcjx(8)+zj*yee*x3*ve**2
  526 continue
  530 continue
  540 continue
cdir$ ivdep
      do 545 m=5,1,-1
        wcx(j,m+5)=wcx(j,m+5)-wcx(j,m)
  545 continue
      wvc=wvc+wcx(j,1)+wcx(j,2)+wcx(j,3)+wcx(j,4)+wcx(j,5)
      wvcs=wvcs+wcx(j,6)+wcx(j,7)+wcx(j,8)+wcx(j,9)+wcx(j,10)
      w3va=w3va+w3x(j,1,1)+w3x(j,2,1)+w3x(j,3,1)+w3x(j,4,1)
      w3v1=w3v1+w3x(j,1,2)+w3x(j,2,2)+w3x(j,3,2)+w3x(j,4,2)
  550 continue
      do 555 m=1,10
  555 wcjx(m)=.5*(wcjx(m)+wckx(m))
      do 560 m=5,1,-1
      wckx(m+5)=wckx(m+5)-wckx(m)
  560 wcjx(m+5)=wcjx(m+5)-wcjx(m)
      wkc=wckx(1)+wckx(2)+wckx(3)+wckx(4)+wckx(5)
      wkcs=wckx(6)+wckx(7)+wckx(8)+wckx(9)+wckx(10)
      wjc=wcjx(1)+wcjx(2)+wcjx(3)+wcjx(4)+wcjx(5)
      wjcs=wcjx(6)+wcjx(7)+wcjx(8)+wcjx(9)+wcjx(10)
      wfc=wfcdd+wfccc
      wfcs=wfcsdd+wfcscc
      wpc=wpcdd+wpccc
      wpcs=wpcsdd+wpcscc
c -------------
c calculate wcb
c -------------
      if (nv.le.6) go to 590
      qe=x**2*dr/(8*s)
      qd=x**2*dr*ksav/12
      wvcb=wcx(1,1)+wcx(1,2)
      wkcb=wckx(1)+wckx(2)
      wjcb=wcjx(1)+wcjx(2)
      do 580 i=1,lg
        x0=gx(i)*rs(i)
        x1=v(i,7)*f(i,1)**2*x0
        x2=2*f(i,7)*v(i,1)*f(i,1)*x0
        x3=-h2m*(f(i,7)*fds(i,1)+f(i,1)*fds(i,7))*x0
        x4=3*v(i,8)*f(i,1)**2*x0
        x5=6*f(i,8)*v(i,1)*f(i,1)*x0
        x6=-3*h2m*(f(i,8)*fds(i,1)+f(i,1)*fds(i,8))*x0
        x7=2*h2m*fp(i,1)*(fp(i,7)+f(i,7)/r(i))*x0
        x8=6*h2m*fp(i,1)*(fp(i,8)+f(i,8)/r(i))*x0
        g1=f(i,1)**2*gx(i)/x
        g2=2*f(i,7)*f(i,1)*gx(i)/x
        g4=3*f(i,1)**2*gx(i)/x
        g5=6*f(i,8)*f(i,1)*gx(i)/x
        do 580 j=1,lg
          ka=iabs(i-j)+1
          kb=min(i+j-1,lg)
          y0=f(j,1)*gx(j)*r(j)
          y1=f(j,7)*y0
          y2=f(j,8)*y0
          y3=2*y1*r(j)
          y4=2*y2*r(j)
          do 580 k=ka,kb
          xx=(rs(i)+rs(j)-rs(k))/(2*r(i)*r(j))
          yy=(rs(i)+rs(k)-rs(j))/(2*r(i)*r(k))
          zz=(rs(j)+rs(k)-rs(i))/(2*r(j)*r(k))
          z1=f(k,1)**2*gx(k)
          z2=qd*(z1-1)*r(k)*xx
          z3=2*qd*z1*sls(k)*r(k)/nu
          z4=qe*z1*(sl(k)*(sldp(k)*r(k)*(xx+yy*zz)+slp(k)*(xx-yy*zz))
     &     -slps(k)*(xx+yy*zz)*r(k))
          z5=4*qe*z1*rllp(k)*yy
          wcx(1,1)=wcx(1,1)+x2*y3*z2
          wcx(7,1)=wcx(7,1)+x1*y3*z2
          wckx(1)=wckx(1)+x3*y3*z2
          wcjx(1)=wcjx(1)+.5*(x3+x7)*y3*z2
          wcx(1,2)=wcx(1,2)+(x2*y3+x5*y4)*(z3+z4)
          wcx(7,2)=wcx(7,2)+x1*(y3*(z3+z4)+y1*z5)
          wcx(8,2)=wcx(8,2)+x4*(y4*(z3+z4)+y2*z5)
          wckx(2)=wckx(2)+(x3*y3+x6*y4)*(z3+z4)
          wcjx(2)=wcjx(2)+.5*((x3+x7)*y3+(x6+x8)*y4)*(z3+z4)
          gnn(i,1)=gnn(i,1)+g2*y3*z2+(g2*y3+g5*y4)*(z3+z4)
          gnn(i,7)=gnn(i,7)+g1*(y3*(z2+z3+z4)+y1*z5)
          gnn(i,8)=gnn(i,8)+g4*(y4*(z3+z4)+y2*z5)
  580 continue
      wvcb=wcx(1,1)+wcx(7,1)+wcx(1,2)+wcx(7,2)+wcx(8,2)-wvcb
      wkcb=wckx(1)+wckx(2)-wkcb
      wjcb=wcjx(1)+wcjx(2)-wjcb
c print ================================================================
  590 if (no.eq.0) go to 600
      write(nlog,1200) wvc,wkc,wfc,wjc,wpc,wvcs,wkcs,wfcs,wjcs,wpcs
      write(nout,1200) wvc,wkc,wfc,wjc,wpc,wvcs,wkcs,wfcs,wjcs,wpcs
 1200 format(/4x,'wvc',5x,'wkc',5x,'wfc',5x,'wjc',5x,'wpc',5x,'wvcs'
     &,4x,'wkcs',4x,'wfcs',4x,'wjcs',4x,'wpcs'/10f8.3)
      if (nv.ge.7) then
        write(nlog,1205) wvcb,wkcb,wjcb
        write(nout,1205) wvcb,wkcb,wjcb
 1205   format(/4x,'wvcb',4x,'wkcb',4x,'wjcb'/3f8.3)
      end if
      write(nout,1210)
 1210 format(/4x,'wc & wcs:dd,de,ee,ddx,cc')
      write(nout,1133)
 1133 format(4x,'c',7x,'t',7x,'s',7x,'st',6x,'tn',6x,'tnt'
     &,5x,'b',7x,'bt',6x,'ke',6x,'jf')
      write(nout,1140) ((wcx(i,j),i=1,8),wckx(j),wcjx(j),j=1,10)
 1140 format(10f8.3)
      write(nout,1220) wfcdd,wfccc,wpcdd,wpccc
 1220 format(/4x,'wfc:dd(ex),cc   wpc:dd(ex),cc'/4f8.3)
c ======================================================================
c ------------------
c calculate wcd,wcds
c ------------------
  600 wvcd=0
      wvcds=0
      wfcddd=0
      wfcdcc=0
      wfcdsd=0
      wfcdsc=0
      wpcddd=0
      wpcdcc=0
      wpcdsd=0
      wpcdsc=0
      do 645 j=1,6,nm
      do 640 l=1+nm,6,nm
      do 635 m=1+nm,l,nm
      qx=2*x
      if (l.eq.m) qx=x
      x1=.5*(ak(j,l,m)*aa(m)+al(j,l,m))
      x2=(ak(j,l,m)*aa(m)+2*al(j,l,m))/3
      x3=(2*ak(j,l,m)*aa(m)+al(j,l,m))/3
      x4=0
      x5=0
      if (j.eq.1) then
        do 610 n=1,4,nm
        x4=x4+.5*(ak(n,l,m)*aa(m)+al(n,l,m))
  610   continue
        x5=aa(l)*af(l)*aa(m)*af(m)
      end if
      if (x1.eq.0..and.x2.eq.0..and.x3.eq.0..and.x4.eq.0.) go to 635
c ---------------
c potential terms
c ---------------
      vd=vc(j,1,1)
      ve=vc(j,1,2)
      vp=vc(j,1,3)
      if (j.eq.1) then
        vd=.5*(vc(l,2,1)+vc(m,2,1))
        ve=.5*(vc(l,2,2)+vc(m,2,2))
        vp=.5*(vc(l,2,3)+vc(m,2,3))
      end if
      do 620 ir=1,lg
      xdd=gx(ir)
      xde=gy(ir)
      xee=gz(ir)
      xcc=gl(ir)
      z=v(ir,j)*f(ir,1)**2*xdd*qx*rs(ir)
      wcdx(j,1)=wcdx(j,1)+.5*z*x1*gdd(ir,l)*gdd(ir,m)*((1+xde)**2+xee)
      wcdx(j,2)=wcdx(j,2)+z*x1*(gdd(ir,l)*gde(ir,m)+gde(ir,l)*gdd(ir,m))
     & *(1+xde)
      wcdx(j,3)=wcdx(j,3)+z*x1*(gde(ir,l)*gde(ir,m)
     & +.5*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))
      wcdx(j,4)=wcdx(j,4)-.5*z*x4*gdd(ir,l)*gdd(ir,m)*xcc**2/nu
      wcdx(j,5)=wcdx(j,5)+z
     & *(gdd(ir,l)*(2*x4*gca(ir,m)-x5*(gca(ir,m)-gcb(ir,m)))
     & +(2*x4*gca(ir,l)-x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))*xcc
      wcdx(j,6)=wcdx(j,6)+.5*z*x1*gdd(ir,l)*gdd(ir,m)
     & *((vd+xde*vp)**2+xee*vp**2)
      wcdx(j,7)=wcdx(j,7)+z*x1*(gdd(ir,l)*gde(ir,m)+gde(ir,l)*gdd(ir,m))
     & *(vd+xde*vp)*ve
      wcdx(j,8)=wcdx(j,8)+z*x1*(gde(ir,l)*gde(ir,m)
     & +.5*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))*ve**2
      wcdx(j,9)=wcdx(j,9)-.5*z*x4*gdd(ir,l)*gdd(ir,m)*xcc**2*ve**2/nu
      wcdx(j,10)=wcdx(j,10)+z
     & *(gdd(ir,l)*(2*x4*gca(ir,m)-x5*(gca(ir,m)-gcb(ir,m)))
     & +(2*x4*gca(ir,l)-x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))*xcc*ve**2
c     -------------------
      gnn(ir,j)=gnn(ir,j)+f(ir,1)**2*xdd*(x1*
     & (.5*gdd(ir,l)*gdd(ir,m)*((vd+xde*vp)**2+xee*vp**2)
     & +(gdd(ir,l)*gde(ir,m)+gde(ir,l)*gdd(ir,m))*(vd+xde*vp)*ve
     & +(gde(ir,l)*gde(ir,m)
     & +.5*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))*ve**2)
     & -.5*x4*gdd(ir,l)*gdd(ir,m)*xcc**2*ve**2/nu
     & +(gdd(ir,l)*(2*x4*gca(ir,m)-x5*(gca(ir,m)-gcb(ir,m)))
     & +(2*x4*gca(ir,l)-x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))
     & *xcc*ve**2)*qx/x
  620 continue
c --------------------
c kinetic energy terms
c --------------------
      if (j.eq.1) then
        do 625 ir=1,lg
        xdd=gx(ir)
        xde=gy(ir)
        xee=gz(ir)
        xcc=gl(ir)
        zk=-h2m*f(ir,1)*fds(ir,1)*xdd*qx*rs(ir)
        wcdkx(1)=wcdkx(1)+.5*zk*x1*gdd(ir,l)*gdd(ir,m)*((1+xde)**2+xee)
        wcdkx(2)=wcdkx(2)+zk*x1*(gdd(ir,l)*gde(ir,m)
     &   +gde(ir,l)*gdd(ir,m))*(1+xde)
        wcdkx(3)=wcdkx(3)+zk*x1*(gde(ir,l)*gde(ir,m)
     &   +.5*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))
        wcdkx(4)=wcdkx(4)-.5*zk*x4*gdd(ir,l)*gdd(ir,m)*xcc**2/nu
        wcdkx(5)=wcdkx(5)+zk
     &   *(gdd(ir,l)*(2*x4*gca(ir,m)-x5*(gca(ir,m)-gcb(ir,m)))
     &   +(2*x4*gca(ir,l)-x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))*xcc
        wcdkx(6)=wcdkx(6)+.5*zk*x1*gdd(ir,l)*gdd(ir,m)
     &   *((vd+xde*vp)**2+xee*vp**2)
        wcdkx(7)=wcdkx(7)+zk*x1*(gdd(ir,l)*gde(ir,m)
     &   +gde(ir,l)*gdd(ir,m))*(vd+xde*vp)*ve
        wcdkx(8)=wcdkx(8)+zk*x1*(gde(ir,l)*gde(ir,m)
     &   +.5*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))*ve**2
        wcdkx(9)=wcdkx(9)-.5*zk*x4*gdd(ir,l)*gdd(ir,m)*xcc**2*ve**2/nu
        wcdkx(10)=wcdkx(10)+zk
     &   *(gdd(ir,l)*(2*x4*gca(ir,m)-x5*(gca(ir,m)-gcb(ir,m)))
     &   +(2*x4*gca(ir,l)-x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))*xcc
     &   *ve**2
        zf=2*h2m*f(ir,1)*fp(ir,1)*slp(ir)*xdd*qx*rs(ir)
        wfcddd=wfcddd+.5*zf*x4*gdd(ir,l)*gdd(ir,m)*xcc/nu
        wfcdcc=wfcdcc-zf
     &   *(gdd(ir,l)*(x4*gca(ir,m)-.5*x5*(gca(ir,m)-gcb(ir,m)))
     &   +(x4*gca(ir,l)-.5*x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))
        wfcdsd=wfcdsd+.5*zf*x4*gdd(ir,l)*gdd(ir,m)*(ve**2-1)*xcc/nu
        wfcdsc=wfcdsc-zf
     &   *(gdd(ir,l)*(x4*gca(ir,m)-.5*x5*(gca(ir,m)-gcb(ir,m)))
     &   +(x4*gca(ir,l)-.5*x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))
     &   *(ve**2-1)
c     -------------------
        zj=h2m*fp(ir,1)**2*xdd*qx*rs(ir)
        wcdjx(1)=wcdjx(1)+.5*zj*x1*gdd(ir,l)*gdd(ir,m)*((1+xde)**2+xee)
        wcdjx(2)=wcdjx(2)+zj*x1*(gdd(ir,l)*gde(ir,m)
     &   +gde(ir,l)*gdd(ir,m))*(1+xde)
        wcdjx(3)=wcdjx(3)+zj*x1*(gde(ir,l)*gde(ir,m)
     &   +.5*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))
        wcdjx(4)=wcdjx(4)-.5*zj*x4*gdd(ir,l)*gdd(ir,m)*xcc**2/nu
        wcdjx(5)=wcdjx(5)+zj
     &   *(gdd(ir,l)*(2*x4*gca(ir,m)-x5*(gca(ir,m)-gcb(ir,m)))
     &   +(2*x4*gca(ir,l)-x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))*xcc
        wcdjx(6)=wcdjx(6)+.5*zj*x1*gdd(ir,l)*gdd(ir,m)
     &   *((vd+xde*vp)**2+xee*vp**2)
        wcdjx(7)=wcdjx(7)+zj*x1*(gdd(ir,l)*gde(ir,m)
     &   +gde(ir,l)*gdd(ir,m))*(vd+xde*vp)*ve
        wcdjx(8)=wcdjx(8)+zj*x1*(gde(ir,l)*gde(ir,m)
     &   +.5*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))*ve**2
        wcdjx(9)=wcdjx(9)-.5*zj*x4*gdd(ir,l)*gdd(ir,m)*xcc**2*ve**2/nu
        wcdjx(10)=wcdjx(10)+zj
     &   *(gdd(ir,l)*(2*x4*gca(ir,m)-x5*(gca(ir,m)-gcb(ir,m)))
     &   +(2*x4*gca(ir,l)-x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))*xcc
     &   *ve**2
        zp=.5*h2m*(f(ir,1)**2*xdd-1)*qx
        wpcddd=wpcddd-.5*zp*x4*gdd(ir,l)*gdd(ir,m)*rlssx(ir)/nu
        wpcdcc=wpcdcc+zp
     &   *(gdd(ir,l)*(x4*gca(ir,m)-.5*x5*(gca(ir,m)-gcb(ir,m)))
     &   +(x4*gca(ir,l)-.5*x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))
     &   *rsdsl(ir)
        wpcdsd=wpcdsd-.5*zp*x4*gdd(ir,l)*gdd(ir,m)*rlssx(ir)/nu
     &   *(ve**2-1)
        wpcdsc=wpcdsc+zp
     &   *(gdd(ir,l)*(x4*gca(ir,m)-.5*x5*(gca(ir,m)-gcb(ir,m)))
     &   +(x4*gca(ir,l)-.5*x5*(gca(ir,l)-gcb(ir,l)))*gdd(ir,m))
     &   *rsdsl(ir)*(ve**2-1)
  625   continue
c ----------------------------
c noncentral correlation terms
c ----------------------------
      else if (j.gt.1) then
        vd=vc(j,2,1)
        ve=vc(j,2,2)
        vp=vc(j,2,3)
        do 630 ir=1,lg
        xdd=gx(ir)
        xde=gy(ir)
        xee=gz(ir)
        z=2*f(ir,1)*v(ir,1)*f(ir,j)*xdd*qx*rs(ir)
        wcdx(1,1)=wcdx(1,1)+.5*z*x1*gdd(ir,l)*gdd(ir,m)
     &   *((1+xde)**2+xee)
        wcdx(1,2)=wcdx(1,2)+z*x1*(gdd(ir,l)*gde(ir,m)+gde(ir,l)
     &   *gdd(ir,m))*(1+xde)
        wcdx(1,3)=wcdx(1,3)+z*(x2*gde(ir,l)*gde(ir,m)
     &   +.5*x3*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))
        wcdx(1,6)=wcdx(1,6)+.5*z*x1*gdd(ir,l)*gdd(ir,m)
     &   *((vd+xde*vp)**2+xee*vp**2)
        wcdx(1,7)=wcdx(1,7)+z*x1*(gdd(ir,l)*gde(ir,m)+gde(ir,l)
     &   *gdd(ir,m))*(vd+xde*vp)*ve
        wcdx(1,8)=wcdx(1,8)+z*(x2*gde(ir,l)*gde(ir,m)
     &   +.5*x3*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))*ve**2
c       --------------
        gnn(ir,1)=gnn(ir,1)+2*f(ir,1)*f(ir,j)*xdd*(x1*
     &   (.5*gdd(ir,l)*gdd(ir,m)*((vd+xde*vp)**2+xee*vp**2)
     &   +(gdd(ir,l)*gde(ir,m)+gde(ir,l)*gdd(ir,m))*(vd+xde*vp)*ve)
     &   +x2*gde(ir,l)*gde(ir,m)*ve**2
     &   +.5*x3*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m))*ve**2)*qx/x
c       --------------
        zk=-h2m*(f(ir,1)*fds(ir,j)+f(ir,j)*fds(ir,1))*xdd*qx*rs(ir)
        wcdkx(1)=wcdkx(1)+.5*zk*x1*gdd(ir,l)*gdd(ir,m)*((1+xde)**2+xee)
        wcdkx(2)=wcdkx(2)+zk*x1*(gdd(ir,l)*gde(ir,m)
     &   +gde(ir,l)*gdd(ir,m))*(1+xde)
        wcdkx(3)=wcdkx(3)+zk*(x2*gde(ir,l)*gde(ir,m)
     &   +.5*x3*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))
        wcdkx(6)=wcdkx(6)+.5*zk*x1*gdd(ir,l)*gdd(ir,m)
     &   *((vd+xde*vp)**2+xee*vp**2)
        wcdkx(7)=wcdkx(7)+zk*x1*(gdd(ir,l)*gde(ir,m)
     &   +gde(ir,l)*gdd(ir,m))*(vd+xde*vp)*ve
        wcdkx(8)=wcdkx(8)+zk*(x2*gde(ir,l)*gde(ir,m)
     &   +.5*x3*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))*ve**2
c       --------------
        zj=2*h2m*fp(ir,1)*fp(ir,j)*xdd*qx*rs(ir)
        wcdjx(1)=wcdjx(1)+.5*zj*x1*gdd(ir,l)*gdd(ir,m)*((1+xde)**2+xee)
        wcdjx(2)=wcdjx(2)+zj*x1*(gdd(ir,l)*gde(ir,m)
     &   +gde(ir,l)*gdd(ir,m))*(1+xde)
        wcdjx(3)=wcdjx(3)+zj*(x2*gde(ir,l)*gde(ir,m)
     &   +.5*x3*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))
        wcdjx(6)=wcdjx(6)+.5*zj*x1*gdd(ir,l)*gdd(ir,m)
     &   *((vd+xde*vp)**2+xee*vp**2)
        wcdjx(7)=wcdjx(7)+zj*x1*(gdd(ir,l)*gde(ir,m)
     &   +gde(ir,l)*gdd(ir,m))*(vd+xde*vp)*ve
        wcdjx(8)=wcdjx(8)+zj*(x2*gde(ir,l)*gde(ir,m)
     &   +.5*x3*(gdd(ir,l)*gee(ir,m)+gee(ir,l)*gdd(ir,m)))*ve**2
  630   continue
      end if
  635 continue
  640 continue
  645 continue
      do 655 j=1,6,nm
      do 650 m=5,1,-1
  650 wcdx(j,m+5)=wcdx(j,m+5)-wcdx(j,m)
      wvcd=wvcd+wcdx(j,1)+wcdx(j,2)+wcdx(j,3)+wcdx(j,4)+wcdx(j,5)
      wvcds=wvcds+wcdx(j,6)+wcdx(j,7)+wcdx(j,8)+wcdx(j,9)+wcdx(j,10)
  655 continue
      do 660 m=1,10
  660 wcdjx(m)=.5*(wcdjx(m)+wcdkx(m))
      do 665 m=5,1,-1
      wcdkx(m+5)=wcdkx(m+5)-wcdkx(m)
  665 wcdjx(m+5)=wcdjx(m+5)-wcdjx(m)
      wkcd=wcdkx(1)+wcdkx(2)+wcdkx(3)+wcdkx(4)+wcdkx(5)
      wkcds=wcdkx(6)+wcdkx(7)+wcdkx(8)+wcdkx(9)+wcdkx(10)
      wjcd=wcdjx(1)+wcdjx(2)+wcdjx(3)+wcdjx(4)+wcdjx(5)
      wjcds=wcdjx(6)+wcdjx(7)+wcdjx(8)+wcdjx(9)+wcdjx(10)
      wfcd=wfcddd+wfcdcc
      wfcds=wfcdsd+wfcdsc
      wpcd=wpcddd+wpcdcc
      wpcds=wpcdsd+wpcdsc
c -------------------------------------
c multiple operator chain contributions
c -------------------------------------
      wvcm=0
      wvcms=0
      wfcmcc=0
      wfcmdd=0
      wfcmsc=0
      wfcmsd=0
      wpcmcc=0
      wpcmdd=0
      wpcmsc=0
      wpcmsd=0
      do 670 ir=1,lg
      z=2*v(ir,1)*f(ir,1)**2*gx(ir)*gl(ir)*x*rs(ir)
      wcmx(1,4)=wcmx(1,4)+z*ghcc(ir,1)
      wcmx(1,9)=wcmx(1,9)+z*ghcc(ir,1)
      zk=-2*h2m*f(ir,1)*fds(ir,1)*gx(ir)*gl(ir)*x*rs(ir)
      wcmkx(4)=wcmkx(4)+zk*ghcc(ir,1)
      wcmkx(9)=wcmkx(9)+zk*ghcc(ir,1)
      zf=2*h2m*f(ir,1)*fp(ir,1)*gx(ir)*slp(ir)*x*rs(ir)
      wfcmcc=wfcmcc-zf*ghcc(ir,1)
      zj=2*h2m*fp(ir,1)**2*gx(ir)*gl(ir)*x*rs(ir)
      wcmjx(4)=wcmjx(4)+zj*ghcc(ir,1)
      wcmjx(9)=wcmjx(9)+zj*ghcc(ir,1)
      zp=.5*h2m*(f(ir,1)**2*gx(ir)-1)*rsdsl(ir)*x
      wpcmcc=wpcmcc+zp*ghcc(ir,1)
      gnn(ir,1)=gnn(ir,1)+2*f(ir,1)**2*gx(ir)*gl(ir)*ghcc(ir,1)
  670 continue
      do 690 j=1+nm,6,nm
      vd=vc(j,1,1)
      ve=vc(j,1,2)
      vp=vc(j,1,3)
      do 680 ir=1,lg
      xdd=gx(ir)
      xde=gy(ir)
      xee=gz(ir)
      xcc=gl(ir)
      z=2*v(ir,j)*f(ir,1)**2*xdd*x*rs(ir)
      wcmx(j,1)=wcmx(j,1)+z*ghdd(ir,j)*((1+xde)**2+xee)
      wcmx(j,2)=wcmx(j,2)+z*ghde(ir,j)*(1+xde)
      wcmx(j,3)=wcmx(j,3)+z*ghed(ir,j)*(1+xde)
      wcmx(j,4)=wcmx(j,4)+z*ghcc(ir,j)*xcc
      wcmx(j,6)=wcmx(j,6)+z*ghdd(ir,j)*((vd+xde*vp)**2+xee*vp**2)
      wcmx(j,7)=wcmx(j,7)+z*ghde(ir,j)*(vd+xde*vp)*ve
      wcmx(j,8)=wcmx(j,8)+z*ghed(ir,j)*(vd+xde*vp)*ve
      wcmx(j,9)=wcmx(j,9)+z*ghcc(ir,j)*xcc*ve**2
      gnn(ir,j)=gnn(ir,j)+2*f(ir,1)**2*xdd*
     & (ghdd(ir,j)*((vd+xde*vp)**2+xee*vp**2)+(ghde(ir,j)
     & +ghed(ir,j))*(vd+xde*vp)*ve+ghcc(ir,j)*xcc*ve**2)
  680 continue
      vd=vc(j,2,1)
      ve=vc(j,2,2)
      vp=vc(j,2,3)
      do 682 ir=1,lg
      xdd=gx(ir)
      xde=gy(ir)
      xee=gz(ir)
      xcc=gl(ir)
      z=4*v(ir,1)*f(ir,1)*f(ir,j)*xdd*x*rs(ir)
      ze=2*af(j)*v(ir,1)*f(ir,1)**2*xdd*x*rs(ir)
      wcmx(1,1)=wcmx(1,1)+z*gfdd(ir,j)*((1+xde)**2+xee)
      wcmx(1,2)=wcmx(1,2)+z*gfde(ir,j)*(1+xde)
      wcmx(1,3)=wcmx(1,3)+z*gfed(ir,j)*(1+xde)
      wcmx(1,4)=wcmx(1,4)+z*gfcc(ir,j)*xcc
      wcmx(1,5)=wcmx(1,5)-ze*ghdd(ir,j)*xcc**2/nu
      wcmx(1,6)=wcmx(1,6)+z*gfdd(ir,j)*((vd+xde*vp)**2+xee*vp**2)
      wcmx(1,7)=wcmx(1,7)+z*gfde(ir,j)*(vd+xde*vp)*ve
      wcmx(1,8)=wcmx(1,8)+z*gfed(ir,j)*(vd+xde*vp)*ve
      wcmx(1,9)=wcmx(1,9)+z*gfcc(ir,j)*xcc*ve**2
      wcmx(1,10)=wcmx(1,10)-ze*ghdd(ir,j)*xcc**2*ve**2/nu
c     --------------
      zk=-2*h2m*(f(ir,1)*fds(ir,j)+f(ir,j)*fds(ir,1))*xdd*x*rs(ir)
      zke=-2*h2m*af(j)*f(ir,1)*fds(ir,1)*xdd*x*rs(ir)
      wcmkx(1)=wcmkx(1)+zk*gfdd(ir,j)*((1+xde)**2+xee)
      wcmkx(2)=wcmkx(2)+zk*gfde(ir,j)*(1+xde)
      wcmkx(3)=wcmkx(3)+zk*gfed(ir,j)*(1+xde)
      wcmkx(4)=wcmkx(4)+zk*gfcc(ir,j)*xcc
      wcmkx(5)=wcmkx(5)-zke*ghdd(ir,j)*xcc**2/nu
      wcmkx(6)=wcmkx(6)+zk*gfdd(ir,j)*((vd+xde*vp)**2+xee*vp**2)
      wcmkx(7)=wcmkx(7)+zk*gfde(ir,j)*(vd+xde*vp)*ve
      wcmkx(8)=wcmkx(8)+zk*gfed(ir,j)*(vd+xde*vp)*ve
      wcmkx(9)=wcmkx(9)+zk*gfcc(ir,j)*xcc*ve**2
      wcmkx(10)=wcmkx(10)-zke*ghdd(ir,j)*xcc**2*ve**2/nu
      zf=h2m*(f(ir,1)*fp(ir,j)+f(ir,j)*fp(ir,1))*slp(ir)*xdd*x*rs(ir)
      zfe=4*h2m*af(j)*f(ir,1)*fp(ir,1)*slp(ir)*xdd*x*rs(ir)
      wfcmcc=wfcmcc-zf*gfcc(ir,j)
      wfcmsc=wfcmsc-zf*gfcc(ir,j)*(ve**2-1)
      wfcmdd=wfcmdd+zfe*ghdd(ir,j)*xcc/nu
      wfcmsd=wfcmsd+zfe*ghdd(ir,j)*xcc*(ve**2-1)/nu
c     --------------
      zj=4*h2m*fp(ir,1)*fp(ir,j)*xdd*x*rs(ir)
      zje=2*h2m*af(j)*fp(ir,1)*fp(ir,1)*xdd*x*rs(ir)
      wcmjx(1)=wcmjx(1)+zj*gfdd(ir,j)*((1+xde)**2+xee)
      wcmjx(2)=wcmjx(2)+zj*gfde(ir,j)*(1+xde)
      wcmjx(3)=wcmjx(3)+zj*gfed(ir,j)*(1+xde)
      wcmjx(4)=wcmjx(4)+zj*gfcc(ir,j)*xcc
      wcmjx(5)=wcmjx(5)-zje*ghdd(ir,j)*xcc**2/nu
      wcmjx(6)=wcmjx(6)+zj*gfdd(ir,j)*((vd+xde*vp)**2+xee*vp**2)
      wcmjx(7)=wcmjx(7)+zj*gfde(ir,j)*(vd+xde*vp)*ve
      wcmjx(8)=wcmjx(8)+zj*gfed(ir,j)*(vd+xde*vp)*ve
      wcmjx(9)=wcmjx(9)+zj*gfcc(ir,j)*xcc*ve**2
      wcmjx(10)=wcmjx(10)-zje*ghdd(ir,j)*xcc**2*ve**2/nu
      zp=.5*h2m*f(ir,1)*f(ir,j)*xdd*x
      zpe=h2m*af(j)*(f(ir,1)**2*xdd-1)*x
      wpcmcc=wpcmcc+zp*gfcc(ir,j)*rsdsl(ir)
      wpcmsc=wpcmsc+zp*gfcc(ir,j)*rsdsl(ir)*(ve**2-1)
      wpcmdd=wpcmdd-zpe*ghdd(ir,j)*rlssx(ir)/nu
      wpcmsd=wpcmsd-zpe*ghdd(ir,j)*rlssx(ir)*(ve**2-1)/nu
      gnn(ir,1)=gnn(ir,1)+4*f(ir,1)*f(ir,j)*xdd
     & *(gfdd(ir,j)*((vd+xde*vp)**2+xee*vp**2)+(gfde(ir,j)
     & +gfed(ir,j))*(vd+xde*vp)*ve+gfcc(ir,j)*xcc*ve**2)
     & -2*af(j)*f(ir,1)**2*xdd*ghdd(ir,j)*xcc**2*ve**2/nu
  682 continue
      do 685 m=5,1,-1
  685 wcmx(j,m+5)=wcmx(j,m+5)-wcmx(j,m)
      wvcm=wvcm+wcmx(j,1)+wcmx(j,2)+wcmx(j,3)+wcmx(j,4)
      wvcms=wvcms+wcmx(j,6)+wcmx(j,7)+wcmx(j,8)+wcmx(j,9)
  690 continue
      do 692 m=1,10
  692 wcmjx(m)=.5*(wcmjx(m)+wcmkx(m))
      do 695 m=5,1,-1
      wcmx(1,m+5)=wcmx(1,m+5)-wcmx(1,m)
      wcmkx(m+5)=wcmkx(m+5)-wcmkx(m)
  695 wcmjx(m+5)=wcmjx(m+5)-wcmjx(m)
      wvcm=wvcm+wcmx(1,1)+wcmx(1,2)+wcmx(1,3)+wcmx(1,4)+wcmx(1,5)
      wvcms=wvcms+wcmx(1,6)+wcmx(1,7)+wcmx(1,8)+wcmx(1,9)+wcmx(1,10)
      wkcm=wcmkx(1)+wcmkx(2)+wcmkx(3)+wcmkx(4)+wcmkx(5)
      wkcms=wcmkx(6)+wcmkx(7)+wcmkx(8)+wcmkx(9)+wcmkx(10)
      wjcm=wcmjx(1)+wcmjx(2)+wcmjx(3)+wcmjx(4)+wcmjx(5)
      wjcms=wcmjx(6)+wcmjx(7)+wcmjx(8)+wcmjx(9)+wcmjx(10)
      wfcm=wfcmcc+wfcmdd
      wfcms=wfcmsc+wfcmsd
      wpcm=wpcmcc+wpcmdd
      wpcms=wpcmsc+wpcmsd
c ------------------
c calculate wcr,wcrs
c ------------------
      wvcr=0
      wfcr=0
      wpcr=0
      do 750 j=1,6,nm
      do 745 i=1,6,nm
      do 740 k=1,i,nm
      x2=acex(i,j,k)
      if (x2.eq.0.) go to 740
      x1=ac(i,j,k)
      qx=4*x
      if (i.eq.k) qx=2*x
      di1dk1=at(i,1)*at(k,1)
      ditdkt=(at(i,5)+at(i,6))*(at(k,5)+at(k,6))
c --------------
c vertex factors
c --------------
      do 730 l=1+nm,6,nm
      dli=ad(l,i)
      dlj=ad(l,j)
      dlk=ad(l,k)
      y1=.25*(dli+dlj+dlk)
      y2=.125*dlj+5*(dli+dlk)/24
      y3=.25*(dli+dlk)
      y4=0
      y5=0
      y6=0
      if (i.eq.1.or.j.eq.1.or.k.eq.1) then
        do 700 n=1,4,nm
  700   y4=y4+ak(l,l,n)*aa(n)*(y1+.125*(ad(i,n)+ad(k,n))*at(j,1))
      end if
      do 705 n=1,4,nm
      do 705 mp=1,6,nm
      y5=y5+.5*(ak(n,i,mp)*ak(j,k,mp)+ak(n,k,mp)*ak(i,j,mp))*aa(mp)
     & *(.25*(dlj+ad(l,n))+.5*ad(l,mp))/x2
  705 y6=y6+.5*(ak(n,i,mp)*ak(j,k,mp)+ak(n,k,mp)*ak(i,j,mp))*aa(mp)
     & *(.125*(dlj+ad(l,n))+.25*ad(l,mp)+y3/3)/x2
      do 710 ir=1,lg
      xdd=gx(ir)
      xde=gy(ir)
      xee=xde**2+gz(ir)
      xcc=gl(ir)**2/nu
      z=f(ir,i)*v(ir,j)*f(ir,k)*qx*rs(ir)
      wcrx(j,1)=wcrx(j,1)+z*x1*(y1*(grfd(ir,l)*(1+2*xde+xee)
     & +grfe(ir,l)*(1+xde)+grfc(ir,l))+y2*(grdd(ir,l)*(1+2*xde+xee)
     & +grde(ir,l)*(1+xde)+grdc(ir,l))+y3*(gred(ir,l)*(1+xde)
     & +gree(ir,l))+y4*(grmd(ir,l)*(1+xde)+grme(ir,l)))*xdd
      wcrx(j,2)=wcrx(j,2)-z*x2*(y5*grfd(ir,l)+y6*grdd(ir,l))*xcc*xdd
      gnn(ir,j)=gnn(ir,j)+f(ir,i)*f(ir,k)*(x1*(y1*(grfd(ir,l)
     & *(1+2*xde+xee)+grfe(ir,l)*(1+xde)+grfc(ir,l))
     & +y2*(grdd(ir,l)*(1+2*xde+xee)+grde(ir,l)*(1+xde)+grde(ir,l))
     & +y3*(gred(ir,l)*(1+xde)+gree(ir,l))+y4*(grmd(ir,l)*(1+xde)
     & +grme(ir,l)))-x2*(y5*grfd(ir,l)+y6*grdd(ir,l))*xcc)*xdd*qx/x
  710 continue
c     -------------------
      if (j.eq.1) then
        do 720 ir=1,lg
        xdd=gx(ir)
        xde=gy(ir)
        xee=xde**2+gz(ir)
        xcc=gl(ir)**2/nu
        zk=-.5*h2m*(f(ir,i)*fds(ir,k)+f(ir,k)*fds(ir,i))*qx*rs(ir)
        wcrkx(1)=wcrkx(1)+zk*x1*(y1*(grfd(ir,l)*(1+2*xde+xee)
     &   +grfe(ir,l)*(1+xde)+grfc(ir,l))+y2*(grdd(ir,l)*(1+2*xde+xee)
     &   +grde(ir,l)*(1+xde)+grdc(ir,l))+y3*(gred(ir,l)*(1+xde)
     &   +gree(ir,l))+y4*(grmd(ir,l)*(1+xde)+grme(ir,l)))*xdd
        wcrkx(2)=wcrkx(2)-zk*x2*(y5*grfd(ir,l)+y6*grdd(ir,l))*xcc*xdd
        zf=h2m*(f(ir,i)*fp(ir,k)+f(ir,k)*fp(ir,i))*slp(ir)*qx*rs(ir)
        wfcr=wfcr+zf*x2*(y5*grfd(ir,l)+y6*grdd(ir,l))*gl(ir)*xdd/nu
        zj=h2m*(fp(ir,i)*fp(ir,k)*rs(ir)+6*f(ir,i)*f(ir,k)*ditdkt)*qx
        wcrjx(1)=wcrjx(1)+zj*x1*(y1*(grfd(ir,l)*(1+2*xde+xee)
     &   +grfe(ir,l)*(1+xde)+grfc(ir,l))+y2*(grdd(ir,l)*(1+2*xde+xee)
     &   +grde(ir,l)*(1+xde)+grdc(ir,l))+y3*(gred(ir,l)*(1+xde)
     &   +gree(ir,l))+y4*(grmd(ir,l)*(1+xde)+grme(ir,l)))*xdd
        wcrjx(2)=wcrjx(2)-zj*x2*(y5*grfd(ir,l)+y6*grdd(ir,l))*xcc*xdd
        wpcr=wpcr-.5*h2m*x2*(f(ir,i)*f(ir,k)*xdd-di1dk1)
     &   *(y5*grfd(ir,l)+y6*grdd(ir,l))*rlssx(ir)*qx/nu
  720   continue
      end if
  730 continue
  740 continue
  745 continue
      wvcr=wvcr+wcrx(j,1)+wcrx(j,2)
  750 continue
      do 760 j=1,2
      wcrjx(j)=.5*(wcrjx(j)+wcrkx(j))
  760 continue
      wkcr=wcrkx(1)+wcrkx(2)
      wjcr=wcrjx(1)+wcrjx(2)
c print ================================================================
      if (no.eq.0) go to 800
      write(nlog,1230) wvcd,wkcd,wfcd,wjcd,wpcd,wvcds,wkcds,wfcds,wjcds
     &,wpcds
      write(nout,1230) wvcd,wkcd,wfcd,wjcd,wpcd,wvcds,wkcds,wfcds,wjcds
     &,wpcds
 1230 format(/4x,'wvcd',4x,'wkcd',4x,'wfcd',4x,'wjcd',4x,'wpcd',3x
     &,'wvcds',3x,'wkcds',3x,'wfcds',3x,'wjcds',3x,'wpcds'/10f8.3)
      write(nout,1240)
 1240 format(/4x,'wcd & wcds:dd,de,ee,ddx,cc')
      write(nout,1133)
      write(nout,1140) ((wcdx(i,j),i=1,8),wcdkx(j),wcdjx(j),j=1,10)
      write(nout,1250) wfcddd,wfcdcc,wpcddd,wpcdcc
 1250 format(/4x,'wfcd:dd(ex),cc  wpcd:dd(ex),cc'/4f8.3)
      write(nlog,1260) wvcm,wkcm,wfcm,wjcm,wpcm,wvcms,wkcms,wfcms,wjcms
     &,wpcms
      write(nout,1260) wvcm,wkcm,wfcm,wjcm,wpcm,wvcms,wkcms,wfcms,wjcms
     &,wpcms
 1260 format(/4x,'wvcm',4x,'wkcm',4x,'wfcm',4x,'wjcm',4x,'wpcm',3x
     &,'wvcms',3x,'wkcms',3x,'wfcms',3x,'wjcms',3x,'wpcms'/10f8.3)
      write(nout,1270)
 1270 format(/4x,'wcm:dd,de,ed,cc,ddx')
      write(nout,1133)
      write(nout,1140) ((wcmx(i,j),i=1,8),wcmkx(j),wcmjx(j),j=1,10)
      write(nout,1275) wfcmdd,wfcmcc,wpcmdd,wpcmcc
 1275 format(/4x,'wfcm:dd(ex),cc  wpcm:dd(ex),cc'/4f8.3)
      write(nlog,1280) wvcr,wkcr,wfcr,wjcr,wpcr
      write(nout,1280) wvcr,wkcr,wfcr,wjcr,wpcr
 1280 format(/4x,'wvcr',4x,'wkcr',4x,'wfcr',4x,'wjcr',4x,'wpcr'/5f8.3)
      write(nout,1290)
 1290 format(/4x,'wcr:d,e')
      write(nout,1133)
      write(nout,1140) ((wcrx(i,j),i=1,8),wcrkx(j),wcrjx(j),j=1,2)
      if (nt.ge.1.and.nt.le.3) then
        d3=dr*float(l3)
        write(nlog,1300) d3,tnia,tnic,tniu,tnix,cut,cut0
        write(nout,1300) d3,tnia,tnic,tniu,tnix,cut,cut0
 1300   format(/4x,'d3',6x,'a',7x,'c',7x,'u',7x,'x',7x,'cut',5x,'cut0'
     &        ,/f8.3,4f8.4,2f8.3)
        write(nlog,1305)
        write(nout,1305)
 1305   format(/4x,'w3v(e2b):df,ef,dg,eg')
        write(nlog,1310)
        write(nout,1310)
 1310   format(4x,'c',7x,'t',7x,'s',7x,'st',6x,'tn',6x,'tnt')
        write(nlog,1150) (w3x(jk,1,1),jk=1,6*4)
        write(nout,1150) (w3x(jk,1,1),jk=1,6*4)
 1150   format(6f8.3)
        if (tnix.ne.0. .or. nt.eq.2 .or. nt.eq.3) then
          write(nlog,1315)
          write(nout,1315)
 1315     format(/4x,'w3v(e2b):df,ef,dg,eg - s-wave')
          write(nlog,1310)
          write(nout,1310)
          write(nlog,1150) (w3x(jk,1,2),jk=1,6*4)
          write(nout,1150) (w3x(jk,1,2),jk=1,6*4)
        end if
      end if
      if (nt.ge.4) then
        write(nlog,1320) tniu,tnia,tnic
        write(nout,1320) tniu,tnia,tnic
 1320   format(/4x,'gam1',4x,'gam2',4x,'gam3',/f8.3,f8.0,f8.1)
      end if
  800 ev6=ev6+wvc+wvcs+wvcd+wvcds+wvcm+wvcms+wvcr
      evb=evb+wvcb
      ek6=ek6+wkc+wkcs+wkcd+wkcds+wkcm+wkcms+wkcr
      ekb=ekb+wkcb
      ef6=ef6+wfc+wfcs+wfcd+wfcds+wfcm+wfcms+wfcr
      ej6=ej6+wjc+wjcs+wjcd+wjcds+wjcm+wjcms+wjcr
      ejb=ejb+wjcb
      ep6=ep6+wpc+wpcs+wpcd+wpcds+wpcm+wpcms+wpcr
      return
      end
