c *id* setpot pot ******************************************************
c subroutine for potential
c subroutine setpot should be called once to initialize some constants,
c obtain values for hbar**2/m, and identify potentials with L.S or L**2
c entry pot is then called once for each position at which the potential
c is desired.
c ----------------------------------------------------------------------
c arguments for setpot
c lpotin: 1=malfliet-tjon V
c (lpot)  2=reid v8
c         3=urbana v14
c         4=argonne v8 (v8' reduction of v14)
c         5=argonne v9 (less simplified v14)
c         6=argonne v14
c
c         7=argonne v18-csbl (experimental cd/csb-large)
c         8=argonne v18-csbs (experimental cd/csb-small)
c         9=argonne v18
c        10=argonne v8'
c        11=argonne v6'
c        12=argonne v4'
c        13=argonne v2'
c        14=argonne v1'
c        15=argonne vx'
c        16=argonne v9'(1D)
c        17=argonne v9'(Da)
c
c        21=argonne v18 v1.9
c        22=argonne v8' v1.9
c        23=argonne v18 v1.7
c        24=argonne v8' v1.7
c
c        26=argonne v18p  (p**2 terms)
c        27=argonne v18pq (p**2 terms & l**2 tensor)
c
c        28=super-soft core(c) v14
c        29=super-soft core(c) v8' modified
c        30=paris
c xmn:    parameter for nucleon mass variation
c gam:    parameter for pion mass variation in OPE (Argonne only)
c rho:    parameter for pion mass variation in TPE-s (Argonne only)
c chi:    parameter for pion mass variation in TPE-L (Argonne only)
c omg:    parameter for heavy-meson mass variation (Argonne only)
c ftp:    parameter for intermediate coupling variation (Argonne only)
c h2m:    returned value of .5*hbar**2(1/mp+1/mn)
c h2mcsb: returned value of .5*hbar**2*(1/mp-1/mn)
c ----------------------------------------------------------------------
c arguments for pot
c lr:   potential returned is v(r)*r**lr
c rr:   relative position r in fm
c vv:   potential in MeV (22 component array)
c vp:   p**2 potential terms (12 component array)
c vw:   subsidiary potential terms in MeV (10 component array)
c ----------------------------------------------------------------------
c order of operators l in vv(l):
c l:    1=1                              2=t1.t2
c       3=s1.s2                          4=(s1.s2)(t1.t2)
c       5=S12 [=3(s1.r)(s2.r)-s1.s2]     6=S12(t1.t2)
c       7=L.S [=L.(s1+s2)/2]             8=L.S(t1.t2)
c       9=L**2                          10=L**2(t1.t2)
c      11=L**2(s1.s2)                   12=L**2(s1.s2)(t1.t2)
c      13=(L.S)**2                      14=(L.S)**2(t1.t2)
c      15=T12 [=3*t1z*t2z-t1.t2]        16=(s1.s2)T12
c      17=S12*T12                       18=L.S*T12
c      19=t1z+t2z                       20=(s1.s2)(t1z+t2z)
c      21=S12(t1z+t2z)                  22=L.S(t1z+t2z)
c where s1=sigma_1, t1=tau_1, t1z=tau_1(z), etc.
c for Paris, p**2 terms appear in place of L**2
c ----------------------------------------------------------------------
c order of operators l in vw(l):
c l:    1=Coulomb with form factor 1-exp(-x)*(1+x*(33+x*(9+x))/48))
c       2=OPEP (s1.s2)(t1.t2)            3=OPEP S12(t1.t2)
c       4=OPEP (s1.s2)T12                5=OPEP S12*T12
c       6=Del2 (t1z+t2z)                 7=Del2 L.A(t1z-t2z)
c       8=OPEP L.(s1xs2)(t1xt2)z         9=OREP L.A(t1z-t2z)
c      10=OREP L.(s1xs2)(t1xt2)z                
c ----------------------------------------------------------------------
      subroutine setpot(lpotin,xmn,gam,rho,chi,omg,ftp,h2m,h2mcsb)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      logical lpotls,lpotll
      common /logpot/ lpotls,lpotll
      logical llpotls(30),llpotll(30)
      dimension vv(22),vp(12),vw(10)
      real*8 mpi0,mpic,mpi,mpis,mpiq,mp,mn,mu0,muc,mu,mus,muq
      real*8 krho,komg,lam,lamp,lamw,mrho,momg,murho,muomg
      dimension pm0(12),pm1(12)
     &,pgva01(12),pgva11(12),pgva00(12),pgva10(12)
     &,pgvb01(12),pgvb11(12),pgvb00(12),pgvb10(12)
     &,pgvls1(12),pgvt1(12),pgvso21(12)
     &,pgvls0(12),pgvt0(12),pgvso20(12)
      common /mypot/ lpot,ftpec,pimass1,pimass2,pimass3,
     & wsrange,tnr,mp,mn
c      save lpot,ftpec,pimass1,pimass2,pimass3,wsrange,tnr
c      save mp,mn
c -------------------
c statement functions
c -------------------
      yc(t)=exp(-t)/x
      yt(t)=(1+3/t+3/t**2)*exp(-t)/x
      yls(t)=-(1+t)*exp(-t)/x**3
      yl2(t)=(1+2/t)*exp(-t)/x**3
      pc(t)=exp(-t)/t
      pt(t)=(1+3/t+3/t**2)
      pls(t)=(1+1/t)/t
      pso2(t)=(1+3/t+3/t**2)/t**2

      llpotls(1) =.false.
      llpotls(2) =.true.
      llpotls(3) =.true.
      llpotls(4) =.true.
      llpotls(5) =.true.
      llpotls(6) =.true.
      llpotls(7) =.true.
      llpotls(8) =.true.
      llpotls(9) =.true.
      llpotls(10) =.true.
      llpotls(11) =.false.
      llpotls(12) =.false.
      llpotls(13) =.false.
      llpotls(14) =.false.
      llpotls(15) =.false.
      llpotls(16) =.true.
      llpotls(17) =.true.
      llpotls(18) =.true.
      llpotls(19) =.true.
      llpotls(20) =.true.
      llpotls(21) =.false.
      llpotls(22) =.false.
      llpotls(23) =.false.
      llpotls(24) =.false.
      llpotls(25) =.false.
      llpotls(26) =.true.
      llpotls(27) =.true.
      llpotls(28) =.true.
      llpotls(29) =.true.
      llpotls(30) =.true.

      llpotll(1) =.false.
      llpotll(2) =.false.
      llpotll(3) =.true.
      llpotll(4) =.false.
      llpotll(5) =.true.
      llpotll(6) =.true.
      llpotll(7) =.true.
      llpotll(8) =.true.
      llpotll(9) =.true.
      llpotll(10) =.false.
      llpotll(11) =.false.
      llpotll(12) =.false.
      llpotll(13) =.false.
      llpotll(14) =.false.
      llpotll(15) =.false.
      llpotll(16) =.true.
      llpotll(17) =.true.
      llpotll(18) =.true.
      llpotll(19) =.true.
      llpotll(20) =.true.
      llpotll(21) =.true.
      llpotll(22) =.false.
      llpotll(23) =.true.
      llpotll(24) =.false.
      llpotll(25) =.true.
      llpotll(26) =.true.
      llpotll(27) =.true.
      llpotll(28) =.true.
      llpotll(29) =.false.
      llpotll(30) =.true.
      small=1e-4
      vsmall=1e-10
c -------------
c set hbar**2/m
c -------------
      h2m=41.47
      h2mcsb=0.
      hc=197.327053
      mp=938.27231*xmn
      mn=939.56563*xmn
      lpot=lpotin
      if (lpot.ge.4 .and. lpot.le.6) h2m=197.33**2/(938.9*xmn)
      if (lpot.ge.7 .and. lpot.le.27) then
        h2m=.5*hc**2*(1./mp+1./mn)
        h2mcsb=.5*hc**2*(1./mp-1./mn)
      end if
c ----------------------------------------------------------------------
c identify potentials that have L.S and/or terms quadratic in L or p
c ----------------------------------------------------------------------
      lpotls=llpotls(lpot)
      lpotll=llpotll(lpot)
c --------------------
c initialize LFP model
c --------------------
      tnr=1.
      if (lpot.eq.3) tnr=exp(-gam*rho)
c ------------------------------
c initialize pion mass variation
c ------------------------------
      if (lpot.ge.4 .and. lpot.le.27) then
        pimass1=gam
        pimass2=rho
        pimass3=chi
        wsrange=1.-(2./3.)*(omg-1.)
        ftpec=ftp
      end if
c --------------------
c paris initialization
c --------------------
      pm1(1)=.684026
      pm1(2)=1.6
      pm1(3)=2.3
      pm1(4)=3.0
      pm1(5)=3.7
      pm1(6)=4.4
      pm1(7)=5.1
      pm1(8)=5.8
      pm1(9)=6.5
      pm1(10)=8.2
      pm1(11)=9.9
      pm1(12)=11.3

      pm0(1)=.699536
      pm0(2)=1.6
      pm0(3)=2.3
      pm0(4)=3.0
      pm0(5)=3.7
      pm0(6)=4.4
      pm0(7)=5.1
      pm0(8)=5.8
      pm0(9)=6.5
      pm0(10)=8.2
      pm0(11)=9.9
      pm0(12)=11.3
      pgva01(1)=-10.077427
      pgva01(2)=-120.49564
      pgva01(3)=-212.36460
      pgva01(4)=-8717.4198
      pgva01(5)= 54383.377
      pgva01(6)=-213421.47
      pgva01(7)= 494583.57
      pgva01(8)=-667153.34
      pgva01(9)= 529575.98
      pgva01(10)=-137034.12
      pgva01(11)=-346971.94
      pgva01(12)=         0

      pgva11(1)= 3.3591422
      pgva11(2)=-86.479568
      pgva11(3)=-465.93111
      pgva11(4)= 1867.3085
      pgva11(5)= 3850.9213
      pgva11(6)=-19674.338
      pgva11(7)= 123231.40
      pgva11(8)=-314493.61
      pgva11(9)= 242424.40
      pgva11(10)= 166904.04
      pgva11(11)=-485343.64
      pgva11(12)=         0

      pgvb01(1)= .0026851393
      pgvb01(2)=.051092455
      pgvb01(3)=-.84264258
      pgvb01(4)=14.736312
      pgvb01(5)=-145.21993
      pgvb01(6)= 841.58389
      pgvb01(7)=-2786.1170
      pgvb01(8)= 5056.4510
      pgvb01(9)=-3367.4205
      pgvb01(10)=-1784.5529
      pgvb01(11)= 5354.8266
      pgvb01(12)=         0

      pgvb11(1)=-.00089504644
      pgvb11(2)=.037488481
      pgvb11(3)=-.89373089
      pgvb11(4)=14.123475
      pgvb11(5)=-146.60152
      pgvb11(6)= 841.91462
      pgvb11(7)=-2839.4273
      pgvb11(8)= 5265.3427
      pgvb11(9)=-3500.0430
      pgvb11(10)=-2487.9479
      pgvb11(11)= 7306.8121
      pgvb11(12)=         0

      pgvls1(1)=         0
      pgvls1(2)=-426.00359
      pgvls1(3)= 26279.517
      pgvls1(4)=-575570.33
      pgvls1(5)= 6003393.4
      pgvls1(6)=-34519443.
      pgvls1(7)=113554590.
      pgvls1(8)=-207292090.
      pgvls1(9)=171315480.
      pgvls1(10)=-86418222.
      pgvls1(11)=         0
      pgvls1(12)=         0

      pgvt1(1)= 3.3591422
      pgvt1(2)=-.85945824
      pgvt1(3)=-104.76340
      pgvt1(4)= 1262.9465
      pgvt1(5)=-18881.061
      pgvt1(6)= 106132.46
      pgvt1(7)=-332119.10
      pgvt1(8)= 555857.62
      pgvt1(9)=-349166.64
      pgvt1(10)=-119450.13
      pgvt1(11)=         0
      pgvt1(12)=         0

      pgvso21(1)=        0
      pgvso21(2)=-.52218640
      pgvso21(3)= 186.44558
      pgvso21(4)=-3709.1115
      pgvso21(5)= 55913.117
      pgvso21(6)=-369985.60
      pgvso21(7)= 1453754.3
      pgvso21(8)=-3135247.1
      pgvso21(9)= 2433908.1
      pgvso21(10)=         0
      pgvso21(11)=         0
      pgvso21(12)=         0

      pgva00(1)= 32.290874
      pgva00(2)=-82.465631
      pgva00(3)= 1232.9384
      pgva00(4)=-16859.879
      pgva00(5)= 172926.83
      pgva00(6)=-768352.77
      pgva00(7)= 2189047.5
      pgva00(8)=-3844728.7
      pgva00(9)= 2799055.9
      pgva00(10)= 502518.28
      pgva00(11)=-2600612.4
      pgva00(12)=         0

      pgva10(1)=-10.763625
      pgva10(2)=-42.973669
      pgva10(3)=-718.56844
      pgva10(4)= 4246.9120
      pgva10(5)=-34574.024
      pgva10(6)= 126711.69
      pgva10(7)=-274168.41
      pgva10(8)= 529607.24
      pgva10(9)=-366067.13
      pgva10(10)=-223036.73
      pgva10(11)= 406838.33
      pgva10(12)=         0

      pgvb00(1)=-.0085980096
      pgvb00(2)=.026814385
      pgvb00(3)=-1.3280693
      pgvb00(4)=10.324289
      pgvb00(5)=-115.27067
      pgvb00(6)= 694.56175
      pgvb00(7)=-2387.9335
      pgvb00(8)= 4238.8011
      pgvb00(9)=-2452.1604
      pgvb00(10)=-1951.2821
      pgvb00(11)= 4180.1160
      pgvb00(12)=         0

      pgvb10(1)= .0028660032
      pgvb10(2)=-.00081798046
      pgvb10(3)=-.53314560
      pgvb10(4)=.83162030
      pgvb10(5)=-31.192395
      pgvb10(6)= 300.41384
      pgvb10(7)=-1241.5067
      pgvb10(8)= 2476.2241
      pgvb10(9)=-1304.3030
      pgvb10(10)=-2149.6577
      pgvb10(11)= 4099.6917
      pgvb10(12)=         0

      pgvls0(1)=         0
      pgvls0(2)=-66.176421
      pgvls0(3)= 2890.3688
      pgvls0(4)=-62592.400
      pgvls0(5)= 691461.41
      pgvls0(6)=-4096914.6
      pgvls0(7)= 14032093.
      pgvls0(8)=-26827468.
      pgvls0(9)= 23511442.
      pgvls0(10)=-14688461.
      pgvls0(11)=         0
      pgvls0(12)=         0

      pgvt0(1)=-10.763625
      pgvt0(2)=-.46818029
      pgvt0(3)= 60.147739
      pgvt0(4)= 352.56941
      pgvt0(5)= 514.32170
      pgvt0(6)= 11637.302
      pgvt0(7)=-44595.415
      pgvt0(8)= 69211.738
      pgvt0(9)=-48127.668
      pgvt0(10)= 7051.4008
      pgvt0(11)=         0
      pgvt0(12)=         0

      pgvso20(1)=        0
      pgvso20(2)=-.62851020
      pgvso20(3)=-76.290197
      pgvso20(4)=-788.27581
      pgvso20(5)=-6490.4798
      pgvso20(6)= 5473.4378
      pgvso20(7)=-32941.912
      pgvso20(8)= 249491.32
      pgvso20(9)=-16012.956
      pgvso20(10)=         0
      pgvso20(11)=         0
      pgvso20(12)=         0
      if (lpot.eq.30) then
      suma01=0
      suma11=0
      suma00=0
      suma10=0
      sumb01=0
      sumb11=0
      sumb00=0
      sumb10=0
      do 1 i=1,11
        suma01=suma01+pgva01(i)/pm1(i)
        suma11=suma11+pgva11(i)/pm1(i)
        suma00=suma00+pgva00(i)/pm0(i)
        suma10=suma10+pgva10(i)/pm0(i)
        sumb01=sumb01+pgvb01(i)/pm1(i)
        sumb11=sumb11+pgvb11(i)/pm1(i)
        sumb00=sumb00+pgvb00(i)/pm0(i)
        sumb10=sumb10+pgvb10(i)/pm0(i)
    1 continue
      pgva01(12)=-suma01*pm1(12)
      pgva11(12)=-suma11*pm1(12)
      pgva00(12)=-suma00*pm0(12)
      pgva10(12)=-suma10*pm0(12)
      pgvb01(12)=-sumb01*pm1(12)
      pgvb11(12)=-sumb11*pm1(12)
      pgvb00(12)=-sumb00*pm0(12)
      pgvb10(12)=-sumb10*pm0(12)
      sumt1=0
      sumt0=0
      sumls1=0
      sumls0=0
      tumt1=0
      tumt0=0
      tumls1=0
      tumls0=0
      do 2 i=1,10
        sumt1=sumt1+pgvt1(i)/pm1(i)
        sumt0=sumt0+pgvt0(i)/pm0(i)
        sumls1=sumls1+pgvls1(i)/pm1(i)
        sumls0=sumls0+pgvls0(i)/pm0(i)
        tumt1=tumt1+pgvt1(i)/pm1(i)**3
        tumt0=tumt0+pgvt0(i)/pm0(i)**3
        tumls1=tumls1+pgvls1(i)/pm1(i)**3
        tumls0=tumls0+pgvls0(i)/pm0(i)**3
    2 continue
      pgvt1(11)=pm1(11)**3*(pm1(12)**2*tumt1-sumt1)
     &          /(pm1(11)**2-pm1(12)**2)
      pgvt1(12)=pm1(12)**3*(pm1(11)**2*tumt1-sumt1)
     &          /(pm1(12)**2-pm1(11)**2)
      pgvt0(11)=pm0(11)**3*(pm0(12)**2*tumt0-sumt0)
     &          /(pm0(11)**2-pm0(12)**2)
      pgvt0(12)=pm0(12)**3*(pm0(11)**2*tumt0-sumt0)
     &          /(pm0(12)**2-pm0(11)**2)
      pgvls1(11)=pm1(11)**3*(pm1(12)**2*tumls1-sumls1)
     &          /(pm1(11)**2-pm1(12)**2)
      pgvls1(12)=pm1(12)**3*(pm1(11)**2*tumls1-sumls1)
     &          /(pm1(12)**2-pm1(11)**2)
      pgvls0(11)=pm0(11)**3*(pm0(12)**2*tumls0-sumls0)
     &          /(pm0(11)**2-pm0(12)**2)
      pgvls0(12)=pm0(12)**3*(pm0(11)**2*tumls0-sumls0)
     &          /(pm0(12)**2-pm0(11)**2)
      sumso21=0
      sumso20=0
      tumso21=0
      tumso20=0
      uumso21=0
      uumso20=0
      do 3 i=1,9
        sumso21=sumso21+pgvso21(i)/pm1(i)
        sumso20=sumso20+pgvso20(i)/pm0(i)
        tumso21=tumso21+pgvso21(i)/pm1(i)**3
        tumso20=tumso20+pgvso20(i)/pm0(i)**3
        uumso21=uumso21+pgvso21(i)/pm1(i)**5
        uumso20=uumso20+pgvso20(i)/pm0(i)**5
    3 continue
      pgvso21(10)=pm1(10)**5*(-pm1(11)**2*pm1(12)**2*uumso21
     &                        +(pm1(11)**2+pm1(12)**2)*tumso21-sumso21)
     &           /((pm1(12)**2-pm1(10)**2)*(pm1(11)**2-pm1(10)**2))
      pgvso21(11)=pm1(11)**5*(-pm1(12)**2*pm1(10)**2*uumso21
     &                        +(pm1(12)**2+pm1(10)**2)*tumso21-sumso21)
     &           /((pm1(10)**2-pm1(11)**2)*(pm1(12)**2-pm1(11)**2))
      pgvso21(12)=pm1(12)**5*(-pm1(10)**2*pm1(11)**2*uumso21
     &                        +(pm1(10)**2+pm1(11)**2)*tumso21-sumso21)
     &           /((pm1(11)**2-pm1(12)**2)*(pm1(10)**2-pm1(12)**2))
      pgvso20(10)=pm0(10)**5*(-pm0(11)**2*pm0(12)**2*uumso20
     &                        +(pm0(11)**2+pm0(12)**2)*tumso20-sumso20)
     &           /((pm0(12)**2-pm0(10)**2)*(pm0(11)**2-pm0(10)**2))
      pgvso20(11)=pm0(11)**5*(-pm0(12)**2*pm0(10)**2*uumso20
     &                        +(pm0(12)**2+pm0(10)**2)*tumso20-sumso20)
     &           /((pm0(10)**2-pm0(11)**2)*(pm0(12)**2-pm0(11)**2))
      pgvso20(12)=pm0(12)**5*(-pm0(10)**2*pm0(11)**2*uumso20
     &                        +(pm0(10)**2+pm0(11)**2)*tumso20-sumso20)
     &           /((pm0(11)**2-pm0(12)**2)*(pm0(10)**2-pm0(12)**2))
      end if
      return
      end subroutine
c **********************************************************************
c      entry pot(lr,rr,vv,vp,vw)
c      subroutine setpot(lpotin,xmn,gam,rho,chi,omg,ftp,h2m,h2mcsb)
      subroutine pot(lr,rr,vv,vp,vw)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c      logical lpotls,lpotll
      common /logpot/ lpotls,lpotll
      logical llpotls(30),llpotll(30)
      dimension vv(22),vp(12),vw(10)
      real*8 mpi0,mpic,mpi,mpis,mpiq,mp,mn,mu0,muc,mu,mus,muq
      real*8 krho,komg,lam,lamp,lamw,mrho,momg,murho,muomg
      dimension pm0(12),pm1(12)
     &,pgva01(12),pgva11(12),pgva00(12),pgva10(12)
     &,pgvb01(12),pgvb11(12),pgvb00(12),pgvb10(12)
     &,pgvls1(12),pgvt1(12),pgvso21(12)
     &,pgvls0(12),pgvt0(12),pgvso20(12)
      common /mypot/ lpot,ftpec,pimass1,pimass2,pimass3,
     & wsrange,tnr,mp,mn
c      save lpot,ftpec,pimass1,pimass2,pimass3,wsrange,tnr
c      save mp,mn
c -------------------
c statement functions
c -------------------
      yc(t)=exp(-t)/x
      yt(t)=(1+3/t+3/t**2)*exp(-t)/x
      yls(t)=-(1+t)*exp(-t)/x**3
      yl2(t)=(1+2/t)*exp(-t)/x**3
      pc(t)=exp(-t)/t
      pt(t)=(1+3/t+3/t**2)
      pls(t)=(1+1/t)/t
      pso2(t)=(1+3/t+3/t**2)/t**2
c --------------------
c paris initialization
c --------------------
      llpotls(1) =.false.
      llpotls(2) =.true.
      llpotls(3) =.true.
      llpotls(4) =.true.
      llpotls(5) =.true.
      llpotls(6) =.true.
      llpotls(7) =.true.
      llpotls(8) =.true.
      llpotls(9) =.true.
      llpotls(10) =.true.
      llpotls(11) =.false.
      llpotls(12) =.false.
      llpotls(13) =.false.
      llpotls(14) =.false.
      llpotls(15) =.false.
      llpotls(16) =.true.
      llpotls(17) =.true.
      llpotls(18) =.true.
      llpotls(19) =.true.
      llpotls(20) =.true.
      llpotls(21) =.false.
      llpotls(22) =.false.
      llpotls(23) =.false.
      llpotls(24) =.false.
      llpotls(25) =.false.
      llpotls(26) =.true.
      llpotls(27) =.true.
      llpotls(28) =.true.
      llpotls(29) =.true.
      llpotls(30) =.true.

      llpotll(1) =.false.
      llpotll(2) =.false.
      llpotll(3) =.true.
      llpotll(4) =.false.
      llpotll(5) =.true.
      llpotll(6) =.true.
      llpotll(7) =.true.
      llpotll(8) =.true.
      llpotll(9) =.true.
      llpotll(10) =.false.
      llpotll(11) =.false.
      llpotll(12) =.false.
      llpotll(13) =.false.
      llpotll(14) =.false.
      llpotll(15) =.false.
      llpotll(16) =.true.
      llpotll(17) =.true.
      llpotll(18) =.true.
      llpotll(19) =.true.
      llpotll(20) =.true.
      llpotll(21) =.true.
      llpotll(22) =.false.
      llpotll(23) =.true.
      llpotll(24) =.false.
      llpotll(25) =.true.
      llpotll(26) =.true.
      llpotll(27) =.true.
      llpotll(28) =.true.
      llpotll(29) =.false.
      llpotll(30) =.true.
      small=1e-4
      vsmall=1e-10

      pm1(1)=.684026
      pm1(2)=1.6
      pm1(3)=2.3
      pm1(4)=3.0
      pm1(5)=3.7
      pm1(6)=4.4
      pm1(7)=5.1
      pm1(8)=5.8
      pm1(9)=6.5
      pm1(10)=8.2
      pm1(11)=9.9
      pm1(12)=11.3

      pm0(1)=.699536
      pm0(2)=1.6
      pm0(3)=2.3
      pm0(4)=3.0
      pm0(5)=3.7
      pm0(6)=4.4
      pm0(7)=5.1
      pm0(8)=5.8
      pm0(9)=6.5
      pm0(10)=8.2
      pm0(11)=9.9
      pm0(12)=11.3
      pgva01(1)=-10.077427
      pgva01(2)=-120.49564
      pgva01(3)=-212.36460
      pgva01(4)=-8717.4198
      pgva01(5)= 54383.377
      pgva01(6)=-213421.47
      pgva01(7)= 494583.57
      pgva01(8)=-667153.34
      pgva01(9)= 529575.98
      pgva01(10)=-137034.12
      pgva01(11)=-346971.94
      pgva01(12)=         0

      pgva11(1)= 3.3591422
      pgva11(2)=-86.479568
      pgva11(3)=-465.93111
      pgva11(4)= 1867.3085
      pgva11(5)= 3850.9213
      pgva11(6)=-19674.338
      pgva11(7)= 123231.40
      pgva11(8)=-314493.61
      pgva11(9)= 242424.40
      pgva11(10)= 166904.04
      pgva11(11)=-485343.64
      pgva11(12)=         0

      pgvb01(1)= .0026851393
      pgvb01(2)=.051092455
      pgvb01(3)=-.84264258
      pgvb01(4)=14.736312
      pgvb01(5)=-145.21993
      pgvb01(6)= 841.58389
      pgvb01(7)=-2786.1170
      pgvb01(8)= 5056.4510
      pgvb01(9)=-3367.4205
      pgvb01(10)=-1784.5529
      pgvb01(11)= 5354.8266
      pgvb01(12)=         0

      pgvb11(1)=-.00089504644
      pgvb11(2)=.037488481
      pgvb11(3)=-.89373089
      pgvb11(4)=14.123475
      pgvb11(5)=-146.60152
      pgvb11(6)= 841.91462
      pgvb11(7)=-2839.4273
      pgvb11(8)= 5265.3427
      pgvb11(9)=-3500.0430
      pgvb11(10)=-2487.9479
      pgvb11(11)= 7306.8121
      pgvb11(12)=         0

      pgvls1(1)=         0
      pgvls1(2)=-426.00359
      pgvls1(3)= 26279.517
      pgvls1(4)=-575570.33
      pgvls1(5)= 6003393.4
      pgvls1(6)=-34519443.
      pgvls1(7)=113554590.
      pgvls1(8)=-207292090.
      pgvls1(9)=171315480.
      pgvls1(10)=-86418222.
      pgvls1(11)=         0
      pgvls1(12)=         0

      pgvt1(1)= 3.3591422
      pgvt1(2)=-.85945824
      pgvt1(3)=-104.76340
      pgvt1(4)= 1262.9465
      pgvt1(5)=-18881.061
      pgvt1(6)= 106132.46
      pgvt1(7)=-332119.10
      pgvt1(8)= 555857.62
      pgvt1(9)=-349166.64
      pgvt1(10)=-119450.13
      pgvt1(11)=         0
      pgvt1(12)=         0

      pgvso21(1)=        0
      pgvso21(2)=-.52218640
      pgvso21(3)= 186.44558
      pgvso21(4)=-3709.1115
      pgvso21(5)= 55913.117
      pgvso21(6)=-369985.60
      pgvso21(7)= 1453754.3
      pgvso21(8)=-3135247.1
      pgvso21(9)= 2433908.1
      pgvso21(10)=         0
      pgvso21(11)=         0
      pgvso21(12)=         0

      pgva00(1)= 32.290874
      pgva00(2)=-82.465631
      pgva00(3)= 1232.9384
      pgva00(4)=-16859.879
      pgva00(5)= 172926.83
      pgva00(6)=-768352.77
      pgva00(7)= 2189047.5
      pgva00(8)=-3844728.7
      pgva00(9)= 2799055.9
      pgva00(10)= 502518.28
      pgva00(11)=-2600612.4
      pgva00(12)=         0

      pgva10(1)=-10.763625
      pgva10(2)=-42.973669
      pgva10(3)=-718.56844
      pgva10(4)= 4246.9120
      pgva10(5)=-34574.024
      pgva10(6)= 126711.69
      pgva10(7)=-274168.41
      pgva10(8)= 529607.24
      pgva10(9)=-366067.13
      pgva10(10)=-223036.73
      pgva10(11)= 406838.33
      pgva10(12)=         0

      pgvb00(1)=-.0085980096
      pgvb00(2)=.026814385
      pgvb00(3)=-1.3280693
      pgvb00(4)=10.324289
      pgvb00(5)=-115.27067
      pgvb00(6)= 694.56175
      pgvb00(7)=-2387.9335
      pgvb00(8)= 4238.8011
      pgvb00(9)=-2452.1604
      pgvb00(10)=-1951.2821
      pgvb00(11)= 4180.1160
      pgvb00(12)=         0

      pgvb10(1)= .0028660032
      pgvb10(2)=-.00081798046
      pgvb10(3)=-.53314560
      pgvb10(4)=.83162030
      pgvb10(5)=-31.192395
      pgvb10(6)= 300.41384
      pgvb10(7)=-1241.5067
      pgvb10(8)= 2476.2241
      pgvb10(9)=-1304.3030
      pgvb10(10)=-2149.6577
      pgvb10(11)= 4099.6917
      pgvb10(12)=         0

      pgvls0(1)=         0
      pgvls0(2)=-66.176421
      pgvls0(3)= 2890.3688
      pgvls0(4)=-62592.400
      pgvls0(5)= 691461.41
      pgvls0(6)=-4096914.6
      pgvls0(7)= 14032093.
      pgvls0(8)=-26827468.
      pgvls0(9)= 23511442.
      pgvls0(10)=-14688461.
      pgvls0(11)=         0
      pgvls0(12)=         0

      pgvt0(1)=-10.763625
      pgvt0(2)=-.46818029
      pgvt0(3)= 60.147739
      pgvt0(4)= 352.56941
      pgvt0(5)= 514.32170
      pgvt0(6)= 11637.302
      pgvt0(7)=-44595.415
      pgvt0(8)= 69211.738
      pgvt0(9)=-48127.668
      pgvt0(10)= 7051.4008
      pgvt0(11)=         0
      pgvt0(12)=         0

      pgvso20(1)=        0
      pgvso20(2)=-.62851020
      pgvso20(3)=-76.290197
      pgvso20(4)=-788.27581
      pgvso20(5)=-6490.4798
      pgvso20(6)= 5473.4378
      pgvso20(7)=-32941.912
      pgvso20(8)= 249491.32
      pgvso20(9)=-16012.956
      pgvso20(10)=         0
      pgvso20(11)=         0
      pgvso20(12)=         0
      rrsave=rr
      if (rr.le.small) rr=small
      vv(:)=0
      vp(:)=0
      vw(:)=0
      x=4.27*rr
      ff=1-exp(-x)*(1+x*(33+x*(9+x))/48)

      vw(1)=1.4399652*ff/rr
      rr=rrsave
      go to (10,20,30,40,40,40,50,50,50,50
     &      ,50,50,50,50,50,50,50,50,50,50
     &      ,50,50,50,50,50,60,60,70,70,80), lpot
c -------------
c malfliet-tjon
c -------------
   10 if (rr.le.vsmall) rr=vsmall
      hc=197.3
      p01=(7.39*exp(-3.11*rr)-2.64*exp(-1.55*rr))*hc/rr
      p10=(7.39*exp(-3.11*rr)-3.22*exp(-1.55*rr))*hc/rr
      p00=0
      p11=0
c ------------------
c V version of force
c ------------------
      vv(1)=.5*(p01+p10)
c ----------------------
c I-III version of force
c ----------------------
c     vv(1)=.25*(3*p01+p10)
c     vv(2)=.25*(  p01-p10)
c -------------------
c v4 version of force
c -------------------
c     vv(1)=.0625*(9*p11+3*p10+3*p01+p00)
c     vv(2)=.0625*(3*p11-3*p10  +p01-p00)
c     vv(3)=.0625*(3*p11  +p10-3*p01-p00)
c     vv(4)=.0625*(  p11  -p10  -p01+p00)
      rr=rrsave
      go to 200
c -------
c reid v8
c -------
   20 if (rr.le.vsmall) rr=vsmall
      u=.7
      x=u*rr
      y1=yc(x)
      y2=yc(2*x)
      y3=yc(3*x)
      y4=yc(4*x)
      y6=yc(6*x)
      y7=yc(7*x)
      yr=yt(x)-(12/x+3/x**2)*y4
      if (rr.le.10*vsmall) yr=23.5/x
      hr=10.463
      vv(1)=-19.874*y2+135.21*y3-1432.3*y4+4196.4*y6+1215.8*y7
      vv(2)=19.874*y2-135.21*y3+319.52*y4-1082.3*y6+405.3*y7
      vv(3)=46.241*y2-135.21*y3-64.78*y4+1398.8*y6-1215.8*y7
      vv(4)=(hr/3)*y1-46.241*y2+135.21*y3+244.06*y4-360.76*y6-405.3*y7
      vv(5)=-26.194*y3+87.943*y4-418.38*y6
      vv(6)=(hr/3)*yr-8.731*y3-87.943*y4+418.38*y6
      vv(7)=177.23*y4-2233.9*y6
      vv(8)=-177.23*y4+159.75*y6
      vw(2)=(hr/3)*y1
      vw(3)=(hr/3)*yr
      rr=rrsave
      go to 200
c ----------
c urbana v14
c ----------
   30 u=.7
      cpi=2
      x=u*rr
      if (rr.le.small) then
        ypi=cpi*rr/u
        tpi=3*cpi**2*rr/u**3
      else
        rcut=1-exp(-cpi*rr*rr)
        ypi=yc(x)*rcut
        tpi=yt(x)*rcut**2
      end if
      tpi2=tpi*tpi*tnr
      ypi=10.463*ypi/3
      tpi=10.463*tpi/3
      ws=1/(1+exp((rr-.5)/.2))
      wp=1/(1+exp((rr-.36)/.17))
      p11=  -4.32  *tpi2+2145.*ws+  ypi
      pt1=  -0.18  *tpi2         +  tpi
      pls1=             -2200.*wp
      pl211=            -  20.*ws
      pls21=             147.5*ws
c following is old urbana v14 (incorrect deuteron)
c     p10=  -6.8009*tpi2+2400.*ws-3*ypi
c following is new urbana v14 (correct deuteron)
      p10=  -6.7983*tpi2+2400.*ws-3*ypi
      pt0=   0.75  *tpi2         -3*tpi
      pls0=                80.*ws
      pl210=              380.*ws
      pls20=-0.2   *tpi2- 230.*ws
      p01=  -6.255 *tpi2+2000.*ws-3*ypi
      pl201=               49.*ws
      p00= -13.2   *tpi2+8700.*ws+9*ypi
      pl200= 0.6   *tpi2- 500.*ws
      vv(1)=.0625*(9*p11+3*p10+3*p01+p00)
      vv(2)=.0625*(3*p11-3*p10  +p01-p00)
      vv(3)=.0625*(3*p11  +p10-3*p01-p00)
      vv(4)=.0625*(  p11  -p10  -p01+p00)
      vv(5)=.25*(3*pt1+pt0)
      vv(6)=.25*(  pt1-pt0)
      vv(7)=.25*(3*pls1+pls0)
      vv(8)=.25*(  pls1-pls0)
      vv(9)= .0625*(9*pl211+3*pl210+3*pl201+pl200)
      vv(10)=.0625*(3*pl211-3*pl210+  pl201-pl200)
      vv(11)=.0625*(3*pl211+  pl210-3*pl201-pl200)
      vv(12)=.0625*(  pl211-  pl210-  pl201+pl200)
      vv(13)=.25*(3*pls21+pls20)
      vv(14)=.25*(  pls21-pls20)
      vw(2)=ypi
      vw(3)=tpi
      go to 200
c ---------------------------
c argonne v14 and derivatives
c lpot = 4 -> v8'
c        5 -> v9'
c        6 -> v14
c ---------------------------
   40 u=138.03/197.33
      u1=pimass1*u
      u2=pimass2*u
      u3=pimass3*u
      cpi=2
      x1=u1*rr
      x2=u2*rr
      x3=u3*rr
      if (rr.le.small) then
        ypi=cpi*rr/u1
        tpi=3*cpi**2*rr/u1**3
        ypib=cpi*rr/u2
        tpib=3*cpi**2*rr/u2**3
        ypid=cpi*rr/u3
        tpid=3*cpi**2*rr/u3**3
      else
        rcut=1-exp(-cpi*rr*rr)
        ypi=rcut*exp(-x1)/x1
        tpi=rcut*(1+3/x1+3/x1**2)*ypi
        ypib=rcut*exp(-x2)/x2
        tpib=rcut*(1+3/x2+3/x2**2)*ypib
        ypid=rcut*exp(-x3)/x3
        tpid=rcut*(1+3/x3+3/x3**2)*ypid
      end if
      pifac=pimass1**3
      ypi=pifac*3.72681*ypi
      tpi=pifac*3.72681*tpi
      tpi2=ftpec*(pimass2**3*tpib)**2
      tpi3=ftpec*(pimass3**3*tpid)**2
      rws=wsrange*.5
      aws=wsrange*.2
      ws=1./(1+exp((rr-rws)/aws))
      p11=  -2.63  *tpi2+1179.*ws+  ypi
      pt1=  -0.91  *tpi2+ 406.*ws+  tpi
      pls1=  0.61  *tpi3- 879.*ws
      pl211=-0.12  *tpi3-   2.*ws
      pls21=-0.54  *tpi3+ 536.*ws
      p10=  -6.5572*tpi2+2700.*ws-3*ypi
      pt0=   2.10  *tpi2- 783.*ws-3*tpi
      pls0=  0.42  *tpi3- 242.*ws
      pl210= 0.48  *tpi3+ 110.*ws
      pls20=-0.55  *tpi3-  44.*ws
      p01=  -8.1188*tpi2+2800.*ws-3*ypi
      pl201= 0.05  *tpi3+  63.*ws
      p00=  -9.12  *tpi2+5874.*ws+9*ypi
      pl200= 0.62  *tpi3- 363.*ws
c -----------
      if (lpot.ne.6) then
        pl2av=0.
c fix average of D waves
c       if (lpot.eq.5) pl2av=.5*(pl210+pl201)+pls20/3
        if (lpot.eq.5) pl2av=1.00*pl201+.00*(pl210+2*pls20/3)
c fix 1D2 partial wave
c       if (lpot.eq.5) pl2av=2*pl201/3
        p00=p00+2*(pl200-pl2av)
        pls0=pls0-2*(pl210-pl2av)-3*pls20
        p11=p11+2*(pl211-pl2av)+4*pls21/3
        pt1=pt1-5*pls21/12
        pls1=pls1-.5*pls21
      end if
c -----------
      vv(1)=.0625*(9*p11+3*p10+3*p01+p00)
      vv(2)=.0625*(3*p11-3*p10  +p01-p00)
      vv(3)=.0625*(3*p11  +p10-3*p01-p00)
      vv(4)=.0625*(  p11  -p10  -p01+p00)
      vv(5)=.25*(3*pt1+pt0)
      vv(6)=.25*(  pt1-pt0)
      vv(7)=.25*(3*pls1+pls0)
      vv(8)=.25*(  pls1-pls0)
      if (lpot.eq.4) go to 45
      if (lpot.eq.5) then
        vv(9)=pl2av
        go to 45
      end if
      vv(9)= .0625*(9*pl211+3*pl210+3*pl201+pl200)
      vv(10)=.0625*(3*pl211-3*pl210+  pl201-pl200)
      vv(11)=.0625*(3*pl211+  pl210-3*pl201-pl200)
      vv(12)=.0625*(  pl211-  pl210-  pl201+pl200)
      vv(13)=.25*(3*pls21+pls20)
      vv(14)=.25*(  pls21-pls20)
   45 vw(2)=ypi
      vw(3)=tpi
      go to 200
c ---------------------------
c argonne v18 and derivatives
c lpot = 7 -> v18-csbl
c        8 -> v18-csbs
c        9 -> v18
c       10 -> v8'
c       11 -> v6'
c       12 -> v4'
c       13 -> v2'
c       14 -> v1'
c       15 -> vx'
c       16 -> v9'(1D)
c       17 -> v9'(Da)
c       21 -> v18 v1.9
c       22 -> v8' v1.9
c       23 -> v18 v1.7
c       24 -> v8' v1.7
c ---------------------------
   50 hc=197.327053
      mpi0=134.9739
      mpic=139.5675
      mpi=(mpi0+2.*mpic)/3.
      mpi0=pimass1*mpi0
      mpic=pimass1*mpic
      mpis=pimass2*mpi
      mpiq=pimass3*mpi
      mu=mpi/hc
      mu0=mpi0/hc
      muc=mpic/hc
      mus=mpis/hc
      muq=mpiq/hc
      fsq=.075
      if (lpot.le.17) then
        cpi=2.1
        rws=wsrange*.5
        aws=.2*wsrange
      else if (lpot.eq.21 .or. lpot.eq.22) then
        cpi=1.9
        rws=wsrange*.525
        aws=.21*wsrange
      else if (lpot.eq.23 .or. lpot.eq.24) then
        cpi=1.7
        rws=wsrange*.55
        aws=.22*wsrange
      end if
      aiws=1./aws
      x=mu*rr
      xs=mus*rr
      xq=muq*rr
      x0=mu0*rr
      xc=muc*rr
      if (rr.le.small) then
        tpis=3*cpi**2*rr/mus**3
        tpiq=3*cpi**2*rr/muq**3
        ypi0=(mpi0/mpic)**2*(mpi0/3)*cpi*rr/mu0
        ypic=(mpic/3)*cpi*rr/muc
        tpi0=3*cpi*ypi0/mu0**2
        tpic=3*cpi*ypic/muc**2
      else
        rcut=1-exp(-cpi*rr*rr)
        tpis=rcut**2*(1+3/xs+3/xs**2)*exp(-xs)/xs
        tpiq=rcut**2*(1+3/xq+3/xq**2)*exp(-xq)/xq
        ypi0=(mpi0/mpic)**2*(mpi0/3)*exp(-x0)*rcut/x0
        ypic=(mpic/3)*exp(-xc)*rcut/xc
        tpi0=(1+(3+3/x0)/x0)*ypi0*rcut
        tpic=(1+(3+3/xc)/xc)*ypic*rcut
      end if
      pifac=pimass1**2
      ypi0=pifac*fsq*ypi0
      ypic=pifac*fsq*ypic
      tpi0=pifac*fsq*tpi0
      tpic=pifac*fsq*tpic
      tpi2=ftpec*(pimass2**3*tpis)**2
      tpi3=ftpec*(pimass3**3*tpiq)**2
      ews=exp((rr-rws)*aiws)
      ws=1/(1+ews)
      ews0=exp(-rws*aiws)
      ws0=1/(1+ews0)
      fws=aiws*ews0*ws0
      wsp=(1+fws*rr)*ws
      wsz=-aiws*ews*ws**2
      wspp=fws*ws+(1+fws*rr)*wsz
      wszz=-2*aiws*ews*ws*wsz-aiws**2*ews*ws**2
      wspdp=2*fws*wsz+(1+fws*rr)*wszz
      if (rr.le.small) then
        wspds=wspdp-2*ews0*aiws**2/(1+ews0)**2
      else
        wspds=wspdp+2*wspp/rr
      end if
      wsx=ws*x
      wsx2=wsx*x
      dypi00=(mpi0/mpic)**2*(mpi0/3)*cpi/mu0
      dypic0=(mpic/3)*cpi/muc
      ypi0p=ypi0-fsq*dypi00*ws*rr/ws0
      ypicp=ypic-fsq*dypic0*ws*rr/ws0
      ypi=(ypi0+2*ypic)/3
      tpi=(tpi0+2*tpic)/3
      ypibar=(ypi0-ypic)/3
      tpibar=(tpi0-tpic)/3
c final version 11/1/93
c nn potential added 11/15/93
c absolutely final version 3/29/94
c totally absolutely final version 6/28/94
      if (lpot.le.17) then
      p11pp=  -7.62701*tpi2+1815.4920*wsp+1847.8059*wsx2+ypi0p
      p11np=  -7.62701*tpi2+1813.5315*wsp+1847.8059*wsx2-ypi0p+2*ypicp
      p11nn=  -7.62701*tpi2+1811.5710*wsp+1847.8059*wsx2+ypi0p
      pt1pp=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0
      pt1np=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2-tpi0+2*tpic
      pt1nn=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0
      p01pp= -11.27028*tpi2+3346.6874*wsp-3*ypi0p
      p01np= -10.66788*tpi2+3126.5542*wsp-3*(-ypi0p+2*ypicp)
      p01nn= -11.27028*tpi2+3342.7664*wsp-3*ypi0p
c location for CD-CSB test version :
c keeps p11, pt1, pt1cs, p01, p01cd the same
c makes p11cd = p01cd/9, 
c       pt1cd = 5*p11cd/7 for tpi2 ; wsx adjusted to preserve phase shift
c while p11cs, p01cs are set with extra class iii-iv terms below
      if (lpot.eq.7 .or. lpot.eq.8) then
        p11pp=  -7.64930*tpi2+1821.6119*wsp+1847.8059*wsx2+ypi0p
        p11np=  -7.58243*tpi2+1797.3707*wsp+1847.8059*wsx2-ypi0p+2*ypicp
        p11nn=  -7.64930*tpi2+1821.6119*wsp+1847.8059*wsx2+ypi0p
        pt1pp=   1.06391*tpi2 -178.6199*wsx -811.2040*wsx2+tpi0
        pt1np=   1.11173*tpi2 -213.0449*wsx -811.2040*wsx2-tpi0+2*tpic
        pt1nn=   1.06391*tpi2 -178.6199*wsx -811.2040*wsx2+tpi0
        p01pp= -11.27028*tpi2+3344.7269*wsp-3*ypi0p
        p01np= -10.66788*tpi2+3126.5542*wsp-3*(-ypi0p+2*ypicp)
        p01nn= -11.27028*tpi2+3344.7269*wsp-3*ypi0p
      end if
c CD test complete
      pls1=    -.62697*tpi3 -570.5571*wsp +819.1222*wsx2
      pl211=    .06709*tpi3 +342.0669*wsp -615.2339*wsx2
      pls21=    .74129*tpi3   +9.3418*wsp -376.4384*wsx2
      p10=    -8.62770*tpi2+2605.2682*wsp +441.9733*wsx2-ypi0p-2*ypicp
      pt0=    1.485601*tpi2-1126.8359*wsx +370.1324*wsx2-tpi0-2*tpic
      pls0=     .10180*tpi3  +86.0658*wsp -356.5175*wsx2
      pl210=   -.13201*tpi3 +253.4350*wsp   -1.0076*wsx2
      pls20=    .07357*tpi3 -217.5791*wsp  +18.3935*wsx2
      pl201=    .12472*tpi3  +16.7780*wsp
      p00=    -2.09971*tpi2+1204.4301*wsp-3*(-ypi0p-2*ypicp)
      pl200=   -.31452*tpi3 +217.4559*wsp
c location for v1.9
      else if (lpot.eq.21 .or. lpot.eq.22) then
        p11pp=  -8.25333*tpi2+1629.4894*wsp+1007.3079*wsx2+ypi0p
        p11np=  -8.25333*tpi2+1627.8623*wsp+1007.3079*wsx2-ypi0p+2*ypicp
        p11nn=  -8.25333*tpi2+1626.2352*wsp+1007.3079*wsx2+ypi0p
        pt1pp=   1.22738*tpi2 -331.5020*wsx -415.4240*wsx2+tpi0
        pt1np=   1.22738*tpi2 -331.5020*wsx -415.4240*wsx2-tpi0+2*tpic
        pt1nn=   1.22738*tpi2 -331.5020*wsx -415.4240*wsx2+tpi0
        pls1=   -1.24596*tpi2 -438.1866*wsp +881.8829*wsx2
        pl211=   0.17268*tpi2 +210.3707*wsp -418.9703*wsx2
        pls21=   0.68968*tpi2  -44.1763*wsp -154.7568*wsx2
        p10=   -10.62968*tpi2+2297.1952*wsp +503.6560*wsx2-ypi0p-2*ypicp
        pt0=     1.43163*tpi2 -932.7628*wsx +415.7518*wsx2-tpi0-2*tpic
        pls0=    0.31692*tpi2   +5.9540*wsp -261.4438*wsx2
        pl210=   0.20369*tpi2 +164.8268*wsp -133.2324*wsx2
        pls20=  -0.08370*tpi2 -162.9074*wsp  +92.1321*wsx2
        p01pp= -11.24918*tpi2+2446.4156*wsp -278.5780*wsx2-3*ypi0p
        p01np= -10.57598*tpi2+2273.0877*wsp -289.7548*wsx2
     &                                          -3*(-ypi0p+2*ypicp)
        p01nn= -11.24918*tpi2+2443.1614*wsp -278.5780*wsx2-3*ypi0p
        pl201=   0.12622*tpi2  +23.4445*wsp  -13.5987*wsx2
        p00=    -2.14060*tpi2+1000.2218*wsp -167.8362*wsx2
     &                                          -3*(-ypi0p-2*ypicp)
        pl200=  -0.30287*tpi2 +178.4343*wsp  -40.3099*wsx2
c v1.9 complete
c location for v1.7
      else if (lpot.eq.23 .or. lpot.eq.24) then
        p11pp=  -9.10461*tpi2+1383.9447*wsp +635.3480*wsx2+ypi0p
        p11np=  -9.10461*tpi2+1382.5743*wsp +635.3480*wsx2-ypi0p+2*ypicp
        p11nn=  -9.10461*tpi2+1381.2039*wsp +635.3480*wsx2+ypi0p
        pt1pp=   1.41302*tpi2 -353.5571*wsx -220.1474*wsx2+tpi0
        pt1np=   1.41302*tpi2 -353.5571*wsx -220.1474*wsx2-tpi0+2*tpic
        pt1nn=   1.41302*tpi2 -353.5571*wsx -220.1474*wsx2+tpi0
        pls1=   -2.32067*tpi2 -294.5597*wsp +905.8846*wsx2
        pl211=   0.31011*tpi2 +124.9726*wsp -295.2155*wsx2
        pls21=   0.56878*tpi2  -53.4883*wsp  -35.2622*wsx2
        p10=   -12.16378*tpi2+1841.4730*wsp +490.2208*wsx2-ypi0p-2*ypicp
        pt0=     1.38999*tpi2 -748.1056*wsx +373.6946*wsx2-tpi0-2*tpic
        pls0=    0.38755*tpi2  -17.6948*wsp -168.8546*wsx2
        pl210=   0.61367*tpi2 +108.7064*wsp -206.9220*wsx2
        pls20=  -0.35600*tpi2 -118.0195*wsp +148.1324*wsx2
        p01pp= -11.28159*tpi2+1735.7847*wsp -290.3768*wsx2-3*ypi0p
        p01np= -11.62496*tpi2+1724.8767*wsp -116.8103*wsx2
     &                                          -3*(-ypi0p+2*ypicp)
        p01nn= -11.28159*tpi2+1733.0438*wsp -290.3768*wsx2-3*ypi0p
        pl201=   0.12387*tpi2  +24.9311*wsp  -14.7805*wsx2
        p00=    -1.94411*tpi2 +808.3897*wsp -277.7773*wsx2
     &                                          -3*(-ypi0p-2*ypicp)
        pl200=  -0.31130*tpi2 +139.5016*wsp  -40.4310*wsx2
c v1.7 complete
      end if
      p11=(p11pp+p11nn+p11np)/3
      p11cd=(.5*(p11pp+p11nn)-p11np)/6
      p11cs=(p11pp-p11nn)/4
      pt1=(pt1pp+pt1nn+pt1np)/3
      pt1cd=(.5*(pt1pp+pt1nn)-pt1np)/6
      pt1cs=(pt1pp-pt1nn)/4
      p01=(p01pp+p01nn+p01np)/3
      p01cd=(.5*(p01pp+p01nn)-p01np)/6
      p01cs=(p01pp-p01nn)/4
c -----------
      if (lpot.ge.10 .and. lpot.ne.21 .and. lpot.ne.23) then
        pl2av=0.
c fix v9' L**2 term to 1D2
        if (lpot.eq.16) pl2av=1.00*pl201+0.00*(pl210+2*pls20/3)
c fix v9' L**2 term to average of 1D2 and 3DJ
        if (lpot.eq.17) pl2av=0.50*pl201+0.50*(pl210+2*pls20/3)
        p00=p00+2*(pl200-pl2av)
        pls0=pls0-2*(pl210-pl2av)-3*pls20
        p11=p11+2*(pl211-pl2av)+4*pls21/3
        pt1=pt1-5*pls21/12
        pls1=pls1-.5*pls21
c fix deuteron in v6' case
        if (lpot.ge.11 .and. lpot.le.15) p10=p10-.3*pls0
c fix deuteron in v4' case
        if (lpot.ge.12 .and. lpot.le.15) p10=p10+.8735*pt0
c project only vc and vt in v2' case
        if (lpot.eq.13) then
          vv(1)=.25*(3*p01+p10)
          vv(2)=.25*(  p01-p10)
          go to 200
c average 1s & 3s in v1' case
        else if (lpot.eq.14) then
          vv(1)=.5*(p01+p10)
          go to 200
c combination for vx' case
        else if (lpot.eq.15) then
          vv(1)=.0625*(9*p11+3*p10+3*p01+p00)
          vv(2)=.0125*(9*p11-5*p10-5*p01+p00)
          vv(3)=vv(2)
          vv(4)=vv(2)
          go to 200
        end if
      end if
c -----------
      vv(1)=.0625*(9*p11+3*p10+3*p01+p00)
      vv(2)=.0625*(3*p11-3*p10  +p01-p00)
      vv(3)=.0625*(3*p11  +p10-3*p01-p00)
      vv(4)=.0625*(  p11  -p10  -p01+p00)
      if (lpot.eq.12) go to 200
      vv(5)=.25*(3*pt1+pt0)
      vv(6)=.25*(  pt1-pt0)
      if (lpot.eq.11) go to 200
      vv(7)=.25*(3*pls1+pls0)
      vv(8)=.25*(  pls1-pls0)
      if (lpot.eq.10 .or. lpot.eq.22 .or. lpot.eq.24) go to 55
      if (lpot.eq.16 .or. lpot.eq.17) then
        vv(9)=pl2av
        go to 55
      end if
      vv(9)= .0625*(9*pl211+3*pl210+3*pl201+pl200)
      vv(10)=.0625*(3*pl211-3*pl210+  pl201-pl200)
      vv(11)=.0625*(3*pl211+  pl210-3*pl201-pl200)
      vv(12)=.0625*(  pl211-  pl210-  pl201+pl200)
      vv(13)=.25*(3*pls21+pls20)
      vv(14)=.25*(  pls21-pls20)
   54 vv(15)=.25*(3*p11cd+p01cd)
      vv(16)=.25*(  p11cd-p01cd)
      vv(17)=pt1cd
      vv(19)=p01cs
   55 vw(2)=ypi
      vw(3)=tpi
c delta function added to ypi!
c     xlm=900./mpi
c     xlm2=xlm*xlm
c     rmu=mu*rr
c     rlm=xlm*rmu
c     ermu=exp(-rmu)/rmu
c     erlm=exp(-rlm)/rmu
c     vw(2)=fsq*(mpi/3)
c    &     *(ermu-xlm2*erlm-.5*xlm*(xlm2-1.)*(1.-2./rlm)*erlm*rmu)
c     vw(3)=fsq*(mpi/3)
c    &     *((1.+3./rmu+3./rmu**2)*ermu-xlm2*(1.+3./rlm+3./rlm**2)*erlm
c    &       -.5*xlm*(xlm2-1.)*(1.+1./rlm)*erlm*rmu)
c     b=4.27
c     br=b*rr
c     fdelta=b**3*(1+br+br**2/3)*exp(-br)/16
c     vw(2)=ypi-fsq*(mpi/3)*fdelta/mu**3
c delta function added to ypi!
      vw(4)=ypibar
      vw(5)=tpibar
c extra class iii-iv csb
      if (lpot.eq.7 .or. lpot.eq.8) then
c -----------------------
c following is for case 1
        if (lpot.eq.7) then
          vcsb=1.50875
          krho=6.1
          komg=0.14
c following is for case 2
        else if (lpot.eq.8) then
          vcsb=1.11160
          krho=3.7
          komg=-.12
        end if
c -----------------------
        vv(19)=vcsb*wsp
        if (rr.ge.vsmall) then
        deltam=(mn-mp)/(mn+mp)
        crho=2.44
        grho=0.55
        fmsqi=vcsb*(hc/(mn+mp))**2
        vw(6)=fmsqi*(2+komg+krho)*wspds
        vv(20)=2*fmsqi*(1+komg)*(1+krho)*wspds/3
        vv(21)=-fmsqi*(1+komg)*(1+krho)*(wspds-3*wspp/rr)/3
        vv(22)=4*fmsqi*(2+komg+krho)*wspp/rr
        vw(7)=4*fmsqi*(krho-komg)*wspp/rr
        mrho=770.
        murho=mrho/hc
        xrho=murho*rr
        rcut=1-exp(-cpi*rr*rr)
        zpic=mpic*((1+1/xc)/xc)*exp(-xc)*rcut**1.5/xc
        zcut=1-exp(-crho*rr*rr)
        zrho=mrho*((1+1/xrho)/xrho)*exp(-xrho)*zcut**1.5/xrho
        vw(8)=2*fsq*deltam*zpic
        vw(9)=2*grho*deltam*(mrho/(mn+mp))**2*zrho
        vw(10)=2*(1+krho)**2*vw(9)
        end if
      end if
c williams thomas miller version
c     lam=1300./hc
c     rmu=mu*rr
c     rlam=lam*rr
c     ermu=exp(-rmu)
c     erlam=exp(-rlam)
c     z1p=-mu*(1+1/rmu)*ermu/rmu
c    &    +lam*(1+1/rlam)*erlam/rmu
c    &    +(lam**2-mu**2)*erlam/(2*mu)
c     vw(11)=2*fsq*deltam*(-hc*z1p/rmu)
c     lamp=1400./hc
c     rmurho=murho*rr
c     rlamp=lamp*rr
c     ermurho=exp(-rmurho)
c     erlamp=exp(-rlamp)
c     z1pp=-murho*(1+1/rmurho)*ermurho/rmurho
c    &    +lamp*(1+1/rlamp)*erlamp/rmurho
c    &    +(lamp**2-murho**2)*erlamp/(2*murho)
c     vw(12)=grho*(mrho/(mn+mp))**2*deltam*(-hc*z1pp/rmurho)
c     vw(13)=2*(1+krho)**2*vw(7)
c     gomg=8.1
c     momg=783.
c     muomg=momg/hc
c     lamw=lamp*momg/mrho
c     rmuomg=muomg*rr
c     rlamw=lamw*rr
c     ermuomg=exp(-rmuomg)
c     erlamw=exp(-rlamw)
c     z1pw=-muomg*(1+1/rmuomg)*ermuomg/rmuomg
c    &    +lamw*(1+1/rlamw)*erlamw/rmuomg
c    &    +(lamw**2-muomg**2)*erlamw/(2*muomg)
c     fpgw=krho*sqrt(grho/gomg)*gomg
c     hwp=-3400.
c     vw(14)=fpgw*(2/(mn+mp)**2)*(hwp/(momg**2-mrho**2))
c    &     *(-hc*z1pp*mrho**2/rmurho+hc*z1pw*momg**2/rmuomg)
      go to 200
c ----------------------------------------
c argonne v18p  (p**2 terms)
c argonne v18pq (p**2 terms & l**2 tensor)
c lpot = 26 -> v18p
c        27 -> v18pq
c ----------------------------------------
   60 hc=197.327053
      mpi0=134.9739
      mpic=139.5675
      mpi=(mpi0+2*mpic)/3
      mu0=mpi0/hc
      muc=mpic/hc
      mu=mpi/hc
      fsq=.075
      cpi=2.1
      rws=.5
      aiws=5.
      x=mu*rr
      x0=mu0*rr
      xc=muc*rr
      if (rr.le.small) then
        tpi=3*cpi**2*rr/mu**3
        ypi0=(mpi0/mpic)**2*(mpi0/3)*cpi*rr/mu0
        ypic=(mpic/3)*cpi*rr/muc
        tpi0=3*cpi*ypi0/mu0**2
        tpic=3*cpi*ypic/muc**2
      else
        expcut=exp(-cpi*rr*rr)
        rcut=1-expcut
        rcutp=2*cpi*rr*expcut
        rcutdp=2*cpi*expcut*(1-2*cpi*rr**2)
        ypi0=(mpi0/mpic)**2*(mpi0/3)*exp(-x0)*rcut/x0
        ypic=(mpic/3)*exp(-xc)*rcut/xc
        tpi0=(1+(3+3/x0)/x0)*ypi0*rcut
        tpic=(1+(3+3/xc)/xc)*ypic*rcut
        y=exp(-x)/x
        yp=-mu*(1+1/x)*y
        ydp=mu**2*(1+2*(1+1/x)/x)*y
        t=1+3*(1+1/x)/x
        tp=-mu*3*(1+2/x)/x**2
        tdp=mu**2*6*(1+3/x)/x**3
        tpi=t*y*rcut**2
        tpip=(tp*y+t*yp)*rcut**2+2*t*y*rcut*rcutp
        tpidp=(tdp*y+2*tp*yp+t*ydp)*rcut**2+4*(tp*y+t*yp)*rcut*rcutp
     &       +2*t*y*(rcut*rcutdp+rcutp**2)
      end if
      ypi0=fsq*ypi0
      ypic=fsq*ypic
      tpi0=fsq*tpi0
      tpic=fsq*tpic
      tpi2=tpi*tpi
      tpi2p=2*tpi*tpip
      tpi2dp=2*(tpi*tpidp+tpip**2)
      expws=exp((rr-rws)*aiws)
      ws=1./(1.+expws)
      wsp=-aiws*expws*ws**2
      wsdp=aiws*wsp+2*wsp**2/ws
      ws0=1./(1+exp(-rws*aiws))
      wsp0=-aiws*exp(-rws*aiws)*ws0**2
      wsm=ws*(1.-wsp0*x/(mu*ws0))
      wsmp=wsp-wsp0*(x*wsp+mu*ws)/(mu*ws0)
      wsmdp=wsdp-wsp0*(x*wsdp+2*mu*wsp)/(mu*ws0)
      wsx=ws*x
      wsx2=wsx*x
      wsx2p=x**2*wsp+2*mu*x*ws
      wsx2dp=x**2*wsdp+4*mu*x*wsp+2*mu**2*ws
      dypi00=(mpi0/mpic)**2*(mpi0/3)*cpi/mu0
      dypic0=(mpic/3)*cpi/muc
      ypi0p=ypi0-fsq*dypi00*ws*rr/ws0
      ypicp=ypic-fsq*dypic0*ws*rr/ws0
      ypi=(ypi0+2*ypic)/3
      tpi=(tpi0+2*tpic)/3
      ypibar=(ypi0-ypic)/3
      tpibar=(tpi0-tpic)/3
c ----------------------------------------
c first version 2002.07.18
c second version 2003.01.25 (deuteron bug fixed)
c plain p**2 version 2003.02.13
c inconsistent pl2p0x and pl2dp0x corrected 2003.08.04
c ----------------------------------------
      if (lpot.eq.26) then
      p11pp= -11.059567*tpi2+2828.6889*wsm+3178.7405*wsx2+ypi0p
      p11np= -11.059567*tpi2+2826.8983*wsm+3178.7405*wsx2-ypi0p+2*ypicp
      p11nn= -11.059567*tpi2+2825.1077*wsm+3178.7405*wsx2+ypi0p
      pt1pp=   1.514600*tpi2 -376.1341*wsx -955.2246*wsx2+tpi0
      pt1np=   1.514600*tpi2 -376.1341*wsx -955.2246*wsx2-tpi0+2*tpic
      pt1nn=   1.514600*tpi2 -376.1341*wsx -955.2246*wsx2+tpi0
      pls1=   -0.732326*tpi2 -576.5707*wsm +925.8808*wsx2
      pl211=  -0.420808*tpi2 +209.2711*wsm  -51.6679*wsx2
      pls21=   1.023793*tpi2  -18.8040*wsm -544.7574*wsx2
      p10=    -8.061451*tpi2+1938.3975*wsm+1895.2537*wsx2-ypi0p-2*ypicp
      pt0=     1.072145*tpi2-1220.5266*wsx +872.6085*wsx2-tpi0-2*tpic
      pls0=    0.328041*tpi2  -74.7348*wsm -271.1134*wsx2
      pl210=   0.226986*tpi2  -19.8530*wsm +104.6043*wsx2
      pls20=  -0.091272*tpi2  -67.3591*wsm  -81.6822*wsx2
c ----------------------------------------
      else if (lpot.eq.27) then
      p11pp=  -9.882847*tpi2+2589.7742*wsm+2952.3910*wsx2+ypi0p
      p11np=  -9.882847*tpi2+2587.9836*wsm+2952.3910*wsx2-ypi0p+2*ypicp
      p11nn=  -9.882847*tpi2+2586.1930*wsm+2952.3910*wsx2+ypi0p
      pt1pp=   1.420069*tpi2 -453.5357*wsx -837.3820*wsx2+tpi0
      pt1np=   1.420069*tpi2 -453.5357*wsx -837.3820*wsx2-tpi0+2*tpic
      pt1nn=   1.420069*tpi2 -453.5357*wsx -837.3820*wsx2+tpi0
      pls1=   -1.749197*tpi2 -493.8470*wsm+1533.0637*wsx2
      pl211=  -0.008159*tpi2 +132.7694*wsm -169.8510*wsx2
      pls21=   0.135181*tpi2  -17.7975*wsm  -46.2542*wsx2
      p10=    -8.351808*tpi2+2325.5929*wsm +957.8091*wsx2-ypi0p-2*ypicp
      pt0=     1.327862*tpi2-1170.8528*wsx +580.5596*wsx2-tpi0-2*tpic
      pls0=    0.060223*tpi2  +58.3208*wsm -126.0235*wsx2
      pl210=  -0.023577*tpi2   +1.8164*wsm +127.2921*wsx2
      pls20=   0.000589*tpi2  -25.1123*wsm   -4.6897*wsx2
      end if
c ----------------------------------------
      p01pp= -10.518030*tpi2+2836.0715*wsm +651.1945*wsx2-3*ypi0p
      p01np= -10.812190*tpi2+2816.4190*wsm+1002.5300*wsx2
     &                                               -3*(-ypi0p+2*ypicp)
      p01nn= -10.518030*tpi2+2832.4903*wsm +651.1945*wsx2-3*ypi0p
      pl201=   0.134747*tpi2   -9.4691*wsm
      p00=    -4.739629*tpi2+1121.2225*wsm+2764.3395*wsx2
     &                                               -3*(-ypi0p-2*ypicp)
      pl200=  -0.227084*tpi2 +166.5629*wsm
c ----------------------------------------
      p11=(p11pp+p11nn+p11np)/3
      p11cd=(.5*(p11pp+p11nn)-p11np)/6
      p11cs=(p11pp-p11nn)/4
      pt1=(pt1pp+pt1nn+pt1np)/3
      pt1cd=(.5*(pt1pp+pt1nn)-pt1np)/6
      pt1cs=(pt1pp-pt1nn)/4
      p01=(p01pp+p01nn+p01np)/3
      p01cd=(.5*(p01pp+p01nn)-p01np)/6
      p01cs=(p01pp-p01nn)/4
      vv(1)=.0625*(9*p11+3*p10+3*p01+p00)
      vv(2)=.0625*(3*p11-3*p10  +p01-p00)
      vv(3)=.0625*(3*p11  +p10-3*p01-p00)
      vv(4)=.0625*(  p11  -p10  -p01+p00)
      vv(5)=.25*(3*pt1+pt0)
      vv(6)=.25*(  pt1-pt0)
      vv(7)=.25*(3*pls1+pls0)
      vv(8)=.25*(  pls1-pls0)
      vv(9)= .0625*(9*pl211+3*pl210+3*pl201+pl200)
      vv(10)=.0625*(3*pl211-3*pl210+  pl201-pl200)
      vv(11)=.0625*(3*pl211+  pl210-3*pl201-pl200)
      vv(12)=.0625*(  pl211-  pl210-  pl201+pl200)
      vv(13)=.25*(3*pls21+pls20)
      vv(14)=.25*(  pls21-pls20)
      vv(15)=.25*(3*p11cd+p01cd)
      vv(16)=.25*(  p11cd-p01cd)
      vv(17)=pt1cd
      vv(19)=p01cs
      vw(2)=ypi
      vw(3)=tpi
      vw(4)=ypibar
      vw(5)=tpibar
c -----------------------
c vpsq & derivative terms
c -----------------------
      if (lpot.eq.26) then
      pl2p11=   -0.420808*tpi2p  +209.2711*wsmp   -51.6679*wsx2p
      pl2p10=    0.226986*tpi2p   -19.8530*wsmp  +104.6043*wsx2p
      pl2dp11=  -0.420808*tpi2dp +209.2711*wsmdp  -51.6679*wsx2dp
      pl2dp10=   0.226986*tpi2dp  -19.8530*wsmdp +104.6043*wsx2dp
c -----------------------
      else if (lpot.eq.27) then
      pl2p11=   -0.008159*tpi2p  +132.7694*wsmp  -169.8510*wsx2p
      pl2p10=   -0.023577*tpi2p    +1.8164*wsmp  +127.2921*wsx2p
      pl2dp11=  -0.008159*tpi2dp +132.7694*wsmdp -169.8510*wsx2dp
      pl2dp10=  -0.023577*tpi2dp   +1.8164*wsmdp +127.2921*wsx2dp
      end if
c -----------------------
      pl2p01=    0.134747*tpi2p    -9.4691*wsmp
      pl2p00=   -0.227084*tpi2p  +166.5629*wsmp
      pl2dp01=   0.134747*tpi2dp   -9.4691*wsmdp
      pl2dp00=  -0.227084*tpi2dp +166.5629*wsmdp
c ----------------------------------------
      psq11=.5*rr**2*pl211
      psq10=.5*rr**2*pl210
      psq01=.5*rr**2*pl201
      psq00=.5*rr**2*pl200
      psqp11=rr*pl211+.5*rr**2*pl2p11
      psqp10=rr*pl210+.5*rr**2*pl2p10
      psqp01=rr*pl201+.5*rr**2*pl2p01
      psqp00=rr*pl200+.5*rr**2*pl2p00
      psqdp11=pl211+2.*rr*pl2p11+.5*rr**2*pl2dp11
      psqdp10=pl210+2.*rr*pl2p10+.5*rr**2*pl2dp10
      psqdp01=pl201+2.*rr*pl2p01+.5*rr**2*pl2dp01
      psqdp00=pl200+2.*rr*pl2p00+.5*rr**2*pl2dp00
      vp(1)=.0625*(9*psq11+3*psq10+3*psq01+psq00)
      vp(2)=.0625*(3*psq11-3*psq10  +psq01-psq00)
      vp(3)=.0625*(3*psq11  +psq10-3*psq01-psq00)
      vp(4)=.0625*(  psq11  -psq10  -psq01+psq00)
      vp(5)=.0625*(9*psqp11+3*psqp10+3*psqp01+psqp00)
      vp(6)=.0625*(3*psqp11-3*psqp10+  psqp01-psqp00)
      vp(7)=.0625*(3*psqp11+  psqp10-3*psqp01-psqp00)
      vp(8)=.0625*(  psqp11-  psqp10-  psqp01+psqp00)
      vp(9)= .0625*(9*psqdp11+3*psqdp10+3*psqdp01+psqdp00)
      vp(10)=.0625*(3*psqdp11-3*psqdp10+  psqdp01-psqdp00)
      vp(11)=.0625*(3*psqdp11+  psqdp10-3*psqdp01-psqdp00)
      vp(12)=.0625*(  psqdp11-  psqdp10-  psqdp01+psqdp00)
      go to 200
c ----------------------
c super-soft core(c) v14
c lpot = 28 -> v14
c        29 -> mod v8'
c ----------------------
   70 if (rr.le.vsmall) rr=vsmall
      x=.7*rr
      rr4=rr**4
      rc4=1-exp(-rr4)
      rc6=1-exp(-rr**6)
      hr=10.463
      p11=144.83*exp(-rr4/.88787**2)
     &   +(-241.34*yc(3.3788*x)+(hr/3)*yc(x))*rc4
      p10=215.32*exp(-rr4/.85807**2)
     &   +(-883.6*yc(3.5042*x)-hr*yc(x))*rc4
      p01=375.*exp(-rr4/.47552**2)
     &   +(-1001.6*yc(3.6071*x)-hr*yc(x))*rc4
      p00=75.653*exp(-rr4/3.0000**2)
     &   +(-286.26*yc(2.0254*x)+3*hr*yc(x))*rc4
      pt1=36.*exp(-rr4/1.0805**2)
     &   +(-110.*yt(3.9529*x)+(hr/3)*yt(x))*rc6
      pt0=-58.951*exp(-rr4/1.3171**2)
     &   +(395.18*yt(4.3098*x)-hr*yt(x))*rc6
      pls1=(520.*yls(5.661*x)-54.85*yls(4.0141*x))*rc6
      pls0=(-40.466*yls(5.768*x)-40.408*yls(4.0676*x))*rc6
      pl211=(6.65*yl2(1.965*x)-.959*yl2(x))*rc6
      pl210=(17.626*yl2(2.6463*x)-.35261*yl2(x))*rc6
      pl201=(14.*yl2(2.5*x)-.35*yl2(x))*rc6
      pl200=(15.633*yl2(2.01*x)+.72581*yl2(x))*rc6
      pq0=-3.9904*yl2(2.4583*x)*rc6
c fix v8'
      if (lpot.eq.29) then
         p00=p00+2*pl200
         pls0=pls0-2*pl210-10*pq0
         p11=p11+2*pl211
c ------------------------
c option for modified v8'
c ------------------------
         p11=p11-111*yc(3.3788*x)*rc4
c ------------------------
      end if
      vv(1)=.0625*(9*p11+3*p10+3*p01+p00)
      vv(2)=.0625*(3*p11-3*p10+  p01-p00)
      vv(3)=.0625*(3*p11+  p10-3*p01-p00)
      vv(4)=.0625*(  p11-  p10-  p01+p00)
      vv(5)=.25*(3*pt1+pt0)
      vv(6)=.25*(  pt1-pt0)
      vv(7)=.25*(3*pls1+pls0)+.75*pq0
      vv(8)=.25*(  pls1-pls0)-.75*pq0
      if (lpot.eq.29) go to 200
      vv(9)= .0625*(9*pl211+3*pl210+3*pl201+pl200)-.75*pq0
      vv(10)=.0625*(3*pl211-3*pl210+  pl201-pl200)+.75*pq0
      vv(11)=.0625*(3*pl211+  pl210-3*pl201-pl200)-.25*pq0
      vv(12)=.0625*(  pl211-  pl210-  pl201+pl200)+.25*pq0
      vv(13)=1.5*pq0
      vv(14)=-1.5*pq0
      vw(2)=(hr/3)*yc(x)
      vw(3)=(hr/3)*yt(x)
      rr=rrsave
      go to 200
c -----
c paris
c -----
   80 p01=0
      p11=0
      p00=0
      p10=0
      pt1=0
      pt0=0
      pls1=0
      pls0=0
      pl201=0
      pl211=0
      pl200=0
      pl210=0
      pls21=0
      pls20=0
      do 105 i=1,12
        if (rr.le.small) then
          rm1=rr*pm1(i)
          rm0=rr*pm0(i)
          vc1=-1+.5*rm1
          vc0=-1+.5*rm0
          vt1=.125*rm1
          vt0=.125*rm0
          vls1=1/3.-.125*rm1
          vls0=1/3.-.125*rm0
          vso21=-1/15.+rm1/48
          vso20=-1/15.+rm0/48
        else
          vc1=pc(pm1(i)*rr)
          vc0=pc(pm0(i)*rr)
          vt1=pt(pm1(i)*rr)*vc1
          vt0=pt(pm0(i)*rr)*vc0
          vls1=pls(pm1(i)*rr)*vc1
          vls0=pls(pm0(i)*rr)*vc0
          vso21=pso2(pm1(i)*rr)*vc1
          vso20=pso2(pm0(i)*rr)*vc0
        end if
        p01=p01+pgva01(i)*vc1
        p11=p11+pgva11(i)*vc1
        p00=p00+pgva00(i)*vc0
        p10=p10+pgva10(i)*vc0
        pt1=pt1+pgvt1(i)*vt1
        pt0=pt0+pgvt0(i)*vt0
        pls1=pls1+pgvls1(i)*vls1
        pls0=pls0+pgvls0(i)*vls0
        pl201=pl201+pgvb01(i)*vc1
        pl211=pl211+pgvb11(i)*vc1
        pl200=pl200+pgvb00(i)*vc0
        pl210=pl210+pgvb10(i)*vc0
        pls21=pls21+pgvso21(i)*vso21
        pls20=pls20+pgvso20(i)*vso20
  105 continue
      vv(1)=.0625*(9*p11+3*p10+3*p01+p00)
      vv(2)=.0625*(3*p11-3*p10+  p01-p00)
      vv(3)=.0625*(3*p11+  p10-3*p01-p00)
      vv(4)=.0625*(  p11-  p10-  p01+p00)
      vv(5)=.25*(3*pt1+pt0)
      vv(6)=.25*(  pt1-pt0)
      vv(7)=.25*(3*pls1+pls0)
      vv(8)=.25*(  pls1-pls0)
      vv(9)= .0625*(9*pl211+3*pl210+3*pl201+pl200)
      vv(10)=.0625*(3*pl211-3*pl210+  pl201-pl200)
      vv(11)=.0625*(3*pl211+  pl210-3*pl201-pl200)
      vv(12)=.0625*(  pl211-  pl210-  pl201+pl200)
      vv(13)=.25*(3*pls21+pls20)
      vv(14)=.25*(  pls21-pls20)
c --------------------
  200 if (lr.ge.1) then
        do 300 l=1,18
          vv(l)=vv(l)*rr**lr
  300   continue
        vw(1)=vw(1)*rr**lr
        vw(2)=vw(2)*rr**lr
        vw(3)=vw(3)*rr**lr
        vw(4)=vw(4)*rr**lr
        vw(5)=vw(5)*rr**lr
      end if
      return
      end
c *id* empot ***********************************************************
c subroutine for electromagnetic potential
c arguments for empot
c lemp: 0=Coulomb w/ff only
c       1=full electromagnetic potential
c xmn:    parameter for nucleon mass variation
c ems:  parameter for electromagnetic alpha variation
c rr:   position r in fm
c vem:  returned potential in MeV (14 component array)
c ----------------------------------------------------------------------
c order of operators in vem(l)
c l:    1=C1    (pp)          2=DF    (pp)          3=C2      (pp)
c       4=VP    (pp)                                5=C1      (np)
c       6=s1.s2 (pp)          7=s1.s2 (nn)          8=s1.s2   (np)
c       9=S12   (pp)         10=S12   (nn)         11=S12     (np)
c      12=L.S   (pp)         13=L.S   (nn)         14=L.S     (np)
c C1 = one-photon-exchange Coulomb with form factor
c C2 = two-photon-exchange Coulomb
c DF = Darwin-Foldy
c VP = vacuum polarization (short-range approximation)
c all other terms from magnetic moment (MM) interactions
c ----------------------------------------------------------------------
      subroutine empot(lemp,xmn,ems,rr,vem)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      dimension vem(14)
      real*8 kr,me,mp,mn,mr,mun,mup
      small=1e-4
      vsmall=1e-10
      alpha=ems/137.03599
      hc=197.327053
      b=4.27
      rrsave=rr
      if (rr.le.small) rr=small
      br=b*rr
      fcoul=1-(1+11*br/16+3*br**2/16+br**3/48)*exp(-br)
      vem(1)=alpha*hc*fcoul/rr
      if (lemp.ge.1) then
        pi=acos(-1.)
        me=0.510999
        mp=938.27231*xmn
        mn=939.56563*xmn
        mr=mp*mn/(mp+mn)
        mup=2.79285
        mun=-1.91304
        gamma=0.577216
        beta=.0189
        kr=me*rr/hc
        ftensor=1-(1+br+br**2/2+br**3/6+br**4/24+br**5/144)*exp(-br)
        fspinor=1-(1+br+br**2/2+7*br**3/48+br**4/48)*exp(-br)
        fdelta=b**3*(1+br+br**2/3)*exp(-br)/16
        fnp=b**2*(15*br+15*br**2+6*br**3+br**4)*exp(-br)/384
        fivp=-gamma-5./6.+abs(log(kr))+6*pi*kr/8
        vem(2)=-alpha*hc**3*fdelta/(4*mp**2)
        vem(3)=-vem(1)**2/mp
        vem(4)=2*alpha*vem(1)*fivp/(3*pi)
        vem(5)=alpha*hc*beta*fnp/rr
        vem(6)=-alpha*hc**3*mup**2*fdelta/(6*mp**2)
        vem(7)=-alpha*hc**3*mun**2*fdelta/(6*mn**2)
        vem(8)=-alpha*hc**3*mup*mun*fdelta/(6*mn*mp)
        vem(9)=-alpha*hc**3*mup**2*ftensor/(4*mp**2*rr**3)
        vem(10)=-alpha*hc**3*mun**2*ftensor/(4*mn**2*rr**3)
        vem(11)=-alpha*hc**3*mup*mun*ftensor/(4*mp*mn*rr**3)
        vem(12)=-alpha*hc**3*(4*mup-1)*fspinor/(2*mp**2*rr**3)
        vem(13)=0
        vem(14)=-alpha*hc**3*mun*fspinor/(2*mn*mr*rr**3)
      else
        do 10 i=2,14
          vem(i)=0
   10   continue
      end if
      rr=rrsave
      return
      end
