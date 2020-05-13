      parameter (nr =2000) ! nr must be even)
      parameter (nz =1001)
      parameter (nx =201)
      parameter (ng =51)
      parameter (ne =15)
      parameter (nlm=18)
      parameter (ncf=2)
      parameter (hr =0.01d0)
      parameter (xx0=1.25d0)

      parameter (xf=200.d0)
      parameter (pi =dacos(-1.d0 ))
      parameter (pi2=pi*pi)
      parameter (pi3=pi*pi2)

      parameter (hbarc=197.32697d0)
      parameter (ampi0=134.9766d0/hbarc)
      parameter (ampic=139.5702d0/hbarc)
      parameter (amp  =( ampi0+2.d0*ampic )/3.d0)
      parameter (amd  =293.1d0/hbarc)

      parameter (ap02=ampi0**2)
      parameter (apc2=ampic**2)
      parameter (amp2=amp**2)
      parameter (amd2=amd**2)

      parameter (gan0=  1.29d0)
      parameter (gan2=gan0**2)
      parameter (gan4=gan2**2)
      parameter (fpp0=184.80d0/hbarc)

      parameter (ct2 =( gan0/fpp0 )**2)
      parameter (ct4 =ct2**2/( 24.d0*pi2 ))
      parameter (ct48=ct2**2/(  8.d0*pi3))

      parameter (dm2=- 5.d0+ 4.d0/gan0**2+1.d0/gan0**4)
      parameter (dk2=-23.d0+10.d0/gan0**2+1.d0/gan0**4)

      parameter (aftt =ct2/(12.d0*pi))
      parameter (aft  =ct48*amp)

      parameter (afcst=aftt*amp2)
      parameter (afcs =4.d0*aft)
      parameter (afct=aft/gan0**4)
