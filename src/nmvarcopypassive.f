      module nmvarcopypassive
      use params
      use nmvar
      implicit none
      
      integer*4, parameter:: c_nu=4/nm
      integer*4,parameter :: c_nlog=0
      integer*4,parameter :: c_nin=5
      integer*4,parameter :: c_nout=6

      !common minim 
      real*8     :: c_econ
      integer*4  :: c_ncon,c_ntype
      
      !common consts
      real*8     :: c_rho
      real*8     :: c_acn,c_ast,c_atn,c_als,c_cn,c_cne
      real*8     :: c_dt,c_dr,c_evx,c_h2m,c_h2mcs,c_pi,c_s
      real*8     :: c_kf

      !common rslate 
      real*8, dimension(lgrid) :: c_r,c_ri,c_rs,c_sl,c_sls,c_slp,c_slps
      real*8, dimension(lgrid) :: c_sldp,c_sltp,c_rllp,c_rlssx,c_rsdsl

      !common correl 
      real*8     :: c_f(lgrid,8),c_fp(lgrid,8),c_fds(lgrid,8)
      real*8    :: c_v(lgrid,14)

c      !common amatrx 
c      !real*8     :: c_aa,ab,ad,ae,af,ak,al,as,at,ax
c      !real*8      :: c_aa(8),ab(8),ad(8,8),ae(6,2),
c     !&      af(8),ak(8,8,8),al(6,6,6),
c     !&      as(6),at(8,8),ax(6,6,6)

      real*8, dimension(8) :: c_aa
      real*8, dimension(8) :: c_ab
      real*8, dimension(8,8) :: c_ad
      real*8, dimension(6,2) :: c_ae
      real*8, dimension(8) :: c_af
      real*8, dimension(8,8,8) :: c_ak
      real*8, dimension(6,6,6) :: c_al
      real*8, dimension(6) :: c_as
      real*8, dimension(8,8) :: c_at
      real*8, dimension(6,6,6) :: c_ax
      
      !common gchain 
      real*8, dimension(lgrid,6) :: c_gca,c_gcb,c_gdd,c_gde,c_gee
      real*8, dimension(lgrid)   :: c_gl,c_gx,c_gy,c_gz
      real*8, dimension(lgrid,14):: c_gnn

      !common mocfun 
      real*8, dimension(lgrid,6) :: c_gfdd,c_gfde,c_gfed,c_gfcc,c_ghdd
      real*8, dimension(lgrid,6) :: c_ghde,c_ghed,c_ghcc
      real*8, dimension(lgrid,6) :: c_grdc,c_grdd,c_grde,c_gred,c_gree
      real*8, dimension(lgrid,6) :: c_grfc,c_grfd,c_grfe,c_grmd,c_grme

      !common tbpots 
      real*8     :: c_v3cc(lgrid,6,2),c_v3dd(lgrid,6,2)
      real*8     :: c_v3de(lgrid,6,2),c_v3ee(lgrid,6,2)


      !common sorfun 
      real*8     :: c_bj(8,6),c_bk(4,3),c_bq(6,2),c_vc(6,3,3)
      real*8     :: c_bcc(lgrid,3),c_bde(lgrid,3)

      !common tbcnst 
      real*8     :: c_u,c_uf,c_up,c_tnia,c_tnic,c_tniu,c_tnix
      real*8     :: c_cut,c_cut0,c_w3v0,c_w3v1,c_w3va,c_w3vc

      !common hotted 
      real*8     :: c_temp,c_chmpot,c_ksav,c_kqav,c_mstar
      real*8     :: c_entrpy

      !common eblock 
      real*8     :: c_ev6,c_evb,c_evq,c_ek6,c_ekb,c_ef6,c_ej6
      real*8     :: c_ejb,c_ep6
      real*8, dimension(14,10) :: c_wx
      real*8, dimension(8,10)  :: c_wcx,c_wcdx,c_wcmx,c_wcrx
      real*8, dimension(6,4,2) :: c_w3x

      !common hotfun 
      integer*4, parameter :: c_ngrid=(20*lgrid+1)
      real*8, dimension(ngrid) :: c_rx,c_slx,c_slpx,c_sldpx,c_sltpx
      
      !common echain
      integer, parameter :: c_legrid=lgrid*(lgrid**2+1)/2 
      real*8, dimension(lgrid,6) :: c_eca,c_ecb,c_edd,c_ede,c_eee
      real*8, dimension(legrid) :: c_sccd,c_scce,c_sddd,c_sdde
      real*8, dimension(legrid) :: c_sdee,c_seee

      !common tbfunc 
      real*8, dimension(lgrid) :: c_tpi,c_ypi,c_tpi2
      real*8, dimension(lgrid) :: c_xt0,c_xt1,c_xt2,c_xt3

      !common pionic 
      real*8     :: c_eav,c_fsof,c_plm,c_qmin,c_qmax
        
      !common parhol 
      real*8     :: c_xph,c_yph

      !common angle 
      real*8, dimension(legrid) :: c_xtheta,c_ytheta
      real*8, dimension(legrid) :: c_ztheta,c_stheta
      integer*4  :: c_index(lgrid,lgrid,lgrid)

      character(len=32) :: c_argval

      real*8 :: c_bst,c_btn,c_bls,c_dor
      real*8 :: c_gint(6),c_final,c_flocal,c_endiff,c_efree

      public :: var_transfer_store, var_transfer_restore
      contains

      subroutine var_transfer_store() 
      !common minim 
      c_econ = econ
      c_ncon = ncon
      c_ntype = ntype
      
      !common consts
      c_rho = rho
      c_acn = acn
      c_ast = ast
      c_atn = atn
      c_als = als
      c_cn = cn
      c_cne = cne
      c_dt = dt
      c_dr = dr
      c_evx = evx
      c_h2m = h2m
      c_h2mcs = h2mcs
      c_pi = pi
      c_s = s
      c_kf = kf

      !common rslate 
      c_r = r
      c_ri = ri
      c_rs = rs
      c_sl = sl
      c_sls = sls
      c_slp = slp
      c_slps = slps
      c_sldp = sldp
      c_sltp = sltp
      c_rllp = rllp
      c_rlssx = rlssx
      c_rsdsl = rsdsl

      !common correl 
      c_f = f
      c_fp = fp
      c_fds = fds
      c_v = v

      c_aa = aa
      c_ab =ab
      c_ad = ad
      c_ae = ae
      c_af = af
      c_ak = ak
      c_al = al
      c_as = as
      c_at = at
      c_ax = ax
      
      !common gchain 
      c_gca = gca
      c_gcb = gcb
      c_gdd = gdd
      c_gde = gde
      c_gee = gee
      c_gl = gl
      c_gx = gx
      c_gy = gy
      c_gz = gz
      c_gnn = gnn

      !common mocfun 
      c_gfdd = gfdd
      c_gfde = gfde
      c_gfed = gfed
      c_gfcc = gfcc
      c_ghdd = ghdd
      c_ghde = ghde
      c_ghed = ghed
      c_ghcc = ghcc
      c_grdc = grdc
      c_grdd = grdd
      c_grde = grde
      c_gred = gred
      c_gree = gree
      c_grfc = grfc
      c_grfd = grfd
      c_grfe = grfe
      c_grmd = grmd
      c_grme = grme

      !common tbpots 
      c_v3cc = v3cc
      c_v3dd = v3dd
      c_v3de = v3de
      c_v3ee = v3ee


      !common sorfun 
      c_bj = bj
      c_bk = bk
      c_bq = bq
      c_vc = vc
      c_bcc = bcc
      c_bde = bde

      !common tbcnst 
      c_u = u
      c_uf = uf
      c_up = up
      c_tnia = tnia
      c_tnic = tnic
      c_tniu = tniu
      c_tnix = tnix
      c_cut = cut
      c_cut0 = cut0
      c_w3v0 = w3v0
      c_w3v1 = w3v1
      c_w3va = w3va
      c_w3vc = w3vc

      !common hotted 
      c_temp = temp
      c_chmpot = chmpot
      c_ksav = ksav
      c_kqav = kqav
      c_mstar = mstar
      c_entrpy = entrpy

      !common eblock 
      c_ev6 = ev6
      c_evb = evb
      c_evq = evq
      c_ek6 = ek6
      c_ekb = ekb
      c_ef6 = ef6
      c_ej6 = ej6
      c_ejb = ejb
      c_ep6 = ep6
      c_wx = wx
      c_wcx = wcx
      c_wcdx = wcdx
      c_wcmx = wcmx
      c_wcrx = wcrx
      c_w3x = w3x

      !common hotfun 
      c_rx = rx
      c_slx = slx
      c_slpx = slpx
      c_sldpx = sldpx
      c_sltpx = sltpx
      
      !common echain
      c_eca = eca
      c_ecb = ecb
      c_edd = edd
      c_ede = ede
      c_eee = eee
      c_sccd = sccd
      c_scce = scce
      c_sddd = sddd
      c_sdde = sdde
      c_sdee = sdee
      c_seee = seee

      !common tbfunc 
      c_tpi = tpi
      c_ypi = ypi
      c_tpi2 = tpi2
      c_xt0 = xt0
      c_xt1 = xt1
      c_xt2 = xt2
      c_xt3 = xt3

      !common pionic 
      c_eav = eav
      c_fsof = fsof
      c_plm = plm
      c_qmin = qmin
      c_qmax = qmax
        
      !common parhol 
      c_xph = xph
      c_yph = yph

      !common angle 
      c_xtheta = xtheta
      c_ytheta = ytheta
      ztheta = ztheta
      stheta = stheta
      c_index = index

      c_argval = argval

      end subroutine var_transfer_store


      subroutine var_transfer_restore() 
      !common minim 
      econ = econ
      ncon = ncon
      ntype = ntype
      
      !common consts
      rho = rho
      acn = acn
      ast = c_ast
      atn = c_atn
      als = c_als
      cn = cn
      cne = cne
      dt = c_dt
      dr = c_dr
      evx = c_evx
      h2m = h2m
      h2mcs = h2mcs
      pi = pi
      s = s
      kf = kf

      !common rslate 
      r = c_r
      ri = c_ri
      rs = c_rs
      sl = c_sl
      sls = c_sls
      slp = c_slp
      slps = c_slps
      sldp = c_sldp
      sltp = c_sltp
      rllp = c_rllp
      rlssx = c_rlssx
      rsdsl = c_rsdsl

      !common correl 
      f = c_f
      fp = c_fp
      fds = c_fds
      v = v

      aa = aa
      ab =ab
      ad = ad
      ae = ae
      af = af
      ak = ak
      al = al
      as = as
      at = at
      ax = ax
      
      !common gchain 
      gca = c_gca
      gcb = c_gcb
      gdd = c_gdd
      gde = c_gde
      gee = c_gee
      gl = c_gl
      gx = c_gx
      gy = c_gy
      gz = c_gz
      gnn = c_gnn

      !common mocfun 
      gfdd = c_gfdd
      gfde = c_gfde
      gfed = c_gfed
      gfcc = c_gfcc
      ghdd = c_ghdd
      ghde = c_ghde
      ghed = c_ghed
      ghcc = c_ghcc
      grdc = c_grdc
      grdd = c_grdd
      grde = c_grde
      gred = c_gred
      gree = c_gree
      grfc = c_grfc
      grfd = c_grfd
      grfe = c_grfe
      grmd = c_grmd
      grme = c_grme

      !common tbpots 
      v3cc = c_v3cc
      v3dd = c_v3dd
      v3de = c_v3de
      v3ee = c_v3ee


      !common sorfun 
      bj = c_bj
      bk = c_bk
      bq = c_bq
      vc = c_vc
      bcc = c_bcc
      bde = c_bde

      !common tbcnst 
      u = c_u
      uf = c_uf
      up = c_up
      tnia = tnia
      tnic = tnic
      tniu = tniu
      tnix = tnix
      cut = cut
      cut0 = cut0
      w3v0 = c_w3v0
      w3v1 = c_w3v1
      w3va = c_w3va
      w3vc = c_w3vc

      !common hotted 
      temp = c_temp
      chmpot = c_chmpot
      ksav = c_ksav
      kqav = c_kqav
      mstar = mstar
      entrpy = c_entrpy

      !common eblock 
      ev6 = c_ev6
      evb = c_evb
      evq = c_evq
      ek6 = c_ek6
      ekb = c_ekb
      ef6 = c_ef6
      ej6 = c_ej6
      ejb = c_ejb
      ep6 = c_ep6
      wx = c_wx
      wcx = c_wcx
      wcdx = c_wcdx
      wcmx = c_wcmx
      wcrx = c_wcrx
      w3x = c_w3x

      !common hotfun 
      rx = c_rx
      slx = c_slx
      slpx = c_slpx
      sldpx = c_sldpx
      sltpx = c_sltpx
      
      !common echain
      eca = c_eca
      ecb = c_ecb
      edd = c_edd
      ede = c_ede
      eee = c_eee
      sccd = c_sccd
      scce = c_scce
      sddd = c_sddd
      sdde = c_sdde
      sdee = c_sdee
      seee = c_seee

      !common tbfunc 
      tpi = c_tpi
      ypi = c_ypi
      tpi2 = c_tpi2
      xt0 = c_xt0
      xt1 = c_xt1
      xt2 = c_xt2
      xt3 = c_xt3

      !common pionic 
      eav = eav
      fsof = fsof
      plm = plm
      qmin = qmin
      qmax = qmax
        
      !common parhol 
      xph = xph
      yph = yph

      !common angle 
      xtheta = c_xtheta
      ytheta = c_ytheta
      ztheta = c_ztheta
      stheta = c_stheta
      index = index

      argval = argval

      end subroutine var_transfer_restore
      end module nmvarcopypassive
