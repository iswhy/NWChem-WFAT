c===============================================================
c===============================================================
      subroutine wfat_2a1e_task_nopw()

      implicit none
      integer rtdb
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      

c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting nat failed', 94, GEOM_ERR)

      
c     ==== Obtaining info on beta ====
      call wfat_get_beta(rtdb, bt_start, bt_last, nbt, dbt)
      dbtr = dbt * deg2rad


c     ==== Obtaining info on gamma ====
      call wfat_get_gamma(rtdb, gm_start, gm_last, ngm, dgm)
      dgmr = dgm * deg2rad

      
      nintg = nat * nqn
      if (.not. ga_create(mt_dcpl, nintg, 1, 'intgl', 0, 0, g_il))
     &     call errquit('WFAT: Cannot create intgl global array.',
     &     77, GA_ERR)

      
c     ==== Calculate and print structure factor and ion yield ====
      do ig = 1, ngm
         gm = gm_start + dgm*(ig-1)
         do ib = 1, nbt
            bt = bt_start + dbt*(ib-1)
            write(un, '(2es12.4)', advance='no') bt, gm
            
            call ga_zero(g_il)
            call wfat_2a1e_intgl_nopw(rtdb, geom, nqn, nid, mid, bt, gm,
     &           gq, zi, celm, sphi, mch, nr, nt, ke, orbid, vcf, erg, 
     &           chmmt0, eps, wfat1, chpos, chmmt, g_il)

            call ga_get(g_il, 1, nintg, 1, 1, dcpl_mb(k_in), nintg)
            call wfat_2a1e_srfc_nopw(bt, gm, erg, chmmt, nat, nqn,
     &           dcpl_mb(k_in), dcpl_mb(k_gn))

            do i = 1, nqn
               i_gn = k_gn + i - 1
               write(un, '(5x, 3es20.8)', advance='no')
     &              real(dcpl_mb(i_gn)), imag(dcpl_mb(i_gn)),
     &              (abs(dcpl_mb(i_gn)))**2
            enddo
            write(un, '("")')
            
         enddo
         write(un, '("")')
      enddo


      if (.not. ga_destroy(g_il))
     &     call errquit('WFAT: Cannot destroy intgl.', 72, GA_ERR)
      
      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_2a1e_intgl_nopw(rtdb, geom, nqn, nid, mid, bt, gm,
     &     gq, zi, celm, sphi, mch, nr, nt, ke, orbid, vcf, erg, chmmt0,
     &     eps, wfat1,
     &     chpos, chmmt, g_il)

      implicit none
#include "util.fh"
#include "rtdb.fh"
#include "stdio.fh"      
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "geom.fh"
      integer rtdb                 ! [input],
      integer geom                 ! [input],
      integer nqn                  ! [input],
      integer nid(nqn)             ! [input],
      integer mid(nqn)             ! [input],
      double precision bt          ! [input], beta orientation angle in degrees.
      double precision gm          ! [input], gamma orientation angle in degrees.
      integer gq                   ! [input],
      double precision zi          ! [input],
      integer celm                 ! [input],
      integer sphi                 ! [input],
      integer mch                  ! [input],
      integer nr                   ! [input],
      integer nt                   ! [input],
      double precision ke          ! [input],
      integer orbid                ! [input],
      double precision vcf(nr*nt)  ! [input],
      double precision erg         ! [input],
      double precision chmmt0(3)   ! [input],
      double precision eps         ! [input],
      logical wfat1                ! [input],
clater      integer nlgl                 ! [input], only used when wfat1=.true.
clater      integer nr_m                 ! [input], only used when wfat1=.true.
clater      integer nt_m                 ! [input], only used when wfat1=.true.
clater      integer nr_p                 ! [input], only used when wfat1=.true.
clater      integer nt_p                 ! [input], only used when wfat1=.true.
clater      double precision cfx(1:*)    ! [input], only used when wfat1=.true.
clater      double precision cfy(1:*)    ! [input], only used when wfat1=.true.
clater      double precision cfz(nr*nt)  ! [input], only used when wfat1=.true.
      double precision chpos(3)    ! [output],
      double precision chmmt(3)    ! [output],
      integer g_il                 ! [output],
clater      integer g_j                  ! [output], only calculated when wfat1=.true.
clater      integer g_k(2)               ! [output], only calculated when wfat1=.true.

      character*12 dt_str
      character*16 tag
      logical ismask_g, ismask_j, ismask_k, t1
      integer ndid, natom, max_npt, nscr, ni0, ni_, nil
      integer nij0, nij_, nij, nik1_0, nik1_, nik1, nik2_0, nik2_, nik2
      integer natomp, nptp, pr_chk(10), g_nptx, nptx, jpt, pr_count, kk
      integer ic
      integer icent, npt, ipt_i, ipt_f, ipt, ipt0, i_l, i_r, i_tu0
      integer l_pt, l_wt, l_sc, l_i0, l_i_, l_il, l_icp, l_ipp, l_tu0
      integer l_ij0, l_ij_, l_ij, l_ik1_0, l_ik1_, l_ik1, l_ik2_0
      integer l_ik2_, l_ik2, l_dr, l_dt
      integer k_pt, k_wt, k_sc, k_i0, k_i_, k_il, k_icp, k_ipp, k_tu0
      integer k_ij0, k_ij_, k_ij, k_ik1_0, k_ik1_, k_ik1, k_ik2_0
      integer k_ik2_, k_ik2, k_dr, k_dt
      double precision mask_g(2), mask_j(2), mask_k(2), npt_tot, t_now
      double precision t_last, t_start, celf, atp(3), znuc, v_c(3)
      double precision v_m(3), vc_psi, wtcelf, dt_prog, dt_all, normch
      double precision iondip(3)
      double complex czero
      integer nch, nch1, imf, nsh
      parameter(nch=1, nch1=3, nsh=1, imf=0)
      parameter(czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)))
      

      call ga_sync()
      ndid = ga_nodeid()


c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)


      ismask_g = .false.
      if (rtdb_get(rtdb, 'wfat:mask_g', mt_dbl, 2, mask_g))
     &     ismask_g = .true.
clater      ismask_j = .false.
clater      if (wfat1 .and. rtdb_get(rtdb, 'wfat:mask_j', mt_dbl, 2, mask_j))
clater     &     ismask_j = .true.
clater      ismask_k = .false.
clater      if (wfat1 .and. rtdb_get(rtdb, 'wfat:mask_k', mt_dbl, 2, mask_k))
clater     &     ismask_k = .true.

      if (ismask_g) write(luout, '("Masking g")')
clater      if (ismask_j) write(luout, '("Masking j")')
clater      if (ismask_k) write(luout, '("Masking k")')
      

c     ==== Obtain properties of the ionized channel ====
c     nch1 = the number of channels that contains the correction terms (needed for WFAT(1)),
c            it is equal to 3 because the ionized channel has three correction terms 
c            associated to the three cartesian components of the dipole operator.
      call wfat_2a1e_orbcent(geom, chmmt0, chpos, chmmt)
      iondip = (/ 0.0d0, 0.0d0, 0.0d0 /)
      call wfat_print_orb_property(nch, orbid, 1, erg, chpos, chmmt,
     &     iondip)


c     ==== Some checks in wfat1 is used ====
      if (wfat1) then
clater         ovl1 = sum(vcf*cfz)
clater         if (.not. ma_alloc_get(MT_DBL, nr*nr, 'dr_dvr', l_dr, k_dr))
clater     &        call errquit('WFAT: Cannot allocate dr_dvr', 22, MA_ERR)
clater         if (.not. ma_alloc_get(MT_DBL, nt*nt, 'dt_dvr', l_dt, k_dt))
clater     &        call errquit('WFAT: Cannot allocate dt_dvr', 25, MA_ERR)
clater         call wfat_2a1e_dpo('z', mch, 0, nr, nt, nr, nt, sphi, sphi, ke,
clater     &        vcf, cfz, .false., dbl_mb(k_dr), .false., dbl_mb(k_dt), 
clater     &        .false., polb_)
clater         polb_ = -2 * polb_
clater         write(luout, '("<wfn0|wfn_z>, pol_zz = ", 2es21.10)')
clater     &        ovl1, polb_
clater         if (.not. ma_free_heap(l_dr))
clater     &        call errquit('WFAT: Cannot deallocate dr_dvr', 1, MA_ERR)
clater         if (.not. ma_free_heap(l_dt))
clater     &        call errquit('WFAT: Cannot deallocate dt_dvr', 1, MA_ERR)
      endif

      
c     ==== Obtaining grid ranges for the current process ====
      call wfat_alloc_quadgrid(gq, max_npt, l_pt, k_pt, l_wt, k_wt)
      call wfat_alloc_gridrange(natom, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange(geom, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natomp, nptp)

      
c     ==== Allocate WFAT integral, its integrands, and some other ====
      if (wfat1) then
clater         nscr = 
clater     &        nr +
clater     &        nt +
clater     &        nr*nt +
clater     &        nr_m +
clater     &        nt_m +
clater     &        nr_m*nt_m +
clater     &        nr_p +
clater     &        nt_p +
clater     &        nr_p*nt_p
      else
         nscr = 
     &        nr +
     &        nt +
     &        nr*nt
      endif
      if (.not. MA_alloc_get(MT_DBL, nscr, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: scratch allocation failed', 59, MA_ERR)
      call wfat_alloc_intgl_nopw(nqn, nch, natom, ni0, ni_, nil, l_i0,
     &     k_i0, l_i_, k_i_, l_il, k_il)
      if (.not. MA_alloc_get(MT_DBL, nqn*nptp, 'tau0', l_tu0, k_tu0))
     &     call errquit('WFAT: tau0 allocation failed', 59, MA_ERR)
      call wfat_collect_tau0(geom, nch, nqn, erg, chpos, nid, mid, bt,
     &     gm, natomp, nptp, int_mb(k_icp), int_mb(k_ipp), gq, zi,
     &     dcpl_mb(k_tu0))
            

c     ==== Calculate WFAT(1) Omega ====
      if (wfat1) then
clater         call wfat_alloc_intgl_wfat1_nopw(nqn, nch, nch1, natom,
clater     &        nptp, nij0, nij_, nij,
clater     &        l_ij0, k_ij0, l_ij_, k_ij_, l_ij, k_ij,
clater     &        nik1_0, nik1_, nik1, l_ik1_0, k_ik1_0,
clater     &        l_ik1_, k_ik1_, l_ik1, k_ik1,
clater     &        nik2_0, nik2_, nik2, l_ik2_0, k_ik2_0,
clater     &        l_ik2_, k_ik2_, l_ik2, k_ik2)
clater         if (.not. MA_alloc_get(MT_DBL, nptp, 'tau1', l_tu1, k_tu1))
clater     &     call errquit('WFAT: tau1 allocation failed', 59, MA_ERR)
clater         call wfat_collect_tau1()
      endif      


      if (ndid == 0) then
         write(luout, '(2/)')
         call util_print_centered(luout, 'Computation the WFAT '//
     &        'integrals', 30, .true.)
         write(luout, '("")')
         write(luout, 901)
         write(luout, '("")')
      endif
 901  format(
     &     5x, 'The calculation of the WFAT integrals has started. ',
     &     'If this is a ''start'' type job, this part can take ', /,
     &     5x, 'up to several minutes, so please wait. The progress ',
     &     'of this calculation will be shown below.')


      
c     ==== BEGIN QUADRATURE SUMMATION ====
c     ==== Loop over atoms ====
      call ga_sync()
      call wfat_get_intgl_chkpt(nptp, pr_chk)
      call wfat_init_progress(nptp, g_nptx, npt_tot)
      if (ndid == 0) write(luout, '(5x, "The total number of '//
     &     'quadrature grid points is ", i9, ".")') int(npt_tot)

      t_last = util_cpusec()
      t_start = t_last
      jpt = 0
      pr_count = 0
      kk = 1
      dcpl_mb(k_i_ : k_i_+ni_-1) = czero
      if (wfat1) then
clater         dcpl_mb(k_ij_ : k_ij_+nij_-1) = czero
clater         dcpl_mb(k_ik1_ : k_ik1_+nik1_-1) = czero
clater         dcpl_mb(k_ik2_ : k_ik2_+nik2_-1) = czero
      endif
      
      do ic = 1, natomp
         icent = int_mb(k_icp + ic - 1)
         if (.not. geom_cent_get(geom, icent, tag, atp, znuc))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)

         call wfat_get_gridpt(int(znuc), gq, max_npt, npt,
     &                        dbl_mb(k_pt : k_pt+3*max_npt-1),
     &                        dbl_mb(k_wt : k_wt+max_npt-1))

         ipt_i = int_mb(k_ipp + (ic-1)*2)
         ipt_f = int_mb(k_ipp + 1 + (ic-1)*2)
c        ==== Loop over all quadrature points around the current atom ====
         do ipt = ipt_i, ipt_f
            jpt = jpt + 1
            ipt0 = k_pt + 3*(ipt-1)
            v_c(1:3) = dbl_mb(ipt0 : ipt0+2)
            call geom_cent2cent(geom, icent, v_c, imf, v_m, t1)

            
c           ==== Calculate vc|psi> ====
            if (wfat1) then
clater               call wfat_2a1e_vcpsi_ord1(geom, sphi, mch, nr, nt, nr_m,
clater     &              nt_m, nr_p, nt_p, ke, eps, zi, vcf, cfx, cfy, cfz, 
clater     &              chpos, v_m, dbl_mb(k_sc), vc_psi, vc_cpsi)
            else
               call wfat_2a1e_vcpsi(geom, sphi, mch, nr, nt, ke, eps, 
     &              zi, vcf, chpos, v_m, dbl_mb(k_sc), vc_psi)
            endif

            
c           ==== Calculate the integrand of I^{(n,m)}_{l,m'} ====
            i_tu0 = k_tu0 + (jpt-1)*nsh*nqn*nch
            call wfat_integrand(nch, nqn, nsh, dcpl_mb(i_tu0), vc_psi, 
     &           dcpl_mb(k_i0))
            if (isnan( abs(dcpl_mb(k_i0)) ))
     &           write(luout, '("isnan", i8, 3es14.6, 6x, 3es21.12E3)')
     &           ic, v_m, dcpl_mb(k_tu0), vc_psi
            
c           ==== Calculate the integrand of J^{(n,m)x_i}_{l,m'} ====
c           =====           and K^{(n,m)}_{l,m'}                ====
            if (wfat1) then
clater               i_om1 = k_om1 + (jpt-1)*nsh*nqn*nch
clater               i_om2 = k_om2 + (jpt-1)*nsh*nqn*nch
clater               call wfat_integrand_wfat1(nch, nqn, nsh, dcpl_mb(k_om),
clater     &              dcpl_mb(i_om1), dcpl_mb(i_om2), vc_cpsi, vc_psi, 
clater     &              dcpl_mb(k_ij0), dcpl_mb(k_ik1_0), dcpl_mb(k_ik2_0))
            endif

c           ==== Apply masks to the integrands ====
            if (ismask_g) call wfat_mask_intd(geom, natom, nch, erg, 
     &           nqn, nsh, v_m, mask_g, dcpl_mb(k_i0))
clater            if (ismask_j .and. wfat1) call wfat_mask_intd(geom, natom, 
clater     &           3*nch, erg, nqn, nsh, v_m, mask_j, dcpl_mb(k_ij0))
clater            if (ismask_k .and. wfat1) then
clater               call wfat_mask_intd(geom, natom, nch, erg, nqn, nsh,
clater     &              v_m, mask_k, dcpl_mb(k_ik1_0))
clater               call wfat_mask_intd(geom, natom, nch, erg, nqn, nsh,
clater     &              v_m, mask_k, dcpl_mb(k_ik2_0))
clater            endif

c           ==== Quadrature summation ====
            call wfat_celf_1pt(geom, icent, v_c, celm, celf)
            wtcelf = celf * dbl_mb(k_wt+ipt-1)

            i_l = k_i_ + (icent-1) * nsh*nqn*nch
            i_r = i_l + nsh*nqn*nch - 1
            dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                         dcpl_mb(k_i0:k_i0+ni0-1) *
     &                         wtcelf                        !! 1)
            
            if (wfat1) then
clater               i_l = k_ij_ + (icent-1) * nsh*nqn*nch1
clater               i_r = i_l + nsh*nqn*nch1 - 1
clater               dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
clater     &                            dcpl_mb(k_ij0:k_ij0+nij0-1) *
clater     &                            wtcelf
clater
clater               i_l = k_ik1_ + (icent-1) * nsh*nqn*nch
clater               i_r = i_l + nsh*nqn*nch - 1
clater               dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
clater     &                            dcpl_mb(k_ik1_0:k_ik1_0+nik1_0-1) *
clater     &                            wtcelf
clater
clater               i_l = k_ik2_ + (icent-1) * nsh*nqn*nch
clater               i_r = i_l + nsh*nqn*nch - 1
clater               dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
clater     &                            dcpl_mb(k_ik2_0:k_ik2_0+nik2_0-1) *
clater     &                            wtcelf
            endif
            
c            jpt = jpt + 1
c            write(luout, '("check: ", i8, 9es15.6)')
c     &           jpt, v_m, vc_psi, dcpl_mb(k_i0), dcpl_mb(k_i_), wtcelf

            
c           ==== Print progress ====
            pr_count = pr_count + 1
            if (pr_count == pr_chk(kk)) then
               t_now = util_cpusec()
               dt_prog = t_now - t_last
               t_last = t_now
               call wfat_intgl_progress(kk, pr_count, dt_prog, npt_tot, 
     &              g_nptx, nptx)
               kk = kk + 1
            endif
            
         enddo
      enddo
      dt_all = t_now - t_start
      write(dt_str, '(f12.2)') dt_all
      if (ndid == 0) write(luout, '(5x, "The computation of the WFAT '//
     &     'integrals takes ", a, " s.")') trim(adjustl(dt_str))

      
c     ==== Reorder elements of intgl ====
      call wfat_reorder_intgl(nsh, natom, nqn, nch,
     &     dcpl_mb(k_i_ : k_i_+ni_-1), dcpl_mb(k_il : k_il+nil-1))
clater      if (wfat1)
clater     &     call wfat_reorder_intgl_wfat1(nsh, natom, nqn, nch,
clater     &     dcpl_mb(k_ij_), dcpl_mb(k_ij),
clater     &     dcpl_mb(k_ik1_), dcpl_mb(k_ik1),
clater     &     dcpl_mb(k_ik2_), dcpl_mb(k_ik2))
      
      
c     ==== Accumulate to global intgl (g_il) ====
      call ga_acc(g_il, 1, nil, 1, 1, dcpl_mb(k_il), nil, 1.0d0)
      call ga_sync()
      call ga_norm1(g_il, normch)      

clater      if (wfat1) then
clater         call ga_acc(g_j, 1, nij, 1, 1, dcpl_mb(k_ij), nij, 1.0d0)
claterc         call ga_zero(g_j)
clater         call ga_sync()
clater         call ga_norm1(g_j, normch)
clater         if (ndid==0) write(luout, '("intgl_j norm = ", es14.6)') normch
clater
clater         call ga_acc(g_k(1), 1, nik1, 1, 1, dcpl_mb(k_ik1), nik1, 1.0d0)
claterc         call ga_zero(g_k(1))
clater         call ga_sync()
clater         call ga_norm1(g_k(1), normch)
clater         if (ndid==0) write(luout, '("intgl_k1 norm = ", es14.6)')
clater     &        normch
clater
clater         call ga_acc(g_k(2), 1, nik2, 1, 1, dcpl_mb(k_ik2), nik2, 1.0d0)
claterc         call ga_zero(g_k(2))
clater         call ga_sync()
clater         call ga_norm1(g_k(2), normch)
clater         if (ndid==0) write(luout, '("intgl_k2 norm = ", es14.6)')
clater     &        normch
clater      endif      


c     ==== Deallocations ====
      call wfat_dealloc_quadgrid(l_pt, l_wt)
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      call wfat_dealloc_intgl_nopw(l_i0, l_i_, l_il)
clater      if (wfat1)
clater     &     call wfat_dealloc_intgl_wfat1(l_ij0, l_ij_, l_ij,
clater     &     l_om1, l_ik1_0, l_ik1_, l_ik1, l_om2, l_ik2_0, l_ik2_, l_ik2)
      call wfat_terminate_progress(g_nptx)

      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_alloc_intgl_nopw(nqn, nch, nat, ni0, ni_, nil,
     &     l_i0, k_i0, l_i_, k_i_, l_il, k_il)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nqn                       ! [input],
      integer nch                       ! [input],
      integer nat                       ! [input],
      integer ni0, ni_, nil             ! [output],
      integer l_i0, k_i0                ! [output],
      integer l_i_, k_i_                ! [output],
      integer l_il, k_il                ! [output],


      ni0 = nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, ni0, 'intgl0', l_i0, k_i0))
     &     call errquit('WFAT: intgl0 allocation failed', 59, MA_ERR)
      ni_ = nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, ni_, 'intgl_', l_i_, k_i_))
     &     call errquit('WFAT: intgl_ allocation failed', 59, MA_ERR)
      nil = nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nil, 'intgl', l_il, k_il))
     &     call errquit('WFAT: intgl allocation failed', 59, MA_ERR)

      end
c===============================================================
      

c=======================================================
c=======================================================
      subroutine wfat_dealloc_intgl_nopw(l_i0, l_i_, l_il)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_i0, l_i_, l_il     ! [input]
      
      if (.not. MA_free_heap(l_i0))
     &     call errquit('WFAT: intgl0 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_i_))
     &     call errquit('WFAT: intgl_ deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_il))
     &     call errquit('WFAT: intgl deallocation failed', 61, MA_ERR)

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_2a1e_srfc_nopw(bt, gm, erg, chmmt, nat, nqn,
     &     intgl, gn)

      implicit none
      double precision bt              ! [input], beta orientation angle in degrees.
      double precision gm              ! [input], gamma orientation angle in degrees.
      double precision erg
      double precision chmmt(3)
      integer nat
      integer nqn
      double complex intgl(nat, nqn)
      double complex gn(nqn)

      integer i
      double precision kp, mu(3), protmat(3,3)

      kp = sqrt(2 * abs(erg))
      call wfat_passrot(bt, gm, protmat)
      mu = matmul(transpose(protmat), -chmmt)

      do i = 1, nqn
         gn(i) = sum(intgl(:,i)) * exp(-kp*mu(3))
      enddo

      end
c=======================================================
      
