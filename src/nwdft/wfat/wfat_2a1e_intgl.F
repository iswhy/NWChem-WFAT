c===============================================================
c===============================================================
      subroutine wfat_2a1e_intgl(rtdb, geom, nqn, nid, mid, lmax_,
     &     gq, zi, celm, mch, nr, nt, ke, orbid, vcf, erg, chmmt0, eps, 
     &     wfat1, chpos, chmmt, g_il)

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
      integer lmax_                ! [input],
      integer gq                   ! [input],
      double precision zi          ! [input],
      integer celm                 ! [input],
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
      double precision chpos(3)    ! [output],
      double precision chmmt(3)    ! [output],
      integer g_il                 ! [output],

      character*12 dt_str
      character*16 tag
      logical ismask_g, t1
      integer ndid, natom, max_npt, nscr, nsh, nom, ni0, ni_, nil
      integer natomp, nptp, pr_chk(10), g_nptx, nptx, jpt, pr_count, kk
      integer ic
      integer icent, npt, ipt_i, ipt_f, ipt, ipt0, i_l, i_r
      integer l_pt, l_wt, l_sc, l_om, l_i0, l_i_, l_il, l_icp, l_ipp
      integer k_pt, k_wt, k_sc, k_om, k_i0, k_i_, k_il, k_icp, k_ipp
      double precision mask_g(2), npt_tot, t_now, t_last, t_start, celf
      double precision atp(3), znuc, v_c(3), v_m(3), vc_psi, wtcelf
      double precision dt_prog, dt_all, normch, iondip(3)
      double complex czero
      integer nch, imf
      parameter(nch=1, imf=0)
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


c     ==== Obtain properties of each channel ====
c     nch = the number of the ionized channels.
c     nch1 = the number of channels that contains the correction terms (needed for WFAT(1)),
c            it is equal to 3*nch because each ionized channel has three correction terms 
c            associated to the three cartesian components of the dipole operator.
clater      if (wfat1) then
clater         nch1 = 3*nch
clater      else
clater         nch1 = 0
clater      endif
cOLD      ncf = nbas * (nch + nch1)
cOLD      if (.not. MA_alloc_get(MT_DBL, ncf, 'vch', l_cf, k_cf))
cOLD     &     call errquit('WFAT: vch allocation failed', 59, MA_ERR)
      call wfat_2a1e_orbcent(geom, chmmt0, chpos, chmmt)
      iondip = (/ 0.0d0, 0.0d0, 0.0d0 /)
      call wfat_print_orb_property(nch, orbid, 1, erg, chpos, chmmt,
     &     iondip)
clater      if (wfat1) then
clater         call wfat_perturbnew_cwfn(nch, spinid, orbid,
clater     &        dbl_mb(k_cf+nbas*nch), pb)
clater         call wfat_check_wfat1ovl(nch, dbl_mb(k_cf),
clater     &        dbl_mb(k_cf+nbas*nch))
clater      endif
cOLD      if (.not. MA_alloc_get(MT_INT, nch+nch1, 'spinid1', l_si, k_si))
cOLD     &     call errquit('WFAT: spinid1 allocation failed', 59, MA_ERR)
cOLD      if (.not. MA_alloc_get(MT_DBL, 3*(nch+nch1), 'chpos1',
cOLD     &     l_pos, k_pos))
cOLD     &     call errquit('WFAT: chpos1 allocation failed', 59, MA_ERR)
cOLD      call wfat_get_wfat1_sppos(wfat1, nch, spinid, chpos, int_mb(k_si),
cOLD     &     dbl_mb(k_pos))      


c     ==== Allocate grid-related arrays ====
      call wfat_alloc_quadgrid(gq, max_npt, l_pt, k_pt, l_wt, k_wt)


c     ==== Allocate WFAT integral, its integrands, and some other ====
      nscr = nr + nt + nr*nt
      if (.not. MA_alloc_get(MT_DBL, nscr, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: scratch allocation failed', 59, MA_ERR)
      call wfat_alloc_intgl(lmax_, nqn, nch, natom, nsh, nom, ni0, ni_,
     &     nil, l_om, k_om, l_i0, k_i0, l_i_, k_i_, l_il, k_il)


c     ==== Obtaining grid ranges for the current process ====
      call wfat_alloc_gridrange(natom, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange(geom, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natomp, nptp)


c     ==== Calculate WFAT(1) Omega ====
clater      if (wfat1) then
clater         call wfat_alloc_intgl_wfat1(lmax_, nqn, nch, nch1, natom,
clater     &        nptp, nsh, nij0, nij_, nij,
clater     &        l_ij0, k_ij0, l_ij_, k_ij_, l_ij, k_ij,
clater     &        nom1, nik1_0, nik1_, nik1, l_om1, k_om1, l_ik1_0, k_ik1_0,
clater     &        l_ik1_, k_ik1_, l_ik1, k_ik1,
clater     &        nom2, nik2_0, nik2_, nik2, l_om2, k_om2, l_ik2_0, k_ik2_0,
clater     &        l_ik2_, k_ik2_, l_ik2, k_ik2)
clater         call wfat_omega_wfat1(geom, nch, nqn, erg, chpos, nid,
clater     &        mid, lmax_, natomp, nptp, int_mb(k_icp), int_mb(k_ipp),
clater     &        gq, nlgl, zion, dcpl_mb(k_om1), dcpl_mb(k_om2))
clater      endif

      
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
clater      if (wfat1) then
clater         dcpl_mb(k_ij_ : k_ij_+nij_-1) = czero
clater         dcpl_mb(k_ik1_ : k_ik1_+nik1_-1) = czero
clater         dcpl_mb(k_ik2_ : k_ik2_+nik2_-1) = czero
clater      endif
      
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
            ipt0 = k_pt + 3*(ipt-1)
            v_c(1:3) = dbl_mb(ipt0 : ipt0+2)
            call geom_cent2cent(geom, icent, v_c, imf, v_m, t1)

            
c           ==== Calculate vc|psi> ====
            if (wfat1) then
               call wfat_2a1e_vcpsi_ord1()
            else
               call wfat_2a1e_vcpsi(geom, mch, nr, nt, ke, eps, zi, vcf,
     &              chpos, v_m, dbl_mb(k_sc), vc_psi)
            endif
            
c           ==== Calculate Omega^(n,m)_(l,m') functions ====
            call wfat_omega_1pt(geom, v_m, nch, erg,
     &           chpos, nqn, nid, mid, lmax_, int(zi),
     &           dcpl_mb(k_om:k_om+nom-1))
            
c           ==== Calculate the integrand of I^{(n,m)}_{l,m'} ====
            call wfat_integrand(nch, nqn, nsh, dcpl_mb(k_om),
     &           vc_psi, dcpl_mb(k_i0))
            if (isnan( abs(dcpl_mb(k_i0)) ))
     &           write(luout, '("isnan", i8, 3es14.6, 6x, 3es21.12E3)')
     &           ic, v_m, dcpl_mb(k_om), vc_psi
            
            

clater            if (wfat1) then
clater               jpt = jpt + 1
clater               i_om1 = k_om1 + (jpt-1)*nsh*nqn*nch
clater               i_om2 = k_om2 + (jpt-1)*nsh*nqn*nch
claterc              ==== Calculate the integrand of J^{(n,m)x_i}_{l,m'} ====
claterc              =====           and K^{(n,m)}_{l,m'}                ====
clater               call wfat_integrand_wfat1(nch, nqn, nsh, dcpl_mb(k_om),
clater     &              dcpl_mb(i_om1), dcpl_mb(i_om2), dbl_mb(k_vc+nch),
clater     &              dbl_mb(k_vc), dcpl_mb(k_ij0), dcpl_mb(k_ik1_0),
clater     &              dcpl_mb(k_ik2_0))
clater            endif

            
            if (ismask_g) call wfat_mask_intd(geom, natom, nch, erg, 
     &           nqn, nsh, v_m, mask_g, dcpl_mb(k_i0))
clater            if (ismask_j .and. wfat1) call wfat_mask_intd(geom, 
clater     &           natom, 3*nch, erg, nqn, nsh, v_m, mask_j,
clater     &           dcpl_mb(k_ij0))
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

clater            if (wfat1) then
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
clater            endif

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


c     ==== Deallocations ====
      call wfat_dealloc_quadgrid(l_pt, l_wt)
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      call wfat_dealloc_intgl(l_om, l_i0, l_i_, l_il)
      call wfat_terminate_progress(g_nptx)
      
      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_orbcent(geom, orbmmt0, orbcent, orbmmt)

      implicit none
#include "geom.fh"
#include "errquit.fh"
      integer geom                 ! [input],
      double precision orbmmt0(3)  ! [input],
      double precision orbcent(3)  ! [output],
      double precision orbmmt(3)   ! [output],
      
      character*16 t
      integer i, nc
      double precision r(3), q, totzn


      if (.not. geom_ncent(geom, nc))
     &     call errquit('WFAT: Cannot get ncenters.', 40, GEOM_ERR)
      
      orbcent = 0.0d0
      do i = 1, nc
         if (.not. geom_cent_get(geom, i, t, r, q))
     &        call errquit('WFAT: geom center info failed', 88,
     &        GEOM_ERR)
         orbcent = orbcent + q * r
      enddo

      if (.not. geom_nuc_charge(geom, totzn))
     &     call errquit('WFAT: Getting totzn failed', 5, GEOM_ERR)
      orbcent = orbcent / totzn
      orbmmt = orbmmt0 - orbcent
      
      end
c===============================================================
      
      
