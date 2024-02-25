c===============================================================
c===============================================================
      subroutine wfat_2a1e_intgl(rtdb, geom, nqn, nid, mid, lmax_,
     &     gq, zi, celm, sphi, mch, nr, nt, ke, vcf, erg, chmmt, chpos,
     &     pottype, nel, potparam, wfat1, nlgl, nr_m, nt_m, nr_p, nt_p, 
     &     cfx, cfy, cfz, fedvr, nfe, ndvr, maxn, xq, wq, ismask, 
     &     mask_param, g_il, g_j, g_k)

      implicit none
#include "util.fh"
#include "rtdb.fh"
#include "stdio.fh"      
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "geom.fh"
      integer rtdb                    ! [input],
      integer geom                    ! [input],
      integer nqn                     ! [input],
      integer nid(nqn)                ! [input],
      integer mid(nqn)                ! [input],
      integer lmax_                   ! [input],
      integer gq                      ! [input],
      double precision zi             ! [input],
      integer celm                    ! [input],
      integer sphi                    ! [input],
      integer mch                     ! [input],
      integer nr                      ! [input],
      integer nt                      ! [input],
      double precision ke             ! [input],
      double precision vcf(nr*nt)     ! [input],
      double precision erg            ! [input],
      double precision chmmt(3)       ! [input],
      double precision chpos(3)       ! [input],
      character*3 pottype             ! [input]
      integer nel                     ! not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)    ! [input], 1)
      logical wfat1                   ! [input],
      integer nlgl                    ! [input], only used when wfat1=.true.
      integer nr_m                    ! [input], only used when wfat1=.true.
      integer nt_m                    ! [input], only used when wfat1=.true.
      integer nr_p                    ! [input], only used when wfat1=.true.
      integer nt_p                    ! [input], only used when wfat1=.true.
      double precision cfx(1:*)       ! [input], only used when wfat1=.true.
      double precision cfy(1:*)       ! [input], only used when wfat1=.true.
cOLD      double precision cfz(nr*nt) ! [input], only used when wfat1=.true.
      double precision cfz(1:*)       ! [input], only used when wfat1=.true.
      logical fedvr                   ! [input],
      integer nfe                     ! [input],
      integer ndvr(*)                 ! [input], size is nfe
      integer maxn                    ! [input],
      double precision xq(*)          ! [input], size is maxn*nfe
      double precision wq(*)          ! [input], size is maxn*nfe
      logical ismask                  ! [input]
      double precision mask_param(2)  ! [input]. Only used when ismask=.true.. Must always be supplied regardless.
      integer g_il                    ! [output],
      integer g_j                     ! [output], only calculated when wfat1=.true.
      integer g_k(2)                  ! [output], only calculated when wfat1=.true.

      character*12 dt_str
      character*16 tag
      logical t1
      integer ndid, natom, max_npt, nscr, nsh, nom, ni0, ni_, nil
      integer nij0, nij_, nij, nik1_0, nik1_, nik1, nik2_0, nik2_, nik2
      integer natomp, nptp, pr_chk(10), g_nptx, nptx, jpt, pr_count, kk
      integer ic, nom1, nom2, i_om1, i_om2, nnnr, nnnt
      integer icent, npt, ipt_i, ipt_f, ipt, ipt0, i_l, i_r
      integer l_pt, l_wt, l_sc, l_om, l_i0, l_i_, l_il, l_icp, l_ipp
      integer l_ij0, l_ij_, l_ij, l_ik1_0, l_ik1_, l_ik1, l_ik2_0
      integer l_ik2_, l_ik2, l_om1, l_om2, l_dr, l_dt
      integer k_pt, k_wt, k_sc, k_om, k_i0, k_i_, k_il, k_icp, k_ipp
      integer k_ij0, k_ij_, k_ij, k_ik1_0, k_ik1_, k_ik1, k_ik2_0
      integer k_ik2_, k_ik2, k_om1, k_om2, k_dr, k_dt
      double precision npt_tot, t_now
      double precision t_last, t_start, celf, atp(3), znuc, v_c(3)
      double precision v_m(3), vc_psi, wtcelf, dt_prog, dt_all, normch
      double precision vc_cpsi(3), ovl1
      double precision polb_
      double complex czero
      integer nch, nch1, imf
      parameter(nch=1, nch1=3, imf=0)
      parameter(czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)))


      call ga_sync()
      ndid = ga_nodeid()


c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
      
      
c     ==== Some tests ====
      call wfat_sometests(rtdb, zi, nch, erg, lmax_, chmmt)


c     ==== Some checks in wfat1 is used ====
      if (wfat1) then
cc         ovl1 = sum(vcf*cfz(1:nr*nt))
cc         nnnr = max(nr*nr, nr_m*nr_m, nr_p*nr_p)
cc         nnnt = max(nt*nt, nt_m*nt_m, nt_p*nt_p)
cc         if (.not. ma_alloc_get(MT_DBL, nnnr, 'dr_dvr', l_dr, k_dr))
cc     &        call errquit('WFAT: Cannot allocate dr_dvr', 22, MA_ERR)
cc         if (.not. ma_alloc_get(MT_DBL, nnnt, 'dt_dvr', l_dt, k_dt))
cc     &        call errquit('WFAT: Cannot allocate dt_dvr', 25, MA_ERR)
cc         call wfat_2a1e_dpo('z', mch, 0, nr, nt, nr, nt, sphi, sphi, ke,
cc     &        vcf, cfz, .false., dbl_mb(k_dr), .false., dbl_mb(k_dt), 
cc     &        .false., polb_)
cc         polb_ = -2 * polb_
cc         write(luout, '("<wfn0|wfn_z>, pol_zz = ", 2es21.10)')
cc     &        ovl1, polb_
cc         call wfat_2a1e_dpo('x', mch+1, -1, nr, nt, nr_p, nt_p, sphi, 
cc     &        sphi, ke, vcf, cfx, .false., dbl_mb(k_dr), .false., 
cc     &        dbl_mb(k_dt), .false., polb_)
cc         polb_ = -2 * polb_
cc         write(luout, '("<wfn0|wfn_x>, pol_xx = ", 2es21.10)')
cc     &        ovl1, polb_
cc         
cc         if (.not. ma_free_heap(l_dr))
cc     &        call errquit('WFAT: Cannot deallocate dr_dvr', 1, MA_ERR)
cc         if (.not. ma_free_heap(l_dt))
cc     &        call errquit('WFAT: Cannot deallocate dt_dvr', 1, MA_ERR)
      endif


c     ==== Allocate grid-related arrays ====
      call wfat_alloc_quadgrid(gq, max_npt, l_pt, k_pt, l_wt, k_wt)


c     ==== Allocate WFAT integral, its integrands, and some other ====
      if (wfat1) then
         if (fedvr) then
            write(luout, '(" nnn = ", 3i4)') nr, nt_m, nt_p
            nscr = maxn**2 + maxn +
     &           (nr + nt + nr*nt) +
     &           (nt_p + nr*nt_p)
            if (mch >= 1) nscr = nscr + (nt_m + nr*nt_m)
         else
            nscr = 
     &           (nr + nt + nr*nt) +
     &           (nr_p + nt_p + nr_p*nt_p)
            if (mch >= 1) nscr = nscr + (nr_m + nt_m + nr_m*nt_m)
         endif
      else
         nscr = 
     &        nr +
     &        nt +
     &        nr*nt
         if (fedvr) nscr = nscr + maxn**2 + maxn
      endif
      if (.not. MA_alloc_get(MT_DBL, nscr, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: scratch allocation failed', 59, MA_ERR)
      call wfat_alloc_intgl(lmax_, nqn, nch, natom, nsh, nom, ni0, ni_,
     &     nil, l_om, k_om, l_i0, k_i0, l_i_, k_i_, l_il, k_il)


c     ==== Obtaining grid ranges for the current process ====
      call wfat_alloc_gridrange(natom, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange(geom, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natomp, nptp)


c     ==== Calculate WFAT(1) Omega ====
      if (wfat1) then
c        nch1 = the number of channels that contains the correction terms (needed for
c               WFAT(1)), it is equal to 3 because the ionized channel has three correction
c               terms associated to the three cartesian components of the dipole operator.
         call wfat_alloc_intgl_wfat1(lmax_, nqn, nch, nch1, natom,
     &        nptp, nsh, nij0, nij_, nij,
     &        l_ij0, k_ij0, l_ij_, k_ij_, l_ij, k_ij,
     &        nom1, nik1_0, nik1_, nik1, l_om1, k_om1, l_ik1_0, k_ik1_0,
     &        l_ik1_, k_ik1_, l_ik1, k_ik1,
     &        nom2, nik2_0, nik2_, nik2, l_om2, k_om2, l_ik2_0, k_ik2_0,
     &        l_ik2_, k_ik2_, l_ik2, k_ik2)
         call wfat_omega_wfat1(geom, nch, nqn, erg, chpos, nid,
     &        mid, lmax_, natomp, nptp, int_mb(k_icp), int_mb(k_ipp),
     &        gq, nlgl, zi, dcpl_mb(k_om1), dcpl_mb(k_om2))
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
         dcpl_mb(k_ij_ : k_ij_+nij_-1) = czero
         dcpl_mb(k_ik1_ : k_ik1_+nik1_-1) = czero
         dcpl_mb(k_ik2_ : k_ik2_+nik2_-1) = czero
      endif

      
c     ==== Loop over atoms ====
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
c           v_c = electron coordinate w.r.t. the current nuclei.
c           v_m = electron coordinate in the molecular frame, which is also the
c                 basis center.

            
c           ==== Calculate vc|psi> ====
            if (wfat1) then
               if (fedvr) then
                  call wfat_2a1efd_vcpsi_ord1(geom, sphi, mch, nfe, 
     &                 ndvr, maxn, xq, wq, nr, nt, nt_m, nt_p, pottype, 
     &                 nel, potparam, zi, vcf, cfx, cfy, cfz, chpos, 
     &                 v_m, dbl_mb(k_sc), vc_psi, vc_cpsi)
               else
                  call wfat_2a1e_vcpsi_ord1(geom, sphi, mch, nr, nt,
     &                 nr_m, nt_m, nr_p, nt_p, ke, pottype, nel, 
     &                 potparam, zi, vcf, cfx, cfy, cfz, chpos, v_m, 
     &                 dbl_mb(k_sc), vc_psi, vc_cpsi)
               endif
            else
               if (fedvr) then
                  call wfat_2a1efd_vcpsi(geom, sphi, mch, nfe, ndvr, 
     &                 maxn, xq, wq, nr, nt, pottype, nel, potparam, zi,
     &                 vcf, chpos, v_m, dbl_mb(k_sc), vc_psi)
               else
                  call wfat_2a1e_vcpsi(geom, sphi, mch, nr, nt, ke, 
     &                 pottype, nel, potparam, zi, vcf, chpos, v_m,
     &                 dbl_mb(k_sc), vc_psi)
               endif
            endif
            
c           ==== Calculate Omega^(n,m)_(l,m') functions ====
            call wfat_omega_1pt(geom, v_m, nch, erg, chpos, nqn, nid, 
     &           mid, lmax_, int(zi), dcpl_mb(k_om:k_om+nom-1))
            
c           ==== Calculate the integrand of I^{(n,m)}_{l,m'} ====
            call wfat_integrand(nch, nqn, nsh, dcpl_mb(k_om), vc_psi, 
     &           dcpl_mb(k_i0))
            if (isnan( abs(dcpl_mb(k_i0)) ))
     &           write(luout, '("isnan", i8, 3es14.6, 6x, 3es21.12E3)')
     &           ic, v_m, dcpl_mb(k_om), vc_psi
            
c           ==== Calculate the integrand of J^{(n,m)x_i}_{l,m'} ====
c           =====           and K^{(n,m)}_{l,m'}                ====
            if (wfat1) then
               i_om1 = k_om1 + (jpt-1)*nsh*nqn*nch
               i_om2 = k_om2 + (jpt-1)*nsh*nqn*nch
               call wfat_integrand_wfat1(nch, nqn, nsh, dcpl_mb(k_om),
     &              dcpl_mb(i_om1), dcpl_mb(i_om2), vc_cpsi, vc_psi, 
     &              dcpl_mb(k_ij0), dcpl_mb(k_ik1_0), dcpl_mb(k_ik2_0))
            endif

c           ==== Apply masks to the integrands ====
            if (ismask) call wfat_mask_intd(geom, natom, nch, 
     &           nqn, nsh, v_m, mask_param, dcpl_mb(k_i0))
            if (ismask .and. wfat1) call wfat_mask_intd(geom, natom, 
     &           3*nch, nqn, nsh, v_m, mask_param, dcpl_mb(k_ij0))
            if (ismask .and. wfat1) then
               call wfat_mask_intd(geom, natom, nch, nqn, nsh,
     &              v_m, mask_param, dcpl_mb(k_ik1_0))
               call wfat_mask_intd(geom, natom, nch, nqn, nsh,
     &              v_m, mask_param, dcpl_mb(k_ik2_0))
            endif

c           ==== Quadrature summation ====
            call wfat_celf_1pt(geom, icent, v_c, celm, celf)
            wtcelf = celf * dbl_mb(k_wt+ipt-1)

            i_l = k_i_ + (icent-1) * nsh*nqn*nch
            i_r = i_l + nsh*nqn*nch - 1
            dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                         dcpl_mb(k_i0:k_i0+ni0-1) *
     &                         wtcelf                        !! 1)
            
            if (wfat1) then
               i_l = k_ij_ + (icent-1) * nsh*nqn*nch1
               i_r = i_l + nsh*nqn*nch1 - 1
               dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                            dcpl_mb(k_ij0:k_ij0+nij0-1) *
     &                            wtcelf

               i_l = k_ik1_ + (icent-1) * nsh*nqn*nch
               i_r = i_l + nsh*nqn*nch - 1
               dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                            dcpl_mb(k_ik1_0:k_ik1_0+nik1_0-1) *
     &                            wtcelf

               i_l = k_ik2_ + (icent-1) * nsh*nqn*nch
               i_r = i_l + nsh*nqn*nch - 1
               dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                            dcpl_mb(k_ik2_0:k_ik2_0+nik2_0-1) *
     &                            wtcelf
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
      if (wfat1)
     &     call wfat_reorder_intgl_wfat1(nsh, natom, nqn, nch,
     &     dcpl_mb(k_ij_), dcpl_mb(k_ij),
     &     dcpl_mb(k_ik1_), dcpl_mb(k_ik1),
     &     dcpl_mb(k_ik2_), dcpl_mb(k_ik2))
      

c     ==== Accumulate to global intgl (g_il) ====
      call ga_acc(g_il, 1, nil, 1, 1, dcpl_mb(k_il), nil, 1.0d0)
      call ga_sync()
      call ga_norm1(g_il, normch)      

      if (wfat1) then
         call ga_acc(g_j, 1, nij, 1, 1, dcpl_mb(k_ij), nij, 1.0d0)
c         call ga_zero(g_j)
         call ga_sync()
         call ga_norm1(g_j, normch)

         call ga_acc(g_k(1), 1, nik1, 1, 1, dcpl_mb(k_ik1), nik1, 1.0d0)
c         call ga_zero(g_k(1))
         call ga_sync()
         call ga_norm1(g_k(1), normch)

         call ga_acc(g_k(2), 1, nik2, 1, 1, dcpl_mb(k_ik2), nik2, 1.0d0)
c         call ga_zero(g_k(2))
         call ga_sync()
         call ga_norm1(g_k(2), normch)
      endif
      

c     ==== Deallocations ====
      call wfat_dealloc_quadgrid(l_pt, l_wt)
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      call wfat_dealloc_intgl(l_om, l_i0, l_i_, l_il)
      if (wfat1)
     &     call wfat_dealloc_intgl_wfat1(l_ij0, l_ij_, l_ij,
     &     l_om1, l_ik1_0, l_ik1_, l_ik1, l_om2, l_ik2_0, l_ik2_, l_ik2)
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
      double precision orbmmt0(3)  ! [input]. Unshifted electronic dipole moment.
      double precision orbcent(3)  ! [output],
      double precision orbmmt(3)   ! [output]. Shifted electronic dipole moment.
      
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
