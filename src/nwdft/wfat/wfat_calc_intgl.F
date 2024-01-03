c=======================================================
c=======================================================
c==== Imam ====
      subroutine wfat_calc_intgl(rtdb, geom, bas, nch, spinid, nqn, nid,
     &                           mid, lmax_, gqual, zi, celm, smh, 
     &                           res_mode, wfat1, nlgl, nb, ns, rdm, 
     &                           vch, erg, chpos, chmmt, iondip, g_il, 
     &                           g_j, g_k)
      
      implicit none
#include "util.fh"
#include "rtdb.fh"
#include "stdio.fh"      
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "errquit.fh"
#include "geom.fh"
      integer rtdb                    ! [input]
      integer geom                    ! [input]
      integer bas                     ! [input]
      integer nch                     ! [input]
      integer spinid(nch)             ! [input], The spin of each column of vch.
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      integer lmax_                   ! [input]
      integer gqual                   ! [input]
      double precision zi             ! [input]
      integer celm                    ! [input]
      double precision smh            ! [input]
      logical res_mode                ! [input]
      logical wfat1                   ! [input]
      integer nlgl                    ! [input], only used when wfat1=.true.
      integer nb                      ! [input]
      integer ns                      ! [input]
      double precision rdm(nb,nb,ns)  ! 
      double precision vch(nb,*)      ! [input], the minimum number of columns is nch if wfat1 is .false., otherwise 4*nch.
      double precision erg(nch)       ! [input]
      double precision chpos(3, nch)  ! [input]
      double precision chmmt(3, nch)  ! [input]
      double precision iondip(3, nch) ! [input]
      integer g_il                    ! [output], global array handle for intgl. For structure factors.
      integer g_j                     ! [output], only calculated when wfat1=.true.
      integer g_k(2)                  ! [output], only calculated when wfat1=.true.

      character*12 dt_str
      character*16 tag
      logical t1, ismask_gjk, ismask_g, ismask_j, ismask_k
      integer natom, nrd, max_npt, nvc, nsh, nom, ni0, ni_
      integer nch1, npt, natomp, nptp, ncpos, nscr, nil, nmo
      integer nij0, nij_, nij
      integer nom1, nik1_0, nik1_, nik1, nom2, nik2_0, nik2_, nik2
      integer i, j, ic, icent, ipt, ipt_i, ipt_f
      integer jpt, ich, iqn, ish, i_vc, ndid
      integer l_pt, l_wt, l_vc, l_om, l_i0
      integer l_i_, l_sc, l_il, l_icp, l_ipp, l_si, l_pos
      integer k_pt, k_wt, k_vc, k_om, k_i0
      integer k_i_, k_sc, k_il, k_icp, k_ipp, k_si, k_pos
      integer l_ij0, l_ij_, l_ij, l_om1, l_ik1_0, l_ik1_, l_ik1, l_om2
      integer l_ik2_0, l_ik2_, l_ik2
      integer k_ij0, k_ij_, k_ij, k_om1, k_ik1_0, k_ik1_, k_ik1, k_om2
      integer k_ik2_0, k_ik2_, k_ik2
      integer is, io, ipt0, i0l, i0r, oml, omr, i_l, i_r, imf, i_rd
      integer g_rd(2), i_mo, i_om1, i_om2, spinid1(nch+3*nch), pr_count
      integer pr_chk(10), kk, g_nptx, g_npt_tot
      double precision atp(3), znuc, v_c(3), v_m(3), mask_gjk(2) 
      double precision mask_g(2), mask_j(2), mask_k(2), normch, celf
      double precision wtcelf, numnel, t_now, t_last, t_start, dt_all
      double precision dt_prog, nptx, npt_tot, kp
      double complex czero
      
      parameter (imf = 0)
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )



      call ga_sync()
      ndid = ga_nodeid()

      
c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)


      ismask_gjk = .false.
      if (rtdb_get(rtdb, 'wfat:mask_gjk', mt_dbl, 2, mask_gjk))
     &     ismask_gjk = .true.
      
      ismask_g = .false.
      if (rtdb_get(rtdb, 'wfat:mask_g', mt_dbl, 2, mask_g))
     &     ismask_g = .true.
      ismask_j = .false.
      if (wfat1 .and. rtdb_get(rtdb, 'wfat:mask_j', mt_dbl, 2, mask_j))
     &     ismask_j = .true.
      ismask_k = .false.
      if (wfat1 .and. rtdb_get(rtdb, 'wfat:mask_k', mt_dbl, 2, mask_k))
     &     ismask_k = .true.

      if (ismask_gjk) write(luout, '("Masking gjk")')
      if (ismask_g) write(luout, '("Masking g")')
      if (ismask_j) write(luout, '("Masking j")')
      if (ismask_k) write(luout, '("Masking k")')

      
      
c     ==== Obtain properties of each channel ====
c     nch = the number of the ionized channels.
c     nch1 = the number of channels that contains the correction terms (needed for WFAT(1)),
c            it is equal to 3*nch because each ionized channel has three correction terms 
c            associated to the three cartesian components of the dipole operator.
      if (wfat1) then
         nch1 = 3*nch
      else
         nch1 = 0
      endif

      if (wfat1 .and. ndid == 0) then
         call wfat_check_wfat1ovl(nch, vch(:,1:nch), 
     &        vch(:,nch+1:4*nch))
      endif
      if (ga_nodeid() == 0) then
         write(luout, '("GTO coefficients:")')
         call wfat_dprintmat(nb, nch+nch1, 'f', vch)
      endif
      call ga_sync()
      if (.not. MA_alloc_get(MT_INT, nch+nch1, 'spinid1', l_si, k_si))
     &     call errquit('WFAT: spinid1 allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, 3*(nch+nch1), 'chpos1',
     &     l_pos, k_pos))
     &     call errquit('WFAT: chpos1 allocation failed', 59, MA_ERR)
      call wfat_get_wfat1_sppos(wfat1, nch, spinid, chpos, int_mb(k_si),
     &     dbl_mb(k_pos))
c     NOTES:
c     k_pos points to a conceptually 2D array where the first nch columns are
c     identical to the columns of chpos. If wfat1=.true., the next 3*nch columns
c     are such that every group of three columns are identical to the corresponding
c     column among the first nch columns. For example, columns within the second
c     group, i.e. the (nch+4)-th up to (nch+6)-th columns, are duplicates of the
c     second column.
      
      
c     ==== Some analyses and tests (optional) ====
      call wfat_sometests(rtdb, zi, nch, erg, lmax_, chmmt)
      call wfat_mulliken(rtdb, geom, bas, nb, ns, rdm, nch, vch, spinid,
     &                   chpos)


c     ==== Allocate grid-related arrays ====
      call wfat_alloc_quadgrid(gqual, max_npt, l_pt, k_pt, l_wt, k_wt)


c     ==== Allocate WFAT integral, its integrands, and some other ====
      nscr = 21*nb + 4*ns + nb * (3*nb+1+ns) + 4*(nch+nch1)
      nvc = nch + nch1
      if (.not. MA_alloc_get(MT_DBL, nscr, 'vscr', l_sc, k_sc))
     &     call errquit('WFAT: vscr allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nvc, 'vcpsi', l_vc, k_vc))
     &     call errquit('WFAT: vcpsi allocation failed', 59, MA_ERR)
      call wfat_alloc_intgl(lmax_, nqn, nch, natom, nsh, nom, ni0, ni_,
     &     nil, l_om, k_om, l_i0, k_i0, l_i_, k_i_, l_il, k_il)
      

c     ==== Obtaining grid ranges for the current process ====
      call wfat_alloc_gridrange(natom, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange(geom, ndid, gqual, int_mb(k_icp),
     &                        int_mb(k_ipp), natomp, nptp)
c      write(luout, '("node:", i3, "  atom1,pt1:", i3, i7, '//
c     &     '"  atom2,pt2:", i3, i7, "  npoint:", i7)')
c     &     ndid, icent_i, ipt_i, icent_f, ipt_f, nptp


c     ==== Calculate WFAT(1) Omega ====
      if (wfat1) then
         call wfat_alloc_intgl_wfat1(lmax_, nqn, nch, nch1, natom,
     &        nptp, nsh, nij0, nij_, nij,
     &        l_ij0, k_ij0, l_ij_, k_ij_, l_ij, k_ij,
     &        nom1, nik1_0, nik1_, nik1, l_om1, k_om1, l_ik1_0, k_ik1_0,
     &        l_ik1_, k_ik1_, l_ik1, k_ik1,
     &        nom2, nik2_0, nik2_, nik2, l_om2, k_om2, l_ik2_0, k_ik2_0,
     &        l_ik2_, k_ik2_, l_ik2, k_ik2)
         call wfat_omega_wfat1(geom, nch, nqn, erg, chpos, nid,
     &        mid, lmax_, natomp, nptp, int_mb(k_icp), int_mb(k_ipp),
     &        gqual, nlgl, zi, dcpl_mb(k_om1), dcpl_mb(k_om2))
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
         dcpl_mb(k_ij_ : k_ij_+nij_-1) = czero
         dcpl_mb(k_ik1_ : k_ik1_+nik1_-1) = czero
         dcpl_mb(k_ik2_ : k_ik2_+nik2_-1) = czero
      endif
      
      do ic = 1, natomp
         icent = int_mb(k_icp + ic - 1)
         if (.not. geom_cent_get(geom, icent, tag, atp, znuc))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)

         call wfat_get_gridpt(int(znuc), gqual, max_npt, npt,
     &                        dbl_mb(k_pt : k_pt+3*max_npt-1),
     &                        dbl_mb(k_wt : k_wt+max_npt-1))

         ipt_i = int_mb(k_ipp + (ic-1)*2)
         ipt_f = int_mb(k_ipp + 1 + (ic-1)*2)
c        ==== Loop over all quadrature points around the current atom ====
         do ipt = ipt_i, ipt_f
            jpt = jpt + 1             ! Needed for wfat1.
            ipt0 = k_pt + 3*(ipt-1)
            v_c(1:3) = dbl_mb(ipt0 : ipt0+2)
            call geom_cent2cent(geom, icent, v_c, imf, v_m, t1)
cc            write(luout, '("        v_m = ", 3e14.5)') v_m

            
c           ==== Calculate vc|psi> ====
            call wfat_dft_vcpsi_1pt(geom, bas, nb, ns, nch+nch1, vch, 
     &           int_mb(k_si), dbl_mb(k_pos), v_m, rdm, zi, smh,
     &           dbl_mb(k_sc:k_sc+nscr-1), ismask_gjk, erg(1), 
     &           mask_gjk(1), mask_gjk(2), dbl_mb(k_vc:k_vc+nvc-1))

            
c           ==== Calculate Omega^(n,m)_(l,m') functions ====
            call wfat_omega_1pt(geom, v_m, nch, erg,
     &           chpos, nqn, nid, mid, lmax_, int(zi),
     &           dcpl_mb(k_om:k_om+nom-1))
            
c           ==== Calculate the integrand of I^{(n,m)}_{l,m'} ====
            call wfat_integrand(nch, nqn, nsh, dcpl_mb(k_om),
     &           dbl_mb(k_vc), dcpl_mb(k_i0))
            if (isnan( abs(dcpl_mb(k_i0)) ))
     &           write(luout, '("isnan", i8, 3es14.6, 6x, 3es21.12E3)')
     &           ic, v_m, dcpl_mb(k_om), dbl_mb(k_vc)
            
            
            if (wfat1) then
cfordebug               write(luout, '("wfat1: ", i8, 4es15.6)')
cfordebug     &              jpt, dbl_mb(k_vc), dbl_mb(k_vc+1), dbl_mb(k_vc+2),
cfordebug     &              dbl_mb(k_vc+3)
               
               i_om1 = k_om1 + (jpt-1)*nsh*nqn*nch
               i_om2 = k_om2 + (jpt-1)*nsh*nqn*nch
c              ==== Calculate the integrand of J^{(n,m)x_i}_{l,m'} ====
c              =====           and K^{(n,m)}_{l,m'}                ====
               call wfat_integrand_wfat1(nch, nqn, nsh, dcpl_mb(k_om),
     &              dcpl_mb(i_om1), dcpl_mb(i_om2), dbl_mb(k_vc+nch),
     &              dbl_mb(k_vc), dcpl_mb(k_ij0), dcpl_mb(k_ik1_0),
     &              dcpl_mb(k_ik2_0))
            endif

            
            if (ismask_g) call wfat_mask_intd(geom, natom, nch,
     &           nqn, nsh, v_m, mask_g, dcpl_mb(k_i0))
            if (ismask_j .and. wfat1) call wfat_mask_intd(geom, 
     &           natom, 3*nch, nqn, nsh, v_m, mask_j,
     &           dcpl_mb(k_ij0))
            if (ismask_k .and. wfat1) then
               call wfat_mask_intd(geom, natom, nch, nqn, nsh,
     &              v_m, mask_k, dcpl_mb(k_ik1_0))
               call wfat_mask_intd(geom, natom, nch, nqn, nsh,
     &              v_m, mask_k, dcpl_mb(k_ik2_0))
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

cfordebug            if (ga_nodeid() == 0)
cfordebug     &           write(luout, '("@grid: ", i8, 2es20.10, '//
cfordebug     &           '3es20.10)')
cfordebug     &           jpt, dcpl_mb(k_om), dbl_mb(k_vc), celf,
cfordebug     &           dbl_mb(k_wt+ipt-1)
            

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
      

cforcheck      write(luout, '("sum CHAN_INT = ", 2e14.5)')
cforcheck     &     sum(dcpl_mb(k_il : k_il+nil-1))
cforcheck      write(luout, '("max of real part = ", e14.5)')
cforcheck     &     maxval(abs(real(dcpl_mb(k_il : k_il+nil-1))))
cforcheck      write(luout, '("max of imag part = ", e14.5)')
cforcheck     &     maxval(abs(imag(dcpl_mb(k_il : k_il+nil-1))))

      
c     ==== Accumulate to global intgl (g_il) ====
      call ga_acc(g_il, 1, nil, 1, 1, dcpl_mb(k_il), nil, 1.0d0)
      call ga_sync()
      call ga_norm1(g_il, normch)
cfordebug      if (ndid==0) write(luout, '("intgl norm = ", es14.6)') normch

      if (wfat1) then
         call ga_acc(g_j, 1, nij, 1, 1, dcpl_mb(k_ij), nij, 1.0d0)
c         call ga_zero(g_j)
         call ga_sync()
         call ga_norm1(g_j, normch)
         if (ndid==0) write(luout, '("intgl_j norm = ", es14.6)') normch

         call ga_acc(g_k(1), 1, nik1, 1, 1, dcpl_mb(k_ik1), nik1, 1.0d0)
c         call ga_zero(g_k(1))
         call ga_sync()
         call ga_norm1(g_k(1), normch)
         if (ndid==0) write(luout, '("intgl_k1 norm = ", es14.6)')
     &        normch

         call ga_acc(g_k(2), 1, nik2, 1, 1, dcpl_mb(k_ik2), nik2, 1.0d0)
c         call ga_zero(g_k(2))
         call ga_sync()
         call ga_norm1(g_k(2), normch)
         if (ndid==0) write(luout, '("intgl_k2 norm = ", es14.6)')
     &        normch
      endif

      
c     ==== Deallocations ====
      call wfat_dealloc_quadgrid(l_pt, l_wt)
      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: vscr deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_vc))
     &     call errquit('WFAT: vc|psi> deallocation failed', 61, MA_ERR)
      call wfat_dealloc_intgl(l_om, l_i0, l_i_, l_il)
      if (wfat1)
     &     call wfat_dealloc_intgl_wfat1(l_ij0, l_ij_, l_ij,
     &     l_om1, l_ik1_0, l_ik1_, l_ik1, l_om2, l_ik2_0, l_ik2_, l_ik2)
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      if (.not. MA_free_heap(l_si))
     &     call errquit('WFAT: spinid1 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_pos))
     &     call errquit('WFAT: chpos1 deallocation failed', 61, MA_ERR)
      call wfat_terminate_progress(g_nptx)
      
      
      end subroutine wfat_calc_intgl
!! 1) The quadrature summation is performed simultaneously for all channels (orbitals),
!!    each of which is given by 3rd index of wf_intgl. It's possible to do this because
!!    the cell function, the quadrature points and weights are the same for all channels.
!!    Note that it may not be possible to perform this simultaneous summation if each of
!!    the channels is assigned their own set of quadrature points. In the current
!!    implementation, all channels share the same set of quadrature points because these
!!    points (and weights) are atom-specific, not channel-specific.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_integrand(nch, nqn, nsh, omg, vc_psi, intgl0)

      implicit none
      integer nch
      integer nqn
      integer nsh
      double complex omg(nsh, nqn, nch)
      double precision vc_psi(nch)
      double complex intgl0(nsh, nqn, nch)

      integer i

      do i = 1, nch
         intgl0(:,:,i) = conjg(omg(:,:,i)) * vc_psi(i)
      enddo

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_integrand_wfat1(nch, nqn, nsh, omg, omg1, omg2,
     &           vxi, vc_psi, intd_j, intd_k1, intd_k2)

      implicit none
#include "stdio.fh"
      integer nch
      integer nqn
      integer nsh
      double complex omg(nsh, nqn, nch)
      double complex omg1(nsh, nqn, nch)
      double complex omg2(nsh, nqn, nch)
      double precision vxi(3, nch)
      double precision vc_psi(nch)
      double complex intd_j(nsh, nqn, 3, nch)
      double complex intd_k1(nsh, nqn, nch)
      double complex intd_k2(nsh, nqn, nch)

      integer i, j

      do i = 1, nch
         intd_k1(:,:,i) = conjg(omg1(:,:,i)) * vc_psi(i)
         intd_k2(:,:,i) = conjg(omg2(:,:,i)) * vc_psi(i)
         do j = 1, 3
            intd_j(:,:,j,i) = conjg(omg(:,:,i)) * vxi(j,i)
         enddo
      enddo

      
cfordebug      write(luout, '(4x, "intgls : ", 5es15.7)')
cfordebug     &     abs(intd_j(2,1,1,1)),
cfordebug     &     abs(intd_j(2,1,2,1)),
cfordebug     &     abs(intd_j(3,1,3,1)),
cfordebug     &     abs(intd_k1(1,1,1)),
cfordebug     &     abs(intd_k2(1,1,1))
      
      end
c=======================================================

      

      
cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_mask_gjk(geom, natom, nch, erg, nqn, nsh, v_m,
cOLD     &     a, s, ef, intd_g, intd_j, intd_k1, intd_k2)
cOLD
cOLD      implicit none
cOLD      integer geom
cOLD      integer natom                           ! [input], the total number of atoms that corresponds to the geometry handle 'geom'.
cOLD      integer nch
cOLD      double precision erg(nch)
cOLD      integer nqn
cOLD      integer nsh
cOLD      double precision v_m(3)                 ! [input], the Cartesian coordinate of the point in molecular-frame
cOLD      double precision a
cOLD      double precision s
cOLD      double precision ef
cOLD      double complex intd_g(nsh, nqn, nch)
cOLD      double complex intd_j(nsh, nqn, 3, nch)
cOLD      double complex intd_k1(nsh, nqn, nch)
cOLD      double complex intd_k2(nsh, nqn, nch)
cOLD
cOLD      logical t1
cOLD      integer i, j, imf
cOLD      double precision v_c(3), r(natom), r_, mask
cOLD      parameter( imf=0 )
cOLD      
cOLD      do i = 1, nch
cOLD         r = 0.0d0
cOLD         do j = 1, natom
cOLD            call geom_cent2cent(geom, imf, v_m, j, v_c, t1)
cOLD            r(j) = sqrt( sum(v_c**2) )
cOLD         enddo
cOLD         r_ = minval(r)
cOLD         mask = 0.5d0 * erfc(a*(r_-s))
cOLD         intd_g(:,:,i) = intd_g(:,:,i) * mask
cOLD         intd_j(:,:,:,i) = intd_j(:,:,:,i) * mask
cOLD         intd_k1(:,:,i) = intd_k1(:,:,i) * mask
cOLD         intd_k2(:,:,i) = intd_k2(:,:,i) * mask
cOLD      enddo
cOLD
cOLD      end
cOLDc=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_reorder_intgl(nsh, ncent, nqn, nch, intgl_, intgl)

      implicit none
      integer nsh
      integer ncent
      integer nqn
      integer nch
      double complex intgl_(nsh, nqn, nch, ncent)
      double complex intgl(nsh, ncent, nqn, nch)

      integer ich, iqn, icent, i_lm
      
      
      do ich = 1, nch
         do iqn = 1, nqn
            do icent = 1, ncent
               do i_lm = 1, nsh
                  intgl(i_lm, icent, iqn, ich) =
     &                     intgl_(i_lm, iqn, ich, icent)
               enddo
            enddo
         enddo
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_reorder_intgl_wfat1(nsh, ncent, nqn, nch,
     &           intl_j_, intl_j, intl_k1_, intl_k1, intl_k2_, intl_k2)

      implicit none
      integer nsh
      integer ncent
      integer nqn
      integer nch
      double complex intl_j_(nsh, nqn, 3, nch, ncent)
      double complex intl_j(nsh, ncent, nqn, 3, nch)
      double complex intl_k1_(nsh, nqn, nch, ncent)
      double complex intl_k1(nsh, ncent, nqn, nch)
      double complex intl_k2_(nsh, nqn, nch, ncent)
      double complex intl_k2(nsh, ncent, nqn, nch)

      integer ich, i, iqn, icent, i_lm
      
      
      do ich = 1, nch
         do i = 1, 3
            do iqn = 1, nqn
               do icent = 1, ncent
                  do i_lm = 1, nsh
                     intl_j(i_lm, icent, iqn, i, ich) =
     &                    intl_j_(i_lm, iqn, i, ich, icent)
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ich = 1, nch
         do iqn = 1, nqn
            do icent = 1, ncent
               do i_lm = 1, nsh
                  intl_k1(i_lm, icent, iqn, ich) =
     &                 intl_k1_(i_lm, iqn, ich, icent)
                  intl_k2(i_lm, icent, iqn, ich) =
     &                 intl_k2_(i_lm, iqn, ich, icent)
               enddo
            enddo
         enddo
      enddo

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_calc_asysphere(geom, rasfac, asc, ras)

      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "stdio.fh"
#include "global.fh"      
      integer geom                 ! [input], geometry handle.
      double precision rasfac      ! [input], factor for determining ras below
      double precision asc(3)      ! [output], coordinate of the center of asymptotic sphere, in Bohr.
      double precision ras         ! [output], radius of asymptotic sphere, in Bohr

      character*16 tag
      integer nat, i
      double precision zn, znt, atp(3), d0, d


      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)

      asc = 0.0d0
      znt = 0.0d0
      do i = 1, nat
         if (.not. geom_cent_get(geom, i, tag, atp, zn))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
         asc = asc + zn*atp
         znt = znt + zn
      enddo
      asc = asc / znt


      d = 0.0d0
      do i = 1, nat
         if (.not. geom_cent_get(geom, i, tag, atp, zn))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
         d0 = sqrt(sum((atp-asc)**2))
         if (d0 > d) d = d0
      enddo
      d = d + 1.0d0          ! 1)
      ras = rasfac * d

      if (ga_nodeid() == 0) write(luout, '(" mol. rad : ", es14.6)') d
      end
c     NOTES:
C       1) d is at least 1.0 Bohr. In atom, it is exactly 1.0 Bohr.
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_alloc_gridrange(nat, l_icp, k_icp, l_ipp, k_ipp)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nat
      integer l_icp, k_icp
      integer l_ipp, k_ipp
      
      
      if (.not. ma_alloc_get(MT_INT, nat, 'icent_p', l_icp, k_icp))
     &     call errquit('WFAT: icent_p allocation failed.', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_INT, 2*nat, 'ipt_p', l_ipp, k_ipp))
     &     call errquit('WFAT: ipt_p allocation failed.', 1, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dealloc_gridrange(l_icp, l_ipp)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_icp, l_ipp

      if(.not. MA_free_heap(l_icp))
     &     call errquit('WFAT: icent_p deallocation failed', 61, MA_ERR)
      if(.not. MA_free_heap(l_ipp))
     &     call errquit('WFAT: ipt_p deallocation failed', 61, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_alloc_quadgrid(gq, max_npt, l_pt, k_pt,
     &     l_wt, k_wt)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer gq          ! [input],
      integer max_npt     ! [output],
      integer l_pt, k_pt  ! [output],
      integer l_wt, k_wt  ! [output],


      call wfat_get_gridmaxpt(gq, max_npt)
      if (.not. MA_alloc_get(MT_DBL, 3*max_npt, 'gpt', l_pt, k_pt))
     &     call errquit('WFAT: gpt allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, max_npt, 'gwt', l_wt, k_wt))
     &     call errquit('WFAT: gwt allocation failed', 59, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dealloc_quadgrid(l_pt, l_wt)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_pt, l_wt       ! [input]

      if (.not. MA_free_heap(l_pt))
     &     call errquit('WFAT: gpt deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_wt))
     &     call errquit('WFAT: gwt deallocation failed', 61, MA_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_alloc_intgl(lmax_, nqn, nch, nat, nsh, nom, ni0, 
     &     ni_, nil, l_om, k_om, l_i0, k_i0, l_i_, k_i_, l_il, k_il)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer lmax_                     ! [input],
      integer nqn                       ! [input],
      integer nch                       ! [input],
      integer nat                       ! [input],
      integer nsh, nom, ni0, ni_, nil   ! [output],
      integer l_om, k_om                ! [output],
      integer l_i0, k_i0                ! [output],
      integer l_i_, k_i_                ! [output],
      integer l_il, k_il                ! [output],


      nsh = (1 + lmax_)**2
      nom = nsh * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nom, 'omega', l_om, k_om))
     &     call errquit('WFAT: omega allocation failed', 59, MA_ERR)
      ni0 = nom
      if (.not. MA_alloc_get(MT_DCPL, ni0, 'intgl0', l_i0, k_i0))
     &     call errquit('WFAT: intgl0 allocation failed', 59, MA_ERR)
      ni_ = nsh * nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, ni_, 'intgl_', l_i_, k_i_))
     &     call errquit('WFAT: intgl_ allocation failed', 59, MA_ERR)
      nil = nsh * nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nil, 'intgl', l_il, k_il))
     &     call errquit('WFAT: intgl allocation failed', 59, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_alloc_intgl_wfat1(lmax_, nqn, nch, nch1, nat, 
     &     nptp, nsh, nij0, nij_, nij,
     &     l_ij0, k_ij0, l_ij_, k_ij_, l_ij, k_ij,
     &     nom1, nik1_0, nik1_, nik1, l_om1, k_om1, l_ik1_0, k_ik1_0, 
     &     l_ik1_, k_ik1_, l_ik1, k_ik1,
     &     nom2, nik2_0, nik2_, nik2, l_om2, k_om2, l_ik2_0, k_ik2_0, 
     &     l_ik2_, k_ik2_, l_ik2, k_ik2)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer lmax_                       ! [input],
      integer nqn                         ! [input],
      integer nch                         ! [input],
      integer nch1                        ! [input],
      integer nat                         ! [input], the total number of atoms in the system.
      integer nptp                        ! [input], the number of points in the current process.
      integer nsh, nij0, nij_, nij        ! [output],
      integer l_ij0, k_ij0                ! [output],
      integer l_ij_, k_ij_                ! [output],
      integer l_ij, k_ij                  ! [output],
      integer nom1, nik1_0, nik1_, nik1   ! [output],
      integer l_om1, k_om1                ! [output],
      integer l_ik1_0, k_ik1_0            ! [output],
      integer l_ik1_, k_ik1_              ! [output],
      integer l_ik1, k_ik1                ! [output],
      integer nom2, nik2_0, nik2_, nik2   ! [output],
      integer l_om2, k_om2                ! [output],
      integer l_ik2_0, k_ik2_0            ! [output],
      integer l_ik2_, k_ik2_              ! [output],
      integer l_ik2, k_ik2                ! [output],

      
      nsh = (1 + lmax_)**2

c     ==== J ====
      nij0 = nsh * nqn * nch1
      if (.not. MA_alloc_get(MT_DCPL, nij0, 'intd_j0',
     &     l_ij0, k_ij0))
     &     call errquit('WFAT: intd_j0 allocation failed', 59, MA_ERR)

      nij_ = nsh * nqn * nch1 * nat
      if (.not. MA_alloc_get(MT_DCPL, nij_, 'intl_j_',
     &     l_ij_, k_ij_))
     &     call errquit('WFAT: intl_j_ allocation failed', 59, MA_ERR)

      nij = nsh * nqn * nch1 * nat
      if (.not. MA_alloc_get(MT_DCPL, nij, 'intl_j', l_ij, k_ij))
     &     call errquit('WFAT: intl_j allocation failed', 59, MA_ERR)


c     ==== Omega1 and K1 ====
      nom1 = nsh * nqn * nch * nptp
      if (.not. MA_alloc_get(MT_DCPL, nom1, 'omg1', l_om1, k_om1))
     &     call errquit('WFAT: omg1 allocation failed', 59, MA_ERR)
      
      nik1_0 = nsh * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nik1_0, 'intd_k1_0',
     &     l_ik1_0, k_ik1_0))
     &     call errquit('WFAT: intd_k1_0 allocation failed', 59, MA_ERR)

      nik1_ = nsh * nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nik1_, 'intl_k1_',
     &     l_ik1_, k_ik1_))
     &     call errquit('WFAT: intl_k1_ allocation failed', 59, MA_ERR)

      nik1 = nsh * nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nik1, 'intl_k1', l_ik1, k_ik1))
     &     call errquit('WFAT: intl_k1 allocation failed', 59, MA_ERR)

      
c     ==== Omega2 and K2 ====
      nom2 = nsh * nqn * nch * nptp
      if (.not. MA_alloc_get(MT_DCPL, nom2, 'omg2', l_om2, k_om2))
     &     call errquit('WFAT: omg2 allocation failed', 59, MA_ERR)
      
      nik2_0 = nsh * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nik2_0, 'intd_k2_0',
     &     l_ik2_0, k_ik2_0))
     &     call errquit('WFAT: intd_k2_0 allocation failed', 59, MA_ERR)

      nik2_ = nsh * nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nik2_, 'intl_k2_',
     &     l_ik2_, k_ik2_))
     &     call errquit('WFAT: intl_k2_ allocation failed', 59, MA_ERR)

      nik2 = nsh * nqn * nch * nat
      if (.not. MA_alloc_get(MT_DCPL, nik2, 'intl_k2', l_ik2, k_ik2))
     &     call errquit('WFAT: intl_k2 allocation failed', 59, MA_ERR)
      
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dealloc_intgl(l_om, l_i0, l_i_, l_il)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_om, l_i0, l_i_, l_il     ! [input]

      
      if (.not. MA_free_heap(l_om))
     &     call errquit('WFAT: omega deallocation failed', 61, MA_ERR)
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
      subroutine wfat_dealloc_intgl_wfat1(l_ij0, l_ij_, l_ij,
     &     l_om1, l_ik1_0, l_ik1_, l_ik1, l_om2, l_ik2_0, l_ik2_, l_ik2)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_ij0, l_ij_, l_ij     ! [input]
      integer l_om1, l_ik1_0, l_ik1_, l_ik1
      integer l_om2, l_ik2_0, l_ik2_, l_ik2

      if (.not. MA_free_heap(l_ij0))
     &     call errquit('WFAT: intd_j0 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_ij_))
     &     call errquit('WFAT: intl_j_ deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_ij))
     &     call errquit('WFAT: intl_j deallocation failed', 61, MA_ERR)

      if (.not. MA_free_heap(l_om1))
     &     call errquit('WFAT: omg1 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_ik1_0))
     &     call errquit('WFAT: intd_k1_0 deallocation failed', 61,
     &     MA_ERR)
      if (.not. MA_free_heap(l_ik1_))
     &     call errquit('WFAT: intl_k1_ deallocation failed', 61,
     &     MA_ERR)
      if (.not. MA_free_heap(l_ik1))
     &     call errquit('WFAT: intl_k1 deallocation failed', 61, MA_ERR)

      if (.not. MA_free_heap(l_om2))
     &     call errquit('WFAT: omg2 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_ik2_0))
     &     call errquit('WFAT: intd_k2_0 deallocation failed', 61,
     &     MA_ERR)
      if (.not. MA_free_heap(l_ik2_))
     &     call errquit('WFAT: intl_k2_ deallocation failed', 61,
     &     MA_ERR)
      if (.not. MA_free_heap(l_ik2))
     &     call errquit('WFAT: intl_k2 deallocation failed', 61, MA_ERR)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_check_wfat1ovl(nch, vch, vch1)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "wfat_qm_param.fh"
      integer nch
      double precision vch(nbas0, nch)
      double precision vch1(nbas0, 3, nch)

      integer i, j, l_ov, k_ov
      double precision val(3)


      write(luout, '("Overlaps of the 0-th order ionized orbital '//
     &     'and their first order corrections")')

      if (.not. ma_alloc_get(mt_dbl, nbas0**2, 'ovl', l_ov, k_ov))
     &     call errquit('WFAT: Cannot allocate ovl ', 74, MA_ERR)
      call ga_get(gwf_ovl, 1, nbas0, 1, nbas0, dbl_mb(k_ov), nbas0)


      do i = 1, nch
         do j = 1, 3
            call wfat_mat_element(nbas0, vch1(:,j,i), vch(:,i),
     &           dbl_mb(k_ov), val(j))
         enddo
         write(luout, '("Channel : ", i5)') i
         write(luout, '("  <v_x|v> = ", es21.12)') val(1)
         write(luout, '("  <v_y|v> = ", es21.12)') val(2)
         write(luout, '("  <v_z|v> = ", es21.12)') val(3)
         write(luout, '(" ")')
      enddo

      if (.not. ma_free_heap(l_ov))
     &     call errquit('WFAT: Cannot deallocate ovl', 74, MA_ERR )

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_intgl_chkpt(npt, chk)

      implicit none
      integer npt
      integer chk(10)

      integer i

      do i = 1, 10
         chk(i) = int( i/10.0d0 * npt )
      enddo
      chk(10) = npt

      end
c=======================================================
      
      
      
      
      

      
c!!=======================================================
c!!=======================================================
csubroutine wfat_check_big_intgl(icent, x_c, y_c, z_c, ii, chpos, &
c                                x_ch, y_ch, z_ch, iorb, &
c                                veffpsi, vspsi, vcpsi, max_omega, cellf, max_ingtl)
c
cc  use, intrinsic :: iso_c_binding
c  implicit none
c  integer icent, ii, iorb
c  double precision x_c, y_c, z_c, chpos(3), x_ch, y_ch, z_ch, &
c                                veffpsi, vspsi, vcpsi, max_omega, cellf, max_ingtl
c
c  
c  print *, " "
c  write(*,'("Checking for super big values:")')
c  write(*,'("  icent = ", i3)') icent
c  write(*,'("  Point and index = ", 3e15.6, i8)') x_c, y_c, z_c, ii
c  write(*,'("  Orbital center = ", 3e15.6)') chpos
c  write(*,'("  Distance from orb. center = ", e15.6)') sqrt(x_ch**2 + y_ch**2 + z_ch**2)
c  write(*,'("  orbital no. = ", i3)')  iorb
c  write(*,'("  veffpsi = ", e15.6)') veffpsi
c  write(*,'("  vspsi = ", e15.6)')   vspsi
c  write(*,'("  vcpsi = ", e15.6)')   vcpsi
c  write(*,'("  omega = ", e15.6)')  max_omega
c  write(*,'("  cell fun. = ", e15.6)') cellf
c  write(*,'("  max. abs. wf_intgl0 = ", e15.6)') max_ingtl
c  stop
c  
cend subroutine wfat_check_big_intgl
c!!=======================================================
