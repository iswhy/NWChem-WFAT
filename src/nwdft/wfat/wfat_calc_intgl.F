c=======================================================
c=======================================================
c==== Imam ====
      subroutine wfat_calc_intgl(rtdb, nch, orbid, spinid, nqn, nid, 
     &                           mid, lmax_, gqual, smh, res_mode, 
     &                           erg, chpos, chmmt, iondip, g_il)

      implicit none
#include "rtdb.fh"
#include "stdio.fh"      
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "errquit.fh"
#include "geom.fh"
#include "wfat.fh"
#include "wfat_qm_param.fh"
      integer rtdb                    ! [input]
      integer nch                     ! [input]
      integer orbid(nch)              ! [input]
      integer spinid(nch)             ! [input]
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      integer lmax_                   ! [input]
      integer gqual                   ! [input]
      double precision smh            ! [input]
      logical res_mode                ! [input]
      double precision erg(nch)       ! [output]
      double precision chpos(3, nch)  ! [output]
      double precision chmmt(3, nch)  ! [output]
      double precision iondip(3, nch) ! [output]
      integer g_il                    ! [output], global array handle for intgl. For structure factors.
      integer g_j                     ! [output]
      integer g_k(2)                  ! [output]
      
      character*16 tag
      character*9 rdm_name(2)
      logical t1, wfat1
      integer natom, nbas, ncf, nrd, max_npt, nvc, nsh, nom, ni0, ni_
      integer npt, natomp, nptp, ncpos, nscr, nil, nmo, nrda, nintj0
      integer nintj_, nintj
      integer i, j, ic, icent, ipt, ipt_i, ipt_f
      integer jpt, ich, iqn, ish, i_vc, ndid
      integer l_cf, l_rd, l_pt, l_wt, l_vc, l_om, l_i0
      integer l_i_, l_sc, l_il, l_icp, l_ipp
      integer l_mo, l_rda
      integer k_cf, k_rd, k_pt, k_wt, k_vc, k_om, k_i0
      integer k_i_, k_sc, k_il, k_icp, k_ipp
      integer k_mo, k_rda
      integer l_intj0, l_intj_, l_intj
      integer k_intj0, k_intj_, k_intj
      integer is, io, ipt0, i0l, i0r, oml, omr, i_l, i_r, imf, i_rd
      integer g_rd(2), i_mo
      double precision atp(3), znuc, v_c(3), v_m(3)
      double precision normch, celf, wtcelf, numnel
      double complex czero
      
      parameter (imf = 0)
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )



      call ga_sync()
      ndid = ga_nodeid()

      
c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom_hd, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
c     ==== Get number of bases ====
      if (.not. bas_numbf(bas_hd, nbas))
     &     call errquit('WFAT: Getting nbas failed in ', 53, UERR)

      
c     ==== Obtain properties of each channel ====
      if (.not. rtdb_get(rtdb, 'wfat:wfat1', mt_log, 1, wfat1))
     &     call errquit('WFAT: Cannot read wfat1 from RTDB', 27,
     &     RTDB_ERR)
c     nch = the number of the ionized channels.
c     nch1 = the number of channels that contains the correction terms (needed for WFAT(1)),
c            it is equal to 3*nch because each ionized channel has three correction terms 
c            associated to the three cartesian components of the dipole operator.
      if (wfat1) then
         nch1 = 3*nch
      else
         nch1 = 0
      endif
      ncf = nbas * (nch + nch1)
      if (.not. MA_alloc_get(MT_DBL, ncf, 'vch', l_cf, k_cf))
     &     call errquit('WFAT: vch allocation failed', 59, MA_ERR)
      call wfat_orb_property(rtdb, nch, orbid, spinid, dbl_mb(k_cf), 
     &                       erg, chpos, chmmt, iondip)
      if (wfat1)
     &     call wfat_perturbnew_cwfn(nch, spinid, orbid,
     &     dbl_mb(k_cf+nbas*nch))
      

c     ==== Get the reduced density matrix ====
      nrd = nbas * nbas * nspin0
      if (.not. MA_alloc_get(MT_DBL, nrd, 'rdm', l_rd, k_rd))
     &     call errquit('WFAT: rdm allocation failed', 59, MA_ERR)
      call wfat_get_rdm(dbl_mb(k_rd:k_rd+nrd-1))

      
c     ==== Some analyses and tests (optional) ====
      call wfat_sometests(zion, nch, erg, lmax_, chmmt)

      call wfat_mulliken(rtdb, geom_hd, bas_hd, nbas, nspin0,
     &                   dbl_mb(k_rd), nch, dbl_mb(k_cf), spinid, 
     &                   chpos)


c     ==== Allocate grid-related arrays ====
      call wfat_alloc_quadgrid(gqual, max_npt, l_pt, k_pt, l_wt, k_wt)


c     ==== Allocate WFAT integral, its integrands, and some other ====
      nscr = 20*nbas + 4*nspin0 + nbas * (3*nbas+1+nspin0) +
     &       4*(nch+nch1)
      nvc = nch + nch1
      if (.not. MA_alloc_get(MT_DBL, nscr, 'vscr', l_sc, k_sc))
     &     call errquit('WFAT: vscr allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nvc, 'vcpsi', l_vc, k_vc))
     &     call errquit('WFAT: vcpsi allocation failed', 59, MA_ERR)
      call wfat_alloc_intgl(lmax_, nqn, nch, natom, nsh, nom, ni0, ni_,
     &     nil, l_om, k_om, l_i0, k_i0, l_i_, k_i_, l_il, k_il)
      if (wfat1) call wfat_alloc_intgl_wfat1(lmax_, nqn, nch1, nat, nsh,
     &     nintj0, nintj_, nintj, l_intj0, k_intj0, l_intj_, k_intj_,
     &     l_intj, k_intj)              add K integrals here
      

c     ==== Obtaining grid ranges for the current process ====
      call wfat_alloc_gridrange(natom, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange(geom_hd, ndid, gqual, int_mb(k_icp),
     &                        int_mb(k_ipp), natomp, nptp)
c      write(luout, '("node:", i3, "  atom1,pt1:", i3, i7, '//
c     &     '"  atom2,pt2:", i3, i7, "  npoint:", i7)')
c     &     ndid, icent_i, ipt_i, icent_f, ipt_f, nptp

      
      if (ndid == 0) write(luout, 901)
 901  format(/,
     &     3x, '*** INFO ***', /,
     &     5x, 'Begin calculating WFAT integrals. If this is a ',
     &     '''start'' type job, this part can take up to ', /,
     &     5x, 'several minutes. Please wait.')
      
      
c     ==== BEGIN QUADRATURE SUMMATION ====
c     ==== Loop over atoms ====
      call ga_sync()
      jpt = 0
      dcpl_mb(k_i_ : k_i_+ni_-1) = czero
      if (wfat1) dcpl_mb(k_intj_ : k_intj_+nintj_-1) = czero
      do ic = 1, natomp
         icent = int_mb(k_icp + ic - 1)
         if (.not. geom_cent_get(geom_hd, icent, tag, atp, znuc))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)

         call wfat_get_gridpt(int(znuc), gqual, max_npt, npt,
     &                        dbl_mb(k_pt : k_pt+3*max_npt-1),
     &                        dbl_mb(k_wt : k_wt+max_npt-1))

         ipt_i = int_mb(k_ipp + (ic-1)*2)
         ipt_f = int_mb(k_ipp + 1 + (ic-1)*2)
c        ==== Loop over all quadrature points around the current atom ====
         do ipt = ipt_i, ipt_f
            ipt0 = k_pt + 3*(ipt-1)
            v_c(1:3) = dbl_mb(ipt0 : ipt0+2)
            call geom_cent2cent(geom_hd, icent, v_c, imf, v_m, t1)
cc            write(luout, '("        v_m = ", 3e14.5)') v_m

            
c           ==== Calculate vc|psi> ====
            call wfat_dft_vcpsi_1pt(geom_hd, bas_hd, nbas, nspin0, 
     &           nch+nch1, dbl_mb(k_cf:k_cf+ncf-1), spinid, chpos, v_m, 
     &           dbl_mb(k_rd:k_rd+nrd-1), zion, smh,
     &           dbl_mb(k_sc:k_sc+nscr-1),
     &           dbl_mb(k_vc:k_vc+nvc-1))

            
c           ==== Calculate Omega^(n,m)_(l,m') functions ====
            call wfat_omega_1pt(geom_hd, v_m, nch, erg,
     &           chpos, nqn, nid, mid, lmax_, int(zion),
     &           dcpl_mb(k_om:k_om+nom-1))
            
c           ==== Calculate the integrand of I^{(n,m)}_{l,m'} ====
            call wfat_integrand(nch, nqn, nsh, dcpl_mb(k_om),
     &           dbl_mb(k_vc), dcpl_mb(k_i0))

            if (wfat1) then
c              ==== Calculate the integrand of J^{(n,m)x_i}_{l,m'} ====
c              =====             K^{(n,m)}_{l,m'}                  ====
               call wfat_integrand_wfat1(nch, nqn, nsh, dcpl_mb(k_om),
     &              dbl_mb(k_vc+nch), dcpl_mb(k_intj0))       add k integrals here
            endif
            
            
c           ==== Quadrature summation ====
            call wfat_celf_1pt(geom_hd, icent, v_c, wf_msm, celf)
            wtcelf = celf * dbl_mb(k_wt+ipt-1)

            i_l = k_i_ + (icent-1) * nsh*nqn*nch
            i_r = i_l + nsh*nqn*nch - 1
            dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                         dcpl_mb(k_i0:k_i0+ni0-1) *
     &                         wtcelf                        !! 1)

            if (wfat1) then
               i_l = k_intj_ + (icent-1) * nsh*nqn*nch1
               i_r = i_l + nsh*nqn*nch1 - 1
               dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                            dcpl_mb(k_intj0:k_intj0+nintj0-1) *
     &                            wtcelf
            endif

            jpt = jpt + 1            
         enddo
      enddo
      

c     ==== Reorder elements of intgl ====      
      call wfat_reorder_intgl(nsh, natom, nqn, nch,
     &     dcpl_mb(k_i_ : k_i_+ni_-1), dcpl_mb(k_il : k_il+nil-1))
      if (wfat1)
     &     call wfat_reorder_intgl_wfat1(nsh, natom, nqn, nch,
     &     dcpl_mb(k_intj_ : k_intj_+nintj_-1),
     &     dcpl_mb(k_intj : k_intj+nintj-1))     add K integrals here
      
      

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
      if (ndid==0) write(luout, '("intgl norm = ", es14.6)') normch

      if (wfat1) then
         call ga_acc(g_j, 1, nintj, 1, 1, dcpl_mb(k_intj), nintj, 1.0d0)
         call ga_sync()

         call ga_norm1(g_j, normch)
         if (ndid==0) write(luout, '("intgl_j norm = ", es14.6)') normch

         do i = 1, 2
            call ga_acc(g_k(i))   resolve here
            call ga_sync()

            call ga_norm1(g_k(i), normch)
            if (ndid==0) write(luout, '("intgl_k norm = ", es14.6)')
     &           normch
         enddo
      endif

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: vch deallocation failed',
     &     61, MA_ERR)
      if (.not. MA_free_heap(l_rd))
     &     call errquit('WFAT: rdm deallocation failed',
     &     61, MA_ERR)
      if (.not. ma_free_heap(l_nu))
     &     call errquit('WFAT: Cannot deallocate vnucvin.', 1, MA_ERR)
      if (.not. ma_free_heap(l_hf))
     &     call errquit('WFAT: Cannot deallocate vhfvin.', 1, MA_ERR)
      if (.not. ma_free_heap(l_xc))
     &     call errquit('WFAT: Cannot deallocate vxcvin.', 1, MA_ERR)
      if (.not. ma_free_heap(l_vs))
     &     call errquit('WFAT: Cannot deallocate vsvin.', 1, MA_ERR)
      call wfat_dealloc_quadgrid(l_pt, l_wt)
      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: vscr deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_vc))
     &     call errquit('WFAT: vc|psi> deallocation failed', 61, MA_ERR)
      call wfat_dealloc_intgl(l_om, l_i0, l_i_, l_il)
      if (wfat1)
     &     call wfat_dealloc_intgl_wfat1(l_intj0, l_intj_, l_intj) add K integrals here
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      if (exasym) then
         if(.not. MA_free_heap(l_mo))
     &        call errquit('WFAT: MO deallocation failed', 61, MA_ERR)
         if(.not. MA_free_heap(l_rda))
     &        call errquit('WFAT: RDMas deallocation failed',
     &        61, MA_ERR)
      endif
      
      
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
      subroutine wfat_integrand_wfat1(nch, nqn, nsh, omg, vxi, int_j)

      implicit none
      integer nch
      integer nqn
      integer nsh
      double complex omg(nsh, nqn, nch)
      double precision vxi(3, nch)
      double complex int_j(nsh, nqn, 3, nch)

      integer i, j

      do i = 1, nch
         do j = 1, 3
            int_j(:,:,j,i) = conjg(omg(:,:,i)) * vxi(j,i)
         enddo
      enddo

      end
c=======================================================
      

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
     &           intj_, intj)

      implicit none
      integer nsh
      integer ncent
      integer nqn
      integer nch
      double complex intj_(nsh, nqn, 3, nch, ncent)
      double complex intj(nsh, ncent, nqn, 3, nch)

      integer ich, i, iqn, icent, i_lm
      
      
      do ich = 1, nch
         do i = 1, 3
            do iqn = 1, nqn
               do icent = 1, ncent
                  do i_lm = 1, nsh
                     intj(i_lm, icent, iqn, i, ich) =
     &                    intj_(i_lm, iqn, i, ich, icent)
                  enddo
               enddo
            enddo
         enddo
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_store_intgl(rtdb, nsh, ncent, nqn, nch,
     &           erg, chpos, chmmt, iondip, g_il)

      implicit none
#include "rtdb.fh"      
#include "stdio.fh"
#include "wfat.fh"
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb                    ! [input]
      integer nsh                     ! [input]
      integer ncent                   ! [input]
      integer nqn                     ! [input]
      integer nch                     ! [input]
      double precision erg(nch)       ! [input]
      double precision chpos(3,nch)   ! [input]
      double precision chmmt(3,nch)   ! [input]
      double precision iondip(3,nch)  ! [input]
      integer g_il                    ! [input]

      character*4096 fname
      logical rtdb_stat, isthere
      integer ndid, l_il_, k_il_
      integer i1, i2, j, ich, iqn, icent, i_lm
      integer nil_
      

      ndid = ga_nodeid()
      if (ndid /= 0)
     &     call errquit('WFAT: The subroutine ''wfat_store_intgl'' '//
     &     'cannot be called by process other than process 0.',
     &     1, GA_ERR)
      rtdb_stat = rtdb_parallel(.false.)


      if (.not. rtdb_put(rtdb, 'wfat:erg', mt_dbl, nch, erg))
     &     call errquit('WFAT: Failed writing erg to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:chpos', mt_dbl, 3*nch, chpos))
     &     call errquit('WFAT: Failed writing chpos to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:chmmt', mt_dbl, 3*nch, chmmt))
     &     call errquit('WFAT: Failed writing chmmt to RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'wfat:iondip', mt_dbl, 3*nch, iondip))
     &     call errquit('WFAT: Failed writing iondip to RTDB.',
     &     1, RTDB_ERR)
      

      nil_ = nsh * ncent * nqn
      if (.not. ma_alloc_get(mt_dcpl, nil_, 'intgl_', l_il_, k_il_))
     &     call errquit('WFAT: Cannot allocate intgl_.', 1, MA_ERR)
      

      if (.not. rtdb_cget(rtdb, 'wfat:intgl_name', 1, fname))
     &     call errquit('WFAT: Cannot read wfat:intgl_name from RTDB.',
     &     1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      if (isthere) then
         write(luout, '(/, 3x, "*** ATTENTION ***")')
         write(luout, '(5x, "A WFAT integral file named '''//
     &        trim(adjustl(fname))//''' has already existed. This '//
     &        'file will be overwritten.")')
      endif
      open(WF_INTGL_UNIT, file=fname, action='write',
     &     form='unformatted')

      
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Writing WFAT integrals to '//
     &     trim(adjustl(fname))//' ... ")')
      write(luout, '(7x, "No. of spherical harmonics = ", i4)') nsh
      write(luout, '(7x, "No. of atoms = ", i4)') ncent
      write(luout, '(7x, "No. of parabolic quantum numbers = ", i4)')
     &     nqn
      write(luout, '(7x, "No. of ionized channels = ", i4)') nch
      

      i1 = -nil_ + 1
      i2 = i1 + nil_ - 1
      do ich = 1, nch
         i1 = i1 + nil_
         i2 = i2 + nil_
         call ga_get(g_il, i1, i2, 1, 1, dcpl_mb(k_il_), nil_)

         j = 0
         do iqn = 1, nqn
            do icent = 1, ncent
               do i_lm = 1, nsh
                  j = j + 1
                  write(WF_INTGL_UNIT) dcpl_mb(k_il_+j-1)
               enddo
            enddo
         enddo
      enddo

      close(WF_INTGL_UNIT)

      if (.not. ma_free_heap(l_il_))
     &     call errquit('WFAT: Cannot deallocate intgl_.', 1, MA_ERR)
      rtdb_stat = rtdb_parallel(.true.)

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_load_intgl(rtdb, nsh, ncent, nqn, nch,
     &           erg, chpos, chmmt, iondip, g_il)

      implicit none
#include "rtdb.fh"      
#include "stdio.fh"
#include "wfat.fh"
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb                    ! [input]
      integer nsh                     ! [input]
      integer ncent                   ! [input]
      integer nqn                     ! [input]
      integer nch                     ! [input]
      double precision erg(nch)       ! [output]
      double precision chpos(3,nch)   ! [output]
      double precision chmmt(3,nch)   ! [output]
      double precision iondip(3,nch)  ! [output]
      integer g_il                    ! [output]

      character*4096 fname
      logical rtdb_stat, isthere
      integer ndid, l_il_, k_il_
      integer i1, i2, j, ich, iqn, icent, i_lm
      integer nil_

      
      ndid = ga_nodeid()
      if (ndid /= 0)
     &     call errquit('WFAT: The subroutine ''wfat_load_intgl'' '//
     &     'cannot be called by process other than process 0.',
     &     1, GA_ERR)
      rtdb_stat = rtdb_parallel(.false.)


      if (.not. rtdb_get(rtdb, 'wfat:erg', mt_dbl, nch, erg))
     &     call errquit('WFAT: Failed reading erg from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:chpos', mt_dbl, 3*nch, chpos))
     &     call errquit('WFAT: Failed reading chpos from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:chmmt', mt_dbl, 3*nch, chmmt))
     &     call errquit('WFAT: Failed reading chmmt from RTDB.',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:iondip', mt_dbl, 3*nch, iondip))
     &     call errquit('WFAT: Failed reading iondip from RTDB.',
     &     1, RTDB_ERR)
      

      nil_ = nsh * ncent * nqn
      if (.not. ma_alloc_get(mt_dcpl, nil_, 'intgl_', l_il_, k_il_))
     &     call errquit('WFAT: Cannot allocate intgl_.', 1, MA_ERR)
      

      if (.not. rtdb_cget(rtdb, 'wfat:intgl_name', 1, fname))
     &     call errquit('WFAT: Cannot read wfat:intgl_name from RTDB.',
     &     1, RTDB_ERR)
      inquire(file=fname, exist=isthere)
      if (.not. isthere) then
         call errquit('WFAT: Cannot find the WFAT integral file '''//
     &        trim(adjustl(fname))//'''. This file is needed when '//
     &        'using ''restart'' as the startup mode. If it is not '//
     &        'specified in the input file, make sure it is located '//
     &        'in the current directory and is named according to '//
     &        '''<file_prefix>.owigl'' where file_prefix is the file '//
     &        'prefix of NWChem files for the present calculation. '//
     &        'Otherwise, it can be given through the input '//
     &        'file using the ''intgl'' input parameter '//
     &        'following the format ''intgl "<string value>"'' '//
     &        '(do not skip the double quotes).',
     &        0, UERR)
      endif
      open(WF_INTGL_UNIT, file=fname, action='read',
     &     form='unformatted')


      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Reading WFAT integrals from '//
     &     trim(adjustl(fname))//' ")')
      write(luout, '(7x, "No. of spherical harmonics = ", i4)') nsh
      write(luout, '(7x, "No. of atoms = ", i4)') ncent
      write(luout, '(7x, "No. of parabolic quantum numbers = ", i4)')
     &     nqn
      write(luout, '(7x, "No. of ionized channels = ", i4)') nch
      

      i1 = -nil_ + 1
      i2 = i1 + nil_ - 1
      do ich = 1, nch
         j = 0
         do iqn = 1, nqn
            do icent = 1, ncent
               do i_lm = 1, nsh
                  j = j + 1
                  read(WF_INTGL_UNIT) dcpl_mb(k_il_+j-1)
               enddo
            enddo
         enddo

         i1 = i1 + nil_
         i2 = i2 + nil_
         call ga_put(g_il, i1, i2, 1, 1, dcpl_mb(k_il_), nil_)
      enddo

      close(WF_INTGL_UNIT)

      if (.not. ma_free_heap(l_il_))
     &     call errquit('WFAT: Cannot deallocate intgl_.', 1, MA_ERR)

      rtdb_stat = rtdb_parallel(.true.)

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
      subroutine wfat_alloc_intgl_wfat1(lmax_, nqn, nch1, nat, nsh,
     &     nintj0, nintj_, nintj, l_intj0, k_intj0, l_intj_, k_intj_,
     &     l_intj, k_intj)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer lmax_                         ! [input],
      integer nqn                           ! [input],
      integer nch1                          ! [input],
      integer nat                           ! [input],
      integer nsh, nintj0, nintj_, nintj    ! [output],
      integer l_intj0, k_intj0              ! [output],
      integer l_intj_, k_intj_              ! [output],
      integer l_intj, k_intj                ! [output],


      nsh = (1 + lmax_)**2

      nintj0 = nsh * nqn * nch1
      if (.not. MA_alloc_get(MT_DCPL, nintj0, 'int_j0',
     &     l_intj0, k_intj0))
     &     call errquit('WFAT: int_j0 allocation failed', 59, MA_ERR)

      nintj_ = nsh * nqn * nch1 * nat
      if (.not. MA_alloc_get(MT_DCPL, nintj_, 'int_j_',
     &     l_intj_, k_intj_))
     &     call errquit('WFAT: int_j_ allocation failed', 59, MA_ERR)

      nintj = nsh * nqn * nch1 * nat
      if (.not. MA_alloc_get(MT_DCPL, nintj, 'int_j', l_intj, k_intj))
     &     call errquit('WFAT: int_j allocation failed', 59, MA_ERR)

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
      subroutine wfat_dealloc_intgl_wfat1(l_intj0, l_intj_, l_intj)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_intj0, l_intj_, l_intj     ! [input]


      if (.not. MA_free_heap(l_intj0))
     &     call errquit('WFAT: int_j0 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_intj_))
     &     call errquit('WFAT: int_j_ deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_intj))
     &     call errquit('WFAT: intj deallocation failed', 61, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_sometests(zi, nch, erg, lmax_, chmmt)

      implicit none
#include "rtdb.fh"      
#include "global.fh"
#include "mafdecls.fh"
      double precision zi
      integer nch
      double precision erg(nch)
      integer lmax_
      double precision chmmt(3, nch)

      integer ndid
      integer l_chgtest(2), l_r1test(2), l_qtest(2)
      integer l_phi1test(2), l_tautest(2), l_otest(2)
      double precision d_chgutest(2)


      ndid = ga_nodeid()

      if (rtdb_get(rtdb, 'wfat:test_chg', mt_int, 2, l_chgtest)) then
         if (ndid==0) call wfat_test_chg(l_chgtest(1), l_chgtest(2),
     &        zi, erg(1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_chgu', mt_dbl, 2, d_chgutest)) then
         if (ndid==0) call wfat_test_chgu(d_chgutest(1), d_chgutest(2))
      endif
      if (rtdb_get(rtdb, 'wfat:test_r1', mt_int, 2, l_r1test)) then
         if (ndid==0) call wfat_test_r1(l_r1test(1), l_r1test(2),
     &        erg(1), zi, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_q', mt_int, 2, l_qtest)) then
         if (ndid==0) call wfat_test_q(l_qtest(1), l_qtest(2),
     &        erg(1), zi, lmax_, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_phi1', mt_int, 2, l_phi1test)) then
         if (ndid==0) call wfat_test_phi1xi(l_phi1test(1), 
     &        l_phi1test(2), erg(1), zi, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_tau', mt_int, 2, l_tautest)) then
         if (ndid==0) call wfat_test_tau(l_tautest(1), l_tautest(2),
     &        erg(1), zi, lmax_, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test2_tau', mt_int, 2, l_tautest)) then
         if (ndid==0) call wfat_test2_tau(l_tautest(1), l_tautest(2),
     &        erg(1), zi, lmax_, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_o', mt_int, 2, l_otest)) then
         if (ndid==0) call wfat_test_o(l_otest(1), l_otest(2),
     &        erg(1), zi)
      endif


      call ga_sync()
      
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
