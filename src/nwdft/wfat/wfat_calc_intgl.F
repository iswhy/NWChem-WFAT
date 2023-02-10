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
      integer g_il                    ! [output], global array handle for intgl.

      character*16 tag
      character*9 rdm_name(2)
      logical t1, exasym
      integer natom, nbas, ncf, nrd, max_npt, nvc, nsh, nom, ni0, ni_
      integer npt, natomp, nptp, ncpos, nscr, nil, nmo, nrda
      integer nnu, nhf, nxc, nvs
      integer i, j, ic, icent, ipt, ipt_i, ipt_f
      integer jpt, ich, iqn, ish, i_vc, ndid
      integer l_cf, l_cf1, l_rd, l_pt, l_wt, l_vc, l_om, l_i0
      integer l_i_, l_sc, l_il, l_icp, l_ipp, l_nu, l_hf, l_xc, l_vs
      integer l_mo, l_rda
      integer k_cf, k_cf1, k_rd, k_pt, k_wt, k_vc, k_om, k_i0
      integer k_i_, k_sc, k_il, k_icp, k_ipp, k_nu, k_hf, k_xc, k_vs
      integer k_mo, k_rda
      integer is, io, ipt0, i0l, i0r, oml, omr, i_l, i_r, imf, i_rd
      integer g_rd(2), i_mo, l_chgtest(2), l_r1test(2), l_qtest(2)
      integer l_tautest(2)
      double precision atp(3), znuc, v_c(3), v_m(3), asc(3), rasfac
      double precision r_mas, ras, exc(2), normch, celf, numnel
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

      
      !!==== Obtain properties of each channel ====!!
      ncf = nbas * nch
      if (.not. MA_alloc_get(MT_DBL, ncf, 'vch', l_cf, k_cf))
     &     call errquit('WFAT: vch allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, ncf, 'v.vch', l_cf1, k_cf1))
     &     call errquit('WFAT: veff.vch allocation failed', 59, MA_ERR)
      call wfat_orb_property(rtdb, nch, orbid, spinid, dbl_mb(k_cf), 
     &                       erg, chpos, chmmt, iondip)


c     ==== Get the reduced density matrix ====
      nrd = nbas * nbas * nspin0
      if (.not. MA_alloc_get(MT_DBL, nrd, 'rdm', l_rd, k_rd))
     &     call errquit('WFAT: rdm allocation failed', 59, MA_ERR)
      call wfat_get_rdm(dbl_mb(k_rd:k_rd+nrd-1))

      
c     ==== Some tests ====
      if (rtdb_get(rtdb, 'wfat:test_chg', mt_int, 2, l_chgtest)) then
         if (ndid==0) call wfat_test_chg(l_chgtest(1), l_chgtest(2),
     &        zion, erg(1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_r1', mt_int, 2, l_r1test)) then
         if (ndid==0) call wfat_test_r1(l_r1test(1), l_r1test(2),
     &        erg(1), zion, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_q', mt_int, 2, l_qtest)) then
         if (ndid==0) call wfat_test_q(l_qtest(1), l_qtest(2),
     &        erg(1), zion, lmax_, chmmt(:,1))
      endif
      if (rtdb_get(rtdb, 'wfat:test_tau', mt_int, 2, l_tautest)) then
         if (ndid==0) call wfat_test_tau(l_tautest(1), l_tautest(2),
     &        erg(1), zion, lmax_, chmmt(:,1))
      endif
      call ga_sync()
      



      

      
c     ==== Calculate v_c|orb> ====
      nnu = nbas * nch
      if (.not. ma_alloc_get(mt_dbl, nnu, 'vnucvin', l_nu, k_nu))
     &     call errquit('WFAT: Cannot allocate vnucvin.', 1, MA_ERR)
      nhf = nbas * nch
      if (.not. ma_alloc_get(mt_dbl, nhf, 'vhfvin', l_hf, k_hf))
     &     call errquit('WFAT: Cannot allocate vhfvin.', 1, MA_ERR)
      nxc = nbas * nch
      if (.not. ma_alloc_get(mt_dbl, nxc, 'vxcvin', l_xc, k_xc))
     &     call errquit('WFAT: Cannot allocate vxcvin.', 1, MA_ERR)
      nvs = nbas * nch
      if (.not. ma_alloc_get(mt_dbl, nvs, 'vsvin', l_vs, k_vs))
     &     call errquit('WFAT: Cannot allocate vsvin.', 1, MA_ERR)
c      call wfat_dft_vcpsi_getcoef(dbl_mb(k_nu),
c     &                   dbl_mb(k_hf),
c     &                   dbl_mb(k_xc),
c     &                   dbl_mb(k_vs))




      call wfat_mulliken(rtdb, geom_hd, bas_hd, nbas, nspin0,
     &                   dbl_mb(k_rd), nch, dbl_mb(k_cf), spinid, 
     &                   chpos)
      
      
            





      
c     ==== Allocate grid-related arrays ====
      call wfat_alloc_quadgrid(gqual, max_npt, l_pt, k_pt, l_wt, k_wt)


c     ==== Allocate WFAT integral, its integrands, and some other ====
      if (.not. rtdb_get(rtdb, 'wfat:exasym', mt_log, 1, exasym))
     &     exasym = .false.
      if (exasym) then
         nscr = 20*nbas + 4*nspin0 + nbas * (3*nbas+1+nspin0) + 4*nch +
     &          nbas + nbas*nch
      else
         nscr = 20*nbas + 4*nspin0 + nbas * (3*nbas+1+nspin0) + 4*nch
      endif
      if (.not. MA_alloc_get(MT_DBL, nscr, 'vscr', l_sc, k_sc))
     &     call errquit('WFAT: vscr allocation failed', 59, MA_ERR)
      nvc = nch
      if (.not. MA_alloc_get(MT_DBL, nvc, 'vcpsi', l_vc, k_vc))
     &     call errquit('WFAT: vcpsi allocation failed', 59, MA_ERR)

      call wfat_alloc_intgl(lmax_, nqn, nch, natom, nsh, nom, ni0, ni_,
     &     nil, l_om, k_om, l_i0, k_i0, l_i_, k_i_, l_il, k_il)
      
      if (.not. rtdb_get(rtdb, 'wfat:rasfac', mt_dbl, 1, rasfac))
     &     rasfac = 9.0d0
      call wfat_calc_asysphere(geom_hd, rasfac, asc, ras)
      if (exasym) then
         nmo = nbas*nbas*nspin0
         if (.not. MA_alloc_get(mt_dbl, nmo, 'MO', l_mo, k_mo))
     &        call errquit('WFAT: Cannot allocate MO', 1, MA_ERR)
         do is = 1, nspin0
            i_mo = k_mo + (is-1)*nbas**2
            call ga_get(g_mov(is), 1, nbas, 1, nbas, dbl_mb(i_mo),
     &           nbas)
         enddo
         nrda = nbas*nbas*nspin0
         if (.not. MA_alloc_get(mt_dbl, nrda, 'RDMas', l_rda, k_rda))
     &        call errquit('WFAT: Cannot allocate RDMas', 1, MA_ERR)
         if (ndid == 0) then
            write(luout, '(" Coordinate of asymp. center : ", 3es14.5)')
     &           asc
            write(luout, '(" Radius of asymp. sphere : ", es14.5)') ras
         endif
      endif
      

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
            r_mas = sqrt(sum((v_m - asc)**2))
            if (exasym .and. r_mas > ras) then
               call wfat_calc_rdmas(nbas, nspin0, 
     &              dbl_mb(k_mo:k_mo+nmo-1), v_m, asc, ras,
     &              dbl_mb(k_rda:k_rda+nrda-1))
               call wfat_dft_vcasy_1pt(geom_hd, bas_hd, nbas, nspin0,
     &              nch, dbl_mb(k_cf:k_cf+ncf-1), erg, spinid, chpos,
     &              v_m, asc, ras, dbl_mb(k_rd:k_rd+nrd-1),
     &              dbl_mb(k_rda:k_rda+nrda-1), zion, smh,
     &              dbl_mb(k_sc:k_sc+nscr-1), dbl_mb(k_vc:k_vc+nvc-1))
            else
               call wfat_dft_vcpsi_1pt(geom_hd, bas_hd, nbas, nspin0, 
     &              nch, dbl_mb(k_cf:k_cf+ncf-1), spinid, chpos, v_m, 
     &              dbl_mb(k_rd:k_rd+nrd-1), zion, smh,
     &              dbl_mb(k_sc:k_sc+nscr-1),
     &              dbl_mb(k_vc:k_vc+nvc-1))
            endif
            

c            call wfat_dft_vcpsi_fromcoef(geom_hd, bas_hd, nbas, nch,
c     &           dbl_mb(k_nu), dbl_mb(k_hf), dbl_mb(k_xc), dbl_mb(k_vs),
c     &           v_m, dbl_mb(k_sc), dbl_mb(k_vc))
            
            
c           ==== Calculate Omega^(n,m)_(l,m') functions ====
            call wfat_omega_1pt(geom_hd, v_m, nch, erg,
     &           chpos, nqn, nid, mid, lmax_, int(zion),
     &           dcpl_mb(k_om:k_om+nom-1))
            
c           ==== Calculate the integrand of I^(n,m)_(l,m') ====
            do i = 1, nch
               i0l = k_i0 + (i-1)*nsh*nqn
               i0r = i0l + nsh*nqn - 1
               oml = k_om + (i-1)*nsh*nqn
               omr = oml + nsh*nqn - 1
               i_vc = k_vc + i - 1
csee               if (abs(dbl_mb(i_vc)) > 1.0d-12) then
                  dcpl_mb(i0l:i0r) = conjg(dcpl_mb(oml:omr)) *
     &                               dbl_mb(i_vc)
csee               else
csee                  dcpl_mb(i0l:i0r) = 0.0d0
csee               endif
                              
            enddo

c           ==== Quadrature summation ====
            call wfat_celf_1pt(geom_hd, icent, v_c, wf_msm, celf)
            i_l = k_i_ + (icent-1) * nsh*nqn*nch
            i_r = i_l + nsh*nqn*nch - 1
csee            if (abs(celf) > 1.0d-12) then               
            dcpl_mb(i_l:i_r) = dcpl_mb(i_l:i_r) +
     &                            dcpl_mb(k_i0:k_i0+ni0-1) *
     &                            celf *
     &                            dbl_mb(k_wt+ipt-1)     !! 1)
csee            endif



cc            if (ndid==0)
cc     &           write(luout, '("(id, icent,ipt) = ", 3i7, '//
cc     &           '" | v_m = ", 3es14.6, '//
cc     &           '" | (vcpsi,omg) = ", 2es14.6, '//
cc     &           '" | sum intgl = ", 2es14.6)')
cc     &           ndid, icent, ipt, v_m, dbl_mb(k_vc),
cc     &           maxval(abs(dcpl_mb(k_om:k_om+nsh*nqn-1))),
cc     &           sum(dcpl_mb(i_l:i_r))
            
            
            
            jpt = jpt + 1
cc            write(luout, '("        sum vcpsi = ", e14.5)')
cc     &           sum(dbl_mb(k_vc:k_vc+nvc-1))
cc            write(luout, '("        sum intgl0 = ", 2e14.5)')
cc     &           sum(dcpl_mb(k_i0:k_i0+ni0-1))
cc            write(luout, '("        sum intgl_ = ", 2e14.5)')
cc     &           sum(dcpl_mb(i_l:i_r))
cc            write(luout, '(i3, i7, i3, i7, " sum omega = ", 2e14.5)')
cc     &           ndid, jpt, icent, ipt, sum(dcpl_mb(k_om:k_om+nom-1))
cc            if (abs(imag(sum(dcpl_mb(k_i0:k_i0+ni0-1)))) > 1.0d-8) then
ccc               call errquit('WFAT: Big imaginary part.', 2, UERR)
cc            endif
            
         enddo
      enddo
      

c     ==== Reorder elements of intgl ====      
      call wfat_reorder_intgl(nsh, natom, nqn, nch,
     &     dcpl_mb(k_i_ : k_i_+ni_-1), dcpl_mb(k_il : k_il+nil-1))

      

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
!!    Note that it may not be possible to perform this simultaneous summation if the
!!    channels may or may not be assigned its own quadrature points. In the current
!!    implementation, all channels don't have quadrature points, so that the quadrature
!!    points (and weights) solely come from the atoms.
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
