c=================================================
c=================================================
      subroutine mewfat_calc_intgl(rtdb, nqn, nid, mid, lmax_, gqual,
     &     smh, rmode, gserg0, gserg1, ionerg, edipm0, edipm1,
     &     chpos, chmmt, dipm1, g_il)

      implicit none
#include "stdio.fh"
#include "wfat.fh"
#include "mewfat_qm_param.fh"      
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb                ! [input],
      integer nqn                 ! [input],
      integer nid(nqn)            ! [input],
      integer mid(nqn)            ! [input],
      integer lmax_               ! [input],
      integer gqual               ! [input],
      double precision smh        ! [input],
      logical rmode               ! [input],
      double precision gserg0     ! [output], neutral GS energy
      double precision gserg1     ! [output], cation GS energy
      double precision ionerg     ! [output],
      double precision edipm0(3)  ! [output], unshifted neutral e. dipole moment
      double precision edipm1(3)  ! [output], unshifted cation e. dipole moment
      double precision chpos(3)   ! [output],
      double precision chmmt(3)   ! [output],
      double precision dipm1(3)   ! [output],
      integer g_il                ! [output],

      integer ndid
      integer nocmo0, nocmo1, npcf, nrcf, nqcf, nscf, natomp, nptp
      integer cx, cy, max_npt, nsh, nom, ni0, ni_, nil
      integer l_mo0, l_mo1, l_pcf, l_rcf, l_qcf, l_scf, l_icp, l_ipp
      integer l_pt, l_wt, l_om, l_i0, l_i_, l_il
      integer k_mo0, k_mo1, k_pcf, k_rcf, k_qcf, k_scf, k_icp, k_ipp
      integer k_pt, k_wt, k_om, k_i0, k_i_, k_il
      double precision normch


      call ga_sync()
      ndid = ga_nodeid()
      

c     ==== Get channel properties ====
      call mewfat_channel_property(rtdb, geom_hd, bas_hd, nbas, gqual,
     &     wf_msm, nel0, nel1, g_mov0, g_mov1, charge0, charge1, 
     &     multip0, multip1, gserg0, gserg1, ionerg, edipm0, edipm1, 
     &     chpos, chmmt, dipm1)


c     ==== Calculate determinant coefficients and get MO's ====
      call mewfat_alloc_ocmo(nbas, nel0, nel1, nocmo0, nocmo1,
     &     l_mo0, k_mo0, l_mo1, k_mo1)
      call mewfat_alloc_detcoef(spinx, spiny, nel0, cx, cy, npcf, nrcf, 
     &     nqcf, nscf, l_pcf, k_pcf, l_rcf, k_rcf, l_qcf, k_qcf, 
     &     l_scf, k_scf)
      call mewfat_calc_detcoef(geom_hd, bas_hd, nbas, spinx, spiny,
     &     nel0, nel1, cx, cy, g_mov0, g_mov1, oskel, dbl_mb(k_mo0),
     &     dbl_mb(k_mo1), dbl_mb(k_pcf), dbl_mb(k_rcf), dbl_mb(k_qcf),
     &     dbl_mb(k_scf))
      if (ndid == 0)
     &     call mewfat_print_detcoef(spinx, nel0, cx, cy, dbl_mb(k_pcf),
     &     dbl_mb(k_rcf), dbl_mb(k_qcf), dbl_mb(k_scf))


c     ==== Obtaining grid ranges for the current process ====
      call wfat_alloc_gridrange(natm, l_icp, k_icp, l_ipp, k_ipp)
      call wfat_get_gridrange(geom_hd, ndid, gqual, int_mb(k_icp),
     &     int_mb(k_ipp), natomp, nptp)


c     ==== Allocate grid-related arrays ====
      call wfat_alloc_quadgrid(gqual, max_npt, l_pt, k_pt, l_wt, k_wt)
      
      
c     ==== Perform quadrature summation ====
      call wfat_alloc_intgl(lmax_, nqn, 1, natm, nsh, nom, ni0, ni_,
     &     nil, l_om, k_om, l_i0, k_i0, l_i_, k_i_, l_il, k_il)
      call mewfat_calc_intgl_qsum(rtdb, geom_hd, bas_hd, nbas, natomp, 
     &     int_mb(k_icp), int_mb(k_ipp), gqual, spinx, spiny, 
     &     nel0, nel1, dbl_mb(k_mo0), dbl_mb(k_mo1), 
     &     dbl_mb(k_pcf), dbl_mb(k_rcf), dbl_mb(k_qcf), dbl_mb(k_scf), 
     &     zion, chpos, ionerg, max_npt, dbl_mb(k_pt), dbl_mb(k_wt), 
     &     wf_msm, nsh, nqn, nid, mid, lmax_, natm, smh,
     &     dcpl_mb(k_om), dcpl_mb(k_i_), dcpl_mb(k_i0), dcpl_mb(k_il))


c     ==== Accumulate to global intgl (g_il) ====
      call ga_acc(g_il, 1, nil, 1, 1, dcpl_mb(k_il), nil, 1.0d0)
      call ga_sync()
      
      call ga_norm1(g_il, normch)
cfordebug      if (ndid==0) write(luout, '("intgl norm = ", es14.6)') normch


c     ==== Deallocations ====
      call mewfat_dealloc_ocmo(l_mo0, l_mo1)
      call mewfat_dealloc_detcoef(l_pcf, l_rcf, l_qcf, l_scf)
      call wfat_dealloc_gridrange(l_icp, l_ipp)
      call wfat_dealloc_quadgrid(l_pt, l_wt)
      call wfat_dealloc_intgl(l_om, l_i0, l_i_, l_il)

      end
c=================================================
      

c=================================================
c=================================================
      subroutine mewfat_calc_intgl_qsum(rtdb, geom, bas, nb, natp, icp,
     &     ipp, gq, sx, sy, n0, n1, ocmo0, ocmo1, pcf, rcf, qcf, scf, 
     &     zi, chpos, ionerg, max_npt, qpt, qwt, mcel, nsh, nqn, nid,
     &     mid, lmax_, nat, smh, omega, intg_, intg0, intg)

      implicit none
#include "util.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
      integer rtdb                                    ! [input], RTDB handle.
      integer geom                                    ! [input], geometry handle.
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of basis functions in bas.
      integer natp                                    ! [input], no. of centers for the current process returned by wfat_get_gridrange.
      integer icp(natp)                               ! [input], centers for the current process returned by wfat_get_gridrange.
      integer ipp(2,natp)                             ! [input], range of points of each centers for the current process returned by wfat_get_gridrange.
      integer gq                                      ! [input], grid quality.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                    ! [input], P det. coefficients returned by mewfat_calc_detcoef. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient returned by mewfat_calc_detcoef.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient returned by mewfat_calc_detcoef. The 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients returned by mewfat_calc_detcoef. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input], charge of the parent ion.
      double precision chpos(3)                       ! [input], center position of the channel returned by wfat_get_gridrange.
      double precision ionerg                         ! [input], ionization energy of the channel returned by wfat_get_gridrange.
      integer max_npt                                 ! [input], max. no. of points for the given gq as returned by wfat_get_gridmaxpt.
      double precision qpt(3, max_npt)                ! [output/scratch], array Quadrature points
      double precision qwt(max_npt)                   ! [output/scratch], array Quadrature weights
      integer mcel                                    ! [input], the m parameter related to the cell function. Usually taken as wf_msm from wfat.fh.
      integer nsh                                     ! [input], no. of spherical harmonics.
      integer nqn                                     ! [input], no. of parabolic quantum number pairs.
      integer nid(nqn)                                ! [input], n parabolic quantum numbers.
      integer mid(nqn)                                ! [input], m parabolic quantum numbers.
      integer lmax_                                   ! [input], max. angular momentum quantum number.
      integer nat                                     ! [input], no. of atoms in the molecule.
      double precision smh                            ! [input], smoothing parameter.
      double complex omega(nsh, nqn)                  ! [output/scratch],
      double complex intg_(nsh, nqn, nat)             ! [output/scratch],
      double complex intg0(nsh, nqn)                  ! [output/scratch],
      double complex intg(nsh, nqn, nat)              ! [output], 

      character*5 si_type
      character*12 dt_str
      character*16 tag
      logical t1, rs, isdft
      integer ndid, nsc, l_sc, k_sc, l_rd0, k_rd0, npt, ic, icent
      integer i, imf, nsc_rs, nptp, pr_chk(10), pr_count, kk
      integer g_nptx, g_npt_tot
      double precision nptx, npt_tot, t_now, t_last, t_start, dt_all
      double precision atp(3), zn, v_m(3), vcpsi, celf, dt_prog
      double complex czero

      parameter (imf = 0)
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )

      logical nwxc_has_cam, nwxc_has_xc
      external nwxc_has_cam, nwxc_has_xc


      ndid = ga_nodeid()
      call ga_sync()

      
c     ==== Calculate neutral RDM in case of DFT ====
      rs = nwxc_has_cam()
      isdft = nwxc_has_xc()
      if (isdft) then
         if (.not. ma_alloc_get(mt_dbl, nb*nb*2, 'rdm0', l_rd0, k_rd0))
     &        call errquit('MEWFAT: Cannot allocate rdm0', 1, MA_ERR)
         call mewfat_calc_rdm0(nb, n0, ocmo0, dbl_mb(k_rd0))

         if (.not. rtdb_cget(rtdb, 'mewfat:si', 1, si_type))
     &        call errquit('MEWFAT: Failed reading SI from RTDB.', 1,
     &        RTDB_ERR)
      endif

      
c     ==== Scratch array ====
      if (rs) then
         nsc_rs = nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0)
      else
         nsc_rs = 0
      endif
      nsc = 20*nb + n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) +
     &      sum(n1*n0) + nsc_rs                ! 1)
      if (.not. ma_alloc_get(mt_dbl, nsc, 'scratch', l_sc, k_sc))
     &     call errquit('MEWFAT: Cannot allocate scratch', 1, MA_ERR)
            

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
      nptp = 0
      do ic = 1, natp
         nptp = nptp + ipp(2,ic) - ipp(1,ic) + 1
      enddo
      call wfat_get_intgl_chkpt(nptp, pr_chk)
      call wfat_init_progress(nptp, g_nptx, npt_tot)
      if (ndid == 0) write(luout, '(5x, "The total number of '//
     &     'quadrature grid points is ", i9, ".")') int(npt_tot)

      t_last = util_cpusec()
      t_start = t_last
      pr_count = 0
      kk = 1
      intg_ = czero              ! intg_(nsh,nqn,nat)     nch=1
      do ic = 1, natp
         icent = icp(ic)        ! icp(natp)
         if (.not. geom_cent_get(geom, icent, tag, atp, zn))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
         call wfat_get_gridpt(int(zn), gq, max_npt, npt, qpt, qwt)   ! qpt(3,max_npt)   qwt(max_npt)
         
c        ==== Loop over all quadrature points around the current atom ====
         do i = ipp(1,ic), ipp(2,ic)           ! ipp(2,natp)
            call geom_cent2cent(geom, icent, qpt(:,i), imf, v_m, t1)

c           ==== Calculate vc|psi> ====
            if (isdft) then
               call mewfat_dft_vcpsi_1pt(geom, bas, nb, sx, sy, n0, n1,
     &              v_m, ocmo0, ocmo1, dbl_mb(k_rd0), pcf, rcf, qcf, 
     &              scf, zi, chpos, smh, rs, si_type, dbl_mb(k_sc),
     &              vcpsi)
            else
               call mewfat_hf_vcpsi_1pt(geom, bas, nb, sx, sy, n0, n1,
     &              v_m, ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, chpos, 
     &              smh, rs, dbl_mb(k_sc), vcpsi)
            endif
            
c           ==== Calculate Omega^(n,m)_(l,m') functions ====
            call wfat_omega_1pt(geom, v_m, 1, ionerg, chpos, nqn,
     &           nid, mid, lmax_, int(zi), omega)      ! omega(nsh,nqn)   nch=1
            
c           ==== Calculate the integrand of I^(n,m)_(l,m') ====
            intg0 = conjg(omega) * vcpsi         ! intg0(nsh,nqn)   nch=1
cfordebug            if (ndid==0) then
cfordebug               write(luout, '("s1", i3, i8, 2es14.5, es14.5, es14.5)')
cfordebug     &              ic, i, sum(omega), vcpsi, vcpsi/rcf
cfordebug               write(luout, '(" ")')
cfordebug            endif

            
c           ==== Quadrature summation ====
            call wfat_celf_1pt(geom, icent, qpt(:,i), mcel, celf)
            intg_(:,:,icent) = intg_(:,:,icent) + intg0 * celf * qwt(i)

            
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
      call wfat_reorder_intgl(nsh, nat, nqn, 1, intg_, intg)    ! intg(nsh,nqn,nat)    nch=1
      call ga_sync()


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_sc))
     &     call errquit('MEWFAT: Cannot deallocate scratch', 1, MA_ERR)
      if (isdft) then
         if (.not. ma_free_heap(l_rd0))
     &        call errquit('MEWFAT: Cannot deallocate rdm0', 8, MA_ERR)
      endif
      call wfat_terminate_progress(g_nptx)
      
      end
c     NOTES:
c     1) This is the min. value for the length of scratch. This value is based on the 
c        info next to scr argument declaration in mewfat_potential.F.
c=================================================


c=================================================
c=================================================
      subroutine mewfat_alloc_ocmo(nb, n0, n1, nocmo0, nocmo1,
     &     l_mo0, k_mo0, l_mo1, k_mo1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nb              ! [input],
      integer n0(2)           ! [input],
      integer n1(2)           ! [input],
      integer nocmo0          ! [output],
      integer nocmo1          ! [output],
      integer l_mo0, k_mo0    ! [output],
      integer l_mo1, k_mo1    ! [output],

      
      nocmo0 = nb * sum(n0)
      if (.not. ma_alloc_get(mt_dbl, nocmo0, 'mo0', l_mo0, k_mo0))
     &     call errquit('MEWFAT: Cannot allocate mo0', 1, MA_ERR)
      nocmo1 = nb * sum(n1)
      if (.not. ma_alloc_get(mt_dbl, nocmo1, 'mo1', l_mo1, k_mo1))
     &     call errquit('MEWFAT: Cannot allocate mo1', 1, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_dealloc_ocmo(l_mo0, l_mo1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_mo0, l_mo1       ! [input]
      
      if (.not. MA_free_heap(l_mo0))
     &     call errquit('WFAT: mo0 deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_mo1))
     &     call errquit('WFAT: mo1 deallocation failed', 61, MA_ERR)

      end
c=================================================
      

c=================================================
c=================================================
      subroutine mewfat_alloc_detcoef(sx, sy, n0, cx, cy, npcf, nrcf, 
     &     nqcf, nscf, l_pcf, k_pcf, l_rcf, k_rcf, l_qcf, k_qcf,
     &     l_scf, k_scf)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer sx              ! [input],
      integer sy              ! [input],
      integer n0(2)           ! [input], no. of electrons in the neutral.
      integer cx              ! [output],
      integer cy              ! [output],
      integer npcf            ! [output], 
      integer nrcf            ! [output], 
      integer nqcf            ! [output], 
      integer nscf            ! [output], 
      integer l_pcf, k_pcf    ! [output], 
      integer l_rcf, k_rcf    ! [output], 
      integer l_qcf, k_qcf    ! [output], 
      integer l_scf, k_scf    ! [output], 


      npcf = n0(sx)
      if (.not. ma_alloc_get(mt_dbl, npcf, 'pcf', l_pcf, k_pcf))
     &     call errquit('MEWFAT: Cannot allocate pcf', 1, MA_ERR)
      nrcf = 1
      if (.not. ma_alloc_get(mt_dbl, nrcf, 'rcf', l_rcf, k_rcf))
     &     call errquit('MEWFAT: Cannot allocate pcf', 1, MA_ERR)
      cx = max(1, n0(sx)-1)
      nqcf = n0(sx) * n0(sx) * cx
      if (.not. ma_alloc_get(mt_dbl, nqcf, 'qcf', l_qcf, k_qcf))
     &     call errquit('MEWFAT: Cannot allocate qcf', 1, MA_ERR)
      cy = max(1, n0(sy))
      nscf = cy * cy
      if (.not. ma_alloc_get(mt_dbl, nscf, 'scf', l_scf, k_scf))
     &     call errquit('MEWFAT: Cannot allocate pcf', 1, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_dealloc_detcoef(l_pcf, l_rcf, l_qcf, l_scf)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_pcf, l_rcf, l_qcf, l_scf


      if (.not. MA_free_heap(l_pcf))
     &     call errquit('WFAT: pcf deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_rcf))
     &     call errquit('WFAT: rcf deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_qcf))
     &     call errquit('WFAT: qcf deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_scf))
     &     call errquit('WFAT: scf deallocation failed', 61, MA_ERR)

      end
c=================================================
