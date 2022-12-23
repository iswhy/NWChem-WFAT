c=======================================================
c=======================================================
      subroutine wfat_test_r1(n, m, erg, zq, chmmt)

      implicit none
#include "mafdecls.fh"
#include "stdio.fh"      
#include "errquit.fh"
#include "wfat.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      double precision chmmt(3)

      character*255 fname
      logical first_notsmalleta
      integer l_pt, k_pt, l_wt, k_wt, l_r1, k_r1
      integer pmax, nlgl, npt, i, id
      double precision kappa, a, b, et_last, de, et, del, z, chgm, chgu
      double precision mu(3), g, d, f, r0_1, o0_1, r0_2, o0_2, r1_
      double precision protmat(3,3)
      double complex zz, chgm_
      double precision factorial
      external factorial


c     ==== Get the passive rotation matrix ===
      call wfat_passrot(-0.0d0, -0.0d0, protmat)
      mu = matmul(transpose(protmat), -chmmt)
      
      
      pmax = 6
      nlgl = 10
      npt = 401
      kappa = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      

      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_pt', l_pt, k_pt))
     &     call errquit('WFAT: lgl_pt allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_wt', l_wt, k_wt))
     &     call errquit('WFAT: lgl_wt allocation failed', 85, MA_ERR)
      call wfat_lgl(-1.0d0, 1.0d0, nlgl-1, dbl_mb(k_pt), dbl_mb(k_wt))
      write(luout, '("Legendre-Gauss-Lobatto nodes and weights:")')
      do i = 1, nlgl
         write(luout, '(i5, 2es15.6)')
     &        i, dbl_mb(k_pt+i-1), dbl_mb(k_wt+i-1)
      enddo

      if (.not. MA_alloc_get(MT_DBL, npt, 'r1', l_r1, k_r1))
     &     call errquit('WFAT: r1 allocation failed', 86, MA_ERR)
      

      et_last = 1.0d0
      de = et_last / (npt-1)
      first_notsmalleta = .true.
      
      call util_file_name('test.wr1', .false., .false., fname)
      open(WF_R1_UNIT, file=fname, action='write')
      do i = 1, npt
         et = 0.0d0 + (i-1)*de
         if (et < 1.0d-2+de) then
            call wfat_r1_smalleta(et, n, m, erg, zq, mu(3), pmax,
     &           del, dbl_mb(k_r1+i-1))
         else
            if (first_notsmalleta) then
               z = kappa * (et-de)
               zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
               call wfat_cchg(a, b, zz, chgm_)
               chgm = real(chgm_)
               call wfat_chgu(a, b, z, chgu, id)

               call wfat_gammaf(a, g)
               d = g / factorial(abs(m))
               f = z**(b/2) * exp(-z/2)
               r0_1 = d * f * chgm
               o0_1 = f * chgu
               first_notsmalleta = .false.
            endif
            
            call wfat_r1(et-de, et, r0_1, o0_1, n, m, erg, zq, mu(3),
     &           nlgl, dbl_mb(k_pt), dbl_mb(k_wt), r0_2, o0_2, r1_)
            dbl_mb(k_r1+i-1) = dbl_mb(k_r1+i-2) + r1_
            r0_1 = r0_2
            o0_1 = o0_2
            del = 0.0d0
         endif

         write(WF_R1_UNIT, '(i5, 3es15.6)') i, et, dbl_mb(k_r1+i-1), del
      enddo
      close(WF_R1_UNIT)


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate lgl_pt', 49, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate lgl_wt', 49, MA_ERR)
      if (.not. ma_free_heap(l_r1))
     &     call errquit('WFAT: Cannot deallocate r1', 49, MA_ERR)

      end
c=======================================================
