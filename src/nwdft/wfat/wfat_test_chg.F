c=======================================================
c=======================================================
      subroutine wfat_test_chg(lmin_, lmax_, zq, erg)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      integer lmin_, lmax_
      double precision zq, erg

      integer l_m, k_m, l_u, k_u
      integer i, l, nl, nr
      parameter (nr=5)
      double precision kappa, lmb, a, aa, b, bb, r(nr)


      write(luout, '("TESTING FOR M(a,b,z) and U(a,b,z)")')

      nl = lmax_ - lmin_ + 1
      if (.not. MA_alloc_get(MT_DBL, nl, 'chgf', l_m, k_m))
     &     call errquit('WFAT: chgf allocation failed', 30, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nl, 'chgu', l_u, k_u))
     &     call errquit('WFAT: chgu allocation failed', 31, MA_ERR)

      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0
      write(luout, '(1x)')
      write(luout, '(" erg = ", es14.6)') erg
      write(luout, '(" kappa = ", es14.6)') kappa
      r = (/ 0.5d0, 1.0d0, 4.8d0, 8.3d0, 24.7d0 /)

      do i = 1, nr
         write(luout, '(1x)')
         write(luout, '("------------------------------------")')
         write(luout, '("z = ", es14.6)') lmb * r(i)
         
         call wfat_chgf_1pt(r(i), lmin_, lmax_, a, b, lmb, dbl_mb(k_m))
         call wfat_chgu_1pt(r(i), lmin_, lmax_, a, b, lmb, dbl_mb(k_u))
         
         write(luout, '(1x)')
         write(luout, '(a14, a4, 2a14)') "a", "b", "M(a,b,z)",
     &        "U(a,b,z)"
         do l = lmin_, lmax_
            aa = l + a
            bb = 2*l + b
            write(luout, '(es14.6, i4, 2es14.6)') aa, int(bb),
     &           dbl_mb(k_m+l-lmin_), dbl_mb(k_u+l-lmin_)
         enddo

         write(luout, '("------------------------------------")')
      enddo


      if (.not. MA_free_heap(l_m))
     &     call errquit('WFAT: Cannot deallocate chgf', 27, MA_ERR)
      if (.not. MA_free_heap(l_u))
     &     call errquit('WFAT: Cannot deallocate chgu', 28, MA_ERR)


      end
c=======================================================
      
