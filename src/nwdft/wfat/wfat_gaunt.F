!==================================================
!==================================================
      subroutine wfat_gaunt_lscan(l1, m1, l2, m2, gn)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer l1
      integer m1
      integer l2
      integer m2
      double precision gn(1:*)

      integer l_min, n, n0, l_3j, k_3j, l_3j0, k_3j0, i, l3, i0
      double precision f

      
c     / l1  l2      l3   \
c     \ m1  m2  -(m1+m2) /
c     where l3=[l_min,l1+l2].
      l_min = max(abs(l1-l2), abs(m1+m2))
      n = (l1+l2) - l_min + 1
      if (.not. MA_alloc_get(MT_DBL, n, '3j', l_3j, k_3j))
     &     call errquit('WFAT: 3j allocation failed', 59, MA_ERR)
      call wfat_3j_jscan(2*l1, 2*m1, 2*l2, 2*m2, dbl_mb(k_3j))

      
c     / l1  l2  l3 \
c     \  0   0   0 /
c     where l3 = [|l1-l2|,l1+l2].
      n0 = 2*min(l1,l2) + 1
      if (.not. MA_alloc_get(MT_DBL, n0, '3j_000', l_3j0, k_3j0))
     &     call errquit('WFAT: 3j_000 allocation failed', 59, MA_ERR)
      call wfat_3j_jscan000(l1, l2, dbl_mb(k_3j0))
cOLD      do i = 1, n0
cOLD         write(luout, '("000 ", i5, f12.8)')
cOLD     &        abs(l1-l2)+i-1, dbl_mb(k_3j0+i-1)
cOLD      enddo

      
      gn(1:n) = 0.0d0
      do i = 1, n
         l3 = l_min + i - 1
         i0 = i + l_min - abs(l1-l2)
         if (mod(l1+l2+l3,2) == 0) then
            call wfat_gaunt_f(l1, m1, l2, m2, l3, f)
            gn(i) = dbl_mb(k_3j0+i0-1) * dbl_mb(k_3j+i-1)
            gn(i) = 2 * (-1)**(m1+m2) * sqrt(f) * gn(i)
         endif
cc         write(luout, '(i5, es16.8, 2f12.8)')
cc     &        l3, f, dbl_mb(k_3j0+i0-1), dbl_mb(k_3j+i-1)
      enddo


      if (.not. MA_free_heap(l_3j))
     &     call errquit('WFAT: 3j deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_3j0))
     &     call errquit('WFAT: 3j000 deallocation failed', 61, MA_ERR)
      end
!==================================================


!==================================================
!==================================================
      subroutine wfat_gaunt_f(l1, m1, l2, m2, l3, f)

      implicit none
#include "stdio.fh"
      integer l1
      integer m1
      integer l2
      integer m2
      integer l3
      double precision f

      integer i, j, nv
      double precision v(6)


      v(1) = l1 + m1
      v(2) = l2 + m2
      v(3) = l3 + m1 + m2
      v(4) = l1 - m1
      v(5) = l2 - m2
      v(6) = l3 - m1 - m2
      do i = 1, 6
         if (abs(v(i)) < 1.0d-12) v(i) = 1.0d0
      enddo
      
      nv = idnint(maxval(v))
      f = 1.0d0
      do i = 1, nv
c         write(luout, '("v = ", 6f12.8, 5x, es16.8)') v, f
         f = f * v(1) * v(2) * v(3) / v(4) / v(5) / v(6)
         do j = 1, 6
            v(j) = v(j) - 1
            if (idnint(v(j)) < 1) v(j) = 1.0d0
         enddo
      enddo

      end
!==================================================
      
      
