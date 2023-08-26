c=======================================================
c=======================================================
      subroutine wfat_assocleg_quad(a, b, m, n, pt, wt)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      double precision a         ![input], left interval boundary.
      double precision b         ![input], right interval boundary.
      integer m                  ![input], the m parameter in P_{l,m}(x).
      integer n                  ![input], the no. of quadrature nodes.
      double precision pt(n)     ![output], quadrature nodes.
      double precision wt(n)     ![output], quadrature weights.

      integer i, j, j_, l_dm, k_dm, l_nd, k_nd, i_dm, i_dm_, lmax
      double precision wt0, factorial
      external factorial
      
      
      if (m < 0)
     &     call errquit('WFAT: The input m of wfat_assocleg_quad '//
     &     'cannot be negative.', 42, UERR)

      lmax = m + n      ! We will solve for the zeros of P_{lmax,m}(x).
      

c     ==== Construct the dipole moment matrix ====
      if (.not. MA_alloc_get(MT_DBL, n**2, 'DM', l_dm, k_dm))
     &     call errquit('WFAT: DM allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, n, 'node', l_nd, k_nd))
     &     call errquit('WFAT: node allocation failed', 59, MA_ERR)
      dbl_mb(k_dm : k_dm+n**2-1) = 0.0d0
      do j = 1, n
         do i = 1, n
            if (j == i + 1) then
               j_ = m + i
               i_dm = k_dm + (j-1)*n + (i-1)
               dbl_mb(i_dm) = sqrt((j_**2 - m**2)/(4*j_**2 - 1.0d0))   ! See eq. (3.24) to (3.26) in Ref. [1].
               i_dm_ = k_dm + (i-1)*n + (j-1)
               dbl_mb(i_dm_) = dbl_mb(i_dm)
            endif
         enddo
      enddo
      call wfat_dsyev_wrap(n, dbl_mb(k_dm), dbl_mb(k_nd))

      
c     ==== Calculate the shifted nodes and weights ====
      do i = 1, n
         i_dm = k_dm + (i-1)*n          ! The first row of the i-th column.
         wt0 = 2**(2*m+1) * (factorial(m))**2/factorial(2*m+1) *
     &        (dbl_mb(i_dm))**2                                     ! eq. (3.45) in Ref. [1].
         pt(i) = (b-a)/2 * (dbl_mb(k_nd+i-1) + (a+b)/(b-a))         ! The shifted node.
         wt(i) = (b-a)/2 * wt0                                      ! The shifted weight.
      enddo


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_dm))
     &     call errquit('WFAT: Cannot deallocate DM', 32, MA_ERR)
      if (.not. ma_free_heap(l_nd))
     &     call errquit('WFAT: Cannot deallocate node', 32, MA_ERR)
      
      end
c     REFERENCES:
c     [1] Spectral Methods - Algorithm, Analysis and Applications by Shen et al.
c=======================================================
