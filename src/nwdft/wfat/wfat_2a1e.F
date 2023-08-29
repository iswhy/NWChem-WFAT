c==========================================================
c==========================================================
      subroutine wfat_lapl_r_laguerre(n, m, ke, lapr)

      implicit none
#include "errquit.fh"
      integer n
      integer m
      double precision ke
      double precision lapr(n,n)

      integer i, j, delkron
      double precision mm, v, normi, normj, laguerre_norm_i
      external delkron, laguerre_norm_i


      if (n < 1) call errquit('WFAT: The input argument n of '//
     &     'wfat_lapl_r_laguerre must be such that n >= 1', 58, UERR)
      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_lapl_r_laguerre must be such that m >= 0', 58, UERR)
      if (k <= 0.0d0) call errquit('WFAT: The input argument k of '//
     &     'wfat_lapl_r_laguerre must be such that k > 0.0', 58, UERR)
      
      mm = dble(m)
      
      do i = 1, n
         normi = sqrt(laguerre_norm_i(i-1,m))
         do j = 1, i
            normj = sqrt(laguerre_norm_i(j-1,m))
            call wfat_laguerre_intg1(mm, mm, mm+2, j-1, i-1, v)
            lapr(j,i) = 1/ke * (0.25d0*v +
     &           i*(i+mm/2) * delkron(i,j-1) -
     &           (i+mm/2)**2 * delkron(i,j) -
     &           (i-1+mm)*(i-1) * delkron(i,j) +
     &           (i-1+mm)*(i-1+mm/2) * delkron(i-1,j)) / normi / normj
            lapr(i,j) = lapr(j,i)
         enddo
      enddo

      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_eval_fbrt(m, nt, theta, tfun)

      implicit none
#include "errquit.fh"
      integer m
      integer nt                    ![input], the number of theta FBR basis states.
      double precision theta        ![input], polar angle in degrees to evaluate tfun.
      double precision tfun(1:nt)   ![output], the theta FBR basis states.

      integer i, l
      double precision x, th, assocleg_norm, pi_const
      parameter(pi_const=3.14159265359d0)
      external assocleg_norm


      if (theta < 0.0d0 .or. theta > 180.0d0)
     &     call errquit('WFAT: The input argument theta of '//
     &     'wfat_eval_fbrt must be such that 0 <= theta <= 180.',
     &     80, UERR)

      th = theta / (180.0d0/pi_const)
      x = cos(th)
      call wfat_assocleg_fixm(abs(m), abs(m)+nt-1, x, tfun)

      do i = 1, nt
         l = abs(m) + i - 1
         tfun(i) = tfun(i) / sqrt(assocleg_norm(l, abs(m)))
      enddo

      end
c==========================================================
