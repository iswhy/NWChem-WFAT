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
      if (ke <= 0.0d0) call errquit('WFAT: The input argument ke of '//
     &     'wfat_lapl_r_laguerre must be such that ke > 0.0', 58, UERR)
      
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
      subroutine wfat_eval_fbrr(m, nr, ke, r, rfun)
c     DESCRIPTION:
c     Evaluates the radial FBR basis based on the associated Laguerre functions.
      implicit none
#include "errquit.fh"
      integer m
      integer nr
      double precision ke
      double precision r
      double precision rfun(nr)

      integer i
      double precision mm, laguerre_norm_i
      external laguerre_norm_i

      if (r < 0.0d0)
     &     call errquit('WFAT: The input argument r of '//
     &     'wfat_eval_fbrr must be such that r >= 0.', 80, UERR)
      
      mm = dble(m)
      call wfat_lagfun(1, nr-1, abs(mm), ke*r, rfun)
c     At this point, rfun(1:nr) contains the generalized Laguerre polynomials
c     L_n^{|m|} from n=0 to n=nr-1.

      do i = 1, nr
         rfun(i) = 1/sqrt(laguerre_norm_i(i-1, abs(m))) *
     &        (ke*r)**(abs(mm)/2) * exp(-ke*r/2) * rfun(i)
      enddo

      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_eval_fbrt(m, nt, theta, tfun)
c     DESCRIPTION:
c     Evaluates the theta FBR basis based on the associated Legendre functions.
      implicit none
#include "errquit.fh"
      integer m                     ![input], the magnetic quantum number.
      integer nt                    ![input], the number of theta FBR basis states.
      double precision theta        ![input], polar angle in degrees to evaluate tfun.
      double precision tfun(1:nt)   ![output], the theta FBR basis states.

      integer i, l
      double precision x, th, assocleg_norm, deg2rad
      parameter(deg2rad=0.01745329251994d0)
      external assocleg_norm


      if (theta < 0.0d0 .or. theta > 180.0d0)
     &     call errquit('WFAT: The input argument theta of '//
     &     'wfat_eval_fbrt must be such that 0 <= theta <= 180.',
     &     80, UERR)

      th = theta * deg2rad
      x = cos(th)
      call wfat_assocleg_fixm(abs(m), abs(m)+nt-1, x, tfun)

      do i = 1, nt
         l = abs(m) + i - 1
         tfun(i) = tfun(i) / sqrt(assocleg_norm(l, abs(m)))
      enddo

      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_dvr_fbr_r(m, nr, ke, u)
c     DESCRIPTION:
c     Calculates the transformation matrix from radial DVR basis to radial
c     FBR basis.
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr
      double precision ke
      double precision u(nr,nr)

      integer l_pt, k_pt, l_wt, k_wt, l_fb, k_fb, j, i
      double precision mm, w, r, fbrr
      
      
      if (.not. MA_alloc_get(MT_DBL, nr, 'pt', l_pt, k_pt))
     &     call errquit('WFAT: pt allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nr, 'wt', l_wt, k_wt))
     &     call errquit('WFAT: wt allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nr, 'fbr', l_fb, k_fb))
     &     call errquit('WFAT: fbr allocation failed', 59, MA_ERR)


c     ==== Get the quadrature nodes and weights ====
      call wfat_lf_quadrature_rule(nr, abs(mm), dbl_mb(k_pt),
     &     dbl_mb(k_wt))

      
c     ==== Compute the transformation matrix elements ====
      mm = dble(m)
      do j = 1, nr      ! Loop over DVR bases.
         w = dbl_mb(k_wt+j-1)
         r = dbl_mb(k_pt+j-1) / ke
         call wfat_eval_fbrr(m, nr, ke, r, dbl_mb(k_fb))
         do i = 1, nr   ! Loop over FBR bases.
            fbrr = dbl_mb(k_fb+i-1)
            u(i,j) = sqrt(w) * (ke*r)**(-abs(mm)/2) * exp(ke*r/2) * fbrr
         enddo
      enddo


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate pt', 32, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate wt', 32, MA_ERR)
      if (.not. ma_free_heap(l_fb))
     &     call errquit('WFAT: Cannot deallocate fb', 32, MA_ERR)
      
      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_dvr_fbr_t(m, nt, u)
c     DESCRIPTION:
c     Calculates the transformation matrix from theta DVR basis to theta
c     FBR basis.
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nt
      double precision u(nt,nt)

      integer j, i, l_pt, k_pt, l_wt, k_wt, l_fb, k_fb
      double precision mm, w, x, fbrt, deg2rad
      parameter (deg2rad=0.01745329251994d0)


      if (.not. MA_alloc_get(MT_DBL, nt, 'pt', l_pt, k_pt))
     &     call errquit('WFAT: pt allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nt, 'wt', l_wt, k_wt))
     &     call errquit('WFAT: wt allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nt, 'fbr', l_fb, k_fb))
     &     call errquit('WFAT: fbr allocation failed', 59, MA_ERR)

      
c     ==== Get the quadrature nodes and weights ====
      call wfat_assocleg_quad(-1.0d0, 1.0d0, abs(m), nt, dbl_mb(k_pt),
     &     dbl_mb(k_wt))


c     ==== Compute the transformation matrix elements ====
      mm = dble(m)
      do j = 1, nt              ! Loop over FBR bases.
         w = dbl_mb(k_wt+j-1)
         x = dbl_mb(k_pt+j-1)
         call wfat_eval_fbrt(m, nt, acos(x)/deg2rad, dbl_mb(k_fb))
         do i = 1, nt           ! Loop over DVR bases.
            fbrt = dbl_mb(k_fb+i-1)
            u(i,j) = sqrt(w) * (1-x**2)**(-abs(mm)/2) * fbrt
         enddo
      enddo

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate pt', 32, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate wt', 32, MA_ERR)
      if (.not. ma_free_heap(l_fb))
     &     call errquit('WFAT: Cannot deallocate fb', 32, MA_ERR)
      
      end
c==========================================================
