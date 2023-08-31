c==========================================================
c==========================================================
      subroutine wfat_2a1e_solve(nr, nt, m, q1, q2, z1, z2, maxiter)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer nr
      integer nt
      integer m
      double precision q1
      double precision q2
      double precision z1
      double precision z2
      integer maxiter
      
      integer l_lr, l_ur, l_ir, l_gr, l_lt, l_ut, l_it, l_gt, l_hm, l_eg
      integer l_x, l_t, l_rw, l_tw
      integer k_lr, k_ur, k_ir, k_gr, k_lt, k_ut, k_it, k_gt, k_hm, k_eg
      integer k_x, k_t, k_rw, k_tw
      integer j, iter
      double precision ke, e0


c     ==== Radial part ====
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'lapr', l_lr, k_lr))
     &     call errquit('WFAT: Cannot allocate lapr', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'u_r', l_ur, k_ur))
     &     call errquit('WFAT: Cannot allocate u_r', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'i_r', l_ir, k_ir))
     &     call errquit('WFAT: Cannot allocate i_r', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'g_r', l_gr, k_gr))
     &     call errquit('WFAT: Cannot allocate g_r', 1, MA_ERR)
      
      call wfat_laplr_fbrr(nr, abs(m), dbl_mb(k_lr))
      call wfat_dvr_fbr_r(m, nr, dbl_mb(k_ur))
c     transpose(U_r) . T_r . U_r
      call wfat_represent_sym(nr, nr, dbl_mb(k_ur), dbl_mb(k_lr),
     &     dbl_mb(k_gr))
      call wfat_get_identity(nr, dbl_mb(k_ir))

c     ==== Get the radial quadrature nodes ====
      if (.not. ma_alloc_get(MT_DBL, nr, 'xquad', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad', 1, MA_ERR)
      call wfat_lf_quadrature_rule(nr, abs(dble(m)), dbl_mb(k_x),
     &     dbl_mb(k_rw))
      
      
c     ==== Theta part ====
      if (.not. ma_alloc_get(MT_DBL, nt, 'lapt', l_lt, k_lt))
     &     call errquit('WFAT: Cannot allocate lapt', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'u_t', l_ut, k_ut))
     &     call errquit('WFAT: Cannot allocate u_t', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'i_t', l_it, k_it))
     &     call errquit('WFAT: Cannot allocate i_t', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'g_t', l_gt, k_gt))
     &     call errquit('WFAT: Cannot allocate g_t', 1, MA_ERR)
      
      do j = 1, nt
         dbl_mb(k_lt+j-1) = dble((abs(m)+j-1) * (abs(m)+j))
      enddo
      call wfat_dvr_fbr_t(m, nt, dbl_mb(k_ut))
c     transpose(U_t) . T_t . U_t
      call wfat_represent_diag(nt, nt, dbl_mb(k_ut), dbl_mb(k_lt),
     &     dbl_mb(k_gt))
      call wfat_get_identity(nt, dbl_mb(k_it))

c     ==== Get the polar quadrature nodes ====
      if (.not. ma_alloc_get(MT_DBL, nt, 'tquad', l_t, k_t))
     &     call errquit('WFAT: Cannot allocate tquad', 1, MA_ERR)
      call wfat_assocleg_quad(-1.0d0, 1.0d0, abs(m), nt, dbl_mb(k_t),
     &     dbl_mb(k_tw))


c     ==== Hamiltonian construction ====
      if (.not. ma_alloc_get(MT_DBL, (nr*nt)**2, 'hm', l_hm, k_hm))
     &     call errquit('WFAT: Cannot allocate hm', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nt, 'erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate erg', 1, MA_ERR)

      e0 = 0.0d0
      ke = sqrt(2 * 1.106d0)
      do iter = 1, maxiter
         call wfat_2a1e_hamiltonian(nr, nt, dbl_mb(k_x), dbl_mb(k_t),
     &        dbl_mb(k_gr), dbl_mb(k_it), dbl_mb(k_ir), dbl_mb(k_gt), 
     &        ke, q1, q2, z1, z2, dbl_mb(k_hm))
      
c        ==== Hamiltonian diagonalization ====
         call wfat_dsyev_wrap(nr*nt, dbl_mb(k_hm), dbl_mb(k_eg))
         write(luout, '("energy, delta energy = ", i5, 2f15.8)')
     &        iter, dbl_mb(k_eg), dbl_mb(k_eg) - e0
         e0 = dbl_mb(k_eg)
         ke = sqrt( 2 * abs(dbl_mb(k_eg)) )
      enddo

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_lr))
     &     call errquit('WFAT: Cannot deallocate lapr', 1, MA_ERR)
      if (.not. ma_free_heap(l_ur))
     &     call errquit('WFAT: Cannot deallocate u_r', 1, MA_ERR)
      if (.not. ma_free_heap(l_ir))
     &     call errquit('WFAT: Cannot deallocate i_r', 1, MA_ERR)
      if (.not. ma_free_heap(l_gr))
     &     call errquit('WFAT: Cannot deallocate g_r', 1, MA_ERR)
      if (.not. ma_free_heap(l_x))
     &     call errquit('WFAT: Cannot deallocate xquad', 1, MA_ERR)

      if (.not. ma_free_heap(l_lt))
     &     call errquit('WFAT: Cannot deallocate lapt', 1, MA_ERR)
      if (.not. ma_free_heap(l_ut))
     &     call errquit('WFAT: Cannot deallocate u_t', 1, MA_ERR)
      if (.not. ma_free_heap(l_it))
     &     call errquit('WFAT: Cannot deallocate i_t', 1, MA_ERR)
      if (.not. ma_free_heap(l_gt))
     &     call errquit('WFAT: Cannot deallocate g_t', 1, MA_ERR)
      if (.not. ma_free_heap(l_t))
     &     call errquit('WFAT: Cannot deallocate tquad', 1, MA_ERR)

      if (.not. ma_free_heap(l_hm))
     &     call errquit('WFAT: Cannot deallocate hm', 1, MA_ERR)
      if (.not. ma_free_heap(l_eg))
     &     call errquit('WFAT: Cannot deallocate erg', 1, MA_ERR)
      
      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_2a1e_hamiltonian(nr, nt, xquad, tquad, gr, it, ir,
     &     gt, ke, q1, q2, z1, z2, hm)

      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"
      integer nr
      integer nt
      double precision xquad(nr)     ! [input], scaled radial quadrature nodes.
      double precision tquad(nt)     ! [input], theta quadrature nodes in degrees.
      double precision gr(nr, nr)
      double precision it(nt, nt)
      double precision ir(nr, nr)
      double precision gt(nt, nt)
      double precision ke
      double precision q1
      double precision q2
      double precision z1
      double precision z2
      double precision hm(nr*nt, nr*nt)

      integer i, j, n, i_tm1, i_tm2, l_tm1, k_tm1, l_tm2, k_tm2, i_r
      integer i_t
      double precision r, t, vdist1(2), vdist2(2), dist1, dist2, pot
      double precision deg2rad
      parameter(deg2rad=0.01745329251994d0)


      if (.not. MA_alloc_get(MT_DBL, (nr*nt)**2, 'tmp1', l_tm1, k_tm1))
     &     call errquit('WFAT: tmp1 allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, (nr*nt)**2, 'tmp2', l_tm2, k_tm2))
     &     call errquit('WFAT: tmp2 allocation failed', 59, MA_ERR)
      
      call wfat_kron_prod(nr, nr, nt, nt, gr, it, dbl_mb(k_tm1))
      call wfat_kron_prod(nr, nr, nt, nt, ir, gt, dbl_mb(k_tm2))

c     ==== Kinetic part ====
      n = nr * nt
      do j = 1, n
         do i = 1, n
            i_tm1 = k_tm1 + (j-1)*n + (i-1)
            i_tm2 = k_tm2 + (j-1)*n + (i-1)
            hm(i,j) = -0.5d0 * (dbl_mb(i_tm1) + dbl_mb(i_tm2)) / ke
         enddo
      enddo

c     ==== Nuclear attraction part ====
      i = 0
      do i_r = 1, nr
         r = xquad(i_r) / ke
         do i_t = 1, nt
            t = tquad(i_t) * deg2rad
            i = i + 1

            vdist1 = (/ r*sin(t), r*cos(t) - z1 /)
            vdist2 = (/ r*sin(t), r*cos(t) - z2 /)
            dist1 = sqrt(sum(vdist1 * vdist1))
            dist2 = sqrt(sum(vdist2 * vdist2))
            pot = - (q1/dist1 + q2/dist2)
            hm(i,i) = hm(i,i) + pot
         enddo
      enddo
      
                  
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_tm1))
     &     call errquit('WFAT: Cannot deallocate tmp1', 32, MA_ERR)
      if (.not. ma_free_heap(l_tm2))
     &     call errquit('WFAT: Cannot deallocate tmp2', 32, MA_ERR)
      
      end
c==========================================================

      
c==========================================================
c==========================================================
      subroutine wfat_laplr_fbrr(n, m, lapr)
c     DESCRIPTION:
C     Calculates the matrix element of the radial derivative part of the Laplacian
c     operator in the radial FBR basis based the generalized Laugerre polynomials.
c        lapr_{ij} = \int dx x^2 R_i^{|m|}(x) T(R_j^{|m|}(x))
c     where
c        T(f(x)) = 1/x \frac{d^2}{dx^2} (xf(x)) ,
c        R_i^{|m|}(x) = 1/\sqrt{h_{i-1,m}} x^{|m|/2} \exp(-x/2) L_{i-1}^{|m|}(x) , and
c        h_{i,m} = \Gamma(|m|+1) \binom{|m|+i}{i} .
      implicit none
#include "errquit.fh"
      integer n
      integer m
      double precision lapr(n,n)

      integer i, j, delkron
      double precision mm, v, normi, normj, laguerre_norm_i
      external delkron, laguerre_norm_i


      if (n < 1) call errquit('WFAT: The input argument n of '//
     &     'wfat_laplr_fbrr must be such that n >= 1', 58, UERR)
      if (m < 0) call errquit('WFAT: The input argument m of '//
     &     'wfat_laplr_fbrr must be such that m >= 0', 58, UERR)

      mm = dble(m)
      do i = 1, n
         normi = sqrt(laguerre_norm_i(i-1,m))
         do j = 1, i
            normj = sqrt(laguerre_norm_i(j-1,m))
            call wfat_laguerre_intg1(mm, mm, mm+2, j-1, i-1, v)
            lapr(j,i) = (0.25d0*v +
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
      subroutine wfat_eval_fbrr(m, nr, x, rfun)
c     DESCRIPTION:
c     Evaluates the radial FBR basis based on the associated Laguerre functions.
      implicit none
#include "errquit.fh"
      integer m
      integer nr
      double precision x
      double precision rfun(nr)

      integer i
      double precision mm, laguerre_norm_i
      external laguerre_norm_i

      if (x < 0.0d0)
     &     call errquit('WFAT: The input argument r of '//
     &     'wfat_eval_fbrr must be such that r >= 0.', 80, UERR)
      
      mm = dble(m)
      call wfat_lagfun(1, nr-1, abs(mm), x, rfun)
c     At this point, rfun(1:nr) contains the generalized Laguerre polynomials
c     L_n^{|m|} from n=0 to n=nr-1.

      do i = 1, nr
         rfun(i) = 1/sqrt(laguerre_norm_i(i-1, abs(m))) *
     &        x**(abs(mm)/2) * exp(-x/2) * rfun(i)
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
      subroutine wfat_dvr_fbr_r(m, nr, u)
c     DESCRIPTION:
c     Calculates the transformation matrix from radial DVR basis to radial
c     FBR basis.
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr
      double precision u(nr,nr)

      integer l_pt, k_pt, l_wt, k_wt, l_fb, k_fb, j, i
      double precision mm, w, x, fbrr
      

      mm = dble(m)
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
      do j = 1, nr      ! Loop over DVR bases.
         w = dbl_mb(k_wt+j-1)
         x = dbl_mb(k_pt+j-1)
         call wfat_eval_fbrr(m, nr, x, dbl_mb(k_fb))
         do i = 1, nr   ! Loop over FBR bases.
            fbrr = dbl_mb(k_fb+i-1)
            u(i,j) = sqrt(w) * x**(-abs(mm)/2) * exp(x/2) * fbrr
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
