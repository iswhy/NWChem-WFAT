c=======================================================
c=======================================================
      subroutine wfat_q_particular(nr, r, n, m, lmax, erg, zq, nlgl,
     &                             q_p1, q_p2)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer nr
      double precision r(nr)    ! [input]
      integer n
      integer m
      integer lmax
      double precision erg
      double precision zq
      integer nlgl
      double precision q_p1(abs(m):lmax, nr)       ! [output], orientation-independent
      double precision q_p2(abs(m):lmax, nr)       ! [output], orientation-dependent

      integer i, l, l_pt, l_wt, k_pt, k_wt
      double precision kappa, lmb, a, b, f, g, h
      double precision g_i0(2,abs(m):lmax), h_i0(2,abs(m):lmax)
      double precision g_i(2,abs(m):lmax), h_i(2,abs(m):lmax)
      double precision gg(2,abs(m):lmax), hh(2,abs(m):lmax)
      double precision chgm(abs(m):lmax), chgu(abs(m):lmax)

      
      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0

      
      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_pt', l_pt, k_pt))
     &     call errquit('WFAT: lgl_pt allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_wt', l_wt, k_wt))
     &     call errquit('WFAT: lgl_wt allocation failed', 85, MA_ERR)
      call wfat_lgl(-1.0d0, 1.0d0, nlgl-1, dbl_mb(k_pt), dbl_mb(k_wt))
      
      
      do i = 1, nr             ! Loop over the increasingly-ordered radial points
         if (r(i) < 1.0d-13) then
            do l = abs(m), lmax
               q_p1(l,i) = 0.0d0
               q_p2(l,i) = 0.0d0
               g_i0(1:2,l) = 0.0d0
               h_i0(1:2,l) = 0.0d0
            enddo
            cycle
         endif

         call wfat_chgm_1pt(r(i), abs(m), lmax, a, b, lmb, chgm)
         call wfat_chgu_1pt(r(i), abs(m), lmax, a, b, lmb, chgu)

         call wfat_calc_gh_i(r(i-1), r(i), n, m, lmax, erg, zq, nlgl,
     &        dbl_mb(k_pt), dbl_mb(k_wt), gg, hh)
         g_i = g_i0 + gg
         h_i = h_i0 + hh
         
         do l = abs(m), lmax
            f = (kappa*r(i))**l * exp(-kappa*r(i))
            g = f * chgm(l)
            h = f * chgu(l)
            q_p1(l, i) = h * g_i(1,l) - g * h_i(1,l)
            q_p2(l, i) = h * g_i(2,l) - g * h_i(2,l)
         enddo

         g_i0 = g_i
         h_i0 = h_i
      enddo


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate lgl_pt', 49, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate lgl_wt', 49, MA_ERR)
      
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_gh_i(r1, r2, n, m, lmax, erg, zq, nlgl,
     &                          pt, wt, gg, hh)

      implicit none
#include "stdio.fh"      
      double precision r1
      double precision r2
      integer n
      integer m
      integer lmax
      double precision erg
      double precision zq
      integer nlgl
      double precision pt(nlgl)
      double precision wt(nlgl)
      double precision gg(2,abs(m):lmax)
      double precision hh(2,abs(m):lmax)

      integer nfem, ngr, id, i, j, l
      double precision pt_4(4), wt_4(4), pt_6(6), wt_6(6)
      double precision dr, dfem, kappa, lmb, a, b, r1_, r2_, r, wt_, f
      double precision g, h, c, c_m, c_p, k_m, k_p, y, ind_oi, ind_od
      double precision chgm(abs(m):lmax+1), chgu(abs(m):lmax)
      double precision min_dfem
      parameter( min_dfem=2.0d0 )


      pt_4 = (/ -1.0d0,
     &          -0.44721359549996d0,
     &           0.44721359549996d0,
     &           1.0d0 /)
      wt_4 = (/ 0.16666666666667d0,
     &          0.83333333333333d0,
     &          0.83333333333333d0,
     &          0.16666666666667d0 /)
      pt_6 = (/ -1.0d0,
     &          -0.76505532392947d0,
     &          -0.28523151648065d0,
     &           0.28523151648065d0,
     &           0.76505532392947d0,
     &           1.0d0 /)
      wt_6 = (/ 0.06666666666667d0,
     &          0.37847495629785d0,
     &          0.55485837703549d0,
     &          0.55485837703549d0,
     &          0.37847495629785d0,
     &          0.06666666666667d0 /)

      
c     ==== Determine the no. of finite elements and grids ====
c     dr = the total interval of r to be integrated
c     dfem = the finite element interval
      dr = r2 - r1
      if (dr < 1.0d-12) then
         nfem = 1
         dfem = dr
         ngr = 2
         id = 0
      else if (dr >= 1.0d-12 .and. dr < 1.0d-3) then
         nfem = 1
         dfem = dr
         ngr = 4
         id = 1
      else if (dr >= 1.0d-3 .and. dr < 1.0d-1) then
         nfem = 1
         dfem = dr
         ngr = 6
         id = 2
      else
         nfem = int(max(dr,min_dfem+1.0d-13) / min_dfem)
         dfem = dr / nfem
         ngr = nlgl
         id = 3
      endif

      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0

      gg = 0.0d0
      hh = 0.0d0
c     ==== Loop over finite elements ====
c     r1_ = lower bound of a FEM
c     r2_ = upper bound of a FEM
c     r = shifted quadrature grid
c     wt_ = shifted quadrature weight
      do i = 1, nfem
         r1_ = r1 + (i-1)*dfem
         r2_ = r1_ + dfem

c        ==== Loop over quadrature grids ====
         do j = 1, ngr
c           ==== Determine the shifted grids and weights ====
            if (id == 0) then
               r = r1_ + (j-1)*dfem
               wt_ = (r2_ - r1_)/2
            else if (id == 1) then
               r = 0.5d0*((r2_-r1_)*pt_4(j) + r2_ + r1_)
               wt_ = (r2_ - r1_)/2 * wt_4(j)
            else if (id == 2) then
               r = 0.5d0*((r2_-r1_)*pt_6(j) + r2_ + r1_)
               wt_ = (r2_ - r1_)/2 * wt_6(j)
            else if (id == 3) then
               r = 0.5d0*((r2_-r1_)*pt(j) + r2_ + r1_)
               wt_ = (r2_ - r1_)/2 * wt(j)
            endif
            if (r < 1.0d-13) cycle           ! When r is zero, the integrands are equal to zero.
            

c           ==== Calculate the necessary Hypergeometric functions ====
            call wfat_chgm_1pt(r, abs(m), lmax+1, a, b, lmb, chgm)     ! 1)
            call wfat_chgu_1pt(r, abs(m), lmax, a, b, lmb, chgu)     

            do l = abs(m), lmax
               f = (kappa*r)**l * exp(-kappa*r)
               g = f * chgm(l)
               h = f * chgu(l)

               call wfat_calc_radnorm(l, n, m, zq, kappa, c)
               call wfat_calc_radnorm(l+1, n, m, zq, kappa, c_p)
               k_p = sqrt(dble((l+1)**2-m**2) / (4*(l+1)**2-1))
               y = k_p * c_p*(kappa*r)**(l+1)*exp(-kappa*r) * chgm(l+1)
               if (l > abs(m)) then
                  call wfat_calc_radnorm(l-1, n, m, zq, kappa, c_m)
                  k_m = sqrt(dble(l**2-m**2) / (4*l**2-1))
                  y = y + k_m *
     &                c_m*(kappa*r)**(l-1)*exp(-kappa*r) * chgm(l-1)
               endif

c              ==== The g_l^\nu part ====
               ind_oi = g * r**3 * y
               ind_od = g * r**2 * c*g
               gg(1,l) = gg(1,l) + wt_*ind_oi
               gg(2,l) = gg(2,l) + wt_*ind_od

c              ==== The h_l^\nu part ====
               ind_oi = h * r**3 * y
               ind_od = h * r**2 * c*g
               hh(1,l) = hh(1,l) + wt_*ind_oi
               hh(2,l) = hh(2,l) + wt_*ind_od
            enddo

         enddo
      enddo
      

      end
c     NOTES:
c     1) lmax+1 instead of just lmax because the integrands involve the evaluation
c        of R_{l+1}^\nu, see below.
c=======================================================
