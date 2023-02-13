c=======================================================
c=======================================================
      subroutine wfat_calc_tau0(n, m, npt, r, th, ph, erg, zq, tau0)

      implicit none
#include "errquit.fh"
      integer n                     ! [input], 
      integer m                     ! [input], 
      integer npt                   ! [input], 
      double precision r(npt)       ! [input], r coordinate in a.u.
      double precision th(npt)      ! [input], theta coordinate in radians.
      double precision ph(npt)      ! [input], phi coordinate in radians.
      double precision erg          ! [input], 
      double precision zq           ! [input], 
      double complex tau0(npt)      ! [output],

      integer i
      double precision kappa, a, b, g, d, eta, r0, xi, phi0, z
      double precision beta0, pi_const, factorial
      double complex zz, chgm_, imnum
      external factorial
      parameter ( pi_const = 3.14159265359d0 )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )
      

      kappa = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))
      

c     ==== Calculate tau_\nu^{(0)} ====
      do i = 1,npt
c        == R_\nu^{(0)} ===
         eta = r(i) * (1 - cos(th(i)))
         z = kappa * eta
         zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
         call wfat_cchg(a, b, zz, chgm_)
         r0 = d * z**(b/2) * exp(-z/2) * real(chgm_)         

c        == phi_\nu^{(0)} ==
         xi = r(i) * (1 + cos(th(i)))         
         call wfat_calc_phi0xi(xi, n, m, erg, phi0)

         tau0 = r0 * phi0 * exp(imnum*m*ph(i)) / sqrt(2*pi_const*eta)
      enddo

      beta0 = zq - kappa*(n + (abs(m)+1)/2.0d0)
      tau0 = -2*kappa**(beta0/kappa-1) * tau0
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_tau0_lexp(n, m, npt, r, th, ph, erg, zq,
     &           lmax, tau0)

      implicit none
#include "stdio.fh"      
#include "errquit.fh"
#include "mafdecls.fh"
      integer n                     ! [input],
      integer m                     ! [input],
      integer npt                   ! [input],
      double precision r(npt)       ! [input], r coordinate in a.u.
      double precision th(npt)      ! [input], theta coordinate in radians.
      double precision ph(npt)      ! [input], phi coordinate in radians.
      double precision erg          ! [input],
      double precision zq           ! [input],
      integer lmax                  ! [input],
      double complex tau0(npt)      ! [output],

      integer i, nsh, nl, l, i_rf, i_sh
      integer l_sh, k_sh, l_hf, k_hf, l_rf, k_rf
      double precision kappa, lmb, a, b
      double complex xx

      
      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0

      
c     ==== Allocations ====
      nsh = (lmax + 1)**2
      if (.not. MA_alloc_get(MT_DCPL, nsh, 'sphm', l_sh, k_sh))
     &     call errquit('WFAT: sphm allocation failed', 59, MA_ERR)
      nl = lmax - abs(m) + 1
      if (.not. MA_alloc_get(MT_DBL, nl, 'chgm', l_hf, k_hf))
     &     call errquit('WFAT: chgm allocation failed', 29, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nl, 'radf', l_rf, k_rf))
     &     call errquit('WFAT: radf allocation failed', 29, MA_ERR)
      

c     ==== Calculate tau_\nu^{(0)} ====
      do i = 1, npt
         call wfat_chgm_1pt(r(i), abs(m), lmax, a, b, lmb, dbl_mb(k_hf))
         call wfat_radf_1pt(n, m, int(zq), lmax, kappa, r(i), 
     &        dbl_mb(k_hf), dbl_mb(k_rf))
         call wfat_sphm_1pt(lmax, th(i), ph(i), dcpl_mb(k_sh))
         
         tau0(i) = 0.0d0
         do l = abs(m), lmax
            i_rf = k_rf + l - abs(m)
            i_sh = k_sh + (l**2 + l + m + 1) - 1
            xx = dbl_mb(i_rf) * dcpl_mb(i_sh)
            tau0(i) = tau0(i) + xx

            write(luout, '("l, R_l, Y_lm = ", i5, 3es15.6)')
     &           l, dbl_mb(i_rf), dcpl_mb(i_sh)
         enddo
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_sh))
     &     call errquit('WFAT: sphm deallocation failed', 71, MA_ERR)
      if (.not. MA_free_heap(l_hf))
     &     call errquit('WFAT: chgm deallocation failed', 71, MA_ERR)
      if (.not. MA_free_heap(l_rf))
     &     call errquit('WFAT: radf deallocation failed', 71, MA_ERR)

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_calc_tau1(n, m, npt, r, th, ph, erg, zq, bt, gm,
     &           chmmt, tau1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer n                     ! [input], 
      integer m                     ! [input], 
      integer npt                   ! [input], 
      double precision r(npt)       ! [input], r coordinate in a.u.
      double precision th(npt)      ! [input], theta coordinate in radians.
      double precision ph(npt)      ! [input], phi coordinate in radians.
      double precision erg          ! [input], 
      double precision zq           ! [input], 
      double precision bt           ! [input], beta orientation angle in degrees. 
      double precision gm           ! [input], gamma orientation angle in degrees.
      double precision chmmt(3)     ! [input], 
      double complex tau1(npt)      ! [output], 

      integer i, l_r1, k_r1, l_et, k_et
      double precision kappa, a, b, g, d, eta, r0, r1, xi, phi0, phi1, z
      double precision beta0, pi_const, factorial
      double complex zz, chgm_, imnum
      external factorial
      parameter ( pi_const = 3.14159265359d0 )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )


      if (npt /= 1)
     &     call errquit('WFAT: For the moment, the argument npt of '//
     &     'wfat_calc_tau1 must be equal to one', 27, UERR)
cforlater      call wfat_sort_asc()     resolve here


      kappa = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))


c     ==== Obtain eta ====
      if (.not. MA_alloc_get(MT_DBL, npt, 'r1', l_r1, k_r1))
     &     call errquit('WFAT: r1 allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'eta', l_et, k_et))
     &     call errquit('WFAT: eta allocation failed', 84, MA_ERR)
      do i = 1, npt
         dbl_mb(k_et+i-1) = r(i) * (1 - cos(th(i)))
      enddo


c     ==== Calculate R_\nu^{(1)} ====
      call wfat_r1(npt, dbl_mb(k_et), n, m, erg, zq, bt, gm, chmmt,
     &     dbl_mb(k_r1))


c     ==== Calculate tau_\nu^{(1)} ====
      do i = 1, npt
c        == R_\nu^{(0)} ===
         eta = dbl_mb(k_et+i-1)
         z = kappa * eta
         zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
         call wfat_cchg(a, b, zz, chgm_)
         r0 = d * z**(b/2) * exp(-z/2) * real(chgm_)

c        == R_\nu^{(1)} ==
         r1 = dbl_mb(k_r1+i-1)

c        == phi_\nu^{(0)} and phi_\nu^{(1)} ==
         xi = r(i) * (1 + cos(th(i)))
         call wfat_calc_phi0xi(xi, n, m, erg, phi0)
         call wfat_calc_phi1xi(xi, n, m, erg, bt, gm, chmmt, phi1)

         tau1(i) = (r0 * phi1 + r1 * phi0) * exp(imnum*m*ph(i)) /
     &             sqrt(2*pi_const*eta)
      enddo
      beta0 = zq - kappa*(n + (abs(m)+1)/2.0d0)
      tau1 = -2*kappa**(beta0/kappa-1) * tau1

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_r1))
     &     call errquit('WFAT: r1 deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_et))
     &     call errquit('WFAT: eta deallocation failed', 46, MA_ERR)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_tau1_lexp(n, m, npt, r, th, ph, erg, zq,
     &           bt, gm, chmmt, lmax, tau1)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer n                     ! [input], 
      integer m                     ! [input], 
      integer npt                   ! [input], 
      double precision r(npt)       ! [input], r coordinate in a.u.
      double precision th(npt)      ! [input], theta coordinate in radians.
      double precision ph(npt)      ! [input], phi coordinate in radians.
      double precision erg          ! [input], 
      double precision zq           ! [input], 
      double precision bt           ! [input], beta orientation angle in degrees. 
      double precision gm           ! [input], gamma orientation angle in degrees.
      double precision chmmt(3)     ! [input], 
      integer lmax                  ! [input], 
      double complex tau1(npt)      ! [output],

      integer nlgl, nsh, nl, i, l, i_q1, i_q2, i_sh
      integer l_sh, k_sh, l_q1, k_q1, l_q2, k_q2
      double precision mu(3), protmat(3,3)
      double complex xx
      parameter(nlgl=8)


      if (npt /= 1)
     &     call errquit('WFAT: For the moment, the argument npt of '//
     &     'wfat_calc_tau1_lexp must be equal to one', 27, UERR)
cforlater      call wfat_sort_asc()     resolve here

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(-bt, -gm, protmat)
      mu = matmul(transpose(protmat), -chmmt)      


c     ==== Calculate Q_l's ====
      nsh = (lmax + 1)**2
      if (.not. MA_alloc_get(MT_DCPL, nsh, 'sphm', l_sh, k_sh))
     &     call errquit('WFAT: sphm allocation failed', 59, MA_ERR)
      nl = lmax - abs(m) + 1
      if (.not. MA_alloc_get(MT_DBL, nl*npt, 'q1', l_q1, k_q1))
     &     call errquit('WFAT: q1 allocation failed', 29, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nl*npt, 'q2', l_q2, k_q2))
     &     call errquit('WFAT: q2 allocation failed', 29, MA_ERR)
      call wfat_calc_q(npt, r, n, m, lmax, nlgl, erg, zq,
     &     dbl_mb(k_q1), dbl_mb(k_q2))

      
c     ==== Calculate tau_\nu^{(1)} ====
      do i = 1, npt
         call wfat_sphm_1pt(lmax, th(i), ph(i), dcpl_mb(k_sh))

         tau1(i) = 0.0d0
         do l = abs(m), lmax
            i_q1 = k_q1 + (i-1)*nl + l - abs(m)
            i_q2 = k_q2 + (i-1)*nl + l - abs(m)
            i_sh = k_sh + (l**2 + l + m + 1) - 1
            xx = (dbl_mb(i_q1) + mu(3)*dbl_mb(i_q2)) * dcpl_mb(i_sh)
            tau1(i) = tau1(i) + xx

            write(luout, '("l, Q_l, Y_lm = ", i5, 3es15.6)')
     &           l, (dbl_mb(i_q1) + mu(3)*dbl_mb(i_q2)), dcpl_mb(i_sh)
         enddo
      enddo
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_sh))
     &     call errquit('WFAT: sphm deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_q1))
     &     call errquit('WFAT: q1 deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_q2))
     &     call errquit('WFAT: q2 deallocation failed', 46, MA_ERR)
      end
c=======================================================
