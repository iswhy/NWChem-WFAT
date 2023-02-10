c=======================================================
c=======================================================
      subroutine wfat_test_tau(n, m, erg, zq, lmax, chmmt)

      implicit none
c#include "wfat.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      integer lmax
      double precision chmmt(3)

      integer npt, i, i_r, i_th, i_ph, i_t0, i_t0l, i_t1, i_t1l
      integer l_r, l_th, l_ph, l_t0, l_t0l, l_t1, l_t1l
      integer k_r, k_th, k_ph, k_t0, k_t0l, k_t1, k_t1l
      double precision kappa, eta, xi, phi
      parameter(npt=1)


      kappa = sqrt(2 * abs(erg))
      
      
c     ==== Allocate r, theta, and phi ====
      if (.not. MA_alloc_get(MT_DBL, npt, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'theta', l_th, k_th))
     &     call errquit('WFAT: theta allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'phi', l_ph, k_ph))
     &     call errquit('WFAT: phi allocation failed', 84, MA_ERR)
      dbl_mb(k_r) = 0.5d0
      dbl_mb(k_th) = 1.2d0
      dbl_mb(k_ph) = 2.0d0

      
c     ==== Allocate tau's ====
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau0', l_t0, k_t0))
     &     call errquit('WFAT: tau0 allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau0_l', l_t0l, k_t0l))
     &     call errquit('WFAT: tau0_l allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau1', l_t1, k_t1))
     &     call errquit('WFAT: tau1 allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, npt, 'tau1_l', l_t1l, k_t1l))
     &     call errquit('WFAT: tau1_l allocation failed', 84, MA_ERR)


c     ==== Compare tau_\nu^{(0)} ====
      call wfat_calc_tau0(n, m, npt, dbl_mb(k_r), dbl_mb(k_th),
     &     dbl_mb(k_ph), erg, zq, dcpl_mb(k_t0))
      call wfat_calc_tau0_lexp(n, m, npt, dbl_mb(k_r), dbl_mb(k_th),
     &     dbl_mb(k_ph), erg, zq, lmax, dcpl_mb(k_t0l))


c     ==== Compare tau_\nu^{(1)} ====
      call wfat_calc_tau1(n, m, npt, dbl_mb(k_r), dbl_mb(k_th),
     &     dbl_mb(k_ph), erg, zq, bt, gm, chmmt, dcpl_mb(k_t1))
      call wfat_calc_tau1_lexp(n, m, npt, dbl_mb(k_r), dbl_mb(k_th),
     &     dbl_mb(k_ph), erg, zq, bt, gm, chmmt, lmax, dcpl_mb(k_t1l))


      do i = 1, npt
c        == Print r, theta, phi ==
         i_r = k_r + i - 1
         i_th = k_th + i - 1
         i_ph = k_ph + i - 1
         write(luout, '(3es15.6)', advance='no')
     &        dbl_mb(i_r), dbl_mb(i_th), dbl_mb(i_ph)

c        == Print eta, xi, phi ==
         eta = dbl_mb(i_r) * (1 - cos( dbl_mb(i_th) ))
         xi = dbl_mb(i_r) * (1 + cos( dbl_mb(i_th) ))
         phi = dbl_mb(i_ph)
         write(luout, '(3es15.6)', advance='no') eta, xi, phi

c        == Print kappa*r, kappa*eta ==
         write(luout, '(2es15.6)', advance='no')
     &        kappa*dbl_mb(i_r), kappa*eta

c        == Print tau0 and tau1 ==
         i_t0 = k_t0 + i - 1
         i_t0l = k_t0l + i - 1
         i_t1 = k_t1 + i - 1
         i_t1l = k_t1l + i - 1
         write(luout, '(2(3x, 2es21.12))', advance='no')
     &        dcpl_mb(i_t0), dcpl_mb(i_t0l)
         write(luout, '(6x, " ")', advance='no')
         write(luout, '(2(3x, 2es21.12))', advance='no')
     &        dcpl_mb(i_t1), dcpl_mb(i_t1l)

         write(luout, '(" ")')
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_r))
     &     call errquit('WFAT: r deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_th))
     &     call errquit('WFAT: theta deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_ph))
     &     call errquit('WFAT: phi deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t0))
     &     call errquit('WFAT: tau0 deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t0))
     &     call errquit('WFAT: tau0_l deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t1))
     &     call errquit('WFAT: tau1 deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t1))
     &     call errquit('WFAT: tau1_l deallocation failed', 46, MA_ERR)
      
      end
c=======================================================
      
