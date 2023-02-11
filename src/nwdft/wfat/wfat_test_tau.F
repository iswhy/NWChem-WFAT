c=======================================================
c=======================================================
      subroutine wfat_test_tau(n, m, erg, zq, lmax, chmmt)

      implicit none
#include "wfat.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      integer lmax
      double precision chmmt(3)

      character*255 fname
      integer npt, i, i_r, i_th, i_ph, i_t0, i_t0l, i_t1, i_t1l
      integer l_r, l_th, l_ph, l_t0, l_t0l, l_t1, l_t1l
      integer k_r, k_th, k_ph, k_t0, k_t0l, k_t1, k_t1l
      double precision kappa, eta, xi, phi, bt, gm
      parameter(npt=1, bt=0.0d0, gm=0.0d0)



      write(luout, '("TESTING FOR TAU")')
      kappa = sqrt(2 * abs(erg))
      
      
c     ==== Allocate r, theta, and phi ====
      if (.not. MA_alloc_get(MT_DBL, npt, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'theta', l_th, k_th))
     &     call errquit('WFAT: theta allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, npt, 'phi', l_ph, k_ph))
     &     call errquit('WFAT: phi allocation failed', 84, MA_ERR)
      dbl_mb(k_r) = 0.5d0
      dbl_mb(k_th) = 0.8d0
      dbl_mb(k_ph) = 0.0d0

      
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


c     ==== Open the test result file and print the header ====
      call util_file_name('test.wtau', .false., .false., fname)
      open(WF_TAUTEST_UNIT, file=fname, action='write')
      
      write(WF_TAUTEST_UNIT, '("n = ", i5)') n
      write(WF_TAUTEST_UNIT, '("m = ", i5)') m
      write(WF_TAUTEST_UNIT, '("energy = ", es15.6)') erg
      write(WF_TAUTEST_UNIT, '("ion charge = ", es15.6)') zq
      write(WF_TAUTEST_UNIT, '("lmax = ", i5)') lmax
      write(WF_TAUTEST_UNIT, '("orbital moment (x y z) = ", 3es15.6)')
     &     chmmt
      write(WF_TAUTEST_UNIT, '("beta, gamma (deg) = ", 2es15.6)') bt, gm

      
      call wfat_test_tau_header()      
      do i = 1, npt
c        == Print r, theta, phi ==
         i_r = k_r + i - 1
         i_th = k_th + i - 1
         i_ph = k_ph + i - 1
         write(WF_TAUTEST_UNIT, '(3es15.6)', advance='no')
     &        dbl_mb(i_r), dbl_mb(i_th), dbl_mb(i_ph)

c        == Print eta, xi, phi ==
         eta = dbl_mb(i_r) * (1 - cos( dbl_mb(i_th) ))
         xi = dbl_mb(i_r) * (1 + cos( dbl_mb(i_th) ))
         phi = dbl_mb(i_ph)
         write(WF_TAUTEST_UNIT, '(3es15.6)', advance='no') eta, xi, phi

c        == Print kappa*r, kappa*eta ==
         write(WF_TAUTEST_UNIT, '(2es15.6)', advance='no')
     &        kappa*dbl_mb(i_r), kappa*eta

c        == Print tau0 and tau1 ==
         i_t0 = k_t0 + i - 1
         i_t0l = k_t0l + i - 1
         i_t1 = k_t1 + i - 1
         i_t1l = k_t1l + i - 1
         write(WF_TAUTEST_UNIT, '(2(3x, 2es21.12))', advance='no')
     &        dcpl_mb(i_t0), dcpl_mb(i_t0l)
         write(WF_TAUTEST_UNIT, '(6x, " ")', advance='no')
         write(WF_TAUTEST_UNIT, '(2(3x, 2es21.12))', advance='no')
     &        dcpl_mb(i_t1), dcpl_mb(i_t1l)

         write(WF_TAUTEST_UNIT, '(" ")')
      enddo

      close(WF_TAUTEST_UNIT)
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_r))
     &     call errquit('WFAT: r deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_th))
     &     call errquit('WFAT: theta deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_ph))
     &     call errquit('WFAT: phi deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t0))
     &     call errquit('WFAT: tau0 deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t0l))
     &     call errquit('WFAT: tau0_l deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t1))
     &     call errquit('WFAT: tau1 deallocation failed', 46, MA_ERR)
      if (.not. MA_free_heap(l_t1l))
     &     call errquit('WFAT: tau1_l deallocation failed', 46, MA_ERR)
      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_test_tau_header()

      implicit none
#include "wfat.fh"
#include "errquit.fh"

      logical isopen


c     ==== Check if the unit is opened ====
      inquire(unit=WF_TAUTEST_UNIT, opened=isopen)
      if (.not. isopen)
     &     call errquit('WFAT: WF_TAUTEST_UNIT must first be opened '//
     &     'before calling wfat_test_tau_header', 1, UERR)
      

c     ==== Column number ====
      write(WF_TAUTEST_UNIT, '("#", i14, 2i15)', advance='no') 1, 2, 3     ! r, theta, phi
      write(WF_TAUTEST_UNIT, '(3i15)', advance='no') 4, 5, 6               ! eta, xi, phi
      write(WF_TAUTEST_UNIT, '(2i15)', advance='no') 7, 8                  ! kappa*r, kappa*eta
      write(WF_TAUTEST_UNIT, '(2(3x, 2i21))', advance='no') 9, 10          ! tau0, tau0_pw
      write(WF_TAUTEST_UNIT, '(6x, " ")', advance='no')
      write(WF_TAUTEST_UNIT, '(2(3x, 2i21))', advance='no') 11, 12         ! tau1, tau1_pw
      write(WF_TAUTEST_UNIT, '(" ")')
      

c     ==== Column title ====
      write(WF_TAUTEST_UNIT, '("#", a14, 2a15)', advance='no')
     &     'r', 'theta', 'phi'
      write(WF_TAUTEST_UNIT, '(3a15)', advance='no')
     &     'eta', 'xi', 'phi'
      write(WF_TAUTEST_UNIT, '(2a15)', advance='no')
     &     'kappa*r', 'kappa*eta'
      write(WF_TAUTEST_UNIT, '(2(3x, 2a21))', advance='no')
     &     'real tau0', 'imag tau0', 'real tau0_pw', 'imag tau0_pw'
      write(WF_TAUTEST_UNIT, '(6x, " ")', advance='no')
      write(WF_TAUTEST_UNIT, '(2(3x, 2a21))', advance='no')
     &     'real tau1', 'imag tau1', 'real tau1_pw', 'imag tau1_pw'
      write(WF_TAUTEST_UNIT, '(" ")')

      end
c=======================================================
      
