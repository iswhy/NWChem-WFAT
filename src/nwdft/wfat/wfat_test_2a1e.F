c==========================================================
c==========================================================
      subroutine wfat_test_2a1e(m, nr, nt)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer nr
      integer nt

      integer i, maxiter, l_cf, k_cf, l_eg, k_eg
      double precision q1, q2, z1, z2, ke


      ke = sqrt(2 * 0.5d0)
      write(luout, '("Initial ke = ", f18.8, " /Bohr")') ke
      call wfat_test_dvr_l0(nr, ke)
      call wfat_test_fbr_l0(nr, ke)
      
      write(luout, '(a)') 'Testing 2a1e ...'

      if (.not. ma_alloc_get(MT_DBL, (nr*nt)**2, 'hm', l_cf, k_cf))
     &     call errquit('WFAT: Cannot allocate hm', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nt, 'erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate erg', 1, MA_ERR)

c     ==== Adapt ke ====
      q1 = 1.0d0
      q2 = 1.0d0
      z1 = 1.0d0
      z2 = -1.0d0
      maxiter = 10
      call wfat_2a1e_solve(m, nr, nt, q1, q2, z1, z2, ke, .true.,
     &     maxiter, dbl_mb(k_cf), dbl_mb(k_eg))

c     ==== Print eigenvalues ====
      write(luout, '("Output ke = ", f18.8, " /Bohr")') ke
      write(luout, '("Eigenvalues for m = ", i3, " and with ke '
     &     //'adaptation:")') m
      do i = 1, nr*nt
         if (dbl_mb(k_eg+i-1) < 0.0d0)
     &        write(luout, '(i4, f18.8)') i, dbl_mb(k_eg+i-1)
      enddo
      write(luout, '("Eigenvalue range = ", f18.8, " H")')
     &     dbl_mb(k_eg+nr*nt-1) - dbl_mb(k_eg)


c     ==== Don't adapt ke ====
      q1 = 1.0d0
      q2 = 1.0d0
      z1 = 1.0d0
      z2 = -1.0d0
      maxiter = 10
      call wfat_2a1e_solve(0, nr, nt, q1, q2, z1, z2, ke, .false.,
     &     0, dbl_mb(k_cf), dbl_mb(k_eg))
      write(luout, '("The lowest eigenvalue for m = 0 and without ke '//
     &     'adaptation = ", f18.8, " H")') dbl_mb(k_eg)

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Cannot deallocate hm', 1, MA_ERR)
      if (.not. ma_free_heap(l_eg))
     &     call errquit('WFAT: Cannot deallocate erg', 1, MA_ERR)
      
      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_test_dvr_l0(nr, ke)
c     DESCRIPTION:
c     This subroutine tests the computation of the eigenvalues of the hydrogen
c     Hamiltonian matrix in the DVR bases.
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer nr
      double precision ke

      integer l_lr, l_ur, l_hm, l_x, l_rw, l_eg, l_po0, l_po
      integer k_lr, k_ur, k_hm, k_x, k_rw, k_eg, k_po0, k_po
      integer i, j, i_hm, i_po0
      double precision x, pot, normi, normj, v, laguerre_norm_i, mm
      external laguerre_norm_i


      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'lapr', l_lr, k_lr))
     &     call errquit('WFAT: Cannot allocate lapr', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'u_r', l_ur, k_ur))
     &     call errquit('WFAT: Cannot allocate u_r', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'g_r', l_hm, k_hm))
     &     call errquit('WFAT: Cannot allocate g_r', 1, MA_ERR)


      call wfat_laplr_fbrr(nr, 0, dbl_mb(k_lr))
      call wfat_dvr_fbr_r(0, nr, dbl_mb(k_ur))
c     T_r^{DVR} = transpose(U_r) . T_r^{FBR} . U_r
      call wfat_represent_sym(nr, nr, dbl_mb(k_ur), dbl_mb(k_lr),
     &     dbl_mb(k_hm))
      dbl_mb(k_hm:k_hm+nr**2-1) = -0.5d0 * dbl_mb(k_hm:k_hm+nr**2-1)
     &     / ke
c      write(luout, '("Kinetic energy matrix:")')
c      call wfat_dprintmat(nr, nr, 'f', dbl_mb(k_hm))


      if (.not. ma_alloc_get(MT_DBL, nr, 'xquad', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr, 'xweight', l_rw, k_rw))
     &     call errquit('WFAT: Cannot allocate xweight', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr, 'erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate erg', 1, MA_ERR)
      

      call wfat_lf_quadrature_rule(nr, 0+2.0d0, dbl_mb(k_x),
     &     dbl_mb(k_rw))
      do i = 1, nr
         x = dbl_mb(k_x+i-1)
         pot = - 1/x / ke**2
         i_hm = k_hm + (i-1)*nr + (i-1)
         dbl_mb(i_hm) = dbl_mb(i_hm) + pot
      enddo


c     ================================================================
c     ==== The block of code below computes the potential matrix  ====
c     ==== elements analytically instead of numerically, that is, ====
c     ==== it does not use the special property of the DVR bases  ====
c     ====  in which the matrix of a local operator is diagonal.  ====
c     ====  Remove the comment below but comment out the entire   ====
c     ====  loop above if you want to use the analytic potential  ====
c     ==== matrix. Using the analytic potential matrix will give  ====
c     ====  a more accurate result which should be identical to   ====
c     ==== the results obtained by diagonalizing in the FBR basis ====
c     ==== (which is done by the wfat_test_fbr_l0 routine) at the ====
c     ====          expense of a longer simulation time.          ====
c     ================================================================
cc      if (.not. ma_alloc_get(MT_DBL, nr**2, 'pot0', l_po0, k_po0))
cc     &     call errquit('WFAT: Cannot allocate pot0', 1, MA_ERR)
cc      if (.not. ma_alloc_get(MT_DBL, nr**2, 'pot', l_po, k_po))
cc     &     call errquit('WFAT: Cannot allocate pot', 1, MA_ERR)
cc      mm = dble(0)
cc      do j = 1, nr
cc         normj = sqrt(laguerre_norm_i(j-1,0+2))
cc         do i = 1, nr
cc            normi = sqrt(laguerre_norm_i(i-1,0+2))
cc            call wfat_laguerre_intg1(mm+2, mm+2, mm+1, i-1, j-1, v)
cc            i_po0 = k_po0 + (j-1)*nr + (i-1)
cc            dbl_mb(i_po0) = - v / normi / normj / ke**2
cc         enddo
cc      enddo
cc      call wfat_represent_sym(nr, nr, dbl_mb(k_ur), dbl_mb(k_po0),
cc     &     dbl_mb(k_po))
cc
cc      dbl_mb(k_hm:k_hm+nr**2-1) = dbl_mb(k_hm:k_hm+nr**2-1) +
cc     &                            dbl_mb(k_po:k_po+nr**2-1)
      
      
      call wfat_dsyev_wrap(nr, dbl_mb(k_hm), dbl_mb(k_eg))
      write(luout, '("Hydrogen eigenvalues for m = 0 using DVR :")')
      do i = 1, nr
         if (dbl_mb(k_eg+i-1) < 0.0d0)
     &        write(luout, '(i4, f18.8)') i, dbl_mb(k_eg+i-1)
      enddo


      if (.not. ma_free_heap(l_lr))
     &     call errquit('WFAT: Cannot deallocate lapr', 1, MA_ERR)
      if (.not. ma_free_heap(l_ur))
     &     call errquit('WFAT: Cannot deallocate u_r', 1, MA_ERR)
      if (.not. ma_free_heap(l_hm))
     &     call errquit('WFAT: Cannot deallocate hm', 1, MA_ERR)
      if (.not. ma_free_heap(l_x))
     &     call errquit('WFAT: Cannot deallocate xquad', 1, MA_ERR)
      if (.not. ma_free_heap(l_rw))
     &     call errquit('WFAT: Cannot deallocate xweight', 1, MA_ERR)
      if (.not. ma_free_heap(l_eg))
     &     call errquit('WFAT: Cannot deallocate erg', 1, MA_ERR)

      end
c==========================================================


c==========================================================
c==========================================================
      subroutine wfat_test_fbr_l0(nr, ke)
c     DESCRIPTION:
c     This subroutine tests the computation of the eigenvalues of the hydrogen
c     Hamiltonian matrix in the FBR bases.
      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
      integer nr
      double precision ke

      integer l_hm, k_hm, l_eg, k_eg, i, j, i_hm
      double precision mm, normi, normj, pot, v, laguerre_norm_i
      external laguerre_norm_i


      ke = sqrt(2 * 0.5d0)
      mm = dble(0)
      
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'hm', l_hm, k_hm))
     &     call errquit('WFAT: Cannot allocate hm', 1, MA_ERR)

      
      call wfat_laplr_fbrr(nr, 0, dbl_mb(k_hm))
      dbl_mb(k_hm:k_hm+nr**2-1) = -0.5d0 * dbl_mb(k_hm:k_hm+nr**2-1)
     &     / ke

      do j = 1, nr
         normj = sqrt(laguerre_norm_i(j-1,0+2))
         do i = 1, nr
            normi = sqrt(laguerre_norm_i(i-1,0+2))
            call wfat_laguerre_intg1(mm+2, mm+2, mm+1, i-1, j-1, v)
            pot = - v / normi / normj / ke**2
            i_hm = k_hm + (j-1)*nr + (i-1)
            dbl_mb(i_hm) = dbl_mb(i_hm) + pot
         enddo
      enddo


      if (.not. ma_alloc_get(MT_DBL, nr, 'erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate erg', 1, MA_ERR)
      call wfat_dsyev_wrap(nr, dbl_mb(k_hm), dbl_mb(k_eg))
      write(luout, '("Hydrogen eigenvalues with m = 0 using FBR :")')
      do i = 1, nr
         if (dbl_mb(k_eg+i-1) < 0.0d0)
     &        write(luout, '(i4, f18.8)') i, dbl_mb(k_eg+i-1)
      enddo


c     ==== Deallocations ====
      if (.not. ma_free_heap(l_hm))
     &     call errquit('WFAT: Cannot deallocate hm', 32, MA_ERR)
      
      end
c==========================================================
      
      
