c=======================================================
c=======================================================
      subroutine wfat_test_r1(n, m, erg, zq, chmmt)

      implicit none
#include "mafdecls.fh"
#include "stdio.fh"      
#include "errquit.fh"
#include "wfat.fh"
      integer n
      integer m
      double precision erg
      double precision zq
      double precision chmmt(3)

      character*255 fname
      logical first_notsmalleta
      integer l_pt, k_pt, l_wt, k_wt, l_r1, k_r1, l_cr, k_cr
      integer pmax, nlgl, npt, i, j, id
      double precision kappa, a, b, et_first, et_last, de, et, del, z
      double precision chgm, chgu, mu(3), g, d, f, r0_1, o0_1, rr, r0
      double precision r0_2, o0_2, r1_, protmat(3,3), beta1, beta0
      double precision cr_1_ref, f2h, f1h, f0h, f1h_, f2h_, dif2, lhs
      double precision rhs, epsi
      double complex zz, chgm_
      double precision factorial
      external factorial


c     ==== Get the passive rotation matrix ===
      call wfat_passrot(-0.0d0, -0.0d0, protmat)
      mu = matmul(transpose(protmat), -chmmt)
      
      
      pmax = 10
      nlgl = 9
      npt = 51
      kappa = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta0 = zq - kappa*(n + (abs(m)+1)/2.0d0)
      beta1 = -mu(3)/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))
      et_first = 0.0d0
      et_last = 10.0d0
      de = (et_last-et_first) / (npt-1)
      
      write(luout, '("TESTING FOR R^{(1)}")')
      write(luout, '("pmax = ", i4)') pmax
      write(luout, '("nlgl = ", i4)') nlgl
      write(luout, '("npt = ", i4)') npt
      write(luout, '("(n,m) = ", 2i4)') n, m
      write(luout, '("a, b = ", 2es15.6)') a, b
      write(luout, '("beta1 = ", es15.6)') beta1
      write(luout, '("gamma(a) = ", es15.6)') g
      write(luout, '("et_0, et_n = ", 2es15.6)') et_first, et_last

      

      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_pt', l_pt, k_pt))
     &     call errquit('WFAT: lgl_pt allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_wt', l_wt, k_wt))
     &     call errquit('WFAT: lgl_wt allocation failed', 85, MA_ERR)
      call wfat_lgl(-1.0d0, 1.0d0, nlgl-1, dbl_mb(k_pt), dbl_mb(k_wt))
      write(luout, '("Legendre-Gauss-Lobatto nodes and weights:")')
      do i = 1, nlgl
         write(luout, '(i5, 2es15.6)')
     &        i, dbl_mb(k_pt+i-1), dbl_mb(k_wt+i-1)
      enddo
      if (.not. MA_alloc_get(MT_DBL, npt, 'r1', l_r1, k_r1))
     &     call errquit('WFAT: r1 allocation failed', 86, MA_ERR)      

      
      cr_1_ref = -g/kappa/factorial(abs(m)+1) * beta1
      if (.not. MA_alloc_get(MT_DBL, pmax, 'r1_coef', l_cr, k_cr))
     &     call errquit('WFAT: r1_coef allocation failed', 62, MA_ERR)
      call wfat_r1_coeff(n, m, erg, zq, mu(3), pmax, dbl_mb(k_cr))
      write(luout, '("Coefficients of R^{(1)} :")')
      do j = 1, pmax
         if (j <= 1) then
            write(luout, '(i3, 2es15.6)') j, dbl_mb(k_cr+j-1), cr_1_ref
         else
            write(luout, '(i3, es15.6)') j, dbl_mb(k_cr+j-1)
         endif
      enddo


      r1_ = 0.0d0
      first_notsmalleta = .true.      
      call util_file_name('test.wr1', .false., .false., fname)
      open(WF_R1_UNIT, file=fname, action='write')
      write(WF_R1_UNIT, '(a5, 4a15)')
     &        'No.', 'eta', 'R^(1)', 'error', 'r1_'
      do i = 1, npt
         et = et_first + (i-1)*de
         if (et < 1.0d0) then
            z = kappa * et
            rr = 0.0d0
            do j = 1, pmax
               rr = rr + dbl_mb(k_cr+j-1) * z**j
            enddo
            dbl_mb(k_r1+i-1) = z**(b/2) * exp(-z/2) * rr
            del = z**(b/2) * exp(-z/2) * dbl_mb(k_cr+pmax-1) * z**pmax
            
c            call wfat_r1_smalleta(et, n, m, erg, zq, mu(3), pmax,
c     &           del, dbl_mb(k_r1+i-1))
         else
            if (first_notsmalleta) then
               z = kappa * (et-de)
               zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
               call wfat_cchg(a, b, zz, chgm_)
               chgm = real(chgm_)
               call wfat_chgu(a, b, z, chgu, id)

               
               f = z**(b/2) * exp(-z/2)
               r0_1 = d * f * chgm
               o0_1 = f * chgu
               first_notsmalleta = .false.
            endif

cOLD            call wfat_r1_trapz(et-de, et, n, m, erg, zq, mu(3), nlgl,
cOLD     &           r1_)


            
            resolve here


            
cOLD            call wfat_r1(et-de, et, r0_1, o0_1, n, m, erg, zq, mu(3),
cOLD     &           nlgl, dbl_mb(k_pt), dbl_mb(k_wt), r0_2, o0_2, r1_)
            dbl_mb(k_r1+i-1) = dbl_mb(k_r1+i-2) + r1_
            r0_1 = r0_2
            o0_1 = o0_2
            del = 0.0d0
         endif

         write(WF_R1_UNIT, '(i5, 4es15.6)')
     &        i, et, dbl_mb(k_r1+i-1), del, r1_
      enddo
      close(WF_R1_UNIT)


c     ==== Check differential equation ====
      call util_file_name('test.wdr1', .false., .false., fname)
      open(WF_DIFR1_UNIT, file=fname, action='write')
      write(WF_DIFR1_UNIT, '(a5, 5a15)')
     &     'No.', 'eta', 'LHS', 'RHS', 'Epsilon', 'R_\nu^(0)'
      do i = 1, npt-2
         if (i <= 2) then
            dif2 = 0.0d0
         else
            f2h = dbl_mb(k_r1+(i+2)-1)
            f1h = dbl_mb(k_r1+(i+1)-1)
            f0h = dbl_mb(k_r1+(i+0)-1)
            f1h_ = dbl_mb(k_r1+(i-1)-1)
            f2h_ = dbl_mb(k_r1+(i-2)-1)
            dif2 = (-f2h + 16*f1h - 30*f0h + 16*f1h_ - f2h_)
     &             / (12*de**2)
         endif

         et = et_first + (i-1)*de
         z = kappa * et
         zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
         call wfat_cchg(a, b, zz, chgm_)
         chgm = real(chgm_)
         r0 = d * z**(b/2) * exp(-z/2) * chgm
               
         lhs = dif2 +
     &        ((1-m**2)/(4*et**2) + beta0/et + erg/2)*dbl_mb(k_r1+i-1)
         rhs = -(beta1/et - mu(3)/2 + et/4) * r0
         epsi = lhs - rhs

         write(WF_DIFR1_UNIT, '(i5, 5es15.6)') i, et, lhs, rhs, epsi, r0
      enddo
      close(WF_DIFR1_UNIT)
      
         
      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_pt))
     &     call errquit('WFAT: Cannot deallocate lgl_pt', 49, MA_ERR)
      if (.not. ma_free_heap(l_wt))
     &     call errquit('WFAT: Cannot deallocate lgl_wt', 49, MA_ERR)
      if (.not. ma_free_heap(l_r1))
     &     call errquit('WFAT: Cannot deallocate r1', 49, MA_ERR)
      if (.not. ma_free_heap(l_cr))
     &     call errquit('WFAT: Cannot deallocate r1_coef', 49, MA_ERR)

      end
c=======================================================
