c=======================================================
c=======================================================
      subroutine wfat_r1_smalleta(et, n, m, erg, zq, mu, pmax,
     &           last, r1)

      implicit none
#include "errquit.fh"      
      double precision et     ! [input]
      integer n               ! [input]
      integer m               ! [input]
      double precision erg    ! [input]
      double precision zq     ! [input]
      double precision mu     ! [input]
      integer pmax            ! [input]
      double precision last   ! [output]
      double precision r1     ! [output]

      logical l_i, l_j, l_k
      integer i, j, k, k_, ii, jj, kk, s1, s2, ubound, sg
      double precision kappa, a, b, z, beta1, ga, gm, mk, w_, t_, u_, v_
      double precision w, t, u, v, dk, ps1, ps2, ps3
      double precision r1_a, r1_b, last_a, last_b
      double precision factorial, rise_fact
      external factorial, rise_fact
      

      if (pmax < 1) call errquit('Argument pmax of '//
     &     'wfat_r1_smalleta must be larger than 1.', 28, UERR)
      if (n < 0) call errquit('Argument n of '//
     &     'wfat_r1_smalleta must be zero or positive.', 28, UERR)
      
      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      z = kappa * et
      beta1 = -mu/2/kappa * (2*n+1+abs(m))
     &        -1/4/kappa**2 * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      call wfat_gammaf(a, ga)


c     ==== Finite k part ====
      r1_b = 0.0d0
      last_b = 0.0d0
      if (abs(m) > 0) then
         ubound = pmax - 1
         do i = 0, ubound
            do j = 0, ubound
               do k = abs(m), 1, -1
                  do k_ = 0, ubound
                     s2 = i + j - k + k_
                     ii = abs(m) + s2 + 1
                     jj = abs(m) + s2 + 2
                     kk = abs(m) + s2 + 3
                     l_i = (ii >= 1 .and. ii <= pmax)
                     l_j = (jj >= 2 .and. jj <= pmax)
                     l_k = (kk >= 3 .and. kk <= pmax)

                     if (l_i .or. l_j .or. l_k) then
                        w_ = (-1)**i/factorial(i)
     &                       *( rise_fact(a,j)/factorial(j)
     &                       /rise_fact(b,j) )
     &                       *( rise_fact(a,k_)/factorial(k_)
     &                       /rise_fact(b,k_) )
     &                       *( factorial(k-1)*rise_fact(1-a+k,abs(m)-k)
     &                       /factorial(abs(m)-k) )
     &                       *(k+k_) 
                     endif

                     if (l_i) then
                        sg = abs(m) + 1 + i + j
                        mk = 1.0d0 / (sg+k_) / (sg-k)
                        t_ = w_ * beta1 * kappa * mk
                        r1_b = r1_b + t_ * z**ii
                        if (ii == pmax) last_b = last_b + t_ * z**ii
                     endif

                     if (l_j) then
                        sg = abs(m) + 2 + i + j
                        mk = 1.0d0 / (sg+k_) / (sg-k)
                        u_ = -w_ * mu/2 * mk
                        r1_b = r1_b + u_ * z**jj
                        if (jj == pmax) last_b = last_b + u_ * z**jj
                     endif

                     if (l_k) then
                        sg = abs(m) + 3 + i + j
                        mk = 1.0d0 / (sg+k_) / (sg-k)
                        v_ = w_ / (4*kappa) * mk
                        r1_b = r1_b + v_ * z**kk
                        if (kk == pmax) last_b = last_b + v_ * z**kk
                     endif

                  enddo
               enddo
            enddo
         enddo
         r1_b = r1_b / ga
         last_b = last_b / ga
      endif
      

c     ==== Infinite k part ====
      ubound = pmax - abs(m) - 1
      r1_a = 0.0d0
      last_a = 0.0d0
      do i = 0, ubound
         do j = 0, ubound
            do k = 0, ubound
               do k_ = 0, ubound
                  s1 = i + j + k + k_
                  ii = abs(m) + s1 + 1
                  jj = abs(m) + s1 + 2
                  kk = abs(m) + s1 + 3
                  l_i = (ii >= abs(m)+1 .and. ii <= pmax)
                  l_j = (jj >= abs(m)+2 .and. jj <= pmax)
                  l_k = (kk >= abs(m)+3 .and. kk <= pmax)
                  
                  if (l_i .or. l_j .or. l_k) then
                     w = (-1)**i/factorial(i)
     &                    *( rise_fact(a,j)/rise_fact(b,j)
     &                    /factorial(j) )
     &                    *( rise_fact(a,k)/rise_fact(b,k)
     &                    /factorial(k) )
     &                    *( rise_fact(a,k_)/rise_fact(b,k_)
     &                    /factorial(k_) )
                     call wfat_psif(k+a, ps1)
                     call wfat_psif(dble(k+1), ps2)
                     call wfat_psif(dble(k+1+abs(m)), ps3)
                     dk = ps1 - ps2 - ps3
                  endif
                  
                  if (l_i) then
                     sg = abs(m) + 1 + i + j
                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
                     t = w * beta1 * kappa * mk
                     r1_a = r1_a + t * z**ii
                     if (ii == pmax) last_a = last_a + t * z**ii
                  endif
                  
                  if (l_j) then
                     sg = abs(m) + 2 + i + j
                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
                     u = -w * mu/2 * mk
                     r1_a = r1_a + u * z**jj
                     if (jj == pmax) last_a = last_a + u * z**jj
                  endif

                  if (l_k) then
                     sg = abs(m) + 3 + i + j
                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
                     v = w / (4*kappa) * mk
                     r1_a = r1_a + v * z**kk
                     if (kk == pmax) last_a = last_a + v * z**kk
                  endif
               enddo
            enddo
         enddo
      enddo
      call wfat_gammaf(a-abs(m), gm)
      r1_a = r1_a * (-1)**(abs(m)+1) / factorial(abs(m)) / gm
      last_a = last_a * (-1)**(abs(m)+1) / factorial(abs(m)) / gm


c     ==== The final value of R_\nu^{(1)}(eta) ====
      r1 = (ga/kappa/factorial(abs(m)))**2 * z**(b/2) * exp(-z/2) *
     &     (r1_a - r1_b)
      last = (ga/kappa/factorial(abs(m)))**2 * z**(b/2) * exp(-z/2) *
     &       (last_a - last_b)

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_r1(et1, et2, r0_1, o0_1, n, m, erg, zq,
     &                   mu, npt, pt, wt, r0_2, o0_2, r1)
c     r0_1 = g/f_m * z1^(b/2) * exp(-z1/2) * M(a, b, z1)
c     r0_2 = g/f_m * z2^(b/2) * exp(-z2/2) * M(a, b, z2)
c     o0_1 = z1^(b/2) * exp(-z1/2) * U(a, b, z1)
c     o0_2 = z2^(b/2) * exp(-z2/2) * U(a, b, z2)
c     where
c     z1 = kappa * et1
c     z2 = kappa * et2
c     a = 1 + |m| + n - zq/kappa
c     b = |m| + 1
c     kappa = sqrt(2 * abs(erg))
c     g = the value of gamma function when its argument is equal to a.
c     f_m = the factorial of |m|.
c     
c     M(a,b,z) is calculated by wfat_cchg subroutine, while U(a,b,z) is calculated
c     by wfat_chgu subroutine.

      implicit none
#include "errquit.fh"            
      double precision et1         ! [input], 
      double precision et2         ! [input], 
      double precision r0_1        ! [input], 
      double precision o0_1        ! [input], 
      integer n                    ! [input], 
      integer m                    ! [input], 
      double precision erg         ! [input], 
      double precision zq          ! [input], 
      double precision mu          ! [input], 
      integer npt                  ! [input], 
      double precision pt(1:npt)   ! [input], 
      double precision wt(1:npt)   ! [input], 
      double precision r0_2        ! [output], 
      double precision o0_2        ! [output], 
      double precision r1          ! [output], 

      integer i, id
      double precision kappa, a, b, beta1, g, d, z2, f2, f, wt_, intd
      double precision r0, o0, et, z, chgm, chgu, chgm2, chgu2
      double complex zz, chgm_
      double precision factorial
      external factorial

      
      if (et1 > et2)
     &     call errquit('Argument et1 and et2 of wfat_r1 '//
     &     'must be such that et1 < et2.', 27, UERR)
      if (et2 < 1.0d-2)
     &     call errquit('Argument et2 of wfat_r1 must be '//
     &     'larger than or equal to 0.01.', 28, UERR)
      if (n < 0) call errquit('Argument n of '//
     &     'wfat_r1 must be zero or positive.', 29, UERR)

      
      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta1 = -mu/2/kappa * (2*n+1+abs(m))
     &        -1/4/kappa**2 * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))
      z2 = kappa * et2
      zz = cmplx(z2, 0.0d0, kind=kind(1.0d0))
      call wfat_cchg(a, b, zz, chgm_)
      chgm2 = real(chgm_)
      call wfat_chgu(a, b, z2, chgu2, id)
            
      f2 = z2**(b/2) * exp(-z2/2)
      r0_2 = d * f2 * chgm2
      o0_2 = f2 * chgu2
      

      r1 = 0.0d0
      do i = 1, npt
         if (i == 1) then
            r0 = r0_1
            o0 = o0_1
         else if (i == npt) then
            r0 = r0_2
            o0 = o0_2
         else
            et = 0.5d0*((et2-et1)*pt(i) + et2 + et1)
            z = kappa * et
            zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
            call wfat_cchg(a, b, zz, chgm_)
            chgm = real(chgm_)
            call wfat_chgu(a, b, z, chgu, id)
   
            f = z**(b/2) * exp(-z/2)
            r0 = d * f * chgm
            o0 = f * chgu
         endif

         wt_ = (et2 - et1)/2 * wt(i)
         intd = (r0*o0_2 - r0_2*o0) * (beta1/et - mu/2 + et/4) *
     &          r0
         r1 = r1 + wt_*intd
      enddo
      r1 = r1 / kappa
      
      end
c=======================================================
      
