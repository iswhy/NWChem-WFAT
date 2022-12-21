c=======================================================
c=======================================================
      subroutine wfat_firstord_r_smalleta(et2, n, m, erg, zq, mu, pmax,
     &     last, r1)

      implicit none

      double precision a, b
      

      if (pmax < abs(m)+1) call errquit()
      
      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      z = kappa * et2


      call wfat_gammaf(a, ga)
      r1_b = 0.0d0
      last_b = 0.0d0
      if (abs(m) > 0) then
         do i = 0, limit
            do j = 0, limit
               do k = abs(m), 1, -1
                  do k_ = 0, limit

                  enddo
               enddo
            enddo
         enddo
         r1_b = r1_b / ga
         last_b = last_b / ga
      endif
      

      limit = pmax-abs(m)-1
      r1_a = 0.0d0
      last_a = 0.0d0
      do i = 0, limit
         do j = 0, limit
            do k = 0, limit
               do k_ = 0, limit
                  s1 = i + j + k + k_
                  ii = abs(m) + s1 + 1
                  jj = abs(m) + s1 + 2
                  kk = abs(m) + s1 + 3
                  l_i = (ii >= abs(m)+1 .and. ii <= pmax)
                  l_j = (jj >= abs(m)+1 .and. jj <= pmax)
                  l_k = (kk >= abs(m)+1 .and. kk <= pmax)
                  
                  if (l_i .or. l_j .or. l_k) then
                     w = (-1)**i * 
     &                   rise_fact(a,j)*rise_fact(a,k)*rise_fact(a,k_) /
     &                   rise_fact(b,j)/rise_fact(b,k)/rise_fact(b,k_) /
     &                   factorial(i)/factorial(j)/factorial(k)/
     &                   factorial(k_)
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

      r1 = (ga/kappa/factorial(abs(m)))**2 * z**(b/2) * exp(-z/2) *
     &     (r1_a - r1_b)
      last = (ga/kappa/factorial(abs(m)))**2 * z**(b/2) * exp(-z/2) *
     &       (last_a - last_b)

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_firstord_r(et1, et2, r0_1, o0_1, n, m, erg, zq,
     &                           mu, npt, pt, wt, r0_2, o0_2, r1)
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
      integer npt
      double precision pt(1:npt)
      double precision wt(1:npt)

      double precision a, b

      if (et1 < 0.0d0)  resolvethis
      call errquit
      
      if (et2 < 1.0d-5)  resolvethis
      call errquit

      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta1 = 
      

      call wfat_gammaf(a, g)
      d = g/factorial(abs(m))
      z2 = kappa * et2
      f2 = z2**(b/2) * exp(-z2/2)
      r0_2 = d * f2 * chgm2
      o0_2 = f2 * chgu2
      

c      call wfat_lgl(et1, et2, npt-1, dbl_mb(k_pt), dbl_mb(k_wt))
      intl = 0.0d0
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
         intl = intl + wt_*intd
      enddo
      

      

      end
c=======================================================
      
