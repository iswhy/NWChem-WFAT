c=======================================================
c=======================================================
      subroutine wfat_r1(neta, eta, n, m, erg, zq, bt, gm, chmmt, r1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer neta
      double precision eta(neta)    ! [input], increasingly ordered eta.
      integer n
      integer m
      double precision erg
      double precision zq
      double precision bt           ! [input], beta orientation angle in degree.
      double precision gm           ! [input], gamma orientation angle in degree.
      double precision chmmt(3)
      double precision r1(neta)

      logical first_notsmalleta
      integer i, j, pmax, l_cr, k_cr
      double precision b, protmat(3,3), mu(3), kappa, z, rr
      double precision zbound
      parameter (zbound=1.0d0, pmax=40)
      

c     ==== Get the passive rotation matrix ===
      call wfat_passrot(-bt, -gm, protmat)
      mu = matmul(transpose(protmat), -chmmt)
      

      kappa = sqrt(2 * abs(erg))
cforlater      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
cforlater      beta0 = zq - kappa*(n + (abs(m)+1)/2.0d0)
cforlater      beta1 = -mu(3)/2/kappa * (2*n+1+abs(m))
cforlater     &        -1.0d0/4/kappa**2
cforlater     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
cforlater      call wfat_gammaf(a, g)
cforlater      d = g / factorial(abs(m))


cforlater      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_pt', l_pt, k_pt))
cforlater     &     call errquit('WFAT: lgl_pt allocation failed', 84, MA_ERR)
cforlater      if (.not. MA_alloc_get(MT_DBL, nlgl, 'lgl_wt', l_wt, k_wt))
cforlater     &     call errquit('WFAT: lgl_wt allocation failed', 85, MA_ERR)
cforlater      call wfat_lgl(-1.0d0, 1.0d0, nlgl-1, dbl_mb(k_pt), dbl_mb(k_wt))

      
c     ==== Coefficients of A1-A2 ====
      if (.not. MA_alloc_get(MT_DBL, pmax, 'r1_coef', l_cr, k_cr))
     &     call errquit('WFAT: r1_coef allocation failed', 62, MA_ERR)
      call wfat_r1_coeff(n, m, erg, zq, mu(3), pmax, dbl_mb(k_cr))

      
      first_notsmalleta = .true.
c     ==== Loop over the increasingly-ordered radial points ====
      do i = 1, neta
         z = kappa * eta(i)
         if (z < zbound) then
            rr = 0.0d0
            do j = 1, pmax
               rr = rr + dbl_mb(k_cr+j-1) * z**j
            enddo
            r1(i) = z**(b/2) * exp(-z/2) * rr

         else
            call errquit('WFAT: For the moment, z cannot exceed 1.0 '//
     &           'in wfat_r1.', 49, UERR)
         endif
      enddo

      
      if (.not. ma_free_heap(l_cr))
     &     call errquit('WFAT: Cannot deallocate r1_coef', 49, MA_ERR)
      
      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_r1_coeff(n, m, erg, zq, mu, pmax, cr1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer n                      ! [input]
      integer m                      ! [input]
      double precision erg           ! [input]
      double precision zq            ! [input]
      double precision mu            ! [input]
      integer pmax                   ! [input]
      double precision cr1(1:pmax)   ! [output]

      logical l_i, l_j, l_k
      integer l_a, k_a, l_b, k_b
      integer i, j, k, k_, ii, jj, kk, k3, s1, s2, ubound, sg
      double precision kappa, a, b, z, beta1, ga, gm, mk, w_, t_, u_, v_
      double precision w, t, u, v, dk, ps1, ps2, ps3
      double precision factorial, rise_fact
      external factorial, rise_fact


      if (pmax < 1) call errquit('Argument pmax of '//
     &     'wfat_r1_smalleta must be larger than 1.', 28, UERR)
      if (n < 0) call errquit('Argument n of '//
     &     'wfat_r1_smalleta must be zero or positive.', 28, UERR)

      cr1 = 0.0d0
      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta1 = -mu/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      call wfat_gammaf(a, ga)


      if (.not. MA_alloc_get(MT_DBL, pmax, 'part_a', l_a, k_a))
     &     call errquit('WFAT: part_a allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, pmax, 'part_b', l_b, k_b))
     &     call errquit('WFAT: part_b allocation failed', 85, MA_ERR)
      

c     ==== Finite k part ====
      dbl_mb(k_b:k_b+pmax-1) = 0.0d0
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
                        k3 = k_b + ii - 1
                        dbl_mb(k3) = dbl_mb(k3) + t_
                     endif

                     if (l_j) then
                        sg = abs(m) + 2 + i + j
                        mk = 1.0d0 / (sg+k_) / (sg-k)
                        u_ = -w_ * mu/2 * mk
                        k3 = k_b + jj - 1
                        dbl_mb(k3) = dbl_mb(k3) + u_
                     endif

                     if (l_k) then
                        sg = abs(m) + 3 + i + j
                        mk = 1.0d0 / (sg+k_) / (sg-k)
                        v_ = w_ / (4*kappa) * mk
                        k3 = k_b + kk - 1
                        dbl_mb(k3) = dbl_mb(k3) + v_
                     endif

                  enddo
               enddo
            enddo
         enddo
         dbl_mb(k_b:k_b+pmax-1) = dbl_mb(k_b:k_b+pmax-1) / ga
      endif
      

c     ==== Infinite k part ====
      ubound = pmax - abs(m) - 1
      dbl_mb(k_a:k_a+pmax-1) = 0.0d0
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
                     k3 = k_a + ii - 1
                     dbl_mb(k3) = dbl_mb(k3) + t
                  endif
                  
                  if (l_j) then
                     sg = abs(m) + 2 + i + j
                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
                     u = -w * mu/2 * mk
                     k3 = k_a + jj - 1
                     dbl_mb(k3) = dbl_mb(k3) + u
                  endif

                  if (l_k) then
                     sg = abs(m) + 3 + i + j
                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
                     v = w / (4*kappa) * mk
                     k3 = k_a + kk - 1
                     dbl_mb(k3) = dbl_mb(k3) + v
                  endif
                  
               enddo
            enddo
         enddo
      enddo
      call wfat_gammaf(a-abs(m), gm)
      dbl_mb(k_a:k_a+pmax-1) = dbl_mb(k_a:k_a+pmax-1)
     &     * (-1)**(abs(m)+1) / factorial(abs(m)) / gm


c     ==== The final value of R_\nu^{(1)}(eta) ====
      cr1 = (ga/kappa/factorial(abs(m)))**2 * 
     &     (dbl_mb(k_a:k_a+pmax-1) - dbl_mb(k_b:k_b+pmax-1))


      write(luout, '("R1 coefficients in wfat_r1")')
      do j = 1, pmax
         write(luout, '(i3, 3es15.6)') j, cr1(j),
     &        dbl_mb(k_a+j-1), dbl_mb(k_b+j-1)
      enddo

      
c     z**(b/2) * exp(-z/2) *

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_a))
     &     call errquit('WFAT: Cannot deallocate part_a', 49, MA_ERR)
      if (.not. ma_free_heap(l_b))
     &     call errquit('WFAT: Cannot deallocate part_b', 49, MA_ERR)
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_o0_i_coeff(n, m, erg, zq, mu, pmax, o0i, o0il)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n                      ! [input]
      integer m                      ! [input]
      double precision erg           ! [input]
      double precision zq            ! [input]
      double precision mu            ! [input]
      integer pmax                   ! [input]
      double precision o0i(1:pmax)   ! [output]
      double precision o0il(1:pmax)  ! [output]

      logical l_i, l_j, l_k
      integer ubound, i, j, k, ii, jj, kk, s1, s2
      double precision kappa, a, b, beta1, v1, v2, ga, gm, w_, t_, u_
      double precision v_, w, ps1, ps2, ps3, dk, mk, t, u, v
      double precision factorial, rise_fact
      external factorial, rise_fact


      if (pmax < 1)
     &     call errquit('WFAT: Argument pmax of wfat_o0_i_coeff '//
     &     'cannot be smaller than 1.', 48, UERR)
      if (n < 0) call errquit('Argument n of '//
     &     'wfat_o0_i_coeff must be zero or positive.', 28, UERR)


      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta1 = -mu/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)

      v1 = 1.0d0 / kappa**2 / factorial(abs(m))
      call wfat_gammaf(a, ga)
      call wfat_gammaf(a-abs(m), gm)
      v2 = (-1)**(abs(m)+1) * ga / (kappa*factorial(m))**2 / gm

      
c     ==== Finite k part ====
      o0i = 0.0d0
      if (abs(m) > 0) then
         ubound = pmax - 1
         do i = 0, ubound
            do j = 0, ubound
               do k = abs(m), 1, -1
                  s2 = i + j - k
                  ii = abs(m) + s2 + 1
                  jj = abs(m) + s2 + 2
                  kk = abs(m) + s2 + 3
                  l_i = (ii >= 1 .and. ii <= pmax)
                  l_j = (jj >= 2 .and. jj <= pmax)
                  l_k = (kk >= 3 .and. kk <= pmax)

                  if (l_i .or. l_j .or. l_k) then
                     w_ = v1 * (-1)**i/factorial(i)
     &                    *( rise_fact(a,j)/factorial(j)
     &                    /rise_fact(b,j) )
     &                    *( factorial(k-1)*rise_fact(1-a+k,abs(m)-k)
     &                    /factorial(abs(m)-k) )
                  endif

                  if (l_i) then
                     t_ = w_ * beta1 * kappa / ii
                     o0i(ii) = o0i(ii) + t_
                  endif

                  if (l_j) then
                     u_ = -w_ * mu/2 / jj
                     o0i(jj) = o0i(jj) + u_
                  endif

                  if (l_k) then
                     v_ = w_ / (4*kappa) / kk
                     o0i(kk) = o0i(kk) + v_
                  endif
               enddo
            enddo
         enddo
      endif      


c     ==== Infinite k part ====
      ubound = pmax - abs(m) - 1
      do i = 0, ubound
         do j = 0, ubound
            do k = 0, ubound
               s1 = i + j + k
               ii = abs(m) + s1 + 1
               jj = abs(m) + s1 + 2
               kk = abs(m) + s1 + 3
               l_i = (ii >= abs(m)+1 .and. ii <= pmax)
               l_j = (jj >= abs(m)+2 .and. jj <= pmax)
               l_k = (kk >= abs(m)+3 .and. kk <= pmax)

               if (l_i .or. l_j .or. l_k) then
                  w = v2 * (-1)**i/factorial(i)
     &                *( rise_fact(a,j)/rise_fact(b,j)
     &                /factorial(j) )
     &                *( rise_fact(a,k)/rise_fact(b,k)
     &                /factorial(k) )
                  call wfat_psif(k+a, ps1)
                  call wfat_psif(dble(k+1), ps2)
                  call wfat_psif(dble(k+1+abs(m)), ps3)
                  dk = ps1 - ps2 - ps3
               endif
               
               if (l_i) then
                  mk = dk/ii - 1.0d0/ii**2
                  t = w * beta1 * kappa * mk
                  o0i(ii) = o0i(ii) + t

                  t = w * beta1 * kappa / ii
                  o0il(ii) = o0il(ii) + t
               endif
               
               if (l_j) then
                  mk = dk/jj - 1.0d0/jj**2
                  u = -w * mu/2 * mk
                  o0i(jj) = o0i(jj) + u

                  u = -w * mu/2 * / jj
                  o0il(jj) = o0il(jj) + u
               endif

               if (l_k) then
                  mk = dk/kk - 1.0d0/kk**2
                  v = w / (4*kappa) * mk
                  o0i(kk) = o0i(kk) + v

                  v = w / (4*kappa) / kk
                  o0il(kk) = o0il(kk) + v
               endif
            enddo
         enddo
      enddo
            

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_r0_i_coeff(n, m, erg, zq, mu, pmax, r0i)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer n                      ! [input]
      integer m                      ! [input]
      double precision erg           ! [input]
      double precision zq            ! [input]
      double precision mu            ! [input]
      integer pmax                   ! [input]
      double precision r0i(1:pmax)

      logical l_i, l_j, l_k
      integer i, j, k_, s1, ii, jj, kk, ubound
      double precision kappa, a, beta1, ga, v, w, t, u, v
      double precision factorial, rise_fact
      external factorial, rise_fact


      if (n < 0) call errquit('Argument n of '//
     &     'wfat_r0_i_coeff must be zero or positive.', 28, UERR)
      

      kappa = sqrt(2*abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      beta1 = -mu/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      call wfat_gammaf(a, ga)
      v = ( ga / kappa / factorial(abs(m)) )**2


c     ==== Infinite k part ====
      r0i = 0.0d0
      ubound = pmax - abs(m) - 1
      do i = 0, ubound
         do j = 0, ubound
            do k_ = 0, ubound
               s1 = i + j + k_
               ii = abs(m) + s1 + 1
               jj = abs(m) + s1 + 2
               kk = abs(m) + s1 + 3
               l_i = (ii >= abs(m)+1 .and. ii <= pmax)
               l_j = (jj >= abs(m)+2 .and. jj <= pmax)
               l_k = (kk >= abs(m)+3 .and. kk <= pmax)

               if (l_i .or. l_j .or. l_k) then
                  w = v * (-1)**i/factorial(i)
     &                *( rise_fact(a,j)/rise_fact(b,j)
     &                /factorial(j) )
     &                *( rise_fact(a,k_)/rise_fact(b,k_)
     &                /factorial(k_) )
               endif
               
               if (l_i) then
                  t = w * beta1 * kappa / ii
                  r0i(ii) = r0i(ii) + t
               endif
               
               if (l_j) then
                  u = -w * mu/2 / jj
                  r0i(jj) = r0i(jj) + u
               endif

               if (l_k) then
                  v = w / (4*kappa) / kk
                  r0i(kk) = r0i(kk) + v
               endif
            enddo
         enddo
      enddo      

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_r0o0_i(eta1, eta2, n, m, erg, zq, mu, nlgl,
     &           pt, wt, rr, oo)

c     RR = \frac{1}{\kappa} \int_{\eta_1}^{\eta2} d\eta' (R_\nu^{(0)}(\eta'))^2
c          (\beta_\nu^{(1)}/\eta' - \mu_z/2 + \eta'/4)
c     OO = \frac{1}{\kappa} \int_{\eta_1}^{\eta2} d\eta' O_\nu^{(0)}(\eta')
c          (\beta_\nu^{(1)}/\eta' - \mu_z/2 + \eta'/4) R_\nu^{(0)}(\eta')
      
      implicit none
#include "stdio.fh"
      double precision eta1
      double precision eta2
      integer n                      ! [input]
      integer m                      ! [input]
      double precision erg           ! [input]
      double precision zq            ! [input]
      double precision mu            ! [input]
      integer nlgl
      double precision pt(nlgl)
      double precision wt(nlgl)
      double precision rr
      double precision oo

      integer nfem, ngr, id, i, j, md
      double precision pt_4(4), wt_4(4), pt_6(6), wt_6(6)
      double precision kappa, a, b, beta1, g, d, de, dfem, eta1_, eta2_
      double precision eta, wt_, ee, z, r0, o0, ind, min_dfem
      double complex zz, chgm_
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


      kappa = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kappa
      b = dble(abs(m) + 1)
      beta1 = -mu/2/kappa * (2*n+1+abs(m))
     &        -1.0d0/4/kappa**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))
      

c     ==== Determine the no. of finite elements and grids ====
c     de = the total interval of r to be integrated
c     dfem = the finite element interval
      de = eta2 - eta1
      if (de < 1.0d-12) then
         nfem = 1
         dfem = de
         ngr = 2
         id = 0
      else if (de >= 1.0d-12 .and. de < 1.0d-3) then
         nfem = 1
         dfem = de
         ngr = 4
         id = 1
      else if (de >= 1.0d-3 .and. de < 1.0d-1) then
         nfem = 1
         dfem = de
         ngr = 6
         id = 2
      else
         nfem = int(max(de,min_dfem+1.0d-13) / min_dfem)
         dfem = de / nfem
         ngr = nlgl
         id = 3
      endif


      rr = 0.0d0
      oo = 0.0d0
c     ==== Loop over finite elements ====
c     eta1_ = lower bound of a FEM
c     eta2_ = upper bound of a FEM
c     eta = shifted quadrature grid
c     wt_ = shifted quadrature weight
      do i = 1, nfem
         eta1_ = eta1 + (i-1)*dfem
         eta2_ = eta1_ + dfem

c        ==== Loop over quadrature grids ====
         do j = 1, ngr
c           ==== Determine the shifted grids and weights ====
            if (id == 0) then
               eta = eta1_ + (j-1)*dfem
               wt_ = (eta2_ - eta1_)/2
            else if (id == 1) then
               eta = 0.5d0*((eta2_-eta1_)*pt_4(j) + eta2_ + eta1_)
               wt_ = (eta2_ - eta1_)/2 * wt_4(j)
            else if (id == 2) then
               eta = 0.5d0*((eta2_-eta1_)*pt_6(j) + eta2_ + eta1_)
               wt_ = (eta2_ - eta1_)/2 * wt_6(j)
            else if (id == 3) then
               eta = 0.5d0*((eta2_-eta1_)*pt(j) + eta2_ + eta1_)
               wt_ = (eta2_ - eta1_)/2 * wt(j)
            endif
            if (eta < 1.0d-13) cycle           ! When eta is zero, the integrands are equal to zero.
            ee = (beta1/eta - mu/2 + eta/4)
            

c           ==== Calculate the necessary Hypergeometric functions ====
            z = kappa * eta
            zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
            call wfat_cchg(a, b, zz, chgm_)
            r0 = d * z**(b/2) * exp(-z/2) * real(chgm_)
            call wfat_chgu(a, b, z, chgu, md)
            o0 = z**(b/2) * exp(-z/2) * chgu
            

c           ==== The R^{(0)} part ====
            ind = r0 * r0 * ee
            rr = rr + wt_*ind

c           ==== The O^{(0)} part ====
            ind = o0 * r0 * ee
            oo = oo + wt_*ind
         enddo
      enddo

      rr = rr / kappa
      oo = oo / kappa
      

      end
c=======================================================
      
      
      
cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_r1_BAD(et1, et2, r0_1, o0_1, n, m, erg, zq,
cOLD     &                   mu, npt, pt, wt, r0_2, o0_2, r1)
cOLDc     r0_1 = g/f_m * z1^(b/2) * exp(-z1/2) * M(a, b, z1)
cOLDc     r0_2 = g/f_m * z2^(b/2) * exp(-z2/2) * M(a, b, z2)
cOLDc     o0_1 = z1^(b/2) * exp(-z1/2) * U(a, b, z1)
cOLDc     o0_2 = z2^(b/2) * exp(-z2/2) * U(a, b, z2)
cOLDc     where
cOLDc     z1 = kappa * et1
cOLDc     z2 = kappa * et2
cOLDc     a = 1 + |m| + n - zq/kappa
cOLDc     b = |m| + 1
cOLDc     kappa = sqrt(2 * abs(erg))
cOLDc     g = the value of gamma function when its argument is equal to a.
cOLDc     f_m = the factorial of |m|.
cOLDc     
cOLDc     M(a,b,z) is calculated by wfat_cchg subroutine, while U(a,b,z) is calculated
cOLDc     by wfat_chgu subroutine.
cOLD
cOLD      implicit none
cOLD#include "errquit.fh"
cOLD#include "stdio.fh"
cOLD      double precision et1         ! [input], 
cOLD      double precision et2         ! [input], 
cOLD      double precision r0_1        ! [input], 
cOLD      double precision o0_1        ! [input], 
cOLD      integer n                    ! [input], 
cOLD      integer m                    ! [input], 
cOLD      double precision erg         ! [input], 
cOLD      double precision zq          ! [input], 
cOLD      double precision mu          ! [input], 
cOLD      integer npt                  ! [input], 
cOLD      double precision pt(1:npt)   ! [input], 
cOLD      double precision wt(1:npt)   ! [input], 
cOLD      double precision r0_2        ! [output], 
cOLD      double precision o0_2        ! [output], 
cOLD      double precision r1          ! [output], 
cOLD
cOLD      integer i, id
cOLD      double precision kappa, a, b, beta1, g, d, z2, f2, f, wt_, intd
cOLD      double precision r0, o0, et, z, chgm, chgu, chgm2, chgu2
cOLD      double complex zz, chgm_
cOLD      double precision factorial
cOLD      external factorial
cOLD
cOLD      
cOLD      if (et1 > et2)
cOLD     &     call errquit('Argument et1 and et2 of wfat_r1 '//
cOLD     &     'must be such that et1 < et2.', 27, UERR)
cOLD      if (et2 < 1.0d-2)
cOLD     &     call errquit('Argument et2 of wfat_r1 must be '//
cOLD     &     'larger than or equal to 0.01.', 28, UERR)
cOLD      if (n < 0) call errquit('Argument n of '//
cOLD     &     'wfat_r1 must be zero or positive.', 29, UERR)
cOLD
cOLD      
cOLD      kappa = sqrt(2*abs(erg))
cOLD      a = 1 + abs(m) + n - zq/kappa
cOLD      b = dble(abs(m) + 1)
cOLD      beta1 = -mu/2/kappa * (2*n+1+abs(m))
cOLD     &        -1.0d0/4/kappa**2
cOLD     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
cOLD      
cOLD      call wfat_gammaf(a, g)
cOLD      d = g / factorial(abs(m))
cOLD      z2 = kappa * et2
cOLD      zz = cmplx(z2, 0.0d0, kind=kind(1.0d0))
cOLD      call wfat_cchg(a, b, zz, chgm_)
cOLD      chgm2 = real(chgm_)
cOLD      call wfat_chgu(a, b, z2, chgu2, id)
cOLD            
cOLD      f2 = z2**(b/2) * exp(-z2/2)
cOLD      r0_2 = d * f2 * chgm2
cOLD      o0_2 = f2 * chgu2
cOLD      
cOLD
cOLD      if (et2 < 5.0d0) then
cOLD         write(luout, '(2es15.6)') et1, et2
cOLD         write(luout, '(4x, a5, 5a15)')
cOLD     &        'No.', 'eta', 'R^(0)', 'O^(0)', 'weight', 'integrand'
cOLD      endif
cOLD      
cOLD      r1 = 0.0d0
cOLD      do i = 1, npt
cOLD         if (i == 1) then
cOLD            et = et1
cOLD            r0 = r0_1
cOLD            o0 = o0_1
cOLD         else if (i == npt) then
cOLD            et = et2
cOLD            r0 = r0_2
cOLD            o0 = o0_2
cOLD         else
cOLD            et = 0.5d0*((et2-et1)*pt(i) + et2 + et1)
cOLD            z = kappa * et
cOLD            zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
cOLD            call wfat_cchg(a, b, zz, chgm_)
cOLD            chgm = real(chgm_)
cOLD            call wfat_chgu(a, b, z, chgu, id)
cOLD   
cOLD            f = z**(b/2) * exp(-z/2)
cOLD            r0 = d * f * chgm
cOLD            o0 = f * chgu
cOLD         endif
cOLD
cOLD         wt_ = (et2 - et1)/2 * wt(i)
cOLD         intd = (r0*o0_2 - r0_2*o0) * (beta1/et - mu/2 + et/4) *
cOLD     &          r0
cOLD         r1 = r1 + wt_*intd
cOLD
cOLD         if (et2 < 5.0d0) then
cOLD            write(luout, '(4x, i5, 5es15.6)')
cOLD     &           i, et, r0, o0, wt_, intd
cOLD         endif
cOLD         
cOLD      enddo
cOLD      r1 = r1 / kappa
cOLD      
cOLD      end
cOLDc=======================================================
cOLD
cOLD
cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_r1_trapz(et1, et2, n, m, erg, zq,
cOLD     &     mu, npt, r1)
cOLD
cOLD      implicit none
cOLD#include "errquit.fh"
cOLD#include "stdio.fh"
cOLD      double precision et1         ! [input], 
cOLD      double precision et2         ! [input], 
cOLD      integer n                    ! [input], 
cOLD      integer m                    ! [input], 
cOLD      double precision erg         ! [input], 
cOLD      double precision zq          ! [input], 
cOLD      double precision mu          ! [input], 
cOLD      integer npt                  ! [input], 
cOLD      double precision r1          ! [output], 
cOLD
cOLD      integer i, id
cOLD      double precision kappa, a, b, beta1, g, d, z2, f2, f, intd
cOLD      double precision de, r0_2, o0_2, r0, o0, et, z, chgm, chgu
cOLD      double precision chgm2, chgu2, wt
cOLD      double complex zz, chgm_
cOLD      double precision factorial
cOLD      external factorial
cOLD
cOLD      
cOLD      if (et1 > et2)
cOLD     &     call errquit('Argument et1 and et2 of wfat_r1 '//
cOLD     &     'must be such that et1 < et2.', 27, UERR)
cOLD      if (et2 < 1.0d-2)
cOLD     &     call errquit('Argument et2 of wfat_r1 must be '//
cOLD     &     'larger than or equal to 0.01.', 28, UERR)
cOLD      if (n < 0) call errquit('Argument n of '//
cOLD     &     'wfat_r1 must be zero or positive.', 29, UERR)
cOLD
cOLD      
cOLD      kappa = sqrt(2*abs(erg))
cOLD      a = 1 + abs(m) + n - zq/kappa
cOLD      b = dble(abs(m) + 1)
cOLD      beta1 = -mu/2/kappa * (2*n+1+abs(m))
cOLD     &        -1.0d0/4/kappa**2
cOLD     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
cOLD      
cOLD      call wfat_gammaf(a, g)
cOLD      d = g / factorial(abs(m))
cOLD      z2 = kappa * et2
cOLD      zz = cmplx(z2, 0.0d0, kind=kind(1.0d0))
cOLD      call wfat_cchg(a, b, zz, chgm_)
cOLD      chgm2 = real(chgm_)
cOLD      call wfat_chgu(a, b, z2, chgu2, id)
cOLD            
cOLD      f2 = z2**(b/2) * exp(-z2/2)
cOLD      r0_2 = d * f2 * chgm2
cOLD      o0_2 = f2 * chgu2
cOLD
cOLD
cOLD      de = (et2-et1)/(npt-1)
cOLD
cOLD      r1 = 0.0d0
cOLD      do i = 1, npt
cOLD         et = et1 + (i-1)*de
cOLD         z = kappa * et
cOLD         zz = cmplx(z, 0.0d0, kind=kind(1.0d0))
cOLD         call wfat_cchg(a, b, zz, chgm_)
cOLD         chgm = real(chgm_)
cOLD         call wfat_chgu(a, b, z, chgu, id)
cOLD         f = z**(b/2) * exp(-z/2)
cOLD         r0 = d * f * chgm
cOLD         o0 = f * chgu
cOLD
cOLD         intd = (r0*o0_2 - r0_2*o0) * (beta1/et - mu/2 + et/4) * r0
cOLD         wt = de
cOLD         if (i==1 .or. i==npt) wt = de/2
cOLD         
cOLD         r1 = r1 + wt * intd
cOLD      enddo
cOLD
cOLD      r1 = r1 / kappa
cOLD      
cOLD      end
cOLDc=======================================================
      

cOLDc=======================================================
cOLDc=======================================================
cOLD      subroutine wfat_r1_smalleta(et, n, m, erg, zq, mu, pmax,
cOLD     &           last, r1)
cOLD
cOLD      implicit none
cOLD#include "errquit.fh"      
cOLD      double precision et     ! [input]
cOLD      integer n               ! [input]
cOLD      integer m               ! [input]
cOLD      double precision erg    ! [input]
cOLD      double precision zq     ! [input]
cOLD      double precision mu     ! [input]
cOLD      integer pmax            ! [input]
cOLD      double precision last   ! [output]
cOLD      double precision r1     ! [output]
cOLD
cOLD      logical l_i, l_j, l_k
cOLD      integer i, j, k, k_, ii, jj, kk, s1, s2, ubound, sg
cOLD      double precision kappa, a, b, z, beta1, ga, gm, mk, w_, t_, u_, v_
cOLD      double precision w, t, u, v, dk, ps1, ps2, ps3
cOLD      double precision r1_a, r1_b, last_a, last_b
cOLD      double precision factorial, rise_fact
cOLD      external factorial, rise_fact
cOLD      
cOLD
cOLD      if (pmax < 1) call errquit('Argument pmax of '//
cOLD     &     'wfat_r1_smalleta must be larger than 1.', 28, UERR)
cOLD      if (n < 0) call errquit('Argument n of '//
cOLD     &     'wfat_r1_smalleta must be zero or positive.', 28, UERR)
cOLD      
cOLD      kappa = sqrt(2*abs(erg))
cOLD      a = 1 + abs(m) + n - zq/kappa
cOLD      b = dble(abs(m) + 1)
cOLD      z = kappa * et
cOLD      beta1 = -mu/2/kappa * (2*n+1+abs(m))
cOLD     &        -1.0d0/4/kappa**2
cOLD     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
cOLD      call wfat_gammaf(a, ga)
cOLD
cOLD
cOLDc     ==== Finite k part ====
cOLD      r1_b = 0.0d0
cOLD      last_b = 0.0d0
cOLD      if (abs(m) > 0) then
cOLD         ubound = pmax - 1
cOLD         do i = 0, ubound
cOLD            do j = 0, ubound
cOLD               do k = abs(m), 1, -1
cOLD                  do k_ = 0, ubound
cOLD                     s2 = i + j - k + k_
cOLD                     ii = abs(m) + s2 + 1
cOLD                     jj = abs(m) + s2 + 2
cOLD                     kk = abs(m) + s2 + 3
cOLD                     l_i = (ii >= 1 .and. ii <= pmax)
cOLD                     l_j = (jj >= 2 .and. jj <= pmax)
cOLD                     l_k = (kk >= 3 .and. kk <= pmax)
cOLD
cOLD                     if (l_i .or. l_j .or. l_k) then
cOLD                        w_ = (-1)**i/factorial(i)
cOLD     &                       *( rise_fact(a,j)/factorial(j)
cOLD     &                       /rise_fact(b,j) )
cOLD     &                       *( rise_fact(a,k_)/factorial(k_)
cOLD     &                       /rise_fact(b,k_) )
cOLD     &                       *( factorial(k-1)*rise_fact(1-a+k,abs(m)-k)
cOLD     &                       /factorial(abs(m)-k) )
cOLD     &                       *(k+k_) 
cOLD                     endif
cOLD
cOLD                     if (l_i) then
cOLD                        sg = abs(m) + 1 + i + j
cOLD                        mk = 1.0d0 / (sg+k_) / (sg-k)
cOLD                        t_ = w_ * beta1 * kappa * mk
cOLD                        r1_b = r1_b + t_ * z**ii
cOLD                        if (ii == pmax) last_b = last_b + t_ * z**ii
cOLD                     endif
cOLD
cOLD                     if (l_j) then
cOLD                        sg = abs(m) + 2 + i + j
cOLD                        mk = 1.0d0 / (sg+k_) / (sg-k)
cOLD                        u_ = -w_ * mu/2 * mk
cOLD                        r1_b = r1_b + u_ * z**jj
cOLD                        if (jj == pmax) last_b = last_b + u_ * z**jj
cOLD                     endif
cOLD
cOLD                     if (l_k) then
cOLD                        sg = abs(m) + 3 + i + j
cOLD                        mk = 1.0d0 / (sg+k_) / (sg-k)
cOLD                        v_ = w_ / (4*kappa) * mk
cOLD                        r1_b = r1_b + v_ * z**kk
cOLD                        if (kk == pmax) last_b = last_b + v_ * z**kk
cOLD                     endif
cOLD
cOLD                  enddo
cOLD               enddo
cOLD            enddo
cOLD         enddo
cOLD         r1_b = r1_b / ga
cOLD         last_b = last_b / ga
cOLD      endif
cOLD      
cOLD
cOLDc     ==== Infinite k part ====
cOLD      ubound = pmax - abs(m) - 1
cOLD      r1_a = 0.0d0
cOLD      last_a = 0.0d0
cOLD      do i = 0, ubound
cOLD         do j = 0, ubound
cOLD            do k = 0, ubound
cOLD               do k_ = 0, ubound
cOLD                  s1 = i + j + k + k_
cOLD                  ii = abs(m) + s1 + 1
cOLD                  jj = abs(m) + s1 + 2
cOLD                  kk = abs(m) + s1 + 3
cOLD                  l_i = (ii >= abs(m)+1 .and. ii <= pmax)
cOLD                  l_j = (jj >= abs(m)+2 .and. jj <= pmax)
cOLD                  l_k = (kk >= abs(m)+3 .and. kk <= pmax)
cOLD                  
cOLD                  if (l_i .or. l_j .or. l_k) then
cOLD                     w = (-1)**i/factorial(i)
cOLD     &                    *( rise_fact(a,j)/rise_fact(b,j)
cOLD     &                    /factorial(j) )
cOLD     &                    *( rise_fact(a,k)/rise_fact(b,k)
cOLD     &                    /factorial(k) )
cOLD     &                    *( rise_fact(a,k_)/rise_fact(b,k_)
cOLD     &                    /factorial(k_) )
cOLD                     call wfat_psif(k+a, ps1)
cOLD                     call wfat_psif(dble(k+1), ps2)
cOLD                     call wfat_psif(dble(k+1+abs(m)), ps3)
cOLD                     dk = ps1 - ps2 - ps3
cOLD                  endif
cOLD                  
cOLD                  if (l_i) then
cOLD                     sg = abs(m) + 1 + i + j
cOLD                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
cOLD                     t = w * beta1 * kappa * mk
cOLD                     r1_a = r1_a + t * z**ii
cOLD                     if (ii == pmax) last_a = last_a + t * z**ii
cOLD                  endif
cOLD                  
cOLD                  if (l_j) then
cOLD                     sg = abs(m) + 2 + i + j
cOLD                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
cOLD                     u = -w * mu/2 * mk
cOLD                     r1_a = r1_a + u * z**jj
cOLD                     if (jj == pmax) last_a = last_a + u * z**jj
cOLD                  endif
cOLD
cOLD                  if (l_k) then
cOLD                     sg = abs(m) + 3 + i + j
cOLD                     mk = dble(k-k_)*dk/(sg+k_)/(sg+k) + 1.0d0/(sg+k)**2
cOLD                     v = w / (4*kappa) * mk
cOLD                     r1_a = r1_a + v * z**kk
cOLD                     if (kk == pmax) last_a = last_a + v * z**kk
cOLD                  endif
cOLD               enddo
cOLD            enddo
cOLD         enddo
cOLD      enddo
cOLD      call wfat_gammaf(a-abs(m), gm)
cOLD      r1_a = r1_a * (-1)**(abs(m)+1) / factorial(abs(m)) / gm
cOLD      last_a = last_a * (-1)**(abs(m)+1) / factorial(abs(m)) / gm
cOLD
cOLD
cOLDc     ==== The final value of R_\nu^{(1)}(eta) ====
cOLD      r1 = (ga/kappa/factorial(abs(m)))**2 * z**(b/2) * exp(-z/2) *
cOLD     &     (r1_a - r1_b)
cOLD      last = (ga/kappa/factorial(abs(m)))**2 * z**(b/2) * exp(-z/2) *
cOLD     &       (last_a - last_b)
cOLD
cOLD      end
cOLDc=======================================================      
