c=======================================================
c=======================================================
      subroutine wfat_calc_o(n, m, erg, zq, o1, o2)

      implicit none
      integer n
      integer m
      double precision erg
      double precision zq
      double precision o1
      double precision o2

      double precision kp, gmm, bt0, o1i, o2i

      
      kp = sqrt(2 * abs(erg))
      gmm = (1.0d0 - m**2) / 4.0d0
      bt0 = zq - kp*(n + (abs(m)+1)/2.0d0)

      
      call wfat_calc_o_intl(, o1i, o2i)

      
      o1 = 1/(4*kp**3)*(1 + gmm**2/2 - gmm*(4*bt0/kp-(bt0/kp)**2)
     &     + 4*bt0/kp + (7.0d0/2)*(bt0/kp)**2 - 4*(bt0/kp)**3 
     &     + (bt0/kp)**4/2) + o1i
      o2 = -1.0d0/(2*kp**2) + o2i

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intl()

      implicit none



      eb = 1.0d0/kappa        ! so that zb = 1.0d0 inside wfat_calc_o_intl_part* below.

      call wfat_calc_o_intl_part1(, eb, , , )
      call wfat_calc_o_intl_part2(, eb, , , )

      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intl_part1(n, m, erg, zq, eb, o_intl)

      implicit none
      integer n
      integer m
      double precision erg
      double precision zq
      double precision eb
      double precision o_intl(2)

      integer i, j, pmax
      double precision kp, zb, gmm, bt1_1, bt1_2
      double precision part_pow(2), part_dir(2)
      parameter (pmax=40)
      double precision o0i(pmax,2), o0il(pmax,2)


      kp = sqrt(2 * abs(erg))
      zb = kp * eb
      bt1_1 = -1.0d0/4/kp**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      bt1_2 = -1.0d0/2/kp * (2*n+1+abs(m))
      gmm = (1.0d0 - m**2) / 4.0d0
      

c     ==== Do the power expansion part ====
      call wfat_o0_i_coeff_nomu(n, m, erg, zq, pmax, o0i, o0il)
      do j = 1, 2
         part_pow(j) = 0.0d0
         do i = 1, pmax
            part_pow(j) = part_pow(j) +
     &                    (o0i(i,j) + log(zb)*o0il(i,j)) * zb**i
         enddo
      enddo

      
c     ==== Do the direct integration part ====
      part_dir(1) = 1/kp * ( eb**2/8 + bt0*eb/(2*kp**2) +
     &              (bt1_1 + gmm/(2*kp**2) + (3*bt0**2)/(2*kp**4)) *
     &              log(eb+1) )
      part_dir(2) = 1/kp * ( eb/2 - (bt1_2 - bt0/kp**2) * log(eb+1) )

      
      o_intl = part_pow + part_dir

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intl_part2( eb )

      implicit none


      o1i = 0.0d0
      o2i = 0.0d0
      do i = 1, ngr

         eta = pt(i) / kp
         call wfat_calc_o_intd(eta, , e1, e2, )

         o1i = o1i + wt * e1
         o2i = o2i + wt * e2
      enddo

      o1i = o1i / kp**2
      o2i = o2i / kp**2

c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intd(eta, n, m, erg, zq, e1, e2)

      implicit none
#include "errquit.fh"
      double precision eta
      integer n
      integer m
      double precision erg
      double precision zq
      double precision e1
      double precision e2

      integer j, pmax, md
      double precision kp, a, b, g, d, bt1_1, bt1_2, chgu, gm, v1, v2
      double precision ps1, ps2, ps3, dj, u0, chg_mu
      double precision factorial, rise_fact
      double complex zz, chgm
      parameter (pmax=40)
      external factorial, rise_fact


      if (kp*eta < 1.0d-14)
     &     call errquit('WFAT: When m=0, the argument eta of '//
     &     'wfat_calc_o_intd must be such that kp*eta >= 1.0E-14',
     &     28, UERR)
      if (eta < 0.0d0)
     &     call errquit('WFAT: The argument eta of wfat_calc_o_intd '//
     &     'must not be negative', 36, UERR)


      kp = sqrt(2 * abs(erg))
      a = 1 + abs(m) + n - zq/kp
      b = dble(abs(m) + 1)
      call wfat_gammaf(a, g)
      d = g / factorial(abs(m))      
      bt1_1 = -1.0d0/4/kp**2
     &        * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      bt1_2 = -1.0d0/2/kp * (2*n+1+abs(m))


c     u = (kp*eta)^|m| * U(a,b,kp*eta)
      if (kp*eta >= 0.0d0 .and. kp*eta < 1.0d0) then
         chgu = 0.0d0
         call wfat_gammaf(a-abs(m), gm)
         v1 = (-1)**b/factorial(abs(m))/gm
         do j = 0, pmax
            call wfat_psif(a+j, ps1)
            call wfat_psif(dble(1+j), ps2)
            call wfat_psif(dble(1+j+abs(m)), ps3)
            dj = ps1 - ps2 - ps3
            u0 = v1 * rise_fact(a,j) / rise_fact(b,j) / factorial(j)
     &           * (kp*eta)**(abs(m)+j) * (log(kp*eta) + dj)
            chgu = chgu + u0
         enddo

         if (abs(m) > 0) then
            v2 = 1.0d0 / g
            do j = 1, abs(m)
               u0 = v2 * factorial(j-1) * rise_fact(1-a+j,abs(m)-j)
     &              / factorial(abs(m)-j) * (kp*eta)**(abs(m)-j)
               chgu = chgu + u0
            enddo
         endif
      else
         call wfat_chgu(a, b, kp*eta, chgu, md)
         chgu = chgu * (kp*eta)**(abs(m))
      endif

      zz = cmplx(kp*eta, 0.0d0, kind=kind(1.0d0))
      call wfat_cchg(a, b, zz, chgm)
      chg_mu = real(chgm) * chgu

      e1 = d * (bt1_1*kp + kp*eta**2/4) * chg_mu
     &     - exp(kp*eta) * (eta/4 + bt0/2/kp**2 + (bt1_1 
     &     + gmm/2/kp**2 + (3*bt0**2)/(2*kp**4)) / (eta+1))

      e2 = d * (bt1_2*kp - kp*eta/2) * chg_mu
     &     + exp(kp*eta) * (1.0d0/2 - (bt1_2 - bt0/kp**2) / (eta+1))

      end
c=======================================================
