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



      o1i = 0.0d0
      o2i = 0.0d0
      do i = 1, ngr

         
         call wfat_calc_o_intd(e1, e2, )

         o1i = o1i + wt * e1
         o2i = o2i + wt * e2
      enddo

      o1i = o1i / kappa
      o2i = o2i / kappa

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_o_intd(z, , , e1, e2)

      implicit none



      bt1_1 =
      bt1_2 = 


      eta = z / kappa
      if (z < 1.0d0) then
         do i = 0, pmax

         enddo

         
      else

         m =
         u = 
         chg_mu = m * u   ! chg_mu = z^(abs(m))*M(a,b,z)*U(a,b,z)

      endif

      f1 = d * (bt1_1*kappa + kappa*eta**2/4) * chg_mu
      e1 = f1
     &     - exp(z) * (eta/4 + bt0/2/kappa**2 + (bt1_1 + gmm/2/kappa**2
     &     + (3*bt0**2)/(2*kappa**4)) / (eta+1))

      f2 = d * (bt1_2*kappa - z/2) * chg_mu
      e2 = f2
     &     + exp(z) * (1.0d0/2 - (bt1_2 - bt0/(kappa**2))
     &     / (eta+1))
      
      

      end
c=======================================================
      
