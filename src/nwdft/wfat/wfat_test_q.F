c=======================================================
c=======================================================
      subroutine wfat_test_q()

      implicit none

      double precision protmat(3,3), mu(3)
      

      nr =
      n =
      m =
      lmax =
      nlgl = 

      r_first =
      r_last =

      kappa = sqrt(2 * abs(erg))
      lmb = 2 * kappa
      a = 1.0d0 - zq/kappa
      b = 2.0d0


c     ==== Get the passive rotation matrix ===
      call wfat_passrot(-0.0d0, -0.0d0, protmat)
      mu = matmul(transpose(protmat), -chmmt)      
      

      if (.not. MA_alloc_get(MT_DBL, nr, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 52, MA_ERR)
      dr = (r_last - r_first) / (nr-1)
      do i = 1, nr
         dbl_mb(k_r+i-1) = r_first + (i-1)*dr
      enddo

      nl = lmax - abs(m) + 1
      if (.not. MA_alloc_get(MT_DBL, nr, 'q_p1', l_qp1, k_qp1))
     &     call errquit('WFAT: q_p1 allocation failed', 52, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nr, 'q_p2', l_qp2, k_qp2))
     &     call errquit('WFAT: q_p2 allocation failed', 52, MA_ERR)
            

      call wfat_q_particular(nr, dbl_mb(k_r), n, m, lmax, erg, zq,
     &     nlgl, dbl_mb(k_qp1), dbl_mb(k_qp2))
      

      call util_file_name('qtest.wr1', .false., .false., fname)
      open(WF_QTEST_UNIT, file=fname, action='write')

      write(WF_QTEST_UNIT, '("#", i5, i15)', advance='no') 1, 2
      do i = 1, lmax-abs(m)+1
         j = 5*i - 2
         write(WF_QTEST_UNIT, '(6x, 5i15)', advance='no')
     &        j, j+1, j+2, j+3, j+4
      enddo
      write(WF_QTEST_UNIT, '(" ")')
      
      write(WF_QTEST_UNIT, '("#", a5, a15)', advance='no') 'No.', 'r'
      do l = abs(m), lmax
         write(WF_QTEST_UNIT, '(6x, 5a15)', advance='no')
     &        'LHS', 'RHS', 'epsilon', 'Q_l', 'R_l'
      enddo
      write(WF_QTEST_UNIT, '(" ")')
      

      do i = 1, nr-2
         r = dbl_mb(k_r+i-1)
         write(WF_QTEST_UNIT, '(1x, i5, es15.6)', advance='no') i, r
         call wfat_chgm_1pt(r, abs(m), lmax+1, a, b, lmb, chgm)

         do l = abs(m), lmax
c           ==== Calculate the derivatives ====
            i_q1_2h = k_qp1 + ((i+2)-1)*nl + l - abs(m)
            i_q1_1h = k_qp1 + ((i+1)-1)*nl + l - abs(m)
            i_q1_0h = k_qp1 + ((i+0)-1)*nl + l - abs(m)
            i_q1_1h_ = k_qp1 + ((i-1)-1)*nl + l - abs(m)
            i_q1_2h_ = k_qp1 + ((i-2)-1)*nl + l - abs(m)

            i_q2_2h = k_qp2 + ((i+2)-1)*nl + l - abs(m)
            i_q2_1h = k_qp2 + ((i+1)-1)*nl + l - abs(m)
            i_q2_0h = k_qp2 + ((i+0)-1)*nl + l - abs(m)
            i_q2_1h_ = k_qp2 + ((i-1)-1)*nl + l - abs(m)
            i_q2_2h_ = k_qp2 + ((i-2)-1)*nl + l - abs(m)

            f2h = dbl_mb(i_q1_2h) + mu(3)*dbl_mb(i_q2_2h)
            f1h = dbl_mb(i_q1_1h) + mu(3)*dbl_mb(i_q2_1h)
            f0h = dbl_mb(i_q1_0h) + mu(3)*dbl_mb(i_q2_0h)
            f1h_ = dbl_mb(i_q1_1h_) + mu(3)*dbl_mb(i_q2_1h_)
            f2h_ = dbl_mb(i_q1_2h_) + mu(3)*dbl_mb(i_q2_2h_)

            if (i == 1) then
               dif2 = 0.0d0
               dif1 = 0.0d0
            else if (i == 2) then
               dif2 = (f1h - 2*f0h + f1h_) / (dr**2)
               dif1 = (f1h - f1h_) / (2*dr)
            else
               dif2 = (-f2h + 16*f1h - 30*f0h + 16*f1h_ - f2h_)
     &                / (12*dr**2)
               dif1 = (-f2h + 8*f1h - 8*f1h_ + f2h_) / (12*dr)
            endif


c           ==== Calculate the terms on the RHS ====
            if (l == abs(m)) then
               r0_m = 0.0d0
            else
               call wfat_calc_radnorm(l-1, n, m, zq, kappa, c_m)
               k_m = sqrt(dble(l**2-m**2) / (4*l**2-1))
               r0_m = c_m * (kappa*r)**(l-1)*exp(-kappa*r) * chgm(l-1)
            endif

            call wfat_calc_radnorm(l+1, n, m, zq, kappa, c_p)
            k_p = sqrt(dble((l+1)**2-m**2) / (4*(l+1)**2-1))
            r0_p = c_p * (kappa*r)**(l+1)*exp(-kappa*r) * chgm(l+1)

            call wfat_calc_radnorm(l, n, m, zq, kappa, c)
            r0 = c * (kappa*r)**l*exp(-kappa*r) * chgm(l)


c           ==== Calculate the LHS and RHS ====
            lhs = r**2*dif2 + 2*r*dif1 +
     &            (-l*(l+1) + 2*zq*r + 2*erg*r**2) * f0h
            rhs = 2*k_m*r**3*r0_m + 2*k_p*r**3*r0_p + 2*mu(3)*r**2*r0
            epsi = lhs - rhs

            write(WF_QTEST_UNIT, '(6x, 5es15.6)', advance='no')
     &           lhs, rhs, epsi, f0h, r0
         enddo

         write(WF_QTEST_UNIT, '(" ")')
      enddo

      

      end
c=======================================================
      
