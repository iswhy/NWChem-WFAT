c=======================================================
c=======================================================
      subroutine wfat_omega_wfat1(geom, nch, nqn, erg, chpos, nid, mid,
     &     lmax, nat, nptp, centid, ptlim, gq, nlgl, zq, omg1, omg2)

      implicit none
      integer geom
      integer nch
      integer nqn
      double precision erg(nch)
      double precision chpos(3,nch)
      integer nid(nqn)
      integer mid(nqn)
      integer lmax
      integer nat
      integer nptp
      integer centid(nat)
      integer ptlim(2,nat)
      integer gq
      integer nlgl
      double precision zq
      double complex omg1((1+lmax)**2, nqn, nch, nptp)
      double complex omg2((1+lmax)**2, nqn, nch, nptp)

      integer i

         

      do i = 1, nch
c        ==== Get the grid points with the current channel center ====
         call wfat_get_gridpt_cents(geom, nat, centid, ptlim, 
     &        chpos(:,i), gq, dbl_mb(k_r), dbl_mb(k_t), dbl_mb(k_p))
         call wfat_dsortx(dbl_mb(k_r), 1, nptp, int_mb(k_id))
         
c        ==== Calculate the spherical harmonics for all (l,m) ====
c        ====   where l=[0,lmax] and m=[-l,l] at each point   ====
         do ii = 1, nptp
            i_sp = k_sp + (ii-1)*(1+lmax)**2
            call wfat_sphm_1pt(lmax, dbl_mb(k_t+ii-1), dbl_mb(k_p+ii-1),
     &           dcpl_mb(i_sp))        check here
         enddo

         do j = 1, nqn
            n = nid(j)
            m = mid(j)

c           ==== Calculate Q_{(l,m)}^\nu at all points above ====
            call wfat_calc_q(nptp, dbl_mb(k_r), n, m, lmax, nlgl, 
     &           erg(i), zq, dbl_mb(k_q1), dbl_mb(k_q2))       check here
            call wfat_reorder_q(lmax, m, nptp, int_mb(k_id),
     &           dbl_mb(k_q1), dbl_mb(k_q2))

c           ==== Calculate WFAT(1) Omega for all (l,m) and ====
c           ====             all points above              ====
            call wfat_omega_wfat1_nm(dbl_mb(k_q1), dbl_mb(k_q2),
     &           dcpl_mb(k_sp), omg1(:,j,i,:), omg2(:,j,i,:))
         enddo
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_omega_wfat1_nm(nptp, lmax, m, q_1, q_2, sphm,
     &     omg1, omg2)

c     Calculate omg1 and omg2 for all (l,m) and all points at a particular channel 
c     and parabolic quantum number.
      implicit none
      integer nptp
      integer lmax
      integer m
      double precision q_1(abs(m):lmax, nptp)       ! [input], orientation-independent
      double precision q_2(abs(m):lmax, nptp)       ! [input], orientation-dependent
      double complex sphm((1+lmax)**2, nptp)        ! [input]
      double complex omg1((1+lmax)**2, nptp)        ! [output]
      double complex omg2((1+lmax)**2, nptp)        ! [output]
      
      integer ii, l, m_, i_sh
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )

      omg1 = czero
      omg2 = czero
      do ii = 1, nptp
         do l = abs(m), lmax
            do m_ = -l, l
               i_sh = l**2 + l + m + 1
               omg1(ish, ii) = q_1(l, ii) * sphm(ish, ii)
               omg2(ish, ii) = q_2(l, ii) * sphm(ish, ii)
            enddo
         enddo
      enddo

      end
c=======================================================
      
