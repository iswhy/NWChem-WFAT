c============================================================
c============================================================
      subroutine wfat_calc_wfat1_rate(bt, gm, nat, nqn, nch, nid, mid,
     &     ef, erg, zi, gn, jn, kn, on, chmmt, pb, ffc, p_rt, rt)

      implicit none
#include "stdio.fh"
      double precision bt       ! [input],
      double precision gm                    ! [input],
      integer nat                            ! [input],
      integer nqn                            ! [input],
      integer nch                            ! [input],
      integer nid(nqn)                       ! [input],
      integer mid(nqn)                       ! [input],
      double precision ef                    ! [input],
      double precision erg(nch)              ! [input],
      double precision zi                    ! [input],
      double complex gn(nat, nqn, nch)       ! [input],
      double complex jn(nat, nqn, nch)       ! [input],
      double complex kn(nat, nqn, nch)       ! [input],
      double precision on(2, nqn, nch)       ! [input],
      double precision chmmt(3,nch)          ! [input],
      double precision pb(3,3,nch)           ! [input],
      double precision ffc(nqn, nch)         ! [output],
      double precision p_rt(nqn,nch)         ! [output],
      double precision rt                    ! [output],
      
      integer i, iqn, n, m
      double precision on_, protmat(3,3), mu(3), lpbz, kp, gmm
      double precision bt0, bt1_1, bt1_2, bt1, aa, bb_
      double complex gn_, an_

      
c      - fix dn for negative m.  
      

c     ==== Get the passive rotation matrix ===
      call wfat_passrot(bt, gm, protmat)
      

      do i = 1, nch
         mu = matmul(transpose(protmat), -chmmt(:,i))
         call wfat_calc_lfpolb(protmat, pb(:,:,i), lpbz)
         kp = sqrt( 2 * abs(erg(i)) )
         
         do iqn = 1, nqn
            n = nid(iqn)
            m = mid(iqn)

            gmm = (1.0d0 - m**2) / 4.0d0
            bt0 = zi - kp*(n + (abs(m)+1)/2.0d0)
            bt1_1 = -1.0d0/4/kp**2
     &              * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
            bt1_2 = -1.0d0/2/kp * (2*n+1+abs(m))
            bt1 = bt1_1 + mu(3)*bt1_2
            aa = -2*bt1/kp - (gmm-2*mu(3)*bt0)/kp**3 - 3*bt0**2/kp**5
            bb_ = -kp*lpbz - mu(3)**2/kp + mu(3)/kp**2 +
     &            4*mu(3)*bt0/kp**3 - (10+18*gmm+3*gmm**2)/24/kp**3 -
     &            (9-6*gmm)*bt0/4/kp**4 - (49+2*gmm)*bt0**2/8/kp**5 +
     &            3*bt0**3/2/kp**6 - bt0**4/8/kp**7

            write(luout, '("A_nu = ", f15.7)') aa
            write(luout, '("B_nu = ", f15.7)') bb_

            
            on_ = on(1,iqn,i) + mu(3)*on(2,iqn,i)
            call wfat_field_factor(ef, erg(i), zi, n, m, ffc(iqn, i))

            gn_ = sum(gn(:,iqn,i))
            an_ = sum(jn(:,iqn,i)) + sum(kn(:,iqn,i))
            p_rt(iqn,i) = ffc(iqn,i) *
     &           (abs(gn_)**2 * 
     &           (1+aa*ef*log(ef/4/kp**2)+bb_*ef+2*ef*on_) +
     &           2*ef*real(conjg(gn_)*an_))
         enddo
      enddo

      rt = sum(p_rt)
      
      end
c============================================================
