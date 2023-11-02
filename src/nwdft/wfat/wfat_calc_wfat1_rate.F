c============================================================
c============================================================
      subroutine wfat_calc_wfat1_rate(bt, gm, nat, nqn, nch, nid, mid,
     &     ef, erg, zi, gn, jn, kn, on, chmmt, pb, justa, ffc, p_rt, rt)

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
      logical justa                          ! [input],
      double precision ffc(nqn, nch)         ! [output],
      double precision p_rt(nqn,nch)         ! [output],
      double precision rt                    ! [output],
      
      integer i, iqn, n, m
      double precision on_, protmat(3,3), mu(3), lpbz, kp, gmm
      double precision bt0, bt1_1, bt1_2, bt1, aa, bb_, ffc_, bb2
      double complex gn_, an_, an2

      
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

            
cfordebug            write(luout, '("B_nu = ", f15.7)') bb_

            on_ = on(1,iqn,i) + mu(3)*on(2,iqn,i)
            call wfat_field_factor(ef, erg(i), zi, n, m, ffc(iqn, i))
            if (justa) then
               ffc_ = 1.0d0
            else
               ffc_ = ffc(iqn,i)
            endif

            gn_ = sum(gn(:,iqn,i))
            an_ = sum(jn(:,iqn,i)) + sum(kn(:,iqn,i))
            p_rt(iqn,i) = ffc_ *
     &           (abs(gn_)**2 * 
     &           (1+aa*ef*log(ef/4/kp**2)+bb_*ef+2*ef*on_) +
     &           2*ef*real(conjg(gn_)*an_))



c           ==== for debug ====
            an2 = on_ + an_/gn_
            bb2 = bb_ + 2 * real(an2)


            write(luout, '("A_nu = ", 6f15.7)')
     &           bt, gm, aa, lpbz, pb(3,3,1), mu(3)
            write(luout, '("   an = ", 2f15.7)') an2
         enddo
      enddo

      
      rt = sum(p_rt)


c     ==== for debug ====
      rt = bb2
      
      end
c============================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_wfat1_yield(bt, gm, zi, nlas, ftlas, alas, 
     &                 wlas, plas, tlas, dlas, nt, tmax, extf, nch, erg,
     &                 nqn, nid, mid, nat, gn, jn, kn, gnn, jnn, knn,
     &                 on, chmmt, pb, yield)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "wfat_field_param.fh"
      double precision bt           ! [input],
      double precision gm           ! [input],
      double precision zi
      integer nlas
      character(len=FLD_TP_LEN) ftlas(nlas)       ! Function type of envelope, 'gaussian' or 'hann'         
      double precision alas(nlas)   ! Peak electric field, in a.u. of e. field.                              
      double precision wlas(nlas)   ! angular frequency, in a.u. (Hartree).                                      
      double precision plas(nlas)   ! Carrier-envelope phase, rad                          
      double precision tlas(nlas)   ! Delay (temporal position of the envelope peak), a.u. of time
      double precision dlas(nlas)   ! Field FWHM, a.u. of time                                       
      integer nt                    
      double precision tmax         ! Time window width, a.u. of time
      logical extf
      integer nch
      double precision erg(nch)
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      integer nat
      double complex gn(nat, nqn, nch)       ! [input],
      double complex jn(nat, nqn, nch)       ! [input],
      double complex kn(nat, nqn, nch)       ! [input],
      double complex gnn(nat, nqn, nch)      ! [input],
      double complex jnn(nat, nqn, nch)      ! [input],
      double complex knn(nat, nqn, nch)      ! [input],
      double precision on(2, nqn, nch)       ! [input],
      double precision chmmt(3,nch)          ! [input],
      double precision pb(3,3,nch)           ! [input],
      double precision yield
      
      integer i, j, k, npi, l_pi, k_pi, nfc, l_fc, k_fc
      double precision t0, dt, t, las, rate
      parameter (t0=0.0d0)

      
      npi = nqn * nch
      if (.not. ma_alloc_get(mt_dbl, npi, 'partial_ion', l_pi, k_pi))
     &     call errquit('WFAT: Cannot allocate partial_ion', 1, MA_ERR)
      nfc = nqn * nch
      if (.not. MA_alloc_get(MT_DBL, nfc, 'ffc', l_fc, k_fc))
     &     call errquit('WFAT: ffc allocation failed', 2, MA_ERR)
      
      
      dt = (tmax-t0)/(nt-1)

      yield = 0.0d0
      do i = 1, nt
         t = t0 + (i-1)*dt
         if (extf) then
            call wfat_eval_extlaser(i, las)
cdebug            las = 0.5d0 * las
cdebug            write(luout, '("laser: ", i5, es14.6)') i, las
         else
            call wfat_eval_laser(nlas, ftlas, alas, wlas, plas, tlas, 
     &           dlas, t, las)
         endif
         
         dbl_mb(k_pi:k_pi+npi-1) = 0.0d0
         if (las >= 1.0d-13) then
            call wfat_calc_wfat1_rate(bt, gm, nat, nqn, nch, nid,
     &           mid, las, erg, zi, gn, jn, kn, on, chmmt, pb, .false.,
     &           dbl_mb(k_fc), dbl_mb(k_pi), rate)
         else if (las <= -1.0d-13) then
            call wfat_calc_wfat1_rate(180.0d0-bt, 180.0d0+gm, nat, nqn,
     &           nch, nid, mid, abs(las), erg, zi, gnn, jnn, knn, on, 
     &           chmmt, pb, .false., dbl_mb(k_fc), dbl_mb(k_pi), rate)
         else
            dbl_mb(k_fc : k_fc+nfc-1) = 0.0d0
            dbl_mb(k_pi : k_pi+npi-1) = 0.0d0
            rate = 0.0d0
         endif

         yield = yield + rate * dt
      enddo

      yield = 1.0d0 - exp(-yield)

      if (.not. ma_free_heap(l_pi))
     &     call errquit('WFAT: Cannot deallocate partial_ion',
     &     1, MA_ERR)
      if (.not. MA_free_heap(l_fc))
     &     call errquit('WFAT: ffc deallocation failed',
     &     3, MA_ERR)
      
      end
c=======================================================
      
