c=======================================================
c=======================================================
      subroutine wfat_omega_wfat1(geom, nch, nqn, erg, chpos, nid, mid,
     &     lmax, nat, nptp, centid, ptlim, gq, nlgl, zq, omg1, omg2)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                                         ! [input],
      integer nch                                          ! [input],
      integer nqn                                          ! [input],
      double precision erg(nch)                            ! [input],
      double precision chpos(3,nch)                        ! [input],
      integer nid(nqn)                                     ! [input],
      integer mid(nqn)                                     ! [input],
      integer lmax                                         ! [input],
      integer nat                                          ! [input],
      integer nptp                                         ! [input],
      integer centid(nat)                                  ! [input],
      integer ptlim(2,nat)                                 ! [input],
      integer gq                                           ! [input],
      integer nlgl                                         ! [input],
      double precision zq                                  ! [input],
      double complex omg1((1+lmax)**2, nqn, nch, nptp)     ! [output],
      double complex omg2((1+lmax)**2, nqn, nch, nptp)     ! [output],

      integer l_r, l_t, l_p, l_sp, l_id, l_q1, l_q2
      integer k_r, k_t, k_p, k_sp, k_id, k_q1, k_q2
      integer nsh, i, ii, j, n, m, nl, i_sp, k


      if (.not. MA_alloc_get(MT_DBL, nptp, 'r', l_r, k_r))
     &     call errquit('WFAT: r allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp, 't', l_t, k_t))
     &     call errquit('WFAT: t allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nptp, 'p', l_p, k_p))
     &     call errquit('WFAT: p allocation failed', 59, MA_ERR)
      nsh = (lmax + 1)**2 * nptp
      if (.not. MA_alloc_get(MT_DCPL, nsh, 'sphm', l_sp, k_sp))
     &     call errquit('WFAT: sphm allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_INT, nptp, 'indx', l_id, k_id))
     &     call errquit('WFAT: indx allocation failed', 59, MA_ERR)
      

      do i = 1, nch
c        ==== Get the grid points with the current channel center ====
         call wfat_get_gridpt_cents(geom, nat, centid, ptlim, 
     &        chpos(:,i), gq, dbl_mb(k_r), dbl_mb(k_t), dbl_mb(k_p))
         write(luout, '("max, min r = ", 2es21.12)')
     &        maxval(dbl_mb(k_r:k_r+nptp-1)),
     &        minval(dbl_mb(k_r:k_r+nptp-1))
         
         write(luout, '("Sorting for channel : ", i4, " .....")',
     &        advance='no') i
         do k = 1, nptp
            int_mb(k_id+k-1) = k
         enddo
         call wfat_dsortx(dbl_mb(k_r), 1, nptp, int_mb(k_id))
         write(luout, '(" [Finished]")')
         write(luout, '("herex1")')
         
c        ==== Calculate the spherical harmonics for all (l,m) ====
c        ====   where l=[0,lmax] and m=[-l,l] at each point   ====
         do ii = 1, nptp
            i_sp = k_sp + (ii-1)*(1+lmax)**2
            call wfat_sphm_1pt(lmax, dbl_mb(k_t+ii-1), dbl_mb(k_p+ii-1),
     &           dcpl_mb(i_sp))
         enddo
         write(luout, '("herex2")')

         do j = 1, nqn
            n = nid(j)
            m = mid(j)
            nl = lmax - abs(m) + 1
            if (.not. MA_alloc_get(MT_DBL, nl*nptp, 'q1', l_q1, k_q1))
     &           call errquit('WFAT: q1 allocation failed', 94, MA_ERR)
            if (.not. MA_alloc_get(MT_DBL, nl*nptp, 'q2', l_q2, k_q2))
     &           call errquit('WFAT: q2 allocation failed', 94, MA_ERR)
            write(luout, '("herex3")')

            
c           ==== Calculate Q_{(l,m)}^\nu at all points above ====
            call wfat_calc_q(nptp, dbl_mb(k_r), n, m, lmax, nlgl, 
     &           erg(i), zq, dbl_mb(k_q1), dbl_mb(k_q2))
            write(luout, '("herex4")')
            write(luout, '("max, min id outs = ", 2i7)')
     &           maxval(int_mb(k_id:k_id+nptp-1)),
     &           minval(int_mb(k_id:k_id+nptp-1))
            
            call wfat_reorder_q(lmax, m, nptp, int_mb(k_id),
     &           dbl_mb(k_q1), dbl_mb(k_q2))
            write(luout, '("herex5")')

c           ==== Calculate WFAT(1) Omega for all (l,m) and ====
c           ====             all points above              ====
            call wfat_omega_wfat1_nm(nptp, lmax, m, dbl_mb(k_q1), 
     &           dbl_mb(k_q2), dcpl_mb(k_sp), omg1(:,j,i,:),
     &           omg2(:,j,i,:))
            write(luout, '("herex6")')

            
            if (.not. MA_free_heap(l_q1)) call errquit('WFAT: q1 '//
     &           'deallocation failed', 32, MA_ERR)
            if (.not. MA_free_heap(l_q2)) call errquit('WFAT: q2 '//
     &           'deallocation failed', 32, MA_ERR)
         enddo
      enddo

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_r))
     &     call errquit('WFAT: r deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_t))
     &     call errquit('WFAT: t deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_p))
     &     call errquit('WFAT: p deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_sp))
     &     call errquit('WFAT: sphm deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_id))
     &     call errquit('WFAT: indx deallocation failed', 61, MA_ERR)
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
               omg1(i_sh, ii) = q_1(l, ii) * sphm(i_sh, ii)
               omg2(i_sh, ii) = q_2(l, ii) * sphm(i_sh, ii)
            enddo
         enddo
      enddo

      end
c=======================================================
      
