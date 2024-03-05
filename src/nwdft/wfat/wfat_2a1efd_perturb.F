c===============================================================
c===============================================================
      subroutine wfat_2a1efd_perturb(geom, mch, sphi, orbid, nfe, ndvr, 
     &     maxn, xq, wq, nr, nt, nt_m, nt_p, pottype, nel, potparam, 
     &     g_cf, erg, cfx, cfy, cfz, polb)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                         ! [input],
      integer mch               ! [input],
      integer sphi                         ! [input],
      integer orbid                        ! [input],
      integer nfe
      integer ndvr(nfe)
      integer maxn
      double precision xq(maxn,nfe)
      double precision wq(maxn,nfe)
      integer nr                           ! [input],
      integer nt                           ! [input],
      integer nt_m                         ! [input],
      integer nt_p                         ! [input],
      character*3 pottype
      integer nel                          ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)
      integer g_cf                         ! [input], GA of size nr*nt x nr*nt
      double precision erg(nr*nt)          ! [input],
      double precision cfx(1:*)            ! [output], The minimum size is nr*nt_p if m=0, if m>0, the minimum size is nr*nt_m + nr*nt_p.
      double precision cfy(1:*)            ! [output], The minimum size is nr*nt_p if m=0, if m>0, the minimum size is nr*nt_m + nr*nt_p.
      double precision cfz(nr*nt)          ! [output],
      double precision polb(3)             ! [output],

      integer g_cfm, l_egm, k_egm, g_cfp, k_egp, l_egp, l_x, k_x
      integer l_cf, k_cf, n_m, n_p
      
c     ==== Solve the m-1 eigenvectors ====
      if (mch >= 1) then
         n_m = nr * nt_m
         if (.not. ga_create(mt_dbl, n_m, n_m, 'cf-', n_m, 0, g_cfm))
     &        call errquit('WFAT: Cannot create cf- global array.',
     &        1, GA_ERR)
         if (.not. ma_alloc_get(MT_DBL, n_m, 'erg-', l_egm, k_egm))
     &        call errquit('WFAT: Cannot allocate erg-', 1, MA_ERR)
         call wfat_2a1efd_solve(mch-1, nfe, ndvr, maxn, xq, wq, nt_m, 
     &        geom, pottype, nel, potparam, g_cfm, dbl_mb(k_egm))
      endif

c     ==== Solve the m+1 eigenvectors ====
      n_p = nr * nt_p
      if (.not. ga_create(mt_dbl, n_p, n_p, 'cf+', n_p, 0, g_cfp))
     &     call errquit('WFAT: Cannot create cf+ global array.',
     &     1, GA_ERR)
      if (.not. ma_alloc_get(MT_DBL, n_p, 'erg+', l_egp, k_egp))
     &     call errquit('WFAT: Cannot allocate erg+', 1, MA_ERR)
      call wfat_2a1efd_solve(mch+1, nfe, ndvr, maxn, xq, wq, nt_p, geom,
     &     pottype, nel, potparam, g_cfp, dbl_mb(k_egp))

c     ==== Obtain the coefficients of the x, y, z ====
c     ====    corrections to the wavefunction     ====
      if (.not. ma_alloc_get(MT_DBL, nr, 'xquad_', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad_', 55, MA_ERR)
      call wfat_fedvr_basnodes(nfe, ndvr, maxn, xq, dbl_mb(k_x))
      call wfat_2a1efd_perturbcf(geom, mch, sphi, orbid, nr, nt, 
     &     nt_m, nt_p, dbl_mb(k_x), g_cf, g_cfm, g_cfp,
     &     erg, dbl_mb(k_egm), dbl_mb(k_egp), cfx, cfy, cfz)

c     ==== Obtain the polarizability ====
      if (.not. ma_alloc_get(MT_DBL, nr*nt, 'coef', l_cf, k_cf))
     &     call errquit('WFAT: Cannot allocate coef.', 67, MA_ERR)
      call ga_get(g_cf, 1, nr*nt, orbid, orbid, dbl_mb(k_cf), nr*nt)
      call wfat_2a1efd_polb(mch, sphi, nr, nt, nt_m, nt_p, dbl_mb(k_x),
     &     dbl_mb(k_cf), cfx, cfy, cfz, polb)
      write(luout, '("Polarizability : ", 3f15.6)') polb

c     ==== Deallocations ====
      if (mch >= 1) then
         if (.not. ga_destroy(g_cfm))
     &        call errquit('WFAT: Cannot destroy cf-.', 45, GA_ERR)
         if (.not. ma_free_heap(l_egm))
     &        call errquit('WFAT: Cannot deallocate erg-', 41, MA_ERR)
      endif
      if (.not. ga_destroy(g_cfp))
     &     call errquit('WFAT: Cannot destroy cf+.', 45, GA_ERR)
      if (.not. ma_free_heap(l_egp))
     &     call errquit('WFAT: Cannot deallocate erg+', 41, MA_ERR)
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Cannot deallocate coef', 41, MA_ERR)

      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_2a1efd_perturbcf(geom, m, sphi, orbid, nr,
     &     nt, nt_m, nt_p, rq, g_cf, g_cfm, g_cfp, erg, erg_m, erg_p, 
     &     cfx, cfy, cfz)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer m
      integer sphi
      integer orbid
      integer nr
      integer nt
      integer nt_m
      integer nt_p
      double precision rq(nr)
      integer g_cf             ! [input], GA of size nr*nt x nr*nt
      integer g_cfm            ! [input], GA of size nr*nt_m x nr*nt_m
      integer g_cfp            ! [input], GA of size nr*nt_p x nr*nt_p
      double precision erg(nr*nt)
      double precision erg_m(nr*nt_m)
      double precision erg_p(nr*nt_p)
      double precision cfx(1:*)
      double precision cfy(1:*)
      double precision cfz(nr*nt)

      integer j0, n, n_m, n_p, i, j, i_, sphi_, ndp
      integer l_dtm, l_dtp, l_dt, l_dp, l_cf, l_tm
      integer k_dtm, k_dtp, k_dt, k_dp, k_cf, k_tm
      double precision dp, dpo


      j0 = orbid
      n = nr * nt
      if (.not. ma_alloc_get(MT_DBL, nr*nt, 'coef', l_cf, k_cf))
     &     call errquit('WFAT: Cannot allocate coef.', 42, MA_ERR)
      call ga_get(g_cf, 1, nr*nt, j0, j0, dbl_mb(k_cf), nr*nt)
      
      
c     ==== m - 1 ====
      if (m >= 1) then
         n_m = nr * nt_m   
         if (.not. ma_alloc_get(MT_DBL, nt_m*nt, 'dt_dvr-', l_dtm,
     &        k_dtm))
     &        call errquit('WFAT: Cannot allocate dt_dvr-', 1, MA_ERR)
         call wfat_dipole_dvrt(m, -1, nt_m, nt, dbl_mb(k_dtm))
      else
         n_m = 0
      endif

      
c     ==== m + 1 ====
      n_p = nr * nt_p
      if (.not. ma_alloc_get(MT_DBL, nt_p*nt, 'dt_dvr+', l_dtp, k_dtp))
     &     call errquit('WFAT: Cannot allocate dt_dvr+', 1, MA_ERR)
      call wfat_dipole_dvrt(m, 1, nt_p, nt, dbl_mb(k_dtp))

      
c     ==== Allocate temp vector to hold dipole column ====
      ndp = max(n, n_m, n_p)
      if (.not. ma_alloc_get(MT_DBL, ndp, 'dp', l_dp, k_dp))
     &     call errquit('WFAT: Cannot allocate dp', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, ndp, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: Cannot allocate tmp', 1, MA_ERR)
      

c     ==== x component of the corrected wave function ====
c     == m-1 sub-component ==
      if (m >= 1) then
         call wfat_2a1efd_dpo_col('x', m-1, m, sphi, sphi, nr, nt_m,
     &        nt, rq, dbl_mb(k_dtm), g_cfm, dbl_mb(k_cf), dbl_mb(k_dp))
         use GA operations below
cOLD         do j = 1, n_m
cOLD            dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_m(j))
cOLD         enddo
cOLD         cfx(1:n_m) = matmul(cf_m, dbl_mb(k_dp:k_dp+n_m-1))
         cfx(1:n_m) = 0.0d0
         do j = 1, n_m
            call ga_get(g_cfm, 1, n_m, j, j, dbl_mb(k_tm), n_m)
            dp = dbl_mb(k_dp+j-1) / (erg(j0)-erg_m(j))
            cfx(1:n_m) = cfx(1:n_m) + dp * dbl_mb(k_tm:k_tm+n_m-1)
         enddo
      endif
c     == m+1 sub-component ==
      call wfat_2a1efd_dpo_col('x', m+1, m, sphi, sphi, nr, nt_p,
     &     nt, rq, dbl_mb(k_dtp), g_cfp, dbl_mb(k_cf), dbl_mb(k_dp))
cOLD      do j = 1, n_p
cOLD         dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_p(j))
cOLD      enddo
cOLD      cfx(n_m+1:n_m+n_p) = matmul(cf_p, dbl_mb(k_dp:k_dp+n_p-1))
      cfx(n_m+1:n_m+n_p) = 0.0d0
      do j = 1, n_p
         call ga_get(g_cfp, 1, n_p, j, j, dbl_mb(k_tm), n_p)
         dp = dbl_mb(k_dp+j-1) / (erg(j0)-erg_p(j))
         cfx(n_m+1:n_m+n_p) = cfx(n_m+1:n_m+n_p) + dp *
     &                        dbl_mb(k_tm:k_tm+n_p-1)
      enddo


c     ==== y component of the corrected wave function ====
      if (sphi==1) sphi_ = 2
      if (sphi==2) sphi_ = 1
c     == m-1 sub-component ==
      if (m >= 1) then
         call wfat_2a1efd_dpo_col('y', m-1, m, sphi_, sphi, nr, nt_m,
     &        nt, rq, dbl_mb(k_dtm), g_cfm, dbl_mb(k_cf), dbl_mb(k_dp))
cOLD         do j = 1, n_m
cOLD            dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_m(j))
cOLD         enddo
cOLD         cfy(1:n_m) = matmul(cf_m, dbl_mb(k_dp:k_dp+n_m-1))
         cfy(1:n_m) = 0.0d0
         do j = 1, n_m
            call ga_get(g_cfm, 1, n_m, j, j, dbl_mb(k_tm), n_m)
            dp = dbl_mb(k_dp+j-1) / (erg(j0)-erg_m(j))
            cfy(1:n_m) = cfy(1:n_m) + dp * dbl_mb(k_tm:k_tm+n_m-1)
         enddo
      endif
c     == m+1 sub-component ==
      call wfat_2a1efd_dpo_col('y', m+1, m, sphi_, sphi, nr, nt_p, nt, 
     &     rq, dbl_mb(k_dtp), g_cfp, dbl_mb(k_cf), dbl_mb(k_dp))
cOLD      do j = 1, n_p
cOLD         dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_p(j))
cOLD      enddo
cOLD      cfy(n_m+1:n_m+n_p) = matmul(cf_p, dbl_mb(k_dp:k_dp+n_p-1))
      cfy(n_m+1:n_m+n_p) = 0.0d0
      do j = 1, n_p
         call ga_get(g_cfp, 1, n_p, j, j, dbl_mb(k_tm), n_p)
         dp = dbl_mb(k_dp+j-1) / (erg(j0)-erg_p(j))
         cfy(n_m+1:n_m+n_p) = cfy(n_m+1:n_m+n_p) + dp *
     &                        dbl_mb(k_tm:k_tm+n_p-1)
      enddo


c     ==== m,m ====
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'dt_dvr', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt_dvr', 1, MA_ERR)
      call wfat_dipole_dvrt(m, 0, nt, nt, dbl_mb(k_dt))
      

c     ==== zz component of the static polarizability ====
      call wfat_2a1efd_dpo_col('z', m, m, sphi, sphi, nr, nt, nt,
     &     rq, dbl_mb(k_dt), g_cf, dbl_mb(k_cf), dbl_mb(k_dp))
cOLD      do j = 1, n
cOLD         if (j /= j0) then
cOLD            dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg(j))
cOLD         else
cOLD            dbl_mb(k_dp+j-1) = 0.0d0
cOLD         endif
cOLD      enddo
cOLD      cfz(1:n) = matmul(cf, dbl_mb(k_dp:k_dp+n-1))
      cfz(1:n) = 0.0d0
      do j = 1, n
         if (j /= j0) then
            call ga_get(g_cf, 1, n, j, j, dbl_mb(k_tm), n)
            dp = dbl_mb(k_dp+j-1) / (erg(j0)-erg(j))
            cfz(1:n) = cfz(1:n) + dp * dbl_mb(k_tm:k_tm+n-1)
         endif
      enddo
      

c     ==== Deallocations ====
      if (m >= 1) then
         if (.not. ma_free_heap(l_dtm))
     &        call errquit('WFAT: Cannot deallocate dt_dvr-', 1, MA_ERR)
      endif
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Cannot deallocate coef', 83, MA_ERR)
      if (.not. ma_free_heap(l_dtp))
     &     call errquit('WFAT: Cannot deallocate dt_dvr+', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt_dvr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dp))
     &     call errquit('WFAT: Cannot deallocate dp', 1, MA_ERR)
      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp', 1, MA_ERR)
      
      end
c===============================================================
      
