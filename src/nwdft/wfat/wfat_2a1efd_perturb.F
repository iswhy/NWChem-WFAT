c===============================================================
c===============================================================
      subroutine wfat_2a1efd_perturb(geom, mch, sphi, orbid, nfe, ndvr, 
     &     maxn, xq, wq, nr, nt, nt_m, nt_p, pottype, nel, potparam, cf,
     &     erg, cfx, cfy, cfz, polb)

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
      double precision cf(nr*nt, nr*nt)    ! [input],
      double precision erg(nr*nt)          ! [input],
      double precision cfx(1:*)            ! [output], The minimum size is nr_p*nt_p if m=0, if m>0, the minimum size is nr_m*nt_m + nr_p*nt_p.
      double precision cfy(1:*)            ! [output], The minimum size is nr_p*nt_p if m=0, if m>0, the minimum size is nr_m*nt_m + nr_p*nt_p.
      double precision cfz(nr*nt)          ! [output],
      double precision polb(3)             ! [output],

      integer l_cfm, l_egm, l_cfp, l_egp, l_x
      integer k_cfm, k_egm, k_cfp, k_egp, k_x
      integer n_m, n_p
      
c     ==== Solve the m-1 eigenvectors ====
      if (mch >= 1) then
         n_m = nr * nt_m
         if (.not. ma_alloc_get(MT_DBL, n_m**2, 'cf-', l_cfm, k_cfm))
     &        call errquit('WFAT: Cannot allocate cf-', 1, MA_ERR)
         if (.not. ma_alloc_get(MT_DBL, n_m, 'erg-', l_egm, k_egm))
     &        call errquit('WFAT: Cannot allocate erg-', 1, MA_ERR)
         call wfat_2a1efd_solve(mch-1, nfe, ndvr, maxn, xq, wq, nt_m, 
     &        geom, pottype, nel, potparam, dbl_mb(k_cfm),
     &        dbl_mb(k_egm))
      endif

c     ==== Solve the m+1 eigenvectors ====
      n_p = nr * nt_p
      if (.not. ma_alloc_get(MT_DBL, n_p**2, 'cf+', l_cfp, k_cfp))
     &     call errquit('WFAT: Cannot allocate cf+', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, n_p, 'erg+', l_egp, k_egp))
     &     call errquit('WFAT: Cannot allocate erg+', 1, MA_ERR)
      call wfat_2a1efd_solve(mch+1, nfe, ndvr, maxn, xq, wq, nt_p, geom,
     &     pottype, nel, potparam, dbl_mb(k_cfp), dbl_mb(k_egp))

c     ==== Obtain the coefficients of the x, y, z ====
c     ====    corrections to the wavefunction     ====
      if (.not. ma_alloc_get(MT_DBL, nr, 'xquad_', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad_', 55, MA_ERR)
      call wfat_fedvr_basnodes(nfe, ndvr, maxn, xq, dbl_mb(k_x))
      call wfat_2a1efd_perturbcf(geom, mch, sphi, orbid, nr, nt, 
     &     nt_m, nt_p, dbl_mb(k_x), cf, dbl_mb(k_cfm), dbl_mb(k_cfp),
     &     erg, dbl_mb(k_egm), dbl_mb(k_egp), cfx, cfy, cfz)

c     ==== Obtain the polarizability ====
      call wfat_2a1efd_polb(mch, sphi, nr, nt, nt_m, nt_p, dbl_mb(k_x),
     &     cf(:,orbid), cfx, cfy, cfz, polb)
      write(luout, '("Polarizability : ", 3f15.6)') polb

c     ==== Deallocations ====
      if (mch >= 1) then
         if (.not. ma_free_heap(l_cfm))
     &        call errquit('WFAT: Cannot deallocate cf-', 41, MA_ERR)
         if (.not. ma_free_heap(l_egm))
     &        call errquit('WFAT: Cannot deallocate erg-', 41, MA_ERR)
      endif
      if (.not. ma_free_heap(l_cfp))
     &     call errquit('WFAT: Cannot deallocate cf+', 41, MA_ERR)
      if (.not. ma_free_heap(l_egp))
     &     call errquit('WFAT: Cannot deallocate erg+', 41, MA_ERR)

      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_2a1efd_perturbcf(geom, m, sphi, orbid, nr,
     &     nt, nt_m, nt_p, rq, cf, cf_m, cf_p, erg, erg_m, erg_p, cfx, 
     &     cfy, cfz)

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
      double precision cf(nr*nt, nr*nt)
      double precision cf_m(nr*nt_m, nr*nt_m)
      double precision cf_p(nr*nt_p, nr*nt_p)
      double precision erg(nr*nt)
      double precision erg_m(nr*nt_m)
      double precision erg_p(nr*nt_p)
      double precision cfx(1:*)
      double precision cfy(1:*)
      double precision cfz(nr*nt)

      integer j0, n, n_m, n_p, i, j, i_, sphi_, ndp
      integer l_dtm, l_dtp, l_dt, l_dp
      integer k_dtm, k_dtp, k_dt, k_dp
      double precision dpo


      j0 = orbid
      n = nr * nt
      
      
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
      

c     ==== x component of the corrected wave function ====
c     == m-1 sub-component ==
      if (m >= 1) then
         call wfat_2a1efd_dpo_col('x', m-1, m, sphi, sphi, nr, nt_m,
     &        nt, rq, dbl_mb(k_dtm), cf_m, cf(:,j0), dbl_mb(k_dp))
         do j = 1, n_m
            dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_m(j))
         enddo
         cfx(1:n_m) = matmul(cf_m, dbl_mb(k_dp:k_dp+n_m-1))
      endif
c     == m+1 sub-component ==
      call wfat_2a1efd_dpo_col('x', m+1, m, sphi, sphi, nr, nt_p,
     &     nt, rq, dbl_mb(k_dtp), cf_p, cf(:,j0), dbl_mb(k_dp))
      do j = 1, n_p
         dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_p(j))
      enddo
      cfx(n_m+1:n_m+n_p) = matmul(cf_p, dbl_mb(k_dp:k_dp+n_p-1))


c     ==== y component of the corrected wave function ====
      if (sphi==1) sphi_ = 2
      if (sphi==2) sphi_ = 1
c     == m-1 sub-component ==
      if (m >= 1) then
         call wfat_2a1efd_dpo_col('y', m-1, m, sphi_, sphi, nr, nt_m,
     &        nt, rq, dbl_mb(k_dtm), cf_m, cf(:,j0), dbl_mb(k_dp))
         do j = 1, n_m
            dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_m(j))
         enddo
         cfy(1:n_m) = matmul(cf_m, dbl_mb(k_dp:k_dp+n_m-1))
      endif
c     == m+1 sub-component ==
      call wfat_2a1efd_dpo_col('y', m+1, m, sphi_, sphi, nr, nt_p, nt, 
     &     rq, dbl_mb(k_dtp), cf_p, cf(:,j0), dbl_mb(k_dp))
      do j = 1, n_p
         dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_p(j))
      enddo
      cfy(n_m+1:n_m+n_p) = matmul(cf_p, dbl_mb(k_dp:k_dp+n_p-1))


c     ==== m,m ====
      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'dt_dvr', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt_dvr', 1, MA_ERR)
      call wfat_dipole_dvrt(m, 0, nt, nt, dbl_mb(k_dt))
      

c     ==== zz component of the static polarizability ====
      call wfat_2a1efd_dpo_col('z', m, m, sphi, sphi, nr, nt, nt,
     &     rq, dbl_mb(k_dt), cf, cf(:,j0), dbl_mb(k_dp))
      do j = 1, n
         if (j /= j0) then
            dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg(j))
         else
            dbl_mb(k_dp+j-1) = 0.0d0
         endif
      enddo
      cfz(1:n) = matmul(cf, dbl_mb(k_dp:k_dp+n-1))
      

c     ==== Deallocations ====
      if (m >= 1) then
         if (.not. ma_free_heap(l_dtm))
     &        call errquit('WFAT: Cannot deallocate dt_dvr-', 1, MA_ERR)
      endif
      if (.not. ma_free_heap(l_dtp))
     &     call errquit('WFAT: Cannot deallocate dt_dvr+', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt_dvr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dp))
     &     call errquit('WFAT: Cannot deallocate dp', 1, MA_ERR)
      
      end
c===============================================================
      
