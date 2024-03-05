c===============================================================
c===============================================================
      subroutine wfat_2a1e_perturb(geom, m, sphi, orbid, pottype, nel, 
     &     potparam, ke, mult, nr, nt, nr_m, nt_m, nr_p, nt_p, g_cf, 
     &     erg, cfx, cfy, cfz, polb)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom                         ! [input],
      integer m                            ! [input],
      integer sphi                         ! [input],
      integer orbid                        ! [input],
      character*3 pottype
      integer nel                    ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)
      double precision ke                  ! [input],
      double precision mult                ! [input],
      integer nr                           ! [input],
      integer nt                           ! [input],
      integer nr_m                         ! [input],
      integer nt_m                         ! [input],
      integer nr_p                         ! [input],
      integer nt_p                         ! [input],
      integer g_cf                         ! [input], GA of size nr*nt x nr*nt
      double precision erg(nr*nt)          ! [input],
      double precision cfx(1:*)            ! [output], The minimum size is nr_p*nt_p if m=0, if m>0, the minimum size is nr_m*nt_m + nr_p*nt_p.
      double precision cfy(1:*)            ! [output], The minimum size is nr_p*nt_p if m=0, if m>0, the minimum size is nr_m*nt_m + nr_p*nt_p.
      double precision cfz(nr*nt)          ! [output],
      double precision polb(3)             ! [output],

      integer g_cfm, g_cfp, l_egm, k_egm, l_egp, k_egp
      integer n_m, n_p
      
c     ==== Solve the m-1 eigenvectors ====
      if (m >= 1) then
         n_m = nr_m * nt_m
         if (.not. ga_create(mt_dbl, n_m, n_m, 'cf-', n_m, 0, g_cfm))
     &        call errquit('WFAT: Cannot create cf- global array.',
     &        7, GA_ERR)
         if (.not. ma_alloc_get(MT_DBL, n_m, 'erg-', l_egm, k_egm))
     &        call errquit('WFAT: Cannot allocate erg-', 1, MA_ERR)
         call wfat_2a1e_solve(m-1, nr_m, nt_m, geom, ke, mult, .false., 
     &        1, 1, pottype, nel, potparam, g_cfm, dbl_mb(k_egm))
      endif

c     ==== Solve the m+1 eigenvectors ====
      n_p = nr_p * nt_p
      if (.not. ga_create(mt_dbl, n_p, n_p, 'cf+', n_p, 0, g_cfp))
     &     call errquit('WFAT: Cannot create cf+ global array.',
     &     7, GA_ERR)
      if (.not. ma_alloc_get(MT_DBL, n_p, 'erg+', l_egp, k_egp))
     &     call errquit('WFAT: Cannot allocate erg+', 1, MA_ERR)
      call wfat_2a1e_solve(m+1, nr_p, nt_p, geom, ke, mult, .false., 1,
     &     1, pottype, nel, potparam, g_cfp, dbl_mb(k_egp))
      
c     ==== Obtain the coefficients of the x, y, z ====
c     ====    corrections to the wavefunction     ====
      call wfat_2a1e_perturbcf(geom, m, sphi, orbid, ke, nr,
     &     nt, nr_m, nt_m, nr_p, nt_p, cf, g_cfm, g_cfp,
     &     erg, dbl_mb(k_egm), dbl_mb(k_egp), cfx, cfy, cfz)
      
c     ==== Obtain the polarizability ====
      call wfat_2a1e_polb(m, sphi, orbid, nr_m, nt_m, nr, nt, nr_p, 
     &     nt_p, ke, dbl_mb(k_egm), erg, dbl_mb(k_egp), dbl_mb(k_cfm),
     &     cf, dbl_mb(k_cfp), .false., polb)
      write(luout, '("Polarizability : ", 3f15.6)') polb

c     ==== Deallocations ====
      if (m >= 1) then
         if (.not. ga_destroy(g_cfm))
     &        call errquit('WFAT: Cannot destroy cf-.', 22, GA_ERR)
         if (.not. ma_free_heap(l_egm))
     &        call errquit('WFAT: Cannot deallocate erg-', 1, MA_ERR)
      endif
      if (.not. ga_destroy(g_cfp))
     &     call errquit('WFAT: Cannot destroy cf+.', 22, GA_ERR)
      if (.not. ma_free_heap(l_egp))
     &     call errquit('WFAT: Cannot deallocate erg+', 1, MA_ERR)

      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_2a1e_perturbcf(geom, m, sphi, orbid, ke, nr,
     &     nt, nr_m, nt_m, nr_p, nt_p, g_cf, g_cfm, g_cfp, erg, erg_m,
     &     erg_p, cfx, cfy, cfz)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer m
      integer sphi
      integer orbid
      double precision ke
      integer nr
      integer nt
      integer nr_m
      integer nt_m
      integer nr_p
      integer nt_p
      integer g_cf             ! [input], GA of size nr*nt x nr*nt
      integer g_cfm            ! [input], GA of size nr_m*nt_m x nr_m*nt_m
      integer g_cfp            ! [input], GA of size nr_p*nt_p x nr_p*nt_p
cOLD      double precision cf(nr*nt, nr*nt)
cOLD      double precision cf_m(nr_m*nt_m, nr_m*nt_m)
cOLD      double precision cf_p(nr_p*nt_p, nr_p*nt_p)
      double precision erg(nr*nt)
      double precision erg_m(nr_m*nt_m)
      double precision erg_p(nr_p*nt_p)
      double precision cfx(1:*)
      double precision cfy(1:*)
      double precision cfz(nr*nt)

      integer j0, n, n_m, n_p, i, j, i_, sphi_, ndp
      integer l_drm, l_dtm, l_drp, l_dtp, l_dr, l_dt, l_dp, l_cf
      integer k_drm, k_dtm, k_drp, k_dtp, k_dr, k_dt, k_dp, k_cf
      double precision dpo


      j0 = orbid
      n = nr * nt
      if (.not. ma_alloc_get(MT_DBL, nr*nt, 'coef', l_cf, k_cf))
     &     call errquit('WFAT: Cannot allocate coef.', 14, MA_ERR)
      call ga_get(g_cf, 1, nr*nt, j0, j0, dbl_mb(k_cf), nr*nt)
      
      
c     ==== m - 1 ====
      if (m >= 1) then
         n_m = nr_m * nt_m
         if (.not. ma_alloc_get(MT_DBL, nr_m*nr, 'dr_dvr-', l_drm, 
     &        k_drm))
     &        call errquit('WFAT: Cannot allocate dr_dvr-', 1, MA_ERR)
         call wfat_dipole_dvrr(m-1, nr_m, m, nr, ke, dbl_mb(k_drm))
   
         if (.not. ma_alloc_get(MT_DBL, nt_m*nt, 'dt_dvr-', l_dtm,
     &        k_dtm))
     &        call errquit('WFAT: Cannot allocate dt_dvr-', 1, MA_ERR)
         call wfat_dipole_dvrt(m, -1, nt_m, nt, dbl_mb(k_dtm))
      else
         n_m = 0
      endif

      
c     ==== m + 1 ====
      n_p = nr_p * nt_p
      if (.not. ma_alloc_get(MT_DBL, nr_p*nr, 'dr_dvr+', l_drp, k_drp))
     &     call errquit('WFAT: Cannot allocate dr_dvr+', 1, MA_ERR)
      call wfat_dipole_dvrr(m+1, nr_p, m, nr, ke, dbl_mb(k_drp))

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
         call wfat_2a1e_dpo_col('x', ke, m-1, m, sphi, sphi, nr_m, nt_m,
     &        nr, nt, dbl_mb(k_drm), dbl_mb(k_dtm), g_cfm, dbl_mb(k_cf),
     &        dbl_mb(k_dp))
         do j = 1, n_m
            dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_m(j))
         enddo
         cfx(1:n_m) = matmul(cf_m, dbl_mb(k_dp:k_dp+n_m-1))
      endif
c     == m+1 sub-component ==
      call wfat_2a1e_dpo_col('x', ke, m+1, m, sphi, sphi, nr_p, nt_p,
     &     nr, nt, dbl_mb(k_drp), dbl_mb(k_dtp), g_cfp, dbl_mb(k_cf),
     &     dbl_mb(k_dp))
      do j = 1, n_p
         dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_p(j))
      enddo
      cfx(n_m+1:n_m+n_p) = matmul(cf_p, dbl_mb(k_dp:k_dp+n_p-1))


c     ==== y component of the corrected wave function ====
      if (sphi==1) sphi_ = 2
      if (sphi==2) sphi_ = 1
c     == m-1 sub-component ==
      if (m >= 1) then
         call wfat_2a1e_dpo_col('y', ke, m-1, m, sphi_, sphi, nr_m, 
     &        nt_m, nr, nt, dbl_mb(k_drm), dbl_mb(k_dtm), g_cfm, 
     &        dbl_mb(k_cf), dbl_mb(k_dp))
         do j = 1, n_m
            dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_m(j))
         enddo
         cfy(1:n_m) = matmul(cf_m, dbl_mb(k_dp:k_dp+n_m-1))
      endif
c     == m+1 sub-component ==
      call wfat_2a1e_dpo_col('y', ke, m+1, m, sphi_, sphi, nr_p, nt_p, 
     &     nr, nt, dbl_mb(k_drp), dbl_mb(k_dtp), g_cfp, dbl_mb(k_cf),
     &     dbl_mb(k_dp))
      do j = 1, n_p
         dbl_mb(k_dp+j-1) = dbl_mb(k_dp+j-1) / (erg(j0)-erg_p(j))
      enddo
      cfy(n_m+1:n_m+n_p) = matmul(cf_p, dbl_mb(k_dp:k_dp+n_p-1))


c     ==== m,m ====
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'dr_dvr', l_dr, k_dr))
     &     call errquit('WFAT: Cannot allocate dr_dvr', 1, MA_ERR)
      call wfat_dipole_dvrr(m, nr, m, nr, ke, dbl_mb(k_dr))

      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'dt_dvr', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt_dvr', 1, MA_ERR)
      call wfat_dipole_dvrt(m, 0, nt, nt, dbl_mb(k_dt))
      

c     ==== zz component of the static polarizability ====
      call wfat_2a1e_dpo_col('z', ke, m, m, sphi, sphi, nr, nt, nr, nt,
     &     dbl_mb(k_dr), dbl_mb(k_dt), g_cf, dbl_mb(k_cf), dbl_mb(k_dp))
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
         if (.not. ma_free_heap(l_drm))
     &        call errquit('WFAT: Cannot deallocate dr_dvr-', 1, MA_ERR)
         if (.not. ma_free_heap(l_dtm))
     &        call errquit('WFAT: Cannot deallocate dt_dvr-', 1, MA_ERR)
      endif
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Cannot deallocate coef', 77, MA_ERR)
      if (.not. ma_free_heap(l_drp))
     &     call errquit('WFAT: Cannot deallocate dr_dvr+', 1, MA_ERR)
      if (.not. ma_free_heap(l_dtp))
     &     call errquit('WFAT: Cannot deallocate dt_dvr+', 1, MA_ERR)
      if (.not. ma_free_heap(l_dr))
     &     call errquit('WFAT: Cannot deallocate dr_dvr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt_dvr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dp))
     &     call errquit('WFAT: Cannot deallocate dp', 1, MA_ERR)
      
      end
c===============================================================


      
      
