c===============================================================
c===============================================================
      subroutine wfat_2a1e_perturb(geom, m, sphi, orbid, eps, ke, nr,
     &     nt, nr_m, nt_m, nr_p, nt_p, cf, erg, cfx, cfy, cfz, polb)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer m
      integer sphi
      integer orbid
      double precision eps
      double precision ke
      integer nr
      integer nt
      integer nr_m
      integer nt_m
      integer nr_p
      integer nt_p
      double precision cf(nr*nt, nr*nt)
      double precision erg(nr*nt)
      double precision cfx(1:*)
      double precision cfy(1:*)
      double precision cfz(nr*nt)
      double precision polb(3)

      integer l_cfm, l_egm, l_cfp, l_egp
      integer k_cfm, k_egm, k_cfp, k_egp
      integer n_m, n_p
      

      if (m >= 1) then
         n_m = nr_m * nt_m
         if (.not. ma_alloc_get(MT_DBL, n_m**2, 'cf-', l_cfm, k_cfm))
     &        call errquit('WFAT: Cannot allocate cf-', 1, MA_ERR)
         if (.not. ma_alloc_get(MT_DBL, n_m, 'erg-', l_egm, k_egm))
     &        call errquit('WFAT: Cannot allocate erg-', 1, MA_ERR)
         call wfat_2a1e_solve(m-1, nr_m, nt_m, geom, ke, .false., 1, 1,
     &        eps, dbl_mb(k_cfm), dbl_mb(k_egm))
      endif


      n_p = nr_p * nt_p
      if (.not. ma_alloc_get(MT_DBL, n_p**2, 'cf+', l_cfp, k_cfp))
     &     call errquit('WFAT: Cannot allocate cf+', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, n_p, 'erg+', l_egp, k_egp))
     &     call errquit('WFAT: Cannot allocate erg+', 1, MA_ERR)
      call wfat_2a1e_solve(m+1, nr_p, nt_p, geom, ke, .false., 1, 1,
     &     eps, dbl_mb(k_cfp), dbl_mb(k_egp))
      

      call wfat_2a1e_perturbcf(geom, m, sphi, orbid, eps, ke, nr,
     &     nt, nr_m, nt_m, nr_p, nt_p, cf, dbl_mb(k_cfm), dbl_mb(k_cfp),
     &      erg, dbl_mb(k_egm), dbl_mb(k_egp), cfx, cfy, cfz)

      call wfat_2a1e_polb(m, sphi, orbid, nr_m, nt_m, nr, nt, nr_p, 
     &     nt_p, ke, dbl_mb(k_egm), erg, dbl_mb(k_egp), dbl_mb(k_cfm),
     &     cf, dbl_mb(k_cfp), .false., polb)


      if (m >= 1) then
         if (.not. ma_free_heap(l_cfm))
     &        call errquit('WFAT: Cannot deallocate cf-', 1, MA_ERR)
         if (.not. ma_free_heap(l_egm))
     &        call errquit('WFAT: Cannot deallocate erg-', 1, MA_ERR)
      endif
      if (.not. ma_free_heap(l_cfp))
     &     call errquit('WFAT: Cannot deallocate cf+', 1, MA_ERR)
      if (.not. ma_free_heap(l_egp))
     &     call errquit('WFAT: Cannot deallocate erg+', 1, MA_ERR)

      end
c===============================================================

      
c===============================================================
c===============================================================
      subroutine wfat_2a1e_perturbcf(geom, m, sphi, orbid, eps, ke, nr,
     &     nt, nr_m, nt_m, nr_p, nt_p, cf, cf_m, cf_p, erg, erg_m,
     &     erg_p, cfx, cfy, cfz)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer m
      integer sphi
      integer orbid
      double precision eps
      double precision ke
      integer nr
      integer nt
      integer nr_m
      integer nt_m
      integer nr_p
      integer nt_p
      double precision cf(nr*nt, nr*nt)
      double precision cf_m(nr_m*nt_m, nr_m*nt_m)
      double precision cf_p(nr_p*nt_p, nr_p*nt_p)
      double precision erg(nr*nt)
      double precision erg_m(nr_m*nt_m)
      double precision erg_p(nr_p*nt_p)
      double precision cfx(1:*)
      double precision cfy(1:*)
      double precision cfz(nr*nt)

      integer j0, n, n_m, n_p, i, j, i_, sphi_
      integer l_drm, l_dtm, l_drp, l_dtp, l_dr, l_dt
      integer k_drm, k_dtm, k_drp, k_dtp, k_dr, k_dt
      double precision dpo


      j0 = orbid
      n = nr * nt

      
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


c     ==== x component of the corrected wave function ====
c     == m-1 sub-component ==
      if (n_m >= 1) then
         
         do i = 1, n_m   ! Loop over eigenstates of the m-1 spectrum.
            cfx(i) = 0.0d0
            do j = 1, n_m    ! Loop over bases of the m-1 spectrum.
               call wfat_2a1e_dpo('x', m, -1, nr_m, nt_m, nr, nt, sphi, 
     &              sphi, ke, cf_m(:,j), cf(:,j0), .true., 
     &              dbl_mb(k_drm), .true., dbl_mb(k_dtm), .false., dpo)
               cfx(i) = cfx(i) + cf_m(i,j) * dpo / (erg(j0)-erg_m(j))
            enddo
         enddo
         
      endif
c     == m+1 sub-component ==
      do i = n_m+1, n_m+n_p    ! Loop over eigenstates of the m+1 spectrum.
         i_ = i - n_m
         cfx(i) = 0.0d0
         do j = 1, n_p             ! Loop over bases of the m+1 spectrum.
            call wfat_2a1e_dpo('x', m, +1, nr_p, nt_p, nr, nt, sphi, 
     &           sphi, ke, cf_p(:,j), cf(:,j0), .true., 
     &           dbl_mb(k_drp), .true., dbl_mb(k_dtp), .false., dpo)
            cfx(i) = cfx(i) + cf_p(i_,j) * dpo / (erg(j0)-erg_p(j))
         enddo
      enddo


c     ==== y component of the corrected wave function ====
      if (sphi==1) sphi_ = 2
      if (sphi==2) sphi_ = 1
c     == m-1 sub-component ==
      if (n_m >= 1) then
         
         do i = 1, n_m   ! Loop over eigenstates of the m-1 spectrum.
            cfy(i) = 0.0d0
            do j = 1, n_m    ! Loop over bases of the m-1 spectrum.
               call wfat_2a1e_dpo('y', m, -1, nr_m, nt_m, nr, nt, sphi_,
     &              sphi, ke, cf_m(:,j), cf(:,j0), .true., 
     &              dbl_mb(k_drm), .true., dbl_mb(k_dtm), .false., dpo)
               cfy(i) = cfy(i) + cf_m(i,j) * dpo / (erg(j0)-erg_m(j))
            enddo
         enddo
         
      endif
c     == m+1 sub-component ==
      do i = n_m+1, n_m+n_p    ! Loop over eigenstates of the m+1 spectrum.
         i_ = i - n_m
         cfy(i) = 0.0d0
         do j = 1, n_p             ! Loop over bases of the m+1 spectrum.
            call wfat_2a1e_dpo('y', m, +1, nr_p, nt_p, nr, nt, sphi_, 
     &           sphi, ke, cf_p(:,j), cf(:,j0), .true., 
     &           dbl_mb(k_drp), .true., dbl_mb(k_dtp), .false., dpo)
            cfy(i) = cfy(i) + cf_p(i_,j) * dpo / (erg(j0)-erg_p(j))
         enddo
      enddo


c     ==== m,m ====
      if (.not. ma_alloc_get(MT_DBL, nr*nr, 'dr_dvr', l_dr, k_dr))
     &     call errquit('WFAT: Cannot allocate dr_dvr', 1, MA_ERR)
      call wfat_dipole_dvrr(m, nr, m, nr, ke, dbl_mb(k_dr))

      if (.not. ma_alloc_get(MT_DBL, nt*nt, 'dt_dvr', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt_dvr', 1, MA_ERR)
      call wfat_dipole_dvrt(m, 0, nt, nt, dbl_mb(k_dt))
      

c     ==== zz component of the static polarizability ====
      do i = 1, n   ! Loop over eigenstates of the m spectrum.
         cfz(i) = 0.0d0
         do j = 1, n    ! Loop over bases of the m spectrum.
            if (j /= j0) then
               call wfat_2a1e_dpo('z', m, 0, nr, nt, nr, nt, sphi,
     &              sphi, ke, cf(:,j), cf(:,j0), .true., dbl_mb(k_dr), 
     &              .true., dbl_mb(k_dt), .false., dpo)
               cfz(i) = cfz(i) + cf(i,j0) * dpo / (erg(j0)-erg(j))
            endif
         enddo
      enddo
      

c     ==== Deallocations ====
      if (m >= 1) then
         if (.not. ma_free_heap(l_drm))
     &        call errquit('WFAT: Cannot deallocate dr_dvr-', 1, MA_ERR)
         if (.not. ma_free_heap(l_dtm))
     &        call errquit('WFAT: Cannot deallocate dt_dvr-', 1, MA_ERR)
      endif
      if (.not. ma_free_heap(l_drp))
     &     call errquit('WFAT: Cannot deallocate dr_dvr+', 1, MA_ERR)
      if (.not. ma_free_heap(l_dtp))
     &     call errquit('WFAT: Cannot deallocate dt_dvr+', 1, MA_ERR)
      if (.not. ma_free_heap(l_dr))
     &     call errquit('WFAT: Cannot deallocate dr_dvr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt_dvr', 1, MA_ERR)
      
      end
c===============================================================
      
