c===============================================================
c===============================================================
      subroutine wfat_2a1e_target(geom, m, sphi, target, eps, 
     &     nr, nt, init_ke, maxiter, cf, erg, dipole, ke)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer m
      integer sphi
      integer target
      double precision eps
      integer nr
      integer nt
      double precision init_ke
      integer maxiter
      double precision cf(nr*nt, nr*nt)
      double precision erg(nr*nt)
      double precision dipole(3)
      double precision ke

      integer l_dr, k_dr, l_dt, k_dt

      
      ke = init_ke
      call wfat_2a1e_solve(m, nr, nt, geom, ke, .true., maxiter, target,
     &     eps, cf, erg)
c      call wfat_2a1e_solve(m, nr, nt, geom, ke, .false., maxiter, 
c     &     target, eps, cf, erg)


c     ==== Obtain the electronic dipole moment ====
      dipole = 0.0d0
      if (.not. ma_alloc_get(MT_DBL, (nr*nr)**2, 'dr', l_dr, k_dr))
     &     call errquit('WFAT: Cannot allocate dr', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, (nt*nt)**2, 'dt', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt', 1, MA_ERR)
      call wfat_2a1e_dpo('z', m, 0, nr, nt, nr, nt, sphi, sphi, ke,
     &     cf(:,target), cf(:,target), .false., dbl_mb(k_dr), .false.,
     &     dbl_mb(k_dt), .false., dipole(3))

      
c     ==== Analysis ====
      call wfat_analyze_target(nr, nt, erg, dipole)

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_dr))
     &     call errquit('WFAT: Cannot deallocate dr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt', 1, MA_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_analyze_target(nr, nt, erg, dipole)

      implicit none
#include "stdio.fh"
      integer nr
      integer nt
      double precision erg(nr*nt)
      double precision dipole(3)

      integer i

      write(luout, '("A few lowest energies : ")')
      do i = 1, min(nr*nt, 10)
         write(luout, '(i4, f16.8)') i, erg(i)
      enddo

      write(luout, '("")')
      write(luout, '("Dipole moment : ")')
      write(luout, '(3f12.6)') dipole
      write(luout, '("")')

      end
c===============================================================
