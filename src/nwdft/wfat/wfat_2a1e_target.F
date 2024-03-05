c===============================================================
c===============================================================
      subroutine wfat_2a1e_target(geom, m, sphi, target, pottype, nel, 
     &     potparam, nr, nt, init_ke, mult, maxiter, g_cf, erg, dipole,
     &     ke)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer geom
      integer m
      integer sphi
      integer target
      character*3 pottype
      integer nel                    ! [input], not used when pottype = 'cou', but still has to be supplied a value.
      double precision potparam(2)
      integer nr
      integer nt
      double precision init_ke
      double precision mult
      integer maxiter
      integer g_cf                   ! [input], GA of size nr*nt x nr*nt
      double precision erg(nr*nt)
      double precision dipole(3)
      double precision ke

      integer l_dr, k_dr, l_dt, k_dt, l_cf, k_cf, ndid

      
      ndid = ga_nodeid()

c     ==== Solve the eigenproblem for GS ====
      ke = init_ke
      call wfat_2a1e_solve(m, nr, nt, geom, ke, mult, .true., maxiter, 
     &     target, pottype, nel, potparam, g_cf, erg)


c     ==== Electronic dipole moment ====
      dipole = 0.0d0
      if (.not. ma_alloc_get(MT_DBL, (nr*nr)**2, 'dr', l_dr, k_dr))
     &     call errquit('WFAT: Cannot allocate dr', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, (nt*nt)**2, 'dt', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt', 1, MA_ERR)

      if (.not. ma_alloc_get(mt_dbl, nr*nt, 'coef', l_cf, k_cf))
     &     call errquit('WFAT: Could not allocate coef.', 34, MA_ERR)
      call ga_get(g_cf, 1, nr*nt, target, target, dbl_mb(k_cf), nr*nt)
      call wfat_2a1e_dpo('z', m, 0, nr, nt, nr, nt, sphi, sphi, ke,
     &     dbl_mb(k_cf), dbl_mb(k_cf), .false., dbl_mb(k_dr), .false.,
     &     dbl_mb(k_dt), .false., dipole(3))

      
c     ==== Analysis ====
      if (ndid == 0)
     &     call wfat_analyze_target(nr, nt, erg, dipole)

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_dr))
     &     call errquit('WFAT: Cannot deallocate dr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt', 1, MA_ERR)
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Could not deallocate coef', 7, MA_ERR)

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
