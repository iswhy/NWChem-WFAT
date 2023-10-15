c===============================================================
c===============================================================
      subroutine wfat_2a1e_prep_system(rtdb, mch, sphi, spinid, orbid, 
     &     eps, nr, nt, init_ke, maxiter)

      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb               ! [input],
      integer mch                ! [output],
      integer sphi               ! [output],
      integer spinid             ! [output],
      integer orbid              ! [output],
      double precision eps       ! [output],
      integer nr                 ! [output],
      integer nt                 ! [output],
      double precision init_ke   ! [output],
      integer maxiter            ! [output],


c     ==== The m (magnetic) quantum number ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_m', mt_int, 1, mch))
     &     call errquit('WFAT: Failed reading 2a1e_m from rtdb',
     &     1, RTDB_ERR)

c     ==== The k quantum number (the xz-plane reflection symmetry) ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_k', mt_int, 1, sphi))
     &     call errquit('WFAT: Failed reading 2a1e_k from rtdb',
     &     1, RTDB_ERR)
      
c     ==== The eigenstate ID within the above (m,k) ====
      if (.not. rtdb_get(rtdb, 'wfat:orbid', mt_int, 1, orbid))
     &     call errquit('WFAT: Failed reading orbid from rtdb',
     &     18, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:spinid', mt_int, 1, spinid))
     &     call errquit('WFAT: Failed reading spinid from rtdb',
     &     18, INPUT_ERR)

c     ==== The soft Coulomb parameter ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_eps', mt_dbl, 1, eps))
     &     call errquit('WFAT: Failed reading 2a1e_eps from rtdb',
     &     1, RTDB_ERR)

c     ==== The number of radial DVR basis ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_nr', mt_int, 1, nr))
     &     call errquit('WFAT: Failed reading 2a1e_nr from rtdb',
     &     1, RTDB_ERR)

c     ==== The number of polar DVR basis ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_nt', mt_int, 1, nt))
     &     call errquit('WFAT: Failed reading 2a1e_nt from rtdb',
     &     1, RTDB_ERR)

c     ==== Initial basis exponent ===
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_inke', mt_dbl, 1, init_ke))
     &     call errquit('WFAT: Failed reading 2a1e_inke from rtdb',
     &     1, RTDB_ERR)

c     ==== Maximum iteration ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_maxit', mt_int, 1, maxiter))
     &     call errquit('WFAT: Failed reading 2a1e_maxit from rtdb',
     &     1, RTDB_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_prep_cf(nr, nt, l_cf, k_cf, l_eg, k_eg)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nr           ! [input]
      integer nt           ! [input]
      integer l_cf, k_cf   ! [output]
      integer l_eg, k_eg   ! [output]

      if (.not. ma_alloc_get(MT_DBL, (nr*nt)**2, '2a1e_coef',
     &     l_cf, k_cf))
     &     call errquit('WFAT: Cannot allocate coef', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nt, '2a1e_erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate erg', 1, MA_ERR)

      end
c===============================================================
      
      
c===============================================================
c===============================================================
      subroutine wfat_2a1e_erase_cf(l_cf, l_eg)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer l_cf, l_eg

      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: 2a1e_coef deallocation failed', 3,
     &     MA_ERR)
      if (.not. MA_free_heap(l_eg))
     &     call errquit('WFAT: 2a1e_erg deallocation failed', 3, MA_ERR)

      end
c===============================================================
      
