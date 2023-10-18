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
     &     call errquit('WFAT: Cannot allocate 2a1e_coef', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nt, '2a1e_erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate 2a1e_erg', 1, MA_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_prep_cf1(wfat1, rtdb, m, nr, nt, nr_m, nt_m, 
     &     nr_p, nt_p, l_cfx, k_cfx, l_cfy, k_cfy, l_cfz, k_cfz)

      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      logical wfat1
      integer rtdb
      integer m               
      integer nr, nt
      integer nr_m, nt_m
      integer nr_p, nt_p
      integer l_cfx, k_cfx
      integer l_cfy, k_cfy
      integer l_cfz, k_cfz    

      integer n_m, n_p, nz


      if (wfat1) then
c        ==== The number of radial DVR basis for m-1 spectrum ====
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_nr_m', mt_int, 1, nr_m))
     &        call errquit('WFAT: Failed reading 2a1e_nr_m from rtdb',
     &        1, RTDB_ERR)
         
c        ==== The number of polar DVR basis for m-1 spectrum ====
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_nt_m', mt_int, 1, nt_m))
     &        call errquit('WFAT: Failed reading 2a1e_nt_m from rtdb',
     &        1, RTDB_ERR)
         
c        ==== The number of radial DVR basis for m+1 spectrum ====
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_nr_p', mt_int, 1, nr_p))
     &        call errquit('WFAT: Failed reading 2a1e_nr_p from rtdb',
     &        1, RTDB_ERR)
         
c        ==== The number of polar DVR basis for m+1 spectrum ====
         if (.not. rtdb_get(rtdb, 'wfat:2a1e_nt_p', mt_int, 1, nt_p))
     &        call errquit('WFAT: Failed reading 2a1e_nt_p from rtdb',
     &        1, RTDB_ERR)
      else
         nr_m = 1
         nt_m = 1
         nr_p = 1
         nt_p = 1
      endif
      
      
c     ==== Allocations ====
      if (wfat1) then
         if (m >= 1) then
            n_m = nr_m * nt_m
         else
            n_m = 0
         endif
         n_p = nr_p * nt_p
         nz = nr * nt
      else
         n_m = 1
         n_p = 1
         nz = 1
      endif

      if (.not. ma_alloc_get(MT_DBL, n_m+n_p, '2a1e_coefx',
     &     l_cfx, k_cfx))
     &     call errquit('WFAT: Cannot allocate 2a1e_coefx', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, n_m+n_p, '2a1e_coefy',
     &     l_cfy, k_cfy))
     &     call errquit('WFAT: Cannot allocate 2a1e_coefy', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nz, '2a1e_coefz',
     &     l_cfz, k_cfz))
     &     call errquit('WFAT: Cannot allocate 2a1e_coefz', 1, MA_ERR)

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


c===============================================================
c===============================================================
      subroutine wfat_2a1e_erase_cf1(l_cfx, l_cfy, l_cfz)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer l_cfx, l_cfy, l_cfz

      if (.not. MA_free_heap(l_cfx))
     &     call errquit('WFAT: 2a1e_coefx deallocation failed', 3,
     &     MA_ERR)
      if (.not. MA_free_heap(l_cfy))
     &     call errquit('WFAT: 2a1e_coefy deallocation failed', 3,
     &     MA_ERR)
      if (.not. MA_free_heap(l_cfz))
     &     call errquit('WFAT: 2a1e_coefz deallocation failed', 3,
     &     MA_ERR)

      end
c===============================================================
      
      
