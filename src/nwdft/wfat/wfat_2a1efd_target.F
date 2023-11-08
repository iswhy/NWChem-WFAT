c==========================================================
c==========================================================
      subroutine wfat_2a1efd_target(rtdb, geom, mch, nt, eps)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"      
#include "errquit.fh"
#include "rtdb.fh"
      integer rtdb
      integer geom
      integer mch
      integer nt
      double precision eps

      integer l_n, l_fn, l_cf, l_eg
      integer k_n, k_fn, k_cf, k_eg
      integer ni, nr


      write(luout, '("herex0")')
      if (.not. rtdb_get(rtdb, 'wfat:2a1efd_nfe', mt_int, 1, ni))
     &     call errquit('WFAT: Failed reading 2a1e_nfe from rtdb',
     &     1, RTDB_ERR)

      write(luout, '("herex1")')
      if (.not. ma_alloc_get(MT_INT, ni, 'ndvr', l_n, k_n))
     &     call errquit('WFAT: Cannot allocate ndvr', 3, MA_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:2a1efd_ndvr', mt_int, ni, 
     &     int_mb(k_n))) call errquit('WFAT: Failed reading 2a1e_ndvr'//
     &     ' from rtdb', 1, RTDB_ERR)

      write(luout, '("herex2")')
      if (.not. ma_alloc_get(MT_DBL, ni, 'fenodes', l_fn, k_fn))
     &     call errquit('WFAT: Cannot allocate fenodes', 3, MA_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:2a1efd_fenodes', mt_dbl, ni, 
     &     dbl_mb(k_fn))) call errquit('WFAT: Failed reading '//
     &     '2a1e_fenodes from rtdb', 1, RTDB_ERR)
      

      write(luout, '("herex3")')
      nr = (sum(int_mb(k_n:k_n+ni-1)) - ni*2) + (ni - 1)
      if (.not. ma_alloc_get(MT_DBL, (nr*nt)**2, 'coeff', l_cf, k_cf))
     &     call errquit('WFAT: Cannot allocate coeff', 3, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nt, 'erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate erg', 3, MA_ERR)


      write(luout, '("herex4")')
      call wfat_2a1efd_solve(mch, ni, int_mb(k_n), dbl_mb(k_fn), nt, 
     &     geom, eps, dbl_mb(k_cf), dbl_mb(k_eg))
      

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_n))
     &     call errquit('WFAT: Cannot deallocate ndvr', 1, MA_ERR)
      if (.not. ma_free_heap(l_fn))
     &     call errquit('WFAT: Cannot deallocate fenodes', 1, MA_ERR)
      if (.not. ma_free_heap(l_cf))
     &     call errquit('WFAT: Cannot deallocate coeff', 1, MA_ERR)
      if (.not. ma_free_heap(l_eg))
     &     call errquit('WFAT: Cannot deallocate erg', 1, MA_ERR)

      end
c==========================================================
      
