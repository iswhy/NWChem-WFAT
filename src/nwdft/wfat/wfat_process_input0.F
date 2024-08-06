c=======================================================
c=======================================================
      subroutine wfat_process_input0(rtdb)

      implicit none
      integer rtdb

      call wfat_process_quadgrid(rtdb)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_process_quadgrid(rtdb)

      implicit none
#include "wfat_const.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb

      logical r
      integer ozn
      double precision thr

      
c     ==== omegaz ====
      if (.not. rtdb_get(rtdb, 'wfat:inp:ozn', mt_int, 1, ozn))
     &     ozn = DEF_OZN
      call wfat_ival_atleast(ozn, 1, 'WFAT', '''omegaz''', 'i', r)
      if (.not. rtdb_put(rtdb, 'wfat:ozn', mt_int, 1, ozn))
     &     call errquit('WFAT: Failed writing ozn to RTDB',
     &     1, RTDB_ERR)


c     ==== cgridthr ====
      if (.not. rtdb_get(rtdb, 'wfat:inp:cgrid_thr', mt_dbl, 1, thr))
     &     thr = DEF_CGRID_THR
      call wfat_dval_atleast(thr, 0.0d0, 'WFAT', '''cgridthr''', 'i', r)
      if (.not. rtdb_put(rtdb, 'wfat:cgrid_thr', mt_dbl, 1, thr))
     &     call errquit('WFAT: Failed writing cgrid_thr to RTDB',
     &     1, RTDB_ERR)      
      
      
      end
c=======================================================
      
