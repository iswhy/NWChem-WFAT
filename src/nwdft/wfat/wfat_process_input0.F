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
#include "stdio.fh"
#include "mafdecls.fh"
      integer rtdb

      logical r
      integer ozn, gs
      double precision thr, b2a

      parameter (b2a=0.52917721067d0)


c     ==== gridset ====
      if (.not. rtdb_get(rtdb, 'wfat:inp:gqual', mt_int, 1, gs))
     &     gs = def_gqual
      call wfat_ival_inrange(gs, 1, 4, 'WFAT', '''gridset''', 'i', r)
      if (.not. rtdb_put(rtdb, 'wfat:gqual', mt_int, 1, gs))
     &     call errquit('WFAT: Failed writing gqual to RTDB',
     &     1, RTDB_ERR)

      
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
c     DEF_CGRID_THR and wfat:inp:cgrid_thr are in angstrom.
      call wfat_dval_atleast(thr, 0.0d0, 'g', 'WFAT', '''cgridthr''',
     &     'e', 'i', r)
      if (.not. rtdb_put(rtdb, 'wfat:cgrid_thr', mt_dbl, 1, thr/b2a))
     &     call errquit('WFAT: Failed writing cgrid_thr to RTDB',
     &     1, RTDB_ERR)      

      end
c=======================================================
      
