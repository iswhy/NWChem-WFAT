c===============================================================
c===============================================================
      subroutine wfat_2a1e_init(rtdb)

      implicit none
#include "geom.fh"
#include "errquit.fh"
#include "wfat_qm_param.fh"
      integer rtdb

      call wfat_2a1e_geom(rtdb, geom_hd)
      if (.not. geom_ncent(geom_hd, natm))
     &     call errquit('WFAT: Getting natm failed', 89, GEOM_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_geom(rtdb, geom)

      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "geom.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input]
      integer geom              ! [output]

      character*16 tag(2)
      integer nc
      double precision ch(2), cd(3,2), r
      

c     ==== Create the geometry object for 2a1e ====
      if (.not. geom_create(geom, 'wfat_geom'))
     &     call errquit('WFAT: Error in creating geometry object.', 71,
     &     GEOM_ERR)

      
c     ==== Read some quantities from RDB ====
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_nc', mt_int, 1, nc))
     &     call errquit('WFAT: Failed reading 2a1e_nc from RTDB.', 1,
     &     RTDB_ERR)
      if (nc <= 0 .or. nc >=3)
     &     call errquit('WFAT: The number of 2a1e centers must be '//
     &     'either 1 or 2.', 93, UERR)
      if (.not. rtdb_get(rtdb, 'wfat:2a1e_q', mt_dbl, nc, ch(1:nc)))
     &     call errquit('WFAT: Failed reading 2a1e_q from RTDB.',
     &     1, RTDB_ERR)
      if (nc == 2) then
         if (.not. rtdb_get(rtdb, 'wfat:internuc', mt_dbl, 1, r))
     &        call errquit('WFAT: Failed reading internuc from RTDB.',
     &        1, RTDB_ERR)
         cd = 0.0d0
         cd(3,1) = ch(2)/sum(ch) * r
         cd(3,2) = -ch(1)/sum(ch) * r
      else if (nc == 1) then
         cd = 0.0d0
      endif

      
c     ==== Set the attributes of 2a1e geom using ====
c     ====       the above RTDB quantities       ====
      tag = (/ 'a1', 'a2' /)
      if (.not. geom_cart_set(geom, nc, tag(1:nc), cd(:,1:nc),
     &     ch(1:nc)))
     &     call errquit('WFAT: Error in initializing geometry object.',
     &     71, GEOM_ERR)

      end
c===============================================================
      
      
