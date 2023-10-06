c===============================================================
c===============================================================
      subroutine wfat_2a1e_init(rtdb)

      implicit none
      


      call wfat_2a1e_geom(rtdb, geom_hd)



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

      character*255 tag(2)
      integer nc
      double precision ch(2), cd(3,2), r
      
      
      if (.not. geom_create(geom, 'wfat_geom'))
     &     call errquit('WFAT: Error in creating geometry object.', 71,
     &     GEOM_ERR)

      
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
      

      tag = (/ 'a1', 'a2' /)
      if (.not. geom_cart_set(geom, nc, tag(1:nc), cd(:,1:nc),
     &     ch(1:nc)))
     &     call errquit('WFAT: Error in initializing geometry object.',
     &     71, GEOM_ERR)

      end
c===============================================================
      
      
