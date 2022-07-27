c=======================================================
c=======================================================
      subroutine mewfat_print_info(rtdb)

      implicit none
#include "global.fh"      
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mewfat_qm_param.fh"
#include "geom.fh"
#include "bas.fh"
      integer rtdb

      logical rtdb_stat
      integer ndid
      

      ndid = ga_nodeid()
      if (ndid /= 0) return

      rtdb_stat = rtdb_parallel(.false.)


c     ==== WFAT module title ====
      write(luout, '(4/)')
      call util_print_centered(luout, 'Weak-Field Asymptotic Theory '//
     &     'Module', 40, .true.)

      
c     ==== WFAT level info ====
      call wfat_print_level(rtdb)


c     ==== Quantum system info ====
      call mewfat_print_qmsys()

      
c     ==== Geometry info ====
      write(luout, '(/)')
      if (.not. geom_print(geom_hd))
     &     call errquit('WFAT: Cannot print geometry info', 1, GEOM_ERR)
      

c     ==== Basis info ====
      write(luout, '(/)')
      if (.not. bas_print(bas_hd))
     &     call errquit('WFAT: Cannot print basis info', 1, BASIS_ERR)

      
c     ==== Printing info ====
      call wfat_print_printsw(rtdb)      
      

c     ==== Orientation angles ====
      call wfat_print_angles(rtdb)
      
            
c     ==== Remaining WFAT info ====
      call wfat_print_others(rtdb)

      
      rtdb_stat = rtdb_parallel(.true.)
      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine mewfat_print_qmsys()

      implicit none
#include "global.fh"
#include "mewfat_qm_param.fh"
#include "stdio.fh"
      integer ndid


      ndid = ga_nodeid()
      
      if (ndid == 0)
     &     write(luout, 900) 'neutral',
     &     trim(movecs0), trim(title0), trim(scftype0), trim(basname0),
     &     nbf0, ns_ion0, nmo0, nel0, multip0, charge0

      if (ndid == 0)
     &     write(luout, 900) 'cation',
     &     trim(movecs1), trim(title1), trim(scftype1), trim(basname1),
     &     nbf1, ns_ion1, nmo1, nel1, multip1, charge1

      if (ndid == 0) then
         write(luout, '(3x, "spinx = ", i3)') spinx
         write(luout, '(3x, "spiny = ", i3)') spiny
      endif
      
 900  format(/,
     &     3x, 'Info of the MO''s of ', a, ' :', /,
     &     6x, 'File name = ', a, /,
     &     6x, 'Job title = ', a, /,
     &     6x, 'SCF type = ', a, /,
     &     6x, 'Basis set name = ', a, /,
     &     6x, 'No. of bases = ', i6, /,
     &     6x, 'No. of spin channels = ', i6, /,
     &     6x, 'No. of MO''s in each channel = ', 2i6, /,
     &     6x, 'No. electrons in each channel = ', 2i6, /,
     &     6x, 'Spin multiplicity = ', i6, /,
     &     6x, 'Charge state = ', f7.3, /)
      
      end
c=======================================================
      
