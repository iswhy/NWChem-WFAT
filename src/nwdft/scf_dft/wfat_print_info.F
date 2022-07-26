c=======================================================
c=======================================================
      subroutine wfat_print_info(rtdb)

      implicit none
#include "global.fh"      
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "wfat_qm_param.fh"
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
      call wfat_print_qmsys(rtdb)
            
      
c     ==== Geometry info ====
      write(luout, '(/)')
      if (.not. geom_print(geom_hd))
     &     call errquit('WFAT: Cannot print geometry info', 1, GEOM_ERR)
      

c     ==== Basis info ====
      write(luout, '(/)')
      if (.not. bas_print(bas_hd))
     &     call errquit('WFAT: Cannot print basis info', 1, BASIS_ERR)

      
c     ==== Channel info ====
      call wfat_print_channel(rtdb, nel)
      

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
      subroutine wfat_print_level(rtdb)

      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb
      
      character*2 wlevel

      if (.not. rtdb_cget(rtdb,'wfat:level', 1, wlevel)) 
     &     call errquit('WFAT: Failed reading wfat level',0, RTDB_ERR)
      write(luout, '(2/)')
      call util_print_centered(luout, 'Level of calculation',
     &     40, .true.)
      if (wlevel == 'oe') then
         call util_print_centered(luout, 'One-electron WFAT',
     &        40, .false.)
      else if (wlevel == 'me') then
         call util_print_centered(luout, 'Many-electron WFAT',
     &        40, .false.)
      else
         call errquit('WFAT: WFAT level unknown.', 1, UERR)
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_qmsys(rtdb)

      implicit none
#include "wfat_qm_param.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb

      character*255 fname
      

      write(luout, '(2/)')
      call util_print_centered(luout, 'Previous SCF info', 20, .true.)
      write(luout, '(1x)')
      write(luout, '(5x, a22, " = ", i4)') "no. of spins", nspin0
      write(luout, '(5x, a22, " = ", i4)') "alpha electrons", nel(1)
      if (nspin0 == 1) then
         write(luout, '(5x, a22, " = ", i4)') "beta electrons", nel(1)
      else
         write(luout, '(5x, a22, " = ", i4)') "beta electrons", nel(2)
      endif
      write(luout, '(5x, a22, " = ", f9.4)') "post-ionization charge",
     &     zion
      if (.not. rtdb_cget(rtdb, 'wfat:movec_name', 1, fname))
     &     call errquit('WFAT: Cannot read wfat:movec_name from RTDB.',
     &     1, RTDB_ERR)
      write(luout, '(5x, a22, " = ", a)')
     &     "Molecular orbital file", trim(adjustl(fname))

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_channel(rtdb, n0)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      integer n0(2)

      character*5 spinnm(2), i_str
      character*7 io_str
      character*12 hl
      integer l_sp, l_or0, l_or
      integer k_sp, k_or0, k_or
      integer nch, i, is, io

      
      spinnm = (/ 'alpha', ' beta' /)      

      write(luout, '(2/)')
      call util_print_centered(luout, 'Ionized Channels', 20, .true.)
      write(luout, '(1x)')

      call wfat_prep_channels(rtdb, n0, nch, l_or0, k_or0, l_sp, k_sp,
     &     l_or, k_or)
      write(luout, '(5x, "No. of channels = ", i6)') nch

      do i = 1, nch
         write(i_str, '(i5)') i
         is = int_mb(k_sp + i - 1)
         io = int_mb(k_or0 + i - 1)
         if (io == 0) then
            hl = 'HOMO'
         else if (io < 0) then
            write(io_str, '(i7)') abs(io)
            hl = "HOMO-"//adjustl(io_str)
         else if (io == 1) then
            hl = 'LUMO'
         else if (io > 1) then
            write(io_str, '(i7)') abs(io-1)
            hl = "LUMO+"//adjustl(io_str)
         endif
         write(luout, '(5x, "Channel ", a5, " : ", a5, "/", a12)')
     &     adjustl(i_str), spinnm(is), adjustl(hl) 
      enddo

      call wfat_erase_channels(l_or0, l_sp, l_or)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_printsw(rtdb)

      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      
      character*3 ispr      
      logical pr_srfc, pr_rate, pr_yild, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom
      

      write(luout, '(2/)')
      call util_print_centered(luout, 'Printing Information', 20,
     &     .true.)
      write(luout, '(1x)')

      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)

      call wfat_tf_to_yn(pr_srfc, .true., ispr)
      write(luout, '(5x, "Print structure factor : ", a3)') ispr
      call wfat_tf_to_yn(pr_rate, .true., ispr)
      write(luout, '(5x, "Print ion rate : ", a3)') ispr
      call wfat_tf_to_yn(pr_yild, .true., ispr)
      write(luout, '(5x, "Print ion yield : ", a3)') ispr
      call wfat_tf_to_yn(pr_avgamma, .true., ispr)
      write(luout, '(5x, "Print gamma-averaged ion yield : ", a3)') ispr
      call wfat_tf_to_yn(pr_avbeta, .true., ispr)
      write(luout, '(5x, "Print beta-averaged ion yield : ", a3)') ispr
      call wfat_tf_to_yn(pr_tmd, .true., ispr)
      write(luout, '(5x, "Print transverse momentum distribution '//
     &     '(TMD) : ", a3)') ispr
      call wfat_tf_to_yn(pr_tmd_atom, .true., ispr)
      write(luout, '(5x, "Print atom-resolved TMD : ", a3)') ispr

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_angles(rtdb)

      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      
      logical pr_srfc, pr_rate, pr_yild, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom
      integer nbt, ngm
      double precision bt_start, bt_last, gm_start, gm_last
      double precision bt_tmd, gm_tmd
      

      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)


      write(luout, '(2/)')
      call util_print_centered(luout, 'Orientation angle '//
     &     'information', 20, .true.)
      write(luout, '(1x)')
      if (pr_srfc .or. pr_rate) then
c        ==== Beta ====
         if (.not. rtdb_get(rtdb, 'wfat:beta_start', mt_dbl, 1,
     &        bt_start))
     &        call errquit('WFAT: Failed reading beta_start from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "First beta (deg) ", bt_start
         if (.not. rtdb_get(rtdb, 'wfat:beta_last', mt_dbl, 1,
     &        bt_last))
     &        call errquit('WFAT: Failed reading beta_last from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "Last beta (deg) ", bt_last
         if (.not. rtdb_get(rtdb, 'wfat:nbeta', mt_int, 1, nbt))
     &        call errquit('WFAT: Failed reading nbeta from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", i7)')
     &        "Number of beta points  ", nbt

c        ==== Gamma ====
         if (.not. rtdb_get(rtdb, 'wfat:gamma_start', mt_dbl, 1,
     &        gm_start))
     &        call errquit('WFAT: Failed reading gamma_start from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "First gamma (deg) ", gm_start
         if (.not. rtdb_get(rtdb, 'wfat:gamma_last', mt_dbl, 1,
     &        gm_last))
     &        call errquit('WFAT: Failed reading gamma_last from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "Last gamma (deg) ", gm_last
         if (.not. rtdb_get(rtdb, 'wfat:ngamma', mt_int, 1, ngm))
     &        call errquit('WFAT: Failed reading ngamma from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", i7)')
     &        "Number of gamma points ", ngm
      endif
      
      if (pr_tmd .or. pr_tmd_atom) then         
         if (.not. rtdb_get(rtdb, 'wfat:beta_tmd', mt_dbl, 1, bt_tmd))
     &        call errquit('WFAT: Failed reading bt_tmd from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "Beta (deg) for TMD ", bt_tmd
         if (.not. rtdb_get(rtdb, 'wfat:gamma_tmd', mt_dbl, 1, gm_tmd))
     &        call errquit('WFAT: Failed reading gm_tmd from rtdb',
     &        1, INPUT_ERR)
         write(luout, '(5x, a23, " = ", f10.6)')
     &        "gamma (deg) for TMD ", gm_tmd
      endif

      write(luout, '(/, 5x, "Note :")')
      write(luout, 601)
 601  format(7x,
     &     "To get an idea of how the molecule is oriented with "     ,
     &     "the above values of orientation ",   /,
     &     7x, 
     &     "angles, imagine the following rotations. First imagine "  ,
     &     "a pair of coinciding axes, the ", /,
     &     7x, 
     &     "(space-fixed) lab axes and the (body-fixed) internal "    ,
     &     "axes. Then rotate the molecule around ", /,
     &     7x, 
     &     "the lab z axis (or internal z, whichever is fine at "     ,
     &     "this point) counterclockwise as seen ",        /,
     &     7x, 
     &     "from lab +z axis by an ammount gamma. Thus, at this "  ,
     &     "point, gamma is the angle subtended ",  /,
     &     7x, 
     &     "by the lab and internal y axes. Next, rotate the ",
     &     "molecule around the lab y axis counter ",/,
     &     7x, 
     &     "clockwise as seen from lab +y axis by an ammount beta ",
     &     "so that beta is the angle between ", /,
     &     7x,
     &     "the lab and internal z axes. These rotations yield the ",
     &     "molecule's final orientation.")

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_others(rtdb)

      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      
      character*255 fname
      logical pr_srfc, pr_rate, pr_yild, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom, autopos
      integer nktmd, lmax, gqual
      double precision fieldmag, omgpos(3)
      

      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)


      write(luout, '(2/)')
      call util_print_centered(luout, 'Other WFAT inputs', 20, .true.)
      write(luout, '(1x)')

cOLDC     ==== Field magnitude ====
cOLD      if (pr_rate .or. pr_tmd .or. pr_tmd_atom) then
cOLD         if (.not. rtdb_get(rtdb, 'wfat:fieldmag', mt_dbl, 1, fieldmag))
cOLD     &        call errquit('WFAT: Failed reading fieldmag from rtdb',
cOLD     &        1, INPUT_ERR)
cOLD         write(luout, '(5x, a28, " = ", es10.3, " a.u.")')
cOLD     &        "Field strength ", fieldmag
cOLD      endif

c     ==== No. of points for TMD ====
      if (pr_tmd .or. pr_tmd_atom) then
         if (.not. rtdb_get(rtdb, 'wfat:nktmd', mt_int, 1, nktmd))
     &        call errquit('WFAT: Failed reading nktmd from rtdb',
     &        1, RTDB_ERR)
         write(luout, '(5x, a28, " = ", i7)')
     &        "No. of x momentum points ", nktmd
         write(luout, '(5x, a28, " = ", i7)')
     &        "No. of y momentum points ", nktmd
      endif

c     ==== Max. angular momentum ====
      if (.not. rtdb_get(rtdb, 'wfat:lmax', mt_int, 1, lmax))
     &     call errquit('WFAT: Failed reading lmax from rtdb',
     &     1, INPUT_ERR)
      write(luout, '(5x, a28, " = ", i3)')
     &     "Max. angular momentum ", lmax

c     ==== Grid quality ====
      if (.not. rtdb_get(rtdb, 'wfat:gqual', mt_int, 1, gqual))
     &     call errquit('WFAT: Failed reading gqual from rtdb',
     &     1, INPUT_ERR)
      write(luout, '(5x, a28, " = ", i3)') "Grid set ", gqual

c     ==== Integral file ====
      if (.not. rtdb_cget(rtdb, 'wfat:intgl_name', 1, fname))
     &     call errquit('WFAT: Cannot read wfat:intgl_name from RTDB',
     &     1, RTDB_ERR)
      write(luout, '(5x, a28, " = ", a)')
     &     "Integral file ", trim(adjustl(fname))

c     ==== Auto/manual positioning of origin ====
      if (.not. rtdb_get(rtdb, 'wfat:autopos', mt_log, 1, autopos))
     &     call errquit('WFAT: Cannot read wfat:autopos from RTDB.',
     &     1, RTDB_ERR)
      if (autopos) then
         write(luout, '(5x, a28, " = ", a)') "Channel pos. ", 'auto'
      else
         if (.not. rtdb_get(rtdb, 'wfat:omgpos', mt_dbl, 3, omgpos))
     &     call errquit('WFAT: Cannot read wfat:omgpos from RTDB.',
     &     1, RTDB_ERR)
         write(luout, '(5x, a28, " = ", 3es14.5)')
     &        "Channel pos. ", omgpos(1:3)

         write(luout, '(/, 7x, "*** ATTENTION ***")')
         write(luout, 602)
      endif
      
 602  format(9x,
     &     "Manual specification for the input 'omgpos' such "      ,
     &     "as in the present run should be used with caution. ",    /,
     &     9x, 
     &     "This feature is only intended for research purposes. "    ,
     &     "Improper choice of the values for 'omgpos' can lead ", /,
     &     9x, 
     &     "to a behavior of the results that is not in agreement, "  ,
     &     "even qualitatively, with results from other methods. ",  /,
     &     9x, 
     &     "To let NWChem pick the right values for the user, "       ,
     &     "use 'auto' as the value for 'omgpos'. ")

      end
c=======================================================
      
