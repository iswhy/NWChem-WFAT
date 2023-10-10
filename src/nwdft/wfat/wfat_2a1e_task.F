c===============================================================
c===============================================================
      logical function task_wfat_2a1e(rtdb)

      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "wfat_qm_param.fh"
#include "wfat.fh"
#include "errquit.fh"
      integer rtdb

      logical res_mode, wfat1, rtdb_stat
      integer l_cf, l_eg, l_npa, l_mpa, l_in
      integer k_cf, k_eg, k_npa, k_mpa, k_in
      integer mch, sphi, orbid, nr, nt, maxiter, i_cf, i_eg, lmax, gqual
      integer nqn, nsh, nintg, g_il, ndid
      double precision eps, init_ke, out_ke, edip(3), chpos(3), chmmt(3)

      integer nch
      parameter(nch=1)


      ndid = ga_nodeid()
      

c     ==============================
c     ==== Initiate WFAT module ====
c     ==============================
      call wfat_2a1e_init(rtdb)


      res_mode = .false.

c     ==== Task title ====
      call ga_sync()
      call wfat_prep_tasktitle('OE-2a1e')


c     ==== Extract the general parameters ====
      call wfat_prep_general(rtdb, gqual, lmax)

      
c     ==== Extract the parameters of the quantum system ====
      call wfat_2a1e_prep_system(rtdb, mch, sphi, orbid, eps, nr, nt, 
     &     init_ke, maxiter)


c     ==== Allocate the coefficients and energies array ====
      call wfat_2a1e_prep_cf(nr, nt, l_cf, k_cf, l_eg, k_eg)


c     ==== Obtaining parabolic quantum numbers ====
      call wfat_prep_paraqn(rtdb, nqn, l_npa, k_npa, l_mpa, k_mpa)
      

c     ==== Solve the Hamiltonian for the target ====
c     ====   eigenstate (the ionized channel)   ====
      call wfat_2a1e_target(geom_hd, mch, sphi, orbid, eps, nr, nt,
     &     init_ke, maxiter, dbl_mb(k_cf), dbl_mb(k_eg), out_ke)


c     ==== Calculate (or load) WFAT integrals ====
      nsh = (1 + lmax)**2
      nintg = nsh * natm * nqn * nch
      if (.not. ga_create(mt_dcpl, nintg, 1, 'intgl', 0, 0, g_il))
     &     call errquit('WFAT: Cannot create intgl global array.',
     &     1, GA_ERR)
      call ga_zero(g_il)
      

c     ==== Compute the OEWFAT integrals ====
      i_cf = k_cf + (orbid-1)*nr*nt
      i_eg = k_eg + (orbid-1)
      call wfat_2a1e_intgl(rtdb, geom_hd, nqn, int_mb(k_npa), 
     &     int_mb(k_mpa), lmax, gqual, zion, wf_msm, mch, nr, nt, 
     &     out_ke, orbid, dbl_mb(i_cf), dbl_mb(i_eg), edip, eps, wfat1, 
     &     chpos, chmmt, g_il)
      

c     ==== Ionization calculations and printings ====
      if (ndid == 0) then
         rtdb_stat = rtdb_parallel(.false.)
         if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &        call errquit('WFAT: Cannot allocate intgl.', 41, MA_ERR)
         call ga_get(g_il, 1, nintg, 1, 1, dcpl_mb(k_in), nintg)
         

         if (wfat1) then
         else
            call wfat_print_wfat0_ints(nch, nqn, natm, lmax, 
     &           int_mb(k_npa), int_mb(k_mpa), dcpl_mb(k_in))
         endif

         if (.not. MA_free_heap(l_in))
     &        call errquit('WFAT: intgl deallocation failed', 3, MA_ERR)
         rtdb_stat = rtdb_parallel(.true.)
      endif
      
      
c     ==== Deallocations ====
      call wfat_erase_paraqn(l_npa, l_mpa)
      call wfat_2a1e_erase_cf(l_cf, l_eg)
      if (.not. ga_destroy(g_il))
     &     call errquit('WFAT: Cannot destroy intgl.', 45, GA_ERR)


      task_wfat_2a1e = .true.
      call wfat_close_task('OE-2a1e')
      
      end
c===============================================================
      
