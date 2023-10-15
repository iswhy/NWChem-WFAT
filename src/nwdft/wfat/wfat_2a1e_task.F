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

      character*2 td
      logical pr_srfc, pr_yild, pr_rate, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom
      logical res_mode, wfat1, rtdb_stat
      integer l_cf, l_eg, l_npa, l_mpa, l_in, l_dpo, l_inj, l_ink1
      integer l_ink2, l_pb
      integer k_cf, k_eg, k_npa, k_mpa, k_in, k_dpo, k_inj, k_ink1
      integer k_ink2, k_pb
      integer mch, sphi, spinid, orbid, nr, nt, maxiter, i_cf, i_eg
      integer lmax, gqual, nqn, nsh, nintg, g_il, ndid, nlgl, ndpo
      integer nintj, nintk, g_j, g_k(2), npb
      double precision eps, init_ke, out_ke, edip(3), chpos(3), chmmt(3)
      double precision fieldmag, pb(3,3)

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
      call wfat_2a1e_prep_system(rtdb, mch, sphi, spinid, orbid, eps, 
     &     nr, nt, init_ke, maxiter)


c     ==== Allocate the coefficients and energies array ====
      call wfat_2a1e_prep_cf(nr, nt, l_cf, k_cf, l_eg, k_eg)


c     ==== Obtaining parabolic quantum numbers ====
      call wfat_prep_paraqn(rtdb, nqn, l_npa, k_npa, l_mpa, k_mpa)


c     ==== Obtaining WFAT(1) information ====
      call wfat_prep_wfat1(rtdb, wfat1, nlgl)
      

c     ==== Obtaining print switches ====
      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)


c     ==== Whether using TD field or no ====
      if (.not. rtdb_cget(rtdb, 'wfat:tdfield', 1, td)) td = 'no'
      if (td=='td') then
clater         call wfat_prep_tdfield(rtdb, nlas, l_al, k_al, las_env,
clater     &        l_wl, k_wl, l_pl, k_pl, l_tl, k_tl, l_dl, k_dl, ntime,
clater     &        tmax, extf)
      else if (td=='ti') then
         call wfat_prep_field(rtdb, fieldmag)
      else if (td=='no') then
         fieldmag = 0.0d0
      endif
      
      
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
      

c     ==== Preparing WFAT(1) integrals ====
      call wfat_prep_wfat1_intgl(wfat1, lmax, natm, nqn, nch, nintj, 
     &     g_j, nintk, g_k, npb, l_pb, k_pb)

      
c     ==== Compute the OEWFAT integrals ====
      i_cf = k_cf + (orbid-1)*nr*nt
      i_eg = k_eg + (orbid-1)
      call wfat_2a1e_intgl(rtdb, geom_hd, nqn, int_mb(k_npa), 
     &     int_mb(k_mpa), lmax, gqual, zion, wf_msm, sphi, mch, nr, nt, 
     &     out_ke, orbid, dbl_mb(i_cf), dbl_mb(i_eg), edip, eps, wfat1, 
     &     chpos, chmmt, g_il)
      

c     ==== Ionization calculations and printings ====
      if (ndid == 0) then
         rtdb_stat = rtdb_parallel(.false.)
         if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &        call errquit('WFAT: Cannot allocate intgl.', 41, MA_ERR)
         call ga_get(g_il, 1, nintg, 1, 1, dcpl_mb(k_in), nintg)


c        ==== When WFAT(1) is in use ====
         call wfat_prep_wfat1_ti(wfat1, nr*nt, 2, nintj, nintk,
     &        ndpo, l_dpo, k_dpo, l_inj, k_inj, l_ink1, k_ink1,
     &        l_ink2, k_ink2)
clater         if (wfat1) then
clater            call ga_get(g_j, 1, nintj, 1, 1, dcpl_mb(k_inj), nintj)
clater            call ga_get(g_k(1), 1, nintk, 1, 1, dcpl_mb(k_ink1), nintk)
clater            call ga_get(g_k(2), 1, nintk, 1, 1, dcpl_mb(k_ink2), nintk)
clater            ctr = 0.0d0
clater            call wfat_calc_dpmo(ctr, dbl_mb(k_dpo))
clater         endif
         

c        ==== Printing the WFAT integrals ====
         if (wfat1) then
         else
            call wfat_print_wfat0_ints(nch, nqn, natm, lmax, 
     &           int_mb(k_npa), int_mb(k_mpa), dcpl_mb(k_in))
         endif


c        ==================================================
c        ==== Calculate and print important quantities ====
c        ==================================================
         if (td == 'td') then
         else
            if (pr_srfc .or. pr_rate)
     &           call wfat_calc_ti(rtdb, geom_hd, pr_srfc, pr_rate,
     &           fieldmag, nch, orbid, spinid, dbl_mb(i_eg), chmmt, 
     &           nqn, int_mb(k_npa), int_mb(k_mpa), lmax, zion, 
     &           dcpl_mb(k_in), wfat1, nr*nt, 2, orbid, dbl_mb(k_dpo), 
     &           pb, chpos, dcpl_mb(k_inj), dcpl_mb(k_ink1), 
     &           dcpl_mb(k_ink2))

clater            if (pr_tmd .or. pr_tmd_atom) 
clater     &           call wfat_calcprint_tmd(rtdb, geom_hd, pr_tmd, 
clater     &           pr_tmd_atom, fieldmag, nch, int_mb(k_or0), 
clater     &           int_mb(k_sp), dbl_mb(k_eg), dbl_mb(k_cmm), nqn, 
clater     &           int_mb(k_npa), int_mb(k_mpa), lmax, dcpl_mb(k_in),
clater     &           res_mode)
         endif

         
         if (.not. MA_free_heap(l_in))
     &        call errquit('WFAT: intgl deallocation failed', 3, MA_ERR)
         rtdb_stat = rtdb_parallel(.true.)
      endif
      call ga_sync()
      
      
c     ==== Deallocations ====
      call wfat_erase_paraqn(l_npa, l_mpa)
      call wfat_2a1e_erase_cf(l_cf, l_eg)
      if (.not. ga_destroy(g_il))
     &     call errquit('WFAT: Cannot destroy intgl.', 45, GA_ERR)


      task_wfat_2a1e = .true.
      call wfat_close_task('OE-2a1e')
      
      end
c===============================================================
      
