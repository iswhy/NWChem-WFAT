c=================================================
c=================================================
      logical function task_mewfat(rtdb)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "global.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "mewfat_qm_param.fh"
#include "wfat_field_param.fh"
      integer rtdb

      character(len=FLD_TP_LEN) las_env(max_nfield)
      character*4096 fname
      logical res_mode, td, rtdb_stat, os, oc, or, extf, isthere
      logical pr_srfc, pr_yild, pr_rate, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom, wfat1
      integer gqual, lmax, nlgl
      integer ndid, nqn, nlas, ntime, nsh, nintg
      integer nerg, ncpos, ncmm, niond
      integer l_npa, l_mpa, l_al, l_ftl, l_wl, l_pl, l_tl, l_dl
      integer l_eg, l_ps, l_cmm, l_idp, l_in
      integer g_il
      integer k_npa, k_mpa, k_al, k_ftl, k_wl, k_pl, k_tl, k_dl
      integer k_eg, k_ps, k_cmm, k_idp, k_in
      double precision tmax, fieldmag, smh
cne
      
      ndid = ga_nodeid()
      task_mewfat = .false.

      if (ndid==0) write(luout, '("WFAT: Moving to wfat directory")')
      
      
c     ==============================
c     ==== Initiate MEWFAT module ====
c     ==============================
      call mewfat_init(rtdb)


c     ==== Task title ====
      call ga_sync()
      call wfat_prep_tasktitle('ME')
      

c     ==== Reading general parameters ====
      call wfat_prep_general(rtdb, gqual, lmax)
      

c     ==== Checking the numerical density ====



c     ==== Obtaining parabolic quantum numbers ====
      call wfat_prep_paraqn(rtdb, nqn, l_npa, k_npa, l_mpa, k_mpa)


c     ==== Obtaining WFAT(1) information ====
      call wfat_prep_wfat1(rtdb, wfat1, nlgl)
      if (wfat1) call errquit('MEWFAT: At the moment, MEWFAT cannot '//
     &     'be used with the first order correction.', 62, UERR)
      

c     ==== Check if research mode is enabled ====
      call wfat_prep_resmode(rtdb, res_mode)
      

c     ==== Obtaining print switches ====
      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)
      if (pr_tmd .or. pr_tmd_atom)
     &     call errquit('MEWFAT: At the moment, TMD calculation is '//
     &     'not yet implemented within MEWFAT', 1, UERR)


c     ==== When TD field is true ====
      if (.not. rtdb_get(rtdb, 'wfat:tdfield', mt_log, 1, td))
     &     td = .false.
      if (td) then
         call wfat_prep_tdfield(rtdb, nlas, l_al, k_al, las_env,
     &        l_wl, k_wl, l_pl, k_pl, l_tl, k_tl, l_dl, k_dl, ntime,
     &        tmax, extf)
      else
         call wfat_prep_field(rtdb, fieldmag)
      endif


c     ==== Allocate channel properties ====
      call wfat_prep_channelprop(1, nerg, ncpos, ncmm, niond,
     &     l_eg, k_eg, l_ps, k_ps, l_cmm, k_cmm, l_idp, k_idp)
      

c      call wfat_test_dtrm(rtdb, geom_hd, bas_hd, nbas0, nspin0, nel, 
c     &                   g_mov, charge0, multip0, dbl_mb(k_moerg), tol2e, 
c     &                   oskel, ixcfit, iolgc)


c     ==== Print Dyson orbital ====
      call mewfat_print_dyson(rtdb, geom_hd, bas_hd, nel0, spinx,
     &     g_mov0, g_mov1, oskel)
      

c     ==== Calculate (or load) WFAT integrals ====
      nsh = (1 + lmax)**2
      nintg = nsh * natm * nqn
      if (.not. ga_create(mt_dcpl, nintg, 1, 'intgl', 0, 0, g_il))
     &     call errquit('WFAT: Cannot create intgl global array.',
     &     1, GA_ERR)
      call ga_zero(g_il)

      smh = 0.0d0
      call util_get_rtdb_state(rtdb, os, oc, or)

      if (.not. rtdb_cget(rtdb, 'wfat:intgl_name', 1, fname))
     &     call errquit('MEWFAT: Cannot read wfat:intgl_name from '//
     &     'RTDB.', 1, RTDB_ERR)
      inquire(file=fname, exist=isthere)      
      if (os .and. .not. or) then
cNEWRESTART      if (.not. isthere .or. os) then
c        ==== Calculate and store WFAT integrals on 'start' startup mode ====
         if (isthere) then
            if (ndid == 0)
     &           write(luout, '("WARNING: An integral file exist. '//
     &           'This file will be overwritten.")')
         endif
         
         call mewfat_store_intgl(rtdb, nsh, natm, nqn, int_mb(k_npa), 
     &        int_mb(k_mpa), lmax, gqual, smh, res_mode, dbl_mb(k_eg), 
     &        dbl_mb(k_ps), dbl_mb(k_cmm), dbl_mb(k_idp), g_il)
         call ga_sync()
      else if (.not. os .and. or) then
cNEWRESTART      else
c        ==== Binary-load MEWFAT integrals on 'restart' startup mode ====
         if (ndid == 0) then
            call mewfat_load_intgl(rtdb, nsh, natm, nqn, 1,
     &           spinx, dbl_mb(k_eg), dbl_mb(k_ps), dbl_mb(k_cmm),
     &           dbl_mb(k_idp), g_il)
         endif
         call ga_sync()
         
         call ga_brdcst(msg_wfat_eg, dbl_mb(k_eg),
     &                  ma_sizeof(MT_DBL, 1, MT_BYTE), 0)
         call ga_brdcst(msg_wfat_cmm, dbl_mb(k_cmm),
     &                  ma_sizeof(MT_DBL, 3*1, MT_BYTE), 0)
      else
         call errquit('MEWFAT: MEWFAT module requires the startup '//
     &        'directive to be either ''start'' or ''restart''.',
     &        1, UERR)
      endif


c     ==== Printings ====
      if (ndid == 0) then
         rtdb_stat = rtdb_parallel(.false.)
         if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &        call errquit('WFAT: Cannot allocate intgl.', 2, MA_ERR)
         call ga_get(g_il, 1, nintg, 1, 1, dcpl_mb(k_in), nintg)


c        ============================================
c        ==== Calculate and print key quantities ====
c        ============================================
         if (td) then
            if (pr_srfc .or. pr_yild .or. pr_avbeta .or. pr_avgamma)
     &           call mewfat_calc_td(rtdb, geom_hd, nel0, natm, pr_srfc,
     &           pr_yild, pr_avbeta, pr_avgamma, nlas, las_env, 
     &           dbl_mb(k_al), dbl_mb(k_wl), dbl_mb(k_pl), dbl_mb(k_tl),
     &           dbl_mb(k_dl), ntime, tmax, extf, 1,
     &           spinx, dbl_mb(k_eg), dbl_mb(k_cmm), nqn,
     &           int_mb(k_npa), int_mb(k_mpa), lmax, zion,
     &           dcpl_mb(k_in))
         else
            if (pr_srfc .or. pr_rate)
     &           call mewfat_calc_ti(rtdb, geom_hd, nel0, natm, pr_srfc,
     &           pr_rate, fieldmag, 1, spinx, dbl_mb(k_eg), 
     &           dbl_mb(k_cmm), nqn, int_mb(k_npa), int_mb(k_mpa),
     &           lmax, zion, dcpl_mb(k_in))
         endif
         
         
         if (.not. MA_free_heap(l_in))
     &        call errquit('WFAT: intgl deallocation failed', 3, MA_ERR)
         rtdb_stat = rtdb_parallel(.true.)
      endif
      

c     ==== Deallocations ====
      call wfat_erase_paraqn(l_npa, l_mpa)
      if (td) call wfat_erase_tdfield(l_al, l_wl, l_pl, l_tl, l_dl)
      call wfat_erase_channelprop(l_eg, l_ps, l_cmm, l_idp)
      if (.not. ga_destroy(g_il))
     &     call errquit('WFAT: Cannot destroy intgl.', 1, GA_ERR)


c     ==== Clean up and terminate task ====
      call wfat_close_task('ME')
      call mewfat_terminate(rtdb)
      

      task_mewfat = .true.
      
      
      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_calc_td(rtdb, geom, n0, nat, pr_srfc, pr_yild,
     &           pr_avbeta, pr_avgamma, nlas, ftlas, alas, wlas,
     &           plas, tlas, dlas, nt, tmax, extf, nch, spinid, erg,
     &           chmmt, nqn, nid, mid, lmax, zi, intgl)

c     DESCRIPTION:
      
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "wfat_field_param.fh"
      integer rtdb              ! [input]
      integer geom
      integer n0(2)             ! [input], no. of electrons in the neutral.
      integer nat
      logical pr_srfc                 ! [input], if .true. then print structure factor.
      logical pr_yild                 ! [input], if .true. then print ion yield.
      logical pr_avbeta               ! [input], if .true. then print beta-averaged ion yield.
      logical pr_avgamma              ! [input], if .true. then print gamma-averaged ion yield.
      integer nlas
      character(len=FLD_TP_LEN) ftlas(nlas)
      double precision alas(nlas)   ! Peak intensity, TW/cm^2                              
      double precision wlas(nlas)   ! Wavelength, nm                                       
      double precision plas(nlas)   ! Carrier-envelope phase, rad                          
      double precision tlas(nlas)   ! Delay (temporal position of the envelope peak), fs   
      double precision dlas(nlas)   ! Field FWHM, fs                                       
      double precision nt
      double precision tmax
      logical extf
      integer nch                     ! [input], number of ionized channels.
      integer spinid(nch)             ! [input], spin value (1 or 2) of ionized channels.
      double precision erg(nch)       ! [input]
      double precision chmmt(3, nch)  ! [input]
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      integer lmax                    ! [input]
      double precision zi             ! [input]
      double complex intgl(1:*)       ! [input]

      integer nsh, nintg, dummy_int, l_in, k_in
      double precision nall
      

      nsh = (1 + lmax)**2
      nintg = nsh * nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &     call errquit('MEWFAT: Cannot allocate intgl.', 2, MA_ERR)
      nall = dble(sum(n0))
      dcpl_mb(k_in : k_in+nintg-1) = sqrt(nall) * intgl(1:nintg)


      dummy_int = 0
c     ATTENTION: Update the arguments of wfat_calc_td below!!!
      call wfat_calc_td(rtdb, geom, pr_srfc, pr_yild, 
     &     pr_avbeta, pr_avgamma, nlas, ftlas, alas, wlas, 
     &     plas, tlas, dlas, nt, tmax, extf, nch, dummy_int, spinid, 
     &     erg, chmmt, nqn, nid, mid, lmax, zi, dcpl_mb(k_in))

      if (.not. ma_free_heap(l_in))
     &     call errquit('MEWFAT: Cannot deallocate intgl', 31, MA_ERR)
      
      end
c=================================================


c===============================================================
c===============================================================
      subroutine mewfat_calc_ti(rtdb, geom, n0, nat, pr_srfc, pr_rate, 
     &           fieldm, nch, spinid, erg, chmmt, nqn, nid, mid,
     &           lmax, zi, intgl)

c     DESCRIPTION:
      
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input]
      integer geom
      integer n0(2)             ! [input], no. of electrons in the neutral.
      integer nat
      logical pr_srfc                 ! [input], if .true. then print structure factor.
      logical pr_rate                 ! [input], if .true. then print ion yield.
      double precision fieldm
      integer nch                     ! [input], number of ionized channels.
      integer spinid(nch)             ! [input], spin value (1 or 2) of ionized channels.
      double precision erg(nch)       ! [input]
      double precision chmmt(3, nch)  ! [input]
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      integer lmax                    ! [input]
      double precision zi             ! [input]
      double complex intgl(1:*) ! [input]

      integer nsh, nintg, dummy_int, l_in, k_in
      double precision nall
      

      nsh = (1 + lmax)**2
      nintg = nsh * nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nintg, 'intgl', l_in, k_in))
     &     call errquit('MEWFAT: Cannot allocate intgl.', 2, MA_ERR)
      nall = dble(sum(n0))
      dcpl_mb(k_in : k_in+nintg-1) = sqrt(nall) * intgl(1:nintg)
      

      dummy_int = 0
      call wfat_calc_ti(rtdb, geom, pr_srfc, pr_rate, fieldm, 
     &     1, dummy_int, spinid, erg, chmmt, nqn, nid, mid, lmax,
     &     zi, dcpl_mb(k_in))


      if (.not. ma_free_heap(l_in))
     &     call errquit('MEWFAT: Cannot deallocate intgl', 83, MA_ERR)
      
      end
c===============================================================
