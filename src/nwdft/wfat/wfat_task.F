c===============================================================
c===============================================================
      logical function task_wfat(rtdb)

      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "global.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "wfat_qm_param.fh"
#include "wfat_field_param.fh"
      integer rtdb

      character(len=FLD_TP_LEN) las_env(max_nfield)
      logical pr_srfc, pr_yild, pr_rate, pr_avbeta, pr_avgamma
      logical pr_tmd, pr_tmd_atom
      logical rtdb_stat, wfat1, res_mode, os, oc, or, int_nosave
      logical dd, td, extf
      character*32 theory
      character*2 prtb_type
      integer ndid, gqual, i, lmax, nlgl, is, nhp
      integer l_nd, l_or0, l_sp, l_or, l_npa, l_mpa, l_eg
      integer l_ps, l_cmm, l_idp, l_in, l_ftl, l_al, l_wl, l_pl, l_tl
      integer l_dl, l_pb
      integer k_nd, k_or0, k_sp, k_or, k_npa, k_mpa, k_eg
      integer k_ps, k_cmm, k_idp, k_in, k_ftl, k_al, k_wl, k_pl, k_tl
      integer k_dl, k_pb
      integer g_il, g_j, g_k(2), nlas, ntime
      integer nch, nqn, nerg, ncpos, ncmm, niond, nsh, nintg, nintj
      integer nintk, prtb_ord
      double precision fieldmag, smh, prtb_bt, prtb_gm
      double precision tmax

      logical task_mewfat
      external task_mewfat

      

      
      ndid = ga_nodeid()

      
c     ==============================
c     ==== Initiate WFAT module ====
c     ==============================
      call wfat_init(rtdb)
      

      res_mode = .false.

      
c     ==== Check whether the current theory is DFT ====
      if(.not. rtdb_cget(rtdb,'wfat:theory',1,theory))
     &     call errquit('WFAT: Error reading previous task theory '//
     &     'from RTDB', 1, INPUT_ERR)
      if (trim(adjustl(theory)) /= 'RHF' .and.
     &    trim(adjustl(theory)) /= 'UHF') then
         if (ndid==0) write(luout, '("WFAT: Error because previous '//
     &        'SCF type is ", a, " instead of RHF.")')
     &        trim(theory)
         call errquit('      WFAT task operation is specified while '//
     &        'the theory is not DFT. At the moment WFAT can only be '//
     &        'used in conjunction with single determinant methods '//
     &        'to obtain quantities such as molecular orbitals and '//
     &        'effective potential. Thus terminating NWChem ...',
     &        1, INPUT_ERR)
      endif
  

      
c     ==== Task title ====
      call ga_sync()
      call wfat_prep_tasktitle('OE')


c     ==== Reading general parameters ====
      call wfat_prep_general(rtdb, gqual, lmax)

      
c     ==== Checking the numerical density ====
      if (.not. MA_alloc_get(MT_DBL, 2, 'numden', l_nd, k_nd))
     &     call errquit('WFAT: numden allocation failed', 59, MA_ERR)
      call dfill(2, 0.0d0, dbl_mb(k_nd), 1)
      call wfat_numdens(bas_hd, geom_hd, nspin0, gqual,
     &        dbl_mb(k_nd : k_nd+1))

      
c     ==== Obtaining ionization channels ====
      call wfat_prep_channels(rtdb, nel, nch, l_or0, k_or0, l_sp, k_sp,
     &     l_or, k_or)


c     ==== Obtaining parabolic quantum numbers ====
      call wfat_prep_paraqn(rtdb, nqn, l_npa, k_npa, l_mpa, k_mpa)
      

c     ==== Obtaining WFAT(1) information ====
      call wfat_prep_wfat1(rtdb, wfat1, nlgl)

      
c     ==== Check if research mode is enabled ====
      call wfat_prep_resmode(rtdb, res_mode)


c     ==== Obtaining print switches ====
      call wfat_prep_printsw(rtdb, pr_srfc, pr_rate, pr_yild, pr_avbeta,
     &     pr_avgamma, pr_tmd, pr_tmd_atom)

      
c     ==== Whether using TD field or no ====
      if (.not. rtdb_get(rtdb, 'wfat:tdfield', mt_log, 1, td))
     &     td = .false.
      if (td) then
         call wfat_prep_tdfield(rtdb, nlas, l_al, k_al, las_env,
     &        l_wl, k_wl, l_pl, k_pl, l_tl, k_tl, l_dl, k_dl, ntime,
     &        tmax, extf)
      else
         call wfat_prep_field(rtdb, fieldmag)
      endif


c     ==== Perturbation theory analysis ====
      if (.not. td) call wfat_calc_perturb(rtdb, fieldmag, pr_rate, 
     &     pr_avbeta, pr_avgamma, pr_tmd, pr_tmd_atom, res_mode,
     &     nch, int_mb(k_or), int_mb(k_sp))


cfordebug     ==== Some tests ====
cfordebug      call wfat_test_numerfexch(rtdb)
      

c     ==== Allocate channel properties ====
      call wfat_prep_channelprop(nch, nerg, ncpos, ncmm, niond,
     &     l_eg, k_eg, l_ps, k_ps, l_cmm, k_cmm, l_idp, k_idp)
      

c     ==== Calculate (or load) WFAT integrals ====
      nsh = (1 + lmax)**2
      nintg = nsh * natm * nqn * nch
      if (.not. ga_create(mt_dcpl, nintg, 1, 'intgl', 0, 0, g_il))
     &     call errquit('WFAT: Cannot create intgl global array.',
     &     1, GA_ERR)
      call ga_zero(g_il)

      
c     ==== Preparing WFAT(1) integrals ====
      call wfat_prep_wfat1_intgl(wfat1, lmax, natm, nqn, nch, g_j, g_k,
     &     l_pb, k_pb)
      
      
      smh = 0.0d0
      call util_get_rtdb_state(rtdb, os, oc, or)
      if (os .and. .not. or) then
c        ==== Calculate WFAT integrals on 'start' startup mode ====
         call wfat_calc_intgl(rtdb, nch, int_mb(k_or), int_mb(k_sp),
     &                        nqn, int_mb(k_npa), int_mb(k_mpa),
     &                        lmax, gqual, smh, res_mode, wfat1, nlgl,
     &                        dbl_mb(k_eg), dbl_mb(k_ps), dbl_mb(k_cmm),
     &                        dbl_mb(k_idp), dbl_mb(k_pb), 
     &                        g_il, g_j, g_k)
         
c        ==== Binary-print the integrals freshly calculated ====
         if (.not. rtdb_get(rtdb, 'wfat:intgl_nosave', mt_log, 1,
     &        int_nosave))
     &        call errquit('WFAT: Cannot read intgl_nosave from RTDB',
     &        1, RTDB_ERR)    
         if (ndid == 0 .and. .not. int_nosave) then
            call wfat_store_intgl(rtdb, nsh, natm, nqn, nch,
     &           dbl_mb(k_eg), dbl_mb(k_ps), dbl_mb(k_cmm),
     &           dbl_mb(k_idp), g_il)
         else if (ndid == 0 .and. int_nosave) then
            write(luout, '("ATTENTION: WFAT integrals are not saved.")')
         endif
         call ga_sync()
      else if (.not. os .and. or) then
c        ==== Binary-load WFAT integrals on 'restart' startup mode ====
         if (ndid == 0) then
            call wfat_load_intgl(rtdb, nsh, natm, nqn, nch,
     &           dbl_mb(k_eg), dbl_mb(k_ps), dbl_mb(k_cmm),
     &           dbl_mb(k_idp), g_il)
            call wfat_print_orb_property(nch, int_mb(k_or), 
     &           int_mb(k_sp), dbl_mb(k_eg), dbl_mb(k_ps),
     &           dbl_mb(k_cmm), dbl_mb(k_idp))
         endif
         call ga_sync()
         
         call ga_brdcst(msg_wfat_eg, dbl_mb(k_eg),
     &                  ma_sizeof(MT_DBL, nch, MT_BYTE), 0)
         call ga_brdcst(msg_wfat_cmm, dbl_mb(k_cmm),
     &                  ma_sizeof(MT_DBL, 3*nch, MT_BYTE), 0)
      else
         call errquit('WFAT: WFAT module requires the startup '//
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
     &           call wfat_calc_td(rtdb, geom_hd, pr_srfc, pr_yild, 
     &           pr_avbeta, pr_avgamma, nlas, las_env,
     &           dbl_mb(k_al), dbl_mb(k_wl), dbl_mb(k_pl), dbl_mb(k_tl),
     &           dbl_mb(k_dl), ntime, tmax, extf, nch, int_mb(k_or0), 
     &           int_mb(k_sp), dbl_mb(k_eg), dbl_mb(k_cmm), nqn,
     &           int_mb(k_npa), int_mb(k_mpa), lmax, zion,
     &           dcpl_mb(k_in))
         else
            if (pr_srfc .or. pr_rate)
     &           call wfat_calc_ti(rtdb, geom_hd, pr_srfc, pr_rate,
     &           fieldmag, nch, int_mb(k_or0), int_mb(k_sp), 
     &           dbl_mb(k_eg), dbl_mb(k_cmm), nqn, int_mb(k_npa),
     &           int_mb(k_mpa), lmax, zion, dcpl_mb(k_in))

            if (pr_tmd .or. pr_tmd_atom) 
     &           call wfat_calcprint_tmd(rtdb, geom_hd, pr_tmd, 
     &           pr_tmd_atom, fieldmag, nch, int_mb(k_or0), 
     &           int_mb(k_sp), dbl_mb(k_eg), dbl_mb(k_cmm), nqn, 
     &           int_mb(k_npa), int_mb(k_mpa), lmax, dcpl_mb(k_in),
     &           res_mode)
         endif
         
         
         if (.not. MA_free_heap(l_in))
     &        call errquit('WFAT: intgl deallocation failed', 3, MA_ERR)
         rtdb_stat = rtdb_parallel(.true.)
      endif
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_nd))
     &     call errquit('WFAT: numden deallocation failed',
     &     3, MA_ERR)
      call wfat_erase_channels(l_or0, l_sp, l_or)
      call wfat_erase_paraqn(l_npa, l_mpa)
      if (td) 
     &     call wfat_erase_tdfield(l_al, l_wl, l_pl, l_tl, l_dl)
      call wfat_erase_channelprop(l_eg, l_ps, l_cmm, l_idp)
      if (.not. ga_destroy(g_il))
     &     call errquit('WFAT: Cannot destroy intgl.', 1, GA_ERR)
      call wfat_erase_wfat1_intgl(g_j, g_k, l_pb)


      call wfat_terminate(rtdb)


cforCHECK      if (.not. rtdb_print(rtdb, .false.))
cforCHECK     &     call errquit('WFAT: Cannot print RTDB.', 1, RTDB_ERR)

      
      task_wfat = .true.
      call wfat_close_task('OE')
      
      end
c===============================================================


cOLD c===============================================================
cOLD c===============================================================
cOLD       subroutine wfat_calcprint_srfc_iony(rtdb, geom, pr_srfc, pr_iony, 
cOLD      &           pr_avbeta, pr_avgamma, td, nlas, ftlas, alas, wlas, 
cOLD      &           plas, tlas, dlas, nt, tmax, nch, orbid0, spinid, erg, 
cOLD      &           chmmt, nqn, nid, mid, lmax, zi, intgl)
cOLD 
cOLD c     DESCRIPTION:
cOLD c     If 'td' is true, calculate and print ion yield. If 'td' is false, calculate and print
cOLD c     ionization rate.
cOLD       
cOLD       implicit none
cOLD #include "global.fh"      
cOLD #include "errquit.fh"
cOLD #include "stdio.fh"
cOLD #include "mafdecls.fh"
cOLD #include "rtdb.fh"
cOLD #include "geom.fh"
cOLD ccOLD #include "wfat_qm_param.fh"      
cOLD       integer rtdb              ! [input]
cOLD       integer geom
cOLD       logical pr_srfc                 ! [input], if .true. then print structure factor.
cOLD       logical pr_iony                 ! [input], if .true. then print ion yield.
cOLD       logical pr_avbeta               ! [input], if .true. then print beta-averaged ion yield.
cOLD       logical pr_avgamma              ! [input], if .true. then print gamma-averaged ion yield.
cOLD       logical td                    ! if td=.false., put any scalar for ftlas, wlas, plas, tlas, and dlas defined below. alas, on the other hand, must contain the field strength.
cOLD       integer nlas
cOLD       character*5 ftlas(nlas)       ! Function type of envelope, 'gauss' or 'hann'         
cOLD       double precision alas(nlas)   ! Peak intensity, TW/cm^2                              
cOLD       double precision wlas(nlas)   ! Wavelength, nm                                       
cOLD       double precision plas(nlas)   ! Carrier-envelope phase, rad                          
cOLD       double precision tlas(nlas)   ! Delay (temporal position of the envelope peak), fs   
cOLD       double precision dlas(nlas)   ! Field FWHM, fs                                       
cOLD       double precision nt
cOLD       double precision tmax
cOLD       integer nch                     ! [input], number of ionized channels.
cOLD       integer orbid0(nch)             ! [input], orbital ID of ionized channels.
cOLD       integer spinid(nch)             ! [input], spin value (1 or 2) of ionized channels.
cOLD       double precision erg(nch)       ! [input]
cOLD       double precision chmmt(3, nch)  ! [input]
cOLD       integer nqn                     ! [input]
cOLD       integer nid(nqn)                ! [input]
cOLD       integer mid(nqn)                ! [input]
cOLD       integer lmax                    ! [input]
cOLD       double precision zi             ! [input]
cOLD       double complex intgl(1:*)       ! [input]
cOLD 
cOLD       logical gap
cOLD       integer nsf, nfc, niy, nbt, ngm, nat
cOLD       integer ib, ig, ndid
cOLD       integer l_sf, l_sfn, l_fc, l_iy, l_aiy
cOLD       integer k_sf, k_sfn, k_fc, k_iy, k_aiy
cOLD       double precision bt_start, bt_last, dbt, dbtr, bt
cOLD       double precision gm_start, gm_last, dgm, dgmr, gm
cOLD       double precision ionyd, a_ionyd
cOLD       double precision deg2rad, pi_const
cOLD 
cOLD       parameter (deg2rad = 0.01745329251994d0)
cOLD       parameter (pi_const = 3.14159265359d0)
cOLD 
cOLD 
cOLD       ndid = ga_nodeid()
cOLD       if (ndid /= 0) return
cOLD 
cOLD 
cOLD 
cOLD c     ==== Obtaining info on beta ====
cOLD       if (.not. rtdb_get(rtdb, 'wfat:beta_start', mt_dbl, 1,
cOLD      &     bt_start))
cOLD      &     call errquit('WFAT: Failed reading beta_start from rtdb',
cOLD      &     1, INPUT_ERR)
cOLD       if (.not. rtdb_get(rtdb, 'wfat:beta_last', mt_dbl, 1,
cOLD      &     bt_last))
cOLD      &     call errquit('WFAT: Failed reading beta_last from rtdb',
cOLD      &     1, INPUT_ERR)
cOLD       if (.not. rtdb_get(rtdb, 'wfat:nbeta', mt_int, 1, nbt))
cOLD      &     call errquit('WFAT: Failed reading nbeta from rtdb',
cOLD      &     1, INPUT_ERR)
cOLD 
cOLD       if (nbt > 1) then
cOLD          dbt = (bt_last - bt_start) / (nbt - 1)
cOLD       else
cOLD          dbt = 0.0d0
cOLD       endif
cOLD       dbtr = dbt * deg2rad
cOLD 
cOLD 
cOLD c     ==== Obtaining info on gamma ====
cOLD       if (.not. rtdb_get(rtdb, 'wfat:gamma_start', mt_dbl, 1,
cOLD      &     gm_start))
cOLD      &     call errquit('WFAT: Failed reading gamma_start from rtdb',
cOLD      &     1, INPUT_ERR)
cOLD       if (.not. rtdb_get(rtdb, 'wfat:gamma_last', mt_dbl, 1,
cOLD      &     gm_last))
cOLD      &     call errquit('WFAT: Failed reading gamma_last from rtdb',
cOLD      &     1, INPUT_ERR)
cOLD       if (.not. rtdb_get(rtdb, 'wfat:ngamma', mt_int, 1, ngm))
cOLD      &     call errquit('WFAT: Failed reading ngamma from rtdb',
cOLD      &     1, INPUT_ERR)
cOLD 
cOLD       if (ngm > 1) then
cOLD          dgm = (gm_last - gm_start) / (ngm - 1)
cOLD       else
cOLD          dgm = 0.0d0
cOLD       endif
cOLD       dgmr = dgm * deg2rad
cOLD 
cOLD       
cOLD c     ==== Allocation of structure factors and ion yields ====
cOLD       if (.not. geom_ncent(geom, nat))
cOLD      &     call errquit('WFAT: Getting natom failed', 31, GEOM_ERR)
cOLD       nsf = nat * nqn * nch
cOLD       if (.not. MA_alloc_get(MT_DCPL, nsf, 'srfc', l_sf, k_sf))
cOLD      &     call errquit('WFAT: srfc allocation failed', 2, MA_ERR)
cOLD       if (td) then
cOLD          if (.not. MA_alloc_get(MT_DCPL, nsf, 'srfcn', l_sfn, k_sfn))
cOLD      &        call errquit('WFAT: srfc allocation failed', 2, MA_ERR)
cOLD       endif
cOLD       nfc = nqn * nch
cOLD       if (.not. MA_alloc_get(MT_DBL, nfc, 'ffc', l_fc, k_fc))
cOLD      &     call errquit('WFAT: ffc allocation failed', 2, MA_ERR)
cOLD       niy = nqn * nch
cOLD       if (.not. MA_alloc_get(MT_DBL, niy, 'p_ionyd', l_iy, k_iy))
cOLD      &     call errquit('WFAT: p_ionyd allocation failed', 2, MA_ERR)
cOLD       if (pr_avbeta .or. pr_avgamma) then
cOLD          if (.not. MA_alloc_get(MT_DBL, niy, 'p_aionyd', l_aiy, k_aiy))
cOLD      &        call errquit('WFAT: p_aionyd allocation failed',
cOLD      &        2, MA_ERR)
cOLD       endif
cOLD       
cOLD 
cOLD       
cOLD c     ==== Initializing the printings ====
cOLD       if (ndid == 0 .and. pr_srfc)
cOLD      &     call wfat_initprint_structuref(nch, orbid0,
cOLD      &          spinid, nqn, nid, mid) 
cOLD       if (ndid == 0 .and. pr_iony)
cOLD      &     call wfat_initprint_ionyield(nch, orbid0, 
cOLD      &          spinid, nqn, nid, mid)
cOLD       if (ndid == 0 .and. pr_avbeta) 
cOLD      &     call wfat_initprint_avbeta(nch, orbid0,
cOLD      &          spinid, nqn, nid, mid)
cOLD       
cOLD       
cOLD c     ==== Calculate and print structure factor and ion yield ====
cOLD       do ig = 1, ngm
cOLD          gm = gm_start + dgm*(ig-1)
cOLD          if (pr_iony .and. pr_avbeta) then
cOLD             dbl_mb(k_aiy:k_aiy+niy-1) = 0.0d0
cOLD             a_ionyd = 0.0d0
cOLD          endif
cOLD          
cOLD          do ib = 1, nbt
cOLD             bt = bt_start + dbt*(ib-1)
cOLD 
cOLD             if (ib /= nbt) then
cOLD                gap = .false.
cOLD             else
cOLD                gap = .true.
cOLD             endif
cOLD             
cOLD             if (pr_srfc .or. pr_iony .or. pr_avbeta) then
cOLD c              ==== Calculate structure factor ====
cOLD                call wfat_calc_structuref(geom, bt, gm, nch, erg,
cOLD      &              nqn, mid, lmax, chmmt, intgl, dcpl_mb(k_sf))
cOLD                if (td)
cOLD      &              call wfat_calc_structuref(geom, 180.0d0-bt, 
cOLD      &              180.0d0+gm, nch, erg, nqn, mid, lmax, chmmt, intgl,
cOLD      &              dcpl_mb(k_sfn))
cOLD                
cOLD                
cOLD c              ==== Print structure factor ====
cOLD                if (ndid == 0 .and. pr_srfc)
cOLD      &              call wfat_print_structuref(nat, bt, gm, nch, nqn,
cOLD      &                                         dcpl_mb(k_sf), gap)
cOLD             endif
cOLD 
cOLD             if (pr_iony .or. pr_avbeta) then
cOLD c              ==== Calculate ion yield ====
cOLD                if (td) then
cOLD                   call wfat_calc_tdionyield(nat, zi, nlas, ftlas, alas, 
cOLD      &                 wlas, plas, tlas, dlas, nt, tmax, nch, erg, nqn, 
cOLD      &                 nid, mid, dcpl_mb(k_sf), dcpl_mb(k_sfn), 
cOLD      &                 dbl_mb(k_fc), dbl_mb(k_iy), ionyd)
cOLD                else
cOLD                   call wfat_calc_ionyield(nat, alas(1), nch, erg,
cOLD      &                 nqn, nid, mid, zi, dcpl_mb(k_sf), dbl_mb(k_fc),
cOLD      &                 dbl_mb(k_iy), ionyd)
cOLD                endif
cOLD                
cOLD 
cOLD c              ==== Print ion yield ====
cOLD                if (ndid == 0 .and. pr_iony)
cOLD      &              call wfat_print_ionyield(bt, gm, nch, nqn, 
cOLD      &                                       dbl_mb(k_iy), ionyd, gap)
cOLD 
cOLD c              ==== Calculate beta-averaged ion yield ====
cOLD                if (pr_avbeta) then
cOLD                   if (.not. td)
cOLD      &                 dbl_mb(k_aiy:k_aiy+niy-1) = 
cOLD      &                 dbl_mb(k_aiy:k_aiy+niy-1) +
cOLD      &                 dbl_mb(k_iy:k_iy+niy-1) * dbtr *
cOLD      &                 sin(bt*deg2rad) / 2
cOLD                   a_ionyd = a_ionyd + ionyd * dbtr * sin(bt*deg2rad) / 2
cOLD                endif
cOLD             endif
cOLD             
cOLD          enddo
cOLD 
cOLD c        ==== Print beta-averaged ion yield ====
cOLD          if (ndid == 0 .and. pr_avbeta) then
cOLD             if (td) then
cOLD cc               call wfat_print_avbeta_td(gm, nch, nqn, a_ionyd)       write this subroutine
cOLD             else
cOLD                call wfat_print_avbeta(gm, nch, nqn, dbl_mb(k_aiy),
cOLD      &              a_ionyd)
cOLD             endif
cOLD          endif
cOLD       enddo
cOLD 
cOLD 
cOLD c     ==== Initializing the printings ====
cOLD       if (ndid == 0 .and. pr_avgamma) 
cOLD      &     call wfat_initprint_avgamma(nch, orbid0,
cOLD      &          spinid, nqn, nid, mid)
cOLD 
cOLD       
cOLD c     ==== Calculate gamma-averaged ion yield ====
cOLD       if (pr_avgamma) then
cOLD 
cOLD          do ib = 1, nbt
cOLD             bt = bt_start + dbt*(ib-1)            
cOLD             dbl_mb(k_aiy:k_aiy+niy-1) = 0.0d0
cOLD             a_ionyd = 0.0d0
cOLD 
cOLD             do ig = 1, ngm
cOLD                gm = gm_start + dgm*(ig-1)
cOLD 
cOLD c              ==== Calculate structure factor ====
cOLD                call wfat_calc_structuref(geom, bt, gm, nch, erg,
cOLD      &              nqn, mid, lmax, chmmt, intgl, dcpl_mb(k_sf))
cOLD 
cOLD                if (td)
cOLD      &              call wfat_calc_structuref(geom, 180.0d0-bt, 
cOLD      &              180.0d0+gm, nch, erg, nqn, mid, lmax, chmmt, intgl,
cOLD      &              dcpl_mb(k_sfn))
cOLD 
cOLD c              ==== Calculate ion yield ====
cOLD                if (td) then
cOLD                   call wfat_calc_tdionyield(nat, zi, nlas, ftlas, alas, 
cOLD      &                 wlas, plas, tlas, dlas, nt, tmax, nch, erg, nqn, 
cOLD      &                 nid, mid, dcpl_mb(k_sf), dcpl_mb(k_sfn), 
cOLD      &                 dbl_mb(k_fc), dbl_mb(k_iy), ionyd)
cOLD                else
cOLD                   call wfat_calc_ionyield(nat, alas(1), nch, erg,
cOLD      &                 nqn, nid, mid, zi, dcpl_mb(k_sf), dbl_mb(k_fc), 
cOLD      &                 dbl_mb(k_iy), ionyd)
cOLD                endif
cOLD 
cOLD                if (.not. td)
cOLD      &              dbl_mb(k_aiy:k_aiy+niy-1) =
cOLD      &              dbl_mb(k_aiy:k_aiy+niy-1) +
cOLD      &              dbl_mb(k_iy:k_iy+niy-1) * dgmr /
cOLD      &              2.0d0 / pi_const
cOLD                a_ionyd = a_ionyd + ionyd * dgmr / 2.0d0 / pi_const
cOLD                
cOLD             enddo
cOLD 
cOLD             if (ndid == 0) then
cOLD                if (td) then
cOLD cc                  call wfat_print_avgamma_td(bt, nch, nqn, a_ionyd)       write this subroutine
cOLD                else
cOLD                   call wfat_print_avgamma(bt, nch, nqn, dbl_mb(k_aiy),
cOLD      &                                    a_ionyd)            
cOLD                endif
cOLD             endif
cOLD             
cOLD          enddo
cOLD       endif
cOLD       
cOLD 
cOLD c     ==== Terminating printings ====      
cOLD       if (ndid==0 .and. pr_srfc) call wfat_endprint_structuref()
cOLD       if (ndid==0 .and. pr_iony) call wfat_endprint_ionyield()
cOLD       if (ndid==0 .and. pr_avbeta) call wfat_endprint_avbeta()
cOLD       if (ndid==0 .and. pr_avgamma) call wfat_endprint_avgamma()
cOLD 
cOLD       
cOLD c     ==== Deallocations ====
cOLD       if (.not. MA_free_heap(l_sf))
cOLD      &     call errquit('WFAT: srfc deallocation failed',
cOLD      &     3, MA_ERR)
cOLD       if (td) then
cOLD          if (.not. MA_free_heap(l_sfn))
cOLD      &        call errquit('WFAT: srfcn deallocation failed',
cOLD      &        3, MA_ERR)
cOLD       endif
cOLD       if (.not. MA_free_heap(l_fc))
cOLD      &     call errquit('WFAT: ffc deallocation failed',
cOLD      &     3, MA_ERR)
cOLD       if (.not. MA_free_heap(l_iy))
cOLD      &     call errquit('WFAT: p_ionyd deallocation failed',
cOLD      &     3, MA_ERR)
cOLD       if (pr_avgamma .or. pr_avbeta) then
cOLD          if (.not. MA_free_heap(l_aiy))
cOLD      &        call errquit('WFAT: p_aionyd deallocation failed',
cOLD      &        3, MA_ERR)
cOLD       endif      
cOLD 
cOLD       end
cOLD c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_calcprint_tmd(rtdb, geom, pr_tmd, pr_tmd_atom, ef,
     &           nch, orbid0, spinid, erg, chmmt, nqn, nid, mid, lmax,
     &           intgl, rmode)
      
      implicit none
#include "global.fh"      
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
ccOLD #include "wfat_qm_param.fh"            
      integer rtdb
      integer geom
      logical pr_tmd
      logical pr_tmd_atom
      double precision ef
      integer nch
      integer orbid0(nch)
      integer spinid(nch)
      double precision erg(nch)
      double precision chmmt(3, nch)
      integer nqn
      integer nid(nqn)
      integer mid(nqn)
      integer lmax
      double complex intgl(1:*)
      logical rmode

      logical gap
      integer ix, iy, ndid, nat
      integer nsf, ntmd, nktmd, ntms
      integer l_sf, l_tm, l_tms
      integer k_sf, k_tm, k_tms
      double precision kx, ky, kmax, dk, bt_tmd, gm_tmd


      ndid = ga_nodeid()
      if (ndid /= 0) return

      
c     ==== Obtaining beta and gamma ====
      if (.not. rtdb_get(rtdb, 'wfat:beta_tmd', mt_dbl, 1, bt_tmd))
     &     call errquit('WFAT: Failed reading bt_tmd from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:gamma_tmd', mt_dbl, 1, gm_tmd))
     &     call errquit('WFAT: Failed reading gm_tmd from rtdb',
     &     1, INPUT_ERR)

      
c     ==== Allocation of structure factors and TMD ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting natom failed', 77, GEOM_ERR)
      nsf = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nsf, 'srfc', l_sf, k_sf))
     &     call errquit('WFAT: srfc allocation failed', 2, MA_ERR)
      call wfat_calc_structuref(geom, bt_tmd, gm_tmd, nch, erg,
     &                          nqn, mid, lmax, chmmt, intgl,
     &                          dcpl_mb(k_sf))
      ntmd = nat * nch
      if (.not. MA_alloc_get(MT_DCPL, ntmd, 'tmd', l_tm, k_tm))
     &     call errquit('WFAT: tmd allocation failed', 2, MA_ERR)
      ntms = nat * nch
      if (.not. MA_alloc_get(MT_DBL, ntmd, 'tmd_ksum', l_tms, k_tms))
     &     call errquit('WFAT: tmd_ksum allocation failed', 2, MA_ERR)

      
c     ==== Calculating max. momentum to calculate the TMD at ====
      call wfat_get_kmax(ef, nch, erg, nqn, nid, mid, kmax)
      if (.not. rtdb_get(rtdb, 'wfat:nktmd', mt_int, 1, nktmd))
     &     call errquit('WFAT: Failed reading nktmd from rtdb',
     &     1, RTDB_ERR)
      dk = 2 * kmax / (nktmd - 1)
      write(luout, '(/, 3x, "*** INFO ***")')
      write(luout, '(5x, "Max. transverse momentum = ", e14.6, '//
     &      '" rad/Bohr")') kmax


c     ==== Initializing the printings ====      
      if (ndid==0 .and. pr_tmd) call wfat_initprint_tmd(nch, orbid0,
     &                                  spinid)
      if (ndid==0 .and. pr_tmd_atom) call wfat_initprint_tmdar(nch, 
     &                                       nat, orbid0, spinid)


c     ==== Calculate the TMDs ====
      dbl_mb(k_tms:k_tms+ntms-1) = 0.0d0
      ky = -kmax
      do iy = 1, nktmd
         kx = -kmax
         do ix = 1, nktmd

            if (ix /= nktmd) then
               gap = .false.
            else
               gap = .true.
            endif
            call wfat_calc_tmd(kx, ky, nch, erg, nqn, nid, mid,
     &                         ef, dcpl_mb(k_sf), dcpl_mb(k_tm))
            
            if (ndid == 0) then
               if (pr_tmd) call wfat_print_tmd(kx, ky, nch,
     &                          dcpl_mb(k_tm), gap)
               if (pr_tmd_atom) call wfat_print_tmdar(kx, ky, nch,
     &                               nat, dcpl_mb(k_tm), gap)
            endif

            if (pr_tmd_atom .and. rmode)
     &           dbl_mb(k_tms:k_tms+ntms-1) = dbl_mb(k_tms:k_tms+ntms-1)
     &                              + abs(dcpl_mb(k_tm:k_tm+ntmd-1))

            kx = kx + dk
         enddo
         ky = ky + dk
      enddo


      if (pr_tmd_atom .and. rmode .and. ndid==0) then
         dbl_mb(k_tms:k_tms+ntms-1) = dbl_mb(k_tms:k_tms+ntms-1)*dk**2
         call wfat_print_tms(geom, nch, nat, dbl_mb(k_tms))
      endif

      
c     ==== Terminating printings ====
      if (ndid==0 .and. pr_tmd) call wfat_endprint_tmd()
      if (ndid==0 .and. pr_tmd_atom) call wfat_endprint_tmdar()
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_sf))
     &     call errquit('WFAT: srfc deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_tm))
     &     call errquit('WFAT: tmd deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_tms))
     &     call errquit('WFAT: tmd_ksum deallocation failed',
     &     3, MA_ERR)



      
      end
c===============================================================
