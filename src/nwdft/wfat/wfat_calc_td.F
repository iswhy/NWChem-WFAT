c===============================================================
c===============================================================
      subroutine wfat_calc_td(rtdb, geom, pr_srfc, pr_yild, 
     &           pr_avbeta, pr_avgamma, nlas, ftlas, alas, wlas, 
     &           plas, tlas, dlas, nt, tmax, extf, nch, orbid0, spinid,
     &           erg, chmmt, nqn, nid, mid, lmax, zi, intgl, wfat1, nb,
     &           ns, orbid, dpmo, pb, chpos, intgl_j, intgl_k1,
     &           intgl_k2)
c     DESCRIPTION:
      
      implicit none
#include "global.fh"      
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "wfat_field_param.fh"
      integer rtdb              ! [input]
      integer geom
      logical pr_srfc                 ! [input], if .true. then print structure factor.
      logical pr_yild                 ! [input], if .true. then print ion yield.
      logical pr_avbeta               ! [input], if .true. then print beta-averaged ion yield.
      logical pr_avgamma              ! [input], if .true. then print gamma-averaged ion yield.
      integer nlas
      character(len=FLD_TP_LEN) ftlas(nlas)
cOLD      character*5 ftlas(nlas)       ! Function type of envelope, 'gauss' or 'hann'         
      double precision alas(nlas)   ! Peak intensity, TW/cm^2                              
      double precision wlas(nlas)   ! Wavelength, nm                                       
      double precision plas(nlas)   ! Carrier-envelope phase, rad                          
      double precision tlas(nlas)   ! Delay (temporal position of the envelope peak), fs   
      double precision dlas(nlas)   ! Field FWHM, fs                                       
      double precision nt
      double precision tmax
      logical extf
      integer nch                     ! [input], number of ionized channels.
      integer orbid0(nch)             ! [input], orbital ID of ionized channels.
      integer spinid(nch)             ! [input], spin value (1 or 2) of ionized channels.
      double precision erg(nch)       ! [input]
      double precision chmmt(3, nch)  ! [input]
      integer nqn                     ! [input]
      integer nid(nqn)                ! [input]
      integer mid(nqn)                ! [input]
      integer lmax                    ! [input]
      double precision zi             ! [input]
      double complex intgl(1:*)       ! [input]
      logical wfat1                       ! [input]
      integer nb                          ! [input], only used when wfat1=.true.
      integer ns                          ! [input], only used when wfat1=.true.
      integer orbid(1:*)                  ! [input], only used when wfat1=.true.
      double precision dpmo(1:*)          ! [input], only used when wfat1=.true. The minimum size should be nb*nb*ns*3.
      double precision pb(1:*)            ! [input], only used when wfat1=.true. The minimum size should be 3*3*nch.
      double precision chpos(1:*)         ! [input], only used when wfat1=.true.
      double complex intgl_j(1:*)         ! [input], only used when wfat1=.true.
      double complex intgl_k1(1:*)        ! [input], only used when wfat1=.true.
      double complex intgl_k2(1:*)        ! [input], only used when wfat1=.true.
      
      logical gap, sumbtgm
      integer nsf, nbt, ngm, nat, ndid, ndgi, ndgn, nsc, nsx, nmmd, npbd
      integer n_jn, n_kn, n_on
      integer ib, ig, i, iqn, i_on
      integer l_gn, l_gnn, l_dgi, l_dgn, l_sc, l_sx, l_mmd, l_pbd, l_jn
      integer l_kn, l_on, l_jnn, l_knn
      integer k_gn, k_gnn, k_dgi, k_dgn, k_sc, k_sx, k_mmd, k_pbd, k_jn
      integer k_kn, k_on, k_jnn, k_knn
      double precision bt_start, bt_last, dbt, dbtr, bt
      double precision gm_start, gm_last, dgm, dgmr, gm
      double precision yield, a_yield, yieldbt, yieldbtgm, f, mthr, dthr
      double precision deg2rad, pi_const

      parameter (deg2rad = 0.01745329251994d0)
      parameter (pi_const = 3.14159265359d0)


      ndid = ga_nodeid()
      if (ndid /= 0) return


c     ==== Obtaining info on beta ====
      call wfat_get_beta(rtdb, bt_start, bt_last, nbt, dbt)
      dbtr = dbt * deg2rad


c     ==== Obtaining info on gamma ====
      call wfat_get_gamma(rtdb, gm_start, gm_last, ngm, dgm)
      dgmr = dgm * deg2rad

      
c     ==== Allocation of structure factors ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting natom failed', 31, GEOM_ERR)
      nsf = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nsf, 'srfc', l_gn, k_gn))
     &     call errquit('WFAT: srfc allocation failed', 2, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, nsf, 'srfcn', l_gnn, k_gnn))
     &     call errquit('WFAT: srfc allocation failed', 2, MA_ERR)
      

c     ==== When WFAT1 is active ====
      if (wfat1) then
         call wfat_prepare_an(nb, nat, nch, nqn, ndgi, l_dgi, k_dgi,
     &        ndgn, l_dgn, k_dgn, nsc, l_sc, k_sc, nsx, l_sx, k_sx, 
     &        nmmd, l_mmd, k_mmd, npbd, l_pbd, k_pbd,
     &        n_jn, l_jn, k_jn, n_kn, l_kn, k_kn, n_on, l_on, k_on)
         if (.not. MA_alloc_get(MT_DCPL, n_jn, 'jnn', l_jnn, k_jnn))
     &        call errquit('WFAT: jnn allocation failed', 2, MA_ERR)
         if (.not. MA_alloc_get(MT_DCPL, n_kn, 'knn', l_knn, k_knn))
     &        call errquit('WFAT: knn allocation failed', 2, MA_ERR)
         
         call wfat_detect_degen0(nch, orbid, spinid, int_mb(k_dgi),
     &        int_mb(k_dgn))
         if (.not. rtdb_get(rtdb, 'wfat:mthr', mt_dbl, 1, mthr))
     &        call errquit('WFAT: Failed reading mthr from rtdb',
     &        1, RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'wfat:dthr', mt_dbl, 1, dthr))
     &        call errquit('WFAT: Failed reading dthr from rtdb',
     &        1, RTDB_ERR)
         do i = 1, nch
            do iqn = 1, nqn
               i_on = k_on + 2*nqn*(i-1) + 2*(iqn-1)
               call wfat_calc_o(nid(iqn), mid(iqn), erg(i), zi,
     &              dbl_mb(i_on))
               write(luout, '("o_nu : ", 2i4, 2es21.12)')
     &              i, iqn, dbl_mb(i_on), dbl_mb(i_on+1)
            enddo
         enddo
      endif

      
c     ==== Initializing the printings ====
      if (ndid == 0 .and. pr_srfc) then
         call wfat_initprint_structuref(nch, orbid0, spinid, nqn,
     &        nid, mid)
         if (wfat1) call wfat_initprint_srfc1(nch, orbid0, spinid, nqn,
     &        nid, mid)
      endif
      if (ndid == 0 .and. pr_yild)
     &     call wfat_initprint_yield()
      if (ndid == 0 .and. pr_avbeta) 
     &     call wfat_initprint_avbeta()

      
      sumbtgm = .false.
      if (rtdb_get(rtdb, 'wfat:sumbtgm', mt_log, 1, sumbtgm)) then
         if (sumbtgm .and. .not. pr_yild) call errquit('WFAT: When '//
     &        'summation of yields over angles is requested, '//
     &        'print_yield must be used.', 1, UERR)
         write(luout, '("INFO : Summation of yields over angles '//
     &        'is requested.")')
      endif
      
      
c     ==== Calculate and print structure factor and ion yield ====
      yieldbtgm = 0.0d0
      do ig = 1, ngm
         gm = gm_start + dgm*(ig-1)
         if (pr_yild .and. pr_avbeta) then
            a_yield = 0.0d0
         endif

         yieldbt = 0.0d0
         do ib = 1, nbt
            bt = bt_start + dbt*(ib-1)

            if (ib /= nbt) then
               gap = .false.
            else
               gap = .true.
            endif
            
            if (pr_srfc .or. pr_yild .or. pr_avbeta) then
c              ==== Calculate structure factor ====
               if (wfat1) then
                  call wfat_get_srfc1(bt, gm, nat, nb, ns,
     &                 nch, nqn, lmax, orbid, spinid, int_mb(k_dgi),
     &                 int_mb(k_dgn), mid, dpmo, erg, chmmt, pb, chpos, 
     &                 mthr, dthr, intgl, intgl_j, intgl_k1, intgl_k2, 
     &                 dbl_mb(k_sc), dcpl_mb(k_sx), dbl_mb(k_mmd), 
     &                 dbl_mb(k_pbd), dcpl_mb(k_gn), dcpl_mb(k_jn),
     &                 dcpl_mb(k_kn))
                  call wfat_get_srfc1(180.0d0-bt, 180.0d0+gm, nat, nb,
     &                 ns, nch, nqn, lmax, orbid, spinid, int_mb(k_dgi),
     &                 int_mb(k_dgn), mid, dpmo, erg, chmmt, pb, chpos, 
     &                 mthr, dthr, intgl, intgl_j, intgl_k1, intgl_k2, 
     &                 dbl_mb(k_sc), dcpl_mb(k_sx), dbl_mb(k_mmd), 
     &                 dbl_mb(k_pbd), dcpl_mb(k_gnn), dcpl_mb(k_jnn),
     &                 dcpl_mb(k_knn))
               else
                  call wfat_get_srfc0(geom, bt, gm, nch, erg,
     &                 nqn, mid, lmax, chmmt, intgl, dcpl_mb(k_gn))
                  call wfat_get_srfc0(geom, 180.0d0-bt, 
     &                 180.0d0+gm, nch, erg, nqn, mid, lmax, chmmt, 
     &                 intgl, dcpl_mb(k_gnn))
               endif
               
               
c              ==== Print structure factor ====
               if (ndid == 0 .and. pr_srfc) then
                  call wfat_print_structuref(nat, bt, gm, nch, nqn,
     &                 dcpl_mb(k_gn), gap)
                  if (wfat1) call wfat_print_srfc1(nat, bt, gm, nch, 
     &                 nqn, dcpl_mb(k_jn), dcpl_mb(k_kn), gap)
               endif
               
            endif

            if (pr_yild .or. pr_avbeta) then
c              ==== Calculate ion yield ====
               if (wfat1) then
                  call wfat_calc_wfat1_yield(bt, gm, zi, nlas, ftlas, 
     &                 alas, wlas, plas, tlas, dlas, nt, tmax, extf, 
     &                 nch, erg, nqn, nid, mid, nat, dcpl_mb(k_gn), 
     &                 dcpl_mb(k_jn), dcpl_mb(k_kn), dcpl_mb(k_gnn),
     &                 dcpl_mb(k_jnn), dcpl_mb(k_knn), dbl_mb(k_on),
     &                 dbl_mb(k_mmd), dbl_mb(k_pbd), yield)
                  
               else
                  call wfat_calc_yield(nat, zi, nlas, ftlas, alas, wlas,
     &                 plas, tlas, dlas, nt, tmax, extf, nch, erg, nqn,
     &                 nid, mid, dcpl_mb(k_gn), dcpl_mb(k_gnn), yield)
               endif
               
c              ==== Print ion yield ====
               if (ndid == 0 .and. pr_yild)
     &              call wfat_print_yield(bt, gm, yield, gap)

c              ==== Calculate beta-averaged ion yield ====
               if (pr_avbeta) then
                  a_yield = a_yield + yield * dbtr * sin(bt*deg2rad) / 2
               endif
            endif

c           ==== Sum the yield over the range of beta ====
            if (sumbtgm) then
               f = 1.0d0
               if (ib==1 .or. ib==nbt) f = 0.5d0
               yieldbt = yieldbt + f * yield * sin(bt*deg2rad) * dbtr
            endif            
         enddo

c        ==== Print beta-averaged ion yield ====
         if (ndid == 0 .and. pr_avbeta) then
            call wfat_print_avbeta(gm, nch, a_yield)
         endif

c        ==== Sum the yield over the range of gamma ====
         if (sumbtgm) then
            f = 1.0d0
            if (ig==1 .or. ig==ngm) f = 0.5d0
            yieldbtgm = yieldbtgm + f * yieldbt * dgmr
         endif         
      enddo


c     ==== Calculate the integral of yields over 4pi ====
c     ====   solid angle using Lebedev quadrature    ====
      if (sumbtgm) then
         write(luout, '("INFO : Sum of the yield over the specified "'//
     &        '"interval is ", es14.6)') yieldbtgm
         write(luout, '("INFO : beta interval   ", 2f11.4, i5)')
     &        bt_start, bt_last, nbt
         write(luout, '("INFO : gamma interval  ", 2f11.4, i5)')
     &        gm_start, gm_last, ngm

         call wfat_sumyield_lebedev(geom, nlas, ftlas, alas, wlas, 
     &        plas, tlas, dlas, nt, tmax, extf, nch,
     &        erg, chmmt, nqn, nid, mid, lmax, zi, intgl)
      endif
      

c     ==== Initializing the printings ====
      if (ndid == 0 .and. pr_avgamma) 
     &     call wfat_initprint_avgamma()

      
c     ==== Calculate gamma-averaged ion yield ====
      if (pr_avgamma) then
         
         do ib = 1, nbt
            bt = bt_start + dbt*(ib-1)            
            a_yield = 0.0d0

            do ig = 1, ngm
               gm = gm_start + dgm*(ig-1)

c              ==== Calculate structure factor ====
               if (wfat1) then
                  call wfat_get_srfc1(bt, gm, nat, nb, ns,
     &                 nch, nqn, lmax, orbid, spinid, int_mb(k_dgi),
     &                 int_mb(k_dgn), mid, dpmo, erg, chmmt, pb, chpos, 
     &                 mthr, dthr, intgl, intgl_j, intgl_k1, intgl_k2, 
     &                 dbl_mb(k_sc), dcpl_mb(k_sx), dbl_mb(k_mmd), 
     &                 dbl_mb(k_pbd), dcpl_mb(k_gn), dcpl_mb(k_jn),
     &                 dcpl_mb(k_kn))
                  call wfat_get_srfc1(180.0d0-bt, 180.0d0+gm, nat, nb,
     &                 ns, nch, nqn, lmax, orbid, spinid, int_mb(k_dgi),
     &                 int_mb(k_dgn), mid, dpmo, erg, chmmt, pb, chpos, 
     &                 mthr, dthr, intgl, intgl_j, intgl_k1, intgl_k2, 
     &                 dbl_mb(k_sc), dcpl_mb(k_sx), dbl_mb(k_mmd), 
     &                 dbl_mb(k_pbd), dcpl_mb(k_gnn), dcpl_mb(k_jnn),
     &                 dcpl_mb(k_knn))
               else
                  call wfat_get_srfc0(geom, bt, gm, nch, erg,
     &                 nqn, mid, lmax, chmmt, intgl, dcpl_mb(k_gn))
                  call wfat_get_srfc0(geom, 180.0d0-bt, 
     &                 180.0d0+gm, nch, erg, nqn, mid, lmax, chmmt, 
     &                 intgl, dcpl_mb(k_gnn))
               endif

c              ==== Calculate ion yield ====
               if (wfat1) then
                  call wfat_calc_wfat1_yield(bt, gm, zi, nlas, ftlas, 
     &                 alas, wlas, plas, tlas, dlas, nt, tmax, extf, 
     &                 nch, erg, nqn, nid, mid, nat, dcpl_mb(k_gn), 
     &                 dcpl_mb(k_jn), dcpl_mb(k_kn), dcpl_mb(k_gnn),
     &                 dcpl_mb(k_jnn), dcpl_mb(k_knn), dbl_mb(k_on),
     &                 dbl_mb(k_mmd), dbl_mb(k_pbd), yield)
               else
                  call wfat_calc_yield(nat, zi, nlas, ftlas, alas, 
     &                 wlas, plas, tlas, dlas, nt, tmax, extf, nch, erg,
     &                 nqn, nid, mid, dcpl_mb(k_gn), dcpl_mb(k_gnn),
     &                 yield)
               endif

               a_yield = a_yield + yield * dgmr / 2.0d0 / pi_const
            enddo

            if (ndid == 0) call wfat_print_avgamma(bt, a_yield)
         enddo
      endif
      

c     ==== Terminating printings ====      
      if (ndid==0 .and. pr_srfc) then
         call wfat_endprint_structuref()
         if (wfat1) call wfat_endprint_srfc1()
      endif
      if (ndid==0 .and. pr_yild) call wfat_endprint_yield()
      if (ndid==0 .and. pr_avbeta) call wfat_endprint_avbeta()
      if (ndid==0 .and. pr_avgamma) call wfat_endprint_avgamma()

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_gn))
     &     call errquit('WFAT: srfc deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_gnn))
     &     call errquit('WFAT: srfcn deallocation failed', 3, MA_ERR)
      if (wfat1) then
         call wfat_erase_an(l_dgi, l_dgn, l_sc, l_sx, l_mmd, 
     &        l_pbd, l_jn, l_kn, l_on)
         if (.not. MA_free_heap(l_jnn))
     &        call errquit('WFAT: jnn deallocation failed', 3, MA_ERR)
         if (.not. MA_free_heap(l_knn))
     &        call errquit('WFAT: knn deallocation failed', 3, MA_ERR)
      endif
      
      end
c===============================================================
