c===============================================================
c===============================================================
      subroutine wfat_calc_td(rtdb, geom, pr_srfc, pr_yild, 
     &           pr_avbeta, pr_avgamma, nlas, ftlas, alas, wlas, 
     &           plas, tlas, dlas, nt, tmax, extf, nch, orbid0, spinid,
     &           erg, chmmt, nqn, nid, mid, lmax, zi, intgl)

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

      logical gap
      integer nsf, nbt, ngm, nat
      integer ib, ig, ndid
      integer l_sf, l_sfn
      integer k_sf, k_sfn
      double precision bt_start, bt_last, dbt, dbtr, bt
      double precision gm_start, gm_last, dgm, dgmr, gm
      double precision yield, a_yield
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

      
c     ==== Allocation of structure factors and ion yields ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting natom failed', 31, GEOM_ERR)
      nsf = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nsf, 'srfc', l_sf, k_sf))
     &     call errquit('WFAT: srfc allocation failed', 2, MA_ERR)
      if (.not. MA_alloc_get(MT_DCPL, nsf, 'srfcn', l_sfn, k_sfn))
     &     call errquit('WFAT: srfc allocation failed', 2, MA_ERR)
      

      
c     ==== Initializing the printings ====
      if (ndid == 0 .and. pr_srfc)
     &     call wfat_initprint_structuref(nch, orbid0,
     &          spinid, nqn, nid, mid) 
      if (ndid == 0 .and. pr_yild)
     &     call wfat_initprint_yield()
      if (ndid == 0 .and. pr_avbeta) 
     &     call wfat_initprint_avbeta()
      
      
c     ==== Calculate and print structure factor and ion yield ====
      do ig = 1, ngm
         gm = gm_start + dgm*(ig-1)
         if (pr_yild .and. pr_avbeta) then
            a_yield = 0.0d0
         endif
         
         do ib = 1, nbt
            bt = bt_start + dbt*(ib-1)

            if (ib /= nbt) then
               gap = .false.
            else
               gap = .true.
            endif
            
            if (pr_srfc .or. pr_yild .or. pr_avbeta) then
c              ==== Calculate structure factor ====
               call wfat_calc_structuref(geom, bt, gm, nch, erg,
     &              nqn, mid, lmax, chmmt, intgl, dcpl_mb(k_sf))
               call wfat_calc_structuref(geom, 180.0d0-bt, 180.0d0+gm, 
     &              nch, erg, nqn, mid, lmax, chmmt, intgl,
     &              dcpl_mb(k_sfn))
               
c              ==== Print structure factor ====
               if (ndid == 0 .and. pr_srfc)
     &              call wfat_print_structuref(nat, bt, gm, nch, nqn,
     &                                         dcpl_mb(k_sf), gap)
            endif

            if (pr_yild .or. pr_avbeta) then
c              ==== Calculate ion yield ====
               call wfat_calc_yield(nat, zi, nlas, ftlas, alas, wlas, 
     &              plas, tlas, dlas, nt, tmax, extf, nch, erg, nqn, 
     &              nid, mid, dcpl_mb(k_sf), dcpl_mb(k_sfn), yield)
               
c              ==== Print ion yield ====
               if (ndid == 0 .and. pr_yild)
     &              call wfat_print_yield(bt, gm, yield, gap)

c              ==== Calculate beta-averaged ion yield ====
               if (pr_avbeta) then
                  a_yield = a_yield + yield * dbtr * sin(bt*deg2rad) / 2
               endif
            endif
            
         enddo

c        ==== Print beta-averaged ion yield ====
         if (ndid == 0 .and. pr_avbeta) then
            call wfat_print_avbeta(gm, nch, a_yield)
         endif
      enddo


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
               call wfat_calc_structuref(geom, bt, gm, nch, erg,
     &              nqn, mid, lmax, chmmt, intgl, dcpl_mb(k_sf))
               call wfat_calc_structuref(geom, 180.0d0-bt, 180.0d0+gm,
     &              nch, erg, nqn, mid, lmax, chmmt, intgl,
     &              dcpl_mb(k_sfn))

c              ==== Calculate ion yield ====
               call wfat_calc_yield(nat, zi, nlas, ftlas, alas, 
     &              wlas, plas, tlas, dlas, nt, tmax, extf, nch, erg, 
     &              nqn, nid, mid, dcpl_mb(k_sf), dcpl_mb(k_sfn), yield)

               a_yield = a_yield + yield * dgmr / 2.0d0 / pi_const
            enddo

            if (ndid == 0) call wfat_print_avgamma(bt, a_yield)
         enddo
      endif
      

c     ==== Terminating printings ====      
      if (ndid==0 .and. pr_srfc) call wfat_endprint_structuref()
      if (ndid==0 .and. pr_yild) call wfat_endprint_yield()
      if (ndid==0 .and. pr_avbeta) call wfat_endprint_avbeta()
      if (ndid==0 .and. pr_avgamma) call wfat_endprint_avgamma()

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_sf))
     &     call errquit('WFAT: srfc deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_sfn))
     &     call errquit('WFAT: srfcn deallocation failed',
     &     3, MA_ERR)
      
      
      end
c===============================================================
