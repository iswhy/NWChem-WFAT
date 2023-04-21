c===============================================================
c===============================================================
      subroutine wfat_calc_ti(rtdb, geom, pr_srfc, pr_rate, fieldm,
     &           nch, orbid0, spinid, erg, chmmt, nqn, nid, mid, lmax, 
     &           zi, intgl, wfat1, nb, ns, orbid, dpmo, pb, intgl_j,
     &           intgl_k1, intgl_k2)

c     DESCRIPTION:
      
      implicit none
#include "global.fh"      
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
      integer rtdb                    ! [input]
      integer geom                    ! [input]
      logical pr_srfc                 ! [input], if .true. then print structure factor.
      logical pr_rate                 ! [input], if .true. then print ion yield.
      double precision fieldm
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
      integer orbid(nch)                  ! [input], only used when wfat1=.true.
      double precision dpmo(nb,nb,ns,3)   ! [input], only used when wfat1=.true.
      double precision pb(3,3,nch)        ! [input], only used when wfat1=.true.
      double complex intgl_j(1:*)         ! [input], only used when wfat1=.true.
      double complex intgl_k1(1:*)        ! [input], only used when wfat1=.true.
      double complex intgl_k2(1:*)        ! [input], only used when wfat1=.true.
      
      logical gap
      integer nsf, nfc, nrt, nbt, ngm, nat, ndgi, ndgn, nsc, nsx, nmmd
      integer nlpb, n_an, ib, ig, ndid
      integer l_sf, l_fc, l_rt, l_dgi, l_dgn, l_sc, l_sx, l_mmd, l_lpb
      integer l_an
      integer k_sf, k_fc, k_rt, k_dgi, k_dgn, k_sc, k_sx, k_mmd, k_lpb
      integer k_an
      double precision bt_start, bt_last, dbt, dbtr, bt
      double precision gm_start, gm_last, dgm, dgmr, gm
      double precision rate
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
      nfc = nqn * nch
      if (.not. MA_alloc_get(MT_DBL, nfc, 'ffc', l_fc, k_fc))
     &     call errquit('WFAT: ffc allocation failed', 2, MA_ERR)
      nrt = nqn * nch
      if (.not. MA_alloc_get(MT_DBL, nrt, 'p_rate', l_rt, k_rt))
     &     call errquit('WFAT: p_rate allocation failed', 2, MA_ERR)


c     ==== When WFAT1 is active ====
      if (wfat1) then
         call wfat_prepare_an(nb, nat, nch, nqn, ndgi, l_dgi, k_dgi,
     &        ndgn, l_dgn, k_dgn, nsc, l_sc, k_sc, nsx, l_sx, k_sx, 
     &        nmmd, l_mmd, k_mmd, nlpb, l_lpb, k_lpb, n_an, l_an, k_an)
         call wfat_detect_degen0(nch, orbid, spinid, int_mb(k_dgi),
     &        int_mb(k_dgn))
      endif
      
      
c     ==== Initializing the printings ====
      if (ndid == 0 .and. pr_srfc)
     &     call wfat_initprint_structuref(nch, orbid0, spinid, nqn,
     &          nid, mid) 
      if (ndid == 0 .and. pr_rate)
     &     call wfat_initprint_rate(nch, orbid0, spinid, nqn, nid, mid)
      
      
c     ==== Calculate and print structure factor and ion yield ====
      do ig = 1, ngm
         gm = gm_start + dgm*(ig-1)
         
         do ib = 1, nbt
            bt = bt_start + dbt*(ib-1)

            if (ib /= nbt) then
               gap = .false.
            else
               gap = .true.
            endif
            
            if (pr_srfc .or. pr_rate) then
c              ==== Calculate structure factor ====
               call wfat_calc_structuref(geom, bt, gm, nch, erg,
     &              nqn, mid, lmax, chmmt, intgl, dcpl_mb(k_sf))

               if (wfat1) call wfat_calc_an(bt, gm, nat, nb, ns, nch,
     &              nqn, lmax, orbid, spinid, int_mb(k_dgi), 
     &              int_mb(k_dgn), mid, dpmo, chmmt, pb, intgl_j, 
     &              intgl_k1, intgl_k2, dbl_mb(k_sc), dcpl_mb(k_sx), 
     &              dbl_mb(k_mmd), dbl_mb(k_lpb), dcpl_mb(k_sf),
     &              dcpl_mb(k_an))

c              ==== Print structure factor ====
               if (ndid == 0 .and. pr_srfc)
     &              call wfat_print_structuref(nat, bt, gm, nch, nqn,
     &              dcpl_mb(k_sf), gap)

cforlater               if (wfat1) call wfat_print_an() resolve here
                  
            endif

            if (pr_rate) then
c              ==== Calculate ionization rate ====
               call wfat_calc_rate(nat, fieldm, nch, erg,
     &              nqn, nid, mid, zi, dcpl_mb(k_sf), dbl_mb(k_fc),
     &              dbl_mb(k_rt), rate)
               
c              ==== Print ionization rate ====
               if (ndid == 0 .and. pr_rate)
     &              call wfat_print_rate(bt, gm, nch, nqn, 
     &              dbl_mb(k_rt), rate, gap)
            endif
            
         enddo
      enddo
      

c     ==== Terminating printings ====      
      if (ndid==0 .and. pr_srfc) call wfat_endprint_structuref()
      if (ndid==0 .and. pr_rate) call wfat_endprint_rate()

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(l_sf))
     &     call errquit('WFAT: srfc deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_fc))
     &     call errquit('WFAT: ffc deallocation failed',
     &     3, MA_ERR)
      if (.not. MA_free_heap(l_rt))
     &     call errquit('WFAT: p_rate deallocation failed',
     &     3, MA_ERR)
      if (wfat1) call wfat_erase_an(l_dgi, l_dgn, l_sc, l_sx, l_mmd, 
     &           l_lpb, l_an)

      end
c===============================================================


c============================================================
c============================================================
      subroutine wfat_detect_degen0(nch, orbid, spinid, dgid0, dgn0)

      implicit none
#include "wfat_qm_param.fh"
#include "mafdecls.fh"
      integer nch          ! [input],
      integer orbid(nch)   ! [input],
      integer spinid(nch)  ! [input],
      integer dgid0(nch)   ! [output], the degenerate-ID of the ionized channels whose MO ID and spin ID are specified in orbid and spinid, respectively.
      integer dgn0(nch)    ! [output], the degree of degeneracy of the ionized channels.

      integer i, io, is

c     dgid0 is a subset of int_mb(k_dgid:k_dgid+nbas0*nspin0-1)
c     dgn0 is a subset of int_mb(k_dgn:k_dgn+nbas0*nspin0-1)
      do i = 1, nch
         io = orbid(i)
         is = spinid(i)
         dgid0(i) = int_mb(k_dgid + (is-1)*nbas0 + (io-1))
         dgn0(i) = int_mb(k_dgn + (is-1)*nbas0 + (io-1))
      enddo

      end
c============================================================


c===============================================================
c===============================================================
      subroutine wfat_get_beta(rtdb, bt_start, bt_last, nbt, dbt)

      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb
      double precision bt_start
      double precision bt_last
      integer nbt
      double precision dbt


      if (.not. rtdb_get(rtdb, 'wfat:beta_start', mt_dbl, 1, bt_start))
     &     call errquit('WFAT: Failed reading beta_start from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:beta_last', mt_dbl, 1, bt_last))
     &     call errquit('WFAT: Failed reading beta_last from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:nbeta', mt_int, 1, nbt))
     &     call errquit('WFAT: Failed reading nbeta from rtdb',
     &     1, INPUT_ERR)

      if (nbt > 1) then
         dbt = (bt_last - bt_start) / (nbt - 1)
      else
         dbt = 0.0d0
      endif

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_get_gamma(rtdb, gm_start, gm_last, ngm, dgm)

      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
      integer rtdb
      double precision gm_start
      double precision gm_last
      integer ngm
      double precision dgm
      

      if (.not. rtdb_get(rtdb, 'wfat:gamma_start', mt_dbl, 1, gm_start))
     &     call errquit('WFAT: Failed reading gamma_start from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:gamma_last', mt_dbl, 1, gm_last))
     &     call errquit('WFAT: Failed reading gamma_last from rtdb',
     &     1, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:ngamma', mt_int, 1, ngm))
     &     call errquit('WFAT: Failed reading ngamma from rtdb',
     &     1, INPUT_ERR)

      if (ngm > 1) then
         dgm = (gm_last - gm_start) / (ngm - 1)
      else
         dgm = 0.0d0
      endif
      
      end
c===============================================================
      
      
c===============================================================
c===============================================================
      subroutine wfat_prepare_an(nb, nat, nch, nqn, ndgi, l_dgi, k_dgi,
     &           ndgn, l_dgn, k_dgn, nsc, l_sc, k_sc, nsx, l_sx, k_sx,
     &           nmmd, l_mmd, k_mmd, nlpb, l_lpb, k_lpb,
     &           n_an, l_an, k_an)
      
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nb             ! [input], 
      integer nat            ! [input], 
      integer nch            ! [input], 
      integer nqn            ! [input], 
      integer ndgi           ! [output],
      integer l_dgi, k_dgi   ! [output],
      integer ndgn           ! [output],
      integer l_dgn, k_dgn   ! [output],
      integer nsc            ! [output],
      integer l_sc, k_sc     ! [output],
      integer nsx            ! [output],
      integer l_sx, k_sx     ! [output],
      integer nmmd           ! [output],
      integer l_mmd, k_mmd   ! [output],
      integer nlpb           ! [output],
      integer l_lpb, k_lpb   ! [output],
      integer n_an           ! [output],
      integer l_an, k_an     ! [output],

      ndgi = nch
      if (.not. MA_alloc_get(MT_INT, ndgi, 'dgid0', l_dgi, k_dgi))
     &     call errquit('WFAT: dgid0 allocation failed', 2, MA_ERR)
      ndgn = nch
      if (.not. MA_alloc_get(MT_INT, ndgn, 'dgn0', l_dgn, k_dgn))
     &     call errquit('WFAT: dgn0 allocation failed', 2, MA_ERR)
      nsc = 2*nb**2 + nb
      if (.not. MA_alloc_get(MT_DBL, nsc, 'scratch', l_sc, k_sc))
     &     call errquit('WFAT: scratch allocation failed', 2, MA_ERR)
      nsx = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nsx, 'scr-cpx', l_sx, k_sx))
     &     call errquit('WFAT: scr-cpx allocation failed', 2, MA_ERR)
      nmmd = 3 * nch
      if (.not. MA_alloc_get(MT_DBL, nmmd, 'chmmtd', l_mmd, k_mmd))
     &     call errquit('WFAT: chmmtd allocation failed', 2, MA_ERR)
      nlpb = nch
      if (.not. MA_alloc_get(MT_DBL, nlpb, 'lfpb', l_lpb, k_lpb))
     &     call errquit('WFAT: lfpb allocation failed', 2, MA_ERR)
      n_an = nat * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, n_an, 'an', l_an, k_an))
     &     call errquit('WFAT: an allocation failed', 2, MA_ERR)

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_erase_an(l_dgi, l_dgn, l_sc, l_sx, l_mmd, l_lpb,
     &                         l_an)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer l_dgi
      integer l_dgn
      integer l_sc
      integer l_sx
      integer l_mmd
      integer l_lpb
      integer l_an

      if (.not. MA_free_heap(l_dgi))
     &     call errquit('WFAT: dgid0 deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_dgn))
     &     call errquit('WFAT: dgn0 deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_sc))
     &     call errquit('WFAT: scratch deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_sx))
     &     call errquit('WFAT: scr-cpx deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_mmd))
     &     call errquit('WFAT: chmmtd deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_lpb))
     &     call errquit('WFAT: lfpb deallocation failed', 3, MA_ERR)
      if (.not. MA_free_heap(l_an))
     &     call errquit('WFAT: an deallocation failed', 3, MA_ERR)

      end
c===============================================================
      
      
