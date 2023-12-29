c=======================================================
c=======================================================
      subroutine wfat_calc_ti_nopw()

      implicit none


c     ==== Obtaining info on beta ====
      call wfat_get_beta(rtdb, bt_start, bt_last, nbt, dbt)
      dbtr = dbt * deg2rad

c     ==== Obtaining info on gamma ====
      call wfat_get_gamma(rtdb, gm_start, gm_last, ngm, dgm)
      dgmr = dgm * deg2rad

      call wfat_get_gridrange(geom, ndid, gq, int_mb(k_icp),
     &     int_mb(k_ipp), natp, nptp)

      call wfat_collect_mfgrid(geom, gq, natp, nptp, int_mb(k_icp),
     &     int_mb(k_ipp), dbl_mb(k_pt), dbl_mb(k_wt))
      
      call wfat_collect_vcpsi(geom, bas, nb, ns, nch, nptp, vch,
     &     spinid, chpos, dbl_mb(k_pt), dbl_mb(k_rd), zi, smh,
     &     ismask_e, ergmask, amsk, bmsk, dbl_mb(k_vc))
      
      call wfat_collect_celf(geom, natp, nptp, int_mb(k_icp), 
     &     int_mb(k_ipp), celm, dbl_mb(k_pt), dbl_mb(k_cl))

      if (ismask)
     &     call wfat_collect_maskf(geom, natom, nptp, dbl_mb(k_pt),
     &     mask_param, dbl_mb(k_ma))
      

      do ig = 1, ngm
         gm = gm_start + dgm*(ig-1)

         do ib = 1, nbt
            bt = bt_start + dbt*(ib-1)

            if (pr_srfc .or. pr_rate) then
               call wfat_collect_tau0(geom, nch, nqn, chpos, erg, nid, 
     &              mid, bt, gm, nptp, zi, dbl_mb(k_pt), dcpl_mb(k_tu0))
               if (wfat1) then
c                  call wfat_collect_tau1()
               endif            
                                             

               if (wfat1) then
c                  call wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
c     &                 natom, int_mb(k_icp), int_mb(k_ipp), 
c     &                 dbl_mb(k_vc), dcpl_mb(k_tu0), dbl_mb(k_cl), 
c     &                 dbl_mb(k_wt), ismask, dbl_mb(k_ma),
c     &                 dcpl_mb(k_acg))
c                  call wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
c     &                 natom, int_mb(k_icp), int_mb(k_ipp), 
c     &                 dbl_mb(k_cvc), dcpl_mb(k_tu0), dbl_mb(k_cl), 
c     &                 dbl_mb(k_wt), ismask, dbl_mb(k_ma),
c     &                 dcpl_mb(k_acj))
c                  call wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
c     &                 natom, int_mb(k_icp), int_mb(k_ipp), 
c     &                 dbl_mb(k_vc), dcpl_mb(k_tu1), dbl_mb(k_cl), 
c     &                 dbl_mb(k_wt), ismask, dbl_mb(k_ma),
c     &                 dcpl_mb(k_ack))
c                  call wfat_get_srfc1(bt, gm, nat, nch, nqn, 0, mid, 
c     &                 erg, chmmt, dcpl_mb(k_acg), dcpl_mb(k_acj),
c     &                 dcpl_mb(k_ack), dcpl_mb(k_dum), nsub, maxn, dgs, 
c     &                 ndgs, dthr, dps, pbs, .true., dbl_mb(k_sc),
c     &                 dcpl_mb(k_sx), dbl_mb(k_pbd), dcpl_mb(k_sg),
c     &                 dcpl_mb(k_sj), dcpl_mb(k_sk))
               else
                  call wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
     &                 natom, int_mb(k_icp), int_mb(k_ipp), 
     &                 dbl_mb(k_vc), dcpl_mb(k_tu0), dbl_mb(k_cl), 
     &                 dbl_mb(k_wt), ismask, dbl_mb(k_ma),
     &                 dcpl_mb(k_acg))
                  call wfat_get_srfc0(geom, bt, gm, nch, erg, nqn, mid,
     &                 0, chmmt, dcpl_mb(k_acg), .true., dcpl_mb(k_sg))
               endif


c              ==== Print structure factor ====
               if (ndid == 0 .and. pr_srfc) then
                  call wfat_print_structuref(nat, bt, gm, nch, nqn,
     &                 dcpl_mb(k_sg), gap)
c                  if (wfat1) call wfat_print_srfc1(nat, bt, gm, nch, 
c     &                 nqn, dcpl_mb(k_sj), dcpl_mb(k_sk), gap)
               endif
            endif

            if (pr_rate) then
c              ==== Calculate ionization rate ====
               if (wfat1) then
c                  call wfat_calc_wfat1_rate(bt, gm, nat, nqn, nch, nid,
c     &                 mid, fieldm, erg, zi, dcpl_mb(k_sg),
c     &                 dcpl_mb(k_sj), dcpl_mb(k_sk), dbl_mb(k_on), 
c     &                 chmmt, dbl_mb(k_pbd), uffc, 
c     &                 dbl_mb(k_fc), dbl_mb(k_rt), rate)
               else
                  call wfat_calc_rate(nat, fieldm, nch, erg, nqn, nid,
     &                 mid, zi, dcpl_mb(k_sg), uffc, dbl_mb(k_fc),
     &                 dbl_mb(k_rt), rate)
               endif
               
c              ==== Print ionization rate ====
               if (ndid == 0 .and. pr_rate)
     &              call wfat_print_rate(bt, gm, nch, nqn, 
     &              dbl_mb(k_rt), rate, gap)
            endif
            
         enddo
      enddo
      

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_acf_nopw(bt, gm, nch, nqn, natp, nptp, 
     &     natom, centid, ptlim, vc_psi, tau0, celf, wt, ismask, maskf,
     &     acf)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      double precision bt                    ! [input], orientation angle in degree.
      double precision gm                    ! [input], orientation angle in degree.
      integer nch
      integer nqn
      integer natp
      integer nptp
      integer natom
      integer centid(natp)
      integer ptlim(2,natp)
      double precision vc_psi(nptp,nch)
      double complex tau0(nptp,nqn,nch)
      double precision celf(nptp)
      double precision wt(nptp)
      logical ismask
      double precision maskf(*)      ! [input], if ismask is .true., its dimension must be at least nptp.
      double complex acf(natom,nqn,nch)   ! [output]

      integer nac, nac_p, l_acp, k_acp, g_ac
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      

      nac = natom * nqn * nch
      if (.not. ga_create(mt_dcpl, nac, 1, 'acf', 0, 0, g_ac))
     &     call errquit('WFAT: Cannot create acf global array.',
     &     1, GA_ERR)
      call ga_zero(g_ac)

c     ==== Calculate the asymptotic coefficients for the ====
c     ====    atoms contained in the current process     ====
      nac_p = natp * nqn * nch
      if (.not. MA_alloc_get(MT_DCPL, nac_p, 'acf_p', l_acp, k_acp))
     &     call errquit('WFAT: acf_p allocation failed', 59, MA_ERR)
      call wfat_calc_acf_nopw_p(nch, nqn, natp, nptp, ptlim, 
     &     vc_psi, tau0, celf, wt, ismask, maskf, dcpl_mb(k_acp))

c     ==== Transfer the values of the asymptotic coefficients in ====
c     ====   the current process to their respective positions   ====
c     ====          in the array that holds all atoms            ====
      call wfat_assign_acf_nopw_p(nch, nqn, natp, natom, centid,
     &     dcpl_mb(k_acp), acf)

c     ==== Accumulate the asymptotic coefficients across all processes ====
      call ga_sync()
      call ga_acc(g_ac, 1, nac, 1, 1, acf, nac, 1.0d0)      
      call ga_sync()
      acf = czero
      call ga_get(g_ac, 1, nac, 1, 1, acf, nac)

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_acp))
     &     call errquit('WFAT: acf_p deallocation failed', 61, MA_ERR)
      if (.not. ga_destroy(g_ac))
     &     call errquit('WFAT: Cannot destroy acf.', 1, GA_ERR)
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_assign_acf_nopw_p(nch, nqn, natp, natom, centid,
     &     acf_p, acf)

      implicit none
      integer nch
      integer nqn
      integer natp
      integer natom
      integer centid(natp)
      double complex acf_p(natp,nqn,nch)
      double complex acf(natom,nqn,nch)

      integer i, j, icp, ic

      do i = 1, nch
         do j = 1, nqn
            do icp = 1, natp
               ic = centid(icp)
               acf(ic,j,i) = acf_p(icp,j,i)
            enddo
         enddo
      enddo

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_calc_acf_nopw_p(nch, nqn, natp, nptp, ptlim, 
     &     vc_psi, tau0, celf, wt, ismask, maskf, acf)

      implicit none
      integer nch
      integer nqn
      integer natp
      integer nptp
      integer ptlim(2,natp)
      double precision vc_psi(nptp,nch)
      double complex tau0(nptp,nqn,nch)
      double precision celf(nptp)
      double precision wt(nptp)
      logical ismask
      double precision maskf(*)          ! [input], if ismask is .true., its dimension must be at least nptp.
      double complex acf(natp,nqn,nch)  ! [output]

      integer i, j, ic, k, kk
      double precision mfun
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )


      acf = czero
      do i = 1, nch
         do j = 1, nqn

            kk = 0
            do ic = 1, natp
               do k = ptlim(1,ic), ptlim(2,ic)
                  kk = kk + 1
                  if (ismask) then
                     mfun = maskf(kk)
                  else
                     mfun = 1.0d0
                  endif
                  
                  acf(ic,j,i) = acf(ic,j,i) +
     &                          conjg(tau0(kk,j,i)) * vc_psi(kk,i) * 
     &                          mfun * celf(kk)
               enddo
            enddo
            
         enddo
      enddo

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_collect_celf(geom, natp, nptp, centid, ptlim,
     &     celm, v_m, celf)

      implicit none
#include "errquit.fh"
      integer geom
      integer natp
      integer nptp
      integer centid(natp)
      integer ptlim(2,natp)
      integer celm
      double precision v_m(3,nptp)
      double precision celf(nptp)

      logical t1
      integer icp, ic, i, j
      double precision v_c(3)

      j = 0
      do icp = 1, natp
         ic = centid(icp)
         do i = ptlim(1,icp), ptlim(2,icp)
            j = j + 1
            call geom_cent2cent(geom, 0, v_m(:,j), ic, v_c, t1)
            call wfat_celf_1pt(geom, ic, v_c, celm, celf(j))
         enddo
      enddo

      if (j /= nptp) call errquit('WFAT: The requirement j=nptp is '//
     &     'not satisfied. This is a bug, check the wfat_collect_celf'//
     &     ' subroutine file.', 99, UERR)
      
      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_collect_maskf(geom, natom, nptp, v_m, mask_param,
     &     maskf)

      implicit none
      integer geom
      integer natom
      integer nptp
      double precision v_m(3,nptp)
      double precision mask_param(2)
      double precision maskf(nptp)

      integer i
      double complex maskf_

      maskf_ = cmplx(1.0d0, 0.0d0, kind=kind(1.0d0))
      do i = 1, nptp
         call wfat_mask_intd(geom, natom, 1, 1, 1, v_m(:,i), 
     &        mask_param, maskf_)
         maskf(i) = real(maskf_)
      enddo

      end
c=======================================================
      
