c=======================================================
c=======================================================
c==== Imam ====
      subroutine wfat_eval_basis_1pt(ghan, bhan, maxder, vpt,
     &                bval, d1bval, d2bval, d3bval)

      implicit none
#include "wfat_qm_param.fh"      
#include "stdio.fh"      
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "errquit.fh"      
      integer ghan                  ! [input] geometry handle
      integer bhan                  ! [input] basis handle
      integer maxder                ! [input] highest derivative order requested
      double precision vpt(3)       ! [input] Cartesian coordinate of the evaluation point
      double precision bval(1:*)    ! [output] the value of basis at vpt
      double precision d1bval(3,*)  ! [output] the value of the 1st derivs of basis at vpt
      double precision d2bval(6,*)  ! [output] the value of the 2nd derivs of basis at vpt
      double precision d3bval(10,*) ! [output] the value of the 3rd derivs of basis at vpt
      
      integer ncont, natom, nbas, npoint
      integer lds, lac, laq, ltag, lrq, lnz, lbcent_info, lbcset_info
      integer ids, iac, iaq, itag, irq, inz, ibcent_info, ibcset_info
      integer ic
      double precision dx, dy, dz

      logical slas
      double precision rsl, kp
      
      parameter (npoint = 1)


c     ==== Filling docset ====
      if (.not. bas_numcont(bhan, ncont))
     &     call errquit('WFAT: Getting ncont failed', 54, BASIS_ERR)
      if (.not. MA_alloc_get(mt_log, ncont, 'docset', lds, ids))
     &     call errquit('WFAT: Cannot allocate docset', ncont, MA_ERR)
      log_mb(ids : ids+ncont-1) = .true.
      
c     ==== Get no. of atoms ====
      if (.not. geom_ncent(ghan, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)

c     ==== Get atomic coordinates and charges ====
      if (.not. MA_alloc_get(mt_dbl, 3*natom, 'a_coord', lac, iac))
     &     call errquit('WFAT: cannot allocate atom coord.',0, MA_ERR)
      if (.not. MA_alloc_get(mt_dbl, natom, 'a_charge', laq, iaq))
     &     call errquit('WFAT: cannot allocate atom charge',0, MA_ERR)
      if (.not. MA_alloc_get(mt_Byte, natom*16, 'a_tag', ltag, itag))
     &     call errquit('WFAT: cannot allocate atom tag',0, MA_ERR)
      if (.not. geom_cart_get(ghan, natom, byte_mb(itag),
     &     dbl_mb(iac), dbl_mb(iaq)))
     &     call errquit('WFAT: Fail to get atom coord and charges',
     &                  0, GEOM_ERR)

      
c     ==== Get the square of point-atom distances ====
      if(.not. MA_alloc_get(mt_dbl, natom, 'rq', lrq, irq))
     &     call errquit("WFAT: Cannot allocate rq",0, MA_ERR)
      do ic = 1, natom
         dx = dbl_mb(iac+3*(ic-1)+0) - vpt(1)
         dy = dbl_mb(iac+3*(ic-1)+1) - vpt(2)
         dz = dbl_mb(iac+3*(ic-1)+2) - vpt(3)
         dbl_mb(irq + npoint*(ic-1)) = dx*dx + dy*dy + dz*dz
      enddo      

c     ==== Filling niz ====
      if(.not. MA_alloc_get(MT_int, natom, 'niz', lnz, inz))
     &     call errquit("WFAT: Cannot allocate niz",0, MA_ERR)
      int_mb(inz : inz+natom-1) = 1

c     ==== Filling bas_cent_info and bas_cset_info ====
      if (.not. MA_alloc_get(mt_int, 3*natom, "bas_cent_info",
     &     lbcent_info, ibcent_info))
     &     call errquit("WFAT: cannot allocate bas_cent_info",0,
     &     MA_ERR)
      if (.not. MA_alloc_get(mt_int, 6*ncont, 'bas_cset_info',
     &     lbcset_info, ibcset_info))
     &     call errquit("WFAT: cannot allocate bas_cset_info",0,
     &     MA_ERR)
      call xc_make_basis_info(bhan, int_mb(ibcent_info),
     &     int_mb(ibcset_info), natom)      

      
c     ==== Calculating basis at the input point ====
      if (.not. bas_numbf(bhan, nbas))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)
c      call xc_eval_basis(bhan, maxder, 
c     &                   bval(1:nbas),
c     &                   d1bval(1:3, 1:nbas),
c     &                   d2bval(1:6, 1:nbas),
c     &                   d3bval(1:10, 1:nbas),      
c     &                   dbl_mb(irq),
c     &                   vpt, dbl_mb(iac), npoint, natom,
c     &                   int_mb(inz), log_mb(ids),
c     &                   int_mb(ibcent_info),
c     &                   int_mb(ibcset_info))

      slas = .false.
      rsl = 5.0d0
      kp = sqrt(2*abs(dbl_mb(k_moerg)))
      call wfat_eval_basis_slater(bhan, maxder, 
     &                   bval(1:nbas),
     &                   d1bval(1:3, 1:nbas),
     &                   d2bval(1:6, 1:nbas),
     &                   d3bval(1:10, 1:nbas),      
     &                   dbl_mb(irq),
     &                   vpt, dbl_mb(iac), npoint, natom,
     &                   int_mb(inz), log_mb(ids),
     &                   int_mb(ibcent_info),
     &                   int_mb(ibcset_info), slas, kp, rsl)

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(lds))
     &     call errquit('WFAT: Cannot deallocate docset', ncont, MA_ERR)
      if (.not. MA_free_heap(lac))
     &     call errquit('WFAT: cannot deallocate atom coord.',0, MA_ERR)
      if (.not. MA_free_heap(laq))
     &     call errquit('WFAT: cannot deallocate atom charge',0, MA_ERR)
      if (.not. MA_free_heap(ltag))
     &     call errquit('WFAT: cannot deallocate atom tag',0, MA_ERR)
      if (.not. MA_free_heap(lrq))
     &     call errquit("WFAT: cannot deallocate rq",0, MA_ERR)
      if (.not. MA_free_heap(lnz))
     &     call errquit("WFAT: cannot deallocate niz",0, MA_ERR)
      if (.not. MA_free_heap(lbcent_info))
     &     call errquit("WFAT: cannot deallocate bcent_info",0, MA_ERR)
      if (.not. MA_free_heap(lbcset_info))
     &     call errquit("WFAT: cannot deallocate bcset_info",0, MA_ERR)

      end
c=======================================================      


c=======================================================
c=======================================================
c==== Imam ====
      subroutine wfat_eval_basis_npt(ghan, bhan, maxder, npt, vpt,
     &                bval, d1bval, d2bval, d3bval)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "errquit.fh"      
      integer ghan                      ! [input] geometry handle
      integer bhan                      ! [input] basis handle
      integer maxder                    ! [input] highest derivative order requested
      integer npt                       ! [input] no. of points
      double precision vpt(3,npt)       ! [input] Cartesian coordinate of the evaluation point
      double precision bval(npt,*)      ! [output] the value of basis at vpt
      double precision d1bval(npt,3,*)  ! [output] the value of the 1st derivs of basis at vpt
      double precision d2bval(npt,6,*)  ! [output] the value of the 2nd derivs of basis at vpt
      double precision d3bval(npt,10,*) ! [output] the value of the 3rd derivs of basis at vpt

      character*16 atag
      integer ncont, natom, nbas
      integer lds, lac, lrq, lnz, lbcent_info, lbcset_info
      integer ids, iac, irq, inz, ibcent_info, ibcset_info
      integer ic, iac1, iac2, i
      double precision dr(3), acord(3), ach
      


c     ==== Filling docset ====
      if (.not. bas_numcont(bhan, ncont))
     &     call errquit('WFAT: Getting ncont failed', 54, BASIS_ERR)
      if (.not. MA_alloc_get(mt_log, ncont, 'docset', lds, ids))
     &     call errquit('WFAT: Cannot allocate docset', ncont, MA_ERR)
      log_mb(ids : ids+ncont-1) = .true.

      
c     ==== Get no. of atoms ====
      if (.not. geom_ncent(ghan, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)

      
c     ==== Get the square of point-atom distances ====
      if (.not. MA_alloc_get(mt_dbl, 3*natom, 'a_coord', lac, iac))
     &     call errquit('WFAT: cannot allocate atom coord.',0, MA_ERR)
      if (.not. MA_alloc_get(mt_dbl, npt*natom, 'rq', lrq, irq))
     &     call errquit("WFAT: Cannot allocate rq",0, MA_ERR)
      do ic = 1, natom
         if (.not. geom_cent_get(ghan, ic, atag, acord, ach))
     &        call errquit('WFAT: Cannot get center info', 1, GEOM_ERR)
         iac1 = iac + (ic-1)*3
         iac2 = iac1 + 2
         dbl_mb(iac1:iac2) = acord
         do i = 1, npt
            dr = acord - vpt(1:3, i)
            dbl_mb(irq + npt*(ic-1) + i - 1) = sum(dr*dr)
         enddo
      enddo
      

c     ==== Filling niz ====
      if(.not. MA_alloc_get(MT_int, natom, 'niz', lnz, inz))
     &     call errquit("WFAT: Cannot allocate niz",0, MA_ERR)
      int_mb(inz : inz+natom-1) = 1

      
c     ==== Filling bas_cent_info and bas_cset_info ====
      if (.not. MA_alloc_get(mt_int, 3*natom, "bas_cent_info",
     &     lbcent_info, ibcent_info))
     &     call errquit("WFAT: cannot allocate bas_cent_info",0,
     &     MA_ERR)
      if (.not. MA_alloc_get(mt_int, 6*ncont, 'bas_cset_info',
     &     lbcset_info, ibcset_info))
     &     call errquit("WFAT: cannot allocate bas_cset_info",0,
     &     MA_ERR)
      call xc_make_basis_info(bhan, int_mb(ibcent_info),
     &     int_mb(ibcset_info), natom)      

      
c     ==== Calculating basis at the input point ====
      if (.not. bas_numbf(bhan, nbas))
     &     call errquit('WFAT: Getting nbas failed', 53, BASIS_ERR)
      call xc_eval_basis(bhan, maxder, 
     &                   bval(1:npt, 1:nbas),
     &                   d1bval(1:npt, 1:3, 1:nbas),
     &                   d2bval(1:npt, 1:6, 1:nbas),
     &                   d3bval(1:npt, 1:10, 1:nbas),      
     &                   dbl_mb(irq),
     &                   vpt, dbl_mb(iac), npt, natom,
     &                   int_mb(inz), log_mb(ids),
     &                   int_mb(ibcent_info),
     &                   int_mb(ibcset_info))

c      call xc_eval_basis(bhan, maxder, 
c     &                   bval,
c     &                   d1bval,
c     &                   d2bval,
c     &                   d3bval,      
c     &                   dbl_mb(irq),
c     &                   vpt, dbl_mb(iac), npt, natom,
c     &                   int_mb(inz), log_mb(ids),
c     &                   int_mb(ibcent_info),
c     &                   int_mb(ibcset_info))

      
c     ==== Deallocations ====
      if (.not. MA_free_heap(lds))
     &     call errquit('WFAT: Cannot deallocate docset', ncont, MA_ERR)
      if (.not. MA_free_heap(lac))
     &     call errquit('WFAT: cannot deallocate atom coord.',0, MA_ERR)
      if (.not. MA_free_heap(lrq))
     &     call errquit("WFAT: cannot deallocate rq",0, MA_ERR)
      if (.not. MA_free_heap(lnz))
     &     call errquit("WFAT: cannot deallocate niz",0, MA_ERR)
      if (.not. MA_free_heap(lbcent_info))
     &     call errquit("WFAT: cannot deallocate bcent_info",0, MA_ERR)
      if (.not. MA_free_heap(lbcset_info))
     &     call errquit("WFAT: cannot deallocate bcset_info",0, MA_ERR)

      end
c=======================================================      
      
