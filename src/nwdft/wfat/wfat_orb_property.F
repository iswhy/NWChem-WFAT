c=======================================================
c=======================================================
      subroutine wfat_orb_property(rtdb, nch, orbid, spinid, 
     &           vcf, erg, orbcent, orbmmt, iondip)

      implicit none
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "wfat_qm_param.fh"
#include "global.fh"
      integer rtdb
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision vcf(nbas0, nch)
      double precision erg(nch)
      double precision orbcent(3, nch)
      double precision orbmmt(3, nch)
      double precision iondip(3, nch)

      logical ismandip
      
            
      call wfat_get_coef(nch, orbid, spinid, vcf, erg)
      if (.not. rtdb_get(rtdb, 'wfat:man_dip', mt_log, 1, ismandip))
     &     ismandip = .false.

      if (ismandip) then
         call wfat_get_man_orbcent(rtdb, nch, orbcent, orbmmt, iondip)
      else
         call wfat_get_orbcent(rtdb, nch, orbid, spinid, 
     &        orbcent, orbmmt, iondip)
      endif

      
      if (ga_nodeid() == 0) then
         call wfat_print_orb_property(nch, orbid, spinid, erg,
     &        orbcent, orbmmt, iondip)
      endif      
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_print_orb_property(nch, orbid, spinid, erg,
     &           orbcent, orbmmt, iondip)

      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision erg(nch)
      double precision orbcent(3,nch)
      double precision orbmmt(3,nch)
      double precision iondip(3,nch)

      integer i, io, is
      

      write(luout, '(2/)')
      call util_print_centered(luout,
     &     'Properties of the user-requested ionization channels',
     &     40, .true.)
      write(luout, '(" ")')
      
      do i = 1, nch
         io = orbid(i)
         is = spinid(i)
         write(luout, '(5x, "Channel : ", i6)') i
         write(luout, '(7x, a17, " = ", i1)')
     &        "ionized spin     ", is
         write(luout, '(7x, a17, " = ", i6)')
     &        "orbital no.      ", io
         write(luout, '(7x, a17, " = ", es14.5, " Hartree")')
     &        "ionization energy", erg(i)
         write(luout, '(7x, a36, " = ", 3es14.5)')
     &        "channel position (xyz, Bohr)        ", orbcent(:, i)
         write(luout, '(7x, a36, " = ", 3es14.5)') 
     &        "unshifted channel moment (xyz, Bohr)",
     &        orbmmt(:, i) + orbcent(:, i)
         write(luout, '(7x, a36, " = ", 3es14.5)') 
     &        "shifted channel moment (xyz, Bohr)  ", orbmmt(:, i)
         write(luout, '(7x, a36, " = ", 3es14.5)')
     &        "shifted ion dipole (xyz, e.Bohr)    ", iondip(:, i)
      enddo

      write(luout, '(" ")')
      write(luout, '(" ")')
      
      end
c=======================================================
      

c=======================================================
c=======================================================            
      subroutine wfat_get_coef(nch, orbid, spinid, vcf, erg)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "wfat_qm_param.fh"      
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision vcf(nbas0,nch)
      double precision erg(nch)

      integer is, io, i, ieg


      do i = 1, nch
         io = orbid(i)
         is = spinid(i)

         if (is > nspin0)
     &        call errquit('WFAT: An element of argument spinid of '//
     &        'wfat_get_coef exceeds the number of spin channels.',
     &        1, UERR)
         call ga_get(g_mov(is), 1, nbas0, io, io, vcf(1,i), nbas0)

         ieg = k_moerg + (is-1)*nbas0 + io - 1
         erg(i) = dbl_mb(ieg)
      enddo

      end
c     Notes:
c       1) Note that this call to movecs_read(...) assumes that the basis set and thus
c          the number of basis functions is equal to nbas0 (argument of this routine).
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_get_channel_coef(nch, orbid, spinid, chcf)
c This routine is OLD.
      implicit none
#include "errquit.fh"
#include "stdio.fh"      
#include "global.fh"
c     #include "cscf.fh"
#include "wfat_qm_param.fh"
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision chcf(nbas0, nch)

      integer ndid
      integer i, is, io, it, nb, nn


      call ga_inquire(g_mov(1), it, nb, nn)
      if (nbas0 /= nb) then
         call errquit('WFAT: Argument nbas0 of wfat_get_channel_coef '//
     &        'must be equal to the number of basis.', 1, UERR)
      endif
      if (maxval(orbid) > nn .or. minval(orbid) < 1) then
         write(luout, '("In wfat_get_channel_coef, norb = ", i5)') nn
         call errquit('WFAT: The elements of argument array orbid '//
     &        'must be positive and less than norb (see above).',
     &        2, UERR)
      endif
      
      
      do i = 1, nch
         is = spinid(i)
         io = orbid(i)
         call ga_get(g_mov(is), 1, nbas0, io, io, chcf(1,i), nbas0)
      enddo

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_get_orbcent(rtdb, nch, orbid, spinid, 
     &           orbcent, orbmmt, iondip)

      implicit none
#include "stdio.fh"
#include "geom.fh"
#include "global.fh"
#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "wfat_qm_param.fh"
      integer rtdb
      integer nch
      integer orbid(nch)
      integer spinid(nch)
      double precision orbcent(3, nch)
      double precision orbmmt(3, nch)
      double precision iondip(3, nch)

      character*8 kst
      logical autopos
      integer ndid, nmp, NMMT, ncf, nom0, niondip
      integer l_mp, l_cf, l_om0, l_md
      integer k_mp, k_cf, k_om0, k_md
      integer is, icf, id, iq, imp, iom0, i, io, k, i_md
      parameter (NMMT = 4)
      double precision vorig(3), ndip0(NMMT), edip0(3), totdip0(3)
      double precision totzn, omgpos(3), ndip(3), edip(3), totdip(3)
      double precision totdip_(3)


      ndid = ga_nodeid()
      vorig = 0.0d0

c     ==== Nuclear moment relative to origin ====
      call geom_nuc_mpole(geom_hd, vorig, 1, ndip0, NMMT)

c     ==== Dipole matrices relative to origin ====
      nmp = nbas0 * nbas0 * NMMT     !! 2)
      if (.not. MA_alloc_get(MT_DBL, nmp, 'mpmat', l_mp, k_mp))
     &     call errquit('WFAT: mpmat allocation failed', 1, MA_ERR)
      call int_mpole_all(bas_hd, nbas0, 1, vorig, dbl_mb(k_mp))
      

c     ==== Get orb. coeffs to calculate electronic moments ====
      ncf = nbas0 * nbas0 * nspin0
      if (.not. MA_alloc_get(MT_DBL, ncf, 'vcf', l_cf, k_cf))
     &     call errquit('WFAT: vcf allocation failed', 1, MA_ERR)
      do is = 1, nspin0
         do io = 1, nbas0
            icf = k_cf + (is-1)*nbas0*nbas0 + (io-1)*nbas0
            call ga_get(g_mov(is), 1, nbas0, io, io,
     &                  dbl_mb(icf), nbas0)
         enddo
      enddo

c     ==== Electronic moments relative to origin ====
      nom0 = nbas0 * nspin0 * 3     !! 1)
      if (.not. MA_alloc_get(MT_DBL, nom0, 'orbmmt0', l_om0, k_om0))
     &     call errquit('WFAT: orbmmt0 allocation failed', 1, MA_ERR)

      do id = 1, 3
         iq = id + 1
         imp = k_mp + (iq-1)*nbas0*nbas0
         iom0 = k_om0 + (id-1)*nbas0*nspin0
         call wfat_ne_expval_1eopr(nbas0, nspin0, nel, dbl_mb(k_cf),
     &        dbl_mb(imp), dbl_mb(iom0), edip0(id))
      enddo
      
c     ==== Total dipole relative to origin ====
      totdip0 = ndip0(NMMT-2:NMMT) - edip0


c     ==== Orb. center and orb. moments relative to orb. center ====
      if (.not. rtdb_get(rtdb, 'wfat:autopos', mt_log, 1, autopos))
     &     call errquit('WFAT: Cannot get autopos from RTDB',
     &     1, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'wfat:omgpos', mt_dbl, 3, omgpos))
     &     call errquit('WFAT: Cannot get omgpos from RTDB',
     &     1, RTDB_ERR)
      if (.not. geom_nuc_charge(geom_hd, totzn))
     &     call errquit('WFAT: Getting totzn failed', 5, UERR)
      do i = 1, nch
         io = orbid(i)
         is = spinid(i)
         do id = 1, 3
c           ==== Compute the orb. center ====
            iom0 = k_om0 + nbas0*nspin0*(id-1) + nbas0*(is-1) + io - 1
            if (autopos) then
               orbcent(id,i) = (totdip0(id) + dbl_mb(iom0)) /
     &                         (totzn - nelec0 + 1)
            else
               orbcent(id,i) = omgpos(id)
            endif

c           ==== Compute orb. moments relative to the orb. center ====
            orbmmt(id,i) = dbl_mb(iom0) - orbcent(id,i)
         enddo
c        ==== Ion dip. moment (better be vanishing) ====
         ndip = ndip0(2:4) - totzn * orbcent(:,i)
         edip = edip0(1:3) - nelec0 * orbcent(:,i)
         totdip = ndip - edip
         iondip(:,i) = totdip + orbmmt(:,i)
      enddo


cOLD c     ==== For manually chosen dipole or moment ====
cOLD       if (rtdb_get(rtdb, 'wfat:niondip', mt_int, 1, niondip)) then
cOLD          if (niondip > nch) call errquit('WFAT: niondip cannot be '//
cOLD      &        'bigger than nch', 1, UERR)
cOLD          if (.not. ma_alloc_get(mt_dbl, niondip*3, 'iondip',
cOLD      &        l_md, k_md))
cOLD      &        call errquit('WFAT: Cannot allocate iondip', 1, MA_ERR)
cOLD 
cOLD          totdip_ = ndip0(2:4) - edip0(1:3)
cOLD          if (ndid==0) write(luout,
cOLD      &        '(" Neutral dipole = ", 3e14.5, " Bohr")') totdip_
cOLD          do k = 1, niondip
cOLD             i_md = k_md + (k-1)*3
cOLD             write(kst, '(i8)') k-1
cOLD             
cOLD cc type1            if (.not. rtdb_get(rtdb, 'wfat:iondip'//trim(adjustl(kst)),
cOLD cc type1     &           mt_dbl, 3, dbl_mb(i_md)))
cOLD cc type1     &           call errquit('WFAT: RTDB entry wfat:iondip'//
cOLD cc type1     &           trim(adjustl(kst))//' is missing', 1, UERR)
cOLD cc type1            iondip(:,k) = dbl_mb(i_md:i_md+2)
cOLD cc type1            orbcent(:,k) = iondip(:,k) / (totzn - nelec0 + 1)
cOLD cc type1            orbmmt(:,k) = -(totdip_ - iondip(:,k)) - orbcent(:,k)
cOLD cc type1            iondip(:,k) = iondip(:,k) -
cOLD cc type1     &                    orbcent(:,k) * (totzn - nelec0 + 1)      ! LHS should be zero.
cOLD 
cOLD             if (.not. rtdb_get(rtdb, 'wfat:orbmmt'//trim(adjustl(kst)),
cOLD      &           mt_dbl, 3, dbl_mb(i_md)))
cOLD      &           call errquit('WFAT: RTDB entry wfat:orbmmt'//
cOLD      &           trim(adjustl(kst))//' is missing', 1, UERR)
cOLD             orbmmt(:,k) = dbl_mb(i_md:i_md+2) - orbcent(:,k)
cOLD          enddo
cOLD 
cOLD          if (.not. ma_free_heap(l_md))
cOLD      &        call errquit('WFAT: Cannot deallocate iondip', 1, MA_ERR)
cOLD       endif
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_mp))
     &     call errquit('WFAT: mpmat deallocation failed', 1, MA_ERR)
      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: vcf deallocation failed', 1, MA_ERR)
      if (.not. MA_free_heap(l_om0))
     &     call errquit('WFAT: orbmmt0 deallocation failed', 1, MA_ERR)

      end
c     Notes:
c       1) orbmmt0 is conceptually a 3D array of size [nbas0, nspin0, 3].
c       2) mpmat is conceptually a 3D array of size [nbas0, nbas0, NMMT].
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_man_orbcent(rtdb, nch, orbcent, orbmmt,
     &     iondip)

      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      integer nch
      double precision orbcent(3, nch)
      double precision orbmmt(3, nch)
      double precision iondip(3, nch)

      integer i

      
      do i = 1, nch
         if (.not. rtdb_get(rtdb, 'wfat:man_orbcent', mt_dbl, 3,
     &        orbcent(:,i)))
     &        call errquit('WFAT: Cannot get man_orbcent from RTDB',
     &        1, RTDB_ERR)
         if (.not. rtdb_get(rtdb, 'wfat:man_orbmmt', mt_dbl, 3,
     &        orbmmt(:,i)))
     &        call errquit('WFAT: Cannot get man_orbmmt from RTDB',
     &        1, RTDB_ERR)
         iondip(:,i) = (/ 1.0d3, 1.0d3, 1.0d3 /)
      enddo

      end
c=======================================================
