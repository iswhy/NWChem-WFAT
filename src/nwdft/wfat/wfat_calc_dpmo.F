c=======================================================
c=======================================================
      subroutine wfat_calc_dpmo(ctr, dpmo)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "wfat_qm_param.fh"
#include "mafdecls.fh"
      double precision ctr(3)                          ! [input],
      double precision dpmo(nbas0, nbas0, nspin0, 3)   ! [output],

      integer NMMT
      parameter(NMMT=4)
      integer nmp, ncf, is, i, iq, l_mp, k_mp, i_mp, l_cf, k_cf, i_cf


c     ==== Dipole matrices in GTO basis relative to ctr ====
      nmp = nbas0 * nbas0 * NMMT
      if (.not. MA_alloc_get(MT_DBL, nmp, 'mpmat', l_mp, k_mp))
     &     call errquit('WFAT: mpmat allocation failed', 1,
     &     MA_ERR)
      call int_mpole_all(bas_hd, nbas0, 1, ctr, dbl_mb(k_mp))


c     ==== Get MO coeffs ====
      ncf = nbas0 * nbas0 * nspin0
      if (.not. MA_alloc_get(MT_DBL, ncf, 'vcf', l_cf, k_cf))
     &     call errquit('WFAT: vcf allocation failed', 1, MA_ERR)
      do is = 1, nspin0
         i_cf = k_cf + (is-1)*nbas0*nbas0
         call ga_get(g_mov(is), 1, nbas0, 1, nbas0, dbl_mb(i_cf), 
     &        nbas0)
      enddo
      

c     ==== Dipole matrices in MO basis relative to ctr ====
      do i = 1, 3
         iq = i + 1
         i_mp = k_mp + (iq-1)*nbas0*nbas0
         do is = 1, nspin0
            i_cf = k_cf + (is-1)*nbas0*nbas0
            call wfat_represent_sym(nbas0, nbas0, dbl_mb(i_cf),
     &           dbl_mb(i_mp), dpmo(:,:,is,i))
         enddo
      enddo


c     ==== Deallocations ====
      if (.not. MA_free_heap(l_mp))
     &     call errquit('WFAT: mpmat deallocation failed', 38, MA_ERR)
      if (.not. MA_free_heap(l_cf))
     &     call errquit('WFAT: vcf deallocation failed', 38, MA_ERR)
      
      end
c=======================================================
      
