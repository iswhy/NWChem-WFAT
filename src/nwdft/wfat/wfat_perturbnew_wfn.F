c=======================================================
c=======================================================
      subroutine wfat_perturbnew_cwfn(nch, spinid, orbid, vch1, pb)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"      
#include "wfat_qm_param.fh"
#include "global.fh"
      integer nch
      integer spinid(nch)
      integer orbid(nch)
      double precision vch1(nbas0, 3, nch)
      double precision pb(3,3,nch)

      integer i, is, imo, io, k, j, m, id0, idm, ieg0, iegm, i_dpo
      integer l_mo, k_mo, l_dpo, k_dpo
      double precision derg, ctr(3)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked


c     ==== Calculate dipole matrix in the MO basis ====      
      if (.not. MA_alloc_get(MT_DBL, nbas0**2*nspin0*3, 'dpmo',
     &     l_dpo, k_dpo)) 
     &     call errquit('WFAT: dpmo allocation failed', 19, MA_ERR)
      ctr = 0.0d0                                     !! 1)
      call wfat_calc_dpmo(ctr, dbl_mb(k_dpo))

      
c     ==== Compute the static polarizabilities in molecular frame ====
      call wfat_calc_polb(nch, orbid, spinid, dbl_mb(k_dpo), pb)


c     ==== MO coefficients ====
      if (.not. MA_alloc_get(MT_DBL, nbas0**2*nspin0, 'mo',
     &     l_mo, k_mo)) 
     &     call errquit('WFAT: MO allocation failed', 19, MA_ERR)
      do is = 1, nspin0
         imo = k_mo + (is-1)*nbas0*nbas0
         call ga_get(g_mov(is), 1, nbas0, 1, nbas0, dbl_mb(imo), nbas0)
      enddo
      
         
c     ==== Calculate the GTO coefficients of ====
c     ==== the first order correction terms  ====
      do i = 1, nch             ! Loop over ionized channels
         io = orbid(i)
         is = spinid(i)
         do k = 1, 3            ! Loop over dipole components

            do j = 1, nbas0     ! Loop over GTO bases
               vch1(j,k,i) = 0.0d0
               do m = 1, nbas0  ! Sum over MO's
                  id0 = int_mb(k_dgid + (is-1)*nbas0 + (io-1))
                  idm = int_mb(k_dgid + (is-1)*nbas0 + (m-1))
                  if (id0 /= idm) then
                     ieg0 = k_moerg + (is-1)*nbas0 + (io-1)
                     iegm = k_moerg + (is-1)*nbas0 + (m-1)
                     derg = dbl_mb(ieg0) - dbl_mb(iegm)
                  
                     imo = k_mo + (is-1)*nbas0*nbas0 + (m-1)*nbas0 +
     &                     (j-1)
                     i_dpo = k_dpo + nbas0**2*nspin0*(k-1) +
     &                       nbas0**2*(is-1) + nbas0*(io-1) + (m-1)
                     vch1(j,k,i) = vch1(j,k,i) +
     &                             dbl_mb(imo) * dbl_mb(i_dpo) / derg
                  endif
               enddo
            enddo
            
         enddo
      enddo
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_mo))
     &     call errquit('WFAT: mo deallocation failed', 20, MA_ERR)
      if (.not. MA_free_heap(l_dpo))
     &     call errquit('WFAT: dpmo deallocation failed', 20, MA_ERR)

      end
c     NOTES:
c     1) The center of the dipole moments does not matter for its use in computing
c        the first order perturbation correction because the center only affects the
c        diagonal elements of the dipole matrix in MO basis while for the aforementioned
c        purpose, we only need the non-diagonal elements.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_polb(nch, orbid, spinid, dpmo, pb)

      implicit none
#include "wfat_qm_param.fh"
#include "mafdecls.fh"
      integer nch                                      ! [input],
      integer orbid(nch)                               ! [input],
      integer spinid(nch)                              ! [input],
      double precision dpmo(nbas0, nbas0, nspin0, 3)   ! [input],
      double precision pb(3,3,nch)                     ! [output],

      integer i, io, is, j2, j1, m, id0, idm, ieg0, iegm
      double precision derg, pb0


      do i = 1, nch
         io = orbid(i)
         is = spinid(i)
         do j2 = 1, 3
            do j1 = 1, 3

               pb(j1,j2,i) = 0.0d0
               do m = 1, nbas0
                  id0 = int_mb(k_dgid + (is-1)*nbas0 + (io-1))
                  idm = int_mb(k_dgid + (is-1)*nbas0 + (m-1))
                  if (id0 /= idm) then
                     ieg0 = k_moerg + (is-1)*nbas0 + (io-1)
                     iegm = k_moerg + (is-1)*nbas0 + (m-1)
                     derg = dbl_mb(ieg0) - dbl_mb(iegm)

                     pb0 = dpmo(m,io,is,j1) * dpmo(m,io,is,j2) / derg
                     pb(j1,j2,i) = pb(j1,j2,i) + pb0
                  endif
               enddo

            enddo
         enddo
      enddo

      end
c=======================================================
      
      
