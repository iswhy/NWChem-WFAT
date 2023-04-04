c=======================================================
c=======================================================
      subroutine wfat_perturbnew_cwfn()

      implicit none
      integer nch
      integer spinid(nch)
      integer orbid(nch)
      double precision vch1(nbas0, 3, nch)

      double precision tmp1(nbas0)

      
c     ==== Calculate dipole matrix in the GTO basis ====
      do i = 1, 3
         g_dp(i) = ga_create_atom_blocked(geom_hd, bas_hd, dpnm(i))
         call ga_zero(g_dp(i))
      enddo
      call int_dip_ga(bas_hd, bas_hd, g_dp(1), g_dp(2), g_dp(3))


c     ==== Calculate dipole matrix in the MO basis ====
      do i = 1, 3
         do is = 1, nspin0
            if (.not. ga_duplicate(g_dp(i), g_dpmo(is,i), dpmonm(is,i)))
     &           call errquit('WFAT: Cannot create '//dpmonm(is,i)//'.',
     &           1, GA_ERR)
            call ga_zero(g_dpmo(is,i))
            call wfat_ga_vrepresent(g_mov(is), g_dp(i), g_dpmo(is,i))
         enddo
      enddo


      if (.not. MA_alloc_get(MT_DBL, nbas0**2*nspin0, 'mo',
     &     l_mo, k_mo)) 
     &     call errquit('WFAT: MO allocation failed', 19, MA_ERR)
      do is = 1, nspin0
         imo = k_mo + (is-1)*nbas0*nbas0
         call ga_get(g_mov(is), 1, nbas0, 1, nbas0, dbl_mb(imo), nbas0)
      endif
      
         
c     ==== Calculate the GTO coefficients of ====
c     ==== the first order correction terms  ====
      do i = 1, nch             ! Loop over ionized channels
         io = orbid(i)
         is = spinid(i)
         do k = 1, 3            ! Loop over dipole components
            call ga_get(g_dpmo(is,k), 1, nbas0, io, io, tmp1(1), nbas0)

            do j = 1, nbas0     ! Loop over GTO bases
               vch1(j,k,i) = 0.0d0
               do m = 1, nbas0  ! Sum over MO's
                  ieg0 = k_moerg + (is-1)*nbas0 + io - 1
                  iegm = k_moerg + (is-1)*nbas0 +  m - 1
                  derg = dbl_mb(ieg0) - dbl_mb(iegm)
               
                  if (abs(derg) > dgnthr) then
                     imo = k_mo + (is-1)*nbas0*nbas0 + (m-1)*nbas0 +
     &                     (j-1)
                     vch1(j,k,i) = vch1(j,k,i) +
     &                             dbl_mb(imo) * tmp1(m) / derg
                  endif
               enddo
            enddo
            
         enddo
      enddo
      

c     ==== Deallocations ====
      if (.not. MA_free_heap(l_mo))
     &     call errquit('WFAT: mo deallocation failed', 20, MA_ERR)
      do i = 1, 3
         if (.not. ga_destroy(g_dp(i)))
     &        call errquit('WFAT: Cannot destroy '//dpnm(i)//'',
     &        1, GA_ERR)         
         do is = 1, nspin0
            if (.not. ga_destroy(g_dp(i)))
     &           call errquit('WFAT: Cannot destroy '//dpmonm(i,is)//'',
     &           1, GA_ERR)
         enddo
      enddo

      end
c=======================================================
      
