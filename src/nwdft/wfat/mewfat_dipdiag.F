c===============================================================
c===============================================================
      subroutine mewfat_dipdiag(rtdb, nchar)

      implicit none
#include "inp.fh"      
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
      integer rtdb
      integer geom
      integer bas
      integer nchar

      character(len=nchar) s
      integer i1, i2, ndets
      parameter(max_ndets=10)

      integer g_mo(max_ndets)


c      write(luout, '("n dipdiag = ", i4)') nchar

      if (.not. rtdb_cget(rtdb, 'mewfat:dipdiag', 1, s))
     &     call errquit('MEWFAT: Cannot read dipdiag', 32, RTDB_ERR)

      i1 = 0
      l = inp_strtok(s, ' ', i1, i2)
      status = movecs_read_header(s(i1:i2), title, basname, scftype,
     &     nb, ns_, nmo, 2)
      ns = 2

      if (.not. ma_alloc_get(mt_dbl, nb, 'occ', l_oc, k_oc))
     &     call errquit('MEWFAT: Allocating occ failed', 14, MA_ERR)
      if (.not. ma_alloc_get(mt_dbl, nb, 'erg', l_eg, k_eg))
     &     call errquit('MEWFAT: Allocating erg failed', 15, MA_ERR)



      i = 0
      i1 = 0
      do while (inp_strtok(s, ' ', i1, i2))
         i = i + 1
         write(luout, '(a)') s(i1:i2)

         if (i > max_ndets) call errquit('MEWFAT: No. of degenerate '//
     &        'Slater determinants exceeds its maximum.', 78, UERR)

         write(si, '(i2.2)') i
         do is = 1, ns
            if (is == 1) then
               nm = 'mo-'//si//'-a'
               if (.not. ga_create(mt_dbl, nb, nb, nm, nb, 0,
     &              g_mo(is,i)))
     &              call errquit("MEWFAT: couldnt create "//nm,
     &              0, GA_ERR)
            else if (is == 2) then
               nm = 'mo_'//si//'-b'
               if (.not. ga_duplicate(g_mo(1,i), g_mo(2,i), nm))
     &              call errquit('MEWFAT: Couldnt create '//nm,
     &              0, GA_ERR)
            endif

            call ga_zero(g_mo(is,i))
         enddo


         do is = 1, ns
            if (ns_ == 2 .or. is == 1) then
               if (.not. movecs_read(s(i1:i2), is, dbl_mb(k_oc),
     $              dbl_mb(k_eg), g_mo(is,i)))
     $              call errquit("WFAT: Failed to read movecs data",
     &              0, UERR)
               if (i == 1) nel(is) = sum(dbl_mb(k_oc:k_oc+nb-1))
            else if (ns_ == 1 .and. is == 2) then
               call ga_copy(g_mo(1,i), g_mo(2,i))
               if (i == 1) then
                  nel(1) = nel(1)/2
                  nel(2) = nel(1)
               endif
            endif
         enddo
                  
      enddo

      ndets = i

      
      if (.not. MA_alloc_get(MT_DBL, nb*nb*2, 'movec1', l_mo1, k_mo1))
     &     call errquit('WFAT: movec1 allocation failed', 84, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nb*nb*2, 'movec2', l_mo2, k_mo2))
     &     call errquit('WFAT: movec2 allocation failed', 84, MA_ERR)
      
      
c     ==== Dipole matrices relative to origin ====
      nmp = nb * nb * NMMT     !! 2)
      if (.not. MA_alloc_get(MT_DBL, nmp, 'mpmat', l_mp, k_mp))
     &     call errquit('WFAT: mpmat allocation failed', 41, MA_ERR)
      orig = 0.0d0
      call int_mpole_all(bas, nb, 1, orig, dbl_mb(k_mp))

      
      do is = 1, 2
         g_ovl(is) = ga_create_atom_blocked(geom, bas, 'MO_overlap')
      enddo
      if (.not. ma_alloc_get(mt_dbl, nb*nb*2, 'moov', l_ovl, k_ovl))
     &     call errquit('WFAT: Cannot allocate moov', 1, MA_ERR)
      
      do j = 1, ndets         ! columns of dipole matrix
         do i = 1, j             ! rows of dipole matrix
c           ==== Overlap matrix between the two MO sets ====
            do is = 1, 2
               call ga_zero(g_ovl(is))
            enddo
            call mewfat_calc_moovl(geom, bas, oskel, g_mo(:,j), 
     &           g_mo(:,i), g_ovl)
            do is = 1, 2
               i_ovl = k_ovl + (is-1)*nb*nb
               call ga_get(g_ovl(is), 1, nb, 1, nb, dbl_mb(i_ovl), nb)

               i_mo1 = k_mo1 + (is-1)*nb*nb
               call ga_get(g_mo(is,i), 1, nb, 1, nb, dbl_mb(i_mo1), nb)
               i_mo2 = k_mo2 + (is-1)*nb*nb
               call ga_get(g_mo(is,j), 1, nb, 1, nb, dbl_mb(i_mo2), nb)

               call mewfat_dipslater()
            enddo

            

            
         enddo
      enddo

      
      
      do is = 1, 2
         if (.not. ga_destroy(g_ovl(is)))
     &        call errquit('WFAT: Cannot destroy MO_overlap.',
     &        1, GA_ERR)
      enddo
      if (.not. ma_free_heap(l_ovl))
     &     call errquit('WFAT: Cannot deallocate moov', 64, MA_ERR)
      


      do i = 1, ndets
         do is = 1, ns
            call ga_destroy(g_mo(is,i))
         enddo
      enddo
      
      end
c===============================================================
