c============================================================
c============================================================
      subroutine wfat_print_integrals(nch, nqn, nat, lmax, nid, mid,
     &     wfat1, intl, intl_j, intl_k1, intl_k2)

      implicit none
#include "stdio.fh"
      integer nch
      integer nqn
      integer nat
      integer lmax
      integer nid(nqn)
      integer mid(nqn)
      logical wfat1
      double complex intl((1+lmax)**2, nat, nqn, nch)       ! [input],
      double complex intl_j((1+lmax)**2, nat, nqn, 3, nch)  ! [input],   only used when wfat1=.true.
      double complex intl_k1((1+lmax)**2, nat, nqn, nch)    ! [input],   only used when wfat1=.true.
      double complex intl_k2((1+lmax)**2, nat, nqn, nch)    ! [input],   only used when wfat1=.true.

      integer i, iqn, l, m, i_sh, j


      write(luout, '("WFAT Integrals")')

      write(luout, '(5a7)', advance='no') 'chan.', 'n', 'm_', 'l', 'm'
      write(luout, '(3x, 2a10)', advance='no')
     &     'I (Re.)', 'I (Im.)'

      if (wfat1) then
         write(luout, ('(3x, 2a10)'), advance='no')
     &        'J_x (Re.)', 'J_x (Im.)'
         write(luout, ('(3x, 2a10)'), advance='no') 
     &        'J_y (Re.)', 'J_y (Im.)'
         write(luout, ('(3x, 2a10)'), advance='no') 
     &        'J_z (Re.)', 'J_z (Im.)'
         write(luout, ('(3x, 2a10)'), advance='no') 
     &        'K_1 (Re.)', 'K_1 (Im.)'
         write(luout, ('(3x, 2a10)'), advance='no') 
     &        'K_2 (Re.)', 'K_2 (Im.)'
      endif
      write(luout, '(" ")')
      
      
      do i = 1, nch
         do iqn = 1, nqn

            do l = abs(mid(iqn)), lmax
               do m = -l, l
                  i_sh = l**2 + l + m + 1
                  write(luout, '(5i7)', advance='no')
     &                 i, nid(iqn), mid(iqn), l, m
                  write(luout, '(3x, 2f10.6)', advance='no')
     &                 sum(intl(i_sh, :, iqn, i))

                  if (wfat1) then
                     do j = 1, 3
                        write(luout, '(3x, 2f10.6)', advance='no')
     &                       sum(intl_j(i_sh, :, iqn, j, i))
                     enddo
                     write(luout, '(3x, 2f10.6)', advance='no')
     &                    sum(intl_k1(i_sh, :, iqn, i))
                     write(luout, '(3x, 2f10.6)', advance='no')
     &                    sum(intl_k2(i_sh, :, iqn, i))
                  endif

                  write(luout, '(" ")')        ! new line
               enddo
               write(luout, '(" ")')    ! gap between different l's
            enddo
            write(luout, '(" ")')       ! gap between different (i,iqn)
            write(luout, '(" ")')       ! gap between different (i,iqn)

         enddo
      enddo

      end
c============================================================
