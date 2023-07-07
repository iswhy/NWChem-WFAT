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
      double complex intl_j(1:*)     ! [input], only used when wfat1=.true.
      double complex intl_k1(1:*)    ! [input], only used when wfat1=.true.
      double complex intl_k2(1:*)    ! [input], only used when wfat1=.true.

      integer i, iqn, l, m, i_sh, j


      write(luout, '("WFAT Integrals")')      
      
      do i = 1, nch
         do iqn = 1, nqn

            write(luout, '(5x, "Channel = ", i3, ",", 4x, '//
     &           '"para. quantum numbers = ", 2i3)') i, nid(iqn),
     &           mid(iqn)

            write(luout, '(2a7)', advance='no') 'l', 'm'
            write(luout, '(3x, 2a10)', advance='no')
     &           'I (Re.)', 'I (Im.)'
            
            if (wfat1) then
               write(luout, ('(3x, 2a10)'), advance='no')
     &              'J_x (Re.)', 'J_x (Im.)'
               write(luout, ('(3x, 2a10)'), advance='no') 
     &              'J_y (Re.)', 'J_y (Im.)'
               write(luout, ('(3x, 2a10)'), advance='no') 
     &              'J_z (Re.)', 'J_z (Im.)'
               write(luout, ('(3x, 2a10)'), advance='no') 
     &              'K_1 (Re.)', 'K_1 (Im.)'
               write(luout, ('(3x, 2a10)'), advance='no') 
     &              'K_2 (Re.)', 'K_2 (Im.)'
            endif
            write(luout, '("")')


            do l = abs(mid(iqn)), lmax
               do m = -l, l
                  i_sh = l**2 + l + m + 1
                  write(luout, '(2i7)', advance='no') l, m
                  write(luout, '(3x, 2f10.6)', advance='no')
     &                 sum(intl(i_sh, :, iqn, i))

                  if (wfat1)
     &                 call wfat_print_wfat1_ints(nch, nqn, nat, lmax,
     &                 i_sh, iqn, i, intl_j, intl_k1, intl_k2)

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


c============================================================
c============================================================
      subroutine wfat_print_wfat1_ints(nch, nqn, nat, lmax, i_sh, i_qn,
     &     i_ch, intl_j, intl_k1, intl_k2)

      implicit none
#include "stdio.fh"
      integer nch
      integer nqn
      integer nat
      integer lmax
      integer i_sh
      integer i_qn
      integer i_ch
      double complex intl_j((1+lmax)**2, nat, nqn, 3, nch)  ! [input],   only used when wfat1=.true.
      double complex intl_k1((1+lmax)**2, nat, nqn, nch)    ! [input],   only used when wfat1=.true.
      double complex intl_k2((1+lmax)**2, nat, nqn, nch)    ! [input],   only used when wfat1=.true.

      integer j

      do j = 1, 3
         write(luout, '(3x, 2f10.6)', advance='no')
     &        sum(intl_j(i_sh, :, i_qn, j, i_ch))
      enddo
      write(luout, '(3x, 2f10.6)', advance='no')
     &     sum(intl_k1(i_sh, :, i_qn, i_ch))
      write(luout, '(3x, 2f10.6)', advance='no')
     &     sum(intl_k2(i_sh, :, i_qn, i_ch))

      end
c============================================================
      
