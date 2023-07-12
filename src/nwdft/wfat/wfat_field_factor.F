c=======================================================
c=======================================================
      subroutine wfat_field_factor(fmagn, orberg_, zion,
     &                             nchan, mchan, fifc)
      
      implicit none
#include "stdio.fh"
      integer nchan, mchan
      double precision zion, fmagn, orberg_
      double precision fifc

      double precision kp

      if (fmagn > 1.0d-13) then
         kp = sqrt(2 * abs(orberg_))
         fifc = kp/2 *
     &        (4*kp**2/fmagn)**(2*zion/kp - 2*nchan - abs(mchan) - 1) *
     &        exp(-2*kp**3/3/fmagn)
      else
         fifc = 0.0d0
      endif

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_field_factor(ef, nch, nqn, erg, nid, mid,
     &     zi)

      implicit none
#include "stdio.fh"
      double precision ef
      integer nch
      integer nqn
      double precision erg(nch)
      integer nid(nqn)
      integer mid(nqn)
      double precision zi

      character*5 i_str
      character*15 e_str
      integer i, j, n, m
      double precision ff

      write(luout, '(2/)')
      call util_print_centered(luout, 'Field factors', 20,
     &     .true.)
      write(luout, '("")')

      write(luout, '(3x, "Field strength = ", f10.6, " a.u.")') ef
      write(luout, '("")')
      
      do i = 1, nch
         write(i_str, '(i5)') i
         write(e_str, '(f15.6)') erg(i)
         write(luout, '(3x, "Channel = ", a, ",", 2x, "energy = ", a)')
     &        trim(adjustl(i_str)), trim(adjustl(e_str))

         write(luout, '(5x, 2a4, a21)') 'n', 'm', 'field factor'
         call wfat_print_repeat(luout, 29, 5, '-')
         write(luout, '("")')
         do j = 1, nqn
            n = nid(j)
            m = mid(j)
            call wfat_field_factor(ef, erg(i), zi, n, m, ff)
            write(luout, '(5x, 2i4, es21.12)') n, m, ff
         enddo
         write(luout, '("")')
      enddo

      end
c=======================================================
      
