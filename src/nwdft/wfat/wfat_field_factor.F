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
      
      kp = sqrt(2 * abs(orberg_))
      fifc = kp/2 *
     &       (4*kp**2/fmagn)**(2*zion/kp - 2*nchan - abs(mchan) - 1) *
     &       exp(-2*kp**3/3/fmagn)

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

      integer i, j, n, m
      double precision ff

      write(luout, '("Field factors : ")')
      do i = 1, nch
         write(luout, '(2x, "Channel : ", i4)') i
         do j = 1, nqn
            n = nid(j)
            m = mid(j)
            call wfat_field_factor(ef, erg(i), zi, n, m, ff)
            write(luout, '(5x, 2i4, es21.12)') n, m, ff
         enddo
         write(luout, '(" ")')
      enddo

      end
c=======================================================
      
