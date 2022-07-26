c=======================================================
c=======================================================
      subroutine wfat_field_factor(fmagn, orberg_, zion,
     &                             nchan, mchan, fifc)
      
      implicit none
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
