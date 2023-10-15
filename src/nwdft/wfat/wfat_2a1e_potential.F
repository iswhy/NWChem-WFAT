c===============================================================
c===============================================================
      subroutine wfat_2a1e_potential(geom, vpt, eps, pot)

      implicit none
#include "errquit.fh"
#include "geom.fh"
      integer geom
      double precision vpt(3)
      double precision eps
      double precision pot

      character*16 t
      integer i, nc
      double precision dr(3), atp(3), q, dist
      

      if (.not. geom_ncent(geom, nc))
     &     call errquit('WFAT: Cannot get ncenters.', 28, GEOM_ERR)

      pot = 0.0d0
      do i = 1, nc
         if (.not. geom_cent_get(geom, i, t, atp, q))
     &        call errquit('WFAT: geom center info failed.', 45,
     &        GEOM_ERR)
         dr = vpt - atp
         dist = sqrt(sum(dr*dr) + eps)
         pot = pot + (-q)/dist
      enddo

      end
c===============================================================
      
      
c===============================================================
c===============================================================
      subroutine wfat_2a1e_vcpsi(geom, sphi, m, nr, nt, ke, eps, zi, 
     &     vcf, chpos, vpt, scr, vc_psi)

      implicit none
#include "stdio.fh"
      integer geom
      integer sphi
      integer m
      integer nr
      integer nt
      double precision ke
      double precision eps
      double precision zi
      double precision vcf(nr*nt)
      double precision chpos(3)
      double precision vpt(3)
      double precision scr(nr+nt+nr*nt)
      double precision vc_psi

      integer pr, pt, pb
      double precision r, t, theta, deg2rad, pot, dr(3), s, a, b, kp, sc
      double precision pi_const, phi
      parameter(deg2rad=0.01745329251994d0, pi_const=3.14159265359d0)


c     pr = pointer to radial DVR basis
c     pt = pointer to polar DVR basis
c     pb = pointer to the global basis
      pr = 1
      pt = 1 + nr
      pb = 1 + nr + nt

      
c     The radial and polar DVR bases are defined in the molecular frame.
c     This is also the frame where the nuclei's positions are defined.
c     Since vpt contains the coordinate of the evaluation point in the
c     molecular frame, r and theta can directly be extracted from vpt
c     without prior shifting.
      a = 1.0d0
      b = ke * 10.0d0
      kp = ke

c     ==== Radial part ====
      r = sqrt(sum(vpt*vpt))
      call wfat_eval_dvrr_mask(m, nr, ke, r, a, b, zi, kp, scr(pr))

c     ==== Theta part ====
      t = vpt(3) / r
      theta = acos(t) / deg2rad
      if (theta > 180.0d0 + 1.0e-12)  theta = 180.0d0
      if (theta <   0.0d0 - 1.0e-12)  theta =   0.0d0
      call wfat_eval_dvrt(m, nt, theta, scr(pt))

c     ==== The (r,\theta) part of the wave function ====
      call wfat_kron_prod(nr, 1, nt, 1, scr(pr), scr(pt), scr(pb))

c     ==== Phi part ====
      phi = atan2(vpt(2), vpt(1))
      if (sphi == 1) then
         if (m == 0) then
            sc = 1.0d0 / sqrt(2*pi_const)
         else
            sc = cos(abs(m)*phi) / sqrt(pi_const)
         endif
      else if (sphi == 2) then
         sc = sin(abs(m)*phi) / sqrt(pi_const)
      endif
      
      call wfat_2a1e_potential(geom, vpt, eps, pot)
      dr = vpt - chpos
      s = zi / sqrt(sum(dr*dr))
      vc_psi = -(pot + s) * sum(vcf(:) * scr(pb:pb+nr*nt-1)) * sc

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_vcpsi_ord1()

      implicit none


      end
c===============================================================
      
