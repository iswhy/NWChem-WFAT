c===============================================================
c===============================================================
      subroutine wfat_2a1efd_vcpsi(geom, sphi, m, nfe, ndvr, maxn, xq, 
     &     wq, nr, nt, eps, zi, vcf, chpos, vpt, scr, vc_psi)

      implicit none
#include "stdio.fh"
      integer geom
      integer sphi
      integer m
      integer nfe
      integer ndvr(nfe)
      integer maxn
      double precision xq(maxn,nfe)
      double precision wq(maxn,nfe)
      integer nr         ! nr=sum(ndvr)-nfe*2+nfe-1
      integer nt
      double precision eps
      double precision zi
      double precision vcf(nr*nt)
      double precision chpos(3)
      double precision vpt(3)
      double precision scr(maxn**2+maxn+nr+nt+nr*nt)
      double precision vc_psi

      character*1 sp
      integer ptm, pr, pt, pb, i
      double precision r, t, theta, deg2rad, pot, dr(3), s, sc
      double precision phi
      parameter(deg2rad=0.01745329251994d0)

      
c     ptm = pointer to temp array for wfat_eval0_fedvr
c     pr = pointer to radial DVR basis
c     pt = pointer to polar DVR basis
c     pb = pointer to the global basis
      ptm = 1
      pr = ptm + maxn**2 + maxn
      pt = pr + nr
      pb = pt + nt

      
c     The radial and polar DVR bases are defined in the molecular frame.
c     This is also the frame where the nuclei's positions are defined.
c     Since vpt contains the coordinate of the evaluation point in the
c     molecular frame, r and theta can directly be extracted from vpt
c     without prior shifting.

c     ==== Radial part ====
      r = sqrt(sum(vpt*vpt))
      call wfat_eval0_fedvr(nfe, ndvr, maxn, xq, wq, nr, r, scr(ptm),
     &     scr(pr))
      scr(pr:pr+nr-1) = scr(pr:pr+nr-1) / r
      
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
      if (sphi == 1) sp = 'c'
      if (sphi == 2) sp = 's'
      call wfat_2a1e_eval_phi(sp, m, phi, sc)
      


      do i = 1, nr*nt
         if (abs(vcf(i)) < 1.0d-9) vcf(i) = 0.0d0
      enddo
      
      call wfat_2a1e_potential(geom, vpt, eps, pot)
      dr = vpt - chpos
      s = zi / sqrt(sum(dr*dr))
      vc_psi = -(pot + s) * sum(vcf(:) * scr(pb:pb+nr*nt-1)) * sc

      end
c===============================================================
