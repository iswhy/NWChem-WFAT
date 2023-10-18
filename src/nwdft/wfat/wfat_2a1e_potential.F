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

      character*1 sp
      integer pr, pt, pb
      double precision r, t, theta, deg2rad, pot, dr(3), s, a, b, kp, sc
      double precision phi
      parameter(deg2rad=0.01745329251994d0)


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
      if (sphi == 1) sp = 'c'
      if (sphi == 2) sp = 's'
      call wfat_2a1e_eval_phi(sp, m, phi, sc)
      
      
      call wfat_2a1e_potential(geom, vpt, eps, pot)
      dr = vpt - chpos
      s = zi / sqrt(sum(dr*dr))
      vc_psi = -(pot + s) * sum(vcf(:) * scr(pb:pb+nr*nt-1)) * sc

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_2a1e_vcpsi_ord1(geom, sphi, m, nr, nt, nr_m, nt_m,
     &     nr_p, nt_p, ke, eps, zi, vcf, cfx, cfy, cfz, chpos, vpt, scr,
     &     vc_psi, vc_cpsi)

      implicit none
      integer geom                    ! [input]
      integer sphi                    ! [input]
      integer m                       ! [input]
      integer nr                      ! [input]
      integer nt                      ! [input]
      integer nr_m                    ! [input]
      integer nt_m                    ! [input]
      integer nr_p                    ! [input]
      integer nt_p                    ! [input]
      double precision ke             ! [input]
      double precision eps            ! [input]
      double precision zi             ! [input]
      double precision vcf(nr*nt)     ! [input]
      double precision cfx(1:*)       ! [input]
      double precision cfy(1:*)       ! [input]
      double precision cfz(nr*nt)     ! [input]
      double precision chpos(3)       ! [input]
      double precision vpt(3)         ! [input]
      double precision scr(1:*)       ! [output],   1)
      double precision vc_psi         ! [output]
      double precision vc_cpsi(3)     ! [output]

      character*1 sp, sp_
      integer pr, pt, pb, pr_m, pt_m, pb_m, pr_p, pt_p, pb_p, n_m, n_p
      double precision r, t, theta, deg2rad, pot, dr(3), s, a, b, kp, sc
      double precision phi, scx_m, scx_p, scy_m, scy_p, basx_m, basx_p
      double precision basy_m, basy_p, basz
      parameter(deg2rad=0.01745329251994d0)

c     1) If m = 0, the minimum size of scr is nr + nt + nr*nt + nr_p + nt_p + nr_p*nt_p.
c        If m >= 1, the minimum size of scr is the above value for m = 0 plus 
c        nr_m + nt_m + nr_m*nt_m.

      
      if (m >= 1) then
         n_m = nr_m * nt_m
      else if (m == 0) then
         n_m = 0
      endif
      n_p = nr_p * nt_p


c     pr = pointer to radial DVR basis of the m spectrum
c     pt = pointer to polar DVR basis of the m spectrum
c     pb = pointer to the global basis of the m spectrum
      pr = 1
      pt = pr + nr
      pb = pt + nt
      
c     pr_m = pointer to radial DVR basis of the m-1 spectrum
c     pt_m = pointer to polar DVR basis of the m-1 spectrum
c     pb_m = pointer to the global basis of the m-1 spectrum
      pr_m = pb + nr*nt
      pt_m = pr_m + nr_m
      pb_m = pt_m + nt_m

c     pr_p = pointer to radial DVR basis of the m+1 spectrum
c     pt_p = pointer to polar DVR basis of the m+1 spectrum
c     pb_p = pointer to the global basis of the m+1 spectrum
      pr_p = pb_m + nr_m*nt_m
      pt_p = pr_p + nr_p
      pb_p = pt_p + nt_p

      
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
      if (m >= 1)
     &     call wfat_eval_dvrr_mask(m-1, nr_m, ke, r, a, b, zi, kp,
     &     scr(pr_m))
      call wfat_eval_dvrr_mask(m+1, nr_p, ke, r, a, b, zi, kp,
     &     scr(pr_p))
      

c     ==== Theta part ====
      t = vpt(3) / r
      theta = acos(t) / deg2rad
      if (theta > 180.0d0 + 1.0e-12)  theta = 180.0d0
      if (theta <   0.0d0 - 1.0e-12)  theta =   0.0d0
      call wfat_eval_dvrt(m, nt, theta, scr(pt))
      if (m >= 1)
     &     call wfat_eval_dvrt(m-1, nt_m, theta, scr(pt_m))
      call wfat_eval_dvrt(m+1, nt_p, theta, scr(pt_p))

      
c     ==== The (r,\theta) part of the wave function ====
      call wfat_kron_prod(nr, 1, nt, 1, scr(pr), scr(pt), scr(pb))
      if (m >= 1)
     &     call wfat_kron_prod(nr_m, 1, nt_m, 1, scr(pr_m), scr(pt_m),
     &     scr(pb_m))
      call wfat_kron_prod(nr_p, 1, nt_p, 1, scr(pr_p), scr(pt_p),
     &     scr(pb_p))

      
c     ==== Phi part ====
      phi = atan2(vpt(2), vpt(1))
      if (sphi == 1) then
         sp = 'c'
         sp_ = 's'
      else if (sphi == 2) then
         sp = 's'
         sp_ = 'c'
      endif
      call wfat_2a1e_eval_phi(sp, m, phi, sc)
      if (m >= 1) call wfat_2a1e_eval_phi(sp, m-1, phi, scx_m)
      call wfat_2a1e_eval_phi(sp, m+1, phi, scx_p)
      if (m >= 1) call wfat_2a1e_eval_phi(sp_, m-1, phi, scy_m)
      call wfat_2a1e_eval_phi(sp_, m+1, phi, scy_p)
      

c     ==== The potential term ====
      call wfat_2a1e_potential(geom, vpt, eps, pot)
      dr = vpt - chpos
      s = zi / sqrt(sum(dr*dr))


c     ==== 0th order wave function ====
      vc_psi = -(pot + s) * sum(vcf(:) * scr(pb:pb+nr*nt-1)) * sc

      
c     ==== 1st order wave function corrections (x,y,z) ====
      if (m >= 1) then
         basx_m = sum(cfx(1:n_m) * scr(pb_m:pb_m+n_m-1)) * scx_m
         basy_m = sum(cfy(1:n_m) * scr(pb_m:pb_m+n_m-1)) * scy_m
      else if (m == 0) then
         basx_m = 0.0d0
         basy_m = 0.0d0
      endif
      basx_p = sum(cfx(n_m+1:n_m+n_p) * scr(pb_p:pb_p+n_p-1)) * scx_p
      basy_p = sum(cfy(n_m+1:n_m+n_p) * scr(pb_p:pb_p+n_p-1)) * scy_p
      basz = sum(cfz(:) * scr(pb:pb+nr*nt-1)) * sc
      vc_cpsi(1) = -(pot + s) * (basx_m + basx_p)
      vc_cpsi(2) = -(pot + s) * (basy_m + basy_p)
      vc_cpsi(3) = -(pot + s) * basz
      
      end
c===============================================================
      
