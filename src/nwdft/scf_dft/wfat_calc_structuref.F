c=======================================================
c=======================================================
c==== Imam ====      
      subroutine wfat_calc_structuref(geom, bt, gm, nch, erg, nqn, mid, 
     &                                lmax_, chmmt, intgl, srfc)
      
      implicit none
c     #include "cdft.fh"
ccOLD #include "wfat_qm_param.fh"
#include "errquit.fh"
#include "geom.fh"
#include "stdio.fh"
      integer geom
      integer nch, nqn, mid(nqn), lmax_
      double precision bt, gm                    ! [input], orientation angle in degree.
      double precision erg(nch), chmmt(3,nch)
      double complex intgl(1:*)    !! 1)
      double complex srfc(1:*)     !! 2)
      
      integer ich, iqn, l, m, i_lm, icent, i, j, nat, nsh
      double precision btr, gmr, wgf, mu(3), kappa, protmat(3,3)
      double precision deg2rad
      double complex g_nu0, g_nu, czero, imnum
      parameter ( deg2rad = 0.01745329251994d0 )
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )

      
c     ==== Get number of atoms ====
      if (.not. geom_ncent(geom, nat))
     &     call errquit('WFAT: Getting nat failed', 77, GEOM_ERR)

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(-bt, -gm, protmat)
      btr = bt * deg2rad
      gmr = gm * deg2rad

      
c     ==== Begin calculating structure factor ====
      nsh = (1 + lmax_)**2

      
      do ich = 1, nch
         mu = matmul(transpose(protmat), -chmmt(:, ich))
         kappa = sqrt(2 * abs(erg(ich)))
         do iqn = 1, nqn
      
            do icent = 1, nat
               g_nu = czero
               do l = abs(mid(iqn)), lmax_
                  do m = -l, l
                     call wfat_wignerf(l, mid(iqn), m, btr, wgf)
c check                  call wfat_wignerf(l, mid(iqn), m, -btr, wgf)
c check                  if (abs(btr) < 1.0d-13) write(*, '("wigner:", e15.6, i3, i3, i3, e15.6)') btr, l, mid(iqn), m, wgf
                     
                     i_lm = l**2 + l + m + 1
                     i = nsh*nat*nqn*(ich-1) + nsh*nat*(iqn-1) +
     &                   nsh*(icent-1) + i_lm
                     g_nu0 = intgl(i) * wgf * exp(-imnum * m * gmr)
c check                  g_nu0 = intgl(i_lm, icent, iqn, ich) * wgf * exp(-imnum * m * -gmr)
                     g_nu = g_nu + g_nu0
                  enddo
               enddo
               j = 1 + nat*nqn*(ich-1) + nat*(iqn-1) + (icent-1)
               srfc(j) = exp(-kappa * mu(3)) * g_nu
            enddo
            
         enddo
      enddo
         
      end
c     Notes:
c       1) intgl is conceptually a 4D array with size nsh * natom * nqn * nch where
c          nsh=(1+lmax_)^2. That is, the first dimension is for spherical harmoncis,
c          second dimension for atoms, third dimension for parabolic quantum numbers, and
c          the last dimension for ionized channels.
c       2) srfc, the structure factor, is conceptually a 3D array with size natom*nqn*nch.
c          That is, the first dimension is for atoms, second dimension for parabolic 
c          quantum numbers, and last dimension for ionized channels.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_passrot(bt0, gm0, protmat)

      implicit none
      double precision bt0, gm0
      double precision protmat(1:3, 1:3)
      double precision bt, gm

      double precision prot_y(3,3), prot_z(3,3), deg2rad
      parameter (deg2rad=0.01745329251994d0)
      
      
      bt = bt0 * deg2rad
      gm = gm0 * deg2rad
      
      prot_y(1:3, 1) = (/  cos(bt),   0.0d0, -sin(bt) /)
      prot_y(1:3, 2) = (/    0.0d0,   1.0d0,    0.0d0 /)
      prot_y(1:3, 3) = (/  sin(bt),   0.0d0,  cos(bt) /)
      
      prot_z(1:3, 1) = (/  cos(gm), sin(gm),    0.0d0 /)
      prot_z(1:3, 2) = (/ -sin(gm), cos(gm),    0.0d0 /)
      prot_z(1:3, 3) = (/    0.0d0,   0.0d0,    1.0d0 /)  
      
      protmat = matmul(prot_z, prot_y)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_wignerf(j, m1, m2, beta, wgf)

      !!  d_{m1,m2}^j (beta)
      !!  beta is in radian
      implicit none
#include "errquit.fh"      
      integer j, m1, m2
      double precision beta
      double precision wgf
      
      integer k, kmax
      double precision nm1, nm2, nm3, nm4, dm1, dm2, dm3 ,dm4
      double precision cc, ss, gg, factorial
      external factorial
      
      
      if (abs(m1) > j .or. abs(m2) > j) then
         call errquit('WFAT: Arguments m1 and m2 of wfat_wignerf '//
     &        'must be such that abs(m1)>j and abs(m2)>j', 1,
     &        UERR)
      endif
      
      wgf = 0.0d0
      kmax = max(j-m1, j+m2, m2-m1)
      do k = 0, kmax
         if (j-m1 < k .or. j+m2 < k .or. k < m2-m1) then
            cycle
         endif
         
         nm1 = factorial(j + m2)
         nm2 = factorial(j - m2)
         nm3 = factorial(j + m1)
         nm4 = factorial(j - m1)
      
         dm1 = factorial(j + m2 - k)
         dm2 = factorial(k)
         dm3 = factorial(j - m1 - k)
         dm4 = factorial(k + m1 - m2)
      
         cc = (cos(beta/2))**(2*j - 2*k + m2 - m1)
         ss = (sin(beta/2))**(2*k - m2 + m1)
         
         gg = (-1)**(k-m2+m1) * sqrt(nm1 * nm2 * nm3 * nm4) /
     &        (dm1 * dm2 * dm3 * dm4) * cc * ss
      
         wgf = wgf + gg
      enddo

      end
c=======================================================
      
