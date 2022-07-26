cc========================================================
cc========================================================
c
c      subroutine wfat_get_gridnseg(znuc, gqual, nseg)
c
c      implicit none
c      integer gqual
c      integer nseg
c      
c
c      
c      if (gqual == 1) then
c         nseg = 
c      else if (gqual == 2) then
c         nseg = 
c      else if (gqual == 3) then
c         nseg = 
c      else
c         call errquit('WFAT: Bad value of gqual', 52, GEOM_ERR)
c      endif
c      
c      end
cc========================================================


c========================================================
c========================================================
      subroutine wfat_get_gridmaxpt(gq, max_npt)

c     DESCRIPTION:
c       This subroutine returns the biggest possible number of quadrature points
c       around an atom for a given grid quality. As of the writing of this code, 
c       this biggest number of points is associated to atoms having nuclear 
c       charge larger than 19. This is because the numbers of radial and angular 
c       quadrature points increase with nuclear charge for a given gq.
      implicit none
      integer gq         ! [input],  grid quality level
      integer max_npt    ! [output], max no. of points for a given gq
      
#include "wfat.fh"
      integer zn, nrad, nseg0, nlang(nseg)
      double precision atr0(nseg)
      parameter (zn = 54)

      call wfat_get_gridnrad(zn, gq, nrad)
      call wfat_get_gridseg(zn, gq, nseg0, atr0(1:nseg), nlang(1:nseg))

      max_npt = nrad * nlang(nseg0)

      end
c========================================================

      
c========================================================
c========================================================
      subroutine wfat_get_gridpt(zn, gq, max_npt, npt, gpt, gwt)

      implicit none
#include "stdio.fh"
#include "mafdecls.fh"
c     wfat.fh contains: m_emc, nseg      
#include "wfat.fh"
#include "errquit.fh"      
#include "global.fh"      
      integer zn                        ! [input]
      integer gq                        ! [input]
      integer max_npt                   ! [input]
      integer npt                       ! [output]
      double precision gpt(3,max_npt)   ! [output] Quadrature points associated with a nucleus
                                        ! having nuclear charge zn. These quadrature points
                                        ! are in the zn nucleus frame of reference.
      double precision gwt(max_npt)     ! [output] Quadrature weights associated with a nucleus
                                        ! having nuclear charge zn.

      integer ndid, nmax, nseg0, nlang0, nlang(nseg), nrad
      integer l_ax, k_ax, l_ay, k_ay, l_az, k_az, l_aw, k_aw
      integer ipt, isg, iq, ia
      double precision dq, q, q0, qn, r, rwt, atr, atr0(nseg)
      double precision t, p
      double precision pi_const
      
      parameter (q0=0.0d0, qn=1.0d0, pi_const=3.14159265359d0)

      
      ndid = ga_nodeid()
      

c     ==== Get radii and no. of angular grids for all segments  ====
      call wfat_get_gridseg(zn, gq, nseg0, atr0(1:nseg), nlang(1:nseg))
      nmax = maxval(nlang)
      if (ndid==0) write(luout, '("herem0", i6)') nmax

      if (.not. MA_alloc_get(MT_DBL, nmax, 'angx', l_ax, k_ax))
     &     call errquit('WFAT: angx allocation failed', 59, MA_ERR)
      if (ndid==0) write(luout, '("herem1")')
      if (.not. MA_alloc_get(MT_DBL, nmax, 'angy', l_ay, k_ay))
     &     call errquit('WFAT: angy allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nmax, 'angz', l_az, k_az))
     &     call errquit('WFAT: angz allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, nmax, 'awt', l_aw, k_aw))
     &     call errquit('WFAT: ang. weight allocation failed', 59,
     &     MA_ERR)
      
      
c     ==== Get atom radius and no. of radial points ====
      call wfat_get_atomrad(zn, atr)
      call wfat_get_gridnrad(zn, gq, nrad)
      if (ndid==0) write(luout, '("herem2")')

c     ==== Calculate quadrature points and weights ====
      dq = (qn - q0)/(nrad - 1)
      ipt = 0
      isg = 0
      do iq = 1, nrad-2
         q = q0 + dq * iq
         
c        == Get radial quadrature point ====
         call wfat_grid_q2r_em(q, m_emc, atr, r)

c        == Get radial quadrature weight ==
         call wfat_grid_rw_em(q, m_emc, atr, nrad, rwt)
         
c        == Determine segment number ==
         if (isg < nseg0) then
            if (r > atr0(isg+1)*atr) then
               isg = isg + 1
               
c              == Get angular quadrature points and weights ==
               nlang0 = nlang(isg)
               call wfat_grid_lebedev(nlang0,                        !! 1)
     &                                dbl_mb(k_ax:k_ax+nlang0-1),
     &                                dbl_mb(k_ay:k_ay+nlang0-1),
     &                                dbl_mb(k_az:k_az+nlang0-1),
     &                                dbl_mb(k_aw:k_aw+nlang0-1) )
            endif
         endif

ccheck         write(luout, '(i7, ") r = ", e13.5)') iq, r
ccheck         write(luout, '(7x, "  nlang = ", i5)') nlang0

c        == Calculate Cartesian grid points and weights ==
         do ia = 1, nlang0
            ipt = ipt + 1
            
            t = acos(dbl_mb(k_az + ia - 1))
            p = atan(dbl_mb(k_ay+ia-1), dbl_mb(k_ax+ia-1))
            
            gpt(1,ipt) = r * sin(t) * cos(p)
            gpt(2,ipt) = r * sin(t) * sin(p)
            gpt(3,ipt) = r * cos(t)
            
            gwt(ipt) = rwt * 4 * pi_const * dbl_mb(k_aw+ia-1)
         enddo
      enddo
      if (ndid==0) write(luout, '("herem3")')

      npt = ipt
      
      if (isg /= nseg0) then
ccheck         write(luout, '("isg = ", i3)') isg
ccheck         write(luout, '("nseg0 = ", i3)') nseg0
         call errquit('WFAT: Wrong last value of isg - '//
     &                'isg should be equal to nseg0', 59, UERR)
      endif
      
      if (npt > max_npt) then
         call errquit('WFAT: max_npt too small', 59, UERR)
      endif
      if (ndid==0) write(luout, '("herem4")')

ccheck      write(luout, '("Deallocating...")')
      
      if (.not. MA_free_heap(l_ax))
     &     call errquit('WFAT: angx deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_ay))
     &     call errquit('WFAT: angy deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_az))
     &     call errquit('WFAT: angz deallocation failed', 61, MA_ERR)
      if (.not. MA_free_heap(l_aw))
     &     call errquit('WFAT: ang. weight deallocation failed', 61,
     &                  MA_ERR)

ccheck      write(luout, '("End of wfat_get_gridpt")')
      end
c     Notes:
c       1) The angular points are only calculated once for each radial segment.
c========================================================


c========================================================
c========================================================
      subroutine wfat_get_gridnpt(zn, gq, npt)

      implicit none
#include "stdio.fh"
#include "wfat.fh"
#include "errquit.fh"            
      integer zn                        ! [input]
      integer gq                        ! [input]
      integer npt                       ! [output]

      integer nseg0, nlang0, nlang(nseg), nrad
      integer ipt, isg, iq, ia
      double precision dq, q, q0, qn, r, atr, atr0(nseg)

      parameter (q0=0.0d0, qn=1.0d0)


c     ==== Get radii and no. of angular grids for all segments  ====
      call wfat_get_gridseg(zn, gq, nseg0, atr0(1:nseg), nlang(1:nseg))

      
c     ==== Get atom radius and no. of radial points ====
      call wfat_get_atomrad(zn, atr)
      call wfat_get_gridnrad(zn, gq, nrad)

      
c     ==== Calculate the number of quadrature points ====
      dq = (qn - q0)/(nrad - 1)
      ipt = 0
      isg = 0
      do iq = 1, nrad-2
         q = q0 + dq * iq
         
c        == Get radial quadrature point ====
         call wfat_grid_q2r_em(q, m_emc, atr, r)
         
c        == Determine segment number ==
         if (isg < nseg0) then
            if (r > atr0(isg+1)*atr) then
               isg = isg + 1               
               nlang0 = nlang(isg)
            endif
         endif

c        == Calculate Cartesian grid points and weights ==
         ipt = ipt + nlang0
      enddo

      npt = ipt
      
      if (isg /= nseg0) then
ccheck         write(luout, '("isg = ", i3)') isg
ccheck         write(luout, '("nseg0 = ", i3)') nseg0
         call errquit('WFAT: Wrong last value of isg - '//
     &                'isg should be equal to nseg0', 59, UERR)
      endif
      
      end
c========================================================      

      
c========================================================
c========================================================
      subroutine wfat_get_gridseg(zn, gq, nseg0, atr0_, nlang_)

      implicit none
#include "errquit.fh"      
#include "wfat.fh"
      integer zn                  ! [input]
      integer gq                  ! [input]
      integer nseg0               ! [output]
      double precision atr0_(*)   ! [output]
      integer nlang_(*)           ! [output]

      double precision atr0(nseg)
      integer nlang(nseg)

      
c     ==== H - He ====
      if (zn >= 1 .and. zn <=2) then
         nseg0 = 5
         atr0 = (/ 0.00d0,
     &             0.25d0,
     &             0.50d0,
     &             1.00d0,
     &             4.50d0 /)
         if (gq == 1) then
            nlang = (/ 110, 434, 350, 350, 194 /)
         else if (gq == 2) then
            nlang = (/ 194, 590, 434, 434, 266 /)
         else if (gq == 3) then
            nlang = (/ 302, 974, 770, 770, 350 /)
         else if (gq == 4) then
            nlang = (/ 590, 1730, 1202, 1202, 770 /)
         endif
c     ==== Li - Ne ====
      else if (zn >= 3 .and. zn <= 10) then
         nseg0 = 5
         atr0 = (/ 0.0000d0,
     &             0.1667d0,
     &             0.5000d0,
     &             0.9000d0,
     &             3.5000d0 /)
         if (gq == 1) then
            nlang = (/  146, 434, 350,  350, 230 /)
         else if (gq == 2) then
            nlang = (/  302, 770, 590,  590, 302 /)
         else if (gq == 3) then
            nlang = (/  590, 1454, 1202, 1202, 590 /)
         else if (gq == 4) then
            nlang = (/  974, 2354, 1730, 1730, 974 /)
         endif
c     ==== Na - Ar ====
      else if (zn >= 11 .and. zn <= 18) then
         nseg0 = 5
         atr0 = (/ 0.0000d0,
     &             0.1000d0,
     &             0.4000d0,
     &             0.8000d0,
     &             2.5000d0 /)
         if (gq == 1) then
            nlang = (/  302,  590, 434,  434, 350 /)
         else if (gq == 2) then
            nlang = (/  590,  974, 770,  770, 434 /)
         else if (gq == 3) then
            nlang = (/  974, 1730, 1454, 1454, 590 /)
         else if (gq == 4) then
            nlang = (/  1454, 3074, 2354, 2354, 1454 /)
         endif
c     ==== K and beyond ====
      else if (zn >= 19) then
         nseg0 = 1
         atr0 = (/ 0.0000d0,
     &             0.0000d0,
     &             0.0000d0,
     &             0.0000d0,
     &             0.0000d0 /)
         if (gq == 1) then
            nlang = (/  590, 0, 0, 0, 0 /)
         else if (gq == 2) then
            nlang = (/  974, 0, 0, 0, 0 /)
         else if (gq == 3) then
            nlang = (/ 1454, 0, 0, 0, 0 /)
         else if (gq == 4) then
            nlang = (/ 3074, 0, 0, 0, 0 /)
         endif
c     ==== Anything else ====
      else
         call errquit('WFAT: Error in determining atr0 and nlang',
     &        52, UERR)
      endif

      atr0_(1:nseg) = atr0
      nlang_(1:nseg) = nlang
      
      end
c========================================================


cc OLD c========================================================
cc OLD c========================================================
cc OLD       subroutine wfat_get_gridseg(zn, gq, nseg0, atr0_, nlang_)
cc OLD 
cc OLD       implicit none
cc OLD #include "errquit.fh"      
cc OLD #include "wfat.fh"
cc OLD       integer zn                  ! [input]
cc OLD       integer gq                  ! [input]
cc OLD       integer nseg0               ! [output]
cc OLD       double precision atr0_(*)   ! [output]
cc OLD       integer nlang_(*)           ! [output]
cc OLD 
cc OLD       double precision atr0(nseg)
cc OLD       integer nlang(nseg)
cc OLD 
cc OLD       
cc OLD c     ==== H - He ====
cc OLD       if (zn >= 1 .and. zn <=2) then
cc OLD          nseg0 = 5
cc OLD          atr0 = (/ 0.00d0,
cc OLD      &             0.25d0,
cc OLD      &             0.50d0,
cc OLD      &             1.00d0,
cc OLD      &             4.50d0 /)
cc OLD          if (gq == 1) then
cc OLD             nlang = (/ 38,  74, 170, 350, 194 /)
cc OLD          else if (gq == 2) then
cc OLD             nlang = (/ 50,  86, 230, 434, 266 /)
cc OLD          else if (gq == 3) then
cc OLD             nlang = (/ 86, 146, 434, 770, 350 /)
cc OLD          endif
cc OLD c     ==== Li - Ne ====
cc OLD       else if (zn >= 3 .and. zn <= 10) then
cc OLD          nseg0 = 5
cc OLD          atr0 = (/ 0.0000d0,
cc OLD      &             0.1667d0,
cc OLD      &             0.5000d0,
cc OLD      &             0.9000d0,
cc OLD      &             3.5000d0 /)
cc OLD          if (gq == 1) then
cc OLD             nlang = (/  50,  86, 266,  350, 230 /)
cc OLD          else if (gq == 2) then
cc OLD             nlang = (/  74, 110, 302,  590, 302 /)
cc OLD          else if (gq == 3) then
cc OLD             nlang = (/ 110, 230, 590, 1202, 590 /)
cc OLD          endif
cc OLD c     ==== Na - Ar ====
cc OLD       else if (zn >= 11 .and. zn <= 18) then
cc OLD          nseg0 = 5
cc OLD          atr0 = (/ 0.0000d0,
cc OLD      &             0.1000d0,
cc OLD      &             0.4000d0,
cc OLD      &             0.8000d0,
cc OLD      &             2.5000d0 /)
cc OLD          if (gq == 1) then
cc OLD             nlang = (/   74, 170, 350,  434, 350 /)
cc OLD          else if (gq == 2) then
cc OLD             nlang = (/  110, 230, 590,  770, 434 /)
cc OLD          else if (gq == 3) then
cc OLD             nlang = (/  230, 350, 770, 1454, 590 /)
cc OLD          endif
cc OLD c     ==== K and beyond ====
cc OLD       else if (zn >= 19) then
cc OLD          nseg0 = 1
cc OLD          atr0 = (/ 0.0000d0,
cc OLD      &             0.0000d0,
cc OLD      &             0.0000d0,
cc OLD      &             0.0000d0,
cc OLD      &             0.0000d0 /)
cc OLD          if (gq == 1) then
cc OLD             nlang = (/  590, 0, 0, 0, 0 /)
cc OLD          else if (gq == 2) then
cc OLD             nlang = (/  974, 0, 0, 0, 0 /)
cc OLD          else if (gq == 3) then
cc OLD             nlang = (/ 1454, 0, 0, 0, 0 /)
cc OLD          endif
cc OLD c     ==== Anything else ====
cc OLD       else
cc OLD          call errquit('WFAT: Error in determining atr0 and nlang',
cc OLD      &        52, UERR)
cc OLD       endif
cc OLD 
cc OLD       atr0_(1:nseg) = atr0
cc OLD       nlang_(1:nseg) = nlang
cc OLD       
cc OLD       end
cc OLD c========================================================      

      
c========================================================
c========================================================
      subroutine wfat_get_atomrad(zn, atr)

c     Atom radii defined in Table 1 of Ref. [1]. These values might have
c     been arbitrarily interpolated.
      implicit none
#include "errquit.fh"      
      integer zn
      double precision atr

      double precision ang2bohr
      parameter (ang2bohr = 1.8897261255d0)
      

c     ==== 1st row ====
      if (zn ==  1) atr = 0.35d0
      if (zn ==  2) atr = 0.30d0 
                               
c     ==== 2nd row ====        
      if (zn ==  3) atr = 1.45d0
      if (zn ==  4) atr = 1.05d0
      if (zn ==  5) atr = 0.80d0 
      if (zn ==  6) atr = 0.70d0
      if (zn ==  7) atr = 0.65d0
      if (zn ==  8) atr = 0.60d0
      if (zn ==  9) atr = 0.50d0
      if (zn == 10) atr = 0.40d0 
                               
c     ==== 3rd row ====        
      if (zn == 11) atr = 1.80d0
      if (zn == 12) atr = 1.50d0
      if (zn == 13) atr = 1.25d0
      if (zn == 14) atr = 1.10d0
      if (zn == 15) atr = 1.10d0 
      if (zn == 16) atr = 1.00d0
      if (zn == 17) atr = 1.00d0
      if (zn == 18) atr = 0.90d0 
                               
c     ==== 4th row ====        
      if (zn == 19) atr = 2.20d0
      if (zn == 20) atr = 1.80d0
      if (zn == 21) atr = 1.60d0 
      if (zn == 22) atr = 1.40d0
      if (zn == 23) atr = 1.40d0 
      if (zn == 24) atr = 1.40d0  
      if (zn == 25) atr = 1.40d0
      if (zn == 26) atr = 1.40d0
      if (zn == 27) atr = 1.35d0
      if (zn == 28) atr = 1.35d0
      if (zn == 29) atr = 1.35d0
      if (zn == 30) atr = 1.35d0
      if (zn == 31) atr = 1.25d0 
      if (zn == 32) atr = 1.25d0 
      if (zn == 33) atr = 1.15d0
      if (zn == 34) atr = 1.15d0
      if (zn == 35) atr = 1.15d0
      if (zn == 36) atr = 1.08d0 

c     ==== 5th row ====         
      if (zn == 37) atr = 2.35d0
      if (zn == 38) atr = 2.00d0
      if (zn == 39) atr = 2.00d0 
      if (zn == 40) atr = 1.95d0 
      if (zn == 41) atr = 1.90d0 
      if (zn == 42) atr = 1.85d0   
      if (zn == 43) atr = 1.80d0 
      if (zn == 44) atr = 1.75d0 
      if (zn == 45) atr = 1.70d0 
      if (zn == 46) atr = 1.65d0 
      if (zn == 47) atr = 1.60d0
      if (zn == 48) atr = 1.55d0
      if (zn == 49) atr = 1.50d0 
      if (zn == 50) atr = 1.45d0
      if (zn == 51) atr = 1.45d0
      if (zn == 52) atr = 1.40d0
      if (zn == 53) atr = 1.40d0
      if (zn == 54) atr = 1.40d0 

      if (zn < 1 .or. zn > 54) then
         call errquit('WFAT: The choice of zn in wfat_get_atomrad '//
     &        'is not available.', 1, UERR)
      endif
      
      atr = atr * ang2bohr
      
      end
c References:
c   [1] Chem. Phys. Lett. 1993, 209, 506.      
c========================================================


c========================================================
c========================================================
      subroutine wfat_get_gridnrad(zn, gq, nrad)

      implicit none
#include "errquit.fh"
      integer zn         ! [input]
      integer gq         ! [input]
      integer nrad       ! [output]


      if (zn >= 1 .and. zn <= 2) then
c     ==== 1st row ====
         if (gq == 1) nrad =  35
         if (gq == 2) nrad =  55
         if (gq == 3) nrad =  90
         if (gq == 4) nrad =  180
      else if (zn >= 3 .and. zn <= 10) then
c     ==== 2nd row ====
         if (gq == 1) nrad =  49
         if (gq == 2) nrad =  70
         if (gq == 3) nrad = 100
         if (gq == 4) nrad = 220
      else if (zn >= 11 .and. zn <= 18) then
c     ==== 3rd row ====
         if (gq == 1) nrad =  88
         if (gq == 2) nrad = 123
         if (gq == 3) nrad = 135
         if (gq == 4) nrad = 290
      else if (zn >= 19 .and. zn <= 36) then
c     ==== 4th row ====
         if (gq == 1) nrad = 112
         if (gq == 2) nrad = 130
         if (gq == 3) nrad = 160
         if (gq == 4) nrad = 350
      else if (zn >= 37 .and. zn <= 54) then
c     ==== 5th row ====
         if (gq == 1) nrad = 123
         if (gq == 2) nrad = 141
         if (gq == 3) nrad = 205
         if (gq == 4) nrad = 440
      else
         call errquit('WFAT: The choice of zn in wfat_get_gridnrad '//
     &        'is not available.', 60, UERR)
      endif

      end
c     References:
c       [1] NWChem User Guide pg. 132-133.
c========================================================


c========================================================
c========================================================
      subroutine wfat_grid_q2r_em(q, m, alpha, r)

      implicit none
      integer m               ! [input]
      double precision alpha  ! [input]
      double precision q      ! [input]
      double precision r      ! [output]
    
      r = alpha * q**m * (1-q)**(-m)
      
      end
c========================================================


c=================================================
c=================================================
      subroutine wfat_grid_rw_em(q, m, alpha, nr, w)

      implicit none
      integer m                ! [input]
      integer nr               ! [input]
      double precision alpha   ! [input]
      double precision q       ! [input]
      double precision w       ! [output]
      
      w = m * alpha**3 * q**(3*m-1) / (nr-1) / (1-q)**(3*m+1)
      
      end
c=================================================


c=================================================
c=================================================
      subroutine wfat_get_gridrange(geomhan, id, gq, ic, ipt,
     &     natom_proc, npt_proc)
c     DESCRIPTION :
c     This subroutine returns the range of centers and the quadrature
c     points associated with each center for the current process. The
c     way the quadrature points are distributed into the process are
c     illustrated in the following example. Assume that the molecule
c     consists of four atoms, hence four centers, and where the 
c     quadrature points around these centers are given in the table
c     below.
cc    ----------------------------------
cc    center1  center2  center3  center4
cc    ----------------------------------
cc        pt1      pt1      pt1      pt1
cc        pt2      pt2      pt2      pt2
cc        pt3      pt3      pt3      pt3
cc        pt4      pt4      pt4         
cc        pt5      pt5      pt5         
cc        pt6      pt6      pt6         
cc        pt7      pt7      pt7         
cc                 pt8      pt8         
cc                 pt9      pt9
cc                         pt10
cc                         pt11
c     Hence, in total there are 30 quadrature points (npt=30). 
c     Further, assume that the job is run on four nodes (nnodes=4). 
c     Then, the number of points, npt_proc, for each node is calculated
c     as,
c        - npt_proc = floor(30/4) + 1    for the first nrem processes
c        - npt_proc = floor(30/4)        for the rest of the processes
c     where nrem is the remainder between npt and nnodes, which for
c     this example gives nrem=2. Therefore, the distribution of ic,
c     ipt, natom_proc, and npt_proc among the processes are
c     
cc    process0 :
cc       natom_proc = 2
cc       npt_proc = 8
cc       ic =  <   1   2 >
cc       ipt = / pt1 pt1 \
cc             \ pt7 pt1 /
cc    process1 :
cc       natom_proc = 1
cc       npt_proc = 8 
cc       ic =  <   2 >
cc       ipt = / pt2 \
cc             \ pt9 /
cc    process2 :
cc       natom_proc = 1
cc       npt_proc = 7
cc       ic =  <   3 >
cc       ipt = / pt1 \
cc             \ pt7 /
cc    process3 :
cc       natom_proc = 2
cc       npt_proc = 7
cc       ic =  <    3   4 >
cc       ipt = /  pt8 pt1 \
cc             \ pt11 pt3 /
      implicit none
#include "geom.fh"
#include "global.fh"
#include "errquit.fh"      
      integer geomhan         ! [input], geometry handle.
      integer id              ! [input], node ID of the current process.
      integer gq              ! [input], the chosen grid quality.
      integer ic(1:*)         ! [output], the first natom_proc elements contains the centers for the current process.
      integer ipt(2,1:*)      ! [output], the first natom_proc columns contains the starting (1st row) and ending (2nd row) points of the corresponding center for the current process.
      integer natom_proc      ! [output], no. of centers for the current process.
      integer npt_proc        ! [output], total no. of points for the current process.

      character*16 tag
      integer nnodes, npt, npt0, natom, nrem, npt_proc0
      integer ic0, ipt0, inode, ia_proc, ipt_proc
      double precision atp(3), znuc


      nnodes = ga_nnodes()

      if (.not. geom_ncent(geomhan, natom))
     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)

      
c     ==== Calculate total number of grid points ====
      npt = 0
      do ic0 = 1, natom
         if (.not. geom_cent_get(geomhan, ic0, tag, atp, znuc))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
         call wfat_get_gridnpt(int(znuc), gq, npt0)
         npt = npt + npt0
      enddo
      nrem = mod(npt, nnodes)
      npt_proc0 = int(npt/nnodes)

      
      inode = 0
      ipt_proc = 0
      ia_proc = 0
c     ==== Loop over atoms ====
      do ic0 = 1, natom
         if (.not. geom_cent_get(geomhan, ic0, tag, atp, znuc))
     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
         call wfat_get_gridnpt(int(znuc), gq, npt0)

c        ==== Loop over grid points in an atom ====
         do ipt0 = 1, npt0
            ipt_proc = ipt_proc + 1

c           ==== Determine the number of grid points in node id ====
            if (inode <= nrem-1) then
               npt_proc = npt_proc0 + 1
            else
               npt_proc = npt_proc0
            endif

            if (inode == id) then
c              ==== Determine initial grid indices ====               
               if (ipt0 == 1 .or. ipt_proc == 1) then
                  ia_proc = ia_proc + 1
                  ic(ia_proc) = ic0
                  ipt(1,ia_proc) = ipt0
               endif
               
c              ==== Determine final grid indices ====
               if (ipt0 == npt0 .or. ipt_proc == npt_proc) then
                  ipt(2,ia_proc) = ipt0
               endif

c              ==== Exit point ====
               if (ipt_proc == npt_proc) then
                  natom_proc = ia_proc
                  return
               endif
            endif

c           ==== Increment node id and reset grid id ====
            if (ipt_proc == npt_proc) then
               inode = inode + 1
               ipt_proc = 0
            endif

         enddo
      enddo
      

      end
c=================================================


cOLDc=================================================
cOLDc=================================================
cOLD      subroutine wfat_get_gridrange(geomhan, id, gq, ic1, ipt1,
cOLD     &                              ic2, ipt2, npt_proc)
cOLD
cOLD      implicit none
cOLD#include "geom.fh"
cOLD#include "global.fh"
cOLD#include "errquit.fh"      
cOLD      integer geomhan
cOLD      integer id
cOLD      integer gq
cOLD      integer ic1
cOLD      integer ipt1
cOLD      integer ic2
cOLD      integer ipt2
cOLD      integer npt_proc
cOLD
cOLD      character*16 tag
cOLD      integer nnodes, npt, npt0, natom, nrem, npt_proc0
cOLD      integer i, inode, ipt_proc, ipt
cOLD      double precision atp(3), znuc
cOLD
cOLD
cOLD      nnodes = ga_nnodes()
cOLD
cOLD      if (.not. geom_ncent(geomhan, natom))
cOLD     &     call errquit('WFAT: Getting natom failed', 51, GEOM_ERR)
cOLD
cOLD      
cOLDc     ==== Calculate total number of grid points ====
cOLD      npt = 0
cOLD      do i = 1, natom
cOLD         if (.not. geom_cent_get(geomhan, i, tag, atp, znuc))
cOLD     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
cOLD         call wfat_get_gridnpt(int(znuc), gq, npt0)
cOLD         npt = npt + npt0
cOLD      enddo
cOLD      nrem = mod(npt, nnodes)
cOLD      npt_proc0 = int(npt/nnodes)
cOLD
cOLD      
cOLD      inode = 0
cOLD      ipt_proc = 0
cOLDc     ==== Loop over atoms ====
cOLD      do i = 1, natom
cOLD         if (.not. geom_cent_get(geomhan, i, tag, atp, znuc))
cOLD     &        call errquit('WFAT: Getting znuc failed', 52, GEOM_ERR)
cOLD         call wfat_get_gridnpt(int(znuc), gq, npt0)
cOLD
cOLDc        ==== Loop over grid points in an atom ====
cOLD         do ipt = 1, npt0
cOLD            ipt_proc = ipt_proc + 1
cOLD
cOLDc           ==== Determine the number of grid points in node id ====
cOLD            if (inode <= nrem-1) then
cOLD               npt_proc = npt_proc0 + 1
cOLD            else
cOLD               npt_proc = npt_proc0
cOLD            endif
cOLD
cOLDc           ==== Determine initial grid indices ====
cOLD            if (inode == id .and. ipt_proc == 1) then
cOLD               ic1 = i
cOLD               ipt1 = ipt
cOLD            endif
cOLD
cOLDc           ==== Determine final grid indices ====
cOLD            if (inode == id .and. ipt_proc == npt_proc) then
cOLD               ic2 = i
cOLD               ipt2 = ipt
cOLD               return
cOLD            endif
cOLD
cOLDc           ==== Increment node id and reset grid id ====
cOLD            if (ipt_proc == npt_proc) then
cOLD               inode = inode + 1
cOLD               ipt_proc = 0
cOLD            endif
cOLD
cOLD         enddo
cOLD      enddo
cOLD      
cOLD
cOLD      end
cOLDc=================================================            
