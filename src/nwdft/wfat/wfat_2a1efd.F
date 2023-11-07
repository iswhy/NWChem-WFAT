c=======================================================
c=======================================================
      subroutine wfat_2a1efd_solve()

      implicit none

      




      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_2a1efd_kin(ni, n, maxn, xq, wq, nb, t)
c     Computes the kinetic energy operator matrix in the FEDVR basis.
      implicit none
#include "errquit.fh"
      integer ni
      integer n(ni)
      integer maxn             ! maxn = maxval(n)
      double precision xq(maxn,ni)
      double precision wq(maxn,ni)
      integer nb               ! nb=sum(n)-ni*2+ni-1
      double precision t(nb,nb)

      logical brg1, brg2, exclude1, exclude2
      integer nb_, j1, i1, m1, j2, i2, m2


      nb_ = (sum(n) - ni*2) + (ni - 1)   ! 1)
      if (nb /= nb_) call errquit('WFAT: Inconsistency found among '//
     &     'the value of the arguments nb, n, and ni of '//
     &     'wfat_2a1efd_kin.', 95, UERR)
      

c     ==== Loop over ket ====
      j2 = 0
      do i2 = 1, ni
         do m2 = 1, n(i2)
c           == Bridge or regular? ==
            brg2 = .false.
            if (m2 == n(i2)) brg2 = .true.

c           == Exclude (1,1) (the 1st basis), (ni,n(ni)) (the last ==
c           == basis), and (:,1) (have been absorbed into bridges) ==
            exclude2 = .false.
            if ((i2==1 .and. m2==1) .or. (i2==ni .and. m2==n(ni))
     &           .or. m2==1) exclude2 = .true.
            if (.not. exclude2) j2 = j2 + 1

c           ==== Loop over bra ====
            j1 = 0
            do i1 = 1, ni
               do m1 = 1, n(i1)
c                 == Bridge or regular? ==
                  brg1 = .false.
                  if (m1 == n(i1)) brg1 = .true.

c                 == Exclude (1,1) (the 1st basis), (ni,n(ni)) (the last ==
c                 == basis), and (:,1) (have been absorbed into bridges) ==
                  exclude1 = .false.
                  if ((i1==1 .and. m1==1) .or. (i1==ni .and. m1==n(ni))
     &                 .or. m1==1) exclude1 = .true.
                  if (.not. exclude1) j1 = j1 + 1

c                 == Calculate the kinetic energy matrix element ==
                  if (.not. (exclude1 .or. exclude2)) then
                     call wfat_2a1efd_kin_elm(brg1, i1, m1, brg2, i2, 
     &                    m2, ni, n, maxn, xq, wq, t(j1,j2))
                  endif

               enddo
            enddo

         enddo
      enddo

      if (j1 /= nb)
     &     call errquit('WFAT: j1 /= nb, this is a bug.', 40, UERR)
      if (j2 /= nb)
     &     call errquit('WFAT: j2 /= nb, this is a bug.', 40, UERR)
      
      end
c     NOTES:
c     1) sum(n)-ni*2 is the number of regular FEDVR bases while ni-1 is the number of bridge
c        FEDVR bases.
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_2a1efd_kin_elm(brg1, i1, m1, brg2, i2, m2, ni, n,
     &     maxn, xq, wq, t)

      implicit none
#include "errquit.fh"
      logical brg1
      integer i1
      integer m1        ! 2 <= m1 <= n(i1)-1. Ignored if brg1=.true..
      logical brg2
      integer i2
      integer m2        ! 2 <= m2 <= n(i2)-1, Ignored if brg2=.true..
      integer ni
      integer n(ni)
      integer maxn             ! maxn = maxval(n)
      double precision xq(maxn,ni)
      double precision wq(maxn,ni)
      double precision t

      integer i, m, n1, n1_, n2, n2_
      double precision df1, df2
      

      n1 = n(i1)
      n1_ = n(i1+1)
      n2 = n(i2)
      n2_ = n(i2+1)

      if ((m1 == 1 .or. m1 == n1) .and. .not. brg1)
     &     call errquit('WFAT: The argument m1 of wfat_2a1efd_kin_elm'//
     &     ' must be such that 2 <= m1 <= n(i1)-1.', 17, UERR)
      if ((m2 == 1 .or. m2 == n2) .and. .not. brg2)
     &     call errquit('WFAT: The argument m2 of wfat_2a1efd_kin_elm'//
     &     ' must be such that 2 <= m2 <= n(i2)-1.', 17, UERR)
      

c     ==== Perform quadrature summation ====
      t = 0.0d0
      do i = 1, ni     ! Loop over FEs.
         do m = 1, n(i)       ! Loop over quadrature points in the i-th FE.
            
c           ==== Evaluate Bra at the (i,m)-th point ====
            if (brg1) then
               if (i == i1 .or. i == i1+1) then
                  call wfat_eval1_fedvr_bridge(i1, n1, xq(:,i1), 
     &                 wq(n1,i1), n1_, xq(:,i1+1), wq(1,i1+1),
     &                 i, m, df1)
               else
                  df1 = 0.0d0
               endif
            else
               if (i == i1) then
                  call wfat_eval1_fedvr(m1, n(i1), xq(:,i1), 
     &                 wq(m1,i1), m, df1)
               else
                  df1 = 0.0d0
               endif
            endif

c           ==== Evaluate Ket at the (i,m)-th point ====
            if (brg2) then
               if (i == i2 .or. i == i2+1) then
                  call wfat_eval1_fedvr_bridge(i2, n2, xq(:,i2), 
     &                 wq(n2,i2), n2_, xq(:,i2+1), wq(1,i2+1),
     &                 i, m, df2)
               else
                  df2 = 0.0d0
               endif
            else
               if (i == i2) then
                  call wfat_eval1_fedvr(m2, n(i2), xq(:,i2), 
     &                 wq(m2,i2), m, df2)
               else
                  df2 = 0.0d0
               endif
            endif
            
c           ==== Quadrature sum ====
            t = t + wq(m,i) * df1 * df2
         enddo
      enddo

      t = t / 2

      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_eval1_fedvr(m, n, x, w, mm, f)
c     DESCRIPTION:
c     Evaluates a regular (non-bridge) FEDVR basis at one of the FEDVR quadrature points
c     WITHIN THE SAME FE. The value of the regular FEDVR basis in the other FEs is zero.
c     The FEDVR basis to be evaluated is characterized by m, n, x, and w, while the
c     evaluation point is specified by mm.
      implicit none
#include "errquit.fh"
      integer m             ! [input]. The DVR index of the FEDVR basis in the i-th FE.
      integer n             ! [input]. The number of DVR points in the i-th FE.
      double precision x(*) ! [input]. The DVR quadrature points in the i-th FE. Needed when m /= mm. When needed, its minimum size should size should be n.
      double precision w    ! [input]. The (i,m)-th quadrature weight.
      integer mm            ! [input]. The DVR index of the evaluation point in the i-th FE.
      double precision f    ! [input]. The value of the (i,m)-th FEDVR basis at the (i,mm)-th quadrature point.

      double precision f0

      if (m < 1 .or. m > n)
     &     call errquit('WFAT: The argument m of '//
     &     'wfat_eval1_lobatto_dm must be in between 1 and n.', 1, UERR)
      if (mm < 1 .or. mm > n)
     &     call errquit('WFAT: The argument mm of '//
     &     'wfat_eval1_lobatto_dm must be in between 1 and n.', 1, UERR)

      if (m /= mm) then
         call wfat_eval1_lobatto_dm(m, n, mm, x, f0)
      else
         call wfat_eval1_lobatto_em(m, n, w, f0)
      endif
      f = f0 / sqrt(w)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval1_fedvr_bridge(i, n1, x1, w1, n2, x2, w2, ii,
     &     mm, f)
c     DESCRIPTION:
c     Evaluates a bridge FEDVR basis at one of the FEDVR quadrature points WITHIN THE SAME
C     OR THE NEXT FE. The value of the bridge FEDVR basis outside of these two FEs is zero.
c     The given value of i represents the bridge function formed by the (i,n1)- and (i+1,1)-th
c     Lobatto shape functions, i.e. the bridge function is the one that bridges the i-th and
c     the (i+1)-th FEs. The bridge FEDVR basis to be evaluated is characterized by i,
c     n1, x1, w1, n2, x2, and w2, while the evaluation point is specified by ii and mm.
      implicit none
#include "errquit.fh"
      integer i
      integer n1             ! [input]. The number of DVR points in the i-th FE.
      double precision x1(*) ! [input]. The DVR quadrature points in the i-th FE. Needed when m /= mm. When needed, its minimum size should be n1.
      double precision w1    ! [input]. The (i,n1)-th quadrature weight.
      integer n2             ! [input]. The number of DVR points in the (i+1)-th FE.
      double precision x2(*) ! [input]. The DVR quadrature points in the (i+1)-th FE. Needed when m /= mm. When needed, its minimum size should be n2.
      double precision w2    ! [input]. The (i+1,1)-th quadrature weight.
      integer ii             ! [input]. The FE index of the evaluation point.
      integer mm             ! [input]. The DVR index of the evaluation point.
      double precision f     ! [input]. The value of the bridge FEDVR basis at the (i,mm)-th quadrature point.

      double precision f1, f2


      if (ii == i) then
         if (mm < 1 .or. mm > n1)
     &        call errquit('WFAT: Since ii=i, the argument mm of  '//
     &        'wfat_eval1_fedvr_bridge must be in between '//
     &        '1 and n1.', 1, UERR)
      else if (ii == i+1) then
         if (mm < 1 .or. mm > n2)
     &        call errquit('WFAT: Since ii=i+1, the argument mm of  '//
     &        'wfat_eval1_fedvr_bridge must be in between '//
     &        '1 and n2.', 1, UERR)
      endif
         
      
c     ==== i side ====
      if (ii == i) then
         if (mm /= n1) then
            call wfat_eval1_lobatto_dm(n1, n1, mm, x1, f1)
         else
            call wfat_eval1_lobatto_em(n1, n1, w1, f1)
         endif
         f2 = 0.0d0
c     ==== i+1 side ====
      else if (ii == i+1) then
         f1 = 0.0d0
         if (mm /= 1) then
            call wfat_eval1_lobatto_dm(1, n2, mm, x2, f2)
         else
            call wfat_eval1_lobatto_em(1, n2, w2, f2)
         endif
      else
         f1 = 0.0d0
         f2 = 0.0d0
      endif
      

      f = (f1 + f2) / sqrt(w1 + w2)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_eval1_lobatto_dm(m, n, mm, x, f)
c     Different m.
      implicit none
#include "errquit.fh"
      integer m
      integer n
      integer mm
      double precision x(n)
      double precision f

      integer k
      double precision a, b

      if (m < 1 .or. m > n)
     &     call errquit('WFAT: The argument m of '//
     &     'wfat_eval1_lobatto_dm must be in between 1 and n.', 1, UERR)
      if (mm < 1 .or. mm > n)
     &     call errquit('WFAT: The argument mm of '//
     &     'wfat_eval1_lobatto_dm must be in between 1 and n.', 1, UERR)
      
      a = 1.0d0 / (x(m) - x(mm))

      b = 1.0d0
      do k = 1, n
         if (k /= m .and. k /= mm) b = b * (x(mm)-x(k)) / (x(m)-x(k))
      enddo
      f = a * b

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_eval1_lobatto_em(m, n, w, f)
c     Equal m.
      implicit none
#include "errquit.fh"
      integer m
      integer n
      double precision w
      double precision f

      integer delkron
      external delkron

      if (m < 1 .or. m > n)
     &     call errquit('WFAT: The argument m of '//
     &     'wfat_eval1_lobatto_em must be in between 1 and n.', 1, UERR)
      
      f = 1.0d0 / (2*w) * (delkron(m,n) - delkron(m,1))
      
      end
c=======================================================
