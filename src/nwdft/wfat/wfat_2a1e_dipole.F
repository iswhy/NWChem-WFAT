!========================================================
!========================================================
      subroutine wfat_dipole_dvrr(m1, n1, m2, n2, ke, dr)
cc     Calculates D, the matrix of the radial part of the dipole matrix in DVR
cc     basis. The elements of D are calculated as
cc     $$
cc       D_{ij} = < R_i^{m1} | r | R_j^{m2} >
cc     $$
cc     where the integration is
cc     $$
cc       \int_0^\infty dr r^2 ...
cc     $$
cc     and $R_i^{m1}$ is the radial DVR function.
      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer m1
      integer n1
      integer m2
      integer n2
      double precision ke
      double precision dr(n1,n2)

      integer l_u1, k_u1, l_u2, k_u2, l_dr0, k_dr0, l_tm, k_tm


      if (m1 < 0)
     &     call errquit('WFAT: The argument m1 of wfat_dipole_dvrr '//
     &     'must be zero or a positive integer.', 85, UERR)
      if (m2 < 0)
     &     call errquit('WFAT: The argument m2 of wfat_dipole_dvrr '//
     &     'must be zero or a positive integer.', 85, UERR)
      if (ke <= 0.0d0)
     &     call errquit('WFAT: The argument ke of wfat_dipole_dvrr '//
     &     'must be a positive number.', 62, UERR)
      

c     ==== Transformation matrix for m1 ====
      if (.not. ma_alloc_get(MT_DBL, n1*n1, 'u1', l_u1, k_u1))
     &     call errquit('WFAT: Cannot allocate u1', 1, MA_ERR)
      call wfat_dvr_fbr_r(m1, n1, dbl_mb(k_u1))

c     ==== Transformation matrix for m2 ====
      if (.not. ma_alloc_get(MT_DBL, n2*n2, 'u2', l_u2, k_u2))
     &     call errquit('WFAT: Cannot allocate u2', 1, MA_ERR)
      call wfat_dvr_fbr_r(m2, n2, dbl_mb(k_u2))

c     ==== The radial part of dipole matrix in FBR basis ====
      if (.not. ma_alloc_get(MT_DBL, n1*n2, 'dr_fbr', l_dr0, k_dr0))
     &     call errquit('WFAT: Cannot allocate dr_fbr', 1, MA_ERR)
      call wfat_dipole_fbrr(m1, n1, m2, n2, ke, dbl_mb(k_dr0))

c     ==== The radial part of dipole matrix in DVR basis ====
      if (.not. ma_alloc_get(MT_DBL, n1*n2, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: Cannot allocate tmp', 1, MA_ERR)
      call wfat_dgemm_wrap('n', 'n', n1, n2, n2, n2, 1.0d0, 
     &     dbl_mb(k_dr0), dbl_mb(k_u2), 0.0d0, dbl_mb(k_tm))
      call wfat_dgemm_wrap('t', 'n', n1, n1, n1, n2, 1.0d0, 
     &     dbl_mb(k_u1), dbl_mb(k_tm), 0.0d0, dr)

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_u1))
     &     call errquit('WFAT: Cannot deallocate u1', 1, MA_ERR)
      if (.not. ma_free_heap(l_u2))
     &     call errquit('WFAT: Cannot deallocate u2', 1, MA_ERR)
      if (.not. ma_free_heap(l_dr0))
     &     call errquit('WFAT: Cannot deallocate dr_fbr', 1, MA_ERR)
      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp', 1, MA_ERR)
      
      end
!========================================================

      
!========================================================
!========================================================
      subroutine wfat_dipole_fbrr(m1, n1, m2, n2, ke, dr)
cc     Calculates D, the matrix of the radial part of the dipole matrix in FBR
cc     basis. The elements of D are calculated as
cc     $$
cc       D_{ij} = < R_i^{m1} | r | R_j^{m2} >
cc     $$
cc     where the integration is
cc     $$
cc       \int_0^\infty dr r^2 ...
cc     $$
cc     and $R_i^{m1}$ is the radial FBR function.
      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer m1
      integer n1
      integer m2
      integer n2
      double precision ke
      double precision dr(n1,n2)

      integer i, j
      double precision a, b, c, normi, normj, laguerre_norm_i
      external laguerre_norm_i


      if (m1 < 0)
     &     call errquit('WFAT: The argument m1 of wfat_dipole_fbrr '//
     &     'must be zero or a positive integer.', 85, UERR)
      if (m2 < 0)
     &     call errquit('WFAT: The argument m2 of wfat_dipole_fbrr '//
     &     'must be zero or a positive integer.', 85, UERR)
      if (ke <= 0.0d0)
     &     call errquit('WFAT: The argument ke of wfat_dipole_fbrr '//
     &     'must be a positive number.', 62, UERR)
      

      a = dble(m1+2)
      b = dble(m2+2)
      c = (m1+m2)/2.0d0 + 3

      do j = 1, n2
         normj = laguerre_norm_i(j-1, m2+2)
         do i = 1, n1
            normi = laguerre_norm_i(i-1, m1+2)
            
            call wfat_laguerre_intg1(a, b, c, i-1, j-1, dr(i,j))
            dr(i,j) = dr(i,j) / sqrt(normi*normj) / ke
         enddo
      enddo

      end
!========================================================


cLOW_ACC !========================================================
cLOW_ACC !========================================================
cLOW_ACC       subroutine wfat_dipole_dvrr_num(m1, n1, m2, n2, ke, dr)
cLOW_ACC 
cLOW_ACC       implicit none
cLOW_ACC #include "mafdecls.fh"      
cLOW_ACC #include "errquit.fh"
cLOW_ACC       integer m1
cLOW_ACC       integer n1
cLOW_ACC       integer m2
cLOW_ACC       integer n2
cLOW_ACC       double precision ke
cLOW_ACC       double precision dr(n1,n2)
cLOW_ACC 
cLOW_ACC       integer l_x, l_rw, l_dv
cLOW_ACC       integer k_x, k_rw, k_dv
cLOW_ACC       integer j, i
cLOW_ACC       double precision x, w, dvi
cLOW_ACC 
cLOW_ACC 
cLOW_ACC       if (.not. ma_alloc_get(MT_DBL, n2, 'xquad', l_x, k_x))
cLOW_ACC      &     call errquit('WFAT: Cannot allocate xquad', 1, MA_ERR)
cLOW_ACC       if (.not. ma_alloc_get(MT_DBL, n2, 'xweight', l_rw, k_rw))
cLOW_ACC      &     call errquit('WFAT: Cannot allocate xweight', 1, MA_ERR)
cLOW_ACC       call wfat_lf_quadrature_rule(n2, dble(abs(m2))+2, dbl_mb(k_x),
cLOW_ACC      &     dbl_mb(k_rw))
cLOW_ACC 
cLOW_ACC       if (.not. ma_alloc_get(MT_DBL, n1, 'dvr_i', l_dv, k_dv))
cLOW_ACC      &     call errquit('WFAT: Cannot allocate dvr_i', 1, MA_ERR)
cLOW_ACC 
cLOW_ACC       
cLOW_ACC       do j = 1, n2
cLOW_ACC          x = dbl_mb(k_x+j-1)
cLOW_ACC          w = dbl_mb(k_rw+j-1)
cLOW_ACC          call wfat_eval_dvrr(m1, n1, ke, x/ke, dbl_mb(k_dv))
cLOW_ACC          do i = 1, n1
cLOW_ACC             dvi = dbl_mb(k_dv+i-1)
cLOW_ACC             dr(i,j) = sqrt(w)/ke/sqrt(ke) * x**(-abs(m2)/2.0d0) * 
cLOW_ACC      &           exp(x/2) * dvi * x/ke
cLOW_ACC          enddo
cLOW_ACC       enddo
cLOW_ACC 
cLOW_ACC 
cLOW_ACC       if (.not. ma_free_heap(l_x))
cLOW_ACC      &     call errquit('WFAT: Cannot deallocate xquad', 1, MA_ERR)
cLOW_ACC       if (.not. ma_free_heap(l_rw))
cLOW_ACC      &     call errquit('WFAT: Cannot deallocate xweight', 1, MA_ERR)
cLOW_ACC       if (.not. ma_free_heap(l_dv))
cLOW_ACC      &     call errquit('WFAT: Cannot deallocate dvr_i', 1, MA_ERR)
cLOW_ACC       
cLOW_ACC       end
cLOW_ACC !========================================================
      

cfordebug!========================================================
cfordebug!========================================================
cfordebug      subroutine wfat_dipole_dvrr_num2(m1, n1, m2, n2, ke, dr)
cfordebug
cfordebug      implicit none
cfordebug#include "stdio.fh"
cfordebug#include "mafdecls.fh"      
cfordebug#include "errquit.fh"
cfordebug      integer m1
cfordebug      integer n1
cfordebug      integer m2
cfordebug      integer n2
cfordebug      double precision ke
cfordebug      double precision dr(n1,n2)
cfordebug
cfordebug      integer i, j, k, l_tm1, k_tm1, l_tm2, k_tm2, i_tm0, i_tmn, i_tm
cfordebug      integer j_tm, ntpz
cfordebug      double precision x0, xn, dx, x, s, s1, sn
cfordebug
cfordebug
cfordebug      ntpz = 10000
cfordebug      if (.not. ma_alloc_get(MT_DBL, ntpz*n1, 'tmp1', l_tm1, k_tm1))
cfordebug     &     call errquit('WFAT: Cannot allocate tmp1', 1, MA_ERR)
cfordebug      if (.not. ma_alloc_get(MT_DBL, ntpz*n2, 'tmp2', l_tm2, k_tm2))
cfordebug     &     call errquit('WFAT: Cannot allocate tmp2', 1, MA_ERR)
cfordebug
cfordebug
cfordebugc     ==== Check the overlap of the DVR bases ====
cfordebugc     == Construct the trapezoid points ==
cfordebug      x0 = 0.0d0
cfordebug      xn = 100.0d0
cfordebug      dx = (xn-x0) / (ntpz-1)
cfordebug      do i = 1, ntpz
cfordebug         x = x0 + (i-1)*dx
cfordebug         
cfordebug         i_tm0 = k_tm1 + (1-1)*ntpz + (i-1)      ! The first DVR at x.
cfordebug         i_tmn = k_tm1 + (n1-1)*ntpz + (i-1)     ! The last DVR at x.
cfordebug         call wfat_eval_dvrr(m1, n1, ke, x/ke, dbl_mb(i_tm0:i_tmn:ntpz))   ! 1)
cfordebug
cfordebug         i_tm0 = k_tm2 + (1-1)*ntpz + (i-1)      ! The first DVR at x.
cfordebug         i_tmn = k_tm2 + (n2-1)*ntpz + (i-1)     ! The last DVR at x.
cfordebug         call wfat_eval_dvrr(m2, n2, ke, x/ke, dbl_mb(i_tm0:i_tmn:ntpz))   ! 1)
cfordebug      enddo
cfordebugc     Notes:
cfordebugc     1) The row of tmp are the trapezoid points, the columns are the DVR functions.
cfordebug
cfordebug
cfordebug      do j = 1, n2
cfordebug         do i = 1, n1
cfordebugc           ==== Trapezoid rule ====
cfordebug            s = 0.0d0
cfordebug            do k = 2, ntpz-1
cfordebug               x = x0 + (k-1)*dx
cfordebug               
cfordebug               i_tm = k_tm1 + (i-1)*ntpz + (k-1)
cfordebug               j_tm = k_tm2 + (j-1)*ntpz + (k-1)
cfordebug               s = s + dbl_mb(i_tm) * (x/ke)**3 * dbl_mb(j_tm)
cfordebug            enddo
cfordebug            s1 = dbl_mb(k_tm1+(i-1)*ntpz) * dbl_mb(k_tm2+(j-1)*ntpz) *
cfordebug     &           (x0/ke)**3
cfordebug            sn = dbl_mb(k_tm1+i*ntpz-1) * dbl_mb(k_tm2+j*ntpz-1) *
cfordebug     &           (xn/ke)**3
cfordebug            
cfordebug            dr(i,j) = dx/ke * (0.5d0*s1 + s + 0.5d0*sn)
cfordebug         enddo
cfordebug      enddo
cfordebug
cfordebug
cfordebug      if (.not. ma_free_heap(l_tm1))
cfordebug     &     call errquit('WFAT: Cannot deallocate tmp1', 1, MA_ERR)
cfordebug      if (.not. ma_free_heap(l_tm2))
cfordebug     &     call errquit('WFAT: Cannot deallocate tmp2', 1, MA_ERR)
cfordebug      end
cfordebug!========================================================
      
      
!========================================================
!========================================================
      subroutine wfat_dipole_dvrr_num_OLD(m, n, ke, dr)

      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer n
      double precision ke
      double precision dr(n,n)

      integer i, l_x, k_x, l_rw, k_rw

      if (.not. ma_alloc_get(MT_DBL, n, 'xquad', l_x, k_x))
     &     call errquit('WFAT: Cannot allocate xquad', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, n, 'xweight', l_rw, k_rw))
     &     call errquit('WFAT: Cannot allocate xweight', 1, MA_ERR)
      call wfat_lf_quadrature_rule(n, abs(dble(m))+2, dbl_mb(k_x),
     &     dbl_mb(k_rw))

      dr = 0.0d0
      do i = 1, n
         dr(i,i) = dbl_mb(k_x+i-1) / ke
      enddo

      if (.not. ma_free_heap(l_x))
     &     call errquit('WFAT: Cannot deallocate xquad', 1, MA_ERR)
      if (.not. ma_free_heap(l_rw))
     &     call errquit('WFAT: Cannot deallocate xweight', 1, MA_ERR)
      
      end
!========================================================


!========================================================
!========================================================
      subroutine wfat_dipole_dvrt(m, p, n1, n2, dt)
cc     Calculates D, the matrix of the polar part of the dipole matrix in DVR
cc     basis. The elements of D depend on the value of p:
cc     1) if p = -1, then
cc          D_{ij} = < T_i^{m-1} | \sin\theta | T_j^m >,
cc     2) if p = 0, then
cc          D_{ij} = < T_i^m | \cos\theta | T_j^m >, or
cc     3) if p = 1, then
cc          D_{ij} = < T_i^{m+1} | \sin\theta | T_j^m >,
cc     where the size of D is n1 by n2, and T_i^m is the polar DVR function.
      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer m
      integer p      ! -1, 0, +1
      integer n1     ! m-1, m, m+1
      integer n2     ! m
      double precision dt(n1,n2)

      integer l_u1, k_u1, l_u2, k_u2, l_dt0, k_dt0, l_tm, k_tm


      if (m < 0)
     &     call errquit('WFAT: The argument m of wfat_dipole_dvrt '//
     &     'must be zero or a positive integer.', 85, UERR)
      if (m == 0 .and. p == -1)
     &     call errquit('WFAT: If the argument m of wfat_dipole_dvrt '//
     &     'has the value of zero, then the allowed values of p is 0 '//
     &     'or 1', 62, UERR)
      if (p /= -1 .and. p /= 0 .and. p /= 1)
     &     call errquit('WFAT: The argument p of wfat_dipole_dvrt '//
     &     'must be equal to either -1, 0, or 1.', 59, UERR)

c     ==== Transformation matrix for m+p ====
      if (.not. ma_alloc_get(MT_DBL, n1*n1, 'u1', l_u1, k_u1))
     &     call errquit('WFAT: Cannot allocate u1', 1, MA_ERR)
      call wfat_dvr_fbr_t(m+p, n1, dbl_mb(k_u1))

c     ==== Transformation matrix for m ====
      if (.not. ma_alloc_get(MT_DBL, n2*n2, 'u2', l_u2, k_u2))
     &     call errquit('WFAT: Cannot allocate u2', 1, MA_ERR)
      call wfat_dvr_fbr_t(m, n2, dbl_mb(k_u2))

c     ==== The polar part of dipole matrix in FBR basis ====
      if (.not. ma_alloc_get(MT_DBL, n1*n2, 'dt_fbr', l_dt0, k_dt0))
     &     call errquit('WFAT: Cannot allocate dt_fbr', 1, MA_ERR)
      call wfat_dipole_fbrt(m, p, n1, n2, dbl_mb(k_dt0))

c     ==== The polar part of dipole matrix in DVR basis ====
      if (.not. ma_alloc_get(MT_DBL, n1*n2, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: Cannot allocate tmp', 1, MA_ERR)
      call wfat_dgemm_wrap('n', 'n', n1, n2, n2, n2, 1.0d0, 
     &     dbl_mb(k_dt0), dbl_mb(k_u2), 0.0d0, dbl_mb(k_tm))
      call wfat_dgemm_wrap('t', 'n', n1, n1, n1, n2, 1.0d0, 
     &     dbl_mb(k_u1), dbl_mb(k_tm), 0.0d0, dt)

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_u1))
     &     call errquit('WFAT: Cannot deallocate u1', 1, MA_ERR)
      if (.not. ma_free_heap(l_u2))
     &     call errquit('WFAT: Cannot deallocate u2', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt0))
     &     call errquit('WFAT: Cannot deallocate dt_fbr', 1, MA_ERR)
      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp', 1, MA_ERR)
      
      end
!========================================================

      
!========================================================
!========================================================
      subroutine wfat_dipole_fbrt(m, p, n1, n2, dt)
cc     Calculates D, the matrix of the polar part of the dipole matrix in FBR
cc     basis. The elements of D depend on the value of p:
cc     1) if p = -1, then
cc          D_{ij} = < T_i^{m-1} | \sin\theta | T_j^m >,
cc     2) if p = 0, then
cc          D_{ij} = < T_i^m | \cos\theta | T_j^m >, or
cc     3) if p = 1, then
cc          D_{ij} = < T_i^{m+1} | \sin\theta | T_j^m >,
cc     where the size of D is n1 by n2, and T_i^m is the polar FBR function.
      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer m
      integer p      ! -1, 0, +1
      integer n1     ! m-1, m, m+1
      integer n2     ! m
      double precision dt(n1,n2)

      integer f, l_min, l_max, nl, l_gn, k_gn, j, l2, i1, i2
      integer i_gn1, i_gn2
      double precision nn, assocleg_norm
      external assocleg_norm
      

      if (m < 0)
     &     call errquit('WFAT: The argument m of wfat_dipole_fbrt '//
     &     'must be zero or a positive integer.', 85, UERR)
      if (m == 0 .and. p == -1)
     &     call errquit('WFAT: If the argument m of wfat_dipole_fbrt '//
     &     'has the value of zero, then the allowed values of p is 0 '//
     &     'or 1', 62, UERR)
      if (p /= -1 .and. p /= 0 .and. p /= 1)
     &     call errquit('WFAT: The argument p of wfat_dipole_fbrt '//
     &     'must be equal to either -1, 0, or 1.', 59, UERR)


      l_max = (m+n2-1) + 1
      l_min = abs(m+p)
      nl =  l_max - l_min + 1
      if (.not. MA_alloc_get(MT_DBL, nl, 'gaunt', l_gn, k_gn))
     &     call errquit('WFAT: gaunt allocation failed', 59, MA_ERR)
      nn = assocleg_norm(1,p)
      

c     ==== Calculate < P_{l1,m+p} | P_{1,p} | P_{l2,m} > ====
      dt = 0.0d0
      do j = 1, n2
         l2 = m + j - 1

c        Obtain < P_{l2,m} | P_{1,p} | P_{l1,m+p} > for l1=[l_min,l2+1]
c        where l_min = max(|l2-1|,|m+p|)
         call wfat_gauntp_lscan(l2, m, 1, p, dbl_mb(k_gn))
         l_min = max(abs(l2-1), abs(m+p))
         l_max = l2 + 1
         nl = l_max - l_min + 1
         
         i1 = l_min - (m+p) + 1
         i2 = i1 + nl - 1
         if (i2 > n2) then
            i2 = n2
            nl = i2 - i1 + 1
         endif
         i_gn1 = k_gn
         i_gn2 = i_gn1 + nl - 1

         if (p == -1) then
            f = 2
         else if (p == 0) then
            f = 1
         else if (p == 1) then
            f = -1
         endif

c        < P_{l1,m+p} | P_{1,p} | P_{l2,m} > = < P_{l2,m} | P_{1,p} | P_{l1,m+p} >
         dt(i1:i2, j) = f * sqrt(nn) * dbl_mb(i_gn1:i_gn2)
      enddo
      
      
      if (.not. MA_free_heap(l_gn))
     &     call errquit('WFAT: gaunt deallocation failed', 61, MA_ERR)
c      call errquit('stop here', 9, UERR)
      
      end
!========================================================


!========================================================
!========================================================
      subroutine wfat_dipole_dvrt_num(m, p, n1, n2, dt)

      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer p
      integer n1
      integer n2
      double precision dt(n1,n2)

      integer i, j, l_t, k_t, l_tw, k_tw, l_dv, k_dv
      double precision x, w, v, dvi, deg2rad
      parameter(deg2rad=0.01745329251994d0)


      if (p /= -1 .and. p /= 0 .and. p /= 1)
     &     call errquit('WFAT: The argument p of wfat_dipole_dvrt_num'//
     &     ' must be equal to either -1, 0, or 1.', 59, UERR)

      
      if (.not. ma_alloc_get(MT_DBL, n2, 'tquad', l_t, k_t))
     &     call errquit('WFAT: Cannot allocate tquad', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, n2, 'tweight', l_tw, k_tw))
     &     call errquit('WFAT: Cannot allocate tweight', 1, MA_ERR)
      call wfat_assocleg_quad(-1.0d0, 1.0d0, abs(m), n2, dbl_mb(k_t),
     &     dbl_mb(k_tw))

      if (.not. ma_alloc_get(MT_DBL, n1, 'dvr_i', l_dv, k_dv))
     &     call errquit('WFAT: Cannot allocate dvr_i', 1, MA_ERR)

      
      do j = 1, n2
         x = dbl_mb(k_t+j-1)
         w = dbl_mb(k_tw+j-1)
         call wfat_eval_dvrt(m+p, n1, acos(x)/deg2rad, dbl_mb(k_dv))
         do i = 1, n1
            dvi = dbl_mb(k_dv+i-1)
            if (abs(p) == 1) then
               v = sqrt(1-x**2)
            else if (p == 0) then
               v = x
            endif
            
            dt(i,j) = sqrt(w) * (1-x**2)**(-abs(m)/2.0d0) * dvi * v
         enddo
      enddo


      if (.not. ma_free_heap(l_t))
     &     call errquit('WFAT: Cannot deallocate tquad', 1, MA_ERR)
      if (.not. ma_free_heap(l_tw))
     &     call errquit('WFAT: Cannot deallocate tweight', 1, MA_ERR)
      if (.not. ma_free_heap(l_dv))
     &     call errquit('WFAT: Cannot deallocate dvr_i', 1, MA_ERR)
      
      end
!========================================================
      

!========================================================
!========================================================
      subroutine wfat_dipole_dvrt_num_OLD(m, n, dt)

      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"
      integer m
      integer n
      double precision dt(n,n)

      integer i, l_t, k_t, l_tw, k_tw

      if (.not. ma_alloc_get(MT_DBL, n, 'tquad', l_t, k_t))
     &     call errquit('WFAT: Cannot allocate tquad', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, n, 'tweight', l_tw, k_tw))
     &     call errquit('WFAT: Cannot allocate tweight', 1, MA_ERR)
      call wfat_assocleg_quad(-1.0d0, 1.0d0, abs(m), n, dbl_mb(k_t),
     &     dbl_mb(k_tw))

      dt = 0.0d0
      do i = 1, n
         dt(i,i) = dbl_mb(k_t+i-1)
      enddo
      
      if (.not. ma_free_heap(l_t))
     &     call errquit('WFAT: Cannot deallocate tquad', 1, MA_ERR)
      if (.not. ma_free_heap(l_tw))
     &     call errquit('WFAT: Cannot deallocate tweight', 1, MA_ERR)
      
      end
!========================================================


!========================================================
!========================================================
      subroutine wfat_2a1e_polb(m, sphi, orbid, nr1, nt1, nr2, nt2,
     &     nr3, nt3, ke, erg1, erg2, erg3, cf1, cf2, cf3, numt, polb)

      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      integer m
      integer sphi         ! [input], sphi=1 for \cos(m\phi) and sphi=2 for \sin(m\phi)
      integer orbid
      integer nr1, nt1
      integer nr2, nt2
      integer nr3, nt3
      double precision ke
      double precision erg1(1:*)                ! [input], unused when m=0
      double precision erg2(nr2*nt2)
      double precision erg3(nr3*nt3) 
      double precision cf1(1:*)                 ! [input], m-1, unused when m=0
      double precision cf2(nr2*nt2, nr2*nt2)    ! [input], m
      double precision cf3(nr3*nt3, nr3*nt3)    ! [input], m+1
      logical numt
      double precision polb(3)  ! [output],

      integer l_dr1, l_dt1, l_dr3, l_dt3, l_dr, l_dt
      integer k_dr1, k_dt1, k_dr3, k_dt3, k_dr, k_dt
      integer io, io0, io_, sphi_
      double precision dpo


      io0 = orbid

      
c     ==== m-1,m ====
      if (m >= 1) then
         if (.not. ma_alloc_get(MT_DBL, nr1*nr2, 'dr_dvr1', l_dr1, 
     &        k_dr1))
     &        call errquit('WFAT: Cannot allocate dr_dvr1', 1, MA_ERR)
         call wfat_dipole_dvrr(m-1, nr1, m, nr2, ke, dbl_mb(k_dr1))
   
         if (.not. ma_alloc_get(MT_DBL, nt1*nt2, 'dt_dvr1', l_dt1,
     &        k_dt1))
     &        call errquit('WFAT: Cannot allocate dt_dvr1', 1, MA_ERR)
         if (.not. numt) then
            call wfat_dipole_dvrt(m, -1, nt1, nt2, dbl_mb(k_dt1))
         else
            call wfat_dipole_dvrt_num(m, -1, nt1, nt2, dbl_mb(k_dt1))
         endif
      endif

      
c     ==== m+1,m ====
      if (.not. ma_alloc_get(MT_DBL, nr3*nr2, 'dr_dvr3', l_dr3, k_dr3))
     &     call errquit('WFAT: Cannot allocate dr_dvr3', 1, MA_ERR)
      call wfat_dipole_dvrr(m+1, nr3, m, nr2, ke, dbl_mb(k_dr3))

      if (.not. ma_alloc_get(MT_DBL, nt3*nt2, 'dt_dvr3', l_dt3, k_dt3))
     &     call errquit('WFAT: Cannot allocate dt_dvr3', 1, MA_ERR)
      if (.not. numt) then
         call wfat_dipole_dvrt(m, 1, nt3, nt2, dbl_mb(k_dt3))
      else
         call wfat_dipole_dvrt_num(m, 1, nt3, nt2, dbl_mb(k_dt3))
      endif
      

c     ==== xx component of the static polarizability ====
      polb(1) = 0.0d0
      if (m >= 1) then
         do io = 1, nr1*nt1
            io_ = 1 + (io-1)*nr1*nt1
            call wfat_2a1e_dpo('x', m, -1, nr1, nt1, nr2, nt2, sphi, 
     &           sphi, ke, cf1(io_), cf2(:,io0), .True., 
     &           dbl_mb(k_dr1), .True., dbl_mb(k_dt1), numt, dpo)
            polb(1) = polb(1) + dpo**2/(erg2(io0)-erg1(io))
         enddo
      endif
      do io = 1, nr3*nt3
         call wfat_2a1e_dpo('x', m, 1, nr3, nt3, nr2, nt2, sphi, 
     &        sphi, ke, cf3(:,io), cf2(:,io0), .True., dbl_mb(k_dr3),
     &        .True., dbl_mb(k_dt3), numt, dpo)
         polb(1) = polb(1) + dpo**2/(erg2(io0)-erg3(io))
      enddo


c     ==== yy component of the static polarizability ====
      polb(2) = 0.0d0
      if (sphi==1) sphi_ = 2
      if (sphi==2) sphi_ = 1
      if (m >= 1) then
         do io = 1, nr1*nt1
            io_ = 1 + (io-1)*nr1*nt1
            call wfat_2a1e_dpo('y', m, -1, nr1, nt1, nr2, nt2, sphi_, 
     &           sphi, ke, cf1(io_), cf2(:,io0), .True., 
     &           dbl_mb(k_dr1), .True., dbl_mb(k_dt1), numt, dpo)
            polb(2) = polb(2) + dpo**2/(erg2(io0)-erg1(io))
         enddo
      endif
      do io = 1, nr3*nt3
         call wfat_2a1e_dpo('y', m, 1, nr3, nt3, nr2, nt2, sphi_, 
     &        sphi, ke, cf3(:,io), cf2(:,io0), .True., dbl_mb(k_dr3),
     &        .True., dbl_mb(k_dt3), numt, dpo)
         polb(2) = polb(2) + dpo**2/(erg2(io0)-erg3(io))
      enddo      
      

c     ==== m,m ====
      if (.not. ma_alloc_get(MT_DBL, nr2*nr2, 'dr_dvr', l_dr, k_dr))
     &     call errquit('WFAT: Cannot allocate dr_dvr', 1, MA_ERR)
      call wfat_dipole_dvrr(m, nr2, m, nr2, ke, dbl_mb(k_dr))

      if (.not. ma_alloc_get(MT_DBL, nt2*nt2, 'dt_dvr', l_dt, k_dt))
     &     call errquit('WFAT: Cannot allocate dt_dvr', 1, MA_ERR)
      if (.not. numt) then
         call wfat_dipole_dvrt(m, 0, nt2, nt2, dbl_mb(k_dt))
      else
         call wfat_dipole_dvrt_num(m, 0, nt2, nt2, dbl_mb(k_dt))
      endif

      
c     ==== zz component of the static polarizability ====
      polb(3) = 0.0d0
      do io = 1, nr2*nt2
         if (io /= io0) then
            call wfat_2a1e_dpo('z', m, 0, nr2, nt2, nr2, nt2, sphi, 
     &           sphi, ke, cf2(:,io), cf2(:,io0), .True., dbl_mb(k_dr),
     &           .True., dbl_mb(k_dt), numt, dpo)
            polb(3) = polb(3) + dpo**2/(erg2(io0)-erg2(io))
         endif
      enddo

      polb = -2 * polb

      
c     ==== Deallocations ====
      if (m >= 1) then
         if (.not. ma_free_heap(l_dr1))
     &        call errquit('WFAT: Cannot deallocate dr_dvr1', 1, MA_ERR)
         if (.not. ma_free_heap(l_dt1))
     &        call errquit('WFAT: Cannot deallocate dt_dvr1', 1, MA_ERR)
      endif
      if (.not. ma_free_heap(l_dr3))
     &     call errquit('WFAT: Cannot deallocate dr_dvr3', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt3))
     &     call errquit('WFAT: Cannot deallocate dt_dvr3', 1, MA_ERR)
      if (.not. ma_free_heap(l_dr))
     &     call errquit('WFAT: Cannot deallocate dr_dvr', 1, MA_ERR)
      if (.not. ma_free_heap(l_dt))
     &     call errquit('WFAT: Cannot deallocate dt_dvr', 1, MA_ERR)

      end
!========================================================
      

!========================================================
!========================================================
      subroutine wfat_2a1e_dpo(comp, m, p, nr1, nt1, nr2, nt2, k1, k2,
     &     ke, cf1, cf2, dr_in, dr, dt_in, dt, numt, dpo)
cc     Calculates D, the matrix of the polar part of the dipole matrix in FBR
cc     basis. The elements of D depend on the value of p:
cc          D_{ij} = < u_{i',m+p,k'} | comp | u_{imk} >,
cc          u_{i01}(r,\theta,\phi) = f_i(r,\theta)/\sqrt(2*pi)
cc          u_{i02}(r,\theta,\phi) = 0
cc          u_{im1}(r,\theta,\phi) = f_i(r,\theta) cos(m\phi)/\sqrt(pi)
cc          u_{im2}(r,\theta,\phi) = f_i(r,\theta) sin(m\phi)/\sqrt(pi)
cc     where comp = x,y,z and the size of D is n1 by n2, and T_i^m is the polar FBR function.
      implicit none
#include "stdio.fh"      
#include "mafdecls.fh"
#include "errquit.fh"
      character*1 comp
      integer m
      integer p
      integer nr1, nt1
      integer nr2, nt2
      integer k1
      integer k2
      double precision ke
      double precision cf1(nr1*nt1)    ! [input], m-1 or m+1
      double precision cf2(nr2*nt2)    ! [input], m
      logical dr_in
      double precision dr(1:*)         ! [input], required only if dr_in=.true., regardless of dr_in, its minimum size should be nr1*nr2.
      logical dt_in
      double precision dt(1:*)         ! [input], required only if dt_in=.true., regardless of dt_in, its minimum size should be nt1*nt2.
      logical numt                     ! [input], meaningful only if dt_in=.true.
      double precision dpo             ! [output],

      integer l_tm, k_tm
      integer i1, i2, ic1, ic2, i_dr, delkron
      double precision s0, cc
      external delkron


      if (comp == 'x' .and. (abs(p) /= 1 .or. k1 /= k2)) then
         dpo = 0.0d0
         return
      endif
      if (comp == 'y' .and. (abs(p) /= 1 .or. k1 == k2)) then
         dpo = 0.0d0
         return
      endif
      if (comp == 'z' .and. (p /= 0 .or. k1 /= k2)) then
         dpo = 0.0d0
         return
      endif
      

      if (m < 0)
     &     call errquit('WFAT: The argument m of wfat_2a1e_dpo '//
     &     'must be zero or a positive integer.', 85, UERR)
      if (m == 0 .and. p == -1)
     &     call errquit('WFAT: If the argument m of wfat_2a1e_dpo '//
     &     'has the value of zero, then the allowed values of p is 0 '//
     &     'or 1.', 62, UERR)
      if (k1 /= 1 .and. k1 /= 2)
     &     call errquit('WFAT: The argument k1 of wfat_2a1e_dpo '//
     &     'must be equal to either 1 or 2.', 59, UERR)
      if (k2 /= 1 .and. k2 /= 2)
     &     call errquit('WFAT: The argument k2 of wfat_2a1e_dpo '//
     &     'must be equal to either 1 or 2.', 59, UERR)
      if (comp /= 'x' .and. comp /= 'y' .and. comp /= 'z')
     &     call errquit('WFAT: The argument comp of wfat_2a1e_dpo '//
     &     'must be equal to either x, y, or z.', 59, UERR)
      if (p == 0 .and. nr1 /= nr2)
     &     call errquit('WFAT: If the argument p of wfat_2a1e_dpo is '//
     &     'equal to zero, then nr1 must be equal to nr2.', 60, UERR)
      if (p == 0 .and. nt1 /= nt2)
     &     call errquit('WFAT: If the argument p of wfat_2a1e_dpo is '//
     &     'equal to zero, then nt1 must be equal to nt2.', 60, UERR)

      
      if (.not. dr_in) then
         call wfat_dipole_dvrr(m+p, nr1, m, nr2, ke, dr)
      endif

      if (.not. dt_in) then
         if (numt) then
            call wfat_dipole_dvrt_num(m, p, nt1, nt2, dt)
         else
            call wfat_dipole_dvrt(m, p, nt1, nt2, dt)
         endif
      endif
      

      if (.not. ma_alloc_get(MT_DBL, nt1, 'tmp', l_tm, k_tm))
     &     call errquit('WFAT: Cannot allocate tmp', 1, MA_ERR)
      dpo = 0.0d0
      do i2 = 1, nr2
         ic2 = 1 + (i2-1)*nt2
         do i1 = 1, nr1
            ic1 = 1 + (i1-1)*nt1

            call wfat_dgemm_wrap('n', 'n', nt1, nt2, nt2, 1, 1.0d0, 
     &           dt, cf2(ic2), 0.0d0, dbl_mb(k_tm))
            s0 = sum(cf1(ic1:ic1+nt1-1) * dbl_mb(k_tm:k_tm+nt1-1))

            i_dr = 1 + (i2-1)*nr1 + (i1-1)
            dpo = dpo + s0 * dr(i_dr)
         enddo
      enddo

      call wfat_2a1e_dpo_phi(comp, k1, m+p, k2, m, cc)
      dpo = cc * dpo
      

c     ==== Deallocations ====
      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp', 1, MA_ERR)      
      
      end
!========================================================


!========================================================
!========================================================
      subroutine wfat_2a1e_dpo_phi(comp, k1, m1, k2, m2, cc)

      implicit none
      character*1 comp
      integer k1
      integer m1
      integer k2
      integer m2
      double precision cc

      double precision norm1, norm2
      integer delkron
      external delkron


c     ==== < cos/sin(m1\phi) | cos(\phi) | cos/sin(m2\phi) > / \pi ====
      if (comp == 'x') then
         if (k1 == k2) then
            if (k1 == 1) then
               cc = 0.5d0 * (delkron(m1,m2+1) + delkron(m1,m2-1) +
     &              delkron(m1,1-m2) + delkron(m1,-1-m2))
            else if (k1 == 2) then
               cc = 0.5d0 * (delkron(m1,m2+1) + delkron(m1,m2-1) -
     &              delkron(m1,1-m2) - delkron(m1,-1-m2))
            endif
         else
            cc = 0.0d0
         endif
c     ==== < cos/sin(m1\phi) | sin(\phi) | cos/sin(m2\phi) > / \pi ====
      else if (comp == 'y') then
         if (k1 /= k2) then
            if (k1 == 1) then
               cc = 0.5d0 * (delkron(m1,1-m2) - delkron(m1,-m2-1) -
     &              delkron(m1,m2+1) + delkron(m1,m2-1))
            else if (k1 == 2) then
               cc = 0.5d0 * (delkron(m1,m2+1) + delkron(m1,1-m2) -
     &              delkron(m1,m2-1) - delkron(m1,-m2-1))
            endif
         else
            cc = 0.0d0
         endif
c     ==== < cos/sin(m1\phi) | cos/sin(m2\phi) > / \pi ====
      else if (comp == 'z') then
         if (k1 == k2) then
            if (m1 == m2 .and. m1 == 0) then
               cc = 2.0d0
            else if (m1 == m2 .and. abs(m1) /= 0) then
               cc = 1.0d0
            else
               cc = 0.0d0
            endif
         else
            cc = 0.0d0
         endif
      endif


c     ==== Normalization_factor / \pi ====
      if (m1 == 0) then
         norm1 = sqrt(2.0d0)
      else
         norm1 = 1.0d0
      endif
      if (m2 == 0) then
         norm2 = sqrt(2.0d0)
      else
         norm2 = 1.0d0
      endif


      cc = cc / norm1 / norm2


      end
!========================================================
