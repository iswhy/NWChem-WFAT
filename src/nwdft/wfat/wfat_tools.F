c=======================================================
c=======================================================
      subroutine wfat_dtranspose(m, n, a, b)

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"      
      integer m, n
      double precision a(m,n)
      double precision b(n,m)

      integer l_tm, k_tm
      integer i, i_tm


      if (.not. ma_alloc_get(mt_dbl, n*m, 'tmp', l_tm, k_tm))      ! 1)
     &     call errquit('WFAT: Cannot allocate tmp in wfat_dtranspose',
     &     1, MA_ERR)
      do i = 1, m
         i_tm = k_tm + (i-1)*n
         dbl_mb(i_tm:i_tm+n-1) = a(i,1:n)
      enddo

      do i = 1, m
         i_tm = k_tm + (i-1)*n
         b(1:n,i) = dbl_mb(i_tm:i_tm+n-1)
      enddo
      if (.not. ma_free_heap(l_tm))
     &     call errquit('WFAT: Cannot deallocate tmp in '//
     &     'wfat_dtranspose', 1, MA_ERR)

      end
      
c     NOTES:
c       1) Temporary array is needed to anticipate a situation where a and b 
c          point to the same memory address, i.e. when a and b in the calling
c          procedure are the same array.
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_dsyev_wrap(n, a, w)

      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"
      integer n                 ! [input] order of matrix a 
      double precision a(n,n)   ! [in/output] On input, matrix to be diagonalized. On output, the eigenvectors.
      double precision w(n)     ! [output] eigenvalues

      integer lwork, info
      integer l_wo, k_wo
      

c     ==== Workspace query ====
      if (.not. ma_alloc_get(MT_DBL, 1, 'work', l_wo, k_wo))
     &     call errquit('WFAT: Cannot allocate work for workspace '//
     &     'query step of wfat_dsyev_wrap.', 1, UERR)

      call dsyev('v', 'u', n, a, n, w, dbl_mb(k_wo), -1, info)
      lwork = int( dbl_mb(k_wo) )

      if (.not. ma_free_heap(l_wo))
     &     call errquit('WFAT: Cannot deallocate work in workspace '//
     &     'query step of wfat_dsyev_wrap.', 1, UERR)


c     ==== Do the diagonalization ====
      if (.not. ma_alloc_get(MT_DBL, lwork, 'work', l_wo, k_wo))
     &     call errquit('WFAT: Cannot allocate work for the '//
     &     'solving step of wfat_dsyev_wrap.', 1, UERR)
      
      call dsyev('v', 'u', n, a, n, w, dbl_mb(k_wo), lwork, info)

      if (info /= 0)
     &     call errquit('WFAT: Error when solving for the solution '//
     &     'using wfat_dsyev_wrap.', 1, UERR)
      
      if (.not. ma_free_heap(l_wo))
     &     call errquit('WFAT: Cannot deallocate work in the '//
     &     'solving step of wfat_dsyev_wrap.', 1, UERR)

      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_dsymv_wrap(n, alpha, a, x, beta, y)

      implicit none
      integer n
      double precision alpha
      double precision a(n,n)
      double precision x(n)
      double precision beta
      double precision y(n)


      call dsymv('u', n, alpha, a, n, x, 1, beta, y, 1)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dgemm_wrap(ta, tb, ma, na, mb, nb, alpha, a, b,
     &     beta, c)

      implicit none
#include "errquit.fh"
      character*1 ta
      character*1 tb
      integer ma
      integer na
      integer mb
      integer nb
      double precision alpha
      double precision a(ma,na)
      double precision b(mb,nb)
      double precision beta
      double precision c(1:*)
      
      integer m, n, k, lda, ldb, ldc
      

      if (ta /= 'n' .and. ta /= 't')
     &     call errquit('wfat_dgemm_wrap: Unknown value of argument ta',
     &     1, UERR)
      if (tb /= 'n' .and. tb /= 't')
     &     call errquit('wfat_dgemm_wrap: Unknown value of argument tb',
     &     1, UERR)

c     m = no. of rows of op(A)
c     n = no. of columns of op(B)
c     k = no. of columns of op(A) = no. of rows of op(B)
      
      if (ta == 'n') then
         m = ma
         k = na
      else
         m = na
         k = ma
      endif

      if (tb == 'n') then
         n = nb
         if (k /= mb)
     &        call errquit('wfat_dgemm_wrap: Inconsistent array sizes',
     &        1, UERR)
      else
         n = mb
         if (k /= nb)
     &        call errquit('wfat_dgemm_wrap: Inconsistent array sizes',
     &        1, UERR)
      endif

      
      if (ta == 'n') then
         lda = max(1,m)
      else
         lda = max(1,k)
      endif

      if (tb == 'n') then
         ldb = max(1,k)
      else
         ldb = max(1,n)
      endif
      
      ldc = max(1,m)

      
      call dgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta,
     &     c, ldc)

c     The codes above are written such that the following conditions are always satisfied:
c       - When ta='n', LDA=ma and ka=k=na. When ta='t', LDA=ma and ka=m=na.
c       - When tb='n', LDB=mb and kb=n=nb. When tb='t', LDB=mb and kb=k=nb.
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_dsymm_wrap(m, n, alpha, a, b, beta, c)
c     c = alpha*a*b + beta*c
      implicit none
      integer m                   ! [input], size of a and no. of rows of b
      integer n                   ! [input], no. of columns of b
      double precision alpha      ! [input], 
      double precision a(m,m)     ! [input], 
      double precision b(m,n)     ! [input], 
      double precision beta       ! [input], 
      double precision c(m,n)     ! [output],   


      call dsymm('l', 'u', m, n, alpha, a, m, b, m, beta, c, m)


      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_dsysv_wrap(n, nrhs, a, b)

      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"      
      integer n
      integer nrhs
      double precision a(n,n)
      double precision b(n,nrhs)

      integer l_wo, k_wo
      integer lwork, ipiv(n), info
      

c     ==== Workspace query ====
      if (.not. ma_alloc_get(MT_DBL, 1, 'work', l_wo, k_wo))
     &     call errquit('WFAT: Cannot allocate work for workspace '//
     &     'query step of wfat_dsysv_wrap.', 1, UERR)

      call dsysv('u', n, nrhs, a, n, ipiv, b, n, dbl_mb(k_wo), -1, 
     &     info)

c      if (.not. info /= 0)
c     &     call errquit('WFAT: Error during workspace query in '//
c     &     'wfat_dsysv_wrap.', 1, UERR)
      lwork = int( dbl_mb(k_wo) )

      if (.not. ma_free_heap(l_wo))
     &     call errquit('WFAT: Cannot deallocate work in workspace '//
     &     'query step of wfat_dsysv_wrap.', 1, UERR)

      
c     ==== Solving the matrix equation ====
      if (.not. ma_alloc_get(MT_DBL, lwork, 'work', l_wo, k_wo))
     &     call errquit('WFAT: Cannot allocate work for the '//
     &     'solving step of wfat_dsysv_wrap.', 1, UERR)
      
      call dsysv('u', n, nrhs, a, n, ipiv, b, n, dbl_mb(k_wo), lwork,
     &     info)

      if (info /= 0)
     &     call errquit('WFAT: Error when solving for the solution '//
     &     'using wfat_dsysv_wrap.', 1, UERR)
      
      if (.not. ma_free_heap(l_wo))
     &     call errquit('WFAT: Cannot deallocate work in the '//
     &     'solving step of wfat_dsysv_wrap.', 1, UERR)

      end
c=======================================================            

      
c=======================================================
c=======================================================
      subroutine wfat_trace(n, a, tr)

      implicit none
      integer n
      double precision a(n,n)
      double precision tr

      integer i
      
      tr = 0.0d0
      do i = 1, n
         tr = tr + a(i,i)
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_project_mat(m, n, v, p)
c     DESCRIPTION :
c     Return the matrix of projection operator defined by
c        p = \sum_{i=1}^n  v(:,i) * v(:,i)^T
c     where the * sign indicates matrix product.
      
      implicit none
      integer m
      integer n
      double precision v(m,n)
      double precision p(m,m)


      p = matmul(v, transpose(v))

      end
c=======================================================      
      

c=======================================================
c=======================================================
      subroutine wfat_ne_expval_1eopr(nb, ns, nel, vcf, opr,
     &           val1e, valne)

      implicit none
      integer nb
      integer ns
      integer nel(ns)
      double precision vcf(nb, nb, ns)
      double precision opr(nb, nb)
      double precision val1e(nb, ns)
      double precision valne

      integer is, io
      

      valne = 0.0d0
      do is = 1, ns
         do io = 1, nb
            call wfat_1e_expval_1eopr(nb, vcf(:,io,is), opr,
     &                            val1e(io,is))
            if (io <= nel(is)) valne = valne + val1e(io,is)
         enddo
      enddo
      
      if (ns == 1) valne = 2 * valne

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_1e_expval_1eopr(nb, vcf, opr, val1e)

      implicit none
      integer nb
      double precision vcf(nb)
      double precision opr(nb, nb)
      double precision val1e

      integer i, j, f
      

      val1e = 0.0d0
      do j = 1, nb
         do i = 1, j
            if (i /= j) then
               f = 2
            else
               f = 1
            endif

            val1e = val1e + f * vcf(j) * vcf(i) * opr(i,j)
         enddo
      enddo
      

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_mat_element(n, u, v, o, val)
c   val = <u|o|v>
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer n
      double precision u(n)
      double precision v(n)
      double precision o(n,n)
      double precision val

      integer l_sc, k_sc


      if (.not. ma_alloc_get(mt_dbl, n, 'scr', l_sc, k_sc))
     &     call errquit('WFAT: Cannot allocate scr in wfat_mat_element',
     &     1, MA_ERR)
      
      call wfat_dsymv_wrap(n, 1.0d0, o, v, 0.0d0, dbl_mb(k_sc))
      val = dot_product(u, dbl_mb(k_sc:k_sc+n-1))

      if (.not. ma_free_heap(l_sc))
     &     call errquit('WFAT: Cannot deallocate scr in '//
     &     'wfat_mat_element', 1, MA_ERR )

      end
c=======================================================


c===============================================================
c===============================================================
      subroutine wfat_ga_mat_element(g_v, g_w, g_a, d)
c  d = < g_v | g_a | g_w >
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"      
      integer g_v
      integer g_w
      integer g_a
      double precision d

      logical l
      integer typev, typew, typea
      integer m_v, n_v, m_w, n_w, m_a, n_a
      integer g_tmp


      call ga_inquire(g_v, typev, m_v, n_v)
      call ga_inquire(g_w, typew, m_w, n_w)
      call ga_inquire(g_a, typea, m_a, n_a)


      if (typev == typew) then
         if (typew /= typea) then
            call errquit('WFAT: Type of argument ''g_w'' and of '//
     &           'argument ''g_a'' of ''wfat_ga_vrepresent'' must be '//
     &           'the same.', 1, UERR)
         endif
      else
         call errquit('WFAT: Type of argument ''g_v'' and of '//
     &        'argument ''g_w'' of ''wfat_ga_vrepresent'' must be '//
     &        'the same.', 1, UERR)
      endif
      

      if (n_w /= 1) then
         call errquit('WFAT: Argument ''g_w'' of '//
     &        '''wfat_ga_mat_element'' must be a column vector.',
     &        1, UERR)
      endif
      if (.not. (m_v == 1 .or. n_v == 1)) then
         call errquit('WFAT: Argument ''g_v'' of '//
     &        '''wfat_ga_mat_element'' must be a 1D array, '//
     &        'i.e. a column or row vector.', 1, UERR)
      endif
      
      if (n_a /= m_w) then
         call errquit('WFAT: In ''wfat_ga_mat_element'', the number '//
     &    'of columns of ''g_a'' must be equal to the number of '//
     &    'rows of ''g_w''.', 1, UERR)
      endif

      if (m_v*n_v /= m_a) then
         call errquit('WFAT: In ''wfat_ga_mat_element'', the number '//
     &        'of elements of ''g_v'' must be equal to the number of '//
     &        'rows of ''g_a''.', 1, UERR)
      endif
      

         

      if (m_a == n_a) then
         l = ga_duplicate(g_w, g_tmp, 'tmp')
      else
         l = ga_create(mt_dbl, m_a, 1, 'tmp', 0, 0, g_tmp)
      endif
      if (.not. l)
     &     call errquit('wfat_ga_mat_element: Cannot create tmp.',
     &     1, GA_ERR)
      call ga_zero(g_tmp)

      
c     ==== g_a * g_w ====
      call ga_dgemm('n', 'n', m_a, n_w, m_w, 1.0d0, g_a, g_w, 0.0d0,
     &     g_tmp)

c     ==== g_v * (g_a * g_w) ====
      d = ga_ddot(g_v, g_tmp)

      
      if (.not. ga_destroy(g_tmp))
     &     call errquit('wfat_ga_vrepresent: Cannot destroy tmp.',
     &     1, GA_ERR)
      
      end
c===============================================================
      

c===============================================================
c===============================================================
      subroutine wfat_ga_vrepresent(g_v, g_a, g_b)
c   g_b = g_v^T . g_a . g_v
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"      
      integer g_v     ! [input]
      integer g_a     ! [input]
      integer g_b     ! [output]

      integer typev, m_v, n_v, typea, m_a, n_a, typeb, m_b, n_b
      integer g_tmp
      
      
      call ga_inquire(g_v, typev, m_v, n_v)
      call ga_inquire(g_a, typea, m_a, n_a)
      call ga_inquire(g_b, typeb, m_b, n_b)

      if (typea == typeb) then
         if (typeb /= typev) then
            call errquit('WFAT: Type of argument g_b and of '//
     &           'argument g_v of wfat_ga_vrepresent must be the same.',
     &           1, UERR)
         endif
      else
         call errquit('WFAT: Type of argument g_a and of argument '//
     &        'g_b of wfat_ga_vrepresent must be the same.', 1, UERR)
      endif
      
      if (n_a /= m_v) then
         call errquit('WFAT: The no. of columns of the matrix '//
     &        'represented by argument g_a of wfat_ga_vrepresent '//
     &        'must be equal the no. of rows of the matrix '//
     &        'represented by g_v', 1, GA_ERR)
      endif
      if (m_a /= n_a) then
         call errquit('WFAT: The matrix represented by argument g_a '//
     &        'of wfat_ga_vrepresent must be square', 1, GA_ERR)
      endif
      if (m_b /= n_b) then
         call errquit('WFAT: The matrix represented by argument g_b '//
     &        'of wfat_ga_vrepresent must be square', 1, GA_ERR)
      endif
      if (n_b /= n_v) then
         call errquit('WFAT: The no. of columns of the matrix '//
     &        'represented by argument g_b of wfat_ga_vrepresent '//
     &        'must be equal the no. of columns of the matrix '//
     &        'represented by g_v', 1, GA_ERR)
      endif
      

      if (.not. ga_create(mt_dbl, m_a, n_v, 'tmp', 0, 0, g_tmp))
     &     call errquit('wfat_ga_vrepresent: Cannot create tmp.',
     &     1, GA_ERR)
      call ga_zero(g_tmp)
      
      call ga_dgemm('n', 'n', m_a, n_v, n_a, 1.0d0, g_a, g_v, 0.0d0,
     &     g_tmp)
      call ga_dgemm('t', 'n', m_b, n_b, m_a, 1.0d0, g_v, g_tmp, 0.0d0,
     &     g_b)

      if (.not. ga_destroy(g_tmp))
     &     call errquit('wfat_ga_vrepresent: Cannot destroy tmp.',
     &     1, GA_ERR)

      end
c===============================================================      

            
c=======================================================
c=======================================================
      subroutine wfat_get_i2cart(ic, ic1, ic2)
c     Based on subroutine xc_drhonuc in $NWCHEM_TOP/src/nwdft/xc/xc_rhogen.F
      
      implicit none
#include "errquit.fh"
      integer ic
      integer ic1, ic2


      if (ic == 1) then
         ic1 = 1
         ic2 = 1
      else if (ic == 2) then
         ic1 = 1
         ic2 = 2
      else if (ic == 3) then 
         ic1 = 1
         ic2 = 3
      else if (ic == 4) then 
         ic1 = 2
         ic2 = 2
      else if (ic == 5) then 
         ic1 = 2
         ic2 = 3
      else if (ic == 6) then 
         ic1 = 3
         ic2 = 3
      else
         call errquit('WFAT: The chosen value for argument ic '//
     &        'of wfat_get_i2cart is unavailable.', 1, UERR)
      endif




cOLD      if (ic == 1) then
cOLD         ic1 = 1
cOLD         ic2 = 1
cOLD      else if (ic == 2) then
cOLD         ic1 = 2
cOLD         ic2 = 2
cOLD      else if (ic == 3) then 
cOLD         ic1 = 3
cOLD         ic2 = 3
cOLD      else if (ic == 4) then 
cOLD         ic1 = 1
cOLD         ic2 = 2
cOLD      else if (ic == 5) then 
cOLD         ic1 = 2
cOLD         ic2 = 3
cOLD      else if (ic == 6) then 
cOLD         ic1 = 1
cOLD         ic2 = 3
cOLD      else
cOLD         call errquit('WFAT: The chosen value for argument ic '//
cOLD     &        'of wfat_get_i2cart is unavailable.', 1, UERR)
cOLD      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_cart2i(ic1, ic2, ic)

      implicit none
#include "errquit.fh"
      integer ic1, ic2
      integer ic      



      if (ic1 == 1 .and. ic2 == 1) then
         ic = 1          
      else if (ic1 == 2 .and. ic2 == 2) then
         ic = 4          
      else if (ic1 == 3 .and. ic2 == 3) then
         ic = 6           
      else if ((ic1 == 1 .and. ic2 == 2) .or.
     &         (ic1 == 2 .and. ic2 == 1)) then
         ic = 2           
      else if ((ic1 == 2 .and. ic2 == 3) .or.
     &         (ic1 == 3 .and. ic2 == 2)) then
         ic = 5           
      else if ((ic1 == 1 .and. ic2 == 3) .or.
     &         (ic1 == 3 .and. ic2 == 1)) then
         ic = 3 
      else
         call errquit('WFAT: The chosen values for arguments ic1 '//
     &        'and ic2 of wfat_get_cart2i are unavailable.', 1, UERR)
      endif




      
cOLD      if (ic1 == 1 .and. ic2 == 1) then
cOLD         ic = 1          
cOLD      else if (ic1 == 2 .and. ic2 == 2) then
cOLD         ic = 2          
cOLD      else if (ic1 == 3 .and. ic2 == 3) then
cOLD         ic = 3           
cOLD      else if ((ic1 == 1 .and. ic2 == 2) .or.
cOLD     &         (ic1 == 2 .and. ic2 == 1)) then
cOLD         ic = 4           
cOLD      else if ((ic1 == 2 .and. ic2 == 3) .or.
cOLD     &         (ic1 == 3 .and. ic2 == 2)) then
cOLD         ic = 5           
cOLD      else if ((ic1 == 1 .and. ic2 == 3) .or.
cOLD     &         (ic1 == 3 .and. ic2 == 1)) then
cOLD         ic = 6 
cOLD      else
cOLD         call errquit('WFAT: The chosen values for arguments ic1 '//
cOLD     &        'and ic2 of wfat_get_cart2i are unavailable.', 1, UERR)
cOLD      endif

      
      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_get_digits(val, dig)

      implicit none
      integer val
      integer dig

      integer aval, i, rem, nmax
      parameter (nmax=10)
      

      aval = abs(val)

      do i = 1, nmax
         rem = mod(aval, 10**i)
         if (rem == aval) then
            dig = i
            exit
         endif
      enddo

      if (val < 0) dig = dig + 1

      end
c=======================================================      


c=======================================================
c=======================================================
      subroutine wfat_dprintmat(m, n, e, a)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
      integer m
      integer n
      character*1 e
      double precision a(m,n)

      character*6 nn
      character*30 efmt
      character*7 efmt0
      integer i


      if (n <= 0 .or. n >= 1000000) then
         call errquit('WFAT: The value of argument n of '//
     &        'wfat_dprintmat must be positive and less than '//
     &        '1000000.', 1, UERR)
      endif

      write(nn, '(i6)') n
      if (e == 'f') then
         efmt0 = "(i10)"
         efmt = "(i5, "//trim(adjustl(nn))//"f10.5)"
      else if (e == 'e') then
         efmt0 = "(i14)"
         efmt = "(i5, "//trim(adjustl(nn))//"e14.5)"
      else
         call errquit('WFAT: The choice of argument e of '//
     &                'wfat_dprintmat is unavailable.', 1, UERR)
      endif

      write(luout, '(5x)', advance='no')
      call flush(luout)
      do i = 1, n
         write(luout, efmt0, advance='no') i
         call flush(luout)
      enddo
      write(luout, '(1x)')
      
      do i = 1, m
         write(luout, efmt) i, a(i, 1:n)
      enddo

      write(luout, '(" ")') 
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_print_repeat(iunit, n, os, c)

      implicit none
#include "errquit.fh"      
      integer iunit, n, os
      character*1 c

      integer i


      if (os > 100)
     &     call errquit('WFAT: Argument os of wfat_print_repeat '//
     &     'cannot be larger than 100.', 1, UERR)

      call util_flush(iunit)
      do i = 1, os
         write(iunit, '(" ")', advance='no')
         call util_flush(iunit)
      enddo

      do i = 1, n
         write(iunit, '(a1)', advance='no') c
         call util_flush(iunit)
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_tf_to_yn(l, e, a)

      implicit none
      logical l       ! [input], input logical to be converted.
      logical e       ! [input], .true. for 'yes/no', .false. for 'y/n'.
      character*(*) a ! [output], result of conversion.


      if (l) then
         if (e)       a = 'yes'
         if (.not. e) a = 'y'
      else
         if (e)       a = 'no'
         if (.not. e) a = 'n'
      endif

      end
c=======================================================            


c=======================================================
c=======================================================
      subroutine wfat_def_intglname(rtdb, lvl, fname)

      implicit none
#include "errquit.fh"
      integer rtdb
      character*2 lvl
      character*(*) fname


      if (lvl == 'oe') then
         call util_file_prefix('owigl', fname)
      else if (lvl == 'me') then
         call util_file_prefix('mwigl', fname)
      else
         call errquit('WFAT: Unknown value of argument ''lvl'' '//
     &        'of ''wfat_def_intglname''', 1, UERR)
      endif
      
      call wfat_prepend_rtdbpath(rtdb, fname)

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_def_movecname(lvl, rtdb, fname)

      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      character*(*) lvl
      integer rtdb
      character*(*) fname


      if (lvl == 'oe') then
         call util_file_prefix('movecs', fname)
      else if (lvl == 'me0') then
         call util_file_prefix('ion0.movecs', fname)
      else if (lvl == 'me1') then
         call util_file_prefix('ion1.movecs', fname)
      else
         call errquit('WFAT: Argument ''lvl'' of wfat_def_movecname '//
     &        'is unknown.', 1, UERR)
      endif
      call wfat_prepend_rtdbpath(rtdb, fname)

      end
c=======================================================      

      
c=======================================================
c=======================================================
      subroutine wfat_prepend_rtdbpath(rtdb, fname)

      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"      
      integer rtdb              ! [input]
      character*(*) fname       ! [input/output]

      character*256 rtdb_name
      integer rtlen, rtlen0, nmlen0, i, pos



      if(.not. rtdb_getfname(rtdb, rtdb_name)) 
     &     call errquit('WFAT: Cannot get RTDB file name', 1, RTDB_ERR)
      
      rtlen = inp_strlen(rtdb_name)
      pos = 0
      do i = 1, rtlen
         if (rtdb_name(i:i) == '/' ) pos = i
      enddo

      rtlen0 = inp_strlen(trim(adjustl(rtdb_name(1:pos))))
      nmlen0 = inp_strlen(trim(adjustl(fname)))
      if (rtlen0 + nmlen0 > len(fname))
     &     call errquit('WFAT: The argument ''fname'' of '//
     &     'wfat_prepend_rtdbpath subroutine is too small to hold '//
     &     'the output string where the path of RTDB will have '//
     &     'been prepended.', 1, UERR)
      fname = trim(adjustl(rtdb_name(1:pos)))//adjustl(fname)
      
cforcheck      write(luout, '("Vector file : ", a256)') adjustl(fname)

      end
c=======================================================


c=======================================================
c=======================================================
      double precision function rise_fact(x, n)
      implicit none
      double precision x
      integer n
      
      integer i
      double precision f

      
      
      if (n == 0) then
         rise_fact = 1.0d0
      else
         f = 1.0d0
         do i = 1, n
            f = f * (x+i-1)
         enddo
         rise_fact = f
      endif
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_legendre(ut,N,x,L)
      implicit none
#include "errquit.fh"      
      character*1 ut
      integer N
      double precision x
      double precision L
      
      integer i
      double precision norm, L0, L1

      
      if (ut == 'y') then
         norm = sqrt((2*N+1.0d0)/2)
      else if (ut == 'n') then
         norm = 1.0d0
      else
         call errquit('WFAT: The argument ut for wfat_legendre '//
     &        'must be either ''y'' or ''n''.', 18, UERR)
      end if
      
      
      L0 = 1.0d0
      L1 = x
      
      if (N == 0) then
         L = L0
         return
      else if (N == 1) then
         L = L1
         return
      end if
      
c     Calculate the N-th order Legendre polynomial using the recursive formula
c     N.L_N(x) = (2N-1).x.L_(N-1)(x) - (N-1).L_(N-2)(x)
      do i = 2,N
         L = ((2*i-1)*x*L1 - (i-1)*L0)/i
         L0 = L1
         L1 = L
      end do
      
      L = norm*L

      end
c=======================================================
      
      
c=======================================================
c=======================================================
      subroutine wfat_lgl(a, b, n, node, weig)
c     This routine returns Legendre-Gauss-Lobatto nodes and weights.
c     a = left end node
c     b = right end node
c     n+1 = number of nodes including the end nodes
c     node = quadrature nodes inlcuding end nodes
c     weig = quadrature weights including those associated with the end nodes.
      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"
      integer n              ! The number of nodes and weights are N+1
      double precision a, b
      double precision node(0:n), weig(0:n)
      
      integer m, l_dm, k_dm, l_nd, k_nd, i, j, k, k_
      double precision lgd

      
      m = n - 2
      if (n < 2)
     &     call errquit('Argument n of wfat_lgl must be bigger than '//
     &     'or equal to 2', 73, UERR)

      
c     ==== Construct the dipole moment matrix ====
      if (.not. MA_alloc_get(MT_DBL, (m+1)**2, 'DM', l_dm, k_dm))
     &     call errquit('WFAT: DM allocation failed', 59, MA_ERR)
      if (.not. MA_alloc_get(MT_DBL, m+1, 'node', l_nd, k_nd))
     &     call errquit('WFAT: node allocation failed', 59, MA_ERR)
      dbl_mb(k_dm : k_dm+(m+1)**2-1) = 0.0d0
      do i = 0,m
         do j = 0,m
            if (j == i+1) then
               k = i*(m+1) + j + 1
               dbl_mb(k_dm+k-1) = sqrt(1.0d0 * j*(j+2)/(2*j+1)/(2*j+3))
               k_ = j*(m+1) + i + 1
               dbl_mb(k_dm+k_-1) = dbl_mb(k_dm+k-1)
            end if
         end do
      end do
cfordebug      call wfat_dprintmat(m+1, m+1, 'e', dbl_mb(k_dm))
      call wfat_dsyev_wrap(m+1, dbl_mb(k_dm), dbl_mb(k_nd))


c     ==== Calculate shifted nodes and weights ====
      node(1:n-1) = dbl_mb(k_nd:k_nd+(m+1)-1)         !node0(0:M)
      node(0) = -1.0d0
      node(n) = 1.0d0
      do i = 0,n
         call wfat_legendre('n', n, node(i), lgd)
         weig(i) = 2.0d0/n/(n+1)/lgd**2              ! Calculate weights
         node(i) = 0.5d0*((b-a)*node(i) + b + a)     ! Calculate shifted nodes   
         weig(i) = (b-a)/2.0d0*weig(i)               ! Calculate shifted weights
      end do

      
c     ==== Deallocations ====
      if (.not. ma_free_heap(l_dm))
     &     call errquit('WFAT: Cannot deallocate DM', 49, MA_ERR)
      if (.not. ma_free_heap(l_nd))
     &     call errquit('WFAT: Cannot deallocate node', 49, MA_ERR)
      end
c=======================================================
