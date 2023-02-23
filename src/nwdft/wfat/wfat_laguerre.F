c=======================================================
c=======================================================
      subroutine wfat_lf_quadrature_rule ( n, alpha, x, w )
!*****************************************************************************80
!
!! LF_QUADRATURE_RULE: Gauss-Laguerre quadrature rule for Lf(n,alpha,x);
!
!  Discussion:
!
!    The integral:
!
!      integral ( 0 <= x < +oo ) exp ( - x ) * x^alpha * f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of the X factor.
!    ALPHA must be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
      implicit none
      integer ( kind = 4 ) n
      real ( kind = 8 ) alpha
      real ( kind = 8 ) x(n)
      real ( kind = 8 ) w(n)
      
      real ( kind = 8 ) bj(n)
      integer ( kind = 4 ) i
      real ( kind = 8 ) i_r8
      real ( kind = 8 ) wfat_r8_gamma
      real ( kind = 8 ) zemu
      external wfat_r8_gamma
      
c
c     Define the zero-th moment.
c
      zemu = wfat_r8_gamma ( alpha + 1.0D+00 )
c   
c     Define the Jacobi matrix.
c   
      do i = 1, n
         i_r8 = real ( i, kind = 8 )
         bj(i) = i_r8 * ( i_r8 + alpha )
      end do
      bj(1:n) = sqrt ( bj(1:n) )
    
      do i = 1, n
         i_r8 = real ( i, kind = 8 )
         x(i) = 2.0D+00 * i_r8 - 1.0D+00 + alpha
      end do
    
      w(1) = sqrt ( zemu )
      w(2:n) = 0.0D+00
c   
c     Diagonalize the Jacobi matrix.
c   
      call wfat_imtqlx ( n, x, bj, w )
      w(1:n) = w(1:n)**2
    
      return
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_imtqlx ( n, d, e, z )
!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine.
!
!    It has been modified to produce the product Q' * Z, where Z is an input
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.
!    The changes consist (essentially) of applying the orthogonal
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
      implicit none
    
      integer ( kind = 4 ) n
    
      real ( kind = 8 ) b
      real ( kind = 8 ) c
      real ( kind = 8 ) d(n)
      real ( kind = 8 ) e(n)
      real ( kind = 8 ) f
      real ( kind = 8 ) g
      integer ( kind = 4 ) i
      integer ( kind = 4 ) ii
      integer ( kind = 4 ), parameter :: itn = 30
      integer ( kind = 4 ) j
      integer ( kind = 4 ) k
      integer ( kind = 4 ) l
      integer ( kind = 4 ) m
      integer ( kind = 4 ) mml
      real ( kind = 8 ) p
      real ( kind = 8 ) prec
      real ( kind = 8 ) r
      real ( kind = 8 ) s
      real ( kind = 8 ) z(n)
    
      prec = epsilon ( prec )
    
      if ( n == 1 ) then
        return
      end if
    
      e(n) = 0.0D+00
    
      do l = 1, n
         j = 0
         do
            do m = l, n

               if ( m == n ) then
                  exit
               end if
    
               if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) +
     &            abs ( d(m+1) ) ) ) then
                  exit
               end if

            end do
            p = d(l)
    
            if ( m == l ) then
               exit
            end if
    
            if ( itn <= j ) then
               write ( *, '(a)' ) ' '
               write ( *, '(a)' ) 'IMTQLX - Fatal error!'
               write ( *, '(a)' ) '  Iteration limit exceeded.'
               write ( *, '(a,i8)' ) '  J = ', j
               write ( *, '(a,i8)' ) '  L = ', l
               write ( *, '(a,i8)' ) '  M = ', m
               write ( *, '(a,i8)' ) '  N = ', n
               stop
            end if
    
            j = j + 1
            g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
            r =  sqrt ( g * g + 1.0D+00 )
            g = d(m) - p + e(l) / ( g + sign ( r, g ) )
            s = 1.0D+00
            c = 1.0D+00
            p = 0.0D+00
            mml = m - l
      
            do ii = 1, mml
               i = m - ii
               f = s * e(i)
               b = c * e(i)
    
               if ( abs ( g ) <= abs ( f ) ) then
                  c = g / f
                  r =  sqrt ( c * c + 1.0D+00 )
                  e(i+1) = f * r
                  s = 1.0D+00 / r
                  c = c * s
               else
                  s = f / g
                  r =  sqrt ( s * s + 1.0D+00 )
                  e(i+1) = g * r
                  c = 1.0D+00 / r
                  s = s * c
               end if
    
               g = d(i+1) - p
               r = ( d(i) - g ) * s + 2.0D+00 * c * b
               p = s * r
               d(i+1) = g + p
               g = c * r - b
               f = z(i+1)
               z(i+1) = s * z(i) + c * f
               z(i) = c * z(i) - s * f
            end do
    
            d(l) = d(l) - p
            e(l) = g
            e(m) = 0.0D+00
         end do
      end do
c    !
c    !  Sorting.
c    !
      do ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
    
         do j = ii, n
            if ( d(j) < p ) then
               k = j
               p = d(j)
            end if
         end do
    
         if ( k /= i ) then
            d(k) = d(i)
            d(i) = p
            p = z(i)
            z(i) = z(k)
            z(k) = p
         end if
      end do
    
      return
      end      
c=======================================================


c=======================================================
c=======================================================
      function wfat_r8_gamma ( x )
!*****************************************************************************80
!
!! WFAT_R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) WFAT_R8_GAMMA, the value of the function.
!
      implicit none
!
!  Coefficients for minimax approximation over (12, INF).
!
      real ( kind = 8 ), dimension ( 7 ) :: c = (/ 
     & -1.910444077728D-03, 
     &  8.4171387781295D-04, 
     & -5.952379913043012D-04, 
     &  7.93650793500350248D-04, 
     & -2.777777777777681622553D-03, 
     &  8.333333333333333331554247D-02, 
     &  5.7083835261D-03 /)
      real ( kind = 8 ), parameter :: eps = 2.22D-16
      real ( kind = 8 ) fact
      integer ( kind = 4 ) i
      integer ( kind = 4 ) n
      real ( kind = 8 ), dimension ( 8 ) :: p = (/ 
     &  -1.71618513886549492533811D+00, 
     &   2.47656508055759199108314D+01, 
     &  -3.79804256470945635097577D+02, 
     &   6.29331155312818442661052D+02, 
     &   8.66966202790413211295064D+02, 
     &  -3.14512729688483675254357D+04, 
     &  -3.61444134186911729807069D+04, 
     &   6.64561438202405440627855D+04 /)
      logical parity
      real ( kind = 8 ), parameter ::
     &     pi = 3.1415926535897932384626434D+00
      real ( kind = 8 ), dimension ( 8 ) :: q = (/ 
     &  -3.08402300119738975254353D+01, 
     &   3.15350626979604161529144D+02, 
     &  -1.01515636749021914166146D+03, 
     &  -3.10777167157231109440444D+03, 
     &   2.25381184209801510330112D+04, 
     &   4.75584627752788110767815D+03, 
     &  -1.34659959864969306392456D+05, 
     &  -1.15132259675553483497211D+05 /)
      real ( kind = 8 ) wfat_r8_gamma
      real ( kind = 8 ) res
      real ( kind = 8 ), parameter ::
     &     sqrtpi = 0.9189385332046727417803297D+00
      real ( kind = 8 ) sum
      real ( kind = 8 ) x
      real ( kind = 8 ), parameter :: xbig = 171.624D+00
      real ( kind = 8 ) xden
      real ( kind = 8 ), parameter :: xinf = 1.0D+30
      real ( kind = 8 ), parameter :: xminin = 2.23D-308
      real ( kind = 8 ) xnum
      real ( kind = 8 ) y
      real ( kind = 8 ) y1
      real ( kind = 8 ) ysq
      real ( kind = 8 ) z
    
      parity = .false.
      fact = 1.0D+00
      n = 0
      y = x
c    !
c    !  Argument is negative.
c    !
      if ( y <= 0.0D+00 ) then
         y = - x
         y1 = aint ( y )
         res = y - y1
    
         if ( res /= 0.0D+00 ) then
            if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
               parity = .true.
            end if
            fact = - pi / sin ( pi * res )
            y = y + 1.0D+00
         else
            res = xinf
            wfat_r8_gamma = res
            return
         end if
    
      end if
c    !
c    !  Argument is positive.
c    !
      if ( y < eps ) then
c    !
c    !  Argument < EPS.
c    !
         if ( xminin <= y ) then
            res = 1.0D+00 / y
         else
            res = xinf
            wfat_r8_gamma = res
            return
         end if
      else if ( y < 12.0D+00 ) then
         y1 = y
c    !
c    !  0.0 < argument < 1.0.
c    !
         if ( y < 1.0D+00 ) then
            z = y
            y = y + 1.0D+00
c    !
c    !  1.0 < argument < 12.0.
c    !  Reduce argument if necessary.
c    !
         else
            n = int ( y ) - 1
            y = y - real ( n, kind = 8 )
            z = y - 1.0D+00
         end if
c    !
c    !  Evaluate approximation for 1.0 < argument < 2.0.
c    !
         xnum = 0.0D+00
         xden = 1.0D+00
         do i = 1, 8
            xnum = ( xnum + p(i) ) * z
            xden = xden * z + q(i)
         end do

         res = xnum / xden + 1.0D+00
c    !
c    !  Adjust result for case  0.0 < argument < 1.0.
c    !
         if ( y1 < y ) then
            res = res / y1
c    !
c    !  Adjust result for case 2.0 < argument < 12.0.
c    !
         else if ( y < y1 ) then
            do i = 1, n
               res = res * y
               y = y + 1.0D+00
            end do
         end if
           
      else
c    !
c    !  Evaluate for 12.0 <= argument.
c    !
         if ( y <= xbig ) then
            ysq = y * y
            sum = c(7)
            do i = 1, 6
               sum = sum / ysq + c(i)
            end do
            sum = sum / y - y + sqrtpi
            sum = sum + ( y - 0.5D+00 ) * log ( y )
            res = exp ( sum )
         else
            res = xinf
            wfat_r8_gamma = res
            return
         end if
      end if
c    !
c    !  Final adjustments and return.
c    !
      if ( parity ) then
         res = - res
      end if
    
      if ( fact /= 1.0D+00 ) then
         res = fact / res
      end if
    
      wfat_r8_gamma = res
    
      return
      end      
c=======================================================
