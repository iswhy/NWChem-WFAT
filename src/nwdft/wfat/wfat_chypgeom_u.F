subroutine chgubi ( a, b, x, hu, id )

c*********************************************************************72
c
cc CHGUBI computes the confluent hypergeometric function with integer argument B.
c
c  Discussion:
c
c    This procedure computes the confluent hypergeometric function
c    U(a,b,x) with integer b ( b = Ò1,Ò2,... )
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    31 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, B, parameters.
c
c    Input, double precision X, the argument.
c
c    Output, double precision HU, the value of U(a,b,x).
c
c    Output, integer ID, the estimated number of significant digits.
c
      implicit none

      double precision a
      double precision a0
      double precision a1
      double precision a2
      double precision b
      double precision da1
      double precision da2
      double precision db1
      double precision db2
      double precision el
      double precision ga
      double precision ga1
      double precision h0
      double precision hm1
      double precision hm2
      double precision hm3
      double precision hmax
      double precision hmin
      double precision hu
      double precision hu1
      double precision hu2
      double precision hw
      integer id
      integer id1
      integer id2
      integer j 
      integer k
      integer m
      integer n
      double precision ps
      double precision r
      double precision rn
      double precision rn1
      double precision s0
      double precision s1
      double precision s2
      double precision sa
      double precision sb
      double precision ua
      double precision ub
      double precision x

      id = -100
      el = 0.5772156649015329D+00
      n = int ( abs ( b - 1 ) )
      rn1 = 1.0D+00
      rn = 1.0D+00
      do j = 1, n
        rn = rn * j
        if ( j .eq. n - 1 ) then
          rn1 = rn
        end if
      end do

      call psi ( a, ps )
      call gamma ( a, ga )

      if ( 0.0D+00 .lt. b ) then
        a0 = a
        a1 = a - n
        a2 = a1
        call gamma ( a1, ga1 )
        ua = ( - 1 ) ** ( n - 1 ) / ( rn * ga1 )
        ub = rn1 / ga * x ** ( - n )
      else
        a0 = a + n
        a1 = a0
        a2 = a
        call gamma ( a1, ga1 )
        ua = ( - 1 ) ** ( n - 1 ) / ( rn * ga ) * x ** n
        ub = rn1 / ga1
      end if

      hm1 = 1.0D+00
      r = 1.0D+00
      hmax = 0.0D+00
      hmin = 1.0D+300

      do k = 1, 150
        r = r * ( a0 + k - 1.0D+00 ) * x / ( ( n + k ) * k )
        hm1 = hm1 + r
        hu1 = abs ( hm1 )
        hmax = max ( hmax, hu1 )
        hmin = min ( hmin, hu1 )
        if ( abs ( hm1 - h0 ) .lt. abs ( hm1 ) * 1.0D-15 ) then
          go to 10
        end if
        h0 = hm1
      end do

10    continue

      da1 = log10 ( hmax )
      if ( hmin .ne. 0.0D+00 ) then
        da2 = log10 ( hmin )
      end if
      id = 15 - int ( abs ( da1 - da2 ) )
      hm1 = hm1 * log ( x )
      s0 = 0.0D+00
      do m = 1, n
        if ( 0.0D+00 .le. b ) then
          s0 = s0 - 1.0D+00 / m
        else
          s0 = s0 + ( 1.0D+00 - a ) / ( m * ( a + m - 1.0D+00 ) )
        end if
      end do
      hm2 = ps + 2.0D+00 * el + s0
      r = 1.0D+00
      hmax = 0.0D+00
      hmin = 1.0D+300
      do k = 1, 150
        s1 = 0.0D+00
        s2 = 0.0D+00
        if ( 0.0D+00 .lt. b ) then
          do m = 1, k
            s1 = s1 - ( m + 2.0D+00 * a - 2.0D+00 ) 
     &        / ( m * ( m + a - 1.0D+00 ) )
          end do
          do m = 1, n
            s2 = s2 + 1.0D+00 / ( k + m )
          end do
        else
          do m = 1, k + n
            s1 = s1 + ( 1.0D+00 - a ) / ( m * ( m + a - 1.0D+00 ) )
          end do
          do m = 1, k
            s2 = s2 + 1.0D+00 / m
          end do
        end if
        hw = 2.0D+00 * el + ps + s1 - s2
        r = r * ( a0 + k - 1.0D+00 ) * x / ( ( n + k ) * k )
        hm2 = hm2 + r * hw
        hu2 = abs ( hm2 )
        hmax = max ( hmax, hu2 )
        hmin = min ( hmin, hu2 )

        if ( abs ( ( hm2 - h0 ) / hm2 ) .lt. 1.0D-15 ) then
          go to 20
        end if
        h0 = hm2
      end do

20    continue

      db1 = log10 ( hmax )
      if ( hmin .ne. 0.0D+00 ) then
        db2 = log10 ( hmin )
      end if
      id1 = 15 - int ( abs ( db1 - db2 ) )
      id = min ( id, id1 )

      if ( n .eq. 0 ) then
        hm3 = 0.0D+00
      else
        hm3 = 1.0D+00
      end if

      r = 1.0D+00
      do k = 1, n - 1
        r = r * ( a2 + k - 1.0D+00 ) / ( ( k - n ) * k ) * x
        hm3 = hm3 + r
      end do

      sa = ua * ( hm1 + hm2 )
      sb = ub * hm3
      hu = sa + sb

      if ( sa .ne. 0.0D+00 ) then
        id1 = int ( log10 ( abs ( sa ) ) )
      end if

      if ( hu .ne. 0.0D+00 ) then
        id2 = int ( log10 ( abs ( hu ) ) )
      end if

      if ( sa * sb .lt. 0.0D+00 ) then
        id = id - abs ( id1 - id2 )
      end if

      return
      end
c =================
