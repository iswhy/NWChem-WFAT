c=======================================================
c=======================================================
      subroutine wfat_chgu_1pt(r, lmin_, lmax_, a_, b_, lmb, chgu)

      implicit none
#include "stdio.fh"      
      integer lmin_, lmax_
      double precision r, a_, b_, lmb
      double precision chgu(lmin_:lmax_)
      
      integer l, id
      double precision a, aa, b, b0, bb, z, chgu_1, chgu_2, chgu_3


c     ==== Prepare initial values for recurrence calculation ====
      z = lmb * r      
      a = a_ + (lmax_ + 1)
      b = b_ + 2*(lmax_ + 1)
      
      b0 = b - 0
c      call wfat_chgubi(a, b0, z, chgu_3, id)
      call wfat_chgu(a, b0, z, chgu_3, id)           ! 1)
      
      b0 = b - 1
c      call wfat_chgubi(a, b0, z, chgu_2, id)
      call wfat_chgu(a, b0, z, chgu_2, id)           ! 1)
      
      b0 = b - 2
c      call wfat_chgubi(a, b0, z, chgu_1, id)
      call wfat_chgu(a, b0, z, chgu_1, id)           ! 1)


c     ==== Calculate the confluent hypergeometric function by recurrence ====
      if (abs(chgu_3) < 5.0d2) then
         call wfat_chgu_1pt_recur(lmin_, lmax_, a_, b_, z,
     &                            chgu_1, chgu_2, chgu_3, chgu)
      else                   ! 2)
         do l = lmin_, lmax_
            aa = l + a_
            bb = 2*l + b_
c            call wfat_chgubi(aa, bb, z, chgu(l), id)
            call wfat_chgu(aa, bb, z, chgu(l), id)           ! 1)
         enddo
      endif
      

c     ==== Scan for too big values ====
      do l = lmin_, lmax_
c           if (chgf(l) > 1.0d+20 .or. isnan(chgf(l)) ) then
         if (isnan(chgu(l)) ) then        
c              print *, "rmax = ", r0, l, chgf(l)
            chgu(l) = 0.0d0
         endif
      enddo
      
      end
c     NOTES:
c     1) wfat_chgu is used instead of wfat_chgubi that is more specifically
c        used for integer b because the latter produced wrong values for
c        large arguments as is evidenced by calling wfat_test_chg.
c     2) When the value of the U(a,b,x) function is large, don't use
c        recurrence, use direct calculation instead. This is because it has
c        been found that recurrence method accumulates error faster when 
c        the values of U(a,b,x) are large.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_chgu_1pt_recur(m, ll, a0, b0, z, chgu_1_, chgu_2_,
     &                               chgu_3_, chgu)
c      !! DESCRIPTION:
c      !!   chgu_1_ = U(ll+1+a0, 2ll+0+b0, z)
c      !!   chgu_2_ = U(ll+1+a0, 2ll+1+b0, z)
c      !!   chgu_3_ = U(ll+1+a0, 2ll+2+b0, z)
c      !!   chgu    = U(l+a0,    2l+b0,    z)  with l = m, m+1, ..., ll
c      !!   m >= 0
      
      implicit none
#include "errquit.fh"
      integer m, ll
      double precision a0, b0, z
      double precision chgu_1_, chgu_2_, chgu_3_
      double precision chgu(m:ll)
      
      integer l
      double precision a, b, chgu_1, chgu_2, chgu_3, chgu_4


      if (m < 0) then
         call errquit('WFAT: Argument m of wfat_chgu_1pt_recur '//
     &                'cannot be negative.', 1, UERR)
      endif
      
            
      do l = ll, m, -1
         a = a0 + l
         b = b0 + 2*l + 1       
      
         !!==== chgu_4 ====!! Use eq. (13.4.18) of Ref. [1]
         chgu_4 = z*chgu_3_ - (b-a-1)*chgu_2_
      
         !!==== chgu_3 ====!! Use eq. (13.4.18) of Ref. [1]
         b = b - 1
         chgu_3 = z*chgu_2_ - (b-a-1)*chgu_1_
         chgu(l) = chgu_3
      
         !!==== Avoid unnecessary calculations ====!!
         if (l == m) exit
      
         !!==== chgu_2 ====!! Use eq. (13.4.16) of Ref. [1]
         b = b - 1
         chgu_2 = (z*chgu_4 - (z+b)*chgu_3) / (a-b)

         !!==== chgu_1 ====!! Use eq. (13.4.16) of Ref. [1]
         b = b - 1
         chgu_1 = (z*chgu_3 - (z+b)*chgu_2) / (a-b)

         !!==== For next iteration ====!!
         chgu_3_ = chgu_3
         chgu_2_ = chgu_2
         chgu_1_ = chgu_1
      enddo
            
      end
c     References:
c        [1] Handbook of Mathematical Functions by Abramowitz and Stegun.
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_chgu(a, b, x, hu, md)

      implicit none
#include "stdio.fh"
      double precision a
      double precision b
      double precision x
      double precision hu
      integer md


      if (x < 3.0d0 .and. abs(int(b)-b) <= 1.0d-14 ) then
         call wfat_chgubi_smallx(a, b, x, hu, md)
      else
         call wfat_chgu_szjj(a, b, x, hu, md)
      endif

      end
c=======================================================
      

c=======================================================
c=======================================================
      subroutine wfat_chgu_szjj ( a, b, x, hu, md )
c*********************************************************************72
c
cc CHGU computes the confluent hypergeometric function U(a,b,x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    27 July 2012
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
c    Output, double precision HU, U(a,b,x).
c
c    Output, integer MD, the method code.
c
      implicit none

      double precision a
      double precision a00
      double precision aa
      double precision b
      double precision b00
      logical bl1
      logical bl2
      logical bl3
      logical bn
      double precision hu
      double precision hu1
      integer id
      integer id1
      logical il1
      logical il2
      logical il3
      integer md
      double precision x

      aa = a - b + 1.0D+00
      il1 = a .eq. int ( a ) .and. a .le. 0.0D+00
      il2 = aa .eq. int ( aa ) .and. aa .le. 0.0D+00
      il3 = abs ( a * ( a - b + 1.0D+00 ) ) / x .le. 2.0D+00
      bl1 = x .le. 5.0D+00 .or. ( x .le. 10.0D+00 .and. a .le. 2.0D+00 )
      bl2 = ( 5.0D+00 .lt. x .and. x .le. 12.5D+00 ) .and. 
     &  ( 1.0D+00 .le. a .and. a + 4.0D+00 .le. b )
      bl3 = 12.5D+00 .lt. x .and. 5.0D+00 .le. a .and. 
     &  a + 5.0D+00 .le. b
      bn = b .eq. int ( b ) .and. b .ne. 0.0D+00
      id1 = -100

      if ( b .ne. int ( b ) ) then
        call wfat_chgus ( a, b, x, hu, id1 )
        md = 1
        if ( 6 .le. id1 ) then
          return
        end if
        hu1 = hu
      end if

      if ( il1 .or. il2 .or. il3 ) then
        call wfat_chgul ( a, b, x, hu, id )
        md = 2
        if ( 6 .le. id ) then
          return
        end if
        if ( id .lt. id1 ) then
          md = 1
          id = id1
          hu = hu1
        end if
      end if

      if ( 0.0D+00 .le. a ) then
        if ( bn .and. ( bl1 .or. bl2 .or. bl3 ) ) then
          call wfat_chgubi ( a, b, x, hu, id )
          md = 3
        else
          call wfat_chguit ( a, b, x, hu, id )
          md = 4
        end if
      else
        if ( b .le. a ) then
          a00 = a
          b00 = b
          a = a - b + 1.0D+00
          b = 2.0D+00 - b
          call wfat_chguit ( a, b, x, hu, id )
          hu = x ** ( 1.0D+00 - b00 ) * hu
          a = a00
          b = b00
          md = 4
        else if ( bn .and. ( .not. il1 ) ) then
          call wfat_chgubi ( a, b, x, hu, id )
          md = 3
        end if
      end if

      if ( id .lt. 6 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CHGU - Warning!'
        write ( *, '(a)' ) '  Accurate results were not obtained.'
      end if

      return
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_chguit ( a, b, x, hu, id )

c*********************************************************************72
c
cc CHGUIT computes the hypergeometric function using Gauss-Legendre integration.
c
c  Discussion:
c
c    This procedure computes the hypergeometric function U(a,b,x) by
c    using Gaussian-Legendre integration (n = 60)
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
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
c    Output, double precision HU, U(a,b,z).
c
c    Output, integer ID, the estimated number of significant digits.
c
      implicit none

      double precision a
      double precision a1
      double precision b
      double precision b1
      double precision c
      double precision d
      double precision f1
      double precision f2
      double precision g
      double precision ga
      double precision hu
      double precision hu0
      double precision hu1
      double precision hu2
      integer id
      integer j
      integer k
      integer m
      double precision s
      double precision t(30)
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision w(30)
      double precision x

      save t
      save w

      data t /  0.259597723012478D-01, 0.778093339495366D-01,
     &          0.129449135396945D+00, 0.180739964873425D+00,
     &          0.231543551376029D+00, 0.281722937423262D+00,
     &          0.331142848268448D+00, 0.379670056576798D+00,
     &          0.427173741583078D+00, 0.473525841761707D+00,
     &          0.518601400058570D+00, 0.562278900753945D+00,
     &          0.604440597048510D+00, 0.644972828489477D+00,
     &          0.683766327381356D+00, 0.720716513355730D+00,
     &          0.755723775306586D+00, 0.788693739932264D+00,
     &          0.819537526162146D+00, 0.848171984785930D+00,
     &          0.874519922646898D+00, 0.898510310810046D+00,
     &          0.920078476177628D+00, 0.939166276116423D+00,
     &          0.955722255839996D+00, 0.969701788765053D+00,
     &          0.981067201752598D+00, 0.989787895222222D+00,
     &          0.995840525118838D+00, 0.999210123227436D+00 /

      data w /  0.519078776312206D-01, 0.517679431749102D-01,
     &          0.514884515009810D-01, 0.510701560698557D-01,
     &          0.505141845325094D-01, 0.498220356905502D-01,
     &          0.489955754557568D-01, 0.480370318199712D-01,
     &          0.469489888489122D-01, 0.457343797161145D-01,
     &          0.443964787957872D-01, 0.429388928359356D-01,
     &          0.413655512355848D-01, 0.396806954523808D-01,
     &          0.378888675692434D-01, 0.359948980510845D-01,
     &          0.340038927249464D-01, 0.319212190192963D-01,
     &          0.297524915007890D-01, 0.275035567499248D-01,
     &          0.251804776215213D-01, 0.227895169439978D-01,
     &          0.203371207294572D-01, 0.178299010142074D-01,
     &          0.152746185967848D-01, 0.126781664768159D-01,
     &          0.100475571822880D-01, 0.738993116334531D-02,
     &          0.471272992695363D-02, 0.202681196887362D-02 /

      id = 7
      a1 = a - 1.0D+00
      b1 = b - a - 1.0D+00
      c = 12.0D+00 / x

      do m = 10, 100, 5

        hu1 = 0.0D+00
        g = 0.5D+00 * c / m
        d = g
        do j = 1, m
          s = 0.0D+00
          do k = 1, 30
            t1 = d + g * t(k)
            t2 = d - g * t(k)
            f1 = exp ( - x * t1 ) * t1 ** a1 * ( 1.0D+00 + t1 ) ** b1
            f2 = exp ( - x * t2 ) * t2 ** a1 * ( 1.0D+00 + t2 ) ** b1
            s = s + w(k) * ( f1 + f2 )
          end do
          hu1 = hu1 + s * g
          d = d + 2.0D+00 * g
        end do

        if ( abs ( 1.0D+00 - hu0 / hu1 ) .lt. 1.0D-07 ) then
          go to 10
        end if

        hu0 = hu1

      end do

10    continue

      call wfat_gammaf ( a, ga )
      hu1 = hu1 / ga

      do m = 2, 10, 2
        hu2 = 0.0D+00
        g = 0.5D+00 / m
        d = g
        do j = 1, m
          s = 0.0D+00
          do k = 1, 30
            t1 = d + g * t(k)
            t2 = d - g * t(k)
            t3 = c / ( 1.0D+00 - t1 )
            t4 = c / ( 1.0D+00 - t2 ) 
            f1 = t3 * t3 / c * exp ( - x * t3 ) * t3 ** a1 
     &        * ( 1.0D+00 + t3 ) ** b1
            f2 = t4 * t4 / c * exp ( - x * t4 ) * t4 ** a1 
     &        * ( 1.0D+00 + t4 ) ** b1
            s = s + w(k) * ( f1 + f2 )
          end do
          hu2 = hu2 + s * g
          d = d + 2.0D+00 * g
        end do

        if ( abs ( 1.0D+00 - hu0 / hu2 ) .lt. 1.0D-07 ) then
          go to 20
        end if

        hu0 = hu2

      end do

20    continue

      call wfat_gammaf ( a, ga )
      hu2 = hu2 / ga
      hu = hu1 + hu2

      return
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_chgul ( a, b, x, hu, id )

c*********************************************************************72
c
cc CHGUL: confluent hypergeometric function U(a,b,x) for large argument X.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
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
      double precision aa
      double precision b
      double precision hu
      integer id
      logical il1
      logical il2
      integer k
      integer nm
      double precision r
      double precision ra
      double precision r0
      double precision x

      id = -100
      aa = a - b + 1.0D+00
      il1 = ( a .eq. int ( a ) ) .and. ( a .le. 0.0D+00 )
      il2 = ( aa .eq. int ( aa ) ) .and. ( aa .le. 0.0D+00 )

      if ( il1 .or. il2 ) then

        if ( il1 ) then
          nm = int ( abs ( a ) )
        end if

        if ( il2 ) then
          nm = int ( abs ( aa ) )
        end if

        hu = 1.0D+00
        r = 1.0D+00
        do k = 1, nm
          r = - r * ( a + k - 1.0D+00 ) * ( a - b + k ) / ( k * x )
          hu = hu + r
        end do
        hu = x ** ( - a ) * hu
        id = 10

      else

        hu = 1.0D+00
        r = 1.0D+00
        do k = 1, 25
          r = - r * ( a + k - 1.0D+00 ) * ( a - b + k ) / ( k * x )
          ra = abs ( r )
          if ( ( 5 .lt. k .and. r0 .le. ra ) .or. ra .lt. 1.0D-15 ) then
            go to 10
          end if
          r0 = ra
          hu = hu + r
        end do

10      continue

        id = int ( abs ( log10 ( ra ) ) )
        hu = x ** ( - a ) * hu

      end if

      return
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_chgus ( a, b, x, hu, id )
c*********************************************************************72
c
cc CHGUS: confluent hypergeometric function U(a,b,x) for small argument X.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    27 July 2012
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
c    Output, double precision HU, U(a,b,x).
c
c    Output, integer ID, the estimated number of significant digits.
c
      implicit none

      double precision a
      double precision b
      double precision d1
      double precision d2
      double precision ga
      double precision gab
      double precision gb
      double precision gb2
      double precision h0
      double precision hmax
      double precision hmin
      double precision hu
      double precision hu0
      double precision hua
      integer id
      integer j
      double precision pi
      double precision r1
      double precision r2
      double precision x
      double precision xg1
      double precision xg2

      id = -100
      pi = 3.141592653589793D+00
      call wfat_gammaf ( a, ga )
      call wfat_gammaf ( b, gb )
      xg1 = 1.0D+00 + a - b
      call wfat_gammaf ( xg1, gab )
      xg2 = 2.0D+00 - b
      call wfat_gammaf ( xg2, gb2 )
      hu0 = pi / sin ( pi * b )
      r1 = hu0 / ( gab * gb )
      r2 = hu0 * x ** ( 1.0D+00 - b ) / ( ga * gb2 )
      hu = r1 - r2
      hmax = 0.0D+00
      hmin = 1.0D+300
      do j = 1, 150
        r1 = r1 * ( a + j - 1.0D+00 ) / ( j * ( b + j - 1.0D+00 ) ) * x
        r2 = r2 * ( a - b + j ) / ( j * ( 1.0D+00 - b + j ) ) * x
        hu = hu + r1 - r2
        hua = abs ( hu )
        hmax = max ( hmax, hua )
        hmin = min ( hmin, hua )
        if ( abs ( hu - h0 ) .lt. abs ( hu ) * 1.0D-15 ) then
          go to 10
        end if
        h0 = hu
      end do

10    continue

      d1 = log10 ( hmax )
      if ( hmin .ne. 0.0D+00 ) then
        d2 = log10 ( hmin )
      end if
      id = 15 - int ( abs ( d1 - d2 ) )

      return
      end
c=======================================================

      
c=======================================================
c=======================================================
      subroutine wfat_chgubi ( a, b, x, hu, id )

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

      call wfat_psif ( a, ps )
      call wfat_gammaf ( a, ga )

      if ( 0.0D+00 .lt. b ) then
        a0 = a
        a1 = a - n
        a2 = a1
        call wfat_gammaf ( a1, ga1 )
        ua = ( - 1 ) ** ( n - 1 ) / ( rn * ga1 )
        ub = rn1 / ga * x ** ( - n )
      else
        a0 = a + n
        a1 = a0
        a2 = a
        call wfat_gammaf ( a1, ga1 )
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
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_chgubi_smallx(a, b, x, hu, kk)

      implicit none
#include "errquit.fh"      
      double precision a
      double precision b
      double precision x
      double precision hu
      integer kk

      integer n, k
      double precision ga, hu1, ln_x, gn, c, hu2, ps1, ps2, ps3, p, gk 
      double precision xx_k, factorial, rise_fact
      external factorial, rise_fact
      integer kmax
      parameter(kmax=200)


      if (abs(int(b)-b) > 1.0d-14)
     &     call errquit('WFAT: Argument b of wfat_chgubi_smallx must '//
     &     'be integer (but with double precision type).', 63, UERR)

      
      n = int(b) - 1
      call wfat_gammaf(a, ga)

      hu1 = 0.0d0
      if (b > 1) then
         do k = 1, n
            hu1 = hu1 + factorial(k-1)*rise_fact(1-a+k, n-k)
     &            /factorial(n-k) * x**(-k)
         enddo
         hu1 = hu1 / ga
      endif


      ln_x = log(x)
      call wfat_gammaf(a-n, gn)
      c = (-1)**(n+1) / ga / gn        ! c = the constant pre-factor for hu2 part
      hu2 = 0.0d0
      do k = 1, kmax
         call wfat_psif(a+k, ps1)
         call wfat_psif(dble(1+k), ps2)
         call wfat_psif(dble(1+k+n), ps3)
         p = ln_x + ps1 - ps2 - ps3
         call wfat_gammaf(a+k, gk)
         xx_k = c * gk / factorial(k) / factorial(n+k) * p * x**k
         hu2 = hu2 + xx_k
         
         if (abs(xx_k/(hu1+hu2)) < 1.0d-12) then
            kk = 100000 + k
            exit
         endif
      enddo
      hu = hu1 + hu2

      end
c=======================================================
      
