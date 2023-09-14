!==================================================
!==================================================
      subroutine wfat_3j_jscan(j1, m1, j2, m2, j3, m3, cf)

      implicit none
      double precision tmp1(j3-abs(m3)+1)
      double precision tmp2(j2-abs(m2)+1)
      double precision cf(j1-abs(m1)+1,j2-abs(m2)+1,j3-abs(m3)+1)


      do j = 1, j3-abs(m3)+1
         do i = 1, j2-abs(m2)+1

            call wfat_3j_jscan0(tmp2(i,j), , , , cf(:,i,j))

         enddo
      enddo


      end
!==================================================


!==================================================
!==================================================
      subroutine wfat_3j_jscan0(tj1, tm1, tj2, tm2, cf)
c     / j1  j2      j3   \
c     \ m1  m2  -(m1+m2) /
c     for j3 = [j_min, j1+j2], where
c     j_min = max(|j1-j2|,|m1+m2|)
      implicit none
      integer tj1
      integer tm1
      integer tj2
      integer tm2
      double precision cf(1:*)

      double precision j1, j2, j3, m1, m2, j_min

      
      j1 = tj1 / 2.0d0
      j2 = tj2 / 2.0d0
      m1 = tm1 / 2.0d0
      m2 = tm2 / 2.0d0

      
c     ==== Determine the type of the first 3-j symbol ====
      j_min = max(abs(j1-j2), abs(m1+m2))
      if ( idnint(2*j_min) == abs(tj1-tj2) ) then
c        / j1  j2   |j1-j2| \
c        \ m1  m2  -(m1+m2) /
         call wfat_3j_special_type1(tj1, tm1, tj2, tm2, cf(1))
c        / j1  j2  |j1-j2|+1 \
c        \ m1  m2   -(m1+m2) /
         call wfat_3j_special_type3(tj1, tm1, tj2, tm2, cf(2))
      else if ( idnint(2*j_min) == abs(tm1+tm2) ) then
c        / j1  j2   |m1+m2| \
c        \ m1  m2  -(m1+m2) /
         call wfat_3j_special_type2(tj1, tm1, tj2, tm2, -(tm1+tm2),
     &        cf(1))
c        / j1  j2   |m1+m2|+1 \
c        \ m1  m2    -(m1+m2) /
         call wfat_3j_special_type4(tj1, tm1, tj2, tm2, -(tm1+tm2),
     &        cf(2))
      endif
      

c     ==== Begin the recurrence ====
      do i = 3, idnint((j1+j2) - j_min + 1)
         j3 = j_min + i - 1

         call wfat_3j_a(j3, , a1)
         call wfat_3j_a(j3-1, , a2)
         call wfat_3j_b(j3-1, , b)
         cf(i) = (-b*cf(i-1) - j3*a2*cf(i-2)) / (j3-1) / a1
      enddo

       
      end
!==================================================


!==================================================
!==================================================
      subroutine wfat_3j_special_type1(tj1, tm1, tj2, tm2, cf)
c     Calculate the 3-j symbol:
c     / j1  j2   |j1-j2| \
c     \ m1  m2  -(m1+m2) /
      implicit none 
      integer tj1    ! [input], 2*j1
      integer tm1    ! [input], 2*m1
      integer tj2    ! [input], 2*j2
      integer tm2    ! [input], 2*m2
      double precision cf     ! [output], the 3-j symbol

      integer i, n, j
      double precision f, j1, m1, j2, m2, j3, v(9), j1_, m1_, j2_, m2_


      if (abs(tj1-tj2) < abs(tm1+tm2)) then
         cf = 0.0d0
         return
      endif
      

c     Calculate
c        <(j1,-m1)(j2,-m2)|(|j1-j2|,-m1-m2)>
c     if j1>=j2, and
c        <(j2,-m2)(j1,-m1)|(|j1-j2|,-m1-m2)>
c     if j1<j2.
      if (tj1 >= tj2) then
         j1 = tj1 / 2.0d0
         m1 = tm1 / 2.0d0
         j2 = tj2 / 2.0d0
         m2 = tm2 / 2.0d0
      else
         j1 = tj2 / 2.0d0
         m1 = tm2 / 2.0d0
         j2 = tj1 / 2.0d0
         m2 = tm1 / 2.0d0
      endif
      j3 = abs(j1-j2)
      
      v(1) = j1 + m1
      v(2) = j1 - m1
      v(3) = 2*j2
      v(4) = 2*j1 - 2*j2 + 1
      v(5) = 2*j1 + 1
      v(6) = j2 + m2
      v(7) = j2 - m2
      v(8) = j1 - j2 + m1 + m2
      v(9) = j1 - j2 - m1 - m2

      
c     ==== Sec. 8.5, eq. (13) in Ref. [1] ====
      n = idnint(maxval(v))
      f = 1.0d0
      do i = 1, n
         f = f * v(1) * v(2) * v(3) * v(4) / v(5) / v(6) / v(7) / v(8)
     &        / v(9)          
         do j = 1, 9
            v(j) = v(j) - 1
            if (idnint(v(j)) < 1) v(j) = 1.0d0
         enddo
      enddo
      cf = (-1)**idnint(j2+m2) * sqrt(f)


c     ==== Sec. 8.1, eq. (11) in Ref. [1] ====
      j1_ = tj1 / 2.0d0
      m1_ = tm1 / 2.0d0
      j2_ = tj2 / 2.0d0
      m2_ = tm2 / 2.0d0
      if (tj1 >= tj2) then
         cf = (-1)**idnint(j3-m1_-m2_+tj1) / sqrt(2*j3+1) * cf
      else
         cf = (-1)**idnint(j1_+3*j2_+2*j3-m1_-m2_) / sqrt(2*j3+1) * cf
      endif
      
      end
c     REFERENCES:
c     [1] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum
c         Theory of Angular Momentum (World Scientific, Singapore, 1988).
!==================================================

      
!==================================================
!==================================================
      subroutine wfat_3j_special_type2(tj1, tm1, tj2, tm2, tm3, cf)
c     Calculate the 3-j symbol:
c     / j1  j2  |m3| \  or  / j1  j2  |m3| \
c     \ m1  m2   m3  /      \ m1  m2   -m3 /
      implicit none 
      integer tj1    ! [input], 2*j1
      integer tm1    ! [input], 2*m1
      integer tj2    ! [input], 2*j2
      integer tm2    ! [input], 2*m2
      integer tm3    ! [input], 2*m3
      double precision cf     ! [output], the 3-j symbol

      integer i, n, j
      double precision f, j1, m1, j2, m2, j3, v(9)


      if (tm1+tm2+tm3 /= 0) then
         cf = 0.0d0
         return
      endif


c     Calculate
c        <(j1,-m1)(j2,-m2)|(|m3|,m3)>
c     if m3>=0, and
c        <(j1,m1)(j2,m2)|(|m3|,-m3)>
c     if m3<0.
      j1 = tj1 / 2.0d0
      if (idnint(tm3) >= 0) then
         m1 = -tm1 / 2.0d0
      else
         m1 = tm1 / 2.0d0
      endif
      j2 = tj2 / 2.0d0
      if (idnint(tm3) >= 0) then
         m2 = -tm2 / 2.0d0
      else
         m2 = tm2 / 2.0d0
      endif
      j3 = abs(tm3) / 2.0d0

      
      v(1) = 2*j3 + 1
      v(2) = j1 + j2 - j3
      v(3) = j1 + m1
      v(4) = j2 + m2
      v(5) = j1 + j2 + j3 + 1
      v(6) = j1 - j2 + j3
      v(7) = -j1 + j2 + j3
      v(8) = j1 - m1
      v(9) = j2 - m2

      
c     ==== Sec. 8.5, eq. (35) in Ref. [1] ====
      n = idnint(maxval(v))
      f = 1.0d0
      do i = 1, n
         f = f * v(1) * v(2) * v(3) * v(4) / v(5) / v(6) / v(7) / v(8)
     &        / v(9)          
         do j = 1, 9
            v(j) = v(j) - 1
            if (idnint(v(j)) < 1) v(j) = 1.0d0
         enddo
      enddo
      cf = (-1)**idnint(j1-m1) * sqrt(f)

      
c     ==== Sec. 8.1, eq. (11) in Ref. [1] ====
      if (idnint(tm3) >= 0) then
         cf = (-1)**idnint(2*j3+tj1) / sqrt(2*j3+1) * cf
      else
         cf = (-1)**idnint(3*j1+j2+tm3/2.0d0) / sqrt(2*j3+1) * cf
      endif
      
      end
c     REFERENCES:
c     [1] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum
c         Theory of Angular Momentum (World Scientific, Singapore, 1988).
!==================================================


!==================================================
!==================================================
      subroutine wfat_3j_special_type3(tj1, tm1, tj2, tm2, cf)
c     / j1  j2  |j1-j2|+1 \
c     \ m1  m2   -(m1+m2) /
      implicit none
      integer tj1    ! [input], 2*j1
      integer tm1    ! [input], 2*m1
      integer tj2    ! [input], 2*j2
      integer tm2    ! [input], 2*m2
      double precision cf     ! [output], the 3-j symbol

      integer i, n, j
      double precision f, j1, m1, j2, m2, j3, v0, v(9), j1_, m1_, j2_
      double precision m2_
      

      if ((tm1 == 0 .and. tm2 == 0) .or. 
     &     (abs(tj1-tj2)+2 < abs(tm1+tm2))) then
         cf = 0.0d0
         return
      endif
      

c     Calculate
c        <(j1,-m1)(j2,-m2)|(|j1-j2|+1,-m1-m2)>
c     if j1>=j2, and
c        <(j2,-m2)(j1,-m1)|(|j1-j2|+1,-m1-m2)>
c     if j1<j2.
      if (tj1 >= tj2) then
         j1 = tj1 / 2.0d0
         m1 = tm1 / 2.0d0
         j2 = tj2 / 2.0d0
         m2 = tm2 / 2.0d0
      else
         j1 = tj2 / 2.0d0
         m1 = tm2 / 2.0d0
         j2 = tj1 / 2.0d0
         m2 = tm1 / 2.0d0
      endif
      j3 = abs(j1-j2)+1

      v0 = 2*j1 - 2*j2 + 3
      v(1) = 2*j2 - 1             
      v(2) = 2*j1 - 2*j2 + 1      
      v(3) = j1 + m1              
      v(4) = j1 - m1              
      v(5) = 2*j1 + 2             
      v(6) = j2 + m2              
      v(7) = j2 - m2              
      v(8) = j1 - j2 + m1 + m2 + 1
      v(9) = j1 - j2 - m1 - m2 + 1

      
c     ==== Sec. 8.5, eq. (17) in Ref. [1] ====
      n = idnint(maxval(v))
      f = 1.0d0
      do i = 1, n
         f = f * v(1) * v(2) * v(3) * v(4) / v(5) / v(6) / v(7) / v(8) 
     &        / v(9)
         do j = 1, 9
            v(j) = v(j) - 1
            if (idnint(v(j)) < 1) v(j) = 1.0d0
         enddo
      enddo
      cf = (-1)**idnint(j2+m2+1) * 2*(j1*m2+j2*m1+m2) * sqrt(v0*f)


c     ==== Sec. 8.1, eq. (11) in Ref. [1] ====
      j1_ = tj1 / 2.0d0
      m1_ = tm1 / 2.0d0
      j2_ = tj2 / 2.0d0
      m2_ = tm2 / 2.0d0
      if (tj1 >= tj2) then
         cf = (-1)**idnint(j3-m1_-m2_+tj1) / sqrt(2*j3+1) * cf
      else
         cf = (-1)**idnint(j1_+3*j2_+2*j3-m1_-m2_) / sqrt(2*j3+1) * cf
      endif
      
      end
c     REFERENCES:
c     [1] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum
c         Theory of Angular Momentum (World Scientific, Singapore, 1988).
!==================================================
      
      
!==================================================
!==================================================
      subroutine wfat_3j_special_type4(tj1, tm1, tj2, tm2, tm3, cf)
c     Calculate the 3-j symbol:
c     / j1  j2  |m3|+1 \  or  / j1  j2  |m3|+1 \
c     \ m1  m2    m3   /      \ m1  m2    -m3  /
      implicit none
      integer tj1    ! [input], 2*j1
      integer tm1    ! [input], 2*m1
      integer tj2    ! [input], 2*j2
      integer tm2    ! [input], 2*m2
      integer tm3    ! [input], 2*m3
      double precision cf     ! [output], the 3-j symbol

      integer i, n, j
      double precision f, j1, m1, j2, m2, j3, c, v0, v(9)
      

      if (tm1+tm2+tm3 /= 0) then
         cf = 0.0d0
         return
      endif


c     Calculate
c        <(j1,-m1)(j2,-m2)|(|m3|+1,m3)>
c     if m3>=0, and
c        <(j1,m1)(j2,m2)|(|m3|+1,-m3)>
c     if m3<0.
      j1 = tj1 / 2.0d0
      if (idnint(tm3) >= 0) then
         m1 = -tm1 / 2.0d0
      else
         m1 = tm1 / 2.0d0
      endif
      j2 = tj2 / 2.0d0
      if (idnint(tm3) >= 0) then
         m2 = -tm2 / 2.0d0
      else
         m2 = tm2 / 2.0d0
      endif
      j3 = abs(tm3) / 2.0d0 + 1
      c = j3 - 1

      v0 = 2*c + 3
      v(1) = 2*c + 1
      v(2) = j1 + j2 - c - 1
      v(3) = j1 + m1
      v(4) = j2 + m2
      v(5) = j1 + j2 + c + 2
      v(6) = j1 - j2 + c + 1
      v(7) = -j1 + j2 + c + 1
      v(8) = j1 - m1
      v(9) = j2 - m2

      
c     ==== Sec. 8.5, eq. (44) in Ref. [1] ====
      n = idnint(maxval(v))
      f = 1.0d0
      do i = 1, n
         f = f * v(1) * v(2) * v(3) * v(4) / v(5) / v(6) / v(7) / v(8)
     &        / v(9)          
         do j = 1, 9
            v(j) = v(j) - 1
            if (idnint(v(j)) < 1) v(j) = 1.0d0
         enddo
      enddo
      cf = (-1)**idnint(j1-m1) * ((j2-m2)*(j2+m2+1) - (j1-m1)*(j1+m1+1))
     &     * sqrt(v0*f)


c     ==== Sec. 8.1, eq. (11) in Ref. [1] ====
      if (idnint(tm3) >= 0) then
         cf = (-1)**idnint(tm3+1+tj1) / sqrt(2*j3+1) * cf
      else
         cf = (-1)**idnint(3*j1+j2+tm3/2.0d0) / sqrt(2*j3+1) * cf
      endif
      
      end
c     REFERENCES:
c     [1] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum
c         Theory of Angular Momentum (World Scientific, Singapore, 1988).
!==================================================
