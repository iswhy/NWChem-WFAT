c done till here 2/7/2014 11:30am
      SUBROUTINE beta_1_3_7_2_1(d_a,k_a_offset,d_b,k_b_offset,d_c,k_c_of
     &fset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i3 ( h4 p7 )_vtrb + = 1 * Sum ( h6 p5 ) * trb ( p5 h6 )_trb * v ( h4 h6 p5 p7 )_v
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER d_c
      INTEGER k_c_offset
      INTEGER NXTVAL
      INTEGER next
      INTEGER nprocs
      INTEGER count
      INTEGER h4b
      INTEGER p7b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p5b
      INTEGER h6b
      INTEGER p5b_1
      INTEGER h6b_1
      INTEGER h4b_2
      INTEGER h6b_2
      INTEGER p7b_2
      INTEGER p5b_2
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER dimb_sort
      INTEGER dimb
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_b_sort
      INTEGER k_b_sort
      INTEGER l_b
      INTEGER k_b
      INTEGER l_c
      INTEGER k_c
      EXTERNAL NXTVAL
      nprocs = GA_NNODES()
      count = 0
      next = NXTVAL(nprocs)
      DO h4b = 1,noab
      DO p7b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h4b-1)+int_mb(k_spin+p7b-1
     &).ne.4)) THEN
      IF (int_mb(k_spin+h4b-1) .eq. int_mb(k_spin+p7b-1)) THEN
      IF (ieor(int_mb(k_sym+h4b-1),int_mb(k_sym+p7b-1)) .eq. ieor(irrep_
     &v,irrep_trb)) THEN
      dimc = int_mb(k_range+h4b-1) * int_mb(k_range+p7b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('beta_1_3_7_2_1',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      DO p5b = noab+1,noab+nvab
      DO h6b = 1,noab
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h6b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h6b-1)) .eq. irrep_trb) 
     &THEN
      CALL TCE_RESTRICTED_2(p5b,h6b,p5b_1,h6b_1)
      CALL TCE_RESTRICTED_4(h4b,h6b,p7b,p5b,h4b_2,h6b_2,p7b_2,p5b_2)
      dim_common = int_mb(k_range+p5b-1) * int_mb(k_range+h6b-1)
      dima_sort = 1
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+h4b-1) * int_mb(k_range+p7b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('beta_1_3_7_2_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &beta_1_3_7_2_1',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h6b_1
     & - 1 + noab * (p5b_1 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+h6b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('beta_1_3_7_2_1',3,MA_ERR
     &)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('beta_1_3_7_2_1',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &beta_1_3_7_2_1',5,MA_ERR)
      IF ((h6b .lt. h4b) .and. (p5b .le. p7b)) THEN
      if(.not.intorb) then
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p7b_2
     & - 1 + (noab+nvab) * (p5b_2 - 1 + (noab+nvab) * (h4b_2 - 1 + (noab
     &+nvab) * (h6b_2 - 1)))))
      else
      CALL GET_HASH_BLOCK_I(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),
     &(p7b_2
     & - 1 + (noab+nvab) * (p5b_2 - 1 + (noab+nvab) * (h4b_2 - 1 + (noab
     &+nvab) * (h6b_2 - 1)))),p7b_2,p5b_2,h4b_2,h6b_2)
      end if
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h6b-1)
     &,int_mb(k_range+h4b-1),int_mb(k_range+p5b-1),int_mb(k_range+p7b-1)
     &,4,2,1,3,-1.0d0)
      END IF
      IF ((h6b .lt. h4b) .and. (p7b .lt. p5b)) THEN
      if(.not.intorb) then
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p5b_2
     & - 1 + (noab+nvab) * (p7b_2 - 1 + (noab+nvab) * (h4b_2 - 1 + (noab
     &+nvab) * (h6b_2 - 1)))))
      else
      CALL GET_HASH_BLOCK_I(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),
     &(p5b_2
     & - 1 + (noab+nvab) * (p7b_2 - 1 + (noab+nvab) * (h4b_2 - 1 + (noab
     &+nvab) * (h6b_2 - 1)))),p5b_2,p7b_2,h4b_2,h6b_2)
      end if
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h6b-1)
     &,int_mb(k_range+h4b-1),int_mb(k_range+p7b-1),int_mb(k_range+p5b-1)
     &,3,2,1,4,1.0d0)
      END IF
      IF ((h4b .le. h6b) .and. (p5b .le. p7b)) THEN
      if(.not.intorb) then
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p7b_2
     & - 1 + (noab+nvab) * (p5b_2 - 1 + (noab+nvab) * (h6b_2 - 1 + (noab
     &+nvab) * (h4b_2 - 1)))))
      else
      CALL GET_HASH_BLOCK_I(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),
     &(p7b_2
     & - 1 + (noab+nvab) * (p5b_2 - 1 + (noab+nvab) * (h6b_2 - 1 + (noab
     &+nvab) * (h4b_2 - 1)))),p7b_2,p5b_2,h6b_2,h4b_2)
      end if
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h4b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+p5b-1),int_mb(k_range+p7b-1)
     &,4,1,2,3,1.0d0)
      END IF
      IF ((h4b .le. h6b) .and. (p7b .lt. p5b)) THEN
      if(.not.intorb) then
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p5b_2
     & - 1 + (noab+nvab) * (p7b_2 - 1 + (noab+nvab) * (h6b_2 - 1 + (noab
     &+nvab) * (h4b_2 - 1)))))
      else
      CALL GET_HASH_BLOCK_I(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),
     &(p5b_2
     & - 1 + (noab+nvab) * (p7b_2 - 1 + (noab+nvab) * (h6b_2 - 1 + (noab
     &+nvab) * (h4b_2 - 1)))),p5b_2,p7b_2,h6b_2,h4b_2)
      end if
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h4b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+p7b-1),int_mb(k_range+p5b-1)
     &,3,1,2,4,-1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('beta_1_3_7_2_1',6,MA_ERR
     &)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('beta_1_3_7_2_1',7,M
     &A_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('beta_1_3_7_2_1',8,M
     &A_ERR)
      END IF
      END IF
      END IF
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &beta_1_3_7_2_1',9,MA_ERR)
      CALL TCE_SORT_2(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+p7b-1)
     &,int_mb(k_range+h4b-1),2,1,1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(p7b -
     & noab - 1 + nvab * (h4b - 1)))
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('beta_1_3_7_2_1',10,MA_ER
     &R)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('beta_1_3_7_2_1',11,
     &MA_ERR)
      END IF
      END IF
      END IF
      next = NXTVAL(nprocs)
      END IF
      count = count + 1
      END DO
      END DO
      next = NXTVAL(-nprocs)
      call GA_SYNC()
      RETURN
      END
