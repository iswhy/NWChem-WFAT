      SUBROUTINE ccsdt_y_tr1_7_1(d_a,k_a_offset,d_b,k_b_offset,d_c,k_c_o
     &ffset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i1 ( h2 h7 h9 p10 )_ytr + = 1 * Sum ( p3 ) * tr ( p3 h9 )_tr * y ( h2 h7 p3 p10 )_y
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
      INTEGER h2b
      INTEGER h7b
      INTEGER h9b
      INTEGER p10b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p3b
      INTEGER p3b_1
      INTEGER h9b_1
      INTEGER h2b_2
      INTEGER h7b_2
      INTEGER p10b_2
      INTEGER p3b_2
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
      DO h2b = 1,noab
      DO h7b = 1,noab
      DO h9b = 1,noab
      DO p10b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h2b-1)+int_mb(k_spin+h7b-1
     &)+int_mb(k_spin+h9b-1)+int_mb(k_spin+p10b-1).ne.8)) THEN
      IF (int_mb(k_spin+h2b-1)+int_mb(k_spin+h7b-1) .eq. int_mb(k_spin+h
     &9b-1)+int_mb(k_spin+p10b-1)) THEN
      IF (ieor(int_mb(k_sym+h2b-1),ieor(int_mb(k_sym+h7b-1),ieor(int_mb(
     &k_sym+h9b-1),int_mb(k_sym+p10b-1)))) .eq. ieor(irrep_y,irrep_tr)) 
     &THEN
      dimc = int_mb(k_range+h2b-1) * int_mb(k_range+h7b-1) * int_mb(k_ra
     &nge+h9b-1) * int_mb(k_range+p10b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('ccsdt_y_tr1_7_1',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      DO p3b = noab+1,noab+nvab
      IF (int_mb(k_spin+p3b-1) .eq. int_mb(k_spin+h9b-1)) THEN
      IF (ieor(int_mb(k_sym+p3b-1),int_mb(k_sym+h9b-1)) .eq. irrep_tr) T
     &HEN
      CALL TCE_RESTRICTED_2(p3b,h9b,p3b_1,h9b_1)
      CALL TCE_RESTRICTED_4(h2b,h7b,p10b,p3b,h2b_2,h7b_2,p10b_2,p3b_2)
      dim_common = int_mb(k_range+p3b-1)
      dima_sort = int_mb(k_range+h9b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+h2b-1) * int_mb(k_range+h7b-1) * int_mb
     &(k_range+p10b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsdt_y_tr1_7_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsdt_y_tr1_7_1',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h9b_1
     & - 1 + noab * (p3b_1 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p3b-1)
     &,int_mb(k_range+h9b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsdt_y_tr1_7_1',3,MA_ER
     &R)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsdt_y_tr1_7_1',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &ccsdt_y_tr1_7_1',5,MA_ERR)
      IF ((h7b .lt. h2b) .and. (p3b .le. p10b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p10b_
     &2 - noab - 1 + nvab * (p3b_2 - noab - 1 + nvab * (h2b_2 - 1 + noab
     & * (h7b_2 - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h7b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p3b-1),int_mb(k_range+p10b-1
     &),4,1,2,3,-1.0d0)
      END IF
      IF ((h7b .lt. h2b) .and. (p10b .lt. p3b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p3b_2
     & - noab - 1 + nvab * (p10b_2 - noab - 1 + nvab * (h2b_2 - 1 + noab
     & * (h7b_2 - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h7b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p10b-1),int_mb(k_range+p3b-1
     &),3,1,2,4,1.0d0)
      END IF
      IF ((h2b .le. h7b) .and. (p3b .le. p10b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p10b_
     &2 - noab - 1 + nvab * (p3b_2 - noab - 1 + nvab * (h7b_2 - 1 + noab
     & * (h2b_2 - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h2b-1)
     &,int_mb(k_range+h7b-1),int_mb(k_range+p3b-1),int_mb(k_range+p10b-1
     &),4,2,1,3,1.0d0)
      END IF
      IF ((h2b .le. h7b) .and. (p10b .lt. p3b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p3b_2
     & - noab - 1 + nvab * (p10b_2 - noab - 1 + nvab * (h7b_2 - 1 + noab
     & * (h2b_2 - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h2b-1)
     &,int_mb(k_range+h7b-1),int_mb(k_range+p10b-1),int_mb(k_range+p3b-1
     &),3,2,1,4,-1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('ccsdt_y_tr1_7_1',6,MA_ER
     &R)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('ccsdt_y_tr1_7_1',7,
     &MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('ccsdt_y_tr1_7_1',8,
     &MA_ERR)
      END IF
      END IF
      END IF
      END DO
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &ccsdt_y_tr1_7_1',9,MA_ERR)
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+p10b-1
     &),int_mb(k_range+h7b-1),int_mb(k_range+h2b-1),int_mb(k_range+h9b-1
     &),3,2,4,1,1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(p10b 
     &- noab - 1 + nvab * (h9b - 1 + noab * (h7b - 1 + noab * (h2b - 1))
     &)))
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('ccsdt_y_tr1_7_1',10,MA_E
     &RR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('ccsdt_y_tr1_7_1',11
     &,MA_ERR)
      END IF
      END IF
      END IF
      next = NXTVAL(nprocs)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      next = NXTVAL(-nprocs)
      call GA_SYNC()
      RETURN
      END
