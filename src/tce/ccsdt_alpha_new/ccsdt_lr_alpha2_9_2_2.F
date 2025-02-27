      SUBROUTINE ccsdt_lr_alpha2_9_2_2(d_a,k_a_offset,d_b,k_b_offset,d_c
     &,k_c_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i2 ( h7 h8 h14 h15 )_ytrb + = 1/3 * Sum ( h4 p1 p2 p3 ) * trb ( p1 p2 p3 h4 h14 h15 )_trb * y ( h4 h7 h8 p1 p2 p3 )_y
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
      INTEGER h7b
      INTEGER h8b
      INTEGER h14b
      INTEGER h15b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p1b
      INTEGER p2b
      INTEGER p3b
      INTEGER h4b
      INTEGER p1b_1
      INTEGER p2b_1
      INTEGER p3b_1
      INTEGER h14b_1
      INTEGER h15b_1
      INTEGER h4b_1
      INTEGER h7b_2
      INTEGER h8b_2
      INTEGER h4b_2
      INTEGER p1b_2
      INTEGER p2b_2
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
      INTEGER nsuperp(3)
      INTEGER isuperp
      INTEGER l_c
      INTEGER k_c
      DOUBLE PRECISION FACTORIAL
      EXTERNAL NXTVAL
      EXTERNAL FACTORIAL
      nprocs = GA_NNODES()
      count = 0
      next = NXTVAL(nprocs)
      DO h7b = 1,noab
      DO h8b = h7b,noab
      DO h14b = 1,noab
      DO h15b = h14b,noab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1
     &)+int_mb(k_spin+h14b-1)+int_mb(k_spin+h15b-1).ne.8)) THEN
      IF (int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1) .eq. int_mb(k_spin+h
     &14b-1)+int_mb(k_spin+h15b-1)) THEN
      IF (ieor(int_mb(k_sym+h7b-1),ieor(int_mb(k_sym+h8b-1),ieor(int_mb(
     &k_sym+h14b-1),int_mb(k_sym+h15b-1)))) .eq. ieor(irrep_y,irrep_trb)
     &) THEN
      dimc = int_mb(k_range+h7b-1) * int_mb(k_range+h8b-1) * int_mb(k_ra
     &nge+h14b-1) * int_mb(k_range+h15b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('ccsdt_lr_alpha2_9_2_2',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
      DO h4b = 1,noab
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)
     & .eq. int_mb(k_spin+h14b-1)+int_mb(k_spin+h15b-1)+int_mb(k_spin+h4
     &b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
     &k_sym+p3b-1),ieor(int_mb(k_sym+h14b-1),ieor(int_mb(k_sym+h15b-1),i
     &nt_mb(k_sym+h4b-1)))))) .eq. irrep_trb) THEN
      CALL TCE_RESTRICTED_6(p1b,p2b,p3b,h14b,h15b,h4b,p1b_1,p2b_1,p3b_1,
     &h14b_1,h15b_1,h4b_1)
      CALL TCE_RESTRICTED_6(h7b,h8b,h4b,p1b,p2b,p3b,h7b_2,h8b_2,h4b_2,p1
     &b_2,p2b_2,p3b_2)
      dim_common = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) * int_m
     &b(k_range+p3b-1) * int_mb(k_range+h4b-1)
      dima_sort = int_mb(k_range+h14b-1) * int_mb(k_range+h15b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+h7b-1) * int_mb(k_range+h8b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsdt_lr_alpha2_9_2_2',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsdt_lr_alpha2_9_2_2',2,MA_ERR)
      IF ((h4b .le. h14b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h15b_
     &1 - 1 + noab * (h14b_1 - 1 + noab * (h4b_1 - 1 + noab * (p3b_1 - n
     &oab - 1 + nvab * (p2b_1 - noab - 1 + nvab * (p1b_1 - noab - 1)))))
     &))
      CALL TCE_SORT_6(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p1b-1)
     &,int_mb(k_range+p2b-1),int_mb(k_range+p3b-1),int_mb(k_range+h4b-1)
     &,int_mb(k_range+h14b-1),int_mb(k_range+h15b-1),6,5,4,3,2,1,1.0d0)
      END IF
      IF ((h14b .lt. h4b) .and. (h4b .le. h15b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h15b_
     &1 - 1 + noab * (h4b_1 - 1 + noab * (h14b_1 - 1 + noab * (p3b_1 - n
     &oab - 1 + nvab * (p2b_1 - noab - 1 + nvab * (p1b_1 - noab - 1)))))
     &))
      CALL TCE_SORT_6(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p1b-1)
     &,int_mb(k_range+p2b-1),int_mb(k_range+p3b-1),int_mb(k_range+h14b-1
     &),int_mb(k_range+h4b-1),int_mb(k_range+h15b-1),6,4,5,3,2,1,-1.0d0)
      END IF
      IF ((h15b .lt. h4b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h4b_1
     & - 1 + noab * (h15b_1 - 1 + noab * (h14b_1 - 1 + noab * (p3b_1 - n
     &oab - 1 + nvab * (p2b_1 - noab - 1 + nvab * (p1b_1 - noab - 1)))))
     &))
      CALL TCE_SORT_6(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p1b-1)
     &,int_mb(k_range+p2b-1),int_mb(k_range+p3b-1),int_mb(k_range+h14b-1
     &),int_mb(k_range+h15b-1),int_mb(k_range+h4b-1),5,4,6,3,2,1,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsdt_lr_alpha2_9_2_2',3
     &,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsdt_lr_alpha2_9_2_2',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &ccsdt_lr_alpha2_9_2_2',5,MA_ERR)
      IF ((h4b .le. h7b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p3b_2
     & - noab - 1 + nvab * (p2b_2 - noab - 1 + nvab * (p1b_2 - noab - 1 
     &+ nvab * (h8b_2 - 1 + noab * (h7b_2 - 1 + noab * (h4b_2 - 1)))))))
      CALL TCE_SORT_6(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h4b-1)
     &,int_mb(k_range+h7b-1),int_mb(k_range+h8b-1),int_mb(k_range+p1b-1)
     &,int_mb(k_range+p2b-1),int_mb(k_range+p3b-1),3,2,1,6,5,4,1.0d0)
      END IF
      IF ((h7b .lt. h4b) .and. (h4b .le. h8b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p3b_2
     & - noab - 1 + nvab * (p2b_2 - noab - 1 + nvab * (p1b_2 - noab - 1 
     &+ nvab * (h8b_2 - 1 + noab * (h4b_2 - 1 + noab * (h7b_2 - 1)))))))
      CALL TCE_SORT_6(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h7b-1)
     &,int_mb(k_range+h4b-1),int_mb(k_range+h8b-1),int_mb(k_range+p1b-1)
     &,int_mb(k_range+p2b-1),int_mb(k_range+p3b-1),3,1,2,6,5,4,-1.0d0)
      END IF
      IF ((h8b .lt. h4b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p3b_2
     & - noab - 1 + nvab * (p2b_2 - noab - 1 + nvab * (p1b_2 - noab - 1 
     &+ nvab * (h4b_2 - 1 + noab * (h8b_2 - 1 + noab * (h7b_2 - 1)))))))
      CALL TCE_SORT_6(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h7b-1)
     &,int_mb(k_range+h8b-1),int_mb(k_range+h4b-1),int_mb(k_range+p1b-1)
     &,int_mb(k_range+p2b-1),int_mb(k_range+p3b-1),2,1,3,6,5,4,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('ccsdt_lr_alpha2_9_2_2',6
     &,MA_ERR)
      nsuperp(1) = 1
      nsuperp(2) = 1
      nsuperp(3) = 1
      isuperp = 1
      IF (p1b .eq. p2b) THEN
      nsuperp(isuperp) = nsuperp(isuperp) + 1
      ELSE
      isuperp = isuperp + 1
      END IF
      IF (p2b .eq. p3b) THEN
      nsuperp(isuperp) = nsuperp(isuperp) + 1
      ELSE
      isuperp = isuperp + 1
      END IF
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,6.0d0/FACTORIAL(
     &nsuperp(1))/FACTORIAL(nsuperp(2))/FACTORIAL(nsuperp(3)),dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('ccsdt_lr_alpha2_9_2
     &_2',7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('ccsdt_lr_alpha2_9_2
     &_2',8,MA_ERR)
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &ccsdt_lr_alpha2_9_2_2',9,MA_ERR)
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h8b-1)
     &,int_mb(k_range+h7b-1),int_mb(k_range+h15b-1),int_mb(k_range+h14b-
     &1),2,1,4,3,1.0d0/3.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h15b 
     &- 1 + noab * (h14b - 1 + noab * (h8b - 1 + noab * (h7b - 1)))))
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('ccsdt_lr_alpha2_9_2_2',1
     &0,MA_ERR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('ccsdt_lr_alpha2_9_2
     &_2',11,MA_ERR)
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
