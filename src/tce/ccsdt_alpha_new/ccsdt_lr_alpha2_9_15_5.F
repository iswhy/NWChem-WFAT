      SUBROUTINE ccsdt_lr_alpha2_9_15_5(d_a,k_a_offset,d_b,k_b_offset,d_
     &c,k_c_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i2 ( h2 p12 h14 h15 )_yttrb + = 8 * Sum ( h7 p5 ) * trb ( p5 p12 h7 h15 )_trb * i3 ( h2 h7 h14 p5 )_yt
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
      INTEGER p12b
      INTEGER h2b
      INTEGER h15b
      INTEGER h14b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p5b
      INTEGER h7b
      INTEGER p12b_1
      INTEGER p5b_1
      INTEGER h15b_1
      INTEGER h7b_1
      INTEGER h2b_2
      INTEGER h7b_2
      INTEGER h14b_2
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
      DO p12b = noab+1,noab+nvab
      DO h2b = 1,noab
      DO h15b = 1,noab
      DO h14b = 1,noab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h2b-1)+int_mb(k_spin+p12b-
     &1)+int_mb(k_spin+h14b-1)+int_mb(k_spin+h15b-1).ne.8)) THEN
      IF (int_mb(k_spin+h2b-1)+int_mb(k_spin+p12b-1) .eq. int_mb(k_spin+
     &h14b-1)+int_mb(k_spin+h15b-1)) THEN
      IF (ieor(int_mb(k_sym+h2b-1),ieor(int_mb(k_sym+p12b-1),ieor(int_mb
     &(k_sym+h14b-1),int_mb(k_sym+h15b-1)))) .eq. ieor(irrep_y,ieor(irre
     &p_t,irrep_trb))) THEN
      dimc = int_mb(k_range+h2b-1) * int_mb(k_range+p12b-1) * int_mb(k_r
     &ange+h14b-1) * int_mb(k_range+h15b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('ccsdt_lr_alpha2_9_15_5',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      DO p5b = noab+1,noab+nvab
      DO h7b = 1,noab
      IF (int_mb(k_spin+p12b-1)+int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+
     &h15b-1)+int_mb(k_spin+h7b-1)) THEN
      IF (ieor(int_mb(k_sym+p12b-1),ieor(int_mb(k_sym+p5b-1),ieor(int_mb
     &(k_sym+h15b-1),int_mb(k_sym+h7b-1)))) .eq. irrep_trb) THEN
      CALL TCE_RESTRICTED_4(p12b,p5b,h15b,h7b,p12b_1,p5b_1,h15b_1,h7b_1)
      CALL TCE_RESTRICTED_4(h2b,h7b,h14b,p5b,h2b_2,h7b_2,h14b_2,p5b_2)
      dim_common = int_mb(k_range+p5b-1) * int_mb(k_range+h7b-1)
      dima_sort = int_mb(k_range+p12b-1) * int_mb(k_range+h15b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+h2b-1) * int_mb(k_range+h14b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsdt_lr_alpha2_9_15_5',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsdt_lr_alpha2_9_15_5',2,MA_ERR)
      IF ((p5b .le. p12b) .and. (h7b .le. h15b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h15b_
     &1 - 1 + noab * (h7b_1 - 1 + noab * (p12b_1 - noab - 1 + nvab * (p5
     &b_1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+p12b-1),int_mb(k_range+h7b-1),int_mb(k_range+h15b-
     &1),4,2,3,1,1.0d0)
      END IF
      IF ((p5b .le. p12b) .and. (h15b .lt. h7b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h7b_1
     & - 1 + noab * (h15b_1 - 1 + noab * (p12b_1 - noab - 1 + nvab * (p5
     &b_1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+p12b-1),int_mb(k_range+h15b-1),int_mb(k_range+h7b-
     &1),3,2,4,1,-1.0d0)
      END IF
      IF ((p12b .lt. p5b) .and. (h7b .le. h15b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h15b_
     &1 - 1 + noab * (h7b_1 - 1 + noab * (p5b_1 - noab - 1 + nvab * (p12
     &b_1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p12b-1
     &),int_mb(k_range+p5b-1),int_mb(k_range+h7b-1),int_mb(k_range+h15b-
     &1),4,1,3,2,-1.0d0)
      END IF
      IF ((p12b .lt. p5b) .and. (h15b .lt. h7b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h7b_1
     & - 1 + noab * (h15b_1 - 1 + noab * (p5b_1 - noab - 1 + nvab * (p12
     &b_1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p12b-1
     &),int_mb(k_range+p5b-1),int_mb(k_range+h15b-1),int_mb(k_range+h7b-
     &1),3,1,4,2,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsdt_lr_alpha2_9_15_5',
     &3,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsdt_lr_alpha2_9_15_5',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &ccsdt_lr_alpha2_9_15_5',5,MA_ERR)
      IF ((h7b .lt. h2b) .and. (h14b .le. p5b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p5b_2
     & - noab - 1 + nvab * (h14b_2 - 1 + noab * (h2b_2 - 1 + noab * (h7b
     &_2 - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h7b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+h14b-1),int_mb(k_range+p5b-1
     &),3,2,1,4,-1.0d0)
      END IF
      IF ((h2b .le. h7b) .and. (h14b .le. p5b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p5b_2
     & - noab - 1 + nvab * (h14b_2 - 1 + noab * (h7b_2 - 1 + noab * (h2b
     &_2 - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h2b-1)
     &,int_mb(k_range+h7b-1),int_mb(k_range+h14b-1),int_mb(k_range+p5b-1
     &),3,1,2,4,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('ccsdt_lr_alpha2_9_15_5',
     &6,MA_ERR)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('ccsdt_lr_alpha2_9_1
     &5_5',7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('ccsdt_lr_alpha2_9_1
     &5_5',8,MA_ERR)
      END IF
      END IF
      END IF
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &ccsdt_lr_alpha2_9_15_5',9,MA_ERR)
      IF ((h14b .le. h15b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h14b-1
     &),int_mb(k_range+h2b-1),int_mb(k_range+h15b-1),int_mb(k_range+p12b
     &-1),2,4,1,3,4.0d0/1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h15b 
     &- 1 + noab * (h14b - 1 + noab * (p12b - noab - 1 + nvab * (h2b - 1
     &)))))
      END IF
      IF ((h15b .le. h14b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h14b-1
     &),int_mb(k_range+h2b-1),int_mb(k_range+h15b-1),int_mb(k_range+p12b
     &-1),2,4,3,1,-4.0d0/1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h14b 
     &- 1 + noab * (h15b - 1 + noab * (p12b - noab - 1 + nvab * (h2b - 1
     &)))))
      END IF
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('ccsdt_lr_alpha2_9_15_5',
     &10,MA_ERR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('ccsdt_lr_alpha2_9_1
     &5_5',11,MA_ERR)
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
