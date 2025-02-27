      SUBROUTINE ccsdt_lr_alpha2_4_3(d_a,k_a_offset,d_b,k_b_offset,d_c,k
     &_c_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i1 ( h6 h7 h11 h12 )_ytratrb + = 1 * Sum ( p1 ) * trb ( p1 h11 )_trb * i2 ( h6 h7 h12 p1 )_ytra
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
      INTEGER h6b
      INTEGER h7b
      INTEGER h11b
      INTEGER h12b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p1b
      INTEGER p1b_1
      INTEGER h11b_1
      INTEGER h6b_2
      INTEGER h7b_2
      INTEGER h12b_2
      INTEGER p1b_2
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
      DO h6b = 1,noab
      DO h7b = h6b,noab
      DO h11b = 1,noab
      DO h12b = 1,noab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1
     &)+int_mb(k_spin+h11b-1)+int_mb(k_spin+h12b-1).ne.8)) THEN
      IF (int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1) .eq. int_mb(k_spin+h
     &11b-1)+int_mb(k_spin+h12b-1)) THEN
      IF (ieor(int_mb(k_sym+h6b-1),ieor(int_mb(k_sym+h7b-1),ieor(int_mb(
     &k_sym+h11b-1),int_mb(k_sym+h12b-1)))) .eq. ieor(irrep_y,ieor(irrep
     &_tra,irrep_trb))) THEN
      dimc = int_mb(k_range+h6b-1) * int_mb(k_range+h7b-1) * int_mb(k_ra
     &nge+h11b-1) * int_mb(k_range+h12b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('ccsdt_lr_alpha2_4_3',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      DO p1b = noab+1,noab+nvab
      IF (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h11b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h11b-1)) .eq. irrep_trb)
     & THEN
      CALL TCE_RESTRICTED_2(p1b,h11b,p1b_1,h11b_1)
      CALL TCE_RESTRICTED_4(h6b,h7b,h12b,p1b,h6b_2,h7b_2,h12b_2,p1b_2)
      dim_common = int_mb(k_range+p1b-1)
      dima_sort = int_mb(k_range+h11b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+h6b-1) * int_mb(k_range+h7b-1) * int_mb
     &(k_range+h12b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsdt_lr_alpha2_4_3',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsdt_lr_alpha2_4_3',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h11b_
     &1 - 1 + noab * (p1b_1 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p1b-1)
     &,int_mb(k_range+h11b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsdt_lr_alpha2_4_3',3,M
     &A_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsdt_lr_alpha2_4_3',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &ccsdt_lr_alpha2_4_3',5,MA_ERR)
      IF ((h12b .le. p1b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p1b_2
     & - noab - 1 + nvab * (h12b_2 - 1 + noab * (h7b_2 - 1 + noab * (h6b
     &_2 - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h6b-1)
     &,int_mb(k_range+h7b-1),int_mb(k_range+h12b-1),int_mb(k_range+p1b-1
     &),3,2,1,4,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('ccsdt_lr_alpha2_4_3',6,M
     &A_ERR)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('ccsdt_lr_alpha2_4_3
     &',7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('ccsdt_lr_alpha2_4_3
     &',8,MA_ERR)
      END IF
      END IF
      END IF
      END DO
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &ccsdt_lr_alpha2_4_3',9,MA_ERR)
      IF ((h11b .le. h12b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h12b-1
     &),int_mb(k_range+h7b-1),int_mb(k_range+h6b-1),int_mb(k_range+h11b-
     &1),3,2,4,1,1.0d0/2.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h12b 
     &- 1 + noab * (h11b - 1 + noab * (h7b - 1 + noab * (h6b - 1)))))
      END IF
      IF ((h12b .le. h11b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h12b-1
     &),int_mb(k_range+h7b-1),int_mb(k_range+h6b-1),int_mb(k_range+h11b-
     &1),3,2,1,4,-1.0d0/2.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h11b 
     &- 1 + noab * (h12b - 1 + noab * (h7b - 1 + noab * (h6b - 1)))))
      END IF
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('ccsdt_lr_alpha2_4_3',10,
     &MA_ERR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('ccsdt_lr_alpha2_4_3
     &',11,MA_ERR)
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
