*
* $Id$
*
      

*    ************************************
*    *					*
*    *		v_bwexc		        *
*    *					*
*    ************************************
*
      
      subroutine v_bwexc(gga,n2ft3d,ispin,dn,
     >                   x_parameter,c_parameter,
     >                   xcp,xce)
      implicit none
      integer gga
      integer n2ft3d
      integer  ispin
      real*8  dn(n2ft3d,2)
      real*8  x_parameter,c_parameter
      real*8  xcp(n2ft3d,2),xce(n2ft3d)


#include "bafdecls.fh"
#include "errquit.fh"
#include "nwxc.fh"
#include "util.fh"


c     **** local variables ****      
      logical value, use_nwxc
      integer nx,ny,nz
      real*8  scal1
      integer rho(2),grx(2),gry(2),grz(2)
      integer agr(3),fn(2),fdn(3),tmp(2),rhog(2)

      integer rhoup(2),grupx(2),grupy(2),grupz(2)
      integer rhodn(2),grdnx(2),grdny(2),grdnz(2)
      integer          grallx(2),grally(2),grallz(2)
      integer          grad(2),gtmp(2)
      integer xagr(2),xfn(2),xfdn(2)
      integer k
      double precision dncut
      parameter (dncut = 1.0d-30)
      double precision dumtau

*     **** external functions ****
      integer  G_indx
      external G_indx

      call nwpw_timing_start(4)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)

      use_nwxc = .false.
      if (util_module_avail("nwxc")) then
         use_nwxc = nwxc_is_on()
      endif
                
*     **********************************
*     ***** restricted calculation *****
*     **********************************
      if (ispin.eq.1) then

c        ***** tempory variables needed rho,grx,gry,grz *****
c        *****                          agr,fn,fdn      *****
        value = BA_push_get(mt_dbl,n2ft3d,'rho', rho(2), rho(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grx',grx(2),grx(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'gry',gry(2),gry(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grz',grz(2),grz(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'agr',agr(2),agr(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'fn',fn(2),fn(1))
        rhog(1) = fn(1)
        value = value.and.
     >        BA_push_get(mt_dbl, 2*n2ft3d,'fdn',fdn(2),fdn(1))
        tmp(1) = fdn(1)
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      !call dcopy(n2ft3d,0.0d0,0,dbl_mb(rho(1)),1)
      !call dcopy(n2ft3d,0.0d0,0,dbl_mb(agr(1)),1)
      call Parallel_shared_vector_zero(.true.,n2ft3d,dbl_mb(agr(1)))
      call Parallel_shared_vector_zero(.true.,n2ft3d,dbl_mb(rho(1)))



c        ***** calculate rho tmp1=rho(g) ****
         call D3dB_rr_Sum(1,dn(1,1),dn(1,1),dbl_mb(rho(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rho(1)))
         call D3dB_r_SMul(1,scal1,dbl_mb(rho(1)),dbl_mb(rhog(1)))
         call D3dB_rc_fft3f(1,dbl_mb(rhog(1)))
         call mask_C(0,dbl_mb(rhog(1)))

c        ***** calculated  grup= grad n ****
         call D3dB_ic_Mul(1,dbl_mb(G_indx(1)),
     >                      dbl_mb(rhog(1)),
     >                      dbl_mb(grx(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(2)),
     >                      dbl_mb(rhog(1)),
     >                      dbl_mb(gry(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(3)),
     >                      dbl_mb(rhog(1)),
     >                      dbl_mb(grz(1)))

         call D3dB_cr_fft3b(1,dbl_mb(grx(1)))
         call D3dB_cr_fft3b(1,dbl_mb(gry(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grz(1)))


c        ***** calculate agr = |grad n| ****
         call D3dB_rr_Sqr(1,dbl_mb(grx(1)),
     >                      dbl_mb(agr(1)))
         call D3dB_rr_Sqr(1,dbl_mb(gry(1)),
     >                      dbl_mb(tmp(1)))

         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(1)))

         call D3dB_rr_Sqr(1,dbl_mb(grz(1)),
     >                      dbl_mb(tmp(1)))

         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(1)))

         if (use_nwxc) then
!$OMP MASTER
           call nwxc_eval_df(1,n2ft3d,dbl_mb(rho(1)),dbl_mb(agr(1)),
     >                       dumtau,xce,
     >                       dbl_mb(fn(1)),dbl_mb(fdn(1)),dumtau)
!$OMP END MASTER
!$OMP BARRIER
c
c          Combine (df/d|grad a|) with (df/d(grad a|grad b))
c
           call D3dB_r_SMul(1,0.5d0,dbl_mb(fdn(1)+n2ft3d),
     >                              dbl_mb(fdn(1)))
c
c          Calculate energy density from energy
c
!$OMP DO
           do k = 1, n2ft3d
             xce(k) = xce(k)/(dbl_mb(rho(1)-1+k)+dncut)
           enddo
!$OMP END DO
           call D3dB_rr_Sqrt1(1,dbl_mb(agr(1)))
           call D3dB_rr_Mul2(1,dbl_mb(agr(1)),dbl_mb(fdn(1)))
         else
           call D3dB_rr_Sqrt1(1,dbl_mb(agr(1)))

c        ***** calculate xce,fn=df/dn, and fdn=df/d|grad n|  ****
         if (gga.eq.10) then
         call gen_PBE96_BW_restricted(n2ft3d,
     >                                dbl_mb(rho(1)),
     >                                dbl_mb(agr(1)),
     >                                x_parameter,c_parameter,
     >                                xce,
     >                                dbl_mb(fn(1)),
     >                                dbl_mb(fdn(1)))
         else if (gga.eq.11) then
         call gen_BLYP_BW_restricted(n2ft3d,
     >                                dbl_mb(rho(1)),
     >                                dbl_mb(agr(1)),
     >                                x_parameter,c_parameter,
     >                                xce,
     >                                dbl_mb(fn(1)),
     >                                dbl_mb(fdn(1)))

         else if (gga.eq.12) then
         call gen_revPBE_BW_restricted(n2ft3d,
     >                                dbl_mb(rho(1)),
     >                                dbl_mb(agr(1)),
     >                                x_parameter,c_parameter,
     >                                xce,
     >                                dbl_mb(fn(1)),
     >                                dbl_mb(fdn(1)))
         else if (gga.eq.13) then
         call gen_PBEsol_BW_restricted(n2ft3d,
     >                                dbl_mb(rho(1)),
     >                                dbl_mb(agr(1)),
     >                                x_parameter,c_parameter,
     >                                xce,
     >                                dbl_mb(fn(1)),
     >                                dbl_mb(fdn(1)))
         else if (gga.eq.14) then
         call gen_HSE_BW_restricted(n2ft3d,
     >                              dbl_mb(rho(1)),
     >                              dbl_mb(agr(1)),
     >                              x_parameter,c_parameter,
     >                              xce,
     >                              dbl_mb(fn(1)),
     >                              dbl_mb(fdn(1)))
         else if (gga.eq.15) then
         call gen_B3LYP_BW_restricted(n2ft3d,
     >                              dbl_mb(rho(1)),
     >                              dbl_mb(agr(1)),
     >                              x_parameter,c_parameter,
     >                              xce,
     >                              dbl_mb(fn(1)),
     >                              dbl_mb(fdn(1)))
         else
          call errquit('bad gga',0,0)
         end if
         endif

c        ***** calculate df/d|grad n| *(grad n)/|grad n| ****
c         call D3dB_rr_Divide(1,dbl_mb(grx(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(grx(1)))
c         call D3dB_rr_Divide(1,dbl_mb(gry(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(gry(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grz(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(grz(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(grx(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(gry(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(grz(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grx(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(grx(1)))
c         call D3dB_rr_Mul(1,dbl_mb(gry(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(gry(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grz(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(grz(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(grx(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(gry(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(grz(1)))

c         call D3dB_r_SMul(1,scal1,dbl_mb(grx(1)),
c     >                            dbl_mb(grx(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(gry(1)),
c     >                            dbl_mb(gry(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grz(1)),
c     >                            dbl_mb(grz(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grx(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(gry(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grz(1)))

         call D3dB_r_Zero_Ends(1,dbl_mb(grx(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(gry(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grz(1)))

         call D3dB_rc_fft3f(1,dbl_mb(grx(1)))
         call D3dB_rc_fft3f(1,dbl_mb(gry(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grz(1)))
          
c         call D3dB_ic_Mul(1,dbl_mb(G_indx(1)),
c     >                      dbl_mb(grx(1)),
c     >                      dbl_mb(grx(1)))
c         call D3dB_ic_Mul(1,dbl_mb(G_indx(2)),
c     >                      dbl_mb(gry(1)),
c     >                      dbl_mb(gry(1)))
c         call D3dB_ic_Mul(1,dbl_mb(G_indx(3)),
c     >                      dbl_mb(grz(1)),
c     >                      dbl_mb(grz(1)))
         call D3dB_ic_Mul2(1,dbl_mb(G_indx(1)),dbl_mb(grx(1)))
         call D3dB_ic_Mul2(1,dbl_mb(G_indx(2)),dbl_mb(gry(1)))
         call D3dB_ic_Mul2(1,dbl_mb(G_indx(3)),dbl_mb(grz(1)))


         call D3dB_cc_Sum(1,dbl_mb(grx(1)),
     >                      dbl_mb(gry(1)),
     >                      dbl_mb(fdn(1)))

c         call D3dB_cc_Sum(1,dbl_mb(grz(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(fdn(1)))
         call D3dB_cc_Sum2(1,dbl_mb(grz(1)),dbl_mb(fdn(1)))

         call mask_C(0,dbl_mb(fdn(1)))

         call D3dB_cr_fft3b(1,dbl_mb(fdn(1)))
         call D3dB_rr_Minus(1,dbl_mb(fn(1)),
     >                        dbl_mb(fdn(1)),
     >                        xcp(1,1))
         call D3dB_r_Zero_Ends(1,xcp(1,1))

*        **** deallocate temporary memory ****
         value = BA_pop_stack(fdn(2))
         value = value.and.BA_pop_stack(fn(2))
         value = value.and.BA_pop_stack(agr(2))
         value = value.and.BA_pop_stack(grz(2))
         value = value.and.BA_pop_stack(gry(2))
         value = value.and.BA_pop_stack(grx(2))
         value = value.and.BA_pop_stack(rho(2))
         if (.not. value) call errquit('cannot pop stack memory',0,
     &       MA_ERR)



*     *******************************************************
*     ***** unrestricted calculation                    *****
*     *******************************************************
      else

c        ***** tempory variables needed rho,grx,gry,grz *****
c        *****                          agr,fn,fdn      *****
        value = BA_push_get(mt_dbl,2*n2ft3d,'rhoup', rhoup(2), rhoup(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grupx',grupx(2),grupx(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grupy',grupy(2),grupy(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grupz',grupz(2),grupz(1))

        value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'rhodn', rhodn(2), rhodn(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grdnx',grdnx(2),grdnx(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grdny',grdny(2),grdny(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grdnz',grdnz(2),grdnz(1))

        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grallx',grallx(2),grallx(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grally',grally(2),grally(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grallz',grallz(2),grallz(1))

        value = value.and.
     >        BA_push_get(mt_dbl, 3*n2ft3d,'xagr',xagr(2),xagr(1))
        agr(1) = xagr(1)
        agr(2) = xagr(1) +   n2ft3d
        agr(3) = xagr(1) + 2*n2ft3d
        value = value.and.
     >        BA_push_get(mt_dbl, 3*n2ft3d,'grad',grad(2),grad(1))
        value = value.and.
     >        BA_push_get(mt_dbl, 3*n2ft3d,'tmp',gtmp(2),gtmp(1))
        value = value.and.
     >        BA_push_get(mt_dbl, 2*n2ft3d,'xfn',xfn(2),xfn(1))
        fn(1) = xfn(1)
        fn(2) = xfn(1)+n2ft3d
        value = value.and.
     >        BA_push_get(mt_dbl, 3*n2ft3d,'xfdn',xfdn(2),xfdn(1))
        fdn(1) = xfdn(1)
        fdn(2) = xfdn(1) +   n2ft3d
        fdn(3) = xfdn(1) + 2*n2ft3d
        tmp(1) = xfdn(1)
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
c      call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhoup(1)),1)
c      call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhodn(1)),1)
c      call dcopy(3*n2ft3d,0.0d0,0,dbl_mb(xagr(1)),1)
      call Parallel_shared_vector_zero(.false.,3*n2ft3d,dbl_mb(xagr(1)))
      call Parallel_shared_vector_zero(.false.,n2ft3d,dbl_mb(rhodn(1)))
      call Parallel_shared_vector_zero(.true.,n2ft3d,dbl_mb(rhoup(1)))

c        ***** calculate rhoup  ****
         call D3dB_r_Copy(1,dn(1,1),dbl_mb(rhoup(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rhoup(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(rhoup(1)),
c     >                            dbl_mb(rhoup(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(rhoup(1)))
!MATHIAS OPENMP OK
         call D3dB_rc_fft3f(1,dbl_mb(rhoup(1)))
         call mask_C(0,dbl_mb(rhoup(1)))

c        ***** calculate   grup= grad nup ****
         call D3dB_ic_Mul(1,dbl_mb(G_indx(1)),
     >                      dbl_mb(rhoup(1)),
     >                      dbl_mb(grupx(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(2)),
     >                      dbl_mb(rhoup(1)),
     >                      dbl_mb(grupy(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(3)),
     >                      dbl_mb(rhoup(1)),
     >                      dbl_mb(grupz(1)))
!MATHIAS OPENMP OK
         call D3dB_cr_fft3b(1,dbl_mb(grupx(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grupy(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grupz(1)))

c        ***** calculate agrup = |grad nup| ****
         call D3dB_rr_Sqr(1,dbl_mb(grupx(1)),
     >                      dbl_mb(agr(1)))
         call D3dB_rr_Sqr(1,dbl_mb(grupy(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(1)),
c     >                      dbl_mb(agr(1)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(1)))

         call D3dB_rr_Sqr(1,dbl_mb(grupz(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(1)),
c     >                      dbl_mb(agr(1)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(1)))

c         call D3dB_rr_Sqrt(1,dbl_mb(agr(1)),
c     >                       dbl_mb(agr(1)))     
         call D3dB_rr_Sqrt1(1,dbl_mb(agr(1)))

c        ***** calculate rhodn  ****
         call D3dB_r_Copy(1,dn(1,2),dbl_mb(rhodn(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rhodn(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(rhodn(1)),
c     >                            dbl_mb(rhodn(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(rhodn(1)))
!MATHIAS OPENMP OK
         call D3dB_rc_fft3f(1,dbl_mb(rhodn(1)))
         call mask_C(0,dbl_mb(rhodn(1)))

c        ***** calculate   grdn = grad ndn ****
         call D3dB_ic_Mul(1,dbl_mb(G_indx(1)),
     >                      dbl_mb(rhodn(1)),
     >                      dbl_mb(grdnx(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(2)),
     >                      dbl_mb(rhodn(1)),
     >                      dbl_mb(grdny(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(3)),
     >                      dbl_mb(rhodn(1)),
     >                      dbl_mb(grdnz(1)))
!MATHIAS OPENMP OK
         call D3dB_cr_fft3b(1,dbl_mb(grdnx(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grdny(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grdnz(1)))

c        ***** calculate agrdn = |grad ndn| ****
         call D3dB_rr_Sqr(1,dbl_mb(grdnx(1)),
     >                      dbl_mb(agr(2)))
         call D3dB_rr_Sqr(1,dbl_mb(grdny(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(2)),
c     >                      dbl_mb(agr(2)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(2)))

         call D3dB_rr_Sqr(1,dbl_mb(grdnz(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(2)),
c     >                      dbl_mb(agr(2)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(2)))

c         call D3dB_rr_Sqrt(1,dbl_mb(agr(2)),
c     >                       dbl_mb(agr(2)))     
         call D3dB_rr_Sqrt1(1,dbl_mb(agr(2)))

c        ***** calculate agr = |grad nup +grad ndn| ****
         call D3dB_rr_Sum(1,dbl_mb(grupx(1)),
     >                      dbl_mb(grdnx(1)),
     >                      dbl_mb(grallx(1)))
         call D3dB_rr_Sum(1,dbl_mb(grupy(1)),
     >                      dbl_mb(grdny(1)),
     >                      dbl_mb(grally(1)))
         call D3dB_rr_Sum(1,dbl_mb(grupz(1)),
     >                      dbl_mb(grdnz(1)),
     >                      dbl_mb(grallz(1)))
         call D3dB_rr_Sqr(1,dbl_mb(grallx(1)),
     >                      dbl_mb(agr(3)))

         call D3dB_rr_Sqr(1,dbl_mb(grally(1)),
     >                      dbl_mb(tmp(1)))
c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(3)),
c     >                      dbl_mb(agr(3)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(3)))

         call D3dB_rr_Sqr(1,dbl_mb(grallz(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(3)),
c     >                      dbl_mb(agr(3)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(3)))

c         call D3dB_rr_Sqrt(1,dbl_mb(agr(3)),
c     >                       dbl_mb(agr(3)))
         call D3dB_rr_Sqrt1(1,dbl_mb(agr(3)))

 
         if (use_nwxc) then
c
c          Copy |grad nup|->grad(1), |grad n|->grad(2), and |grad ndn|->grad(3)
c
           call D3dB_r_Copy(1,dbl_mb(agr(1)),dbl_mb(grad(1)))
           call D3dB_r_Copy(1,dbl_mb(agr(2)),dbl_mb(grad(1)+2*n2ft3d))
           call D3dB_r_Copy(1,dbl_mb(agr(3)),dbl_mb(grad(1)+n2ft3d))
c
c          Replace |grad x| with |grad x|^2
c
           call D3dB_rr_Sqr1(1,dbl_mb(grad(1)))
           call D3dB_rr_Sqr1(1,dbl_mb(grad(1)+2*n2ft3d))
           call D3dB_rr_Sqr1(1,dbl_mb(grad(1)+n2ft3d))
c
c          Replace |grad n|^2 with |(grad nup|grad ndn)|^2
c
           call D3dB_rr_Sub2(1,dbl_mb(grad(1)),dbl_mb(grad(1)+n2ft3d))
           call D3dB_rr_Sub2(1,dbl_mb(grad(1)+2*n2ft3d),
     +                         dbl_mb(grad(1)+n2ft3d))
           call D3dB_r_SMul1(1,0.5d0,dbl_mb(grad(1)+n2ft3d))
c
c          Evaluate the functional
c
!$OMP MASTER
           call nwxc_eval_df(2,n2ft3d,dn,dbl_mb(grad(1)),
     >                       dumtau,xce,
     >                       dbl_mb(fn(1)),dbl_mb(fdn(1)),dumtau)
!$OMP END MASTER
!$OMP BARRIER
c
c          Replace f with the energy density f/n
c
!$OMP DO
           do k = 1, n2ft3d
             xce(k) = xce(k)/(dn(k,1)+dn(k,2)+dncut)
           end do
!$OMP END DO
c          
c          Replace (df/d|grad nup|^2) with (df/d|grad nup|)
c
           call D3dB_rr_Mul2(1,dbl_mb(agr(1)),dbl_mb(fdn(1)))
           call D3dB_r_SMul1(1,2.0d0,dbl_mb(fdn(1)))
c
c          Replace (df/d|grad ndn|^2) with (df/d|grad ndn|)
c
           call D3dB_rr_Mul(1,dbl_mb(agr(2)),dbl_mb(fdn(3)),
     +                      dbl_mb(gtmp(1)+n2ft3d))
           call D3dB_r_SMul1(1,2.0d0,dbl_mb(gtmp(1)+n2ft3d))
c
c          Replace (df/d|(grad nup|grad ndn)|^2) with (df/d|grad n|)
c
           call D3dB_rr_Mul(1,dbl_mb(agr(3)),dbl_mb(fdn(2)),
     +                      dbl_mb(gtmp(1)+2*n2ft3d))
c
c          Put the results back into fdn
c
c           call dcopy(n2ft3d,dbl_mb(gtmp(1)+n2ft3d),1,
c     +                       dbl_mb(fdn(2)),1) 
c           call dcopy(n2ft3d,dbl_mb(gtmp(1)+2*n2ft3d),1,
c     +                       dbl_mb(fdn(3)),1) 
           call Parallel_shared_vector_copy(.false.,n2ft3d,
     >                                      dbl_mb(gtmp(1)+n2ft3d),
     >                                      dbl_mb(fdn(2)))
           call Parallel_shared_vector_copy(.true.,n2ft3d,
     >                                      dbl_mb(gtmp(1)+2*n2ft3d),
     >                                      dbl_mb(fdn(3)))
c
         else
c        ***** calculate xce,fn=(df/dnup,df/dndn) and               ****
c        *****  fdn=(dfx/d|grad nup|,dfx/d|grad ndn|,dfc/d|grad n|) ****
         if (gga.eq.10) then
         call gen_PBE96_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >                                  x_parameter,c_parameter,
     >                                  xce,
     >                                  dbl_mb(fn(1)),
     >                                  dbl_mb(fdn(1)))

         else if (gga.eq.11) then
         call gen_BLYP_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >                                  x_parameter,c_parameter,
     >                                  xce,
     >                                  dbl_mb(fn(1)),
     >                                  dbl_mb(fdn(1)))
         else if (gga.eq.12) then
         call gen_revPBE_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >                                  x_parameter,c_parameter,
     >                                  xce,
     >                                  dbl_mb(fn(1)),
     >                                  dbl_mb(fdn(1)))
         else if (gga.eq.13) then
         call gen_PBEsol_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >                                  x_parameter,c_parameter,
     >                                  xce,
     >                                  dbl_mb(fn(1)),
     >                                  dbl_mb(fdn(1)))
         else if (gga.eq.14) then
         call gen_HSE_BW_unrestricted(n2ft3d,dn,
     >                                dbl_mb(agr(1)),
     >                                x_parameter,c_parameter,
     >                                xce,
     >                                dbl_mb(fn(1)),
     >                                dbl_mb(fdn(1)))
         else if (gga.eq.15) then
         call gen_B3LYP_BW_unrestricted(n2ft3d,dn,
     >                                dbl_mb(agr(1)),
     >                                x_parameter,c_parameter,
     >                                xce,
     >                                dbl_mb(fn(1)),
     >                                dbl_mb(fdn(1)))
         else
          call errquit('bad gga',0,0)
         end if
         end if

         
*        **** calculate df/d|grad nup|* (grad nup)/|grad nup|  ****
*        **** calculate df/d|grad ndn|* (grad ndn)/|grad ndn|  ****
*        **** calculate df/d|grad n|  * (grad n)/|grad n|  ****
c         call D3dB_rr_Divide(1,dbl_mb(grupx(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(grupx(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grupy(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(grupy(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grupz(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(grupz(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grdnx(1)),
c     >                         dbl_mb(agr(2)),
c     >                         dbl_mb(grdnx(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grdny(1)),
c     >                         dbl_mb(agr(2)),
c     >                         dbl_mb(grdny(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grdnz(1)),
c     >                         dbl_mb(agr(2)),
c     >                         dbl_mb(grdnz(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grallx(1)),
c     >                         dbl_mb(agr(3)),
c     >                         dbl_mb(grallx(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grally(1)),
c     >                         dbl_mb(agr(3)),
c     >                         dbl_mb(grally(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grallz(1)),
c     >                         dbl_mb(agr(3)),
c     >                         dbl_mb(grallz(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(grupx(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(grupy(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(grupz(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(2)),dbl_mb(grdnx(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(2)),dbl_mb(grdny(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(2)),dbl_mb(grdnz(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(3)),dbl_mb(grallx(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(3)),dbl_mb(grally(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(3)),dbl_mb(grallz(1)))

c         call D3dB_rr_Mul(1,dbl_mb(grupx(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(grupx(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grupy(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(grupy(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grupz(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(grupz(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grdnx(1)),
c     >                      dbl_mb(fdn(2)),
c     >                      dbl_mb(grdnx(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grdny(1)),
c     >                      dbl_mb(fdn(2)),
c     >                      dbl_mb(grdny(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grdnz(1)),
c     >                      dbl_mb(fdn(2)),
c     >                      dbl_mb(grdnz(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grallx(1)),
c     >                      dbl_mb(fdn(3)),
c     >                      dbl_mb(grallx(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grally(1)),
c     >                      dbl_mb(fdn(3)),
c     >                      dbl_mb(grally(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grallz(1)),
c     >                      dbl_mb(fdn(3)),
c     >                      dbl_mb(grallz(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(grupx(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(grupy(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(grupz(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(2)),dbl_mb(grdnx(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(2)),dbl_mb(grdny(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(2)),dbl_mb(grdnz(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(3)),dbl_mb(grallx(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(3)),dbl_mb(grally(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(3)),dbl_mb(grallz(1)))


*        **** calculate (df/d|grad nup|* (grad nup)/|grad nup|)  ****
*        ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****
*        **** calculate (df/d|grad ndn|* (grad ndn)/|grad ndn|)  ****
*        ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****
c         call D3dB_rr_Sum(1,dbl_mb(grupx(1)),
c     >                      dbl_mb(grallx(1)),
c     >                      dbl_mb(grupx(1)))
c         call D3dB_rr_Sum(1,dbl_mb(grupy(1)),
c     >                      dbl_mb(grally(1)),
c     >                      dbl_mb(grupy(1)))
c         call D3dB_rr_Sum(1,dbl_mb(grupz(1)),
c     >                      dbl_mb(grallz(1)),
c     >                      dbl_mb(grupz(1)))
c         call D3dB_rr_Sum(1,dbl_mb(grdnx(1)),
c     >                      dbl_mb(grallx(1)),
c     >                      dbl_mb(grdnx(1)))
c         call D3dB_rr_Sum(1,dbl_mb(grdny(1)),
c     >                      dbl_mb(grally(1)),
c     >                      dbl_mb(grdny(1)))
c         call D3dB_rr_Sum(1,dbl_mb(grdnz(1)),
c     >                      dbl_mb(grallz(1)),
c     >                      dbl_mb(grdnz(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grallx(1)),dbl_mb(grupx(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grally(1)),dbl_mb(grupy(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grallz(1)),dbl_mb(grupz(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grallx(1)),dbl_mb(grdnx(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grally(1)),dbl_mb(grdny(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grallz(1)),dbl_mb(grdnz(1)))


c         call D3dB_r_SMul(1,scal1,dbl_mb(grupx(1)),
c     >                            dbl_mb(grupx(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grupy(1)),
c     >                            dbl_mb(grupy(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grupz(1)),
c     >                            dbl_mb(grupz(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grdnx(1)),
c     >                            dbl_mb(grdnx(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grdny(1)),
c     >                            dbl_mb(grdny(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grdnz(1)),
c     >                            dbl_mb(grdnz(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grupx(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grupy(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grupz(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grdnx(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grdny(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grdnz(1)))

         call D3dB_r_Zero_Ends(1,dbl_mb(grupx(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grupy(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grupz(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grdnx(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grdny(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grdnz(1)))

*        **** put sums in k-space ***
!MATHIAS OPENMP OK
         call D3dB_rc_fft3f(1,dbl_mb(grupx(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grupy(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grupz(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grdnx(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grdny(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grdnz(1)))

*        **** multiply sums by G vector ***
c         call D3dB_ic_Mul(1,dbl_mb(G_indx(1)),
c     >                      dbl_mb(grupx(1)),
c     >                      dbl_mb(grupx(1)))
c         call D3dB_ic_Mul(1,dbl_mb(G_indx(2)),
c     >                      dbl_mb(grupy(1)),
c     >                      dbl_mb(grupy(1)))
c         call D3dB_ic_Mul(1,dbl_mb(G_indx(3)),
c     >                      dbl_mb(grupz(1)),
c     >                      dbl_mb(grupz(1)))
c         call D3dB_ic_Mul(1,dbl_mb(G_indx(1)),
c     >                      dbl_mb(grdnx(1)),
c     >                      dbl_mb(grdnx(1)))
c         call D3dB_ic_Mul(1,dbl_mb(G_indx(2)),
c     >                      dbl_mb(grdny(1)),
c     >                      dbl_mb(grdny(1)))
c         call D3dB_ic_Mul(1,dbl_mb(G_indx(3)),
c     >                      dbl_mb(grdnz(1)),
c     >                      dbl_mb(grdnz(1)))
         call D3dB_ic_Mul2(1,dbl_mb(G_indx(1)),dbl_mb(grupx(1)))
         call D3dB_ic_Mul2(1,dbl_mb(G_indx(2)),dbl_mb(grupy(1)))
         call D3dB_ic_Mul2(1,dbl_mb(G_indx(3)),dbl_mb(grupz(1)))
         call D3dB_ic_Mul2(1,dbl_mb(G_indx(1)),dbl_mb(grdnx(1)))
         call D3dB_ic_Mul2(1,dbl_mb(G_indx(2)),dbl_mb(grdny(1)))
         call D3dB_ic_Mul2(1,dbl_mb(G_indx(3)),dbl_mb(grdnz(1)))

*        **** add up two dot products ****
         call D3dB_cc_Sum(1,dbl_mb(grupx(1)),
     >                      dbl_mb(grupy(1)),
     >                      dbl_mb(fdn(1)))

c         call D3dB_cc_Sum(1,dbl_mb(grupz(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(fdn(1)))
         call D3dB_cc_Sum2(1,dbl_mb(grupz(1)),dbl_mb(fdn(1)))

         call D3dB_cc_Sum(1,dbl_mb(grdnx(1)),
     >                      dbl_mb(grdny(1)),
     >                      dbl_mb(fdn(2)))

c         call D3dB_cc_Sum(1,dbl_mb(grdnz(1)),
c     >                      dbl_mb(fdn(2)),
c     >                      dbl_mb(fdn(2)))
         call D3dB_cc_Sum2(1,dbl_mb(grdnz(1)),dbl_mb(fdn(2)))

*        **** put back in r-space and subtract from df/dnup,df/dndn ****
         call mask_C(0,dbl_mb(fdn(1)))
         call mask_C(0,dbl_mb(fdn(2)))
!MATHIAS OPENMP OK
         call D3dB_cr_fft3b(1,dbl_mb(fdn(1)))
         call D3dB_cr_fft3b(1,dbl_mb(fdn(2)))
         call D3dB_rr_Minus(1,dbl_mb(fn(1)),
     >                        dbl_mb(fdn(1)),
     >                        xcp(1,1))
         call D3dB_rr_Minus(1,dbl_mb(fn(2)),
     >                        dbl_mb(fdn(2)),
     >                        xcp(1,2))
         call D3dB_r_Zero_Ends(1,xcp(1,1))
         call D3dB_r_Zero_Ends(1,xcp(1,2))


*        **** deallocate temporary memory ****
         value = BA_pop_stack(xfdn(2))
         value = value.and.BA_pop_stack(xfn(2))
         value = value.and.BA_pop_stack(gtmp(2))
         value = value.and.BA_pop_stack(grad(2))
         value = value.and.BA_pop_stack(xagr(2))

         value = value.and.BA_pop_stack(grallz(2))
         value = value.and.BA_pop_stack(grally(2))
         value = value.and.BA_pop_stack(grallx(2))
         value = value.and.BA_pop_stack(grdnz(2))
         value = value.and.BA_pop_stack(grdny(2))
         value = value.and.BA_pop_stack(grdnx(2))
         value = value.and.BA_pop_stack(rhodn(2))
         value = value.and.BA_pop_stack(grupz(2))
         value = value.and.BA_pop_stack(grupy(2))
         value = value.and.BA_pop_stack(grupx(2))
         value = value.and.BA_pop_stack(rhoup(2))
         if (.not. value) call errquit('cannot pop stack memory',0,
     &       MA_ERR)


      end if
      
      call nwpw_timing_end(4)
        
      return
      end



 
