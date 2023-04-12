c============================================================
c============================================================
      subroutine wfat_calc_an()

      implicit none
      integer nat
      integer nb
      integer ns
      integer nch
      integer nqn
      integer lmax
      double complex intl_j((1+lmax)**2, nat, nqn, 3, nch)  ! [input],
      double complex intl_k1((1+lmax)**2, nat, nqn, nch)    ! [input],
      double complex intl_k2((1+lmax)**2, nat, nqn, nch)    ! [input],
      double complex srfc(nat, nqn, nch)            ! [in/output],
      double complex an(nat, nqn, nch)              ! [output],

      integer dgid0(nch), dgn0(nch), mfid(nch)
      double precision protmat(3,3), cu(3)
      

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(-bt, -gm, protmat)
      cu = protmat(:,3)
      

      
c     Diagonalize the dipole interaction Hamiltonian within any degenerate
c     manifolds that are ionized.
      chmmtd = chmmt
      call wfat_detect_degen0(nch, orbid, spinid, dgid0, dgn0)

      mfid = 0      
      do i = 1, nch
c        ==== Get the degenerate manifold ID ====
         call wfat_get_degen_mfold(nch, orbid, spinid, dgid0, dgn0, i,
     &        mfid, orbid0)
         
         dgtest = (dgn0(i) > 1)
         isdiag = .true.
         if (dgtest) then
c           ==== Construct the H_int within the degenerate MO's ====
            subht = 0.0d0
            do k = 1, 3
               call wfat_submatrix(nb, nb, dpmo(:,:,k), dgn0(i), 
     &              dgn0(i), orbid0, orbid0, subdpmo(:,:,k))
               subht = subht + cu(k)*subdpmo(:,:,k)
            enddo
            call wfat_isdiag(dgn0(i), subht, diagthr, isdiag)

            
            if (.not. isdiag) then
c              ==== Diagonalize H_int in this manifold ====
               call wfat_dsyev_wrap(dgn0(i), subht, ev)

c              ==== Construct the new structure factors ====
               call wfat_lincomb_srfc(nat, nch, nqn, dgn0(i), subht, 
     &              mfid, dcpl_mb(k_sc), srfc)

               calculate new dipole moments
            endif
         endif

         calculate polarizabilities

         do j = 1, dgn0(i)
            j_ = mfid(j)
            call wfat_calc_an0(, , , intl_j(:,:,:,:,j_),
     &           intl_k1(:,:,:,j_), intl_k2(:,:,:,j_), an(:,:,j_))
         enddo
         
         if (dgtest) then
            if (.not. isdiag)
     &           call wfat_lincomb_srfc(nat, nch, nqn, dgn0(i), subht,
     &           mfid, dcpl_mb(k_sc), an)
         endif
      enddo
      
          
      chmmtd = 

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_calc_an0(nat, nqn, mid, lmax, bt, gm, mmt,
     &     intl_j, intl_k1, intl_k2, an)

      implicit none
      integer nat
      integer nqn
      integer mid(nqn)
      integer lmax
      double precision bt
      double precision gm
      double precision mmt(3)
      double complex intl_j((1+lmax)**2, nat, nqn, 3)  ! [input],
      double complex intl_k1((1+lmax)**2, nat, nqn)    ! [input],
      double complex intl_k2((1+lmax)**2, nat, nqn)    ! [input],
      double complex an(nat, nqn)              ! [output],

      integer iqn, ic, l, m, ish, i
      double precision protmat(3,3), cu(3), btr, gmr, mu(3), deg2rad
      double precision wgf
      double complex aa, bb, xx, czero, imnum
      parameter ( deg2rad = 0.01745329251994d0 )
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(-bt, -gm, protmat)
      cu = protmat(:,3)      
      btr = bt * deg2rad
      gmr = gm * deg2rad
      

c     IMPORTANT:
c     At this point, if ich corresponds to a degenerate state, mu and alpha must
c     be associated with the interaction-Hamiltonian-diagonalizing eigenstates,
c     which are a linear combination of the original degenerate states.
      mu = matmul(transpose(protmat), -mmt)
      do iqn = 1, nqn
         do ic = 1, nat

            an(ic, iqn) = czero
            do l = abs(mid(iqn)), lmax
               do m = -l, l
                  ish = l**2 + l + m + 1
                  
                  aa = intl_k1(ish,ic,iqn) +
     &                 mu(3)*intl_k2(ish,ic,iqn)

                  bb = czero
                  do i = 1, 3
                     bb = bb + cu(i)*intl_j(ish,ic,iqn,i)
                  enddo

                  call wfat_wignerf(l, mid(iqn), m, btr, wgf)
                  xx = (aa + bb) * wgf * exp(-imnum * m * gmr)

                  an(ic, iqn) = an(ic, iqn) + xx
               enddo
            enddo
            
         enddo
      enddo

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_detect_degen0(nch, orbid, spinid, dgid0, dgn0)

      implicit none
#include "wfat_qm_param.fh"
      integer nch          ! [input],
      integer orbid(nch)   ! [input],
      integer spinid(nch)  ! [input],
      integer dgid0(nch)   ! [output], the ID of the degenerate MO's whose MO ID and spin ID are specified in orbid and spinid, respectively.
      integer dgn0(nch)    ! [output], the degree of degeneracy of the above MO's.

      integer i, io, is

      do i = 1, nch
         io = orbid(i)
         is = spinid(i)
         dgid0(i) = int_mb(k_dgid + (is-1)*nbas0 + (io-1))
         dgn0(i) = int_mb(k_dgn + (is-1)*nbas0 + (io-1))
      enddo

      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_get_degen_mfold(nch, orbid, spinid, dgid0, dgn0, 
     &     ich, mfid, orbid0)

      implicit none
#include "errquit.fh"
      integer nch          ! [input],
      integer orbid(nch)   ! [input],
      integer spinid(nch)  ! [input],
      integer dgid0(nch)   ! [input], this array is produced by wfat_detect_degen0.
      integer dgn0(nch)    ! [input], this array is produced by wfat_detect_degen0.
      integer ich          ! [input], an integer whose value satisfies 1 <= ich <= nch.
      integer mfid(1:*)    ! [output], 1 <= mfid(i) <= nch, with i runs from 1 to dgn0(ich).
      integer orbid0(1:*)  ! [output], 1 <= orbid0(i) <= nbas0, with i runs from 1 to dgn0(ich).
      
      integer k, jch, j, j_, jo


      if (ich > nch) call errquit('WFAT: The argument ich of '//
     &     'wfat_get_degen_mfold must not be larger than nch', 34, UERR)
      
      k = 0
      do jch = 1, nch
         if (dgid0(jch) == dgid0(ich) .and. spinid(jch) == spinid(ich))
     &        then
            k = k + 1
            mfid(k) = jch
         endif
      enddo

      if (k /= dgn0(ich)) call errquit('WFAT: An inconsistency is '//
     &     'detected when identifying degenerate manifold. This is '//
     &     'most likely caused by not all degenerate MO''s within '//
     &     'a given degenerate manifold being included in the '//
     &     'ionized channels.', 23, UERR)

      do j = 1, dgn0(ich)
         j_ = mfid(j)
         jo = orbid(j_)
         orbid0(j) = jo
      enddo
      
      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_lincomb_srfc(nat, nch, nqn, n, cf, cid, scr, s)

      implicit none
      integer nat
      integer nch
      integer nqn
      integer n
      double precision cf(n,n)
      integer cid(n)
      double complex scr(nat, nqn, n)
      double complex s(nat, nqn, nch)

      integer j, j_, i, i_


      do j = 1, n
         scr(:,:,j) = czero
         do i = 1, n
            i_ = cid(i)
            scr(:,:,j) = scr(:,:,j) + cf(i,j)*s(:,:,i_)
         enddo
      enddo

c     ==== Replace the old structure factors with the new ones ====
      do j = 1, n
         j_ = cid(j)
         s(:,:,j_) = scr(:,:,j)
      enddo

      end
c============================================================
