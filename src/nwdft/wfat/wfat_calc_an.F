c============================================================
c============================================================
      subroutine wfat_calc_an(bt, gm, nat, nb, ns, nch, nqn, lmax,
     &     orbid, spinid, dgid0, dgn0, mid, dpmo, chmmt, pb, intl_j, 
     &     intl_k1, intl_k2, scr, scrx, chmmtd, lfpb, srfc, an)

      implicit none
#include "errquit.fh"
      double precision bt     ! [input], beta in degrees
      double precision gm     ! [input], gamma in degrees
      integer nat
      integer nb
      integer ns
      integer nch
      integer nqn
      integer lmax
      integer orbid(nch)
      integer spinid(nch)
      integer dgid0(nch)                                    ! [input],
      integer dgn0(nch)                                     ! [input],
      integer mid(nqn)                                      ! [input],
      double precision dpmo(nb,nb,ns,3)                     ! [input],
      double precision chmmt(3,nch)                         ! [input],
      double precision pb(3,3,nch)                          ! [input],
      double complex intl_j((1+lmax)**2, nat, nqn, 3, nch)  ! [input],
      double complex intl_k1((1+lmax)**2, nat, nqn, nch)    ! [input],
      double complex intl_k2((1+lmax)**2, nat, nqn, nch)    ! [input],
      double precision scr(1:*)                             ! [output], minimum size is 2*nb**2 + nb
      double complex scrx(1:*)                              ! [output], minimum size is nat*nqn*nch
      double precision chmmtd(3,nch)                        ! [output],
      double precision lfpb(nch)                            ! [output],
      double complex srfc(nat, nqn, nch)                    ! [in/output],
      double complex an(nat, nqn, nch)                      ! [output],

      logical iscalc(nch), isdiag, dgtest
      integer nht, ndp, nev, nlc, p_ht0, p_htn, p_dp0, p_dpn, p_ev0
      integer p_evn, p_lc0, p_lcn, i, is, nn, k, j, j_, mfid(nch)
      integer orbid0(nch)
      double precision protmat(3,3), cu(3), diagthr

      parameter(diagthr=1.0d-12)

c     ==== Double precision scratch ====
      nht = nb**2
      ndp = nb**2
      nev = nb
      
      p_ht0 = 1                 ! p_ht0 = pointer to the (dynamic) submatrix of H_int
      p_htn = p_ht0 + nht - 1
      p_dp0 = p_htn + 1         ! p_dp0 = pointer to the (dynamic) submatrix of dpmo
      p_dpn = p_dp0 + ndp - 1
      p_ev0 = p_dpn + 1         ! p_ev0 = pointer to the (dynamic) eigenvalue of the submatrix of H_int
      p_evn = p_ev0 + nev - 1

c     Double complex scratch ====
      nlc = nat * nqn * nch
      
      p_lc0 = 1
      p_lcn = p_lc0 + nlc - 1

      
c     ==== Get the passive rotation matrix ===
      call wfat_passrot(-bt, -gm, protmat)
      cu = protmat(:,3)
      

c     Diagonalize the dipole interaction Hamiltonian within any
c     degenerate manifolds that are ionized.
      chmmtd = chmmt

      iscalc = .false.
      mfid = 0      
      do i = 1, nch
         is = spinid(i)
         nn = dgn0(i)
         if (.not. iscalc(i)) then
            
c           ==== Get the degenerate manifold ID ====
            call wfat_get_degen_mfold(nch, orbid, spinid, dgid0, dgn0,
     &           i, mfid, orbid0)
            
            dgtest = (nn > 1)
            isdiag = .true.
            if (dgtest) then
c              ==== Construct the H_int within the degenerate MO's ====
               scr(p_ht0:p_ht0+nn**2-1) = 0.0d0
               do k = 1, 3
                  call wfat_submatrix(nb, nb, dpmo(:,:,is,k), nn, 
     &                 nn, orbid0, orbid0, scr(p_dp0))
                  scr(p_ht0:p_ht0+nn**2-1) = scr(p_ht0:p_ht0+nn**2-1) +
     &                                    cu(k)*scr(p_dp0:p_dp0+nn**2-1)
               enddo
               call wfat_isdiag(nn, scr(p_ht0), diagthr, isdiag)
         
               
               if (.not. isdiag) then
                  call errquit('WFAT: For the moment, degenerate '//
     &                 'ionized channels must diagonalize the '//
     &                 'interaction Hamiltonian', 48, UERR)
                  
c                 ==== Diagonalize H_int in this manifold ====
                  call wfat_dsyev_wrap(nn, scr(p_ht0), scr(p_ev0))
         
c                 ==== Construct the new structure factors ====
                  call wfat_lincomb_srfc(nat, nch, nqn, nn,
     &                 scr(p_ht0), mfid, scrx(p_lc0), srfc)
         
c                 ==== Update first moments of ionized orbitals ====
cforlater                  call wfat_update_orbmmt(chmmtd)
               endif
            endif

c           ==== Calculate a_\nu and polarizabilities ====
            do j = 1, nn
               j_ = mfid(j)
               call wfat_calc_an0(nat, nqn, mid, lmax, bt, gm, 
     &              chmmtd(:,j_), intl_j(:,:,:,:,j_),
     &              intl_k1(:,:,:,j_), intl_k2(:,:,:,j_), an(:,:,j_))
               call wfat_calc_lfpolb(protmat, pb(:,:,j_), lfpb(j_))
            enddo

c           ==== Update a_\nu and polarizabilities ====
            if (dgtest) then
               if (.not. isdiag) then
                  call wfat_lincomb_srfc(nat, nch, nqn, nn,
     &                 scr(p_ht0), mfid, scrx(p_lc0), an)
cforlater                  call wfat_update_polb()
               endif
            endif
         
            do j = 1, nn
               j_ = mfid(j)
               iscalc(j_) = .true.
            enddo

         endif
      enddo
      
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
      subroutine wfat_get_degen_mfold(nch, orbid, spinid, dgid0, dgn0, 
     &     ich, mfid, orbid0)
c     Given ich, :
c     1) mfid will contain the ID's of the ionized channels that are
c        in the same degenerate manifold and has the same spin as ich (it
c        includes ich too) , and
c     2) orbid0 will contain the ID's of the MO's that are in the same
c        degenerate manifold and has the same spin as the MO ich corresponds 
c        to (it includes the MO corresponding to ich too).
c     Hence, for a given ich, mfid and orbid0 constitute two equivalent
c     representations of the same set of degenerate MO's.

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
      integer cid(n)                   ! [input], 1 <= cid(i) <= nch
      double complex scr(nat, nqn, n)  ! [output],
      double complex s(nat, nqn, nch)  ! [in/output],

      integer j, j_, i, i_
      double complex czero
      parameter ( czero = cmplx(0.0d0, 0.0d0, kind=kind(1.0d0)) )

      
c     ==== Construct the linear combinations ====
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


c============================================================
c============================================================
      subroutine wfat_update_orbmmt()

      implicit none
cforlater      double precision scr(n,n)
cforlater      double precision mmtd()      ! [output], New orbital moments in molecular frame and relative to the corresponding orbital center.



      end
c============================================================


c============================================================
c============================================================
      subroutine wfat_calc_lfpolb(protmat, pb, lfpb)

      implicit none
      double precision protmat(3,3)
      double precision pb(3,3)
      double precision lfpb

      call wfat_1e_expval_1eopr(3, protmat(:,3), pb, lfpb)
      lfpb = -2 * lfpb
      
      end
c============================================================
      
      
