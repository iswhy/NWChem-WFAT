c=======================================================
c=======================================================
      subroutine wfat_calc_polb(nch, orbid, spinid, dpmo, pb)
c     DESCRIPTION:
c     Computes the molecular frame static polarizability tensor for the MOs whose index
c     are contained in orbid. The tensor $\alpha$ is calculated as
c     $$
c     alpha_{ij} = \sum_n \langle u_m | x_i | u_n \rangle \langle u_n | x_j | u_m \rangle
c                  \frac{1}{e_m - e_n} .
c     $$
c     where u_n is the MO for which the polarizability tensor is computed, u_m is the MOs
c     (including u_n), and e_n is the orbital energy. i = 1,2,3.
      implicit none
#include "wfat_qm_param.fh"
#include "mafdecls.fh"
      integer nch                                      ! [input],
      integer orbid(nch)                               ! [input],
      integer spinid(nch)                              ! [input],
      double precision dpmo(nbas0, nbas0, nspin0, 3)   ! [input],
      double precision pb(3,3,nch)                     ! [output],

      integer i, io, is, j2, j1, m, id0, idm, ieg0, iegm
      double precision derg, pb0


      do i = 1, nch      ! Loop over ionization channels.
         io = orbid(i) 
         is = spinid(i)
         do j2 = 1, 3         ! Loop over Cartesian coordinates.
            do j1 = 1, 3         ! Loop over Cartesian coordinates.

               pb(j1,j2,i) = 0.0d0
               do m = 1, nbas0      ! Loop over the MOs.
                  id0 = int_mb(k_modgid + (is-1)*nbas0 + (io-1))    ! Index of the reference MO.
                  idm = int_mb(k_modgid + (is-1)*nbas0 + (m-1))     ! Index of the scanned MOs.
                  if (id0 /= idm) then
                     ieg0 = k_moerg + (is-1)*nbas0 + (io-1)
                     iegm = k_moerg + (is-1)*nbas0 + (m-1)
                     derg = dbl_mb(ieg0) - dbl_mb(iegm)

                     pb0 = dpmo(m,io,is,j1) * dpmo(m,io,is,j2) / derg
                     pb(j1,j2,i) = pb(j1,j2,i) + pb0
                  endif
               enddo

            enddo
         enddo
      enddo

      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_calc_polb_eco(n, orbid0, ispin, dpmo, pb)
c     DESCRIPTION :
c     The more economic version of the wfat_calc_polb subroutine.
      implicit none
#include "wfat_qm_param.fh"
#include "mafdecls.fh"
      integer n                          ! [input],
      integer orbid0(n)                  ! [input],
      integer ispin                      ! [input],
      double precision dpmo(nbas0, n, 3) ! [input],
      double precision pb(3,3,n)         ! [output],

      integer i, io, is, j2, j1, m, id0, idm, ieg0, iegm
      double precision derg, pb0


      is = ispin
      do i = 1, n
         io = orbid0(i)
         do j2 = 1, 3
            do j1 = 1, 3

               pb(j1,j2,i) = 0.0d0
               do m = 1, nbas0
                  id0 = int_mb(k_modgid + (is-1)*nbas0 + (io-1))
                  idm = int_mb(k_modgid + (is-1)*nbas0 + (m-1))
                  if (id0 /= idm) then
                     ieg0 = k_moerg + (is-1)*nbas0 + (io-1)
                     iegm = k_moerg + (is-1)*nbas0 + (m-1)
                     derg = dbl_mb(ieg0) - dbl_mb(iegm)

                     pb0 = dpmo(m,i,j1) * dpmo(m,i,j2) / derg
                     pb(j1,j2,i) = pb(j1,j2,i) + pb0
                  endif
               enddo

            enddo
         enddo
      enddo

      end
c=======================================================
      

c===============================================================
c===============================================================
      subroutine wfat_polb_dgn_sub(nch, orbid, nsub, maxn, dgs, dgss,
     &     ndgs, nb, ns, dpmo, erg, pbs)
c     DESCRIPTION:
c     nsub = The number of degenerate subspaces, a nondegenerate state is considered a
c            degenerate state in a subspace of dimension 1.
c     maxn = The largest dimension of all degenerate subspaces, i.e. maxn=maxval(ndgs).
c     dgs = The channel ID arranged into degenerate subspaces (columns). Each subspace
c           (column) must only have states of the same spin. However, note that dgs has
c           no information about the spins of its columns. This information is not needed
c           because dgs contains the indices of the channels which already encode the spin
c           of that channel. The range of its elements is 1 <= dgs(i,j) <= nch. Its
c           minimum size ismaxn*nsub.
c     dgss = The spin of the degenerate subspaces.
c     ndgs = The dimension of each degenerate subspace. Its minimum size is nsub.
c     pbs = The static polarizability matrix in molecular frame in each degenerate
c           subspace. The 6 indices in the 3rd dimension correspond to xx, yy, zz, xy, xz,
c           yz. Its minimum size is maxn*maxn*6*nsub.

      implicit none
#include "errquit.fh"
      integer nch                             ! [input]
      integer orbid(nch)                      ! [input]
      integer nsub                            ! [input]
      integer maxn                            ! [input]
      integer dgs(maxn,nsub)                  ! [input]
      integer dgss(nsub)                      ! [input]
      integer ndgs(nsub)                      ! [input]
      integer nb                              ! [input]
      integer ns                              ! [input]
      double precision dpmo(nb,nb,ns,3)       ! [input]
      double precision erg(nb,ns)             ! [input]
      double precision pbs(maxn,maxn,6,nsub)  ! [output]

      integer i, is, nn, j, jj, idv(maxn), k, k1, k2, j1, j2
      integer cartid(2,6)


      if (maxn < maxval(ndgs)) call errquit('WFAT: Incompatibility '//
     &     'between maxn and ndgs detected.', 55, UERR)

c     ==== Cartesian ID translation ====
      cartid(:,1) = (/ 1, 1 /)
      cartid(:,2) = (/ 2, 2 /)
      cartid(:,3) = (/ 3, 3 /)
      cartid(:,4) = (/ 1, 2 /)
      cartid(:,5) = (/ 1, 3 /)
      cartid(:,6) = (/ 2, 3 /)

      
      do i = 1, nsub    ! Loop over degenerate subspaces.
         is = dgss(i)
         nn = ndgs(i)

c        ==== Construct the degenerate subspace orbital ID array ====
         do j = 1, nn        ! Loop over states in the current degenerate subspace.
            jj = dgs(j,i)         ! 1 <= jj <= nch
            idv(j) = orbid(jj)    ! 1 <= idv(j) <= nb
         enddo

c        ==== Compute the static polarizabilities ====
c        ====   within this degenerate subspace   ====
         do k = 1, 6         ! Loop over Cartesian components of polarizability tensor.
            k1 = cartid(1,k)
            k2 = cartid(2,k)
            do j2 = 1, nn
               do j1 = 1, nn
                  call wfat_calc_polb0(nn, idv, nb, 
     &                 dpmo(:,idv(j1),is,k1), dpmo(:,idv(j2),is,k2),
     &                 erg(:,is), pbs(j1,j2,k,i))
               enddo
            enddo
         enddo

      enddo

      end
c===============================================================


c===============================================================
c===============================================================
      subroutine wfat_calc_polb0(ndg, dgnid, nb, dp1, dp2, erg, pb)
c     Computes the following quantity:
c     $$
c        p = -2 \sum_{i \notin DG}^N d_1[i] * d_2[i] / (e0-e[i])
c     $$
c     where DG contains list of integers to be excluded from the sum and e0
c     is equal to e[i'] with i' being one of the integers in the list DG. This
c     implies that the value of e[i] for any i in DG must be (numerically)
c     identical. This subroutine can be used to compute the matrix elements of
c     a Cartesian component of the static polarizability tensor operator
c     where in this case $d_1$ and $d_2$ are arrays that hold
c     $\rangle \phi_{i'} | x_1 | \psi \rangle$ and
c     $\rangle \psi | x_2 | \phi_{j'} \rangle$, respectively (the different $\psi$'s
c     are represented by the indices of these arrays). The sum runs over
c     $\psi$'s. This subroutine returns the $(i'j')$-th element of the $x_1 x_2$
c     component of the static polarizability tensor operator within the degenerate
c     subspace defined by the list DG. $\phi_{i'}$ and $\phi_{j'}$ must have their
c     indices contained in DG.
      
c     ndg = The lenght of dgnid.
c     dgnid = The list $DG$, the range of its elements is 1 <= dgnid(i) <= nb.
c     nb = The number of bases, $N$.
c     dp1 = The array holding the first array of dipole moment matrix elements, $d_1$.
c     dp2 = The array holding the second array of dipole moment matrix elements, $d_2$.
c     erg = The energy of each $\psi$ used in dp1 and dp2, $e$
c     pb = $p$.
      implicit none
      integer ndg               ! [input]
      integer dgnid(ndg)        ! [input]
      integer nb                ! [input]
      double precision dp1(nb)  ! [input]
      double precision dp2(nb)  ! [input]
      double precision erg(nb)  ! [input]
      double precision pb       ! [output]

      logical isdgn, wfat_check_invec
      integer io0, io
      double precision pb0
      external wfat_check_invec

      io0 = dgnid(1)
      pb = 0.0d0
      do io = 1, nb
         isdgn = wfat_check_invec(io, ndg, dgnid)
         if (.not. isdgn) then
            pb0 = dp1(io)*dp2(io) / (erg(io0)-erg(io))
            pb = pb + pb0
         endif
      enddo
      pb = -2 * pb

      end
c===============================================================
