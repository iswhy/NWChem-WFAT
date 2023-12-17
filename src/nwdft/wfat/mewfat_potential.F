c=======================================================
c=======================================================
      subroutine mewfat_hf_vcpsi_1pt(geom, bas, nb, sx, sy, n0, n1,
     &     vpt, ocmo0, ocmo1, pcf, rcf, qcf, scf, zi, ccent, smh, rs, 
     &     scr, vcpsi)

      implicit none
#include "global.fh"
#include "stdio.fh"      
#include "errquit.fh"
      integer geom                                    ! [input], geometry handle.
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of bases in bas.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      double precision vpt(3)                         ! [input], evaluation point.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision pcf(n0(sx))                    ! [input], P det. coefficients. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input],
      double precision ccent(3)                       ! [input],
      double precision smh                            ! [input],
      logical rs                                      ! [input],
      double precision scr(1:*)                       ! [input], scratch array. Minimum number of elements should be 20*nb + n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0).
      double precision vcpsi                          ! [output],

      integer nb0, pb0, pb0n, nb1, pb1, pb1n, nb2, pb2, pb2n
      integer nb3, pb3, pb3n, nsc, psc, pscn, pscn_, ndid
      double precision vhfpsi, v1psi

      ndid = ga_nodeid()


c     ==== pb0 = pointer to evaluated bases ====
      nb0 = nb
      pb0 = 1
      pb0n = pb0 + nb0 - 1

c     ==== pb1 = pointer to evaluated 1st derivative bases ====
      nb1 = 3 * nb
      pb1 = pb0n + 1
      pb1n = pb1 + nb1 - 1
      
c     ==== pb2 = pointer to evaluated 2nd derivative bases ====
      nb2 = 6 * nb
      pb2 = pb1n + 1
      pb2n = pb2 + nb2 - 1

c     ==== pb3 = pointer to evaluated 3rd derivative bases ====
      nb3 = 10 * nb
      pb3 = pb2n + 1
      pb3n = pb3 + nb3 - 1

c     ==== psc = pointer to scratch array ====
      nsc = n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0)
      psc = pb3n + 1
      pscn = psc + nsc - 1
            

c     ==== Calculate basis and its derivatives ====
      call errquit('Use wfat_eval_basis_radf instead of '//
     &     'wfat_eval_basis_1pt below', 22, UERR)
      call wfat_eval_basis_1pt(geom, bas, 0, vpt,
     &                         scr(pb0 : pb0n), scr(pb1 : pb1n),
     &                         scr(pb2 : pb2n), scr(pb3 : pb3n))

      
c     ==== Calculate vhf_psi ====
cOLD      rs = .false.
      call mewfat_hf_veepsi_1ptbas(bas, nb, n0, n1, sx, sy, vpt,
     &     scr(pb0:pb0n), ocmo0, ocmo1, pcf, rcf, qcf, scf, 
     &     scr(psc:pscn), rs, 'none', vhfpsi)

      
c     ==== Calculate v1_psi = vnuc_psi + vs_psi ====
      pscn_ = psc + n0(sx) - 1
      call mewfat_v1psi_1ptbas(geom, nb, sx, n0, vpt, scr(pb0:pb0n),
     &     ocmo0, pcf, rcf, zi, ccent, smh, scr(psc:pscn_), v1psi)


c     ==== Calculate vc_psi ====
      vcpsi = vhfpsi + v1psi

cfordebug      if (ndid==0) write(luout, '("    v= ", 4es14.5)')
cfordebug     &     vhfpsi, vhfpsi/rcf, v1psi, v1psi/rcf
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine mewfat_dft_vcpsi_1pt(geom, bas, nb, sx, sy, n0, n1,
     &     vpt, ocmo0, ocmo1, rdm0, pcf, rcf, qcf, scf, zi, ccent, smh,
     &     rs, si_type, scr, vcpsi)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
      integer geom                                    ! [input], geometry handle.
      integer bas                                     ! [input], basis handle.
      integer nb                                      ! [input], no. of bases in bas.
      integer sx                                      ! [input], sx=1 or 2. The ionized spin channel.
      integer sy                                      ! [input], if sx=1 then sy=2 and vice versa. The complementary channel.
      integer n0(2)                                   ! [input], no. of electrons in neutral.
      integer n1(2)                                   ! [input], no. of electrons in cation.
      double precision vpt(3)                         ! [input], evaluation point.
      double precision ocmo0(nb, n0(1)+n0(2))         ! [input], occupied neutral MO coeffs. The first n0(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision ocmo1(nb, n1(1)+n1(2))         ! [input], occupied cation MO coeffs. The first n1(1) columns hold spin-a coeffs., the rest holds spin-b coeffs.
      double precision rdm0(nb,nb,2)                  ! [input], RDM of the neutral (that corresponds to ocmo0).
      double precision pcf(n0(sx))                    ! [input], P det. coefficients. Its elements index neutral MO's.
      double precision rcf                            ! [input], R det. coefficient.
      double precision qcf(n0(sx), n0(sx), n0(sx)-1)  ! [input], Q det. coefficient. 1st and 2nd dimension index neutral MO's, while the 3rd index cation MO's.
      double precision scf(n0(sy), n0(sy))            ! [input], S det. coefficients. Its rows index neutral MO's while its columns index cation MO's.
      double precision zi                             ! [input],
      double precision ccent(3)                       ! [input],
      double precision smh                            ! [input],
      logical rs                                      ! [input],
      character*5 si_type                             ! [input], self-interaction type.
      double precision scr(1:*)                       ! [input], scratch array. Minimum number of elements should be 20*nb + n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0).
      double precision vcpsi                          ! [output],

      character*4 xclvl
      integer nb0, pb0, pb0n, nb1, pb1, pb1n, nb2, pb2, pb2n
      integer nb3, pb3, pb3n, nsc, psc, pscn, pscn_, dbas, nspin
      double precision vhfpsi, vxcpsi, v1psi

      logical xc_chktau
      external xc_chktau
      

c     ==== pb0 = pointer to evaluated bases ====
      nb0 = nb
      pb0 = 1
      pb0n = pb0 + nb0 - 1

c     ==== pb1 = pointer to evaluated 1st derivative bases ====
      nb1 = 3 * nb
      pb1 = pb0n + 1
      pb1n = pb1 + nb1 - 1
      
c     ==== pb2 = pointer to evaluated 2nd derivative bases ====
      nb2 = 6 * nb
      pb2 = pb1n + 1
      pb2n = pb2 + nb2 - 1

c     ==== pb3 = pointer to evaluated 3rd derivative bases ====
      nb3 = 10 * nb
      pb3 = pb2n + 1
      pb3n = pb3 + nb3 - 1

c     ==== psc = pointer to scratch array ====
      nsc = n0(sx) + nb**2 + nb*max(maxval(n0),maxval(n1)) + sum(n1*n0)    change to include rs
      psc = pb3n + 1
      pscn = psc + nsc - 1

      
c     ==== Determine DFT XC level ====
      if (xc_chktau()) then
         call errquit('MEWFAT: At the moment, mGGA functional is not '//
     &                'yet supported in WFAT.', 11, UERR)
      endif
      call wfat_get_xclvl(dbas, xclvl)
      

c     ==== Calculate basis and its derivatives ====
      call errquit('Use wfat_eval_basis_radf instead of '//
     &     'wfat_eval_basis_1pt below', 24, UERR)
      call wfat_eval_basis_1pt(geom, bas, dbas, vpt,
     &                         scr(pb0 : pb0n), scr(pb1 : pb1n),
     &                         scr(pb2 : pb2n), scr(pb3 : pb3n))

      
c     ==== Calculate v_psi ====
      call mewfat_hf_veepsi_1ptbas(bas, nb, n0, n1, sx, sy, vpt,
     &     scr(pb0:pb0n), ocmo0, ocmo1, pcf, rcf, qcf, scf, 
     &     scr(psc:pscn), rs, si_type, vhfpsi)


c     ==== Calculate vxc_psi ====
      pscn_ = psc + n0(sx) - 1
cOLD      call mewfat_vxcpsi_1ptbas_OLD(nb, sx, n0, scr(pb0), scr(pb1), 
cOLD     &     scr(pb2), ocmo0, rdm0, pcf, rcf, smh, xclvl, scr(psc:pscn_),
cOLD     &     vxcpsi)
      if (n0(1) == n0(2)) then
         nspin = 1
      else
         nspin = 2
      endif
      call mewfat_vxcpsi_1ptbas(nb, nspin, sx, n0, scr(pb0), 
     &     scr(pb1), scr(pb2), ocmo0, rdm0, pcf, rcf, smh, xclvl, 
     &     scr(psc:pscn_), vxcpsi)
      
      
c     ==== Calculate v1_psi = vnuc_psi + vs_psi ====
      pscn_ = psc + n0(sx) - 1
      call mewfat_v1psi_1ptbas(geom, nb, sx, n0, vpt, scr(pb0:pb0n),
     &     ocmo0, pcf, rcf, zi, ccent, smh, scr(psc:pscn_), v1psi)


c     ==== Calculate vc_psi ====
      vcpsi = vhfpsi + vxcpsi + v1psi

      
      end
c=======================================================
