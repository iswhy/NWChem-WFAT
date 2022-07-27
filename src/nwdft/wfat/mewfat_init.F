c=================================================
c=================================================
      subroutine mewfat_init(rtdb)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"      
#include "geom.fh"
#include "bas.fh"      
#include "mewfat_qm_param.fh"
#include "wfat_rs.fh"      
      integer rtdb

      character*5 si_type
      character*20 basis_name, trans_name
      logical use_nwxc, rtdb_stat, status, reorder, isdft
      integer ndid, basnormid, noio
      double precision nel0_(2), nel1_(2), zn, rscent0(3)

      logical int_normalize, nwxc_has_xc, nwxc_has_cam
      external int_normalize, nwxc_has_xc, nwxc_has_cam
      
      
      ndid = ga_nodeid()

c     ==== Check existence, format, and ====
c     ====     necessities of inputs    ====
      if (ndid == 0) then
         rtdb_stat = rtdb_parallel(.false.)
         call wfat_process_input(rtdb)
         rtdb_stat = rtdb_parallel(.true.)
      endif
      call ga_sync()

cdebug      if (.not. rtdb_print(rtdb, .false.))
cdebug     &     call errquit('WFAT: Cannot print RTDB.', 1, RTDB_ERR)
      call ga_sync()


c     ==== Geometry ====
      if (.not. geom_create(geom_hd, 'mewfat_geom'))
     &     call errquit('MEWFAT: Error in creating geometry object.',
     &     1, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom_hd, 'geometry'))
     &     call errquit('MEWFAT: Error in loading geometry from rtdb',
     &     2, GEOM_ERR)
      if (.not. geom_ncent(geom_hd, natm))
     &     call errquit('MEWFAT: Getting natm failed', 51, GEOM_ERR)
c     NOTES: Since there is only one RTDB for a given MEWFAT job, the geometries for
c     the neutral and cation are implicitly assumed to be identical. No checks of
c     consistency are possible, hence the user is responsible for ensuring that the
c     geometries are identical in the input file.


c     ==== Basis ====
      if (.not. bas_create(bas_hd, 'mewfat_bas'))
     &     call errquit('WFAT: Error in creating basis object.',
     &     1, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom_hd, bas_hd, 'ao basis'))
     &     call errquit('WFAT: Error in loading basis from rtdb',
     &     3, BASIS_ERR)
      if(.not. bas_name(bas_hd, basis_name, trans_name))
     &     call errquit('WFAT: Error in getting basis name.',
     &     1, BASIS_ERR)
      if (.not. bas_numbf(bas_hd, nbas))
     &     call errquit('WFAT: Getting nbas failed in ', 53, UERR)

      if (.not. int_normalize(rtdb, bas_hd))
     &     call errquit('WFAT: Error in normalizing contraction '//
     &     'coefficients.', 1, BASIS_ERR)
      if (.not. bas_norm_get(bas_hd, basnormid))
     &     call errquit('WFAT: Error getting basis norm ID',
     &     1, BASIS_ERR)
c     NOTES: Since there is only one RTDB for a given MEWFAT job, the bases for
c     the neutral and cation are implicitly assumed to be identical. No checks of
c     consistency are possible, hence the user is responsible for ensuring that the
c     bases are identical in the input file.
      
      
c     ==== Neutral molecular orbitals ====
      call mewfat_alloc_mo(rtdb, 0, movecs0, title0, basname0, scftype0,
     &           nbf0, ns_ion0, nmo0, me_nocc0, me_nerg0, l_moocc0, 
     &           k_moocc0, l_moerg0, k_moerg0, g_mov0, nel0_)
      

c     ==== Cation molecular orbitals ====
      call mewfat_alloc_mo(rtdb, 1, movecs1, title1, basname1, scftype1,
     &           nbf1, ns_ion1, nmo1, me_nocc1, me_nerg1, l_moocc1, 
     &           k_moocc1, l_moerg1, k_moerg1, g_mov1, nel1_)


c     ==== No. of bases consistency checks ====
      if (nbf0 /= nbf1 .or. nbf0 /= nbas)
     &     call errquit('MEWFAT: No. of bases for neutral and cation '//
     &     'does not match', 1, UERR)
      

c     ==== Determine no. of electrons ====
      if (.not. rtdb_get(rtdb, 'mewfat:spinx', mt_int, 1, spinx))
     &     call errquit('MEWFAT: Failed reading spinx from RTDB.', 1,
     &     RTDB_ERR)
      if (spinx == 1) then
         spiny = 2
      else if (spinx == 2) then
         spiny = 1
      else
cmovetosanitycheck         call errquit('MEWFAT: Unknown value for the argument ''spinx'' '//
cmovetosanitycheck     &        'of mewfat_calc_pcf.', 1, UERR)
      endif
      nel0 = int(nel0_)
      nel1(spinx) = nel0(spinx) - 1
      nel1(spiny) = nel0(spiny)

      write(luout, '("nel:", 2i4, 3x, 2f8.3, 3x, 2i4)')
     &     nel1, nel1_, int(nel1_)
      
      if (nel1(1) /= int(nel1_(1)) .or. nel1(2) /= int(nel1_(2))) then
         write(luout, 700) nel1, int(nel1_)
cforlater         call errquit('MEWFAT: Inconsistency in the no. of electrons '//
cforlater     &        'in the cation.', 1, UERR)
      endif

           
c     ==== Calculate parent charge ====
      if (.not. geom_nuc_charge(geom_hd, zn))
     &     call errquit('WFAT: Cannot get nuclear charge', 1, GEOM_ERR)
      zion = zn - sum(nel0) + 1


c     ==== Calculate charge state and multiplicity ====
      charge0 = zn - sum(nel0)
      charge1 = zn - sum(nel1)
      multip0 = nel0(1) - nel0(2) + 1
      multip1 = nel1(1) - nel1(2) + 1
      
      
c     ==== Molecular integrals ====
      call int_init(rtdb, 1, bas_hd)
      if (.not. rtdb_get(rtdb, 'dft:tol2e', mt_dbl, 1, tol2e))
     &     call errquit('WFAT: Failed reading tol2e from RTDB.', 1,
     &     RTDB_ERR)
      call scf_get_fock_param(rtdb, tol2e)
      call fock_force_direct(rtdb)          ! Not sure why but making this call speed up walltime.
      call schwarz_init(geom_hd, bas_hd)

      
c     ==== Reorder neutral MO's ====
      mo0_swapped = .false.
      if (.not. rtdb_get(rtdb, 'dft:skeleton', mt_log, 1, oskel))
     &     call errquit('WFAT: Failed reading oskel from RTDB.', 1,
     &     RTDB_ERR)
cOLD      if (.not. rtdb_get(rtdb, 'mewfat:reord', mt_log, 1, reorder))
cOLD     &     reorder = .false.
cOLD      if (reorder)
      call mewfat_reorder_neutmo(geom_hd, bas_hd, nbas, nel0, 
     &     spinx, g_mov0, g_mov1, oskel, dbl_mb(k_moerg0), 
     &     dbl_mb(k_moocc0), mo0_swapped)          ! Should only be called after calling int_init(...).
      

c     ==== XC functionals ====
      call nwxc_rtdb_load(rtdb, 'dft', use_nwxc)
      if (.not. use_nwxc)
     &     call errquit('WFAT: NWXC library for the density '//
     &     'functionals must be used when WFAT task is requested.',
     &     1, UERR)


c     ==== DFT self-interaction ====
      isdft = nwxc_has_xc()
      if (isdft) then
         if (.not. rtdb_cget(rtdb, 'mewfat:si', 1, si_type))
     &        call errquit('MEWFAT: Failed reading SI from RTDB.', 1,
     &        RTDB_ERR)
         if (ndid == 0)
     &        write(luout, '(3x, "Self-interaction type : ", a)')
     &        si_type
      endif
      

c     ==== Use of CAM functional ====
      if (nwxc_has_cam()) then
         if (.not. geom_create(rsgeom, rsgeom_name))
     &        call errquit('WFAT: Cannot create rsgeom', 1, GEOM_ERR)
         rscent0 = (/ 0.0d0, 5.0d0, 0.0d0 /)
         if (.not. geom_cart_set(rsgeom, 1, rstag, rscent0, rschg))
     &        call errquit('WFAT: Cannot set rsgeom', 1, GEOM_ERR)
         if (.not. geom_rtdb_store(rtdb, rsgeom, rsgeom_name))
     &        call errquit('WFAT: Cannot store rsgeom', 1, GEOM_ERR)
         if (.not. bas_create(rsbas, rsbas_name))
     &        call errquit('WFAT: Cannot create rsbas', 1, BASIS_ERR)
         if (.not. bas_rtdb_load(rtdb, rsgeom, rsbas, rsbas_name))
     &        call errquit('WFAT: Cannot load rsbas', 1, BASIS_ERR)
      endif


c     ==== Miscellaneous ====
      if (.not. rtdb_get(rtdb, 'dft:noio', mt_int, 1, noio))
     &     call errquit('WFAT: Failed reading noio from RTDB.', 1,
     &     RTDB_ERR)
      iolgc = .true.
      if (noio.eq.1) iolgc = .false.
      ixcfit = 0


c     ==== Printing important information ====
      call mewfat_print_info(rtdb)
      

 700  format(/,
     &     3x, 'The no. of electrons in the cation deduced from ',
     &     'neutral''s calculation and user''s choice ', /,
     &     3x, 'for the ionized electron spin do not match those ',
     &     'from cation''s calculation.', /,
     &     6x, '- from user''s choice : ', 2i4, /,
     &     6x, '- from cation''s calculation : ', 2i4) 
           

cc      if (scftype0 /= 'dft' .or. scftype1 /= 'dft')
cc     &     call errquit('MEWFAT: The SCF type used to obtain the MO''s '//
cc     &     'of both neutral and cation must be ''dft''.', 1, UERR)
      

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_alloc_mo(rtdb, c, movecs, title, basname, 
     &           scftype, nb, ns, nmo, me_nocc, me_nerg, l_oc, k_oc,
     &           l_eg, k_eg, g_mo, nel)
      
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb             ! [input]
      integer c                ! [input]
      character*(*) movecs     ! [output]
      character*(*) title      ! [output]
      character*(*) basname    ! [output]
      character*(*) scftype    ! [output]
      integer nb               ! [output]
      integer ns               ! [output]
      integer nmo(2)           ! [output]
      integer me_nocc          ! [output]
      integer me_nerg          ! [output]
      integer l_oc             ! [output]
      integer k_oc             ! [output]
      integer l_eg             ! [output]
      integer k_eg             ! [output]
      integer g_mo(2)          ! [output]
      double precision nel(2)  ! [output]
      
      character*1 cs
      character*7 cs_nm
      character*2 spin_name
      character*12 mov_name 
      logical status
      integer is, ioc, ieg, ns_

      logical movecs_read_header, movecs_read
      external movecs_read_header, movecs_read
      

c     ==== Determine whether handling neutral or cation ====
      if (c == 0) then
         cs = '0'
         cs_nm = 'neutral'
      else if (c == 1) then
         cs = '1'
         cs_nm = 'cation'
      else
         call errquit('MEWFAT: Unknwon value of argument ''c'' of '//
     &        'the subroutine ''mewfat_alloc_mo''.', 1, UERR)
      endif

      
c     ==== Getting movecs file names ====
      if (.not. rtdb_cget(rtdb, 'mewfat:movec'//cs//'_name',
     &     1, movecs))
     &     call errquit('MEWFAT: Cannot read movec'//cs//'_name '//
     &     'from RTDB', 1, RTDB_ERR)


c     ==== Reading movecs header ====
      status = movecs_read_header(movecs, title, basname, scftype,
     &         nb, ns_, nmo, 2)
      if (.not. status)
     &     call errquit('MEWFAT: Error reading info of '//cs_nm//' '//
     &     'MO''s', 1, UERR)
      if (ns_ == 2) then       
         ns = ns_
      else if (ns_ == 1) then     ! Take care for spin channels in the case restricted MO's.
         ns = 2
         nmo(2) = nmo(1)
      else
         call errquit('MEWFAT: No. of spin channels is not 1 nor 2',
     &        1, UERR)
      endif

      me_nocc = nb * ns
      if (.not. ma_alloc_get(mt_dbl, me_nocc, 'me_moocc'//cs,
     &     l_oc, k_oc))
     &     call errquit('MEWFAT: Allocating occnum failed', 1, MA_ERR)
      me_nerg = nb * ns
      if (.not. ma_alloc_get(mt_dbl, me_nerg, 'me_moerg'//cs,
     &     l_eg, k_eg))
     &     call errquit('MEWFAT: Allocating erg failed', 2, MA_ERR)


c     ==== Creating MO coefficients global arrrays ====
      do is = 1, ns
         if (is == 1) then
            spin_name = '_a'
            mov_name = 'mo_'//trim(cs_nm)//spin_name
            if (.not. ga_create(mt_dbl, nb, nb, mov_name, nb, 0,
     &           g_mo(is)))
     &           call errquit("WFAT: couldnt create "//mov_name,
     &           0, GA_ERR)
         else if (is == 2) then
            spin_name = '_b'
            mov_name = 'mo_'//trim(cs_nm)//spin_name
            if (.not. ga_duplicate(g_mo(1), g_mo(2), mov_name))
     &           call errquit('MEWFAT: Couldnt create '//mov_name,
     &           0, GA_ERR)
         endif
         
         call ga_zero(g_mo(is))
      enddo      


c     ====Compute numbers of electrons ====
      nel = 0.0d0
      do is = 1, ns
         ioc = k_oc + (is-1) * nb
         ieg = k_eg + (is-1) * nb
         if (ns_ == 2 .or. is == 1) then
            if (.not. movecs_read(movecs, is, dbl_mb(ioc),    !! 1)
     $                dbl_mb(ieg), g_mo(is)))
     $           call errquit("WFAT: Failed to read movecs data",
     &           0, UERR)
            nel(is) = sum(dbl_mb(ioc:ioc+nb-1))
         else if (ns_ == 1 .and. is == 2) then
            call ga_copy(g_mo(1), g_mo(2))
            dbl_mb(ieg:ieg+nb-1) = dbl_mb(k_eg:k_eg+nb-1)
            dbl_mb(k_oc:k_oc+nb-1) = dbl_mb(k_oc:k_oc+nb-1) / 2
            dbl_mb(ioc:ioc+nb-1) = dbl_mb(k_oc:k_oc+nb-1)
            nel(1) = nel(1)/2
            nel(2) = nel(1)
         endif
      enddo
      
      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_reorder_neutmo(geom, bas, nb, n0, sx,
     &           g_mo0, g_mo1, oskel, moerg, moocc, swapped)
      
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer geom                  ! [input],  
      integer bas                   ! [input], 
      integer nb                    ! [input],
      integer n0(2)                 ! [input],
      integer sx                    ! [input], 
      integer g_mo0(2)              ! [in/output], 
      integer g_mo1(2)              ! [input], 
      logical oskel                 ! [input], 
      double precision moerg(nb,2)  ! [in/output], 
      double precision moocc(nb,2)  ! [in/output], 
      logical swapped               ! [output], 

      integer is, i, j, nn, nsc, npcf, i_ovl
      integer l_ovl, l_iv, l_sc, l_pc, l_j
      integer k_ovl, k_iv, k_sc, k_pc, k_j
      integer g_ovl(2)
      double precision erg_j, occ_j, max_p, pcf

      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      

c     ==== Calculate overlap between MO's ====
      do is = 1, 2
         g_ovl(is) = ga_create_atom_blocked(geom, bas, 'MO_overlap')
         call ga_zero(g_ovl(is))
      enddo
      call mewfat_calc_moovl(geom, bas, oskel, g_mo0, g_mo1, g_ovl)
      if (.not. ma_alloc_get(mt_dbl, nb*nb*2, 'moov', l_ovl, k_ovl))
     &     call errquit('WFAT: Cannot allocate moov', 1, MA_ERR)
      do is = 1, 2
         i_ovl = k_ovl + (is-1)*nb*nb
         call ga_get(g_ovl(is), 1, nb, 1, nb, dbl_mb(i_ovl), nb)
      enddo
      

c     ==== Calculate P det. coefficients and reorder ====
      nn = n0(sx)-1
      if (.not. ma_alloc_get(mt_int, nn, 'ipvt', l_iv, k_iv)) 
     &     call errquit('MEWFAT: failed to allocate ipvt', 0, MA_ERR)
      nsc = nn*(nn+1)         ! The '+1' is for the 'work' space.
      if (.not. ma_alloc_get(mt_dbl, nsc, 'scratch', l_sc, k_sc))
     &     call errquit('MEWFAT: Cannot allocate mod_ovl.', 1, MA_ERR)
      npcf = n0(sx)
      if (.not. ma_alloc_get(mt_dbl, npcf, 'pcf', l_pc, k_pc))
     &     call errquit('MEWFAT: Cannot allocate pcf.', 29, MA_ERR)
      if (n0(sx) > 1) then
         call mewfat_calc_pcf(sx, n0, nb, dbl_mb(k_ovl), dbl_mb(k_sc),
     &        int_mb(k_iv), dbl_mb(k_pc))
      endif
      j = 1
      max_p = 0.0d0
      do i = 1, n0(sx)
         pcf = dbl_mb(k_pc+i-1)
         if (ga_nodeid()==0)
     &        write(luout, '("pcf: ", i4, es14.6)') i, pcf
         if (abs(pcf) > max_p) then
            max_p = abs(pcf)
            j = i
         endif
      enddo
      if (ga_nodeid()==0) write(luout, 301) j
 301  format(3x, 'Max P det. coefficient occurs for neutral occupied ',
     &     'MO no. ', i4, '.')
      

c     ==== Reordering ====
      if (.not. ma_alloc_get(mt_dbl, nb, 'ocmo_j', l_j, k_j))
     &     call errquit('MEWFAT: Cannot allocate ocmo_j.', 52, MA_ERR)
      swapped = .false.
      if (j < n0(sx)) then
         if (ga_nodeid()==0) write(luout, '("REORDERING . . .")')
c        ==== Swap coefficients ====
         call ga_get(g_mo0(sx), 1, nb, j, j, dbl_mb(k_j), nb)
         call ga_copy_patch('n', g_mo0(sx), 1, nb, n0(sx), n0(sx),
     &                      g_mo0(sx), 1, nb, j, j)
         call ga_put(g_mo0(sx), 1, nb, n0(sx), n0(sx), dbl_mb(k_j), nb)

c        ==== Swap energies ====
         erg_j = moerg(j,sx)
         moerg(j,sx) = moerg(n0(sx),sx)
         moerg(n0(sx),sx) = erg_j

c        ==== Swap occupation no. ====
         occ_j = moocc(j,sx)
         moocc(j,sx) = moocc(n0(sx),sx)
         moocc(n0(sx),sx) = occ_j

         swapped = .true.
      else if (j > n0(sx)) then
         call errquit('MEWFAT: Error in mewfat_reorder_neutmo',
     &        37, UERR)
      endif


c     ==== Deallocations ====
      do is = 1, 2
         if (.not. ga_destroy(g_ovl(is)))
     &        call errquit('MEWFAT: Cannot destroy MO_overlap.',
     &        47, GA_ERR)
      enddo      
      if (.not. ma_free_heap(l_ovl))
     &     call errquit('MEWFAT: Cannot deallocate moov', 81, MA_ERR)
      if (.not. ma_free_heap(l_iv))
     &     call errquit('MEWFAT: Cannot deallocate ipvt', 81, MA_ERR)
      if (.not. ma_free_heap(l_sc))
     &     call errquit('MEWFAT: Cannot deallocate mod_ovl', 81, MA_ERR)
      if (.not. ma_free_heap(l_pc))
     &     call errquit('MEWFAT: Cannot deallocate pcf', 81, MA_ERR)
      if (.not. ma_free_heap(l_j))
     &     call errquit('MEWFAT: Cannot deallocate ocmo_j', 81, MA_ERR)

      end
c=================================================


c=================================================
c=================================================
      subroutine mewfat_terminate(rtdb)

      implicit none
#include "rtdb.fh"      
#include "errquit.fh"
#include "mafdecls.fh"      
#include "mewfat_qm_param.fh"
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
      integer rtdb
      
      character*2 sp
      character*12 fn
      integer is


c     ==== Neutral deallocations ====
      if (.not. ma_free_heap(l_moocc0))
     &     call errquit('WFAT: Deallocating occnum failed', 3, MA_ERR)
      if (.not. ma_free_heap(l_moerg0))
     &     call errquit('WFAT: Deallocating erg failed', 3, MA_ERR)


c     ==== Cation deallocations ====
      if (.not. ma_free_heap(l_moocc1))
     &     call errquit('WFAT: Deallocating occnum failed', 3, MA_ERR)
      if (.not. ma_free_heap(l_moerg1))
     &     call errquit('WFAT: Deallocating erg failed', 3, MA_ERR)
      
      
c     ==== Geometry ====
      if (.not. geom_destroy(geom_hd))
     &     call errquit('MEWFAT: Error in destroying geometry object.',
     &     1, GEOM_ERR)

      
c     ==== Basis ====
      if (.not. bas_destroy(bas_hd))
     &     call errquit('MEWFAT: Error in destroying basis object.',
     &     1, BASIS_ERR)      

      
c     ==== Molecular orbitals ====
      do is = 1, 2
         if (is == 1) then
            sp = '_a'
         else if (is == 2) then
            sp = '_b'
         endif
   
         fn = 'mo_neutral'//sp
         if (.not. ga_destroy(g_mov0(is)))
     $        call errquit ("MEWFAT: couldnt destroy "//fn, 0, GA_ERR)
         fn = 'mo_cation'//sp
         if (.not. ga_destroy(g_mov1(is)))
     $        call errquit ("MEWFAT: couldnt destroy "//fn, 0, GA_ERR)
      enddo

      
c     ==== Molecular integrals ====
      call int_terminate()

      end
c=================================================
      
