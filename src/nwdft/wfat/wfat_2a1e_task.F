c===============================================================
c===============================================================
      logical function task_wfat_2a1e(rtdb)

      implicit none



c     ==============================
c     ==== Initiate WFAT module ====
c     ==============================
      call wfat_2a1e_init(rtdb)


      res_mode = .false.

c     ==== Task title ====
      call ga_sync()
      call wfat_prep_tasktitle('OE-2a1e')



      call wfat_2a1e_prep_system(q1, q2, z1, z2, eps, sphi, m, orbid, 
     &     maxiter, init_ke, nr, nt)
      call wfat_2a1e_prep_cf(nr, nt, l_cf, k_cf, l_eg, k_eg)

      if (.not. ma_alloc_get(MT_DBL, (nr*nt)**2, 'coef', l_cf, k_cf))
     &     call errquit('WFAT: Cannot allocate coef', 1, MA_ERR)
      if (.not. ma_alloc_get(MT_DBL, nr*nt, 'erg', l_eg, k_eg))
     &     call errquit('WFAT: Cannot allocate erg', 1, MA_ERR)


      call wfat_2a1e_target()



      call wfat_2a1e_intgl()
      
      end
c===============================================================
      
