c=======================================================
c=======================================================
      subroutine wfat_input(rtdb)

      implicit none
#include "stdio.fh"      
#include "global.fh"      
      integer rtdb



      if (ga_nodeid() .ne. 0) return      

      write(luout, '(/)')
      call util_print_centered(
     &     luout, '-------------------------------------------------',
     &     40, .false.)
      call util_print_centered(
     &     luout, '<<<<<< Begin reading inputs in WFAT block >>>>>>', 
     &     40, .true.)
      
      
c     ==== Read WFAT inputs from input file ====
      call wfat_read_input(rtdb)



      write(luout, '(/)')
      call util_print_centered(
     &     luout, '-------------------------------------------------',
     &     40, .false.)
      call util_print_centered(
     &     luout, '<<<<<< Finish reading inputs in WFAT block >>>>>>',
     &     40, .true.)
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_read_input(rtdb)

      implicit none
#include "errquit.fh"
#include "wfat_input.fh"
#include "inp.fh"      
#include "stdio.fh"
#include "rtdb.fh"      
      integer rtdb
      
      logical tt
      character*50 field
      character*255 ss
      integer i1, i2, k

      logical wfat_inp_paraqn, wfat_inp_movec, wfat_inp_omgpos
      logical wfat_inp_field, wfat_inp_perturb, wfat_inp_dt
      logical wfat_inp_tmax
      external wfat_inp_paraqn, wfat_inp_movec, wfat_inp_omgpos
      external wfat_inp_field, wfat_inp_perturb, wfat_inp_dt
      external wfat_inp_tmax

      
      WFAT_PRINT_OPTIONS = (/ 'srfc', 'rate', 'yield', 'avbeta', 
     &                        'avgamma', 'tmd', 'tmd_atom' /)


c     ==== Initial status of WFAT inputs ====
      print_exist = .false.
      print_stat = .false.

      nfield = 0
      field_exist = .false.
      field_stat = .false.

      dt_exist = .false.
      dt_stat = .false.
      
      tmax_exist = .false.
      tmax_stat = .false.
      
      beta_exist = .false.
      beta_stat = .false.

      gamma_exist = .false.
      gamma_stat = .false.
      
      orbid_exist = .false.
      orbid_stat = .false.

      paraqn_exist = .false.
      paraqn_stat = .false.
      
      nktmd_exist = .false.
      nktmd_stat = .false.

      beta_tmd_exist = .false.
      beta_tmd_stat = .false.

      gamma_tmd_exist = .false.
      gamma_tmd_stat = .false.

      lmax_exist = .false.
      lmax_stat = .false.

      gqual_exist = .false.
      gqual_stat = .false.

      omgpos_exist = .false.
      omgpos_stat = .false.

      sitype_exist = .false.
      sitype_stat = .false.

      dgnthr_exist = .false.
      dgnthr_stat =  .false.
      
      intgl_exist = .false.
      intgl_stat = .false.

      movec_exist = .false.
      movec_stat = .false.
      movec1_exist = .false.
      movec1_stat = .false.

      wfat1_exist = .false.
      wfat1_stat = .false.

      perturb_exist = .false.
      perturb_stat = .false.
      
      prtb_type_exist = .false.
      prtb_type_stat = .false.
      
      prtb_ord_exist = .false.
      prtb_ord_stat = .false.
      
      prtb_bt_exist = .false.
      prtb_bt_stat = .false.
      
      prtb_gm_exist = .false.
      prtb_gm_stat = .false.      

      prtb_dgnthr_exist = .false.
      prtb_dgnthr_stat = .false. 
      
      prtb_offthr_exist = .false.
      prtb_offthr_stat = .false. 
      
      prtb_e1thr_exist = .false. 
      prtb_e1thr_stat = .false.
      
      prtb_itrthr_exist = .false.
      prtb_itrthr_stat = .false. 
      
      prtb_maxitr_exist = .false.
      prtb_maxitr_stat = .false. 
      
      prtb_ordthr_exist = .false.
      prtb_ordthr_stat = .false. 
      
      prtb_maxord_exist = .false.
      prtb_maxord_stat = .false.

      print_srfc = .false.
      print_rate = .false.
      print_yild = .false.
      print_avbeta = .false.
      print_avgamma = .false.
      print_tmd = .false.
      print_tmd_atom = .false.


c     ========================================================
c     ==================== INPUTS PARSING ====================
c     ========================================================

      do while (inp_read()) 
         if (inp_a(field)) then

            write(luout, '("input field : ", a)') field

c           ==== ABSOLUTELY MANDATORY INPUTS ====
            if (inp_compare(.false., 'print', field)) then
               print_exist = .true.
               
               if (print_exist)
     &              write(luout, '("Found print directive")')
               
               if (.not. inp_line(ss))
     &              call errquit('WFAT: Error in reading the '//
     &              '''print'' line.', 1, INPUT_ERR)
               i1 = 0       ! i1 must be initiated to zero before first call of inp_strtok
               tt = inp_strtok(ss, ' ', i1, i2) !! The 'print' keyword
               
               do while (inp_strtok(ss, ' ', i1, i2))        !! The rest of 'print' line
                  tt = inp_match(WFAT_NPRINT, .false.,
     &                 ss(i1:i2), WFAT_PRINT_OPTIONS, k)
                  
                  if (WFAT_PRINT_OPTIONS(k) == 'srfc') then
                     print_srfc = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'rate') then
                     print_rate = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'yield') then
                     print_yild = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'avbeta') then
                     print_avbeta = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'avgamma') then
                     print_avgamma = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'tmd') then
                     print_tmd = .true.
                  else if (WFAT_PRINT_OPTIONS(k) == 'tmd_atom') then
                     print_tmd_atom = .true.
                  else
                     print_stat = .false.
                     exit
                  endif
                  print_stat = .true.
               enddo

c           ==== CONDITIONALLY MANDATORY INPUTS ====    
cOLD            else if (inp_compare(.false., 'fieldmag', field)) then
cOLD               fieldmag_exist = .true.
cOLD               fieldmag_stat = inp_f(fieldmag)

            else if (inp_compare(.false., 'field', field)) then
               nfield = nfield + 1
               field_exist = .true.
               field_stat = wfat_inp_field(nfield, fld_tp, 
     &              fld_fr, fld_fr0, fld_fr_un, 
     &              fld_ce, fld_ce0, fld_ce_un, 
     &              fld_wi, fld_wi0, fld_wi_un,
     &              fld_mx, fld_mx0, fld_mx_un,
     &              fld_cp, fld_cp0, fld_cp_un)
               if (is_extfield) then
                  dt_exist = .true.
                  dt_stat = .true.
                  tmax_exist = .true.
                  tmax_stat = .true.
               endif

            else if (inp_compare(.false., 'dt', field)) then
               dt_exist = .true.
               dt_stat = wfat_inp_dt(wfat_dt, wfat_dt0, dt_un)

            else if (inp_compare(.false., 'tmax', field)) then
               tmax_exist = .true.
               tmax_stat = wfat_inp_tmax(wfat_tmax, wfat_tmax0, tmax_un)
               
            else if (inp_compare(.false., 'beta', field)) then
               beta_exist = .true.
               beta_start_stat = inp_f(wfat_beta_start)
               beta_n_stat = inp_i(wfat_nbeta)
               if (.not. beta_n_stat .or. wfat_nbeta==1) then
                  wfat_nbeta = 1
                  wfat_beta_last = wfat_beta_start
                  beta_last_stat = .true.
               else
                  beta_last_stat = inp_f(wfat_beta_last)
               endif
               beta_stat = beta_start_stat 
            else if (inp_compare(.false., 'gamma', field)) then
               gamma_exist = .true.
               gamma_start_stat = inp_f(wfat_gamma_start)
               gamma_n_stat = inp_i(wfat_ngamma)
               if (.not. gamma_n_stat .or. wfat_ngamma==1) then
                  wfat_ngamma = 1
                  wfat_gamma_last = wfat_gamma_start
                  gamma_last_stat = .true.
               else
                  gamma_last_stat = inp_f(wfat_gamma_last)
               endif
               gamma_stat = gamma_start_stat
            else if (inp_compare(.false., 'beta_tmd', field)) then
               beta_tmd_exist = .true.
               beta_tmd_stat = inp_f(beta_tmd)
            else if (inp_compare(.false., 'gamma_tmd', field)) then
               gamma_tmd_exist = .true.
               gamma_tmd_stat = inp_f(gamma_tmd)
               
c           ==== CONDITIONALLY OPTIONAL INPUTS ====
            else if (inp_compare(.false., 'nktmd', field)) then
               nktmd_exist = .true.
               nktmd_stat = inp_i(nktmd)
               
c           ==== OPTIONAL INPUTS ====
            else if (inp_compare(.false., 'movec_id', field)) then
               orbid_exist = .true.
               orbid_stat = wfat_inp_movec(MAX_N_ORBID, wfat_orbid0,
     &                      wfat_spinid, n_orbid)
            else if(inp_compare(.false., 'paraqn', field)) then
               paraqn_exist = .true.
               paraqn_stat = wfat_inp_paraqn(MAX_N_NPARA, npara, mpara,
     &                       n_npara)
            else if (inp_compare(.false., 'lmax', field)) then
               lmax_exist = .true.
               lmax_stat = inp_i(wfat_lmax)
            else if (inp_compare(.false., 'gridqual', field)) then
               gqual_exist = .true.
               gqual_stat = inp_i(wfat_gqual)
            else if (inp_compare(.false., 'omgpos', field)) then
               omgpos_exist = .true.
               omgpos_stat = wfat_inp_omgpos(autopos, omgpos)
            else if (inp_compare(.false., 'si_type', field)) then
               sitype_exist = .true.
               sitype_stat = inp_a(si_type)
            else if (inp_compare(.false., 'degenthr', field)) then
               dgnthr_exist = .true.
               dgnthr_stat = inp_f(dgnthr)
            else if (inp_compare(.false., 'intgl', field)) then
               intgl_exist = .true.
               intgl_stat = inp_a(wfat_intgl_name)
            else if (inp_compare(.false., 'vectors', field)) then
               movec_exist = .true.
               movec_stat = inp_a(wfat_movec_name)
               if (inp_a(wfat_movec1_name)) then
                  movec1_exist = .true.
                  movec1_stat = .true.
               endif
               call inp_prev_field()
            else if (inp_compare(.false., 'wfat1', field)) then
               wfat1_exist = .true.
               wfat1_stat = .true.
               wfat1 = .true.
            else if (inp_compare(.false., 'perturb', field)) then
               perturb_exist = .true.
               perturb_stat = wfat_inp_perturb(
     &              prtb_type,   prtb_type_exist,   prtb_type_stat,
     &              prtb_ord,    prtb_ord_exist,    prtb_ord_stat,
     &              prtb_bt,     prtb_bt_exist,     prtb_bt_stat,
     &              prtb_gm,     prtb_gm_exist,     prtb_gm_stat,
     &              prtb_dgnthr, prtb_dgnthr_exist, prtb_dgnthr_stat,
     &              prtb_offthr, prtb_offthr_exist, prtb_offthr_stat,
     &              prtb_e1thr,  prtb_e1thr_exist,  prtb_e1thr_stat,
     &              prtb_itrthr, prtb_itrthr_exist, prtb_itrthr_stat,
     &              prtb_maxitr, prtb_maxitr_exist, prtb_maxitr_stat,
     &              prtb_ordthr, prtb_ordthr_exist, prtb_ordthr_stat,
     &              prtb_maxord, prtb_maxord_exist, prtb_maxord_stat)

c           ==== End ====
            else if (inp_compare(.false., 'end', field)) then
               exit
            else
               call errquit('WFAT: The specified keyword '''
     &              //trim(adjustl(field))//''' is unavailable.', 1,
     &              INPUT_ERR)
            endif
         else
            call errquit('WFAT: The first keyword of an input line '//
     &           'must be a valid character type.',
     &           1, INPUT_ERR)
         endif
      enddo


c     ==== Obtain the theory of previous task ====
      if (rtdb_cget(rtdb,'scf:scftype', 1, wfat_theory)) then
         write(luout, '("Previous SCF type : ", a)') trim(wfat_theory)
      else
         call errquit('task: no task input for theory?',0, RTDB_ERR)
      endif
      if (.not. rtdb_cput(rtdb, 'wfat:theory', 1, trim(wfat_theory)))
     &     call errquit('WFAT: Error writing theory to rtdb',
     &     1, RTDB_ERR)
      

c     ==== Check whether the current theory is DFT ====
      if (trim(adjustl(wfat_theory)) /= 'RHF' .and.
     &    trim(adjustl(wfat_theory)) /= 'UHF' ) then
         write(luout, '("WFAT: WARNING!")')
         write(luout, 100)
      endif
 100  format('      WFAT input directive is found while the SCF type ',
     &       'is not RHF. At the moment WFAT can only be used in ',
     &       'conjunction with single determinant methods to obtain ',
     &       'quantities such as molecular orbitals and effective ',
     &       'potential. All WFAT related inputs will thus be ignored.')
      
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_field(i, tp, fr, fr0, fr_un,
     &     ce, ce0, ce_un, wi, wi0, wi_un, mx, mx0, mx_un,
     &     cp, cp0, cp_un)

c     field
c        type < dc | gaussian | hann | external <"laser_file".laser> >
c        frequency <double> [ au | ev | thz | nm | bohr ]
c        center <double> [ fs | au ]
c        width <double> [ fs | au ]
c        max <double> [ au | tw/cm2 ]
c        cep <double> [ rad | deg ]
c     end
      implicit none
#include "wfat_input.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      integer i                                        ! [input],
      character(len=FLD_TP_LEN) tp(max_nfield)         ! [output],
      double precision fr(max_nfield)                  ! [output],
      double precision fr0(max_nfield)                 ! [output],
      character(len=FLD_FR_UN_LEN) fr_un(max_nfield)   ! [output],
      double precision ce(max_nfield)                  ! [output],
      double precision ce0(max_nfield)                 ! [output],
      character(len=FLD_CE_UN_LEN) ce_un(max_nfield)   ! [output],
      double precision wi(max_nfield)                  ! [output],
      double precision wi0(max_nfield)                 ! [output],
      character(len=FLD_WI_UN_LEN) wi_un(max_nfield)   ! [output],
      double precision mx(max_nfield)                  ! [output],
      double precision mx0(max_nfield)                 ! [output],
      character(len=FLD_MX_UN_LEN) mx_un(max_nfield)   ! [output],
      double precision cp(max_nfield)                  ! [output],
      double precision cp0(max_nfield)                 ! [output],
      character(len=FLD_CP_UN_LEN) cp_un(max_nfield)   ! [output],

      character*50 field
      character*2 max_str
      logical tp_stat, fr_stat, ce_stat, wi_stat, mx_stat, cp_stat
      logical tp_exist, fr_exist, ce_exist, wi_exist, mx_exist, cp_exist
      integer k

      logical wfat_inp_externfield
      double precision wfat_convert_freq, wfat_convert_time,
     &     wfat_convert_efield, wfat_convert_angle
      external wfat_convert_freq, wfat_convert_time,
     &     wfat_convert_efield, wfat_convert_angle, wfat_inp_externfield


      wfat_inp_field = .false.
      is_extfield = .false.
      tp_exist = .false.
      fr_exist = .false.
      ce_exist = .false.
      wi_exist = .false.
      mx_exist = .false.
      cp_exist = .false.


      TP_OPTS = (/ 'dc', 'gaussian', 'hann', 'external' /)
      FR_UNIT = (/ 'au', 'ev', 'thz', 'nm', 'bohr' /)
      CE_UNIT = (/ 'fs', 'au' /)
      WI_UNIT = (/ 'fs', 'au' /)
      MX_UNIT = (/ 'au', 'tw/cm2' /)
      CP_UNIT = (/ 'rad', 'deg' /)
      

      write(max_str, '(i2)') max_nfield
      if (i > max_nfield)
     &     call errquit('WFAT: No. of fields exceeds the maximum, '//
     &     'which is '//trim(adjustl(max_str)), 1, INPUT_ERR)
     

      do while (inp_read())
         if (inp_a(field)) then

            if (inp_compare(.false., 'type', field)) then
               tp_exist = .true.
               tp_stat = inp_a(tp(i))
               if (.not. inp_match(TP_NOPTS, .false., tp(i), TP_OPTS,
     &              k)) call errquit('WFAT: The choice of field type '//
     &              'is unknown', 1, INPUT_ERR)
               if (tp(i) == 'external') then
                  if (i > 1) call errquit('WFAT: When there is a '//
     &                 'field with the type of ''external'', that '//
     &                 'field is the only field that can be specified'//
     &                 ' in the input file.', 1, INPUT_ERR)
                  is_extfield = wfat_inp_externfield()
               endif
               
            else if (inp_compare(.false., 'frequency', field)) then
               fr_exist = .true.
               fr_stat = inp_f(fr0(i))
               fr_un(i) = 'au'
               if (inp_a(field)) fr_un(i) = field
               fr(i) = wfat_convert_freq(fr_un(i), fr0(i), 'au')
               
            else if (inp_compare(.false., 'center', field)) then
               ce_exist = .true.
               ce_stat = inp_f(ce0(i))
               ce_un(i) = 'au'
               if (inp_a(field)) ce_un(i) = field
               ce(i) = wfat_convert_time(ce_un(i), ce0(i), 'au')
               
            else if (inp_compare(.false., 'width', field)) then
               wi_exist = .true.
               wi_stat = inp_f(wi0(i))
               wi_un(i) = 'au'
               if (inp_a(field)) wi_un(i) = field
               wi(i) = wfat_convert_time(wi_un(i), wi0(i), 'au')
               
            else if (inp_compare(.false., 'max', field)) then
               mx_exist = .true.
               mx_stat = inp_f(mx0(i))
               mx_un(i) = 'au'
               if (inp_a(field)) mx_un(i) = field
               mx(i) = wfat_convert_efield(mx_un(i), mx0(i), 'au')
               
            else if (inp_compare(.false., 'cep', field)) then
               cp_exist = .true.
               cp_stat = inp_f(cp0(i))
               cp_un(i) = 'rad'
               if (inp_a(field)) cp_un(i) = field
               cp(i) = wfat_convert_angle(cp_un(i), cp0(i), 'rad')

            else if (inp_compare(.false., 'end', field)) then
               exit
            else
               call errquit('WFAT: Unknown input directive '''//
     &              trim(adjustl(field))//''' is found under the '//
     &              '''field'' subdirective', 1, INPUT_ERR)
            endif

         endif
      enddo
      

c     ==== Checking unit choice availability ====
      if (.not. tp_exist) then
         call errquit('WFAT: The ''field'' subdirective '//
     &        '''type'' is missing', 1, INPUT_ERR)
      endif
      
      if (tp(i) /= 'external') then
         if (mx_exist) then
            if (.not. inp_match(MX_NUNIT, .false., mx_un(i),
     &           MX_UNIT, k))
     &           call errquit('WFAT: The unit choice for ''field'' '//
     &           'subdirective ''max'' is unavailable', 1, INPUT_ERR)
         else
            call errquit('WFAT: The ''field'' subdirective '//
     &           '''max'' is missing', 1, INPUT_ERR)
         endif
      endif

      if (tp(i) /= 'dc' .and. tp(i) /= 'external') then
         if (fr_exist) then
            if (.not. inp_match(FR_NUNIT, .false., fr_un(i), FR_UNIT,
     &           k)) call errquit('WFAT: The unit choice for ''field'' '
     &           //'subdirective ''frequency'' is unavailable',
     &           1, INPUT_ERR)
         else
            call errquit('WFAT: The ''field'' subdirective '//
     &           '''frequency'' is missing when it is necessary',
     &           1, INPUT_ERR)
         endif

         if (ce_exist) then
            if (.not. inp_match(CE_NUNIT, .false., ce_un(i), CE_UNIT, 
     &           k)) call errquit('WFAT: The unit choice for ''field'' '
     &           //'subdirective ''center'' is unavailable',
     &           1, INPUT_ERR)
         else
            call errquit('WFAT: The ''field'' subdirective '//
     &           '''center'' is missing when it is necessary',
     &           1, INPUT_ERR)
         endif

         if (wi_exist) then
            if (.not. inp_match(WI_NUNIT, .false., wi_un(i), WI_UNIT, 
     &           k)) call errquit('WFAT: The unit choice for ''field'' '
     &           //'subdirective ''width'' is unavailable',
     &           1, INPUT_ERR)
         else
            call errquit('WFAT: The ''field'' subdirective '//
     &           '''width'' is missing when it is necessary',
     &           1, INPUT_ERR)
         endif

         if (cp_exist) then
            if (.not. inp_match(CP_NUNIT, .false., cp_un(i), CP_UNIT, 
     &           k)) call errquit('WFAT: The unit choice for ''field'' '
     &           //'subdirective ''cep'' is unavailable', 1, INPUT_ERR)
         else
            call errquit('WFAT: The ''field'' subdirective '//
     &           '''cep'' is missing when it is necessary',
     &           1, INPUT_ERR)
         endif
      endif
      

c     ==== When field type is 'dc' ====
      if (tp(i) == 'dc' .or. tp(i) == 'external') then
         fr0(i) = 0.0d0
         ce0(i) = 0.0d0
         wi0(i) = 0.0d0
         cp0(i) = 0.0d0
         fr(i) = 0.0d0
         ce(i) = 0.0d0
         wi(i) = 0.0d0
         cp(i) = 0.0d0
      endif
      if (tp(i) == 'external') then
         mx(i) = 0.0d0
         mx0(i) = 0.0d0
      endif
      

      if (tp(i) == 'dc') then
         wfat_inp_field = tp_stat .and. mx_stat
      else if (tp(i)=='gaussian' .or. tp(i)=='hann') then
         wfat_inp_field = tp_stat .and. fr_stat .and. ce_stat .and.
     &                    wi_stat .and. mx_stat .and. cp_stat
      else if (tp(i) == 'external') then
         wfat_inp_field = is_extfield
      else
         wfat_inp_field = .false.
      endif
      
      end
c=======================================================

      
c=======================================================
c=======================================================
      logical function wfat_inp_externfield()
c     <string t_unit>
c     <string f_unit>
c     <integer data_start>  <integer data_end>
c     data:
c     <double t(1)>     <double f(1)>
c     <double t(2)>     <double f(2)>
c     ...               ...          
     
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "inp.fh"
#include "wfat_input.fh"
      character(len=EXT_TM_UN_LEN) t_unit
      character(len=EXT_FL_UN_LEN) f_unit
      character*5 ss
      integer k, h, i, j, ios
      double precision time, field

      double precision wfat_convert_time, wfat_convert_efield
      external wfat_convert_time, wfat_convert_efield


      wfat_inp_externfield = .false.
      TMEXT_UNIT = (/ 'fs', 'au' /)
      FLEXT_UNIT = (/ 'au', 'tw/cm2' /)
      fld_t_ext(:) = 0.0d0
      fld_f_ext(:) = 0.0d0

      
      if (inp_a(extfield_name)) then
         open(h, file=extfield_name, action='read')

c        ==== Time data unit ====
         read(h, '(a)') t_unit
         if (.not. inp_match(TMEXT_NUNIT, .false., t_unit, TMEXT_UNIT, 
     &        k)) call errquit('WFAT: The unit choice for time data of '
     &        //'the external field is unavailable', 1, INPUT_ERR)

c        ==== Field data unit ====
         read(h, '(a)') f_unit
         if (.not. inp_match(FLEXT_NUNIT, .false., f_unit, FLEXT_UNIT,
     &        k))
     &        call errquit('WFAT: The unit choice for field data of '//
     &        'the external field is unavailable', 1, INPUT_ERR)

c        ==== Data range ====
         read(h, *) ext_start, ext_end
         ntime_ext = ext_end - ext_start + 1
         if (ntime_ext > MAX_EXTERNFIELD) then
            write(luout, 100) ntime_ext, MAX_EXTERNFIELD
            call errquit('WFAT: Error inside the external field file.',
     &           1, INPUT_ERR)
         endif

c        ==== Field and time data ====
         read(h, '(a)') ss
         if (ss == 'data:') then
            j = 0
            do i = 1, MAX_EXTERNFIELD
               if (i >= ext_start .and. i <= ext_end) then
                  read(h, *, iostat=ios) time, field
                  if (ios /= 0) call errquit('WFAT: Bad data within '//
     &                 'the ''data:'' scope in the file containing'//
     &                 'the external field.', 1, INPUT_ERR)
                  j = j + 1
                  fld_t_ext(j) = wfat_convert_time(t_unit, time, 'au')
                  fld_f_ext(j) = wfat_convert_efield(f_unit, field,
     &                 'au')
               else if (i < ext_start) then
                  read(h, *)
               else if (i > ext_end) then
                  exit
               endif
            enddo
         else
            call errquit('WFAT: The keyword ''data:'' is not found '//
     &           'in the external field file.', 1, INPUT_ERR)
         endif
         
         close(h)

c        ==== dt and tmax ====
         dt_ext = fld_t_ext(2) - fld_t_ext(1)
         write(luout, '("dt_ext = ", es14.6)') dt_ext
         tmax_ext = fld_t_ext(ntime_ext) - fld_t_ext(1)
      else
         call errquit('WFAT: Filename of the file containing the '//
     &        'external field is incorrectly or not specified.',
     &        1, INPUT_ERR)
      endif


      write(luout, '(/, 3x, "*** ATTENTION ***")')
      write(luout, 101)

      wfat_inp_externfield = .true.

      
 100  format('ERROR: The no. of data points specified for the ',
     &     'external field, which is ', i6, /,
     &     7x, ', exceeds the maximum, which is ', i6, '.')

 101  format(
     &     5x, 'The field has been taken from an external file.',
     &     'In this case, the user is ', /,
     &     5x, 'responsible in ensuring that the time data values ',
     &     'are equally spaced.')
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_dt(dt, dt0, un)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "wfat_input.fh"
      double precision dt, dt0
      character(len=DT_UN_LEN) un

      integer k
      double precision wfat_convert_time
      external wfat_convert_time

      
      DT_UNIT = (/ 'au', 'fs' /)

      wfat_inp_dt = inp_f(dt0)
      if (inp_a(un)) then
         dt = wfat_convert_time(un, dt0, 'au')
         if (.not. inp_match(DT_NUNIT, .false., un, DT_UNIT, k)) 
     &        call errquit('WFAT: The unit choice for the directive '//
     &        '''dt'' is unavailable', 1, INPUT_ERR)
      else
         dt = dt0
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_tmax(tmax, tmax0, un)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "wfat_input.fh"
      double precision tmax, tmax0
      character(len=TMAX_UN_LEN) un

      integer k
      double precision wfat_convert_time
      external wfat_convert_time

      
      TMAX_UNIT = (/ 'au', 'fs' /)

      wfat_inp_tmax = inp_f(tmax0)
      if (inp_a(un)) then
         tmax = wfat_convert_time(un, tmax0, 'au')
         if (.not. inp_match(TMAX_NUNIT, .false., un, TMAX_UNIT, k)) 
     &        call errquit('WFAT: The unit choice for the directive '//
     &        '''tmax'' is unavailable', 1, INPUT_ERR)
      else
         tmax = tmax0
      endif
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_movec(nelm, orbid0, spinid0, nor0)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "wfat_input.fh"
#include "inp.fh"      
      integer nelm
      integer orbid0(nelm)
      integer spinid0(nelm)
      integer nor0

      character*4 nelm_str
      character*1 afield
      integer i, j, ior, ior1, ior2, is, io
      integer ifield
      integer orbid0_(nelm), spinid0_(nelm)

      logical wfat_check_invec
      external wfat_check_invec



      wfat_inp_movec = .false.

      ior = 1
      is = 0
c     ==== Loop for reading spins ====
      do i = 1, 2
         
         if (.not. inp_a(afield)) then
            if (is .ne. 0) then
               exit
            else
               call errquit('WFAT: No molecular orbitals to be '//
     &              'ionized have been specified. Specify at '//
     &              'least one orbital using the '//
     &              'format ''x 0 -1 -2 ...'' where ''x'' is either '//
     &              '''a'' or ''b''. The trailing list of integers '//
     &              'give the spatial orbitals having the x spin. At '//
     &              'the moment, only negative integers '//
     &              '(corresponding to occupied orbitals) are ',
     &              'supported.',1, INPUT_ERR)
            endif
         endif

         is = is + 1
         if (afield == 'a') then
            spinid0_(is) = 1
         else if (afield == 'b') then
            spinid0_(is) = 2
         else
            call errquit('WFAT: The spin specification for '//
     &           'input keyword ''movec_id'' can only take '//
     &           'either ''a'' or ''b''.', 1, INPUT_ERR)
         endif
         
         io = 0
c        ==== Loop for reading m quantum numbers ====
         do j = 1, MAX_MOVEC_FIELDS
            if (inp_i(ifield)) then
               io = io + 1
               orbid0_(io) = ifield
            else
               if (io == 0)
     &              call errquit('WFAT: One spin channel has not '//
     &              'been given any orbitals.', 1, INPUT_ERR)
               exit
            endif
         enddo

         ior1 = ior
         ior2 = ior1 + io - 1
         nor0 = ior2

c     ==== Check duplicate in the spins ====
         if (is > 1) then
            if (wfat_check_invec(spinid0_(is), is-1, 
     &           spinid0_(1:is-1))) 
     &           call errquit('WFAT: Both of the specified spin '//            
     &           'channels cannot be identical.', 1, INPUT_ERR)
         endif
         spinid0(ior1:ior2) = spinid0_(is)

c        ==== Check duplicate in orbital ====
         if (io > 1) then
            if (wfat_check_invec(orbid0_(io), io-1, 
     &           orbid0_(1:io-1)))
     &           call errquit('WFAT: Two of the specified spatial '//
     &           'orbitals that correspond to '//afield//
     &           '-spin are identical. '//
     &           'Each of the spatial orbitals for any given '//
     &           'spin channel must occur only once.',
     &           1, INPUT_ERR)
         endif
         orbid0(ior1:ior2) = orbid0_(1:io)
         

         ior = ior + io
         if (ior > nelm) then
            write(nelm_str, '(i4)')  nelm
            call errquit('WFAT: The number of the specified '//
     &           'molecular orbitals exceeds limit, which '//
     &           'is equal to '//nelm_str//'.', 1, INPUT_ERR)
         endif

      enddo

      wfat_inp_movec = .true.



      end
c=======================================================

      
c=======================================================
c=======================================================
      logical function wfat_inp_paraqn(nelm, n0, m0, nqn0)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "wfat_input.fh"
#include "inp.fh"      
      integer nelm
      integer n0(nelm)
      integer m0(nelm)
      integer nqn0

      character*4 nelm_str, n_str
      character*5 afield
      logical cond1, cond2
      integer i, j, iqn, iqn1, iqn2, i_n, i_m
      integer ifield
      integer n0_(nelm), m0_(nelm)


      logical wfat_check_invec
      external wfat_check_invec



      wfat_inp_paraqn = .false.

      iqn = 1
      i_n = 0
c     ==== Loop for reading n quantum numbers ====
      do i = 1, MAX_PARAQN_NFIELDS
         
         if (.not. inp_a(afield)) then
            if (i_n .ne. 0) then
               exit
            else
               call errquit('WFAT: No n-parabolic quantum '//
     &              'number has been specified. Specify at '//
     &              'least one n-quantum number using the '//
     &              'format ''nx'' where ''x'' is zero or '//
     &              'positive integer.',1, INPUT_ERR)
            endif
         endif

         cond1 = (afield(1:1) == 'n')
         cond2 = (verify(trim(adjustl(afield(2:NMAX_NFIELD))),
     &           '0123456789') == 0)
         if (cond1 .and. cond2) then
            i_n = i_n + 1
            read(afield(2:NMAX_NFIELD), NFIELD_FORM) n0_(i_n)

            i_m = 0
c           ==== Loop for reading m quantum numbers ====
            do j = 1, MAX_PARAQN_MFIELDS
               if (inp_i(ifield)) then
                  i_m = i_m + 1
                  m0_(i_m) = ifield
               else
                  if (i_m == 0)
     &                 call errquit('WFAT: One n-parabolic quantum '//
     &                 'number has no m-parabolic quantum numbers.',
     &                 1, INPUT_ERR)
                  exit
               endif
            enddo

            iqn1 = iqn
            iqn2 = iqn1 + i_m - 1
            nqn0 = iqn2

c     ==== Check duplicate in n ====
            if (i_n > 1) then
               if (wfat_check_invec(n0_(i_n), i_n-1, 
     &              n0_(1:i_n-1)))
     &              call errquit('WFAT: Two of the specified n- '//
     &              'parabolic quantum numbers are identical. '//
     &              'Each of the n-quantum numbers must occur '//
     &              'only once.', 1, INPUT_ERR)
            endif
            n0(iqn1:iqn2) = n0_(i_n)

c     ==== Check duplicate in m ====
            if (i_m > 1) then
               if (wfat_check_invec(m0_(i_m), i_m-1, 
     &              m0_(1:i_m-1))) then
                  write(n_str, '(i4)') n0(i_n)
                  call errquit('WFAT: Two of the specified m-'//
     &                 'parabolic quantum numbers that correspond to '//
     &                 'n='//trim(adjustl(n_str))//' are identical. '//
     &                 'Each of the m-quantum numbers for any given '//
     &                 'n-quantum number must occur only once.',
     &                 1, INPUT_ERR)
               endif
            endif
            m0(iqn1:iqn2) = m0_(1:i_m)

            
            iqn = iqn + i_m
            if (iqn > nelm) then
               write(nelm_str, '(i4)')  nelm
               call errquit('WFAT: The number of the specified '//
     &              'parabolic quantum numbers exceeds limit, which '//
     &              'is equal to '//nelm_str//'.', 1, INPUT_ERR)
            endif
         else
            call errquit('WFAT: Wrong n-parabolic quantum number '//
     &           'input format. The correct format is nx where x is '//
     &           'zero or positive integer.', 1, INPUT_ERR)
         endif

      enddo

      wfat_inp_paraqn = .true.
      
      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_omgpos(au, omgpos)

      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "inp.fh"      
      logical au
      double precision omgpos(3)

      character*4 field
      logical l
      integer i


      wfat_inp_omgpos = .false.


      if (inp_f(omgpos(1))) then
         au = .false.
         do i = 2, 3
            if (.not. inp_f(omgpos(i))) then
               call errquit('WFAT: Input directive ''omgpos'' '//
     &           'expects either a string ''auto'' or a sequence of '//
     &           'three decimal numbers.', 1, INPUT_ERR)
            endif
         enddo
      else
         l = inp_a(field)
         if (inp_compare(.false., 'auto', field)) then
            au = .true.
            omgpos = 0.0d0
         else            
            call errquit('WFAT: Input directive ''omgpos'' expects '//
     &           'either a string ''auto'' or a sequence of three '//
     &           'decimal numbers.', 1, INPUT_ERR)
         endif
      endif
      

      wfat_inp_omgpos = .true.

      end
c=======================================================


c=======================================================
c=======================================================
      logical function wfat_inp_perturb(
     &                            type,   type_exist,   type_stat,
     &                            ord,    ord_exist,    ord_stat,
     &                            bt,     bt_exist,     bt_stat,
     &                            gm,     gm_exist,     gm_stat,
     &                            dgnthr, dgnthr_exist, dgnthr_stat,
     &                            offthr, offthr_exist, offthr_stat, 
     &                            e1thr,  e1thr_exist,  e1thr_stat, 
     &                            itrthr, itrthr_exist, itrthr_stat,  
     &                            maxitr, maxitr_exist, maxitr_stat,  
     &                            ordthr, ordthr_exist, ordthr_stat,
     &                            maxord, maxord_exist, maxord_stat)

      implicit none
#include "inp.fh"
#include "errquit.fh"      
      character*2 type                   ! [output]
      logical type_exist, type_stat      ! [output]
      integer ord                        ! [output]
      logical ord_exist, ord_stat        ! [output]
      double precision bt                ! [output]
      logical bt_exist, bt_stat          ! [output]
      double precision gm                ! [output]
      logical gm_exist, gm_stat          ! [output]
      double precision dgnthr            ! [output]
      logical dgnthr_exist, dgnthr_stat  ! [output]
      double precision offthr            ! [output]
      logical offthr_exist, offthr_stat  ! [output]
      double precision e1thr             ! [output]
      logical e1thr_exist, e1thr_stat    ! [output]
      double precision itrthr            ! [output]
      logical itrthr_exist, itrthr_stat  ! [output]
      integer maxitr                     ! [output]
      logical maxitr_exist, maxitr_stat  ! [output]
      double precision ordthr            ! [output]
      logical ordthr_exist, ordthr_stat  ! [output]
      integer maxord                     ! [output]
      logical maxord_exist, maxord_stat  ! [output]

      
      character*50 field
      integer nfields, i
      

      type_exist = .false.
      type_stat = .false.

      ord_exist = .false.
      ord_stat = .false.

      bt_exist = .false.
      bt_stat = .false.

      gm_exist = .false.
      gm_stat = .false.

      dgnthr_exist = .false.
      dgnthr_stat = .false. 

      offthr_exist = .false.
      offthr_stat = .false. 

      e1thr_exist = .false. 
      e1thr_stat = .false.

      itrthr_exist = .false.
      itrthr_stat = .false. 

      maxitr_exist = .false.
      maxitr_stat = .false. 

      ordthr_exist = .false.
      ordthr_stat = .false. 

      maxord_exist = .false.
      maxord_stat = .false. 


      nfields = inp_n_field()
      do i = 2, nfields
         if (inp_a(field)) then

            if (inp_compare(.false., 'bw', field)) then
               type = 'bw'
               type_exist = .true.
               type_stat = .true.
            else if (inp_compare(.false., 'rs', field)) then
               type = 'rs'
               type_exist = .true.
               type_stat = .true.
            else if (inp_compare(.false., 'order', field)) then
               ord_exist = .true.
               ord_stat = inp_i(ord)
            else if (inp_compare(.false., 'beta', field)) then
               bt_exist = .true.
               bt_stat = inp_f(bt)
            else if (inp_compare(.false., 'gamma', field)) then
               gm_exist = .true.
               gm_stat = inp_f(gm)

            else if (inp_compare(.false., 'dgnthr', field)) then
               dgnthr_exist = .true.
               dgnthr_stat = inp_f(dgnthr)
            else if (inp_compare(.false., 'offthr', field)) then
               offthr_exist = .true.
               offthr_stat = inp_f(offthr)
            else if (inp_compare(.false., 'e1thr', field)) then
               e1thr_exist = .true.
               e1thr_stat = inp_f(e1thr)
            else if (inp_compare(.false., 'itrthr', field)) then
               itrthr_exist = .true.
               itrthr_stat = inp_f(itrthr)
            else if (inp_compare(.false., 'maxitr', field)) then
               maxitr_exist = .true.
               maxitr_stat = inp_i(maxitr)
            else if (inp_compare(.false., 'ordthr', field)) then
               ordthr_exist = .true.
               ordthr_stat = inp_f(ordthr)
            else if (inp_compare(.false., 'maxord', field)) then
               maxord_exist = .true.
               maxord_stat = inp_i(maxord)
            
            else
               call errquit('WFAT: The specified keyword '''
     &              //trim(adjustl(field))//''' inside ''perturb'' '//
     &              'directive is unavailable.', 1, INPUT_ERR)
            endif
cOLD         else
cOLD            call errquit('WFAT: One parameter for directive '//
cOLD     &           '''perturb'' is not preceded by a string.', 
cOLD     &           1, INPUT_ERR)
         endif
      enddo


      if (type /= 'bw' .and. type /= 'rs') then        ! Just for safety.
         type_stat = .false.
      endif

      
      wfat_inp_perturb = .true.

      end
c=======================================================                  
      
      
cc=======================================================
cc=======================================================
c      subroutine wfat_process_input()
c
c      implicit none
c
c
cc     ==== beta ====
c      if (wfat_nbeta > 1) then
c         dbt = (wfat_beta_last-wfat_beta_start) / (wfat_nbeta-1)
c         do i = 1, wfat_nbeta
c            wfat_beta(i) = wfat_beta_start + (i-1)*dbt
c         enddo
c      else
c         wfat_beta(1) = wfat_beta_start
c      endif
c
cc     ==== gamma ====
c      if (wfat_ngamma > 1) then
c         dgm = (wfat_gamma_last-wfat_gamma_start) / (wfat_ngamma-1)
c         do i = 1, wfat_ngamma
c            wfat_gamma(i) = wfat_gamma_start + (i-1)*dgm
c         enddo
c      else
c         wfat_gamma(1) = wfat_gamma_start
c      endif
c
cc     ==== orbid ====
c      
c      
c      
c
c      end
cc=======================================================
     
