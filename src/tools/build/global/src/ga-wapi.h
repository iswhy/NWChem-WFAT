#ifndef WAPI_H_
#define WAPI_H_

#include <stdio.h>

#include "gacommon.h"
#include "typesf2c.h"


/* Routines from base.c */
extern logical wnga_allocate(Integer g_a);
extern logical wnga_compare_distr(Integer g_a, Integer g_b);
extern logical wnga_create(Integer type, Integer ndim,
                           Integer *dims, char* name,
                           Integer *chunk, Integer *g_a);
extern logical wnga_create_config(Integer type, Integer ndim,
                                  Integer *dims, char* name,
                                  Integer *chunk, Integer p_handle,
                                  Integer *g_a);
extern logical wnga_create_ghosts(Integer type, Integer ndim,
                                  Integer *dims, Integer *width, char* name,
                                  Integer *chunk, Integer *g_a);
extern logical wnga_create_ghosts_irreg(Integer type, Integer ndim,
                                        Integer *dims, Integer *width,
                                        char* name,
                                        Integer *map, Integer *block,
                                        Integer *g_a);
extern logical wnga_create_ghosts_irreg_config(Integer type, Integer ndim,
                                               Integer *dims, Integer *width,
                                               char* name,
                                               Integer *map, Integer *block,
                                               Integer p_handle, Integer *g_a);
extern logical wnga_create_ghosts_config(Integer type, Integer ndim,
                                         Integer *dims, Integer *width,
                                         char* name,
                                         Integer *chunk, Integer p_handle,
                                         Integer *g_a);
extern logical wnga_create_irreg(Integer type, Integer ndim,
                                 Integer *dims, char* name,
                                 Integer *map, Integer *block, Integer *g_a);
extern logical wnga_create_irreg_config(Integer type, Integer ndim,
                                        Integer *dims, char* name,
                                        Integer *map,
                                        Integer *block, Integer p_handle,
                                        Integer *g_a);
extern Integer wnga_create_handle();
extern logical wnga_create_mutexes(Integer num);
extern logical wnga_destroy(Integer g_a);
extern logical wnga_destroy_mutexes();
extern void wnga_distribution(Integer g_a, Integer proc, Integer *lo, Integer *hi);
extern logical wnga_duplicate(Integer g_a, Integer *g_b, char *array_name);
extern void wnga_fill(Integer g_a, void* val);
extern void wnga_get_block_info(Integer g_a, Integer *num_blocks,
                                Integer *block_dims);
extern logical wnga_get_debug();
extern Integer wnga_get_dimension(Integer g_a);
extern void wnga_get_proc_grid(Integer g_a, Integer *dims);
extern void wnga_get_proc_index(Integer g_a, Integer iproc, Integer *index);
extern logical wnga_has_ghosts(Integer g_a);
extern void wnga_initialize();
extern int  wnga_initialized();
extern void wnga_initialize_ltd(Integer limit);
extern void wnga_inquire(Integer g_a, Integer *type, Integer *ndim, Integer *dims);
extern void wnga_inquire_type(Integer g_a, Integer *type);
extern Integer wnga_inquire_memory();
extern void wnga_inquire_name(Integer g_a, char **array_name);
extern logical wnga_is_mirrored(Integer g_a);
extern void wnga_list_nodeid(Integer *list, Integer nprocs);
extern logical wnga_locate(Integer g_a, Integer *subscript, Integer *owner);
extern Integer wnga_locate_num_blocks(Integer g_a, Integer *lo, Integer *hi);
extern logical wnga_locate_nnodes(Integer g_a, Integer *lo, Integer *hi, Integer *np);
extern logical wnga_locate_region(Integer g_a, Integer *lo, Integer *hi,
                                  Integer *map, Integer *proclist, Integer *np);
extern void wnga_lock(Integer mutex);
extern Integer wnga_ndim(Integer g_a);
extern void wnga_mask_sync(Integer begin, Integer end);
extern Integer wnga_memory_avail();
extern logical wnga_memory_limited();
extern void wnga_merge_distr_patch(Integer g_a, Integer *alo, Integer *ahi,
                                   Integer g_b, Integer *blo, Integer *bhi);
extern void wnga_merge_mirrored(Integer g_a);
extern void wnga_nblock(Integer g_a, Integer *nblock);

extern Integer wnga_nnodes();
extern Integer wnga_nodeid();
extern Integer wnga_pgroup_absolute_id(Integer grp, Integer pid);
extern Integer wnga_pgroup_create(Integer *list, Integer count);
extern logical wnga_pgroup_destroy(Integer grp);
extern Integer wnga_pgroup_get_default();
extern Integer wnga_pgroup_get_mirror();
extern Integer wnga_pgroup_get_world();
extern void wnga_pgroup_set_default(Integer grp);
extern Integer wnga_pgroup_split(Integer grp, Integer grp_num);
extern Integer wnga_pgroup_split_irreg(Integer grp, Integer mycolor);
extern Integer wnga_pgroup_nnodes(Integer grp);
extern Integer wnga_pgroup_nodeid(Integer grp);
extern void wnga_proc_topology(Integer g_a, Integer proc, Integer* subscript);
extern void wnga_randomize(Integer g_a, void* val);
extern Integer wnga_get_pgroup(Integer g_a);
extern Integer wnga_get_pgroup_size(Integer grp_id);
extern void wnga_set_array_name(Integer g_a, char *array_name);
extern void wnga_set_block_cyclic(Integer g_a, Integer *dims);
extern void wnga_set_block_cyclic_proc_grid(Integer g_a, Integer *dims, Integer *proc_grid);
extern void wnga_set_chunk(Integer g_a, Integer *chunk);
extern void wnga_set_data(Integer g_a, Integer ndim, Integer *dims, Integer type);
extern void wnga_set_debug(logical flag);
extern void wnga_set_ghosts(Integer g_a, Integer *width);
extern void wnga_set_irreg_distr(Integer g_a, Integer *map, Integer *block);
extern void wnga_set_irreg_flag(Integer g_a, logical flag);
extern void wnga_set_memory_limit(Integer mem_limit);
extern void wnga_set_pgroup(Integer g_a, Integer p_handle);
extern void wnga_set_restricted(Integer g_a, Integer *list, Integer size);
extern void wnga_set_restricted_range(Integer g_a, Integer lo_proc, Integer hi_proc);
extern void wnga_terminate();
extern Integer wnga_total_blocks(Integer g_a);
extern void wnga_unlock(Integer mutex);
extern logical wnga_uses_ma();
extern logical wnga_uses_proc_grid(Integer g_a);
extern logical wnga_valid_handle(Integer g_a);
extern Integer wnga_verify_handle(Integer g_a);
extern void wnga_check_handle(Integer g_a, char *string);

/* Routines from onesided.c */
extern void wnga_acc(Integer g_a, Integer *lo, Integer *hi, void *buf,
                     Integer *ld, void *alpha);
extern void wnga_access_idx(Integer g_a, Integer *lo, Integer *hi,
                            AccessIndex *index, Integer *ld);
extern void wnga_access_ptr(Integer g_a, Integer *lo, Integer *hi, void *ptr,
                            Integer *ld);
extern void wnga_access_block_idx(Integer g_a, Integer idx,
                                  AccessIndex* index, Integer *ld);
extern void wnga_access_block_ptr(Integer g_a, Integer idx, void* ptr,
                                  Integer *ld);
extern void wnga_access_block_grid_idx(Integer g_a, Integer* subscript,
                                       AccessIndex *index, Integer *ld);
extern void wnga_access_block_grid_ptr(Integer g_a, Integer *index, void* ptr,
                                       Integer *ld);
extern void wnga_access_block_segment_idx(Integer g_a, Integer proc,
                                          AccessIndex* index, Integer *len);
extern void wnga_access_block_segment_ptr(Integer g_a, Integer proc,
                                          void* ptr, Integer *len);
extern void wnga_alloc_gatscat_buf(Integer nelems);
extern void wnga_fence();
extern void wnga_free_gatscat_buf();
extern void wnga_gather2d(Integer g_a, void *v, Integer *i, Integer *j,
                          Integer nv);
extern void wnga_gather(Integer g_a, void* v, void *subscript,
                        Integer c_flag, Integer nv);
extern void wnga_get(Integer g_a, Integer *lo, Integer *hi,
                     void *buf, Integer *ld);
extern void wnga_init_fence();
extern void wnga_nbacc(Integer g_a, Integer *lo, Integer *hi, void *buf,
                       Integer *ld, void *alpha, Integer *nbhndl);
extern void wnga_nbget(Integer g_a, Integer *lo, Integer *hi, void *buf,
                       Integer *ld, Integer *nbhandle);
extern void wnga_nbput(Integer g_a, Integer *lo, Integer *hi, void *buf,
                       Integer *ld, Integer *nbhandle);
extern void wnga_nbput_notify(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld, Integer g_b, Integer *ecoords, void *bufn, Integer *nbhandle);
extern void wnga_nbwait_notify(Integer *nbhandle);
extern Integer wnga_nbtest(Integer *nbhandle);
extern void wnga_nbwait(Integer *nbhandle);
extern void wnga_put(Integer g_a, Integer *lo, Integer *hi, void *buf,
                     Integer *ld);
extern void wnga_pgroup_sync(Integer grp_id);
extern Integer wnga_read_inc(Integer g_a, Integer *subscript, Integer inc);
extern void wnga_release(Integer g_a, Integer *lo, Integer *hi);
extern void wnga_release_block(Integer g_a, Integer iblock);
extern void wnga_release_block_grid(Integer g_a, Integer *index);
extern void wnga_release_block_segment(Integer g_a, Integer iproc);
extern void wnga_release_update(Integer g_a, Integer *lo, Integer *hi);
extern void wnga_release_update_block(Integer g_a, Integer iblock);
extern void wnga_release_update_block_grid(Integer g_a, Integer *index);
extern void wnga_release_update_block_segment(Integer g_a, Integer iproc);
extern void wnga_scatter2d(Integer g_a, void *v, Integer *i, Integer *j, Integer nv);
extern void wnga_scatter(Integer g_a, void *v, void *subscript,
                         Integer c_flag, Integer nv);
extern void wnga_scatter_acc2d(Integer g_a, void *v, Integer *i, Integer *j,
                               Integer nv, void *alpha);
extern void wnga_scatter_acc(Integer g_a, void* v, void *subscript,
                             Integer c_flag, Integer nv, void *alpha);
extern void wnga_strided_acc(Integer g_a, Integer *lo, Integer *hi, Integer *skip,
                             void *buf, Integer *ld, void *alpha);
extern void wnga_strided_get(Integer g_a, Integer *lo, Integer *hi, Integer *skip,
                             void *buf, Integer *ld);
extern void wnga_strided_put(Integer g_a, Integer *lo, Integer *hi, Integer *skip,
                             void *buf, Integer *ld);
extern void wnga_sync();
extern DoublePrecision wnga_wtime();

/* Routines from datatypes.c */
extern Integer wnga_type_f2c(Integer type);
extern Integer wnga_type_c2f(Integer type);

/* Routines from collect.c */
extern void wnga_msg_brdcst(Integer type, void *buffer, Integer len, Integer root);
extern void wnga_brdcst(Integer type, void *buf, Integer len, Integer originator);
extern void wnga_pgroup_brdcst(Integer grp_id, Integer type, void *buf, Integer len, Integer originator);
extern void wnga_msg_sync();
extern void wnga_msg_pgroup_sync(Integer grp_id);
extern void wnga_pgroup_gop(Integer p_grp, Integer type, void *x, Integer n, char *op);
extern void wnga_gop(Integer type, void *x, Integer n, char *op);

/* Routines from elem_alg.c */
extern void wnga_abs_value_patch(Integer g_a, Integer *lo, Integer *hi);
extern void wnga_recip_patch(Integer g_a, Integer *lo, Integer *hi);
extern void wnga_add_constant_patch(Integer g_a, Integer *lo, Integer *hi, void *alpha);
extern void wnga_abs_value(Integer g_a);
extern void wnga_add_constant(Integer g_a, void *alpha);
extern void wnga_recip(Integer g_a);
extern void wnga_elem_multiply(Integer g_a, Integer g_b, Integer g_c);
extern void wnga_elem_divide(Integer g_a, Integer g_b, Integer g_c);
extern void wnga_elem_maximum(Integer g_a, Integer g_b, Integer g_c);
extern void wnga_elem_minimum(Integer g_a, Integer g_b, Integer g_c);
extern void wnga_elem_multiply_patch(Integer g_a,Integer *alo,Integer *ahi,Integer g_b,Integer *blo,Integer *bhi,Integer g_c,Integer *clo,Integer *chi);
extern void wnga_elem_divide_patch(Integer g_a,Integer *alo,Integer *ahi,Integer g_b,Integer *blo,Integer *bhi,Integer g_c,Integer *clo,Integer *chi);
extern void wnga_elem_maximum_patch(Integer g_a,Integer *alo,Integer *ahi,Integer g_b,Integer *blo,Integer *bhi,Integer g_c,Integer *clo,Integer *chi);
extern void wnga_elem_minimum_patch(Integer g_a,Integer *alo,Integer *ahi,Integer g_b,Integer *blo,Integer *bhi,Integer g_c,Integer *clo,Integer *chi);
extern void wnga_elem_step_divide_patch(Integer g_a,Integer *alo,Integer *ahi, Integer g_b,Integer *blo,Integer *bhi,Integer g_c, Integer *clo,Integer *chi);
extern void wnga_elem_stepb_divide_patch(Integer g_a,Integer *alo,Integer *ahi, Integer g_b,Integer *blo,Integer *bhi,Integer g_c, Integer *clo,Integer *chi);
extern void wnga_step_mask_patch(Integer g_a,Integer *alo,Integer *ahi, Integer g_b,Integer *blo,Integer *bhi,Integer g_c, Integer *clo,Integer *chi);
extern void wnga_step_bound_info_patch(Integer g_xx, Integer *xxlo, Integer *xxhi, Integer g_vv, Integer *vvlo, Integer *vvhi, Integer g_xxll, Integer *xxlllo, Integer *xxllhi, Integer g_xxuu, Integer *xxuulo, Integer *xxuuhi, void *boundmin, void* wolfemin, void *boundmax);
extern void wnga_step_max_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, void *result);
extern void wnga_step_max(Integer g_a, Integer g_b, void *retval);
extern void wnga_step_bound_info(Integer g_xx, Integer g_vv, Integer g_xxll, Integer g_xxuu, void *boundmin, void *wolfemin, void *boundmax);

/* Routines from ga_solve_seq.c */
extern void wnga_lu_solve_seq(char *trans, Integer g_a, Integer g_b);

/* Routines from global.util.c */
extern void wnga_print_stats();
extern void wnga_error(char *string, Integer icode);
extern Integer wnga_cluster_nodeid();
extern Integer wnga_cluster_nprocs(Integer node);
extern Integer wnga_cluster_procid(Integer node, Integer loc_proc_id);
extern Integer wnga_cluster_nnodes();
extern Integer wnga_cluster_proc_nodeid(Integer proc);
extern void wnga_print_file(FILE *file, Integer g_a);
extern void wnga_print(Integer g_a);
extern void wnga_print_patch_file2d(FILE *file, Integer g_a, Integer ilo, Integer ihi, Integer jlo, Integer jhi, Integer pretty);
extern void wnga_print_patch2d(Integer g_a, Integer ilo, Integer ihi, Integer jlo, Integer jhi, Integer pretty);
extern void wnga_print_patch_file(FILE *file, Integer g_a, Integer *lo, Integer *hi, Integer pretty);
extern void wnga_print_patch(Integer g_a, Integer *lo, Integer *hi, Integer pretty);
extern void wnga_print_distribution(int fstyle, Integer g_a);
extern void wnga_summarize(Integer verbose);

/* Routines from ghosts.c */
extern void wnga_access_ghost_ptr(Integer g_a, Integer dims[], void* ptr, Integer ld[]);
extern void wnga_access_ghost_element(Integer g_a, AccessIndex* index, Integer subscript[], Integer ld[]);
extern void wnga_access_ghost_element_ptr(Integer g_a, void *ptr, Integer subscript[], Integer ld[]);
extern void wnga_access_ghosts(Integer g_a, Integer dims[], AccessIndex* index, Integer ld[]);
extern void wnga_release_ghost_element(Integer g_a, Integer subscript[]);
extern void wnga_release_update_ghost_element(Integer g_a, Integer subscript[]);
extern void wnga_release_ghosts(Integer g_a);
extern void wnga_release_update_ghosts(Integer g_a);
extern void wnga_get_ghost_block(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld);
extern void wnga_update1_ghosts(Integer g_a);
extern logical wnga_update2_ghosts(Integer g_a);
extern logical wnga_update3_ghosts(Integer g_a);
extern logical wnga_set_update4_info(Integer g_a);
extern logical wnga_update4_ghosts(Integer g_a);
extern logical wnga_update44_ghosts(Integer g_a);
extern logical wnga_update55_ghosts(Integer g_a);
extern logical wnga_update_ghost_dir(Integer g_a, Integer pdim, Integer pdir, logical pflag);
extern logical wnga_update5_ghosts(Integer g_a);
extern logical wnga_set_update5_info(Integer g_a);
extern void wnga_update_ghosts(Integer g_a);
extern void wnga_update_ghosts_nb(Integer g_a, Integer *nbhandle);
extern logical wnga_update6_ghosts(Integer g_a);
extern logical wnga_update7_ghosts(Integer g_a);
extern void wnga_ghost_barrier();
extern void wnga_nbget_ghost_dir(Integer g_a, Integer *mask, Integer *nbhandle);
extern logical wnga_set_ghost_info(Integer g_a);
extern void wnga_set_ghost_corner_flag(Integer g_a, logical flag);

/* Routines from global.nalg.c */
extern void wnga_zero(Integer g_a);
extern void wnga_copy(Integer g_a, Integer g_b);
extern void wnga_dot(int type, Integer g_a, Integer g_b, void *value);
extern void wnga_scale(Integer g_a, void* alpha);
extern void wnga_add(void *alpha, Integer g_a, void* beta, Integer g_b, Integer g_c);
extern void wnga_transpose(Integer g_a, Integer g_b);

/* Routines from global.npatch.c */
extern void wnga_copy_patch(char *trans, Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi);
extern void wnga_zero_patch(Integer g_a, Integer *lo, Integer *hi);
extern logical wnga_patch_intersect(Integer *lo, Integer *hi, Integer *lop, Integer *hip, Integer ndim);
extern logical wnga_comp_patch(Integer andim, Integer *alo, Integer *ahi, Integer bndim, Integer *blo, Integer *bhi);
extern void wnga_dot_patch(Integer g_a, char *t_a, Integer *alo, Integer *ahi, Integer g_b, char *t_b, Integer *blo, Integer *bhi, void *retval);
extern void wnga_fill_patch(Integer g_a, Integer *lo, Integer *hi, void* val);
extern void wnga_scale_patch(Integer g_a, Integer *lo, Integer *hi, void *alpha);
extern void wnga_add_patch(void *alpha, Integer g_a, Integer *alo, Integer *ahi, void *beta, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi);

/* Routines from select.c */

extern void wnga_select_elem(Integer g_a, char* op, void* val, Integer *subscript);

/* Routines from ga_malloc.c */

extern Integer wnga_memory_avail_type(Integer datatype);

/* Routines from sparse.c */

extern void wnga_patch_enum(Integer g_a, Integer lo, Integer hi, void* start, void* stride);
extern void wnga_scan_copy(Integer g_a, Integer g_b, Integer g_sbit, Integer lo, Integer hi);
extern void wnga_scan_add(Integer g_a, Integer g_b, Integer g_sbit, Integer lo, Integer hi, Integer excl);
extern void wnga_pack(Integer g_a, Integer g_b, Integer g_sbit, Integer lo, Integer hi, Integer* icount);
extern void wnga_unpack(Integer g_a, Integer g_b, Integer g_sbit, Integer lo, Integer hi, Integer* icount);
extern logical wnga_create_bin_range(Integer g_bin, Integer g_cnt, Integer g_off, Integer *g_range);
extern void wnga_bin_sorter(Integer g_bin, Integer g_cnt, Integer g_off);
extern void wnga_bin_index(Integer g_bin, Integer g_cnt, Integer g_off, Integer *values, Integer *subs, Integer n, Integer sortit);

/* Routines from matrix.c */

extern void wnga_median_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi, Integer g_m, Integer *mlo, Integer *mhi);
extern void wnga_median(Integer g_a, Integer g_b, Integer g_c, Integer g_m);
extern void wnga_norm_infinity(Integer g_a, double *nm);
extern void wnga_norm1(Integer g_a, double *nm);
extern void wnga_get_diag(Integer g_a, Integer g_v);
extern void wnga_add_diagonal(Integer g_a, Integer g_v);
extern void wnga_set_diagonal(Integer g_a, Integer g_v);
extern void wnga_shift_diagonal(Integer g_a, void *c);
extern void wnga_zero_diagonal(Integer g_a);
extern void wnga_scale_rows(Integer g_a, Integer g_v);
extern void wnga_scale_cols(Integer g_a, Integer g_v);

/* Routines from ga_symmetr.c */

extern void wnga_symmetrize(Integer g_a);

/* Routines from global.periodic.c */

extern void wnga_periodic(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld, void *alpha, Integer op_code);

/* Routines from matmul.c */

extern void wnga_matmul(char *transa, char *transb, void *alpha, void *beta, Integer g_a, Integer ailo, Integer aihi, Integer ajlo, Integer ajhi, Integer g_b, Integer bilo, Integer bihi, Integer bjlo, Integer bjhi, Integer g_c, Integer cilo, Integer cihi, Integer cjlo, Integer cjhi);
extern void wnga_matmul_mirrored(char *transa, char *transb, void *alpha, void *beta, Integer g_a, Integer ailo, Integer aihi, Integer ajlo, Integer ajhi, Integer g_b, Integer bilo, Integer bihi, Integer bjlo, Integer bjhi, Integer g_c, Integer cilo, Integer cihi, Integer cjlo, Integer cjhi);
extern void wnga_matmul_patch(char *transa, char *transb, void *alpha, void *beta, Integer g_a, Integer alo[], Integer ahi[], Integer g_b, Integer blo[], Integer bhi[], Integer g_c, Integer clo[], Integer chi[]);

/* Routines from ga_diag_seqc.c */

extern void wnga_diag_seq(Integer g_a, Integer g_s, Integer g_v, DoublePrecision *eval);
extern void wnga_diag_std_seq(Integer g_a, Integer g_v, DoublePrecision *eval);

/* Routines from peigstubs.c */

extern void wnga_diag(Integer g_a, Integer g_s, Integer g_v, DoublePrecision *eval);
extern void wnga_diag_std(Integer g_a, Integer g_v, DoublePrecision *eval);
extern void wnga_diag_reuse(Integer reuse, Integer g_a, Integer g_s, Integer g_v, DoublePrecision *eval);

/* Routines from sclstubs.c */

extern void wnga_lu_solve_alt(Integer tran, Integer g_a, Integer g_b);
extern void wnga_lu_solve(char *tran, Integer g_a, Integer g_b);
extern Integer wnga_llt_solve(Integer g_a, Integer g_b);
extern Integer wnga_solve(Integer g_a, Integer g_b);
extern Integer wnga_spd_invert(Integer g_a);

/* Routines from DP.c */

extern void wnga_copy_patch_dp(char *t_a, Integer g_a, Integer ailo, Integer aihi, Integer ajlo, Integer ajhi, Integer g_b, Integer bilo, Integer bihi, Integer bjlo, Integer bjhi);
extern DoublePrecision wnga_ddot_patch_dp(Integer g_a, char *t_a, Integer ailo, Integer aihi, Integer ajlo, Integer ajhi, Integer g_b, char *t_b, Integer bilo, Integer bihi, Integer bjlo, Integer bjhi);

/* Routines from ga_trace.c */

extern double wnga_timer();

/*Routines for types from base.c*/

extern int wnga_register_type(size_t size);
extern int wnga_deregister_type(int type);

/*Routines for field-wise GA operations*/

extern void wnga_get_field(Integer g_a, Integer *lo, Integer *hi, Integer foff, Integer fsize, void *buf, Integer *ld);
extern void wnga_nbget_field(Integer g_a, Integer *lo, Integer *hi, Integer foff, Integer fsize,void *buf, Integer *ld, Integer *nbhandle);
extern void wnga_nbput_field(Integer g_a, Integer *lo, Integer *hi, Integer foff, Integer fsize,void *buf, Integer *ld, Integer *nbhandle);
extern void wnga_put_field(Integer g_a, Integer *lo, Integer *hi, Integer foff, Integer fsize,void *buf, Integer *ld);

#endif /* WAPI_H_ */
