/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio.h 4609 2017-12-07 07:26:38Z wkliao $
 */

#ifndef H_GIO
#define H_GIO

#include <errno.h>
extern int errno;

/* constants defined in gio_params.F90 */
#define MAX_DATA_FIELDS        256
#define MAX_DATA_BLOCKS         64
#define MAX_STATIC_GRID_FIELDS  32
#define MAX_FILES              100
#define MAX_PROPS               32
#define MAX_DIMS                32
#define MAX_ATTR                32
#define MAX_FLD_PER_FILE        64
#define IO_AGGREGATION           8
#define MAX_INDICES              2

#define POWER2(x) (1 << (x))

int debug;

#define MODE_NONBLOCKING_COLL  0
#define MODE_NONBLOCKING_INDEP 1
#define MODE_BLOCKING_COLL     2
#define MODE_BLOCKING_INDEP    3

#include <gio_defaults.h>
#include <gio_file.h>
#include <gio_prop.h>
#include <gio_data.h>
#include <gio_grid.h>
#include <gio_dim.h>

/*
   call gio_init(iotype,             &
                 level_max,          &
                 km,                 &
                 block_size,         &
                 im,                 &
                 jm,                 &
                 descfile,           &
                 outfileconfig,      &
                 npe_io,             &
                 cdf_output_path,    &  ! default cdf file location
                 cdf_output_freq,    &  ! default cdf frequency
                 cdf_output_nsamples,&  ! default number of samples per file
                 cdf_grid_option,    &  ! default cdf grid option
                 cdf_sep_grid,       &  ! File for separate grid files
                 restart_output,     &
                 restart_interval,   &
                 l_restart_overwrite )
*/
typedef struct {
    /* from initialize_ZGrd() */
    int  block_size;

    int  nlevel;       /* will be set to grid, level, km in gio_init() */

    /* from MODULE gio_sched */
    MPI_Comm  gio_world_comm, gio_io_comm, gio_block_comm;
    int  gio_nprocs;
    int  gio_me;
    int  num_io_processors;
    int  gio_iotype;
    int  gio_do_io;
    int  using_interleaved;
    int  using_direct;
    int  is_schedule_set;
    int  using_collectives;
    int *cell_mask_sizes;
    int *edge_mask_sizes;
    int *crnr_mask_sizes;

    gio_defaults    defaults;
    gio_prop        props;
    gio_grid        grid;
    gio_data        data;
    gio_dim         dim;
    gio_file        file;

    /* from io_netcdfx.F90: counts for cluster bit mask */
    int  clus_cells, clus_edges, clus_crns;
    int *req_ids;
    int  num_reqs;

    /* added to track the number of dumps */
    int  num_dumps;

    /* from module gio_statistics */
    MPI_Info used_info;
    int    num_records_written;
    int    num_records_read;
    double total_time_in_API;
    double time_in_nf_put_var;
    double time_in_nf_put_var_grid;
    double time_in_nf_get_var;
    double time_in_nf_create;
    double time_in_nf_open;
    double time_in_nf_close;
    double time_in_nf_inq_varid;
    double time_in_nf_inq_dimlen;
    double time_in_nf_put_att;
    double time_in_nf_get_att;
    double time_in_nf_def_dim;
    double time_in_nf_def_var;
    double time_in_nf_enddef;
    double time_in_update_time;
    double time_in_API_copy;
    double time_in_avgs;
    double time_in_nf_iput;
    double time_in_nf_iget;
    double time_in_nf_wait;

    int        num_grid_vars;
    int        num_fld_vars;
    MPI_Offset bytes_API_write;
    MPI_Offset bytes_API_read;
    MPI_Offset pnc_write_amount;
    MPI_Offset pnc_read_amount;
} gio_parameters;


void check(int status, char *message, int nodie);
int gio_open(gio_parameters *gio_param, MPI_Comm comm, char *fname, int flags,
                      MPI_Info *used_info, int *fid);
int gio_inq_varid(gio_parameters *gio_param, int fid, char *fieldname, int *fldid);
int gio_open_file(gio_parameters *gio_param, int ifile, char *fname);
void gio_close_file(gio_parameters *gio_param, int ifile);
void gio_sync_file(gio_parameters *gio_param, int fileID);
int total_number_of_cells(gio_parameters *gio_param);
void gio_write(gio_parameters *gio_param, int field_idx, int ifile, void *buffer,
                      MPI_Offset offset, int time_idx, MPI_Offset blen);
void gio_read(gio_parameters *gio_param, int field_idx, int ifile, void *buffer,
                      MPI_Offset offset, int time_idx, MPI_Offset blen);
void gio_update_time_field(gio_parameters *gio_param, int ifile, double *the_time);
int gio_open_restart_read(gio_parameters *gio_param, int ifile, char *fname);
void gio_read_restart_attrs(gio_parameters *gio_param, int ifile);
void gio_read_time_field(gio_parameters *gio_param, int ifile,
                      double *the_time, MPI_Offset *nsteps);

void gio_term_statistics(gio_parameters *gio_param);

#include <gcrm.h>

/* API declartions */
int gio_init(gcrm_parameters *gcrm_param, gio_parameters *gio_param);
void gio_terminate(gio_parameters *gio_param);
int gio_read_file_config(gio_parameters*, const char*);
int gio_read_data_config(gio_parameters*, const char*);
void gio_set_strides(gio_parameters *gio_param);
int gio_init_dims(gio_dim *dim);
void gio_set_dimension(int idx, int num_dims, char **dims_str, gio_data *data,
                       gio_grid *grid, gio_dim *dim);
void gio_set_limits(gio_parameters *gio_param, int idx, int *imin, int *imax);
void gio_check_file_data(gio_parameters *gio_param);
void gio_grid_setup(gio_parameters *gio_param, int **nghbr_tags,
                    double **center, double ***corner, double *grid_area,
                    int panel_id, int gmin[2], int gmax[2]);
void gio_grid_setup_pole(gio_parameters *gio_param, int *nghbr_tags, double *center,
                         double **corner, double grid_area, char pole_id);
void gio_register_dfield(char *field_name, double *darray, int lo[2], int hi[2],
                         int panel_id, gio_parameters *gio_param);
void gio_register_ifield(char *field_name, int    *iarray, int lo[2], int hi[2],
                         int panel_id, gio_parameters *gio_param);
void gio_register_dpole(char *field_name, double *dval, char pole_id, gio_data *data);
void gio_register_ipole(char *field_name, int *dval, char pole_id, gio_data *data);
void gio_register_dlevel(gio_data *data, char *field_name, double *dvals, int length);
void gio_register_index(gio_dim *dim, char *t_index_name, int *iarray,
                        int ival, int length);
void gio_register_restart_dbl(gio_data *data, const char *t_name, double t_data);
void gio_register_restart_int(gio_data *data, const char *t_name, int t_data);

/* from gio_driver.c */
void gio_driver(gio_parameters *gio_param, double model_time);

/* from average.c */
void gio_allct_avg_descriptor(gio_parameters *gio_param, int field_idx);
void gio_free_avg_descriptor(gio_parameters *gio_param);
void gio_accumulate_field_data(gio_parameters *gio_param, int field_idx);
void gio_accumulate_all_fields(gio_parameters *gio_param);
void gio_eval_average(gio_parameters *gio_param,
                      int             field_idx,
                      char            pole,
                      int             iblock,
                      int             time_idx,
                      int             species_idx,
                      float          *buffer);
void gio_reset_average(gio_parameters *gio_param,
                       int             field_idx);

/* from ZGrd_output.c */
void ZGrd_write_output(gio_parameters *gio_param, double time_ZGrd);

/* from gio_sched.c */
void gio_init_schedule(gio_parameters *gio_param);

/* from gio_utility.c */
void gio_convert_time(double  model_time, char *time_str);

/* from nonblock_io.c */
int gio_iput_double(gio_parameters *gio_param,
                    int             fileID,
                    int             fieldID,
                    MPI_Offset     *start,
                    MPI_Offset     *count,
                    double         *buffer);
int gio_iget_double(gio_parameters *gio_param,
                    int             fileID,
                    int             fieldID,
                    MPI_Offset     *start,
                    MPI_Offset     *count,
                    double         *buffer);
int gio_iput_real(gio_parameters *gio_param,
                  int             fileID,
                  int             fieldID,
                  MPI_Offset     *start,
                  MPI_Offset     *count,
                  float          *buffer);
int gio_iget_real(gio_parameters *gio_param,
                  int             fileID,
                  int             fieldID,
                  MPI_Offset     *start,
                  MPI_Offset     *count,
                  float          *buffer);
int gio_iput_int(gio_parameters *gio_param,
                 int             fileID,
                 int             fieldID,
                 MPI_Offset     *start,
                 MPI_Offset     *count,
                 int            *buffer);
int gio_iget_int(gio_parameters *gio_param,
                 int             fileID,
                 int             fieldID,
                 MPI_Offset     *start,
                 MPI_Offset     *count,
                 int            *buffer);
int gio_get_num_w_reqs(gio_parameters *gio_param,
                       int             field_idx);
int gio_is_output_field(gio_parameters *gio_param,
                        int             ifile,
                        int             ifield);
int gio_is_output_grid(gio_parameters *gio_param,
                       int             ifile,
                       int             ifield);
void gio_allocate_req_ids(gio_parameters *gio_param,
                          int             ifile,
                          int             is_grid,
                          int             ifield);
void gio_deallocate_io_buffer(gio_parameters *gio_param,
                              int             ifile,
                              int             is_grid,
                              int             ifield);
void gio_write_verticalgrid(gio_parameters *gio_param,
                            int             ifile,
                            int             field_idx);
void gio_write_field(gio_parameters *gio_param,
                     int             ifile,
                     int             field_idx);
void gio_nonblocking_io_wait(gio_parameters *gio_param,
                             int             ifile);

/* from copy.c */
void gio_copy_to_buffer(gio_parameters *gio_param,
                        int             field_idx,
                        char            pole,
                        int             iblock,
                        int             time_idx,
                        int             spc_idx,
                        void           *buffer);

extern int
gio_put_real(gio_parameters *gio_param, int fid, int fldid, 
             MPI_Offset *start, MPI_Offset *count, float *buffer);

extern void
gio_free_req_ids(gio_parameters *gio_param);

#endif
