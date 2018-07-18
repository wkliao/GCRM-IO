/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "gcrm.h"
#include "gio.h"
#include "util.h"

/*  Evaluate the morton index for an element in an lsize x lsize block
 *  from the standard cartesian i,j indices.
 */
/*----< gio_local_morton_index() >--------------------------------------------*/
static
int gio_local_morton_index(int y,
                           int x,
                           int res)
{
    int k, idx, xbit[20], ybit[20];

    for (k=0; k<res; k++) {
        xbit[k] = x % 2;
        x = (x - xbit[k]) / 2;
        ybit[k] = y % 2;
        y = (y - ybit[k]) / 2;
    }
    idx = 0;
    for (k=res-1; k>=0; k--)
        idx = 4*idx + (2*ybit[k] + xbit[k]);

    return idx;
/* should we return idx as index in C strats with 0
    return (idx + 1);
*/
}

/*----< load_config_info() >--------------------------------------------------*/
static
int load_config_info(gio_parameters *gio_param,
                     char           *config_file,
                     char           *output_config)
{
    /* read desc file */
    gio_read_data_config(gio_param, config_file);

    /* read fcfg file */
    gio_read_file_config(gio_param, output_config);

    return 1;
}

/*----< gio_init() >----------------------------------------------------------*/
/*
In GIO library, gio_init() is defined as the following.
subroutine
gio_init(iotype,              &  ! direct, interleaved, or nonblocking
         res,                 &  ! resolution
         level,               &  ! level
         blocksize,           &  ! blocksize
         istride,             &  ! stride size for i index
         jstride,             &  ! stride size for j index
         config_file,         &  ! config file to read
         output_config,       &  ! config file that specifies output file(s)
         num_ioprocs,         &  ! number of IO procs
         cdf_output_path,     &  ! default cdf file location
         cdf_freq,            &  ! default cdf frequency
         cdf_nsamples,        &  ! default number of samples per file
         cdf_grid_opt,        &  ! default cdf grid option
         cdf_sep_grid_file,   &
         prestart_fname,      &
         prestart_interval,   &
         prestart_overwrite)

In ZGrd_initialize.F90, gio_init() is called in the following.
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
int gio_init(gcrm_parameters *gcrm_param,
             gio_parameters  *gio_param)
{
    /* In gio_init(), there is load_config_info() that reads
     * read IO descriptor parameter file, ZGrd.desc
     * In load_config_info(), there are
     * gio_read_data_config() and gio_read_file_config()
     */
    int i, j;
    double starttime = MPI_Wtime();

    gio_defaults *defaults = &gio_param->defaults;
    gio_data     *data     = &gio_param->data;
    gio_grid     *grid     = &gio_param->grid;
    gio_dim      *dim      = &gio_param->dim;
    gio_file     *file     = &gio_param->file;
    MODULE_ZGrd_output *zgrd_output = &gcrm_param->zgrd_output;

    gio_param->used_info               = MPI_INFO_NULL;
    gio_param->num_records_written     = 0;
    gio_param->num_records_read        = 0;
    gio_param->time_in_nf_put_var      = 0.0;
    gio_param->time_in_nf_put_var_grid = 0.0;
    gio_param->time_in_nf_get_var      = 0.0;
    gio_param->total_time_in_API       = 0.0;
    gio_param->time_in_nf_create       = 0.0;
    gio_param->time_in_nf_open         = 0.0;
    gio_param->time_in_nf_close        = 0.0;
    gio_param->time_in_nf_inq_varid    = 0.0;
    gio_param->time_in_nf_inq_dimlen   = 0.0;
    gio_param->time_in_nf_put_att      = 0.0;
    gio_param->time_in_nf_get_att      = 0.0;
    gio_param->time_in_nf_def_dim      = 0.0;
    gio_param->time_in_nf_def_var      = 0.0;
    gio_param->time_in_nf_enddef       = 0.0;
    gio_param->time_in_update_time     = 0.0;
    gio_param->time_in_API_copy        = 0.0;
    gio_param->time_in_avgs            = 0.0;
    gio_param->time_in_nf_iput         = 0.0;
    gio_param->time_in_nf_iget         = 0.0;
    gio_param->time_in_nf_wait         = 0.0;
    gio_param->bytes_API_write         = 0;
    gio_param->bytes_API_read          = 0;
    gio_param->pnc_write_amount        = 0;
    gio_param->pnc_read_amount         = 0;
    gio_param->num_dumps               = 0;

    gio_param->gio_world_comm = MPI_COMM_WORLD;

    MPI_Comm_size(MPI_COMM_WORLD, &gio_param->gio_nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &gio_param->gio_me);

    for (i=0; i<MAX_FILES; i++) {
        file->gio_files[i].fileID = -1;
        for (j=0; j<MAX_DATA_FIELDS; j++)
            file->gio_files[i].fields[j] = -1;
    }

    /* Setup iotype.  nonblocking is default if no match is found */
    gio_param->gio_iotype           = MODE_NONBLOCKING_COLL;
    gio_param->req_ids               = NULL;
    gio_param->using_interleaved     = 0;
    gio_param->using_direct          = 0;
    if (!strcmp(zgrd_output->iotype, "nonblocking_collective"))
        gio_param->gio_iotype = MODE_NONBLOCKING_COLL;
    else if (!strcmp(zgrd_output->iotype, "nonblocking_independent"))
        gio_param->gio_iotype = MODE_NONBLOCKING_INDEP;
    else if (!strcmp(zgrd_output->iotype, "blocking_collective"))
        gio_param->gio_iotype = MODE_BLOCKING_COLL;
    else if (!strcmp(zgrd_output->iotype, "blocking_independent"))
        gio_param->gio_iotype = MODE_BLOCKING_INDEP;
#if 0
    /* not supported yet */
    else if (!strcmp(zgrd_output->iotype, "interleaved")) {
        gio_param->using_interleaved = 1;
        gio_param->gio_iotype = 0;
    }
    else if (!strcmp(zgrd_output->iotype, "direct")) {
        /* not supported yet */
        gio_param->using_direct = 1;
        gio_param->gio_iotype = 0;
    }
#endif
    else {
        fprintf(stderr,"Error: unsupported I/O method %s\n",zgrd_output->iotype);
        return 0;
    }

    strcpy(file->gio_restart_fname, zgrd_output->restart_output);
    file->gio_restart_interval  = zgrd_output->restart_interval;
    file->gio_restart_overwrite = zgrd_output->l_restart_overwrite;

    gio_param->is_schedule_set   = 0;
    gio_param->using_collectives = 1;

    if (zgrd_output->cdf_output_path[0] != '\0') {
        int len = strlen(zgrd_output->cdf_output_path);
        strcpy(defaults->gio_default_cdf_path, zgrd_output->cdf_output_path);

        if (zgrd_output->cdf_output_path[len-1] != '/') {
            defaults->gio_default_cdf_path[len] = '/';
            defaults->gio_default_cdf_path[len+1]   = '\0';
        }
    }

    if (zgrd_output->cdf_output_freq > 0)
        /* Changed to seconds to avoid roundoff problems */
        defaults->gio_default_frequency = zgrd_output->cdf_output_freq;

    if (zgrd_output->cdf_output_nsamples > 0)
       defaults->gio_default_nsamples = zgrd_output->cdf_output_nsamples;

    if (!strncmp(zgrd_output->cdf_grid_option, "nogrid", 6))
        defaults->gio_grid_option = 1;
    else if (!strncmp(zgrd_output->cdf_grid_option, "all_files", 9))
        defaults->gio_grid_option = 2;
    else if (!strncmp(zgrd_output->cdf_grid_option, "sep_grid", 8)) {
        defaults->gio_grid_option = 3;
        if (zgrd_output->cdf_sep_grid[0] == '\0')
            strcpy(defaults->gio_sep_grid_file, "grid");
        else
            strcpy(defaults->gio_sep_grid_file, zgrd_output->cdf_sep_grid);
    }

    /*  Set basic grid parameters */

    grid->refine_param = gcrm_param->grid_params.level_max;
    grid->nlevel       = gcrm_param->vertical.km;
    grid->grid_length  = POWER2(gcrm_param->grid_params.level_max);

    grid->gio_grid_size    = 10 * grid->grid_length * grid->grid_length + 2;
    grid->gio_corners_size = 2 * (grid->gio_grid_size - 2);
    grid->gio_edges_size   = 3 * (grid->gio_grid_size - 2);

    grid->gio_istride = gcrm_param->grid_params.im;
    grid->gio_jstride = gcrm_param->grid_params.jm;

    gio_param->num_io_processors = gcrm_param->npe_io;
/* This benchmark only support nonblocking I/O method where all processes
 * are the I/O processes.

    if (gio_param->num_io_processors > gio_param->gio_nprocs ||
        gio_param->num_io_processors < 1) {
        printf("Warning: invalid ioproc request (%d) - setting to all processors: %d\n",
               gio_param->num_io_processors, gio_param->gio_nprocs);
        gio_param->num_io_processors = gio_param->gio_nprocs;
    }
*/

    /* blocking and nonblocking I/O methods make all processes do I/O */
    if (gio_param->using_interleaved == 0 && gio_param->using_direct == 0)
        gio_param->num_io_processors = gio_param->gio_nprocs;

    int io_block_size = gio_param->gio_nprocs / gio_param->num_io_processors;
    gio_param->gio_do_io = (gio_param->gio_me % io_block_size) ? 0 : 1;

    /* gcrm_param->grid_params.im is set in init_MODULE_grid_params()
       grid->block_size is used in gio library only */
    grid->block_size      = gcrm_param->grid_params.im - 2;
    grid->num_grid_blocks = 0;

    /* Initialize fields to null values (if necessary) */
    grid->north_pole_flag = 0;
    grid->south_pole_flag = 0;

    /* initialize data descriptors */
    /*   call init_gio_datadesc() */

    for (i=0; i<MAX_DATA_FIELDS; i++) {
        data->gio_descriptors[i].io_buf       = NULL;
        data->gio_descriptors[i].index_choice = -1;
        data->gio_descriptors[i].cell_idx     = -1;
        data->gio_descriptors[i].time_idx     = -1;
        data->gio_descriptors[i].intg_idx     = -1;
    }

    /* initialize arrays for managing extra dimensions */
    for (i=0; i<MAX_INDICES; i++) {
        dim->index_sizes[i] = 0;
        dim->index_names[i][0] = '\0';
        dim->species_sizes[i] = 0;
        dim->species_names[i][0] = '\0';
    }

    /* read desc and fcfg files */
    load_config_info(gio_param, zgrd_output->descfile,
                                zgrd_output->outfileconfig);

    /* determine total number of grid blocks */
    grid->total_grid_blocks = 10 * POWER2(grid->grid_length/grid->block_size);
    gio_param->cell_mask_sizes = calloc_1D_int(grid->total_grid_blocks);
    gio_param->edge_mask_sizes = calloc_1D_int(grid->total_grid_blocks);
    gio_param->crnr_mask_sizes = calloc_1D_int(grid->total_grid_blocks);

    /* allocate 2D array with Morton-ordering index map and initialize it */
    grid->grid_map = (int**) calloc_2D_int(grid->block_size, grid->block_size);

    for (j=0; j<grid->block_size; j++)
        for (i=0; i<grid->block_size; i++)
            grid->grid_map[j][i] = gio_local_morton_index(j, i,
                                       gcrm_param->grid_params.level_max);

    gio_param->total_time_in_API += MPI_Wtime() - starttime;

    return 1;
}

/*----< gio_grid_finalize() >-------------------------------------------------*/
static
void gio_grid_finalize(gio_parameters *gio_param)
{
    int i;
    grid_block *gio_grid_blocks = gio_param->grid.gio_grid_blocks;

    /* the buffers allocation is done in gio_register.c */
    for (i=0; i<gio_param->grid.num_grid_blocks; i++) {
        free_1D_dbl(gio_grid_blocks[i].grid_center_lat);
        free_1D_dbl(gio_grid_blocks[i].grid_center_lon);
        free_2D_dbl(gio_grid_blocks[i].grid_corner_lat);
        free_2D_dbl(gio_grid_blocks[i].grid_corner_lon);
        free_2D_dbl(gio_grid_blocks[i].grid_edge_lat);
        free_2D_dbl(gio_grid_blocks[i].grid_edge_lon);
        free_2D_int(gio_grid_blocks[i].grid_neighbors);
        free_2D_int(gio_grid_blocks[i].grid_corners);
        free_2D_int(gio_grid_blocks[i].grid_edges);
        free_3D_int(gio_grid_blocks[i].grid_endpoints);
    }
    for (i=0; i<gio_param->data.num_fields; i++) {
        if (gio_param->data.gio_descriptors[i].has_levels &&
            gio_param->data.gio_descriptors[i].float_level_data != NULL) {
            free_1D_flt(gio_param->data.gio_descriptors[i].float_level_data);
            gio_param->data.gio_descriptors[i].float_level_data = NULL;
        }
    }
}

/*----< gio_terminate() >-----------------------------------------------------*/
void gio_terminate(gio_parameters *gio_param)
{
    int i;
    double starttime = MPI_Wtime();

    /* Close any open files */
    for (i=0; i<gio_param->file.num_files; i++) {
        if (gio_param->file.gio_files[i].fileID != -1) {
           gio_close_file(gio_param, i);
           gio_param->file.gio_files[i].fileID = -1;
        }
    }

    /* free up memory allocated and used by gio library */
    free_1D_int(gio_param->cell_mask_sizes); /* gio_init() */
    free_1D_int(gio_param->edge_mask_sizes);
    free_1D_int(gio_param->crnr_mask_sizes);
    free_2D_int(gio_param->grid.grid_map);

    gio_free_avg_descriptor(gio_param);
    gio_grid_finalize(gio_param);

    /* generate a report of API timing and other stats */
    // gio_print_timing_info();
    gio_param->total_time_in_API += MPI_Wtime() - starttime;
}

