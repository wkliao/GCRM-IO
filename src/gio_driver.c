/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio_driver.c 4609 2017-12-07 07:26:38Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "gcrm.h"
#include "gio.h"

static void gio_set_mask_sizes(gio_parameters *gio_param);

/*----< gio_driver() >--------------------------------------------------------*/
/* Main driver for writing history and restart files using the gio
 * descriptors.
 * Averages are computed for variables specified as averages.
 * @param model_time the current model time used to decide if it is time
 * to write
 */
void gio_driver(gio_parameters *gio_param,
                double          model_time)
{
    int i, j, now, restart_now, created;
    char time_str[24], filename[512];
    double dbl_freq, local_start;
    static int  wrote_grid=0;

    file_descriptor *gio_files = gio_param->file.gio_files;

#define alwaysclose
    /* if alwaysclose defined, each file will be closed after each write.
     * Since this has known severe performance problems on some systems
     * (hopper), the other option is to close at the restart interval
     * and/or when we change to a new file.
     */

    /* duplicate the if check here to reduce function call overhead */
    if (! gio_param->is_schedule_set)
        gio_init_schedule(gio_param);

    /* accumulate averages first */
    gio_accumulate_all_fields(gio_param);

    /* put barrier and timer here after averaging since averaging should not be 
       included in io times. */
    MPI_Barrier(MPI_COMM_WORLD);
    local_start = MPI_Wtime();

    /* Compute this to optimize open/close by only closing (and reopening)
     * at the restart interval */
    restart_now = -1;
    dbl_freq = gio_files[gio_param->file.restart_idx].frequency;
    if (dbl_freq != 0) restart_now = (int)model_time % (int)dbl_freq;

    /* set the cluster mask counts */
    gio_set_mask_sizes(gio_param);
    gio_param->clus_cells = 0;
    for (i=0; i<gio_param->grid.total_grid_blocks; i++)
        gio_param->clus_cells += gio_param->cell_mask_sizes[i];
    gio_param->clus_edges = gio_param->clus_cells;
    gio_param->clus_crns  = gio_param->clus_cells;

    for (i=0; i<gio_param->file.num_files; i++) {
        /* If descriptor is not completely registered, skip to next file */
        if (! gio_files[i].complete) continue;

        /* If descriptor is for the special separate grid file (there should
           only be one), open the file which populates the grid, then close and
           we are done. */
        if (gio_files[i].grid_ONLY) {
            if (! wrote_grid) {
                sprintf(filename, "%s%s.nc",
                        gio_param->defaults.gio_default_cdf_path,
                        gio_param->defaults.gio_sep_grid_file);
                created = gio_open_file(gio_param, i, filename);
                gio_close_file(gio_param, i);
                wrote_grid = 1;
            }
            continue;
        }

        /* Decide if it is time to write data to this file (each file can have
           a different write frequency) */
        now = 1;
        if (gio_files[i].frequency > 0)
            now = (int)model_time % (int)(gio_files[i].frequency);
        if (now == 0) {  /* 0 means to write now */
            gio_param->num_dumps++;

            if (gio_files[i].samples_written >= gio_files[i].nsamples ||
                gio_files[i].last_time_str[0] == '\0') {
                gio_convert_time(model_time, time_str);
                strcpy(gio_files[i].last_time_str, time_str);
#ifndef alwaysclose
                if (gio_files[i].fileID != -1)
                    /* No need to bcast this since all procs execute and it
                       always resets the fileID */
                    gio_close_file(gio_param, i);
#endif
            } else
                strcpy(time_str, gio_files[i].last_time_str);

            if (gio_files[i].is_clustered_file)
                gio_convert_time(model_time, time_str);

            /* compose the file name */
            sprintf(filename, "%s%s%s.nc", gio_files[i].directory_name,
                              gio_files[i].file_prefix, time_str);

            if (i == gio_param->file.restart_idx &&
                gio_param->file.gio_restart_overwrite)
                sprintf(filename, "%s%s", gio_files[i].directory_name,
                                  gio_files[i].file_prefix);

            /* Clustered */
            if (gio_files[i].is_clustered_file) {
                if (gio_param->clus_cells == 0 ||
                    gio_param->clus_crns  == 0 ||
                    gio_param->clus_edges == 0)
                    continue;
            }

            /* fileID should have been initialized to -1; use it to decide if
               we need to open/create a file, so we can close less frequently
               on systems where on systems where open/close is expensive
               (lustre)
             */
            if (gio_files[i].fileID == -1) {
                /* Must determine filename and create/open the file
                   open will abort if error; file id is set upon success */
                created = gio_open_file(gio_param, i, filename);
                if (gio_param->num_io_processors != gio_param->gio_nprocs)
                    MPI_Bcast(&gio_files[i].fileID, 1, MPI_INT, 0,
                              gio_param->gio_world_comm);
                if (created) gio_files[i].samples_written = 0;
            }

            /* Write non-grid fields, Grid fields have been written during
               the call to gio_open_file(). Available options are writing
               grid fields in a separate, grid-only file or grid fields in
               all files. */

            /* allocate nonblocking I/O request array for all non-grid fields */
            if (gio_param->gio_iotype != MODE_BLOCKING_INDEP)
                gio_allocate_req_ids(gio_param, i, 0, -1);

            /* update time field (grid_only files won't have time field) */
            if (! gio_files[i].grid_ONLY)
                gio_update_time_field(gio_param, i, &model_time);

            for (j=0; j<gio_files[i].nflds; j++) {
                /* check whether we are writing this field */
                if (! gio_is_output_field(gio_param, i, j)) continue;

                gio_write_field(gio_param, i, gio_files[i].fields[j]);

                if (gio_param->gio_iotype == MODE_BLOCKING_COLL)
                    /* TODO: this is to flushing nonblocking calls per field,
                     * not a true blocking collective write. */
                    gio_nonblocking_io_wait(gio_param, i);
            }

            /* wait for all nonblocking requests to complete */
            if (gio_param->gio_iotype == MODE_NONBLOCKING_COLL ||
                gio_param->gio_iotype == MODE_NONBLOCKING_INDEP)
                gio_nonblocking_io_wait(gio_param, i);
            /* free allocated space */
            if (gio_param->using_interleaved == 0 && gio_param->using_direct == 0)
                gio_deallocate_io_buffer(gio_param, i, 0, -1);
            gio_free_req_ids(gio_param);

#ifdef alwaysclose
            gio_close_file(gio_param, i);
#else
            if (gio_files[i].is_clustered_file)
                gio_close_file(gio_param, i);
            else if (i == gio_param->file.restart_idx)
                gio_close_file(gio_param, i);
            else if (restart_now == 0 && model_time != 0)
                gio_close_file(gio_param, i);
#endif
            gio_files[i].samples_written++;

        } /* End if (now == 0) */
    } /* Main loop on files */

    gio_param->total_time_in_API += MPI_Wtime() - local_start;
}

/*----< gio_set_mask_sizes() >------------------------------------------------*/
/* Function to determine if block is being written out based on mask array */
static
void gio_set_mask_sizes(gio_parameters *gio_param)
{
    int i, j, ldx, mask_idx, num_cells, nsd;
    int imin, imax, jmin, jmax, ilen, jlen;
    double *mask_ptr;
    int *tmp_cell, *tmp_edge, *tmp_crnr;
    int total_number_of_cells=0, rank;
    int total_grid_blocks = gio_param->grid.total_grid_blocks;

    data_descriptor *gio_descriptors = gio_param->data.gio_descriptors;
    MPI_Comm_rank(gio_param->gio_io_comm, &rank);

    /* Find index of mask array */
    i = 0;
    while (i < gio_param->data.num_fields)
        if (!strncmp(gio_descriptors[i++].field_name, "clus_mask", 9))
            break;

    if (i == gio_param->data.num_fields) {
        printf("Error: Failure in gio_set_mask_size to find mask array\n");
        ABORT
    }
    mask_idx = i;

    tmp_cell = (int*) tcalloc(total_grid_blocks, sizeof(int));
    tmp_edge = (int*) tcalloc(total_grid_blocks, sizeof(int));
    tmp_crnr = (int*) tcalloc(total_grid_blocks, sizeof(int));

    for (nsd=0; nsd<gio_descriptors[mask_idx].num_data_blocks; nsd++) {
        num_cells = 0;
        mask_ptr = gio_descriptors[mask_idx].gio_data_blocks[nsd].double_data;
        imin = gio_descriptors[mask_idx].gio_data_blocks[nsd].imin;
        imax = gio_descriptors[mask_idx].gio_data_blocks[nsd].imax;
        jmin = gio_descriptors[mask_idx].gio_data_blocks[nsd].jmin;
        jmax = gio_descriptors[mask_idx].gio_data_blocks[nsd].jmax;
        ilen = imax - imin + 1;
        jlen = jmax - jmin + 1;
        for (j=0; j<jlen; j++) {
            for (i=0; i<ilen; i++) {
                if (mask_ptr[j * gio_param->grid.gio_istride + i] != 0.0)
                    num_cells++;
            }
        }
        total_number_of_cells += num_cells;
        if (num_cells != 0) {
            ldx = gio_descriptors[mask_idx].gio_data_blocks[nsd].block_index;
            if (ldx == 0)
              tmp_cell[ldx] = ilen*jlen+2;
            else
              tmp_cell[ldx] = ilen*jlen;

            tmp_edge[ldx] = 3*ilen*jlen;
            tmp_crnr[ldx] = 2*ilen*jlen;
        }
    }

    MPI_Allreduce(tmp_cell, gio_param->cell_mask_sizes, total_grid_blocks,
                  MPI_INT, MPI_SUM, gio_param->gio_world_comm);
    MPI_Allreduce(tmp_edge, gio_param->edge_mask_sizes, total_grid_blocks,
                  MPI_INT, MPI_SUM, gio_param->gio_world_comm);
    MPI_Allreduce(tmp_crnr, gio_param->crnr_mask_sizes, total_grid_blocks,
                  MPI_INT, MPI_SUM, gio_param->gio_world_comm);

    tfree(tmp_cell);
    tfree(tmp_edge);
    tfree(tmp_crnr);
}
