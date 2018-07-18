/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: nonblocking_io.c 4699 2018-07-18 00:36:54Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <pnetcdf.h>

#include <gcrm.h>
#include <gio.h>

static
void gio_write_block(gio_parameters *gio_param, int ifile, int field_idx,
                     int iblock, int buf_off, char pole, MPI_Offset offset,
                     MPI_Offset blen, int idx, int isp);

/*! ---------------------------------------------------------------------------
  ! Purpose: Nonblocking write a double variable
  ! --------------------------------------------------------------------------*/
int gio_iput_double(gio_parameters *gio_param,
                    int             fileID,
                    int             fieldID,
                    MPI_Offset     *start,
                    MPI_Offset     *count,
                    double         *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    err = ncmpi_iput_vara_double(fileID, fieldID, start, count, buffer,
                                 &gio_param->req_ids[gio_param->num_reqs++]);
    gio_param->time_in_nf_iput += MPI_Wtime() - starttime;
    return err;
}

int gio_iget_double(gio_parameters *gio_param,
                    int             fileID,
                    int             fieldID,
                    MPI_Offset     *start,
                    MPI_Offset     *count,
                    double         *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    err = ncmpi_iget_vara_double(fileID, fieldID, start, count, buffer,
                                 &gio_param->req_ids[gio_param->num_reqs++]);
    gio_param->time_in_nf_iget += MPI_Wtime() - starttime;
    return err;
}

/*! ---------------------------------------------------------------------------
  ! Purpose: Nonblocking write a real variable
  ! --------------------------------------------------------------------------*/
int gio_iput_real(gio_parameters *gio_param,
                  int             fileID,
                  int             fieldID,
                  MPI_Offset     *start,
                  MPI_Offset     *count,
                  float          *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    err = ncmpi_iput_vara_float(fileID, fieldID, start, count, buffer,
                                &gio_param->req_ids[gio_param->num_reqs++]);
    gio_param->time_in_nf_iput += MPI_Wtime() - starttime;
    return err;
}

int gio_iget_real(gio_parameters *gio_param,
                  int             fileID,
                  int             fieldID,
                  MPI_Offset     *start,
                  MPI_Offset     *count,
                  float          *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    err = ncmpi_iget_vara_float(fileID, fieldID, start, count, buffer,
                                &gio_param->req_ids[gio_param->num_reqs++]);
    gio_param->time_in_nf_iget += MPI_Wtime() - starttime;
    return err;
}

/*! ---------------------------------------------------------------------------
  ! Purpose: Nonblocking write an integer variable
  ! --------------------------------------------------------------------------*/
int gio_iput_int(gio_parameters *gio_param,
                 int             fileID,
                 int             fieldID,
                 MPI_Offset     *start,
                 MPI_Offset     *count,
                 int            *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    err = ncmpi_iput_vara_int(fileID, fieldID, start, count, buffer,
                              &gio_param->req_ids[gio_param->num_reqs++]);
    gio_param->time_in_nf_iput += MPI_Wtime() - starttime;
    return err;
}

int gio_iget_int(gio_parameters *gio_param,
                 int             fileID,
                 int             fieldID,
                 MPI_Offset     *start,
                 MPI_Offset     *count,
                 int            *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    err = ncmpi_iget_vara_int(fileID, fieldID, start, count, buffer,
                              &gio_param->req_ids[gio_param->num_reqs++]);
    gio_param->time_in_nf_iget += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_get_num_w_reqs() >-------------------------------------------------*/
/* calculate the number of write requests (including all blocks and poles) */
int gio_get_num_w_reqs(gio_parameters *gio_param,
                       int             field_idx)
{
    int idx_min, idx_max, imin4, imax4, num_w_reqs;

    /* make writing level data part of non-blocking I/O */
    if (gio_param->data.gio_descriptors[field_idx].is_level_data)
        return 1;

    num_w_reqs = gio_param->grid.num_grid_blocks;
    if (! gio_param->data.gio_descriptors[field_idx].no_poles) {
        if (gio_param->grid.north_pole_flag) /* this proc has north pole */
            num_w_reqs++;
        if (gio_param->grid.south_pole_flag) /* this proc has south pole */
            num_w_reqs++;
    }
    imin4 = gio_param->data.gio_descriptors[field_idx].imin4;
    imax4 = gio_param->data.gio_descriptors[field_idx].imax4;
    gio_set_limits(gio_param, field_idx, &idx_min, &idx_max);
    num_w_reqs *= (idx_max - idx_min + 1) * (imax4 - imin4 + 1);

    return num_w_reqs;
}

/*---< gio_is_output_field() >------------------------------------------------*/
/* check if we are outputing this field */
int gio_is_output_field(gio_parameters *gio_param,
                        int             ifile,
                        int             ifield) 
{
    if (gio_param->file.gio_files[ifile].is_grid[ifield] ||
        gio_param->file.gio_files[ifile].is_time_field[ifield])
        /* Skip this field ... it is static data */
        return 0;

    if (gio_param->data.gio_descriptors[gio_param->file.gio_files[ifile].fields[ifield]].no_data)
        /* no model data exists for this field */
        return 0;

    return 1;
}

/*---< gio_is_output_grid() >-------------------------------------------------*/
/* check if we are outputing this grid */
int gio_is_output_grid(gio_parameters *gio_param,
                       int             ifile,
                       int             ifield)
{
    int field_idx;

    if (gio_param->file.gio_files[ifile].fields[ifield] == -1) return 0;

    /* field is a known one */
    field_idx = gio_param->file.gio_files[ifile].fields[ifield];

    /* skip dummy variables */
    if (gio_param->data.gio_descriptors[field_idx].no_data) return 0;

    /* skip vertical variables; handled elsewhere */
    // if (! gio_param->gio_do_nonblocking_io && gio_param->data.gio_descriptors[field_idx].has_levels) return 0;
    // if (gio_param->data.gio_descriptors[field_idx].has_levels) return 0;
    if (gio_param->data.gio_descriptors[field_idx].has_levels) return 1;

    /* skip time variable */
    if (gio_param->file.gio_files[ifile].is_time_field[field_idx]) return 0;

    if (gio_param->file.gio_files[ifile].is_grid[ifield]) return 1;

    return 0;
}

/*---< gio_allocate_req_ids() >-----------------------------------------------*/
/* allocate nonblocking request ID array req_ids(:) */
void gio_allocate_req_ids(gio_parameters *gio_param,
                          int             ifile,
                          int             is_grid,
                          int             ifield)
{
    int i, idx_start, idx_end, num_reqs;

    if (ifield == -1) {  /* find the number of requests for all fields */
        idx_start = 0;
        idx_end   = gio_param->file.gio_files[ifile].nflds;
    } else {             /* find the number of requests for this field only */
        idx_start = ifield;
        idx_end   = ifield + 1;
    }

    num_reqs = 1;  /* at least one, the time stamp written by root process */

    /* accumulate the number of requests for all fields */
    for (i=idx_start; i<idx_end; i++) {
        if (  is_grid && ! gio_is_output_grid(gio_param, ifile, i)) continue;
        if (! is_grid && ! gio_is_output_field(gio_param, ifile, i)) continue;
        int field_idx = gio_param->file.gio_files[ifile].fields[i];
        // if (gio_param->data.gio_descriptors[field_idx].is_level_data) continue;
        num_reqs += gio_get_num_w_reqs(gio_param, field_idx);
    }

    /* allocate the request ID array */
    gio_param->req_ids = calloc_1D_int(num_reqs);
    gio_param->num_reqs = 0;
}

/*----< gio_deallocate_io_buffer() >------------------------------------------*/
/* deallocate all grid fields' I/O buffers used in nonblocking I/O */
void gio_deallocate_io_buffer(gio_parameters *gio_param,
                              int             ifile,
                              int             is_grid,
                              int             ifield)
{
    int idx_start, idx_end, i, field_idx;
    file_descriptor *gio_files = gio_param->file.gio_files;

    if (ifield == -1) {  /* for all fields */
        idx_start = 0;
        idx_end   = gio_files[ifile].nflds;
    } else {             /* for this field only */
        idx_start = ifield;
        idx_end   = ifield + 1;
    }

    for (i=idx_start; i<idx_end; i++) {
        data_descriptor *dscrpt;

        if (  is_grid && ! gio_is_output_grid (gio_param, ifile, i)) continue;
        if (! is_grid && ! gio_is_output_field(gio_param, ifile, i)) continue;

        field_idx = gio_files[ifile].fields[i];
        dscrpt    = &gio_param->data.gio_descriptors[field_idx];
        if (dscrpt->is_level_data) continue;

        /* deallocate buffers used by the nonblocking I/O */
        if (dscrpt->io_buf != NULL) {
            tfree(dscrpt->io_buf);
            dscrpt->io_buf = NULL;
        }
    }
}

/*---< gio_write_verticalgrid() >---------------------------------------------*/
/* only root process write the vertical grid variables */
void gio_write_verticalgrid(gio_parameters *gio_param,
                            int             ifile,
                            int             field_idx)
{
    char *fldname, err_msg[128];
    int err, fieldID, dimid;
    MPI_Offset start, count;
    data_descriptor *dscrpt    = &gio_param->data.gio_descriptors[field_idx];
    file_descriptor *gio_files = gio_param->file.gio_files;

    if (gio_param->gio_me > 0) return;

    fldname = dscrpt->output_field_name;

    err = gio_inq_varid(gio_param, gio_files[ifile].fileID, fldname,
                            &fieldID);
    if (err != NC_NOERR) {
        sprintf(err_msg, "getting varid: %s\n",fldname);
        NC_CHECK(err, err_msg, 0)
    }

    err = ncmpi_inq_vardimid(gio_files[ifile].fileID, fieldID, &dimid);
    if (err != NC_NOERR) {
        sprintf(err_msg, "getting dimid: %s\n",fldname);
        NC_CHECK(err, err_msg, 0)
    }

    err = ncmpi_inq_dimlen(gio_files[ifile].fileID, dimid, &count);
    if (err != NC_NOERR) {
        sprintf(err_msg, "getting dimlen: %s\n",fldname);
        NC_CHECK(err, err_msg, 0)
    }
    start = 0;

    /* Since vertical grid is small, let root process write the entire
       field array */
    if (gio_param->gio_iotype == MODE_BLOCKING_INDEP)
        err = gio_put_real(gio_param,
                           gio_files[ifile].fileID,
                           fieldID,
                           &start, &count,
                           dscrpt->float_level_data);
    else
        err = gio_iput_real(gio_param,
                            gio_files[ifile].fileID,
                            fieldID,
                            &start, &count,
                            dscrpt->float_level_data);
    if (err != NC_NOERR) {
        sprintf(err_msg, "gio_iput/put_real: %s\n",fldname);
        NC_CHECK(err, err_msg, 0)
    }

    gio_param->bytes_API_write += count*4;

    /* Update stats specific to this file */
    dscrpt->record_count++;
    dscrpt->bytes_written += count*4;

}

/*---< gio_write_field() >----------------------------------------------------*/
/* Write all blocks of a field */
void gio_write_field(gio_parameters *gio_param,
                     int             ifile,
                     int             field_idx)
{
    int no_poles, psize, allocate_size, one_block_size;
    int i, buf_off, idx, idx_min, idx_max, isp, imin4, imax4;
    MPI_Offset offset, record_size, blen;

    gio_grid        *grid   = &gio_param->grid;
    data_descriptor *dscrpt = &gio_param->data.gio_descriptors[field_idx];

    /* post the nonblocking writes for vertical grid variables */
    if (dscrpt->is_level_data) {
        gio_write_verticalgrid(gio_param, ifile, field_idx);
        return;
    }

    /* calculate the total number of write requests, one write per block */
    psize          = dscrpt->column_size;
    one_block_size = psize * grid->block_size * grid->block_size;

    /* each field has this many writes (extra space for poles) */
    allocate_size = grid->num_grid_blocks * one_block_size + psize * 2;

    gio_set_limits(gio_param, field_idx, &idx_min, &idx_max);
    imin4 = dscrpt->imin4;
    imax4 = dscrpt->imax4;

    allocate_size *= (idx_max - idx_min + 1) * (imax4 - imin4 + 1);
    dscrpt->allocate_size = allocate_size;

    /* allocate write buffers */
         if (dscrpt->data_type == flt_type) allocate_size *= sizeof(float);
    else if (dscrpt->data_type == int_type) allocate_size *= sizeof(int);
    else if (dscrpt->data_type == dbl_type) allocate_size *= sizeof(double);

    dscrpt->io_buf = tmalloc(allocate_size);

    /* Loop over time integration index (if applicable) */

    /* wkliao: Q. when idx > 1, where are the poles stored in the nc variable? */
    no_poles    = dscrpt->no_poles;
    record_size = dscrpt->grid_size * psize;
    buf_off = 0;

    for (isp=imin4; isp<=imax4; isp++) {
        for (idx=idx_min; idx<=idx_max; idx++) {
            /* pack and write pole data */
            blen = psize;
            if (! no_poles) {
                if (grid->north_pole_flag) {
                    /* this proc has north pole */
                    // offset = record_size * (idx - idx_min);
                    offset = 0;
                    gio_write_block(gio_param, ifile, field_idx, 0, buf_off, 'n',
                                    offset, blen, idx, isp);
                    buf_off += psize;
                }
                else
                    gio_write_block(gio_param, ifile, field_idx, 0, 0, 'n',
                                    0, 0, 0, 0);

                if (grid->south_pole_flag) {
                    /* this proc has south pole */
                    // offset = record_size * (idx - idx_min) + psize;
                    offset = psize;
                    gio_write_block(gio_param, ifile, field_idx, 0, buf_off, 's',
                                    offset, blen, idx, isp);
                    buf_off += psize;
                }
                else
                    gio_write_block(gio_param, ifile, field_idx, 0, 0, 's',
                                    0, 0, 0, 0);
            }

            /* pack and write the cell blocks */
            blen = one_block_size;
            for (i=0; i<grid->num_grid_blocks; i++) {
                /* write one block at a time */

                /* avoid integer overflow */
                offset = grid->gio_grid_blocks[i].block_index;
                offset *= one_block_size;
                if (! no_poles) offset += psize * 2;
                // offset = offset + record_size * (idx - 1)
    
                gio_write_block(gio_param, ifile, field_idx, i, buf_off, 'x',
                                offset, blen, idx, isp);
                buf_off += one_block_size;
            }
        }
    }

    gio_reset_average(gio_param, field_idx);
}

/*----< gio_nonblocking_io_wait() >-------------------------------------------*/
/* wait for all nonblocking I/O requests to complete */
void gio_nonblocking_io_wait(gio_parameters *gio_param,
                             int             ifile)
{
    int i, err, *req_sts=NULL;
    double starttime;

    if (gio_param->gio_iotype == MODE_BLOCKING_INDEP) return;

    starttime = MPI_Wtime();
    if (gio_param->num_reqs > 0)
        req_sts = (int*) malloc(gio_param->num_reqs * sizeof(int));

    if (gio_param->gio_iotype == MODE_NONBLOCKING_COLL ||
        gio_param->gio_iotype == MODE_BLOCKING_COLL)
        err = ncmpi_wait_all(gio_param->file.gio_files[ifile].fileID,
                             gio_param->num_reqs,
                             gio_param->req_ids,
                             req_sts);
    else
        err = ncmpi_wait    (gio_param->file.gio_files[ifile].fileID,
                             gio_param->num_reqs,
                             gio_param->req_ids,
                             req_sts);
    if (err != NC_NOERR) {
        printf("ERROR: rank=%d at FILE=%s LINE=%d ncmpi_wait(_all): %s\n",
               gio_param->gio_me, __FILE__, __LINE__, ncmpi_strerror(err));
        return;
    }
    /* check each nonblocking I/O status for errors */
    for (i=0; i<gio_param->num_reqs; i++) {
        if (req_sts[i] != NC_NOERR) {
            printf("ERROR: rank=%d at FILE=%s LINE=%d ncmpi_wait(_all) status[%d]: %s\n",
                   gio_param->gio_me, __FILE__, __LINE__, i, ncmpi_strerror(req_sts[i]));
        }
    }

    if (gio_param->num_reqs > 0)
        free(req_sts);
/*
    if (gio_param->req_ids != NULL)
        tfree(gio_param->req_ids);
    gio_param->req_ids = NULL;
*/

    /* because all num_reqs are committed, reset num_reqs to 0, if no error */
    gio_param->num_reqs = 0;

    gio_param->time_in_nf_wait += MPI_Wtime() - starttime;
}

/*----< gio_free_req_ids() >-------------------------------------------------*/
/* wait for all nonblocking I/O requests to complete */
void gio_free_req_ids(gio_parameters *gio_param)
{
    if (gio_param->gio_iotype == MODE_BLOCKING_INDEP) return;

    if (gio_param->req_ids != NULL)
        tfree(gio_param->req_ids);
    gio_param->req_ids = NULL;
}

/*----< gio_write_block() >---------------------------------------------------*/
/* Post a PnetCDF nonblocking write for a block */
static
void gio_write_block(gio_parameters *gio_param,
                     int             ifile,
                     int             field_idx,
                     int             iblock,
                     int             buf_off,
                     char            pole,
                     MPI_Offset      offset,
                     MPI_Offset      blen,
                     int             idx,
                     int             isp)
{
    data_descriptor *dscrpt = &gio_param->data.gio_descriptors[field_idx];

    if (blen == 0) {
        gio_write(gio_param, field_idx, ifile, NULL, 0, 0, 0);
        return;
    }

    /* copy a block to the I/O buffer and post the non-blocking write call */
    if (dscrpt->data_type == flt_type) {
        float *flt_buf = (float*)dscrpt->io_buf;
        if (dscrpt->average_data)
            gio_eval_average(gio_param, field_idx, pole, iblock, idx, isp,
                             flt_buf+buf_off);
        else
            gio_copy_to_buffer(gio_param, field_idx, pole, iblock, idx, isp,
                               flt_buf+buf_off);

        gio_write(gio_param, field_idx, ifile, flt_buf+buf_off,
                  offset, idx, blen);
    }
    else if (dscrpt->data_type == int_type) {
        int *int_buf = (int*)dscrpt->io_buf;
        gio_copy_to_buffer(gio_param, field_idx, pole, iblock, idx, isp,
                           int_buf+buf_off);

        gio_write(gio_param, field_idx, ifile, int_buf+buf_off,
                  offset, idx, blen);
    }
    else if (dscrpt->data_type == dbl_type) {
        double *dbl_buf = (double*)dscrpt->io_buf;
        gio_copy_to_buffer(gio_param, field_idx, pole, iblock, idx, isp,
                           dbl_buf+buf_off);

        gio_write(gio_param, field_idx, ifile, dbl_buf+buf_off,
                  offset, idx, blen);
    }
}
