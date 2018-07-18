/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: average.c 4607 2017-12-07 07:21:49Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <gcrm.h>
#include <gio.h>
#include <util.h>

/* In average.F90 */

/*----< gio_allct_avg_descriptor() >------------------------------------------*/
/* Allocate memory for descriptors that will contain averaged data and copy
 * pointers from parent descriptors. */
void gio_allct_avg_descriptor(gio_parameters *gio_param,
                              int             field_idx)
{
    int i, buf_size, pole_size;

    data_descriptor *dscrpt = &gio_param->data.gio_descriptors[field_idx];
    data_descriptor *parent = &gio_param->data.gio_descriptors[dscrpt->parent];
    /* This indirect way of allocating and assigning memory is a workaround
       for IBM compilers */

    /* allocate space in averaging buffers for data and initialize buffers */
    dscrpt->num_data_blocks = parent->num_data_blocks;

    for (i=0; i<dscrpt->num_data_blocks; i++) {
        dscrpt->gio_data_blocks[i].imin = parent->gio_data_blocks[i].imin;
        dscrpt->gio_data_blocks[i].imax = parent->gio_data_blocks[i].imax;
        dscrpt->gio_data_blocks[i].jmin = parent->gio_data_blocks[i].jmin;
        dscrpt->gio_data_blocks[i].jmax = parent->gio_data_blocks[i].jmax;

        dscrpt->gio_average_data[i].imin = parent->gio_data_blocks[i].imin;
        dscrpt->gio_average_data[i].imax = parent->gio_data_blocks[i].imax;
        dscrpt->gio_average_data[i].jmin = parent->gio_data_blocks[i].jmin;
        dscrpt->gio_average_data[i].jmax = parent->gio_data_blocks[i].jmax;
    }

    buf_size   = dscrpt->ld1 * dscrpt->ld2 * dscrpt->ld3;
    pole_size  = buf_size;
    buf_size  *= gio_param->grid.gio_istride * gio_param->grid.gio_jstride;
    dscrpt->assigned = parent->assigned;

    for (i=0; i<dscrpt->num_data_blocks; i++) {
        dscrpt->gio_data_blocks[i].double_data =
        parent->gio_data_blocks[i].double_data;

        dscrpt->gio_average_data[i].real_data = calloc_1D_flt(buf_size);
    }

    if (! dscrpt->no_poles) {
        if (gio_param->grid.south_pole_flag) {
            dscrpt->double_south_data = parent->double_south_data;
            dscrpt->avg_south_data    = calloc_1D_flt(pole_size);
        }
        if (gio_param->grid.north_pole_flag) {
            dscrpt->double_north_data = parent->double_north_data;
            dscrpt->avg_north_data    = calloc_1D_flt(pole_size);
        }
    }
    dscrpt->avg_count = 0;
}

/*----< gio_free_avg_descriptor() >-------------------------------------------*/
void gio_free_avg_descriptor(gio_parameters *gio_param)
{
    int i, j;

    for (i=0; i<gio_param->data.num_averages; i++) {
        int field_idx = gio_param->data.avg_field_ids[i];
        data_descriptor *dscrpt = &gio_param->data.gio_descriptors[field_idx];

        for (j=0; j<dscrpt->num_data_blocks; j++) {
            if (dscrpt->gio_average_data[j].real_data == NULL)
                continue;
            free_1D_flt(dscrpt->gio_average_data[j].real_data);
            dscrpt->gio_average_data[j].real_data = NULL;
        }
        if (! dscrpt->no_poles) {
            if (gio_param->grid.south_pole_flag) {
                if (dscrpt->avg_south_data != NULL)
                    free_1D_flt(dscrpt->avg_south_data);
                dscrpt->avg_south_data = NULL;
            }
            if (gio_param->grid.north_pole_flag) {
                if (dscrpt->avg_north_data != NULL)
                    free_1D_flt(dscrpt->avg_north_data);
                dscrpt->avg_north_data = NULL;
            }
        }
        dscrpt->avg_count = 0;
    }
}


#define COPY_AVG {                                                      \
    for (k=imin3; k<=imax3; k+=inc3) {                                  \
        int q = p;                                                      \
        int icnt1 = 0;                                                  \
        for (n=imin1; n<=imax1; n+=inc1) {                              \
            idx = mdx + icnt3 + icnt1;                                  \
            ldx = q;                                                    \
            for (m=imin2; m<=imax2; m+=inc2) {                          \
                avg_data[idx++] = dbl_data[ldx];                        \
                ldx += inc2;                                            \
            }                                                           \
            icnt1 += stride2_stride3;                                   \
            q     += inc1*ld2;                                          \
        }                                                               \
        icnt3 += stride2;                                               \
        p     += inc3*ld1_ld2_istride_jstride;                          \
    }                                                                   \
}

/*----< gio_accumulate_field_data() >-----------------------------------------*/
/* Accumulate current values of field data into averaging buffer */
void gio_accumulate_field_data(gio_parameters *gio_param,
                               int             field_idx)
{
    double *dbl_data;
    float  *avg_data;
    int     ld1, ld2, ld3, stride1, stride2, stride3;
    int     imin1, imax1, inc1, imin2, imax2, inc2, imin3, imax3, inc3;
    int     istride, jstride, ilen, jlen;
    int     i, j, k, m, n, ldx, idx, imin, imax, jmin, jmax, nvals, res;
    int     nb;

    data_descriptor *dscrpt = &gio_param->data.gio_descriptors[field_idx];

    /* Start setting up some general parameters */
    istride = gio_param->grid.gio_istride;
    jstride = gio_param->grid.gio_jstride;

    imin1 = dscrpt->imin1;
    imax1 = dscrpt->imax1;
    inc1  = dscrpt->inc1;
    ld1   = dscrpt->ld1;
    imin2 = dscrpt->imin2;
    imax2 = dscrpt->imax2;
    inc2  = dscrpt->inc2;
    ld2   = dscrpt->ld2;
    imin3 = dscrpt->imin3;
    imax3 = dscrpt->imax3;
    inc3  = dscrpt->inc3;
    ld3   = dscrpt->ld3;
 
    stride1 = dscrpt->stride1;
    stride2 = dscrpt->stride2;
    stride3 = dscrpt->stride3;

    int ld1_ld2                 = ld1*ld2;
    int ld1_ld2_istride         = ld1_ld2*istride;
    int ld1_ld2_istride_jstride = ld1_ld2_istride*jstride;
    int stride2_stride3         = stride2*stride3;
    int stride1_stride2_stride3 = stride1*stride2_stride3;
    int u = (imin2-1) + (imin1-1)*ld2 + (imin3-1)*ld1_ld2_istride_jstride;
    /* Note u is the start index in C, translated from Fortran style */

    /* Accumulate north and south pole data */

    if (! dscrpt->no_poles) {
        if (gio_param->grid.north_pole_flag) {
            dbl_data = dscrpt->double_north_data;
            if (dbl_data == NULL)
                printf("%d: north pole dbl_data not associated for field: %s\n",
                       gio_param->gio_me, dscrpt->field_name);
            avg_data = dscrpt->avg_north_data;
            if (avg_data == NULL)
                printf("%d: north pole avg_data not associated for field: %s\n",
                       gio_param->gio_me, dscrpt->field_name);
            if (avg_data != NULL && dbl_data != NULL) {
                int icnt3=0, p=u, mdx=0;
                COPY_AVG
            }
        }
        if (gio_param->grid.south_pole_flag) {
            dbl_data = dscrpt->double_south_data;
            if (dbl_data == NULL)
                printf("%d: south pole dbl_data not associated for field: %s\n",
                       gio_param->gio_me, dscrpt->field_name);
            avg_data = dscrpt->avg_south_data;
            if (avg_data == NULL)
                printf("%d: south pole avg_data not associated for field: %s\n",
                       gio_param->gio_me, dscrpt->field_name);
            if (avg_data != NULL && dbl_data != NULL) {
                int icnt3=0, p=u, mdx=0;
                COPY_AVG
            }
        }
    }

    /* Accumulate remaining field data */
    nvals = gio_param->grid.block_size;
    res = 0;
    while (nvals > 1) {
        res++;
        nvals /= 2;
    }
    for (nb=0; nb<dscrpt->num_data_blocks; nb++) {
        imin = dscrpt->gio_data_blocks[nb].imin;
        imax = dscrpt->gio_data_blocks[nb].imax;
        jmin = dscrpt->gio_data_blocks[nb].jmin;
        jmax = dscrpt->gio_data_blocks[nb].jmax;
        ilen = imax - imin + 1;
        jlen = jmax - jmin + 1;
        dbl_data = dscrpt->gio_data_blocks[nb].double_data;
        if (dbl_data == NULL) {
            printf("%d: data blocks dbl_data not associated for field: %s\n",
                   gio_param->gio_me, dscrpt->field_name);
            continue;
        }
        avg_data = dscrpt->gio_average_data[nb].real_data;
        if (avg_data == NULL) {
            printf("%d: data blocks avg_data not associated for field: %s\n",
                   gio_param->gio_me, dscrpt->field_name);
            continue;
        }
        int jdx = 0;
        for (j=0; j<jlen; j++) {
            int kdx = jdx;
            for (i=0; i<ilen; i++) {
                int icnt3, mdx, p;
                mdx = gio_param->grid.grid_map[j][i] * stride1_stride2_stride3;
                p = u + kdx;
                icnt3 = 0;
                COPY_AVG
                kdx += ld1_ld2;
            }
            jdx += ld1_ld2_istride;
        }
    }

    /* Increment counter */
    dscrpt->avg_count++;
}

/*----< gio_accumulate_all_fields() >-----------------------------------------*/
/* Accumulate data for all fields */
void gio_accumulate_all_fields(gio_parameters *gio_param)
{
    int i;
    double starttime = MPI_Wtime();

    for (i=0; i<gio_param->data.num_fields; i++)
        if (gio_param->data.gio_descriptors[i].average_data)
            gio_accumulate_field_data(gio_param, i);
    gio_param->time_in_avgs += MPI_Wtime() - starttime;
}

/*----< gio_eval_average() >--------------------------------------------------*/
/* Copy averaging buffer to IO buffer and divide by the number of images in the
   accumulated data to get the average
*/
void gio_eval_average(gio_parameters *gio_param,
                      int             field_idx,
                      char            pole,
                      int             iblock,
                      int             time_idx,
                      int             species_idx,
                      float          *buffer)
{
    float *flt_buf, norm;
    int i, nsize, block_size, record_offset;
    double starttime = MPI_Wtime();

    data_descriptor *dscrpt = &gio_param->data.gio_descriptors[field_idx];

    if (dscrpt->avg_count <= 0) {
        /* TODO: Some kind of error condition */
        return;
    }
    block_size = gio_param->grid.block_size;
    norm = 1.0 / (float)(dscrpt->avg_count);

    /* Figure out how much data is in block.
       Start by setting up some general parameters */
    nsize = dscrpt->column_size;

    /* Figure out offset for variables with time integration index */

    record_offset = time_idx * nsize * block_size * block_size;
    if (pole == 'n' || pole == 'N')
        flt_buf = dscrpt->avg_north_data + record_offset;
    else if (pole == 's' || pole == 'S')
        flt_buf = dscrpt->avg_south_data + record_offset;
    else {
        nsize *= block_size * block_size;
        flt_buf = dscrpt->gio_average_data[iblock].real_data + record_offset;
    }
    for (i=0; i<nsize; i++)
        buffer[i] = flt_buf[i] * norm;

    gio_param->time_in_avgs += MPI_Wtime() - starttime;
}

/*----< gio_reset_average() >-------------------------------------------------*/
/* Reset counters and buffers on averaged fields to zero */
void gio_reset_average(gio_parameters *gio_param,
                       int             field_idx)
{
    int i, nsize, isize, imin, jmin, imax, jmax, ilen, jlen;
    double starttime = MPI_Wtime();

    data_descriptor *dscptr = &gio_param->data.gio_descriptors[field_idx];

    if (dscptr->average_data) {
        dscptr->avg_count = 0;

        /* Figure out how much data is in block.
           Start by setting up some general parameters */
        nsize = dscptr->column_size * sizeof(float);

        if (gio_param->grid.north_pole_flag)
            memset(dscptr->avg_north_data, 0, nsize);
        if (gio_param->grid.south_pole_flag)
            memset(dscptr->avg_south_data, 0, nsize);
        for (i=0; i<dscptr->num_data_blocks; i++) {
            imin = dscptr->gio_data_blocks[i].imin;
            imax = dscptr->gio_data_blocks[i].imax;
            jmin = dscptr->gio_data_blocks[i].jmin;
            jmax = dscptr->gio_data_blocks[i].jmax;
            ilen = imax - imin + 1;
            jlen = jmax - jmin + 1;
            isize = nsize * ilen * jlen;
            memset(dscptr->gio_average_data[i].real_data, 0, isize);
        }
    }
    gio_param->time_in_avgs += MPI_Wtime() - starttime;
}
