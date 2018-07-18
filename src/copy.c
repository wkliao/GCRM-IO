/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: copy.c 1858 2013-06-12 04:49:45Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>   /* access() */
#include <sys/time.h> /* gettimeofday() */
#include <time.h>     /* strftime() */

#include <mpi.h>
#include <pnetcdf.h>

#include "gcrm.h"
#include "gio.h"

/*----< COPY_LOOP() >---------------------------------------------------------*/
#define COPY_LOOP(dest_buf, src_buf) {                                      \
    int k, n, m;                                                            \
    for (k=0; k<klen; k++) {                                                \
        for (n=0; n<nlen; n++) {                                            \
            for (m=0; m<mlen; m++) {                                        \
                dest_buf[idx] = src_buf[ldx]; /* implicit type convert */   \
                ldx += inc2;                                                \
                idx++;                                                      \
            }                                                               \
            ldx = ldx - mlen*inc2 + ld2*inc1;                               \
            idx = idx - mlen + stride2_stride3;                             \
        }                                                                   \
        ldx = ldx - nlen*ld2*inc1 + ld1_ld2_istride_jstride;                \
        idx = idx - nlen*stride2_stride3 + stride2;                         \
    }                                                                       \
}

#define COPY_BUFFER(type, nbuf, sbuf, buf) {                            \
    type *dest = buffer;                                                \
    if (write_north_pole)                                               \
        COPY_LOOP(dest, nbuf)                                           \
    else if (write_south_pole)                                          \
        COPY_LOOP(dest, sbuf)                                           \
    else {                                                              \
        int i, j;                                                       \
        for (j=0; j<jlen; j++) {                                        \
            for (i=0; i<ilen; i++) {                                    \
                int mdx = gio_param->grid.grid_map[j][i];               \
                idx = mdx*stride1_stride2_stride3;                      \
                COPY_LOOP(dest, buf);                                   \
                ldx = ldx - klen*ld1_ld2_istride_jstride + ld1_ld2;     \
            }                                                           \
            ldx = ldx - ilen*ld1_ld2 + ld1_ld2_istride;                 \
        }                                                               \
    }                                                                   \
}

/*----< gio_copy_to_buffer() >------------------------------------------------*/
void gio_copy_to_buffer(gio_parameters *gio_param,
                        int             field_idx,
                        char            pole,
                        int             iblock,
                        int             time_idx,
                        int             spc_idx,
                        void           *buffer)
{
    /* Copy grid data to IO buffer (size of one block) */
    int idx, ldx, res;
    int istride, jstride, ilen, jlen;
    int inc1, inc2, inc3;
    int ld2, stride2;
    int ld1_ld2, ld1_ld2_istride, ld1_ld2_istride_jstride;
    int stride2_stride3, stride1_stride2_stride3;
    int klen, nlen, mlen;
    int write_north_pole, write_south_pole, write_data_block;
    int record_offset, species_offset;
    int nvals;
    int imin1, imax1, imin2, imax2, imin3, imax3, imin4, imax4;
    int ld1, ld3, stride1, stride3;
    int no_poles;
    double starttime = MPI_Wtime();

    gio_dim         *dim    = &gio_param->dim;
    gio_grid        *grid   = &gio_param->grid;
    data_descriptor *dscptr = &gio_param->data.gio_descriptors[field_idx];
    data_block      *blk    = &dscptr->gio_data_blocks[iblock];
    elm_type         data_type = dscptr->data_type;

    no_poles         = dscptr->no_poles;
    write_north_pole = 0;
    write_south_pole = 0;
    write_data_block = 0;

    if (pole == 'N' || pole == 'n') {
        if (grid->north_pole_flag && ! no_poles) {
            write_north_pole = 1;
            if ((data_type == int_type && dscptr->int_north_data    == NULL) ||
                (data_type != int_type && dscptr->double_north_data == NULL)) {
                printf("Failure in gio_copy_to_buffer for %s North pole data not assigned\n",
                       dscptr->field_name);
                ABORT
            }
        }
    }
    else if (pole == 'S' || pole == 's') {
        if (grid->south_pole_flag && ! no_poles) {
            write_south_pole = 1;
            if ((data_type == int_type && dscptr->int_south_data    == NULL) ||
                (data_type != int_type && dscptr->double_south_data == NULL)) {
                printf("Failure in gio_copy_to_buffer for %s South pole data not assigned\n",
                       dscptr->field_name);
                ABORT
            }
        }
    }
    else {
        write_data_block = 1;
        if ((data_type == int_type && blk->integer_data == NULL) ||
            (data_type != int_type && blk->double_data  == NULL)) {
            printf("Failure in gio_copy_to_buffer for %s Data not assigned for block %d\n",
                   dscptr->field_name, iblock);
            ABORT
        }
        ilen = blk->imax - blk->imin + 1;
        jlen = blk->jmax - blk->jmin + 1;

        nvals = grid->block_size;
        res = 0;
        while (nvals > 1) {
            res++;
            nvals /= 2;
        }
    }

    /* Start setting up some general parameters */

    istride = grid->gio_istride;
    jstride = grid->gio_jstride;
    imin1   = dscptr->imin1;
    imax1   = dscptr->imax1;
    inc1    = dscptr->inc1;
    ld1     = dscptr->ld1;
    imin2   = dscptr->imin2;
    imax2   = dscptr->imax2;
    inc2    = dscptr->inc2;
    ld2     = dscptr->ld2;
    imin3   = dscptr->imin3;
    imax3   = dscptr->imax3;
    inc3    = dscptr->inc3;
    ld3     = dscptr->ld3;
    imin4   = dscptr->imin4;
    imax4   = dscptr->imax4;
    stride1 = dscptr->stride1;
    stride2 = dscptr->stride2;
    stride3 = dscptr->stride3;

    klen  = (imax3-imin3+1) / inc3;
    if ((imax3-imin3+1) % inc3) klen++;
    nlen  = (imax1-imin1+1) / inc1;
    if ((imax1-imin1+1) % inc1) nlen++;
    mlen  = (imax2-imin2+1) / inc2;
    if ((imax2-imin2+1) % inc2) mlen++;

    ld1_ld2                 = ld1*ld2;
    ld1_ld2_istride         = ld1_ld2*istride;
    ld1_ld2_istride_jstride = ld1_ld2_istride*jstride;
            stride2_stride3 =         stride2*stride3;
    stride1_stride2_stride3 = stride1*stride2_stride3;

    /* Evaluate record_offset for variables with time integration values */

    record_offset = ld1_ld2_istride_jstride * ld3;
    if (dscptr->index_choice >= 0)
        species_offset = record_offset*dim->index_sizes[dscptr->index_choice];
    else
        species_offset = record_offset;

    idx = 0;
    ldx = imin2-1 + (imin1-1)*ld2 + (imin3-1)*ld1_ld2_istride_jstride
        + time_idx*record_offset + (spc_idx-1)*species_offset;

    if (data_type == int_type)
        COPY_BUFFER(int, dscptr->int_north_data,
                         dscptr->int_south_data,
                         blk->integer_data)
    else if (data_type == flt_type)
        COPY_BUFFER(float, dscptr->double_north_data,
                           dscptr->double_south_data,
                           blk->double_data)
    else if (data_type == dbl_type)
        COPY_BUFFER(double, dscptr->double_north_data,
                            dscptr->double_south_data,
                            blk->double_data)

    gio_param->time_in_API_copy += MPI_Wtime() - starttime;
}



