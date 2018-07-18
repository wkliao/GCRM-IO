/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio_API_statistics.c 1858 2013-06-12 04:49:45Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include "gio.h"

#define ALLREDUCE_SUM(buf,type) MPI_Allreduce(MPI_IN_PLACE, buf, 1, type, MPI_SUM, MPI_COMM_WORLD);

/*----< gio_term_statistics() >-----------------------------------------------*/
void gio_term_statistics(gio_parameters *gio_param)
{
    int rank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double d_tmp[18], d_sum[18];
    d_tmp[0]  = gio_param->total_time_in_API;
    d_tmp[1]  = gio_param->time_in_nf_put_var;
    d_tmp[2]  = gio_param->time_in_nf_put_var_grid;
    d_tmp[3]  = gio_param->time_in_nf_get_var;
    d_tmp[4]  = gio_param->time_in_nf_create;
    d_tmp[5]  = gio_param->time_in_nf_open;
    d_tmp[6]  = gio_param->time_in_nf_close;
    d_tmp[7]  = gio_param->time_in_nf_inq_varid;
    d_tmp[8]  = gio_param->time_in_nf_inq_dimlen;
    d_tmp[9]  = gio_param->time_in_nf_put_att;
    d_tmp[10] = gio_param->time_in_nf_get_att;
    d_tmp[11] = gio_param->time_in_nf_def_dim;
    d_tmp[12] = gio_param->time_in_nf_def_var;
    d_tmp[13] = gio_param->time_in_nf_enddef;
    d_tmp[14] = gio_param->time_in_update_time;
    d_tmp[15] = gio_param->time_in_nf_iput;
    d_tmp[16] = gio_param->time_in_nf_iget;
    d_tmp[17] = gio_param->time_in_nf_wait;
    MPI_Reduce(d_tmp, d_sum, 18, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    gio_param->total_time_in_API       = d_sum[0];
    gio_param->time_in_nf_put_var      = d_sum[1];
    gio_param->time_in_nf_put_var_grid = d_sum[2];
    gio_param->time_in_nf_get_var      = d_sum[3];
    gio_param->time_in_nf_create       = d_sum[4];
    gio_param->time_in_nf_open         = d_sum[5];
    gio_param->time_in_nf_close        = d_sum[6];
    gio_param->time_in_nf_inq_varid    = d_sum[7];
    gio_param->time_in_nf_inq_dimlen   = d_sum[8];
    gio_param->time_in_nf_put_att      = d_sum[9];
    gio_param->time_in_nf_get_att      = d_sum[10];
    gio_param->time_in_nf_def_dim      = d_sum[11];
    gio_param->time_in_nf_def_var      = d_sum[12];
    gio_param->time_in_nf_enddef       = d_sum[13];
    gio_param->time_in_update_time     = d_sum[14];
    gio_param->time_in_nf_iput         = d_sum[15];
    gio_param->time_in_nf_iget         = d_sum[16];
    gio_param->time_in_nf_wait         = d_sum[17];

    MPI_Offset ll_tmp[2], ll_sum[2];
    ll_tmp[0] = gio_param->bytes_API_write;
    ll_tmp[1] = gio_param->bytes_API_read;
    MPI_Reduce(ll_tmp, ll_sum, 2, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
    gio_param->bytes_API_write = ll_sum[0];
    gio_param->bytes_API_read  = ll_sum[1];

    if (rank > 0) return;

    printf("GIO: ----------------------------------------------------------\n");
    printf("GIO statistics: (cumulative across all MPI processes)\n");
#define PRN_TIME(timing) printf("GIO: %-32s %20.4f\n", #timing, gio_param->timing);
    PRN_TIME(total_time_in_API)
    PRN_TIME(time_in_nf_put_var)
    PRN_TIME(time_in_nf_put_var_grid)
    PRN_TIME(time_in_nf_put_att)
    PRN_TIME(time_in_nf_def_dim)
    PRN_TIME(time_in_nf_def_var)
    PRN_TIME(time_in_update_time)
    PRN_TIME(time_in_nf_create)
    PRN_TIME(time_in_nf_open)
    PRN_TIME(time_in_nf_close)
    PRN_TIME(time_in_nf_enddef)
    PRN_TIME(time_in_nf_inq_varid)
    PRN_TIME(time_in_nf_inq_dimlen)
    PRN_TIME(time_in_nf_iput)
    PRN_TIME(time_in_nf_iget)
    PRN_TIME(time_in_nf_wait)
    PRN_TIME(time_in_API_copy)
    PRN_TIME(time_in_avgs)
    printf("GIO: ----------------------------------------------------------\n");

    double wr_amount = gio_param->bytes_API_write;
    printf("GIO: bytes_API_write (Bytes):         %20.0f\n", wr_amount);
    wr_amount /= 1048576;
    printf("GIO: bytes_API_write (MiB):           %20.4f\n", wr_amount);
    wr_amount /= 1024;
    printf("GIO: bytes_API_write (GiB):           %20.4f\n", wr_amount);

    if (gio_param->total_time_in_API > 0.0) {
        double wr_bandwidth = gio_param->bytes_API_write;
        wr_bandwidth /= gio_param->total_time_in_API/nprocs;
        printf("GIO: bandwidth for writes (MiB/sec):  %20.4f\n", wr_bandwidth/1048576);
        printf("GIO: bandwidth for writes (GiB/sec):  %20.4f\n", wr_bandwidth/1073741824);
    }
    else
        printf("GIO: total_time_in_API =              %20.4f\n",0.0);
    printf("GIO: ----------------------------------------------------------\n");
}

