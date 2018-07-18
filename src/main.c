/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "gcrm.h"
#include "gio.h"
#include "grid_params.h"

/*----< print_timing_result() >-----------------------------------------------*/
void print_timing_result(gcrm_parameters *gcrm_param,
                         gio_parameters  *gio_param,
                         int              num_dumps,
                         double           init_timing,
                         double           comp_timing,
                         double           io_timing,
                         double           final_timing)
{
    int rank, nprocs;
    double  maxT[4], d_tmp[4];
    MPI_Offset w_size, r_size, io_amount;
    MODULE_grid_params *grid_params = &gcrm_param->grid_params;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* get the max timing among all processes */
    d_tmp[0] = init_timing;
    d_tmp[1] = comp_timing;
    d_tmp[2] = io_timing;
    d_tmp[3] = final_timing;
    MPI_Reduce(d_tmp, maxT, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    init_timing  = maxT[0];
    comp_timing  = maxT[1];
    io_timing    = maxT[2];
    final_timing = maxT[3];

    /* get the malloc usage among all processes */
    MPI_Offset usage = get_max_mem_alloc();
    MPI_Offset max_usage, min_usage, avg_usage;
    MPI_Reduce(&usage, &max_usage, 1, MPI_OFFSET, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&usage, &min_usage, 1, MPI_OFFSET, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&usage, &avg_usage, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&gio_param->pnc_write_amount, &w_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&gio_param->pnc_read_amount,  &r_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
    io_amount = w_size + r_size;

    /* report performance statistics results */
    if (rank == 0 && io_timing > 0.0) {
        int    i, j, num_files, num_grid_vars=0, num_fld_vars=0;
        char   *phys_str = "disabled";
        double io_size, bandwidth;

        /* count the number of files, grid variables, field variables written */
        num_files = gio_param->file.num_files;
        for (i=0; i<gio_param->file.num_files; i++) {
            if (! gio_param->file.gio_files[i].complete) continue;
            if (  gio_param->file.gio_files[i].grid_ONLY) continue;
            if (  gio_param->file.gio_files[i].frequency == 0) {
                num_files--;
                continue;
            }

            for (j=0; j<gio_param->file.gio_files[i].nflds; j++)
                if (gio_is_output_field(gio_param, i, j))
                    num_fld_vars++;
        }
        for (i=0; i<MAX_DATA_FIELDS; i++)
            if (gio_param->data.gio_descriptors[i].is_grid)
                num_grid_vars++;

        printf("====  %s   released on  RELEASE_DATE ====\n", PACKAGE_STRING);
        printf("----  Run-time parameters ---------------------------------\n");
        printf("----  Number of processes                            = %d\n", nprocs);
        printf("----  level_max  (global horizontal grid resolution) = %d\n", grid_params->level_max);
        printf("----  sbdmn_iota (see grid_params.h)                 = %d\n", grid_params->sbdmn_iota);
        printf("----  level_glbl (see grid_params.h)                 = %d\n", grid_params->level_glbl);
        printf("----  km         (number of vertical layers)         = %d\n", gcrm_param->vertical.km);
        printf("----  cell_max   (global number of cells)            = %d\n", grid_params->cell_max);
        printf("----  im         (local  number of cells along i)    = %d\n", grid_params->im);
        printf("----  jm         (local  number of cells along j)    = %d\n", grid_params->jm);
        printf("----  nsdm_glbl  (global number of blocks)           = %d\n", grid_params->nsdm_glbl);
        printf("----  nsdm       (local  number of blocks)           = %d\n", grid_params->nsdm);
        if (gcrm_param->enable_physics) phys_str = "enabled";
        printf("----  using physics variables                        = %s\n", phys_str);
        printf("----  number of files written                        = %d\n", num_files);
        printf("----  number of grid  variables written              = %d\n", num_grid_vars);
        printf("----  number of field variables written              = %d\n", num_fld_vars);
        printf("----  number of snapshot dumps                       = %d\n", num_dumps);
        printf("----  I/O method                                     = %s\n", gcrm_param->zgrd_output.iotype);
        printf("----  cdf_output_path                                = %s\n", gcrm_param->zgrd_output.cdf_output_path);
        // printf("----  number of snapshot dumps                       = %d\n", gio_param->num_dumps);
        // printf("----  number of snapshot files                       = %d\n", gio_param->file.num_files);
        printf("---------------------------------------------------------------\n");
        printf("Timing results (max among all processes)\n");
        printf("init     time=%17.2f sec\n", init_timing);
        printf("Max comp time=%17.2f sec\n", comp_timing);
        printf("Max I/O  time=%17.2f sec\n", io_timing);
        printf("finalize time=%17.2f sec\n", final_timing);
        printf("---------------------------------------------------------------\n");
        io_size = w_size;
        printf("Write  amount=%17.2f MB      = %14.2f MiB\n",       io_size/1000000.0,      io_size/1048576.0);
        io_size = r_size;
        printf("Read   amount=%17.2f MB      = %14.2f MiB\n",       io_size/1000000.0,      io_size/1048576.0);
        io_size = w_size + r_size;
        printf("I/O    amount=%17.2f MB      = %14.2f MiB\n",       io_size/1000000.0,      io_size/1048576.0);
        printf("             =%17.2f GB      = %14.2f GiB\n",       io_size/1000000000.0,   io_size/1073741824.0);
        bandwidth = io_size / io_timing;
        printf("I/O bandwidth=%17.2f MB/sec  = %14.2f MiB/sec\n", bandwidth/1000000.0,    bandwidth/1048576.0);
        printf("             =%17.2f GB/sec  = %14.2f GiB/sec\n", bandwidth/1000000000.0, bandwidth/1073741824.0);
        printf("---------------------------------------------------------------\n");
        printf("memory MAX usage (among %d procs) = %13.2f MiB\n", nprocs, (double)max_usage/1048576);
        printf("memory MIN usage (among %d procs) = %13.2f MiB\n", nprocs, (double)min_usage/1048576);
        printf("memory AVG usage (among %d procs) = %13.2f MiB\n", nprocs, (double)avg_usage/nprocs/1048576);
    }
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char   *in_fname;
    int     rank, nprocs, num_dumps=0, verbose=1;
    double  timing, init_timing, comp_timing=0, io_timing=0, final_timing=0;

    gcrm_parameters *gcrm_param;
    gio_parameters  *gio_param;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    debug = 0;

    if (argc == 1) /* default input parameter file is zgrd.in */
        in_fname = "zgrd.in";
    else if (argc == 2)
        in_fname = argv[1];
    else {
        if (!rank) printf("Usage: %s [parameter file]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    init_timing = MPI_Wtime();

    /* set malloc tracking counters to zeros */
    reset_mem_alloc();

    /* initialize all parameters set by GCRM, including
     *     reading input parameter file,
     *     memory allocations for all buffers, and
     *     initialize the contents for all grid variables
     */
    /* zero out the entire contents of gcrm_param and gio_param */
    gcrm_param = (gcrm_parameters*) tcalloc(1, sizeof(gcrm_parameters));
    gio_param  = (gio_parameters*)  tcalloc(1, sizeof(gio_parameters));

    gio_param->file.gio_files = (file_descriptor*) tcalloc(MAX_FILES, sizeof(file_descriptor));
    gio_param->data.gio_descriptors = (data_descriptor*) tcalloc(MAX_DATA_FIELDS, sizeof(data_descriptor));
int i; for (i=0; i<MAX_DATA_FIELDS; i++) gio_param->data.gio_descriptors[i].io_buf = NULL;
    gio_param->grid.gio_grid_blocks = (grid_block*) tcalloc(MAX_DATA_BLOCKS, sizeof(grid_block));

    gcrm_param_init(in_fname, gcrm_param, gio_param);

    /* initialize parameters for IO scheduling */
    gio_init_schedule(gio_param);
    init_timing = MPI_Wtime() - init_timing;

    comp_timing = io_timing = 0.0;
    int l_ZGrd = 1;
    int random_seed = 0;
    while (l_ZGrd) { /* the main loop */

        /* run dummy computation to emulate the compute cost */ 
        timing = MPI_Wtime();
        pseudo_computation(gcrm_param, random_seed);
        random_seed++;
        comp_timing += MPI_Wtime() - timing;

        /* advance the simulation time of amount idt_ZGrd[0] */
        gcrm_param->time.itime_ZGrd[0] = gcrm_param->time.idt_ZGrd[0]
                                       * (gcrm_param->time.step_count_ZGrd-1);
        gcrm_param->time.time_ZGrd = gcrm_param->time.itime_ZGrd[0]
                                   / gcrm_param->time.itime_ZGrd[1];

        if ((int)gcrm_param->time.time_ZGrd %
            gio_param->defaults.gio_default_frequency == 0)
            num_dumps++;

        if (gcrm_param->time.itime_ZGrd[0] >= gcrm_param->time.itime_end_ZGrd[0])
            l_ZGrd = 0;

        /* measure I/O cost */
        MPI_Barrier(MPI_COMM_WORLD);
        timing = MPI_Wtime();
        /* ZGrd_write_output() calls gio_driver(time_ZGrd)
         * checkpoint write for analysis data is carried out here */
        ZGrd_write_output(gio_param, gcrm_param->time.time_ZGrd);

#ifdef NOT_YET
        /* checkpoint write (in binary format) for restart */
        if (gcrm_param->l_use_binary_restart)
            restart(gcrm_param);
#endif
        io_timing += MPI_Wtime() - timing;

        gcrm_param->time.step_count_ZGrd++;

        if (verbose && rank == 1)
            printf("while loop num_dumps=%d time_ZGrd=%f gio_default_frequency=%d\n",num_dumps,gcrm_param->time.time_ZGrd,gio_param->defaults.gio_default_frequency);
    }

    /* close all files if there is any left opened */
    timing = MPI_Wtime();
    gio_terminate(gio_param);
    io_timing += MPI_Wtime() - timing;

    final_timing = MPI_Wtime();

    /* free the allocated memory space for all buffers */
    gcrm_param_finalize(gcrm_param, gio_param);

    /* print detailed I/O statistics on stdout */
    if (verbose) gio_term_statistics(gio_param);
    final_timing = MPI_Wtime() - final_timing;

    /* print the MPI-IO hints to stdout */
    if (verbose && gio_param->used_info != MPI_INFO_NULL) {
        if (gio_param->gio_me == 0)
            print_mpi_info(&gio_param->used_info);
        MPI_Info_free(&gio_param->used_info);
    }

    /* print high-level I/O and malloc statistics to stdout */
    print_timing_result(gcrm_param, gio_param, num_dumps,
                        init_timing, comp_timing, io_timing, final_timing);

    tfree(gio_param->file.gio_files);
    tfree(gio_param->data.gio_descriptors);
    tfree(gio_param->grid.gio_grid_blocks);
    tfree(gio_param);
    tfree(gcrm_param);

    /* check and report if any malloc has yet to be freed */
    MPI_Offset mem_alloc = get_mem_alloc();
    if (mem_alloc != 0)
        printf("memory allocated yet to be freed = %lld B\n", mem_alloc);
    check_mem_alloc();

    MPI_Finalize();
    return 0;
}

