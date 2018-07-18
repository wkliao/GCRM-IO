/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_output.c 4603 2017-12-07 07:12:52Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <errno.h>
extern int errno;

#include "ZGrd_output.h"
#include "util.h"


/*----< init_MODULE_ZGrd_output() >-------------------------------------------*/
/* from ZGrd_output.F90
 * Data updates:
 *     zgrd_output->descfile
 *     zgrd_output->outfileconfig
 *     zgrd_output->iotype
 *     zgrd_output->cdf_output_path
 *     zgrd_output->cdf_grid_option
 *     zgrd_output->cdf_sep_grid
 *     zgrd_output->cdf_output_freq
 *     zgrd_output->cdf_output_nsamples
 *     zgrd_output->cdf_restart_path
 *     zgrd_output->l_restart_overwrite
 *     zgrd_output->restart_output[0]
 *     zgrd_output->l_read_restart_tmp
 *     zgrd_output->l_use_binary_restart
 *     zgrd_output->restart_interval
 */
void init_MODULE_ZGrd_output(MODULE_ZGrd_output *zgrd_output)
{
    /* parameters from MODULE ZGrd_output */
    strcpy(zgrd_output->descfile,         "./ZGrd.desc");
    strcpy(zgrd_output->outfileconfig,    "./ZGrd.fcfg");
    strcpy(zgrd_output->iotype,           "nonblocking_collective");
    strcpy(zgrd_output->cdf_output_path,  "./");
    strcpy(zgrd_output->cdf_grid_option,  "all_files");
    strcpy(zgrd_output->cdf_sep_grid,     "grid_");
    zgrd_output->cdf_output_freq          = 600;
    zgrd_output->cdf_output_nsamples      = 8;
    strcpy(zgrd_output->cdf_restart_path, "./");

    zgrd_output->l_restart_overwrite  = 0;
    zgrd_output->restart_output[0]    = '\0';
    zgrd_output->l_read_restart_tmp   = 0;
    zgrd_output->l_use_binary_restart = 1;
    zgrd_output->restart_interval     = 0.0;
}


/*----< io_load_inputs() >----------------------------------------------------*/
/* from ZGrd_output.F90
 * Data updated:
 *     zgrd_output->iotype
 *     zgrd_output->descfile
 *     zgrd_output->outfileconfig
 *     zgrd_output->cdf_output_path
 *     zgrd_output->cdf_output_freq
 *     zgrd_output->cdf_output_nsamples
 *     zgrd_output->cdf_grid_option
 *     zgrd_output->cdf_sep_grid
 *     zgrd_output->l_use_binary_restart
 *     zgrd_output->l_read_restart_tmp
 *     zgrd_output->l_restart_overwrite
 *     zgrd_output->restart_interval
 *     zgrd_output->restart_input
 *     zgrd_output->restart_output
 */
void io_load_inputs(const char         *fname,
                    MODULE_ZGrd_output *zgrd_output)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        /* only root process reads from the file */

        char line[256], *value, *key;
        FILE *fp;

        OPEN_FILE(fname)

        while (fgets(line, 256, fp) != NULL) {
            GET_LINE_TOKEN_PAIR

            /* followings are from a call to io_load_inputs() */
            if (!strcasecmp(key, "iotype"))
                strcpy(zgrd_output->iotype, value);
            else if (!strcasecmp(key, "io_desc_file"))
                strcpy(zgrd_output->descfile, value);
            else if (!strcasecmp(key, "io_config_file"))
                strcpy(zgrd_output->outfileconfig, value);
            else if (!strcasecmp(key, "cdf_output_path"))
                strcpy(zgrd_output->cdf_output_path, value);
            else if (!strcasecmp(key, "cdf_output_frequency"))
                zgrd_output->cdf_output_freq = atoi(value);
            else if (!strcasecmp(key, "cdf_output_nsamples"))
                zgrd_output->cdf_output_nsamples = atoi(value);
            else if (!strcasecmp(key, "cdf_grid_option"))
                strcpy(zgrd_output->cdf_grid_option, value);
            else if (!strcasecmp(key, "cdf_sep_grid"))
                strcpy(zgrd_output->cdf_sep_grid, value);
            else if (!strcasecmp(key, "use_binary_restart")) {
                zgrd_output->l_use_binary_restart = 1;
                if (!strcasecmp(value, "false"))
                    zgrd_output->l_use_binary_restart = 0;
            }
            else if (!strcasecmp(key, "restart_interval"))
                zgrd_output->restart_interval = atof(value);
            else if (!strcasecmp(key, "read_restart")) {
                zgrd_output->l_read_restart_tmp = 0;
                if (!strcasecmp(value, "true"))
                    zgrd_output->l_read_restart_tmp = 1;
            }
            else if (!strcasecmp(key, "restart_input"))
                strcpy(zgrd_output->restart_input, value);
            else if (!strcasecmp(key, "restart_output"))
                strcpy(zgrd_output->restart_output, value);
            else if (!strcasecmp(key, "restart_overwrite")) {
                zgrd_output->l_restart_overwrite = 0;
                if (!strcasecmp(value, "true"))
                    zgrd_output->l_restart_overwrite = 1;
            }
        }
        fclose(fp);

        char *env_str;
        if ((env_str = getenv("GCRM_OUTPUT_PATH")) != NULL)
            strcpy(zgrd_output->cdf_output_path, env_str);
    }
    MPI_Bcast(zgrd_output, sizeof(MODULE_ZGrd_output), MPI_BYTE, 0, MPI_COMM_WORLD);
}

#include "gio.h"

/*----< ZGrd_write_output() >-------------------------------------------------*/
/* from ZGrd_output.F90 */
void ZGrd_write_output(gio_parameters *gio_param,
                       double          time_ZGrd)
{
    gio_driver(gio_param, time_ZGrd);
}

