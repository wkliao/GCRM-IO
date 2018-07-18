/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_params.c 4603 2017-12-07 07:12:52Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <errno.h>
extern int errno;

#include "ZGrd_params.h"
#include "util.h"


/*----< ZGrd_load_inputs() >--------------------------------------------------*/
/* from ZGrd_params.F90
 * Data updates:
 *    zgrd_params->l_ZGrd
 *    zgrd_params->l_density
 *    zgrd_params->l_tracer
 *    zgrd_params->l_vorticity
 *    zgrd_params->l_divergence
 *    zgrd_params->l_theta
 *    zgrd_params->l_restart
 *    zgrd_params->l_timing
 *    zgrd_params->l_diffusion_div
 *    zgrd_params->l_diffusion_eta
 *    zgrd_params->l_diffusion_tht
 *    zgrd_params->l_eta_enstrophy_damping
 *    zgrd_params->l_verbose
 *    zgrd_params->l_report
 *    zgrd_params->equation_set_select
 *    zgrd_params->ZGrd_strng
 *    zgrd_params->initial_condition_select
 *    zgrd_params->path_output
 *    zgrd_params->path_restart
 *    // zgrd_params->path_data
 *    zgrd_params->file_restart
 *    zgrd_params->sfc_flux_method
 *    zgrd_params->k_diffusion_div
 *    zgrd_params->k_diffusion_tht
 *    zgrd_params->grid_point_path[0]
 *    zgrd_params->grid_point_select
 *    zgrd_params->io_mode
 *    zgrd_params->physics_mode
 *    zgrd_params->verbose_count
 *    zgrd_params->initial_condition_select
 */
void ZGrd_load_inputs(const char         *fname,
                      MODULE_ZGrd_params *zgrd_params)
{
    /* from ZGrd_params.F90 */
    /* set default values for input parameters define in MODULE ZGrd_params */
    zgrd_params->l_ZGrd       = 1; /* the model will continue to run */
    zgrd_params->l_density    = 1; /* timestep the continuity equation */
    zgrd_params->l_tracer     = 0; /* timestep the tracer     equation */
    zgrd_params->l_vorticity  = 1; /* timestep the vorticity  equation */
    zgrd_params->l_divergence = 1; /* timestep the divergence equation */
    zgrd_params->l_theta      = 1; /* timestep the potential temperature equation */

    strcpy(zgrd_params->equation_set_select, "unified"); /* baroclinic, anelastic, unified */
    strcpy(zgrd_params->ZGrd_strng, "ZGrd");
    zgrd_params->l_restart    = 0;
    zgrd_params->l_timing     = 1;

    zgrd_params->initial_condition_select = 6;
    strcpy(zgrd_params->path_output,  "../../output/ZGrd");
    strcpy(zgrd_params->path_restart, "../../output/ZGrd");
    // strcpy(zgrd_params->path_data,    "../../data");
    strcpy(zgrd_params->file_restart, "../../output/ZGrd/restart_ZGrd");

    strcpy(zgrd_params->sfc_flux_method, "Louis"); /* Louis or Deardorff */
    zgrd_params->l_diffusion_div = 1;
    zgrd_params->l_diffusion_eta = 1;
    zgrd_params->l_diffusion_tht = 0;
    zgrd_params->l_eta_enstrophy_damping = 0;
    zgrd_params->k_diffusion_div = 1.6E+17;
    zgrd_params->k_diffusion_tht = 1.6E+17;

    zgrd_params->grid_point_path[0] = '\0';
    strcpy(zgrd_params->grid_point_select, "bisect");
    strcpy(zgrd_params->io_mode, "parallel");
    strcpy(zgrd_params->physics_mode, "FULL PHYSICS");

    zgrd_params->l_verbose     = 0;
    zgrd_params->l_report      = 0;
    zgrd_params->verbose_count = 1;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) { /* only root process reads from the file */
        char line[256], *value, *key;
        FILE *fp;

        OPEN_FILE(fname)

        while (fgets(line, 256, fp) != NULL) {
            GET_LINE_TOKEN_PAIR

            /* followings are from a call to ZGRD_load_inputs() */
            if (!strcasecmp(key, "initial_condition"))
                zgrd_params->initial_condition_select = atoi(value);
            if (!strcasecmp(key, "output_path")) {
                int len = strlen(value);
                strcpy(zgrd_params->path_output, value);
                if (value[len-1] != '/') { /* add "/" at the end of path */
                    value[len] = '/';
                    value[len+1]   = '\0';
                }
            }
/* removed unused data_path as this code does not support reading grids etc 
            else if (!strcasecmp(key, "data_path")) {
                strcpy(zgrd_params->path_data, value);
                sprintf(zgrd_params->grid_point_path,
                        "%s/grid/grid_points/tweaked", zgrd_params->path_data);
            }
*/
            if (!strcasecmp(key, "physics_mode"))
                strcpy(zgrd_params->physics_mode, value);
        }
        fclose(fp);

        char *env_str;
        if ((env_str = getenv("GCRM_OUTPUT_PATH")) != NULL) {
            int len = strlen(env_str);
            strcpy(zgrd_params->path_output, env_str);
            if (zgrd_params->path_output[len-1] != '/') { /* add "/" at the end of path */
                zgrd_params->path_output[len] = '/';
                zgrd_params->path_output[len+1]   = '\0';
            }
        }
    }
    BCAST_INT(zgrd_params->initial_condition_select)
    BCAST_STR(zgrd_params->path_output,  256)
    // BCAST_STR(zgrd_params->grid_point_path, 256)
    BCAST_STR(zgrd_params->physics_mode, 128)
}

