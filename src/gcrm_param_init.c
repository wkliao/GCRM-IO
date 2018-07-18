/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gcrm_param_init.c 4602 2017-12-07 07:07:51Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "gcrm.h"
#include "gio.h"
#include "ZGrd_clusdet.h"
#include "ZGrd_output.h"
#include "ZGrd_params.h"
#include "ZGrd_params_tracer.h"
#include "ZGrd_params_vertical.h"
#include "ZGrd_params_time.h"
#include "ZGrd_vars_diagnostic.h"
#include "ZGrd_vars_prognostic.h"
#include "util.h"

/*----< parallel_load_inputs() >----------------------------------------------*/
/* from parallel/parallel_params.F90
 * Data updated:
 *     gcrm_param->npe_io
 */
static
void parallel_load_inputs(const char      *fname,
                          gcrm_parameters *gcrm_param)
{
    if (gcrm_param->rnk_wrld == 0) { /* only root reads from the file */
        char line[256], *value, *key;
        FILE *fp;

        OPEN_FILE(fname)

        while (fgets(line, 256, fp) != NULL) {
            GET_LINE_TOKEN_PAIR

            /* from parallel_load_inputs() called in ZGrd_main */
            if (!strcasecmp(key, "io_procs")) {
                gcrm_param->npe_io = atoi(value);
                break;
            }
        }
        fclose(fp);
    }
    BCAST_INT(gcrm_param->npe_io)
}

/*----< initialize_ZGrd() >---------------------------------------------------*/
/* initialize_ZGrd() is defined in ZGrd_initialize.F90 */
static
void initialize_ZGrd(const char      *in_fname,
                     gcrm_parameters *gcrm_param,
                     gio_parameters  *gio_param)
{
    MODULE_grid_params *grid_params = &gcrm_param->grid_params;

    gcrm_param->enable_physics=1;

    if (!strncasecmp(gcrm_param->zgrd_params.physics_mode, "disabled", 7))
        gcrm_param->enable_physics = 0;  /* not disabled */

    /* set subdomains. partition global grid and assign pieces to processes
       in grid/grid_subdomain.F90 */
    initialize_subdomain(grid_params, &gcrm_param->grid_subdomain);

    /* set tree data structure. each process builds its portion of the global
       grid.  grid/grid_connectivity.F90 */
    initialize_grid_connectivity(grid_params, "world");

    /* set parallel communication for ghost cell updates of 2D-array data
       structure */
    initialize_wrap(&gcrm_param->wrap_data,
                    grid_params,
                    "world","ZGrd",
                    grid_params->level_max,
                    grid_params->sbdmn[grid_params->level_max].sbdmn_iota,
                    grid_params->sbdmn[grid_params->level_max].lst,
                    0);

    /* set grid metrics */
    initialize_grid_metrics(grid_params,
                            &gcrm_param->grid_subdomain,
                            &gcrm_param->grid_metrics,
                            &gcrm_param->wrap_data,
                            "ZGrd",
                            gcrm_param->zgrd_params.grid_point_select,
                            NULL);

    /* call initialize_params_vertical() in ZGrd_params_vertical.F90 */
    initialize_params_vertical(&gcrm_param->vertical);

    /* call initialize_params_time() in ZGrd_params_time.F90 */
    initialize_params_time(in_fname, &gcrm_param->time,
                           gcrm_param->grid_params.level_max);

    init_MODULE_ZGrd_params_tracer(&gcrm_param->tracer,
                                   gcrm_param->enable_physics);

    /* this must be called after initialize_params_time(), as time.nprog
     * and time.ntend are defined there */
    init_MODULE_ZGrd_vars_prognostic(&gcrm_param->vars_prog,
                                     gcrm_param->grid_params.im,
                                     gcrm_param->grid_params.jm,
                                     gcrm_param->vertical.km,
                                     gcrm_param->time.nprog,
                                     gcrm_param->grid_params.nsdm,
                                     gcrm_param->time.ntend,
                                     gcrm_param->tracer.wtrm,
                                     gcrm_param->tracer.trcm);

    initialize_vars_prognostic(&gcrm_param->vars_prog,
                                     gcrm_param->grid_params.im,
                                     gcrm_param->grid_params.jm,
                                     gcrm_param->vertical.km,
                                     gcrm_param->time.nprog,
                                     gcrm_param->grid_params.nsdm,
                                     gcrm_param->time.ntend,
                                     gcrm_param->tracer.wtrm,
                                     gcrm_param->tracer.trcm);

    init_MODULE_ZGrd_vars_diagnostic(&gcrm_param->vars_diag,
                                     gcrm_param->grid_params.im,
                                     gcrm_param->grid_params.jm,
                                     gcrm_param->vertical.km,
                                     gcrm_param->grid_params.nsdm,
                                     gcrm_param->tracer.wtrm,
                                     gcrm_param->tracer.trcm,
                                     gcrm_param->grid_params.edgm,
                                     gcrm_param->grid_params.crnm);

    initialize_vars_diagnostic(&gcrm_param->vars_diag,
                                     gcrm_param->grid_params.im,
                                     gcrm_param->grid_params.jm,
                                     gcrm_param->vertical.km,
                                     gcrm_param->grid_params.nsdm,
                                     gcrm_param->tracer.wtrm,
                                     gcrm_param->tracer.trcm,
                                     gcrm_param->grid_params.edgm,
                                     gcrm_param->grid_params.crnm);

    if (gcrm_param->enable_physics)
        init_MODULE_phys_vars_diag_out(&gcrm_param->phys_vars,
                                     gcrm_param->grid_params.nsdm,
                                     gcrm_param->vertical.km,
                                     gcrm_param->grid_params.jm,
                                     gcrm_param->grid_params.im);

    /* pass GCRM parameters to GIO */
    gio_init(gcrm_param, gio_param);

    /* initialize_gio_api() is called in ZGrd_register.F90
     * It sets contents of all grid variables and register all grid and
     * non-grid variables to the gio library */
    initialize_gio_api(gio_param,
                       &gcrm_param->tracer,
                       &gcrm_param->grid_params,
                       &gcrm_param->grid_metrics,
                       &gcrm_param->vertical,
                       &gcrm_param->vars_diag,
                       &gcrm_param->vars_prog,
                       &gcrm_param->clusdet,
                       &gcrm_param->time,
                       &gcrm_param->phys_vars,
                       gcrm_param->enable_physics);
}

/*----< pseudo_computation() >------------------------------------------------*/
void pseudo_computation(gcrm_parameters *gcrm_param, int random_seed)
{
    /* set random number seed */
    srandom(gcrm_param->rnk_wrld + random_seed);

    /* just some dummy computation of random number generation */
    randomize_MODULE_ZGrd_clusdet(&gcrm_param->clusdet,
                                  gcrm_param->grid_params.im,
                                  gcrm_param->grid_params.jm,
                                  gcrm_param->vertical.km,
                                  gcrm_param->grid_params.nsdm,
                                  gcrm_param->grid_params.edgm,
                                  gcrm_param->grid_params.crnm,
                                  gcrm_param->npe_wrld);

    randomize_ZGrd_vars_prognostic(&gcrm_param->vars_prog,
                                   gcrm_param->grid_params.im,
                                   gcrm_param->grid_params.jm,
                                   gcrm_param->vertical.km,
                                   gcrm_param->time.nprog,
                                   gcrm_param->grid_params.nsdm,
                                   gcrm_param->time.ntend,
                                   gcrm_param->tracer.wtrm,
                                   gcrm_param->tracer.trcm);

    randomize_ZGrd_vars_diagnostic(&gcrm_param->vars_diag,
                                   gcrm_param->grid_params.im,
                                   gcrm_param->grid_params.jm,
                                   gcrm_param->vertical.km,
                                   gcrm_param->grid_params.nsdm,
                                   gcrm_param->tracer.wtrm,
                                   gcrm_param->tracer.trcm,
                                   gcrm_param->grid_params.edgm,
                                   gcrm_param->grid_params.crnm);

    if (gcrm_param->enable_physics)
        randomize_phys_vars_diag_out(&gcrm_param->phys_vars,
                                     gcrm_param->grid_params.nsdm,
                                     gcrm_param->vertical.km,
                                     gcrm_param->grid_params.jm,
                                     gcrm_param->grid_params.im);
}

/*----< gcrm_param_init() >---------------------------------------------------*/
/* this subroutine is based on ZGrd_main.F90  */
int gcrm_param_init(char                *in_fname,
                    gcrm_parameters     *gcrm_param,
                    gio_parameters      *gio_param)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &gcrm_param->rnk_wrld);
    MPI_Comm_size(MPI_COMM_WORLD, &gcrm_param->npe_wrld);

    /* set random number seed */
    srandom(gcrm_param->rnk_wrld);

    /* All functions called below with prefix name init_MODULE_ correspond to
     * the F90 modules that define array variables with static sizes and
     * dimensions.
     * In this C version, definition for all static array variables have been
     * changed to dynamic memory allocation, and initialized to zeros.
     * Their malloc must be done before their initializations.
     */

    /* File: grid/grid_params.F90, Fortran MODULE grid_params
     * Data updated:
     *     grid_params->level_max
     *     gcrm_param->grid_params->sbdmn_iota
     *     gcrm_param->grid_params->level_glbl
     *     gcrm_param->grid_params->cell_max
     *     gcrm_param->grid_params->nsdm_glbl
     *     gcrm_param->grid_params->nsdm
     *     gcrm_param->grid_params->im
     *     gcrm_param->grid_params->jm
     *     gcrm_param->grid_params->edgm
     *     gcrm_param->grid_params->crnm
     *     gcrm_param->grid_params->tag_nonexistent
     *     gcrm_param->grid_params->grid_node_total
     *     gcrm_param->grid_params->grid_node_size 
     *     gcrm_param->grid_params->grid_node_memory_max
     *     gcrm_param->grid_params->alfalfa
     *     gcrm_param->grid_params->pentagon_latitude
     *     gcrm_param->grid_params->sbdmn
     *     gcrm_param->grid_params->path_next
     *     gcrm_param->grid_params->path_real
     *     gcrm_param->grid_params->path_ghst
     *     gcrm_param->grid_params->nm_lvl     
     *     gcrm_param->grid_params->nm_real_lvl
     *     gcrm_param->grid_params->nm_ghst_lvl
     *     gcrm_param->grid_params->grid[12].dn[4]
     */
    init_MODULE_grid_params(&gcrm_param->grid_params, in_fname);

    /* File: grid/grid_subdomain.F90 Fortran MODULE grid_subdomain
     * Data updated:
     *     gcrm_param->grid_subdomain->level_threshold
     *     gcrm_param->grid_subdomain->cell_min
     *     gcrm_param->grid_subdomain->distribution_pattern
     *     gcrm_param->grid_subdomain->l_sbdmn_pntgn_north
     *     gcrm_param->grid_subdomain->l_sbdmn_pntgn_south
     *     gcrm_param->grid_subdomain->l_sbdmn_north_pole
     *     gcrm_param->grid_subdomain->l_sbdmn_south_pole
     */
    init_MODULE_grid_subdomain(&gcrm_param->grid_subdomain,
                               gcrm_param->grid_params.nsdm);

    /* File: grid/grid_metrics.F90 Fortran MODULE grid_metrics
     * Data malloc-ed
     *     gcrm_param->grid_metrics->l_msk_fac
     *     gcrm_param->grid_metrics->l_msk_crn
     *     gcrm_param->grid_metrics->l_msk_edg
     *     gcrm_param->grid_metrics->tag_glbl
     *     gcrm_param->grid_metrics->tag_glbl_nghbr
     *     gcrm_param->grid_metrics->point
     *     gcrm_param->grid_metrics->corner
     *     gcrm_param->grid_metrics->area
     *     gcrm_param->grid_metrics->area_inv
     *     gcrm_param->grid_metrics->d_point
     *     gcrm_param->grid_metrics->d_point_inv
     *     gcrm_param->grid_metrics->d_edge
     *     gcrm_param->grid_metrics->d_edge_inv
     *     gcrm_param->grid_metrics->point_crn
     *     gcrm_param->grid_metrics->area_crn
     *     gcrm_param->grid_metrics->area_inv_crn
     *     gcrm_param->grid_metrics->area_kite_crn
     *     gcrm_param->grid_metrics->d_point_crn
     *     gcrm_param->grid_metrics->d_point_inv_crn
     *     gcrm_param->grid_metrics->d_edge_crn
     *     gcrm_param->grid_metrics->d_edge_inv_crn
     *     gcrm_param->grid_metrics->rlx_wght_crn
     *     gcrm_param->grid_metrics->point_edg
     *     gcrm_param->grid_metrics->nrm_edg
     *     gcrm_param->grid_metrics->tng_edg
     *     gcrm_param->grid_metrics->area_edg
     *     gcrm_param->grid_metrics->wghts_crn
     *     gcrm_param->grid_metrics->vctr_wghts_crn
     *     gcrm_param->grid_metrics->vctr_wghts_crn2
     *     gcrm_param->grid_metrics->wghts_wind_crn
     *     gcrm_param->grid_metrics->laplacian_wghts
     *     gcrm_param->grid_metrics->laplacian_wghts_crn
     * Data updated:
     *     gcrm_param->grid_metrics->grid_point_mode
     *     gcrm_param->grid_metrics->grid_point_file
    */
    init_MODULE_grid_metrics(&gcrm_param->grid_metrics,
                             gcrm_param->grid_params.im,
                             gcrm_param->grid_params.jm,
                             gcrm_param->grid_params.nsdm,
                             gcrm_param->grid_params.edgm,
                             gcrm_param->grid_params.crnm);

    /* File: ZGrd/ZGrd_output.F90, Fortran MODULE ZGrd_output
     * Data updates:
     *     gcrm_param->zgrd_output->descfile
     *     gcrm_param->zgrd_output->outfileconfig
     *     gcrm_param->zgrd_output->iotype
     *     gcrm_param->zgrd_output->cdf_output_path
     *     gcrm_param->zgrd_output->cdf_grid_option
     *     gcrm_param->zgrd_output->cdf_sep_grid
     *     gcrm_param->zgrd_output->cdf_output_freq
     *     gcrm_param->zgrd_output->cdf_output_nsamples
     *     gcrm_param->zgrd_output->cdf_restart_path
     *     gcrm_param->zgrd_output->l_restart_overwrite
     *     gcrm_param->zgrd_output->restart_output[0]
     *     gcrm_param->zgrd_output->l_read_restart_tmp
     *     gcrm_param->zgrd_output->l_use_binary_restart
     *     gcrm_param->zgrd_output->restart_interval
     */
    init_MODULE_ZGrd_output(&gcrm_param->zgrd_output);

    /* array memory allocation in ZGrd_clusdet.F90
     * Data malloc-ed:
     *     gcrm_param->clusdet->clus_mask
     *     gcrm_param->clusdet->clus_mask_edge
     *     gcrm_param->clusdet->clus_mask_crn
     *     gcrm_param->clusdet->cluster
     *     gcrm_param->clusdet->pairs
     * Data updated:
     *     gcrm_param->clusdet->min_cluster_size
     *     gcrm_param->clusdet->clus_mark_all
     */
    init_MODULE_ZGrd_clusdet(&gcrm_param->clusdet,
                             gcrm_param->grid_params.im,
                             gcrm_param->grid_params.jm,
                             gcrm_param->vertical.km,
                             gcrm_param->grid_params.nsdm,
                             gcrm_param->grid_params.edgm,
                             gcrm_param->grid_params.crnm,
                             gcrm_param->npe_wrld);

    init_MODULE_wrap_data(&gcrm_param->wrap_data);

    /* Now below are the functions called in main() in ZGrd_main.F90 */
    //     call gin_load_inputs() in utilities/ginput.F90

    /* parallel_load_inputs() is in parallel/parallel_params.F90
     * It only reads io_procs from the input parameter file
     * Data updated:
     *     gcrm_param->npe_io
     */
    parallel_load_inputs(in_fname, gcrm_param);

    /* io_load_inputs() is defined in ZGrd_output.F90
     * It reads input parameter file (e.g. zgrd.in)
     * Data updated:
     *     gcrm_param->zgrd_output->iotype
     *     gcrm_param->zgrd_output->descfile
     *     gcrm_param->zgrd_output->outfileconfig
     *     gcrm_param->zgrd_output->cdf_output_path
     *     gcrm_param->zgrd_output->cdf_output_freq
     *     gcrm_param->zgrd_output->cdf_output_nsamples
     *     gcrm_param->zgrd_output->cdf_grid_option
     *     gcrm_param->zgrd_output->cdf_sep_grid
     *     gcrm_param->zgrd_output->l_use_binary_restart
     *     gcrm_param->zgrd_output->l_read_restart_tmp
     *     gcrm_param->zgrd_output->l_restart_overwrite
     *     gcrm_param->zgrd_output->restart_interval
     *     gcrm_param->zgrd_output->restart_input
     *     gcrm_param->zgrd_output->restart_output
     */
    io_load_inputs(in_fname, &gcrm_param->zgrd_output);

    /* ZGrd_load_inputs() is defined in ZGrd_params.F90
     * Data updates:
     *     gcrm_param->zgrd_params->l_ZGrd
     *     gcrm_param->zgrd_params->l_density
     *     gcrm_param->zgrd_params->l_tracer
     *     gcrm_param->zgrd_params->l_vorticity
     *     gcrm_param->zgrd_params->l_divergence
     *     gcrm_param->zgrd_params->l_theta
     *     gcrm_param->zgrd_params->l_restart
     *     gcrm_param->zgrd_params->l_timing
     *     gcrm_param->zgrd_params->l_diffusion_div
     *     gcrm_param->zgrd_params->l_diffusion_eta
     *     gcrm_param->zgrd_params->l_diffusion_tht
     *     gcrm_param->zgrd_params->l_eta_enstrophy_damping
     *     gcrm_param->zgrd_params->l_verbose
     *     gcrm_param->zgrd_params->l_report
     *     gcrm_param->zgrd_params->equation_set_select
     *     gcrm_param->zgrd_params->ZGrd_strng
     *     gcrm_param->zgrd_params->initial_condition_select
     *     gcrm_param->zgrd_params->path_output
     *     gcrm_param->zgrd_params->path_restart
     *     // gcrm_param->zgrd_params->path_data
     *     gcrm_param->zgrd_params->file_restart
     *     gcrm_param->zgrd_params->sfc_flux_method
     *     gcrm_param->zgrd_params->k_diffusion_div
     *     gcrm_param->zgrd_params->k_diffusion_tht
     *     gcrm_param->zgrd_params->grid_point_path[0]
     *     gcrm_param->zgrd_params->grid_point_select
     *     gcrm_param->zgrd_params->io_mode
     *     gcrm_param->zgrd_params->physics_mode
     *     gcrm_param->zgrd_params->verbose_count
     *     gcrm_param->zgrd_params->initial_condition_select
     */
    ZGrd_load_inputs(in_fname, &gcrm_param->zgrd_params);

    /* clusdet_initialize() is defined in ZGrd_clusdet.F90
     * It only reads clusdet_enable from input parameter file.
     * Data updated:
     *     gcrm_param->clusdet->clusdet_enable
     */
    clusdet_initialize(in_fname, &gcrm_param->clusdet);

    /* initialize_ZGrd() is defined in ZGrd_initialize.F90 */
    initialize_ZGrd(in_fname, gcrm_param, gio_param);

    return 1;
}

/*----< gcrm_param_finalize() >-----------------------------------------------*/
int gcrm_param_finalize(gcrm_parameters *gcrm_param,
                        gio_parameters  *gio_param)
{
    if (gcrm_param->enable_physics)
        finalize_MODULE_phys_vars_diag_out(&gcrm_param->phys_vars);

    finalize_MODULE_ZGrd_vars_diagnostic(&gcrm_param->vars_diag);
    finalize_MODULE_ZGrd_vars_prognostic(&gcrm_param->vars_prog);
    finalize_MODULE_wrap_data(&gcrm_param->wrap_data);
    finalize_MODULE_ZGrd_clusdet(&gcrm_param->clusdet);
    finalize_MODULE_grid_metrics(&gcrm_param->grid_metrics);
    finalize_MODULE_grid_subdomain(&gcrm_param->grid_params,
                                   &gcrm_param->grid_subdomain);
    finalize_MODULE_grid_params(&gcrm_param->grid_params);

    return 1;
}
