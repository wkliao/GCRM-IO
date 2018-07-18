/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_params.h 4603 2017-12-07 07:12:52Z wkliao $
 */

#ifndef H_ZGRD_PARAMS
#define H_ZGRD_PARAMS

typedef struct {

    /* from MODULE ZGrd_params -----------------------------------------------*/
    int    l_ZGrd;
    int    l_density;
    int    l_tracer;
    int    l_vorticity;
    int    l_divergence;
    int    l_theta;

    char   equation_set_select[32];
    char   ZGrd_strng[128];

    int    l_restart;
    int    l_timing;

    /* the following 3 are moved to ZGrd_output.h
    int    l_use_binary_restart;
    int    l_read_restart_tmp;
    char   restart_input[512];
     */

    int    initial_condition_select;

    char   path_output[256];
    char   path_restart[256];
    char   path_data[256];
    char   file_restart[256];

    char   sfc_flux_method[256];
    int    l_diffusion_div;
    int    l_diffusion_eta;
    int    l_diffusion_tht;
    int    l_eta_enstrophy_damping;

    double k_diffusion_div;
    double k_diffusion_tht;

    char   grid_point_path[256];
    char   grid_point_select[256];
    char   io_mode[128];
    char   physics_mode[128];

    int    l_verbose;
    int    l_report;
    int    verbose_count;

} MODULE_ZGrd_params;

void ZGrd_load_inputs(const char *fname, MODULE_ZGrd_params *zgrd_param);

#endif
