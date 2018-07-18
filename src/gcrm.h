/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gcrm.h 4376 2017-08-03 21:52:01Z wkliao $
 */

#ifndef H_GCRM
#define H_GCRM

#include "grid_params.h"
#include "grid_subdomain.h"
#include "wrap_data.h"
#include "grid_metrics.h"
#include "util.h"
#include "ZGrd_output.h"
#include "ZGrd_params.h"
#include "ZGrd_params_vertical.h"
#include "ZGrd_params_tracer.h"
#include "ZGrd_params_time.h"
#include "ZGrd_vars_diagnostic.h"
#include "ZGrd_vars_prognostic.h"
#include "ZGrd_clusdet.h"
#include "phys_vars_diag_out.h"

int debug;

typedef struct {

    MODULE_ZGrd_output    zgrd_output;    /* defined in ZGrd_output.h */
    MODULE_ZGrd_params    zgrd_params;    /* defined in ZGrd_params.h */
    MODULE_grid_params    grid_params;    /* defined in grid_params.h */
    MODULE_grid_subdomain grid_subdomain; /* defined in grid_subdomain.h */
    MODULE_wrap_data      wrap_data;      /* defined in wrap_data.h */
    MODULE_grid_metrics   grid_metrics;   /* defined in grid_metrics.h */

    /* defined in MODULE parallel_params -------------------------------------*/
    int npe_wrld; /* number of process elements within the world communicator */

    int npe_io;   /* number of process elements within an io group */
    int rnk_wrld; /* an identifier for the local process within the
                     world communicator */

    int rnk_nonexistent;  /* = -1 an identifier for a non-existent process */

    char rnk_wrld_strng[8]; /* a string identifier for the local process within
                               the world communicator */

    MODULE_ZGrd_params_time time; /* defined in ZGrd_params_time.h */

    /* from subroutine clusdet() in ZGrd_clusdet.F90 */
    int    clusdet_enable;

    /* whether to use physics variables */
    int    enable_physics;

    MODULE_ZGrd_params_vertical vertical;  /* defined in ZGrd_params_vertical.h */
    MODULE_ZGrd_params_tracer   tracer;    /* defined in ZGrd_params_tracer.h */
    MODULE_ZGrd_vars_diagnostic vars_diag; /* defined in ZGrd_vars_diagnostic.h */
    MODULE_ZGrd_vars_prognostic vars_prog; /* defined in ZGrd_vars_prognostic.h */
    MODULE_ZGrd_clusdet         clusdet;   /* defined in ZGrd_clusdet.h */
    MODULE_phys_vars_diag_out   phys_vars; /* defined in phys_vars_diag_out.h */

} gcrm_parameters;



#include "gio.h"
int gcrm_param_init(char *in_fname, gcrm_parameters *gcrm_param,
                    gio_parameters *gio_param);
int gcrm_param_finalize(gcrm_parameters *gcrm_param,
                        gio_parameters *gio_param);
void initialize_gio_api(gio_parameters  *gio_param,
                        MODULE_ZGrd_params_tracer   *tracer,
                        MODULE_grid_params          *grid_params,
                        MODULE_grid_metrics         *grid_metrics,
                        MODULE_ZGrd_params_vertical *vertical,
                        MODULE_ZGrd_vars_diagnostic *vars_diag,
                        MODULE_ZGrd_vars_prognostic *vars_prog,
                        MODULE_ZGrd_clusdet         *clusdet,
                        MODULE_ZGrd_params_time     *time,
                        MODULE_phys_vars_diag_out   *phys_vars,
                        int                          enable_physics);

void pseudo_computation(gcrm_parameters *gcrm_param, int random_seed);

#endif
