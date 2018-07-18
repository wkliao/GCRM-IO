/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_params_tracer.c 1858 2013-06-12 04:49:45Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>

#include "ZGrd_params_tracer.h"

/*----< init_MODULE_ZGrd_params_tracer() >------------------------------------*/
void init_MODULE_ZGrd_params_tracer(MODULE_ZGrd_params_tracer *tracer,
                                    int                        enable_physics)
{
    tracer->l_wtr = 1;
    tracer->l_trc = 0;

    if (enable_physics)
        tracer->wtrm =  6;  /* number of water species */
    else
        tracer->wtrm =  1;  /* number of passive tracers */

    tracer->trcm =  1;

    /* total number of tracers */
    tracer->tracer_total = tracer->wtrm + tracer->trcm;

    /* the indices for the various water species (F2C z-based) */
    tracer->iqwv = 0;
    tracer->iqcw = 1;
    tracer->iqrw = 2;
    tracer->iqci = 3;
    tracer->iqsn = 4;
    tracer->iqgr = 5;
}

