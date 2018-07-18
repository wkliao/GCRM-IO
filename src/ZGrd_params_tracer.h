/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_params_tracer.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_ZGRD_PARAMS_TRACER
#define H_ZGRD_PARAMS_TRACER

typedef struct {

    /* from MODULE ZGrd_params_tracer ----------------------------------------*/
    int    l_wtr, l_trc;
    int    wtrm, /* number of water species */
           trcm; /* number of passive tracers */

    int    tracer_total;  /* total number of tracers */

    /*  the indices for the various water species */
    int    iqwv,  /* water vapor */
           iqcw,  /* cloud water */
           iqrw,  /* rain */
           iqci,  /* cloud ice */
           iqsn,  /* snow */
           iqgr;  /* graupel */

} MODULE_ZGrd_params_tracer;

void init_MODULE_ZGrd_params_tracer(MODULE_ZGrd_params_tracer *tracer,
                                    int enable_physics);

#endif
