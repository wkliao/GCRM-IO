/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_output.h 4603 2017-12-07 07:12:52Z wkliao $
 */

#ifndef H_ZGRD_OUTPUT
#define H_ZGRD_OUTPUT

#if HAVE_CONFIG_H
#include "config.h"
#endif

/*
   MODULE ZGrd_output
*/


typedef struct {
    char   descfile[256];
    char   outfileconfig[256];
    char   iotype[128];
    char   cdf_output_path[128];
    char   cdf_grid_option[128];
    /* int    cdf_grid_option; 2: all_files, 1:nogrid, 3:sep_grid */
    char   cdf_sep_grid[128];
    float  cdf_output_freq;
    int    cdf_output_nsamples;
    char   cdf_restart_path[256];

    int    l_restart_overwrite;
    double restart_interval;
    char   restart_output[512];

    /* the following 3 are from ZGrd_params.h */
    int    l_use_binary_restart;
    int    l_read_restart_tmp;
    char   restart_input[512];

} MODULE_ZGrd_output;


void init_MODULE_ZGrd_output(MODULE_ZGrd_output *zgrd_output);
void io_load_inputs(const char *fname, MODULE_ZGrd_output *zgrd_output);

#endif
