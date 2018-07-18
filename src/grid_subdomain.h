/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: grid_subdomain.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_GRID_SUBDOMAIN
#define H_GRID_SUBDOMAIN

#if HAVE_CONFIG_H
#include "config.h"
#endif

typedef struct {
    int  level_threshold;      /* = 2,  resolution at which one processor owns 
                                        the whole grid */
    int  cell_min;             /* = 16, minimum number of cells per subdomain 
                                        must be one of 4,16,64,256,1024,... */
    int  distribution_pattern; /* = 2   pattern to assign subdomains to procs
                                        at the finest resolution */

    int *l_sbdmn_pntgn_north;  /* [grid_params.nsdm]  .TRUE. iff subdomain contains 
                                                       mid-latitiude pentagon */
    int *l_sbdmn_pntgn_south;  /* [grid_params.nsdm]  .TRUE. iff subdomain contains 
                                                       mid-latitiude pentagon */

    int *l_sbdmn_north_pole;   /* [grid_params.nsdm]  .TRUE. iff subdomain contains 
                                                       pole pentagon */
    int *l_sbdmn_south_pole;   /* [grid_params.nsdm]  .TRUE. iff subdomain contains 
                                                       pole pentagon */
} MODULE_grid_subdomain;

#include "grid_params.h"

void init_MODULE_grid_subdomain(MODULE_grid_subdomain *grid_subdomain, int nsdm);
void finalize_MODULE_grid_subdomain(MODULE_grid_params *grid_params,
                                    MODULE_grid_subdomain *grid_subdomain);
void initialize_subdomain(MODULE_grid_params *grid_params,
                          MODULE_grid_subdomain *grid_subdomain);
int get_big_block_index(int lvl, int iota, int nsd_glbl, int *ix);

#endif
