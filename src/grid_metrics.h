/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: grid_metrics.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_GRID_METRICS
#define H_GRID_METRICS

#if HAVE_CONFIG_H
#include "config.h"
#endif

typedef struct {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  logical masks
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    int  ***l_msk_fac;            /* [nsdm][jm][im] logical mask for faces
                                                    (cell centers) */
    int ****l_msk_crn;            /* [nsdm][jm][im][crnm] logical mask for
                                                    corners */
    int ****l_msk_edg;            /*[nsdm][jm][im][edgm] logical mask for
                                                    edges */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  metrics defined at cell centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    int      ***tag_glbl;         /* [nsdm][jm][im] */
    int     ****tag_glbl_nghbr;   /* [nsdm][jm][im][6] */
    double  ****point;            /* [nsdm][jm][im][3] */
    double *****corner;           /* [nsdm][jm][im][6][3] */
    double   ***area;             /* [nsdm][jm][im] */
    double   ***area_inv;         /* [nsdm][jm][im] */
    double  ****d_point;          /* [nsdm][jm][im][6] */
    double  ****d_point_inv;      /* [nsdm][jm][im][6] */
    double  ****d_edge;           /* [nsdm][jm][im][6] */
    double  ****d_edge_inv;       /* [nsdm][jm][im][6] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  metrics defined at cell corners (abbreviated crn)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    double *****point_crn;        /* [nsdm][jm][im][crnm][3] */
    double  ****area_crn;         /* [nsdm][jm][im][crnm] */
    double  ****area_inv_crn;     /* [nsdm][jm][im][crnm] */
    double *****area_kite_crn;    /* [nsdm][jm][im][crnm][3] */
    double *****d_point_crn;      /* [nsdm][jm][im][crnm][3] */
    double *****d_point_inv_crn;  /* [nsdm][jm][im][crnm][3] */
    double *****d_edge_crn;       /* [nsdm][jm][im][crnm][3] */
    double *****d_edge_inv_crn;   /* [nsdm][jm][im][crnm][3] */
    double *****rlx_wght_crn;     /* [nsdm][jm][im][crnm][0:3] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  metrics defined at cell edges (abbreviated edg)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    double *****point_edg;        /* [nsdm][jm][im][edgm][3] */
    double *****nrm_edg;          /* [nsdm][jm][im][edgm][3] */
    double *****tng_edg;          /* [nsdm][jm][im][edgm][3] */
    double  ****area_edg;         /* [nsdm][jm][im][edgm] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  weights
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    double  *****wghts_crn;       /* [nsdm][jm][im][crnm][3] weights to average cell centers to corners */
    double ******vctr_wghts_crn;  /* [nsdm][jm][im][crnm][2][3] */
    double ******vctr_wghts_crn2; /* [nsdm][jm][im][crnm][2][3] */
/*-----------------------------------------------------------------------
!  for wind at the corners in the GCRM
!-----------------------------------------------------------------------*/
    double ******wghts_wind_crn;      /* [nsdm][jm][im][crnm][3][3] */
    double   ****laplacian_wghts;     /* [nsdm][jm][im][6] */
    double  *****laplacian_wghts_crn; /* [nsdm][jm][im][crnm][3] */

    int grid_point_mode; /* = 2; */
    int grid_point_file; /* = 19; */

    int file_ncid;
    int var_grd_point_xyz_ncid;
} MODULE_grid_metrics;

#include "grid_params.h"
#include "grid_subdomain.h"
#include "grid_metrics.h"
#include "wrap_data.h"

void initialize_grid_metrics(MODULE_grid_params    *grid_params,
                             MODULE_grid_subdomain *grid_subdomain,
                             MODULE_grid_metrics   *grid_metrics,
                             MODULE_wrap_data      *wrap_data,
                             char                  *wrap_name,
                             char                  *grid_point_select,
                             char                  *grid_point_path_tmpry);

void init_MODULE_grid_metrics(MODULE_grid_metrics *grid_metrics,
                              int im,   int jm, int nsdm, int edgm, int crnm);

void finalize_MODULE_grid_metrics(MODULE_grid_metrics *grid_metrics);

#endif
