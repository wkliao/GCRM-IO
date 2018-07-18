/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_clusdet.h 4605 2017-12-07 07:17:10Z wkliao $
 */

#ifndef H_ZGrd_clusdet
#define H_ZGrd_clusdet

#include <stdlib.h>

typedef struct {
    double /* DIMENSION(im, jm, nsdm) */
      ***clus_mask;         /* mask for cluster identification, faces */

    double /* DIMENSION(edgm, im, jm, nsdm) */
      ****clus_mask_edge;    /* mask for cluster identification, edges */

    double /* DIMENSION(2, crnm, im, jm, nsdm) */
      *****clus_mask_crn;     /* mask for cluster identification, corners */

    int nclust;          /* cluster count */

    int /* DIMENSION(im*jm*nsdm*npe_wrld) */
      *cluster;           /* cluster lookup table */

    int npairs;          /* merge pair count */

    int /* DIMENSION(jm*nsdm*npe_wrld) */
      *pairs;             /* merge pairs */

    int Top, Left, Botleft, Me; /* shortcuts for relative cell locations */

    int clus_num_cell, clus_num_edge, clus_num_crn;

    int min_cluster_size; /*  = 3  ignore clusters smaller than this */

    int clus_mark_all;    /*  = 0  set to 1 to output entire node if 
                                   _anything_ in it is interesting */

    int clusdet_enable; /* added from subroutine clusdet() */

} MODULE_ZGrd_clusdet;

void randomize_MODULE_ZGrd_clusdet(MODULE_ZGrd_clusdet *clusdet, int im, int jm,
                              int km, int nsdm, int edgm, int crnm,
                              int npe_wrld);
void init_MODULE_ZGrd_clusdet(MODULE_ZGrd_clusdet *clusdet, int im, int jm,
                              int km, int nsdm, int edgm, int crnm,
                              int npe_wrld);

void clusdet_initialize(const char *fname, MODULE_ZGrd_clusdet *clusdet);

void finalize_MODULE_ZGrd_clusdet(MODULE_ZGrd_clusdet *clusdet);

#endif
