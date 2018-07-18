/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_clusdet.c 4603 2017-12-07 07:12:52Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
extern int errno;

#include <mpi.h>

#include "ZGrd_clusdet.h"
#include "util.h"


/*----< randomize_MODULE_ZGrd_clusdet() >-------------------------------------*/
void randomize_MODULE_ZGrd_clusdet(MODULE_ZGrd_clusdet *clusdet,
                                   int                  im,
                                   int                  jm,
                                   int                  km,
                                   int                  nsdm,
                                   int                  edgm,
                                   int                  crnm,
                                   int                  npe_wrld)
{
    /* DIMENSION(im, jm, nsdm) */
    random_3D_dbl(clusdet->clus_mask,  nsdm, jm, im);

    /* DIMENSION(edgm, im, jm, nsdm) */
    random_4D_dbl(clusdet->clus_mask_edge,  nsdm, jm, im, edgm);

    /* DIMENSION(2, crnm, im, jm, nsdm) */
    random_5D_dbl(clusdet->clus_mask_crn,  nsdm, jm, im, crnm, 2);

    /* DIMENSION(im*jm*nsdm*npe_wrld) */
    random_1D_int(clusdet->cluster,  im*jm*nsdm*npe_wrld);

    /* DIMENSION(jm*nsdm*npe_wrld) */
    random_1D_int(clusdet->pairs,  jm*nsdm*npe_wrld);
}

/*----< init_MODULE_ZGrd_clusdet() >------------------------------------------*/
/* from ZGrd_clusdet.F90
 * Data malloc-ed:
 *     clusdet->clus_mask
 *     clusdet->clus_mask_edge
 *     clusdet->clus_mask_crn
 *     clusdet->cluster
 *     clusdet->pairs
 * Data updated:
 *     clusdet->min_cluster_size
 *     clusdet->clus_mark_all
 */
void init_MODULE_ZGrd_clusdet(MODULE_ZGrd_clusdet *clusdet,
                              int                  im,
                              int                  jm,
                              int                  km,
                              int                  nsdm,
                              int                  edgm,
                              int                  crnm,
                              int                  npe_wrld)
{
    /* DIMENSION(im, jm, nsdm) */
    clusdet->clus_mask = calloc_3D_dbl(nsdm, jm, im);

    /* DIMENSION(edgm, im, jm, nsdm) */
    clusdet->clus_mask_edge = calloc_4D_dbl(nsdm, jm, im, edgm);

    /* DIMENSION(2, crnm, im, jm, nsdm) */
    clusdet->clus_mask_crn = calloc_5D_dbl(nsdm, jm, im, crnm, 2);

    /* DIMENSION(im*jm*nsdm*npe_wrld) */
    clusdet->cluster = calloc_1D_int(im*jm*nsdm*npe_wrld);

    /* DIMENSION(jm*nsdm*npe_wrld) */
    clusdet->pairs = calloc_1D_int(jm*nsdm*npe_wrld);

    clusdet->min_cluster_size = 3;

    clusdet->clus_mark_all = 0;
}

/*----< clusdet_initialize() >------------------------------------------------*/
/* from ZGrd_clusdet.F90
 * Data updated:
 *     clusdet->clusdet_enable
 */
void clusdet_initialize(const char          *fname,
                        MODULE_ZGrd_clusdet *clusdet)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) { /* only root process reads from the file */
        char line[256], *value, *key;
        FILE *fp;

        OPEN_FILE(fname)

        while (fgets(line, 256, fp) != NULL) {
            GET_LINE_TOKEN_PAIR

            /* This actually is called from clusdet() in ZGrd_main */
            if (!strcasecmp(key, "clusdet_enable")) {
                clusdet->clusdet_enable = 0;
                if (!strcasecmp(value, "true"))
                    clusdet->clusdet_enable = 1;
                break;
            }
        }
        fclose(fp);
    }
    BCAST_INT(clusdet->clusdet_enable)
}

/*----< finalize_MODULE_ZGrd_clusdet() >--------------------------------------*/
void finalize_MODULE_ZGrd_clusdet(MODULE_ZGrd_clusdet *clusdet)
{
    free_3D_dbl(clusdet->clus_mask);
    free_4D_dbl(clusdet->clus_mask_edge);
    free_5D_dbl(clusdet->clus_mask_crn);
    free_1D_int(clusdet->cluster);
    free_1D_int(clusdet->pairs);
}


