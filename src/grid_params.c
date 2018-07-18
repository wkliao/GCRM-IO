/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: grid_params.c 4603 2017-12-07 07:12:52Z wkliao $
 */

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* memset() */
#include <mpi.h>

#include "grid_params.h"
#include "gcrm.h"


/*----< infile_read_int() >---------------------------------------------------*/
static
void infile_read_int(const char *fname,
                     const char *in_key,
                     int        *in_value)
{
    char line[256], *value, *key;
    FILE *fp;

    OPEN_FILE(fname)

    while (fgets(line, 256, fp) != NULL) {
        GET_LINE_TOKEN_PAIR

        if (!strcasecmp(key, in_key)) {
            *in_value = atoi(value);
            break;
        }
    }
    fclose(fp);
}

/*----< init_MODULE_grid_params() >------------------------------------------*/
/* data updated:
 *     grid_params->level_max
 *     grid_params->sbdmn_iota
 *     grid_params->level_glbl
 *     grid_params->cell_max
 *     grid_params->nsdm_glbl
 *     grid_params->nsdm
 *     grid_params->im
 *     grid_params->jm
 *     grid_params->edgm
 *     grid_params->crnm
 *     grid_params->tag_nonexistent
 *     grid_params->grid_node_total
 *     grid_params->grid_node_size 
 *     grid_params->grid_node_memory_max
 *     grid_params->alfalfa
 *     grid_params->pentagon_latitude
 *     grid_params->sbdmn
 *     grid_params->path_next
 *     grid_params->path_real
 *     grid_params->path_ghst
 *     grid_params->nm_lvl     
 *     grid_params->nm_real_lvl
 *     grid_params->nm_ghst_lvl
 *     grid_params->grid[12].dn[4]
 */
void init_MODULE_grid_params(MODULE_grid_params *grid_params,
                             const char         *in_fname)
{
    int rank, npe_wrld;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npe_wrld);

    if (rank == 0) /* only root process reads from the file */
        infile_read_int(in_fname, "level_max", &grid_params->level_max);
    MPI_Bcast(&grid_params->level_max, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* In the Fortran GCRM codes, SBDMN_IOTA is set at the configure time,
       as in grid/grid_params.F90:85:      sbdmn_iota  =  SBDMN_IOTA,  &!
       Below is extracted from file "configure".

       SBDMN_IOTA=0
       found=no
       for nblocks in 10 40 160 640 2560 10240 40960 163840 655360 2621440
       do
           val=$((nblocks%NPE_WRLD))
           if test $val -eq 0; then :
         found=yes; break
       else
         SBDMN_IOTA=$((SBDMN_IOTA+1))
       fi
       done
       if test "x$SBDMN_IOTA" = x0; then :
         if test "x$NPE_WRLD" != x1; then :
         SBDMN_IOTA=1
       fi
    */
    int i, j;
    int nblocks[10] = {10, 40, 160, 640, 2560, 10240, 40960, 163840, 655360, 2621440};

    for (i=0; i<10; i++)
        if (nblocks[i] % npe_wrld == 0)
            break;
    if (i == 10) {
        printf("nsdm_glbl must be evenly divisble by NPE_WRLD (currently=%d)\n",npe_wrld);
        printf("nsdm_glbl = {10, 40, 160, 640, 2560, 10240, 40960, 163840, 655360, 2621440}\n");
        printf("nsdm_glbl = 10 * 2**(2 * sbdmn_iota)\n");
        printf("GCRM picks the first number from the above list that is divisible by NPE_WRLD\n");
        printf("Please adjust NPE_WRLD accordingly\n");
        ABORT
    }

    /* sbdmn_iota is the first index of nblocks[] that is divisible by
       npe_wrld */
    grid_params->sbdmn_iota = i;
    // if (i == 0) grid_params->sbdmn_iota = 1;

    grid_params->level_glbl = grid_params->sbdmn_iota + 1;

    int level_max = grid_params->level_max;
    int cell_max2 = POWER2(level_max);
    cell_max2 = cell_max2 * cell_max2;
    grid_params->cell_max = 2 + 10 * cell_max2;
    /* +2 are cells of north and south poles */
    grid_params->nsdm_glbl = 10 * POWER2(2 * grid_params->sbdmn_iota);
    grid_params->nsdm = grid_params->nsdm_glbl / npe_wrld;
    grid_params->im = 2 + POWER2(level_max - grid_params->sbdmn_iota);
    grid_params->jm = grid_params->im;
    /* +2 are ghost cells along i and j directions */

    grid_params->edgm=3;
    grid_params->crnm=2;

    grid_params->tag_nonexistent = -999;

    /* the total number of grid nodes allocated by the local process */
    grid_params->grid_node_total =   0;
    grid_params->grid_node_size  = 300; /* size in bytes of one grid node */
    grid_params->grid_node_memory_max = 500000000;

    grid_params->alfalfa = /* 2 Pi/5 */
         1.2566370614359172953850573533118011536788677597500423;
    grid_params->pentagon_latitude = /* latitude of icosahedron vertex */
         0.4636476090008061162142562314612144020285370542861202;

    // sbdmn_node sbdmn[level_max+1];
    grid_params->sbdmn = (sbdmn_node*) tcalloc((level_max+1), sizeof(sbdmn_node));

    grid_params->path_next = (grid_node**) tmalloc((level_max+1) * sizeof(grid_node*));
    grid_params->path_real = (grid_node**) tmalloc((level_max+1) * sizeof(grid_node*));
    grid_params->path_ghst = (grid_node**) tmalloc((level_max+1) * sizeof(grid_node*));

    grid_params->nm_lvl      = (int*) tmalloc((level_max+1) * sizeof(int));
    grid_params->nm_real_lvl = (int*) tmalloc((level_max+1) * sizeof(int));
    grid_params->nm_ghst_lvl = (int*) tmalloc((level_max+1) * sizeof(int));

    for (i=0; i<12; i++)
        for (j=0; j<4; j++)
            grid_params->grid[i].dn[j] = NULL;
}

/*----< finalize_MODULE_grid_params() >---------------------------------------*/
void finalize_MODULE_grid_params(MODULE_grid_params *grid_params)
{
    if (grid_params->sbdmn     != NULL)   tfree(grid_params->sbdmn);
    if (grid_params->path_next != NULL)   tfree(grid_params->path_next);
    if (grid_params->path_real != NULL)   tfree(grid_params->path_real);
    if (grid_params->path_ghst != NULL)   tfree(grid_params->path_ghst);

    if (grid_params->nm_lvl      != NULL) tfree(grid_params->nm_lvl);
    if (grid_params->nm_real_lvl != NULL) tfree(grid_params->nm_real_lvl);
    if (grid_params->nm_ghst_lvl != NULL) tfree(grid_params->nm_ghst_lvl);

    finalize_grid_connectivity(grid_params);
}

/*----< set_ptr() >----------------------------------------------------------*/
grid_node* set_ptr(grid_node *grid,
                   int        level_max,
                   int        level,
                   int        tag_glbl)  /* 1-based */
{
    int i, *twopowtwo, m, *path, lvl;
    grid_node *ptr;

    /* Fortran twopowtwo(0:level_max) and path(0:level_max) */
    twopowtwo = (int*) malloc((level_max+1) * sizeof(int));
    path      = (int*) calloc((level_max+1),  sizeof(int));

    for (i=0; i<=level_max; i++)
        twopowtwo[i] = POWER2(2*(level_max-i));

    switch (tag_glbl) {
        case 1:
            m = 0;
            // memset(path, 0, (level_max+1)*sizeof(int));
            path[0] = m;
            break;
        case 2:
            m = 1;
            // memset(path, 0, (level_max+1)*sizeof(int));
            path[0] = m;
            break;
        default:
            m = 2+(tag_glbl-3)/twopowtwo[0];
            for (i=0; i<=level_max; i++)  /* Fortran grid[1:12] */
                path[i] = (tag_glbl-grid[m].tag_glbl)/twopowtwo[i] % 4;
            path[0] = m;
            break;
    }
    /* note that values in path[] is 0-based */

    ptr = grid[path[0]].nghbr[0];
    for (lvl=1; lvl<=level; lvl++) {
        if (ptr->dn[path[lvl]] != NULL)  /* Fortran dn[0:3] */
            ptr = ptr->dn[path[lvl]];
        else
            ptr = NULL;
    }

    free(path);
    free(twopowtwo);
    return ptr;
}


/*----< get_index() >--------------------------------------------------------*/
void get_index(grid_node *grid,
               int        level_max,
               int        level,
               int        iota,
               int       *subdomain_list, /*[subdomain_list_len]: 0 based*/
               int        subdomain_list_len,
               int        tag_glbl,
               int       *ix)  /* [3] values are 1-based */
{
    int i, el, nsd_glbl, nsd;
    grid_node *ptr;

    for (i=0; i<3; i++) ix[i] = -1;

    ptr = set_ptr(grid, level_max, level, tag_glbl);

    el = POWER2(level - iota);

    switch (tag_glbl) {
        case 1:
            ix[0] = 2;
            ix[1] = el + 2;
            /* nsd_glbl is global block ID, 0-based */
            nsd_glbl = POWER2(2*iota) - 1;
            break;
        case 2:
            ix[0] = el + 2;
            ix[1] = 2;
            nsd_glbl = (((POWER2(iota)+1)*(POWER2(iota)-1))/3)+9*POWER2(2*iota);
            break;
        default:
            ix[0] = 2 + (ptr->i - 1) % el;
            ix[1] = 2 + (ptr->j - 1) % el;
            nsd_glbl = (tag_glbl-3)/(POWER2(2*(level_max-iota)));
    }

    for (nsd=0; nsd<subdomain_list_len; nsd++) {
        if (subdomain_list[nsd] == nsd_glbl) {
            ix[2] = nsd + 1;
            break;
        }
    }
}
