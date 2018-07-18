/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: grid_subdomain.c 4601 2017-12-07 07:05:48Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "gcrm.h"
#include "gio.h"
#include "util.h"

static void set_sbdmn_map(int ***sbdmn_map, int iota, int lvl_map, int i,
                          int j, int pnl, int *nsd);
static void push_extended_list(extended_list_node **extended_list_head,
                               int nsd_glbl);

/*----< init_MODULE_grid_subdomain() >----------------------------------------*/
/* Data updated:
 *     grid_subdomain->level_threshold
 *     grid_subdomain->cell_min
 *     grid_subdomain->distribution_pattern
 *     grid_subdomain->l_sbdmn_pntgn_north
 *     grid_subdomain->l_sbdmn_pntgn_south
 *     grid_subdomain->l_sbdmn_north_pole
 *     grid_subdomain->l_sbdmn_south_pole
 */
void init_MODULE_grid_subdomain(MODULE_grid_subdomain *grid_subdomain,
                                int                    nsdm)
{
    grid_subdomain->level_threshold      = 2;
    grid_subdomain->cell_min             = 16;
    grid_subdomain->distribution_pattern = 2;

    grid_subdomain->l_sbdmn_pntgn_north = calloc_1D_int(nsdm);
    grid_subdomain->l_sbdmn_pntgn_south = calloc_1D_int(nsdm);

    grid_subdomain->l_sbdmn_north_pole = calloc_1D_int(nsdm);
    grid_subdomain->l_sbdmn_south_pole = calloc_1D_int(nsdm);
}

/*----< finalize_MODULE_grid_subdomain() >------------------------------------*/
void finalize_MODULE_grid_subdomain(MODULE_grid_params    *grid_params,
                                    MODULE_grid_subdomain *grid_subdomain)
{
    int i, level_max=grid_params->level_max;
    sbdmn_node *sbdmn = grid_params->sbdmn;

    for (i=0; i<=level_max; i++) {
        free_1D_int(sbdmn[i].proc);
        while (sbdmn[i].extended_list_head != NULL) {
            extended_list_node *tmpry = sbdmn[i].extended_list_head;
            sbdmn[i].extended_list_head = tmpry->next;
            tfree(tmpry);
        }
        if (sbdmn[i].lst != NULL) free_1D_int(sbdmn[i].lst);
    }

    free_1D_int(grid_subdomain->l_sbdmn_south_pole);
    free_1D_int(grid_subdomain->l_sbdmn_north_pole);
    free_1D_int(grid_subdomain->l_sbdmn_pntgn_south);
    free_1D_int(grid_subdomain->l_sbdmn_pntgn_north);
}

/*----< initialize_subdomain() >----------------------------------------------*/
void initialize_subdomain(MODULE_grid_params    *grid_params,
                          MODULE_grid_subdomain *grid_subdomain)
{
    int l_include, npe_wrld, rnk_wrld;
    int lvl, m, n, nsd, iota, i, j, k, pnl, lvl_map, p, q,
        ix[3], sbdmn_north[5], sbdmn_south[5], sbdmn_offset, rnk;

    int ***sbdmn_map, **map_ix;

    sbdmn_node *sbdmn = grid_params->sbdmn;
    int level_max = grid_params->level_max;

    MPI_Comm_size(MPI_COMM_WORLD, &npe_wrld);
    MPI_Comm_rank(MPI_COMM_WORLD, &rnk_wrld);

/*-----------------------------------------------------------------------
!  check that the gobal grid is built to sufficiently high resolution 
!  to allow subdomain blocks to be properly attached at resolutions
!  greater than level_glbl
!-----------------------------------------------------------------------*/
    if (grid_params->level_glbl < grid_params->sbdmn_iota) {
        printf(" initialize_subdomain : level_glbl < sbdmn_iota\n");
        printf(" sbdmn_iota                  = %d\n",grid_params->sbdmn_iota);
        printf(" level_glbl                  = %d\n",grid_params->level_glbl);
        ABORT
    }
/*-----------------------------------------------------------------------
!  check that gobal number of subdomains is evenly divisible by npe_wrld
!-----------------------------------------------------------------------*/
    if (grid_params->nsdm_glbl % npe_wrld) {
        printf(" initialize_subdomain : the global number of subdomains\n");
        printf(" it not divisble by npe_wrld                           \n");
        printf(" nsdm_glbl                  = %d\n",grid_params->nsdm_glbl);
        printf(" npe_wrld                   = %d\n",npe_wrld);
        ABORT
    }
/*-----------------------------------------------------------------------
!  the finest grid resolution
!-----------------------------------------------------------------------*/
    sbdmn[level_max].cell_max   = grid_params->cell_max;
    sbdmn[level_max].sbdmn_iota = grid_params->sbdmn_iota;
    sbdmn[level_max].nsdm       = grid_params->nsdm;
    sbdmn[level_max].nsdm_glbl  = grid_params->nsdm_glbl;

    sbdmn[level_max].lst = calloc_1D_int(grid_params->nsdm);

    if ((grid_params->cell_max-2)/grid_params->nsdm_glbl < grid_subdomain->cell_min) {
        printf(" initialize_subdomain : subdomain blocks contains less\n");
        printf(" than the minimum number of cells                     \n");
        printf(" level_max                  = %d\n",level_max);
        printf(" sbdmn[level_max].cell_max  = %d\n",sbdmn[level_max].cell_max);
        printf(" cell_min                   = %d\n",grid_subdomain->cell_min);
        printf(" (cell_max-2)/nsdm_glbl     = %d\n",(grid_params->cell_max-2)/grid_params->nsdm_glbl);
        ABORT
    }

    if (grid_subdomain->distribution_pattern == 1) { /* strided  */
        for (i=0; i<grid_params->nsdm_glbl; i+=npe_wrld)
            sbdmn[level_max].lst[i] = rnk_wrld + i;
            /* lst[] are block IDs, 0-based */
    }
    else if (grid_subdomain->distribution_pattern == 2) { /* contiguous */
        for (i=0; i<grid_params->nsdm; i++)
            sbdmn[level_max].lst[i] = grid_params->nsdm * rnk_wrld + i;
            /* lst[] are block IDs, 0-based */
    }
    else {
        printf("Error: grid_subdomain->distribution_pattern undefined\n");
        ABORT
    }
/*-----------------------------------------------------------------------
!  finest to level_threshold (see resolution_level_partition.nb or sbdmn.nb)
!-----------------------------------------------------------------------*/
    for (lvl=level_max-1; lvl>=grid_subdomain->level_threshold+1; lvl--) {
        sbdmn[lvl].cell_max = 2+10*POWER2(2*lvl);

        if (((sbdmn[lvl].cell_max-2)/sbdmn[lvl+1].nsdm_glbl) >= grid_subdomain->cell_min)
            sbdmn[lvl].sbdmn_iota = sbdmn[lvl+1].sbdmn_iota;
        else
            sbdmn[lvl].sbdmn_iota = sbdmn[lvl+1].sbdmn_iota-1;

        if (sbdmn[lvl].sbdmn_iota < 0) {
            printf(" initialize_subdomain : subdomain blocks contains less\n");
            printf(" than the minimum number of cells                     \n");
            printf(" lvl                  = %d\n",lvl);
            printf(" sbdmn[lvl].cell_max  = %d\n",sbdmn[lvl].cell_max);
            printf(" cell_min             = %d\n",grid_subdomain->cell_min);
            ABORT
        }

        sbdmn[lvl].nsdm_glbl = 10*POWER2(2*(sbdmn[lvl].sbdmn_iota));

        if (sbdmn[lvl].sbdmn_iota == sbdmn[lvl+1].sbdmn_iota) {
            sbdmn[lvl].nsdm = sbdmn[lvl+1].nsdm;
            sbdmn[lvl].lst = calloc_1D_int(sbdmn[lvl].nsdm);
            memcpy(sbdmn[lvl].lst, sbdmn[lvl+1].lst, sbdmn[lvl].nsdm * sizeof(int));
        } else {
            sbdmn[lvl].nsdm = 0;
            for (i=0; i<sbdmn[lvl+1].nsdm; i++)
                if ((sbdmn[lvl+1].lst[i]) % 4 == 0)
                    sbdmn[lvl].nsdm++;
            sbdmn[lvl].lst = calloc_1D_int(sbdmn[lvl].nsdm);
            n = 0;
            for (nsd=0; nsd<sbdmn[lvl+1].nsdm; nsd++)
                if ((sbdmn[lvl+1].lst[nsd]) % 4 == 0)
                    /* lst[] are block IDs, 0-based */
                    sbdmn[lvl].lst[n++] = sbdmn[lvl+1].lst[nsd] / 4;
        }
    }
/*-----------------------------------------------------------------------
!  IF (level_max == level_threshold) THEN STOP
!-----------------------------------------------------------------------*/
    if (level_max == grid_subdomain->level_threshold) {
        printf(" initialize_subdomain : level_max == level_threshold\n");
        printf("             all processes must have nsdm unique\n");
        printf("             subdomains on the finest resolution.\n");
        ABORT
    }
/*-----------------------------------------------------------------------
!  level_threshold to coarsest
!-----------------------------------------------------------------------*/
    for (lvl=grid_subdomain->level_threshold; lvl>=0; lvl--) {

        sbdmn[lvl].cell_max = 2+10*POWER2(2*lvl);

        sbdmn[lvl].sbdmn_iota = 0;

        if (rnk_wrld==0) {
            sbdmn[lvl].nsdm_glbl = 10;
            sbdmn[lvl].nsdm      = 10;
            sbdmn[lvl].lst       = calloc_1D_int(sbdmn[lvl].nsdm);
            for (i=0; i<10; i++)
                sbdmn[lvl].lst[i] = i; /* values of lst[] are 0-based */
        } else {
            sbdmn[lvl].nsdm_glbl = 10;
            sbdmn[lvl].nsdm      =  0;
            sbdmn[lvl].lst       = NULL;
            // sbdmn[lvl].lst = calloc_1D_int(sbdmn[lvl].nsdm+1);
        }
    }
/*-----------------------------------------------------------------------
! set other stuff
!-----------------------------------------------------------------------*/
    for (lvl=level_max; lvl>=0; lvl--) {

        iota = sbdmn[lvl].sbdmn_iota;
        /* nsd_north_glbl is subdomain block ID of north pole */
        sbdmn[lvl].nsd_north_glbl=       POWER2(2*iota) - 1;
        /* nsd_south_glbl is subdomain block ID of south pole */
        sbdmn[lvl].nsd_south_glbl=2*(1+7*POWER2(2*iota+1))/3 - 1;

        sbdmn[lvl].im = 2+POWER2(lvl-iota);
        sbdmn[lvl].jm = 2+POWER2(lvl-iota);

        if ((lvl>0) && (lvl-iota<=0)) {
            printf(" initialize_subdomain :: a subdomain block must contain at least four cells.\n");
            printf("                         lvl = %d\n",lvl);
            printf(" sbdmn[lvl].sbdmn_iota = %d\n",sbdmn[lvl].sbdmn_iota);
            ABORT
        }

/* north pole */
        for (n=sbdmn[lvl].nsdm-1; n>=0; n--) /* find the max index fulfill the condition */
            if (sbdmn[lvl].lst[n] == sbdmn[lvl].nsd_north_glbl)
                /* both lst[] and nsd_north_glbl are block IDs, 0-based */
                break;
        if (n >= 0) {
            sbdmn[lvl].l_agent_north = 1;
            sbdmn[lvl].nsd_north = n;   /* nsd_north is 0-based block ID */
        } else {
            sbdmn[lvl].l_agent_north = 0;
            sbdmn[lvl].nsd_north = -1;
        }
/* south pole */
        for (n=sbdmn[lvl].nsdm-1; n>=0; n--)
            if (sbdmn[lvl].lst[n] == sbdmn[lvl].nsd_south_glbl)
                /* both lst[] and nsd_south_glbl are block IDs, 0-based */
                break;
        if (n >= 0) {
            sbdmn[lvl].l_agent_south = 1;
            sbdmn[lvl].nsd_south = n;  /* nsd_south is 0-based block ID */
        } else {
            sbdmn[lvl].l_agent_south = 0;
            sbdmn[lvl].nsd_south = -1;
        }
    }

    grid_params->l_agent_north = sbdmn[level_max].l_agent_north;
    grid_params->l_agent_south = sbdmn[level_max].l_agent_south;

    grid_params->nsd_north = sbdmn[level_max].nsd_north;
    grid_params->nsd_south = sbdmn[level_max].nsd_south;

/*-----------------------------------------------------------------------
!  set l_sbdmn_pntgn_north and l_sbdmn_pntgn_south
!-----------------------------------------------------------------------*/
    for (i=0; i<5; i++)
        sbdmn_north[i] = POWER2(2*grid_params->sbdmn_iota) * (2*i);

    for (nsd=0; nsd<grid_params->nsdm; nsd++) {
        grid_subdomain->l_sbdmn_pntgn_north[nsd] = 0;
        for (i=0; i<5; i++)  /* values in lst[] is 0-based */
            if (sbdmn[level_max].lst[nsd] == sbdmn_north[i])
                break;
        if (i < 5) grid_subdomain->l_sbdmn_pntgn_north[nsd] = 1;
    }

    for (i=0; i<5; i++)
        sbdmn_south[i] = POWER2(2*grid_params->sbdmn_iota) * (2*i+1);

    for (nsd=0; nsd<grid_params->nsdm; nsd++) {
        grid_subdomain->l_sbdmn_pntgn_south[nsd] = 0;
        for (i=0; i<5; i++)  /* values in lst[] is 0-based */
            if (sbdmn[level_max].lst[nsd] == sbdmn_south[i])
                break;
        if (i < 5) grid_subdomain->l_sbdmn_pntgn_south[nsd] = 1;
    }

/*-----------------------------------------------------------------------
!  set l_sbdmn_north_pole and l_sbdmn_south_pole
!-----------------------------------------------------------------------*/
    for (i=0; i<5; i++)
        sbdmn_north[i] = POWER2(2*grid_params->sbdmn_iota) * (2*i+1) - 1;
        /*subdomains adjacent to north pole */

    for (nsd=0; nsd<grid_params->nsdm; nsd++) {
        grid_subdomain->l_sbdmn_north_pole[nsd] = 0;
        for (i=0; i<5; i++)  /* values in lst[] is 0-based */
            if (sbdmn[level_max].lst[nsd] == sbdmn_north[i])
                break;
        if (i < 5) grid_subdomain->l_sbdmn_north_pole[nsd] = 1;
    }

    sbdmn_offset = (POWER2(grid_params->sbdmn_iota)+1)*
                   (POWER2(grid_params->sbdmn_iota)-1)/3;
    for (i=0; i<5; i++)
        sbdmn_south[i] = sbdmn_offset + POWER2(2*grid_params->sbdmn_iota)
                                      * (2*i+1);
        /*subdomains adjacent to south pole */

    for (nsd=0; nsd<grid_params->nsdm; nsd++) {
        grid_subdomain->l_sbdmn_south_pole[nsd] = 0;
        for (i=0; i<5; i++)  /* values in lst[] is 0-based */
            if (sbdmn[level_max].lst[nsd] == sbdmn_south[i])
                break;
        if (i < 5) grid_subdomain->l_sbdmn_south_pole[nsd] = 1;
    }
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  SET EXTENDED SUBDOMAIN LISTS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    for (lvl=0; lvl<=level_max; lvl++) {
/*-----------------------------------------------------------------------
!  set subdomain map
!-----------------------------------------------------------------------*/
        iota = sbdmn[lvl].sbdmn_iota;
        sbdmn_map = calloc_3D_int(10, POWER2(iota)+2, POWER2(iota)+2);

        for (i=0; i<10; i++)
            for (j=0; j<POWER2(iota)+2; j++)
                for (k=0; k<POWER2(iota)+2; k++)
                    sbdmn_map[i][j][k] = -1;

        lvl_map=1;
        i=2;
        j=2;
        nsd=0;
        for (pnl=0; pnl<10; pnl++)  /* sbdmn_map[] is 0 based */
            set_sbdmn_map(sbdmn_map, sbdmn[lvl].sbdmn_iota,
                          lvl_map, i, j, pnl+1, &nsd);
/*-----------------------------------------------------------------------
!  wrap subdomain map
!-----------------------------------------------------------------------*/
        p = q = POWER2(iota)+1;
/*  northern hemisphere panels */
        for (pnl=0; pnl<=8; pnl+=2) {  /* careful on Fortran index MOD(pnl+7,10)+1 */
            for (i=1; i<q; i++)
                sbdmn_map[pnl][i][0] = sbdmn_map[(pnl+8)%10][q-1][p-i];
            for (i=1; i<p; i++)
                sbdmn_map[pnl][0][i] = sbdmn_map[(pnl+9)%10][q-1][i];
            for (i=1; i<q; i++)
                sbdmn_map[pnl][i][p] = sbdmn_map[(pnl+1)%10][i][1];
            for (i=1; i<p; i++)
                sbdmn_map[pnl][q][i] = sbdmn_map[(pnl+2)%10][q-i][1];
        }
/*  southern hemisphere panels */
        for (pnl=1; pnl<=9; pnl+=2) {
            for (i=1; i<q; i++)
                sbdmn_map[pnl][i][0] = sbdmn_map[(pnl+9)%10][i][p-1];
            for (i=1; i<p; i++)
                sbdmn_map[pnl][0][i] = sbdmn_map[(pnl+8)%10][q-i][p-1];
            for (i=1; i<q; i++)
                sbdmn_map[pnl][i][p] = sbdmn_map[(pnl+2)%10][1][p-i];
            for (i=1; i<p; i++)
                sbdmn_map[pnl][q][i] = sbdmn_map[(pnl+1)%10][1][i];
        }
/*-----------------------------------------------------------------------
!  set map index
!-----------------------------------------------------------------------*/
        map_ix = calloc_2D_int(sbdmn[lvl].nsdm_glbl, 3);

        for (j=0; j<sbdmn[lvl].nsdm_glbl; j++)
            for (i=0; i<3; i++)
                map_ix[j][i] = -1;
 
        for (pnl=0; pnl<10; pnl++)
            for (j=1; j<q; j++)
                for (i=1; i<p; i++) {
                    k = sbdmn_map[pnl][j][i];
                    map_ix[k][0] = i;
                    map_ix[k][1] = j;
                    map_ix[k][2] = pnl;
                }
/*-----------------------------------------------------------------------
!  set extended subdomain lists
!-----------------------------------------------------------------------*/
        sbdmn[lvl].extended_list_head = NULL;

        if (sbdmn[lvl].nsdm > 0) {
/*  push the local subdomains to the extended list */
            for (n=0; n<sbdmn[lvl].nsdm; n++)
                push_extended_list(&sbdmn[lvl].extended_list_head, sbdmn[lvl].lst[n]);

/*  push to extended list the local neighbor subdomains */
            int di[8] = {1, 1, 0,-1,-1,-1, 0, 1};
            int dj[8] = {0, 1, 1, 1, 0,-1,-1,-1};

            for (n=0; n<sbdmn[lvl].nsdm; n++) {
                for (i=0; i<3; i++)  /* values in lst[] and map_ix[] are 0-based */
                    ix[i] = map_ix[sbdmn[lvl].lst[n]][i];
                for (m=0; m<8; m++) {
                    if (sbdmn_map[ix[2]][ix[1]+dj[m]][ix[0]+di[m]] != -1)
                        push_extended_list(&sbdmn[lvl].extended_list_head,
                                           sbdmn_map[ix[2]][ix[1]+dj[m]][ix[0]+di[m]]);
                }
            }
/*  push to extended list the north pole subdomains */
            for (m=0; m<5; m++)
                sbdmn_north[m] = POWER2(2*iota)*(2*m+1) - 1;

            l_include = 0;
            for (m=0; m<5; m++)
                for (i=0; i<sbdmn[lvl].nsdm; i++)  /* values in lst[] is 0-based */
                    if (sbdmn[lvl].lst[i] == sbdmn_north[m]) {
                        l_include = 1;
                        m = 5;
                        break;
                    }

            if (l_include)
                for (m=0; m<5; m++)
                    push_extended_list(&sbdmn[lvl].extended_list_head, sbdmn_north[m]);

/*  push to extended list the south pole subdomains */
            sbdmn_offset = ((POWER2(iota)+1)*(POWER2(iota)-1))/3;
            for (m=0; m<5; m++)
                sbdmn_south[m] = sbdmn_offset + POWER2(2*iota)*(2*m+1);

            l_include = 0;
            for (m=0; m<5; m++)
                for (i=0; i<sbdmn[lvl].nsdm; i++)
                    if (sbdmn[lvl].lst[i] == sbdmn_south[m]) {
                        l_include = 1;
                        m = 5;
                        break;
                    }

            if (l_include)
                for (m=0; m<5; m++)
                    push_extended_list(&sbdmn[lvl].extended_list_head, sbdmn_south[m]);

        }
        free_3D_int(sbdmn_map);
        free_2D_int(map_ix);
    }

/*-----------------------------------------------------------------------
!  SET THE PROCESS NUMBER
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  the finest grid resolution
!-----------------------------------------------------------------------*/
    sbdmn[level_max].proc = calloc_1D_int(sbdmn[level_max].nsdm_glbl);
    for (i=0; i<sbdmn[level_max].nsdm_glbl; i++)
        sbdmn[level_max].proc[i] = RNK_NONEXISTENT;

    switch (grid_subdomain->distribution_pattern) {
        case 1: /* strided */
             for (nsd=0; nsd<sbdmn[level_max].nsdm_glbl; nsd++)
                 sbdmn[level_max].proc[nsd] = nsd % npe_wrld;  /* proc[] is 0-based */
             break;
        case 2: /* contiguous */
             for (rnk=0; rnk<npe_wrld; rnk++)
                 for (i=grid_params->nsdm*rnk; i<grid_params->nsdm*(rnk+1); i++)
                     sbdmn[level_max].proc[i] = rnk;
             break;
        default:
             printf("Error: grid_subdomain->distribution_pattern is invalid\n");
             ABORT
             break;
    }
/*-----------------------------------------------------------------------
!  finest to level_threshold (see sbdmn.nb)
!-----------------------------------------------------------------------*/
    for (lvl=level_max-1; lvl>=grid_subdomain->level_threshold+1; lvl--) {
        sbdmn[lvl].proc = calloc_1D_int(sbdmn[lvl].nsdm_glbl);
        for (i=0; i<sbdmn[lvl].nsdm_glbl; i++)
            sbdmn[lvl].proc[i] = RNK_NONEXISTENT;
        if (sbdmn[lvl].sbdmn_iota == sbdmn[lvl+1].sbdmn_iota) {
            for (i=0; i<sbdmn[lvl].nsdm_glbl; i++)
                sbdmn[lvl].proc[i] = sbdmn[lvl+1].proc[i];
        }
        else {
            for (nsd=0; nsd<sbdmn[lvl].nsdm_glbl; nsd++)
                sbdmn[lvl].proc[nsd] = sbdmn[lvl+1].proc[4*nsd];
        }
    }
/*-----------------------------------------------------------------------
!  level_threshold to coarsest
!-----------------------------------------------------------------------*/
    for (lvl=0; lvl<=grid_subdomain->level_threshold; lvl++)
        sbdmn[lvl].proc = calloc_1D_int(sbdmn[lvl].nsdm_glbl);
}

/*----< set_sbdmn_map() >-----------------------------------------------------*/
static
void set_sbdmn_map(int ***sbdmn_map,
                   int    iota,
                   int    lvl_map,
                   int    i,
                   int    j,
                   int    pnl,
                   int   *nsd)
{
    int del;

    del = POWER2(iota - lvl_map);

    if (lvl_map <= iota) {
        set_sbdmn_map(sbdmn_map, iota, lvl_map+1, i    , j    , pnl, nsd);
        set_sbdmn_map(sbdmn_map, iota, lvl_map+1, i+del, j    , pnl, nsd);
        set_sbdmn_map(sbdmn_map, iota, lvl_map+1, i+del, j+del, pnl, nsd);
        set_sbdmn_map(sbdmn_map, iota, lvl_map+1, i    , j+del, pnl, nsd);
    }
    else {
        sbdmn_map[pnl-1][j-1][i-1] = *nsd;
        (*nsd)++;
    }
}

/*----< push_extended_list() >------------------------------------------------*/
static
void push_extended_list(extended_list_node **extended_list_head,
                        int                  nsd_glbl)  /* global block index */
{
    extended_list_node *tmpry;

    /* the input nsd_glbl is 0 based */

    if (*extended_list_head == NULL) {
        *extended_list_head = (extended_list_node*)tmalloc(sizeof(extended_list_node));
        tmpry = *extended_list_head;
        tmpry->nsd_glbl = nsd_glbl;
        tmpry->next = NULL;
    }
    else {
/* search for existing nsd_glbl */
        tmpry = *extended_list_head;
        while (tmpry != NULL) {
            if (tmpry->nsd_glbl == nsd_glbl) return;
            tmpry = tmpry->next;
        }
/* add to the end of the list */
        tmpry = *extended_list_head;
        while (tmpry->next != NULL)
            tmpry = tmpry->next;
        tmpry->next = (extended_list_node*)tmalloc(sizeof(extended_list_node));
        tmpry = tmpry->next;
        tmpry->nsd_glbl = nsd_glbl;
        tmpry->next = NULL;
    }
}

/*----< get_proc() >---------------------------------------------------------*/
int get_proc(sbdmn_node *sbdmn,  /* Fortran sbdmn(0:level_max) */
             int         level,
             int         tag_glbl)

{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PURPOSE : given a grid resolution level and global tag tag_glbl 
!           return the process rnk that owns that grid point
!.......................................................................*/
    int nsd=0;

    if (tag_glbl >= 3)
        nsd = (tag_glbl-3)/(sbdmn[level].im * sbdmn[level].jm);
    else if (tag_glbl == 1)
        nsd = sbdmn[level].nsd_north_glbl;
    else if (tag_glbl == 2)
        nsd = sbdmn[level].nsd_south_glbl;

    /* sbdmn[*].proc = calloc_1D_int(sbdmn[*].nsdm_glbl); */
    return sbdmn[level].proc[nsd];
}

/*----< get_big_block_index() ------------------------------------------------*/
/*
   FUNCTION get_big_block_index (lvl,iota,nsd_glbl) RESULT (ix)
!.......................................................................
! PURPOSE : find the index of nsd_glbl within the big block data structure
! 
!   INPUT : lvl      -> the grid resolution level (r value). typically
!                       this is level_max.
!           iota     -> determines the number of subdomain blocks of the
!                       horizontal domain decomposition. typically this
!                       is sbdmn_iota.
!           nsd_glbl -> global number of subdomain.  typically this is
!                       nsd_glbl
!
!  OUTPUT : ix(3)    -> big block index or panel index.  this number is 
!                       in the range [1,10]
!           ix(1:2)  -> the i and j starting index within the big block.
!                       each of these numbers is in the range [1,2**lvl]
!
! EXAMPLE : suppose the big block global array is called x1 
!           and the local subdomain is called x0. then to copy x0 into the x1
!
!    ix(:) = CALL get_big_block_index (level_max,sbdmn_iota,nsd_glbl)
!    x1(ix(1):ix(1)+im-1,ix(2):ix(2)+jm-1,ix(3)) = x0(2:im-1,2:jm-1)
*/
int get_big_block_index(int lvl,
                        int iota,
                        int nsd_glbl, /* global block ID, 0-based value */
                        int *ix)      /* [3] values are 1-based */
{
    int i;
    int nblkm = POWER2(2 * iota); /* the number of subdomains in one big block */

    ix[2] = (nsd_glbl/nblkm) + 1; /* the global big block index [1,10] */
    ix[0] = ix[1] = 0;

    for (i=1; i<=iota; i++) {
        nsd_glbl %= nblkm;
        nblkm    /= 4;

        switch (nsd_glbl/nblkm) {
            case 0:
                break;
            case 1:
                ix[0] += POWER2(iota-i);
                break;
            case 2:
                ix[0] += POWER2(iota-i);
                ix[1] += POWER2(iota-i);
                break;
            case 3:
                ix[1] += POWER2(iota-i);
                break;
            default: break;
        }
   }

   ix[0] *= POWER2(lvl-iota);
   ix[1] *= POWER2(lvl-iota);

   ix[0]++;
   ix[1]++;

   return 0;
}

