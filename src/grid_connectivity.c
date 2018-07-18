/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: grid_connectivity.c 4425 2017-08-28 03:42:05Z wkliao $
 */

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* memset() */
#include <mpi.h>

#include "grid_params.h"
#include "util.h"

static void connectivity_1(MODULE_grid_params *grid_params, grid_node *ptr);
static void connectivity_2(grid_node *ptr);
static int build(MODULE_grid_params *grid_params, grid_node *ptr, int n,
                 int d_tag_glbl);
static int get_direction(grid_node *ptr1, grid_node *ptr2);
static void set_real(sbdmn_node *sbdmn, int level_max, grid_node *ptr);
static grid_node* spiral_path(int sprl_dpth, grid_node *ptr0, int tag_optional);
static void set_ghst(grid_node *ptr, int level_max);
static void set_proc(grid_node  *ptr, int level_max, sbdmn_node *sbdmn);
static void set_path(int lvl, int *count_next, int *count_real, int *count_ghst,
                     grid_node **ptr_next, grid_node **ptr_real,
                     grid_node **ptr_ghst, grid_node *ptr);
static void set_ix_2D(grid_node *ptr, grid_node *grid, int level_max,
                      sbdmn_node *sbdmn);
static grid_node* twist(grid_node *ptr0, grid_node *ptr1, int ntwist);
static void set_nghbr_lst(grid_node *ptr);
static void free_nghbr_lst(grid_node *ptr);

/*----< initialize_grid_connectivity() >--------------------------------------*/
void initialize_grid_connectivity(MODULE_grid_params *grid_params,
                                  char               *communicator_name)
{
    int i, n, m, lvl, count_next, count_real, count_ghst, tag_locl, level_max;
    grid_node *path_head, *grid;
    grid_node *ptr, *ptr_next, *ptr_real, *ptr_ghst;

    grid = grid_params->grid;
    level_max = grid_params->level_max;
    path_head = (grid_node*) tcalloc(1, sizeof(grid_node));

/*-----------------------------------------------------------------------
!  set level 0
!-----------------------------------------------------------------------*/
    /* should grid[].i and grid[].j start from 0 index ? */
    grid[0].i = 1;
    grid[0].j = 2;
    grid[1].i = 2;
    grid[1].j = 1;
    for (i=2; i<12; i++) {
        grid[i].i = 1;
        grid[i].j = 1;
    }

    for (i=0; i<12; i++) {
        grid[i].level  = 0;
        grid[i].l_real = 0;
        grid[i].l_ghst = 0;

        grid[i].l_north = (i%2 == 0) ? 1 : 0;
        grid[i].l_south = (i%2 == 0) ? 0 : 1;

        grid[i].l_pentagon = 1;

        grid[i].l_pole_north = 0;
        grid[i].l_pole_south = 0;

        grid[i].l_pentagon_north = (i > 1 && i%2 == 0) ? 1 : 0;
        grid[i].l_pentagon_south = (i > 1 && i%2  > 0) ? 1 : 0;

        grid[i].proc = RNK_NONEXISTENT;
    }
    grid[0].l_pole_north = 1;
    grid[1].l_pole_south = 1;
/*-----------------------------------------------------------------------
!  nullify resolution links at level 0
!-----------------------------------------------------------------------*/
    for (n=0; n<12; n++) { /* assume SIZE (grid) == 12 */
        for (m=0; m<4; m++)
            grid[n].dn[m] = NULL;
        grid[n].up = NULL;
    }

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  SET CONNECTIONS AT LEVEL 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!-----------------------------------------------------------------------
!  link the north pole to mid-latitude pentagons
!-----------------------------------------------------------------------*/
    grid[0].nghbr[0] = &grid[ 0];
    grid[0].nghbr[1] = &grid[ 2];
    grid[0].nghbr[2] = &grid[ 4];
    grid[0].nghbr[3] = &grid[ 6];
    grid[0].nghbr[4] = &grid[ 8];
    grid[0].nghbr[5] = &grid[10];
    grid[0].nghbr[6] = NULL;

/*-----------------------------------------------------------------------
!  link the south pole to mid-latitude pentagons
!-----------------------------------------------------------------------*/
    grid[1].nghbr[0] = &grid[ 1];
    grid[1].nghbr[1] = &grid[11];
    grid[1].nghbr[2] = &grid[ 9];
    grid[1].nghbr[3] = &grid[ 7];
    grid[1].nghbr[4] = &grid[ 5];
    grid[1].nghbr[5] = &grid[ 3];
    grid[1].nghbr[6] = NULL;

/*-----------------------------------------------------------------------
!  link the northern hemosphere pentagons
!-----------------------------------------------------------------------*/
    for (n=0; n<=8; n+=2) {
        grid[n+2].nghbr[0] = &grid[n+2];
        grid[n+2].nghbr[1] = &grid[2 + (n+1)%10];
        grid[n+2].nghbr[2] = &grid[2 + (n+2)%10];
        grid[n+2].nghbr[3] = &grid[0];
        grid[n+2].nghbr[4] = NULL;
        grid[n+2].nghbr[5] = &grid[2 + (n+8)%10];
        grid[n+2].nghbr[6] = &grid[2 + (n+9)%10];
    }
/*-----------------------------------------------------------------------
!  link the southern hemosphere pentagons
!-----------------------------------------------------------------------*/
    for (n=0; n<=8; n+=2) {
        grid[n+3].nghbr[0] = &grid[n+3];
        grid[n+3].nghbr[1] = &grid[1];
        grid[n+3].nghbr[2] = &grid[2 + (n+ 3)%10];
        grid[n+3].nghbr[3] = &grid[2 + (n+ 2)%10];
        grid[n+3].nghbr[4] = &grid[2 + (n+10)%10];
        grid[n+3].nghbr[5] = &grid[2 + (n+ 9)%10];
        grid[n+3].nghbr[6] = NULL;
    }
/*-----------------------------------------------------------------------
!  set global tags at level 0
!-----------------------------------------------------------------------*/
    grid[0].tag_glbl = 1;
    grid[1].tag_glbl = 2;
    for (n=2; n<12; n++)
        grid[n].tag_glbl = 3+(POWER2(2*level_max))*(n-2);

/*-----------------------------------------------------------------------
!  set spiral seach tags at level 0
!-----------------------------------------------------------------------*/
    for (n=0; n<12; n++)
        grid[n].tag_sprl = grid_params->tag_nonexistent;

    for (n=0; n<12; n++)
        connectivity_1(grid_params, grid[n].nghbr[0]);

    for (n=2; n<12; n++)
        connectivity_2(grid[n].nghbr[0]);

    for (n=0; n<12; n++)
        set_real(grid_params->sbdmn, level_max, grid[n].nghbr[0]);

    for (n=0; n<12; n++)
        set_ghst(grid[n].nghbr[0], level_max);

/*-----------------------------------------------------------------------
!  make a path through all and real and ghost cells associated
!  with the local process and record the size of local arrays
!-----------------------------------------------------------------------*/
    for (lvl=0; lvl<=level_max; lvl++) {
        path_head->next_next = NULL;
        path_head->next_real = NULL;
        path_head->next_ghst = NULL;

        ptr_next = path_head;
        ptr_real = path_head;
        ptr_ghst = path_head;

        count_next = 0;
        count_real = 0;
        count_ghst = 0;

        for (n=0; n<12; n++)
            set_path(lvl,
                     &count_next,
                     &count_real,
                     &count_ghst,
                     &ptr_next,
                     &ptr_real,
                     &ptr_ghst,
                     grid[n].nghbr[0]);

        grid_params->path_next[lvl] = path_head->next_next;
        grid_params->path_real[lvl] = path_head->next_real;
        grid_params->path_ghst[lvl] = path_head->next_ghst;

        grid_params->nm_lvl     [lvl] = count_next;
        grid_params->nm_real_lvl[lvl] = count_real;
        grid_params->nm_ghst_lvl[lvl] = count_ghst;
    }
    grid_params->nm      = grid_params->nm_lvl     [level_max];
    grid_params->nm_real = grid_params->nm_real_lvl[level_max];
    grid_params->nm_ghst = grid_params->nm_ghst_lvl[level_max];

/*-----------------------------------------------------------------------
!  set the local tags
!-----------------------------------------------------------------------*/
    for (lvl=0; lvl<=level_max; lvl++) {
        tag_locl = 0;
        ptr = grid_params->path_next[lvl];
        while (ptr != NULL) {
            ptr->tag_locl = tag_locl;
            ptr = ptr->next_next;
            tag_locl++;  /* local tag (ID), 0-based */
        }
    }

    for (n=0; n<12; n++)
        set_ix_2D(grid[n].nghbr[0], grid_params->grid, level_max,
                  grid_params->sbdmn);

    for (n=0; n<12; n++)
        set_nghbr_lst(grid[n].nghbr[0]);

    //   CALL communicate_processor (communicator_name,grid(:))

    for (n=0; n<12; n++)
        set_proc(grid[n].nghbr[0], level_max, grid_params->sbdmn);

    // CALL report_grid_node_total (communicator_name)
    tfree(path_head);
}

/*----< free_connectivity_1() >-----------------------------------------------*/
static
void free_connectivity_1(grid_node *ptr)
{
    int n;
    if (ptr == NULL) return;

    for (n=0; n<4; n++) {
        if (ptr->dn[n] != NULL) {
            free_connectivity_1(ptr->dn[n]);
            tfree(ptr->dn[n]);
            ptr->dn[n] = NULL;
        }
    }
}

/*----< finalize_grid_connectivity() >----------------------------------------*/
void finalize_grid_connectivity(MODULE_grid_params *grid_params)
{
    int i;

    for (i=0; i<12; i++)
        free_nghbr_lst(grid_params->grid[i].nghbr[0]);

    for (i=0; i<12; i++) {
        free_connectivity_1(grid_params->grid[i].nghbr[0]);

        if (grid_params->grid[i].l_real) {
            if (grid_params->grid[i].nghbr_lst != NULL)
                free_2D_int(grid_params->grid[i].nghbr_lst);
            grid_params->grid[i].nghbr_lst = NULL;
        }
    }
}

/*----< connectivity_1() >----------------------------------------------------*/
static
void connectivity_1(MODULE_grid_params *grid_params,
                    grid_node          *ptr)
{
    int n, dn_max, m, d_tag_glbl, level_max;

    level_max = grid_params->level_max;

    if (ptr->level < level_max) {

        if (ptr->l_pole_north || ptr->l_pole_south)
            dn_max = 0;
        else
            dn_max = 3;

        d_tag_glbl = POWER2(2*(level_max-(ptr->level+1)));

        for (n=0; n<=dn_max; n++) {
            if (build(grid_params, ptr, n, d_tag_glbl)) {
                ptr->dn[n] = (grid_node*) tcalloc(1, sizeof(grid_node));

                if (ptr->dn[n] == NULL) {
                    printf("Error: malloc grid_node\n");
                    ABORT
                }

                grid_params->grid_node_total++;

                if (grid_params->grid_node_size * grid_params->grid_node_total >
                    grid_params->grid_node_memory_max) {
                    printf(" connectivity_1 :: TOO MUCH GRID NODE MEMORY ALLOCATED\n");
                    printf(" grid_node_total        = %lld",grid_params->grid_node_total);
                    printf(" grid_node_size         = %lld",grid_params->grid_node_size);
                    printf(" grid_node_memory_total = %lld",grid_params->grid_node_size*grid_params->grid_node_total);
                    printf(" grid_node_memory_max   = %lld",grid_params->grid_node_memory_max);
                    ABORT
                }

                ptr->dn[n]->up = ptr;

                ptr->dn[n]->nghbr[0] = ptr->dn[n];
                for (m=1; m<7; m++) ptr->dn[n]->nghbr[m] = NULL;
                for (m=0; m<4; m++) ptr->dn[n]->dn[m]    = NULL;

                ptr->dn[n]->level = ptr->level + 1;

                ptr->dn[n]->tag_glbl = ptr->tag_glbl + d_tag_glbl*n;
                ptr->dn[n]->tag_locl = grid_params->tag_nonexistent;
                ptr->dn[n]->tag_sprl = grid_params->tag_nonexistent;

                ptr->dn[n]->l_pole_north = ptr->l_pole_north;
                ptr->dn[n]->l_pole_south = ptr->l_pole_south;

                ptr->dn[n]->l_north = ptr->l_north;
                ptr->dn[n]->l_south = ptr->l_south;

                if (n == 0) {
                    ptr->dn[n]->l_pentagon       = ptr->l_pentagon;
                    ptr->dn[n]->l_pentagon_north = ptr->l_pentagon_north;
                    ptr->dn[n]->l_pentagon_south = ptr->l_pentagon_south;
                }
                else {
                    ptr->dn[n]->l_pentagon       = 0;
                    ptr->dn[n]->l_pentagon_north = 0;
                    ptr->dn[n]->l_pentagon_south = 0;
                }

                ptr->dn[n]->l_real = 0;
                ptr->dn[n]->l_ghst = 0;

                ptr->dn[n]->i = 2*(ptr->i-1)+1;
                ptr->dn[n]->j = 2*(ptr->j-1)+1;

                if (n == 1)
                    ptr->dn[n]->i++;
                else if (n == 2) {
                    ptr->dn[n]->i++;
                    ptr->dn[n]->j++;
                }
                else if (n == 3)
                   ptr->dn[n]->j++;

                ptr->dn[n]->proc        = RNK_NONEXISTENT;
                ptr->dn[n]->l_point_set = 0;
                ptr->dn[n]->point[0] = -1.0;
                ptr->dn[n]->point[1] = -1.0;
                ptr->dn[n]->point[2] = -1.0;

                ptr->next_next = NULL;
                ptr->next_real = NULL;
                ptr->next_ghst = NULL;
                ptr->next_sprl = NULL;

                connectivity_1(grid_params, ptr->dn[n]);
            }
        }
    }
}


/*----< build() >-------------------------------------------------------------*/
static
int build(MODULE_grid_params *grid_params,
          grid_node          *ptr,
          int                 n,
          int                 d_tag_glbl)
{
    int lvl, tag_range_min, tag_range_max, iota,tag_glbl, level_max;
    extended_list_node *tmpry;

    if (ptr->tag_glbl == 1 || ptr->tag_glbl == 2) return 1;

    level_max = grid_params->level_max;

    if (ptr->level < grid_params->level_glbl ||
        ptr->level == level_max-1)
        return 1;

/*-----------------------------------------------------------------------
!  build from outside
!-----------------------------------------------------------------------*/
    tag_range_min = ptr->tag_glbl+(POWER2(2*(level_max-ptr->level))/4)*(n  );
    tag_range_max = ptr->tag_glbl+(POWER2(2*(level_max-ptr->level))/4)*(n+1);

    for (lvl=ptr->level+1; lvl<=level_max; lvl++) {
        tmpry = grid_params->sbdmn[lvl].extended_list_head;
        while (tmpry != NULL) {
            iota = grid_params->sbdmn[lvl].sbdmn_iota;
            tag_glbl = 3+(POWER2(2*(level_max-iota)))*(tmpry->nsd_glbl);
            /* nsd_glbl is block index, 0-based */
            if (tag_range_min <= tag_glbl && tag_glbl < tag_range_max)
                return 1;
            tmpry = tmpry->next;
        }
    }
/*-----------------------------------------------------------------------
!  build from inside
!-----------------------------------------------------------------------*/
    tag_glbl = ptr->tag_glbl + d_tag_glbl*n;
    for (lvl=ptr->level+1; lvl<=level_max; lvl++) {
        tmpry = grid_params->sbdmn[lvl].extended_list_head;
        while (tmpry != NULL) {
            iota = grid_params->sbdmn[lvl].sbdmn_iota;
            tag_range_min = 3+(POWER2(2*(level_max-iota)))*(tmpry->nsd_glbl);
            tag_range_max = 3+(POWER2(2*(level_max-iota)))*(tmpry->nsd_glbl+1);
            if (tag_range_min <= tag_glbl && tag_glbl < tag_range_max)
                return 1;
            tmpry = tmpry->next;
        }
    }
    return 0;
}

#define LINK(n1, n2, ptr1, ptr2) { \
    (ptr1)->nghbr[n1] = ptr2;      \
    (ptr2)->nghbr[n2] = ptr1;      \
}

/*----< connectivity_2() >----------------------------------------------------*/
static
void connectivity_2(grid_node *ptr)
{
    int n;

/* note the Fortran index of nghbr and dn starts from 0
    grid_node *nghbr[7];  [0:6]
    grid_node *dn[4];     [0:3]
*/

/*-----------------------------------------------------------------------
!  type A links.  link 1, link 2 and link 3
!-----------------------------------------------------------------------*/
    if (ptr->dn[0] != NULL) {
        if (ptr->dn[1] != NULL) LINK(1, 4, ptr->dn[0], ptr->dn[1])
        if (ptr->dn[2] != NULL) LINK(2, 5, ptr->dn[0], ptr->dn[2])
        if (ptr->dn[3] != NULL) LINK(3, 6, ptr->dn[0], ptr->dn[3])
    }
/*-----------------------------------------------------------------------
!  type B links.
!-----------------------------------------------------------------------*/
    if (ptr->dn[1] != NULL) {
        /* link 4 */
        if (ptr->nghbr[1] != NULL) {
            if (ptr->nghbr[1]->dn[0] != NULL) {
                n = get_direction(ptr->nghbr[1], ptr->nghbr[0]);
                LINK(1, n, ptr->dn[1],ptr->nghbr[1]->dn[0])
            }
        }
        /* link 5 */
        if (ptr->l_south && ptr->i == POWER2(ptr->level)) {
             if (ptr->nghbr[2]        != NULL &&
                 ptr->nghbr[2]->dn[1] != NULL)
                 LINK(2, 6, ptr->dn[1], ptr->nghbr[2]->dn[1])
        } else {
             if (ptr->nghbr[1]        != NULL &&
                 ptr->nghbr[1]->dn[3] != NULL)
                 LINK(2, 5, ptr->dn[1], ptr->nghbr[1]->dn[3])
        }
        /* link 6 */
        if (ptr->dn[2] != NULL)
            LINK(3, 6, ptr->dn[1], ptr->dn[2])
    }
/*-----------------------------------------------------------------------
!  type C links.
!-----------------------------------------------------------------------*/
    if (ptr->dn[2] != NULL) {
        /* link 7 */
        if (ptr->l_south && ptr->i == POWER2(ptr->level)) {
            if (ptr->nghbr[2]        != NULL &&
                ptr->nghbr[2]->dn[1] != NULL)
                LINK(1, 5, ptr->dn[2], ptr->nghbr[2]->dn[1])
        } else {
            if (ptr->nghbr[1]        != NULL &&
                ptr->nghbr[1]->dn[3] != NULL)
                LINK(1, 4, ptr->dn[2], ptr->nghbr[1]->dn[3])
        }
        /* link 8 */
        if (ptr->nghbr[2]        != NULL &&
            ptr->nghbr[2]->dn[0] != NULL) {
            n = get_direction(ptr->nghbr[2], ptr->nghbr[0]);
            LINK(2, n, ptr->dn[2], ptr->nghbr[2]->dn[0])
        }
        /* link 9 */
        if (ptr->l_north && ptr->j == POWER2(ptr->level)) {
            if (ptr->nghbr[2]        != NULL &&
                ptr->nghbr[2]->dn[3] != NULL)
                LINK(3, 5, ptr->dn[2], ptr->nghbr[2]->dn[3])
        } else {
            if (ptr->nghbr[3]        != NULL &&
                ptr->nghbr[3]->dn[1] != NULL)
                LINK(3, 6, ptr->dn[2], ptr->nghbr[3]->dn[1])
        }
    }
/*-----------------------------------------------------------------------
!  type D links.
!-----------------------------------------------------------------------*/
    if (ptr->dn[3] != NULL) {
        /* link 10 */
        if (ptr->dn[2] != NULL)
            LINK(1, 4, ptr->dn[3], ptr->dn[2])
        /* link 11 */
        if (ptr->l_north && ptr->j == POWER2(ptr->level)) {
            if (ptr->nghbr[2]        != NULL &&
                ptr->nghbr[2]->dn[3] != NULL)
                LINK(2, 4, ptr->dn[3], ptr->nghbr[2]->dn[3])
        } else {
            if (ptr->nghbr[3]        != NULL &&
                ptr->nghbr[3]->dn[1] != NULL)
                LINK(2, 5, ptr->dn[3], ptr->nghbr[3]->dn[1])
        }
        /* link 12 */
        if (ptr->nghbr[3]        != NULL &&
            ptr->nghbr[3]->dn[0] != NULL) {
            n = get_direction(ptr->nghbr[3], ptr->nghbr[0]);
            LINK(3, n, ptr->dn[3], ptr->nghbr[3]->dn[0])
        }
    }

    for (n=0; n<4; n++)
        if (ptr->dn[n] != NULL)
            connectivity_2(ptr->dn[n]);
}

/*----< get_direction() >-----------------------------------------------------*/
static
int get_direction(grid_node *ptr1,
                  grid_node *ptr2)
{
    int n;
    for (n=1; n<7; n++) { /* grid_node *nghbr[7]; [0:6] */
        if (ptr1->nghbr[n] != NULL &&
            ptr1->nghbr[n]->tag_glbl == ptr2->tag_glbl)
            return n;
    }

    printf(" ERROR get_direction : ptr1->tag_glbl = %d\n",ptr1->tag_glbl);
    printf("                       ptr2->tag_glbl = %d\n",ptr2->tag_glbl);

    for (n=1; n<7; n++) { /* grid_node *nghbr[7]; [0:6] */
        if (ptr1->nghbr[n] != NULL)
            printf("      ptr1->nghbr[%d] = %d\n",n,ptr1->nghbr[n]->tag_glbl);
        else
            printf("      ptr1->nghbr[%d] = NULL\n",n);
    }
    return -1;
}

/*----< set_real() >----------------------------------------------------------*/
static
void set_real(sbdmn_node *sbdmn,     /* Fortran sbdmn(0:level_max) */
              int         level_max,
              grid_node  *ptr)
{
    int iota, nsd_glbl, n;

    switch (ptr->tag_glbl) {
        case(1):
            if (sbdmn[ptr->level].l_agent_north) ptr->l_real = 1;
            break;
        case(2) :
            if (sbdmn[ptr->level].l_agent_south) ptr->l_real = 1;
            break;
        default:
            iota = sbdmn[ptr->level].sbdmn_iota;
            /* nsd_glbl is global block ID, 0-based */
            nsd_glbl = (ptr->tag_glbl-3)/(POWER2(2*(level_max-iota)));
            ptr->l_real = 0;
            for (n=0; n<sbdmn[ptr->level].nsdm; n++)
                /* lst[] is the list of block IDs, 0-based */
                if (sbdmn[ptr->level].lst[n] == nsd_glbl) {
                    ptr->l_real = 1;
                    break;
                }
    }

    for (n=0; n<4; n++)  /* Fortran dn(0:3) */
        if (ptr->dn[n] != NULL)
            set_real(sbdmn, level_max, ptr->dn[n]);
}

/*----< set_ghst() >----------------------------------------------------------*/
static
void set_ghst(grid_node *ptr,
              int        level_max)
{
    int sprl_dpth, n;
    grid_node *ptr1, *ptr_sprl;

    if (ptr->l_real) {

        if (ptr->level == level_max)
            sprl_dpth = 4;
        else
            sprl_dpth = 1;

        ptr_sprl = spiral_path(sprl_dpth, ptr, -999);
        ptr1     = ptr_sprl;
        while (ptr1 != NULL) {
            if (ptr1->l_real == 0)
                ptr1->l_ghst = 1;
            ptr1 = ptr1->next_sprl;
        }
    }

    for (n=0; n<4; n++)
        if (ptr->dn[n] != NULL)
            set_ghst(ptr->dn[n], level_max);
}

/*----< set_proc() >----------------------------------------------------------*/
static
void set_proc(grid_node  *ptr,
              int         level_max,
              sbdmn_node *sbdmn)
{
    int level, tag_glbl, tag_tmpry, nsd, n;

    level    = ptr->level;
    tag_glbl = ptr->tag_glbl;

    if (tag_glbl >= 3) {
        tag_tmpry = (tag_glbl-3)/(POWER2(2*(level_max-level))) + 3;
        nsd = (tag_tmpry-3)/((sbdmn[level].im-2)*(sbdmn[level].jm-2));
    } else {
        if (tag_glbl == 1) nsd = sbdmn[level].nsd_north_glbl;
        if (tag_glbl == 2) nsd = sbdmn[level].nsd_south_glbl;
    }

    /* sbdmn[level].proc[] is allocated with size of sbdmn[level].nsdm_glbl */
    ptr->proc = sbdmn[level].proc[nsd]; /* nsd is subdomain block ID, 0-based */

    for (n=0; n<4; n++)
        if (ptr->dn[n] != NULL)
            set_proc(ptr->dn[n], level_max, sbdmn);
}

/*----< set_path() >----------------------------------------------------------*/
static
void set_path(int         lvl,
              int        *count_next,
              int        *count_real,
              int        *count_ghst,
              grid_node **ptr_next,
              grid_node **ptr_real,
              grid_node **ptr_ghst,
              grid_node  *ptr)
{
    int n;

    if (ptr->level == lvl &&
        (ptr->l_real == 1 || ptr->l_ghst == 1)) {
        (*count_next)++;
        (*ptr_next)->next_next = ptr;
        (*ptr_next)            = ptr;
    }

    if (ptr->level == lvl && ptr->l_real == 1) {
        (*count_real)++;
        (*ptr_real)->next_real = ptr;
        (*ptr_real)            = ptr;
    }

    if (ptr->level == lvl && ptr-> l_ghst == 1) {
        (*count_ghst)++;
        (*ptr_ghst)->next_ghst = ptr;
        (*ptr_ghst)            = ptr;
    }

    for (n=0; n<4; n++)
        if (ptr->dn[n] != NULL)
            set_path(lvl, count_next, count_real, count_ghst,
                     ptr_next, ptr_real, ptr_ghst, ptr->dn[n]);
}

/*----< set_ix_2D() >---------------------------------------------------------*/
static
void set_ix_2D(grid_node  *ptr,
               grid_node  *grid,
               int         level_max,
               sbdmn_node *sbdmn)
{
    int n, lvl;

    lvl = ptr->level;

    get_index(grid, level_max, lvl, sbdmn[lvl].sbdmn_iota,
              sbdmn[lvl].lst, sbdmn[lvl].nsdm, ptr->tag_glbl, ptr->ix_2D);

    for (n=0; n<4; n++)
        if (ptr->dn[n] != NULL)
            set_ix_2D(ptr->dn[n], grid, level_max, sbdmn);
}

/*----< free_nghbr_lst() >----------------------------------------------------*/
static
void free_nghbr_lst(grid_node *ptr)
{
    int n;

    for (n=0; n<4; n++)
        if (ptr->dn[n] != NULL)
            free_nghbr_lst(ptr->dn[n]);

    if (ptr->l_real)
        free_2D_int(ptr->nghbr_lst);
    ptr->nghbr_lst = NULL;
}

/*----< set_nghbr_lst() >-----------------------------------------------------*/
static
void set_nghbr_lst(grid_node *ptr)
{
    int n, ix5[5]={1,2,3,4,5}, ix6[6]={1,2,3,4,5,6};
    grid_node *ptr1, *ptr2;

    if (ptr->l_real) {
        if (ptr->l_pentagon) {

            ptr->nghbr_lst = calloc_2D_int(2, 5);

            /* IF (ptr0%l_pole_north    ) ix5(:) = (/ 1,2,3,4,5 /)
               IF (ptr0%l_pole_south    ) ix5(:) = (/ 1,2,3,4,5 /)
               IF (ptr0%l_pentagon_north) ix5(:) = (/ 1,2,3,5,6 /)
               IF (ptr0%l_pentagon_south) ix5(:) = (/ 1,2,3,4,5 /)
             */
            if (ptr->l_pentagon_north) {
                ix5[3]=5; ix5[4]=6;
            }

            for (n=0; n<5; n++) {
                ptr1 = ptr->nghbr[ix5[n]];
                /* neighbor block lists are block IDs, 0-based */
                ptr->nghbr_lst[0][n] = ptr1->tag_locl;
                ptr2 = twist(ptr1, ptr, 2);
                ptr->nghbr_lst[1][n] = ptr2->tag_locl;
            }
        } else {
            ptr->nghbr_lst = calloc_2D_int(2, 6);

            for (n=0; n<6; n++) {
                ptr1 = ptr->nghbr[ix6[n]];
                ptr->nghbr_lst[0][n] = ptr1->tag_locl;
                ptr2 = twist(ptr1, ptr, 2);
                ptr->nghbr_lst[1][n] = ptr2->tag_locl;
            }
        }
    }

    for (n=0; n<4; n++)
        if (ptr->dn[n] != NULL)
            set_nghbr_lst(ptr->dn[n]);
}


/*----< twist() >-------------------------------------------------------------*/
static
grid_node* twist(grid_node *ptr0,
                 grid_node *ptr1,
                 int        ntwist)
{
    int l_found, n, n1=0, ix5[5]={1,2,3,4,5}, ix6[6]={1,2,3,4,5,6};
    grid_node *ptr2;

    if (ptr0->l_pentagon) {
        /* IF (ptr0%l_pole_north    ) ix5(:) = (/ 1,2,3,4,5 /)
           IF (ptr0%l_pole_south    ) ix5(:) = (/ 1,2,3,4,5 /)
           IF (ptr0%l_pentagon_north) ix5(:) = (/ 1,2,3,5,6 /)
           IF (ptr0%l_pentagon_south) ix5(:) = (/ 1,2,3,4,5 /)
        */
        if (ptr0->l_pentagon_north) {
            ix5[3]=5; ix5[4]=6;
        }

        l_found = 0;
        for (n=0; n<5; n++) {
           if (ptr0->nghbr[ix5[n]]->tag_glbl == ptr1->tag_glbl) {
              l_found = 1;
              n1 = n;
              break;
           }
        }
        cshift_1D(ix5, 5, sizeof(int), ntwist);
        ptr2 = ptr0->nghbr[ix5[n1]];
    } else {
        l_found = 0;
        for (n=0; n<6; n++) {
           if (ptr0->nghbr[ix6[n]]->tag_glbl == ptr1->tag_glbl) {
              l_found = 1;
              n1 = n;
              break;
           }
        }
        cshift_1D(ix6, 6, sizeof(int), ntwist);
        ptr2 = ptr0->nghbr[ix6[n1]];
    }
    return ptr2;
}

/*----< spiral_path() >-------------------------------------------------------*/
grid_node* spiral_path(int        sprl_dpth,
                       grid_node *ptr0,
                       int        tag_optional)
{
    int l_found, tag_sprl, nsprl, n, position, sprl_lngth;
    grid_node *ptr0_sprl, *ptr1_sprl=NULL, *ptr2_sprl;
   
    int sprl_len[20] = {    6,  18,  36,  60,  90, 126, 168, 216, 270, 330,
                          396, 468, 546, 630, 720, 816, 918,1026,1140,1260};

    ptr0_sprl = ptr0;

    if (ptr0->nghbr[1] != NULL) {
        ptr0_sprl->next_sprl = ptr0->nghbr[1];
        ptr1_sprl            = ptr0_sprl->next_sprl;
    } else {
        printf(" spiral_path :: neighbor 1 not ASSOCIATED\n");
        printf(" spiral_path :: ptr0->tag_glbl = %d\n",ptr0->tag_glbl);
        printf(" spiral_path :: ptr0->level    = %d\n",ptr0->level);
        ABORT
    }

/* wkliao: temp solution assume -999 is an invalid value for tag_glbl */
    if (tag_optional != -999)
        tag_sprl = tag_optional;
    else
        tag_sprl = ptr0->tag_glbl;

    ptr0_sprl->tag_sprl = tag_sprl;
    ptr1_sprl->tag_sprl = tag_sprl;

    sprl_lngth = sprl_len[MIN(20, sprl_dpth)-1] - 1;
    if (ptr0->level == 0) sprl_lngth = MIN(sprl_len[0], sprl_lngth);
    if (ptr0->level == 1) sprl_lngth = MIN(sprl_len[1], sprl_lngth);
    if (ptr0->level == 2) sprl_lngth = MIN(sprl_len[2], sprl_lngth);
    if (ptr0->level == 3) sprl_lngth = MIN(sprl_len[3], sprl_lngth);

    for (nsprl=0; nsprl<sprl_lngth; nsprl++) {
   
        position = 0;
        for (n=1; n<7; n++) {
            position++;
            if (ptr1_sprl->nghbr[n] != NULL &&
                ptr1_sprl->nghbr[n]->tag_glbl == ptr0_sprl->tag_glbl)
                break;
        }

        l_found = 0;
        for (n=1; n<7; n++) {
            position = (position+4) % 6 + 1;
            if (ptr1_sprl->nghbr[position] != NULL) {
                ptr2_sprl = ptr1_sprl->nghbr[position];
                if (ptr2_sprl->tag_sprl != tag_sprl) {
                    l_found = 1;
                    break;
                }
            }
        }

        if (! l_found) {
            printf(" spiral_path :: (.NOT.l_found)\n");
            printf(" spiral_path :: ptr0->tag_glbl = %d\n",ptr0->tag_glbl);
            printf(" spiral_path :: ptr0->level    = %d\n",ptr0->level);
            printf(" spiral_path :: sprl_lngth     = %d\n",sprl_lngth);
            ABORT
        } else {
            ptr2_sprl->tag_sprl = tag_sprl;
            ptr1_sprl->next_sprl = ptr2_sprl;
            ptr0_sprl = ptr1_sprl;
            ptr1_sprl = ptr2_sprl;
        }
    }

    ptr2_sprl->next_sprl = NULL;

    return ptr0;
}

