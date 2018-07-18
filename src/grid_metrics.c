/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: grid_metrics.c 4604 2017-12-07 07:15:42Z wkliao $
 */

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* memset() */
#include <math.h>
#include <mpi.h>

#include "physical_params.h"
#include "grid_metrics.h"
#include "grid_subdomain.h"
#include "gcrm.h"
#include "grid_utilities.h"


static void read_grid_point(MODULE_grid_params *grid_params,
                            char *grid_point_select, char *grid_point_path);
static void initialize_grid_metrics_tree(grid_node *ptr);
static void initialize_grid_metrics_face(MODULE_grid_params *grid_params,
                            MODULE_grid_metrics *grid_metrics,
                            MODULE_wrap_data *wrap_data, char *wrap_name);
static void initialize_grid_metrics_corner(MODULE_grid_params *grid_params,
                            MODULE_wrap_data *wrap_data, char *wrap_name,
                            double  ****point, double *****point_crn,
                            double *****d_point_crn, double *****d_point_inv_crn,
                            double *****d_edge_crn, double *****d_edge_inv_crn,
                            double *****rlx_wght_crn);
static void initialize_wghts_crn(MODULE_grid_params  *grid_params,
                            MODULE_grid_subdomain *grid_subdomain,
                            double   ****point, double  *****point_crn,
                            int      ****l_msk_crn, double  *****wghts_crn,
                            double ******vctr_wghts_crn,
                            double ******vctr_wghts_crn2,
                            double  *****area_kite_crn);
static void initialize_grid_metrics_edge(MODULE_grid_params  *grid_params,
                            MODULE_wrap_data    *wrap_data, char *wrap_name,
                            double  ****point, double *****point_crn,
                            double *****point_edg, double  ****area_edg,
                            double *****nrm_edg, double *****tng_edg);
static void set_grid_point(char *grid_point_select, grid_node *ptr,
                            int level_max, int level);

/* MODULE grid_metrics */

/*----< init_MODULE_grid_metrics() >-----------------------------------------*/
/* Data malloc-ed
 *     grid_metrics->l_msk_fac
 *     grid_metrics->l_msk_crn
 *     grid_metrics->l_msk_edg
 *     grid_metrics->tag_glbl
 *     grid_metrics->tag_glbl_nghbr
 *     grid_metrics->point
 *     grid_metrics->corner
 *     grid_metrics->area
 *     grid_metrics->area_inv
 *     grid_metrics->d_point
 *     grid_metrics->d_point_inv
 *     grid_metrics->d_edge
 *     grid_metrics->d_edge_inv
 *     grid_metrics->point_crn
 *     grid_metrics->area_crn
 *     grid_metrics->area_inv_crn
 *     grid_metrics->area_kite_crn
 *     grid_metrics->d_point_crn
 *     grid_metrics->d_point_inv_crn
 *     grid_metrics->d_edge_crn
 *     grid_metrics->d_edge_inv_crn
 *     grid_metrics->rlx_wght_crn
 *     grid_metrics->point_edg
 *     grid_metrics->nrm_edg
 *     grid_metrics->tng_edg
 *     grid_metrics->area_edg
 *     grid_metrics->wghts_crn
 *     grid_metrics->vctr_wghts_crn
 *     grid_metrics->vctr_wghts_crn2
 *     grid_metrics->wghts_wind_crn
 *     grid_metrics->laplacian_wghts
 *     grid_metrics->laplacian_wghts_crn
 * Data updated:
 *     grid_metrics->grid_point_mode
 *     grid_metrics->grid_point_file
 */
void init_MODULE_grid_metrics(MODULE_grid_metrics *grid_metrics,
                              int                  im,
                              int                  jm,
                              int                  nsdm,
                              int                  edgm,
                              int                  crnm)
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  logical masks
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* logical mask for faces (cell centers) */
    grid_metrics->l_msk_fac = calloc_3D_int(nsdm, jm, im);
    /* logical mask for corners */
    grid_metrics->l_msk_crn = calloc_4D_int(nsdm, jm, im, crnm);
    /* logical mask for edges */
    grid_metrics->l_msk_edg = calloc_4D_int(nsdm, jm, im, edgm);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  metrics defined at cell centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    grid_metrics->tag_glbl       = calloc_3D_int(nsdm, jm, im);
    grid_metrics->tag_glbl_nghbr = calloc_4D_int(nsdm, jm, im, 6);
    grid_metrics->point          = calloc_4D_dbl(nsdm, jm, im, 3);
    grid_metrics->corner         = calloc_5D_dbl(nsdm, jm, im, 6, 3);
    grid_metrics->area           = calloc_3D_dbl(nsdm, jm, im);
    grid_metrics->area_inv       = calloc_3D_dbl(nsdm, jm, im);
    grid_metrics->d_point        = calloc_4D_dbl(nsdm, jm, im, 6);
    grid_metrics->d_point_inv    = calloc_4D_dbl(nsdm, jm, im, 6);
    grid_metrics->d_edge         = calloc_4D_dbl(nsdm, jm, im, 6);
    grid_metrics->d_edge_inv     = calloc_4D_dbl(nsdm, jm, im, 6);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  metrics defined at cell corners (abbreviated crn)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    grid_metrics->point_crn       = calloc_5D_dbl(nsdm, jm, im, crnm, 3);
    grid_metrics->area_crn        = calloc_4D_dbl(nsdm, jm, im, crnm);
    grid_metrics->area_inv_crn    = calloc_4D_dbl(nsdm, jm, im, crnm);
    grid_metrics->area_kite_crn   = calloc_5D_dbl(nsdm, jm, im, crnm, 3);
    grid_metrics->d_point_crn     = calloc_5D_dbl(nsdm, jm, im, crnm, 3);
    grid_metrics->d_point_inv_crn = calloc_5D_dbl(nsdm, jm, im, crnm, 3);
    grid_metrics->d_edge_crn      = calloc_5D_dbl(nsdm, jm, im, crnm, 3);
    grid_metrics->d_edge_inv_crn  = calloc_5D_dbl(nsdm, jm, im, crnm, 3);
    grid_metrics->rlx_wght_crn    = calloc_5D_dbl(nsdm, jm, im, crnm, 4);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  metrics defined at cell edges (abbreviated edg)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    grid_metrics->point_edg = calloc_5D_dbl(nsdm, jm, im, edgm, 3);
    grid_metrics->nrm_edg   = calloc_5D_dbl(nsdm, jm, im, edgm, 3);
    grid_metrics->tng_edg   = calloc_5D_dbl(nsdm, jm, im, edgm, 3);
    grid_metrics->area_edg  = calloc_4D_dbl(nsdm, jm, im, edgm);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  weights
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* weights to average cell centers to corners */
    grid_metrics->wghts_crn       = calloc_5D_dbl(nsdm, jm, im, crnm, 3);
    grid_metrics->vctr_wghts_crn  = calloc_6D_dbl(nsdm, jm, im, crnm, 2, 3);
    grid_metrics->vctr_wghts_crn2 = calloc_6D_dbl(nsdm, jm, im, crnm, 2, 3);
/*-----------------------------------------------------------------------
!  for wind at the corners in the GCRM
!-----------------------------------------------------------------------*/
    grid_metrics->wghts_wind_crn      = calloc_6D_dbl(nsdm, jm, im, crnm, 3, 3);
    grid_metrics->laplacian_wghts     = calloc_4D_dbl(nsdm, jm, im, 6);
    grid_metrics->laplacian_wghts_crn = calloc_5D_dbl(nsdm, jm, im, crnm, 3);

    grid_metrics->grid_point_mode = 2;
    grid_metrics->grid_point_file = 19;
}

/*----< finalize_MODULE_grid_metrics() >--------------------------------------*/
void finalize_MODULE_grid_metrics(MODULE_grid_metrics *grid_metrics)
{
    free_3D_int(grid_metrics->l_msk_fac);
    free_4D_int(grid_metrics->l_msk_crn);
    free_4D_int(grid_metrics->l_msk_edg);
    free_3D_int(grid_metrics->tag_glbl);
    free_4D_int(grid_metrics->tag_glbl_nghbr);
    free_4D_dbl(grid_metrics->point);
    free_5D_dbl(grid_metrics->corner);
    free_3D_dbl(grid_metrics->area);
    free_3D_dbl(grid_metrics->area_inv);
    free_4D_dbl(grid_metrics->d_point);
    free_4D_dbl(grid_metrics->d_point_inv);
    free_4D_dbl(grid_metrics->d_edge);
    free_4D_dbl(grid_metrics->d_edge_inv);
    free_5D_dbl(grid_metrics->point_crn);
    free_4D_dbl(grid_metrics->area_crn);
    free_4D_dbl(grid_metrics->area_inv_crn);
    free_5D_dbl(grid_metrics->area_kite_crn);
    free_5D_dbl(grid_metrics->d_point_crn);
    free_5D_dbl(grid_metrics->d_point_inv_crn);
    free_5D_dbl(grid_metrics->d_edge_crn);
    free_5D_dbl(grid_metrics->d_edge_inv_crn);
    free_5D_dbl(grid_metrics->rlx_wght_crn);
    free_5D_dbl(grid_metrics->point_edg);
    free_5D_dbl(grid_metrics->nrm_edg);
    free_5D_dbl(grid_metrics->tng_edg);
    free_4D_dbl(grid_metrics->area_edg);
    free_5D_dbl(grid_metrics->wghts_crn);
    free_6D_dbl(grid_metrics->vctr_wghts_crn);
    free_6D_dbl(grid_metrics->vctr_wghts_crn2);
    free_6D_dbl(grid_metrics->wghts_wind_crn);
    free_4D_dbl(grid_metrics->laplacian_wghts);
    free_5D_dbl(grid_metrics->laplacian_wghts_crn);
}

/*----< initialize_grid_metrics() >-------------------------------------------*/
void initialize_grid_metrics(MODULE_grid_params    *grid_params,
                             MODULE_grid_subdomain *grid_subdomain,
                             MODULE_grid_metrics   *grid_metrics,
                             MODULE_wrap_data      *wrap_data,
                             char                  *wrap_name,
                             char                  *grid_point_select,
                             char                  *grid_point_path_tmpry)
{
    int n, i, j, k, nsd, n1, n2, ix[3];
    char grid_point_path[128];
    grid_node *ptr, *grid;

    int im        = grid_params->im,
        jm        = grid_params->jm,
        nsdm      = grid_params->nsdm,
        level_max = grid_params->level_max;

    grid = grid_params->grid;

/*-----------------------------------------------------------------------
!  read the positions of the grid points
!-----------------------------------------------------------------------*/
    read_grid_point(grid_params, grid_point_select, grid_point_path);

/*-----------------------------------------------------------------------
!  apply a rotation to all grid points
!-----------------------------------------------------------------------*/
/*
    if (0) {
        for (n=0; n<12; n++)
            rotate_grid(grid[n].nghbr[0].p);
    }
*/
/*-----------------------------------------------------------------------
!  set grid metrics within the tree structure
!-----------------------------------------------------------------------*/
    for (n=0; n<12; n++)  /* grid[12] */
        initialize_grid_metrics_tree(grid[n].nghbr[0]);
/*-----------------------------------------------------------------------
!  set grid metrics for local subdomains
!-----------------------------------------------------------------------*/
    initialize_grid_metrics_face(grid_params, grid_metrics, wrap_data, wrap_name);
/*-----------------------------------------------------------------------
!  set logical masks
!-----------------------------------------------------------------------*/
    for (k=0; k<nsdm; k++) {
        for (j=0; j<jm; j++) {
            for (i=0; i<im; i++)
                grid_metrics->l_msk_fac[k][j][i] = 0;
        }
        for (j=1; j<jm-1; j++) {
            for (i=1; i<im-1; i++)
                grid_metrics->l_msk_fac[k][j][i] = 1;
        }
    }
    if (grid_params->l_agent_north) /* nsd_north is block ID, 0-based */
        grid_metrics->l_msk_fac[grid_params->nsd_north][jm-1][1] = 1;
    if (grid_params->l_agent_south) /* nsd_south is block ID, 0-based */
        grid_metrics->l_msk_fac[grid_params->nsd_south][1][im-1] = 1;
/*-----------------------------------------------------------------------
!  set metrics associated with the corner
!-----------------------------------------------------------------------*/
    initialize_grid_metrics_corner(grid_params,
                                   wrap_data,
                                   wrap_name,
                                   grid_metrics->point,
                                   grid_metrics->point_crn,
                                   grid_metrics->d_point_crn,
                                   grid_metrics->d_point_inv_crn,
                                   grid_metrics->d_edge_crn,
                                   grid_metrics->d_edge_inv_crn,
                                   grid_metrics->rlx_wght_crn);

/*-----------------------------------------------------------------------
!  set logical masks and weights associated with corners
!-----------------------------------------------------------------------*/
    initialize_wghts_crn(grid_params,
                         grid_subdomain,
                         grid_metrics->point,
                         grid_metrics->point_crn,
                         grid_metrics->l_msk_crn,
                         grid_metrics->wghts_crn,
                         grid_metrics->vctr_wghts_crn,
                         grid_metrics->vctr_wghts_crn2,
                         grid_metrics->area_kite_crn);
/*-----------------------------------------------------------------------
!  set the EDGE grid point, normal and tangent directions
!-----------------------------------------------------------------------*/
    initialize_grid_metrics_edge(grid_params,
                                 wrap_data,
                                 wrap_name,
                                 grid_metrics->point,
                                 grid_metrics->point_crn,
                                 grid_metrics->point_edg,
                                 grid_metrics->area_edg,
                                 grid_metrics->nrm_edg,
                                 grid_metrics->tng_edg);
/*-----------------------------------------------------------------------
!  write the grid to a file
!-----------------------------------------------------------------------*/
/*
        if (rnk_wrld == 0) printf(" grid_metrics :: write_grid \n");
         write_grid();
*/
/*-----------------------------------------------------------------------
!  set neighbor list
!-----------------------------------------------------------------------*/

    ptr = grid_params->path_real[level_max];

/*-----------------------------------------------------------------------
!  loop over all real grid points (not ghost) of the local process
!-----------------------------------------------------------------------*/
    for (n1=0; n1<grid_params->nm_real; n1++) {
/*-----------------------------------------------------------------------
!  loop over neighboring directions.  the pentagons will have a direction
!  where the pointer is not associated.
!-----------------------------------------------------------------------*/
        get_index(grid, level_max, level_max, grid_params->sbdmn_iota,
                  grid_params->sbdmn[level_max].lst,
                  grid_params->sbdmn[level_max].nsdm, ptr->tag_glbl, ix);
                  /* valuses in ix[] are 1-based */
        i   = ix[0] - 1;
        j   = ix[1] - 1;
        nsd = ix[2] - 1;

        for (n2=1; n2<=6; n2++) {
            n = 0;
            if (ptr->nghbr[n2] != NULL)
                grid_metrics->tag_glbl_nghbr[nsd][j][i][n2-1] =
                     ptr->nghbr[n2]->tag_glbl;
                /* global tag of the neighboring cell */
            else
                n = n2;

            if (n > 0) {
                if (n > 1)
                    grid_metrics->tag_glbl_nghbr[nsd][j][i][n-1] =
                    grid_metrics->tag_glbl_nghbr[nsd][j][i][n-2];
                else
                    grid_metrics->tag_glbl_nghbr[nsd][j][i][0] =
                    grid_metrics->tag_glbl_nghbr[nsd][j][i][5];
            }
        }
        ptr = ptr->next_real;
    }
}

/*----< read_grid_point() >---------------------------------------------------*/
static
void read_grid_point(MODULE_grid_params *grid_params,
                     char               *grid_point_select,
                     char               *grid_point_path)
{
    int n, level;
    grid_node *grid = grid_params->grid;

/*-----------------------------------------------------------------------
!  set grid points
!-----------------------------------------------------------------------*/
    if (!strcmp(grid_point_select, "bisect")) {
        for (n=0; n<12; n++)  /* grid[12] */
            grid[n].l_point_set = 1;

        grid[0].point[0] = 0;
        grid[0].point[1] = 0;
        grid[0].point[2] = 1;
        grid[1].point[0] = 0;
        grid[1].point[1] = 0;
        grid[1].point[2] = -1;

        for (n=0; n<=8; n+=2) {
            double lonlat[2];
            lonlat[0] = grid_params->alfalfa * (-0.5 + n/2);
            lonlat[1] = grid_params->pentagon_latitude;
            lonlat_to_xyz(lonlat, grid[n+2].point);

            lonlat[0] = grid_params->alfalfa * (n/2);
            lonlat[1] = grid_params->pentagon_latitude * -1.0;
            lonlat_to_xyz(lonlat, grid[n+3].point);
        }

        for (level=0; level<grid_params->level_max; level++)
            for (n=0; n<12; n++)  /* grid[12] */
                set_grid_point(grid_point_select, grid[n].nghbr[0],
                               grid_params->level_max, level);
    }
    else {
        printf("initialize_grid_metrics :: cannot determine grid_point_select\n");
        ABORT
    }
}

/*----< initialize_grid_metrics_tree() >-------------------------------------*/
static
void initialize_grid_metrics_tree(grid_node *ptr)
{
    int l_lst[6], l_bad, n, m;
    grid_node *ptr_a, *ptr_b;
    l_bad = 0;

    if (! ptr->nghbr[0]->l_point_set) l_bad = 1;

    for (m=1; m<7; m++) {  /* Fortran nghbr(0:6) */
        if (ptr->nghbr[m] != NULL) {
            l_lst[m-1] = 1;
            if (! ptr->nghbr[m]->l_point_set) l_bad = 1;
        }
        else
            l_lst[m-1] = 0;
    }

    int count = 0;
    for (m=0; m<6; m++)
        if (l_lst[m] == 1)
            count++;

    if (l_bad || count < 5) {
        for (m=0; m<6; m++) /* double corner[6][3]; */
            for (n=0; n<3; n++)
                ptr->corner[m][n] = 0.0;
        ptr->area     = 0.0;
        ptr->area_inv = 0.0;
        for (m=0; m<2; m++) {
            ptr->area_crn[m]     = 0.0;  /* double area_crn[2]; */
            ptr->area_inv_crn[m] = 0.0;  /* double area_inv_crn[2]; */
        }
    }
    else {
/*-----------------------------------------------------------------------
!  set cell corners
!-----------------------------------------------------------------------*/
        for (m=1; m<=6; m++) {  /* Fortran nghbr(0:6) */
            if (ptr->nghbr[(m+5)%6+1] != NULL)
                ptr_a = ptr->nghbr[(m+5)%6+1];
            else
                ptr_a = ptr->nghbr[(m+4)%6+1];

            if (ptr->nghbr[m%6+1] != NULL)
                ptr_b = ptr->nghbr[m%6+1];
            else
                ptr_b = ptr->nghbr[(m+1)%6+1];

            voronoi_corner(ptr->point, ptr_a->point, ptr_b->point,
                           ptr->corner[m-1]);
        }
/*-----------------------------------------------------------------------
!  set FACE areas and the inverse of the areas
!-----------------------------------------------------------------------*/
        if (ptr->l_pentagon) {
            ptr->area = 5.0 * spherical_triangle_area(ptr->point,
                                                      ptr->corner[0],
                                                      ptr->corner[1]);
        }
        else {
            ptr->area = 0.0;
            for (m=1; m<=6; m++)
                ptr->area += spherical_triangle_area(ptr->point,
                                                     ptr->corner[m-1],
                                                     ptr->corner[m%6]);
        }
        ptr->area     *= A * A;
        ptr->area_inv  = 1.0 / ptr->area;
/*-----------------------------------------------------------------------
!  set CORNER areas and the inverse of the CORNER areas
!-----------------------------------------------------------------------*/
        ptr->area_crn[0]     = 0.0;
        ptr->area_crn[1]     = 0.0;
        ptr->area_inv_crn[0] = 0.0;
        ptr->area_inv_crn[1] = 0.0;

        if (ptr->l_pole_north || ptr->l_pole_south) {
            ptr->area_crn[0]     = 0.0;
            ptr->area_crn[1]     = 0.0;
            ptr->area_inv_crn[0] = 0.0;
            ptr->area_inv_crn[1] = 0.0;
        }
        else {
            for (m=1; m<=2; m++) {
                if (ptr->nghbr[m  ] == NULL) break;
                if (ptr->nghbr[m+1] == NULL) break;
                ptr_a = ptr->nghbr[m];
                ptr_b = ptr->nghbr[m+1];
                ptr->area_crn[m-1] = spherical_triangle_area(ptr->point,
                                                             ptr_a->point,
                                                             ptr_b->point);
                ptr->area_crn[m-1]     *= A*A;
                ptr->area_inv_crn[m-1]  = 1.0/ptr->area_crn[m-1];
            }
        }
    }

    for (n=0; n<4; n++) {
        if (ptr->dn[n] != NULL)
            initialize_grid_metrics_tree(ptr->dn[n]);
    }
}

/*----< initialize_grid_metrics_face() >--------------------------------------*/
static
void initialize_grid_metrics_face(MODULE_grid_params  *grid_params,
                                  MODULE_grid_metrics *grid_metrics,
                                  MODULE_wrap_data    *wrap_data,
                                  char                *wrap_name)
{
    int offset, sbdmn_north[5], sbdmn_south[5], nsd, tag_nsd;
    int i, j, k, l, m, n;
    grid_node *ptr0, *ptr1, *grid;

    int im         = grid_params->im,
        jm         = grid_params->jm,
        nsdm       = grid_params->nsdm,
        crnm       = grid_params->crnm,
        sbdmn_iota = grid_params->sbdmn_iota,
        level_max  = grid_params->level_max;

    sbdmn_node *sbdmn = grid_params->sbdmn;
    grid = grid_params->grid;

    offset = ((POWER2(sbdmn_iota)+1)*(POWER2(sbdmn_iota)-1))/3;

    for (m=0; m<5; m++) {
        sbdmn_north[m] = POWER2(2*sbdmn_iota) * (2*m+1) - 1;
        sbdmn_south[m] = POWER2(2*sbdmn_iota) * (2*m+1) + offset;
    }

    for (nsd=0; nsd<nsdm; nsd++) {

        /* sbdmn[].lst[*] are block IDs, 0-based, tags are 1-based */
        tag_nsd = 3+(POWER2(2*(level_max-sbdmn_iota)))*(sbdmn[level_max].lst[nsd]);
        ptr0 = set_ptr(grid, level_max, level_max, tag_nsd);

        for (j=1; j<jm-1; j++) {
            ptr1 = ptr0->nghbr[0];
            for (i=1; i<im-1; i++) {
                grid_metrics->tag_glbl[nsd][j][i] = ptr1->tag_glbl;
                grid_metrics->area    [nsd][j][i] = ptr1->area;
                grid_metrics->area_inv[nsd][j][i] = ptr1->area_inv;

                for (k=0; k<3; k++)
                    grid_metrics->point[nsd][j][i][k] = ptr1->point[k];

                /* grid_metrics->corner[nsdm][jm][im][6][3] */
                for (l=0; l<6; l++)
                    for (k=0; k<3; k++)
                        grid_metrics->corner[nsd][j][i][l][k] = ptr1->corner[l][k];

                for (k=0; k<2; k++) {
                    grid_metrics->area_crn    [nsd][j][i][k] = ptr1->area_crn[k];
                    grid_metrics->area_inv_crn[nsd][j][i][k] = ptr1->area_inv_crn[k];
                }
/*-----------------------------------------------------------------------
!  distances and lengths associated with the hexagonal grid
!-----------------------------------------------------------------------*/
                for (n=0; n<6; n++) {  /* nghbr(0:6) */
                    if (ptr1->nghbr[n+1] != NULL) {
                        grid_metrics->d_point[nsd][j][i][n] =
                             arch_distance(ptr1->point,
                                           ptr1->nghbr[n+1]->point);
                        grid_metrics->d_edge[nsd][j][i][n] =
                             arch_distance(grid_metrics->corner[nsd][j][i][(n+5)%6],
                                           grid_metrics->corner[nsd][j][i][n]);

                        grid_metrics->d_point_inv[nsd][j][i][n] = 1.0 /
                        grid_metrics->d_point    [nsd][j][i][n];
                        grid_metrics->d_edge_inv [nsd][j][i][n] = 1.0 /
                        grid_metrics->d_edge     [nsd][j][i][n];

                        grid_metrics->laplacian_wghts[nsd][j][i][n] =
                        grid_metrics->d_edge         [nsd][j][i][n] /
                        grid_metrics->d_point        [nsd][j][i][n];
                    }
                }
                ptr1 = ptr1->nghbr[1];
            }
            ptr0 = ptr0->nghbr[3];
        }
/*-----------------------------------------------------------------------
!  north pole
!-----------------------------------------------------------------------*/
        for (j=0; j<5; j++)  /* Fortran sbdmn(0:level_max) */
            if (sbdmn_north[j] == sbdmn[level_max].lst[nsd])
                /* sbdmn[].lst[*] are block IDs, 0-based */
                break;

        if (j < 5) {
            ptr0 = set_ptr(grid, level_max, level_max, 1);

            grid_metrics->tag_glbl[nsd][jm-1][1] = ptr0->tag_glbl;
            grid_metrics->area    [nsd][jm-1][1] = ptr0->area;
            grid_metrics->area_inv[nsd][jm-1][1] = ptr0->area_inv;

            for (k=0; k<3; k++)
                grid_metrics->point[nsd][jm-1][1][k] = ptr0->point[k];

            /* grid_metrics->corner[nsdm][jm][im][6][3] */
            for (l=0; l<6; l++)
                for (k=0; k<3; k++)
                    grid_metrics->corner[nsd][jm-1][1][l][k] = ptr0->corner[l][k];

            for (k=0; k<5; k++) {
                grid_metrics->d_point[nsd][jm-1][1][k] =
                     arch_distance(ptr0->point, ptr0->nghbr[1]->point);

                grid_metrics->d_edge[nsd][jm-1][1][k] =
                     arch_distance(grid_metrics->corner[nsd][jm-1][1][0],
                                   grid_metrics->corner[nsd][jm-1][1][1]);

                grid_metrics->d_point_inv[nsd][jm-1][1][k] = 1.0 /
                    grid_metrics->d_point[nsd][jm-1][1][k];

                grid_metrics->d_edge_inv [nsd][jm-1][1][k] = 1.0 /
                    grid_metrics->d_edge [nsd][jm-1][1][k];

                grid_metrics->laplacian_wghts[nsd][jm-1][1][k] =
                        grid_metrics->d_edge [nsd][jm-1][1][k] /
                        grid_metrics->d_point[nsd][jm-1][1][k];
            }
        }
/*-----------------------------------------------------------------------
!  south pole
!-----------------------------------------------------------------------*/
        for (j=0; j<5; j++)
            if (sbdmn_south[j] == sbdmn[level_max].lst[nsd])
                /* sbdmn[].lst[*] are block IDs, 0-based */
                break;

        if (j < 5) {
            ptr0 = set_ptr(grid, level_max, level_max, 2);

            grid_metrics->tag_glbl[nsd][1][im-1] = ptr0->tag_glbl;
            grid_metrics->area    [nsd][1][im-1] = ptr0->area;
            grid_metrics->area_inv[nsd][1][im-1] = ptr0->area_inv;

            for (k=0; k<3; k++)
                grid_metrics->point[nsd][1][im-1][k] = ptr0->point[k];

            for (l=0; l<6; l++)
                for (k=0; k<3; k++)
                    grid_metrics->corner[nsd][1][im-1][l][k] = ptr0->corner[l][k];

            for (k=0; k<5; k++) {
                grid_metrics->d_point[nsd][1][im-1][k] =
                     arch_distance(ptr0->point, ptr0->nghbr[1]->point);

                grid_metrics->d_edge[nsd][1][im-1][k] =
                     arch_distance(grid_metrics->corner[nsd][1][im-1][0],
                                   grid_metrics->corner[nsd][1][im-1][1]);

                grid_metrics->d_point_inv[nsd][1][im-1][k] = 1.0 /
                    grid_metrics->d_point[nsd][1][im-1][k];

                grid_metrics->d_edge_inv [nsd][1][im-1][k] = 1.0 /
                    grid_metrics->d_edge [nsd][1][im-1][k];

                grid_metrics->laplacian_wghts[nsd][1][im-1][k] =
                        grid_metrics->d_edge [nsd][1][im-1][k] /
                        grid_metrics->d_point[nsd][1][im-1][k];
            }
        }
    }
/*-----------------------------------------------------------------------
!  wrap local metrics
!-----------------------------------------------------------------------*/
    int dimlen[6];

    // point = calloc_4D_dbl(nsdm, jm, im, 3);
    dimlen[0] = nsdm;
    dimlen[1] = jm;
    dimlen[2] = im;
    dimlen[3] = 3;
    for (m=0; m<3; m++)
        wrap_face_1lyr(wrap_data, wrap_name, grid_metrics->point, dimlen, m);

    // corner = calloc_5D_dbl(nsdm, jm, im, 6, 3);
    dimlen[3] = 2; /* corner(3,6,im,jm,nsdm) but using (:,1:2,:,:,:) */
    dimlen[4] = 3;
    wrap_vrtx_1lyr(wrap_data, wrap_name, grid_metrics->corner, dimlen);

    // area = calloc_3D_dbl(nsdm, jm, im);
    wrap_face_1lyr(wrap_data, wrap_name, &grid_metrics->area, dimlen, -1);

    // area_inv = calloc_3D_dbl(nsdm, jm, im);
    wrap_face_1lyr(wrap_data, wrap_name, &grid_metrics->area_inv, dimlen, -1);

    // area_crn = calloc_4D_dbl(nsdm, jm, im, crnm);
    dimlen[3] = crnm;
    wrap_vrtx_scalar_1lyr(wrap_data, wrap_name, grid_metrics->area_crn, dimlen);

    // area_inv_crn = calloc_4D_dbl(nsdm, jm, im, crnm);
    wrap_vrtx_scalar_1lyr(wrap_data, wrap_name, grid_metrics->area_inv_crn, dimlen);

    // d_edge = calloc_4D_dbl(nsdm, jm, im, 6);
    dimlen[3] = 3;   /* using only (1:3,:,:,:) */
    wrap_edge_scalar_1lyr(wrap_data, wrap_name, grid_metrics->d_edge, dimlen);

    // d_edge_inv = calloc_4D_dbl(nsdm, jm, im, 6);
    wrap_edge_scalar_1lyr(wrap_data, wrap_name, grid_metrics->d_edge_inv, dimlen);

/*-----------------------------------------------------------------------
!  get laplacian weights for corners
!-----------------------------------------------------------------------*/
    for (nsd=0; nsd<nsdm; nsd++) {
        for (j=1; j<jm-1; j++) {
            for (i=1; i<im-1; i++) {
                grid_metrics->laplacian_wghts_crn[nsd][j][i][0][0] = 1.0 /
                     grid_metrics->laplacian_wghts[nsd][j][i][0];

                grid_metrics->laplacian_wghts_crn[nsd][j][i][0][1] =
                     arch_distance(grid_metrics->point[nsd][j][i+1],
                                   grid_metrics->point[nsd][j+1][i+1]) /
                     arch_distance(grid_metrics->corner[nsd][j][i][0],
                                   grid_metrics->corner[nsd][j][i+1][1]);

                grid_metrics->laplacian_wghts_crn[nsd][j][i][0][2] = 1.0 /
                     grid_metrics->laplacian_wghts[nsd][j][i][1];
                grid_metrics->laplacian_wghts_crn[nsd][j][i][1][0] =
                     grid_metrics->laplacian_wghts_crn[nsd][j][i][0][2];
                grid_metrics->laplacian_wghts_crn[nsd][j][i][1][1] =
                     arch_distance(grid_metrics->point[nsd][j+1][i+1],
                                   grid_metrics->point[nsd][j+1][i]) /
                     arch_distance(grid_metrics->corner[nsd][j][i][1],
                                   grid_metrics->corner[nsd][j+1][i][0]);
                grid_metrics->laplacian_wghts_crn[nsd][j][i][1][2] = 1.0 /
                     grid_metrics->laplacian_wghts[nsd][j][i][2];
            }
        }
    }
/*-----------------------------------------------------------------------
!  set lengths on the earth
!-----------------------------------------------------------------------*/
    for (nsd=0; nsd<nsdm; nsd++) {
        for (j=0; j<jm; j++) {
            for (i=0; i<im; i++) {
                for (k=0; k<6; k++) {
                    grid_metrics->d_point    [nsd][j][i][k] *= A;
                    grid_metrics->d_point_inv[nsd][j][i][k] *= 1.0/A;
                    grid_metrics->d_edge     [nsd][j][i][k] *= A;
                    grid_metrics->d_edge_inv [nsd][j][i][k] *= 1.0/A;
                }
            }
        }
    }
}

/*----< initialize_grid_metrics_corner() >------------------------------------*/
static
void initialize_grid_metrics_corner(MODULE_grid_params  *grid_params,
                                    MODULE_wrap_data    *wrap_data,
                                    char       *wrap_name,
                                    double  ****point,           /* [nsdm][jm][im][3] */
                                    double *****point_crn,       /* [nsdm][jm][im][crnm][3] */
                                    double *****d_point_crn,     /* [nsdm][jm][im][crnm][3] */
                                    double *****d_point_inv_crn, /* [nsdm][jm][im][crnm][3] */
                                    double *****d_edge_crn,      /* [nsdm][jm][im][crnm][3] */
                                    double *****d_edge_inv_crn,  /* [nsdm][jm][im][crnm][3] */
                                    double *****rlx_wght_crn)    /* [nsdm][jm][im][crnm][4] */
{
    int i, j, k, l, nsd, dimlen[6];

    int im         = grid_params->im,
        jm         = grid_params->jm,
        nsdm       = grid_params->nsdm,
        crnm       = grid_params->crnm;
/*-----------------------------------------------------------------------
!  set point_crn
!-----------------------------------------------------------------------*/
    for (nsd=0; nsd<nsdm; nsd++) {
        for (j=1; j<jm-1; j++) {
            for (i=1; i<im-1; i++) {
                voronoi_corner(point[nsd][j  ][i  ],
                               point[nsd][j  ][i+1],
                               point[nsd][j+1][i+1],
                               point_crn[nsd][j][i][0]);
                voronoi_corner(point[nsd][j  ][i  ],
                               point[nsd][j+1][i+1],
                               point[nsd][j+1][i  ],
                               point_crn[nsd][j][i][1]);
            }
        }
    }
    dimlen[0] = nsdm;
    dimlen[1] = jm;
    dimlen[2] = im;
    dimlen[3] = crnm;
    dimlen[4] = 3;
    wrap_vrtx_1lyr(wrap_data, wrap_name, point_crn, dimlen);
/*-----------------------------------------------------------------------
!  set d_point_crn
!-----------------------------------------------------------------------*/
    for (nsd=0; nsd<nsdm; nsd++) {
        for (j=1; j<jm-1; j++) {
            for (i=1; i<im-1; i++) {
                d_point_crn[nsd][j][i][0][0] =
                    arch_distance(point_crn[nsd][j][i][0],
                                  point_crn[nsd][j-1][i][2]);
                d_point_crn[nsd][j][i][0][1] =
                    arch_distance(point_crn[nsd][j][i][0],
                                  point_crn[nsd][j][i+1][1]);
                d_point_crn[nsd][j][i][0][2] =
                    arch_distance(point_crn[nsd][j][i][0],
                                  point_crn[nsd][j][i][1]);

                d_point_crn[nsd][j][i][1][0] =
                    arch_distance(point_crn[nsd][j][i][1],
                                  point_crn[nsd][j][i][0]);
                d_point_crn[nsd][j][i][1][1] =
                    arch_distance(point_crn[nsd][j][i][1],
                                  point_crn[nsd][j+1][i][0]);
                d_point_crn[nsd][j][i][1][2] =
                    arch_distance(point_crn[nsd][j][i][1],
                                  point_crn[nsd][j][i-1][0]);
            }
        }
    }
/*-----------------------------------------------------------------------
!  set d_edge_crn
!-----------------------------------------------------------------------*/
    for (nsd=0; nsd<nsdm; nsd++) {
        for (j=1; j<jm-1; j++) {
            for (i=1; i<im-1; i++) {
                d_edge_crn[nsd][j][i][0][0] =
                    arch_distance(point[nsd][j][i],
                                  point[nsd][j][i+1]);
                d_edge_crn[nsd][j][i][0][1] =
                    arch_distance(point[nsd][j][i+1],
                                  point[nsd][j+1][i+1]);
                d_edge_crn[nsd][j][i][0][2] =
                    arch_distance(point[nsd][j+1][i+1],
                                  point[nsd][j][i]);

                d_edge_crn[nsd][j][i][1][0] =
                    arch_distance(point[nsd][j][i],
                                  point[nsd][j+1][i+1]);
                d_edge_crn[nsd][j][i][1][1] =
                    arch_distance(point[nsd][j+1][i+1],
                                  point[nsd][j+1][i]);
                d_edge_crn[nsd][j][i][1][2] =
                    arch_distance(point[nsd][j+1][i],
                                  point[nsd][j][i]);
            }
        }
    }
/*-----------------------------------------------------------------------
!  set rlx_wght_crn
!-----------------------------------------------------------------------*/
    for (nsd=0; nsd<nsdm; nsd++) {
        for (j=1; j<jm-1; j++) {
            for (i=1; i<im-1; i++) {  /* rlx_wght_crn's last dim is (0:3) */
                rlx_wght_crn[nsd][j][i][0][1] =  d_edge_crn[nsd][j][i][0][0]
                                              / d_point_crn[nsd][j][i][0][0];
                rlx_wght_crn[nsd][j][i][0][2] =  d_edge_crn[nsd][j][i][0][1]
                                              / d_point_crn[nsd][j][i][0][1];
                rlx_wght_crn[nsd][j][i][0][3] =  d_edge_crn[nsd][j][i][0][2]
                                              / d_point_crn[nsd][j][i][0][2];

                double sum = 0.0;
                for (k=0; k<4; k++) sum += rlx_wght_crn[nsd][j][i][0][k];
                rlx_wght_crn[nsd][j][i][0][0] = 1.0/sum;

                rlx_wght_crn[nsd][j][i][1][1] =  d_edge_crn[nsd][j][i][1][0]
                                              / d_point_crn[nsd][j][i][1][0];
                rlx_wght_crn[nsd][j][i][1][2] =  d_edge_crn[nsd][j][i][1][1]
                                              / d_point_crn[nsd][j][i][1][1];
                rlx_wght_crn[nsd][j][i][1][3] =  d_edge_crn[nsd][j][i][1][2]
                                              / d_point_crn[nsd][j][i][1][2];

                sum = 0.0;
                for (k=0; k<4; k++) sum += rlx_wght_crn[nsd][j][i][1][k];
                rlx_wght_crn[nsd][j][i][1][0] = 1.0 / sum;
            }
        }
    }
/*-----------------------------------------------------------------------
!  make as big as the earth
!-----------------------------------------------------------------------*/
    for (nsd=0; nsd<nsdm; nsd++)
        for (j=1; j<jm-1; j++)
            for (i=1; i<im-1; i++)
                for (k=0; k<crnm; k++)
                    for (l=0; l<3; l++) {
                        d_point_crn[nsd][j][i][k][l] *= A;
                        d_edge_crn [nsd][j][i][k][l] *= A;
                    }
/*-----------------------------------------------------------------------
!  form the inverse
!-----------------------------------------------------------------------*/
    for (nsd=0; nsd<nsdm; nsd++)
        for (j=1; j<jm-1; j++)
            for (i=1; i<im-1; i++)
                for (k=0; k<crnm; k++)
                    for (l=0; l<3; l++) {
                        if (d_point_crn[nsd][j][i][k][l] != 0.0)
                            d_point_inv_crn[nsd][j][i][k][l] = 1.0 /
                                d_point_crn[nsd][j][i][k][l];
                        if (d_edge_crn[nsd][j][i][k][l] != 0.0)
                            d_edge_inv_crn[nsd][j][i][k][l] = 1.0 /
                                d_edge_crn[nsd][j][i][k][l];
                    }
}

/*----< initialize_wghts_crn() >---------------------------------------------*/
static
void initialize_wghts_crn(MODULE_grid_params  *grid_params,
                          MODULE_grid_subdomain *grid_subdomain,
                          double   ****point,           /* [nsdm][jm][im][3] */
                          double  *****point_crn,       /* [nsdm][jm][im][crnm][3] */
                          int      ****l_msk_crn,       /* [nsdm][jm][im][crnm] */
                          double  *****wghts_crn,       /* [nsdm][jm][im][crnm][3] */
                          double ******vctr_wghts_crn,  /* [nsdm][jm][im][crnm][2][3] */
                          double ******vctr_wghts_crn2, /* [nsdm][jm][im][crnm][2][3] */
                          double  *****area_kite_crn)   /* [nsdm][jm][im][crnm][3] */
{
    int l_crn[2], n, i, j, k, l, m, nsd;   /* crnm == 2 */
    double p[4][3];  /* Fortran (3,0:3) */

    int im         = grid_params->im,
        jm         = grid_params->jm,
        nsdm       = grid_params->nsdm,
        crnm       = grid_params->crnm;

    for (nsd=0; nsd<nsdm; nsd++) {
        for (j=0; j<jm-1; j++) {
            for (i=0; i<im-1; i++) {
                for (k=0; k<3; k++) {
                    p[0][k] = point[nsd][j  ][i  ][k];
                    p[1][k] = point[nsd][j  ][i+1][k];
                    p[2][k] = point[nsd][j+1][i+1][k];
                    p[3][k] = point[nsd][j+1][i  ][k];
                }

                l_crn[0] = l_crn[1] = 1;

                if (i == 0 && j == 0) {
                    if (grid_subdomain->l_sbdmn_pntgn_south[nsd]) l_crn[0] = 0;
                    if (grid_subdomain->l_sbdmn_pntgn_north[nsd]) l_crn[1] = 0;
                }
/*-----------------------------------------------------------------------
!  logical mask
!-----------------------------------------------------------------------*/
                l_msk_crn[nsd][j][i][0] = l_crn[0];
                l_msk_crn[nsd][j][i][1] = l_crn[1];
/*-----------------------------------------------------------------------
!  weights to area-weight interpolate from cell centers to cell corners
!-----------------------------------------------------------------------*/
                for (n=0; n<2; n++) {
                    if (l_crn[n]) {
                        area_corner_kites(p[0], p[n+1], p[n+2],
                                          area_kite_crn[nsd][j][i][n]);
                        for (k=0; k<3; k++)
                            wghts_crn[nsd][j][i][n][k] =
                                 area_kite_crn[nsd][j][i][n][k] /
                                 spherical_triangle_area(p[0], p[n+1], p[n+2]);
                    }
                }
/*-----------------------------------------------------------------------
!  weights for gradient defined at corners
!-----------------------------------------------------------------------*/
                for (n=0; n<2; n++)
                    if (l_crn[n])
                        set_vctr_wghts_crn(p[0], p[n+1], p[n+2],
                                           vctr_wghts_crn[nsd][j][i][n]);
            } /* loop i */
        } /* loop j */
/*-----------------------------------------------------------------------
!  weights for gradient defined at corners with corner inputs
!-----------------------------------------------------------------------*/
        for (j=1; j<jm-1; j++) {
            for (i=1; i<im-1; i++) {
                set_vctr_wghts_crn(point_crn[nsd][j-1][i-1][1],
                                   point_crn[nsd][j  ][i+1][1],
                                   point_crn[nsd][j  ][i  ][1],
                                   vctr_wghts_crn2[nsd][j][i][0]);
                set_vctr_wghts_crn(point_crn[nsd][j  ][i  ][0],
                                   point_crn[nsd][j+1][i  ][0],
                                   point_crn[nsd][j  ][i-1][0],
                                   vctr_wghts_crn2[nsd][j][i][1]);
            }
        }
/*-----------------------------------------------------------------------
!  north pole
!-----------------------------------------------------------------------*/
        if (grid_subdomain->l_sbdmn_north_pole[nsd]) {

            for (k=0; k<3; k++) {
                p[0][k] = point[nsd][jm-1][0][k];
                p[1][k] = point[nsd][jm-1][1][k];
                p[2][k] = point[nsd][0][im-1][k];
            }

            l_msk_crn[nsd][jm-1][0][0] = 1;
            area_corner_kites(p[0], p[1], p[2],
                              area_kite_crn[nsd][jm-1][0][0]);

            for (k=0; k<3; k++)
                wghts_crn[nsd][jm-1][0][0][k] =
                      area_kite_crn[nsd][jm-1][0][0][k] /
                      spherical_triangle_area(p[0], p[1], p[2]);

            set_vctr_wghts_crn(p[0], p[1], p[2],
                               vctr_wghts_crn[nsd][jm-1][0][0]);

            for (k=0; k<3; k++) {
                p[0][k] = point[nsd][jm-1][1][k];
                p[1][k] = point[nsd][jm-1][2][k];
                p[2][k] = point[nsd][0][im-1][k];
            }

            l_msk_crn[nsd][jm-1][1][0] = 1;
            area_corner_kites(p[0], p[1], p[2],
                              area_kite_crn[nsd][jm-1][1][0]);
            for (k=0; k<3; k++)
                wghts_crn[nsd][jm-1][1][0][k] =
                      area_kite_crn[nsd][jm-1][1][0][k] /
                      spherical_triangle_area(p[0], p[1], p[2]);

            set_vctr_wghts_crn(p[0], p[1], p[2],
                               vctr_wghts_crn[nsd][jm-1][1][0]);
        }
/*-----------------------------------------------------------------------
!  south pole
!-----------------------------------------------------------------------*/
        if (grid_subdomain->l_sbdmn_south_pole[nsd]) {

            for (k=0; k<3; k++) {
                p[0][k] = point[nsd][0   ][im-1][k];
                p[1][k] = point[nsd][jm-1][0   ][k];
                p[2][k] = point[nsd][1   ][im-1][k];
            }

            l_msk_crn[nsd][0][im-1][1] = 1;
            area_corner_kites(p[0], p[1], p[2],
                              area_kite_crn[nsd][0][im-1][1]);
            for (k=0; k<3; k++)
                wghts_crn[nsd][0][im-1][1][k] =
                      area_kite_crn[nsd][0][im-1][1][k] /
                      spherical_triangle_area(p[0], p[1], p[2]);

            set_vctr_wghts_crn(p[0], p[1], p[2],
                               vctr_wghts_crn[nsd][0][im-1][1]);

            for (k=0; k<3; k++) {
                p[0][k] = point[nsd][1   ][im-1][k];
                p[1][k] = point[nsd][jm-1][0   ][k];
                p[2][k] = point[nsd][2   ][im-1][k];
            }

            l_msk_crn[nsd][1][im-1][1] = 1;
            area_corner_kites(p[0], p[1], p[2],
                              area_kite_crn[nsd][1][im-1][1]);
            for (k=0; k<3; k++)
                wghts_crn[nsd][1][im-1][1][k] =
                      area_kite_crn[nsd][1][im-1][1][k] /
                      spherical_triangle_area(p[0], p[1], p[2]);

            set_vctr_wghts_crn(p[0], p[1], p[2],
                               vctr_wghts_crn[nsd][1][im-1][1]);
        }
    }
/*-----------------------------------------------------------------------
!  make as big as the earth
!-----------------------------------------------------------------------*/
    for (nsd=0; nsd<nsdm; nsd++) {
        for (j=0; j<jm; j++) {
            for (i=0; i<im; i++) {
                for (k=0; k<crnm; k++) {
                    for (l=0; l<2; l++) {
                        for (m=0; m<3; m++) {
                            vctr_wghts_crn [nsd][j][i][k][l][m] *= A;
                            vctr_wghts_crn2[nsd][j][i][k][l][m] *= A;
                        }
                    }
                    for (m=0; m<3; m++) {
                        area_kite_crn[nsd][j][i][k][m] *= A*A;
                    }
                }
            }
        }
    }
}

/*----< initialize_grid_metrics_edge() >--------------------------------------*/
static
void initialize_grid_metrics_edge(MODULE_grid_params  *grid_params,
                                  MODULE_wrap_data    *wrap_data,
                                  char       *wrap_name,
                                  double  ****point,     /* [nsdm][jm][im][3] */
                                  double *****point_crn, /* [nsdm][jm][im][crnm][3] */
                                  double *****point_edg, /* [nsdm][jm][im][edgm][3] */
                                  double  ****area_edg,  /* [nsdm][jm][im][edgm] */
                                  double *****nrm_edg,   /* [nsdm][jm][im][edgm][3] */
                                  double *****tng_edg)   /* [nsdm][jm][im][edgm][3] */
{
    int e, i, j, nsd, di[3], dj[3];
    double p_out[3];

    int im         = grid_params->im,
        jm         = grid_params->jm,
        nsdm       = grid_params->nsdm,
        edgm       = grid_params->edgm;

    di[0] = 1; di[1] = 1; di[2] = 0;
    dj[0] = 0; dj[1] = 1; dj[2] = 1;

    for (nsd=0; nsd<nsdm; nsd++) {
        for (j=1; j<jm-1; j++) {
            for (i=1; i<im-1; i++) {
                for (e=0; e<3; e++) {
                    mid_point(    point[nsd][j      ][i      ],
                                  point[nsd][j+dj[e]][i+di[e]],
                              point_edg[nsd][j      ][i      ][e]);
                    tangent_to_sphere(point_edg[nsd][j      ][i      ][e],
                                          point[nsd][j+dj[e]][i+di[e]],
                                        nrm_edg[nsd][j      ][i      ][e]);
                    cross_product(point_edg[nsd][j][i][e],
                                    nrm_edg[nsd][j][i][e],
                                  p_out);
                    unit_vector(p_out, tng_edg[nsd][j][i][e]);
                }
            }
        }
    }

    for (nsd=0; nsd<nsdm; nsd++) {
        for (j=1; j<jm-1; j++) {
            for (i=1; i<im-1; i++) {
                area_edg[nsd][j][i][0] =
                        spherical_triangle_area(point    [nsd][j  ][i  ],
                                                point_crn[nsd][j-1][i  ][1],
                                                point_crn[nsd][j  ][i  ][0]) +
                        spherical_triangle_area(point    [nsd][j  ][i+1],
                                                point_crn[nsd][j  ][i  ][0],
                                                point_crn[nsd][j-1][i  ][1]);
                area_edg[nsd][j][i][1] =
                        spherical_triangle_area(point    [nsd][j  ][i  ],
                                                point_crn[nsd][j  ][i  ][0],
                                                point_crn[nsd][j  ][i  ][1]) +
                        spherical_triangle_area(point    [nsd][j+1][i+1],
                                                point_crn[nsd][j  ][i  ][1],
                                                point_crn[nsd][j  ][i  ][0]);
                area_edg[nsd][j][i][2] =
                        spherical_triangle_area(point    [nsd][j  ][i  ],
                                                point_crn[nsd][j  ][i  ][1],
                                                point_crn[nsd][j  ][i-1][0]) +
                        spherical_triangle_area(point    [nsd][j+1][i  ],
                                                point_crn[nsd][j  ][i-1][0],
                                                point_crn[nsd][j  ][i  ][1]);
            }
        }
    }

    int dimlen[5];
    dimlen[0] = nsdm;
    dimlen[1] = jm;
    dimlen[2] = im;
    dimlen[3] = edgm;
    dimlen[4] = 3;
    wrap_edge_1lyr(wrap_data, wrap_name, point_edg, dimlen);
    wrap_edge_1lyr(wrap_data, wrap_name, nrm_edg, dimlen);
    wrap_edge_1lyr(wrap_data, wrap_name, tng_edg, dimlen);

    wrap_edge_scalar_1lyr(wrap_data, wrap_name, area_edg, dimlen);
}

/*----< set_grid_point() >----------------------------------------------------*/
static
void set_grid_point(char      *grid_point_select,
                    grid_node *ptr,
                    int        level_max,
                    int        level)
{
    int n, m;
    double pert_coeff;

#define PI 3.14159265358979323846264338327950288
    pert_coeff = 0.1 * (2.0 * PI) / (10.0*POWER2(level_max-1));

    if (!strcmp(grid_point_select, "bisect")) {
        if (ptr->level == level) {
            if (ptr->l_point_set) {
                if (ptr->dn[0] != NULL) {
                    ptr->dn[0]->l_point_set = 1;
                    for (n=0; n<3; n++)
                        ptr->dn[0]->point[n] = ptr->point[n];
                }
                for (m=1; m<=3; m++) {
                    if ((ptr->nghbr[m] != NULL) && (ptr->dn[m] != NULL)) {
                        if (ptr->nghbr[m]->l_point_set) {
                            ptr->dn[m]->l_point_set = 1;
                            mid_point(ptr->nghbr[0]->point,
                                      ptr->nghbr[m]->point,
                                      ptr->dn[m]->point);
/*
                            if (0) {
                                double perturbation[3];
                                CALL RANDOM_NUMBER (perturbation)
                                for (n=0; n<3; n++) {
                                    perturbation[n] = 2.0 * (perturbation[n]- 0.5);
                                    ptr->dn[m].p.point[n] +=
                                               pert_coeff * perturbation[n];
                                }
                                unit_vector(ptr->dn[m].p.point);
                            }
*/
                        }
                    }
                }
            }
        }
    }

    for (n=0; n<4; n++)
        if (ptr->dn[n] != NULL)
            set_grid_point(grid_point_select, ptr->dn[n], level_max, level);
}

/* rotate pentagon 1 and pentagon 3 to the equator */
#define alphalph 0.553574358897045251508532730089268520035023822700716

/*----< rotate_grid() >-------------------------------------------------------*/
void rotate_grid(grid_node *ptr)
{
    int n;
    double rho, rotation_matrix_x[3][3], rotation_matrix_y[3][3], tmp[3];

    rho = 0.0 * PI / 4.0;
    rotation_matrix_x[0][0] = 1.0;
    rotation_matrix_x[0][1] = 0.0;
    rotation_matrix_x[0][2] = 0.0;
    rotation_matrix_x[1][0] = 0.0;
    rotation_matrix_x[1][1] =  cos(rho);
    rotation_matrix_x[1][2] = -sin(rho);
    rotation_matrix_x[2][0] = 0.0;
    rotation_matrix_x[2][1] = sin(rho);
    rotation_matrix_x[2][2] = cos(rho);

    mat_vec_mul(3, 3, rotation_matrix_x, ptr->point, tmp);
    unit_vector(tmp, tmp);
    memcpy(ptr->point, tmp, 3*sizeof(double));

    rho = alphalph;
    rotation_matrix_y[0][0] = cos(rho);
    rotation_matrix_y[0][1] = 0.0;
    rotation_matrix_y[0][2] = -sin(rho);
    rotation_matrix_y[1][0] = 0.0;
    rotation_matrix_y[1][1] = 1.0;
    rotation_matrix_y[1][2] = 0.0;
    rotation_matrix_y[2][0] = sin(rho);
    rotation_matrix_y[2][1] = 0.0;
    rotation_matrix_y[2][2] = cos(rho);

    mat_vec_mul(3, 3, rotation_matrix_y, ptr->point, tmp);
    unit_vector(tmp, tmp);
    memcpy(ptr->point, tmp, 3*sizeof(double));

    for (n=0; n<4; n++)
        if (ptr->dn[n] != NULL)
            rotate_grid(ptr->dn[n]);
}

