/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio_register.c 4604 2017-12-07 07:15:42Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "gio.h"
#include "util.h"

/*----< Morton_index() >------------------------------------------------------*/
static
int Morton_index(int *lo,           /* [2] */
                 int *hi,           /* [2] */
                 int  panel_id,
                 int  refine_param,
                 int  block_size)
{
    /* Evaluate Morton-ordered block index */
    int i, i_index, j_index, nblock, lmorton;

    i_index = (lo[0] - 1) / block_size;
    j_index = (lo[1] - 1) / block_size;
    nblock = POWER2(refine_param)/block_size;
    lmorton = 0;
    while (nblock > 1) {
        lmorton++;
        nblock /= 2;
    }
    int *i_binary = (int*) malloc(lmorton * sizeof(int));
    int *j_binary = (int*) malloc(lmorton * sizeof(int));

    for (i=0; i<lmorton; i++) {
        i_binary[i] = i_index % 2;
        j_binary[i] = j_index % 2;
        i_index = (i_index - i_binary[i]) / 2;
        j_index = (j_index - j_binary[i]) / 2;
    }

    i_index = 0;
    for (i=lmorton-1; i>=0; i--)
        i_index = 4*i_index + (2*j_binary[i] + i_binary[i]);

    free(i_binary);
    free(j_binary);

    nblock = POWER2(refine_param)/block_size;

    return panel_id*nblock*nblock + i_index;
}

/*----< gio_register_dfield() >----------------------------------------------*/
void gio_register_dfield(char           *field_name,
                         double         *darray, 
                         int             lo[2],     /* 1-based */
                         int             hi[2],     /* 1-based */
                         int             panel_id,  /* 0-based */
                         gio_parameters *gio_param)
{
    int i, idx, nblk;

    gio_data     *data     = &gio_param->data;
    gio_grid     *grid     = &gio_param->grid;

    data_descriptor *gio_descriptors = data->gio_descriptors;

    /* Find index to appropriate data descriptor */
    for (i=0; i<data->num_fields; i++) {
        if (!strcmp(gio_descriptors[i].field_name, field_name))
            break;
    }
    if (i == data->num_fields) {
        printf("Error: gio_register_dfield - could not find descriptor: %s\n",
               field_name);
        ABORT
    }
    idx = i;

    nblk = gio_descriptors[idx].num_data_blocks;
    gio_descriptors[idx].num_data_blocks++;

    gio_descriptors[idx].gio_data_blocks[nblk].double_data = darray;
    gio_descriptors[idx].assigned = 1;
    gio_descriptors[idx].gio_data_blocks[nblk].imin = lo[0];
    gio_descriptors[idx].gio_data_blocks[nblk].jmin = lo[1];
    gio_descriptors[idx].gio_data_blocks[nblk].imax = hi[0];
    gio_descriptors[idx].gio_data_blocks[nblk].jmax = hi[1];
    gio_descriptors[idx].gio_data_blocks[nblk].panel_id = panel_id;

    /* Evaluate Morton-ordered block index */
    gio_descriptors[idx].gio_data_blocks[nblk].block_index =
        Morton_index(lo, hi, panel_id, grid->refine_param, grid->block_size);
}

/*----< gio_register_ifield() >----------------------------------------------*/
void gio_register_ifield(char           *field_name,
                         int            *iarray, 
                         int             lo[2],
                         int             hi[2],
                         int             panel_id,
                         gio_parameters *gio_param)
{
    int i, idx, nblk;

    gio_data     *data     = &gio_param->data;
    gio_grid     *grid     = &gio_param->grid;

    data_descriptor *gio_descriptors = data->gio_descriptors;

    /* Find index to appropriate data descriptor */
    for (i=0; i<data->num_fields; i++) {
        if (!strcmp(gio_descriptors[i].field_name, field_name))
            break;
    }
    if (i == data->num_fields) {
        printf("Error: gio_register_ifield - could not find descriptor: %s\n",
               field_name);
        ABORT
    }
    idx = i;

    nblk = gio_descriptors[idx].num_data_blocks;
    gio_descriptors[idx].num_data_blocks++;

    gio_descriptors[idx].gio_data_blocks[nblk].integer_data = iarray;
    gio_descriptors[idx].assigned = 1;
    gio_descriptors[idx].gio_data_blocks[nblk].imin = lo[0];
    gio_descriptors[idx].gio_data_blocks[nblk].jmin = lo[1];
    gio_descriptors[idx].gio_data_blocks[nblk].imax = hi[0];
    gio_descriptors[idx].gio_data_blocks[nblk].jmax = hi[1];
    gio_descriptors[idx].gio_data_blocks[nblk].panel_id = panel_id;

    /* Evaluate Morton-ordered block index */
    gio_descriptors[idx].gio_data_blocks[nblk].block_index =
        Morton_index(lo, hi, panel_id, grid->refine_param, grid->block_size);
}

/*----< gio_register_dpole() >-----------------------------------------------*/
void gio_register_dpole(char     *field_name,
                        double   *dval,
                        char      pole_id,
                        gio_data *data)
{
    int i, idx;

    data_descriptor *gio_descriptors = data->gio_descriptors;

    /* Find index to appropriate data descriptor */
    for (i=0; i<data->num_fields; i++) {
        if (!strcmp(gio_descriptors[i].field_name, field_name))
            break;
    }
    if (i == data->num_fields) {
        printf("Error: gio_register_dpole - could not find descriptor: %s\n",
               field_name);
        ABORT
    }
    idx = i;

    if (pole_id == 'N' || pole_id == 'n') {
        gio_descriptors[idx].double_north_data = dval;
        gio_descriptors[idx].assigned = 1;
    }
    else if (pole_id == 'S' || pole_id == 's') {
        gio_descriptors[idx].double_south_data = dval;
        gio_descriptors[idx].assigned = 1;
    }
    else {
        printf("Error: gio_register_dpole - could not find pole_id: %c\n",
               pole_id);
        ABORT
    }
}

/*----< gio_register_ipole() >-----------------------------------------------*/
void gio_register_ipole(char     *field_name,
                        int      *ival,
                        char      pole_id,
                        gio_data *data)
{
    int i, idx;

    data_descriptor *gio_descriptors = data->gio_descriptors;

    /* Find index to appropriate data descriptor */
    for (i=0; i<data->num_fields; i++) {
        if (!strcmp(gio_descriptors[i].field_name, field_name))
            break;
    }
    if (i == data->num_fields) {
        printf("Error: gio_register_ipole - could not find descriptor: %s\n",
               field_name);
        ABORT
    }
    idx = i;

    if (pole_id == 'N' || pole_id == 'n') {
        gio_descriptors[idx].int_north_data = ival;
        gio_descriptors[idx].assigned = 1;
    }
    else if (pole_id == 'S' || pole_id == 's') {
        gio_descriptors[idx].int_south_data = ival;
        gio_descriptors[idx].assigned = 1;
    }
    else {
        printf("Error: gio_register_dpole - could not find pole_id: %c\n",
               pole_id);
        ABORT
    }
}

/*----< gio_register_dlevel() >----------------------------------------------*/
void gio_register_dlevel(gio_data *data,
                         char     *field_name,
                         double   *dvals,
                         int       length)
{
    int i, idx, nlvl;

    data_descriptor *gio_descriptors = data->gio_descriptors;

    /* Find index to appropriate data descriptor */
    for (i=0; i<data->num_fields; i++) {
        if (!strcmp(gio_descriptors[i].field_name, field_name))
            break;
    }
    if (i == data->num_fields) {
        printf("Error: gio_register_dfield - could not find descriptor: %s\n",
               field_name);
        ABORT
    }
    idx = i;

    if (gio_descriptors[idx].has_levels) {
        gio_descriptors[idx].float_level_data = calloc_1D_flt(length);

        /* copy the double precision data into our single precision array */

        for (nlvl=0; nlvl<length; nlvl++)
            gio_descriptors[idx].float_level_data[nlvl] = dvals[nlvl];

        if (strcmp(field_name, "layers") == 0)
            data->layers_length = length;
        else if (strcmp(field_name, "interfaces") == 0)
            data->interfaces_length = length;

        gio_descriptors[idx].is_level_data = 1;
        gio_descriptors[idx].assigned = 1;
    }
    else {
        printf("Error: gio_register_dlevel - field not declared as level data: %s\n",
               field_name);
        ABORT
    }
}


/*----< gio_get_panel_index() >----------------------------------------------*/
/* Find panel ID and indices (i,j) within the panel for a given cell with
 * index tag_label. This is the label coming from the GCRM code.
 *   level     : maximum level for grid (input)
 *   panel_size: 2**(2*level) (save time recomputing this quantity by
 *               adding as an argument (input)
 *   tag_label : GCRM label for cell (input)
 *   i,j       : i,j indices of cell in panel (output)
 *   panel_id  : panel in which cell is located (output) 
 */
static
void gio_get_panel_index(int  level,
                         int  panel_size,
                         int  tag_label,
                         int *i,
                         int *j,
                         int *panel_id)
{
    int ii, n, nmx, ifac;

    /* Assume that this function will not be called for north and south pole
     * (tag_label <= 2) */

    n         = tag_label - 3;
    nmx       = panel_size;
    *panel_id = n/nmx;

    *i = 0;
    *j = 0;
    for (ii=1; ii<=level; ii++) {
        n %= nmx;
        nmx /= 4;

        ifac = POWER2(level-ii);
        switch (n/nmx) {
            case 0:                         break;
            case 1: *i += ifac;             break;
            case 2: *i += ifac; *j += ifac; break;
            case 3:             *j += ifac; break;
            default:                        break;
        }
    }
    (*i)++;
    (*j)++;
}

/*----< gio_morton_index() >-------------------------------------------------*/
static
int gio_morton_index(int x,
                     int y,
                     int panel_id,
                     int lsize,
                     int res)
{
    int k, idx;
    int xbit[20], ybit[20];

    x--;
    y--;
    for (k=0; k<res; k++) {
        xbit[k] = x % 2;
        x = (x - xbit[k]) / 2;
        ybit[k] = y % 2;
        y = (y - ybit[k]) / 2;
    }
    idx = 0;
    for (k=res-1; k>=0; k--)
        idx = 4*idx + (2*ybit[k] + xbit[k]);

    /* return (idx + panel_id*lsize + 3); ! 1 based index */
    return (idx + panel_id*lsize + 2);   /* 0 based index */
}

/*----< gio_grid_setup() >---------------------------------------------------*/
void gio_grid_setup(gio_parameters   *gio_param,
                    int             **nghbr_tags, /* [][6] */
                    double          **center,     /* [][3] */
                    double         ***corner,     /* [][6][3] */
                    double           *grid_area,  /* [gio_istride*gio_jstride] */
                    int               panel_id,
                    int               gmin[2], 
                    int               gmax[2])
{
    /* Need ghost cells for tag array, hence dimensions */
    int idx, i, j, k, ii, jj;
    int ldx, ndx;
    int iii, jjj, panel_size, nvals, res;
    int is_northern;
    int **nghbr_ptr, **crnr_ptr, **edge_ptr;
    int ***end_ptr;
    double x, y, z, xc, yc, zc, x2, y2, z2, rn;
    double starttime = MPI_Wtime();

    gio_grid     *grid     = &gio_param->grid;

    int gio_istride = grid->gio_istride;
    int gio_jstride = grid->gio_jstride;

    /*  Add another grid block to count */
    idx = grid->num_grid_blocks;
    grid->num_grid_blocks++;

    if (panel_id % 2 == 0)
        is_northern = 1;
    else
        is_northern = 0;

    grid_block *gio_grid_blocks = grid->gio_grid_blocks;

    /* Evaluate Morton-ordered block index */
    gio_grid_blocks[idx].block_index = 
        Morton_index(gmin, gmax, panel_id, grid->refine_param, grid->block_size);

/* should the indices be reduced by 1?
gmin[0]--; gmin[1]--;
gmax[0]--; gmax[1]--;
*/
    int ij_stride = gio_istride * gio_jstride;
    /* convert cell centers to lat-long coordinates */

    gio_grid_blocks[idx].grid_center_lat = calloc_1D_dbl(ij_stride);
    gio_grid_blocks[idx].grid_center_lon = calloc_1D_dbl(ij_stride);

    /* computing grid_center_lat ... */
    for (j=gmin[1]; j<=gmax[1]; j++) {
        jj = j - gmin[1];
        for (i=gmin[0]; i<=gmax[0]; i++) {
            ii   = i - gmin[0];
            ldx  = jj*gio_istride + ii;
            x    = center[ldx][0];
            y    = center[ldx][1];
            z    = center[ldx][2];
            if (z >  1.0) z =  1.0;
            if (z < -1.0) z = -1.0;
            if (x != 0.0 || y != 0.0)
                gio_grid_blocks[idx].grid_center_lon[ldx] = atan2(y,x);
            else
                gio_grid_blocks[idx].grid_center_lon[ldx] = 0.0;
            gio_grid_blocks[idx].grid_center_lat[ldx] = asin(z);
        }
    }

    /* convert cell corners to lat-long coordinates */
    gio_grid_blocks[idx].grid_corner_lat = calloc_2D_dbl(ij_stride, 6);
    gio_grid_blocks[idx].grid_corner_lon = calloc_2D_dbl(ij_stride, 6);

    /* computing grid_corner_lat ... */
    for (j=gmin[1]; j<=gmax[1]; j++) {
        jj = j - gmin[1];
        for (i=gmin[0]; i<=gmax[0]; i++) {
            ii = i - gmin[0];
            ldx  = jj*gio_istride + ii;
            for (k=0; k<6; k++) {
                x = corner[ldx][k][0];
                y = corner[ldx][k][1];
                z = corner[ldx][k][2];
                if (z >  1.0) z =  1.0;
                if (z < -1.0) z = -1.0;
                if (x != 0.0 || y != 0.0)
                    gio_grid_blocks[idx].grid_corner_lon[ldx][k] = atan2(y,x);
                else
                    gio_grid_blocks[idx].grid_corner_lon[ldx][k] = 0.0;
                gio_grid_blocks[idx].grid_corner_lat[ldx][k] = asin(z);
            }
        }
    }

    /* Create cell edge center array from corner array */
    gio_grid_blocks[idx].grid_edge_lat = calloc_2D_dbl(ij_stride, 6);
    gio_grid_blocks[idx].grid_edge_lon = calloc_2D_dbl(ij_stride, 6);

    for (j=gmin[1]; j<=gmax[1]; j++) {
        jj = j - gmin[1];
        for (i=gmin[0]; i<=gmax[0]; i++) {
            ii = i - gmin[0];
            ldx  = jj*gio_istride + ii;
            for (k=0; k<6; k++) {
                if (k > 0) {
                    x  = corner[ldx][k-1][0];
                    y  = corner[ldx][k-1][1];
                    z  = corner[ldx][k-1][2];
                    x2 = corner[ldx][k][0];
                    y2 = corner[ldx][k][1];
                    z2 = corner[ldx][k][2];
                } else {
                    x  = corner[ldx][5][0];
                    y  = corner[ldx][5][1];
                    z  = corner[ldx][5][2];
                    x2 = corner[ldx][0][0];
                    y2 = corner[ldx][0][1];
                    z2 = corner[ldx][0][2];
                }
                xc = 0.5*(x + x2);
                yc = 0.5*(y + y2);
                zc = 0.5*(z + z2);
                rn = sqrt(xc*xc+yc*yc+zc*zc);
                xc = xc/rn;
                yc = yc/rn;
                zc = zc/rn;
                if (zc >  1.0) zc =  1.0;
                if (zc < -1.0) zc = -1.0;
                if (xc != 0.0 || yc != 0.0)
                    gio_grid_blocks[idx].grid_edge_lon[ldx][k] = atan2(yc,xc);
                else
                    gio_grid_blocks[idx].grid_edge_lon[ldx][k] = 0.0;
                gio_grid_blocks[idx].grid_edge_lat[ldx][k] = asin(zc);
            }
        }
    }

    /* assign pointer to grid_area */
    gio_grid_blocks[idx].grid_area = grid_area;

    /* construct neighbor list */
    gio_grid_blocks[idx].grid_neighbors = calloc_2D_int(ij_stride, 6);

    nghbr_ptr = gio_grid_blocks[idx].grid_neighbors;
    panel_size = grid->grid_length * grid->grid_length;

    for (j=gmin[1]; j<=gmax[1]; j++) {
        jj = j - gmin[1];
        for (i=gmin[0]; i<=gmax[0]; i++) {
            ii = i - gmin[0];
            ldx  = jj*gio_istride + ii;

            /* Relabel neighbors from GCRM convention to Morton-ordering
               convention */
            for (k=0; k<6; k++) {
                ndx = nghbr_tags[ldx][k];
                if (ndx > 2) {
                    gio_get_panel_index(grid->refine_param, panel_size, ndx,
                                        &iii, &jjj, &panel_id);
                    ndx = gio_morton_index(iii, jjj, panel_id, panel_size,
                                           grid->refine_param);
                } else
                  ndx--;
                nghbr_ptr[ldx][k] = ndx;
            }
        }
    }
    /* Loop over grid cells and construct corner array for all cells on
     * panels */
    gio_grid_blocks[idx].grid_corners = calloc_2D_int(ij_stride, 6);
    crnr_ptr = gio_grid_blocks[idx].grid_corners;

    /*  Check if block has pentagon */
    for (j=gmin[1]; j<=gmax[1]; j++) {
        jj = j - gmin[1];
        for (i=gmin[0]; i<=gmax[0]; i++) {
            ii = i - gmin[0];
            ldx  = jj*gio_istride + ii;
            ndx = gio_morton_index(i, j, panel_id, panel_size,
                                   grid->refine_param);

            if (i == 1 && j == 1 && is_northern) {
                /*
                 *  cell is pentagon in northern hemisphere
                 */
                crnr_ptr[ldx][0] = 2*(ndx-2);
                crnr_ptr[ldx][1] = 2*(ndx-2) + 1;
                ndx = 2*(nghbr_ptr[ldx][4]-2) + 1;
                crnr_ptr[ldx][2] = ndx;
                crnr_ptr[ldx][3] = ndx;
                ndx = 2*(nghbr_ptr[ldx][4]-2);
                crnr_ptr[ldx][4] = ndx;
                ndx = 2*(nghbr_ptr[ldx][5]-2) + 1;
                crnr_ptr[ldx][5] = ndx;
            }
            else if (i == 1 && j == 1 && !is_northern) {
                /*
                 *  cell is pentagon in southern hemisphere
                 */
                crnr_ptr[ldx][0] = 2*(ndx-2);
                crnr_ptr[ldx][1] = 2*(ndx-2) + 1;
                ndx = 2*(nghbr_ptr[ldx][3]-2);
                crnr_ptr[ldx][2] = ndx;
                ndx = 2*(nghbr_ptr[ldx][4]-2) + 1;
                crnr_ptr[ldx][3] = ndx;
                ndx = 2*(nghbr_ptr[ldx][4]-2);
                crnr_ptr[ldx][4] = ndx;
                crnr_ptr[ldx][5] = ndx;
            }
            else if (i == 1 && j != 1 && is_northern) {
                /*
                 * cell is hexagon along northwest edge in northern hemisphere
                 */
                crnr_ptr[ldx][0] = 2*(ndx-2);
                crnr_ptr[ldx][1] = 2*(ndx-2) + 1;
                ndx = 2*(nghbr_ptr[ldx][3]-2) + 1;
                crnr_ptr[ldx][2] = ndx;
                ndx = 2*(nghbr_ptr[ldx][3]-2);
                crnr_ptr[ldx][3] = ndx;
                ndx = 2*(nghbr_ptr[ldx][4]-2) + 1;
                crnr_ptr[ldx][4] = ndx;
                ndx = 2*(nghbr_ptr[ldx][5]-2) + 1;
                crnr_ptr[ldx][5] = ndx;
            }
            else if (i != 1 && j == 1 && !is_northern) {
                /*
                 * cell is hexagon along southwest edge in southern hemisphere
                 */
                crnr_ptr[ldx][0] = 2*(ndx-2);
                crnr_ptr[ldx][1] = 2*(ndx-2) + 1;
                ndx = 2*(nghbr_ptr[ldx][3]-2);
                crnr_ptr[ldx][2] = ndx;
                ndx = 2*(nghbr_ptr[ldx][4]-2);
                crnr_ptr[ldx][3] = ndx;
                ndx = 2*(nghbr_ptr[ldx][5]-2) + 1;
                crnr_ptr[ldx][4] = ndx;
                ndx = 2*(nghbr_ptr[ldx][5]-2);
                crnr_ptr[ldx][5] = ndx;
            }
            else {
                /*
                 * cell is regular hexagon
                 */
                crnr_ptr[ldx][0] = 2*(ndx-2);
                crnr_ptr[ldx][1] = 2*(ndx-2) + 1;
                ndx = 2*(nghbr_ptr[ldx][3]-2);
                crnr_ptr[ldx][2] = ndx;
                ndx = 2*(nghbr_ptr[ldx][4]-2) + 1;
                crnr_ptr[ldx][3] = ndx;
                ndx = 2*(nghbr_ptr[ldx][4]-2);
                crnr_ptr[ldx][4] = ndx;
                ndx = 2*(nghbr_ptr[ldx][5]-2) + 1;
                crnr_ptr[ldx][5] = ndx;
            }
        }
    }

    /* Loop over grid cells and construct edge array for all cells on
     * panels */
    gio_grid_blocks[idx].grid_edges = calloc_2D_int(ij_stride, 6);
    edge_ptr = gio_grid_blocks[idx].grid_edges;

    /*  Check if block has pentagon */
    for (j=gmin[1]; j<=gmax[1]; j++) {
        jj = j - gmin[1];
        for (i=gmin[0]; i<=gmax[0]; i++) {
            ii = i - gmin[0];
            ldx  = jj*gio_istride + ii;
            ndx = gio_morton_index(i, j, panel_id, panel_size,
                                   grid->refine_param);

            if (i == 1 && j == 1 && is_northern) {
                /*
                 *  cell is pentagon in northern hemisphere
                 */
                edge_ptr[ldx][0] = 3*(ndx-2);
                edge_ptr[ldx][1] = 3*(ndx-2) + 1;
                edge_ptr[ldx][2] = 3*(ndx-2) + 2; ndx = 3*(nghbr_ptr[ldx][4]-2) + 1;
                edge_ptr[ldx][3] = ndx;
                edge_ptr[ldx][4] = ndx;           ndx = 3*(nghbr_ptr[ldx][5]-2) + 2;
                edge_ptr[ldx][5] = ndx;
            }
            else if (i == 1 && j == 1 && !is_northern) {
                /*
                 *  cell is pentagon in southern hemisphere
                 */
                edge_ptr[ldx][0] = 3*(ndx-2);
                edge_ptr[ldx][1] = 3*(ndx-2) + 1;
                edge_ptr[ldx][2] = 3*(ndx-2) + 2; ndx = 3*(nghbr_ptr[ldx][3]-2);
                edge_ptr[ldx][3] = ndx;           ndx = 3*(nghbr_ptr[ldx][4]-2) + 1;
                edge_ptr[ldx][4] = ndx;
                edge_ptr[ldx][5] = ndx;
            }
            else if (i == 1 && j != 1 && is_northern) {
                /*
                 * cell is hexagon along northwest edge in northern hemisphere
                 */
                edge_ptr[ldx][0] = 3*(ndx-2);
                edge_ptr[ldx][1] = 3*(ndx-2) + 1;
                edge_ptr[ldx][2] = 3*(ndx-2) + 2; ndx = 3*(nghbr_ptr[ldx][3]-2) + 1;
                edge_ptr[ldx][3] = ndx;           ndx = 3*(nghbr_ptr[ldx][4]-2) + 2;
                edge_ptr[ldx][4] = ndx;           ndx = 3*(nghbr_ptr[ldx][5]-2) + 2;
                edge_ptr[ldx][5] = ndx;
            }
            else if (i != 1 && j == 1 && !is_northern) {
                /*
                 * cell is hexagon along southwest edge in southern hemisphere
                 */
                edge_ptr[ldx][0] = 3*(ndx-2);
                edge_ptr[ldx][1] = 3*(ndx-2) + 1;
                edge_ptr[ldx][2] = 3*(ndx-2) + 2; ndx = 3*(nghbr_ptr[ldx][3]-2);
                edge_ptr[ldx][3] = ndx;           ndx = 3*(nghbr_ptr[ldx][4]-2);
                edge_ptr[ldx][4] = ndx;           ndx = 3*(nghbr_ptr[ldx][5]-2) + 1;
                edge_ptr[ldx][5] = ndx;
            }
            else {
                /*
                 * cell is regular hexagon
                 */
                edge_ptr[ldx][0] = 3*(ndx-2);
                edge_ptr[ldx][1] = 3*(ndx-2) + 1;
                edge_ptr[ldx][2] = 3*(ndx-2) + 2; ndx = 3*(nghbr_ptr[ldx][3]-2);
                edge_ptr[ldx][3] = ndx;           ndx = 3*(nghbr_ptr[ldx][4]-2) + 1;
                edge_ptr[ldx][4] = ndx;           ndx = 3*(nghbr_ptr[ldx][5]-2) + 2;
                edge_ptr[ldx][5] = ndx;
            }
        }
    }

    /* Loop over grid cells and construct array that represents the two corners
     * for each unique edge associated with the panel */
    gio_grid_blocks[idx].grid_endpoints = calloc_3D_int(ij_stride, 3, 2);
    end_ptr = gio_grid_blocks[idx].grid_endpoints;

    nvals = grid->block_size;
    res = 0;
    while (nvals > 1) {
        res++;
        nvals /= 2;
    }

    for (j=gmin[1]; j<=gmax[1]; j++) {
        jj = j - gmin[1];
        for (i=gmin[0]; i<=gmax[0]; i++) {
            ii = i - gmin[0];
            ldx  = jj*gio_istride + ii;
            end_ptr[ldx][0][0] = crnr_ptr[ldx][5];
            end_ptr[ldx][0][1] = crnr_ptr[ldx][0];
            end_ptr[ldx][1][0] = crnr_ptr[ldx][0];
            end_ptr[ldx][1][1] = crnr_ptr[ldx][1];
            end_ptr[ldx][2][0] = crnr_ptr[ldx][1];
            end_ptr[ldx][2][1] = crnr_ptr[ldx][2];
        }
    }

    gio_grid_blocks[idx].panel_id = panel_id;
    gio_grid_blocks[idx].imin = gmin[0];
    gio_grid_blocks[idx].jmin = gmin[1];
    gio_grid_blocks[idx].imax = gmax[0];
    gio_grid_blocks[idx].jmax = gmax[1];

#define REG_DFIELD(name, buf) gio_register_dfield(name, buf, gmin, gmax, panel_id, gio_param);
#define REG_IFIELD(name, buf) gio_register_ifield(name, buf, gmin, gmax, panel_id, gio_param);

    REG_DFIELD("grid_center_lat",    gio_grid_blocks[idx].grid_center_lat)
    REG_DFIELD("grid_center_lon",    gio_grid_blocks[idx].grid_center_lon)
    REG_DFIELD("corner_cell_map_lat",gio_grid_blocks[idx].grid_corner_lat[0])
    REG_DFIELD("corner_cell_map_lon",gio_grid_blocks[idx].grid_corner_lon[0])
    REG_DFIELD("grid_corner_lat",    gio_grid_blocks[idx].grid_corner_lat[0])
    REG_DFIELD("grid_corner_lon",    gio_grid_blocks[idx].grid_corner_lon[0])
    REG_DFIELD("grid_edge_lat",      gio_grid_blocks[idx].grid_edge_lat[0])
    REG_DFIELD("grid_edge_lon",      gio_grid_blocks[idx].grid_edge_lon[0])
    REG_DFIELD("area",               gio_grid_blocks[idx].grid_area)

    REG_IFIELD("cell_neighbors",     gio_grid_blocks[idx].grid_neighbors[0])
    REG_IFIELD("cell_corners",       crnr_ptr[0])
    REG_IFIELD("cell_edges",         edge_ptr[0])
    REG_IFIELD("edge_corners",       end_ptr[0][0])

    gio_param->total_time_in_API += MPI_Wtime() - starttime;
}

/*----< gio_grid_setup_pole() >----------------------------------------------*/
void gio_grid_setup_pole(gio_parameters *gio_param,
                         int            *nghbr_tags, /* [6] */
                         double         *center,     /* [3] */
                         double        **corner,     /* [6][3] */
                         double          grid_area,
                         char            pole_id)
{
    int iii, jjj, k, ndx, panel_id, panel_size;
    int ijstride=0;
    double x, y, z, x2, y2, z2, xc, yc, zc, rn;
    double starttime = MPI_Wtime();

    gio_grid   *grid = &gio_param->grid;
    grid_block *gio_grid_blocks = grid->gio_grid_blocks;

    panel_size = grid->grid_length * grid->grid_length;

    if (pole_id == 'n' || pole_id == 'N') {
        ijstride = (grid->gio_jstride-1) * grid->gio_istride + 1;
        grid->north_pole_flag = 1;
    }
    else if (pole_id == 's' || pole_id == 'S') {
        ijstride = grid->gio_istride *2 - 1;
        grid->south_pole_flag = 1;
    }
    else {
        printf("Failure in grid_setup_pole - could not find pole_id for: %c\n",pole_id);
        ABORT
    }

    if (center[0] != 0.0 || center[1] != 0.0)
        gio_grid_blocks[0].grid_center_lon[ijstride]
                = atan2(center[1], center[0]);
    else
        gio_grid_blocks[0].grid_center_lon[ijstride] = 0.0;

    z = center[2];
    if (z >  1.0) z =  1.0;
    if (z < -1.0) z = -1.0;
    gio_grid_blocks[0].grid_center_lat[ijstride] = asin(z);

    for (k=0; k<6; k++) {
        if (corner[k][0] != 0.0 || corner[k][1] != 0.0)
            gio_grid_blocks[0].grid_corner_lon[ijstride][k]
                  = atan2(corner[k][1], corner[k][0]);
        else
            gio_grid_blocks[0].grid_corner_lon[ijstride][k] = 0.0;

        z = corner[k][2];
        if (z >  1.0) z =  1.0;
        if (z < -1.0) z = -1.0;
        gio_grid_blocks[0].grid_corner_lat[ijstride][k] = asin(z);
    }
    for (k=0; k<6; k++) {
        if (k > 0) {
            x  = corner[k-1][0];
            y  = corner[k-1][1];
            z  = corner[k-1][2];
            x2 = corner[k][0];
            y2 = corner[k][1];
            z2 = corner[k][2];
        } else {
            x  = corner[5][0];
            y  = corner[5][1];
            z  = corner[5][2];
            x2 = corner[0][0];
            y2 = corner[0][1];
            z2 = corner[0][2];
        }
        xc = 0.5 * (x + x2);
        yc = 0.5 * (y + y2);
        zc = 0.5 * (z + z2);
        rn = sqrt(xc*xc + yc*yc + zc*zc);
        xc = xc/rn;
        yc = yc/rn;
        zc = zc/rn;
        if (zc >  1.0) zc =  1.0;
        if (zc < -1.0) zc = -1.0;
        if (xc != 0.0 || yc != 0.0)
            gio_grid_blocks[0].grid_edge_lon[ijstride][k] = atan2(yc,xc);
        else
          gio_grid_blocks[0].grid_edge_lon[ijstride][k] = 0.0;

        gio_grid_blocks[0].grid_edge_lat[ijstride][k] = asin(zc);
    }
    gio_grid_blocks[0].grid_area[ijstride] = grid_area;

    for (k=0; k<6; k++)
        gio_grid_blocks[0].grid_neighbors[ijstride][k] = nghbr_tags[k];

    /* Relabel neighbors from GCRM convention to Morton-ordering
     *  convention */
    for (k=0; k<6; k++) {
        ndx = gio_grid_blocks[0].grid_neighbors[ijstride][k];
        gio_get_panel_index(grid->refine_param, panel_size, ndx,
                            &iii, &jjj, &panel_id);
        ndx = gio_morton_index(iii, jjj, panel_id, panel_size,
                            grid->refine_param);
        gio_grid_blocks[0].grid_neighbors[ijstride][k] = ndx;
    }

    /* Create corner list */
    for (k=0; k<6; k++) {
        ndx = gio_grid_blocks[0].grid_neighbors[ijstride][k];
        if (grid->north_pole_flag)
            gio_grid_blocks[0].grid_corners[ijstride][k] = 2*(ndx-2) + 1;
        else
            gio_grid_blocks[0].grid_corners[ijstride][k] = 2*(ndx-2);
    }
    /* Create edge list */
    for (k=0; k<6; k++) {
        ndx = gio_grid_blocks[0].grid_neighbors[ijstride][k];
        if (grid->north_pole_flag)
            gio_grid_blocks[0].grid_edges[ijstride][k] = 3*(ndx-2) + 2;
        else
            gio_grid_blocks[0].grid_edges[ijstride][k] = 3*(ndx-2);
    }

#define REG_DPOLE(name, buf) gio_register_dpole(name, buf, pole_id, &gio_param->data);
#define REG_IPOLE(name, buf) gio_register_ipole(name, buf, pole_id, &gio_param->data);
    k = ijstride;
    REG_DPOLE("grid_center_lat",     &gio_grid_blocks[0].grid_center_lat[k])
    REG_DPOLE("grid_center_lon",     &gio_grid_blocks[0].grid_center_lon[k])
    REG_DPOLE("corner_cell_map_lat",  gio_grid_blocks[0].grid_corner_lat[k])
    REG_DPOLE("corner_cell_map_lon",  gio_grid_blocks[0].grid_corner_lon[k])
    REG_DPOLE("area",                &gio_grid_blocks[0].grid_area[k])

    REG_IPOLE("cell_neighbors",       gio_grid_blocks[0].grid_neighbors[k])
    REG_IPOLE("cell_corners",         gio_grid_blocks[0].grid_corners[k])
    REG_IPOLE("cell_edges",           gio_grid_blocks[0].grid_edges[k])

    gio_param->total_time_in_API += MPI_Wtime() - starttime;
}

/*----< gio_register_restart_dbl() >-----------------------------------------*/
/* Register data for use in restart file */
void gio_register_restart_dbl(gio_data   *data,
                              const char *t_name,
                              double      t_data)
{
    int idx = data->num_restart_dbl;
    strcpy(data->restart_dbl_names[idx], t_name);
    data->restart_dbl[idx] = t_data;
    data->num_restart_dbl++;
}

/*----< gio_register_restart_int() >-----------------------------------------*/
void gio_register_restart_int(gio_data   *data,
                              const char *t_name,
                              int         t_data)
{
    int idx = data->num_restart_int;
    strcpy(data->restart_int_names[idx], t_name);
    data->restart_int[idx] = t_data;
    data->num_restart_int++;
}

/*----< gio_check_file_data() >-----------------------------------------------*/
/* Check to make sure that all data that should be written out by a file
 * is registered */
void gio_check_file_data(gio_parameters *gio_param)
{
    int i, j, k, ifld, i_ok, r_ok, icnt, is_restart, rcv_rst, ok;

    file_descriptor *gio_files = gio_param->file.gio_files;
    data_descriptor *gio_descriptors = gio_param->data.gio_descriptors;

    for (i=0; i<gio_param->file.num_files; i++) {
        ok = 1;
        is_restart = 0;
        icnt = 0;
        for (j=0; j<gio_files[i].nflds; j++) {
            ifld = gio_files[i].fields[j];
            if (!gio_descriptors[ifld].assigned &&
                strncmp(gio_descriptors[ifld].field_name, "time", 4))
                ok = 0;

            if (gio_descriptors[ifld].assigned ||
                !strncmp(gio_descriptors[ifld].field_name, "time", 4)) {
                gio_files[i].fields[icnt] = ifld;
                for (k=0; k<MAX_DIMS; k++)
                    gio_files[i].field_dim_ids[k][icnt] =
                    gio_files[i].field_dim_ids[k][j];
                gio_files[i].is_grid[icnt] = gio_files[i].is_grid[j];
                gio_files[i].is_time_field[icnt] =
                gio_files[i].is_time_field[j];
                icnt++;
            }
            if (gio_descriptors[ifld].is_restart) is_restart = 1;
        }

        gio_files[i].nflds = icnt;
        if (ok) i_ok = 1;
        else    i_ok = 0;

        MPI_Allreduce(&i_ok, &r_ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        MPI_Allreduce(&is_restart, &rcv_rst, 1, MPI_INT, MPI_PROD,
                      MPI_COMM_WORLD);
        if (r_ok > 0)
            gio_files[i].complete = 1;
        else {
            if (rcv_rst == 0)
                gio_files[i].complete = 0;
            else
                gio_files[i].complete = 1;

            if (gio_param->gio_me == 0 &&
                strcmp(gio_files[i].file_prefix, "restart.nc"))
                printf("Unregistered data for file %s\n",
                       gio_files[i].file_prefix);
        }
    }
}
