/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio_grid.h 4609 2017-12-07 07:26:38Z wkliao $
 */

#ifndef H_GIO_GRID
#define H_GIO_GRID

#include <gio.h>

/* grid block descriptor */
typedef struct {
    double   *grid_center_lat;
    double   *grid_center_lon;
    double  **grid_corner_lat;
    double  **grid_corner_lon;
    double  **grid_edge_lat;
    double  **grid_edge_lon;
    double   *grid_area;
    double  **crnrarea;
    double  **latedgarea;
    double  **lonedgarea;
    int      *tags;
    int     **grid_neighbors;
    int     **grid_corners;
    int     **grid_edges;
    int    ***grid_endpoints;

    /* what portion of geodesic grid does this block correspond to? */
    int       imin, imax, jmin, jmax;
    int       block_index;
    int       panel_id;
} grid_block;

typedef struct {
    /* This module contains information about the portion of grid
     *  owned by this processor
     */

    /* Does this process contain data for one or both poles? */
    int   north_pole_flag;
    int   south_pole_flag;

    /* General characteristics of grid */
    int   nlevel;
    int   refine_param;
    int   grid_length;
    int   block_size;
    int   total_grid_blocks;
    int   gio_grid_size;
    int   gio_corners_size; /* Total number of unique corners (2 * (gio_grid_size-2)) */
    int   gio_edges_size;   /* Total number of unique edges   (3 * (gio_grid_size-3)) */
    int   gio_istride;
    int   gio_jstride;
    int **grid_map;

    /*  number of grid blocks held by processor */
    int   num_grid_blocks;
    grid_block *gio_grid_blocks; /* [MAX_DATA_BLOCKS] */

} gio_grid;

#endif
