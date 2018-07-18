/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: grid_params.h 4603 2017-12-07 07:12:52Z wkliao $
 */

#ifndef H_GRID_PARAMS
#define H_GRID_PARAMS

#if HAVE_CONFIG_H
#include "config.h"
#endif

/*
   MODULE grid_params
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Purpose:
!
!     This module specifies parameters related to the horizontal structure
!     of the model and horizontal grid of the atmosphere.
!
!  Define:
!
!    level_max -> resolution of the grid. This parameter determines the
!                 global horizontal grid resolution. TABLE 1 shows a few.
!                 examples. User can change level_max to change the grid
!                 resolution and hence problem size.
!
!                           ____________________________________
!                           |              number   resolution |  
!                           | level_max   of cells     (km)    |  
!                           | ---------------------------------|
!                           |     5           10242   250.2    | 
!                           |     6           40962   125.1    |  
!                           |     7          163842    62.55   | 
!                           |     8          655362    31.27   | 
!                           |     9         2621442    15.64   |
!                           |    10        10485762     7.82   |
!                           |    11        41943042     3.91   |
!                           |    12       167772162     1.95   |
!                           |    13       671088642     0.977  |
!                           |    14      2684354562     0.487  |
!                           |    15     10737418242     0.244  |
!                 TABLE 1.  |    16     42949672962     0.122  |
!                           ------------------------------------
!
!     sbdmn_iota -> determines the global number of subdomain blocks that
!                   constitute the horizontal domain decomposition. Table 2
!                   shows a few examples. User can change this number to set
!                   the global number of subdomain blocks.
!
!                           _____________________________ 
!                           |             global number |
!                           |             of subdomains |
!                           | sbdmn_iota   (=nsdm_glbl) |
!                           | --------------------------|
!                           |      0             10     |
!                           |      1             40     |
!                           |      2            160     |
!                           |      3            640     |
!                           |      4           2560     |
!                           |      5          10240     |
!                           |      6          40960     |
!                 TABLE 2.  |      7         163840     |
!                           -----------------------------
!
!     level_glbl        -> maximum level to which the tree grid grid data
!                          stucture is globally generated.  at higher 
!                          resolutions the local process generates only its 
!                          portion of the grid.
!
!     cell_max          -> the global number of cells on the finest resolution
!
!     nsdm_glbl         -> the global number of subdomains on the finest
!                          resolution
!
!     nsdm              -> the number of subdomain blocks managed by the local 
!                          process
!
!     im, jm            -> the extent of local arrays in the i-direction
!                          and j-direction, respectively
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#define RNK_NONEXISTENT -1

typedef struct grid_node grid_node;
struct grid_node {
    int l_real;
    int l_ghst;
    int l_north;
    int l_south;
    int l_pentagon;
    int l_pole_north;
    int l_pole_south;
    int l_pentagon_north;
    int l_pentagon_south;
    int l_point_set;

    int level;
    int tag_glbl;
    int tag_locl;
    int tag_sprl;
    int i,j;
    int ix_2D[3];
    int proc;

    int **nghbr_lst;
    double point[3];
    double corner[6][3];
    double area,area_inv;
    double area_crn[2];
    double area_inv_crn[2];
    double f;

    grid_node *nghbr[7]; /* [0:6] */
    grid_node *up;
    grid_node *dn[4];    /* [0:3] */
    grid_node *next_next;
    grid_node *next_real;
    grid_node *next_ghst;
    grid_node *next_sprl;
};

typedef struct extended_list_node extended_list_node;
struct extended_list_node {
    int                 nsd_glbl;  /* global block ID, 0-based */
    extended_list_node *next;
};

typedef struct {
    int l_agent_north;
    int l_agent_south; /* local process manages pole grid point */
    int cell_max;
    int sbdmn_iota;
    int im;
    int jm;
    int nsdm_glbl;       /* number of global blocks */
    int nsdm;            /* number of local  blocks */
    int nsd_north_glbl;  /* global north pole block ID, 0-based */
    int nsd_south_glbl;  /* global south pole block ID, 0-based */
    int nsd_north;       /* local  north pole block ID, 0-based */
    int nsd_south;       /* local  south pole block ID, 0-based */
    int *lst;            /* [nsdm] block IDs, 0-based */
    int *proc;
    extended_list_node *extended_list_head;
} sbdmn_node;

typedef struct {

    int       level_max;    /* grid resolution level (r value) */
    int       sbdmn_iota;   /* number of subdomain blocks of the horizontal
                               domain decomposition */
    int       level_glbl;
    int       cell_max;
    int       nsdm_glbl, nsdm;
    int       im, jm;

    int       edgm, crnm;
    int       tag_nonexistent;

    long long grid_node_total; /* the total number of grid nodes */
    long long grid_node_size;  /* allocated by the local process */
    long long grid_node_memory_max; /* size in bytes of one grid node */

    int       l_agent_north; /* local process manages pole grid point */
    int       l_agent_south;
    int       nsd_north;
    int       nsd_south;

    double    alfalfa;
    double    pentagon_latitude;

    grid_node grid[12];

    // sbdmn_node sbdmn[level_max+1]; Fortran sbdmn(0:level_max)
    sbdmn_node *sbdmn;

    grid_node **path_next; /* [level_max+1] */
    grid_node **path_real; /* [level_max+1] */
    grid_node **path_ghst; /* [level_max+1] */

    int nm;
    int *nm_lvl;      /* [level_max+1] */
    int nm_real;
    int *nm_real_lvl; /* [level_max+1] */
    int nm_ghst;
    int *nm_ghst_lvl; /* [level_max+1] */

} MODULE_grid_params;

void init_MODULE_grid_params(MODULE_grid_params *grid_params, const char *in_fname);
void finalize_MODULE_grid_params(MODULE_grid_params *grid_params);
grid_node* set_ptr(grid_node *grid, int level_max, int level, int tag_glbl);
void initialize_grid_connectivity(MODULE_grid_params *grid_params,
                                  char               *communicator_name);
void finalize_grid_connectivity(MODULE_grid_params *grid_params);

void get_index(grid_node *grid,
               int        level_max,
               int        level,
               int        iota,
               int       *subdomain_list,
               int        subdomain_list_len,
               int        tag_glbl,
               int       *ix);

#endif
