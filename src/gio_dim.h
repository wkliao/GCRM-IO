/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio_dim.h 4599 2017-12-07 07:02:13Z wkliao $
 */

#ifndef H_GIO_DIM
#define H_GIO_DIM

#include <gio.h>

typedef struct {
    int *p;
} index_ptr;

typedef struct {
    int p;
} index_val_ptr;

typedef struct {
    /* Indexed list of dimension names and sizes */
    int  num_file_dims;
    int  num_index_arrays;
    int  num_species_dims;
    char file_dim_names[MAX_DIMS][32];
    int  file_dim_sizes[MAX_DIMS];

    int  index_sizes[MAX_INDICES];
    char index_names[MAX_INDICES][32];
    int  species_sizes[MAX_INDICES];
    char species_names[MAX_INDICES][32];

    index_ptr     index_arrays[MAX_INDICES];
    index_val_ptr   index_vals[MAX_INDICES];

    int  slice_Z, slice_comp;
    int  index_choice;

} gio_dim;

#endif
