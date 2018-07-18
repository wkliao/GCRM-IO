/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio_file.h 4609 2017-12-07 07:26:38Z wkliao $
 */

#ifndef H_GIO_FILE
#define H_GIO_FILE

#include <gio.h>

typedef struct {
    char file_prefix[128];
    char directory_name[128];
    char last_time_str[24];    /* Used to keep track of filenames we have used for this file */
    int  frequency;            /* in seconds */
    int  grid_reduction;
    int  grid_3D_to_2D;
    int  grid_ONLY;
    int  contain_grids;        /* whether this file contains grid variables */
    int  fields[MAX_DATA_FIELDS];  /* indices into gio_descriptors structure */

/* wkl: TODO switch the 2D dimensions from Fortran to C */
    int  field_dim_ids[MAX_DIMS][MAX_DATA_FIELDS]; /* IDs of dimensions for each data field */
    int  is_grid[MAX_DATA_FIELDS];
    int  is_time_field[MAX_DATA_FIELDS]; /* Time field is handled in a different manner */
    int  nflds;               /* Number of fields for this file type */
    int  fileID;              /* Unit number returned by ncmpi_create() or ncmpi_open() ... ordinary integer */
    int  nsamples;            /* Max. number of time samples to put in this file */
    int  samples_written;     /* To keep track of the number of samples we actually write */
    int  samples_read;        /* To keep track of the number of samples we actually read */

    char dimension_names[MAX_DIMS][32];
    int  dimension_sizes[MAX_DIMS];
    int  dimension_ids[MAX_DIMS];
    int  hdims[MAX_DIMS];     /* true for horizontal(distributed) dims */

    int  dimension_num;
    int  complete;             /* all data associated with this file has been registered */
    int  grid_level;           /* grid resolution of all fields in file */
    int  crsn_average;         /* average values instead of decimating if grid is coarsened */

    int is_clustered_file;     /* clustered or full */
} file_descriptor;

typedef struct {
    /* This module contains all the data descriptors for the system */

    /*  number of files defined in configuration file */
    int num_files;
    int restart_idx;

    /*  Some user defined restart parameters */
    double gio_restart_interval;
    char   gio_restart_fname[512];
    int    gio_restart_overwrite;

    /*  Vector of data descriptors */
    file_descriptor *gio_files; /* [MAX_FILES] */

} gio_file;

#endif
