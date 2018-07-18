/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: dimensions.c 4606 2017-12-07 07:21:06Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "gcrm.h"
#include "gio.h"

/*----< gio_init_dims() >-----------------------------------------------------*/
/* Initialize dimension parsing data structures */
int gio_init_dims(gio_dim *dim)
{
    int i;
    dim->num_file_dims    = 0;
    dim->num_index_arrays = 0;
    dim->num_species_dims = 0;

    for (i=0; i<MAX_DIMS; i++) {
      dim->file_dim_names[i][0] = '\0';
      dim->file_dim_sizes[i] = 0;
    }
    return 1;
}

/*--< gio_parse_dimension() >-------------------------------------------------*/
/* High level routine to parse dimension line from descriptor file */
void gio_parse_dimension(gio_parameters *gio_param,
                         int             idx,
                         char           *dims_str)
{
    int i;
    data_descriptor *gio_descriptors = gio_param->data.gio_descriptors;

    /* Initialize descriptor */
    gio_descriptors[idx].num_dims = 0;
    for (i=0; i<MAX_DIMS; i++)
        gio_descriptors[idx].field_dims[i] = 0;

    // gio_set_dimension_flags(gio_param, idx, dims_str);
    // gio_set_dimension_limits(gio_param, idx);
}

/*----< gio_add_file_dim_name() >--------------------------------------------*/
/* Add dimension name to global list of possible dimensions */
static
void gio_add_file_dim_name(gio_dim *dim,
                           char    *name)
{
    int i;

    for (i=0; i<dim->num_file_dims; i++)
        if (!strcmp(name, dim->file_dim_names[i]))
            break;

    if (i == dim->num_file_dims) {
        strcpy(dim->file_dim_names[dim->num_file_dims], name);
        dim->num_file_dims++;
        if (dim->num_file_dims > MAX_DIMS) {
            printf("Number of file dimensions %d exceeds MAX_DIMS\n",
                   dim->num_file_dims);
            ABORT
        }
    }
}

/*----< gio_add_file_dim_size() >---------------------------------------------*/
/* Add dimension size to global list of possible dimensions */
void gio_add_file_dim_size(gio_dim *dim,
                           char    *name,
                           int      size)
{
    int i;
    for (i=0; i<dim->num_file_dims; i++)
        if (!strcmp(name, dim->file_dim_names[i]))
            break;
    if (i < dim->num_file_dims)
        dim->file_dim_sizes[i] = size;
}

#define GIO_ADD_FILE_DIM(dim_str, dim_value, dim_idx) {                       \
    int j;                                                                    \
    descriptors[idx].field_dims[i] = dim_value;                               \
    for (j=0; j<dim->num_file_dims; j++) {                                    \
        if (!strcasecmp(dim->file_dim_names[j], dim_str)) break;              \
    }                                                                         \
    if (j == dim->num_file_dims) {                                            \
        /* new file dimension name is found */                                \
        strcpy(dim->file_dim_names[dim->num_file_dims], dim_str);             \
        dim->file_dim_sizes[dim->num_file_dims] = dim_value;                  \
        dim->num_file_dims++;                                                 \
        if (dim->num_file_dims > MAX_DIMS) {                                  \
            printf("Error: Number of file dimensions exceeds MAX_DIMS\n");    \
            ABORT                                                             \
        }                                                                     \
    }                                                                         \
    strcpy(descriptors[idx].raw_dim_names[dim_idx-1], dim_str);               \
    descriptors[idx].raw_dim_sizes[dim_idx-1] = dim_value;                    \
    gio_add_file_dim_name(dim, dim_str);                                      \
    gio_add_file_dim_size(dim, dim_str, dim_value);                           \
}

/*----< gio_set_dimension() >-----------------------------------------------*/
/* This combines gio_set_dimension_flags() and gio_set_dimension_limits() */
void gio_set_dimension(int        idx,
                       int        num_dims,
                       char     **dims_str, /* [num_dims] */
                       gio_data  *data,
                       gio_grid  *grid,     /* IN */
                       gio_dim   *dim)
{
    int i, num_unknown_indices=0;
    data_descriptor *descriptors = data->gio_descriptors;

    for (i=0; i<MAX_DIMS; i++)
        descriptors[idx].field_dims[i] = 0;

    /* Initialize descriptor (note descriptors has been zero-ed) */
    for (i=0; i<7; i++) {
        descriptors[idx].raw_dim_sizes[i] = 1;
        descriptors[idx].raw_dim_names[i][0] = '\0';
    }

    descriptors[idx].num_dims   = num_dims;

    descriptors[idx].imin1      = 1;
    descriptors[idx].imin2      = 1;
    descriptors[idx].imin3      = 1;
    descriptors[idx].imin4      = 1;
    descriptors[idx].imax1      = 1;
    descriptors[idx].imax2      = 1;
    descriptors[idx].imax3      = 1;
    descriptors[idx].imax4      = 1;
    descriptors[idx].inc1       = 1;
    descriptors[idx].inc2       = 1;
    descriptors[idx].inc3       = 1;
    descriptors[idx].ld1        = 1;
    descriptors[idx].ld2        = 1;
    descriptors[idx].ld3        = 1;
    descriptors[idx].no_poles   = 0;
    descriptors[idx].is_vector  = 0;
    descriptors[idx].has_endpts = 0;
    descriptors[idx].is_static  = 1;
    descriptors[idx].use_layers = 0;
    descriptors[idx].is_corner  = 0;
    descriptors[idx].is_general = 1;
    descriptors[idx].grid_size  = grid->gio_grid_size;

    for (i=0; i<num_dims; i++) {
        if (!strcasecmp(dims_str[i], "cells")) {
            GIO_ADD_FILE_DIM("cells", data->grid_cells_dim, 4)
        }
        else if (!strcasecmp(dims_str[i], "time")) {
            descriptors[idx].is_static = 0;
            GIO_ADD_FILE_DIM("time", data->time_dim, 7)
        }
        else if (!strcasecmp(dims_str[i], "vector")) {
            descriptors[idx].imin2     = 1;
            descriptors[idx].imax2     = 3;
            descriptors[idx].inc2      = 1;
            descriptors[idx].ld2       = 3;
            descriptors[idx].is_vector = 1;
            GIO_ADD_FILE_DIM("vector", data->grid_vectors_dim, 1)
        }
        else if (!strcasecmp(dims_str[i], "i_uniq_edge")) {
            descriptors[idx].imin1     = 1;
            descriptors[idx].imax1     = 3;
            descriptors[idx].inc1      = 1;
            descriptors[idx].ld1       = 3;
            descriptors[idx].no_poles  = 1;
            descriptors[idx].grid_size = grid->gio_edges_size;
            GIO_ADD_FILE_DIM("edges", data->grid_all_uniq_edge_dim, 4)
        }
        else if (!strcasecmp(dims_str[i], "i_uniq_corner")) {
            descriptors[idx].imin1     = 1;
            descriptors[idx].imax1     = 2;
            descriptors[idx].inc1      = 1;
            descriptors[idx].ld1       = 2;
            descriptors[idx].no_poles  = 1;
            descriptors[idx].grid_size = grid->gio_corners_size;
            GIO_ADD_FILE_DIM("corners", data->grid_all_uniq_crnr_dim, 4)
        }
        else if (!strcasecmp(dims_str[i], "cellneighbors")) {
            descriptors[idx].imin1     = 1;
            descriptors[idx].imax1     = 6;
            descriptors[idx].inc1      = 1;
            descriptors[idx].ld1       = 6;
            descriptors[idx].is_corner = 1;
            GIO_ADD_FILE_DIM("cellneighbors", data->grid_cellneighbors_dim, 2)
        }
        else if (!strcasecmp(dims_str[i], "celledges")) {
            descriptors[idx].imin1     = 1;
            descriptors[idx].imax1     = 6;
            descriptors[idx].inc1      = 1;
            descriptors[idx].ld1       = 6;
            descriptors[idx].is_corner = 1;
            GIO_ADD_FILE_DIM("celledges", data->grid_celledges_dim, 2)
        }
        else if (!strcasecmp(dims_str[i], "cellcorners")) {
            descriptors[idx].imin1     = 1;
            descriptors[idx].imax1     = 6;
            descriptors[idx].inc1      = 1;
            descriptors[idx].ld1       = 6;
            descriptors[idx].is_corner = 1;
            GIO_ADD_FILE_DIM("cellcorners", data->grid_cellcorners_dim, 2)
        }
        else if (!strcasecmp(dims_str[i], "i_cell_corners")) {
            descriptors[idx].imin1     = 1;
            descriptors[idx].imax1     = 2;
            descriptors[idx].inc1      = 1;
            descriptors[idx].ld1       = 6;
            descriptors[idx].no_poles  = 1;
            descriptors[idx].grid_size = grid->gio_corners_size;
            GIO_ADD_FILE_DIM("corners", data->grid_all_uniq_crnr_dim, 4)
        }
        else if (!strcasecmp(dims_str[i], "i_cell_edges")) {
            descriptors[idx].imin1     = 1;
            descriptors[idx].imax1     = 3;
            descriptors[idx].inc1      = 1;
            descriptors[idx].ld1       = 3;
            descriptors[idx].no_poles  = 1;
            descriptors[idx].grid_size = grid->gio_edges_size;
            GIO_ADD_FILE_DIM("edges", data->grid_all_uniq_edge_dim, 4)
        }
        else if (!strcasecmp(dims_str[i], "ns_ew")) {
            descriptors[idx].imin2      = 1;
            descriptors[idx].imax2      = 2;
            descriptors[idx].inc2       = 1;
            descriptors[idx].ld2        = 2;
            descriptors[idx].has_endpts = 1;
            GIO_ADD_FILE_DIM("ns_ew", 2, 1)
        }
        else if (!strcasecmp(dims_str[i], "endpoints")) {
            descriptors[idx].imin2      = 1;
            descriptors[idx].imax2      = 2;
            descriptors[idx].inc2       = 1;
            descriptors[idx].ld2        = 2;
            descriptors[idx].has_endpts = 1;
            GIO_ADD_FILE_DIM("endpoints", data->grid_endpoints_dim, 1)
        }
        else if (!strcasecmp(dims_str[i], "layers")) {
            descriptors[idx].imin3      = 1;
            descriptors[idx].imax3      = grid->nlevel;
            descriptors[idx].inc3       = 1;
            descriptors[idx].ld3        = grid->nlevel;
            descriptors[idx].has_levels = 1;
            descriptors[idx].use_layers = 1;
            GIO_ADD_FILE_DIM("layers", data->grid_layers_dim, 3)
        }
        else if (!strcasecmp(dims_str[i], "interfaces")) {
            descriptors[idx].imin3      = 1;
            descriptors[idx].imax3      = grid->nlevel + 1;
            descriptors[idx].inc3       = 1;
            descriptors[idx].ld3        = grid->nlevel + 1;
            descriptors[idx].has_levels = 1;
            descriptors[idx].use_layers = 0;
            GIO_ADD_FILE_DIM("interfaces", data->grid_interfaces_dim, 3)
        }
        else {
            /* Check to see if dimension corresponds to an index array */
            strcpy(descriptors[idx].index_dimension_name[num_unknown_indices],
                   dims_str[i]);
            strcpy(descriptors[idx].raw_dim_names[4+num_unknown_indices],
                   dims_str[i]);
            descriptors[idx].num_unknown_dims = ++num_unknown_indices;
        }
    }

    /* Evaluate maximum data size */
    MPI_Offset size_check = 1;
    for (i=0; i<num_dims; i++) {
        /* Exclude -1 dimension which implies time */
        if (descriptors[idx].field_dims[i] != -1)
           size_check *= descriptors[idx].field_dims[i];
    }

    int bytes = 4;
    if (descriptors[idx].data_type == dbl_type) bytes = 8;
    data->field_size_check = MAX(data->field_size_check, size_check*bytes);
}

/*----< gio_set_file_dims_for_fields() >-------------------------------------*/
/* Determine dimensions that will appear in file once file descriptor values
 * have all been set */
void gio_set_file_dims_for_fields(gio_parameters *gio_param)
{
    int j, idx, icnt, imin, imax, inc;

    gio_data *data = &gio_param->data;
    data_descriptor *gio_descriptors = gio_param->data.gio_descriptors;

    /* loop over descriptors */
    for (idx=0; idx<data->num_fields; idx++) {
        if (gio_descriptors[idx].no_data) continue;

        /* Find out how large dimensions actually are */

        if (gio_descriptors[idx].raw_dim_names[2][0] != '\0') {
            icnt = 0;
            imin = gio_descriptors[idx].imin3;
            imax = gio_descriptors[idx].imax3;
            inc = gio_descriptors[idx].inc3;
            for (j=imin; j<=imax; j+=inc)
                icnt++;
            gio_descriptors[idx].raw_dim_sizes[2] = icnt;
        }
        if (gio_descriptors[idx].raw_dim_names[1][0] != '\0') {
            icnt = 0;
            imin = gio_descriptors[idx].imin1;
            imax = gio_descriptors[idx].imax1;
            inc = gio_descriptors[idx].inc1;
            for (j=imin; j<=imax; j+=inc)
                icnt++;
            gio_descriptors[idx].raw_dim_sizes[1] = icnt;
        }
        if (gio_descriptors[idx].raw_dim_names[0][0] != '\0') {
            icnt = 0;
            imin = gio_descriptors[idx].imin2;
            imax = gio_descriptors[idx].imax2;
            inc = gio_descriptors[idx].inc2;
            for (j=imin; j<=imax; j+=inc)
                icnt++;
            gio_descriptors[idx].raw_dim_sizes[0] = icnt;
        }

        /* Temporarily set values for time integration and species indices */
        if (gio_descriptors[idx].raw_dim_names[4][0] != '\0') {
            if (gio_descriptors[idx].index_cur != -1)
              gio_descriptors[idx].raw_dim_sizes[4] = 1;
            else
              gio_descriptors[idx].raw_dim_sizes[4] = -1;
        }
        if (gio_descriptors[idx].raw_dim_names[5][0] != '\0')
          gio_descriptors[idx].raw_dim_sizes[5] = 1;
    }
}

/*----< gio_compress_file_dims_for_fields() >---------------------------------*/
/* Finish evaluating unknown dimension sizes and eliminate any
 * dimensions of size 1 from descriptor */
static
void gio_compress_file_dims_for_fields(gio_parameters *gio_param)
{
    int idx, icnt, idim, unk_idx;
    int has_time_idx;

    gio_dim      *dim      = &gio_param->dim;
    gio_data     *data     = &gio_param->data;
    data_descriptor *gio_descriptors = gio_param->data.gio_descriptors;

    /* loop over descriptors */
    for (idx=0; idx<data->num_fields; idx++) {
        /* Check to see if time integration index is present */
        has_time_idx = 0;
        if (gio_descriptors[idx].raw_dim_names[4][0] != '\0') {
            for (idim=0; idim<dim->num_index_arrays; idim++) {
                int slen=strlen(dim->index_names[idim]);
                if (!strncmp(gio_descriptors[idx].raw_dim_names[4],
                             dim->index_names[idim], slen)) {
                    has_time_idx = 1;
                    break;
                }
            }
        }

        /* Check to see if species index is present */
        if (has_time_idx && gio_descriptors[idx].num_unknown_dims > 1)
            unk_idx = 2;
        else if (!has_time_idx && gio_descriptors[idx].num_unknown_dims > 0)
            unk_idx = 1;
        else
            unk_idx = 0;

        if (gio_descriptors[idx].raw_dim_names[3+unk_idx][0] != '\0') {
            for (idim=0; idim<dim->num_species_dims; idim++) {
                int slen=strlen(dim->species_names[idim]);
                if (!strncmp(gio_descriptors[idx].raw_dim_names[3+unk_idx],
                             dim->species_names[idim], slen)) {
                    has_time_idx = 1;
                    gio_descriptors[idx].raw_dim_sizes[3+unk_idx] =
                    dim->species_sizes[idim];
                }
            }
        }

        /* Clean up descriptors with no data */
        if (gio_descriptors[idx].no_data) {
            for (idim=0; idim<7; idim++)
                gio_descriptors[idx].raw_dim_names[idim][0] = '\0';
        }
        icnt = 0;
        gio_descriptors[idx].cell_idx = -1;
        gio_descriptors[idx].time_idx = -1;
        for (idim=0; idim<7; idim++) {
            if (gio_descriptors[idx].raw_dim_names[idim][0] != '\0' &&
                gio_descriptors[idx].raw_dim_sizes[idim] != 1) {
                strcpy(gio_descriptors[idx].file_dim_names[icnt],
                       gio_descriptors[idx].raw_dim_names[idim]);
                gio_descriptors[idx].file_dims[icnt] =
                gio_descriptors[idx].raw_dim_sizes[idim];
                if (idim == 3)
                    gio_descriptors[idx].cell_idx = icnt;
                if (idim == 6)
                    gio_descriptors[idx].time_idx = icnt;
                if (gio_descriptors[idx].intg_idx == idim)
                    gio_descriptors[idx].intg_idx = icnt;
                icnt++;
            }
        }
        gio_descriptors[idx].num_file_dims = icnt;
    }
}

/*----< gio_set_file_header_dimensions() >------------------------------------*/
/* Set up dimension fields that appear at top of each file */
void  gio_set_file_header_dimensions(gio_parameters *gio_param)
{
    int i, ifile, ifld, idim, idesc, dimension_num;

    gio_data     *data     = &gio_param->data;
    gio_file     *file     = &gio_param->file;
    data_descriptor *gio_descriptors = data->gio_descriptors;
    file_descriptor *gio_files = file->gio_files;

    for (ifile=0; ifile<file->num_files; ifile++) {
        dimension_num = 0;
        gio_files[ifile].dimension_num = 0;
        for (ifld=0; ifld<gio_files[ifile].nflds; ifld++) {
            idesc = gio_files[ifile].fields[ifld];
            for (i=0; i<MAX_DIMS; i++)
                gio_files[ifile].field_dim_ids[i][ifld] = 0;

            for (idim=0; idim<gio_descriptors[idesc].num_file_dims; idim++) {
                for (i=0; i<gio_files[ifile].dimension_num; i++)
                    if (!strcmp(gio_descriptors[idesc].file_dim_names[idim],
                                gio_files[ifile].dimension_names[i]))
                        break;
                if (i < gio_files[ifile].dimension_num)
                    gio_files[ifile].field_dim_ids[idim][ifld] = i;
                else {
                    strcpy(gio_files[ifile].dimension_names[dimension_num],
                           gio_descriptors[idesc].file_dim_names[idim]);
                    gio_files[ifile].dimension_sizes[dimension_num] =
                      gio_descriptors[idesc].file_dims[idim];
                    gio_files[ifile].field_dim_ids[idim][ifld] = dimension_num;

                    gio_files[ifile].hdims[dimension_num] = 0;
                    if (gio_descriptors[idesc].cell_idx == idim)
                       gio_files[ifile].hdims[dimension_num] = 1;

                    gio_files[ifile].dimension_num = ++dimension_num;
                    if (dimension_num > MAX_DIMS) {
                        printf("Number of file dimensions %d exceeds MAX_DIMS\n",
                               dimension_num);
                        ABORT
                    }
                }
            }
        }
    }
}

/*----< gio_register_index() >-----------------------------------------------*/
/* Register arrays that will be used as dimension variables */
void gio_register_index(gio_dim *dim,
                        char    *t_index_name,
                        int     *iarray,
                        int      ival,
                        int      length)
{
    int idx;

    idx = dim->num_index_arrays;
    strcpy(dim->index_names[idx], t_index_name);
    dim->index_sizes[idx] = length;
    dim->index_arrays[idx].p = iarray;
    dim->index_vals[idx].p = ival;
    dim->num_index_arrays++;
}

/*----< gio_set_limits() >---------------------------------------------------*/
/*  Figure out what limits should be used for messaging loop */
void gio_set_limits(gio_parameters *gio_param,
                    int             idx,
                    int            *imin,
                    int            *imax)
{
    int i, index_size, index_choice, index_cur;

    index_choice = gio_param->data.gio_descriptors[idx].index_choice;
    index_size   = gio_param->dim.index_sizes[index_choice];
    index_cur    = gio_param->data.gio_descriptors[idx].index_cur;

    if (index_choice >= 0) {
        if (index_cur == -1) {
            *imin = 0;
            *imax = index_size-1;
        }
        else {
            for (i=0; i<=index_choice; i++) {
                if (gio_param->dim.index_arrays[index_choice].p[i] == index_cur) {
                    *imin = i;
                    *imax = i;
                    break;
                }
            }
            if (i == index_choice+1) {
                printf("Error: could not find index_cur from dim.index_arrays[]\n");
                ABORT
            }
        }
    }
    else {
        *imin = 0;
        *imax = 0;
    }
}

/*----< gio_set_index_dim_sizes() >-------------------------------------------*/
/* Set dimension sizes for time integration indices, if necessary */
void gio_set_index_dim_sizes(gio_parameters *gio_param)
{
    int idx, jdx, kdx, is_intg_idx;

    gio_data     *data     = &gio_param->data;
    gio_dim      *dim      = &gio_param->dim;
    gio_file     *file     = &gio_param->file;

    data_descriptor *gio_descriptors = gio_param->data.gio_descriptors;

    for (idx=0; idx<data->num_fields; idx++) {
        gio_descriptors[idx].intg_idx = -1;
        if (gio_descriptors[idx].num_unknown_dims > 0) {

            /* Check to if unknown dimension matches one of the time
               integration fields */
            for (jdx=0; jdx<dim->num_index_arrays; jdx++)
                if (!strcmp(gio_descriptors[idx].index_dimension_name[0],
                            dim->index_names[jdx]))
                    break;
            is_intg_idx = (jdx < dim->num_index_arrays) ? 1 : 0;

            /* Find location of time integration index */

            if (is_intg_idx && gio_descriptors[idx].index_cur == -1) {
                for (jdx=0; jdx<gio_descriptors[idx].num_file_dims; jdx++) {
                    if (!strcmp(gio_descriptors[idx].index_dimension_name[0],
                                gio_descriptors[idx].file_dim_names[jdx])) {
                        for (kdx=0; kdx<dim->num_index_arrays; kdx++)
                            if (!strcmp(gio_descriptors[idx].index_dimension_name[0],
                                        dim->index_names[kdx]))
                                break;
                        gio_descriptors[idx].intg_idx = jdx;
                        gio_descriptors[idx].file_dims[jdx] = dim->index_sizes[kdx];
                    }
                }
            }
        }
    }
    for (idx=0; idx<file->num_files; idx++) {
        for (jdx=0; jdx<file->gio_files[idx].dimension_num; jdx++) {
            if (file->gio_files[idx].dimension_sizes[jdx] == -1) {
                for (kdx=0; kdx<dim->num_index_arrays; kdx++)
                    if (!strcmp(file->gio_files[idx].dimension_names[jdx],
                                dim->index_names[kdx]))
                        break;
                file->gio_files[idx].dimension_sizes[jdx] = dim->index_sizes[kdx];
            }
        }
    }
}


/*----< gio_set_unknown_dim_sizes() >----------------------------------------*/
/* Figure out sizes of unknown dimensions */
static
void gio_set_unknown_dim_sizes(gio_parameters *gio_param)
{
    int idx, j, found_time_idx, found_species_idx;

    gio_data *data = &gio_param->data;
    gio_dim  *dim  = &gio_param->dim;

    data_descriptor *gio_descriptors = gio_param->data.gio_descriptors;

    for (idx=0; idx<data->num_fields; idx++) {
        /* Check to see if an unknown dimension is being used and adjust
           descriptor properties accordingly */

        if (gio_descriptors[idx].num_unknown_dims > 0) {
            found_time_idx = 0;
            for (j=0; j<dim->num_index_arrays; j++)
                if (!strcmp(gio_descriptors[idx].index_dimension_name[0],
                            dim->index_names[j]))
                    break;

            if (j < dim->num_index_arrays) {
                dim->index_choice = j;
                found_time_idx = 1;

                strcpy(gio_descriptors[idx].raw_dim_names[4],
                       dim->index_names[dim->index_choice]);
                if (gio_descriptors[idx].index_cur == -1)
                    gio_descriptors[idx].raw_dim_sizes[4] =
                    dim->index_sizes[dim->index_choice];
                else {
                    gio_descriptors[idx].raw_dim_sizes[4] = 1;
                    gio_descriptors[idx].index_cur = 0;
                }
                gio_descriptors[idx].index_choice = dim->index_choice;
                gio_add_file_dim_name(dim, dim->index_names[dim->index_choice]);
                gio_add_file_dim_size(dim, dim->index_names[dim->index_choice],
                                      dim->index_sizes[dim->index_choice]);
            }

            int dim_idx;
            if (found_time_idx && gio_descriptors[idx].num_unknown_dims > 1)
                dim_idx = 1;
            else if (! found_time_idx)
                dim_idx = 0;
            else
                dim_idx = -1;
    
            found_species_idx = 0;
            if (dim_idx >= 0) {
                for (j=0; j<dim->num_species_dims; j++)
                    if (!strcmp(gio_descriptors[j].index_dimension_name[0],
                                dim->species_names[j])) {
                        found_species_idx = 1;
                        break;
                    }
            }

            if (found_species_idx) {
                dim->index_choice = j;
                found_species_idx = 1;
                strcpy(gio_descriptors[idx].raw_dim_names[4+dim_idx],
                       dim->species_names[dim->index_choice]);
                gio_descriptors[idx].raw_dim_sizes[4+dim_idx] =
                dim->species_sizes[dim->index_choice];
                gio_add_file_dim_name(dim, dim->species_names[dim->index_choice]);
                gio_add_file_dim_size(dim, dim->species_names[dim->index_choice],
                                           dim->species_sizes[dim->index_choice]);
            }

            /*  ERROR: invalid dimension name */

            if (!found_time_idx && !found_species_idx) {
                printf("Found invalid dimension name: %s\n",
                       gio_descriptors[idx].index_dimension_name[0]);
                ABORT
            }
        }
    }
}

/*----< gio_set_strides() >---------------------------------------------------*/
/* Evaluate strides and column sizes for all data fields */
void gio_set_strides(gio_parameters *gio_param)
{
    int idx, j;
    int stride1, stride2, stride3;
    int imin1, imin2, imin3, imax1, imax2, imax3;
    int inc1, inc2, inc3;

    gio_set_unknown_dim_sizes(gio_param);
    gio_set_file_dims_for_fields(gio_param);
    gio_compress_file_dims_for_fields(gio_param);
    gio_set_index_dim_sizes(gio_param);

    gio_data *data = &gio_param->data;
    data_descriptor *gio_descriptors = gio_param->data.gio_descriptors;

    for (idx=0; idx<data->num_fields; idx++) {
        /* Set remaining descriptor properties */

        if (!gio_descriptors[idx].no_data) {
            imin1 = gio_descriptors[idx].imin1;
            imax1 = gio_descriptors[idx].imax1;
            inc1  = gio_descriptors[idx].inc1;
            imin2 = gio_descriptors[idx].imin2;
            imax2 = gio_descriptors[idx].imax2;
            inc2  = gio_descriptors[idx].inc2;
            imin3 = gio_descriptors[idx].imin3;
            imax3 = gio_descriptors[idx].imax3;
            inc3  = gio_descriptors[idx].inc3;
            if (!gio_descriptors[idx].has_levels) {
                imin3 = 1;
                imax3 = 1;
                inc3  = 1;
            }
            stride1 = 0;
            for (j=imin1; j<=imax1; j+=inc1)
                stride1++;

            stride2 = 0;
            for (j=imin2; j<=imax2; j+=inc2)
                stride2++;

            stride3 = 0;
            for (j=imin3; j<=imax3; j+=inc3)
                stride3++;
            gio_descriptors[idx].stride1 = stride1;
            gio_descriptors[idx].stride2 = stride2;
            gio_descriptors[idx].stride3 = stride3;
            gio_descriptors[idx].column_size = stride1*stride2*stride3;
        }
    }

    /* Now that field dimensions are okay, fix up file header dimensions */
    gio_set_file_header_dimensions(gio_param);
}

