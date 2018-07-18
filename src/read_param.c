/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#define __USE_GNU
#include <string.h>

#include <mpi.h>

#include "gio.h"
#include "util.h"

/*----< read_io_desc() >-----------------------------------------------------*/
/* this is corresponding to gio_read_data_config() in init_term.F90
 * Data updated:
 *     data->field_size_check
 *     data->num_grid_flds
 *     data->num_fields
 *     props->num_props
 *     props->gio_properties[props->num_props].prop_name
 *     props->gio_properties[props->num_props].prop_value
 *     descriptors[ndesc].field_name
 *     descriptors[ndesc].output_field_name
 *     descriptors[ndesc].num_dims
 *     descriptors[ndesc].is_grid
 *     descriptors[ndesc].is_restart
 *     descriptors[ndesc].no_data
 *     descriptors[ndesc].assigned
 *     descriptors[ndesc].units
 *     descriptors[ndesc].long_name
 *     descriptors[ndesc].standard_name
 *     descriptors[ndesc].data_type
 *     descriptors[ndesc].imin3
 *     descriptors[ndesc].imax3
 *     descriptors[ndesc].inc3
 *     descriptors[ndesc].imin2
 *     descriptors[ndesc].imax2
 *     descriptors[ndesc].inc2
 *     descriptors[ndesc].index_cur
 *     descriptors[ndesc].field_attr_name[num_attrs]
 *     descriptors[ndesc].field_attr_value[num_attrs]
 *     descriptors[ndesc].num_attrs
 *     gio_param->dim.slice_Z
 *     gio_param->dim.slice_comp
*/
static
int read_io_desc(gio_parameters *gio_param,
                 const char     *descfile)
{
    char line[256];
    FILE *fp;
    char *key, *value, *prop, *next_str;
    int ndesc, have_name, have_long_name, have_type, have_dimension;

    gio_prop        *props       = &gio_param->props;
    gio_data        *data        = &gio_param->data;
    data_descriptor *descriptors = data->gio_descriptors;

    /* TODO: check valid value of in_param->refine_param */
    /* TODO: check valid value of in_param->nlevel       */

    /* Keep track of the largest field, its value will be inclreaed in
       gio_set_dimension_limits(). It is used to check against 2^31 and
       decide to use format CDF, CDF-2, or CDF-5 */
    data->field_size_check = 0;

    OPEN_FILE(descfile)

    have_name      = 0;
    have_long_name = 0;
    have_type      = 0;
    have_dimension = 0;

    ndesc            = 0;
    props->num_props = 0;

    /* read data descriptors */
    line[0] = '\0';
    while (ndesc <= MAX_DATA_FIELDS) {

        if (!strncasecmp(line, "name", 4)) {
            /* A new field is detected already from the previous iteration
             * make sure the last descriptor we processed had all the
             * required elements
             */
            if (descriptors[ndesc-1].no_data) {
                if (have_name == 0 || have_long_name == 0 || have_type == 0) {
                    /* ERROR missing a required token */
                    char msg[256];
                    sprintf(msg, "IO Descriptor for %s is missing a required entry:",
                            descriptors[ndesc-1].standard_name);
                    if (have_name == 0)      strcat(msg, " name");
                    if (have_long_name == 0) strcat(msg, " long_name");
                    if (have_type == 0)      strcat(msg, " type");
                    printf("Error: %s\n", msg);
                    ABORT
                }
            }
            else {
                if (have_name == 0 || have_long_name == 0 ||
                    have_type == 0 || have_dimension == 0) {
                    /* ERROR missing a required token */
                    char msg[256];
                    sprintf(msg, "IO Descriptor for %s is missing a required entry:",
                            descriptors[ndesc-1].standard_name);
                    if (have_name == 0)      strcat(msg, " name");
                    if (have_long_name == 0) strcat(msg, " long_name");
                    if (have_type == 0)      strcat(msg, " type");
                    if (have_dimension == 0) strcat(msg, " dimension");
                    printf("Error: %s\n", msg);
                    ABORT
                }
                /* We have the dimension info, so decide how big this
                 * field will be. WE need to know the max field size
                 */
                // field_size_check = MAX(field_size_check, my_size_check*4);
                // wkl: not used ?
            }
        }
        else {
            if (fgets(line, 256, fp) == NULL) break;

            /* line contains a '\n' at the end */
            key = strtok(line, "= \t\n");
            if (key == NULL || key[0] == '#' || key[0] == '!') continue;

            if (ndesc == MAX_DATA_FIELDS) {
                printf("Error: Too many data fields specified. Increase MAX_DATA_FIELDS\n");
                ABORT
            }
        }

        /* Required Token/Values: 
         * name (becomes field_name) MUST BE FIRST token
         * units
         * type
         * long_name
         * dimension
         *
         * All other tokens, as long as they have an associated
         * value will be stored as attribute names/values
         */

        if (!strcasecmp(key, "property")) {
            if ((prop = strtok(NULL, " ")) == NULL) /* property name */
                continue;

            /* skip the leading blanks (space and tab) */
            next_str = prop + strlen(prop) + 1;
            while (*next_str == ' ' || *next_str == '\t') next_str++;

            if ((value = strtok(next_str, "\'\"\n")) == NULL) continue;
            strcpy(props->gio_properties[props->num_props].prop_name,  prop);
            strcpy(props->gio_properties[props->num_props].prop_value, value);
            props->num_props++;
        }
        else if (!strcasecmp(key, "name")) {
            /* reset required token info */
            have_name      = 1;
            have_long_name = 0;
            have_type      = 0;
            have_dimension = 0;

            if ((value = strtok(NULL, " \'\"\t\n")) == NULL)
                continue; /* name value is missing */

            strcpy(descriptors[ndesc].field_name, value);
            /* Set output_field_name to the same thing as field_name. User
             * can over-ride this via the output_field_name keyword */
            strcpy(descriptors[ndesc].output_field_name, value);

            while (fgets(line, 256, fp) != NULL) {
                key = strtok(line, "= \t\n");
                if (key == NULL || key[0] == '#' || key[0] == '!')
                    continue;
                if (!strcasecmp(key, "name")) {
                    break;
                }

                /* skip the leading blanks (space and tab) */
                next_str = key + strlen(key) + 1;
                while (*next_str == ' ' || *next_str == '\t') next_str++;

                if ((value = strtok(next_str, "\'\"\n")) == NULL) {
                    /* key without value, only isgridvar and restart allowed */
                    if (!strcasecmp(key, "isgridvar")) {
                        descriptors[ndesc].is_grid = 1;
                        data->num_grid_flds++;
                    }
                    else if (!strcasecmp(key, "restart"))
                        descriptors[ndesc].is_restart = 1;
                    else
                        printf("Warning: desc key=%s has no value\n",key);
                    continue;
                }

                if (!strcasecmp(key, "no_data")) {
                    if (value[0] == 't' || value[0] == 'T') {
                        descriptors[ndesc].no_data  = 1;
                        descriptors[ndesc].assigned = 1;
                    }
                    /* else default is false */
                }
                else if (!strcasecmp(key, "output_name"))
                    /* User wants an "alternate" field name in the output */
                    strcpy(descriptors[ndesc].output_field_name, value);
                else if (!strcasecmp(key, "units"))
                    /* units element in descriptors structure */
                    strcpy(descriptors[ndesc].units, value);
                else if (!strcasecmp(key, "long_name")) {
                    have_long_name = 1;
                    strcpy(descriptors[ndesc].long_name, value);
                }
                else if (!strcasecmp(key, "standard_name"))
                    strcpy(descriptors[ndesc].standard_name, value);
                else if (!strcasecmp(key, "type")) {
                    if (!strcasecmp(value, "integer"))
                        descriptors[ndesc].data_type = int_type;
                    else if (!strcasecmp(value, "float"))
                        descriptors[ndesc].data_type = flt_type;
                    else if (!strcasecmp(value, "double"))
                        descriptors[ndesc].data_type = dbl_type;
                    else {
                        printf("Error: invalid data type (name=%s type=%s)\n",
                               descriptors[ndesc].field_name,value);
                        ABORT
                    }
                    have_type = 1;
                }
                else if (!strcasecmp(key, "slice")) {
                    if (!have_dimension) {
                        printf("Error: Slice directive MUST come after dimension directive: (name=%s)\n",
                               descriptors[ndesc].field_name);
                        ABORT
                    }
                    key = strtok(value, " \t");
                    /* key must not be NULL */

                    /* Handle the "slice" option
                     * syntax: slice Z N
                     *         slice Z n,m,i (start at n, go to m, increment by i)
                     * or
                     * syntax: slice component N (for fields like wind_crn in GCRM
                     */
                    if (!strcasecmp(key, "z")) {
                        /* slicing the Z (layers or interfaces) component */
                        /* parse the start, end and increment */
                        char *start_val, *end_val, *inc_val;
                        start_val = strtok(NULL, " \t");
                        descriptors[ndesc].imin3 = atoi(start_val);
                        if ((end_val = strtok(NULL, " \t")) == NULL) {
                            descriptors[ndesc].imax3 = atoi(start_val);
                            descriptors[ndesc].inc3 = 1;
                        }
                        else {
                            descriptors[ndesc].imax3 = atoi(end_val);
                            if ((inc_val = strtok(NULL, " \t")) == NULL)
                                descriptors[ndesc].inc3 = 1;
                            else
                                descriptors[ndesc].inc3 = atoi(inc_val);
                        }
                        gio_param->dim.slice_Z    = 1;
                        gio_param->dim.slice_comp = 0;
                    }
                    else if (!strcasecmp(key, "component")) {
                        /* slicing the wind_crn (or similer array) */
                        char *start_val, *end_val, *inc_val;
                        start_val = strtok(NULL, " \t");
                        descriptors[ndesc].imin2 = atoi(start_val);
                        if ((end_val = strtok(NULL, " \t")) == NULL) {
                            descriptors[ndesc].imax2 = atoi(start_val);
                            descriptors[ndesc].inc2 = 1;
                        }
                        else {
                            descriptors[ndesc].imax2 = atoi(end_val);
                            if ((inc_val = strtok(NULL, " \t")) == NULL)
                                descriptors[ndesc].inc2 = 1;
                            else
                                descriptors[ndesc].inc2 = atoi(inc_val);
                        }
                        gio_param->dim.slice_Z    = 0;
                        gio_param->dim.slice_comp = 1;
                    }
                }
                else if (!strcasecmp(key, "index")) {
                    if (!strcasecmp(value, "all"))
                        descriptors[ndesc].index_cur = -1;
                }
                else if (!strcasecmp(key, "dimension")) {
                    int   num_dims=0;
                    char *dims_str[MAX_DIMS], *str;
                    descriptors[ndesc].num_dims = 0;
                    if ((str = strtok(value, " \t")) != NULL) {
                        dims_str[num_dims++] = str;
                        while ((str = strtok(NULL, " \t")) != NULL)
                            dims_str[num_dims++] = str;
                        gio_set_dimension(ndesc, num_dims, dims_str, data,
                                            &gio_param->grid, &gio_param->dim);
                        have_dimension = 1;
                    }
                }
                else {
                    /* Do nothing. For the no_data possibility, we have a dummy token */
                    if (descriptors[ndesc].no_data == 0) {
                        /* This token is not a "required" token, so we will stuff
                         * it into a field_attr_name token, providing a value has
                         * also been given
                         */
                        if (descriptors[ndesc].num_attrs + 1 > MAX_ATTR) {
                            printf("GIO Error: MAX_ATTR exceeded\n");
                            ABORT
                        }
                        int num_attrs = descriptors[ndesc].num_attrs;
                        strcpy(descriptors[ndesc].field_attr_name[num_attrs],  key);
                        strcpy(descriptors[ndesc].field_attr_value[num_attrs], value);
                        descriptors[ndesc].num_attrs++;
                    }
                }
            }
            ndesc++;
            // if (!strcasecmp(key, "name")) continue;
        }
    }
    fclose(fp);
    data->num_fields = ndesc;

    return 1;
}

/*----< gio_read_data_config() >---------------------------------------------*/
/* this read in *.desc file */
int gio_read_data_config(gio_parameters *gio_param,
                         const char     *descfile)
{
    int bc_len;
    gio_prop *prop = &gio_param->props;
    gio_data *data = &gio_param->data;
    gio_grid *grid = &gio_param->grid;

    data->grid_cells_dim         = grid->gio_grid_size;
    data->grid_all_uniq_crnr_dim = grid->gio_corners_size;
    data->grid_all_uniq_edge_dim = grid->gio_edges_size;

    data->grid_vectors_dim       = 3;

    data->grid_interfaces_dim    = grid->nlevel+1;
    data->grid_layers_dim        = grid->nlevel;

    data->time_dim               = -1; /* IO code will change this to NC_UNLIMITED
                                          when we create the dimensions in the file */
    data->grid_celledges_dim     =  6; /* IO code will change this to 6 (celledges) */
    data->grid_cellcorners_dim   =  6; /* IO code will change this to 6 (cellorners) */
    data->grid_cellneighbors_dim =  6; /* IO code will change this to 6 (cellneighbors) */
    data->grid_endpoints_dim     =  2; /* IO code will change this to 2 (endpoints) */
    data->num_grid_flds          =  0; /* Total number of static fields describing grid */

    /* Initialize dimension parsing */
    gio_init_dims(&gio_param->dim);

    if (gio_param->gio_me == 0)
        read_io_desc(gio_param, descfile);

    MPI_Bcast(&prop->num_props, 1, MPI_INT, 0, MPI_COMM_WORLD);
    bc_len = prop->num_props * sizeof(prop_descriptor);
    MPI_Bcast(prop->gio_properties, bc_len, MPI_BYTE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&data->num_fields, 1, MPI_INT, 0, MPI_COMM_WORLD);
    bc_len = data->num_fields * sizeof(data_descriptor);
    MPI_Bcast(data->gio_descriptors, bc_len, MPI_BYTE, 0, MPI_COMM_WORLD);

    /* Bcast below may not be necessary, but only root has non-zero value */
    MPI_Bcast(&data->field_size_check, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&data->num_grid_flds, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&data->num_fields, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gio_param->dim.slice_Z, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gio_param->dim.slice_comp, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* This ensures that all processors know that variable sizes will exceed
     * 32 bit indexing. Required for pnetcdf at a minimum */
    data->have_large_data = 0;
    if (data->field_size_check > 2147483648)
        data->have_large_data = 1;

    /* Finish of setting up dimensions for writing fields to files */

    return 1;
}

/*----< gio_create_avg_descriptor() >------------------------------------------*/
static
int gio_create_avg_descriptor(gio_data *data,
                              int       field_idx_orig)
{
    if (data->num_fields >= MAX_DATA_FIELDS) {
        printf("Averaging is creating too many data fields.  Increate MAX_DATA_FIELDS\n");
        ABORT
    }

    /* copy old descriptor properties to new descriptor and then overwrite
     * new features */
    int nf = data->num_fields;
    data->gio_descriptors[nf] = data->gio_descriptors[field_idx_orig];
    data->gio_descriptors[nf].parent = field_idx_orig;
    data->gio_descriptors[nf].average_data = 1;

    /* create new names for averaged data */
    strcpy(data->gio_descriptors[nf].field_name,
           data->gio_descriptors[field_idx_orig].field_name);
    strcat(data->gio_descriptors[nf].field_name, "_avg");

    strcpy(data->gio_descriptors[nf].output_field_name,
           data->gio_descriptors[field_idx_orig].output_field_name);
    strcat(data->gio_descriptors[nf].output_field_name, "_avg");

    data->num_fields++;

    return 1;
}

/*----< gio_inspect_file() >-------------------------------------------------*/
static
int gio_inspect_file(gio_file *file,
                     gio_data *data,
                     int       ifile,
                     int      *found_grid_flds,  /* [MAX_STATIC_GRID_FIELDS] */
                     int      *grid_descriptors) /* [MAX_STATIC_GRID_FIELDS] */
{
    int inv_map[MAX_DATA_FIELDS];
    int i, icnt;
    int which;

    /* Look for grid fields, and add them if we don't have them
     *
     * First, locate all the grid descriptors
     * by looking for them in the gio_descriptors structure
     */ 

    icnt = 0;
    for (i=0; i<MAX_DATA_FIELDS; i++) {
        inv_map[i] = 0;
        if (data->gio_descriptors[i].is_grid) {
            grid_descriptors[icnt] = i;
            inv_map[i] = icnt;
            icnt++;
            if (icnt >= MAX_STATIC_GRID_FIELDS) {
                printf("Number of static grid fields exceeds MAX_STATIC_GRID_FIELDS parameter %d\n",MAX_STATIC_GRID_FIELDS);
                ABORT
            }
        }
    }

    /* What name is field i? */
    for (i=0; i<MAX_STATIC_GRID_FIELDS; i++)
        found_grid_flds[i] = -1;

    if (ifile < 0) {
        /* We won't find file name in the gio_files() structure because
         * we are using the "sep_grid" option
         * So, just set found_grid_flds to -1 and return
         */
        return 1;
    }

    for (i=0; i<file->gio_files[ifile].nflds; i++) {
        which = file->gio_files[ifile].fields[i];
        if (data->gio_descriptors[which].is_grid) {
            if (icnt >= MAX_STATIC_GRID_FIELDS) {
                printf("Number of static grid fields exceeds MAX_STATIC_GRID_FIELDS parameter %d\n",MAX_STATIC_GRID_FIELDS);
                ABORT
            }
            found_grid_flds[inv_map[i]] = 1;
        }
    }

    return 1;
}

/*----< read_io_fcfg() >-----------------------------------------------------*/
/* this is corresponding to gio_read_file_config() in init_term.F90 */
static
int read_io_fcfg(gio_parameters *gio_param,
                 const char     *outfileconfig) /* IN */
{
    int i;
    gio_defaults *defaults = &gio_param->defaults;
    gio_file     *file     = &gio_param->file;
    gio_data     *data     = &gio_param->data;
    /* only changes are file and data->avg_field_ids[] and data->num_averages */

    data->num_averages = 0;

    for (i=0; i<MAX_FILES; i++) {
        if (defaults->gio_default_cdf_path[0] != '\0')
            /* Set directory_name to the gio_default_path default */
            strcpy(file->gio_files[i].directory_name,
                   defaults->gio_default_cdf_path);
        if (defaults->gio_default_frequency > 0)
            /* Set frequency to the supplied default frequency */
            file->gio_files[i].frequency = defaults->gio_default_frequency;
        if (defaults->gio_default_nsamples > 0)
            /* Set the nsamples to the supplied default */
            file->gio_files[i].nsamples = defaults->gio_default_nsamples;
        file->gio_files[i].fileID     = -1;
        file->gio_files[i].grid_level = gio_param->grid.refine_param;
    }

    char line[256];
    int  num_files = -1;
    int ifield_idx=0, nflds_tmp=0;
    FILE *fp;

    OPEN_FILE(outfileconfig)

    /* Read the file.fcfg until we run out of data or hit
       an error or otherwise bail */
    while (fgets(line, 256, fp) != NULL && num_files < MAX_FILES) {
        /* line contains a '\n' at the end */
        char *value;
        char *key = strtok(line, "= \t\n");
        if (key == NULL || key[0] == '#' || key[0] == '!') continue;

        value = strtok(NULL, " !\'\"\t\n");
        if (value == NULL) {
            printf("IO token %s has no value\n", key);
            ABORT
        }

        /* from parallel_load_inputs() called in ZGrd_main */
        if (!strncasecmp(key, "file_prefix", 11)) {
            num_files++;
            strcpy(file->gio_files[num_files].file_prefix, value);
            file->gio_files[num_files].nflds = 0;

            /*  check to see if this is a clustered field file */
            file->gio_files[num_files].is_clustered_file = 0;
            if (strcasestr(value, "clustered") != NULL)
                file->gio_files[num_files].is_clustered_file = 1;

            ifield_idx = 0; /* Set the field index to 0, will reset when we
                               actually have a field token */
            nflds_tmp  = 0; /* Which field for this file */
        }
        else if (!strncasecmp(key, "grid_only", 9)) {
            /* the grid_ONLY option is used to control
             * writing a separate grid file, so that
             * we only write grid fields to it. 
             *
             * We are using default values for the CDF
             * files now, and this is a little cumbersome
             *
             * PONDER and FIX
             */
            file->gio_files[num_files].grid_ONLY = 0;
            if (value[0] == 't' || value[0] == 'T')
                file->gio_files[num_files].grid_ONLY = 1;
            continue;
        }
        else if (!strncasecmp(key, "base_directory", 14))
            strcpy(file->gio_files[num_files].directory_name, value);
        else if (!strncasecmp(key, "nsamples", 8)) {
            file->gio_files[num_files].nsamples = atoi(value);
            if (file->gio_files[num_files].nsamples <= 0) {
                printf("ERROR: nsamples must be a positive integer. setting nsamples to 1\n");
                file->gio_files[num_files].nsamples = 1;
            }
        }
        /* Parse field specifications */
        else if (!strncasecmp(key, "field", 5)) {
            /* Syntax is
             *
             *  field name-of-field 
             *
             * So, we need to get the field name, make sure the name is OK,
             * then get the dimensions and make sure those dimensions are OK
             *
             * should match a name we know about. Let's check.
             * Make sure value (a field name) exists in data structure
             */
            int nd;
            int idesc_idx  = -1;
            ifield_idx++;

            for (nd=0; nd<MAX_DATA_FIELDS; nd++) {
                char *name = data->gio_descriptors[nd].field_name;

                if (!strcasecmp(value, name)) {
                    /* Check to see if field is being averaged. If it is,
                     * then create a new data descriptor for it */
                    value = strtok(NULL, " \t\n");
                    if (value != NULL && !strncasecmp(value, "avg", 3)) {
                        idesc_idx = data->num_fields; /* index of new element */
                        gio_create_avg_descriptor(data, nd);
                        data->avg_field_ids[data->num_averages++] = idesc_idx;
                    }
                    else {
                        /* Index of matching field in the gio_descriptors structure */
                        idesc_idx = nd;
                    }
                    break; /* Get out of loop nd */
                }
            }

            if (idesc_idx < 0) {
                printf("IO field not specified in config file: %s\n", value);
                ABORT
            }

            /* fields[nflds] is an index into gio_descriptors */
            file->gio_files[num_files].fields[nflds_tmp] = idesc_idx;

            /* grid fields are special. See if the field name we will write 
               is a grid field */
            file->gio_files[num_files].is_grid[nflds_tmp] = 0;

            /* Find out if this is the "time" field */
            if (data->gio_descriptors[idesc_idx].is_grid)
                file->gio_files[num_files].is_grid[nflds_tmp] = 1;
            else if (!strncmp(data->gio_descriptors[idesc_idx].field_name, "time", 4))
                file->gio_files[num_files].is_time_field[nflds_tmp] = 1;
            else
                file->gio_files[num_files].is_grid[nflds_tmp] = 0;

            file->gio_files[num_files].nflds = ++nflds_tmp;
        }
        else if (!strncasecmp(key, "frequency", 9))
            file->gio_files[num_files].frequency = atoi(value);
        else if (!strncasecmp(key, "coarsen", 7)) {
            file->gio_files[num_files].grid_level = atoi(value);

            /* Parse additional token and see if we need to average
               instead of sample data */
            value = strtok(NULL, " \t");
            if (!strncasecmp(value, "average", 7))
                file->gio_files[num_files].crsn_average = 1;
        }
        else {
            printf("IO Invalid token: %s\n",value);
            ABORT
        }
    } /* while */

    fclose(fp);

    /*===================================================================
     * CDF default work:
     *
     * After we have read the fcfg file, we need to make some
     * final "default" checks. 
     *
     * if the user wants the grid data place in all files, then
     * we need to look through all the files specified (xxx)
     * and make sure all the grid fields are included
     *
     * This will also mean that we have to update the dimension info 
     */

    /* A total of 13 grid fields currently */
    int ifile;
    int found_grid_flds[MAX_STATIC_GRID_FIELDS];
    int grid_descriptors[MAX_STATIC_GRID_FIELDS];

    if (defaults->gio_grid_option == 1) {
        /* No grid data, so we don't have to look */
        for (ifile=0; ifile<=num_files; ifile++)
            file->gio_files[ifile].contain_grids = 0;
    }
    else if (defaults->gio_grid_option == 2) {
        /* User wants all grid data in all output files */

        /* num_files is the last index */
        for (ifile=0; ifile<=num_files; ifile++) {
            /* no add grid fields if this file contains coarsened data */
            if (file->gio_files[ifile].grid_level !=
                gio_param->grid.refine_param) {
                file->gio_files[ifile].contain_grids = 0;
                continue;
            }
            file->gio_files[ifile].contain_grids = 1;

            /* Add grid fields. */
            gio_inspect_file(file, data, ifile, found_grid_flds,
                             grid_descriptors);

            /* Look through the found_grid_flds array and
             * see if we need to add that grid field to this file */
            nflds_tmp = file->gio_files[ifile].nflds;

            int grid_fld;
            for (grid_fld=0; grid_fld<data->num_grid_flds; grid_fld++) {
                if (found_grid_flds[grid_fld] < 0) {
                    /* Add grid field to the list of fields for this file */
                    file->gio_files[ifile].fields[nflds_tmp]  =
                          grid_descriptors[grid_fld];
                    file->gio_files[ifile].is_grid[nflds_tmp] = 1;
                    file->gio_files[ifile].nflds              = ++nflds_tmp;
                }
            }
        }
    }
    else if (defaults->gio_grid_option == 3) {
        /* User wants a separate file for grid variables.
         * We will add a file to the list, named as specified in
         * gio_sep_grid_file
         */
        if (num_files == MAX_FILES-1) {
            printf("IO max files exceeded.  Increase MAX_FILES\n");
            ABORT
        }

        /* files that are not grid-only */
        for (ifile=0; ifile<=num_files; ifile++)
            file->gio_files[ifile].contain_grids = 0;

        /* grid-only file */
        num_files++;
        strcpy(file->gio_files[num_files].file_prefix,
               defaults->gio_sep_grid_file);
        file->gio_files[num_files].nflds     = 0;
        file->gio_files[num_files].grid_ONLY = 1;
        file->gio_files[num_files].contain_grids = 1;

        /* Call gio_inspect_file() just to get the indices of the grid
         * descriptors pass the "ifile" argument to gio_inspect_file as
         * "-1" so we will not just return the found_grid_flds as all
         * missing */
        gio_inspect_file(file, data, -1, found_grid_flds, grid_descriptors);

        /* Add the grid fields to this new file */
        for (i=0; i<MAX_DATA_FIELDS; i++)
            file->gio_files[num_files].is_grid[i] = 1;

        file->gio_files[num_files].nflds   = data->num_grid_flds;

        for (i=0; i<data->num_grid_flds; i++)
            file->gio_files[num_files].fields[i] = grid_descriptors[i];
    }
    /*  End of default work */

    /*===================================================================
     * Construct restart file
     */
    num_files++;
    nflds_tmp = 0;
    for (i=0; i<data->num_fields; i++) {
        if (data->gio_descriptors[i].is_restart) {
            file->gio_files[num_files].fields[nflds_tmp] = i;
            if (!strncmp(data->gio_descriptors[nflds_tmp].field_name,
                         "time",4))
                file->gio_files[num_files].is_time_field[nflds_tmp] = 1;
            file->gio_files[num_files].nflds = ++nflds_tmp;
        }
    }

    /* If restart fields found finish file parameters, otherwise get rid
       of file */
    if (nflds_tmp > 0) {
        file->restart_idx = num_files;
        /* TODO trim off file extension */
        strcpy(file->gio_files[num_files].file_prefix,
               file->gio_restart_fname);
        file->gio_files[num_files].grid_ONLY = 0;
        file->gio_files[num_files].nsamples  = 1;
        file->gio_files[num_files].frequency = file->gio_restart_interval;
    }
    else {
        num_files--;
        file->restart_idx = 0;
    }

    /* For each file, we will need to know what dimensions are required
     *
     * We can figure this out using info in the gio_descriptors structure
     */
    int file_dims[MAX_DIMS];
    for (i=0; i<=num_files; i++) {
        int ifld, file_num_dims=0;
        for (ifld=0; ifld<file->gio_files[i].nflds; ifld++) {
            int index = file->gio_files[i].fields[ifld];
            int ndims = data->gio_descriptors[index].num_dims;
            if (data->gio_descriptors[index].no_data)
                /* fields with no data will not have dimensions
                   so don't bother looking */
                continue;

            int nd;
            for (nd=0; nd<ndims; nd++) {
                int idim = data->gio_descriptors[index].field_dims[nd];
                int nd2, match=-1;
                for (nd2=0; nd2<file_num_dims; nd2++) {
                    if (idim == file_dims[nd2]) {
                        match = nd2;
                        break;
                    }
                }
                if (match == -1)
                    file_dims[file_num_dims++] = idim;
            }
        }
    }

    file->num_files = num_files + 1;
/* check the last one, it looks like for restart, but has no file_prefix defined */

    return 1;
}

/*----< gio_read_file_config() >---------------------------------------------*/
/* this read in *.fcfg file */
int gio_read_file_config(gio_parameters *gio_param,
                         const char     *outfileconfig)
{
#if 1
    gio_data *data = &gio_param->data;
    gio_file *file = &gio_param->file;

    if (gio_param->gio_me == 0)
        read_io_fcfg(gio_param, outfileconfig);

    MPI_Bcast(&data->num_averages, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (data->num_averages > 0)
        MPI_Bcast(data->avg_field_ids, data->num_averages, MPI_INT, 0,
                  MPI_COMM_WORLD);

    MPI_Bcast(&data->num_fields, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(data->gio_descriptors, sizeof(data_descriptor)*data->num_fields,
              MPI_BYTE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&file->restart_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&file->num_files, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (file->num_files > 0)
        MPI_Bcast(file->gio_files, sizeof(file_descriptor)*file->num_files,
                  MPI_BYTE, 0, MPI_COMM_WORLD);
#else
    read_io_fcfg(gio_param, outfileconfig);
#endif

    return 1;
}

