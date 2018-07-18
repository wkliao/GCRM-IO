/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: io_netcdfx.c 4609 2017-12-07 07:26:38Z wkliao $
 */

#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>   /* access() */
#include <sys/time.h> /* gettimeofday() */
#include <time.h>     /* strftime() */

#include <mpi.h>
#include <pnetcdf.h>

#include "gio.h"
#include "util.h"


static int gio_create(gio_parameters *gio_param, MPI_Comm comm, char *fname,
                      int *fid);
static int gio_inq_dimlen(gio_parameters *gio_param, int fid, int dimid,
                      MPI_Offset *dimlen);
static int gio_inq_dim(gio_parameters *gio_param, int fid, int dimid,
                      char *dimname, MPI_Offset *dimlen);
static int gio_enddef(gio_parameters *gio_param, int fid);
static int gio_put_att_text(gio_parameters *gio_param, int fid, int fldid,
                      char *attname, char *value);
static int gio_put_att_double(gio_parameters *gio_param, int fid, int fldid,
                      char *attname, double value);
static int gio_get_att_double(gio_parameters *gio_param, int fid, int fldid,
                      char *attname, double *value);
static int gio_put_att_int(gio_parameters *gio_param, int fid, int fldid,
                      char *attname, int value);
static int gio_get_att_int(gio_parameters *gio_param, int fid, int fldid,
                      char *attname, int *value);
static int gio_put_att_intarr(gio_parameters *gio_param, int fid, int fldid,
                      char *attname, int *value, int arrlen);
static int gio_get_att_intarr(gio_parameters *gio_param, int fid, int fldid,
                      char *attname, int *value, MPI_Offset arrlen);
static int gio_def_dim(gio_parameters *gio_param, int fid, char *dimname,
                      MPI_Offset dimsize, int *dimID);
static int gio_def_var(gio_parameters *gio_param, int fid, char *fieldname,
                      nc_type vartype, int numdims, int *dims, int *fldid,
                      int *csizes);
static int gio_put_double(gio_parameters *gio_param, int fid, int fldid,
                      MPI_Offset *start, MPI_Offset *count, double *buffer);
static int gio_get_double(gio_parameters *gio_param, int fid, int fldid,
                      MPI_Offset *start, MPI_Offset *count, double *buffer);
/*
static int gio_put_real(gio_parameters *gio_param, int fid, int fldid,
                      MPI_Offset *start, MPI_Offset *count, float *buffer);
*/
static int gio_get_real(gio_parameters *gio_param, int fid, int fldid,
                      MPI_Offset *start, MPI_Offset *count, float *buffer);
static int gio_put_int(gio_parameters *gio_param, int fid, int fldid,
                      MPI_Offset *start, MPI_Offset *count, int *buffer);
static int gio_get_int(gio_parameters *gio_param, int fid, int fldid,
                      MPI_Offset *start, MPI_Offset *count, int *buffer);
static void gio_create_dims(gio_parameters *gio_param, int nftype);
static void gio_add_field_descriptor(gio_parameters *gio_param, int ifile,
                      int nfield);
static void gio_set_property(gio_parameters *gio_param, gio_prop *props,
                      int fileID, int nprops, MPI_Comm comm);
static void gio_add_verticalgrid(gio_parameters *gio_param, int ifile);
static void gio_write_grid(gio_parameters *gio_param, int ifile);
static void gio_write_restart_attrs(gio_parameters *gio_param, int ifile);
static void gio_coarsen_field(gio_parameters *gio_param, int field_idx,
                      int ifile, MPI_Offset offset, MPI_Offset *start,
                      MPI_Offset *count, void *buffer, MPI_Offset cellsize);;

// #define r0create

/*---< gio_close() >----------------------------------------------------------*/
static
int gio_close(int fid)
{
    return ncmpi_close(fid);
}


/*---< gio_check_file() >-----------------------------------------------------*/
/* Check the existence of a file. Let rank 0 check and broadcast the result.
   This strategy is in case the number of MPI processes is high, having all
   processes checking the file may cost high due to accessing the file
   metadata server concurrently.
*/
static
int gio_check_file(MPI_Comm  comm,
                   char     *fname)
{
    int my_rank, exists;

    MPI_Comm_rank(comm, &my_rank);

    /* if access() is available, use it to check if file already exists */
    exists = 0;
    if (my_rank == 0) { /* root checks if file exists */
        /* remove the file system type prefix name if there is any.
         * For example, path=="lustre:/home/foo/testfile.nc",
         * use "/home/foo/testfile.nc" when calling access()
         */
        char *filename = strchr(fname, ':');
        if (filename == NULL) /* no prefix */
            filename = (char*)fname;
        else
            filename++;

        if (access(filename, F_OK) == 0) exists = 1;
    }
#ifndef r0create
    /* inform all processes whether the file exists */
    MPI_Bcast(&exists, 1, MPI_INT, 0, comm);
#endif
    return exists;
}

/*---< gio_open() >-----------------------------------------------------------*/
int gio_open(gio_parameters *gio_param,
             MPI_Comm        comm,
             char           *fname,
             int             flags,
             MPI_Info       *used_info,
             int            *fid)
{
    int    err;
    double starttime;
    MPI_Info info;

    starttime = MPI_Wtime();

// #define USE_PREV_INFO
/* if USE_PREV_INFO is defined, the MPI info used at the first time when
   opening a file is re-used for the successive file opens */

#ifdef USE_PREV_INFO
    if (*used_info != MPI_INFO_NULL) {
#endif
        MPI_Info_create(&info);
        MPI_Info_set(info, "romio_no_indep_rw", "true");
        MPI_Info_set(info, "romio_ds_write",    "disable");
        MPI_Info_set(info, "romio_cb_write",    "enable");

        /* cb_buffer_size will be ignored by cray's MPI-IO on Lustre, unless
           you also mess with cb_align.
           This info will take effect on other file systems.
        call MPI_Info_set (info, " cb_buffer_size ",   "67108864", err)
         */

        err = ncmpi_open(comm, fname, flags, info, fid);
        MPI_Info_free(&info);
#ifdef USE_PREV_INFO
        ncmpi_get_file_info(*fid, used_info);
    }
    else
        err = ncmpi_open(comm, fname, flags, *used_info, fid);
#else
    if (*used_info == MPI_INFO_NULL)
        ncmpi_get_file_info(*fid, used_info);
#endif
    gio_param->time_in_nf_open += MPI_Wtime() - starttime;
    if (err != NC_NOERR) return err;

    if (gio_param->gio_iotype == MODE_NONBLOCKING_INDEP ||
        gio_param->gio_iotype == MODE_BLOCKING_INDEP)
        err = ncmpi_begin_indep_data(*fid);

    return err;
}

/*---< gio_create() >---------------------------------------------------------*/
static
int gio_create(gio_parameters *gio_param,
               MPI_Comm        comm,
               char           *fname,
               int            *fid)
{
    double starttime;
    int  mode, err;
    char *env_str, *hdr_space, *align_size;
    MPI_Info info;

    starttime = MPI_Wtime();

    /* set the file access mode */
    mode = NC_CLOBBER | NC_64BIT_OFFSET;
    if (gio_param->data.have_large_data)
        mode = NC_CLOBBER | NC_64BIT_DATA;

#if defined(USE_PREV_INFO) && !defined(r0create)
    if (gio_param->used_info == MPI_INFO_NULL) {
#endif
        /* Set some defaults for header buffer and alignment but accept
           env var overrides
           http://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/VariableAlignment
         */
        align_size = "1";       /* no alignment */
        hdr_space = "1048576";  /* give header at least 1 MB */
        if (gio_param->grid.refine_param < 8)
            hdr_space = "2048";   /* not much for low res */
        env_str = getenv("nc_header_align_size");
        if (env_str != NULL) hdr_space = env_str;
        env_str = getenv("nc_var_align_size");
        if (env_str != NULL) align_size = env_str;

        /* These are hints that may allow extra room in the netcdf header */
        MPI_Info_create(&info);
        MPI_Info_set (info, "nc_header_align_size", hdr_space);
        MPI_Info_set (info, "nc_var_align_size",    align_size);
        MPI_Info_set (info, "romio_no_indep_rw",    "true");
        MPI_Info_set (info, "romio_ds_write",       "disable");
        MPI_Info_set (info, "romio_cb_write",       "enable");

        /* cb_buffer_size will be ignored by cray's library unless you also
           mess with cb_align.  Use for WLK's library.
           call MPI_Info_set (info, " cb_buffer_size ",   "67108864", err)
         */
        err = ncmpi_create(comm, fname, mode, info, fid);
        MPI_Info_free(&info);
#if defined(USE_PREV_INFO) && !defined(r0create)
        err = ncmpi_get_file_info(*fid, &gio_param->used_info);
    }
    else
        err = ncmpi_create(comm, fname, mode, gio_param->used_info, fid);
#endif
#ifndef USE_PREV_INFO
    if (gio_param->used_info == MPI_INFO_NULL)
        err = ncmpi_get_file_info(*fid, &gio_param->used_info);
#endif
    gio_param->time_in_nf_create += MPI_Wtime() - starttime;

    return err;
}

/*---< gio_inq_varid() >------------------------------------------------------*/
/* get the varid for a given variable name.  */
int gio_inq_varid(gio_parameters *gio_param,
                  int             fid,
                  char           *fieldname,
                  int            *fldid)
{
    double starttime = MPI_Wtime();
    int err = ncmpi_inq_varid(fid, fieldname, fldid);
    gio_param->time_in_nf_inq_varid += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_inq_dimlen() >-----------------------------------------------------*/
/* get the dimension length for a given variable id. */
static
int gio_inq_dimlen(gio_parameters *gio_param,
                   int             fid,
                   int             dimid,
                   MPI_Offset     *dimlen)
{
    double starttime = MPI_Wtime();
    int err = ncmpi_inq_dimlen(fid, dimid, dimlen);
    gio_param->time_in_nf_inq_dimlen += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_inq_dim() >--------------------------------------------------------*/
/* get the dimension name and length for a variable */
static
int gio_inq_dim(gio_parameters *gio_param,
                int             fid,
                int             dimid,
                char           *dimname,
                MPI_Offset     *dimlen)
{
    double starttime = MPI_Wtime();
    int err = ncmpi_inq_dim (fid, dimid, dimname, dimlen);
    gio_param->time_in_nf_inq_dimlen += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_enddef() >---------------------------------------------------------*/
/* tell io library we are done defining data */
static
int gio_enddef(gio_parameters *gio_param,
               int             fid)
{
    int err;
    double starttime = MPI_Wtime();
    err = ncmpi_enddef(fid);
    gio_param->time_in_nf_enddef += MPI_Wtime() - starttime;
    if (err != NC_NOERR) return err;

    if (gio_param->gio_iotype == MODE_NONBLOCKING_INDEP ||
        gio_param->gio_iotype == MODE_BLOCKING_INDEP)
        err = ncmpi_begin_indep_data(fid);

    return err;
}

/*---< gio_put_att_text() >---------------------------------------------------*/
/* Write a text attribute and time the operations */
static
int gio_put_att_text(gio_parameters *gio_param,
                     int             fid,
                     int             fldid,
                     char           *attname,
                     char           *value)
{
    int err = NC_NOERR;  /* initialize err, in case attlen .eq. 0 */
    double starttime = MPI_Wtime();
    MPI_Offset attlen;

    if ((attlen = strlen(value)) > 0)
        err = ncmpi_put_att_text(fid, fldid, attname, attlen, value);

    gio_param->time_in_nf_put_att += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_put_att_double() >-------------------------------------------------*/
/* Write a double-type attribute and time the operation */
static
int gio_put_att_double(gio_parameters *gio_param,
                       int             fid,
                       int             fldid,
                       char           *attname,
                       double          value)
{
    int    err;
    double starttime = MPI_Wtime();
    err = ncmpi_put_att_double(fid, fldid, attname, NC_DOUBLE, 1, &value);
    gio_param->time_in_nf_put_att += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_get_att_double() >-------------------------------------------------*/
/* get a double-type attribute and time the operation */
static
int gio_get_att_double(gio_parameters *gio_param,
                       int             fid,
                       int             fldid,
                       char           *attname,
                       double         *value)
{
    double starttime = MPI_Wtime();
    int err = ncmpi_get_att_double(fid, fldid, attname, value);
    gio_param->time_in_nf_get_att += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_put_att_int() >----------------------------------------------------*/
/* Write an integer attribute and time the operation */
static
int gio_put_att_int(gio_parameters *gio_param,
                    int             fid,
                    int             fldid,
                    char           *attname,
                    int             value)
{
    double starttime = MPI_Wtime();
    int err = ncmpi_put_att_int(fid, fldid, attname, NC_INT, 1, &value);
    gio_param->time_in_nf_put_att += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_get_att_int() >----------------------------------------------------*/
/* Read an integer attribute and time the operation */
static
int gio_get_att_int(gio_parameters *gio_param,
                    int             fid,
                    int             fldid,
                    char           *attname,
                    int            *value)
{
    double starttime = MPI_Wtime();
    int err = ncmpi_get_att_int(fid, fldid, attname, value);
    gio_param->time_in_nf_get_att += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_put_att_intarr() >-------------------------------------------------*/
/* Write an integer array attribute and time the operations */
static
int gio_put_att_intarr(gio_parameters *gio_param,
                       int             fid,
                       int             fldid,
                       char           *attname,
                       int            *value,
                       int             arrlen)
{
    int err = NC_NOERR;

    if (arrlen > 0) {
        double starttime = MPI_Wtime();
        err = ncmpi_put_att_int(fid, fldid, attname, NC_INT, arrlen, value);
        gio_param->time_in_nf_put_att += MPI_Wtime() - starttime;
    }
    return err;
}

/*---< gio_get_att_intarr() >-------------------------------------------------*/
/* Read an integer attribute and time the operations */
static
int gio_get_att_intarr(gio_parameters *gio_param,
                       int             fid,
                       int             fldid,
                       char           *attname,
                       int            *value,
                       MPI_Offset      arrlen)
{
    double starttime = MPI_Wtime();
    int err = ncmpi_get_att_int(fid, fldid, attname, value);
    gio_param->time_in_nf_get_att += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_def_dim() >--------------------------------------------------------*/
static
int gio_def_dim(gio_parameters *gio_param,
                int             fid,
                char           *dimname,
                MPI_Offset      dimsize,
                int            *dimID)
{
    double starttime = MPI_Wtime();
    int err = ncmpi_def_dim(fid, dimname, dimsize, dimID);
    gio_param->time_in_nf_def_dim += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_def_var() >--------------------------------------------------------*/
static
int gio_def_var(gio_parameters *gio_param,
                int             fid,
                char           *fieldname,
                nc_type         vartype,
                int             numdims,
                int            *dims,
                int            *fldid,
                int            *csizes)
{
    double starttime = MPI_Wtime();
    int err = ncmpi_def_var(fid, fieldname, vartype, numdims, dims, fldid);
    gio_param->time_in_nf_def_var += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_put_double() >-----------------------------------------------------*/
static
int gio_put_double(gio_parameters *gio_param,
                   int             fid,
                   int             fldid,
                   MPI_Offset     *start,
                   MPI_Offset     *count,
                   double         *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    if (gio_param->gio_iotype == MODE_BLOCKING_COLL)
        err = ncmpi_put_vara_double_all(fid, fldid, start, count, buffer);
    else
        err = ncmpi_put_vara_double(fid, fldid, start, count, buffer);
    gio_param->time_in_nf_put_var += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_get_double() >-----------------------------------------------------*/
static
int gio_get_double(gio_parameters *gio_param,
                   int             fid,
                   int             fldid,
                   MPI_Offset     *start,
                   MPI_Offset     *count,
                   double         *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    if (gio_param->gio_iotype == MODE_BLOCKING_COLL)
        err = ncmpi_get_vara_double_all(fid, fldid, start, count, buffer);
    else
        err = ncmpi_get_vara_double(fid, fldid, start, count, buffer);
    gio_param->time_in_nf_get_var += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_put_real() >-------------------------------------------------------*/
int gio_put_real(gio_parameters *gio_param,
                 int             fid,
                 int             fldid,
                 MPI_Offset     *start,
                 MPI_Offset     *count,
                 float          *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    if (gio_param->gio_iotype == MODE_BLOCKING_COLL)
        err = ncmpi_put_vara_float_all(fid, fldid, start, count, buffer);
    else
        err = ncmpi_put_vara_float(fid, fldid, start, count, buffer);
    gio_param->time_in_nf_put_var += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_get_real() >-------------------------------------------------------*/
static
int gio_get_real(gio_parameters *gio_param,
                 int             fid,
                 int             fldid,
                 MPI_Offset     *start,
                 MPI_Offset     *count,
                 float          *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    if (gio_param->gio_iotype == MODE_BLOCKING_COLL)
        err = ncmpi_get_vara_float_all(fid, fldid, start, count, buffer);
    else
        err = ncmpi_get_vara_float(fid, fldid, start, count, buffer);
    gio_param->time_in_nf_get_var += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_put_int() >--------------------------------------------------------*/
static
int gio_put_int(gio_parameters *gio_param,
                int             fid,
                int             fldid,
                MPI_Offset     *start,
                MPI_Offset     *count,
                int            *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    if (gio_param->gio_iotype == MODE_BLOCKING_COLL)
        err = ncmpi_put_vara_int_all(fid, fldid, start, count, buffer);
    else
        err = ncmpi_put_vara_int(fid, fldid, start, count, buffer);
    gio_param->time_in_nf_put_var += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_get_int() >--------------------------------------------------------*/
static
int gio_get_int(gio_parameters *gio_param,
                int             fid,
                int             fldid,
                MPI_Offset     *start,
                MPI_Offset     *count,
                int            *buffer)
{
    int err;
    double starttime = MPI_Wtime();
    if (gio_param->gio_iotype == MODE_BLOCKING_COLL)
        err = ncmpi_get_vara_int_all(fid, fldid, start, count, buffer);
    else
        err = ncmpi_get_vara_int(fid, fldid, start, count, buffer);
    gio_param->time_in_nf_get_var += MPI_Wtime() - starttime;
    return err;
}

/*---< gio_open_file() >------------------------------------------------------*/
/* Open/Create a file as needed
   If the file must be created, then the header is defined
*/
int gio_open_file(gio_parameters *gio_param,
                  int             ifile,
                  char           *fname)
{
    int err, nfield, nprops, fmode, created=0;
    int restart_idx = gio_param->file.restart_idx;
    MPI_Comm comm=gio_param->gio_io_comm;
    file_descriptor *gio_files = &gio_param->file.gio_files[ifile];

    if (gio_param->gio_do_io) {
#ifdef r0create
      comm = MPI_COMM_SELF;
      if (gio_param->gio_me == 0)
#else
      if (1)
#endif
      {
        if ((ifile == restart_idx) || (! gio_check_file(comm, fname))) {
            err = gio_create(gio_param, comm, fname, &gio_files->fileID);
            if (err != NC_NOERR) {
                char err_msg[128];
                sprintf(err_msg, "Unable to create file: %s",fname);
                NC_CHECK(err, err_msg, 1)
            }

            created = 1;
          
            /* define dimensions in file header */
            gio_create_dims(gio_param, ifile);

            /* define variables and their attributes */
            for (nfield=0; nfield<gio_files->nflds; nfield++) {
                if (ifile == restart_idx) {
                    if (gio_files->is_grid[nfield] == 0)
                        gio_add_field_descriptor(gio_param, ifile, nfield);
                }
                else
                    gio_add_field_descriptor(gio_param, ifile, nfield);
            }

            /* define global attribute */
            for (nprops=0; nprops<gio_param->props.num_props; nprops++)
                gio_set_property(gio_param, &gio_param->props,
                                 gio_files->fileID, nprops, comm);

            /* write global time stamp attribute */
            gio_set_property(gio_param, &gio_param->props, gio_files->fileID,
                             -1, comm);

            if (ifile == restart_idx)
                gio_write_restart_attrs(gio_param, ifile);

            err = gio_enddef(gio_param, gio_files->fileID);
            NC_CHECK(err, "Failure in gio_enddef", 0)

            /* we only want to write vertical grid variable from the head
               node and not for restart files */
            if (ifile != restart_idx && gio_param->gio_iotype == MODE_BLOCKING_INDEP)
                gio_add_verticalgrid(gio_param, ifile);

#ifdef r0create
            gio_close_file(gio_param, ifile);
#endif
        }  /* checkfile */
      }    /* Endif for if (gio_param->gio_me == 0) */
    } /* End if i'm an io proc */

#ifdef r0create
    /* Gotta tell everyone since we check this later for grid writing
       This is required only when file create is done by root process */
    MPI_Bcast(&created, 1, MPI_INT, 0, gio_param->gio_world_comm);

    /* IO procs must open the file now.  
       After this point, all the IO processes have an open file
       The non IO processes have "ifile" defined. */
    if (gio_param->gio_do_io)
#else
    if (gio_param->gio_do_io && ! created)
        /* when all procs create but file exists, IO procs must open file now */
#endif
    {
        comm = gio_param->gio_io_comm;
        fmode = NC_WRITE;
        err = gio_open(gio_param, comm, fname, fmode, &gio_param->used_info,
                       &gio_files->fileID);
        if (err != NC_NOERR) {
            char err_msg[128];
            sprintf(err_msg, "Unable to open file: %s", fname);
            NC_CHECK(err, err_msg, 0)
        }
    }

    /* This must be called on all processors; created must have been bcast */
    if (created && ifile != restart_idx &&
        gio_param->file.gio_files[ifile].contain_grids)
        gio_write_grid(gio_param, ifile);

    return created;
}

/*---< gio_close_file() >-----------------------------------------------------*/
void gio_close_file(gio_parameters *gio_param,
                    int             ifile)
{
    if (gio_param->gio_do_io) {  /* IO procs only */
        char err_msg[128], *fname=gio_param->file.gio_files[ifile].file_prefix;
        int err, ncid=gio_param->file.gio_files[ifile].fileID;
        double starttime = MPI_Wtime();
        MPI_Offset size;

        err = ncmpi_inq_put_size(ncid, &size);
        if (err != NC_NOERR) {
            sprintf(err_msg, "ncmpi_inq_put_size(): %s", fname);
            NC_CHECK(err, err_msg, 0)
        }
        gio_param->pnc_write_amount += size;

        err = ncmpi_inq_get_size(ncid, &size);
        if (err != NC_NOERR) {
            sprintf(err_msg, "ncmpi_inq_get_size(): %s", fname);
            NC_CHECK(err, err_msg, 0)
        }
        gio_param->pnc_read_amount += size;

        err = gio_close(ncid);
        if (err != NC_NOERR) {
            sprintf(err_msg, "closing file: %s", fname);
            NC_CHECK(err, err_msg, 0)
        }
        gio_param->time_in_nf_close += MPI_Wtime() - starttime;
    }
    gio_param->file.gio_files[ifile].fileID = -1;
}

/*---< gio_sync_file() >------------------------------------------------------*/
/* Use sync to flush yet avoid open/close cycle.  Both are costly.
   open has proved to be very costly.
*/
void gio_sync_file(gio_parameters *gio_param,
                   int             fileID)
{
    int err;
    double starttime;

    if (! gio_param->gio_do_io) return;  /* IO procs only */

    starttime = MPI_Wtime();
    err = ncmpi_sync(fileID);
    NC_CHECK(err, "synching file.", 0)
    gio_param->time_in_nf_close += MPI_Wtime() - starttime;
}

/*---< total_number_of_cells() >----------------------------------------------*/
int total_number_of_cells(gio_parameters *gio_param)
{
    int i, j, mask_idx, num_cells, nsd, nvalue;
    int imin, imax, jmin, jmax, ilen, jlen;
    double *mask_ptr;

    data_descriptor *dscptrs = gio_param->data.gio_descriptors;

    nvalue = 0;

    /* Find index of mask array */
    mask_idx = 0;
    for (i=0; i<gio_param->data.num_fields; i++)
        if (!strcmp(dscptrs[i].field_name, "clus_mask"))
            break;

    if (i == gio_param->data.num_fields) {
        printf("Failure in gio_set_mask_size to find mask array\n");
        ABORT
    }
    mask_idx = i;

    for (nsd=0; nsd<dscptrs[mask_idx].num_data_blocks; nsd++) {
        num_cells = 0;
        mask_ptr = dscptrs[mask_idx].gio_data_blocks[nsd].double_data;
        imin = dscptrs[mask_idx].gio_data_blocks[nsd].imin;
        imax = dscptrs[mask_idx].gio_data_blocks[nsd].imax;
        jmin = dscptrs[mask_idx].gio_data_blocks[nsd].jmin;
        jmax = dscptrs[mask_idx].gio_data_blocks[nsd].jmax;
        ilen = imax - imin + 1;
        jlen = jmax - jmin + 1;
        for (j=0; j<jlen; j++) {
            for (i=0; i<ilen; i++) {
                if (mask_ptr[j*gio_param->grid.gio_istride + i] != 0.0)
                    num_cells++;
            }
        }
        nvalue += num_cells;
    }
    return nvalue;
}


/*---< gio_write() >----------------------------------------------------------*/
/* Write a buffer which will be one of integer, float or double as 
   determined by the field_idx's data type.
*/
void gio_write(gio_parameters *gio_param,
               int             field_idx,
               int             ifile,
               void           *buffer,
               MPI_Offset      offset,
               int             time_idx,
               MPI_Offset      blen)
{
    char *fieldname;
    int i, fileID, fieldID, bytes, err;
    int numdims, celldim, timedim, indexdim;
    MPI_Offset cellsize, start[4], count[4];
    double starttime = MPI_Wtime();

    /* cluster stuff */
    int block_index, my_rank, offset_sum=offset;
    int cluster_cells, cmask_length;

    if (! gio_param->gio_do_io) return;  /* for io procs only */

    data_descriptor *dscptr    = &gio_param->data.gio_descriptors[field_idx];
    file_descriptor *gio_files = &gio_param->file.gio_files[ifile];

    if (blen == 0) {
#if 0
        if (gio_param->gio_iotype == MODE_NONBLOCKING_COLL ||
            gio_param->gio_iotype == MODE_NONBLOCKING_INDEP ||
            gio_param->gio_iotype == MODE_BLOCKING_INDEP)
            return;

        /* for blocking API, participate collective I/O with zero-sized request */
        for (i=0; i<dscptr->num_file_dims; i++) start[i] = count[i] = 0;
        fieldname = dscptr->output_field_name;
        gio_inq_varid(gio_param, gio_files->fileID, fieldname, &fieldID);
        gio_put_real(gio_param, gio_files->fileID, fieldID, start, count, buffer);
#endif
        return;
    }

    /* all cluster related */
    MPI_Comm_rank(gio_param->gio_io_comm, &my_rank);

    /* for clustered data */
    cmask_length  = 1;
    cluster_cells = 0;
    if (gio_files->is_clustered_file && blen > 0) {
        cluster_cells = total_number_of_cells(gio_param);
        offset_sum    = 0;
        cellsize      = 1;
        for (i=0; i<dscptr->cell_idx; i++)
            cellsize *= dscptr->file_dims[i];

        block_index = offset / blen;
        if (gio_param->clus_cells != 0 ||
            gio_param->clus_edges != 0 ||
            gio_param->clus_crns  != 0 ) {

            // block_index++;
            cmask_length = gio_param->cell_mask_sizes[block_index];
            for (i=0; i<=block_index; i++)
                offset_sum += gio_param->cell_mask_sizes[i];
            offset_sum -= cmask_length;
            offset = offset_sum * cellsize;
        }
        else {
            cmask_length = 0 ;
            offset_sum = offset/cellsize;
        }
    }
    /* end cluster related */

    fieldname = dscptr->output_field_name;
    err = gio_inq_varid(gio_param, gio_files->fileID, fieldname, &fieldID);
    if (err != NC_NOERR) {
        char err_msg[128];
        sprintf(err_msg, "getting varid: %s",fieldname);
        NC_CHECK(err, err_msg, 0)
    }

    /*  scan over dimensions to set up count and start arrays array */
    celldim  = dscptr->cell_idx;
    timedim  = dscptr->time_idx;
    indexdim = dscptr->intg_idx;
    numdims  = dscptr->num_file_dims;

    cellsize = 1;
    for (i=0; i<celldim; i++)
        cellsize *= dscptr->file_dims[i];

    for (i=0; i<numdims; i++) {
        count[i] = dscptr->file_dims[i];
        start[i] = 0;
    }

    if (! gio_files->is_clustered_file)
        offset /= cellsize;
    else
        offset = offset_sum;

    if (celldim >= 0) {  /* F2C: shouldn't it be >= 0? */
        if (gio_files->is_clustered_file) {
            if (cmask_length != 0)
                count[celldim] = blen/cellsize;
            else 
                count[celldim] = 0;
        }
        else
            count[celldim] = blen/cellsize;

//     if (indexdim >= 0)
//         start[celldim] = offset % dscptr->file_dims[celldim];
//     else
           start[celldim] = offset;
    }
    if (indexdim >= 0) {
//      start[indexdim] = offset / dscptr->file_dims[celldim];
        count[indexdim] = 1;
        start[indexdim] = time_idx;
    }
    
    if (timedim >= 0) {
        count[timedim] = 1;
        start[timedim] = gio_files->samples_written;
    }

    if (count[celldim] == 0) {
        count[indexdim] = 0;
        count[timedim]  = 0;
        start[celldim]  = 0;
        start[indexdim] = 0;
        start[timedim]  = 0;
    }

    /* Modify data if cell is coarsened. */
    gio_coarsen_field(gio_param, field_idx, ifile, offset, start, count,
                      buffer, cellsize);

    if (!strcmp(dscptr->field_name, "pressure")) {
        if (gio_files->is_clustered_file)
            printf("r/s/c : %d %lld %lld\n", my_rank, start[celldim],
                                             count[celldim]); 
    }

    bytes = 4;
    if (dscptr->data_type == dbl_type) bytes = 8;

    fileID = gio_files->fileID;
/*
    if (gio_files->is_clustered_file) then
        if (dscptr->data_type == flt_type) then
            err = gio_put_real_async(fileID, fieldID, start, count, buffer)
        else if( isinteger ) then
        else if( isdouble ) then
        endif
*/
    /* F2C dimension re-order */
    MPI_Offset sstart[4], scount[4];
    for (i=0; i<numdims; i++) {
        sstart[numdims-1-i] = start[i];
        scount[numdims-1-i] = count[i];
    }
    for (i=0; i<numdims; i++) {
        start[i] = sstart[i];
        count[i] = scount[i];
    }

    if (dscptr->data_type == flt_type) {
        if (gio_param->gio_iotype == MODE_BLOCKING_INDEP)
            err = gio_put_real(gio_param, fileID, fieldID, start, count, buffer);
        else
            err = gio_iput_real(gio_param, fileID, fieldID, start, count, buffer);
    }
    else if (dscptr->data_type == int_type) {
        if (gio_param->gio_iotype == MODE_BLOCKING_INDEP)
            err = gio_put_int(gio_param, fileID, fieldID, start, count, buffer);
        else
            err = gio_iput_int(gio_param, fileID, fieldID, start, count, buffer);
    }
    else if (dscptr->data_type == dbl_type) {
        if (gio_param->gio_iotype == MODE_BLOCKING_INDEP)
            err = gio_put_double(gio_param, fileID, fieldID, start, count, buffer);
        else
            err = gio_iput_double(gio_param, fileID, fieldID, start, count, buffer);
    }
    else
        NC_CHECK(-1, "Unsupported data type", 0)

    NC_CHECK(err, "writing variable.", 0)

    gio_param->num_records_written++;
    gio_param->bytes_API_write += bytes*blen;

    /* Update stats specific to this file */
    dscptr->record_count++;
    dscptr->bytes_written += bytes*blen;
    dscptr->time_for_writes += MPI_Wtime() - starttime;
}

/*---< gio_read() >-----------------------------------------------------------*/
/* Read a variable
   Reads are done from all processors - no IO aggregation at the application
   level
*/
void gio_read(gio_parameters *gio_param,
              int             field_idx,
              int             ifile,
              void           *buffer,
              MPI_Offset      offset,
              int             time_idx,
              MPI_Offset      blen)
{
    char *fieldname;
    int i, fileID, fieldID, bytes, err;
    int numdims, celldim, timedim, indexdim;
    int isgrid;
    double starttime;
    MPI_Offset myoffset, cellsize, start[4], count[4];

    starttime = MPI_Wtime();

    data_descriptor *dscptr = &gio_param->data.gio_descriptors[field_idx];
    file_descriptor *gio_files = &gio_param->file.gio_files[ifile];

    fieldname = dscptr->output_field_name;
    fileID = gio_files->fileID;

    err = gio_inq_varid(gio_param, fileID, fieldname, &fieldID);
    if (err != NC_NOERR) {
        char err_msg[128];
        sprintf(err_msg, "getting varid: %s",fieldname);
        NC_CHECK(err, err_msg, 0)
    }

    if (blen == 0) {
#if 0
        if (gio_param->gio_iotype == MODE_NONBLOCKING_COLL ||
            gio_param->gio_iotype == MODE_NONBLOCKING_INDEP ||
            gio_param->gio_iotype == MODE_BLOCKING_INDEP)
            return;

        /* for blocking API, participate collective I/O with zero-sized request */
        for (i=0; i<dscptr->num_dims; i++) start[i] = count[i] = 0;
        gio_get_real(gio_param, fileID, fieldID, start, count, buffer);
#endif
        return;
    }

    /* scan over dimensions to set up count and start arrays array */
    celldim  = dscptr->cell_idx;
    timedim  = dscptr->time_idx;
    numdims  = dscptr->num_dims;
    indexdim = dscptr->intg_idx;

    cellsize = 1;
    for (i=0; i<celldim; i++)
        cellsize *= dscptr->file_dims[i];

    for (i=0; i<numdims; i++) {
        count[i] = dscptr->file_dims[i];
        start[i] = 0;
    }

    myoffset = offset/cellsize;
    if (celldim >= 0) {  /* F2C: shouldn't it be >= 0? */
        count[celldim] = blen/cellsize;

       if (indexdim >= 0)
           start[celldim] = myoffset % dscptr->file_dims[celldim];
       else
           start[celldim] = myoffset;
    }

    if (indexdim >= 0) {
        count[indexdim] = 1;
        // start[indexdim] = myoffset/dscptr->file_dims[celldim];
        start[indexdim] = time_idx;
    }

    if (timedim >= 0) {
        count[timedim] = 1;
        start[timedim] = gio_files->samples_read;
    }

    bytes = 4;
    if (dscptr->data_type == dbl_type) bytes = 8;
    gio_param->bytes_API_read += bytes*blen;

    isgrid = gio_files[ifile].is_grid[field_idx];
    MPI_Barrier(MPI_COMM_WORLD);

    if (dscptr->data_type == flt_type) {
        if (gio_param->gio_iotype == MODE_BLOCKING_INDEP)
            err = gio_get_real(gio_param, fileID, fieldID, start, count, buffer);
        else
            err = gio_iget_real(gio_param, fileID, fieldID, start, count, buffer);
    }
    else if (dscptr->data_type == int_type) {
        if (gio_param->gio_iotype == MODE_BLOCKING_INDEP)
            err = gio_get_int(gio_param, fileID, fieldID, start, count, buffer);
        else
            err = gio_iget_int(gio_param, fileID, fieldID, start, count, buffer);
    }
    else if (dscptr->data_type == dbl_type) {
        if (gio_param->gio_iotype == MODE_BLOCKING_INDEP)
            err = gio_get_double(gio_param, fileID, fieldID, start, count, buffer);
        else
            err = gio_iget_double(gio_param, fileID, fieldID, start, count, buffer);
    }
    else
        NC_CHECK(-1, "Unsupported data type", 0)

    NC_CHECK(err, "reading variable.", 0)

    // call gio_stat_reads(bytes*blen)
    gio_param->num_records_read++;
    /* Update stats specific to this file */
    dscptr[field_idx].record_count++;
    dscptr[field_idx].bytes_read += bytes*blen;
    dscptr[field_idx].time_for_reads += MPI_Wtime() - starttime;

    /* All IO procs wait here until all IO is done */
    MPI_Barrier(MPI_COMM_WORLD);
}

/*----< gio_create_dims() >---------------------------------------------------*/
static
void gio_create_dims(gio_parameters *gio_param,
                     int             nftype)
{
    int i, err, dimID, coarsen_factor;
    MPI_Offset dimsize;
    char *tmp_dim_name;

    file_descriptor *gio_files = gio_param->file.gio_files;

    for (i=0; i<gio_files[nftype].dimension_num; i++) {
        tmp_dim_name = gio_files[nftype].dimension_names[i];
        /*  Handle time differently */
        if (!strncmp(tmp_dim_name, "time", 4))
            dimsize = NC_UNLIMITED;
        else {
            dimsize = gio_files[nftype].dimension_sizes[i];

            /* Fix up dimension size for coarsened files, if necessary */
            if (gio_files[nftype].grid_level != gio_param->grid.refine_param &&
                (!strncmp(tmp_dim_name, "cells",   5) ||
                 !strncmp(tmp_dim_name, "corners", 7) ||
                 !strncmp(tmp_dim_name, "edges",   6))) {
                coarsen_factor = gio_param->grid.refine_param
                               - gio_files[nftype].grid_level;
                coarsen_factor = POWER2(2*coarsen_factor);
                if (!strncmp(tmp_dim_name, "cells", 5))
                    dimsize =  (dimsize-2)/coarsen_factor + 2;
                else
                    dimsize =  dimsize/coarsen_factor;
            }
            if (gio_files[nftype].is_clustered_file ) {
                if (gio_param->clus_cells != 0 ||
                    gio_param->clus_crns  != 0 ||
                    gio_param->clus_edges != 0) {

                    if (!strncmp(tmp_dim_name, "cells", 5))
                        dimsize = gio_param->clus_cells;
                    else if (!strncmp(tmp_dim_name, "corners", 7))
                        dimsize = gio_param->clus_crns;
                    else if (!strncmp(tmp_dim_name, "edges", 6))
                        dimsize = gio_param->clus_edges;
                }
            }
        }

        err = gio_def_dim(gio_param, gio_files[nftype].fileID, tmp_dim_name,
                          dimsize, &dimID);
        if (err != NC_NOERR) {
            char err_msg[128];
            sprintf(err_msg, "defining dimension: %s", tmp_dim_name);
            NC_CHECK(err, err_msg, 0)
        }
        gio_files[nftype].dimension_ids[i] = dimID;
    }
}


/*----< gio_add_field_descriptor() >------------------------------------------*/
/* Create field in a specific fileID
   The information on dimensions is held in the
   gio_descriptors[[nfield]].field_dims[] array
   This is a list of indices that point to the information stored in the
   gio_files[ifile].dimension_ids[:] array
*/
static
void gio_add_field_descriptor(gio_parameters *gio_param,
                              int             ifile,
                              int             nfield)
{
    int idesc;  /* Index into gio_descriptors structure
                   corresponding to nfield */
    int fieldID, nd, idx, numdims;
    int the_dims[MAX_DIMS], dimids[MAX_DIMS];  /* dimension IDs  */
    int dim_sizes[MAX_DIMS]; /* dimension values  */
    int nattr, err;
    char *attr_name, *attr_value, *tmp_dim_name, err_msg[128];
    int chunks[MAX_DIMS];
    int coarsen_factor;
    int distdim, bytesize=4;
    nc_type  vartype;

    data_descriptor *dscptr;
    file_descriptor *gio_files = &gio_param->file.gio_files[ifile];

    idesc = gio_files->fields[nfield];
    dscptr = &gio_param->data.gio_descriptors[idesc];

    numdims = dscptr->num_file_dims;
    distdim = 0;

    for (nd=0; nd<numdims; nd++) {
        idx = gio_files->field_dim_ids[nd][nfield];
        the_dims[nd]  = gio_files->dimension_ids[idx];
        dim_sizes[nd] = gio_files->dimension_sizes[idx];
        tmp_dim_name  = gio_files->dimension_names[idx];

        /* Fix up dimension size for coarsened files, if necessary */
        if (gio_files->grid_level != gio_param->grid.refine_param &&
            (!strncmp(tmp_dim_name, "cells",   5) ||
             !strncmp(tmp_dim_name, "corners", 7) ||
             !strncmp(tmp_dim_name, "edges",   6))) {
            coarsen_factor = gio_param->grid.refine_param
                           - gio_files->grid_level;
            coarsen_factor = POWER2(2*coarsen_factor);
            if (!strncmp(tmp_dim_name, "cells", 5))
                dim_sizes[nd] = (dim_sizes[nd]-2)/coarsen_factor + 2;
            else
                dim_sizes[nd] = dim_sizes[nd]/coarsen_factor;
        }
    }
    for (nd=0; nd<numdims; nd++)  /* wkliao: F2C dimension order */
        dimids[nd] = the_dims[numdims-nd-1];

    if (dscptr->data_type == flt_type)
        vartype = NC_FLOAT;
    else if (dscptr->data_type == int_type)
        vartype = NC_INT;
    else if (dscptr->data_type == dbl_type) {
        vartype = NC_DOUBLE;
        bytesize = 8;
    }
    else
        NC_CHECK(-1, "Insupported variable type", 0)

    /* pcell is never set for pnetcdf */
    if (distdim != 0)
       err = gio_def_var(gio_param, gio_files->fileID,
                         dscptr->output_field_name,
                         vartype, numdims, dimids, &fieldID, chunks);
    else
       err = gio_def_var(gio_param, gio_files->fileID,
                         dscptr->output_field_name,
                         vartype, numdims, dimids, &fieldID, NULL);

    if (err != NC_NOERR) {
        sprintf(err_msg, "defining: %s", dscptr->output_field_name);
        NC_CHECK(err, err_msg, 0)
    }

    /* Then define a set of standard field attributes (long_name, units,
       standard_name, etc.) */
    err = gio_put_att_text(gio_param, gio_files->fileID, fieldID,
                           "long_name", dscptr->long_name);
    if (err != NC_NOERR) {
        sprintf(err_msg, "adding attribute: %s", dscptr->long_name);
        NC_CHECK(err, err_msg, 0)
    }
    err = gio_put_att_text(gio_param, gio_files->fileID, fieldID,
                           "units", dscptr->units);
    if (err != NC_NOERR) {
        sprintf(err_msg, "adding attribute: %s", dscptr->units);
        NC_CHECK(err, err_msg, 0)
    }
    err = gio_put_att_text(gio_param, gio_files->fileID, fieldID,
                           "standard_name", dscptr->standard_name);
    if (err != NC_NOERR) {
        sprintf(err_msg, "adding attribute: %s", dscptr->standard_name);
        NC_CHECK(err, err_msg, 0)
    }

    /* Write all the attributes as provided by the config files */
    for (nattr=0; nattr<dscptr->num_attrs; nattr++) {
        attr_name  = dscptr->field_attr_name[nattr];
        attr_value = dscptr->field_attr_value[nattr];
        err = gio_put_att_text(gio_param, gio_files->fileID, fieldID,
                               attr_name, attr_value);
        if (err != NC_NOERR) {
            sprintf(err_msg, "adding attribute: %s", attr_name);
            NC_CHECK(err, err_msg, 0)
        }
    }
}


/*----< gio_update_time_field() >---------------------------------------------*/
void gio_update_time_field(gio_parameters *gio_param,
                           int             ifile,
                           double         *the_time) /* time in seconds */
{
    int err, timevid;
    double starttime;
    MPI_Offset start, count;

    if (!gio_param->gio_do_io) return;

    starttime = MPI_Wtime();
    file_descriptor *gio_files = gio_param->file.gio_files;

    err = gio_inq_varid(gio_param, gio_files[ifile].fileID, "time", &timevid);
    NC_CHECK(err,"getting varid: time", 0)

    count = 1;
    start = gio_files[ifile].samples_written;
    /* Time is tricky because its a replicated scalar
       We write 0 data except on the head node for pnetcdf
       Bug in nc4/hdf prevents use of this strategy for pnetcdf */
    if (gio_param->gio_me != 0) count = 0;

    if (gio_param->gio_iotype == MODE_BLOCKING_INDEP) {
        err = gio_put_double(gio_param, gio_files[ifile].fileID,
                             timevid, &start, &count, the_time);
        NC_CHECK(err, "writing time.", 0)
    } else {
        if (gio_param->gio_me == 0) {
            err = gio_iput_double(gio_param, gio_files[ifile].fileID,
                                  timevid, &start, &count, the_time);
            NC_CHECK(err, "nonblocking writing time.", 0)
        }
    }

    /* Not totally sure if this should be counted as a record */
    gio_param->num_records_written++;

    /* only count the bytes on proc 0 */
    if (gio_param->gio_me == 0) gio_param->bytes_API_write += sizeof(double);

    // gio_stat_writes(bytes);
    gio_param->time_in_update_time += MPI_Wtime() - starttime;
}

/*----< gio_set_property() >--------------------------------------------------*/
/* Add a global attribute (aka "property") to fileID */
static
void gio_set_property(gio_parameters *gio_param,
                      gio_prop       *props,
                      int             fileID,
                      int             nprops,
                      MPI_Comm        comm)
{
    int err;

    if (nprops < 0) {
        /* put history attribute in file */
        char *gio_timestamp = props->gio_timestamp;
        struct timeval tv;
        time_t curtime;

        gettimeofday(&tv, NULL);
        curtime=tv.tv_sec;
        strftime(gio_timestamp, 128, "created on: %m-%d-%Y at %T",
                 localtime(&curtime));

        /* Ensure that all processors have the same date value
           wkliao: a temporary solution */
        if (comm != MPI_COMM_SELF)
            MPI_Bcast(gio_timestamp, 128, MPI_CHAR, 0, comm);

#define USE_FIXED_TIME_STAMP
#ifdef USE_FIXED_TIME_STAMP
        err = gio_put_att_text(gio_param, fileID, NC_GLOBAL, "history",
                               "created on: 08-03-2017 at 16:33:00");
#else
        err = gio_put_att_text(gio_param, fileID, NC_GLOBAL, "history",
                               gio_timestamp);
#endif
        NC_CHECK(err, "adding attribute: history", 0)
    }
    else {
        prop_descriptor *gio_properties = props->gio_properties;

        err = gio_put_att_text(gio_param, fileID, NC_GLOBAL,
                               gio_properties[nprops].prop_name,
                               gio_properties[nprops].prop_value);
        if (err != NC_NOERR) {
            char err_msg[128];
            sprintf(err_msg, "adding attribute: %s",
                    gio_properties[nprops].prop_name);
            NC_CHECK(err, err_msg, 0)
        }
    }
}


/*----< gio_add_verticalgrid() >----------------------------------------------*/
/* Write vertical grid data which exists on every process.
   This data is flagged as has_levels, presumed to be single precision,
*/
static
void gio_add_verticalgrid(gio_parameters *gio_param,
                          int             ifile)
{
    /* note that nonblocking method does not call this function.
       it calls gio_write_verticalgrid() instead */
    char *fldname;
    int i, j, err, fieldID;
    double starttime;
    MPI_Offset start, count;

    starttime = MPI_Wtime();
    gio_data *data = &gio_param->data;
    file_descriptor *gio_files = &gio_param->file.gio_files[ifile];
    data_descriptor *dscptrs = gio_param->data.gio_descriptors;

    for (i=0; i<gio_files->nflds; i++) {
        j = gio_files->fields[i];  /* Index into descriptors array */

        if (dscptrs[j].is_level_data) {

            /* only root writes level variables (layers or interfaces) */
            if (gio_param->gio_me > 0) continue;

            fldname = dscptrs[j].output_field_name;
            if (strcmp(fldname, "layers") == 0)
                count = data->layers_length;
            else if (strcmp(fldname, "interfaces") == 0)
                count = data->interfaces_length;

            /* Hmm if we don't create from proc 0, then all procs are writing
               the same data. how do the underlying libs work with this */
            start = 0;

            err = gio_inq_varid(gio_param, gio_files->fileID, fldname,
                                &fieldID);
            if (err != NC_NOERR) {
                char err_msg[128];
                sprintf(err_msg, "getting varid: %s",fldname);
                NC_CHECK(err, err_msg, 0)
            }

            if (gio_param->gio_iotype == MODE_BLOCKING_INDEP)
                err = gio_put_real(gio_param, gio_files->fileID, fieldID, &start,
                                   &count, dscptrs[j].float_level_data);
            else
                err = gio_iput_real(gio_param, gio_files->fileID, fieldID, &start,
                                    &count, dscptrs[j].float_level_data);
            if (err != NC_NOERR) {
                char err_msg[128];
                sprintf(err_msg, "put real: %s\n",fldname);
                NC_CHECK(err, err_msg, 0)
            }

            gio_param->bytes_API_write += count * sizeof(float);
            gio_param->num_records_written++;
            /* call gio_stat_writes(4*count1d(1)) */
        }
    }
    gio_param->time_in_nf_put_var_grid += MPI_Wtime() - starttime;
}


/*----< gio_write_grid() >----------------------------------------------------*/
/* Write distributed grid information which excludes any grid
   info that is known on all processors.
*/
static
void gio_write_grid(gio_parameters *gio_param,
                    int             ifile)
{
    int nfield;
    double starttime, endtime;

    starttime = MPI_Wtime();
    file_descriptor *gio_files = &gio_param->file.gio_files[ifile];

    /* allocate the nonblocking I/O request array for all the grid variables
       to be written to this file */
    if (gio_param->gio_iotype != MODE_BLOCKING_INDEP)
        gio_allocate_req_ids(gio_param, ifile, 1, -1);

    for (nfield=0; nfield<gio_files->nflds; nfield++) {
        /* check whether we are writing this grid */
        if (! gio_is_output_grid(gio_param, ifile, nfield)) continue;

        // if (gio_param->gio_do_nonblocking_io)
            /* All processes post nonblocking writes */
            gio_write_field(gio_param, ifile, gio_files->fields[nfield]);

/* this benchmark only works for nonblocking I/O option. The call to
   gio_exchng_and_write() is an option available on Fortran's GIO library
        else
            gio_exchng_and_write(gio_files->fields[nfield], ifile);
*/
    }

    /* wait for all nonblocking requests to complete */
    if (gio_param->gio_iotype != MODE_BLOCKING_INDEP)
        gio_nonblocking_io_wait(gio_param, ifile);
    /* free the allocated buffers */
    if (gio_param->using_interleaved == 0 && gio_param->using_direct == 0)
        gio_deallocate_io_buffer(gio_param, ifile, 1, -1);
    gio_free_req_ids(gio_param);

    endtime = MPI_Wtime();
    gio_param->time_in_nf_put_var_grid += endtime - starttime;
}


/*#############################################################################
!# Restart related special code
!#############################################################################*/

/*----< gio_open_restart_read() >---------------------------------------------*/
int gio_open_restart_read(gio_parameters *gio_param,
                          int             ifile,
                          char           *fname)
{
    int err;

    /* Note that all processors read restart files */
    err = gio_open(gio_param, gio_param->gio_world_comm, fname, NC_NOWRITE,
                   &gio_param->used_info,
                   &gio_param->file.gio_files[ifile].fileID);
    if (err != NC_NOERR) {
        char err_msg[128];
        sprintf(err_msg, "opening restart file for read: %s", fname);
        NC_CHECK(err, err_msg, 0)
    }

    gio_read_restart_attrs(gio_param, ifile);

    return 1;
}


/*----< gio_write_restart_attrs() >-------------------------------------------*/
static
void gio_write_restart_attrs(gio_parameters *gio_param,
                             int             ifile)
{
    int i, err;
    char tmp_name[128];

    gio_dim  *dim  = &gio_param->dim;
    gio_data *data = &gio_param->data;

    int fileid = gio_param->file.gio_files[ifile].fileID;
    int fieldid = NC_GLOBAL;
    for (i=0; i<dim->num_index_arrays; i++) {
        err = gio_put_att_int(gio_param, fileid, fieldid, dim->index_names[i],
                              dim->index_sizes[i]);
        NC_CHECK(err, "adding int attribute", 0)
        sprintf(tmp_name, "%s_index", dim->index_names[i]);
        err = gio_put_att_intarr(gio_param, fileid, fieldid, tmp_name,
                                 dim->index_arrays[i].p, dim->index_sizes[i]);
        NC_CHECK(err, "adding int array attribute", 0)
    }
    for (i=0; i<data->num_restart_dbl; i++) {
        double dval = data->restart_dbl[i];
        err = gio_put_att_double(gio_param, fileid, fieldid,
                                 data->restart_dbl_names[i], dval);
        NC_CHECK(err, "adding double attribute", 0)
    }
    for (i=0; i<data->num_restart_int; i++) {
        int ival = data->restart_int[i];
        err = gio_put_att_int(gio_param, fileid, fieldid,
                              data->restart_int_names[i], ival);
        NC_CHECK(err, "adding double attribute", 0)
    }
}

/*----< gio_read_restart_attrs() >--------------------------------------------*/
void gio_read_restart_attrs(gio_parameters *gio_param,
                            int             ifile)
{
    int i, err;
    int fileid, fieldid, ival;
    char tmp_name[128];
    double dval;

    gio_data *data = &gio_param->data;
    gio_dim  *dim  = &gio_param->dim;

    fileid = gio_param->file.gio_files[ifile].fileID;
    fieldid = NC_GLOBAL;

    for (i=0; i<data->num_restart_int; i++) {
        err = gio_get_att_int(gio_param, fileid, fieldid, dim->index_names[i],
                              &dim->index_sizes[i]);
        NC_CHECK(err, "reading integer attribute", 0)

        sprintf(tmp_name, "%s_index", dim->index_names[i]);
        err = gio_get_att_intarr(gio_param, fileid, fieldid, tmp_name,
                                 dim->index_arrays[i].p, dim->index_sizes[i]);
        NC_CHECK(err, "reading integer array attribute", 0)
    }
    for (i=0; i<data->num_restart_dbl; i++) {
        err = gio_get_att_double(gio_param, fileid, fieldid,
                                 data->restart_dbl_names[i], &dval);
        NC_CHECK(err, "reading double attribute", 0)
        data->restart_dbl[i] = dval;
    }
    for (i=0; i<data->num_restart_int; i++) {
        err = gio_get_att_int(gio_param, fileid, fieldid,
                              data->restart_int_names[i], &ival);
        NC_CHECK(err, "reading double attribute", 0)
        data->restart_int[i] = ival;
    }
}


/*----< gio_read_time_field() >-----------------------------------------------*/
/* Read the latest time field */
void gio_read_time_field(gio_parameters *gio_param,
                         int             ifile,
                         double         *the_time,
                         MPI_Offset     *nsteps)
{
    /* the_time (time in seconds) ... only allow 1 time in restarts */

    int err, timevid, time_dimid;
    MPI_Offset count, start;
    double starttime = MPI_Wtime();

    file_descriptor *gio_files = &gio_param->file.gio_files[ifile];

    err = gio_inq_varid(gio_param, gio_files->fileID, "time", &timevid);
    NC_CHECK(err,"getting varid: time", 0)

    /* variable time is a 1D array (UNLIMITED) */
    err = ncmpi_inq_vardimid(gio_files->fileID, timevid, &time_dimid);
    NC_CHECK(err, "getting dimension ID", 0)

    err = gio_inq_dimlen(gio_param, gio_files->fileID, time_dimid, nsteps);
    NC_CHECK(err, "getting dimension len", 0)

    /* Have all the procs read the time. This will likely be slow...
       Possible optimization -  read once and broadcast instead?? */
    start = *nsteps;  /* read the latest time step */
    count = 1;

    err = gio_get_double(gio_param, gio_files->fileID, timevid, &start,
                         &count, the_time);
    NC_CHECK(err, "getting time", 0)

    gio_param->time_in_update_time += MPI_Wtime() - starttime;
}

/*----< gio_coarsen_field() >-------------------------------------------------*/
/* Apply coarsening operation to field, if necessary.  */
static
void gio_coarsen_field(gio_parameters *gio_param,
                       int             field_idx,
                       int             ifile,
                       MPI_Offset      offset,
                       MPI_Offset     *start,
                       MPI_Offset     *count,
                       void           *buffer,
                       MPI_Offset      cellsize)
{
    int j, idx, jdx, iavg;
    int celldim, coarsen_factor;
    int offset_pole, offset_size, iloop_dim;
    int cellnum, base_cell, iloop, icell;
    int average;

    gio_grid        *grid            = &gio_param->grid;
    file_descriptor *gio_files       = &gio_param->file.gio_files[ifile];
    data_descriptor *dscptr = &gio_param->data.gio_descriptors[field_idx];

    celldim = dscptr->cell_idx;

    /* Modify data if cell is coarsened. */
    if (gio_files->grid_level >= gio_param->grid.refine_param)
        return;

    coarsen_factor = gio_param->grid.refine_param - gio_files->grid_level;
    coarsen_factor = POWER2(2*coarsen_factor);
    average = gio_files->crsn_average;
    if (dscptr->data_type == int_type) average = 0;

    /* Check if cell has poles */
    if (dscptr->no_poles)
        offset_pole = 0;
    else {
        if (offset == 0)
            offset_pole = 2;
        else
            offset_pole = 0;
    }
    offset_size = offset_pole*cellsize;
    cellnum = count[celldim];
    idx = 0;
    jdx = 0;

    /* Cell-centered fields */

    if (! dscptr->no_poles) {
        if (average) {
#define AVERAGING(type, buf) {                                           \
    type *buf_ptr = buf;                                                 \
    type *buf_avg = (type*) calloc(cellsize, sizeof(type));              \
    idx = offset_size;                                                   \
    for (icell=offset_pole; icell<cellnum; icell+=coarsen_factor) {      \
        for (iavg=0; iavg<coarsen_factor; iavg++) {                      \
            jdx = (icell+iavg-offset_pole)*cellsize + offset_size;       \
            for (j=0; j<cellsize; j++)                                   \
                buf_avg[j] += buf_ptr[jdx++];                            \
        }                                                                \
        for (j=0; j<cellsize; j++)                                       \
            buf_ptr[idx++] = buf_avg[j]/coarsen_factor;                  \
    }                                                                    \
    free(buf_avg);                                                       \
}
            if (dscptr->data_type == flt_type)
                AVERAGING(float,  buffer)
            else
                AVERAGING(double, buffer)
        }
        else {
#define BUFFER_COPY(type, buf) {                                         \
    type *buf_ptr = buf;                                                 \
    idx = offset_size;                                                   \
    for (icell=offset_pole; icell<cellnum; icell+=coarsen_factor) {      \
        jdx = (icell-offset_pole)*cellsize + offset_size;                \
        memmove(buf_ptr+idx, buf_ptr+jdx, cellsize*sizeof(type));        \
        idx += cellsize;                                                 \
    }                                                                    \
}
            if (dscptr->data_type == flt_type)
                BUFFER_COPY(float,  buffer)
            else if (dscptr->data_type == int_type)
                BUFFER_COPY(int,    buffer)
            else if (dscptr->data_type == dbl_type)
                BUFFER_COPY(double, buffer)
        }
        count[celldim] = (cellnum-offset_pole)/coarsen_factor + offset_pole;
        start[celldim] = start[celldim]-1;
        if (start[celldim] != 0)
            start[celldim] = (start[celldim]-2)/coarsen_factor + 2;
        start[celldim]++;
    }
    else {
        /* Corner or edge-centered fields */
        iloop_dim = 0;
        if (dscptr->grid_size == grid->gio_edges_size) {
            iloop_dim = 3;
            base_cell = cellnum/3;
        }
        else if (dscptr->grid_size == grid->gio_corners_size) {
            iloop_dim = 2;
            base_cell = cellnum/2;
        }
        if (average) {
#define AVERAGING2(type, dtype, buf) {                                        \
    type  *buf_ptr = buf;                                                     \
    type **buf_avg = (type**) calloc_2D_##dtype(iloop_dim, cellsize);         \
    for (icell=0; icell<base_cell; icell+=coarsen_factor) {                   \
        for (iavg=0; iavg<coarsen_factor; iavg++) {                           \
            for (iloop=0; iloop<iloop_dim; iloop++) {                         \
                jdx = iloop*cellsize                                          \
                    + (icell+iavg-2)*iloop_dim*base_cell + offset_size;       \
                for (j=0; j<cellsize; j++)                                    \
                    buf_avg[iloop][j] += buf_ptr[jdx++];                      \
            }                                                                 \
        }                                                                     \
        idx = offset_size;                                                    \
        for (iloop=0; iloop<iloop_dim; iloop++) {                             \
            for (j=0; j<cellsize; j++)                                        \
                buf_ptr[idx++] = buf_avg[iloop][j]/coarsen_factor;            \
        }                                                                     \
    }                                                                         \
    free_2D_##dtype(buf_avg);                                                 \
}

            if (dscptr->data_type == flt_type)
                AVERAGING2(float,  flt, buffer)
            else
                AVERAGING2(double, dbl, buffer)
        }
        else {
#define BUFFER_COPY2(type, buf) {                                           \
    type *buf_ptr = buf;                                                    \
    idx = offset_size;                                                      \
    for (icell=0; icell<base_cell; icell+=coarsen_factor) {                 \
        jdx = icell*iloop_dim*base_cell + offset_size;                      \
        for (iloop=0; iloop<iloop_dim; iloop++) {                           \
            memmove(buf_ptr+idx, buf_ptr+jdx, cellsize*sizeof(type));       \
            jdx += cellsize;                                                \
            idx += cellsize;                                                \
        }                                                                   \
    }                                                                       \
}
            if (dscptr->data_type == flt_type)
                BUFFER_COPY2(float,  buffer)
            else if (dscptr->data_type == int_type)
                BUFFER_COPY2(int,    buffer)
            else if (dscptr->data_type == dbl_type)
                BUFFER_COPY2(double, buffer)
        }
        count[celldim] = cellnum/coarsen_factor;
        start[celldim] = (start[celldim]-1) / coarsen_factor;
    }
}


