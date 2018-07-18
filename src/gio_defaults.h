/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio_defaults.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_GIO_DEFAULTS
#define H_GIO_DEFAULTS

typedef struct {
/*
!  This module contains API defaults to make
!  construction of the atmos.fcfg, gcrm.fcfg much easier.
!
!  The default values are read from atmos.in, gcrm.in, etc. and
!  pased to gio_init via the command line.
!
*/
  char gio_default_cdf_path[256];
  char gio_sep_grid_file[256];
  int  gio_default_frequency;  /* in seconds */
  int  gio_default_nsamples;   /* default nsamples per file */
  int  gio_grid_option;        /* 1=nogrid, 2=allfiles get grid, 3=grid */
} gio_defaults;

#endif
