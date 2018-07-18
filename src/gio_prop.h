/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio_prop.h 4600 2017-12-07 07:03:41Z wkliao $
 */

#ifndef H_GIO_PROP
#define H_GIO_PROP

typedef struct {
    char prop_name[32];
    char prop_value[128];
} prop_descriptor;

typedef struct {
/*
 *  This module contains all the global properties (attribues in netcdf-speak)
 */

    /*  number of properties defined in configuration file */
    int             num_props;
    prop_descriptor gio_properties[MAX_PROPS];

    char gio_datestr[8];
    char gio_timestr[16];
    char gio_timestamp[128];
} gio_prop;

#endif
