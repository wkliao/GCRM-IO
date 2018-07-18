/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio_sched.c 4405 2017-08-19 23:05:34Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <gcrm.h>
#include <gio.h>


/*----< gio_init_schedule() >-------------------------------------------------*/
void gio_init_schedule(gio_parameters *gio_param)
{
    int i;

    if (! gio_param->is_schedule_set) {
        if (gio_param->using_interleaved == 0 && gio_param->using_direct == 0)
            gio_param->gio_io_comm = MPI_COMM_WORLD;  /* all processes do I/O */
        /* this code only performs PnetCDF blocking and nonblocking I/O option */

        /* Need to finish allocating memory for averages */
        for (i=0; i<gio_param->data.num_averages; i++)
            gio_allct_avg_descriptor(gio_param,
                                     gio_param->data.avg_field_ids[i]);

        gio_set_strides(gio_param);
        gio_check_file_data(gio_param);

        gio_param->is_schedule_set = 1;
    }
}
