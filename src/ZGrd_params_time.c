/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_params_time.c 4605 2017-12-07 07:17:10Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <errno.h>
extern int errno;

#include "ZGrd_params_time.h"
#include "util.h"


/*----< set_tendency_weights() >----------------------------------------------*/
/* defined in MODULE ZGrd_params_time in ZGrd_params_time.F90 */
static
void set_tendency_weights(int     step_count,
                          double *weights)    /* [ntend] */
{
    switch (step_count) {
        case 1:  weights[0] = 1.0;
                 weights[1] = 0.0;
                 weights[2] = 0.0;
                 break;
        case 2:  weights[0] = 1.5;
                 weights[1] = 0.5;
                 weights[2] = 0.0;
                 break;
        default: weights[0] =  23.0 / 12.0;
                 weights[1] = -16.0 / 12.0;
                 weights[2] =   5.0 / 12.0;
                 break;
     }
}

/*----< initialize_params_time() >--------------------------------------------*/
void initialize_params_time(const char              *fname,
                            MODULE_ZGrd_params_time *time,
                            int                      level_max)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* parameters from MODULE ZGrd_params_time */
    time->nprog       = 2;
    time->ntend       = 3;
    time->end_days    = 12;
    time->end_hours   = 0;
    time->end_minutes = 0;
    time->end_seconds = 0;
    strcpy(time->time_unit_strng, "h");

    /* the following is set in initialize_params_time() */
    /* time->idt_ZGrd[0] is the amount of time advanced in each simulation
       iteration */
    time->idt_ZGrd[1] = 1;  /* all 2nd elements are set to 1 */
    switch (level_max) {
        case 3: /* 0000000642 */
                time->idt_ZGrd[0] = 900;
                time->step_count_rad= 4;
                time->step_count_sfx= 1;
                break;
        case 4: /* 0000002562 */
                time->idt_ZGrd[0] = 900;
                time->idt_ZGrd[0] = 120;
                time->step_count_rad= 4;
                time->step_count_sfx= 1;
                break;
        case 5: /* 0000010242 */
                time->idt_ZGrd[0] = 600;
                time->idt_ZGrd[0] =  60;
                time->step_count_rad= 6;
                time->step_count_sfx= 1;
                break;
        case 6: /* 0000040962 */
                time->idt_ZGrd[0] = 300;
                time->idt_ZGrd[0] =  30;
                time->step_count_rad=12;
                time->step_count_sfx= 1;
                break;
        case 7: /* 0000163842 */
                time->idt_ZGrd[0] = 180;
                time->idt_ZGrd[0] =  15;
                time->step_count_rad=20;
                time->step_count_sfx= 1;
                break;
        case 8: /* 0000655362 */
                time->idt_ZGrd[0] = 90;
                time->idt_ZGrd[0] =  8;
                time->step_count_rad=20;
                time->step_count_sfx= 1;
                break;
        case 9: /* 0002621442 */
                time->idt_ZGrd[0] = 50;
                time->idt_ZGrd[0] =  4;
                time->step_count_rad=18;
                time->step_count_sfx= 1;
                break;
        case 10: /* 0010485762 */
                time->idt_ZGrd[0] = 24;
                time->idt_ZGrd[0] =  2;
                time->step_count_rad=15;
                time->step_count_sfx= 2;
                break;
        case 11: /* 0041943042 */
                time->idt_ZGrd[0] = 12;
                time->idt_ZGrd[0] =  1;
                time->step_count_rad=20;
                time->step_count_sfx= 5;
                break;
        case 12: /* 0167772162 */
                time->idt_ZGrd[0] = 6;
                time->idt_ZGrd[0] = 1;
                time->idt_ZGrd[1] = 2;
                time->step_count_rad=20;
                time->step_count_sfx= 10;
                break;
        case 13: /* 0671088642 */
                time->idt_ZGrd[0] = 3;
                time->step_count_rad=20;
                time->step_count_sfx= 20;
                break;
        default:
                printf("initialize_params_time :: CANNOT DETERMINE dt_ZGrd\n");
                break;
    }

    /*  KLS Time doesn't matter in this stripped down version so set
        them to 60 for reasonable speed
    time->idt_ZGrd[0] = 60;
    time->idt_ZGrd[1] =  1;
    */

    time->dt_ZGrd = time->idt_ZGrd[0] / time->idt_ZGrd[1];

    int int_msg[4];
    if (rank == 0) { /* only root process reads from the file */
        char line[256], *value, *key;
        FILE *fp;

        OPEN_FILE(fname)

        while (fgets(line, 256, fp) != NULL) {
            GET_LINE_TOKEN_PAIR

            if (!strcasecmp(key, "end_days"))
                time->end_days = atoi(value);
            else if (!strcasecmp(key, "end_hours"))
                time->end_hours = atoi(value);
            else if (!strcasecmp(key, "end_minutes"))
                time->end_minutes = atoi(value);
            else if (!strcasecmp(key, "end_seconds"))
                time->end_seconds = atoi(value);

        }
        fclose(fp);

        int_msg[0] = time->end_days;
        int_msg[1] = time->end_hours;
        int_msg[2] = time->end_minutes;
        int_msg[3] = time->end_seconds;
    }
    MPI_Bcast(int_msg, 4, MPI_INT, 0, MPI_COMM_WORLD);
    time->end_days    = int_msg[0];
    time->end_hours   = int_msg[1];
    time->end_minutes = int_msg[2];
    time->end_seconds = int_msg[3];

    time->itime_start_ZGrd[0] = 0.0;
    time->itime_end_ZGrd[0]   =     time->end_seconds
                                    + 60*(time->end_minutes
                                    + 60*(time->end_hours
                                    + 24* time->end_days));

    /* give all integer time and timestep variables the denominator of
       idt_ZGrd and adjust the numerator */
    time->itime_start_ZGrd[1] = time->idt_ZGrd[1];
    time->itime_start_ZGrd[0] = time->itime_start_ZGrd[0]*time->idt_ZGrd[1];
    time->itime_end_ZGrd[1]   = time->idt_ZGrd[1];
    time->itime_end_ZGrd[0]   = time->itime_end_ZGrd[0]*time->idt_ZGrd[1];
    time->idt_output[1]       = time->idt_ZGrd[1];
    time->idt_output[0]       = time->idt_output[0] * time->idt_ZGrd[1];
    time->idt_restart[1]      = time->idt_ZGrd[1];
    time->idt_restart[0]      = time->idt_restart[0] * time->idt_ZGrd[1];
    time->itime_ZGrd[1]       = time->idt_ZGrd[1];
    time->step_count_ZGrd     = 1;

    int i;
    for (i=0; i<time->nprog; i++)
        time->prog_index[i] = i;
    time->np0 = time->prog_index[0];
    time->np1 = time->prog_index[1];

    for (i=0; i<time->ntend; i++)
        time->tend_index[i] = i;
    time->nm0 = time->tend_index[0];
    time->nm1 = time->tend_index[1];
    time->nm2 = time->tend_index[2];

    set_tendency_weights(time->step_count_ZGrd, time->tend_weights);
}

