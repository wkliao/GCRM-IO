/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_params_time.h 4603 2017-12-07 07:12:52Z wkliao $
 */

#ifndef H_ZGRD_PARAMS_TIME
#define H_ZGRD_PARAMS_TIME

typedef struct {

    double dt_ZGrd; /* timestep of the model */

    int  nprog; /* = 2 number of time levels stored for each prognostic */
    int  ntend; /* = 3 number of time levels stored for each tendency */

    int  prog_index[2]; /* [nprog] */
    int  tend_index[3]; /* [ntend] */
    int  np0, np1, nm0, nm1, nm2;

    double tend_weights[3]; /* [ntend] */

    int step_count_ZGrd; /* counts the number of timesteps performed by the model */
    int step_count_rad;  /* the ratio of the radiation time step to 
                            the basic timestep */
    int step_count_sfx;  /* the ratio of the surface flux time step to  */

    int end_days;
    int end_hours;
    int end_minutes;
    int end_seconds;

    double time_start_ZGrd; /* time at the beginning of the simulation */
    double time_ZGrd;       /* current time in the atmosphere model
                               time_ZGrd = dt_ZGrd*FLOAT (step_count_ZGrd-1) */
    double time_end_ZGrd;   /* time at the end of the simulation [s] */

    /* These variables are specified as quotients, first element is numerator,
       second is denominator */
    double idt_output[2];
    double idt_restart[2];
    double idt_ZGrd[2];
    double itime_ZGrd[2];
    double itime_start_ZGrd[2];
    double itime_end_ZGrd[2];

    char time_unit_strng[128];

} MODULE_ZGrd_params_time;

void initialize_params_time(const char *fname, MODULE_ZGrd_params_time *time, int level_max);

#endif
