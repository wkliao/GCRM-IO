/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: gio_utility.c 1858 2013-06-12 04:49:45Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "gcrm.h"
#include "gio.h"

/* ---------------------------------------------------------------------------
! Convert numeric time to time string for output
! ---------------------------------------------------------------------------*/
void gio_convert_time(double  model_time,
                      char   *time_str)
{
    int seconds_rem;
    int secs_in_yr;
    int seconds_rem2;
    int day, year, month, hh, mm, ss;

    seconds_rem  = model_time;
    secs_in_yr   = 365 * 86400;
    year         = seconds_rem / secs_in_yr + 1; /* NOTE: Copied logic from BUGS5, which assumes no LEAP YEARS */
    year         = year + 1900;                  /* Add 1900 years to the computation */

    seconds_rem  = seconds_rem % secs_in_yr;
    seconds_rem2 = seconds_rem;                  /* Use as an integer now */
    day          = seconds_rem2 / 86400 + 1;

         if (day <=  31) { month =  1; }
    else if (day <=  59) { month =  2; day -=  31; }
    else if (day <=  90) { month =  3; day -=  59; }
    else if (day <= 120) { month =  4; day -=  90; }
    else if (day <= 151) { month =  5; day -= 120; }
    else if (day <= 181) { month =  6; day -= 151; }
    else if (day <= 212) { month =  7; day -= 181; }
    else if (day <= 243) { month =  8; day -= 212; }
    else if (day <= 273) { month =  9; day -= 243; }
    else if (day <= 304) { month = 10; day -= 243; }
    else if (day <= 334) { month = 11; day -= 304; }
    else                 { month = 12; day -= 334; }

    seconds_rem2 = seconds_rem2 % 86400;
    hh           = seconds_rem2 / 3600;
    seconds_rem2 = seconds_rem2 % 3600;
    mm           = seconds_rem2 / 60;
    seconds_rem2 = seconds_rem2 % 60;
    ss           = seconds_rem2;

    /* Construct a date/time string from this info */
    sprintf(time_str, "%d%02d%02d_%02d%02d%02d", year,month,day,hh,mm,ss);
}

