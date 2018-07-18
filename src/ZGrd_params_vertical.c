/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_params_vertical.c 1858 2013-06-12 04:49:45Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "gcrm.h"
#include "gio.h"
#include "ZGrd_params.h"
#include "util.h"

/* from ZGrd_params_vertical.F90 */

/*----< set_vertical_coordinate_levels() >------------------------------------*/
void set_vertical_coordinate_levels(int     km,
                                    double  z_top,
                                    double *z,     /* [km+1] */
                                    double *z_lyr) /* [km] */
{
    int k, select_coordinate_levels = 1;
    double dz_ratio = 5.0;
    double dz_0, dz_low, c2, c1;
    double *z_0, *z_lyr_0;

    z_0     = (double*) malloc((km+1) * sizeof(double));
    z_lyr_0 = (double*) malloc(km     * sizeof(double));

    if (km == 1) {
            z[0] = 0.0;
        z_lyr[0] = z_top / 2.0;
            z[1] = z_top;
    }
    else {
        switch(select_coordinate_levels) {
            case (1):
                for (k=0; k<km+1; k++)
                    z[k] = z_top * k / km;
                for (k=0; k<km; k++)
                    z_lyr[k] = (z[k]+z[k+1]) / 2.0;
                break;
            case (2):
                dz_0   = z_top / km;
                dz_low = dz_0 / dz_ratio;

                c2 = (dz_0 - dz_low) / (dz_0 * (z_top - dz_0));
                c1 = 1.0 - c2 * z_top;

                for (k=0; k<km+1; k++)
                   z_0[k] = dz_0 * k;

                for (k=0; k<km; k++)
                   z_lyr_0[k] = z_0[k] + dz_0 / 2.0;

                for (k=0; k<km+1; k++)
                   z[k] = (c1 + c2 * z_0[k]) * z_0[k];

                for (k=0; k<km; k++)
                  z_lyr[k] = (c1 + c2 * z_lyr_0[k]) * z_lyr_0[k];
                break;
            default: break;
        }
    }
    free(z_0);
    free(z_lyr_0);
}

/*----< initialize_params_vertical() >----------------------------------------*/
void initialize_params_vertical(MODULE_ZGrd_params_vertical *vertical)
{
    int k;

    vertical->km      = 256; /* number of vertical layers in the  atmosphere */
    vertical->z_top   = 18000.0; /* model top [m] */
    vertical->prs_top = 20.0; /* pressure at the model top [Pa] */

    memset(vertical->zeta,           0, 257 * sizeof(double));
    memset(vertical->zeta_lyr,       0, 256 * sizeof(double));
    memset(vertical->d_zeta,         0, 256 * sizeof(double));
    memset(vertical->d_zeta_ifc,     0, 257 * sizeof(double));
    memset(vertical->d_zeta_inv,     0, 256 * sizeof(double));
    memset(vertical->d_zeta_inv_ifc, 0, 257 * sizeof(double));

    set_vertical_coordinate_levels(vertical->km,
                                   vertical->z_top,
                                   vertical->z,
                                   vertical->z_lyr);

    for (k=0; k<vertical->km; k++)
        vertical->dz[k] = vertical->z[k] - vertical->z[k-1];

    vertical->dz_ifc[0] = vertical->z_lyr[0] - vertical->z[0];
    if (vertical->km > 1) {
        for (k=0; k<vertical->km-1; k++)
            vertical->dz_ifc[k] = vertical->z_lyr[k+1] - vertical->z_lyr[k];
    }
    vertical->dz_ifc[vertical->km] = vertical->z[vertical->km]
                                   - vertical->z_lyr[vertical->km-1];

    for (k=0; k<vertical->km; k++)
        vertical->dz_inv[k]     = 1.0 / vertical->dz[k];

    for (k=0; k<vertical->km+1; k++)
        vertical->dz_ifc_inv[k] = 1.0 / vertical->dz_ifc[k];
}

