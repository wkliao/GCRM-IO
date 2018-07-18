/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_params_vertical.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_ZGRD_PARAMS_VERTICAL
#define H_ZGRD_PARAMS_VERTICAL

typedef struct {

    /* from MODULE ZGrd_params_vertical --------------------------------------*/
    int    km;    /* = 256  number of vertical layers in the  atmosphere */
    double z_top; /* = 18000.0 model top [m] */
    double prs_top; /* = 20.0 pressure at the model top [Pa] */

    double              z[257]; /* all [257] are from 0 to 256 in Fortran */
    double          z_lyr[256];
    double             dz[256];
    double         dz_inv[256];
    double         dz_ifc[257];
    double     dz_ifc_inv[257];

    double           zeta[257];
    double       zeta_lyr[256];
    double         d_zeta[256];
    double     d_zeta_ifc[257];
    double     d_zeta_inv[256];
    double d_zeta_inv_ifc[257];

} MODULE_ZGrd_params_vertical;

void initialize_params_vertical(MODULE_ZGrd_params_vertical *vertical);
void set_vertical_coordinate_levels(int km, double z_top, double *z, double *z_lyr);

#endif
