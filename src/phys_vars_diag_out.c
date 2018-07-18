/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: phys_vars_diag_out.c 1858 2013-06-12 04:49:45Z wkliao $
 */

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "phys_vars_diag_out.h"
#include "util.h"


/*----< init_MODULE_phys_vars_diag_out() >------------------------------------*/
void init_MODULE_phys_vars_diag_out(MODULE_phys_vars_diag_out *vars,
                                    int                        nsdm,
                                    int                        km_phys,
                                    int                        jm,
                                    int                        im)
{
    vars->km_phys = km_phys;

    vars->l_accum = 0;  /* if true, physics diagnostic output is accumulated,
                           else only instantaneous value saved at time of
                                output */

    /* microphysics diagnostics */
    vars->pr           = calloc_3D_dbl(nsdm, jm, im);
    vars->prfz         = calloc_3D_dbl(nsdm, jm, im);
    vars->q_latent     = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->src_local_qv = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->src_local_qc = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->src_local_qi = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->src_local_qs = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->src_local_qr = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->src_local_qg = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->src_fall_qs  = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->src_fall_qr  = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->src_fall_qg  = calloc_4D_dbl(nsdm, km_phys, jm, im);
   
    /* radiation diagnostics */
    vars->q_lwrad      = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->q_swrad      = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->flut         = calloc_3D_dbl(nsdm, jm, im);
    vars->flus         = calloc_3D_dbl(nsdm, jm, im);
    vars->flds         = calloc_3D_dbl(nsdm, jm, im);
    vars->fsdt         = calloc_3D_dbl(nsdm, jm, im);
    vars->fsut         = calloc_3D_dbl(nsdm, jm, im);
    vars->fsds         = calloc_3D_dbl(nsdm, jm, im);
    vars->fsus         = calloc_3D_dbl(nsdm, jm, im);
    vars->q_lwrad_cs   = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->q_swrad_cs   = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->flut_cs      = calloc_3D_dbl(nsdm, jm, im);
    vars->flus_cs      = calloc_3D_dbl(nsdm, jm, im);
    vars->flds_cs      = calloc_3D_dbl(nsdm, jm, im);
    vars->fsut_cs      = calloc_3D_dbl(nsdm, jm, im);
    vars->fsds_cs      = calloc_3D_dbl(nsdm, jm, im);
    vars->fsus_cs      = calloc_3D_dbl(nsdm, jm, im);
   
    /* turbulence diagnostics */
    vars->evap         = calloc_3D_dbl(nsdm, jm, im);
    vars->fss          = calloc_3D_dbl(nsdm, jm, im);
    vars->q_tht_sgs    = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->qwv_sgs_tend = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->qcw_sgs_tend = calloc_4D_dbl(nsdm, km_phys, jm, im);
    vars->qci_sgs_tend = calloc_4D_dbl(nsdm, km_phys, jm, im);

    /* below 3 are actually defined in MODULE phys_vars_diagnostic */
    vars->surfaceT     = calloc_3D_dbl(nsdm, jm, im);
    vars->zrough       = calloc_3D_dbl(nsdm, jm, im);
    vars->gwet         = calloc_3D_dbl(nsdm, jm, im);
}

/*----< finalize_MODULE_phys_vars_diag_out() >--------------------------------*/
void finalize_MODULE_phys_vars_diag_out(MODULE_phys_vars_diag_out *vars)
{
    /* microphysics diagnostics */
    free_3D_dbl(vars->pr);
    free_3D_dbl(vars->prfz);
    free_4D_dbl(vars->q_latent);
    free_4D_dbl(vars->src_local_qv);
    free_4D_dbl(vars->src_local_qc);
    free_4D_dbl(vars->src_local_qi);
    free_4D_dbl(vars->src_local_qs);
    free_4D_dbl(vars->src_local_qr);
    free_4D_dbl(vars->src_local_qg);
    free_4D_dbl(vars->src_fall_qs);
    free_4D_dbl(vars->src_fall_qr);
    free_4D_dbl(vars->src_fall_qg);
   
    /* radiation diagnostics */
    free_4D_dbl(vars->q_lwrad);
    free_4D_dbl(vars->q_swrad);
    free_3D_dbl(vars->flut);
    free_3D_dbl(vars->flus);
    free_3D_dbl(vars->flds);
    free_3D_dbl(vars->fsdt);
    free_3D_dbl(vars->fsut);
    free_3D_dbl(vars->fsds);
    free_3D_dbl(vars->fsus);
    free_4D_dbl(vars->q_lwrad_cs);
    free_4D_dbl(vars->q_swrad_cs);
    free_3D_dbl(vars->flut_cs);
    free_3D_dbl(vars->flus_cs);
    free_3D_dbl(vars->flds_cs);
    free_3D_dbl(vars->fsut_cs);
    free_3D_dbl(vars->fsds_cs);
    free_3D_dbl(vars->fsus_cs);
   
    /* turbulence diagnostics */
    free_3D_dbl(vars->evap);
    free_3D_dbl(vars->fss);
    free_4D_dbl(vars->q_tht_sgs);
    free_4D_dbl(vars->qwv_sgs_tend);
    free_4D_dbl(vars->qcw_sgs_tend);
    free_4D_dbl(vars->qci_sgs_tend);

    /* below 3 are actually defined in MODULE phys_vars_diagnostic */
    free_3D_dbl(vars->surfaceT);
    free_3D_dbl(vars->zrough);
    free_3D_dbl(vars->gwet);
}

/*----< randomize_phys_vars_diag_out() >--------------------------------------*/
void randomize_phys_vars_diag_out(MODULE_phys_vars_diag_out *vars,
                                  int                        nsdm,
                                  int                        km_phys,
                                  int                        jm,
                                  int                        im)
{
    /* microphysics diagnostics */
    random_3D_dbl(vars->pr,           nsdm, jm, im);
    random_3D_dbl(vars->prfz,         nsdm, jm, im);
    random_4D_dbl(vars->q_latent,     nsdm, km_phys, jm, im);
    random_4D_dbl(vars->src_local_qv, nsdm, km_phys, jm, im);
    random_4D_dbl(vars->src_local_qc, nsdm, km_phys, jm, im);
    random_4D_dbl(vars->src_local_qi, nsdm, km_phys, jm, im);
    random_4D_dbl(vars->src_local_qs, nsdm, km_phys, jm, im);
    random_4D_dbl(vars->src_local_qr, nsdm, km_phys, jm, im);
    random_4D_dbl(vars->src_local_qg, nsdm, km_phys, jm, im);
    random_4D_dbl(vars->src_fall_qs,  nsdm, km_phys, jm, im);
    random_4D_dbl(vars->src_fall_qr,  nsdm, km_phys, jm, im);
    random_4D_dbl(vars->src_fall_qg,  nsdm, km_phys, jm, im);
   
    /* radiation diagnostics */
    random_4D_dbl(vars->q_lwrad,      nsdm, km_phys, jm, im);
    random_4D_dbl(vars->q_swrad,      nsdm, km_phys, jm, im);
    random_3D_dbl(vars->flut,         nsdm, jm, im);
    random_3D_dbl(vars->flus,         nsdm, jm, im);
    random_3D_dbl(vars->flds,         nsdm, jm, im);
    random_3D_dbl(vars->fsdt,         nsdm, jm, im);
    random_3D_dbl(vars->fsut,         nsdm, jm, im);
    random_3D_dbl(vars->fsds,         nsdm, jm, im);
    random_3D_dbl(vars->fsus,         nsdm, jm, im);
    random_4D_dbl(vars->q_lwrad_cs,   nsdm, km_phys, jm, im);
    random_4D_dbl(vars->q_swrad_cs,   nsdm, km_phys, jm, im);
    random_3D_dbl(vars->flut_cs,      nsdm, jm, im);
    random_3D_dbl(vars->flus_cs,      nsdm, jm, im);
    random_3D_dbl(vars->flds_cs,      nsdm, jm, im);
    random_3D_dbl(vars->fsut_cs,      nsdm, jm, im);
    random_3D_dbl(vars->fsds_cs,      nsdm, jm, im);
    random_3D_dbl(vars->fsus_cs,      nsdm, jm, im);
   
    /* turbulence diagnostics */
    random_3D_dbl(vars->evap,         nsdm, jm, im);
    random_3D_dbl(vars->fss,          nsdm, jm, im);
    random_4D_dbl(vars->q_tht_sgs,    nsdm, km_phys, jm, im);
    random_4D_dbl(vars->qwv_sgs_tend, nsdm, km_phys, jm, im);
    random_4D_dbl(vars->qcw_sgs_tend, nsdm, km_phys, jm, im);
    random_4D_dbl(vars->qci_sgs_tend, nsdm, km_phys, jm, im);

    /* below 3 are actually defined in MODULE phys_vars_diagnostic */
    random_3D_dbl(vars->surfaceT,     nsdm, jm, im);
    random_3D_dbl(vars->zrough,       nsdm, jm, im);
    random_3D_dbl(vars->gwet,         nsdm, jm, im);
}

/*----< init_phys_vars_diag_out() >-------------------------------------------*/
void init_phys_vars_diag_out(MODULE_phys_vars_diag_out *vars,
                             int                        nsdm,
                             int                        km_phys,
                             int                        jm,
                             int                        im)
{
    /* this subroutine resets the accumulation counters and diagnostic
       output holders */
    int nsdm_jm_im         = nsdm * jm * im * sizeof(double);
    int nsdm_km_phys_jm_im = km_phys * nsdm_jm_im;

    /* microphysics diagnostics */
    memset(vars->pr,           0, nsdm_jm_im);
    memset(vars->prfz,         0, nsdm_jm_im);
    memset(vars->q_latent,     0, nsdm_km_phys_jm_im);
    memset(vars->src_local_qv, 0, nsdm_km_phys_jm_im);
    memset(vars->src_local_qc, 0, nsdm_km_phys_jm_im);
    memset(vars->src_local_qi, 0, nsdm_km_phys_jm_im);
    memset(vars->src_local_qs, 0, nsdm_km_phys_jm_im);
    memset(vars->src_local_qr, 0, nsdm_km_phys_jm_im);
    memset(vars->src_local_qg, 0, nsdm_km_phys_jm_im);
    memset(vars->src_fall_qs,  0, nsdm_km_phys_jm_im);
    memset(vars->src_fall_qr,  0, nsdm_km_phys_jm_im);
    memset(vars->src_fall_qg,  0, nsdm_km_phys_jm_im);
   
    /* radiation diagnostics */
    memset(vars->q_lwrad,      0, nsdm_km_phys_jm_im);
    memset(vars->q_swrad,      0, nsdm_km_phys_jm_im);
    memset(vars->flut,         0, nsdm_jm_im);
    memset(vars->flus,         0, nsdm_jm_im);
    memset(vars->flds,         0, nsdm_jm_im);
    memset(vars->fsdt,         0, nsdm_jm_im);
    memset(vars->fsut,         0, nsdm_jm_im);
    memset(vars->fsds,         0, nsdm_jm_im);
    memset(vars->fsus,         0, nsdm_jm_im);
    memset(vars->q_lwrad_cs,   0, nsdm_km_phys_jm_im);
    memset(vars->q_swrad_cs,   0, nsdm_km_phys_jm_im);
    memset(vars->flut_cs,      0, nsdm_jm_im);
    memset(vars->flus_cs,      0, nsdm_jm_im);
    memset(vars->flds_cs,      0, nsdm_jm_im);
    memset(vars->fsut_cs,      0, nsdm_jm_im);
    memset(vars->fsds_cs,      0, nsdm_jm_im);
    memset(vars->fsus_cs,      0, nsdm_jm_im);
   
    /* turbulence diagnostics */
    memset(vars->evap,         0, nsdm_jm_im);
    memset(vars->fss,          0, nsdm_jm_im);
    memset(vars->q_tht_sgs,    0, nsdm_km_phys_jm_im);
    memset(vars->qwv_sgs_tend, 0, nsdm_km_phys_jm_im);
    memset(vars->qcw_sgs_tend, 0, nsdm_km_phys_jm_im);
    memset(vars->qci_sgs_tend, 0, nsdm_km_phys_jm_im);
}

/*----< norm_phys_vars_diag_out() >-------------------------------------------*/
void norm_phys_vars_diag_out(MODULE_phys_vars_diag_out *vars,
                             int                        nsdm,
                             int                        km_phys,
                             int                        jm,
                             int                        im)
{
    /* this subroutine normalizes the diagnostic output variables */

#define FAC3D(buf) for(n=0;n<nsdm;n++) for(j=0;j<jm;j++) for(i=0;i<im;i++) buf[n][j][i] *= fac;
#define FAC4D(buf) for(n=0;n<nsdm;n++) for(k=0;k<km_phys;k++) for(j=0;j<jm;j++) for(i=0;i<im;i++) buf[n][k][j][i] *= fac;
    int i, j, k, n;
    double fac;
    /* microphysics diagnostics */
    if (vars->count_micro_var > 1) {
        fac = 1.0 / vars->count_micro_var;
        FAC3D(vars->pr)
        FAC3D(vars->prfz)
        FAC4D(vars->q_latent)
        FAC4D(vars->src_local_qv)
        FAC4D(vars->src_local_qc)
        FAC4D(vars->src_local_qi)
        FAC4D(vars->src_local_qs)
        FAC4D(vars->src_local_qr)
        FAC4D(vars->src_local_qg)
        FAC4D(vars->src_fall_qs)
        FAC4D(vars->src_fall_qr)
        FAC4D(vars->src_fall_qg)
    }
   
    /* radiation diagnostics */
    if (vars->count_rad_var > 1) {
        fac = 1.0 / vars->count_rad_var;
        FAC4D(vars->q_lwrad)
        FAC4D(vars->q_swrad)
        FAC3D(vars->flut)
        FAC3D(vars->flus)
        FAC3D(vars->flds)
        FAC3D(vars->fsdt)
        FAC3D(vars->fsut)
        FAC3D(vars->fsds)
        FAC3D(vars->fsus)
        FAC4D(vars->q_lwrad_cs)
        FAC4D(vars->q_swrad_cs)
        FAC3D(vars->flut_cs)
        FAC3D(vars->flus_cs)
        FAC3D(vars->flds_cs)
        FAC3D(vars->fsut_cs)
        FAC3D(vars->fsds_cs)
        FAC3D(vars->fsus_cs)
    }

    /* turbulence diagnostics */
    if (vars->count_turb_var > 1) {
        fac = 1.0 / vars->count_turb_var;
        FAC3D(vars->evap)
        FAC3D(vars->fss)
        FAC4D(vars->q_tht_sgs)
        FAC4D(vars->qwv_sgs_tend)
        FAC4D(vars->qcw_sgs_tend)
        FAC4D(vars->qci_sgs_tend)
    }
}

