/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_vars_diagnostic.c 1858 2013-06-12 04:49:45Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>

#include "ZGrd_vars_diagnostic.h"
#include "util.h"

/*----< init_MODULE_ZGrd_vars_diagnostic() >----------------------------------*/
void init_MODULE_ZGrd_vars_diagnostic(MODULE_ZGrd_vars_diagnostic *vars_diag,
                                      int im,
                                      int jm,
                                      int km,
                                      int nsdm,
                                      int wtrm,
                                      int trcm,
                                      int edgm,
                                      int crnm)
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  QUANTITIES DEFINED AT A SINGLE LEVEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (im,jm,     nsdm) */
    vars_diag->f           = calloc_3D_dbl(nsdm, jm, im);
    vars_diag->geopot_surf = calloc_3D_dbl(nsdm, jm, im);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  relative vorticity, streamfunction and velocity potential and other stuff
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (im,jm,  km,nsdm) */
    vars_diag->relative   = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->psi        = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->chi        = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->ke         = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->prs_lyr    = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->exn_lyr    = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->exn_prm    = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->exner_lyr  = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->geopot_lyr = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->tmp_ex     = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->w_lyr      = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->rho_t1     = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->rho_t2     = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->rho_adv_t1 = calloc_4D_dbl(nsdm, km, jm, im);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  variables defined at layer interfaces
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (im,jm,0:km,nsdm) */
    vars_diag->prs     = calloc_4D_dbl(nsdm, km+1, jm, im);
    vars_diag->rho_ifc = calloc_4D_dbl(nsdm, km+1, jm, im);
    vars_diag->exn     = calloc_4D_dbl(nsdm, km+1, jm, im);
    vars_diag->exner   = calloc_4D_dbl(nsdm, km+1, jm, im);
    vars_diag->geopot  = calloc_4D_dbl(nsdm, km+1, jm, im);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  sources/sinks, may be either interface or layer, if layer then level 0
!  not used
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (im,jm,  km,nsdm) */
    vars_diag->Q_heating = calloc_4D_dbl(nsdm, km, jm, im);

    /* double (im,jm,  km,nsdm,wtrm) */
    vars_diag->wtr_lyr_src = calloc_5D_dbl(wtrm, nsdm, km, jm, im);

    /* double (im,jm,  km,nsdm,trcm) */
    vars_diag->trc_lyr_src = calloc_5D_dbl(trcm, nsdm, km, jm, im);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  normal component of horizontal wind at edges and at layer centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (  edgm,im,jm,  km,nsdm) */
    vars_diag->wnd = calloc_5D_dbl(nsdm, km, jm, im, edgm);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  wind at corners on the face grid and at layer centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (2,crnm,im,jm,  km,nsdm) */
    vars_diag->wnd_crn = calloc_6D_dbl(nsdm, km, jm, im, crnm, 2);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  wind at edges on the face grid and at layer centers, normal and tangential,
!     and geographic
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (2,edgm,im,jm,  km,nsdm) */
    vars_diag->wnd_edg    = calloc_6D_dbl(nsdm, km, jm, im, edgm, 2);
    vars_diag->wnd_edg_uv = calloc_6D_dbl(nsdm, km, jm, im, edgm, 2);

    /* double (im,jm,km,nsdm) */
    vars_diag->relative_edgg   = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->divergence_edgg = calloc_4D_dbl(nsdm, km, jm, im);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  wind at edges on the face grid and at layer interfaces
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (  edgm,im,jm,0:km,nsdm) */
    vars_diag->wnd_ifc    = calloc_5D_dbl(nsdm, km+1, jm, im, edgm);
    vars_diag->rhoflx_ifc = calloc_5D_dbl(nsdm, km+1, jm, im, edgm);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  fluxes at edges on the face grid and at layer centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (  edgm,im,jm,  km,nsdm) */
    vars_diag->rhoflx = calloc_5D_dbl(nsdm, km, jm, im, edgm);
    vars_diag->flx    = calloc_5D_dbl(nsdm, km, jm, im, edgm);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  fluxes at edges on the face grid and at layer interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (  edgm,im,jm,0:km,nsdm) */
    vars_diag->flx_ifc = calloc_5D_dbl(nsdm, km+1, jm, im, edgm);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  horizontal diffusivity and viscosity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (       im,jm,  km,nsdm) */
    vars_diag->rkh_h = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->rkm_h = calloc_4D_dbl(nsdm, km, jm, im);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  vertical diffusivity and viscosity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (       im,jm,km-1,nsdm) */
    vars_diag->rkh_v = calloc_4D_dbl(nsdm, km-1, jm, im);
    /* double (       im,jm,km-1,nsdm) */
    vars_diag->rkm_v = calloc_4D_dbl(nsdm, km-1, jm, im);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  vertical diffusive fluxes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (       im,jm,  km,nsdm) */
    vars_diag->tht_vdif = calloc_4D_dbl(nsdm, km, jm, im);
    vars_diag->qwv_vdif = calloc_4D_dbl(nsdm, km, jm, im);

    /* double (       im,jm,0:km,nsdm) */
    vars_diag->fss_vdif = calloc_4D_dbl(nsdm, km+1, jm, im);
    vars_diag->fws_vdif = calloc_4D_dbl(nsdm, km+1, jm, im);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  ventilation mass fluxes for surface fluxes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (  im,jm,nsdm) */
    vars_diag->ven2d = calloc_3D_dbl(nsdm, jm, im);

    /* double (crnm,im,jm,nsdm) */
    vars_diag->ven2d_crn = calloc_4D_dbl(nsdm, jm, im, crnm);
}
 
/*----< randomize_ZGrd_vars_diagnostic() >------------------------------------*/
void randomize_ZGrd_vars_diagnostic(MODULE_ZGrd_vars_diagnostic *vars_diag,
                                    int im,
                                    int jm,
                                    int km,
                                    int nsdm,
                                    int wtrm,
                                    int trcm,
                                    int edgm,
                                    int crnm)
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  QUANTITIES DEFINED AT A SINGLE LEVEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (im,jm,     nsdm) */
    random_3D_dbl(vars_diag->f,            nsdm, jm, im);
    random_3D_dbl(vars_diag->geopot_surf,  nsdm, jm, im);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  relative vorticity, streamfunction and velocity potential and other stuff
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (im,jm,  km,nsdm) */
    random_4D_dbl(vars_diag->relative,    nsdm, km, jm, im);
    random_4D_dbl(vars_diag->psi,         nsdm, km, jm, im);
    random_4D_dbl(vars_diag->chi,         nsdm, km, jm, im);
    random_4D_dbl(vars_diag->ke,          nsdm, km, jm, im);
    random_4D_dbl(vars_diag->prs_lyr,     nsdm, km, jm, im);
    random_4D_dbl(vars_diag->exn_lyr,     nsdm, km, jm, im);
    random_4D_dbl(vars_diag->exn_prm,     nsdm, km, jm, im);
    random_4D_dbl(vars_diag->exner_lyr,   nsdm, km, jm, im);
    random_4D_dbl(vars_diag->geopot_lyr,  nsdm, km, jm, im);
    random_4D_dbl(vars_diag->tmp_ex,      nsdm, km, jm, im);
    random_4D_dbl(vars_diag->w_lyr,       nsdm, km, jm, im);
    random_4D_dbl(vars_diag->rho_t1,      nsdm, km, jm, im);
    random_4D_dbl(vars_diag->rho_t2,      nsdm, km, jm, im);
    random_4D_dbl(vars_diag->rho_adv_t1,  nsdm, km, jm, im);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  variables defined at layer interfaces
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (im,jm,0:km,nsdm) */
    random_4D_dbl(vars_diag->prs,      nsdm, km+1, jm, im);
    random_4D_dbl(vars_diag->rho_ifc,  nsdm, km+1, jm, im);
    random_4D_dbl(vars_diag->exn,      nsdm, km+1, jm, im);
    random_4D_dbl(vars_diag->exner,    nsdm, km+1, jm, im);
    random_4D_dbl(vars_diag->geopot,   nsdm, km+1, jm, im);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  sources/sinks, may be either interface or layer, if layer then level 0
!  not used
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (im,jm,  km,nsdm) */
    random_4D_dbl(vars_diag->Q_heating,  nsdm, km, jm, im);

    /* double (im,jm,  km,nsdm,wtrm) */
    random_5D_dbl(vars_diag->wtr_lyr_src,  wtrm, nsdm, km, jm, im);

    /* double (im,jm,  km,nsdm,trcm) */
    random_5D_dbl(vars_diag->trc_lyr_src,  trcm, nsdm, km, jm, im);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  normal component of horizontal wind at edges and at layer centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (  edgm,im,jm,  km,nsdm) */
    random_5D_dbl(vars_diag->wnd,  nsdm, km, jm, im, edgm);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  wind at corners on the face grid and at layer centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (2,crnm,im,jm,  km,nsdm) */
    random_6D_dbl(vars_diag->wnd_crn,  nsdm, km, jm, im, crnm, 2);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  wind at edges on the face grid and at layer centers, normal and tangential,
!     and geographic
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (2,edgm,im,jm,  km,nsdm) */
    random_6D_dbl(vars_diag->wnd_edg,     nsdm, km, jm, im, edgm, 2);
    random_6D_dbl(vars_diag->wnd_edg_uv,  nsdm, km, jm, im, edgm, 2);

    /* double (im,jm,km,nsdm) */
    random_4D_dbl(vars_diag->relative_edgg,    nsdm, km, jm, im);
    random_4D_dbl(vars_diag->divergence_edgg,  nsdm, km, jm, im);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  wind at edges on the face grid and at layer interfaces
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (  edgm,im,jm,0:km,nsdm) */
    random_5D_dbl(vars_diag->wnd_ifc,     nsdm, km+1, jm, im, edgm);
    random_5D_dbl(vars_diag->rhoflx_ifc,  nsdm, km+1, jm, im, edgm);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  fluxes at edges on the face grid and at layer centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (  edgm,im,jm,  km,nsdm) */
    random_5D_dbl(vars_diag->rhoflx,  nsdm, km, jm, im, edgm);
    random_5D_dbl(vars_diag->flx,     nsdm, km, jm, im, edgm);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  fluxes at edges on the face grid and at layer interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (  edgm,im,jm,0:km,nsdm) */
    random_5D_dbl(vars_diag->flx_ifc,  nsdm, km+1, jm, im, edgm);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  horizontal diffusivity and viscosity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (       im,jm,  km,nsdm) */
    random_4D_dbl(vars_diag->rkh_h,  nsdm, km, jm, im);
    random_4D_dbl(vars_diag->rkm_h,  nsdm, km, jm, im);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  vertical diffusivity and viscosity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (       im,jm,km-1,nsdm) */
    random_4D_dbl(vars_diag->rkh_v,  nsdm, km-1, jm, im);
    /* double (       im,jm,km-1,nsdm) */
    random_4D_dbl(vars_diag->rkm_v,  nsdm, km-1, jm, im);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  vertical diffusive fluxes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (       im,jm,  km,nsdm) */
    random_4D_dbl(vars_diag->tht_vdif,  nsdm, km, jm, im);
    random_4D_dbl(vars_diag->qwv_vdif,  nsdm, km, jm, im);

    /* double (       im,jm,0:km,nsdm) */
    random_4D_dbl(vars_diag->fss_vdif,  nsdm, km+1, jm, im);
    random_4D_dbl(vars_diag->fws_vdif,  nsdm, km+1, jm, im);
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  ventilation mass fluxes for surface fluxes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* double (  im,jm,nsdm) */
    random_3D_dbl(vars_diag->ven2d,  nsdm, jm, im);

    /* double (crnm,im,jm,nsdm) */
    random_4D_dbl(vars_diag->ven2d_crn,  nsdm, jm, im, crnm);
}

/*----< initialize_vars_diagnostic() >----------------------------------------*/
void initialize_vars_diagnostic(MODULE_ZGrd_vars_diagnostic *vars_diag,
                                int im,
                                int jm,
                                int km,
                                int nsdm,
                                int wtrm,
                                int trcm,
                                int edgm,
                                int crnm)
{
    /* Originally, this function initializes all vriables to zeros */
    return;
}

/*----< finalize_MODULE_ZGrd_vars_diagnostic() >------------------------------*/
void finalize_MODULE_ZGrd_vars_diagnostic(MODULE_ZGrd_vars_diagnostic *vars_diag)
{
    free_3D_dbl(vars_diag->f);
    free_3D_dbl(vars_diag->geopot_surf);
    free_4D_dbl(vars_diag->relative);
    free_4D_dbl(vars_diag->psi);
    free_4D_dbl(vars_diag->chi);
    free_4D_dbl(vars_diag->ke);
    free_4D_dbl(vars_diag->prs_lyr);
    free_4D_dbl(vars_diag->exn_lyr);
    free_4D_dbl(vars_diag->exn_prm);
    free_4D_dbl(vars_diag->exner_lyr);
    free_4D_dbl(vars_diag->geopot_lyr);
    free_4D_dbl(vars_diag->tmp_ex);
    free_4D_dbl(vars_diag->w_lyr);
    free_4D_dbl(vars_diag->rho_t1);
    free_4D_dbl(vars_diag->rho_t2);
    free_4D_dbl(vars_diag->rho_adv_t1);
    free_4D_dbl(vars_diag->prs);
    free_4D_dbl(vars_diag->rho_ifc);
    free_4D_dbl(vars_diag->exn);
    free_4D_dbl(vars_diag->exner);
    free_4D_dbl(vars_diag->geopot);
    free_4D_dbl(vars_diag->Q_heating);
    free_5D_dbl(vars_diag->wtr_lyr_src);
    free_5D_dbl(vars_diag->trc_lyr_src);
    free_5D_dbl(vars_diag->wnd);
    free_6D_dbl(vars_diag->wnd_crn);
    free_6D_dbl(vars_diag->wnd_edg);
    free_6D_dbl(vars_diag->wnd_edg_uv);
    free_4D_dbl(vars_diag->relative_edgg);
    free_4D_dbl(vars_diag->divergence_edgg);
    free_5D_dbl(vars_diag->wnd_ifc);
    free_5D_dbl(vars_diag->rhoflx_ifc);
    free_5D_dbl(vars_diag->rhoflx);
    free_5D_dbl(vars_diag->flx);
    free_5D_dbl(vars_diag->flx_ifc);
    free_4D_dbl(vars_diag->rkh_h);
    free_4D_dbl(vars_diag->rkm_h);
    free_4D_dbl(vars_diag->rkh_v);
    free_4D_dbl(vars_diag->rkm_v);
    free_4D_dbl(vars_diag->tht_vdif);
    free_4D_dbl(vars_diag->qwv_vdif);
    free_4D_dbl(vars_diag->fss_vdif);
    free_4D_dbl(vars_diag->fws_vdif);
    free_3D_dbl(vars_diag->ven2d);
    free_4D_dbl(vars_diag->ven2d_crn);
}

