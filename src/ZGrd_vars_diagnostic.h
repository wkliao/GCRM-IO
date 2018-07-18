/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_vars_diagnostic.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_ZGrd_vars_diagnostic
#define H_ZGrd_vars_diagnostic

#include <stdlib.h>
#include "util.h"

typedef struct {
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  QUANTITIES DEFINED AT A SINGLE LEVEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,     nsdm) */
      ***f,            /* Coriolis force (=2.0*½*SIN (lat))       [s^-1] */
      ***geopot_surf;  /* geopotential of the Earth's surface [m^2 s^-2] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  relative vorticity, streamfunction and velocity potential and other stuff
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,  km,nsdm) */
      ****relative,
      ****psi,
      ****chi,
      ****ke,
      ****prs_lyr,
      ****exn_lyr,
      ****exn_prm,
      ****exner_lyr,
      ****geopot_lyr,
      ****tmp_ex,
      ****w_lyr,
      ****rho_t1,
      ****rho_t2,
      ****rho_adv_t1;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  variables defined at layer interfaces
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,0:km,nsdm) */
      ****prs,
      ****rho_ifc,
      ****exn,
      ****exner,
      ****geopot;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  sources/sinks, may be either interface or layer, if layer then level 0 not used
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,  km,nsdm) */
      ****Q_heating;   /* K/s */
   double /* DIMENSION(im,jm,  km,nsdm,wtrm) */
      *****wtr_lyr_src; /* kg/kg/s */
   double /* DIMENSION(im,jm,  km,nsdm,trcm) */
      *****trc_lyr_src; /* kg/kg/s */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  normal component of horizontal wind at edges and at layer centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(  edgm,im,jm,  km,nsdm) */
      *****wnd;  /* [m s^-1] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  wind at corners on the face grid and at layer centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(2,crnm,im,jm,  km,nsdm) */
      ******wnd_crn;  /* horizontal velocity at corners [m s^-1] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  wind at edges on the face grid and at layer centers, normal and tangential,
!     and geographic
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(2,edgm,im,jm,  km,nsdm) */
      ******wnd_edg,    /* horizontal velocity at edges normal/tangential  [m s^-1] */
      ******wnd_edg_uv; /* horizontal velocity at edges ew/ns              [m s^-1] */
   double /*  DIMENSION(im,jm,km,nsdm) */
      ****relative_edgg,
      ****divergence_edgg;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  wind at edges on the face grid and at layer interfaces
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(  edgm,im,jm,0:km,nsdm) */
      *****wnd_ifc,
      *****rhoflx_ifc;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  fluxes at edges on the face grid and at layer centers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(  edgm,im,jm,  km,nsdm) */
      *****rhoflx,
      *****flx;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  fluxes at edges on the face grid and at layer interface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(  edgm,im,jm,0:km,nsdm) */
      *****flx_ifc;
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  horizontal diffusivity and viscosity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(       im,jm,  km,nsdm) */
      ****rkh_h; /* diffusivity */
   double /* DIMENSION(       im,jm,  km,nsdm) */
      ****rkm_h; /* viscosity */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  vertical diffusivity and viscosity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(       im,jm,km-1,nsdm) */
      ****rkh_v, /* diffusivity at layer base */
      ****rkm_v; /* viscosity, at interfaces of layers */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  vertical diffusive fluxes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(       im,jm,  km,nsdm) */
      ****tht_vdif, /* theta flux */
      ****qwv_vdif; /* water vapor flux */
   double /* DIMENSION(       im,jm,0:km,nsdm) */
      ****fss_vdif, /* sensible heat flux diagnostic */
      ****fws_vdif; /* latent heat flux diagnostic */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  ventilation mass fluxes for surface fluxes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(  im,jm,nsdm) */
      ***ven2d; /* for theta and water vapor flux */
   double /* DIMENSION(crnm,im,jm,nsdm) */
      ****ven2d_crn; /* for momentum fluxes (at corners), use curl and div of this */
} MODULE_ZGrd_vars_diagnostic;

/* API declarations */
void init_MODULE_ZGrd_vars_diagnostic(MODULE_ZGrd_vars_diagnostic *vars_diag,
                                      int im, int jm, int km, int nsdm,
                                      int wtrm, int trcm, int edgm, int crnm);
void finalize_MODULE_ZGrd_vars_diagnostic(MODULE_ZGrd_vars_diagnostic *vars_diag);
void initialize_vars_diagnostic(MODULE_ZGrd_vars_diagnostic *vars_diag,
                                int im, int jm, int km, int nsdm,
                                int wtrm, int trcm, int edgm, int crnm);

void randomize_ZGrd_vars_diagnostic(MODULE_ZGrd_vars_diagnostic *vars_diag,
                                    int im,
                                    int jm,
                                    int km,
                                    int nsdm,
                                    int wtrm,
                                    int trcm,
                                    int edgm,
                                    int crnm);

#endif
