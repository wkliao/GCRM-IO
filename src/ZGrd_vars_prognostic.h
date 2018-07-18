/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_vars_prognostic.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_ZGrd_vars_prognostic
#define H_ZGrd_vars_prognostic

#include <stdlib.h>

typedef struct {

   double /* DIMENSION(im,jm,     nprog,nsdm) */
      ****mss;       /* column integrated mass                          [???] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  PROGNOSTICS DEFINED AT CELL CENTERS (FACES) AND LAYER CENTERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,  km,nprog,nsdm) */
      *****rho,      /* density                                     [kg m^-3] */
      *****rho_adv,  /* density advection                      [kg m^-3 s^-1] */
      *****eta,      /* absolute vorticity                             [s^-1] */
      *****div,      /* divergence                                     [s^-1] */
      *****tht_lyr;  /* potential temperature defined in the LZ fashion   [K] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  PROGNOSTICS DEFINED AT CELL CENTERS (FACES) AND LAYER INTERFACES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,0:km,nprog,nsdm) */
      *****tht_ifc,  /* potential temperature defined in the CP fashion   [K] */
      *****w;        /* vertical velocity                            [m s^-1] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  PROGNOSTIC TRACERS DEFINED AT CELL CENTERS (FACES) AND LAYER CENTERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,  km,nprog,nsdm,wtrm) */
      ******wtr_lyr;  /* mixing ratio of water species   [kg kg^-1] */
   double /* DIMENSION(im,jm,  km,nprog,nsdm,trcm) */
      ******trc_lyr;  /* mixing ratio of passive tracers [kg kg^-1] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  TIME TENDENCIES DEFINED AT CELL CENTERS (FACES) AND LAYER CENTERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,  km,ntend,nsdm) */
      *****rho_f,    /* tendency of density                  [kg m^-3 s^-1] */
      *****rho_adv_f,/* tendency of density advection        [kg m^-3 s^-2] */
      *****eta_f,    /* tendency of absolute vorticity       [s^-2] */
      *****div_f,    /* tendency of divergence               [s^-2] */
      *****tht_lyr_f;/* tendency of LZ potential temperature [K s^-1] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  SUM OF TIME TENDENCIES DEFINED AT CELL CENTERS (FACES) AND LAYER CENTERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,  km,      nsdm) */
      ****div_f_sum1, /* sum of Adams-Bashforth tendencies of divergence [s^-2] */
      ****div_f_sum2; /* sum of Forward         tendencies of divergence [s^-2] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  TIME TENDENCIES DEFINED AT CELL CENTERS (FACES) AND LAYER INTERFACES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,0:km,ntend,nsdm) */
      *****tht_ifc_f, /* tendency of CP potential temperature [K s^-1] */
      *****w_f;       /* tendency of vertical velocity        [m s^-2] */
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  TIME TENDENCIES OF PROGNOSTIC TRACERS DEFINED AT CELL CENTERS (FACES) 
!  AND LAYER CENTERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   double /* DIMENSION(im,jm,  km,ntend,nsdm,wtrm) */
      ******wtr_lyr_f;/* mixing ratio of water species   [kg kg^-1 s^-1] */
   double /* DIMENSION(im,jm,  km,ntend,nsdm,trcm) */
      ******trc_lyr_f;/* mixing ratio of passive tracers [kg kg^-1 s^-1] */
} MODULE_ZGrd_vars_prognostic;


/* API declarations */
void init_MODULE_ZGrd_vars_prognostic(MODULE_ZGrd_vars_prognostic *vars_prog,
                                      int im, int jm, int km, int nprog,
                                      int nsdm, int ntend, int wtrm, int trcm);

void finalize_MODULE_ZGrd_vars_prognostic(MODULE_ZGrd_vars_prognostic *vars_prog);

void initialize_vars_prognostic(MODULE_ZGrd_vars_prognostic *vars_prog,
                                int im, int jm, int km, int nprog,
                                int nsdm, int ntend, int wtrm, int trcm);

void randomize_ZGrd_vars_prognostic(MODULE_ZGrd_vars_prognostic *vars_prog,
                                    int                          im,
                                    int                          jm,
                                    int                          km,
                                    int                          nprog,
                                    int                          nsdm,
                                    int                          ntend,
                                    int                          wtrm,
                                    int                          trcm);

#endif
