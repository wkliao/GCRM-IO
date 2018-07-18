/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_vars_prognostic.c 1858 2013-06-12 04:49:45Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>

#include "ZGrd_vars_prognostic.h"
#include "util.h"

/*----< init_MODULE_ZGrd_vars_prognostic() >----------------------------------*/
void init_MODULE_ZGrd_vars_prognostic(MODULE_ZGrd_vars_prognostic *vars_prog,
                                      int                          im,
                                      int                          jm,
                                      int                          km,
                                      int                          nprog,
                                      int                          nsdm,
                                      int                          ntend,
                                      int                          wtrm,
                                      int                          trcm)
{
    /* double  DIMENSION(im,jm,     nprog,nsdm) */
    vars_prog->mss = calloc_4D_dbl(nsdm, nprog, jm, im);

    /* double  DIMENSION(im,jm,  km,nprog,nsdm) */
    vars_prog->rho     = calloc_5D_dbl(nsdm, nprog, km, jm, im);
    vars_prog->rho_adv = calloc_5D_dbl(nsdm, nprog, km, jm, im);
    vars_prog->eta     = calloc_5D_dbl(nsdm, nprog, km, jm, im);
    vars_prog->div     = calloc_5D_dbl(nsdm, nprog, km, jm, im);
    vars_prog->tht_lyr = calloc_5D_dbl(nsdm, nprog, km, jm, im);

    /* double  DIMENSION(im,jm,0:km,nprog,nsdm) */
    vars_prog->tht_ifc = calloc_5D_dbl(nsdm, nprog, km+1, jm, im);
    vars_prog->w       = calloc_5D_dbl(nsdm, nprog, km+1, jm, im);

    /* double  DIMENSION(im,jm,  km,nprog,nsdm,wtrm) */
    vars_prog->wtr_lyr = calloc_6D_dbl(wtrm, nsdm, nprog, km, jm, im);
    /* double  DIMENSION(im,jm,  km,nprog,nsdm,trcm) */
    vars_prog->trc_lyr = calloc_6D_dbl(trcm, nsdm, nprog, km, jm, im);

    /* double  DIMENSION(im,jm,  km,ntend,nsdm) */
    vars_prog->rho_f     = calloc_5D_dbl(nsdm, ntend, km, jm, im);
    vars_prog->rho_adv_f = calloc_5D_dbl(nsdm, ntend, km, jm, im);
    vars_prog->eta_f     = calloc_5D_dbl(nsdm, ntend, km, jm, im);
    vars_prog->div_f     = calloc_5D_dbl(nsdm, ntend, km, jm, im);
    vars_prog->tht_lyr_f = calloc_5D_dbl(nsdm, ntend, km, jm, im);

    /* double  DIMENSION(im,jm,  km,      nsdm) */
    vars_prog->div_f_sum1 = calloc_4D_dbl(nsdm, km, jm, im);
    vars_prog->div_f_sum2 = calloc_4D_dbl(nsdm, km, jm, im);

    /* double  DIMENSION(im,jm,0:km,ntend,nsdm) */
    vars_prog->tht_ifc_f = calloc_5D_dbl(nsdm, ntend, km+1, jm, im);
    vars_prog->w_f       = calloc_5D_dbl(nsdm, ntend, km+1, jm, im);

    /* double  DIMENSION(im,jm,  km,ntend,nsdm,wtrm) */
    vars_prog->wtr_lyr_f = calloc_6D_dbl(wtrm, nsdm, ntend, km, jm, im);
    /* double  DIMENSION(im,jm,  km,ntend,nsdm,trcm) */
    vars_prog->trc_lyr_f = calloc_6D_dbl(trcm, nsdm, ntend, km, jm, im);
}

/*----< randomize_ZGrd_vars_prognostic() >------------------------------------*/
void randomize_ZGrd_vars_prognostic(MODULE_ZGrd_vars_prognostic *vars_prog,
                                    int                          im,
                                    int                          jm,
                                    int                          km,
                                    int                          nprog,
                                    int                          nsdm,
                                    int                          ntend,
                                    int                          wtrm,
                                    int                          trcm)
{
    /* double  DIMENSION(im,jm,     nprog,nsdm) */
    random_4D_dbl(vars_prog->mss, nsdm, nprog, jm, im);

    /* double  DIMENSION(im,jm,  km,nprog,nsdm) */
    random_5D_dbl(vars_prog->rho,      nsdm, nprog, km, jm, im);
    random_5D_dbl(vars_prog->rho_adv,  nsdm, nprog, km, jm, im);
    random_5D_dbl(vars_prog->eta,      nsdm, nprog, km, jm, im);
    random_5D_dbl(vars_prog->div,      nsdm, nprog, km, jm, im);
    random_5D_dbl(vars_prog->tht_lyr,  nsdm, nprog, km, jm, im);

    /* double  DIMENSION(im,jm,0:km,nprog,nsdm) */
    random_5D_dbl(vars_prog->tht_ifc,  nsdm, nprog, km+1, jm, im);
    random_5D_dbl(vars_prog->w,        nsdm, nprog, km+1, jm, im);

    /* double  DIMENSION(im,jm,  km,nprog,nsdm,wtrm) */
    random_6D_dbl(vars_prog->wtr_lyr,  wtrm, nsdm, nprog, km, jm, im);
    /* double  DIMENSION(im,jm,  km,nprog,nsdm,trcm) */
    random_6D_dbl(vars_prog->trc_lyr,  trcm, nsdm, nprog, km, jm, im);

    /* double  DIMENSION(im,jm,  km,ntend,nsdm) */
    random_5D_dbl(vars_prog->rho_f,      nsdm, ntend, km, jm, im);
    random_5D_dbl(vars_prog->rho_adv_f,  nsdm, ntend, km, jm, im);
    random_5D_dbl(vars_prog->eta_f,      nsdm, ntend, km, jm, im);
    random_5D_dbl(vars_prog->div_f,      nsdm, ntend, km, jm, im);
    random_5D_dbl(vars_prog->tht_lyr_f,  nsdm, ntend, km, jm, im);

    /* double  DIMENSION(im,jm,  km,      nsdm) */
    random_4D_dbl(vars_prog->div_f_sum1,  nsdm, km, jm, im);
    random_4D_dbl(vars_prog->div_f_sum2,  nsdm, km, jm, im);

    /* double  DIMENSION(im,jm,0:km,ntend,nsdm) */
    random_5D_dbl(vars_prog->tht_ifc_f,  nsdm, ntend, km+1, jm, im);
    random_5D_dbl(vars_prog->w_f,        nsdm, ntend, km+1, jm, im);

    /* double  DIMENSION(im,jm,  km,ntend,nsdm,wtrm) */
    random_6D_dbl(vars_prog->wtr_lyr_f,  wtrm, nsdm, ntend, km, jm, im);
    /* double  DIMENSION(im,jm,  km,ntend,nsdm,trcm) */
    random_6D_dbl(vars_prog->trc_lyr_f,  trcm, nsdm, ntend, km, jm, im);
}

/*----< initialize_vars_prognostic() >----------------------------------------*/
void initialize_vars_prognostic(MODULE_ZGrd_vars_prognostic *vars_prog,
                                int                          im,
                                int                          jm,
                                int                          km,
                                int                          nprog,
                                int                          nsdm,
                                int                          ntend,
                                int                          wtrm,
                                int                          trcm)
{
    /* Originally, this function initializes all vriables to zeros */
    return;
}

/*----< finalize_MODULE_ZGrd_vars_prognostic() >------------------------------*/
void finalize_MODULE_ZGrd_vars_prognostic(MODULE_ZGrd_vars_prognostic *vars_prog)
{
    free_4D_dbl(vars_prog->mss);
    free_5D_dbl(vars_prog->rho);
    free_5D_dbl(vars_prog->rho_adv);
    free_5D_dbl(vars_prog->eta);
    free_5D_dbl(vars_prog->div);
    free_5D_dbl(vars_prog->tht_lyr);
    free_5D_dbl(vars_prog->tht_ifc);
    free_5D_dbl(vars_prog->w);
    free_6D_dbl(vars_prog->wtr_lyr);
    free_6D_dbl(vars_prog->trc_lyr);
    free_5D_dbl(vars_prog->rho_f);
    free_5D_dbl(vars_prog->rho_adv_f);
    free_5D_dbl(vars_prog->eta_f);
    free_5D_dbl(vars_prog->div_f);
    free_5D_dbl(vars_prog->tht_lyr_f);
    free_4D_dbl(vars_prog->div_f_sum1);
    free_4D_dbl(vars_prog->div_f_sum2);
    free_5D_dbl(vars_prog->tht_ifc_f);
    free_5D_dbl(vars_prog->w_f);
    free_6D_dbl(vars_prog->wtr_lyr_f);
    free_6D_dbl(vars_prog->trc_lyr_f);
}
