/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: ZGrd_register.c 1858 2013-06-12 04:49:45Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "gio.h"
#include "grid_params.h"
#include "grid_metrics.h"
#include "ZGrd_params_vertical.h"
#include "ZGrd_vars_diagnostic.h"
#include "ZGrd_vars_prognostic.h"
#include "ZGrd_clusdet.h"
#include "ZGrd_params_time.h"

#include "util.h"


/*----< initialize_gio_api() >------------------------------------------------*/
/* from ZGrd_register.F90
 * This subroutine is designed to attach data fields to
 * data descriptors inside the GCRM IO API */
void initialize_gio_api(gio_parameters              *gio_param,
                        MODULE_ZGrd_params_tracer   *tracer,
                        MODULE_grid_params          *grid_params,
                        MODULE_grid_metrics         *grid_metrics,
                        MODULE_ZGrd_params_vertical *vertical,
                        MODULE_ZGrd_vars_diagnostic *vars_diag,
                        MODULE_ZGrd_vars_prognostic *vars_prog,
                        MODULE_ZGrd_clusdet         *clusdet,
                        MODULE_ZGrd_params_time     *time,
                        MODULE_phys_vars_diag_out   *phys_vars,
                        int                          enable_physics)
{
    int n, panel_id, nsdm, im, jm, km, xm[3], glo[3], ghi[3];

    nsdm = grid_params->nsdm;
    im   = grid_params->im;
    jm   = grid_params->jm;
    km   = vertical->km;

    int iqwv = tracer->iqwv,
        iqcw = tracer->iqcw,
        iqrw = tracer->iqrw,
        iqci = tracer->iqci,
        iqsn = tracer->iqsn,
        iqgr = tracer->iqgr;

    /* the 4 arrays below are defined in grid/grid_metrics.F90
       They are used to created geodescic grid */
    int    ****tag_glbl_nghbr = grid_metrics->tag_glbl_nghbr;
    double          ****point = grid_metrics->point;
    double        *****corner = grid_metrics->corner;
    double            ***area = grid_metrics->area;

/*==============================================================================
 *  Register grid data
 *=============================================================================*/
    /* nsdm is the number of blocks owned by this process */
    for (n=0; n<nsdm; n++) {

        get_big_block_index(grid_params->level_max,
                            grid_params->sbdmn_iota,
                            grid_params->sbdmn[grid_params->level_max].lst[n],
                            xm);
        /* xm[] returned is 1-based */

        panel_id = xm[2] - 1;
        glo[0]   = xm[0];
        glo[1]   = xm[1];
        ghi[0]   = glo[0] + im - 3;
        ghi[1]   = glo[1] + jm - 3;

        gio_grid_setup(gio_param,
                       tag_glbl_nghbr[n][1]+1,  /* [][6] */
                                point[n][1]+1,  /* [][3] */
                               corner[n][1]+1,  /* [][6][3] */
// wkliao: bug here?  12/23/2012
                                 // area[n][1]+1,  /* [gio_istride*gio_jstride] */
                                 area[n][0],  /* [gio_istride*gio_jstride] */
                       panel_id, glo, ghi);
    }

    int nsd_north = grid_params->nsd_north;  /* north pole block ID */
    if (grid_params->l_agent_north)
        gio_grid_setup_pole(gio_param,
                            tag_glbl_nghbr[nsd_north][jm-1][1],
                                     point[nsd_north][jm-1][1],
                                    corner[nsd_north][jm-1][1],
                                      area[nsd_north][jm-1][1],
                            'n');

    int nsd_south = grid_params->nsd_south;  /* south pole block ID */
    if (grid_params->l_agent_south)
        gio_grid_setup_pole(gio_param,
                            tag_glbl_nghbr[nsd_south][1][im-1],
                                     point[nsd_south][1][im-1],
                                    corner[nsd_south][1][im-1],
                                      area[nsd_south][1][im-1],
                            's');

//==============================================================================
//  Register remaining data
//==============================================================================

//==============================================================================
//  Set exner/relative vorticity
//==============================================================================
#if NOT_YET
!   do k = 1, km
!      tmpry01(:,:,k,:) = exner_lyr(:,:,k,:) + exner_lyr_ref(k)
!      relative(:,:,k,:) = eta(:,:,k,np1,:) - f(:,:,:)
!   end do
!  IF (l_theta_CP) THEN
!    do k = 0, km
!       tmp_ex(:,:,k,:) = exner(:,:,k,:) + exner_ifc_ref(k)
!    end do
!  ELSE
    /* vars_diag->tmp_ex = calloc_4D_dbl(nsdm, km, jm, im); */
    for (k=0; k<km; k++) {
        int j, i;
        for (n=0; n<nsdm; n++)
            for (j=0; j<jm; j++)
                for (i=0; i<im; i++)
                    vars_diag->tmp_ex[n][k][j][i] =
                    vars_diag->exner_lyr[n][k][j][i] +
                    ZGrd_vars_reference.exner_lyr_ref[k];
                    /* in ZGrd, initialize_REF() is never called? */
    }
#endif
//==============================================================================
//  Register remaining data
//==============================================================================

    for (n=0; n<nsdm; n++) {

        get_big_block_index(grid_params->level_max,
                            grid_params->sbdmn_iota,
                            grid_params->sbdmn[grid_params->level_max].lst[n],
                            xm);

        panel_id = xm[2] - 1;
        glo[0]   = xm[0];
        glo[1]   = xm[1];
        ghi[0]   = glo[0] + im - 3;
        ghi[1]   = glo[1] + jm - 3;

        /* Cell data */
#define REG_DIAG(name, buf) gio_register_dfield(name, &vars_diag->buf, glo, ghi, panel_id, gio_param);
#define REG_PROG(name, buf) gio_register_dfield(name, &vars_prog->buf, glo, ghi, panel_id, gio_param);

        REG_DIAG("strm_func",                   psi[n][0][1][1])
        REG_DIAG("vel_pot",                     chi[n][0][1][1])
        REG_DIAG("dbl_psi",                     psi[n][0][1][1])
        REG_DIAG("dbl_chi",                     chi[n][0][1][1])
        REG_PROG("vorticity",                   eta[n][0][0][1][1])
        REG_DIAG("ke",                           ke[n][0][1][1])
        REG_PROG("divergence",                  div[n][0][0][1][1])
        REG_DIAG("divergence_edge", divergence_edgg[n][0][1][1])
        REG_PROG("dbl_eta",                     eta[n][0][0][1][1])
        REG_PROG("dbl_eta_f",                 eta_f[n][0][0][1][1])
        REG_PROG("mass",                        rho[n][0][0][1][1])
        REG_DIAG("dbl_exn_lyr",             exn_lyr[0][0][1][1])
        REG_DIAG("dbl_exn_prm",             exn_prm[0][0][1][1])
        REG_DIAG("dbl_prs_lyr",             prs_lyr[0][0][1][1])
        REG_PROG("dbl_rho",                     rho[n][0][0][1][1])
        REG_PROG("dbl_rho_f",                 rho_f[n][0][0][1][1])
        REG_PROG("dbl_rho_adv",             rho_adv[n][0][0][1][1])
        REG_PROG("dbl_rho_adv_f",         rho_adv_f[n][0][0][1][1])
        REG_PROG("dbl_div",                     div[n][0][0][1][1])
        REG_PROG("dbl_div_f",                 div_f[n][0][0][1][1])
        REG_DIAG("relative",               relative[n][0][1][1])
        REG_DIAG("relative_edge",     relative_edgg[n][0][1][1])
        REG_DIAG("dbl_relative",           relative[n][0][1][1])
        REG_PROG("dbl_mss",                     mss[n][0][1][1])

        // 2D data
        gio_register_dfield("clus_mask", &clusdet->clus_mask[n][1][1],
                            glo, ghi, panel_id, gio_param);

        //  Data at interfaces
        REG_DIAG("pressure",                    prs[n][0][1][1])
        REG_DIAG("dbl_prs_ifc",                 prs[n][0][1][1])
        REG_DIAG("dbl_exn_ifc",                 exn[n][0][1][1])
        REG_PROG("dbl_w",                         w[n][0][0][1][1])
        REG_PROG("dbl_w_f",                     w_f[n][0][0][1][1])
        REG_PROG("w_vert",                        w[n][0][0][1][1])
#ifdef CURRENTLY_DISABLED
        if (l_theta_CP) {
            REG_PROG("temperature_ifc",     tht_ifc[n][0][0][1][1])
            REG_PROG("dbl_tht_ifc",         tht_ifc[n][0][0][1][1])
            REG_PROG("dbl_tht_ifc_f",     tht_ifc_f[n][0][0][1][1])
            REG_DIAG("heat_flux_vdiff_lyr",fss_vdif[n][0][1][1])
            /* ADDED EXNER */
            REG_DIAG("exner_ifc",            tmp_ex[n][0][1][1])
        }
        REG_PROG("tracer1",              trc_ifc[0][n][0][0][1][1])
        REG_PROG("tracer2",              trc_ifc[1][n][0][0][1][1])
        REG_PROG("tracer3",              trc_ifc[2][n][0][0][1][1])
        REG_PROG("tracer4",              trc_ifc[3][n][0][0][1][1])
#endif
        REG_DIAG("geopotential",             geopot[n][0][1][1])
        REG_DIAG("heat_flux_vdiff_ifc",    fss_vdif[n][0][1][1])
// IF (l_theta_LZ) THEN
        REG_PROG("temperature_lyr",         tht_lyr[n][0][0][1][1])
        REG_PROG("dbl_tht_lyr",             tht_lyr[n][0][0][1][1])
        REG_PROG("dbl_tht_lyr_f",         tht_lyr_f[n][0][0][1][1])
//ADDED EXNER
        REG_DIAG("exner_lyr",                tmp_ex[n][0][1][1])
// ENDIF

        if (enable_physics) {
#define REG_PHY(name, buf) gio_register_dfield(name, &phys_vars->buf, glo, ghi, panel_id, gio_param);
        // ADDED PHYSICS
#ifdef OBSOLETE
        if (l_theta_CP) {
            REG_PHY("wtr_flux_vdiff_lyr",     fws_vdif[n][0][1][1])
            // wtr(im,jm,0:km,nprog,nsdm,nwtr)
            REG_PHY("water_vapor_ifc",       wtr[iqwv][n][0][0][1][1])
            REG_PHY("cloud_water_ifc",       wtr[iqcw][n][0][0][1][1])
            REG_PHY("rain_mmr_ifc",          wtr[iqrw][n][0][0][1][1])
            REG_PHY("cloud_ice_ifc",         wtr[iqci][n][0][0][1][1])
            REG_PHY("snow_mmr_ifc",          wtr[iqsn][n][0][0][1][1])
            REG_PHY("graupel_mmr_ifc",       wtr[iqgr][n][0][0][1][1])
            REG_PHY("heating_sw_ifc",          q_swrad[n][0][1][1])
            REG_PHY("heating_lw_ifc",          q_lwrad[n][0][1][1])
            REG_PHY("heating_sw_cs_ifc",    q_swrad_cs[n][0][1][1])
            REG_PHY("heating_lw_cs_ifc",    q_lwrad_cs[n][0][1][1])
            REG_PHY("heating_latent_ifc",     q_latent[n][0][1][1])
            REG_PHY("qwv_tend_micro_ifc", src_local_qv[n][0][1][1])
            REG_PHY("qcw_tend_micro_ifc", src_local_qc[n][0][1][1])
            REG_PHY("qci_tend_micro_ifc", src_local_qi[n][0][1][1])
            REG_PHY("qrw_tend_micro_ifc", src_local_qr[n][0][1][1])
            REG_PHY("qsn_tend_micro_ifc", src_local_qs[n][0][1][1])
            REG_PHY("qgr_tend_micro_ifc", src_local_qg[n][0][1][1])
            REG_PHY("dbl_vapor_ifc",         wtr[iqwv][n][0][0][1][1])
            REG_PHY("dbl_cloudwater_ifc",    wtr[iqcw][n][0][0][1][1])
            REG_PHY("dbl_rainmmr_ifc",       wtr[iqrw][n][0][0][1][1])
            REG_PHY("dbl_cloudice_ifc",      wtr[iqci][n][0][0][1][1])
            REG_PHY("dbl_snowmmr_ifc_f",     wtr[iqsn][n][0][0][1][1])
            REG_PHY("dbl_graupelmmr_ifc_f",  wtr[iqgr][n][0][0][1][1])
        }
#endif
        // if (l_theta_LZ) {
            REG_DIAG("wtr_flux_vdiff_ifc",          fws_vdif[n][0][1][1])
            REG_PROG("water_vapor_lyr",        wtr_lyr[iqwv][n][0][0][1][1])
            REG_PROG("cloud_water_lyr",        wtr_lyr[iqcw][n][0][0][1][1])
            REG_PROG("rain_mmr_lyr",           wtr_lyr[iqrw][n][0][0][1][1])
            REG_PROG("cloud_ice_lyr",          wtr_lyr[iqci][n][0][0][1][1])
            REG_PROG("snow_mmr_lyr",           wtr_lyr[iqsn][n][0][0][1][1])
            REG_PROG("graupel_mmr_lyr",        wtr_lyr[iqgr][n][0][0][1][1])
            REG_PHY("heating_sw_lyr",                q_swrad[n][0][1][1])
            REG_PHY("heating_lw_lyr",                q_lwrad[n][0][1][1])
            REG_PHY("heating_sw_cs_lyr",          q_swrad_cs[n][0][1][1])
            REG_PHY("heating_lw_cs_lyr",          q_lwrad_cs[n][0][1][1])
            REG_PHY("heating_latent_lyr",           q_latent[n][0][1][1])
            REG_PHY("qwv_tend_micro_lyr",       src_local_qv[n][0][1][1])
            REG_PHY("qcw_tend_micro_lyr",       src_local_qc[n][0][1][1])
            REG_PHY("qci_tend_micro_lyr",       src_local_qi[n][0][1][1])
            REG_PHY("qrw_tend_micro_lyr",       src_local_qr[n][0][1][1])
            REG_PHY("qsn_tend_micro_lyr",       src_local_qs[n][0][1][1])
            REG_PHY("qgr_tend_micro_lyr",       src_local_qg[n][0][1][1])
            REG_PROG("dbl_vapor_lyr",        wtr_lyr  [iqwv][n][0][0][1][1])
            REG_PROG("dbl_cloudwater_lyr",   wtr_lyr  [iqcw][n][0][0][1][1])
            REG_PROG("dbl_rainmmr_lyr",      wtr_lyr  [iqrw][n][0][0][1][1])
            REG_PROG("dbl_cloudice_lyr",     wtr_lyr  [iqci][n][0][0][1][1])
            REG_PROG("dbl_snowmmr_lyr",      wtr_lyr  [iqsn][n][0][0][1][1])
            REG_PROG("dbl_graupelmmr_lyr",   wtr_lyr  [iqgr][n][0][0][1][1])
            REG_PROG("dbl_vapor_lyr_f",      wtr_lyr_f[iqwv][n][0][0][1][1])
            REG_PROG("dbl_cloudwater_lyr_f", wtr_lyr_f[iqcw][n][0][0][1][1])
            REG_PROG("dbl_rainmmr_lyr_f",    wtr_lyr_f[iqrw][n][0][0][1][1])
            REG_PROG("dbl_cloudice_lyr_f",   wtr_lyr_f[iqci][n][0][0][1][1])
            REG_PROG("dbl_snowmmr_lyr_f",    wtr_lyr_f[iqsn][n][0][0][1][1])
            REG_PROG("dbl_graupelmmr_lyr_f", wtr_lyr_f[iqgr][n][0][0][1][1])
        // }
        //  one layer fields (sfc or TOA)
        REG_PHY("prec_tot",         pr[n][1][1])
        REG_PHY("prec_frz",       prfz[n][1][1])
        REG_PHY("olr",            flut[n][1][1])
        REG_PHY("swinc",          fsdt[n][1][1])
        // Register the "f" field, which is 2d, static
        REG_PHY("dbl_ts",     surfaceT[n][1][1])
        REG_PHY("dbl_gwet",       gwet[n][1][1])
        REG_PHY("dbl_zrough",   zrough[n][1][1])
        } /* of if (enable_physics) */

// Corner data
        REG_DIAG("u",                       wnd_crn[n][0][1][1][0][0])
        REG_DIAG("v",                       wnd_crn[n][0][1][1][0][0])
//      REG_PROG("u_src",              wind_src_crn[n][0][1][1][0][0])
//      REG_PROG("v_src",              wind_src_crn[n][0][1][1][0][0])
// Edge data
        REG_DIAG("wind",                        wnd[n][0][1][1][0])
//      REG_DIAG("flux",                        flx[n][0][1][1][0])
//      REG_PROG("mass_flux",                mssflx[n][0][1][1][0])
//      REG_DIAG("flux_ifc",                flx_ifc[n][0][1][1][0])
//      REG_PROG("mass_flux_ifc",        mssflx_ifc[n][0][1][1][0])
    } /* end of loop n */

// NOTE: The length values remain the same as above since this
//       data is heavily strided

#define REG_DIAG_N(name, buf) gio_register_dpole(name, &vars_diag->buf, 'n',&gio_param->data);
#define REG_PROG_N(name, buf) gio_register_dpole(name, &vars_prog->buf, 'n',&gio_param->data);
    if (grid_params->l_agent_north) {
        int nsd = grid_params->nsd_north;  /* north pole block ID */
        jm--;  /* jm is not array index in C */

        REG_DIAG_N("strm_func",                   psi[nsd][0][jm][1])
        REG_DIAG_N("vel_pot",                     chi[nsd][0][jm][1])
        REG_DIAG_N("dbl_psi",                     psi[nsd][0][jm][1])
        REG_DIAG_N("dbl_chi",                     chi[nsd][0][jm][1])
        REG_PROG_N("divergence",                  div[nsd][0][0][jm][1])
        REG_DIAG_N("divergence_edge", divergence_edgg[nsd][0][jm][1])
        REG_PROG_N("vorticity",                   eta[nsd][0][0][jm][1])
        REG_DIAG_N("ke",                           ke[nsd][0][jm][1])
        REG_PROG_N("dbl_eta",                     eta[nsd][0][0][jm][1])
        REG_PROG_N("dbl_eta_f",                 eta_f[nsd][0][0][jm][1])
        REG_PROG_N("mass",                        rho[nsd][0][0][jm][1])
        REG_DIAG_N("dbl_exn_lyr",             exn_lyr[nsd][0][jm][1])
        REG_DIAG_N("dbl_exn_prm",             exn_prm[nsd][0][jm][1])
        REG_DIAG_N("dbl_prs_lyr",             prs_lyr[nsd][0][jm][1])
        REG_PROG_N("dbl_rho",                     rho[nsd][0][0][jm][1])
        REG_PROG_N("dbl_rho_f",                 rho_f[nsd][0][0][jm][1])
        REG_PROG_N("dbl_rho_adv",             rho_adv[nsd][0][0][jm][1])
        REG_PROG_N("dbl_rho_adv_f",         rho_adv_f[nsd][0][0][jm][1])
        REG_PROG_N("dbl_div",                     div[nsd][0][0][jm][1])
        REG_PROG_N("dbl_div_f",                 div_f[nsd][0][0][jm][1])
        REG_PROG_N("dbl_mss",                     mss[nsd][0][jm][1])

        gio_register_dpole("clus_mask", &clusdet->clus_mask[nsd][jm][1],
                           'n', &gio_param->data);

        REG_DIAG_N("pressure",                    prs[nsd][0][jm][1])
        REG_DIAG_N("dbl_prs_ifc",                 prs[nsd][0][jm][1])
        REG_DIAG_N("dbl_exn_ifc",                 exn[nsd][0][jm][1])
        REG_PROG_N("dbl_w",                         w[nsd][0][0][jm][1])
        REG_PROG_N("dbl_w_f",                     w_f[nsd][0][0][jm][1])
        REG_PROG_N("w_vert",                        w[nsd][0][0][jm][1])
#ifdef CURRENTLY_DISABLED
        if (l_theta_CP) {
            REG_PROG_N("temperature_ifc",     tht_ifc[nsd][0][0][jm][1])
            REG_PROG_N("dbl_tht_ifc",         tht_ifc[nsd][0][0][jm][1])
            REG_PROG_N("dbl_tht_ifc_f",     tht_ifc_f[nsd][0][0][jm][1])
            REG_DIAG_N("heat_flux_vdiff_lyr",fss_vdif[nsd][0][jm][1])
            /* ADDED EXNER */
            REG_DIAG_N("exner_ifc",            tmp_ex[nsd][0][jm][1])
        }
#endif
// IF (l_theta_LZ) THEN
        REG_PROG_N("temperature_lyr",         tht_lyr[nsd][0][0][jm][1])
        REG_PROG_N("dbl_tht_lyr",             tht_lyr[nsd][0][0][jm][1])
        REG_PROG_N("dbl_tht_lyr_f",         tht_lyr_f[nsd][0][0][jm][1])
//ADDED EXNER
        REG_DIAG_N("exner_lyr",                tmp_ex[nsd][0][jm][1])
// ENDIF
        REG_DIAG_N("relative",               relative[nsd][0][jm][1])
        REG_DIAG_N("relative_edge",     relative_edgg[nsd][0][jm][1])
        REG_DIAG_N("dbl_relative",           relative[nsd][0][jm][1])
//      REG_DIAG_N("tracer1",              trc_ifc[0][nsd][0][0][jm][1])
//      REG_DIAG_N("tracer2",              trc_ifc[1][nsd][0][0][jm][1])
//      REG_DIAG_N("tracer3",              trc_ifc[2][nsd][0][0][jm][1])
//      REG_DIAG_N("tracer4",              trc_ifc[3][nsd][0][0][jm][1])

        REG_DIAG_N("heat_flux_vdiff_ifc",    fss_vdif[nsd][0][jm][1])
        REG_DIAG_N("geopotential",             geopot[nsd][0][jm][1])
//      REG_DIAG_N("q_heating",             Q_heating[nsd][0][jm][1])
//      REG_DIAG_N("pot_tmp_adv",            tht_horz[nsd][0][jm][1])

        if (enable_physics) {
#define REG_PHY_N(name, buf) gio_register_dpole(name, &phys_vars->buf, 'n',&gio_param->data);
        // ADDED PHYSICS
#ifdef OBSOLETE
        if (l_theta_CP) {
            REG_DIAG_N("wtr_flux_vdiff_lyr",     fws_vdif[nsd][0][jm][1])
            REG_PHY_N("water_vapor_ifc",        wtr[iqwv][nsd][0][0][jm][1])
            REG_PHY_N("cloud_water_ifc",        wtr[iqcw][nsd][0][0][jm][1])
            REG_PHY_N("rain_mmr_ifc",           wtr[iqrw][nsd][0][0][jm][1])
            REG_PHY_N("cloud_ice_ifc",          wtr[iqci][nsd][0][0][jm][1])
            REG_PHY_N("snow_mmr_ifc",           wtr[iqsn][nsd][0][0][jm][1])
            REG_PHY_N("graupel_mmr_ifc",        wtr[iqgr][nsd][0][0][jm][1])
            REG_PHY_N("heating_sw_ifc",           q_swrad[nsd][0][jm][1])
            REG_PHY_N("heating_lw_ifc",           q_lwrad[nsd][0][jm][1])
            REG_PHY_N("heating_sw_cs_ifc",     q_swrad_cs[nsd][0][jm][1])
            REG_PHY_N("heating_lw_cs_ifc",     q_lwrad_cs[nsd][0][jm][1])
            REG_PHY_N("heating_latent_ifc",      q_latent[nsd][0][jm][1])
            REG_PHY_N("qwv_tend_micro_ifc",  src_local_qv[nsd][0][jm][1])
            REG_PHY_N("qcw_tend_micro_ifc",  src_local_qc[nsd][0][jm][1])
            REG_PHY_N("qci_tend_micro_ifc",  src_local_qi[nsd][0][jm][1])
            REG_PHY_N("qrw_tend_micro_ifc",  src_local_qr[nsd][0][jm][1])
            REG_PHY_N("qsn_tend_micro_ifc",  src_local_qs[nsd][0][jm][1])
            REG_PHY_N("qgr_tend_micro_ifc",  src_local_qg[nsd][0][jm][1])
            REG_PHY_N("dbl_vapor_ifc",          wtr[iqwv][nsd][0][0][jm][1])
            REG_PHY_N("dbl_cloudwater_ifc",     wtr[iqcw][nsd][0][0][jm][1])
            REG_PHY_N("dbl_rainmmr_ifc",        wtr[iqrw][nsd][0][0][jm][1])
            REG_PHY_N("dbl_cloudice_ifc",       wtr[iqci][nsd][0][0][jm][1])
            REG_PHY_N("dbl_snowmmr_ifc",        wtr[iqsn][nsd][0][0][jm][1])
            REG_PHY_N("dbl_graupelmmr_ifc",     wtr[iqgr][nsd][0][0][jm][1])
            REG_PHY_N("dbl_vapor_ifc_f",      wtr_f[iqwv][nsd][0][0][jm][1])
            REG_PHY_N("dbl_cloudwater_ifc_f", wtr_f[iqcw][nsd][0][0][jm][1])
            REG_PHY_N("dbl_rainmmr_ifc_f",    wtr_f[iqrw][nsd][0][0][jm][1])
            REG_PHY_N("dbl_cloudice_ifc_f",   wtr_f[iqci][nsd][0][0][jm][1])
            REG_PHY_N("dbl_snowmmr_ifc_f",    wtr_f[iqsn][nsd][0][0][jm][1])
            REG_PHY_N("dbl_graupelmmr_ifc_f", wtr_f[iqgr][nsd][0][0][jm][1])
        }
#endif
        // if (l_theta_LZ) {
            REG_DIAG_N("wtr_flux_vdiff_ifc",         fws_vdif[nsd][0][jm][1])
            REG_PROG_N("water_vapor_lyr",       wtr_lyr[iqwv][nsd][0][0][jm][1])
            REG_PROG_N("cloud_water_lyr",       wtr_lyr[iqcw][nsd][0][0][jm][1])
            REG_PROG_N("rain_mmr_lyr",          wtr_lyr[iqrw][nsd][0][0][jm][1])
            REG_PROG_N("cloud_ice_lyr",         wtr_lyr[iqci][nsd][0][0][jm][1])
            REG_PROG_N("snow_mmr_lyr",          wtr_lyr[iqsn][nsd][0][0][jm][1])
            REG_PROG_N("graupel_mmr_lyr",       wtr_lyr[iqgr][nsd][0][0][jm][1])
            REG_PHY_N("heating_sw_lyr",               q_swrad[nsd][0][jm][1])
            REG_PHY_N("heating_lw_lyr",               q_lwrad[nsd][0][jm][1])
            REG_PHY_N("heating_sw_cs_lyr",         q_swrad_cs[nsd][0][jm][1])
            REG_PHY_N("heating_lw_cs_lyr",         q_lwrad_cs[nsd][0][jm][1])
            REG_PHY_N("heating_latent_lyr",          q_latent[nsd][0][jm][1])
            REG_PHY_N("qwv_tend_micro_lyr",      src_local_qv[nsd][0][jm][1])
            REG_PHY_N("qcw_tend_micro_lyr",      src_local_qc[nsd][0][jm][1])
            REG_PHY_N("qci_tend_micro_lyr",      src_local_qi[nsd][0][jm][1])
            REG_PHY_N("qrw_tend_micro_lyr",      src_local_qr[nsd][0][jm][1])
            REG_PHY_N("qsn_tend_micro_lyr",      src_local_qs[nsd][0][jm][1])
            REG_PHY_N("qgr_tend_micro_lyr",      src_local_qg[nsd][0][jm][1])
            REG_PROG_N("dbl_vapor_lyr",         wtr_lyr[iqwv][nsd][0][0][jm][1])
            REG_PROG_N("dbl_cloudwater_lyr",    wtr_lyr[iqcw][nsd][0][0][jm][1])
            REG_PROG_N("dbl_rainmmr_lyr",       wtr_lyr[iqrw][nsd][0][0][jm][1])
            REG_PROG_N("dbl_cloudice_lyr",      wtr_lyr[iqci][nsd][0][0][jm][1])
            REG_PROG_N("dbl_snowmmr_lyr",       wtr_lyr[iqsn][nsd][0][0][jm][1])
            REG_PROG_N("dbl_graupelmmr_lyr",    wtr_lyr[iqgr][nsd][0][0][jm][1])
            REG_PROG_N("dbl_vapor_lyr_f",     wtr_lyr_f[iqwv][nsd][0][0][jm][1])
            REG_PROG_N("dbl_cloudwater_lyr_f",wtr_lyr_f[iqcw][nsd][0][0][jm][1])
            REG_PROG_N("dbl_rainmmr_lyr_f",   wtr_lyr_f[iqrw][nsd][0][0][jm][1])
            REG_PROG_N("dbl_cloudice_lyr_f",  wtr_lyr_f[iqci][nsd][0][0][jm][1])
            REG_PROG_N("dbl_snowmmr_lyr_f",   wtr_lyr_f[iqsn][nsd][0][0][jm][1])
            REG_PROG_N("dbl_graupelmmr_lyr_f",wtr_lyr_f[iqgr][nsd][0][0][jm][1])
        // }

        // one layer fields (sfc or TOA)
        REG_PHY_N("prec_tot",       pr[nsd][jm][1])
        REG_PHY_N("prec_frz",     prfz[nsd][jm][1])
        REG_PHY_N("olr",          flut[nsd][jm][1])
        REG_PHY_N("swinc",        fsdt[nsd][jm][1])
        REG_PHY_N("dbl_ts",   surfaceT[nsd][jm][1])
        REG_PHY_N("dbl_gwet",     gwet[nsd][jm][1])
        REG_PHY_N("dbl_zrough", zrough[nsd][jm][1])
        } /* of if (enable_physics) */
    }

#define REG_DIAG_S(name, buf) gio_register_dpole(name, &vars_diag->buf, 's',&gio_param->data);
#define REG_PROG_S(name, buf) gio_register_dpole(name, &vars_prog->buf, 's',&gio_param->data);
    if (grid_params->l_agent_south) {
        int nsd = grid_params->nsd_south;  /* south pole block ID */
        im--;  /* im is not array index in C */

        REG_DIAG_S("strm_func",                    psi[nsd][0][1][im])
        REG_DIAG_S("vel_pot",                      chi[nsd][0][1][im])
        REG_DIAG_S("dbl_psi",                      psi[nsd][0][1][im])
        REG_DIAG_S("dbl_chi",                      chi[nsd][0][1][im])
        REG_PROG_S("divergence",                   div[nsd][0][0][1][im])
        REG_DIAG_S("divergence_edge",  divergence_edgg[nsd][0][1][im])
        REG_PROG_S("vorticity",                    eta[nsd][0][0][1][im])
        REG_DIAG_S("ke",                            ke[nsd][0][1][im])
        REG_PROG_S("dbl_eta",                      eta[nsd][0][0][1][im])
        REG_PROG_S("dbl_eta_f",                  eta_f[nsd][0][0][1][im])
        REG_PROG_S("mass",                         rho[nsd][0][0][1][im])
        REG_PROG_S("dbl_rho",                      rho[nsd][0][0][1][im])
        REG_PROG_S("dbl_rho_f",                  rho_f[nsd][0][0][1][im])
        REG_PROG_S("dbl_rho_adv",              rho_adv[nsd][0][0][1][im])
        REG_PROG_S("dbl_rho_adv_f",          rho_adv_f[nsd][0][0][1][im])
        REG_DIAG_S("dbl_exn_lyr",              exn_lyr[nsd][0][1][im])
        REG_DIAG_S("dbl_exn_prm",              exn_prm[nsd][0][1][im])
        REG_DIAG_S("dbl_prs_lyr",              prs_lyr[nsd][0][1][im])
        REG_PROG_S("dbl_div",                      div[nsd][0][0][1][im])
        REG_PROG_S("dbl_div_f",                  div_f[nsd][0][0][1][im])
        REG_PROG_S("dbl_mss",                      mss[nsd][0][1][im])

        gio_register_dpole("clus_mask", &clusdet->clus_mask[nsd][1][im],
                           's', &gio_param->data);

        REG_DIAG_S("pressure",                     prs[nsd][0][1][im])
        REG_DIAG_S("dbl_prs_ifc",                  prs[nsd][0][1][im])
        REG_DIAG_S("dbl_exn_ifc",                  exn[nsd][0][1][im])
        REG_PROG_S("dbl_w",                          w[nsd][0][0][1][im])
        REG_PROG_S("dbl_w_f",                      w_f[nsd][0][0][1][im])
        REG_PROG_S("w_vert",                         w[nsd][0][0][1][im])
#ifdef CURRENTLY_DISABLED
        if (l_theta_CP) {
            REG_PROG_S("temperature_ifc",      tht_ifc[nsd][0][0][1][im])
            REG_PROG_S("dbl_tht_ifc",          tht_ifc[nsd][0][0][1][im])
            REG_PROG_S("dbl_tht_ifc_f",      tht_ifc_f[nsd][0][0][1][im])
            REG_DIAG_S("heat_flux_vdiff_lyr", fss_vdif[nsd][0][1][im])
            /* ADDED EXNER */
            REG_DIAG_S("exner_ifc",             tmp_ex[nsd][0][1][im])
        }
#endif
// IF (l_theta_LZ) THEN
        REG_PROG_S("temperature_lyr",          tht_lyr[nsd][0][0][1][im])
        REG_PROG_S("dbl_tht_lyr",              tht_lyr[nsd][0][0][1][im])
        REG_PROG_S("dbl_tht_lyr_f",          tht_lyr_f[nsd][0][0][1][im])
        REG_DIAG_S("heat_flux_vdiff_ifc",     fss_vdif[nsd][0][1][im])
//ADDED EXNER
        REG_DIAG_S("exner_lyr",                 tmp_ex[nsd][0][1][im])
// ENDIF
        REG_DIAG_S("relative",                relative[nsd][0][1][im])
        REG_DIAG_S("relative_edge",      relative_edgg[nsd][0][1][im])
        REG_DIAG_S("dbl_relative",            relative[nsd][0][1][im])
//      REG_DIAG_S("tracer1",               trc_ifc[0][nsd][0][0][1][im])
//      REG_DIAG_S("tracer2",               trc_ifc[1][nsd][0][0][1][im])
//      REG_DIAG_S("tracer3",               trc_ifc[2][nsd][0][0][1][im])
//      REG_DIAG_S("tracer4",               trc_ifc[3][nsd][0][0][1][im])
        REG_DIAG_S("geopotential",              geopot[nsd][0][1][im])
//      REG_DIAG_S("q_heating",              Q_heating[nsd][0][1][im])
//      REG_DIAG_S("pot_tmp_adv",             tht_horz[nsd][0][1][im])

        if (enable_physics) {
#define REG_PHY_S(name, buf) gio_register_dpole(name, &phys_vars->buf, 's',&gio_param->data);
#ifdef OBSOLETE
        // ADDED PHYSICS
        if (l_theta_CP) {
            REG_DIAG_S("wtr_flux_vdiff_ifc",     fws_vdif[nsd][0][1][im])
            REG_PHY_S("water_vapor_ifc",        wtr[iqwv][nsd][0][0][1][im])
            REG_PHY_S("cloud_water_ifc",        wtr[iqcw][nsd][0][0][1][im])
            REG_PHY_S("rain_mmr_ifc",           wtr[iqrw][nsd][0][0][1][im])
            REG_PHY_S("cloud_ice_ifc",          wtr[iqci][nsd][0][0][1][im])
            REG_PHY_S("snow_mmr_ifc",           wtr[iqsn][nsd][0][0][1][im])
            REG_PHY_S("graupel_mmr_ifc",        wtr[iqgr][nsd][0][0][1][im])
            REG_PHY_S("heating_sw_ifc",           q_swrad[nsd][0][1][im])
            REG_PHY_S("heating_lw_ifc",           q_lwrad[nsd][0][1][im])
            REG_PHY_S("heating_sw_cs_ifc",     q_swrad_cs[nsd][0][1][im])
            REG_PHY_S("heating_lw_cs_ifc",     q_lwrad_cs[nsd][0][1][im])
            REG_PHY_S("heating_latent_ifc",      q_latent[nsd][0][1][im])
            REG_PHY_S("qwv_tend_micro_ifc",  src_local_qv[nsd][0][1][im])
            REG_PHY_S("qcw_tend_micro_ifc",  src_local_qc[nsd][0][1][im])
            REG_PHY_S("qci_tend_micro_ifc",  src_local_qi[nsd][0][1][im])
            REG_PHY_S("qrw_tend_micro_ifc",  src_local_qr[nsd][0][1][im])
            REG_PHY_S("qsn_tend_micro_ifc",  src_local_qs[nsd][0][1][im])
            REG_PHY_S("qgr_tend_micro_ifc",  src_local_qg[nsd][0][1][im])
            REG_PHY_S("dbl_vapor_ifc",          wtr[iqwv][nsd][0][0][1][im])
            REG_PHY_S("dbl_cloudwater_ifc",     wtr[iqcw][nsd][0][0][1][im])
            REG_PHY_S("dbl_rainmmr_ifc",        wtr[iqrw][nsd][0][0][1][im])
            REG_PHY_S("dbl_cloudice_ifc",       wtr[iqci][nsd][0][0][1][im])
            REG_PHY_S("dbl_snowmmr_ifc",        wtr[iqsn][nsd][0][0][1][im])
            REG_PHY_S("dbl_graupelmmr_ifc",     wtr[iqgr][nsd][0][0][1][im])
            REG_PHY_S("dbl_vapor_ifc_f",      wtr_f[iqwv][nsd][0][0][1][im])
            REG_PHY_S("dbl_cloudwater_ifc_f", wtr_f[iqcw][nsd][0][0][1][im])
            REG_PHY_S("dbl_rainmmr_ifc_f",    wtr_f[iqrw][nsd][0][0][1][im])
            REG_PHY_S("dbl_cloudice_ifc_f",   wtr_f[iqci][nsd][0][0][1][im])
            REG_PHY_S("dbl_snowmmr_ifc_f",    wtr_f[iqsn][nsd][0][0][1][im])
            REG_PHY_S("dbl_graupelmmr_ifc_f", wtr_f[iqgr][nsd][0][0][1][im])
        }
#endif
        //  if (l_theta_LZ) {
            REG_DIAG_S("wtr_flux_vdiff_ifc",         fws_vdif[nsd][0][1][im])
            REG_PROG_S("water_vapor_lyr",       wtr_lyr[iqwv][nsd][0][0][1][im])
            REG_PROG_S("cloud_water_lyr",       wtr_lyr[iqcw][nsd][0][0][1][im])
            REG_PROG_S("rain_mmr_lyr",          wtr_lyr[iqrw][nsd][0][0][1][im])
            REG_PROG_S("cloud_ice_lyr",         wtr_lyr[iqci][nsd][0][0][1][im])
            REG_PROG_S("snow_mmr_lyr",          wtr_lyr[iqsn][nsd][0][0][1][im])
            REG_PROG_S("graupel_mmr_lyr",       wtr_lyr[iqgr][nsd][0][0][1][im])
            REG_PHY_S("heating_sw_lyr",               q_swrad[nsd][0][1][im])
            REG_PHY_S("heating_lw_lyr",               q_lwrad[nsd][0][1][im])
            REG_PHY_S("heating_sw_cs_lyr",         q_swrad_cs[nsd][0][1][im])
            REG_PHY_S("heating_lw_cs_lyr",         q_lwrad_cs[nsd][0][1][im])
            REG_PHY_S("heating_latent_lyr",          q_latent[nsd][0][1][im])
            REG_PHY_S("qwv_tend_micro_lyr",      src_local_qv[nsd][0][1][im])
            REG_PHY_S("qcw_tend_micro_lyr",      src_local_qc[nsd][0][1][im])
            REG_PHY_S("qci_tend_micro_lyr",      src_local_qi[nsd][0][1][im])
            REG_PHY_S("qrw_tend_micro_lyr",      src_local_qr[nsd][0][1][im])
            REG_PHY_S("qsn_tend_micro_lyr",      src_local_qs[nsd][0][1][im])
            REG_PHY_S("qgr_tend_micro_lyr",      src_local_qg[nsd][0][1][im])
            REG_PROG_S("dbl_vapor_lyr",         wtr_lyr[iqwv][nsd][0][0][1][im])
            REG_PROG_S("dbl_cloudwater_lyr",    wtr_lyr[iqcw][nsd][0][0][1][im])
            REG_PROG_S("dbl_rainmmr_lyr",       wtr_lyr[iqrw][nsd][0][0][1][im])
            REG_PROG_S("dbl_cloudice_lyr",      wtr_lyr[iqci][nsd][0][0][1][im])
            REG_PROG_S("dbl_snowmmr_lyr",       wtr_lyr[iqsn][nsd][0][0][1][im])
            REG_PROG_S("dbl_graupelmmr_lyr",    wtr_lyr[iqgr][nsd][0][0][1][im])
            REG_PROG_S("dbl_vapor_lyr_f",     wtr_lyr_f[iqwv][nsd][0][0][1][im])
            REG_PROG_S("dbl_cloudwater_lyr_f",wtr_lyr_f[iqcw][nsd][0][0][1][im])
            REG_PROG_S("dbl_rainmmr_lyr_f",   wtr_lyr_f[iqrw][nsd][0][0][1][im])
            REG_PROG_S("dbl_cloudice_lyr_f",  wtr_lyr_f[iqci][nsd][0][0][1][im])
            REG_PROG_S("dbl_snowmmr_lyr_f",   wtr_lyr_f[iqsn][nsd][0][0][1][im])
            REG_PROG_S("dbl_graupelmmr_lyr_f",wtr_lyr_f[iqgr][nsd][0][0][1][im])
        // }

        // one layer fields (sfc or TOA)
        REG_PHY_S("prec_tot",       pr[nsd][1][im])
        REG_PHY_S("prec_frz",     prfz[nsd][1][im])
        REG_PHY_S("olr",          flut[nsd][1][im])
        REG_PHY_S("swinc",        fsdt[nsd][1][im])
        REG_PHY_S("dbl_ts",   surfaceT[nsd][1][im])
        REG_PHY_S("dbl_gwet",     gwet[nsd][1][im])
        REG_PHY_S("dbl_zrough", zrough[nsd][1][im])
        } /* of if (enable_physics) */
    }

// Vertical levels
    gio_register_dlevel(&gio_param->data, "layers", vertical->z_lyr, km);
// Vertical interfaces
    gio_register_dlevel(&gio_param->data, "interfaces", vertical->z, km+1);

// Time integration arrays
    gio_register_index(&gio_param->dim, "nprog", time->prog_index,
                       time->nprog, 2);
    gio_register_index(&gio_param->dim, "ntend", time->tend_index,
                       time->ntend, 3);

// Restart parameters
    gio_register_restart_int(&gio_param->data, "km", km);
    gio_register_restart_int(&gio_param->data, "level_max",
                             grid_params->level_max);
    gio_register_restart_dbl(&gio_param->data, "prs_top", vertical->prs_top);
}

