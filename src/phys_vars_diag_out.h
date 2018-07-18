/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: phys_vars_diag_out.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_PHYS_VARS_DIAG_OUT
#define H_PHYS_VARS_DIAG_OUT

typedef struct {
/*=======================================================================
!This module stores the physics variables we wish to save as
!  diagnostic output.
!=======================================================================*/

    /* km_phys is actually defined in MODULE phys_vertical */
    int km_phys;

    int l_accum;

    /* microphysics diagnostics */
    int count_micro_var;
    double 
       ***          pr, /* (im,jm,        nsdm) total precipitation                      (kg/m^2/s) */
       ***        prfz, /* (im,jm,        nsdm) frozen precipitation                     (kg/m^2/s) */
      ****    q_latent, /* (im,jm,km_phys,nsdm) latent heating rate                      (K/d)      */
      ****src_local_qv, /* (im,jm,km_phys,nsdm) local microphysical water vapor tendency (kg/kg/s)  */
      ****src_local_qc, /* (im,jm,km_phys,nsdm) local microphysical cloud water tendency (kg/kg/s)  */
      ****src_local_qi, /* (im,jm,km_phys,nsdm) local microphysical cloud ice tendency   (kg/kg/s)  */
      ****src_local_qs, /* (im,jm,km_phys,nsdm) local microphysical snow tendency        (kg/kg/s)  */
      ****src_local_qr, /* (im,jm,km_phys,nsdm) local microphysical rain tendency        (kg/kg/s)  */
      ****src_local_qg, /* (im,jm,km_phys,nsdm) local microphysical graupel tendency     (kg/kg/s)  */
      **** src_fall_qs, /* (im,jm,km_phys,nsdm) snow tendency due to falling precip      (kg/kg/s)  */
      **** src_fall_qr, /* (im,jm,km_phys,nsdm) rain tendency due to falling precip      (kg/kg/s)  */
      **** src_fall_qg; /* (im,jm,km_phys,nsdm) graupel tendency due to falling precip    kg/kg/s)  */
   
    /* radiation diagnostics */
    int count_rad_var;
    double
      ****   q_lwrad, /* (im,jm,km_phys,nsdm) lw radiative heating rate, all-sky     (K/d)         */
      ****   q_swrad, /* (im,jm,km_phys,nsdm) sw radiative heating rate, all-sky     (K/d)         */
       ***      flut, /* (im,jm,        nsdm) lw flux, toa upward component, all-sky (W/m^2) <olr> */
       ***      flus, /* (im,jm,        nsdm) lw flux, sfc upward component, all-sky (W/m^2)       */
       ***      flds, /* (im,jm,        nsdm) lw flux, sfc dwnwrd component, all-sky (W/m^2)       */
       ***      fsdt, /* (im,jm,        nsdm) sw flux, toa dwnwrd component          (W/m^2)       */
       ***      fsut, /* (im,jm,        nsdm) sw flux, toa upward component, all-sky (W/m^2)       */
       ***      fsds, /* (im,jm,        nsdm) sw flux, sfc dwnwrd component, all-sky (W/m^2)       */
       ***      fsus, /* (im,jm,        nsdm) sw flux, sfc upward component, all-sky (W/m^2)       */
      ****q_lwrad_cs, /* (im,jm,km_phys,nsdm) lw radiative heating rate, clear-sky   (K/d)         */
      ****q_swrad_cs, /* (im,jm,km_phys,nsdm) sw radiative heating rate, clear-sky   (K/d)         */
       ***   flut_cs, /* (im,jm,        nsdm) lw flux, toa upward component, clr-sky (W/m^2)       */
       ***   flus_cs, /* (im,jm,        nsdm) lw flux, sfc upward component, clr-sky (W/m^2)       */
       ***   flds_cs, /* (im,jm,        nsdm) lw flux, sfc dwnwrd component, clr-sky (W/m^2)       */
       ***   fsut_cs, /* (im,jm,        nsdm) sw flux, toa upward component, clr-sky (W/m^2)       */
       ***   fsds_cs, /* (im,jm,        nsdm) sw flux, sfc dwnwrd component, clr-sky (W/m^2)       */
       ***   fsus_cs; /* (im,jm,        nsdm) sw flux, sfc upward component, clr-sky (W/m^2)       */
   
    /* turbulence diagnostics */
    int count_turb_var;
    double
       ***        evap, /* (im,jm,        nsdm) sfc evaporation                    (kg/m^2/s) */
       ***         fss, /* (im,jm,        nsdm) sfc sensible heat flux             (W/m^2)    */
      ****   q_tht_sgs, /* (im,jm,km_phys,nsdm) sgs turbulent heating rate         (K/d)      */
      ****qwv_sgs_tend, /* (im,jm,km_phys,nsdm) sgs turbulent water vapor tendency (kg/kg/s)  */
      ****qcw_sgs_tend, /* (im,jm,km_phys,nsdm) sgs turbulent cloud water tendency (kg/kg/s)  */
      ****qci_sgs_tend; /* (im,jm,km_phys,nsdm) sgs turbulent cloud ice tendency   (kg/kg/s)  */
   

    /* below 3 are actually defined in MODULE phys_vars_diagnostic */
    double
       ***surfaceT,   /* (im,jm,nsdm) lower boundary temperature (K) */
       ***zrough,     /* (im,jm,nsdm) lower boundary roughness length (m) */
       ***gwet;       /* (im,jm,nsdm) lower boundary wetness (-) */

} MODULE_phys_vars_diag_out;

/* API declarations */
void init_MODULE_phys_vars_diag_out(MODULE_phys_vars_diag_out *vars,
                                    int nsdm, int km_phys, int jm, int im);
void finalize_MODULE_phys_vars_diag_out(MODULE_phys_vars_diag_out *vars);
void init_phys_vars_diag_out(MODULE_phys_vars_diag_out *vars, int nsdm, int km_phys, int jm, int im);
void norm_phys_vars_diag_out(MODULE_phys_vars_diag_out *vars, int nsdm, int km_phys, int jm, int im);
void randomize_phys_vars_diag_out(MODULE_phys_vars_diag_out *vars,
                                  int                        nsdm,
                                  int                        km_phys,
                                  int                        jm,
                                  int                        im);

#endif
