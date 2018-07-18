/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: physical_params.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_PHYSICAL_PARAMS
#define H_PHYSICAL_PARAMS

#if HAVE_CONFIG_H
#include "config.h"
#endif

#define      A            6.371229E+06
#define      GRAV         9.80616
#define      OMEGA        7.29211E-5
#define      GAS_CONST_R  287.04
#define      SPEC_HEAT_CV 717.60
#define      SPEC_HEAT_CP 1004.64
#define      KAPPA        GAS_CONST_R/SPEC_HEAT_CP
#define      TICE         273.15
#define      HLTM         2.52E+06
#define      HLTF         3.336E+05
#define      VK           0.40
#define      DELTA        0.608
#define      MWDRY        28.991
#define      P0_SFC       100000.0

#endif
