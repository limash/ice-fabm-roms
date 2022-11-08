/*
** svn $Id: upw_ice.h 1001 2020-01-10 22:41:16Z arango $
*******************************************************************************
** Copyright (c) 2002-2020 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for upw_ice Test.
**
** Application flag:   UPW_ICE
** Input script:       ocean_upw_ice.in
*/

#define UV_ADV
#define UV_COR
#define UV_LDRAG
#define UV_VIS2
#undef  MIX_GEO_UV
#define MIX_S_UV
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define DJ_GRADPS
#define TS_DIF2
#undef  TS_DIF4
#undef  MIX_GEO_TS
#define MIX_S_TS

#define BIOLOGY
#define NEMURO

#define SALINITY
#define SOLVE3D
#define AVERAGES
#define DIAGNOSTICS_TS
#define DIAGNOSTICS_UV

#define ANA_GRID
#define ANA_INITIAL

/* surface forcing */
#ifdef SOLVE3D
# ifndef ATM2OCN_FLUXES
#  define CORE_FORCING
#  define BULK_FLUXES
#  define CCSM_FLUXES
#  undef ARCTIC_MERRA_HACK
# endif
# if defined BULK_FLUXES
#  define LONGWAVE_OUT
#  define SOLAR_SOURCE
#  define EMINUSP
#  define ANA_LRFLUX
#  define ANA_SRFLUX
#  undef ANA_ALBEDO
#  define ALBEDO
#  define ANA_SNOW
#  undef ALBEDO_CLOUD
#  undef ALBEDO_CURVE  /* for water */
#  undef ICE_ALB_EC92  /* for ice */
#  undef ALBEDO_CSIM   /* for ice */
#  undef ALBEDO_FILE  /* for both */
#  undef LONGWAVE
#  define ANA_PAIR
#  define ANA_TAIR
#  define ANA_HUMIDITY
#  define ANA_WINDS
#  define ANA_CLOUD
#  define ANA_RAIN
# endif
#endif
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX

/* ice */
#ifdef SOLVE3D
# undef CICE_MODEL
# ifdef CICE_MODEL
#  define SNOWFALL
#  define SNOW_FROM_RAIN
# endif

# define  ICE_MODEL
# undef NO_SNOW
# ifdef ICE_MODEL
#  define ANA_ICE
#  define SNOWFALL
#  define  ICE_THERMO
#  define  ICE_MK
#  undef  ICE_MOMENTUM
#  undef  ICE_MOM_BULK
#  undef  ICE_EVP
#  undef  ICE_STRENGTH_QUAD
#  define  ICE_ADVECT   /* Note that we need these two for the */
#  define  ICE_SMOLAR   /* timestepping to work correctly.     */
#  undef  ICE_UPWIND
#  define  ICE_BULK_FLUXES
#  define ICE_CONVSNOW
#  undef  MELT_PONDS
#  define ICE_I_O
# endif
#endif

#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define RI_SPLINES
#else
# define ANA_VMIX
#endif

#if defined BIO_FENNEL  || defined ECOSIM || \
    defined NPZD_POWELL || defined NEMURO
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_SRFLUX
#endif

#if defined NEMURO
# define HOLLING_GRAZING
# undef  IVLEV_EXPLICIT
#endif

#ifdef BIO_FENNEL
# define CARBON
# define DENITRIFICATION
# define BIO_SEDIMENT
# define DIAGNOSTICS_BIO
#endif

#ifdef PERFECT_RESTART
# undef  AVERAGES
# undef  DIAGNOSTICS_BIO
# undef  DIAGNOSTICS_TS
# undef  DIAGNOSTICS_UV
# define OUT_DOUBLE
#endif
