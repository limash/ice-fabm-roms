# svn $Id: Linux-ftn.mk 199 2008-07-25 18:57:58Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2008 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for CRAY FTN cross-compiler with Linux
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# CLEAN          Name of cleaning executable after C-preprocessing
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF libary directory
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := ifort
#           FFLAGS := -e I -e m
           FFLAGS := -fPIC
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional
               LD := $(FC)
          LDFLAGS :=
               AR := ar
          ARFLAGS := -r
            MKDIR := mkdir -p
               RM := rm -f
           RANLIB := touch
	     PERL := perl
             TEST := test

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#
# Library locations, can be overridden by environment variables.
#

ifdef USE_NETCDF4
    NETCDF_INCDIR ?= /global/hds/software/cpu/eb1/netCDF-Fortran/4.4.3-intel-2016a/include
    NETCDF_LIBDIR ?= /global/hds/software/cpu/eb1/netCDF/4.4.0-intel-2016a/lib
    HDF5_LIBDIR ?= /global/hds/software/cpu/eb1/HDF5/1.8.16-intel-2016a/lib
else
    NETCDF_INCDIR ?= $(NETCDF_DIR)/include
    NETCDF_LIBDIR ?= $(NETCDF_DIR)/lib
endif

             LIBS := -L$(NETCDF_LIBDIR) -lmpifort -lmpi -lnetcdff
ifdef USE_NETCDF4
             LIBS += -L$(HDF5_LIBDIR) -lhdf5_hl -lhdf5 -lz
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /work/apps/arpack/2002.08.23-pgi/lib
             LIBS += -L$(PARPACK_LIBDIR) 
 endif
    ARPACK_LIBDIR ?= /work/apps/arpack/2002.08.23-pgi/lib
             LIBS += -L$(ARPACK_LIBDIR) 
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
endif

ifdef USE_DEBUG
#           FFLAGS += -G 0
#           FFLAGS += -C -g -Mchkptr -Mchkstk -Ktrap=fp -traceback
           FFLAGS += -g
else
#           FFLAGS += -O 3,aggress
           FFLAGS += -O3
endif

ifdef USE_MCT
       MCT_INCDIR ?= /usr/local/mct/include
       MCT_LIBDIR ?= /usr/local/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif

ifdef USE_ESMF
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) -lesmf -lC
endif

###PWA Inserted 03/04/2017
ifdef USE_FABM
             LIBS += -L$(FABM_LIBDIR) -lfabm
# Remember that this means the compiler will look for 'libfabm.a'
endif
###PWA Inserted 03/04/2017

#
# Use full path of compiler.
#
               FC := $(shell which ${FC})

#
# Set free form format in source files to allow long string for
# local directory and compilation flags inside the code.
#

$(SCRATCH_DIR)/analytical.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/mod_strings.o: FFLAGS := $(MY_FFLAGS) -Mfree

#
# Supress free format in SWAN source files since there are comments
# beyond column 72.
#

ifdef USE_SWAN

$(SCRATCH_DIR)/ocpcre.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/ocpids.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/ocpmix.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanparll.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanser.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swmod1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swmod2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/m_constants.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/m_fileio.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/mod_xnl4v5.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/serv_xnl4v5.o: FFLAGS += -Mfree

endif
