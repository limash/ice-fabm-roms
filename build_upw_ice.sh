#!/bin/bash
#
# svn $Id: build_roms.bash 995 2020-01-10 04:01:28Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2020 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS/TOMS Compiling BASH Script                                       :::
#                                                                       :::
# Script to compile an user application where the application-specific  :::
# files are kept separate from the ROMS source code.                    :::
#                                                                       :::
# Q: How/why does this script work?                                     :::
#                                                                       :::
# A: The ROMS makefile configures user-defined options with a set of    :::
#    flags such as ROMS_APPLICATION. Browse the makefile to see these.  :::
#    If an option in the makefile uses the syntax ?= in setting the     :::
#    default, this means that make will check whether an environment    :::
#    variable by that name is set in the shell that calls make. If so   :::
#    the environment variable value overrides the default (and the      :::
#    user need not maintain separate makefiles, or frequently edit      :::
#    the makefile, to run separate applications).                       :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build_roms.bash [options]                                        :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#                                                                       :::
#    -p macro    Prints any Makefile macro value. For example,          :::
#                                                                       :::
#                  build.bash -p FFLAGS                                 :::
#                                                                       :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
# Notice that sometimes the parallel compilation fail to find MPI       :::
# include file "mpif.h".                                                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

which_MPI=openmpi                             # default, overwritten below

parallel=0
clean=1
dprint=0

MY_CPP_FLAGS=

while [ $# -gt 0 ]
do
  case "$1" in
    -j )
      shift
      parallel=1
      test=`echo $1 | grep '^[0-9]\+$'`
      if [ "$test" != "" ]; then
        NCPUS="-j $1"
        shift
      else
        NCPUS="-j"
      fi
      ;;

    -p )
      shift
      clean=0
      dprint=1
      debug="print-$1"
      shift
      ;;

    -noclean )
      shift
      clean=0
      ;;

    * )
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              omit argument for all avaliable CPUs"
      echo ""
      echo "-p macro    Prints any Makefile macro value"
      echo "              For example:  build.bash -p FFLAGS"
      echo ""
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
      ;;
  esac
done

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions.

export   ROMS_APPLICATION=UPW_ICE

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.

export        MY_ROOT_DIR=${HOME}/src
export     MY_PROJECT_DIR=${PWD}
export MY_COMPILE_DIR=${MY_PROJECT_DIR}

# The path to the user's local current ROMS source code.
#
# If using svn locally, this would be the user's Working Copy Path (WCPATH).
# Note that one advantage of maintaining your source code locally with svn
# is that when working simultaneously on multiple machines (e.g. a local
# workstation, a local cluster and a remote supercomputer) you can checkout
# the latest release and always get an up-to-date customized source on each
# machine. This script is designed to more easily allow for differing paths
# to the code and inputs on differing machines.

export       MY_ROMS_SRC=${MY_PROJECT_DIR}/ROMS_v37_KATE

NEW_SRC_DIR=${MY_ROMS_SRC}/ROMS/Nonlinear/Biology

# Set path of the directory containing makefile configuration (*.mk) files.
# The user has the option to specify a customized version of these files
# in a different directory than the one distributed with the source code,
# ${MY_ROMS_SRC}/Compilers. If this is the case, you need to keep these
# configurations files up-to-date.

 export         COMPILERS=${MY_ROMS_SRC}/Compilers
#export         COMPILERS=${HOME}/Compilers/ROMS

# # Check if we have any modified source files

if [ -s Src_modify ]; then
    cd Src_modify/
#    gotModifiedSource=`ls *.F *.h *.mk`
    gotModifiedSource=`ls *.F *.h *.mk makefile`
###Modified PWA 04/03/2015
    cd ..
fi

# Replace the original files with the modifications
if [ "$gotModifiedSource" != "" ]; then

    # Copy locally modified source to main ROMS directory
    for ModSrc_modify in $gotModifiedSource; do

        # Check where original resides
###PWA Modified 04/03/2015
        #	origFile=`find $MY_ROMS_SRC -name $ModSrc_modify
        if [ "$ModSrc_modify" == "makefile" ]; then
        # For the makefile, consider only the makefile in $MY_ROMS_SRC (PWA)
            echo "Copying makefile"
            echo $ModSrc_modify
            origFile=$MY_ROMS_SRC/makefile
        else
            origFile=`find $MY_ROMS_SRC -name $ModSrc_modify`
            if [ "$origFile" != "" ]; then
                echo "Found original file:"
                echo $origFile
            fi
        fi
###PWA Modified 04/03/2015

        if [ -f "$origFile" ]; then

            # Moving original and copying user-modifed source code
            # first checking if the original already exists with
            # the .orig extension
            if [ ! -f "$origFile.orig" ]; then
                mv $origFile $origFile.orig
                echo "Moving $origFile to $origFile.orig"
            fi

            # Copying from local source directory to repository
            cp Src_modify/$ModSrc_modify $origFile
            echo "Copying Src_modify/$ModSrc_modify to $origFile"

            if [ ! -f USER_MODIFIED_CODE_IN_REPO ]; then

                # Touch file to notify that user modified code has been
                # placed in the repository
                touch USER_MODIFIED_CODE_IN_REPO

            fi
        else

	    ##PWA changes

        # No such file in repository, proceed to copy in new file (PWA)
            cp Src_modify/$ModSrc_modify $NEW_SRC_DIR/.
            echo "Copying Src_modify/$ModSrc_modify to $NEW_SRC_DIR"

            if [ ! -f USER_MODIFIED_CODE_IN_REPO ]; then

                # Touch file to notify that user modified code has been
                # placed in the repository
                touch USER_MODIFIED_CODE_IN_REPO

            fi
	    ##End of PWA changes
        fi
    done
fi
# Removing user modified source code in repository
# KHC - 20110209
# NMK - 2013
rollback() {
    cd $MY_COMPILE_DIR

    if [ -f USER_MODIFIED_CODE_IN_REPO ]; then

    # Find source code files with ".orig"-ending and
    # remove ending
        filelist=`find "$MY_ROMS_SRC" -name *.orig`

        if [ "$filelist" != "" ]; then

            for oldFileName in $filelist; do

            # extract basename
            newFileName=`basename $oldFileName .orig`
            fileDirectory=`dirname $oldFileName`
            mv $oldFileName  $fileDirectory/$newFileName

            echo "Moved $oldFileName  to $fileDirectory/$newFileName"

            done

        else # Empty filelist, no such files in repository

            echo "Did not find any .orig-files in the repository, empty file deleted"

        fi

        # Remove empty file

        rm -f USER_MODIFIED_CODE_IN_REPO

    fi
}
trap 'rollback; exit 99' 0

#--------------------------------------------------------------------------
# Set tunable CPP options.
#--------------------------------------------------------------------------
#
# Sometimes it is desirable to activate one or more CPP options to run
# different variants of the same application without modifying its header
# file. If this is the case, specify each options here using the -D syntax.
# Notice also that you need to use shell's quoting syntax to enclose the
# definition.  Both single or double quotes work. For example,
#
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DAVERAGES"
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DDEBUGGING"
#
# can be used to write time-averaged fields. Notice that you can have as
# many definitions as you want by appending values.

#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -D"

#--------------------------------------------------------------------------
# Compiler options.
#--------------------------------------------------------------------------
#
# Other user defined environmental variables. See the ROMS makefile for
# details on other options the user might want to set here. Be sure to
# leave the switches meant to be off set to an empty string or commented
# out. Any string value (including off) will evaluate to TRUE in
# conditional if-statements.

#export           USE_MPI=on            # distributed-memory parallelism
#export        USE_MPIF90=on            # compile with mpif90 script
#export         which_MPI=mpich         # compile with MPICH library
#export         which_MPI=mpich2        # compile with MPICH2 library
#export         which_MPI=mvapich2      # compile with MVAPICH2 library
#export         which_MPI=openmpi       # compile with OpenMPI library

#export        USE_OpenMP=on            # shared-memory parallelism

 export              FORT=ifort
#export              FORT=gfortran
#export              FORT=pgi

 export         USE_DEBUG=on            # use Fortran debugging flags
 export         USE_LARGE=on            # activate 64-bit compilation
 export       USE_NETCDF4=on            # compile with NetCDF-4 library
#export   USE_PARALLEL_IO=on            # Parallel I/O with NetCDF-4/HDF5


#--------------------------------------------------------------------------
# If Earth System Model (ESM) coupling, set location of ESM component
# libraries and modules. The strategy is to compile and link each ESM
# component separately first and ROMS last since it is driving the
# coupled system. Only the ESM components activated are considered and
# the rest ignored.
#--------------------------------------------------------------------------

export        WRF_SRC_DIR=${HOME}/ocean/repository/WRF

if [ -n "${USE_DEBUG:+1}" ]; then
  export     CICE_LIB_DIR=${MY_PROJECT_DIR}/Build_ciceG
  export   COAMPS_LIB_DIR=${MY_PROJECT_DIR}/Build_coampsG
  export    REGCM_LIB_DIR=${MY_PROJECT_DIR}/Build_regcmG
  export      WAM_LIB_DIR=${MY_PROJECT_DIR}/Build_wamG
# export      WRF_LIB_DIR=${MY_PROJECT_DIR}/Build_wrfG
  export      WRF_LIB_DIR=${WRF_SRC_DIR}
else
  export     CICE_LIB_DIR=${MY_PROJECT_DIR}/Build_cice
  export   COAMPS_LIB_DIR=${MY_PROJECT_DIR}/Build_coamps
  export    REGCM_LIB_DIR=${MY_PROJECT_DIR}/Build_regcm
  export      WAM_LIB_DIR=${MY_PROJECT_DIR}/Build_wam
  export      WRF_LIB_DIR=${MY_PROJECT_DIR}/Build_wrf
# export      WRF_LIB_DIR=${WRF_SRC_DIR}
fi

#--------------------------------------------------------------------------
# If applicable, use my specified library paths.
#--------------------------------------------------------------------------

 export USE_MY_LIBS=no            # use system default library paths
#export USE_MY_LIBS=yes           # use my customized library paths

MY_PATHS=${COMPILERS}/my_build_paths.bash

if [ "${USE_MY_LIBS}" = "yes" ]; then
  source ${MY_PATHS} ${MY_PATHS}
fi

#--------------------------------------------------------------------------
# The rest of this script sets the path to the users header file and
# analytical source files, if any. See the templates in User/Functionals.
#--------------------------------------------------------------------------
#
# If applicable, use the MY_ANALYTICAL_DIR directory to place your
# customized biology model header file (like fennel.h, nemuro.h, ecosim.h,
# etc).

 export     MY_HEADER_DIR=${MY_ROMS_SRC}/Apps/Upw_ice

#export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}

# Put the binary to execute in the following directory.

 export            BINDIR=${MY_PROJECT_DIR}/Case_upw_ice

# Put the f90 files in a project specific Build directory to avoid conflict
# with other projects.

if [ -n "${USE_DEBUG:+1}" ]; then
 export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build_romsG
else
 export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build_roms
fi

###PWA Inserted 04/03/2015
export FABM_INCDIR=/home/schmiak/.local/include
export FABM_LIBDIR=/home/schmiak/.local/lib
#This defines the location of the FABM *.mod files and library,
#for use in the makefile if needed.
###PWA Inserted 04/03/2015

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${MY_ROMS_SRC}

# Stop if activating both MPI and OpenMP at the same time.

if [ -n "${USE_MPI:+1}" ] && [ -n "${USE_OpenMP:+1}" ]; then
  echo "You cannot activate USE_MPI and USE_OpenMP at the same time!"
  exit 1
fi

#--------------------------------------------------------------------------
# Compile.
#--------------------------------------------------------------------------

# Remove build directory.

if [ $clean -eq 1 ]; then
  make clean
fi

# Compile (the binary will go to BINDIR set above).

if [ $dprint -eq 1 ]; then
  make $debug
else
  if [ $parallel -eq 1 ]; then
    make $NCPUS
  else
    make
  fi
fi
