#!/bin/bash

fabm_src="../fabm/src"

mkdir -p Build_fabm
rm -r Build_fabm/* 2> /dev/null

# Compile, build, and install FABM with ROMS as a host
cmake -S $fabm_src -B $PWD/Build_fabm \
    -DFABM_HOST=roms \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_INSTALL_PREFIX=$HOME/.local

cd Build_fabm
make install
cd ..

