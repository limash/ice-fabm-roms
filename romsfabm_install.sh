#!/bin/bash

fabm_src="/home/schmiak/src/fabm/src"

mkdir -p roms-fabm-build
rm -r roms-fabm-build/* 2> /dev/null

# Compile, build, and install FABM with ROMS as a host
cmake -S $fabm_src -B $PWD/roms-fabm-build \
    -DFABM_HOST=roms \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_INSTALL_PREFIX=$PWD/local

cd roms-fabm-build
make install
cd ..

