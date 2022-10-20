#!/bin/bash

python_fabm_src="../fabm/src/drivers/python"

mkdir -p python-fabm-build
rm -r python-fabm-build/* 2> /dev/null

cmake -S $python_fabm_src -B $PWD/python-fabm-build \
    -DCMAKE_Fortran_COMPILER=gfortran 

cd python-fabm-build
make install
cd ..

