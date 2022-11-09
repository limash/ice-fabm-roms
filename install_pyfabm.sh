#!/bin/bash

python_fabm_src="../fabm/src/drivers/python"

mkdir -p Build_pyfabm
rm -r Build_pyfabm/* 2> /dev/null

cmake -S $python_fabm_src -B $PWD/Build_pyfabm \
    -DCMAKE_Fortran_COMPILER=ifort 

cd Build_pyfabm
make install
cd ..

