#!/bin/bash
module purge

#module load oneapi/eng-compiler/2022.01.30.003
#module load oneapi/eng-compiler/2022.01.30.005
#module load oneapi/eng-compiler/2022.01.30.006
module load oneapi/eng-compiler/2022.01.30.007


export HDF5_INC=/localdisk/quanshao/ExaStar/hdf57/include
export HDF5_LIB=/localdisk/quanshao/ExaStar/hdf57/lib64

export THORNADO_DIR=/localdisk/quanshao/ExaStar/thornado-lab
export WEAKLIB_DIR=/localdisk/quanshao/ExaStar/weaklib
export WEAKLIB_TABLES_DIR=/localdisk/quanshao/ExaStar/weaklib-tables
export THORNADO_MACHINE=beacon_intel
#export OMP_NUM_THREADS=1

rm ./ApplicationDriver_beacon_intel
rm make.ifx.log
make clean

module list |& tee -a make.ifx.log

#time make ApplicationDriver USE_OMP_OL=TRUE |& tee -a make.log
#time make ApplicationDriver USE_OMP=TRUE |& tee -a make.ifort.log
#( time make ApplicationDriver USE_GPU=TRUE USE_OMP_OL=TRUE USE_ONEMKL=TRUE ) |& tee -a make.ifx.log
#( time make ApplicationDriver USE_OMP_OL=TRUE USE_GPU=TRUE USE_CUDA=FALSE USE_ONEMKL=TRUE ) |& tee -a make.ifx.log
( time make ApplicationDriver USE_OMP_OL=TRUE USE_GPU=TRUE USE_CUDA=FALSE USE_ONEMKL=TRUE ) |& tee -a make.ifx.log
