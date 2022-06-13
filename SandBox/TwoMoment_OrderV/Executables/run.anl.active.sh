#!/bin/bash
module purge

#module load oneapi/eng-compiler/2022.01.30.003
#module load oneapi/eng-compiler/2022.01.30.005
#module load oneapi/eng-compiler/2022.01.30.006
module load oneapi/eng-compiler/2022.01.30.007

module load iprof

#module load mpi/aurora_mpich/icc-sockets/41.1

export LTTNG_HOME=/localdisk/quanshao
mkdir -p $LTTNG_HOME

export LD_LIBRARY_PATH=/localdisk/quanshao/ExaStar/hdf57/lib64:$LD_LIBRARY_PATH

export THORNADO_DIR=/localdisk/quanshao/ExaStar/thornado-lab
export WEAKLIB_DIR=/localdisk/quanshao/ExaStar/weaklib
export WEAKLIB_TABLES_DIR=/localdisk/quanshao/ExaStar/weaklib-tables
export THORNADO_MACHINE=beacon_intel
#export OMP_NUM_THREADS=1

export LIBOMPTARGET_PLUGIN=LEVEL0
export SYCL_DEVICE_FILTER=LEVEL_ZERO
export LIBOMPTARGET_DEBUG=0
unset EnableWalkerPartition
export ZE_AFFINITY_MASK=0.0
#export EnableWalkerPartition=1
#export EnableImplicitScaling=1
##export ZE_AFFINITY_MASK=0.0
export LIBOMPTARGET_PLUGIN_PROFILE=T
#export OMP_TARGET_OFFLOAD=DISABLED
export OMP_TARGET_OFFLOAD=MANDATORY
#export OMP_TARGET_OFFLOAD=DISABLED
#unset OMP_TARGET_OFFLOAD
#export OMP_NUM_THREADS=1
## The following seems working well for the SineWaveStream app.
export LIBOMPTARGET_LEVEL0_MEMORY_POOL=device,16,32

ulimit -s unlimited
rm output.log

module list |& tee -a output.log

#valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-fds=yes ./ApplicationDriver_beacon_intel  |& tee -a output.val.log
( time iprof ./ApplicationDriver_beacon_intel ) |& tee -a output.log
#( time  gdb ./ApplicationDriver_beacon_intel ) |& tee output.log
#( time  gdb-oneapi ./ApplicationDriver_beacon_intel ) |& tee output.log
