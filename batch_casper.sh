#!/bin/bash
#PBS -N MICM 
#PBS -A NTDD0004
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=300GB:ngpus=1
#PBS -l gpu_type=v100
#PBS -l walltime=01:59:00
#PBS -q casper 
#PBS -j oe
#PBS -k eod

ulimit -s unlimited

cd $PBS_O_WORKDIR

# copy the preprocessor chemistry solver to the src folder
cp test/preprocessor_output/*F90 src/preprocessor_output/

# clean up the build folder
if [ ! -d "./build" ]
then
   mkdir build
fi
cd build
rm -rf *

###################################################################################################
# Set some options on Casper following Matt's suggestions: https://github.com/NCAR/micm/issues/16 #
# Note that using CMAKE options may drop some linkages and lead to a compilation failure          #
###################################################################################################
# unload any modules currently loaded
module purge

# load modules
module load ncarenv/1.3
# module load intel/19.1.1
module load nvhpc/22.2
module load openmpi/4.1.1
module load ncarcompilers/0.5.0
module load cmake/3.22.0

export JSON_FORTRAN_HOME=/glade/scratch/sunjian/temp/json-fortran-8.3.0/build
export NETCDF_HOME=/glade/u/apps/dav/opt/netcdf/4.8.1/nvhpc/22.2

# build a MICM test
#cmake -D ENABLE_UTIL_ONLY=ON ..
#cmake -D ENABLE_UTIL_ONLY=ON -D ENABLE_NETCDF=ON ..
#cmake -D ENABLE_UTIL_ONLY=ON -D ENABLE_NSYS=ON ..
#cmake -D ENABLE_UTIL_ONLY=ON -D ENABLE_OPENACC=OFF ..
cmake -D ENABLE_UTIL_ONLY=ON -D ENABLE_NETCDF=ON -D ENABLE_OPENACC=OFF ..
#cmake -D ENABLE_UTIL_ONLY=ON -D ENABLE_OPENACC=OFF -D CMAKE_BUILD_TYPE=DEBUG ..
make VERBOSE=1       # VERBOSE shows whether the desired flags are applied or not

# run a MICM test
make test

# save the output to a desired directory
outdir="../output"
if [ ! -d $outdir ]
then
   mkdir $outdir
fi
#mv ./Testing/Temporary/LastTest.log $outdir/gpu1000000.log 
