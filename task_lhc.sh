#!/bin/bash

#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 5000

index=$1
main_dir=$2

source /public3/soft/modules/module.sh 
source /public3/home/c9l2000009/hexin/bagelrc
export BAGEL_DIR_QM=~/hexin/local/bagel
export BAGEL_EXE_QM=~/hexin/local/bagel/bin/BAGEL
export PSND_PYTHON=~/lhc/PSiNaD_250310/scripts/
export PSND_SCRIPTS_PATH=~/lhc/PSiNaD_250310/scripts

export OMP_NUM_THREADS=24
export BAGEL_NUM_THREADS=24

~/lhc/PSiNaD_250310/build/psinad -w -d ./$main_dir/$index -p param.json -BGIDX=$index

