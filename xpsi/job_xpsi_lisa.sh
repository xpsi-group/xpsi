#!/bin/bash
#SBATCH -N 4
#SBATCH --tasks-per-node=8
#SBATCH -t 3-00:00:00
#SBATCH -p normal
#SBATCH --job-name=syntd1

echo start of job in directory $SLURM_SUBMIT_DIR
echo number of nodes is $SLURM_JOB_NUM_NODES
echo the allocated nodes are:
echo $SLURM_JOB_NODELIST

module load 2019
module load intel/2019b
module load Python/2.7.15-intel-2019b

## cp -r $HOME/NICER_analyses/J0030_ST_PST $TMPDIR

## cd $TMPDIR/J0030_ST_PST

cd home/tsalmi/xpsi/xpsi_dev/xpsi/

export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
## export LD_LIBRARY_PATH=$HOME/MultiNest_v3.11_CMake/multinest/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/tsalmi/multinest/MultiNest_v3.12_CMake/multinest/lib/:$LD_LIBRARY_PATH
export PATH=/home/tsalmi/gsl/bin:$PATH


## srun python main_run1.py > out_run1 2> err_run1

srun python compute_job.py > out_run1 2> err_run1

## cp run1* out_run1 err_run1 $HOME/NICER_analyses/J0030_ST_PST/.

cp out_run1 err_run1 home/tsalmi/xpsi/xpsi_dev/xpsi/run_test1/.

#end of job file
