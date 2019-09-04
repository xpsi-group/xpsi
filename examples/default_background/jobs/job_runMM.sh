#!/bin/bash
#SBATCH -N 10
#SBATCH --tasks-per-node=32
#SBATCH -t 1-00:00:00
#SBATCH -p broadwell
#SBATCH --job-name=run2_MM

echo start of job in directory $SLURM_SUBMIT_DIR
echo number of nodes is $SLURM_JOB_NUM_NODES
echo the allocated nodes are:
echo $SLURM_JOB_NODELIST

module load intel/2017b
module load python/2.7.9

cp -r $HOME/xpsi_example_twospots/ $TMPDIR

cd $TMPDIR/xpsi_example_twospots/

export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=$HOME/MultiNest_v3.11_CMake/multinest/lib:$LD_LIBRARY_PATH

srun python main_run2_MM.py > out_run2_MM 2> err_run2_MM

cp run2* out_run2_MM err_run2_MM $HOME/xpsi_example_twospots/.
#end of job file
