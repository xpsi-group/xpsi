#!/bin/bash
#SBATCH -N 20
#SBATCH --tasks-per-node=32
#SBATCH -t 06:00:00
#SBATCH -p broadwell
#SBATCH --job-name=run2_MM_r1

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

srun python main_run2_MM_r1.py > out_run2_MM_r1 2> err_run2_MM_r1

cp run2* out_run2_MM_r1 err_run2_MM_r1 $HOME/xpsi_example_twospots/.
#end of job file
