#!/bin/sh
#SBATCH -N 1
#SBATCH --tasks-per-node=48
#SBATCH -t 3-00:00:00
#SBATCH -J synt1df
#SBATCH -o out # STDOUT
#SBATCH -e err # STDERR
#SBATCH --partition=neutron-star
#SBATCH --mem 252000


# ##SBATCH --nodelist=helios-cn024
start=$(date +%s)
echo $start >startsts

#RUN THESE 2 lines IN TERMINAL BEFORE SUBMITTING THE JOB:  
# module load anaconda3  
# source activate xpsi 

module load mpi
export OMPI_MCA_mpi_warn_on_fork=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd home/tsalmi/xpsi/xpsi_dev/xpsi/

export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH
#export PYTHONPATH="/home/tsalmi/anaconda3/envs/xpsi/lib/python2.7/site-packages/xpsi-0.7.5-py2.7-linux-x86_64.egg"

export LD_LIBRARY_PATH=$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib:$LD_LIBRARY_PATH

mpiexec -n 48 python compute_job.py > out_run1f 2> err_run1f
#mpiexec -n 144 python main.py @config.ini --multinest

cp out_run1f err_run1f home/tsalmi/xpsi/xpsi_dev/xpsi/run_test1f/.

end=$(date +%s)
echo $end >endsts
runtime=$((end-start))
echo $runtime >timingsts
#end of job file
