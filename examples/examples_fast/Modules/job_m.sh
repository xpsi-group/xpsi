#!/bin/bash
#SBATCH -N 2
#SBATCH --ntasks-per-node=126
#SBATCH --time=30-00:00:00
#SBATCH --job-name=PT_MCMC_multinest_prior
#SBATCH --output=out_multinest_prior
#SBATCH --error=err_multinest_prior
#SBATCH --mem=0
#SBATCH --partition=neutron-star
#SBATCH --exclude=helios-cn[020,021,023]

module purge
module load gnu12
module load openmpi4
module load gsl

source "$HOME/venv311/xpsi_hack/bin/activate"

export LD_LIBRARY_PATH="$VIRTUAL_ENV/lib:$LD_LIBRARY_PATH"

# Prevent oversubscription
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMPI_MCA_pml=ob1

# Run JAX on CPU
export JAX_PLATFORMS=cpu

start=$(date +%s)

echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "MPI tasks: $SLURM_NTASKS"
echo "Start: $(date)"

mpiexec -n "$SLURM_NTASKS" python main_m.py

####python main.py --model TP_wind1 --nlive 5000 --resume


end=$(date +%s)
runtime=$((end - start))

echo "$runtime" > timingrun_real
echo "$SLURM_JOB_ID" > jobid_real

echo "All runs completed."
echo "Total runtime: $runtime seconds"
