#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=126
#SBATCH --time=00:20:00
#SBATCH --job-name=test_pool
#SBATCH --output=out_test_pool
#SBATCH --error=err_test_pool
#SBATCH --mem=0
#SBATCH --partition=neutron-star
#SBATCH --exclude=helios-cn[020,022,023]

module purge
module load gnu12
module load openmpi4
module load gsl

source "$HOME/venv311/xpsi_hack/bin/activate"
export LD_LIBRARY_PATH="$VIRTUAL_ENV/lib:$LD_LIBRARY_PATH"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export JAX_PLATFORMS=cpu

# Each variant is killed after 3 minutes if the pool hangs.

echo "=== 1. default ==="
timeout 180 mpiexec -n "$SLURM_NTASKS" python test_pool.py
echo "exit code: $?   (124 = timed out, i.e. hung)"

echo; echo "=== 2. MPROBE=0 (no matched-probe receives) ==="
MPROBE=0 timeout 180 mpiexec -n "$SLURM_NTASKS" python test_pool.py
echo "exit code: $?"

echo; echo "=== 3. OMPI_MCA_pml=ob1 (Open MPI's own transport instead of UCX) ==="
OMPI_MCA_pml=ob1 timeout 180 mpiexec -n "$SLURM_NTASKS" python test_pool.py
echo "exit code: $?"
