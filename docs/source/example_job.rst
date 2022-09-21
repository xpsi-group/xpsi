.. _example_job:

Example job
===========

Snellius
--------

The following job script ``job.sh`` is an example job script used by the X-PSI core team for NICER data analysis on the Snellius system (see :ref:`hpcsystems`).

.. code-block:: bash

    #!/bin/bash
    #SBATCH -N 5
    #SBATCH --tasks-per-node=128
    #SBATCH -t 1-00:00:00
    #SBATCH -p thin
    #SBATCH --job-name=run1

    echo start of job in directory $SLURM_SUBMIT_DIR
    echo number of nodes is $SLURM_JOB_NUM_NODES
    echo the allocated nodes are:
    echo $SLURM_JOB_NODELIST

    module load 2021
    module load intel/2021a
    module load Python/2.7.18-GCCcore-10.3.0-bare

    cp -r $HOME/NICER_analyses/J0437/CST_PDT/ "$TMPDIR"

    export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export GOTO_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export LD_PRELOAD=/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_def.so.1:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_avx2.so.1:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_core.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_intel_lp64.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_intel_thread.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/compiler/2021.2.0/linux/compiler/lib/intel64_lin/libiomp5.so
    export LD_LIBRARY_PATH=$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib:$LD_LIBRARY_PATH

    srun python main_run1.py @config.ini > out_run1 2> err_run1

    cp run1* out_run1 err_run1 $HOME/NICER_analyses/J0437/CST_PDT/CST_PDT_outputs/init/.
    #end of job file

A corresponding resume script would look like:

.. code-block:: bash

    #!/bin/bash
    #SBATCH -N 30
    #SBATCH --tasks-per-node=128
    #SBATCH -t 2-00:00:00
    #SBATCH -p thin
    #SBATCH --job-name=run1_r1

    echo start of job in directory $SLURM_SUBMIT_DIR
    echo number of nodes is $SLURM_JOB_NUM_NODES
    echo the allocated nodes are:
    echo $SLURM_JOB_NODELIST

    module load 2021
    module load intel/2021a
    module load Python/2.7.18-GCCcore-10.3.0-bare

    cp -r $HOME/NICER_analyses/J0437/CST_PDT/ "$TMPDIR"

    export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export GOTO_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export LD_PRELOAD=/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_def.so.1:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_avx2.so.1:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_core.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_intel_lp64.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/mkl/2021.2.0/lib/intel64/libmkl_intel_thread.so:/sw/arch/Centos8/EB_production/2021/software/imkl/2021.2.0-iimpi-2021a/compiler/2021.2.0/linux/compiler/lib/intel64_lin/libiomp5.so
    export LD_LIBRARY_PATH=$HOME/multinest/MultiNest_v3.12_CMake/multinest
    export LD_LIBRARY_PATH=$HOME/MultiNest_v3.11_CMake/multinest/lib:$LD_LIBRARY_PATH

    srun python main_run1_resume1.py @config.ini> out_run1_resume1 2> err_run1_resume1

    cp run1* out_run1_resume1 err_run1_resume1 $HOME/NICER_analyses/J0437/CST_PDT/CST_PDT_outputs/resume1/.
    #end of job file

Note how ``srun`` is aware of the MPI World, so there is no need specify the
number of processes to spawn as a flag argument. Also note that the number of
processes (which we set to equal the number of physical cores per node)
specified in the top directives is far higher than for the initial run. This
is because parallelisation efficiency scales with the local rejection fraction
during a nested sampling iteration.

Finally, note that only the root process will generate output for inspection.

Helios
------
