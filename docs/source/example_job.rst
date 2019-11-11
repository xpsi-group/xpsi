.. _example_job:

Example job
===========

The following job script ``job.sh`` is an example from :ref:`R19`, for the
Slurm batch system on the Cartesius system (see :ref:`surfsystems`).

.. code-block:: bash

   #!/bin/bash
   #SBATCH -N 5
   #SBATCH --tasks-per-node=32
   #SBATCH -t 1-00:00:00
   #SBATCH -p broadwell
   #SBATCH --job-name=run1

   echo start of job in directory $SLURM_SUBMIT_DIR
   echo number of nodes is $SLURM_JOB_NUM_NODES
   echo the allocated nodes are:
   echo $SLURM_JOB_NODELIST

   module load intel/2017b
   module load python/2.7.9

   cp -r $HOME/NICER_analyses/J0030_ST_PST $TMPDIR

   cd $TMPDIR/J0030_ST_PST

   export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

   export OMP_NUM_THREADS=1
   export OPENBLAS_NUM_THREADS=1
   export GOTO_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   export LD_LIBRARY_PATH=$HOME/MultiNest_v3.11_CMake/multinest/lib:$LD_LIBRARY_PATH

   srun python main_run1.py > out_run1 2> err_run1

   cp run1* out_run1 err_run1 $HOME/NICER_analyses/J0030_ST_PST/.
   #end of job file

A corresponding resume script would look like:

.. code-block:: bash

    #!/bin/bash
    #SBATCH -N 30
    #SBATCH --tasks-per-node=32
    #SBATCH -t 2-00:00:00
    #SBATCH -p broadwell
    #SBATCH --job-name=run1_r1

    echo start of job in directory $SLURM_SUBMIT_DIR
    echo number of nodes is $SLURM_JOB_NUM_NODES
    echo the allocated nodes are:
    echo $SLURM_JOB_NODELIST

    module load intel/2017b
    module load python/2.7.9

    cp -r $HOME/NICER_analyses/J0030_ST_PST $TMPDIR

    cd $TMPDIR/J0030_ST_PST

    export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export GOTO_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export LD_LIBRARY_PATH=$HOME/MultiNest_v3.11_CMake/multinest/lib:$LD_LIBRARY_PATH

    srun python main_run1_resume1.py > out_run1_resume1 2> err_run1_resume1

    cp run1* out_run1_resume1 err_run1_resume1 $HOME/NICER_analyses/J0030_ST_PST/.
    #end of job file

Note how ``srun`` is aware of the MPI World, so there is no need specify the
number of processes to spawn as a flag argument. Also note that the number of
processes (which we set to equal the number of physical cores per node)
specified in the top directives is far higher than for the initial run. This
is because parallelisation efficiency scales with the local rejection fraction
during a nested sampling iteration.

Finally, note that only the root process will generate output for inspection.
