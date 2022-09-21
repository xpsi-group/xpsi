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

If X-PSI is installed using the conda environment in the Helios system (see :ref:`hpcsystems`), we can use the following type of job script:

.. code-block:: bash

   #!/bin/bash
   #SBATCH -N 2
   #SBATCH --tasks-per-node=126
   #SBATCH -t 3-00:00:00
   #SBATCH --partition=neutron-star
   #SBATCH --job-name=run1
   #SBATCH --o out
   #SBATCH --e err
   #SBATCH --mem 160000

   module purge
   module load anaconda2/2019-10
   conda activate xpsi
   module load openmpi/3.1.6

   export OMP_NUM_THREADS=1
   export OPENBLAS_NUM_THREADS=1
   export GOTO_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   export LD_LIBRARY_PATH=$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib/:$LD_LIBRARY_PATH
   export PATH=$HOME/gsl/bin:$PATH

   srun -mca btl_tcp_if_include ib0 python main_run1.py > out_run1 2> err_run1

   #end of job file

However, note that the above example is not making use of the faster communication between different nodes by using the scratch space called ``hddstore``.

If X-PSI is installed using the ``python --user`` approach, we can use the following type of job script (making now also use of the scratch space):

.. code-block:: bash

   #!/bin/bash
   #SBATCH -N 2
   #SBATCH --tasks-per-node=126
   #SBATCH -t 3-00:00:00
   #SBATCH --partition=neutron-star
   #SBATCH --job-name=run1
   #SBATCH --o out
   #SBATCH --e err
   #SBATCH --mem 160000

   module purge
   module load anaconda2/python2.7.16
   module load openmpi/3.1.6

   export OMP_NUM_THREADS=1
   export OPENBLAS_NUM_THREADS=1
   export GOTO_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages/:$PYTHONPATH
   export LD_LIBRARY_PATH=$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib/:$LD_LIBRARY_PATH
   export PATH=$HOME/gsl/bin:$PATH
   export JOB_DIR=$HOME/NICER_analyses/J0030_ST_PST/STPST_modules

   #Save output only in the primary node
   export OUTPUT_FOLDER=$(mktemp -d -p /hddstore/$USER)
   echo $OUTPUT_FOLDER $SLURMD_NODENAME
   mkdir $OUTPUT_FOLDER/run1
   cd $OUTPUT_FOLDER

   #Copy the input data to be visible for all the nodes (and make sure your paths point to hddstore):
   srun -n $SLURM_JOB_NUM_NODES --ntasks-per-node=1 mkdir -p /hddstore/$USER/data
   sleep 1
   srun -n $SLURM_JOB_NUM_NODES --ntasks-per-node=1 cp -r $HOME/NICER_analyses/data/* /hddstore/$USER/data/
   sleep 1

   mpiexec -n 252 -mca btl_tcp_if_include ib0 python $JOB_DIR/main.py @$JOB_DIR/config.ini > out1 2> err1
   sleep 1

   #Move your output from scratch to storage space.
   mkdir -p /zfs/helios/filer0/$USER/
   cp -r $OUTPUT_FOLDER/* /zfs/helios/filer0/tsalmi/J0740_GAMMA/STU_NICER_hr/

   #Clean the scratch automatically here.
   #But remember to remove manually in each node, if the main program ends by crashing.
   rm -rf $OUTPUT_FOLDER
   srun -n $SLURM_JOB_NUM_NODES --ntasks-per-node=1 rm -rf /hddstore/$USER/data

   #end of job file
