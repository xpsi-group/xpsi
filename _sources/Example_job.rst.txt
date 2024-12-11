.. _example_job:

Example job
===========

For both jobs you will need these 
`auxiliary files <https://zenodo.org/record/7113931>`_ inside the ``model_data/``
directory.

Snellius
--------

The following job script ``job.sh`` is an example job script for analysis on the Snellius system (see :ref:`hpcsystems`).

.. code-block:: bash

    #!/bin/bash
    #SBATCH -N 5
    #SBATCH --tasks-per-node=128
    #SBATCH -t 1-00:00:00
    #SBATCH -p thin
    #SBATCH --job-name=run1
    #SBATCH --mail-user=my_email@gmail.com
    #SBATCH --mail-type=END    

    echo start of job in directory $SLURM_SUBMIT_DIR
    echo number of nodes is $SLURM_JOB_NUM_NODES
    echo the allocated nodes are:
    echo $SLURM_JOB_NODELIST

    module purge

    module load 2024
    module load foss/2024a
    module load SciPy-bundle/2024.05-gfbf-2024a
    module load wrapt/1.16.0-gfbf-2024a
    module load matplotlib/3.9.2-gfbf-2024a

    source $HOME/venvs/xpsi_py3/bin/activate
    
    cp -r $HOME/xpsi/examples/examples_modeling_tutorial/* $TMPDIR/
    mkdir $TMPDIR/run
    cd $TMPDIR/

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export GOTO_NUM_THREADS=1
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/multinest/MultiNest_v3.12_CMake/multinest/lib/

    srun python TestRun_BB.py > out1 2> err1
    cp -r out1 err1 run $HOME/xpsi/examples/examples_modeling_tutorial/.
    #end of job file

This can be submitted with ``sbatch job.sh``. In a corresponding resume script one would typically increase the number of processes. This is because parallelisation efficiency scales with the local rejection fraction during a nested sampling iteration. Make sure that also the previous sample output files are copied to the $TMPDIR directory, and that the resume option for MultiNest is turned on in the main Python file (in ``TestRun_BB.py`` for this example).

Note how ``srun`` is aware of the MPI World, so there is no need specify the
number of processes to spawn as a flag argument.

Finally, note that only the root process will generate output for inspection.

Helios
------

For Helios, we can use the following type of job script:

.. code-block:: bash

    #!/bin/bash
    #SBATCH -N 2
    #SBATCH --tasks-per-node=126
    #SBATCH -t 1-00:00:00
    #SBATCH -J example_run
    #SBATCH -o example_run_%j.out
    #SBATCH -e example_run_%j.err
    #SBATCH --partition=neutron-star
    #SBATCH --mem 160000
    #SBATCH --mail-user=my_email@gmail.com
    #SBATCH --mail-type=END 

    module purge
    module load gnu12
    module load openmpi4
    module load gsl
    source $HOME/venv311/xpsi/bin/activate

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export GOTO_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib:$LD_LIBRARY_PATH
    export PATH=$HOME/gsl/bin:$PATH

    export JOB_DIR=$HOME/xpsi/examples/examples_modeling_tutorial
    export OUTPUT_FOLDER=$(mktemp -d -p /hddstore/$USER)
    echo $OUTPUT_FOLDER $SLURMD_NODENAME
    mkdir $OUTPUT_FOLDER/run
    cd $OUTPUT_FOLDER

    #Copy the input data to be visible for all the nodes (and make sure your paths point to hddstore):
    srun -n $SLURM_JOB_NUM_NODES --ntasks-per-node=1 cp -r $JOB_DIR/model_data/ $OUTPUT_FOLDER 
    sleep 1

    mpiexec -n 252 -mca btl_tcp_if_include ib0 python $JOB_DIR/TestRun_BB.py

    #Move your output from scratch to storage space.
    mkdir -p /zfs/helios/filer0/$USER/
    cp -r $OUTPUT_FOLDER/* /zfs/helios/filer0/$USER/

    #Clean the scratch automatically here.
    #But remember to remove manually in each node, if the main program ends by crashing.
    rm -rf $OUTPUT_FOLDER
    
