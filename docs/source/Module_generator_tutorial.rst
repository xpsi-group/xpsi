.. _module_generator_tutorial:

Module generator tutorial
=========================

To develop models in X-PSI, users are recommended to to create custom derived classes for the different X-PSI model
components (e.g. ``Prior``, ``Photosphere``, etc.) that require further modification. Examples of such `Custom` classes
can be found in :ref:`./Modeling.ipynb`. These classes are preferably to be written as distinct modules within a project
directory rather than in a notebook as is done in :ref:`./Modeling.ipynb`. The modules can then be imported by a `Main`
script, which can be executed using an MPI command in a shell to exploit parallelism for resource-intensive likelihood
evaluations.

The structure of these `Main` and `Custom` modules are largely the same for most user-defined models. It is usually the
scope of these models (e.g. hotspot models, required data and instrument files, prior definitions, etc.) that are
subject to change. Therefore, to simplify the process and avoid the user from having to create these modules every time
for a new model, we have created the program ``module_generator`` which takes in arguments to define the scope of the
model. Below, we provide instructions on how to use it.

Create a symlink to the module generator in your project directory. You can run the help command as shown below to
obtain details regarding its usage.

.. code-block:: bash

    ln -s <path/to/xpsi>/xpsi/module_generator.py ./
    python module_generator.py -h

The output of the above code block is:

.. code-block:: bash

    /=============================================\
    | X-PSI: X-ray Pulse Simulation and Inference |
    |---------------------------------------------|
    |                Version: 2.0.0               |
    |---------------------------------------------|
    |      https://xpsi-group.github.io/xpsi      |
    \=============================================/

    Imported GetDist version: 1.4
    Imported nestcheck version: 0.2.1
    Parsing configuration file...
    usage: module_generator.py [-h] [--generate-config-file]
                               [--telescope TELESCOPE] [--instrument INSTRUMENT]
                               [--source SOURCE] --frequency FREQUENCY
                               [--model MODEL] --hot-region-model
                               {ST,CST,EST,PST,CDT,EDT,PDT}
                               [--antipodal-reflection-symmetry]
                               [--break-hot-region-exchange-degeneracy-with BREAK_HOT_REGION_EXCHANGE_DEGENERACY_WITH]
                               [--is-antiphased IS_ANTIPHASED] [--prefix PREFIX]
                               [--hot-atmosphere-model HOT_ATMOSPHERE_MODEL]
                               [--hot-atmosphere-load]
                               [--elsewhere-atmosphere-model ELSEWHERE_ATMOSPHERE_MODEL]
                               [--elsewhere-atmosphere-load]
                               [--attenuation-model ATTENUATION_MODEL]
                               [--background-model] [--background-shared-instance]
                               [--background-shared-class]
                               [--background-parameters [BACKGROUND_PARAMETERS ...]]
                               [--print-MPI-rank] [--config-path CONFIG_PATH]
                               [--module-directory-path MODULE_DIRECTORY_PATH]
                               [--main-module MAIN_MODULE]
                               [--custom-signal-module CUSTOM_SIGNAL_MODULE]
                               [--custom-instrument-module CUSTOM_INSTRUMENT_MODULE]
                               [--custom-photosphere-module CUSTOM_PHOTOSPHERE_MODULE]
                               [--custom-interstellar-module CUSTOM_INTERSTELLAR_MODULE]
                               [--custom-prior-module CUSTOM_PRIOR_MODULE]
                               [--custom-background-module CUSTOM_BACKGROUND_MODULE]

    Script for automated generation of X-PSI model module set. Usage: python
    module_generator.py [-h] @<generate.ini>

    options:
      -h, --help            show this help message and exit
      --generate-config-file
                            Generate the meta configuration file template.
      --telescope TELESCOPE
                            Telescope name, e.g., NICER. Use argument once per
                            telescope name, and no whitespaces.
      --instrument INSTRUMENT
                            Name of an instrument on-board a telescope, e.g., XTI.
                            Can use one or more instrument names per telescope
                            name, and no whitespaces.
      --source SOURCE       The name of the star, e.g., PSR J0740+6620.
      --frequency FREQUENCY
                            The coordinate spin frequency of the star (Hz).
      --model MODEL         A custom model name, e.g., ST-U + NSX-H, otherwise the
                            model name is constructed from other arguments.
      --hot-region-model {ST,CST,EST,PST,CDT,EDT,PDT}
                            The name of the hot-region model, e.g., ST. Maximum of
                            two argument uses.
      --antipodal-reflection-symmetry
                            Are the two hot regions related via antipodal
                            reflection symmetry? E.g., ST-S.
      --break-hot-region-exchange-degeneracy-with BREAK_HOT_REGION_EXCHANGE_DEGENERACY_WITH
                            Hot region parameter name to break hot-region exchange
                            degeneracy with when there are two hot-regions of the
                            same type that are not antipodally reflection-
                            symmetric, e.g., ST+ST (ST-U). An example is e.g.,
                            "super_temperature".
      --is-antiphased IS_ANTIPHASED
                            Specify whether the hot regions are anti-phased w.r.t
                            to Earth. If True, the cell mesh shifts by pi radians
                            about the stellar rotation axis for pulse integration
                            and therefore the hot region at phase zero is aligned
                            with the meridian on which the observer’s antipode
                            lies.
      --prefix PREFIX       Specify the prefixes for hot region parameter naming.
      --hot-atmosphere-model HOT_ATMOSPHERE_MODEL
                            Name of atmosphere model within hot regions, e.g.,
                            blackbody or NSX-H.
      --hot-atmosphere-load
                            Does a numeric atmosphere table need to be loaded from
                            disk for the hot regions?
      --elsewhere-atmosphere-model ELSEWHERE_ATMOSPHERE_MODEL
                            Name of atmosphere model elsewhere, e.g., blackbody or
                            NSX-H.
      --elsewhere-atmosphere-load
                            Does a numeric atmosphere table need to be loaded from
                            disk for elsewhere?
      --attenuation-model ATTENUATION_MODEL
                            Name of interstellar attenuation model, e.g., tbnew.
      --background-model    Include an incident background component?
      --background-shared-instance
                            Do all instruments share the same background model
                            instance?
      --background-shared-class
                            Do all instrument models share a background class?
      --background-parameters [BACKGROUND_PARAMETERS ...]
                            Background model parameter names.
      --print-MPI-rank      Print MPI rank from main module?
      --config-path CONFIG_PATH
                            If main module is imported, use this argument to
                            specify the relative or absolute path to the
                            configuration file.
      --module-directory-path MODULE_DIRECTORY_PATH
                            Absolute path to directory to write module files to.
      --main-module MAIN_MODULE
                            Name of the main module.
      --custom-signal-module CUSTOM_SIGNAL_MODULE
                            Name of the module containing the CustomSignal
                            subclass.
      --custom-instrument-module CUSTOM_INSTRUMENT_MODULE
                            Name of the module containing the CustomInstrument
                            subclass.
      --custom-photosphere-module CUSTOM_PHOTOSPHERE_MODULE
                            Name of the module containing the CustomPhotosphere
                            subclass.
      --custom-interstellar-module CUSTOM_INTERSTELLAR_MODULE
                            Name of the module containing the CustomInterstellar
                            subclass.
      --custom-prior-module CUSTOM_PRIOR_MODULE
                            Name of the module containing the CustomPrior
                            subclass.
      --custom-background-module CUSTOM_BACKGROUND_MODULE
                            Name of the module containing the CustomBackground
                            subclass(es).

Most of the flags displayed above describe the command-line arguments that the user needs to pass to define the kind of
model to generate. The user can choose to have these arguments written into a ``generate.ini`` file for the module
generator to read, instead of passing them individually from the command line.

The ``generate.ini`` file can be created as follows:

.. code-block:: bash

    python module_generator.py --generate-config-file

The corresponding output is:

.. code-block:: bash

    /=============================================\
    | X-PSI: X-ray Pulse Simulation and Inference |
    |---------------------------------------------|
    |                Version: 2.0.0               |
    |---------------------------------------------|
    |      https://xpsi-group.github.io/xpsi      |
    \=============================================/

    Imported GetDist version: 1.4
    Imported nestcheck version: 0.2.1
    Parsing configuration file...
    Configuration file generated.

Let's take a look at the ``generate.ini`` file created:

.. code-block:: bash

    cat generate.ini

It should look like this:

.. code-block:: bash

    ##----------------------------##
    ## telescope instrument flags ##
    ##----------------------------##
    --telescope=
    #--telescope=
    --instrument=
    #--instrument=


    ##---------------------##
    ## target source flags ##
    ##---------------------##
    --source=
    --frequency=


    ##-------------##
    ## model flags ##
    ##-------------##
    #--model=
    --hot-region-model=
    #--hot-region-model=
    #--antipodal-reflection-symmetry
    #--break-hot-region-exchange-degeneracy-with=super_colatitude
    --is-antiphased=
    #--is-antiphased=
    --prefix=
    #--prefix=
    --hot-atmosphere-model=
    #--hot-atmosphere-load
    --elsewhere-atmosphere-model=
    #--elsewhere-atmosphere-load
    --attenuation-model=


    #--background-model
    --background-shared-instance
    --background-shared-class
    #--background-parameters ## enter one name per line below
    #powerlaw_index
    #powerlaw_normalization





    ##---------------------##
    ## miscellaneous flags ##
    ##---------------------##
    --print-MPI-rank


    ##-------------##
    ## write flags ##
    ##-------------##
    --config-path=
    --module-directory-path=
    --main-module=main
    --custom-signal-module=CustomSignal
    --custom-instrument-module=CustomInstrument
    --custom-photosphere-module=CustomPhotosphere
    --custom-interstellar-module=CustomInterstellar
    --custom-prior-module=CustomPrior
    --custom-background-module=CustomBackground%                                                                                                                                                                            (xpsi_py3) dc1408@sp-byods-145-109-127-41 docs % cat generate.ini
    ##----------------------------##
    ## telescope instrument flags ##
    ##----------------------------##
    --telescope=
    #--telescope=
    --instrument=
    #--instrument=


    ##---------------------##
    ## target source flags ##
    ##---------------------##
    --source=
    --frequency=


    ##-------------##
    ## model flags ##
    ##-------------##
    #--model=
    --hot-region-model=
    #--hot-region-model=
    #--antipodal-reflection-symmetry
    #--break-hot-region-exchange-degeneracy-with=super_colatitude
    --prefix=
    #--prefix=
    --hot-atmosphere-model=
    #--hot-atmosphere-load
    --elsewhere-atmosphere-model=
    #--elsewhere-atmosphere-load
    --attenuation-model=


    #--background-model
    --background-shared-instance
    --background-shared-class
    #--background-parameters ## enter one name per line below
    #powerlaw_index
    #powerlaw_normalization





    ##---------------------##
    ## miscellaneous flags ##
    ##---------------------##
    --print-MPI-rank


    ##-------------##
    ## write flags ##
    ##-------------##
    --config-path=
    --module-directory-path=
    --main-module=main
    --custom-signal-module=CustomSignal
    --custom-instrument-module=CustomInstrument
    --custom-photosphere-module=CustomPhotosphere
    --custom-interstellar-module=CustomInterstellar
    --custom-prior-module=CustomPrior
    --custom-background-module=CustomBackground

We can modify the ``generate.ini`` file as per our need by filling up, commenting and/or removing the arguments provided.
An example of a filled out ``generate.ini`` file is present in ``../examples/examples_module_generator`` which creates a
`CST+PDT` hot-region model for PSR J0030+0451 using `NSX-H` atmosphere model.

The ``generate.ini`` file can then be used to create the required `Main` and `Custom` modules as follows:

.. code-block:: bash

    python module_generator.py @generate.ini

The corresponding output below reflects the arguments passed. Note that any empty arguments that aren't commented out
will take in the default value specified.

.. code-block:: bash

    /=============================================\
    | X-PSI: X-ray Pulse Simulation and Inference |
    |---------------------------------------------|
    |                Version: 2.0.0               |
    |---------------------------------------------|
    |      https://xpsi-group.github.io/xpsi      |
    \=============================================/

    Imported GetDist version: 1.4
    Imported nestcheck version: 0.2.1
    Parsing configuration file...
    --telescope=NICER
    --instrument=XTI
    --source=PSR J0030+0451
    --frequency=205
    --model=CST+PDT
    --hot-region-model=CST
    --hot-region-model=PDT
    --is-antiphased=False
    --is-antiphased=False
    --prefix=p
    --prefix=s
    --hot-atmosphere-model=NSX-H
    --hot-atmosphere-load
    --attenuation-model=tbnew
    --print-MPI-rank
    --config-path=./config.ini
    --module-directory-path=./_auto_modules
    --main-module=main
    --custom-signal-module=CustomSignal
    --custom-instrument-module=CustomInstrument
    --custom-photosphere-module=CustomPhotosphere
    --custom-interstellar-module=CustomInterstellar
    --custom-prior-module=CustomPrior
    Configuration file parsed.

Let's take a look at the files generated.

.. code-block:: bash

    ls _auto_modules

In the output below we can see the `Main` and `Custom` files that have been created in the ``_auto_module`` directory as
prompted in ``generate.ini``.

.. code-block:: bash

    CustomInstrument.py
    CustomInterstellar.py
    CustomPhotosphere.py
    CustomPrior.py
    CustomSignal.py
    __init__.py
    main.py

Now that the necessary modules for the model have been generated, we need to pass command-line arguments to specify the
external files required (e.g. for data, instrument, interstellar extinction, atmosphere model, background, etc.),
provide additional details required to read these files, and specify our prior definitions.

The list of arguments we can pass to the modules and their details can be obtained by running the help command for
``main.py`` as follows:

.. code-block:: bash

     python _auto_modules/main.py -h

Here we do not display the output of the help command given the large number of potential arguments that the modules can
accept. Again, in order to avoid passing individual arguments, we can make a ``config.ini`` file for the modules to
read.

The ``config.ini`` file can be created as follows:

.. code-block:: bash

    python _auto_modules/main.py --generate-config-file

The corresponding output is:

.. code-block:: bash

    /=============================================\
    | X-PSI: X-ray Pulse Simulation and Inference |
    |---------------------------------------------|
    |                Version: 2.0.0               |
    |---------------------------------------------|
    |      https://xpsi-group.github.io/xpsi      |
    \=============================================/

    Imported GetDist version: 1.4
    Imported nestcheck version: 0.2.1
    Parsing configuration file...
    Configuration file generated.

You can confirm that the ``config.ini`` has been generated by running the following command:

.. code-block:: bash

    ls config.ini

You can check the content of the empty ``config.ini`` file by running:

.. code-block:: bash

    cat config.ini

Again, we don't display the empty ``config.ini`` file in here given the large number of potential arguments to pass to
the modules. We can modify the ``config.ini`` file as per our need by filling up, commenting and/or removing the
arguments provided. An example of a filled out ``config.ini`` file is present in ``../examples/examples_module_generator``.
The files specified in the examples config file can be found on `Zenodo <https://zenodo.org/record/7113931#.Y90fHi8w35k>`_.

The modules can then be run as follows:

.. code-block:: bash

    python _auto_modules/main.py @config.ini [--multinest] [-emcee]

The additional flags specify the sampler to be used. Note that any empty arguments that aren't commented out
will take in the default value(s) specified.
