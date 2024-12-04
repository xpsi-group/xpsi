.. _citation:

Citation
--------

If X-PSI proves to be a useful tool for your work, please cite the project
as a software acknowledgement, e.g.:

.. code-block:: latex

    X-PSI (\url{https://github.com/xpsi-group/xpsi.git}, \citet{xpsi})

and cite the `Riley et al. (2023) X-PSI Journal of Open Source Software (JOSS) paper <https://doi.org/10.21105/joss.04977>`_

.. code-block:: latex

	@ARTICLE{xpsi,
       	author = {{Riley}, Thomas E. 	and {Choudhury}, Devarshi 
				     	and {Salmi}, Tuomo 
     					and {Vinciguerra}, Serena 
					and {Kini}, Yves 
					and {Dorsman}, Bas 
					and {Watts}, Anna L. 
					and {Huppenkothen}, D. 
					and {Guillot}, Sebastien },
        title = "{X-PSI: A Python package for neutron star X-ray pulse simulation and inference}",
        journal = "{Journal of Open Source Software}",
	publisher = {The Open Journal},
        year = 2023,
	volume = {8}, 
	number = {82}, 
	pages = {4977},
       	doi  = {10.21105/joss.04977},
	url = {https://doi.org/10.21105/joss.04977}
	}


If you wish to cite published applications of the X-PSI package you can find a list of papers
that have used the software via :ref:`Applications`. 

X-PSI has numerous dependencies that users also benefit from.
Citations to these packages may be found in the software acknowledgements of both the X-PSI 
JOSS paper and our :ref:`Applications` papers, and in the
article text where the relevant software is applied. Links to the software and
associated preprints and journal articles may be found in the bibliographies.
Please follow the guidelines demonstrated in e.g. `Riley et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021ApJ...918L..27R/abstract>`_
to ensure the authors of the software receive the recognition they should.
The dependencies often have their own citation instructions that should be
followed first and foremost, but as an example from `Riley et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021ApJ...918L..27R/abstract>`_. 

.. code-block:: latex

    \newcommand{\project}[1]{\textsl{#1}\xspace}

    \software{Python/C~language~\citep{python2007},
              GNU~Scientific~Library~\citep[GSL;][]{Gough:2009},
              NumPy~\citep{Numpy2011},
              Cython~\citep{cython2011},
              SciPy~\citep{Scipy},
              OpenMP~\citep{openmp},
              MPI~\citep{MPI},
              \project{MPI for Python}~\citep{mpi4py},
              Matplotlib~\citep{Hunter:2007,matplotlibv2},
              IPython~\citep{IPython2007},
              Jupyter~\citep{Kluyver:2016aa},
              \MultiNest~\citep{MultiNest_2009},
              \textsc{PyMultiNest}~\citep{PyMultiNest},
              \project{GetDist}~\citep[][\url{https://github.com/cmbant/getdist}]{Lewis19},
              \project{nestcheck}~\citep{higson2018nestcheck,higson2018sampling,higson2019diagnostic},
              \project{fgivenx}~\citep{fgivenx},
              \XPSI~\texttt{v1.0} (\url{https://github.com/ThomasEdwardRiley/xpsi}; \citealt{xpsi}).

    @article{fgivenx,
        doi = {10.21105/joss.00849},
        year  = {2018},
        month = {Aug},
        publisher = {The Open Journal},
        volume = {3},
        number = {28},
        pages = {849},
        author = {Will Handley},
        title = {fgivenx: Functional Posterior Plotter},
        journal = {The Journal of Open Source Software}
    }

    @article{higson2019diagnostic,
        title={nestcheck: diagnostic tests for nested sampling calculations},
        author={Higson, Edward and Handley, Will and Hobson, Mike and Lasenby, Anthony},
        journal={Monthly Notices of the Royal Astronomical Society},
        year={2019},
        volume={483},
        number={2},
        pages={2044--2056},
        doi={10.1093/mnras/sty3090},
        archivePrefix={arXiv},
        arxivId={1804.06406}
    }

    @article{higson2018sampling,
        title={Sampling Errors in Nested Sampling Parameter Estimation},
        author={Higson, Edward and Handley, Will and Hobson, Mike and Lasenby, Anthony},
        year={2018},
        journal={Bayesian Analysis},
        number={3},
        volume={13},
        pages={873--896},
        doi={10.1214/17-BA1075},
    }

    @article{higson2018nestcheck,
        title={nestcheck: error analysis, diagnostic tests and plots for nested sampling calculations},
        author={Higson, Edward},
        year={2018},
        journal={Journal of Open Source Software},
        number={29},
        pages={916},
        volume={3},
        doi={10.21105/joss.00916},
    }

    @ARTICLE{Lewis19,
        author = {{Lewis}, Antony},
        title = "{GetDist: a Python package for analysing Monte Carlo samples}",
        journal = {arXiv e-prints},
        keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Cosmology and Nongalactic Astrophysics, Physics - Data Analysis, Statistics and Probability},
        year = 2019,
        month = oct,
        eid = {arXiv:1910.13970},
        pages = {arXiv:1910.13970},
        archivePrefix = {arXiv},
        eprint = {1910.13970},
        primaryClass = {astro-ph.IM},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2019arXiv191013970L},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

    @ARTICLE{PyMultiNest,
        author = {{Buchner}, J. and {Georgakakis}, A. and {Nandra}, K. and {Hsu}, L. and {Rangel}, C. and {Brightman}, M. and {Merloni}, A. and {Salvato}, M. and {Donley}, J. and {Kocevski}, D.},
        title = "{X-ray spectral modelling of the AGN obscuring region in the CDFS: Bayesian model selection and catalogue}",
        journal = {\aap},
        archivePrefix = "arXiv",
        eprint = {1402.0004},
        primaryClass = "astro-ph.HE",
        keywords = {accretion, accretion disks, methods: data analysis, methods: statistical, galaxies: nuclei, X-rays: galaxies, galaxies: high-redshift},
        year = 2014,
        month = apr,
        volume = 564,
        eid = {A125},
        pages = {A125},
        doi = {10.1051/0004-6361/201322971},
        adsurl = {http://adsabs.harvard.edu/abs/2014A\%26A...564A.125B},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

    @ARTICLE{MultiNest_2009,
        author = {{Feroz}, F. and {Hobson}, M.~P. and {Bridges}, M.},
        title = "{MULTINEST: an efficient and robust Bayesian inference tool for cosmology and particle physics}",
        journal = {\mnras},
        archivePrefix = "arXiv",
        eprint = {0809.3437},
        keywords = {methods: data analysis , methods: statistical},
        year = 2009,
        month = oct,
        volume = 398,
        pages = {1601-1614},
        doi = {10.1111/j.1365-2966.2009.14548.x},
        adsurl = {http://adsabs.harvard.edu/abs/2009MNRAS.398.1601F},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

    @conference{Kluyver:2016aa,
        Author = {Thomas Kluyver and Benjamin Ragan-Kelley and Fernando P{\'e}rez and Brian Granger and Matthias Bussonnier and Jonathan Frederic and Kyle Kelley and Jessica Hamrick and Jason Grout and Sylvain Corlay and Paul Ivanov and Dami{\'a}n Avila and Safia Abdalla and Carol Willing},
        Booktitle = {Positioning and Power in Academic Publishing: Players, Agents and Agendas},
        Editor = {F. Loizides and B. Schmidt},
        Organization = {IOS Press},
        Pages = {87 - 90},
        Title = {Jupyter Notebooks -- a publishing format for reproducible computational workflows},
        Year = {2016}
    }

    @ARTICLE{IPython2007,
        author={F. {Perez} and B. E. {Granger}},
        journal={Computing in Science Engineering},
        title={IPython: A System for Interactive Scientific Computing},
        year={2007},
        volume={9},
        number={3},
        pages={21-29},
        keywords={data visualisation;natural sciences computing;object-oriented languages;object-oriented programming;parallel programming;software libraries;IPython;interactive scientific computing;comprehensive library;data visualization;distributed computation;parallel computation;Scientific computing;Libraries;Data visualization;Spine;Supercomputers;Hardware;Data analysis;Testing;Production;Parallel processing;Python;computer languages;scientific programming;scientific computing},
        doi={10.1109/MCSE.2007.53},
        ISSN={1521-9615},
        month={May}
    }

    @misc{matplotlibv2,
        author       = {Michael Droettboom and
                        Thomas A Caswell and
                        John Hunter and
                        Eric Firing and
                        Jens Hedegaard Nielsen and
                        Antony Lee and
                        Elliott Sales de Andrade and
                        Nelle Varoquaux and
                        David Stansby and
                        Benjamin Root and
                        Phil Elson and
                        Darren Dale and
                        Jae-Joon Lee and
                        Ryan May and
                        Jouni K. Seppänen and
                        Jody Klymak and
                        Damon McDougall and
                        Andrew Straw and
                        Paul Hobson and
                        cgohlke and
                        Tony S Yu and
                        Eric Ma and
                        Adrien F. Vincent and
                        Steven Silvester and
                        Charlie Moad and
                        Jan Katins and
                        Nikita Kniazev and
                        Tim Hoffmann and
                        Federico Ariza and
                        Peter Würtz},
        title        = {matplotlib/matplotlib v2.2.2},
        month        = mar,
        year         = 2018,
        publisher    = {Zenodo},
        doi          = {10.5281/zenodo.1202077},
    }

    @Article{Hunter:2007,
        Author    = {Hunter, J. D.},
        Title     = {Matplotlib: A 2D graphics environment},
        Journal   = {Computing in Science \& Engineering},
        Volume    = {9},
        Number    = {3},
        Pages     = {90--95},
        abstract  = {Matplotlib is a 2D graphics package used for Python for
        application development, interactive scripting, and publication-quality
        image generation across user interfaces and operating systems.},
        publisher = {IEEE COMPUTER SOC},
        doi       = {10.1109/MCSE.2007.55},
        year      = 2007
    }

    @article{mpi4py,
        author = {Lisandro Dalc\'{i}n and Rodrigo Paz and Mario Storti and Jorge D'El\'{i}a},
        title = {MPI for Python: Performance improvements and MPI-2 extensions},
        journal = {Journal of Parallel and Distributed Computing},
        volume = {68},
        number = {5},
        pages = {655-662},
        year = {2008},
        issn = {0743-7315},
        doi = {10.1016/j.jpdc.2007.09.005},
        keywords = {Message passing, MPI, High-level languages, Parallel Python}
    }

    @techreport{MPI,
        author = {Forum, Message P},
        title = {MPI: A Message-Passing Interface Standard},
        year = {1994},
        url = {https://www.mpi-forum.org/docs/mpi-1.0/mpi-10.ps},
        publisher = {University of Tennessee},
        address = {Knoxville, TN, USA}
    }

    @article{openmp,
        Author = {Dagum, Leonardo and Menon, Ramesh},
        Date-Added = {2014-07-24 11:13:01 +0000},
        Date-Modified = {2014-07-24 11:13:01 +0000},
        Journal = {Computational Science \& Engineering, IEEE},
        Number = {1},
        Pages = {46--55},
        Publisher = {IEEE},
        Title = {OpenMP: an industry standard API for shared-memory programming},
        Volume = {5},
        Year = {1998}
    }

    @misc{Scipy,
        author = {Eric Jones and Travis Oliphant and Pearu Peterson and others},
        title = {{SciPy}: Open source scientific tools for {Python}},
        year = {2001--},
        url = "http://www.scipy.org/",
        note = {[Online; accessed 21.06.2019]}
    }

    @ARTICLE{cython2011,
        author={S. {Behnel} and R. {Bradshaw} and C. {Citro} and L. {Dalcin} and D. S. {Seljebotn} and K. {Smith}},
        journal={Computing in Science Engineering},
        title={Cython: The Best of Both Worlds},
        year={2011},
        volume={13},
        number={2},
        pages={31-39},
        keywords={C language;numerical analysis;Python language extension;Fortran code;numerical loops;Cython language;programming language;Sparse matrices;Runtime;Syntactics;Computer programs;Programming;Python;Cython;numerics;scientific computing},
        doi={10.1109/MCSE.2010.118},
        ISSN={1521-9615},
        month={March}
    }

    @ARTICLE{Numpy2011,
        author={S. {van der Walt} and S. C. {Colbert} and G. {Varoquaux}},
        journal={Computing in Science Engineering},
        title={The NumPy Array: A Structure for Efficient Numerical Computation},
        year={2011},
        volume={13},
        number={2},
        pages={22-30},
        keywords={data structures;high level languages;mathematics computing;numerical analysis;numerical computation;numpy array;numerical data;high level language;Python programming language;Arrays;Numerical analysis;Performance evaluation;Computational efficiency;Finite element methods;Vector quantization;Resource management;Python;NumPy;scientific programming;numerical computations;programming libraries}, 
        doi={10.1109/MCSE.2011.37},
        ISSN={1521-9615},
        month={March}
    }

    @book{Gough:2009,
        author = {Gough, Brian},
        title = {GNU Scientific Library Reference Manual - Third Edition},
        year = {2009},
        isbn = {0954612078, 9780954612078},
        edition = {3rd},
        publisher = {Network Theory Ltd.}
    }

    @ARTICLE{Python2007,
        author={T. E. {Oliphant}},
        journal={Computing in Science Engineering},
        title={Python for Scientific Computing},
        year={2007},
        volume={9},
        number={3},
        pages={10-20},
        keywords={high level languages;Python;scientific computing;steering language;scientific codes;high-level language;Scientific computing;High level languages;Libraries;Writing;Application software;Embedded software;Software standards;Standards development;Internet;Prototypes;Python;computer languages;scientific programming;scientific computing}, 
        doi={10.1109/MCSE.2007.58},
        ISSN={1521-9615},
        month={May}
    }

    @ARTICLE{2021JOSS....6.3001B,
        author = {{Buchner}, Johannes},
        title = "{UltraNest - a robust, general purpose Bayesian inference engine}",
        journal = {The Journal of Open Source Software},
        keywords = {C, Monte Carlo, Python, Nested Sampling, C++, Bayesian inference, Fortran, Bayes factors, Statistics - Computation, Astrophysics - Instrumentation and Methods for Astrophysics},
        year = 2021,
        month = apr,
        volume = {6},
        number = {60},
        eid = {3001},
        pages = {3001},
        doi = {10.21105/joss.03001},
        archivePrefix = {arXiv},
        eprint = {2101.09604},
        primaryClass = {stat.CO},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2021JOSS....6.3001B},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
        
    }


