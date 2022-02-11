.. _example_script:

Example script and modules
==========================

The following model script is an example from :ref:`R19`. One would run the
program from the command line on a desktop via:

.. code-block:: bash

    mpiexec -n 4 python [-m mpi4py] main.py

where ``main.py`` is of course the script name, and where for a machine with
four physical cores, we spawn as many MPI processes.

.. note:: Make sure to set OpenMP environment variables appropriately
          (see, e.g., :ref:`surfsystems`) to avoid hybrid parallelisation by
          libraries linked to NumPy and GSL, for instance.

The optional ``-m`` argument is to run the module :mod:`mpi4py` as main, which
enables elegant termination of all MPI processes if one process encounters
an unhandled exception, for instance; without out this flag, MPI processes
may hang if an unhandled exception is raised in a subset of processes. An
example of this is an exception thrown from a :class:`~.Prior` subclass during
set-up for sampling: only the root (rank zero) process encounters the
exception, whilst all other wait for it and will hang indefinitely. If one
is not monitoring the output during the start of the program to check that
sensible things are happening, this can be missed, which is especially
problematic on a HPC system on which one is budgeted time and thus there is a
serious risk of resource wastage from the perspective of all users of the
system.

The script below is an example like that found on the `Modeling`
page, but without being interleaved with verbose explanations. The model
implemented here is more involved than shown on the `Modeling`
page, and than those defined for the example parameter
estimation problems that can be found in the ``xpsi/examples`` directory of
the repository.

.. note::

    The script and modules below are being updated to use
    the current API (``v0.5.0``).

Main
^^^^

.. code-block:: python

    """ main.py """

    from __future__ import print_function, division

    import numpy as np
    import math

    import xpsi

    from xpsi.global_imports import gravradius

    from CustomInstrument import CustomInstrument
    from CustomInterstellar import CustomInterstellar
    from CustomSignal import CustomSignal
    from CustomPrior import CustomPrior
    from CustomPhotosphere import CustomPhotosphere

    path = 'data/NICER_J0030_PaulRay_fixed_evt_25to299__preprocessed.txt'
    obs_settings = dict(counts=np.loadtxt(path, dtype=np.double),
                        channels=np.arange(25, 300),
                        phases=np.linspace(0.0, 1.0, 33),
                        first=0, last=274,
                        exposure_time=1936864.0)

    data = xpsi.Data(**obs_settings)

    bounds = dict(alpha = (0.5,1.5),
                  beta = (0.0,1.0),
                  gamma = (0.5,1.5))

    NICER = CustomInstrument.from_SWG(bounds = bounds,
                                      values = {},
                                      ARF = 'model_data/ni_xrcall_onaxis_v1.02_arf.txt',
                                      RMF = 'model_data/nicer_upd_d49_matrix.txt',
                                      ratio = 'model_data/crab_ratio_SA80_d49.txt',
                                      max_input=700,
                                      min_input=0,
                                      channel_edges = 'model_data/nicer_upd_energy_bounds.txt')

    interstellar = CustomInterstellar.from_SWG('model_data/interstellar_phot_frac.txt',
                                           bounds = dict(column_density = (0.0,5.0)))

    signal = CustomSignal(data = data,
                          instrument = NICER,
                          interstellar = interstellar,
                          cache = True,
                          workspace_intervals = 1000,
                          epsrel = 1.0e-8,
                          epsilon = 1.0e-3,
                          sigmas = 10.0)

    bounds = dict(mass = (1.0, 3.0),
                  radius = (3.0*gravradius(1.0), 16.0),
                  distance = (0.05, 2.0),
                  cos_inclination = (0.0, math.cos(0.001)))

    spacetime = xpsi.Spacetime(bounds, dict(frequency = 1.0/(4.87e-3)))

    bounds = dict(super_colatitude = (0.001, math.pi - 0.001),
                  super_radius = (0.001, math.pi/2.0 - 0.001),
                  phase_shift = (None, None),
                  super_temperature = (5.1, 6.8))

    primary = xpsi.HotRegion(bounds=bounds,
                                values={},
                                symmetry=True,
                                omit=False,
                                cede=False,
                                concentric=False,
                                sqrt_num_cells=24,
                                min_sqrt_num_cells=10,
                                max_sqrt_num_cells=64,
                                do_fast=False,
                                num_leaves=80,
                                num_rays=200,
                                is_secondary=False,
                                prefix='p')

    # we transform to these geometric parameters, so see CustomPrior instead
    # for inverse sampling setup
    bounds = dict(super_colatitude = (None, None), # see CustomPrior
                    super_radius = (None, None), # see CustomPrior
                    phase_shift = (-0.5, 0.5),
                    super_temperature = (5.1, 6.8),
                    omit_colatitude = (0.0, math.pi),
                    omit_radius = (None, None), # see CustomPrior
                    omit_azimuth = (None, None)) # see CustomPrior

    # overlap of an omission region and
    # and a radiating super region
    secondary = xpsi.HotRegion(bounds=bounds,
                                values={},
                                symmetry=True,
                                omit=True,
                                cede=False,
                                concentric=False,
                                sqrt_num_cells=24,
                                min_sqrt_num_cells=10,
                                max_sqrt_num_cells=64,
                                num_leaves=80,
                                num_rays=200,
                                do_fast=False,
                                is_secondary=True,
                                prefix='s')

    from xpsi import HotRegions

    hot = HotRegions((primary, secondary))

    photosphere = CustomPhotosphere(hot = hot, elsewhere = None,
                                    values=dict(mode_frequency = spacetime['frequency']))

    photosphere.hot_atmosphere = 'model_data/nsx_H_v171019.out'

    star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

    likelihood = xpsi.Likelihood(star = star, signals = signal,
                                 num_energies = 128,
                                 threads = 1,
                                 externally_updated = True)

    prior = CustomPrior()

    likelihood.prior = prior

    p = [1.4033703360094012,
         13.378462458584202,
         0.32897884439908337,
         math.cos(1.004349731136371),
         0.4542555093514883,
         2.1937752730930784,
         0.07916088420116879,
         6.106556223820221,
         0.4768294130316574,
         2.7162985247930496,
         0.32234225478780626,
         6.1173049179880445,
         2.7463301464251777,
         0.2844169651751102,
         -0.048326090505605386,
         1.0335682718716097,
         0.02227107198360202,
         0.8748566319738948,
         0.4604998629950954]

    # source code changes since model was applied, so let's be a
    # bit lenient when checking the likelihood function
    likelihood.check(None, [-36316.354394388654], 1.0e-4,
                     physical_points=[p])

    wrapped_params = [0] * len(likelihood)
    wrapped_params[likelihood.index('s__phase_shift')] = 1

    runtime_params = {'resume': False,
                      'importance_nested_sampling': False,
                      'multimodal': False,
                      'n_clustering_params': None,
                      'outputfiles_basename': './run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1',
                      'n_iter_before_update': 100,
                      'n_live_points': 1000,
                      'sampling_efficiency': 0.3,
                      'const_efficiency_mode': False,
                      'wrapped_params': wrapped_params,
                      'evidence_tolerance': 0.1,
                      'max_iter': -1,
                      'verbose': True}

    xpsi.Sample.nested(likelihood, prior, check_kwargs, **runtime_params)


We proceed to show the custom modules required for the model.

.. todo::

    Write more extensive inline comments for clarity, and clean up where
    applicable.

Photosphere
^^^^^^^^^^^

.. code-block:: python

    """ CustomPhotosphere.py """

    import numpy as np
    import math

    import xpsi

    class CustomPhotosphere(xpsi.Photosphere):
        """ A photosphere extension to preload the numerical atmosphere NSX. """

        @xpsi.Photosphere.hot_atmosphere.setter
        def hot_atmosphere(self, path):
            NSX = np.loadtxt(path, dtype=np.double)
            logT = np.zeros(35)
            logg = np.zeros(11)
            mu = np.zeros(67)
            logE = np.zeros(166)

            reorder_buf = np.zeros((35,11,67,166))

            index = 0
            for i in range(reorder_buf.shape[0]):
                for j in range(reorder_buf.shape[1]):
                    for k in range(reorder_buf.shape[3]):
                       for l in range(reorder_buf.shape[2]):
                            logT[i] = NSX[index,3]
                            logg[j] = NSX[index,4]
                            logE[k] = NSX[index,0]
                            mu[reorder_buf.shape[2] - l - 1] = NSX[index,1]
                            reorder_buf[i,j,reorder_buf.shape[2] - l - 1,k] = 10.0**(NSX[index,2])
                            index += 1

            buf = np.zeros(np.prod(reorder_buf.shape))

            bufdex = 0
            for i in range(reorder_buf.shape[0]):
                for j in range(reorder_buf.shape[1]):
                    for k in range(reorder_buf.shape[2]):
                       for l in range(reorder_buf.shape[3]):
                            buf[bufdex] = reorder_buf[i,j,k,l]; bufdex += 1

            self._hot_atmosphere = (logT, logg, mu, logE, buf)

Instrument
^^^^^^^^^^

.. code-block:: python

    """ CustomInstrument.py """

    from __future__ import print_function, division

    import numpy as np
    import math

    import xpsi
    from xpsi import Parameter

    class CustomInstrument(xpsi.Instrument):
        """ Methods and attributes specific to the NICER instrument.

        Currently tailored to the NICER light-curve SWG model specification.

        """
        def __init__(self, *args):
            """ Set channel edges attribute. """
            super(CustomInstrument, self).__init__(*args)

        def construct_matrix(self):
            """ Implement response matrix parameterisation. """
            matrix = self['alpha']*self.matrix
            matrix[matrix < 0.0] = 0.0
            return matrix

        def __call__(self, signal, *args):
            """ Overwrite. """

            matrix = self.construct_matrix()

            self._cached_signal = np.dot(matrix, signal)

            return self._cached_signal

        @classmethod
        def from_SWG(cls,
                     bounds, values,
                     ARF, RMF,
                     max_input, min_input=0,
                     channel_edges=None):
            """ Constructor which converts files into :class:`numpy.ndarray`s.

            :param str ARF: Path to ARF which is compatible with
                                    :func:`numpy.loadtxt`.

            :param str RMF: Path to RMF which is compatible with
                                    :func:`numpy.loadtxt`.

            :param str channel_edges: Optional path to edges which is compatible with
                                    :func:`numpy.loadtxt`.

            """
            ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
            RMF = np.loadtxt(RMF, dtype=np.double, skiprows=3, usecols=-1)

            if channel_edges:
                channel_edges = np.loadtxt(channel_edges, dtype=np.double, skiprows=3)

            matrix = np.zeros((1501,3980))

            for i in range(3980):
                matrix[:,i] = RMF[i*1501:(i+1)*1501]

            if min_input != 0:
                min_input = int(min_input)

            max_input = int(max_input)

            edges = np.zeros(ARF[min_input:max_input,3].shape[0]+1, dtype=np.double)

            edges[0] = ARF[min_input,1]; edges[1:] = ARF[min_input:max_input,2]

            RSP = np.ascontiguousarray(np.zeros(matrix[25:300,min_input:max_input].shape), dtype=np.double)

            for i in range(RSP.shape[0]):
                RSP[i,:] = matrix[i+25,min_input:max_input] * ARF[min_input:max_input,3] * 49.0/52.0

            channels = np.arange(25, 300)

            alpha = Parameter('alpha',
                              strict_bounds = (0.0,2.0),
                              bounds = bounds.get('alpha', None),
                              doc = 'alpha',
                              symbol = r'$\alpha$',
                              value = values.get('alpha', None))


            return cls(RSP, edges, channels, channel_edges[25:301, -2], alpha)

Interstellar
^^^^^^^^^^^^

.. code-block:: python

    """ CustomInterstellar.py """

    from __future__ import print_function, division

    import numpy as np
    import math

    import xpsi
    from xpsi import Parameter

    from scipy.interpolate import Akima1DInterpolator

    class CustomInterstellar(xpsi.Interstellar):
        """ Apply interstellar attenuation. """

        def __init__(self, energies, attenuation, bounds, values = {}):

            assert len(energies) == len(attenuation), 'Array length mismatch.'

            self._lkp_energies = energies # for lookup
            self._lkp_attenuation = attenuation # for lookup

            N_H = Parameter('column_density',
                            strict_bounds = (0.0,10.0),
                            bounds = bounds.get('column_density', None),
                            doc = 'Units of 10^20 cm^-2.',
                            symbol = r'$N_{\rm H}$',
                            value = values.get('column_density', None))

            super(CustomInterstellar, self).__init__(N_H)

        def attenuation(self, energies):
            """ Interpolate the attenuation coefficients.

            Useful for post-processing.

            """
            return self._interpolate(energies)**(self['column_density']/0.4)

        def _interpolate(self, energies):
            """ Helper. """
            try:
                self._interpolator
            except AttributeError:
                self._interpolator = Akima1DInterpolator(self._lkp_energies,
                                                         self._lkp_attenuation)
                self._interpolator.extrapolate = True

            return self._interpolator(energies)

        @classmethod
        def from_SWG(cls, path, **kwargs):
            """ Load attenuation file from the NICER SWG. """

            temp = np.loadtxt(path, dtype=np.double)

            energies = temp[0:351,0]

            attenuation = temp[0:351,2]

            return cls(energies, attenuation, **kwargs)

Signal
^^^^^^

.. code-block:: python

    """ CustomSignal.py """

    from __future__ import print_function, division

    import numpy as np
    import math

    import xpsi

    from xpsi.likelihoods.default_background_marginalisation import eval_marginal_likelihood
    from xpsi.likelihoods.default_background_marginalisation import precomputation

    class CustomSignal(xpsi.Signal):
        """ A custom calculation of the logarithm of the likelihood.

        We extend the :class:`xpsi.Signal.Signal` class to make it callable.

        We overwrite the body of the __call__ method. The docstring for the
        abstract method is copied.

        """

        def __init__(self, workspace_intervals = 1000, epsabs = 0, epsrel = 1.0e-8,
                     epsilon = 1.0e-3, sigmas = 10.0, support = None, *args, **kwargs):
            """ Perform precomputation. """

            super(CustomSignal, self).__init__(*args, **kwargs)

            try:
                self._precomp = precomputation(self._data.counts.astype(np.int32))
            except AttributeError:
                print('No data... can synthesise data but cannot evaluate a '
                      'likelihood function.')
            else:
                self._workspace_intervals = workspace_intervals
                self._epsabs = epsabs
                self._epsrel = epsrel
                self._epsilon = epsilon
                self._sigmas = sigmas

                if support is not None:
                    self._support = support
                else:
                    self._support = -1.0 * np.ones((self._data.counts.shape[0],2))
                    self._support[:,0] = 0.0

        @property
        def support(self):
            return self._support

        @support.setter
        def support(self, obj):
            self._support = obj

        def __call__(self, phase_shifts, *args, **kwargs):
            self.shifts = np.array(phase_shifts)

            self.loglikelihood, self.expected_counts, self.background_signal, self.background_given_support = \
                    eval_marginal_likelihood(self._data.exposure_time,
                                              self._data.phases,
                                              self._data.counts,
                                              self._signals,
                                              self._phases,
                                              self._shifts,
                                              self._precomp,
                                              self._support,
                                              self._workspace_intervals,
                                              self._epsabs,
                                              self._epsrel,
                                              self._epsilon,
                                              self._sigmas,
                                              kwargs.get('llzero'))

Prior
^^^^^

.. code-block:: python

    """ CustomPrior.py """

    from __future__ import print_function, division

    import numpy as np
    import math
    from scipy.stats import truncnorm

    import xpsi
    from xpsi.global_imports import _G, _csq, _km, _M_s, _2pi
    from xpsi.global_imports import gravradius, inv_gravradius

    from xpsi.cellmesh.mesh_tools import eval_cedeCentreCoords

    from scipy.interpolate import Akima1DInterpolator

    class CustomPrior(xpsi.Prior):
        """ A custom (joint) prior distribution.

        Source: PSR J0030+0451
        Model variant: ST+PST
            Two single-temperature hot regions with unshared parameters
            and different complexity levels.

        Parameter vector: (print the likelihood object)

        * p[0] = (rotationally deformed) gravitational mass (solar masses)
        * p[1] = coordinate equatorial radius (km)
        * p[2] = distance (kpc)
        * p[3] = cos(inclination of Earth to rotational axis)
        * p[4] = primary cap phase shift (cycles); (alias for initial azimuth, periodic)
        * p[5] = primary centre colatitude (radians)
        * p[6] = primary angular radius (radians)
        * p[7] = primary log10(comoving NSX FIH effective temperature [K])
        * p[8] = secondary cap phase shift (cycles)
        * p[9] = secondary centre colatitude (radians)
        * p[10] = secondary angular radius (radians)
        * p[11] = secondary omit colatitude (radians)
        * p[12] = secondary omit angular radius (radians)
        * p[13] = secondary omit azimuth (radians); periodic
        * p[14] = secondary log10(comoving NSX FIH effective temperature [K])
        * p[15] = hydrogen column density (10^20 cm^-2)
        * p[16] = instrument parameter alpha
        * p[17] = instrument parameter beta
        * p[18] = instrument parameter gamma

        """

        __derived_names__ = ['compactness',
                             's__annulus_width',
                             's__transformed_phase',
                             's__f',
                             's__xi',
                             's__super_offset_fraction',
                             's__super_offset_azi']

        a_f = 0.0
        b_f = 2.0
        a_xi = 0.001
        b_xi = math.pi/2.0 - a_xi

        vals = np.linspace(0.0, b_xi, 1000)

        def __init__(self):
            """ Construct mapping from unit interval. """

            self.interpolator = Akima1DInterpolator(self._vector_super_radius_mass(self.vals), self.vals)
            self.interpolator.extrapolate = True

        def __call__(self, p = None):
            """ Evaluate distribution at ``p``.

            :param list p: Model parameter values.

            :return: Logarithm of the distribution evaluated at ``p``.

            """
            temp = super(CustomPrior, self).__call__(p)
            if not np.isfinite(temp):
                return temp

            # based on contemporary EOS theory
            if not self.parameters['radius'] <= 16.0:
                return -np.inf

            ref = self.parameters.star.spacetime # shortcut

            # polar radius at photon sphere for ~static star (static ambient spacetime)
            R_p = 1.0 + ref.epsilon * (-0.788 + 1.030 * ref.zeta)
            if R_p < 1.5 / ref.R_r_s:
                return -np.inf

            # limit polar radius to try to exclude deflections >= \pi radians
            # due to oblateness this does not quite eliminate all configurations
            # with deflections >= \pi radians
            #if R_p < 1.76 / ref.R_r_s:
            #    return -np.inf

            mu = math.sqrt(-1.0 / (3.0 * ref.epsilon * (-0.788 + 1.030 * ref.zeta)))

            # 2-surface cross-section have a single maximum in |z|
            # i.e., an elliptical surface; minor effect on support, if any,
            # for high spin frequenies
            if mu < 1.0:
                return -np.inf

            ref = self.parameters # redefine shortcut

            phi = (0.5 + ref['s__phase_shift']) * _2pi
            phi -= ref['s__omit_azimuth']
            phi = ref['p__phase_shift'] * _2pi - phi

            ang_sep = xpsi.HotRegion.psi(ref['s__super_colatitude'],
                                         phi,
                                         ref['p__super_colatitude'])

            # hot regions cannot overlap
            if ang_sep < ref['p__super_radius'] + ref['s__super_radius']:
                return -np.inf

            return 0.0

        def _I(self, x):
            return x * np.log(self.b_xi/self.a_xi)

        def _II(self, x):
            return 2.0*(x - self.a_xi) - x*np.log(x/self.b_xi)

        def _scalar_super_radius_mass(self, x):
            if x >= self.a_xi:
                mass = self._II(x)
            else:
                mass = self._I(x)

            return mass

        def _vector_super_radius_mass(self, x):
            masses = np.zeros(len(x))

            for i, _ in enumerate(x):
                masses[i] = self._scalar_super_radius_mass(_)

            masses /= (self.b_f - self.a_f)
            masses /= (self.b_xi - self.a_xi)

            return masses

        def _inverse_sample_cede_radius(self, x, psi):
            if psi < self.a_xi:
                return self.a_xi*np.exp(x * np.log(self.b_xi/self.a_xi))
            elif psi >= self.a_xi and x <= 1.0/(1.0 + np.log(self.b_xi/psi)):
                return x*psi*(1.0 + np.log(self.b_xi/psi))
            else:
                return psi*np.exp(x*(1.0 + np.log(self.b_xi/psi)) - 1.0)

        def inverse_sample(self, hypercube = None):
            """ Draw sample uniformly from the distribution via inverse sampling.

            :param hypercube: A pseudorandom point in an n-dimensional hypercube.

            :return: A parameter ``list``.

            """
            to_cache = self.parameters.vector

            if hypercube is None:
                hypercube = np.random.rand(len(self))

            _ = super(CustomPrior, self).inverse_sample(hypercube)

            ref = self.parameters # redefine shortcut

            # draw from flat prior in inclination
            idx = ref.index('cos_inclination')
            a, b = ref.get_param('cos_inclination').bounds
            a = math.acos(a); b = math.acos(b)
            ref['cos_inclination'] = math.cos(b + (a - b) * hypercube[idx])

            idx = ref.index('distance')
            ref['distance'] = truncnorm.ppf(hypercube[idx],
                                            -10.0, 10.0,
                                            loc=0.325, scale=0.009)

            idx = ref.index('p__phase_shift')
            phase = 0.35 + 0.2 * hypercube[idx]
            if phase > 0.5:
                ref['p__phase_shift'] = phase - 1.0
            else:
                ref['p__phase_shift'] = phase

            idx = ref.index('s__phase_shift')
            phase = -0.25 + hypercube[idx]
            if phase > 0.5:
                ref['s__phase_shift'] = phase - 1.0
            else:
                ref['s__phase_shift'] = phase

            idx = ref.index('s__omit_radius')
            ref['s__omit_radius'] = float(self.interpolator(hypercube[idx]))

            idx = ref.index('s__super_radius')
            ref['s__super_radius'] = self._inverse_sample_cede_radius(hypercube[idx],
                                                                      ref['s__omit_radius'])

            idx = ref.index('s__super_colatitude')
            if ref['s__omit_radius'] <= ref['s__super_radius']:
                # temp var
                t = hypercube[idx] * (ref['s__super_radius'] + ref['s__omit_radius'])
            else:
                # temp var
                t = ref['s__omit_radius'] - ref['s__super_radius']
                t += 2.0 * hypercube[idx] * ref['s__super_radius']

            idx = ref.index('s__omit_azimuth')
            # temp var
            u = hypercube[idx] * _2pi

            # function from mesh tools module
            # in this case the ceding region is the "super" region, which
            # cedes to the omission region
            ref['s__super_colatitude'], ref['s__omit_azimuth'] = \
                    eval_cedeCentreCoords(ref['s__omit_colatitude'], t, u)

            ref['s__omit_azimuth'] *= -1.0

            idx = ref.index('alpha')
            ref['alpha'] = truncnorm.ppf(hypercube[idx],
                                         -5.0, 5.0,
                                         loc=1.0, scale=0.1)

            idx = ref.index('gamma')
            ref['gamma'] = truncnorm.ppf(hypercube[idx],
                                         -5.0, 5.0,
                                         loc=1.0, scale=0.1)

            # restore proper cache
            for parameter, cache in zip(self.parameters, to_cache):
                parameter.cached = cache

            return self.parameters.vector # only free parameter values returned

        def transform(self, p, old_API = False):
            """ A transformation for post-processing.

            Note that if you want to use dictionary-like access to values,
            you could make a dictionary, e.g.:

            .. code-block:: python

                ref = dict(zip(self.parameters.names, p))

            and use the ``__getitem__`` functionality of ``ref`` instead of
            numeric indexing.

            """

            p = list(p) # copy

            if old_API:
                idx = self.parameters.index('cos_inclination')
                p[idx] = math.cos(p[idx])

            # used ordered names and values
            ref = dict(zip(self.parameters.names, p))

            # compactness ratio M/R_eq
            p += [gravradius(ref['mass']) / ref['radius']]

            p += [ref['s__super_radius'] - ref['s__omit_radius']]

            if ref['s__phase_shift'] > 0.0:
                p += [ref['s__phase_shift'] - 1.0]
            else:
                p += [ref['s__phase_shift']]

            temp = eval_cedeCentreCoords(-1.0*ref['s__omit_colatitude'],
                                         ref['s__super_colatitude'],
                                         -1.0*ref['s__omit_azimuth'])

            azi = temp[1]

            if azi < 0.0:
                azi += 2.0*math.pi

            p += [ref['s__omit_radius']/ref['s__super_radius'] \
                  if ref['s__omit_radius'] <= ref['s__super_radius'] \
                  else 2.0 - ref['s__super_radius']/ref['s__omit_radius']] # f

            p += [ref['s__super_radius'] if ref['s__omit_radius'] \
                  <= ref['s__super_radius'] else ref['s__omit_radius']] # xi

            p += [temp[0]/(ref['s__super_radius'] + ref['s__omit_radius']) \
                  if ref['s__omit_radius'] <= ref['s__super_radius'] \
                  else (temp[0] - ref['s__omit_radius'] + ref['s__super_radius'])/(2.0*ref['s__super_radius'])] # kappa

            p += [azi/math.pi]

            return p
