.. _example_script:

Example script and modules
==========================

The following model script is an example from :ref:`R19`. One would run the
program from the command line on a desktop via:

.. code-block:: bash

    mpiexec -n 4 python main.py

where ``main.py`` is of course the script name, and where for a machine with
four physical cores, we spawn as many MPI processes.

.. note:: Make sure to set OpenMP environment variables appropriately
          (see, e.g., :ref:`surfsystems`) to avoid hybrid parallelisation by
          libraries linked to NumPy and GSL, for instance.

The script below is an example like that found on the `Model construction`
page, but without being interleaved with verbose explanations. The model
implemented here is more involved than shown on the `Model construction`
page, and than those defined for the example parameter
estimation problems that can be found in the ``xpsi/examples`` directory of
the repository.

.. note::

    The script has been updated to use updated object names in the API
    since ``v0.1``.

Main
^^^^

.. code-block:: python

    """ main.py """

    from __future__ import print_function, division

    import numpy as np
    import math

    import xpsi

    print('Rank reporting: %d' % xpsi._rank)

    from CustomData import CustomData
    from CustomInstrument import CustomInstrument
    from CustomInterstellar import CustomInterstellar
    from CustomPulse import CustomPulse
    from CustomSpacetime import CustomSpacetime
    from CustomPrior import CustomPrior
    from CustomPhotosphere import CustomPhotosphere

    data = CustomData.from_SWG('data/NICER_J0030_PaulRay_fixed_evt_25to299__preprocessed.txt', 1936864.0)

    NICER = CustomInstrument.from_SWG(num_params=3,
                        bounds=[(0.5,1.5),(0.0,1.0),(0.5,1.5)],
                        ARF = 'model_data/ni_xrcall_onaxis_v1.02_arf.txt',
                        RMF = 'model_data/nicer_upd_d49_matrix.txt',
                        ratio = 'model_data/crab_ratio_SA80_d49.txt',
                        max_input=700,
                        min_input=0,
                        chan_edges = 'model_data/nicer_upd_energy_bounds.txt')

    interstellar = CustomInterstellar.from_SWG('model_data/interstellar_phot_frac.txt',
                                               num_params = 1,
                                               bounds = [(0.0, 5.0)])

    pulse = CustomPulse(tag = 'all',
                        num_params = 2,
                        bounds = [(0.35, 0.55), (-0.25,0.75)],
                        data = data,
                        instrument = NICER,
                        interstellar = interstellar,
                        energies_per_interval = 0.25,
                        fast_rel_energies_per_interval = 0.5,
                        default_energy_spacing = 'logspace',
                        adaptive_energies = False,
                        adapt_exponent = None,
                        store = False,
                        workspace_intervals = 1000,
                        epsrel = 1.0e-8,
                        epsilon = 1.0e-3,
                        sigmas = 10.0)

    from xpsi.global_imports import _c, _G, _M_s, _dpr, gravradius

    bounds = [(0.235, 0.415),
              (1.0, 3.0),
              (3.0 * gravradius(1.0), 16.0),
              (0.001, math.pi/2.0)]

    spacetime = CustomSpacetime(num_params = 4, bounds = bounds, S = 1.0/(4.87e-3))

    bounds = [(0.001, math.pi - 0.001),
              (0.001, math.pi/2.0 - 0.001),
              (5.1, 6.8)]

    primary = xpsi.HotRegion(num_params=3, bounds=bounds,
                                symmetry=True,
                                hole=False,
                                cede=False,
                                concentric=False,
                                sqrt_num_cells=24,
                                min_sqrt_num_cells=10,
                                max_sqrt_num_cells=64,
                                do_fast=False,
                                fast_sqrt_num_cells=8,
                                fast_min_sqrt_num_cells=8,
                                fast_max_sqrt_num_cells=16,
                                fast_num_leaves=32,
                                fast_num_rays=100,
                                num_leaves=80,
                                num_rays=200)

    bounds = [(0.001, math.pi - 0.001),
              (0.001, math.pi/2.0 - 0.001),
              (0.001, math.pi - 0.001),
              (0.0, 2.0),
              (0.0, 2.0*math.pi),
              (5.1, 6.8)]

    secondary = xpsi.HotRegion(num_params=6, bounds=bounds,
                                  symmetry=True,
                                  hole=True,
                                  cede=False,
                                  concentric=False,
                                  sqrt_num_cells=24,
                                  min_sqrt_num_cells=10,
                                  max_sqrt_num_cells=64,
                                  do_fast=False,
                                  fast_sqrt_num_cells=8,
                                  fast_min_sqrt_num_cells=8,
                                  fast_max_sqrt_num_cells=16,
                                  fast_num_leaves=32,
                                  fast_num_rays=100,
                                  num_leaves=80,
                                  num_rays=200,
                                  is_secondary=True)

    from xpsi import TwoHotRegions

    hot = TwoHotRegions((primary, secondary))

    photosphere = CustomPhotosphere(num_params = 0, bounds = [],
                                    tag = 'all', hot = hot, elsewhere = None)

    photosphere.hot_atmosphere = 'model_data/nsx_H_v171019.out'

    star = xpsi.Star(spacetime = spacetime, photospheres = photosphere)

    likelihood = xpsi.Likelihood(star = star, pulses = pulse, threads=1)

    prior = CustomPrior(bounds=likelihood.bounds, spacetime=spacetime)

    likelihood.prior = prior

    import time

    p = [0.328978844399083370E+00,
            0.140337033600940120E+01,
            0.133784624585842025E+02,
            0.100434973113637094E+01,
            0.219377527309307840E+01,
            0.791608842011687908E-01,
            0.610655622382022134E+01,
            0.271629852479304956E+01,
            0.322342254787806259E+00,
            0.274633014642517770E+01,
            0.284416965175110226E+00,
            -0.483260905056053860E-01,
            0.611730491798804454E+01,
            0.460499862995095377E+00,
            0.103356827187160971E+01,
            0.222710719836020192E-01,
            0.874856631973894849E+00,
            0.454255509351488285E+00,
            0.476829413031657379E+00]

    t = time.time()
    ll = likelihood(p) # check ll = -36316.354394388654
    print('p: ', ll, time.time() - t)

    runtime_params = {'resume': False,
                      'importance_nested_sampling': False,
                      'multimodal': False,
                      'n_clustering_params': None,
                      'outputfiles_basename': './run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1',
                      'n_iter_before_update': 100,
                      'n_live_points': 1000,
                      'sampling_efficiency': 0.3,
                      'const_efficiency_mode': False,
                      'wrapped_params': [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1],
                      'evidence_tolerance': 0.1,
                      'max_iter': -1,
                      'verbose': True}

    xpsi.Sample.MultiNest(likelihood, prior, **runtime_params)


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
        """ A photosphere extension to preload the numerical atmosphere NSX.

        """

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

Spacetime
^^^^^^^^^

.. code-block:: python

    """ CustomSpacetime.py """
    import numpy as np
    import math

    import xpsi

    class CustomSpacetime(xpsi.Spacetime):
        """ A custom spacetime object.

        For the NICER SWG synthetic data parameter recovery exercise, the coordinate
        rotation frequency of the star is fixed.

        """

        def __init__(self, num_params, bounds, S):
            """
            :param int num_params: The number of spacetime parameters.

            :param float S: The coordinate rotation frequency (Hz).

            """
            super(CustomSpacetime, self).__init__(num_params, bounds)

            try:
                self._S = float(S)
            except TypeError:
                raise TypeError('Coordinate spin frequency must be a ``float``.')
            else:
                self._Omega = 2.0 * math.pi * S

Data
^^^^

.. code-block:: python

    """ CustomData.py """

    from __future__ import print_function

    import numpy as np
    import math

    import xpsi

    class CustomData(xpsi.Data):
        """ Custom data container.

        """
        def __init__(self, first, last, counts, phases, exposure_time):
            """
            :param counts: A :class:`numpy.ndarray` of count numbers. The rows of
                           the array must map to a contiguous subset of instrument
                           output channels, with the zeroth row corresponding to
                           the :attr:`first` channel, and the last row
                           corresponding to the channel :attr:`last` minus one.
                           The columns must map to the phases given by
                           :obj:`phases`.

            :param phases: A :class:`numpy.ndarray` of phase *edges* of intervals
                           in which the *synthetic* photons arrive.

            :param exposure_time: The total exposure time in seconds.

            """
            # Execute parent initialisation code
            super(CustomData, self).__init__(first, last)

            try:
                assert isinstance(counts, np.ndarray)
            except AssertionError:
                raise TypeError('Counts object is not a ``numpy.ndarray``.')
            else:
                self._counts = counts

            try:
                assert self._counts.shape[0] == self._last - self._first
            except AssertionError:
                raise AssertionError('The number of rows must be compatible '
                                     'with the first and last output channel '
                                     'numbers.')

            try:
                assert isinstance(phases, np.ndarray)
            except AssertionError:
                raise TypeError('Phases object is not a ``numpy.ndarray``.')
            else:
                self._phases = phases

            self._exposure_time = exposure_time

        @property
        def exposure_time(self):
            """ Get the total exposure time in seconds. """
            return self._exposure_time

        @property
        def counts(self):
            """ Get the photon count data. """
            return self._counts

        @property
        def phases(self):
            """ Get the phases. """
            return self._phases

        @classmethod
        def from_SWG(cls, path, *args):
            """ Constructor which loads photon data from a .txt file.

            :param str path: Path to .txt file which is converted into a
                             two-dimensional :class:`numpy.ndarray`.

            """
            try:
                data = np.loadtxt(path, dtype=np.double)
            except (OSError, IOError, TypeError, ValueError):
                print('Data file could not be loaded.')
                raise

            first = 0; last = 275

            phases = np.linspace(0.0, 1.0, 33)

            return cls(first, last, data, phases, *args)

Instrument
^^^^^^^^^^

.. code-block:: python

    """ CustomInstrument.py """

    from __future__ import print_function, division

    import numpy as np
    import math

    import xpsi

    class CustomInstrument(xpsi.Instrument):
        """ Methods and attributes specific to the NICER instrument.

        """
        def __init__(self, ratio, PI_channels, chan_edges, *args):
            """ Set channel edges attribute. """
            super(CustomInstrument, self).__init__(*args)
            self._ratio = ratio
            self._PI_channels = PI_channels
            self._chan_edges = chan_edges

            self._modified = self.matrix.copy()
            for i in range(self._modified.shape[0]):
                self._modified[i,:] *= self._ratio[i]

        @property
        def channels(self):
            return self._PI_channels

        @property
        def channel_edges(self):
            """ Get the channel edges. """
            return self._chan_edges

        def _construct_matrix(self, p):
            """ Implement response matrix parameterisation. """
            matrix = p[0]*p[1]*self._modified + (1.0 - p[1])*p[2]*self.matrix

            matrix[matrix < 0.0] = 0.0

            return matrix

        def __call__(self, p, signal, *args):
            """ Overwrite. """

            matrix = self._construct_matrix(p)

            self._folded_signal = np.dot(matrix, signal)

            return self._folded_signal

        @classmethod
        def from_SWG(cls, num_params, bounds,
                     ARF, RMF, ratio, max_input, min_input=0, chan_edges=None,
                     offset_correction=None):
            """ Constructor which converts files into :class:`numpy.ndarray`s.

            :param str ARF: Path to ARF which is compatible with
                                    :func:`numpy.loadtxt`.

            :param str RMF: Path to RMF which is compatible with
                                    :func:`numpy.loadtxt`.

            :param str ratio: Path to channel-by-channel ratio file.

            :param str chan_edges: Optional path to edges which is compatible with
                                    :func:`numpy.loadtxt`.

            """
            try:
                ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
                RMF = np.loadtxt(RMF, dtype=np.double, skiprows=3, usecols=-1)
                ratio = np.loadtxt(ratio, dtype=np.double, skiprows=3)[:,2]
                if chan_edges:
                    chan_edges = np.loadtxt(chan_edges, dtype=np.double, skiprows=3)
            except (OSError, IOError, TypeError, ValueError):
                print('A file could not be loaded.')
                raise

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

            PI_channels = np.arange(25, 300)

            ratios = ratio[:275]
            ratios[:10] = ratio[10]

            return cls(ratios, PI_channels, chan_edges[25:301,-2],
                       num_params, bounds, RSP, edges)

Interstellar
^^^^^^^^^^^^

.. code-block:: python

    """ CustomInterstellar.py """

    from __future__ import division

    import numpy as np
    import math

    from scipy.interpolate import Akima1DInterpolator

    import xpsi

    class CustomInterstellar(xpsi.Interstellar):
        """ Apply interstellar absorption. """

        def __init__(self, absorption, **kwargs):

            super(CustomInterstellar, self).__init__(**kwargs)

            self._supplied = absorption[0:351,:]

            self._energies = np.zeros(700, dtype=np.double)
            self._absorption = np.zeros(700, dtype=np.double)

            for i in range(self._supplied.shape[0]-1):
                att_diff = self._supplied[i+1, 1] - self._supplied[i, 1]
                E_diff = self._supplied[i+1, 0] - self._supplied[i, 0]
                self._absorption[2*i] = self._supplied[i,1] + 0.25*att_diff
                self._absorption[2*i+1] = self._supplied[i,1] + 0.75*att_diff
                self._energies[2*i] = self._supplied[i,0] + 0.25*E_diff
                self._energies[2*i+1] = self._supplied[i,0] + 0.75*E_diff

        @property
        def absorption(self):
            return self._absorption

        def __call__(self, p, channel_range, pulse):

            for i in range(pulse.shape[1]):
                pulse[:,i] *= self._absorption**(p[0]/0.4)

        def _interpolate(self, E):
            try:
                self._interpolator
            except AttributeError:
                self._interpolator = Akima1DInterpolator(self._supplied[:,0],
                                                         self._supplied[:,1])
                self._interpolator.extrapolate = True

            return self._interpolator(E)

        def interp_and_absorb(self, p, E, signal):
            """ Interpolate the absorption coefficients and apply. """

            for i in range(signal.shape[1]):
                signal[:,i] *= self._interpolate(E)**(p[0]/0.4)

        @classmethod
        def from_SWG(cls, path, **kwargs):
            """ Load absorption file from the NICER SWG. """

            temp = np.loadtxt(path, dtype=np.double)

            absorption = temp[:,::2]

            return cls(absorption, **kwargs)


Pulse
^^^^^

.. code-block:: python

    """ CustomPulse.py """

    from __future__ import print_function, division

    import numpy as np
    import math

    import xpsi

    from xpsi.likelihoods.default_background_marginalisation import eval_loglike_phaseIntervals_maximise as eval_loglike_maximise
    from xpsi.likelihoods.default_background_marginalisation import precomputation
    from xpsi.tools import phase_interpolator
    from xpsi.tools.phase_integrator import phase_integrator
    from xpsi.tools.synthesise import synthesise
    from xpsi.global_imports import _kpc

    class CustomPulse(xpsi.Pulse):
        """ A custom calculation of the logarithm of the likelihood.

        We extend the :class:`xpsi.Pulse.Pulse` class to make it callable.

        We overwrite the body of the __call__ method. The docstring for the
        abstract method is copied.

        """

        def __init__(self, workspace_intervals = 1000, epsabs = 0, epsrel = 1.0e-8,
                     epsilon = 1.0e-3, sigmas = 10.0, **kwargs):
            """ Perform precomputation. """

            super(CustomPulse, self).__init__(**kwargs)

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

        def __call__(self, p, *args, **kwargs):
            """

            Parameter vector:

            * p[0] = phase shift primary (alias for initial azimuth/phase of photosphere)
            * p[1] = phase shift secondary

            """
            self.shift = np.array(p)

            self.loglikelihood, self.expected_counts, self.background_signal = \
                    eval_loglike_maximise(self._data.exposure_time,
                                          self._data.phases,
                                          self._data.counts,
                                          self._pulse,
                                          self._phases,
                                          self._shift,
                                          self._precomp,
                                          self._workspace_intervals,
                                          self._epsabs,
                                          self._epsrel,
                                          self._epsilon,
                                          self._sigmas,
                                          kwargs.get('llzero'))

        __call__.__doc__ = xpsi.Pulse.__call__.__doc__ + __call__.__doc__

        def synthesise(self):
            """" Overwrite. """

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

    a_f = 0.0
    b_f = 2.0
    a_xi = 0.001
    b_xi = math.pi/2.0 - a_xi

    class CustomPrior(xpsi.Prior):
        """ A custom (joint) prior distribution.

        Source: PSR J0030+0451
        Model variant: ST+PST

        Parameter vector:

        * p[0] = distance (kpc)
        * p[1] = (rotationally deformed) gravitational mass (solar masses)
        * p[2] = coordinate equatorial radius (km)
        * p[3] = inclination of Earth to rotational axis (radians)
        * p[4] = primary centre colatitude (radians)
        * p[5] = primary angular radius (radians)
        * p[6] = primary log10(comoving NSX FIH effective temperature [K])
        * p[7] = secondary centre colatitude (radians)
        * p[8] = secondary angular radius (radians)
        * p[9] = secondary hole colatitude (radians)
        * p[10] = secondary hole angular radius (radians)
        * p[11] = secondary hole azimuth (radians); periodic
        * p[12] = secondary log10(comoving NSX FIH effective temperature [K])
        * p[13] = hydrogen column density (10^20 cm^-2)
        * p[14] = instrument parameter a
        * p[15] = instrument parameter b
        * p[16] = instrument parameter c
        * p[17] = primary cap phase shift (cycles); (alias for initial azimuth, periodic)
        * p[18] = secondary cap phase shift (cycles)

        Note that the unit hypercube to physical transformation is constructed
        for the phases by inverse sampling a flat prior on [-0.25,0.75].
        There is then no need for a periodic boundary and we need to worry about
        accuracy at the boundary.

        """
        def __init__(self, bounds, spacetime):
            # Execute abstract parent initialiser
            super(CustomPrior, self).__init__(bounds)

            assert isinstance(spacetime, xpsi.Spacetime),\
                    'Invalid type for ambient spacetime object.'

            self._spacetime = spacetime

            vals = np.linspace(0.0, b_xi, 1000)

            self._interpolator = Akima1DInterpolator(self._vector_super_radius_mass(vals), vals)
            self._interpolator.extrapolate = True

        def __call__(self, p):
            """ Evaluate distribution at :obj:`p`.

            :param list p: Model parameters values.

            :return: Logarithm of the distribution evaluated at :obj:`p`.

            """
            i = self._spacetime.num_params
            self._spacetime.update(*p[:i])

            if not self._spacetime.R <= 16.0*_km:
                return -np.inf

            if not 1.5 < self._spacetime.R_r_s:
                return -np.inf

            epsilon = self._spacetime.epsilon
            zeta = self._spacetime.zeta
            mu = math.sqrt(-1.0 / (3.0 * epsilon * (-0.788 + 1.030 * zeta)))

            # 2-surface cross-section have a single maximum in |z|
            # i.e., an elliptical surface
            if mu < 1.0:
                return -np.inf

            # polar radius causality for ~static star (static ambient spacetime)
            R_p = 1.0 + epsilon * (-0.788 + 1.030 * zeta)

            if R_p < 1.5 / self._spacetime.R_r_s:
                return -np.inf

            # hot regions cannot overlap
            theta_p = p[4]
            phi_s = (0.5 + p[18]) * _2pi - p[11]
            phi = p[17] * _2pi - phi_s # include ceding azimuth
            rho_p = p[5]

            theta_s = p[7]
            rho_s = p[8]

            ang_sep = xpsi.HotRegion._psi(theta_s, phi, theta_p)

            if ang_sep < rho_p + rho_s:
                return -np.inf

            return 0.0

        @staticmethod
        def _I(x):
            return x * np.log(b_xi/a_xi)

        @staticmethod
        def _II(x):
            return 2.0*(x - a_xi) - x*np.log(x/b_xi)

        def _scalar_super_radius_mass(self, x):
            if x >= a_xi:
                mass = self._II(x)
            else:
                mass = self._I(x)

            return mass

        def _vector_super_radius_mass(self, x):
            masses = np.zeros(len(x))

            for i, _ in enumerate(x):
                masses[i] = self._scalar_super_radius_mass(_)

            masses /= (b_f - a_f)
            masses /= (b_xi - a_xi)

            return masses

        @staticmethod
        def _inverse_sample_cede_radius(x, psi):
            if psi < a_xi:
                return a_xi*np.exp(x * np.log(b_xi/a_xi))
            elif psi >= a_xi and x <= 1.0/(1.0 + np.log(b_xi/psi)):
                return x*psi*(1.0 + np.log(b_xi/psi))
            else:
                return psi*np.exp(x*(1.0 + np.log(b_xi/psi)) - 1.0)

        def inverse_sample(self, hypercube):
            """ Draw sample uniformly from the distribution via inverse sampling.

            :param hypercube: A pseudorandom point in an n-dimensional hypercube.

            :return: A parameter ``list``.

            """
            p = super(CustomPrior, self).inverse_sample(hypercube)

            # distance
            p[0] = truncnorm.ppf(hypercube[0], -10.0, 10.0, loc=0.325, scale=0.009)

            # instrument parameter a
            p[-5] = truncnorm.ppf(hypercube[-5], -5.0, 5.0, loc=1.0, scale=0.1)

            # instrument parameter c
            p[-3] = truncnorm.ppf(hypercube[-3], -5.0, 5.0, loc=1.0, scale=0.1)

            # hole radius
            p[10] = float(self._interpolator(hypercube[10]))

            # cede radius
            p[8] = self._inverse_sample_cede_radius(hypercube[8], p[10])

            if p[10] <= p[8]:
                p[7] = hypercube[7] * (p[8] + p[10])
            else:
                p[7] = p[10] - p[8] + 2.0*hypercube[7]*p[8]

            p[7], p[11] = eval_cedeCentreCoords(p[9], p[7], p[11])

            p[11] *= -1.0

            if p[-2] > 0.5:
                p[-2] -= 1.0

            if p[-1] > 0.5:
                p[-1] -= 1.0

            return p

        def inverse_sample_and_transform(self, hypercube):
            """ A transformation for post-processing. """

            p = self.transform(self.inverse_sample(hypercube))

            return p

        @staticmethod
        def transform(p):
            """ A transformation for post-processing. """

            if not isinstance(p, list):
                p = list(p)

            p += [gravradius(p[1]) / p[2]]

            p += [p[8] - p[10]]

            if p[18] > 0.0:
                p += [p[18] - 1.0]
            else:
                p += [p[18]]

            temp = eval_cedeCentreCoords(-1.0*p[9], p[7], -1.0*p[11])

            azi = temp[1]

            if azi < 0.0:
                azi += 2.0*math.pi

            p += [p[10]/p[8] if p[10] <= p[8] else 2.0 - p[8]/p[10]] # f

            p += [p[8] if p[10] <= p[8] else p[10]] # xi

            p += [temp[0]/(p[8] + p[10]) if p[10] <= p[8] else (temp[0] - p[10] + p[8])/(2.0*p[8])] # kappa

            p += [azi/math.pi]

            return p
