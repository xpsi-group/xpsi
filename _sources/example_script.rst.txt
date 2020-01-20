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

The script below is an example like that found on the `Modelling`
page, but without being interleaved with verbose explanations. The model
implemented here is more involved than shown on the `Modelling`
page, and than those defined for the example parameter
estimation problems that can be found in the ``xpsi/examples`` directory of
the repository.

.. note::

    The script and modules below are being updated to use
    the current API (``v0.3.2``).

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
    from CustomPrior import CustomPrior
    from CustomPhotosphere import CustomPhotosphere

    from xpsi.global_imports import _c, _G, _M_s, _dpr, gravradius

    data = CustomData.from_SWG('data/NICER_J0030_PaulRay_fixed_evt_25to299__preprocessed.txt', 1936864.0)

    # bounds on instrument model parameters
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

    pulse = CustomPulse(data = data,
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

    spacetime = xpsi.Spacetime.fixed_spin(1.0/(4.87e-3))

    bounds = dict(super_colatitude = (0.001, math.pi - 0.001),
                  super_radius = (0.001, math.pi/2.0 - 0.001),
                  phase_shift = (0.0, 0.2), # defined relative to 0.35 cycles
                  super_temperature = (5.1, 6.8))

    bounds = [(0.35, 0.55), (-0.25,0.75)],

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

    likelihood = xpsi.Likelihood(star = star, pulses = pulse, threads=1,
                                 externally_updated = True)

    prior = CustomPrior()

    likelihood.prior = prior

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

    # see CustomPrior docstring for parameter names
    p = [0.140337033600940120E+01,
            0.133784624585842025E+02,
            0.328978844399083370E+00,
            0.100434973113637094E+01,
            0.454255509351488285E+00,
            0.219377527309307840E+01,
            0.791608842011687908E-01,
            0.610655622382022134E+01,
            0.476829413031657379E+00,
            0.271629852479304956E+01,
            0.322342254787806259E+00,
            0.274633014642517770E+01,
            0.284416965175110226E+00,
            -0.483260905056053860E-01,
            0.611730491798804454E+01,
            0.460499862995095377E+00,
            0.103356827187160971E+01,
            0.222710719836020192E-01,
            0.874856631973894849E+00]

    # let's require that checks pass before starting to sample
    check_kwargs = dict(hypercube_points = None,
                        physical_points = p,
                        loglikelihood_call_vals = [-36316.35439439],
                        rtol_loglike = 1.0e-8)

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
    from xpsi import Parameter

    class CustomInstrument(xpsi.Instrument):
    """ Methods and attributes specific to the NICER instrument.

    Currently tailored to the NICER light-curve SWG model specification.

    """
    def __init__(self, ratio, channels, channel_edges, *args):
        """ Set channel edges attribute. """
        super(CustomInstrument, self).__init__(*args)

        self._ratio = ratio
        self._channels = channels
        self._channel_edges = channel_edges

        self._modified = self.matrix.copy()
        for i in range(self._modified.shape[0]):
            self._modified[i,:] *= self._ratio[i]

    @property
    def channels(self):
        return self._channels

    @property
    def channel_edges(self):
        """ Get the channel edges. """
        return self._channel_edges

    def construct_matrix(self):
        """ Implement response matrix parameterisation. """
        matrix = self['alpha']*self['beta']*self._modified
        matrix += (1.0 - self['beta'])*self['gamma']*self.matrix

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
                 ARF, RMF, ratio,
                 max_input, min_input=0,
                 channel_edges=None):
        """ Constructor which converts files into :class:`numpy.ndarray`s.

        :param str ARF: Path to ARF which is compatible with
                                :func:`numpy.loadtxt`.

        :param str RMF: Path to RMF which is compatible with
                                :func:`numpy.loadtxt`.

        :param str ratio: Path to channel-by-channel ratio file.

        :param str channel_edges: Optional path to edges which is compatible with
                                :func:`numpy.loadtxt`.

        """
        ARF = np.loadtxt(ARF, dtype=np.double, skiprows=3)
        RMF = np.loadtxt(RMF, dtype=np.double, skiprows=3, usecols=-1)
        ratio = np.loadtxt(ratio, dtype=np.double, skiprows=3)[:,2]

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

        ratios = ratio[:275]
        ratios[:10] = ratio[10]

        alpha = Parameter('alpha',
                          strict_bounds = (0.0,2.0),
                          bounds = bounds.get('alpha', None),
                          doc = 'alpha',
                          symbol = r'$\alpha$',
                          value = values.get('alpha', None))

        beta = Parameter('beta',
                          strict_bounds = (0.0,1.0),
                          bounds = bounds.get('beta', None),
                          doc = 'beta',
                          symbol = r'$\beta$',
                          value = values.get('beta', None))

        gamma = Parameter('gamma',
                          strict_bounds = (0.0,2.0),
                          bounds = bounds.get('gamma', None),
                          doc = 'gamma',
                          symbol = r'$\gamma$',
                          value = values.get('gamma', None))

        return cls(ratios, channels, channel_edges,
                   RSP, edges, alpha, beta, gamma)

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

        def __init__(self, absorption, bounds, values = {}):

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

            N_H = Parameter('column_density',
                            strict_bounds = (0.0,10.0),
                            bounds = bounds.get('column_density', None),
                            doc = 'Units of 10^20 cm^-2.',
                            symbol = r'$N_{\rm H}$',
                            value = values.get('column_density', None))

            super(CustomInterstellar, self).__init__(N_H)

        @property
        def absorption(self):
            return self._absorption

        def __call__(self, energies, pulse):

            for i in range(pulse.shape[1]):
                pulse[:,i] *= self._absorption**(self['column_density']/0.4)

        def _interpolate(self, E):
            """ Helper. """
            try:
                self._interpolator
            except AttributeError:
                self._interpolator = Akima1DInterpolator(self._supplied[:,0],
                                                         self._supplied[:,1])
                self._interpolator.extrapolate = True

            return self._interpolator(E)

        def interp_and_absorb(self, E, signal):
            """ Interpolate the absorption coefficients and apply.

            Useful for post-processing.

            """

            for i in range(signal.shape[1]):
                signal[:,i] *= self._interpolate(E)**(self['column_density']/0.4)

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

        def __call__(self, phase_shifts, *args, **kwargs):
            self.shift = np.array(phase_shifts)

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

        Parameter vector: (print the likelihood object)

        * p[0] = (rotationally deformed) gravitational mass (solar masses)
        * p[1] = coordinate equatorial radius (km)
        * p[2] = distance (kpc)
        * p[3] = inclination of Earth to rotational axis (radians)
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

        a_f = 0.0
        b_f = 2.0
        a_xi = 0.001
        b_xi = math.pi/2.0 - a_xi

        vals = np.linspace(0.0, b_xi, 1000)

        interpolator = Akima1DInterpolator(self._vector_super_radius_mass(vals), vals)
        interpolator.extrapolate = True

        def __init__(self):
            """ Nothing to be done. """
            pass

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
            #if R_p < 1.5 / ref.R_r_s:
            #    return -np.inf

            # limit polar radius to try to exclude deflections >= \pi radians
            # due to oblateness this does not quite eliminate all configurations
            # with deflections >= \pi radians
            R_p = 1.0 + ref.epsilon * (-0.788 + 1.030 * ref.zeta)
            if R_p < 1.76 / ref.R_r_s:
                return -np.inf

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

        @staticmethod
        def _I(x):
            return x * np.log(self.b_xi/self.a_xi)

        @staticmethod
        def _II(x):
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

        @staticmethod
        def _inverse_sample_cede_radius(x, psi):
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

            super(CustomPrior, self).inverse_sample(hypercube)

            ref = self.parameters # redefine shortcut

            idx = ref.index('distance')
            ref['distance'] = truncnorm.ppf(hypercube[idx],
                                            -10.0, 10.0,
                                            loc=0.325, scale=0.009)

            idx = ref.index('p__phase_shift')
            ref['p_phase_shift'] = 0.35 + 0.2 * hypercube[idx]
            if ref['p__phase_shift'] > 0.5:
                ref['p__phase_shift'] -= 1.0

            idx = ref.index('s__phase_shift')
            ref['s_phase_shift'] = -0.25 + hypercube[idx]
            if ref['s__phase_shift'] > 0.5:
                ref['s__phase_shift'] -= 1.0

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

        def inverse_sample_and_transform(self, hypercube = None):
            """ A transformation for post-processing. """

            p = self.transform(self.inverse_sample(hypercube))

            return p

        @staticmethod
        def transform(p):
            """ A transformation for post-processing. """

            if not isinstance(p, list):
                p = list(p)

            p += [gravradius(p[0]) / p[1]]

            p += [p[10] - p[12]]

            if p[8] > 0.0:
                p += [p[8] - 1.0]
            else:
                p += [p[8]]

            temp = eval_cedeCentreCoords(-1.0*p[11], p[9], -1.0*p[13])

            azi = temp[1]

            if azi < 0.0:
                azi += 2.0*math.pi

            p += [p[12]/p[10] if p[12] <= p[10] else 2.0 - p[10]/p[12]] # f

            p += [p[10] if p[12] <= p[10] else p[12]] # xi

            p += [temp[0]/(p[10] + p[12]) if p[12] <= p[10] else (temp[0] - p[12] + p[10])/(2.0*p[10])] # kappa

            p += [azi/math.pi]

            return p
