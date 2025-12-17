""" Class for managing likelihood function evaluation.

"""
__all__ = ["Likelihood"]

from xpsi import _rank
from xpsi.utils import make_verbose

from xpsi.global_imports import *

from xpsi.Star import Star
from xpsi.Signal import Signal, LikelihoodError, construct_energy_array
from xpsi.Prior import Prior
from xpsi.ParameterSubspace import ParameterSubspace
from xpsi.EmissionModels import EmissionModel, EmissionModels
from xpsi import HotRegion
from xpsi import Elsewhere

class Likelihood(ParameterSubspace):
    """ A container for all objects related to likelihood evaluation.

    A collective for objects pertaining to a statistical analysis of
    X-ray signals. These objects include X-ray data (sub)sets, the model
    instruments used to acquire these data sets, and the model star and
    model backgrounds.

    :param obj star:
        An instance of :class:`~.Star.Star`. This instance is the model star.

    :param list signals:
        Either:

            * a single :class:`~.Signal.Signal` instance
            * a list of :class:`~.Signal.Signal` instances
            * a list of lists of :class:`~.Signal.Signal` instances

    :param float num_energies:
        Number of energies to compute specific photon flux signals at for
        likelihood evaluation. These energies will be distributed linearly
        in logarithmic space within the union of waveband coverages
        achieved by some set of instruments. Gaps in waveband coverage will
        be skipped.

    :param float fast_rel_num_energies:
        Fraction of the normal number of energies to use in *fast* mode.

    :param int threads:
        The number of ``OpenMP`` threads to spawn for integration. The default
        number of threads used by low-level integration routines is ``1``. The
        number can be increased if the parallelisation paradigm on a cluster
        is changed such that, e.g., per node there is one thread per CPU,
        instead of one ``OpenMPI`` process. It is recommended that
        :obj:`threads` is ``1``; more threads are useful when performing
        integrals at high resolution, but not necessarily when integrating
        *many* times as when sampling, because the MPI parallelisation
        paradigm is invoked.

    :param float llzero:
        The minimum log-likelihood setting for MultiNest. Points whose
        log-likelihood is lower than this value are ignored, which is useful
        for defining complicated priors.

    :param bool externally_updated:
        Update the parameter values upon call to the likelihood object?
        If so, then pass ``False``. A parameter vector then needs to be passed
        to the likelihood object. Otherwise, safely assume that the
        new parameter values are set externally, e.g., in a prior object
        when inverse sampling the prior for nested sampling.

    :param obj prior:
        Instance of subclass of :class:`~.Prior.Prior`.

    :param float max_energy:
        Optional maximum of energy set for signal computation. If no maximum
        is requested (the default), then the maximum is equal to the maximum
        energy from the loaded instrument response models.

    """
    def __init__(self, star, signals,
                 emission_models = None,
                 num_energies = 128,
                 fast_rel_num_energies = 0.25,
                 threads = 1, llzero = -1.0e90,
                 externally_updated = False,
                 prior = None,
                 max_energy = None):

        self.star = star
        self.signals = signals
        self.emission_models = emission_models

        self._do_fast = False

        self._num_energies = num_energies
        self._fast_rel_num_energies = fast_rel_num_energies

        for photosphere, signals in zip(star.photospheres, self._signals):
            try:
                for signal in signals:
                    assert photosphere.prefix == signal.photosphere, \
                        'Each signal subspace must have a photosphere \
                         attribute that matches the identification prefix \
                         of a photosphere object, by convention, and the order \
                         of the list of signal-object lists must match the \
                         order of the list of photosphere objects.'
            except AttributeError:
                pass # quietly assume one photosphere object

            energies = construct_energy_array(num_energies,
                                              list(signals),
                                              max_energy)
            num = int( fast_rel_num_energies * num_energies )
            fast_energies = construct_energy_array(num,
                                                   list(signals),
                                                   max_energy)

            for signal in signals:
                signal.energies = energies
                signal.phases = photosphere.surface.phases_in_cycles

                if photosphere.surface.do_fast:
                    signal.fast_energies = fast_energies
                    signal.fast_phases = photosphere.surface.fast_phases_in_cycles
                    self._do_fast = True

                # Add emission phase arrays
                if self._emission_models is not None:
                    phases = signal.phases
                    for model in self._emission_models:
                        phases.append( model._phases )
                    signal.phases = phases

        self.threads = threads

        self.llzero = llzero

        self.externally_updated = externally_updated

        if prior is not None:
            self.prior = prior

        # merge subspaces
        super(Likelihood, self).__init__(self._star, *(self._signals + [prior]), self._emission_models)

    @property
    def threads(self):
        """ Get the number of threads spawned for integration. """
        return self._threads

    @threads.setter
    def threads(self, n):
        """ Set the number of threads spawned for integration. """
        try:
            self._threads = int(n)
        except TypeError:
            raise TypeError('Thread number must be an integer.')

    @property
    def do_fast(self):
        """ Does fast mode need to be activated? """
        return self._do_fast

    @property
    def star(self):
        """ Get the instance of :class:`~.Star.Star`. """
        return self._star

    @star.setter
    def star(self, obj):
        # only one star object for this version
        if not isinstance(obj, Star):
            raise TypeError('Invalid type for a Star object.')
        else:
            self._star = obj

    @property
    def signals(self):
        """ Get the list of :class:`~.Signal.Signal` instances. """
        if len(self._signals) == 1:
            if len(self._signals[0]) == 1:
                return self._signals[0][0]
        return self._signals

    @signals.setter
    def signals(self, obj):
        # infer how user supplied the signal objects
        if isinstance(obj, Signal):
            self._signals = [[obj,],]
        elif isinstance(obj, list):
            if all(isinstance(o, list) for o in obj):
                for l in obj:
                    if not all(isinstance(o, Signal) for o in l):
                        raise TypeError('Invalid type for a signal object.')

                self._signals = obj
            elif all(isinstance(o, Signal) for o in obj):
                self._signals = [obj,]
            else:
                raise TypeError('Invalid type for a signal object.')

    @property
    def signal(self):
        """ Get the sole signal instance or throw exception. """
        if len(self._signals) > 1 or len(self._signals[0]) > 1:
            raise ValueError('There is more than one signal instance.')
        return self._signals[0][0]

    @property
    def prior(self):
        return self._prior

    @prior.setter
    def prior(self, obj):
        """ Set the prior object if useful, e.g., in nested sampling.

        :param prior:
            A callable instance of :class:`~.Prior.Prior` which is
            useful for nested sampling if inverse prior sampling
            is not possible or not straightforward. If not ``None``
            (default), the prior is called. If the point lies
            exterior to the prior hypervolume, ensure the callable
            returns minus :obj:`numpy.inf`. The purpose is to mix the
            likelihood function and the prior for definition of
            the prior support. Inverse sampling is
            still performed between hard parameter bounds, or on
            any prior factors not handled in the call method. If the
            (conditional) prior to be evaluated is not uniform,
            return the logarithm of the prior density, which is
            combined with the likelihood evaluation. The likelihood
            is not evaluated if the log-prior is infinite.

        """
        if not isinstance(obj, Prior):
            raise TypeError('The prior object is not an instance of the '
                            'Prior class.')

        self._prior = obj
        try:
            self.merge(obj)
        except TypeError:
            pass # assume no prior parameters
        self._prior.parameters = self # a reference to the parameter container

    @prior.deleter
    def prior(self):
        """ Delete the reference to the prior. """

        try:
            del self._prior
        except AttributeError:
            pass # nothing to be done

    @property
    def emission_models(self):
        return self._emission_models

    @emission_models.setter
    def emission_models(self, models):
        if isinstance(models, EmissionModels):
            self._emission_models = models
        elif isinstance(models, EmissionModel):
            self._emission_models = EmissionModels(models)
        elif models is None:
            self._emission_models = None
        else:
            print( 'Warning: emission_models is not an EmissionModels object. No emission models will be used' )
            self._emission_models = None

        # Reference to parameter container into the emission models
        if self._emission_models is not None:
            for model in self._emission_models:
                model.parameters = self # a reference to the parameter container

    @property
    def llzero(self):
        """ Get the minimum log-likelihood setting passed to MultiNest. """
        return self._llzero

    @llzero.setter
    def llzero(self, value):
        """ Set the minimum log-likelihood that is viewed as zero. """
        try:
            self._llzero = float(value)
        except TypeError:
            raise TypeError('Zero log-likelihood number must be a float.')

    @property
    def random_near_llzero(self):
        """ Get the minimum log-likelihood scaled randomly by up to an
            order of magnitude. """
        return self._llzero * (0.1 + 0.9*_np.random.rand(1))

    @property
    def externally_updated(self):
        """ Safely assume parameters are already updated upon call to self? """
        return self._externally_updated

    @externally_updated.setter
    def externally_updated(self, updated):
        if not isinstance(updated, bool):
            raise TypeError('Provide a boolean to define update protocol.')
        self._externally_updated = updated

    @staticmethod
    def _divide(obj, x):
        """ Helper operator to check for compatibility first.

        As an example, if fast mode is activated for some hot region but not
        another, :obj:`obj` would be ``None`` and thus a safeguard is needed.

        """
        if isinstance(obj, _np.ndarray):
            return obj / x
        else:
            return None

    def _driver(self, fast_mode=False, synthesise=False, force_update=False, **kwargs):
        """ Main likelihood evaluation driver routine. """

        self._star.activate_fast_mode(fast_mode)

        star_updated = False
        if self._star.needs_update or force_update: # ignore fast parameters in this version
            try:
                if fast_mode or not self._do_fast:
                    fast_total_counts = None
                else:
                    fast_total_counts = tuple(signal.fast_total_counts for\
                                                        signal in self._signals)

                self._star.update(fast_total_counts, self.threads,force_update=force_update)

            except xpsiError as e:
                if isinstance(e, HotRegion.RayError):
                    print('Warning: HotRegion.RayError raised.')
                elif isinstance(e, Elsewhere.RayError):
                    print('Warning: Elsewhere.RayError raised.')


                return self.random_near_llzero

            for photosphere, signals in zip(self._star.photospheres, self._signals):
                try:
                    if fast_mode:
                        energies = signals[0].fast_energies
                    else:
                        energies = signals[0].energies

                    photosphere.integrate(energies, self.threads)
                except xpsiError as e:
                    try:
                        prefix = ' prefix ' + photosphere.prefix
                    except AttributeError:
                        prefix = ''
                    if isinstance(e, HotRegion.PulseError):
                        print('Warning: HotRegion.PulseError raised for '
                              'photosphere%s.' % prefix)
                    elif isinstance(e, Elsewhere.IntegrationError):
                        print('Warning: Elsewhere.IntegrationError for '
                              'photosphere%s.' % prefix)
                    elif isinstance(e, HotRegion.AtmosError):
                        raise
                    elif isinstance(e, Elsewhere.AtmosError):
                        raise

                    print('Parameter vector: ', super(Likelihood,self).__call__())
                    return self.random_near_llzero

            # Add emission models from outside the photosphere
            if self._emission_models is not None:
                self._emission_models.update(self.threads,force_update=force_update)
                self._emission_models.integrate(signals[0].energies, self.threads)

            star_updated = True

        # register the signals by operating with the instrument response
        for signals, photosphere in zip(self._signals, self._star.photospheres):
            for signal in signals:
                if star_updated or signal.needs_update:
                    
                    # Define which photosphere signal to use
                    if signal.isI:
                        photosphere_signal = photosphere.signal
                    elif signal.isQ or signal.isQn:
                        photosphere_signal =  photosphere.signalQ
                    elif signal.isU or signal.isUn:
                        photosphere_signal =  photosphere.signalU
                    else:
                        raise TypeError('Signal type must be either I, Q, U, Qn, or Un.')
                    
                    # Apply this choice to register appropriate signal
                    signal_to_register = tuple( tuple(self._divide(component,
                                                                self._star.spacetime.d_sq)
                                                            for component in hot_region)
                                                    for hot_region in photosphere_signal)
                    signal.register(signal_to_register, fast_mode=fast_mode, threads=self.threads)

                    # Normalize if required
                    if signal.isQn or signal.isUn:

                        # Save old value
                        polarized_signal = signal.signals

                        # Compute the I component
                        signal_to_nomalize = tuple( tuple(self._divide(component,
                                                                self._star.spacetime.d_sq)
                                                            for component in hot_region)
                                                    for hot_region in photosphere.signal)
                        signal.register(signal_to_nomalize, fast_mode=fast_mode, threads=self.threads)
                        Isignal = signal.signals

                        # Normalize
                        for ihot in range(len(polarized_signal)):
                            signal._signals[ihot]=_np.where(Isignal[ihot]==0.0, 0.0, polarized_signal[ihot]/Isignal[ihot])

                    # Add emission models from outside the photosphere
                    if self._emission_models is not None:
                        signal.register(tuple( tuple( component for component in model)
                                                    for model in self._emission_models.signal),
                                        fast_mode=fast_mode, threads=self.threads, reset=False)

                    reregistered = True
                else:
                    reregistered = False

                if not fast_mode and reregistered:
                    if synthesise:
                        hot = photosphere.surface
                        try:
                            kws = kwargs.pop(signal.prefix)
                        except AttributeError:
                            kws = {}

                        shifts = [h['phase_shift'] for h in hot.objects]
                        signal.shifts = _np.array(shifts)
                        signal.synthesise(threads=self._threads, **kws)
                    else:
                        try:
                            hot = photosphere.surface
                            shifts = [h['phase_shift'] for h in hot.objects]

                            # Add model shifts if needed
                            if self._emission_models is not None:
                                shifts_emission_models = [model['phase_shift'] for model in self._emission_models]
                                shifts = shifts + shifts_emission_models

                            signal.shifts = _np.array(shifts)
                            signal(threads=self._threads, llzero=self._llzero)
                        except LikelihoodError:
                            try:
                                prefix = ' prefix ' + signal.prefix
                            except AttributeError:
                                prefix = ''
                            print('Warning: LikelihoodError raised for '
                                  'signal%s.' % prefix)
                            print('Parameter vector: ', super(Likelihood,self).__call__())
                            return self.random_near_llzero


        return star_updated

    def reinitialise(self):
        """ Reinitialise the likelihood object.

        Useful if some resolution settings in child objects were changed
        (namely, the number of signal phases) that need to be communicated
        to other child objects.

        """

        self.__init__(self._star,
                      self._signals,
                      self._emission_models,
                      self._num_energies,
                      self._fast_rel_num_energies,
                      self._threads,
                      self._llzero)

    def __call__(self, p = None, reinitialise = False, force = False):
        """ Evaluate the logarithm of the joint likelihood over all pulsations.

        :param list p:
            Parameter vector if parameters not updated externally.

        :param optional[bool] reinitialise:
            Call ``self.reinitialise()``?

        :param optional[bool] force:
            Force complete reevaluation even if some parameters are unchanged.
            To faciliate this, all parameter caches are cleared.

        :returns: The logarithm of the likelihood.

        """
        if reinitialise: # for safety if settings have been changed
            self.reinitialise() # do setup again given exisiting object refs
            self.clear_cache() # clear cache and values
        elif force: # no need to reinitialise, just clear cache and values
            self.clear_cache()

        if not self.externally_updated: # do not safely assume already handled
            if p is None: # expected a vector of values instead of nothing
                raise TypeError('Parameter values have not been updated.')
            super(Likelihood, self).__call__(p) # update free parameters

        if self.needs_update or force:
            try:
                logprior = self._prior(p) # pass vector just in case wanted
            except AttributeError:
                pass
            else:
                if not _np.isfinite(logprior):
                    # we need to restore due to premature return
                    super(Likelihood, self).__call__(self.cached)
                    return self.random_near_llzero

            if self._do_fast:
                # perform a low-resolution precomputation to direct cell
                # allocation
                x = self._driver(fast_mode=True,force_update=force)
                if not isinstance(x, bool):
                    super(Likelihood, self).__call__(self.cached) # restore
                    return x
                elif x:
                    x = self._driver(force_update=force)
                    if not isinstance(x, bool):
                        super(Likelihood, self).__call__(self.cached) # restore
                        return x
            else:
                x = self._driver(force_update=force)
                if not isinstance(x, bool):
                    super(Likelihood, self).__call__(self.cached) # restore
                    return x

            # memoization: update parameter value caches
            super(Likelihood, self).__call__(self.vector)

        loglikelihood = 0.0
        for signals in self._signals:
            for signal in signals:
                try:
                    loglikelihood += signal.loglikelihood
                except AttributeError as e:
                    print("ERROR: It looks like X-PSI falsely thought that the signal does not need to be updated and thus skipped an essential part of the calculation. If not sampling, please use ``force=True`` or ``force_update=True`` option for the likelihood evaluation, or if sampling please set ``likelihood.externally_updated = True``")
                    raise

        if loglikelihood <= self.llzero:
            return self.random_near_llzero

        try:
            return loglikelihood + logprior
        except NameError:
            return loglikelihood

    @make_verbose('Checking likelihood and prior evaluation before '
                  'commencing sampling', 'Checks passed')
    def check(self,
              hypercube_points,
              loglikelihood_call_vals,
              rtol_loglike,
              atol_loglike=0.0,
              logprior_call_vals=None,
              rtol_logprior=None,
              atol_logprior=None,
              physical_points=None,
              force_update=False,
              numpy_allclose=False):
        """ Perform checks on the likelihood evaluator and the prior density.

        Can be called from :func:`~.Sample.nested` to execute a check before
        automatically commencing a sampling process.

        :param ndarray[n,m] hypercube_points:
            A set of ``n`` points in the unit hypercube, where ``m`` is
            dimensionality (``self.num_params``) of the sampling space -- i.e.,
            of the hypercube. If you want to pass the physical points instead,
            just pass ``None`` here.

        :param optional(ndarray[n,m]) physical_points:
            A set of ``n`` points in the physical parameter space, where
            ``m`` is dimensionality (``self.num_params``) of the sampling space.
            The ``hypercube_points``, if not ``None``, will be ignored.

        :param optional[bool] force_update:
            Force everything to be re-calculated regardless of what was computed before.
            This can be used to prevent errors in cases when the automatic check 
            for update need is not working as intended.

        :param optional[bool] numpy_allclose:
            Determine whether the allclose function of numpy is used when evaluating
            the closeness of the given and calculated likelihood. By default, a fallback
            implementation is used, which also prints the likelihood values.

        """
        if numpy_allclose:
            from numpy import allclose
        else:
            yield 'Not using ``allclose`` function from NumPy.'
            yield 'Using fallback implementation instead.'

            @make_verbose('Checking closeness of likelihood arrays:',
                          'Closeness evaluated')
            def allclose(a, b, rtol, atol, equal_nan=None):
                """ Fallback based on NumPy v1.17. """
                for _a, _b in zip(a, b):
                    yield '%.10e | %.10e .....' % (_a, _b)
                yield ~((_np.abs(a - b) > atol + rtol*_np.abs(b)).any())

        lls = []
        lps = [] if logprior_call_vals is not None else None


        if physical_points is not None:
            physical_points = _np.array(physical_points)

            try:
                if physical_points.ndim > 2:
                    raise AttributeError("Dimension of physical_points is > 2")

                elif physical_points.ndim == 1:
                    physical_points = physical_points.reshape(1,len(physical_points))

            except AttributeError as e:
                e.message = ('Physical points for likelihood check must form '
                             'a numpy.ndarray with at most two dimensions.')
                raise

            if physical_points.shape[1] != len(self):
                raise IndexError('Vector size does not match number of '
                                 'free parameters of likelihood function.')

            _cached = self.externally_updated

            self.externally_updated = False

            for point in physical_points:
                lls.append(self.__call__(point,force=force_update))
                if lps is not None:
                    lps.append(self._prior(point))

            self.externally_updated = _cached
        else:
            hypercube_points = _np.array(hypercube_points)

            try:
                if hypercube_points.ndim > 2:
                    raise AttributeError

                elif hypercube_points.ndim == 1:
                    hypercube_points = hypercube_points.reshape(1,len(hypercube_points))

            except AttributeError as e:
                e.message = ('Unit hypercube points for likelihood check must '
                             'form a numpy.ndarray with at most two dimensions.')
                raise

            if hypercube_points.shape[1] != len(self):
                raise IndexError('Vector size does not match number of '
                                 'free parameters of likelihood function.')

            for point in hypercube_points:
                phys_point = self._prior.inverse_sample(point)
                lls.append(self.__call__(phys_point,force=force_update))
                if lps is not None:
                    lps.append(self._prior(phys_point))

        try:
            proceed = allclose(_np.array(lls),
                               _np.array(loglikelihood_call_vals),
                               rtol=rtol_loglike, atol=atol_loglike,
                               equal_nan=False)
        except Exception as e:
            raise Exception('Log-likelihood value checks failed on process %i '
                            'with exception value: %s' % (_rank, str(e)))
        else:
            if proceed:
                yield 'Log-likelihood value checks passed on root process.'
            else:
                raise ValueError('Log-likelihood value checks failed on process '
                                 '%i' % _rank)

        if lps is not None:
            try:
                proceed = allclose(_np.array(lps),
                                   _np.array(logprior_call_vals),
                                   rtol=rtol_logprior, atol=atol_logprior,
                                   equal_nan=False)
            except Exception as e:
                raise Exception('Log-prior value checks failed on process %i '
                                'with exception value: %s' % (_rank, str(e)))
            else:
                if proceed:
                    yield 'Log-prior value checks passed on root process.'
                else:
                    raise ValueError('Log-prior value checks failed on '
                                     'process %i' % _rank)

    def synthesise(self, p, reinitialise=False, force=False, **kwargs):
        """ Synthesise pulsation data.

        :param list p: Parameter vector.

        :param optional[bool] reinitialise: Call ``self.reinitialise()``?

        :param optional[bool] force:
            Force complete reevaluation even if some parameters are unchanged.

        :param dict kwargs:
            Keyword arguments propagated to custom signal synthesis methods.
            Examples of such arguments include exposure times or
            required total count numbers (see example notebooks).

        """
        if reinitialise: # for safety if settings have been changed
            self.reinitialise() # do setup again given exisiting object refs
            self.clear_cache() # clear cache and values
        elif force: # no need to reinitialise, just clear cache and values
            self.clear_cache()

        if p is None: # expected a vector of values instead of nothing
            raise TypeError('Parameter values are None.')
        super(Likelihood, self).__call__(p) # update free parameters

        try:
            logprior = self._prior(p) # pass vector just in case wanted
        except AttributeError:
            pass
        else:
            if not _np.isfinite(logprior):
                because_of_1D_bounds = False
                for param in self._prior.parameters:
                    if param.bounds[0] is not None:
                        if not param.bounds[0] <= param.value:
                            because_of_1D_bounds = True
                    if param.bounds[1] is not None:
                        if not param.value <= param.bounds[1]:
                            because_of_1D_bounds = True
                if because_of_1D_bounds:
                    print("Warning: Prior check failed, because at least one of the parameters is not within the hard 1D-bounds. No synthetic data will be produced.")
                else:
                    print('Warning: Prior check failed because a requirement set in CustomPrior has failed. No synthetic data will be produced.')
                # we need to restore due to premature return
                super(Likelihood, self).__call__(self.cached)
                return None

        if self._do_fast:
            # perform a low-resolution precomputation to direct cell
            # allocation
            x = self._driver(fast_mode=True,force_update=force)
            if not isinstance(x, bool):
                super(Likelihood, self).__call__(self.cached) # restore
                return None
            elif x:
                x = self._driver(synthesise=True,force_update=force, **kwargs)
                if not isinstance(x, bool):
                    super(Likelihood, self).__call__(self.cached) # restore
                    return None
        else:
            x = self._driver(synthesise=True,force_update=force, **kwargs)
            if not isinstance(x, bool):
                super(Likelihood, self).__call__(self.cached) # restore
                return None
