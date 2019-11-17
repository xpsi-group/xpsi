""" Class for managing likelihood function evaluation.

"""

from __future__ import division, print_function

__all__ = ["Likelihood", "TagError"]

from . import _rank
from . import make_verbose

from .global_imports import *
from . import global_imports

from .Star import Star
from .Pulse import Pulse, LikelihoodError
from .Background import Background
from .Prior import Prior
from .ParameterSpace import ParameterSpace
from . import HotRegion
from . import TwoHotRegions
from . import Elsewhere

class TagError(xpsiError):
    """
    Raised if an ordered list of photosphere objects and an ordered list
    of pulse objects do have matching identification tags.

    """

class Likelihood(ParameterSpace):
    """ A container for all objects related to likelihood evaluation.

    A collective for objects pertaining to a statistical analysis of
    X-ray pulses. These objects include X-ray data (sub)sets, the model
    instruments used to acquire these data sets, and the model star and
    model backgrounds.

    :param obj star: An instance of :class:`~.Star.Star`. This instance
                     is the model star.

    :param list pulses: A list of :class:`~.Pulse.Pulse`
                        instances.

    :param int threads: The number of ``OpenMP`` threads to spawn for
                        integration. The default number of threads used
                        by low-level integration routines is ``1``. The
                        number can be increased if the parallelisation
                        paradigm on a cluster is changed such that, e.g.,
                        per node there is one thread per CPU, instead of
                        one ``OpenMPI`` process. It is recommended that
                        :obj:`threads` is ``1``; more threads are useful
                        when performing integrals at high resolution, and
                        not necessarily when integrating *many* times as
                        when sampling.

    :param float llzero: The minimum log-likelihood setting for
                         MultiNest. Points whose log-likelihood is lower
                         than this value are ignored, which is useful for
                         defining complicated priors.

    """
    def __init__(self, star, pulses, threads=1, llzero=-1.0e90):
        if not isinstance(star, Star):
            raise TypeError('Invalid type for a Star object.')
        else:
            self._star = star

        try:
            iter(pulses)
        except TypeError:
            if isinstance(pulses, Pulse):
                self._pulses = [pulses]
            else:
                raise TypeError('Invalid type for a pulse object.')
        else:
            for pulse in pulses:
                if not isinstance(pulse, Pulse):
                    raise TypeError('Invalid type for a pulse object.')
            self._pulses = pulses

        self._do_fast = False

        for photosphere, pulse in zip(star.photospheres, self._pulses):
            assert photosphere.tag == pulse.tag, \
                'Photosphere and pulse object lists must have matching \
                 identification tags.'

            pulse.phases = photosphere.hot.phases_in_cycles
            pulse.fast_phases = photosphere.hot.fast_phases_in_cycles
            if photosphere.hot.do_fast: self._do_fast = True

        self._theta = [0.0] * self.num_params

        try:
            assert int(threads) == threads, 'Thread number must be an integer.'
        except TypeError:
            raise TypeError('Thread number must be an integer.')
        else:
            self._threads = int(threads)

        self._llzero = llzero

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
    def num_params(self):
        """ Get the number of dimensions of the parameter space. """

        n = self._star.num_params

        for pulse in self._pulses:
            n += pulse.total_params

        return n

    @property
    def bounds(self):
        """ Get a list of tuples of hard parameter bounds. """

        b = []

        b += self._star.spacetime.bounds

        for photosphere in self._star.photospheres:
            b += photosphere.bounds
            b += photosphere.hot.bounds
            if photosphere.elsewhere is not None:
                b += photosphere.elsewhere.bounds

        for pulse in self._pulses:
            if pulse.interstellar is not None:
                b += pulse.interstellar.bounds
            if pulse.instrument is not None:
                b += pulse.instrument.bounds
            if pulse.background is not None:
                b += pulse.background.bounds
            b += pulse.bounds
        return b

    def __str__(self):
        """ Print the model objects in order, with parameter numbers. """
        print('Printing the model objects in order...\n')
        print('Object [tag]: number of parameters')
        print('--------------------------------------')
        print('Spacetime: %i' % self._star.spacetime.num_params)
        for photosphere in self._star.photospheres:
            print("Photosphere ['%s']: %i"
                    % (photosphere.tag, photosphere.total_params))
            print("   --> Hot ['%s']: %i"
                    % (photosphere.tag, photosphere.hot.num_params))
            try:
                print("   --> Elsewhere ['%s']: %i"
                        % (photosphere.tag, photosphere.elsewhere.num_params))
            except AttributeError:
                pass
        for pulse in self._pulses:
            print("Pulse ['%s']: %i" % (pulse.tag, pulse.total_params))
        print('--------------------------------------')

        return 'Total number of model parameters: %i' % self.num_params

    @property
    def star(self):
        """ Get the instance of :class:`~.Star.Star`. """
        return self._star

    @property
    def pulses(self):
        """ Get the list of :class:`~.Pulse.Pulse` instances. """
        return self._pulses

    @property
    def theta(self):
        """ Get the last parameter vector for joint likelihood evaluation.

        :setter: Stores the last parameter vector so that if only the fast
                 parameters change, recompute only the fast contribution the
                 likelihood.

        :param p: A list of model parameters values.

        """
        return self._theta

    @property
    def prior(self):
        return self._prior

    @prior.setter
    def prior(self, obj):
        """ Set the prior object if useful, e.g., in nested sampling.

        :param prior: A callable instance of :class:`~.Prior.Prior` which is
                      useful for nested sampling if inverse prior sampling
                      is not possible or not straightforward. If not ``None``
                      (default), the prior is called. If the point lies
                      exterior to the prior hypervolume, ensure the callable
                      returns minus :obj:`numpy.inf`. The purpose is to mix the
                      likelihood function and the prior for definition of
                      hard nested sampling boundaries. Inverse sampling is
                      still performed between hard parameter bounds, or on
                      any prior factors not handled in the call method. If the
                      (conditional) prior to be evaluated is not uniform,
                      return the logarithm of the prior density, which is
                      combined with the likelihood evaluation. The likelihood
                      is not evaluated if the log-prior is infinite.

        """
        if isinstance(obj, Prior):
            self._prior = obj
        else:
            raise TypeError('The prior object is not an instance of the '
                            'Prior class.')

    @property
    def llzero(self):
        """ Get the minimum log-likelihood setting passed to MultiNest. """
        return self._llzero

    @property
    def random_near_llzero(self):
        """ Get the minimum log-likelihood scaled randomly by up to an
            order of magnitude. """
        return self._llzero * (0.1 + 0.9*_np.random.rand(1))

    @property
    def less_than_llzero(self):
        """ Get a number less than the minimum log-likelihood threshold. """
        return 1.1 * self._llzero

    @staticmethod
    def _divide(obj, x):
        if isinstance(obj, _np.ndarray):
            return obj / x
        else:
            return None

    def _driver(self, p, fast_mode=False, synthesise=False, **kwargs):
        """ Main likelihood evaluation driver routine. """

        self._star.activate_fast_mode(fast_mode)

        star_updated = False

        s = self._star.num_params
        if (p[1:s] != self._theta[1:s]).any():
            try:
                if fast_mode or not self._do_fast:
                    fast_total_counts = tuple([None]*len(self._pulses))
                else:
                    fast_total_counts = tuple(pulse.fast_total_counts for\
                                                        pulse in self._pulses)

                self._star.update(p[:s], fast_total_counts, self.threads)
            except xpsiError as e:
                if isinstance(e, HotRegion.RayError):
                    print('Warning: HotRegion.RayError raised.')
                elif isinstance(e, Elsewhere.RayError):
                    print('Warning: Elsewhere.RayError raised.')

                print('Parameter vector: ', p)
                return self.random_near_llzero

            for photosphere, pulse in zip(self._star.photospheres, self._pulses):
                try:
                    try:
                        if fast_mode:
                            energies = pulse.fast_energies
                        elif self._do_fast:
                            energies = pulse.energies
                        else:
                            energies = pulse.default_energies
                    except AttributeError:
                        energies = pulse.logspace_energies

                    photosphere.integrate(energies, self.threads)
                except xpsiError as e:
                    if isinstance(e, HotRegion.PulseError):
                        print('Warning: HotRegion.PulseError raised for '
                              'photosphere tag %s.' % photosphere.tag)
                    elif isinstance(e, Elsewhere.IntegrationError):
                        print('Warning: Elsewhere.IntegrationError for '
                              'photosphere tag %s.' % photosphere.tag)

                    print('Parameter vector: ', p)
                    return self.random_near_llzero

            star_updated = True

        # update distance if changed
        if not star_updated and p[0] != self._theta[0]:
            self._star.spacetime.d = p[0]
            star_updated = True

        # fold the pulse through the instrument response
        i = s
        for pulse, photosphere in zip(self._pulses, self._star.photospheres):
            j = i + pulse.total_params - pulse.num_params

            if star_updated or (p[i:j] != self._theta[i:j]).any():
                pulse.fold(tuple(
                                 tuple(self._divide(component,
                                                    self._star.spacetime.d_sq)
                                       for component in hotRegion)
                                 for hotRegion in photosphere.pulse),
                           p[i:j], fast_mode=fast_mode, threads=self.threads)
                refolded = True
            else:
                refolded = False

            if not fast_mode:
                k = j + pulse.num_params
                if refolded or (p[j:k] != self._theta[j:k]).any():
                    if synthesise:
                        pulse.synthesise(p[j:k],
                                         threads=self._threads,
                                         **kwargs)
                    else:
                        try:
                            pulse(p[j:k], *[list(p[:s]) + list(p[i:j])],
                                  threads=self._threads, llzero=self._llzero)
                        except LikelihoodError:
                            print('Warning: LikelihoodError raised for '
                                  'pulse tag %s.' % pulse.tag)
                            print('Parameter vector: ', p)
                            return self.random_near_llzero

            i = j

        return star_updated

    def __call__(self, p, reinitialise=False):
        """ Evaluate the logarithm of the joint likelihood over all pulsations.

        :param list p: Parameter vector.
        :param optional[bool] reinitialise:
            Reinitialise? Useful if some resolution settings in child
            objects were changed (namely, the number of pulse phases) that
            need to be communicated to other child objects.

        :return: The logarithm of the likelihood.

        """
        if reinitialise:
            self.__init__(self._star,
                          self._pulses,
                          self._threads,
                          self._llzero)

        p = _np.array(p)

        if (p != self._theta).any():
            try:
                logprior = self._prior(p)
            except AttributeError:
                pass
            else:
                if not _np.isfinite(logprior):
                    self._theta = p
                    return self.less_than_llzero

            if self._do_fast:
                # perform a low-resolution precomputation to direct cell
                # allocation
                x = self._driver(p, fast_mode=True)
                if not isinstance(x, bool):
                    self._theta = p
                    return x
                elif x:
                    x = self._driver(p)
                    if not isinstance(x, bool):
                        self._theta = p
                        return x
            else:
                x = self._driver(p)
                if not isinstance(x, bool):
                    self._theta = p
                    return x

        # store parameter vector for next iteration
        self._theta = p

        loglikelihood = 0.0
        for pulse in self._pulses:
            loglikelihood += pulse.loglikelihood

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
              physical_points=None):
        """ Perform checks on the likelihood evaluator and the prior
            density function.

        :param hypercube_points:
            :class:`numpy.ndarray` of ``n`` points in the unit hypercube, with
            shape ``(n, m)``, where ``m`` is the dimensionality of the sampling
            space -- i.e., of the hypercube.

        .. todo::
            Write a fallback routine porting directly from NumPy?

        """

        try:
            from _np import allclose
        except ImportError:
            yield 'Cannot import ``allclose`` function from NumPy...'
            yield 'Using fallback implementation...'

            def allclose(a, b, rtol, atol, equal_nan):
                return False

        lls = []
        lps = [] if logprior_call_vals is not None else None

        if physical_points is not None:
            for point in physical_points:
                lls.append(self.__call__(point))
                if lps is not None:
                    lps.append(self._prior(point))
        else:
            for point in hypercube_points:
                phys_point = self._prior.inverse_sample(point)
                lls.append(self.__call__(phys_point))
                if lps is not None:
                    lps.append(self._prior(phys_point))

        try:
            if allclose(_np.array(lls), loglikelihood_call_vals,
                        rtol=rtol_loglike, atol=atol_loglike,
                        equal_nan=False):
                yield 'Log-likelihood value checks passed on root process.'
            else:
                raise ValueError('Log-likelihood value checks failed on process '
                                 '%i' % _rank)
        except Exception as e:
            raise Exception('Log-likelihood value checks failed on process %i'
                            'with exception value: %s' % (_rank, e.value))

        if lps is not None:
            try:
                if allclose(_np.array(lps), logprior_call_vals,
                            rtol=rtol_logprior, atol=atol_logprior,
                            equal_nan=False):
                    yield 'Log-prior value checks passed on root process.'
                else:
                    raise ValueError('Log-prior value checks failed on '
                                     'process %i' % _rank)
            except Exception as e:
                raise Exception('Log-prior value checks failed on process %i '
                                'with exception value: %s' % (_rank, e.value))


    def synthesise(self, p, reinitialise=False, **kwargs):
        """ Synthesise pulsation data.

        :param list p: Parameter vector.

        :param optional[bool] reinitialise:
            Reinitialise? Useful if some resolution settings in child
            objects were changed (namely, the number of pulse phases) that
            need to be communicated to other child objects.

        """
        if reinitialise:
            self.__init__(self._star,
                          self._pulses,
                          self._threads,
                          self._llzero)

        p = _np.array(p)

        if (p != self._theta).any():
            try:
                logprior = self._prior(p)
            except AttributeError:
                pass
            else:
                if not _np.isfinite(logprior):
                    self._theta = p
                    return None

            if self._do_fast:
                # perform a low-resolution precomputation to direct cell
                # allocation
                x = self._driver(p, fast_mode=True)
                if not isinstance(x, bool):
                    self._theta = p
                    return None
                elif x:
                    x = self._driver(p, synthesise=True, **kwargs)
                    if not isinstance(x, bool):
                        self._theta = p
                        return None
            else:
                x = self._driver(p, synthesise=True, **kwargs)
                if not isinstance(x, bool):
                    self._theta = p
                    return None

        self._theta = p
