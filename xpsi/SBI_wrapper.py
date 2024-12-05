import torch
import numpy as np
import xpsi
from xpsi import Likelihood
from xpsi.tools.synthesise import synthesise_exposure
from xpsi.tools.synthesise import synthesise_given_total_count_number
from random import randint


class Custom_SBI_Likelihood(xpsi.Likelihood):

    """
    Custom likelihood function for use with SBI.

    Modifies the `_driver` and `synthesise` methods from the base class 
    to return `model_flux` that is the synthesised signal.
    """

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

            star_updated = True

        # register the signals by operating with the instrument response
        for signals, photosphere in zip(self._signals, self._star.photospheres):
            for signal in signals:
                if star_updated or signal.needs_update:
                    if signal.isI:
                        signal.register(tuple(
                                         tuple(self._divide(component,
                                                      self._star.spacetime.d_sq)
                                               for component in hot_region)
                                         for hot_region in photosphere.signal),
                                    fast_mode=fast_mode, threads=self.threads)
                    elif signal.isQ:
                        signal.register(tuple(
                                         tuple(self._divide(component,
                                                      self._star.spacetime.d_sq)
                                               for component in hot_region)
                                         for hot_region in photosphere.signalQ),
                                    fast_mode=fast_mode, threads=self.threads)
                    elif signal.isU:
                        signal.register(tuple(
                                         tuple(self._divide(component,
                                                      self._star.spacetime.d_sq)
                                               for component in hot_region)
                                         for hot_region in photosphere.signalU),
                                    fast_mode=fast_mode, threads=self.threads)
                    elif signal.isQn:
                        signal.register(tuple(
                                         tuple(self._divide(component,
                                                      self._star.spacetime.d_sq)
                                               for component in hot_region)
                                         for hot_region in photosphere.signalQ),
                                    fast_mode=fast_mode, threads=self.threads)
                        Qsignal = signal.signals
                        signal.register(tuple(
                                         tuple(self._divide(component,
                                                      self._star.spacetime.d_sq)
                                               for component in hot_region)
                                         for hot_region in photosphere.signal),
                                    fast_mode=fast_mode, threads=self.threads)
                        Isignal = signal.signals
                        for ihot in range(len(photosphere.signalQ)):
                            signal._signals[ihot]=np.where(Isignal[ihot]==0.0, 0.0, Qsignal[ihot]/Isignal[ihot])
                    elif signal.isUn:
                        signal.register(tuple(
                                         tuple(self._divide(component,
                                                      self._star.spacetime.d_sq)
                                               for component in hot_region)
                                         for hot_region in photosphere.signalU),
                                    fast_mode=fast_mode, threads=self.threads)
                        Usignal = signal.signals
                        signal.register(tuple(
                                         tuple(self._divide(component,
                                                      self._star.spacetime.d_sq)
                                               for component in hot_region)
                                         for hot_region in photosphere.signal),
                                    fast_mode=fast_mode, threads=self.threads)
                        Isignal = signal.signals
                        for ihot in range(len(photosphere.signalU)):
                            signal._signals[ihot]=np.where(Isignal[ihot]==0.0, 0.0, Usignal[ihot]/Isignal[ihot])
                    else:
                        raise TypeError('Signal type must be either I, Q, U, Qn, or Un.')
                    reregistered = True
                else:
                    reregistered = False

                if not fast_mode and reregistered:
                    if synthesise:
                        hot = photosphere.hot
                        try:
                            kws = kwargs.pop(signal.prefix)
                        except AttributeError:
                            kws = {}

                        shifts = [h['phase_shift'] for h in hot.objects]
                        signal.shifts = np.array(shifts)
                        model_flux = signal.synthesise(threads=self._threads, **kws)
                    else:
                        try:
                            hot = photosphere.hot
                            shifts = [h['phase_shift'] for h in hot.objects]
                            signal.shifts = np.array(shifts)
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

        return star_updated, model_flux

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

        :returns ndarray model_flux:
            2D array of synthesised counts.

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
            if not np.isfinite(logprior):
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
            x, model_flux = self._driver(fast_mode=True)
            if not isinstance(x, bool):
                super(Likelihood, self).__call__(self.cached) # restore
                return model_flux
            elif x:
                x, model_flux = self._driver(synthesise=True, **kwargs)
                if not isinstance(x, bool):
                    super(Likelihood, self).__call__(self.cached) # restore
                    return model_flux
        else:
            x, model_flux = self._driver(synthesise=True, **kwargs)
            if not isinstance(x, bool):
                super(Likelihood, self).__call__(self.cached) # restore
                return model_flux
            
            else:
                return model_flux

class SynthesiseData(xpsi.Data):
    """ Custom data container to enable synthesis.

    :param ndarray channels:
        Instrument channel numbers which must be equal in number to the first
        dimension of the count matrix.

    :param ndarray phases:
        Phases of the phase bins which must be equal in number to the second
        dimension of the count matrix.

    :param int first:
        First channel index number to include in the synthesised data.

    :param int last:
        Last channel index number to include in the synthesised data.
    """
    def __init__(self, channels, phases, first, last):
        self.channels = channels
        self._phases = phases
        try:
            self._first = int(first)
            self._last = int(last)
        except TypeError:
            raise TypeError('The first and last channels must be integers.')
        if self._first >= self._last:
            raise ValueError('The first channel number must be lower than the '
                             'the last channel number.')
        
def synthesise(self,
                exposure_time = None,
                expected_source_counts = None,
                nchannels = None,
                nphases = None,
                seed=0,
                **kwargs):
    
    """
    Synthesise data set.

    :param float exposure_time:
        Exposure time in seconds to scale the expected count rate.
    :param float expected_source_counts:
        Total expected number of source counts.
    :param int nchannels:
        Number of channels in the synthesised data.
    :param int nphases:
        Number of phase bins in the synthesised data.
    :param optional[int] seed:
        Seed for random number generation for Poisson noise in synthesised data.

    :return ndarray synthetic:
        The synthesised data set.
    """
    if nchannels is None or nphases is None:
        raise ValueError('nchannels and nphases must be specified.')
    bkg = np.zeros((nchannels, nphases))

    if exposure_time:
        _, synthetic, _= synthesise_exposure(exposure_time,
                                                self._data.phases,
                                                self._signals,
                                                self._phases,
                                                self._shifts,
                                                np.sum(bkg),
                                                bkg,
                                                gsl_seed=seed)
    elif expected_source_counts:
        _, synthetic, _, _ = synthesise_given_total_count_number(self._data.phases,
                                                expected_source_counts,
                                                self._signals,
                                                self._phases,
                                                self._shifts,
                                                np.sum(bkg),
                                                bkg,
                                                gsl_seed=seed)
        
    else:
        raise ValueError('Must specify either exposure time or expected source counts.')
    
    return synthetic

class xpsi_wrappers:
    """
    Class that wraps the xpsi likelihood and prior into a SBI compatible interface.

    :param xpsi.Prior prior:
        xpsi.Prior instance.
    :param xpsi.Likelihood likelihood:
        xpsi.Likelihood instance.
    :param dict instr_kwargs:
        Instrument keyword arguments for the likelihood synthesise method.
    :param optional[bool] train_using_CNNs:
        Whether to use CNNs for training. Defaults to True.
    """
    # Todo: Enable cuda parallelisation.
    def __init__(self, prior, likelihood, 
                 instr_kwargs={}, 
                 train_using_CNNs=True):
        self.prior = prior
        self.likelihood = likelihood
        self.instr_kwargs = instr_kwargs
        self.train_using_CNNs = train_using_CNNs
    
    def sample(self, sample_shape=torch.Size([])):
        """
        Sample from the prior distribution.

        :param torch.Size sample_shape:
            The shape of the sample. Defaults to torch.Size([]).
        :return optional[torch.Tensor] samples:
            The sampled values. Returns torch.cuda.Tensor is CUDA is available.
        """
        if len(sample_shape) > 0:
            samples = self.prior.draw(sample_shape[0])
            return torch.Tensor(samples[0]).cuda() if torch.cuda.is_available() else torch.Tensor(samples[0])
        else:
            return torch.Tensor(self.prior.draw(1)[0]).cuda() if torch.cuda.is_available() else torch.Tensor(self.prior.draw(1)[0])
    
    def log_prob(self, parameter_vector):
        """
        Compute the log probability of the parameter vector.

        :param torch.Tensor parameter_vector:
            The parameter vector.
        :return torch.Tensor log_probability:
            The log probability of the parameter vector.
        """
        parameter_vector = torch.Tensor(parameter_vector).cuda() if torch.cuda.is_available() else torch.Tensor(parameter_vector)
        log_probability = []
        for i in range(parameter_vector.shape[0]):
            log_probability.append(self.prior(parameter_vector[i,:]))
        return torch.Tensor(log_probability).cuda() if torch.cuda.is_available() else torch.Tensor(log_probability)
    
    def simulator(self, parameter_vector):
        """
        Compute the likelihood of the parameter vector.

        :param torch.Tensor parameter_vector:
            The parameter vector for which to simulate pulse profile.
        :return torch.Tensor model_flux:
            The pulse profile for the input parameter vector.
        """
        self.instr_kwargs['seed'] = randint(0,1000000000000000)
        parameter_vector = np.array(parameter_vector.cpu())
        model_flux = self.likelihood.synthesise(parameter_vector, force=True, 
                                                instr=self.instr_kwargs)

        if self.train_using_CNNs==True:
            model_flux = torch.Tensor(model_flux).cuda() if torch.cuda.is_available() else torch.Tensor(model_flux)
        else:
            model_flux = torch.flatten(torch.Tensor(model_flux)).cuda() if torch.cuda.is_available() else torch.flatten(torch.Tensor(model_flux))

        return model_flux
