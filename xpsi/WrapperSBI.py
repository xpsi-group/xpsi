""" Wrappers exposing the X-PSI likelihood to the ``sbi`` package. """

__all__ = ["CustomLikelihoodSBI", "SynthesiseData", "SBIWrapper"]

from random import randint

import numpy as np
import torch

import xpsi
from xpsi import Likelihood
from xpsi.tools.synthesise import synthesise_exposure
from xpsi.tools.synthesise import synthesise_given_total_count_number


class CustomLikelihoodSBI(xpsi.Likelihood):
    """ Likelihood that returns synthesised signals for use with SBI.

    Overrides :meth:`~xpsi.Likelihood.Likelihood.synthesise` so that, instead
    of writing synthetic data to disk, it returns the synthesised counts
    (``model_flux``) that constitute the SBI training data. The base-class
    :meth:`~xpsi.Likelihood.Likelihood._driver` is reused unchanged; the
    synthesised signal is collected from each signal after driving.
    """

    def __init__( self, *args, **kwargs ):
        super(CustomLikelihoodSBI, self).__init__( *args, **kwargs )

        # Check that there is only one signal for now
        if len(self._signals) > 1 or len(self._signals[0]) > 1:
            raise ValueError('There is more than one signal instance. '
                             'This is currently not supported for SBI, '
                             'please use only one signal.')

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
            The synthesised counts, or ``None`` if the prior check failed.
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
                    print('Warning: Prior check failed, because at least one of '
                          'the parameters is not within the hard 1D-bounds. No '
                          'synthetic data will be produced.')
                else:
                    print('Warning: Prior check failed because a requirement set '
                          'in CustomPrior has failed. No synthetic data will be '
                          'produced.')
                # we need to restore due to premature return
                super(Likelihood, self).__call__(self.cached)
                return None

        # clear any previously synthesised signals so a failed drive cannot
        # return stale data
        self._signals[0][0].synthesised = None   # Only support one instrument for now

        x = self._driver(synthesise=True, force_update=force, **kwargs)
        if not isinstance(x, bool): # driver bailed out; restore cached values
            super(Likelihood, self).__call__(self.cached)

        model_flux = self._signals[0][0].synthesised
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
                             'last channel number.')


def synthesise_for_sbi(self,
                       exposure_time=None,
                       expected_source_counts=None,
                       nchannels=None,
                       nphases=None,
                       seed=0,
                       bkg=None, 
                       **kwargs):
    """ Synthesise a data set.

    Bound onto a :class:`~xpsi.Signal.Signal` subclass (see example notebooks),
    this replaces the file-writing :meth:`~xpsi.Signal.Signal.synthesise` with
    an in-memory version that returns the synthesised counts. The result is
    also stored on the signal as ``self.synthesised`` so the likelihood layer
    can retrieve it.

    :param float exposure_time:
        Exposure time in seconds to scale the expected count rate.
    :param float expected_source_counts:
        Total expected number of source counts.
    :param int nchannels:
        Number of channels in the synthesised data.
    :param int nphases:
        Number of phase bins in the synthesised data.
    :param optional[int] seed:
        Seed for the Poisson noise random number generation.
    :param optional[numpy.ndarray] bkg:
            Background for the synthesized pulse profile. Default set to zero 
            if none is provided. 

    :returns ndarray synthetic:
        The synthesised data set.
    """
    if nchannels is None or nphases is None:
        raise ValueError('nchannels and nphases must be specified.')

    # Set the background to zero if no background is provided 
    if bkg is None: 
        bkg = np.zeros((nchannels, nphases))

    if exposure_time:
        _, synthetic, _ = synthesise_exposure(exposure_time,
                                               self._data.phases,
                                               self._signals,
                                               self._phases,
                                               self._shifts,
                                               np.sum(bkg),
                                               bkg,
                                               gsl_seed=seed)
    elif expected_source_counts:
        _, synthetic, _, _ = synthesise_given_total_count_number(
                                               self._data.phases,
                                               expected_source_counts,
                                               self._signals,
                                               self._phases,
                                               self._shifts,
                                               np.sum(bkg),
                                               bkg,
                                               gsl_seed=seed)
    else:
        raise ValueError('Must specify either exposure time or expected '
                         'source counts.')

    self.synthesised = synthetic
    return synthetic


class SBIWrapper:
    """ Wrap the X-PSI likelihood and prior into an SBI-compatible interface.

    :param xpsi.Prior prior:
        xpsi.Prior instance.
    :param xpsi.Likelihood likelihood:
        xpsi.Likelihood instance.
    :param optional[dict] instr_kwargs:
        Instrument keyword arguments for the likelihood synthesise method.
    :param optional[bool] train_using_CNNs:
        Whether to use CNNs for training. Defaults to True.
    :param optional[str] device:
        Torch device on which to place the tensors. Defaults to ``"cuda"`` when
        available, otherwise ``"cpu"``.
    """
    # Todo: Enable cuda parallelisation.
    def __init__(self, prior, likelihood,
                 instr_kwargs=None,
                 train_using_CNNs=True,
                 device="cuda" if torch.cuda.is_available() else "cpu"):
        self.prior = prior
        self.likelihood = likelihood
        self.instr_kwargs = instr_kwargs if instr_kwargs is not None else {}
        self.train_using_CNNs = train_using_CNNs
        self.device = device

    def sample(self, sample_shape=torch.Size([])):
        """ Sample from the prior distribution.

        :param torch.Size sample_shape:
            The shape of the sample. Defaults to torch.Size([]).
        :returns torch.Tensor:
            The sampled values, on ``self.device``.
        """
        if len(sample_shape) > 0:
            samples = self.prior.draw(sample_shape[0])
            return torch.Tensor(samples[0]).to(self.device)
        return torch.Tensor(self.prior.draw(1)[0]).to(self.device)

    def log_prob(self, parameter_vector):
        """ Compute the log probability of the parameter vector.

        :param torch.Tensor parameter_vector:
            The parameter vector.
        :returns torch.Tensor:
            The log probability of each parameter vector.
        """
        parameter_vector = torch.Tensor(parameter_vector).to(self.device)
        log_probability = [self.prior(parameter_vector[i, :])
                           for i in range(parameter_vector.shape[0])]
        return torch.Tensor(log_probability).to(self.device)

    def simulator(self, parameter_vector):
        """ Simulate the pulse profile for the parameter vector.

        :param torch.Tensor parameter_vector:
            The parameter vector for which to simulate the pulse profile.
        :returns torch.Tensor model_flux:
            The pulse profile for the input parameter vector.
        """
        # The kwarg key passed to synthesise must match the signal's prefix
        # ('instr' in the example notebook), which routes it in Likelihood._driver.
        self.instr_kwargs['seed'] = randint(0, 1000000000000000)
        parameter_vector = np.array(parameter_vector.cpu())
        model_flux = self.likelihood.synthesise(parameter_vector, force=True,
                                                 instr=self.instr_kwargs)

        model_flux = torch.Tensor(model_flux)
        if not self.train_using_CNNs:
            model_flux = torch.flatten(model_flux)
        return model_flux.to(self.device)
