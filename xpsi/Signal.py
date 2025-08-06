__all__ = ["Signal", "LikelihoodError"]

from xpsi.global_imports import *

from xpsi.Data import Data
from xpsi.Instrument import Instrument, ChannelError
from xpsi.Background import Background
from xpsi.Interstellar import Interstellar

from xpsi.tools.energy_integrator import energy_integrator
#from xpsi.tools.energy_interpolator import energy_interpolator
#from xpsi.tools.phase_integrator import phase_integrator

from abc import abstractmethod
from xpsi.Parameter import Parameter
from xpsi.ParameterSubspace import ParameterSubspace
from xpsi.tools.synthesise import synthesise_exposure as _synthesise_expo

from astropy.io import fits
from copy import deepcopy

class LikelihoodError(xpsiError):
    """ Raised if there is a problem with the value of the log-likelihood. """

class Signal(ParameterSubspace):
    """
    A signal is constituted by some X-ray dataset, a model instrument
    with which that data was acquired, a model background, and an object for
    modelling interstellar processes.

    The methods in this class must transform incident specific flux signals
    into a structure congruent to that of the data for the purpose of
    evaluation of the custom likelihood implemented via subclassing.

    :param obj data:
        An instance of :class:`~.Data.Data`.

    :param obj instrument:
        An instance of :class:`~.Instrument.Instrument`.

    :param obj background:
        If not ``None``, an instance of :class:`~.Background.Background`.
        It is assumed if one constructs a model using instances of
        :class:`~.Background.Background` that the background needs to be
        registered by a model instrument. If ``None``, it is still possible
        for one to define and use background parameters in a custom subclass
        of :class:`~.Signal`. In particular, background parameters for some
        model which directly specifies background contribution in units of
        count/s per *output* channels. These background parameters can even
        *be* the counts/s in output channels.

    :param obj interstellar:
        If not ``None``, an instance of :class:`~.Interstellar.Interstellar`.
        To be applied to the incident signal as a callable that modifies the
        signal in place.

    :param str photosphere_prefix:
        The ``str`` prefix of the photosphere object with which this signal
        object is associated.

    :param bool cache:
        Cache intermediary signals during likelihood evalation? When performing
        post-processing, this needs to be activated for full functionality of
        the :mod:`~.xpsi.PostProcessing` module. For likelihood function
        evaluation during sampling, caching should be deactivated because it is
        not used. It might be useful to activate caching also when preparing a
        model for a sampling application, to check the likelihood function
        works as intended.

    :param bool store:
        Deprecated. You can use this or ``cache``, which has the same effect.
        
    :param str stokes:
        Define the type of the signal. Options are Stokes "I" (default), "Q",
        "U", "Qn" (Q/I), and "Un" (U/I).

    """
    def __init__(self,
                 data,
                 instrument,
                 background = None,
                 interstellar = None,
                 support = None,
                 photosphere_prefix = None,
                 cache = False,
                 bounds = None,
                 values = None,
                 stokes = "I",
                 min_channel = 0,
                 max_channel = -1,
                 *args,
                 **kwargs):

        self._isI = False
        self._isQ = False
        self._isU = False
        self._isQn = False
        self._isUn = False        

        if stokes == "I":
            self._isI = True
        elif stokes == "Q":
            self._isQ = True
        elif stokes == "U":
            self._isU = True
        elif stokes == "Qn":
            self._isQn = True
        elif stokes == "Un":
            self._isUn = True            
            
        else:
             raise TypeError('param stokes for likelihood must be either "I", "Q", "U", "Qn", or "Un".')
             
        if not isinstance(data, Data):
            raise TypeError('Invalid type for a data object.')
        else:
            self._data = data
            self._original_data = deepcopy(data)

        if not isinstance(instrument, Instrument):
            raise TypeError('Invalid type for an instrument object.')
        else:
            self._instrument = instrument
            self._original_instrument = deepcopy( instrument )

        # Trimming the data and response so they fit together
        if min_channel != 0 or max_channel != -1:
            try:
                self._data.trim_data(min_channel, max_channel)
            except AttributeError:
                print('WARNING : There is no counts in data object. This is normal if you are trying to synthesise data.'
                      'Otherwise something is very wrong and do not continue')
            self._instrument.trim_response(min_channel, max_channel)

        a, b = data.index_range
        if (len(self._data.channels) != len(self._instrument.channels[a:b])):
            raise ChannelError( 'Size of the channel array declared for event data '
                                'does not match the size of the channel array declared'
                                ' for the loaded instrument response (sub)matrix.')

        if (self._data.channels != self._instrument.channels[a:b]).any():
            raise ChannelError('Channel array declared for event data does not '
                               'match channel array declared for the loaded '
                               'instrument response (sub)matrix.')

        # Check that they come from the same instrument
        if hasattr( self._data , 'instrument' ) and hasattr( self._instrument , 'name' ):
            assert self._data.instrument == self._instrument.name, 'Data and Instrument come from different instruments'

        self._identify_waveband()

        if background is not None:
            if not isinstance(background, Background):
                raise TypeError('Invalid type for a background object.')
            else:
                self._background = background
        else:
            self._background = None

        if support is not None:

            if self._data.counts.shape[0]==support.shape[0]:
                self._support = support
            else:
                raise TypeError("Data spectrum and background support must the have same shape")
        else:
            try :
                self._support = -1.0 * _np.ones((self._data.counts.shape[0],2))
                self._support[:,0] = 0.0
            except AttributeError:
                pass

        if interstellar is not None:
            if not isinstance(interstellar, Interstellar):
                raise TypeError('Invalid type for an interstellar object.')
            else:
                self._interstellar = interstellar
        else:
            self._interstellar = None

        if photosphere_prefix is not None:
            self._photosphere = photosphere_prefix

        cache = kwargs.get('store', cache)
        if not isinstance(cache, bool):
            raise TypeError('Activate or deactivate caching with a boolean.')
        self._cache = cache

        if bounds is None: bounds = {}
        if values is None: values = {}


        doc = """
        The phase shift for the signal, a periodic parameter [cycles].
        """
        phase_bounds = bounds.get('phase_shift', None)
        phase_value = values.get('phase_shift', 0.0 if phase_bounds is None else None)
        if phase_value is None:
            if not phase_bounds or None in phase_bounds:
                raise ValueError('Phase-shift bounds must be specified.')
            elif _np.array([not _np.isfinite(b) for b in phase_bounds]).any():
                raise ValueError('Phase-shift bounds must be finite.')
            elif not (0.0 <= (phase_bounds[1] - phase_bounds[0]) <= 1.0):
                raise ValueError('Phase bounds must be separated by '
                                 'a maximum of one cycle.')

        phase_shift = Parameter('phase_shift',
                                strict_bounds = (-_np.infty, _np.infty),
                                bounds = phase_bounds,
                                doc = doc,
                                symbol = r'$\phi$',
                                value = phase_value)

        # merge the subspaces; order unimportant
        super(Signal, self).__init__(self._instrument,
                                     self._background,
                                     self._interstellar,
                                     phase_shift,
                                     *args, **kwargs)


    @property
    def isI(self):
        """ ... """
        return self._isI
        
    @isI.setter
    def isI(self, b):
        """ ... """
        self._isI = b
        
    @property
    def isQ(self):
        """ ... """
        return self._isQ
        
    @isQ.setter
    def isQ(self, b):
        """ ... """
        self._isQ = b
        
    @property
    def isU(self):
        """ ... """
        return self._isU
        
    @isU.setter
    def isU(self, b):
        """ ... """
        self._isU = b 
        
    @property
    def isQn(self):
        """ ... """
        return self._isQn
        
    @isQn.setter
    def isQn(self, b):
        """ ... """
        self._isQn = b
        
    @property
    def isUn(self):
        """ ... """
        return self._isUn
        
    @isUn.setter
    def isUn(self, b):
        """ ... """
        self._isUn = b         
        
        
    @property
    def background(self):
        """ Get the instance of :class:`~.Background.Background`."""
        return self._background

    @property
    def interstellar(self):
        """ Get the instance of :class:`~.Interstellar.Interstellar`."""
        return self._interstellar

    @property
    def instrument(self):
        """ Get the instance of :class:`~.Instrument.Instrument`."""
        return self._instrument
    
    @property
    def original_data(self):
        """ Get the a copy of the original instance of :class:`~.Data.Data`."""
        return self._original_data

    @property
    def original_instrument(self):
        """ Get the a copy of the original instance of :class:`~.Instrument.Instrument`."""
        return self._original_instrument
    
    @property
    def photosphere(self):
        return self._photosphere

    def _identify_waveband(self):
            """ Bound the waveband for signal integration.

            Constructs an array of energy edges for instrument operation.
            This method thus automatically constructs energy bounds for this
            a particular instrument. At energies between these bounds signals
            are calculated. This requires details about the contiguous
            subset of output channels the photon data spans (in an instance of
            the :class:`~.Data.Data` class) and the redistribution matrix of the
            model instrument (in an instance of the
            :class:`~.Instrument.Instrument` class).

            :raises IndexError:
                If the channel range of the data object is not consistent with
                the instrument object.

            """
            a, b = self._data.index_range

            def search(i, j, k):
                while self._instrument.matrix[i,j] == 0.0:
                    j += k
                return j

            a = search(a, 0, 1)
            b = self._instrument.matrix.shape[1] + search(b-1, -1, -1) + 1

            self._input_interval_range = (a, b)
            self._energy_edges = self._instrument.energy_edges[a:b + 1]
            self._energy_mids = (self._energy_edges[:-1] + self._energy_edges[1:])/2.0

    @property
    def fast_energies(self):
        """ Get coarse array of energies for fast-mode likelihood evals. """
        return self._fast_energies

    @fast_energies.setter
    def fast_energies(self, energies):
        """ Set energies for fast mode."""
        self._fast_energies = energies

    def create_energy_array(self, rel_num_energies=10.0):
        """ Get a (finer) array of energies spanning instrument waveband.

        Useful for getting an appropriately bounded and spaced set of energies
        for signal interpolation.

        :param float rel_num_energies:
            The number of energies desired as a fraction of the number of
            energies implemented for incident signal integration.

        """
        L = self.energy_edges[0]
        R = self.energy_edges[-1]
        energies = _np.logspace(_np.log10(L), _np.log10(R),
                                int(rel_num_energies * len(self.energies)),
                                base=10.0)
        return energies

    @property
    def energy_edges(self):
        """ Get a :class:`numpy.ndarray` of energy edges. """
        return self._energy_edges

    def register(self, signals, fast_mode=False, threads=1):
        """  Register an incident signal by operating with the response matrix.

        A :class:`numpy.ndarray` is stored as an instance attribute containing
        source signal for each *output* channel in units of counts cm^2/s
        (assuming instrument effective area units are cm^2).

        """
        if fast_mode:
            try:
                del self.fast_total_counts
            except AttributeError:
                pass

            for hotRegion in signals:
                fast_total_counts = []

                for component, phases in zip(hotRegion, self.fast_phases):
                    if component is None:
                        fast_total_counts.append(None)
                    else:
                        integrated = energy_integrator(threads,
                                                       component,
                                                       _np.log10(self.fast_energies),
                                                       _np.log10(self._energy_edges))

                        # move interstellar to star?
                        if self._interstellar is not None:
                            self._interstellar(self._energy_mids, integrated)

                        temp = self._instrument(integrated,
                                                self._input_interval_range,
                                                self._data.index_range)

                        fast_total_counts.append(_np.sum(temp))

                self.fast_total_counts = tuple(fast_total_counts)
        else:
            try:
                del self.signals
            except AttributeError:
                pass

            if self.cache:
                try:
                    del self.incident_specific_flux_signals
                except AttributeError:
                    pass

                for hotRegion in signals: # iterate over hot regions
                    signal = None
                    for component in hotRegion: # add other components
                        try:
                            signal += component
                        except TypeError:
                            signal = component.copy()
                    # cache total hot region signal
                    self.incident_specific_flux_signals = signal

                try:
                    del self.incident_flux_signals
                except AttributeError:
                    pass

                try:
                    self.execute_custom_cache_instructions()
                except NotImplementedError:
                    pass # no custom caching targets

            for hotRegion in signals:
                integrated = None
                for component in hotRegion:
                    temp = energy_integrator(threads,
                                             component,
                                             _np.log10(self._energies),
                                             _np.log10(self._energy_edges))
                    try:
                        integrated += temp
                    except TypeError:
                        integrated = temp

                if self.cache:
                    self.incident_flux_signals = integrated.copy()

                if self._interstellar is not None:
                    self._interstellar(self._energy_mids, integrated)

                self.signals = self._instrument(integrated,
                                                self._input_interval_range,
                                                self._data.index_range)

            if self._background is not None:
                try:
                    self._background(self._energy_edges,
                                     self._data.phases)
                except TypeError:
                    print('Error when evaluating the incident background.')
                    raise

                self._background.registered_background = \
                                self._instrument(self._background.incident_background,
                                                 self._input_interval_range,
                                                 self._data.index_range)

    @property
    def num_components(self):
        return len(self._signals)

    @property
    def phases(self):
        return [phases.copy() for phases in self._phases]

    @phases.setter
    def phases(self, obj):
        if not isinstance(obj, list):
            obj = [obj]
        self._phases = obj

    @property
    def fast_phases(self):
        return [phases.copy() for phases in self._fast_phases]

    @fast_phases.setter
    def fast_phases(self, obj):
        if not isinstance(obj, list):
            obj = [obj]
        self._fast_phases = obj

    @property
    def energies(self):
        return self._energies

    @energies.setter
    def energies(self, obj):
            self._energies = obj

    @energies.deleter
    def energies(self):
        del self._energies

    @property
    def fast_total_counts(self):
        return tuple(self._fast_total_counts)

    @fast_total_counts.setter
    def fast_total_counts(self, obj):
        try:
            self._fast_total_counts.append(obj)
        except AttributeError:
            self._fast_total_counts = [obj]

    @fast_total_counts.deleter
    def fast_total_counts(self):
        del self._fast_total_counts

    @property
    def store(self):
        return self._cache

    @store.setter
    def store(self, value):
        if isinstance(value, bool):
            self._cache = value
        else:
            raise ValueError('Signal storage requires boolean activation.')

    @property
    def cache(self):
        return self._cache

    @cache.setter
    def cache(self, value):
        if isinstance(value, bool):
            self._cache = value
        else:
            raise ValueError('Signal storage requires boolean activation.')

    @property
    def data(self):
        """ Get the stored data object. """
        return self._data

    @data.setter
    def data(self, data):
        """ Set the data object. """
        if isinstance(data, Data):
                self._data = data
        else:
            raise TypeError('The data object is of an invalid type.')

    @property
    def signals(self):
        """ Get the stored channel-by-channel signal components. """
        return tuple(signal.copy() for signal in self._signals)

    @signals.setter
    def signals(self, obj):
        try:
            self._signals.append(obj)
        except AttributeError:
            self._signals = [obj]

    @signals.deleter
    def signals(self):
        del self._signals

    @property
    def incident_specific_flux_signals(self):
        """ Get the incident signal components. """
        return tuple(s.copy() for s in self._incident_specific_flux_signals)

    @incident_specific_flux_signals.setter
    def incident_specific_flux_signals(self, obj):
        try:
            self._incident_specific_flux_signals.append(obj)
        except AttributeError:
            self._incident_specific_flux_signals = [obj]

    @incident_specific_flux_signals.deleter
    def incident_specific_flux_signals(self):
        del self._incident_specific_flux_signals

    @property
    def incident_flux_signals(self):
        """ Get the incident flux signal components.

        These signals are integrated over a set of energy intervals spanning
        the instrument waveband.

        """
        return tuple(s.copy() for s in self._incident_flux_signals)

    @incident_flux_signals.setter
    def incident_flux_signals(self, obj):
        try:
            self._incident_flux_signals.append(obj)
        except AttributeError:
            self._incident_flux_signals = [obj]

    @incident_flux_signals.deleter
    def incident_flux_signals(self):
        del self._incident_flux_signals

    @property
    def expected_counts(self):
        return self._expected_counts

    @expected_counts.setter
    def expected_counts(self, obj):
        self._expected_counts = obj

    @expected_counts.deleter
    def expected_counts(self):
        del self._expected_counts

    @property
    def shifts(self):
        """ Returns the hot region phase plus the instrument phase-shift."""
        return self._shifts + self['phase_shift']

    @shifts.setter
    def shifts(self, obj):
        if isinstance(obj, _np.ndarray) and len(obj) == len(self._phases):
            self._shifts = obj
        else:
            raise TypeError('Store phase shift parameters as a 1D ndarray.')

    @shifts.deleter
    def shifts(self):
        del self._shifts

    @property
    def background_signal(self):
        """ Get stored background. """
        return self._background_signal

    @background_signal.setter
    def background_signal(self, obj):
        if isinstance(obj, _np.ndarray):
            self._background_signal = obj

    @background_signal.deleter
    def background_signal(self):
        del self._background_signal

    @property
    def caching_target_names(self):
        """ Just return the names of the caching targets. """
        return self._caching_targets

    @property
    def caching_targets(self):
        """ Get a dictionary of model objects for caching.

        Called by the post-processing module.

        :raises AttributeError:
            If a property is not set in methods of a subclass, or if the
            ``self.store`` property is not ``True``.

        """
        try:
            self._caching_targets
        except AttributeError:
            print('Caching targets not declared.')
            raise

        return {target: getattr(self, target) for target in self._caching_targets}

    @caching_targets.setter
    def caching_targets(self, obj):
        if isinstance(obj, list):
            if all(isinstance(o, _six.string_types) for o in obj):
                if all(hasattr(self, o) for o in obj):
                    self._caching_targets = obj
                    return None

        raise ValueError('Invalid caching targets.')

    def execute_custom_cache_instructions(self):
        """ Subclass and overwrite to specify custom cache objects.

        The default cached objects, when ``cache`` mode is activated, are
        handled in the :meth:`~.Signal.register` method.

        """
        raise NotImplementedError('Cache method not implemented.')

    @property
    def loglikelihood(self):
        """ Return the logarithm of the likelihood.

        :raises AttributeError: If property not set in methods of a subclass.

        """
        return self._loglikelihood

    @loglikelihood.setter
    def loglikelihood(self, ll):
        """ Check and store the logarithm of the likelihood. """

        if _np.isnan(ll):
            raise LikelihoodError('Log-likelihood is ``NaN``.')

        if not _np.isfinite(ll):
            self._loglikelihood = -_np.inf
        else:
            self._loglikelihood = ll

    @abstractmethod
    def __call__(self, **kwargs):
        """ Compute the logarithm of the likelihood and store it as a property.

        The keyword arguments currently communicated by an
        :class:`~.Likelihood.Likelihood` instance are as follows.

        :param int threads:
            Number of ``OpenMP`` threads to use for likelihood evaluation.
            This argument can be ignored if not required.

        :param float llzero:
            The minimum log-likelihood setting for MultiNest. Points whose
            log-likelihood is lower than this value are ignored.

        """

    def synthesise(self, phase_shifts, directory, **kwargs):
        """ Synthesise signal data according to the generative model.

        :param iterable phase_shifts:
            Container of phase shift :class:`~.Parameter.Parameter` instances,
            one per hot region, communicated by the likelihood object from
            the star object. The order is equal to the order of the hot
            region objects stored in ``photosphere.hot.objects``.

        :param str directory: Path to directory in which to write synthetic
                              data. It is recommended that the ``prefix`` of
                              the signal appears in the filename.
        :param int threads:
            Number of ``OpenMP`` threads to use for likelihood evaluation.
            This argument can be ignored if not required.

        """
        raise NotImplementedError('Cannot synthesise data.')
    
    @property
    def bolometric_chi2( self ):
        """ Return the bolometric chi2 value for the current likelihood value. """

        # Get the pulse of data and model
        pulse_data = self.data.counts.sum( axis = 0 )
        pulse_model = self.expected_counts.sum( axis = 0 )

        # Compute the chi squared
        chi2 = _np.sum( (pulse_data - pulse_model)**2 / pulse_model )

        return chi2
    
    ################ Sythesize data example and useful functions ################

    def synthesise_exposure(self,
                exposure_time,
                data_BKG,
                instrument_name,
                backscal_ratio=1.0,
                name='synthetic',
                directory='./',
                seed=None,
                save_background=True,
                **kwargs):
        
            """ Synthesise data set in FITS format that can be used in input of X-PSI.
            Data can be synthesised with or without background.
            
            :param float exposure_time: 
                Exposure time in seconds.

            :param array data_BKG: 
                Array of expected background counts in the background extraction region.
                It can be either integers or floats if needed.
            
            :param str instrument_name: 
                Instrument name to be written to the header. 
                Beware that it must be the same as the one in the background and response files.
                
            :param float backscal_ratio: 
                Ratio of the backscal parameters of the background to the source one. 
                This is used to rescale the background counts to be extracted from a 
                similar area as the source ones. 
                If the input background is considered from the same region as the source, set it to 1.
            
            :param str name: 
                Base name of synthetic data files.
                
            :param str directory: 
                Directory in which to write the synthetic data files.
                
            :param int seed: 
                Random seed for the pseudo-random number generator.
                
            :param bool save_background: 
                If True, the synthetic background is saved in the FITS file.
                If False, it is not saved.
                For NICER for instance, it does not need to be saved as it is unknown.
            """
            
            # Create or find directory
            try:
                if not _os.path.isdir(directory):
                    _os.mkdir(directory)
            except OSError:
                print('Cannot create write directory.')
                raise
    
            # Rescale background
            rescaled_background = data_BKG / backscal_ratio
            expected_background_counts = rescaled_background.sum()
            background_txt = rescaled_background.flatten()

            # Make the phase resolved background (constant for all phases but needed for _synthesise_expo)
            background_syn = _np.zeros((len(self._data.channels),len(self._data.phases)-1))      
            for j in range(len(self._data.phases)-1):
                background_syn[:,j] = background_txt/(len(self._data.phases)-1)
    
            # Synthesise
            self._expected_counts, synthetic, _  = _synthesise_expo(exposure_time,
                                                                self._data.phases,
                                                                self._signals,
                                                                self._phases,
                                                                self._shifts,
                                                                expected_background_counts,
                                                                background_syn,
                                                                gsl_seed=seed)
            
            # Write expected counts
            kwargs = {'channels':self._data.channels,
                    'counts':self._expected_counts,
                    'filename':_os.path.join(directory, name),
                    'instrument':instrument_name,
                    'exposure':exposure_time,
                    'backscale':1.0}
            self._write(**kwargs)

            # Write synthetic counts (Poisson realisation of expected counts)
            kwargs['counts'] = synthetic
            kwargs['filename'] += '_realization'
            self._write(**kwargs)

            # If requested, write background
            if save_background:

                # Write expected background
                kwargs = {'channels':self._data.channels,
                    'counts':data_BKG.flatten(),
                    'filename':_os.path.join(directory, name+'_bkg'),
                    'instrument':instrument_name,
                    'exposure':exposure_time,
                    'backscale':backscal_ratio}
                self._write(**kwargs)

                # Write synthetic background (Poisson realisation of expected background)
                kwargs['counts'] = _np.random.poisson( data_BKG.flatten() )
                kwargs['filename'] += '_realization'
                self._write(**kwargs)


    def _write(self, **kwargs):
        """ Wrapper to write either PHA or EVT file.
        This depends on whether the data is phase-resolved or not. """

        if len(self._data.phases) > 2:
            self._write_EVT(**kwargs)
        else:
            self._write_PHA(**kwargs)


    def _write_EVT(self, 
                channels,
                counts,
                filename,
                instrument,
                exposure, 
                backscale ):
        """ Write an event file in the FITS format. This is used for phase resolved data.
         
        :param array channels: 
            Array of channels for which data is written.
        
        :param array counts: 
            2D Array of counts, matching the channels on the 1st dimension, for which data is written.
        
        :param str filename: 
            Name of the FITS file.
        
        :param str instrument: 
            Instrument name to be written to the header. 
            Beware that it must be the same as the one in the background and response files. 
        
        :param float exposure: 
            Exposure time in seconds to be written to the header.
        
        :param float backscale: 
            Backscale parameter to be written to the header.
        """

        # Create the header
        header = {
            'INSTRUME': instrument, 
            'HDUCLAS1': 'EVENTS',
            'HDUCLAS2': 'none',
            'DATE': '2025-07-31',
            'EXPOSURE': exposure,
            'BACKSCAL': backscale,
            'ANCRFILE': 'none',
            'RESPFILE': 'none'
        }

        # Get the phase array
        phase_borders = _np.linspace( 0., 1., counts.shape[1]+1 )
        phase_centers = (phase_borders[1:] + phase_borders[:-1]) / 2
        PI, PULSE_PHASE = [], []

        # Fill the arrays
        for i in range(counts.shape[0]):
            for j in range(counts.shape[1]):
                for k in range( int(counts[i,j]) ):
                    PI.append( channels[i] )
                    PULSE_PHASE.append( phase_centers[j] )

        # Create the columns
        cols = [
            fits.Column(name='PI', array=_np.array(PI), format='1J'),
            fits.Column(name='PULSE_PHASE', array=_np.array(PULSE_PHASE), format='D'),
        ]

        # Create a primary HDU (header/data unit) for the header info
        primary_hdu = fits.PrimaryHDU()
        bin_table = fits.BinTableHDU.from_columns(cols, name='EVENTS')

        # Set header keywords
        for key, value in header.items():
            bin_table.header[key] = value

        # Create the FITS file (PHA file)
        hdul = fits.HDUList([primary_hdu, bin_table])

        # Save to a PHA file
        hdul.writeto(filename+'.evt', overwrite=True)

    def _write_PHA(self, 
                channels,
                counts,
                filename,
                instrument,
                exposure, 
                backscale ):
        """ Write a spectrum file in the FITS format. This is used for data with no phase information.
         
        :param array channels: 
            Array of channels for which data is written.
        
        :param array counts: 
            1D Array of counts, matching the channels, for which data is written.
        
        :param str filename: 
            Name of the FITS file.
        
        :param str instrument: 
            Instrument name to be written to the header. 
            Beware that it must be the same as the one in the background and response files. 
        
        :param float exposure: 
            Exposure time in seconds to be written to the header.
        
        :param float backscale: 
            Backscale parameter to be written to the header.
        """

        # Create the header
        header = {
            'INSTRUME': instrument, 
            'HDUCLAS1': 'SPECTRUM',
            'HDUCLAS2': 'none',
            'DATE': '2025-07-31',
            'EXPOSURE': exposure,
            'BACKSCAL': backscale,
            'ANCRFILE': 'none',
            'RESPFILE': 'none',
        }

        # Create the columns
        cols = [
            fits.Column(name='CHANNEL', array=channels, format='1J'),
            fits.Column(name='COUNTS', array=counts, format='D'),
        ]

        # Create a primary HDU (header/data unit) for the header info
        primary_hdu = fits.PrimaryHDU()
        bin_table = fits.BinTableHDU.from_columns(cols, name='SPECTRUM')

        # Set header keywords
        for key, value in header.items():
            bin_table.header[key] = value

        # Create the FITS file (PHA file)
        hdul = fits.HDUList([primary_hdu, bin_table])

        # Save to a PHA file
        hdul.writeto(filename+'.pha', overwrite=True)
            


def construct_energy_array(num_energies, signals, max_energy=None):
        """ Construct an array of photon energies for integration.

        :param int num_energies:
            Number of energies, distributed over union of wavebands covered
            by instruments that registered the data signals.

        :param list signals:
            An unordered list of :class:`~.Signal` instances.

        """
        ordered = [] # in waveband upper-limit, highest to lowest
        coverage_gaps = [] # highest to lowest

        # locate coverage gaps if any
        for _ in range(len(signals)):
            # find upper limit in energy from those remaining
            for signal in signals:
                try:
                    MAX
                except NameError:
                    MAX = signal.energy_edges[-1]
                    s = signal

                E = signal.energy_edges[-1]
                if E > MAX:
                    MAX = E
                    s = signal

            ordered.append(s)
            signals.remove(s)

            if len(ordered) > 1:
                for signal in ordered[:-1]:
                    try:
                        MIN
                    except NameError:
                        MIN = signal.energy_edges[0]

                    E = signal.energy_edges[0]
                    if E < MIN:
                        MIN = E

                if MAX < MIN: # MAX from above
                    coverage_gaps.append((MAX, MIN))

            del MAX

        # find global limits
        _signal_max = ordered[0].energy_edges[-1]
        if max_energy is not None and max_energy < _signal_max:
            MAX = max_energy

            # respect maximum energy setting
            _coverage_gaps = []
            for _coverage_gap in coverage_gaps:

                if _coverage_gap[0] < MAX <= _coverage_gap[1]:
                     MAX = _coverage_gap[0]

                if MAX > _coverage_gap[1]:
                    _coverage_gaps.append(_coverage_gap)

            coverage_gaps = _coverage_gaps
        else:
            MAX = _signal_max

        for signal in ordered:
            try:
                MIN
            except NameError:
                MIN = signal.energy_edges[0]

            E = signal.energy_edges[0]
            if E < MIN:
                MIN = E

        interval = _np.log10(MAX) - _np.log10(MIN)

        # account for gaps to conserve total number of energies requested
        for _coverage_gap in coverage_gaps:
            interval -= ( _coverage_gap[1] - _coverage_gap[0] )

        energies = _np.array([])

        # distribute energies over waveband intervals
        for i in range(len(coverage_gaps) + 1):
            if i == 0:
                U = MAX
            else:
                U = coverage_gaps[i-1][0]

            if i == len(coverage_gaps):
                L = MIN
            else:
                L = coverage_gaps[i][1]

            frac = ( _np.log10(U) - _np.log10(L) ) / interval
            num = int( _m.ceil(frac * num_energies) )

            energies = _np.append(energies,
                                  _np.logspace(_np.log10(L), _np.log10(U),
                                  int(num),
                                  base=10.0)[::-1])

        return _np.ascontiguousarray(energies[::-1])
