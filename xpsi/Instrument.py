__all__ = ["Instrument"]

from xpsi.global_imports import *

from xpsi.utils import make_verbose

from xpsi.ParameterSubspace import ParameterSubspace

from xpsi import Parameter
from xpsi.utils import make_verbose

from xpsi.PileupModule import XrayPileup

from astropy.io import fits
from astropy.table import Table

class ResponseError(xpsiError):
    """ Raised if there is a problem with the input response matrix. """

class EdgesError(xpsiError):
    """ Raised if there is a problem with the input energy edges. """

class ChannelError(xpsiError):
    """ Raised if there is a problem with the input channel numbers. """

class Instrument(ParameterSubspace):
    r""" Base class for astronomical X-ray instruments on-board space telescopes.

    The body of the initialiser must not be changed to ensure inter-module
    compatibility, but can be extended if appropriate using a call to
    ``super().__init__``. Specialist constructors can be defined in a subclass
    using the ``@classmethod`` decorator.

    :param ndarray[p,q] matrix:
        A :math:`p \\times q` matrix which is the product of a redistribution
        matrix and effective area vector. The input energy intervals must
        increase along the columns of :attr:`matrix`, and the output channels
        must increase along the rows of :attr:`matrix`. The *units* of the
        elements must be that of an *effective* area (:math:`cm^2`). Generally
        there will be some available calibration product, and deviations from
        this nominal response model will be parametrised. So here load some
        nominal response matrix.

    :param ndarray[q+1] energy_edges:
        Energy edges in keV of the instrument energy intervals which must be
        congruent to the first dimension of the :attr:`matrix`: the number of
        edges must be :math:`q + 1`. The edges must be monotonically increasing.
        These edges will correspond to the nominal response matrix and any
        deviation from this matrix (see above).

    :param ndarray[p] channels:
        Instrument channel numbers which must be equal in number to the number
        of rows of the :attr:`matrix`. The number of channels must therefore be
        :math:`p`. These channels will correspond to the nominal response
        matrix and any deviation from this matrix (see above). In common usage
        patterns, the channel numbers will increase monotonically with row
        number, and usually increment by one (but this is not necessary).

    .. note::

        That these channel numbers are not used to index the loaded instrument
        (sub)matrix. The :attr:`xpsi.Data.index_range` property returns
        bounding row numbers that index the loaded instrument response
        (sub)matrix in order to operate on an incident signal flux. The
        channel array contained in :attr:`xpsi.Data.channels` must be a
        contiguous (ordered) subset of the channel array loaded here.

    .. note::

        The dimensions of the response matrix need not be equal, but it is
        required that the number of input intervals be greater than or equal to
        the number of output channels -- i.e., :math:`p \leq q`. If :math:`p <
        q` then it is implied that subsets of adjacent output channels are
        effectively grouped together.

    :param ndarray[p+1] channel_edges:
        The channel (energy) edges of the instrument, in keV. The array must
        be congruent to the zeroth dimension of the :attr:`matrix`: the number
        of edges must be :math:`p + 1`. The edges must be monotonically
        increasing. These edges will correspond to the nominal response matrix
        and any deviation from this matrix (see above).

    :param tuple args:
        Container of parameter instances.

    :param dict kwargs:
        If you want to prefix parameters of an instance of this instrument
        subspace with an identifier, pass it as keyword argument and it will
        find its way to the base class.

    """
    def __init__(self, matrix, energy_edges, channels, channel_edges=None,
                 *args, **kwargs):

        self.matrix = matrix
        self.energy_edges = energy_edges
        self.channels = channels
        if channel_edges is not None:
            self.channel_edges = channel_edges

        super(Instrument, self).__init__(*args, **kwargs)

    @property
    def matrix(self):
        r""" Get the reference response matrix.

        In common usage patterns there will be some fiducial or nominal
        response matrix that either defines fixed instrument operation or
        is a basis for parametrised deviations. This matrix is usually a
        calibration product distributed by an instrument calibration team.

        A matrix of dimension :math:`p \\times q`. Here :math:`p` must be the
        number of input energy intervals, and :math:`q \geq p` the number of
        output channels.

        .. note::

            The attribute :attr:`matrix` must be assigned, and it must be
            a :class:`numpy.ndarray` for use with :func:`numpy.dot` (even
            if the matrix is sparse to some degree).

        """
        return self._matrix

    @matrix.setter
    def matrix(self, matrix):
        """ Set the matrix. """
        try:
            assert isinstance(matrix, _np.ndarray)
            assert matrix.ndim == 2
            assert (matrix >= 0.0).all()
        except AssertionError:
            raise ResponseError('Input matrix must be a two-dimensional ndarray with all matrix elements that are zero or positive.')
        try:
            for i in range(matrix.shape[0]):
                assert matrix[i,:].any()
        except AssertionError:
            raise ResponseError('Each row of the matrix must contain at least one positive number.')
        try:
            for j in range(matrix.shape[1]):
                assert matrix[:,j].any()
        except AssertionError:
            raise ResponseError('Each column of the matrix must contain at least one positive number.')
        self._matrix = matrix

    def construct_matrix(self):
        """ Construct the response matrix if it is parameterised.

        If customising, do operations to calculate a matrix, and return it.
        You can access parameters (free, fixed, and derived) via the container
        access ``self[<name>]``.

        If the instrument operation is fixed, you might not need to subclass,
        because the default behaviour is to return the nominal response you
        loaded. If for some reason the matrix you loaded is to be modified in
        some fixed manner, possibly as a function of some custom fixed
        parameters that you defined, you would also have to subclass and
        provide the correct implementation of this method.

        """
        return self.matrix

    def __call__(self, signal, irange, orange):
        """ Register an incident signal.

        :param ndarray[m,n] signal:
            An :math:`m \\times n` matrix, where input energy interval
            increments along rows, and phase increases along columns. The
            number of rows, :math:`m`, must equal the number of columns of
            :attr:`matrix`: :math:`m=q`.

        :param array-like irange:
            Indexable object with two elements respectively denoting the
            indices of the first and last *input* intervals. The response
            matrix :attr:`matrix` must be indexable with these numbers, i.e.,
            they must satisfy :math:`indx < q`.

        :param array-like orange:
            Indexable object with two elements respectively denoting the
            indices of the first and last *output* channels. The response
            matrix :attr:`matrix` must be indexable with these numbers, i.e.,
            they must satisfy :math:`indx < p`.

        :return: *ndarray[p,n]* containing the registered signal.

        .. note::

            The product of the most recent operation is stored as the property
            :attr:`cached_signal`.

        """
        matrix = self.construct_matrix()

        self._cached_signal = _np.dot(matrix[orange[0]:orange[1],
                                             irange[0]:irange[1]], signal)

        return self._cached_signal

    @property
    def cached_signal(self):
        """ Get the cached registered signal. """
        return self._cached_signal

    @property
    def energy_edges(self):
        """ Get the energy edges of the instrument, in keV.

        A :class:`numpy.ndarray` of edges of the input energy intervals which
        map to channels defined in the data space.

        """
        return self._energy_edges

    @energy_edges.setter
    def energy_edges(self, energy_edges):
        """ Set the energy edges in keV. """

        if not isinstance(energy_edges, _np.ndarray):
            try:
                energy_edges = _np.array(energy_edges)
            except TypeError:
                raise EdgesError('Energy edges must be in a one-dimensional array of positive increasing values.')

        try:
            assert energy_edges.ndim == 1
            assert (energy_edges >= 0.0).all()
            assert energy_edges.shape[0] == self._matrix.shape[1] + 1
            assert not (energy_edges[1:] <= energy_edges[:-1]).any()
        except AssertionError:
            raise EdgesError('Energy edges must be in a one-dimensional array of positive increasing values, with a '
                             'length equal to number of energy intervals in the matrix + 1.')

        self._energy_edges = energy_edges

    @property
    def channel_edges(self):
        """ Get the channel (energy) edges of the instrument, in keV.

        A :class:`numpy.ndarray` of edges of the registered energy intervals
        labelled as channels defined in the data space. This is relevant when
        there is a detector-by-detector gain scale applied to event data (such
        as for NICER instrument calibration products), meaning that the
        redistribution matrix is effectively shared by detectors and the
        channels across detectors can share an energy scale definition.

        An incident photon of given energy then has a registered-energy
        distribution that generally peaks in the vicinity of the true photon
        energy. The resdistribution matrix will have some energy resolution
        (along with other features such as shelves). With thanks to Paul S.
        Ray for explaining the choice to calibrate in this manner.

        .. note::

            If you made a channel cut that results in a non-contiguous subset
            of channels, you will need to overwrite the setter method because
            the checks will fail.

        """
        return self._channel_edges

    @channel_edges.setter
    def channel_edges(self, channel_edges):
        """ Set the channel (energy) edges in keV. """

        if not isinstance(channel_edges, _np.ndarray):
            try:
                channel_edges = _np.array(channel_edges)
            except TypeError:
                raise EdgesError('Channel edges must be in a one-dimensional array of positive increasing values.')

        try:
            assert channel_edges.ndim == 1
            assert (channel_edges >= 0.0).all()
            assert channel_edges.shape[0] == self._matrix.shape[0] + 1
            assert not (channel_edges[1:] <= channel_edges[:-1]).any()
        except AssertionError:
            raise EdgesError('Channel edges must be in a one-dimensional array of positive increasing values, with a '
                             'length equal to the number of channel intervals in the matrix + 1.')

        self._channel_edges = channel_edges

    @property
    def channels(self):
        """ Get the array of channels corresponding to rows of the matrix.

        The matrix being the loaded instrument response (sub)matrix.

        """
        return self._channels

    @channels.setter
    @make_verbose('Setting channels for loaded instrument response (sub)matrix',
                  'Channels set')
    def channels(self, channel_array):
        if not isinstance(channel_array, _np.ndarray):
            try:
                channel_array = _np.array(channel_array)
            except TypeError:
                raise ChannelError('Channel numbers must be in an array.')

        try:
            assert channel_array.ndim == 1
            assert (channel_array >= 0).all()
            assert channel_array.shape[0] == self._matrix.shape[0]
        except AssertionError:
            raise ChannelError('Channel numbers must be in a one-dimensional array of positive integers (including zero), with a '
                             'length equal to the number of channel in the matrix.')


        if (channel_array[1:] - channel_array[:-1] != 1).any():
            print('WARNING: Channel numbers do not uniformly increment by one.')

        self._channels = channel_array

        yield

    @make_verbose('Trimming instrument response',
                  'Instrument response trimmed')
    def trim_response(self, 
                      min_channel=0,
                      max_channel=-1,
                      tolerance=0.0 ):
        """ Trim the instrument response to the specified channel range.

        :param int min_channel:
            The minimum channel number to include in the trimmed response.

        :param int max_channel:
            The maximum channel number to include in the trimmed response.
            -1 will use the last channel.

        :param float tolerance:
            The tolerance value to use for trimming the instrument response.
            Trimming occurs at the first energy where every channel’s cumulative response is above the (1−tolerance) quantile.
            This allows to trim long tails of the response with little weight.
            If zero, no trimming is performed.

        """
        
        # Make the table of required channels
        assert min_channel >= self.channels[0]
        if max_channel == -1:
            max_channel = self.channels[-1]
        assert max_channel <= self.channels[-1]
        old_channels = self.channels
        new_channels_indexes = [ min_channel <= c <= max_channel for c in self.channels]

        # Adapt the matrix channel wise
        self.matrix = self.matrix[new_channels_indexes]
        self.channels = self.channels[new_channels_indexes]

        # Compute the cumulative of the response for all the channels
        cumsum = self.matrix.cumsum(axis=1)
        for i in range( cumsum.shape[0] ):
            cumsum[i] /= cumsum[i,-1]

        # Extract with the adapted percentile
        under_tolerance = [ _np.any( cumsum[:,i] <= 1 - tolerance ) for i in range( cumsum.shape[1] ) ]
        new_input_indexes = _np.where( under_tolerance )[-1]

        # Re-trim the response
        self.matrix = self.matrix[:,new_input_indexes]

        # Get the edges of energies for both input and channel
        new_energy_edges = [ self.energy_edges[k] for k in new_input_indexes ]
        self.energy_edges = _np.hstack( (new_energy_edges , self.energy_edges[ _np.where( self.energy_edges == new_energy_edges[-1] )[0] + 1 ] ) )
        if hasattr( self , 'channel_edges' ):
            new_channels_edges = [ self.channel_edges[ _np.where(old_channels==chan)[0][0]] for chan in self.channels]
            self.channel_edges = _np.hstack( (new_channels_edges , self.channel_edges[_np.where( old_channels==self.channels[-1])[0] + 1]) )

        # Print if any trimming happens
        if len(old_channels) > len(self.channels):
            print(f'Triming channels of the response matrix because a smaller channel range has been requested.\n '
                  f'Now min_channel={self.channels[0]} and max_channel={self.channels[-1]}')

        if len(under_tolerance) > len(new_input_indexes):
            print(f'Triming energy inputs of the response matrix with tolerance {tolerance}.\n '
                  f'Now min_energy={self.energy_edges[0]} and max_energy={self.energy_edges[-1]}')

        # If ARF and RMF, trim them
        if hasattr( self , 'ARF' ):
            self.RMF = self.RMF[new_channels_indexes][:,new_input_indexes]
            self.ARF = self.ARF[new_input_indexes]

    @classmethod
    @make_verbose('Loading instrument response matrix from OGIP compliant files',
                  'Response matrix loaded')
    def from_ogip_fits(cls,
              RMF_path,
              ARF_path=None,
              min_channel=0,
              max_channel=-1,
              min_input=1,
              max_input=-1,
              datafolder=None,
              **kwargs):
        """ Loading method for Instrument using OGIP defined ARF/RMF or RSP.

        :param str RMF_path:
            The path to the RMF file which should be OGIP compliant. Path to the OGIP compliant RSP file it ARF_path is None. 

        :param str | None ARF_path:
            The path to the ARF file which should be OGIP compliant or None if the RMF_path points to a RSP file. 

        :param int min_channel:
            The minimum channel for which the instrument response is loaded.

        :param int max_channel:
            The maximum channel for which the instrument response is loaded.
            -1 will use the last channel.

        :param int min_input:
            The minimum input energy number for which the instrument response is loaded.

        :param int max_input:
            The maximum input energy number for which the instrument response is loaded.
            -1 will use the last input.

        :param str | None datafolder:
            The path to the folder which contains both ARF and RMF files, if not specified in RMF_path or ARF_path.

        :return:
            Instrument instance with loaded instrument response matrix.
        """

        if datafolder:
            ARF_path = _os.path.join( datafolder, ARF_path ) if ARF_path is not None else None
            RMF_path = _os.path.join( datafolder, RMF_path )

        # Open useful values in ARF/RMF/RSP
        with fits.open( RMF_path ) as RMF_hdul:
            RMF_header = RMF_hdul['MATRIX'].header
        RMF_instr = RMF_header['INSTRUME'] 
        DETCHANS = RMF_header['DETCHANS']
        NUMGRP = RMF_header['NAXIS2']
        TLMIN = RMF_header['TLMIN4']
        TLMAX = RMF_header['TLMAX4']

        # Handle the RSP case
        if ARF_path is not None:
            with fits.open( ARF_path ) as ARF_hdul:
                assert RMF_instr == ARF_hdul['SPECRESP'].header['INSTRUME']

        # Get the values and change the -1 values if requried
        if min_channel == 0:
            min_channel = TLMIN
        if max_channel == -1:
            max_channel = DETCHANS -1
        if max_input == -1:
            max_input = NUMGRP

        # Perform routine checks
        assert min_channel >= TLMIN and max_channel <= TLMAX
        assert min_input >= 0 and max_input <= NUMGRP

        # If everything in order, get the data
        with fits.open( RMF_path ) as RMF_hdul:
            RMF_MATRIX = RMF_hdul['MATRIX'].data
            RMF_EBOUNDS = RMF_hdul['EBOUNDS'].data

        # Get all the channels and input energies
        channels = _np.array( RMF_EBOUNDS['CHANNEL'] )
        inputs = _np.arange( 1, NUMGRP+1 )

        # Get the channels from the data
        RMF = _np.zeros((DETCHANS, NUMGRP))
        for i, (N_GRP, F_CHAN, N_CHAN, RMF_line) in enumerate( zip(RMF_MATRIX['N_GRP'], RMF_MATRIX['F_CHAN'], RMF_MATRIX['N_CHAN'], RMF_MATRIX['MATRIX']) ):

            # Skip if needed
            if N_GRP == 0:
                continue

            # Check the values
            if not isinstance(F_CHAN, _np.ndarray ):
                F_CHAN = [F_CHAN]
                N_CHAN = [N_CHAN]

            # Add the values to the RMF
            n_skip = 0 
            for f_chan, n_chan in zip(F_CHAN,N_CHAN):

                if n_chan == 0:
                    continue

                RMF[f_chan-TLMIN:f_chan+n_chan-TLMIN,i] += RMF_line[n_skip:n_skip+n_chan]
                n_skip += n_chan

        # Get the indexes of channels
        channels_indexes = _np.where( ( channels >= min_channel ) & ( channels <= max_channel ) )[0]
        inputs_indexes = _np.where( ( inputs >= min_input ) & ( inputs <= max_input ) )[0]
        channels = channels[channels_indexes]
        inputs = inputs[inputs_indexes]
        RMF = RMF[channels_indexes][:,inputs_indexes]

        # Make the RSP, depending on the input files
        if ARF_path is None:
            RSP = RMF
            ARF_area = RMF.sum( axis=0 )
        else:
            ARF = Table.read(ARF_path, 'SPECRESP')
            ARF_area = ARF['SPECRESP'][inputs_indexes]
            RSP = RMF * ARF_area

        # Find empty columns and lines
        empty_channels = _np.all(RSP == 0, axis=1)
        empty_inputs = _np.all(RSP == 0, axis=0)
        RSP = RSP[~empty_channels][:,~empty_inputs]
        channels = channels[ ~empty_channels ]
        inputs = inputs[ ~empty_inputs ]

        # Get the edges of energies for both input and channel
        energy_edges = _np.append( RMF_MATRIX['ENERG_LO'][inputs-1], RMF_MATRIX['ENERG_HI'][inputs[-1]-1]).astype(dtype=_np.double)
        energies = (RMF_MATRIX['ENERG_LO']+RMF_MATRIX['ENERG_HI'])/2
        channel_energy_edges = _np.append(RMF_EBOUNDS['E_MIN'][channels-TLMIN],RMF_EBOUNDS['E_MAX'][channels[-1]-TLMIN])

        # Print informations
        if empty_inputs.sum() > 0:
            print(f'Triming the response matrix because it contains rows with only 0 values.\n '
                  f'Now min_energy={energy_edges[0]} and max_energy={energy_edges[-1]}')
        if empty_channels.sum() > 0:
            print(f'Triming the response matrix because it contains columns with only 0 values.\n'
                  f' Now min_channel={channels[0]} and max_channel={channels[-1]}')

        # Make the instrument
        Instrument = cls(RSP,
                         energy_edges,
                         channels,
                         channel_energy_edges,
                         **kwargs)
        
        # Add ARF and RMF for plotting
        Instrument.RMF = RMF[~empty_channels][:,~empty_inputs]
        Instrument.ARF = ARF_area[~empty_inputs]
        Instrument.name = RMF_instr
        Instrument.energies = energies[inputs_indexes][~empty_inputs]

        return Instrument


class InstrumentPileup(Instrument):

    """ Instrument with pile-up computation """

    def __call__(self, signal, *args):
        """ Overwrite. """
        ## Compute spectrum with pileup
        piled_spectrum = self.pileup.analyze(signal,
                                             alpha=self['grade_migration'],
                                             psf_frac=self['psf_fraction'])

        self._cached_signal = piled_spectrum
        return self._cached_signal

    @classmethod
    @make_verbose('Loading instrument response matrix',
                  'Response matrix loaded')
    def from_ogip_fits(cls,
              Data_path,
              bounds=dict(),
              values=dict(),
              **kwargs):
        
        """ Load any instrument response matrix. """
    
        # Load the default instrument class
        Instrument = super().from_ogip_fits(**kwargs)

        ## -------- INITIALIZATION OF THE PILEUP ---------------

        ## Definition of the parameters for the pile-up
        alpha_grade = Parameter('grade_migration',
                    strict_bounds = (0.0,1.0),
                    bounds = bounds.get('grade_migration', None),
                    doc = 'Grade migration factor : probability that the piled event is not rejected as "bad event"',
                    symbol = r'$G_n$',
                    value = values.get('grade_migration', None))
    
        psffrac = Parameter('psf_fraction',
                    strict_bounds = (0.8,1.0),  #min in sherpa set at 0.85
                    bounds = bounds.get('psf_fraction', None),
                    doc = 'fraction of events in the source extraction region to which pileup will be applied',
                    symbol = r'$PSF_frac',
                    value = values.get('psf_fraction', 0.9))
    
        nregions = Parameter('nregions',
                    strict_bounds = (0.0,10.0),
                    bounds = None,     ##this parameter should always be fixed  - should be 1.0 for point sources
                    doc = 'number of regions to which pileup model will be applied independently',
                    symbol = r'$N_{regions}$',
                    value = values.get('nregions', 1.0))

        g0 = Parameter('g0',
                    strict_bounds = (0.0,1.0),
                    bounds = None,     ##this parameter should always be fixed
                    doc = 'grade correction for single photon detection',
                    symbol = r'$g_0$',
                    value = values.get('g0', 1.0))
    
        npiled = Parameter('npiled',
                    strict_bounds = (0.0,100.0),
                    bounds = None,     ##this parameter should always be fixed 
                    doc = 'number of photons considered for pileup in a single frame',
                    symbol = r'$N_{phot}$',
                    value = values.get('npiled', 5)) #or 30 in sherpa

        # Open the headers for values
        with fits.open( Data_path ) as hdul:
            Data_header = hdul['SPECTRUM'].header 

        if kwargs.get( 'datafolder', None ) is not None:
            ARF_path = _os.path.join( kwargs.get( 'datafolder'), kwargs.get( 'ARF_path' ) )
        else:
            ARF_path = kwargs.get( 'ARF_path' )
        with fits.open( ARF_path ) as hdul:
            ARF_header = hdul['SPECRESP'].header 

        ## make sure that EXPTIME and FRACEXPO exist - should be the case if it's Chandra data
        ## EXPTIME is in the _evt2.fits file, might need to be manually copied to the data file
        frametime = Data_header['EXPTIME']
        frac_expo = ARF_header['FRACEXPO']

        frame = Parameter('frame_time',
                    strict_bounds = (0.0,10.0),
                    bounds = None,    
                    doc = 'good exposure time per frame',
                    symbol = r'$\tau$ (s)',
                    value = frametime )

        fracexpo = Parameter('frac_expo',
                    strict_bounds = (0.0,1.0),
                    bounds = None,    
                    doc = 'fraction of the frame exposure time to create effective frame exposure time',
                    symbol = r'$f_{expo}$',
                    value = frac_expo) 

        # Merge the parameters
        Instrument.merge(alpha_grade,psffrac,nregions,g0,npiled,frame,fracexpo)

        ##Initialization of the pileup module
        pileup = XrayPileup(Instrument)
        Instrument.pileup = pileup

        return Instrument
