from xpsi.global_imports import *

from xpsi.Spacetime import Spacetime
from xpsi.HotRegion import HotRegion
from xpsi.Elsewhere import Elsewhere
from xpsi.Everywhere import Everywhere

from xpsi.Parameter import Parameter
from xpsi.ParameterSubspace import ParameterSubspace

class Photosphere(ParameterSubspace):
    """ A photosphere embedded in an ambient Schwarzschild spacetime.

    :param obj hot:
        An instance of :class:`~.HotRegion.HotRegion` (or a derived class).
        This objects represents the hot regions of the surface that in most
        use-cases will be assumed to contain radiating material that is hotter
        than that *elsewhere*.

    :param obj elsewhere:
        An instance of :class:`~.Elsewhere.Elsewhere` (or a derived class).

    :param obj everywhere:
        An instance of :class:`~.Everywhere.Everywhere` (or a derived class).

    .. note::

        You cannot specify the surface radiation field *everywhere* if you
        use hot regions (the latter usage may also include specification of
        the radiation field *elsewhere*).

    :param dict bounds:
        Bounds are supplied for instantiation of a frequency parameter.
        The parameter name ``'mode_frequency'`` must be a key in the
        dictionary unless the parameter is *fixed* or *derived*. If a bound
        is ``None`` that bound is set equal to a strict hard-coded bound.
        If ``None``, lock the coordinate rotation frequency of a mode of
        asymmetry in the photosphere to a fixed frequency, e.g., the stellar
        rotation frequency. If bounds are passed, the frequency is interpreted
        as a free parameter.

    :param dict values:
        Either the fixed value of the mode frequency, a callable if the
        frequency is *derived*, or a value upon initialisation if the
        frequency is free. The dictionary must have a key with name
        ``'mode_frequency'`` if it is *fixed* or *derived*.
        If the asymmetry is locked to the stellar spin, then you need to pass
        the spin frequency. If fixed but different to the spin frequency, this
        value needs to be passed instead. In the hot region base class this
        mode frequency is applied to normalise the ray lags instead of the
        stellar rotation frequency.

    :param boolean stokes:
        A Boolean that determines whether the signals for all the Stokes I, Q,
        and U parameters are calculated. If False, only Stokes I is calculated.

    :param iterable custom:
        A :class:`~.Parameter.Parameter` instance or iterable over such
        instances. Might be useful for calling image plane extensions and
        passing global variables, without having to instantiate
        surface-discretisation classes and without having to handle global
        variable values at compile time or from disk for runtime access.

    .. note::

        In basic modelling patterns the frequency is the spin frequency,
        and thus you only need to explicitly pass the spin as ``value`` whilst
        leaving ``bounds`` to default. If the spin frequency happens to be a
        free parameter (perhaps with informative prior information), then
        pass a callable instead that can be used to get the spin frequency
        dynamically when the derived mode frequency variable is called for.

    """
    required_names = ['mode_frequency']

    def __init__(self,
                 hot = None, elsewhere = None,
                 everywhere = None,
                 bounds = None, values = None,
                 stokes=False,
                 custom = None,
                 **kwargs):

        if everywhere is not None:
            if hot or elsewhere is not None:
                raise ValueError('Cannot use hot region nor elsewhere '
                                 'functionality if constructing the '
                                 'radiation field everywhere.')
            if not isinstance(everywhere, Everywhere):
                raise TypeError('Invalid type for everywhere object.')
            self._everywhere_atmosphere = ()
            
        elif hot is None and elsewhere is None:
            pass # can call image-plane extensions

        else:
            if elsewhere is not None:
                if not isinstance(elsewhere, Elsewhere):
                    raise TypeError('Invalid type for an elsewhere object.')

                if hot is None:
                    raise ValueError('Hot region object(s) must be used in '
                                     'conjuction with an elsewhere object.')

            self._elsewhere_atmosphere = ()
                                              # including derived classes
            if hot is not None and hot is not isinstance(hot, HotRegion):
                if hasattr(hot, 'objects'):
                    for obj in getattr(hot, 'objects'):
                        if not isinstance(obj, HotRegion):
                            raise TypeError('Invalid object for the hot '
                                            'region(s).')
                else:
                    raise TypeError('Invalid object for the hot region(s).')

        self._hot = hot
        self._hot_atmosphere = ()
        self._hot_atmosphere_Q = ()
        self._elsewhere = elsewhere
        self._everywhere = everywhere
        self._stokes = stokes
        if hot is not None:
            self._surface = self._hot
        else:
            self._surface = self._everywhere
            self.surface.objects = [self.surface]

        if bounds is None: bounds = {}
        if values is None: values = {}

        doc = """
        Coordinate frequency of the mode of radiative asymmetry in the
        photosphere that is assumed to generate the pulsed signal [Hz].
        """
        mode_frequency = Parameter('mode_frequency',
                                   strict_bounds = (0.0, 2000.0),
                                   bounds = bounds.get('mode_frequency', None),
                                   doc = doc,
                                   symbol = r'$f_{\rm mode}$',
                                   value = values.get('mode_frequency', None))

        if stokes:
            doc = """
            Spin axis position angle measured from the north counterclock-
            wise to the projection of the rotation axis on the plane of the
            sky [in radians].
            """
            spin_axis_position_angle = Parameter('spin_axis_position_angle',
                                       strict_bounds = (-_np.pi/2.0, _np.pi/2.0),
                                       bounds = bounds.get('spin_axis_position_angle', None),
                                       doc = doc,
                                       symbol = r'$\chi_{0}$',
                                       value = values.get('spin_axis_position_angle', None))

            super(Photosphere, self).__init__(mode_frequency, spin_axis_position_angle,
                                              hot, elsewhere, everywhere,
                                              custom,
                                              **kwargs)
        else:
            super(Photosphere, self).__init__(mode_frequency,
                                              hot, elsewhere, everywhere,
                                              custom,
                                              **kwargs)

    def load_NSX_table( self, path, Tcol=4, gcol=5, mucol=1, Ecol=0, spe_Icol=2 ):
        """
        Loading the nsx atmosphere table provided in path 
        giving the colums in the table corresponding to 
        - the logarithm of local comoving effective temperature logTeff(K) (Tcol, default=4), 
        - the logarithm of effective surface gravity logg(cm s^-2) (gcol, default=5), 
        - the cosine of the angle from the local surface normal mu = cos(theta) (mucol, default=1), 
        - the logarithm of the photon energy log(E/kTeff) (Ecol, default=0) 
        - the one-dimensional buffer of specific intensity log(Inu/Teff^3) (spe_Icol, default=2)
        """

        # Load tables and get sizes
        table = _np.loadtxt(path, dtype=_np.double)
        lenlogT = len( _np.unique(table[:,Tcol]) )
        lenlogg = len( _np.unique(table[:,gcol]) )
        lenmu = len( _np.unique(table[:,mucol]) )
        lenlogE = len( _np.unique(table[:,Ecol]) )

        # Make respective tables
        logT = _np.zeros( lenlogT )
        logg = _np.zeros( lenlogg )
        mu = _np.zeros( lenmu )
        logE = _np.zeros( lenlogE )

        reorder_buf = _np.zeros((lenlogT,
                                lenlogg,
                                lenmu,
                                lenlogE,))

        index = 0
        for i in range(lenlogT):
            for j in range(lenlogg):
                for k in range(lenlogE):
                    for l in range(lenmu):
                        logT[i] = table[index,Tcol]
                        logg[j] = table[index,gcol]
                        logE[k] = table[index,Ecol]
                        mu[reorder_buf.shape[2] - l - 1] = table[index,mucol]
                        reorder_buf[i,j,reorder_buf.shape[2] - l - 1,k] = 10.0**(table[index,spe_Icol])
                        index += 1

        buf = _np.zeros(_np.prod(reorder_buf.shape))

        bufdex = 0
        for i in range(lenlogT):
                for j in range(lenlogg):
                    for k in range(lenmu):
                        for l in range(lenlogE):
                            buf[bufdex] = reorder_buf[i,j,k,l]; bufdex += 1

        return logT, logg, mu, logE, buf

    @property
    def hot_atmosphere(self):
        """ Get the numerical atmosphere buffers for hot regions if used.

        To preload a numerical atmosphere into a buffer, subclass and
        overwrite the setter. The underscore attribute set by the setter
        must be an :math:`n`-tuple whose :math:`n^{th}` element is an
        :math:`(n-1)`-dimensional array flattened into a one-dimensional
        :class:`numpy.ndarray`. The first :math:`n-1`
        elements of the :math:`n`-tuple must each be an ordered one-dimensional
        :class:`numpy.ndarray` of parameter values for the purpose of
        multi-dimensional interpolation in the :math:`n^{th}` buffer. The
        first :math:`n-1` elements must be ordered to match the index
        arithmetic applied to the :math:`n^{th}` buffer. An example would be
        ``self._hot_atmosphere = (logT, logg, mu, logE, buf)``, where:
        ``logT`` is a logarithm of local comoving effective temperature;
        ``logg`` is a logarithm of effective surface gravity;
        ``mu`` is the cosine of the angle from the local surface normal;
        ``logE`` is a logarithm of the photon energy; and
        ``buf`` is a one-dimensional buffer of intensities of size given by
        the product of sizes of the first :math:`n-1` tuple elements.

        It is highly recommended that buffer preloading is used, instead
        of loading from disk in the customisable radiation field extension
        module, to avoid reading from disk for every signal
        (likelihood) evaluation. This can be a non-negligible waste of compute
        resources. By preloading in Python, the memory is allocated and
        references to that memory are not in general deleted until a sampling
        script exits and the kernel stops. The likelihood callback accesses
        the same memory upon each call without I/O.

        """
        return self._hot_atmosphere

    @hot_atmosphere.setter
    def hot_atmosphere(self, path,Tcol=3, gcol=4, mucol=1, Ecol=0, spe_Icol=2):

        if 'nsx' in path:
            # Read and set attributes of NSX model table
            logT, logg, mu, logE, buf = self.load_NSX_table( path ,Tcol, gcol, mucol, Ecol, spe_Icol)
            self._hot_atmosphere = (logT, logg, mu, logE, buf)

        else:
            ## if you want to set another model
            """ Implement if required. """
            raise NotImplementedError('Implement setter if required.')


    @property
    def hot_atmosphere_Q(self):
        """ Get the numerical atmosphere Stokes Q buffers for hot regions if used.

        To preload a numerical atmosphere into a buffer, subclass and
        overwrite the setter. The underscore attribute set by the setter
        must be an :math:`n`-tuple whose :math:`n^{th}` element is an
        :math:`(n-1)`-dimensional array flattened into a one-dimensional
        :class:`numpy.ndarray`. The first :math:`n-1`
        elements of the :math:`n`-tuple must each be an ordered one-dimensional
        :class:`numpy.ndarray` of parameter values for the purpose of
        multi-dimensional interpolation in the :math:`n^{th}` buffer. The
        first :math:`n-1` elements must be ordered to match the index
        arithmetic applied to the :math:`n^{th}` buffer. An example would be
        ``self._hot_atmosphere_Q = (logT, logg, mu, logE, buf)``, where:
        ``logT`` is a logarithm of local comoving effective temperature;
        ``logg`` is a logarithm of effective surface gravity;
        ``mu`` is the cosine of the angle from the local surface normal;
        ``logE`` is a logarithm of the photon energy; and
        ``buf`` is a one-dimensional buffer of intensities of size given by
        the product of sizes of the first :math:`n-1` tuple elements.

        It is highly recommended that buffer preloading is used, instead
        of loading from disk in the customisable radiation field extension
        module, to avoid reading from disk for every signal
        (likelihood) evaluation. This can be a non-negligible waste of compute
        resources. By preloading in Python, the memory is allocated and
        references to that memory are not in general deleted until a sampling
        script exits and the kernel stops. The likelihood callback accesses
        the same memory upon each call without I/O.

        """
        return self._hot_atmosphere_Q

    @hot_atmosphere_Q.setter
    def hot_atmosphere_Q(self, path, Tcol=3, gcol=4, mucol=1, Ecol=0, spe_Icol=2):

        if 'nsx' in path:
            # Read and set attributes of NSX model table
            logT, logg, mu, logE, buf = self.load_NSX_table( path ,Tcol, gcol, mucol, Ecol, spe_Icol)
            self._hot_atmosphere_Q = (logT, logg, mu, logE, buf)

        else:
            """ Implement if required. """
            raise NotImplementedError('Implement setter if required.')


    @property
    def elsewhere_atmosphere(self):
        """ Get the numerical atmosphere buffers for elsewhere if used.

        To preload a numerical atmosphere into a buffer, subclass and
        overwrite the setter. The underscore attribute set by the setter
        must be an :math:`n`-tuple whose :math:`n^{th}` element is an
        :math:`(n-1)`-dimensional array flattened into a one-dimensional
        :class:`numpy.ndarray`. The first :math:`n-1`
        elements of the :math:`n`-tuple must each be an ordered one-dimensional
        :class:`numpy.ndarray` of parameter values for the purpose of
        multi-dimensional interpolation in the :math:`n^{th}` buffer. The
        first :math:`n-1` elements must be ordered to match the index
        arithmetic applied to the :math:`n^{th}` buffer. An example would be
        ``self._hot_atmosphere = (logT, logg, mu, logE, buf)``, where:
        ``logT`` is a logarithm of local comoving effective temperature;
        ``logg`` is a logarithm of effective surface gravity;
        ``mu`` is the cosine of the angle from the local surface normal;
        ``logE`` is a logarithm of the photon energy; and
        ``buf`` is a one-dimensional buffer of intensities of size given by
        the product of sizes of the first :math:`n-1` tuple elements.

        It is highly recommended that buffer preloading is used, instead
        of loading from disk in the customisable radiation field extension
        module, to avoid reading from disk for every signal
        (likelihood) evaluation. This can be a non-negligible waste of compute
        resources. By preloading in Python, the memory is allocated and
        references to that memory are not in general deleted until a sampling
        script exits and the kernel stops. The likelihood callback accesses
        the same memory upon each call without I/O.

        """
        return self._elsewhere_atmosphere

    @elsewhere_atmosphere.setter
    def elsewhere_atmosphere(self, path, Tcol=3, gcol=4, mucol=1, Ecol=0, spe_Icol=2):

        if 'nsx' in path:
            # Read and set attributes of NSX model table
            logT, logg, mu, logE, buf = self.load_NSX_table( path ,Tcol, gcol, mucol, Ecol, spe_Icol)
            self._elsewhere_atmosphere = (logT, logg, mu, logE, buf)

        else:
            """ Implement if required. """
            raise NotImplementedError('Implement setter if required.')
    
    @property
    def everywhere_atmosphere(self):
        """ Get the numerical atmosphere buffers for eerywhere if used.

        To preload a numerical atmosphere into a buffer, subclass and
        overwrite the setter. The underscore attribute set by the setter
        must be an :math:`n`-tuple whose :math:`n^{th}` element is an
        :math:`(n-1)`-dimensional array flattened into a one-dimensional
        :class:`numpy.ndarray`. The first :math:`n-1`
        elements of the :math:`n`-tuple must each be an ordered one-dimensional
        :class:`numpy.ndarray` of parameter values for the purpose of
        multi-dimensional interpolation in the :math:`n^{th}` buffer. The
        first :math:`n-1` elements must be ordered to match the index
        arithmetic applied to the :math:`n^{th}` buffer. An example would be
        ``self._hot_atmosphere = (logT, logg, mu, logE, buf)``, where:
        ``logT`` is a logarithm of local comoving effective temperature;
        ``logg`` is a logarithm of effective surface gravity;
        ``mu`` is the cosine of the angle from the local surface normal;
        ``logE`` is a logarithm of the photon energy; and
        ``buf`` is a one-dimensional buffer of intensities of size given by
        the product of sizes of the first :math:`n-1` tuple elements.

        It is highly recommended that buffer preloading is used, instead
        of loading from disk in the customisable radiation field extension
        module, to avoid reading from disk for every signal
        (likelihood) evaluation. This can be a non-negligible waste of compute
        resources. By preloading in Python, the memory is allocated and
        references to that memory are not in general deleted until a sampling
        script exits and the kernel stops. The likelihood callback accesses
        the same memory upon each call without I/O.

        """
        return self._everywhere_atmosphere

    @everywhere_atmosphere.setter
    def everywhere_atmosphere(self, path, Tcol=3, gcol=4, mucol=1, Ecol=0, bufcol=2):

        if 'nsx' in path:
            # Read and set attributes of NSX model table
            logT, logg, mu, logE, buf = self.load_NSX_table( path ,Tcol, gcol, mucol, Ecol, bufcol)
            self._everywhere_atmosphere = (logT, logg, mu, logE, buf)

        else:
            """ Implement if required. """
            raise NotImplementedError('Implement setter if required.')

    @property
    def hot(self):
        """ Get the instance of :class:`~.HotRegion.HotRegion`. """
        return self._hot

    @property
    def elsewhere(self):
        """ Get the instance of :class:`~.Elsewhere.Elsewhere`. """
        return self._elsewhere

    @property
    def everywhere(self):
        """ Get the instance of :class:`~.Everywhere.Everywhere`. """
        return self._everywhere

    @property
    def surface(self):
        """ Get the instance of :class:`~.HotRegion.HotRegion` or
        :class:`~.Everywhere.Everywhere` depending on which approach is used
        in the modelling. """
        return self._surface

    @property
    def spacetime(self):
        """ Return instance of :class:`~.Spacetime.Spacetime`. """
        return self._spacetime

    @spacetime.setter
    def spacetime(self, obj):
        if not isinstance(obj, Spacetime):
            raise TypeError('Invalid type for spacetime object.')
        # otherwise store a reference to the spacetime object
        self._spacetime = obj

    @property
    def stokes(self):
        """ Get the stokes option. If True, a full Stokes vector is computed and
        stored in signal, signalQ, and signalU. """
        return self._stokes

    def embed(self, fast_total_counts, threads):
        """ Embed the photosphere in an ambient Schwarzschild spacetime.

        In other words, generate a discrete representation of the photospheric
        radiation field and the null mapping from the photosphere to infinity,
        for use in flux integrators called by distant observers.

        """
        if self._everywhere is not None:
            self._everywhere.embed(self._spacetime,
                                   self,
                                   threads)
        else:
            if self._elsewhere is not None:
                self._elsewhere.embed(self._spacetime, threads)

                if self._hot is not None:
                    self._hot.embed(self._spacetime,
                                    self,
                                    fast_total_counts,
                                    threads,
                                    self._elsewhere._compute_cellParamVecs)
            elif self._hot is not None:
                self._hot.embed(self._spacetime,
                                self,
                                fast_total_counts,
                                threads)

    def integrate(self, energies, threads):
        """ Integrate over the photospheric radiation field.

        :param energies:
            A one-dimensional :class:`numpy.ndarray` of energies in keV.

        :param int threads:
            Number of ``OpenMP`` threads to spawn for signal integration.

        :param bool stokes:
            If activated, a full Stokes vector is computed and stored in signal, signalQ, and signalU.

        """
        if self._everywhere is not None:
            if self._stokes:
                raise NotImplementedError('Stokes option for everywhere not implmented yet.')      
            spectrum = self._everywhere.integrate(self._spacetime,
                                                   energies,
                                                   threads,
                                                   self._everywhere_atmosphere)
            if spectrum.ndim == 1:
                self._signal = ((spectrum.reshape(-1,1),),)
            else:
                self._signal = ((spectrum,),)
        else:
            if self._elsewhere is not None:
                spectrum = self._elsewhere.integrate(self._spacetime,
                                                     energies,
                                                     threads,
                                                     *self._elsewhere_atmosphere)

            if self._hot is not None:
                try:
                    else_atm_ext = self._elsewhere.atm_ext
                except:
                    else_atm_ext = None

                if self._stokes:
                    self._signal, self._signalQ, self._signalU  = self._hot.integrate_stokes(self._spacetime,
                                                   energies,
                                                   threads,
                                                   self._hot_atmosphere,
                                                   self._hot_atmosphere_Q,
                                                   self._elsewhere_atmosphere,
                                                   else_atm_ext)
                    if not isinstance(self._signal[0], tuple):
                        self._signal = (self._signal,)
                    if not isinstance(self._signalQ[0], tuple):
                        self._signalQ = (self._signalQ,)
                    if not isinstance(self._signalU[0], tuple):
                        self._signalU = (self._signalU,)
                    #Rotate the Stokes parameters based on position of the spin axis:
                    chi_rad = self["spin_axis_position_angle"]
                    tempQ = [list(x) for x in self._signalQ]
                    tempU = [list(x) for x in self._signalU]
                    for ih in range(0,len(self._signalQ)):
                        for ic in range(0,len(self._signalQ[ih][:])):
                            tempQ[ih][ic] = _np.cos(2.0*chi_rad) * tempQ[ih][ic] - _np.sin(2.0*chi_rad) * tempU[ih][ic]
                            tempU[ih][ic] = _np.sin(2.0*chi_rad) * tempQ[ih][ic] + _np.cos(2.0*chi_rad) * tempU[ih][ic]
                    self._signalQ = tuple(map(tuple, tempQ))
                    self._signalU = tuple(map(tuple, tempU))
                else:
                    self._signal = self._hot.integrate(self._spacetime,
                                                   energies,
                                                   threads,
                                                   self._hot_atmosphere,
                                                   self._elsewhere_atmosphere,
                                                   else_atm_ext)
                    if not isinstance(self._signal[0], tuple):
                        self._signal = (self._signal,)

                # add time-invariant component to first time-dependent component
                if self._elsewhere is not None:
                    for i in range(self._signal[0][0].shape[1]):
                        self._signal[0][0][:,i] += spectrum

    @property
    def signal(self):
        """ Get the stored signal (Stokes I).

        :returns:
            A tuple of tuples of *ndarray[m,n]*.
            Here :math:`m` is the number of energies, and
            :math:`n` is the number of phases. Units are photon/s/keV; the
            distance is a fast parameter so the areal units are not yet
            factored in. If the signal is a spectrum because the signal is
            time-invariant, then :math:`n=1`.

        """
        return self._signal
        
    @property
    def signalQ(self):
        """ Get the stored Stokes Q signal.

        :returns:
            A tuple of tuples of *ndarray[m,n]*.
            Here :math:`m` is the number of energies, and
            :math:`n` is the number of phases. Units are photon/s/keV; the
            distance is a fast parameter so the areal units are not yet
            factored in. If the signal is a spectrum because the signal is
            time-invariant, then :math:`n=1`.

        """
        if not self._stokes:
            raise Exception("Need to set stokes=True for the Photosphere object "
            "to calculate Stokes Q signal.")
        return self._signalQ
        
    @property
    def signalU(self):
        """ Get the stored Stokes U signal.

        :returns:
            A tuple of tuples of *ndarray[m,n]*.
            Here :math:`m` is the number of energies, and
            :math:`n` is the number of phases. Units are photon/s/keV; the
            distance is a fast parameter so the areal units are not yet
            factored in. If the signal is a spectrum because the signal is
            time-invariant, then :math:`n=1`.

        """
        if not self._stokes:
            raise Exception("Need to set stokes=True for the Photosphere object "
            "to calculate Stokes U signal.")
        return self._signalU
