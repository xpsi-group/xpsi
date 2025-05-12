from __future__ import division, print_function
import pandas as pd

from .global_imports import *
from . import global_imports

from .cellmesh.mesh_tools import allocate_cells as _allocate_cells
from .cellmesh.mesh import construct_spot_cellMesh as _construct_spot_cellMesh
from .cellmesh.polar_mesh import construct_polar_cellMesh as _construct_polar_cellMesh
from .cellmesh.rays import compute_rays as _compute_rays

from .Parameter import Parameter, Derive
from .ParameterSubspace import ParameterSubspace
import xpsi
from xpsi.global_imports import _c, _G, _dpr, gravradius, _csq, _km, _2pi
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate

class Temp_Interpolator_shells(ParameterSubspace):

	"""
    Interpolates shell temperature data onto an X-PSI angular grid.

    Supports input from `.npz`, `.csv`, or `.txt` files containing shell-averaged or snapshot data.
    Handles multiple surface regions (hotspots, elsewhere, everywhere) and computes temperature either
    from flux (if `bhac_shell_avg=True`) or directly from provided values.

    Attributes:
        filename (str): Path to input data file.
        coderes (int): Input grid resolution (coderes x coderes).
        xpsi_theta (np.ndarray): X-PSI grid polar angles.
        xpsi_phi (np.ndarray): X-PSI grid azimuthal angles.
        first_spot, second_spot, elsewhere_xpsi, everywhere_xpsi (bool): Region flags.
        bhac_shell_avg (bool): If True, compute temperature from flux; else use raw temperature.

    Methods:
        read_regrid(): Extract θ, φ, flux, and tracer from shell snapshot data (Used for BHAC (https://gitlab.itp.uni-frankfurt.de/BHAC-release/bhac) snapshots).
        read_average_shell(): Load θ, φ, temp from shell-averaged data file (Both for BHAC and general tabulated temperature values).
        temp_interpolation_func(): Interpolate temp onto the X-PSI mesh with boundary handling.
	"""
	required_names = ['filename',
					'first_spot',
					'second_spot',
					'elsewhere_xpsi',
					'everywhere_xpsi',
					'bhac_shell_avg']

	def __init__(self,
				filename ='data',
				coderes = 512, 
				xpsi_theta = None,
				xpsi_phi = None,
				first_spot = False,
				second_spot = False,
				elsewhere_xpsi = False,
				everywhere_xpsi = False,
				bhac_shell_avg = True,
				**kwargs):

		self.filename=filename
		self.num_cells_theta = num_cells_theta
		self.coderes = int(coderes)
		self.num_cells_phi = num_cells_phi
		self.xpsi_theta = xpsi_theta
		self.xpsi_phi = xpsi_phi
		self.first_spot = first_spot
		self.second_spot = second_spot
		self.elsewhere_xpsi = elsewhere_xpsi
		self.everywhere_xpsi = everywhere_xpsi
		self.bhac_shell_avg = bhac_shell_avg

	def read_regrid(self,filename,coderes):
		"""
		Read shell snapshots data at a fixed radius and reorganize it into angular grids 
		(theta, phi) along with flux and tracer values, depending on the selected surface region.

		This function assumes the input file contains a flattened 2D grid of shape 
		(coderes x coderes), with Cartesian coordinates and scalar fields such as 
		temperature and passive tracers. It extracts angular coordinates from the 
		Cartesian data and rearranges the data into symmetric grids for use in 
		X-PSI surface emission modeling.

		Args:
			filename (str): Path to the BHAC .csv file containing shell slice data (columns: X, Y, Z, T_MArt, tr1).
			coderes (int): Resolution of the input data grid. The input file is assumed to represent 
						a square grid of size (coderes x coderes).

		Raises:
			ValueError: If none of the region flags (first_spot, second_spot, elsewhere_xpsi, everywhere_xpsi) are set.

		Returns:
			Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: Arrays for:
				- phicode (2D array): azimuthal coordinate φ
				- thetacode (2D array): polar coordinate θ
				- Fluxcode (2D array): flux or temperature field
				- tracercode (2D array): passive tracer value
			All arrays are of shape (N, coderes), where N depends on the region type.
		"""

		data = pd.read_csv(filename)
		solarmass = 1.989e33
		c = 2.99792458e10
		G = 6.674e-8
		rstar = 10.0 * 1.e5
		rg = rstar/4.0
		M = (rg * c**2/G)/solarmass

		Ledd = 1.26e38 * M #1.68 mass of the NS in solar mass
		Mdot_units = 0.01*Ledd/(0.1 * c * c)
		rho_units = Mdot_units/(rg * rg * c ) #for efficiency = 0.1
		sigma = 5.6704e-5

		x = data['X']
		y = data['Y']
		z = data['Z']
		thetabhac = _np.arccos(z/_np.sqrt(x**2 + y**2 + z**2))
		phibhac = _np.arctan2(y,x) 
		Flux = data['T_MArt']
		k = int(coderes/2)
		tracer = data['tr1']
		num_cells_lt = 0
		
		if (self.first_spot or self.second_spot or self.elsewhere_xpsi or self.everywhere_xpsi) is False:
			raise ValueError('Neither of the 2 hotspot or elsewhere region specified')
		#Generating the required segments
		if self.first_spot:

			index_i = 0
			index_j = int(k - num_cells_lt)

			phicode = _np.zeros((index_j-index_i,coderes))
			thetacode = _np.zeros((index_j-index_i,coderes))
			Fluxcode = _np.zeros((index_j-index_i,coderes))
			tracercode = _np.zeros((index_j-index_i,coderes))

			for i in range(index_i,index_j):
				for j in range(0,k):

					phicode[i,j+k] = phibhac[j+i*coderes]
					thetacode[i,j+k] = thetabhac[j+i*coderes]

					Fluxcode[i,j+k] = Flux[j+i*coderes]
					tracercode[i,j+k] = tracer[j+i*coderes]
				for j in range(k,coderes):

					phicode[i,j-k] = phibhac[j+i*coderes]
					thetacode[i,j-k] = thetabhac[j+i*coderes]

					Fluxcode[i,j-k] = Flux[j+i*coderes]
					tracercode[i,j-k] = tracer[j+i*coderes]
			
		if self.second_spot:	

			index_i = int(k + num_cells_lt)
			index_j = coderes

			phicode = _np.zeros((index_j-index_i,coderes))
			thetacode = _np.zeros((index_j-index_i,coderes))
			Fluxcode = _np.zeros((index_j-index_i,coderes))
			tracercode = _np.zeros((index_j-index_i,coderes))

			for i in range(index_i,index_j):
				for j in range(0,k):

					phicode[i-index_i,j+k] = phibhac[j+i*coderes] 
					thetacode[i-index_i,j+k] = thetabhac[j+i*coderes]

					Fluxcode[i-index_i,j+k] = Flux[j+i*coderes]
					tracercode[i-index_i,j+k] = tracer[j+i*coderes]

				for j in range(k,coderes):

					phicode[i-index_i,j-k] = phibhac[j+i*coderes]
					thetacode[i-index_i,j-k] = thetabhac[j+i*coderes]

					Fluxcode[i-index_i,j-k] = Flux[j+i*coderes]
					tracercode[i-index_i,j-k] = tracer[j+i*coderes]

		if self.elsewhere_xpsi:

			num_cells_lt = 64
			index_i = int(k - num_cells_lt)+1
			index_j = int(k + num_cells_lt)-1

			phicode = _np.zeros((index_j-index_i,coderes))
			thetacode = _np.zeros((index_j-index_i,coderes))
			Fluxcode = _np.zeros((index_j-index_i,coderes))
			tracercode = _np.zeros((index_j-index_i,coderes))

			for i in range(index_i,index_j):
				for j in range(0,k):

					phicode[i-index_i,j+k] = _np.arctan2(y[j+i*coderes],x[j+i*coderes]) 
					thetacode[i-index_i,j+k] = _np.arccos(z[j+i*coderes]/_np.sqrt(x[j+i*coderes]**2 \
						+ y[j+i*coderes]**2 + z[j+i*coderes]**2))

					Fluxcode[i-index_i,j+k] = Flux[j+i*coderes]
					tracercode[i-index_i,j+k] = tracer[j+i*coderes]

				for j in range(k,coderes):

					phicode[i-index_i,j-k] = _np.arctan2(y[j+i*coderes],x[j+i*coderes]) 
					thetacode[i-index_i,j-k] = _np.arccos(z[j+i*coderes]/_np.sqrt(x[j+i*coderes]**2 \
						+ y[j+i*coderes]**2 + z[j+i*coderes]**2))

					Fluxcode[i-index_i,j-k] = Flux[j+i*coderes]
					tracercode[i-index_i,j-k] = tracer[j+i*coderes]

		if self.everywhere_xpsi:

			index_i = 0
			index_j = coderes

			phicode = _np.zeros((index_j-index_i,coderes))
			thetacode = _np.zeros((index_j-index_i,coderes))
			Fluxcode = _np.zeros((index_j-index_i,coderes))
			tracercode = _np.zeros((index_j-index_i,coderes))
			Flux_min = 5.e-9

			for i in range(index_i,index_j):
				for j in range(0,k):

					phicode[i,j+k] = phibhac[j+i*coderes]
					thetacode[i,j+k] = thetabhac[j+i*coderes]

					if (tracer[j+i*coderes] > 0.8):
						Fluxcode[i,j+k] = _np.log10(((_np.abs(Flux[j+i*coderes])*tracer[j+i*coderes] + (1.0 - tracer[j+i*coderes])*Flux_min)*rho_units*c**3/sigma)**(1./4.))
					else:
						Fluxcode[i,j+k] = _np.log10((((1.0 - tracer[j+i*coderes])*Flux_min)*rho_units*c**3/sigma)**(1./4.))
					tracercode[i,j+k] = tracer[j+i*coderes]
				for j in range(k,coderes):

					phicode[i,j-k] = phibhac[j+i*coderes]
					thetacode[i,j-k] = thetabhac[j+i*coderes]

					if (tracer[j+i*coderes] > 0.8):
						Fluxcode[i,j-k] = _np.log10(((_np.abs(Flux[j+i*coderes])*tracer[j+i*coderes] + (1.0 - tracer[j+i*coderes])*Flux_min)*rho_units*c**3/sigma)**(1./4.))
					else:
						Fluxcode[i,j-k] = _np.log10((((1.0 - tracer[j+i*coderes])*Flux_min)*rho_units*c**3/sigma)**(1./4.))
					tracercode[i,j-k] = tracer[j+i*coderes]

		return phicode,thetacode,Fluxcode,tracercode


	def read_average_shell(self,filename,
                        coderes,
                        bhac_shell_avg):
		"""
		Load angular grid data (theta, phi, temp) from a shell .npz, .csv, or .txt file.

		Supports either computing temperature from flux using physical scaling (if `bhac_shell_avg=True`)
		or reading precomputed temperature directly. The data is reshaped to (coderes, coderes),
		with input validation for CSV/TXT formats.

		Args:
			filename (str): Path to .npz, .csv, or .txt file containing BHAC shell slice data.
			coderes (int): Grid resolution per axis; total data points must equal coderes x coderes.
			bhac_shell_avg (bool): If True, compute temperature from flux; otherwise use 'temp' field.

		Returns:
			Tuple[np.ndarray, np.ndarray, np.ndarray]: theta, phi, and log10(temperature) arrays of shape (coderes, coderes).
		"""
		def compute_temperature_from_flux(flux):
			
			solarmass = 1.989e33
			c = 2.99792458e10
			G = 6.674e-8
			rstar = 10.0 * 1.e5
			rg = rstar/4.0                        #rg = G*M*solarmass/c**2
			M = (rg * c**2/G)/solarmass			  #Mdot is 1% of Eddington Limit [Mdot = 0.01 * Ledd/(efficiency * c^2)]
			Ledd = 1.26e38 * M                    #Eddington Luminosity = 1.26e38 (M/M_solarmass) erg/s
			Mdot_units = 0.01*Ledd/(0.1 * c * c)  #Mass accretion rate is 1% of Eddington
			rho_units = Mdot_units/(rg * rg * c ) #for efficiency = 0.1
			sigma = 5.6704e-5
			temp = _np.log10(((_np.abs(flux)*rho_units*c**3/sigma)**(1./4.)))
			return temp

		if filename.endswith('.npz'):
			data = _np.load(filename)
			theta = data['theta'][:,:,0]
			phi = data['phi'][:,:,0]
			if bhac_shell_avg:
				temp = compute_temperature_from_flux(data['flux'])
			else:
				temp = _np.log10(data['temp'])

		elif filename.endswith(('.csv', '.txt')):
			df = pd.read_csv(filename, delim_whitespace=filename.endswith('.txt'))
			required = ['theta', 'phi', 'temp']
			if not all(col in df.columns for col in required):
				raise ValueError(f"Missing required columns: {required}")

			expected_size = coderes * coderes
			actual_size = len(df['theta'])
			if actual_size != expected_size:
				raise ValueError(f"Data size mismatch: expected {expected_size} rows "
								f"(for coderes={coderes}), but got {actual_size}")

			theta = df['theta'].to_numpy().reshape(coderes, coderes)
			phi = df['phi'].to_numpy().reshape(coderes, coderes)
			tempcode = df['temp'].to_numpy().reshape(coderes, coderes)
			temp = _np.log10(tempcode)

		else:
			raise ValueError("Input file must be a .npz, .csv, or .txt file.")

		return theta,phi,temp

	def temp_interpolation_avg(self,
							thetacode,
							phicode,
							Tempcode,
							T_everywhere):
		"""
		Interpolate a temperature map (Tempcode) defined on an angular grid (thetacode, phicode)
		onto the X-PSI model grid (self.xpsi_theta, self.xpsi_phi), using cubic interpolation.

		Edge cells are manually corrected using boundary values from Tempcode, and any interpolated
		temperatures below T_everywhere are clipped to T_everywhere.

		Args:
			thetacode (np.ndarray): 2D array of polar angles (θ) from the input temperature map.
			phicode (np.ndarray): 2D array of azimuthal angles (φ) from the input temperature map.
			Tempcode (np.ndarray): 2D temperature map to be interpolated.
			T_everywhere (float): Minimum allowed temperature after interpolation (floor value).

		Returns:
			np.ndarray: Interpolated temperature map on the X-PSI angular grid (shape matches self.xpsi_theta).
		"""

		lenphi = (_np.shape(self.xpsi_phi)[1])
		lentheta = (_np.shape(self.xpsi_phi)[0])
		Tempxpsi = _np.zeros((lentheta,lenphi))

		Tempxpsi = interpolate.griddata((phicode.ravel(), thetacode.ravel()), Tempcode.ravel(), (self.xpsi_phi, self.xpsi_theta), method='cubic')

		for i in range(0,lentheta):
			for j in range(0,lenphi):
				#Take care of the edge cells
				if (self.xpsi_theta[i,j] < thetacode[0,_np.shape(thetacode)[1]-1]):
					Tempxpsi[i,j] = Tempcode[0,_np.shape(thetacode)[1]-1]
				elif (self.xpsi_theta[i,j] > thetacode[_np.shape(thetacode)[0]-1,_np.shape(thetacode)[1]-1]):
					Tempxpsi[i,j] = Tempcode[_np.shape(thetacode)[0]-1,_np.shape(thetacode)[1]-1]
				elif (self.xpsi_phi[i,j] > phicode[_np.shape(thetacode)[0]-1,_np.shape(phicode)[1]-1] or self.xpsi_phi[i,j] < phicode[_np.shape(thetacode)[0]-1,0]):
					Tempxpsi[i,j] = Tempcode[_np.shape(thetacode)[0]-1,0]
				if (Tempxpsi[i,j] < T_everywhere):
					Tempxpsi[i,j] = T_everywhere
				# if (self.xpsi_theta[i,j] < _np.pi/2):
				# 	Tempxpsi[i,j] = T_everywhere 
		return Tempxpsi
				
	def temp_interpolation_flux(self,
							thetacode,
							phicode,
							Fluxcode,
							tracercode,
							tracer_threshold,
							T_everywhere):

		lenphi = (_np.shape(self.xpsi_phi)[1])
		lentheta = (_np.shape(self.xpsi_phi)[0])
		Fluxxpsi = _np.zeros((lentheta,lenphi))
		Tempxpsi = _np.zeros((lentheta,lenphi))
		Tracerxpsi = _np.zeros((lentheta,lenphi))

		#...........................................................................
		Fluxxpsi = interpolate.griddata((phicode.ravel(), thetacode.ravel()), Fluxcode.ravel(), (self.xpsi_phi, self.xpsi_theta), method='cubic')
		Tracerxpsi = interpolate.griddata((phicode.ravel(), thetacode.ravel()), tracercode.ravel(), (self.xpsi_phi, self.xpsi_theta), method='cubic')	
		for i in range(0,lentheta):
			for j in range(0,lenphi):
				#Take care of the edge cells
				if (self.xpsi_theta[i,j] < thetacode[0,_np.shape(thetacode)[1]-1]):
					Fluxxpsi[i,j] = Fluxcode[0,_np.shape(thetacode)[1]-1]
					Tracerxpsi[i,j] = tracercode[0,_np.shape(thetacode)[1]-1]
				elif (self.xpsi_theta[i,j] > thetacode[_np.shape(thetacode)[0]-1,_np.shape(thetacode)[1]-1]):
					Fluxxpsi[i,j] = Fluxcode[_np.shape(thetacode)[0]-1,_np.shape(thetacode)[1]-1]
					Tracerxpsi[i,j] = tracercode[_np.shape(thetacode)[0]-1,_np.shape(thetacode)[1]-1]
				elif (self.xpsi_phi[i,j] > phicode[_np.shape(thetacode)[0]-1,_np.shape(phicode)[1]-1] or self.xpsi_phi[i,j] < phicode[_np.shape(thetacode)[0]-1,0]):
					Fluxxpsi[i,j] = Fluxcode[_np.shape(thetacode)[0]-1,0]
					Tracerxpsi[i,j] = tracercode[_np.shape(thetacode)[0]-1,0]
				
				Tempxpsi[i,j] = Fluxxpsi[i,j]
				
		return Tempxpsi