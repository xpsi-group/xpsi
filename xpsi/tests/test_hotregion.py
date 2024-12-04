import numpy as np
from xpsi.HotRegion import HotRegion
from xpsi.Spacetime import Spacetime
from xpsi.Photosphere import Photosphere
import pytest

class TestHotRegion(object):
    
    def setup_class(cls):
        # Just to comment, this all relies on SpaceTime and Atmosphere classes
        # working correctly, so those should be unit tested separately.
        cls.super_colatitude_bounds = (None, None)
        cls.super_radius_bounds = (None, None)
        cls.phase_shift_bounds = (0.0, 0.1)
        cls.super_temperature_bounds = (5.1, 6.8)

        cls.hot_bounds = {'super_colatitude': cls.super_colatitude_bounds,
                      'super_radius': cls.super_radius_bounds,
                      'phase_shift': cls.phase_shift_bounds,
                      'super_temperature': cls.super_temperature_bounds,}
	
        cls.hot_values = {'super_radius': np.pi/4,
                          'super_colatitude': np.pi/4,
                          'super_temperature': 6.}

        cls.HotRegion = HotRegion(cls.hot_bounds, cls.hot_values)
        cls.space_values = {'frequency': 401., #Hz
                            'distance': 5., #kpc
                            'mass': 2., #Msun
                            'radius': 10., #km
                            'cos_inclination': 0.5
                            }
        cls.space_bounds = {'frequency': (None, None), #Hz
                            'distance': (None, None), #kpc
                            'mass': (None, None), #Msun
                            'radius': (None, None), #km
                            'cos_inclination': (None, None)
                            }
        
        cls.Spacetime = Spacetime(bounds=cls.space_bounds, values=cls.space_values)
        cls.mock_photosphere = {'mode_frequency': cls.space_values['frequency']}
        cls.Photosphere = Photosphere(hot=cls.HotRegion, bounds={'mode_frequency': (None, None)}, 
                                      values={'mode_frequency': cls.space_values['frequency']})
    
    def test_hotregion_class_initializes(self):
        HotRegion(self.hot_bounds, self.hot_values)
    
    def test_hotregion_breaks_with_bounds_missing(self):
        with pytest.raises(TypeError):
            HotRegion(self.hot_values)
            
    def test_hotregion_breaks_with_values_missing(self):
        with pytest.raises(TypeError):
            HotRegion(self.hot_bounds)
            
    def test_hotregion_breaks_if_phase_shift_bounds_unspecified(self):
        bounds = self.hot_bounds.copy()
        bounds['phase_shift'] = (None, None)
        with pytest.raises(ValueError):
            HotRegion(bounds, self.hot_values)
            
    def test_hotregion_breaks_if_phase_shift_bounds_infinite(self):
        bounds = self.hot_bounds.copy()
        bounds['phase_shift'] = (-np.inf, np.inf)
        with pytest.raises(ValueError):
            HotRegion(bounds, self.hot_values)
            
    def test_hotregion_breaks_if_phase_shift_bounds_more_than_one_cycle(self):
        bounds = self.hot_bounds.copy()
        bounds['phase_shift'] = (0, 1.1)
        with pytest.raises(ValueError):
            HotRegion(bounds, self.hot_values)
            
    def test_hotregion_produces_default_parameters(self):
        self.HotRegion.get_param('phase_shift')
        self.HotRegion.get_param('super_colatitude')
        self.HotRegion.get_param('super_radius')
        self.HotRegion.get_param('super_temperature')
        self.HotRegion.get_param('cede_colatitude')
        self.HotRegion.get_param('cede_radius')
        self.HotRegion.get_param('cede_azimuth')
        self.HotRegion.get_param('omit_colatitude')
        self.HotRegion.get_param('omit_radius')
        self.HotRegion.get_param('omit_azimuth')
        
        
    def test_hotregion_handles_cede(self):
        bounds = self.hot_bounds.copy()
        bounds['cede_temperature'] = (5.1, 6.8)
        HotRegionCede = HotRegion(bounds, self.hot_values, cede = True)
        HotRegionCede.get_param('cede_temperature')
        
    @pytest.mark.parametrize("declaration, split, expected_error", [
    (True, False, None),                 # No error, symmetry=True, split=False
    (False, False, None),                # No error, symmetry=False, split=False
    (False, True, "TypeError"),          # Error, split=True, symmetry=False
    ("invalid", False, "TypeError"),     # Error, invalid type for symmetry
])
    def test_symmetry(self, declaration, split, expected_error):
        hot_region = HotRegion(self.hot_bounds, self.hot_values, split=split)
    
        if expected_error:
            with pytest.raises(eval(expected_error)):
                hot_region.symmetry = declaration
        else:
            hot_region.symmetry = declaration
            assert hot_region._integrator is not None
            assert hot_region._integratorIQU is not None
            
    @pytest.mark.parametrize("num_rays", [50, 50.5, "invalid"])
    def test_set_num_rays(self, num_rays):
        # Test that integer, float and invalid values for num_rays have corect set and get. Raise error if invalid.
        
        fast_num_rays = 50
        if isinstance(num_rays, float | int):
            self.HotRegion.set_num_rays(num_rays, fast_num_rays)
            assert self.HotRegion.num_rays == int(num_rays), f"Expected num_rays to be {num_rays}, got {self.HotRegion.num_rays}"
            assert self.HotRegion._fast_num_rays == fast_num_rays, f"Expected _fast_num_rays to be {fast_num_rays}, got {self.HotRegion._fast_num_rays}"    
        else:
            with pytest.raises(ValueError):
                self.HotRegion.set_num_rays(num_rays, fast_num_rays)
                
    def test_image_order_limit_getter(self):
        # Set a value directly to the private attribute
        self.HotRegion._image_order_limit = 5
    
        # Validate the getter returns the correct value
        assert self.HotRegion.image_order_limit == 5, (
            f"Expected image_order_limit to be 5, but got {self.HotRegion.image_order_limit}"
        )

    @pytest.mark.parametrize("image_order_limit", [10, None, "invalid", 3.5])
    def test_image_order_limit_setter_valid(self, image_order_limit):
        # Test with a positive integer
        if image_order_limit == 10:
            self.HotRegion.image_order_limit = 10
            assert self.HotRegion._image_order_limit == 10, (
                f"Expected _image_order_limit to be 10, but got {self.HotRegion._image_order_limit}"
            )
        elif image_order_limit == None:
            self.HotRegion.image_order_limit = None
            assert self.HotRegion._image_order_limit is None, (
                f"Expected _image_order_limit to be None, but got {self.HotRegion._image_order_limit}"
            )
        elif image_order_limit == "invalid":
            # Test with a non-integer value
            with pytest.raises(TypeError, match="Image order limit must be an positive integer if not None."):
                self.HotRegion.image_order_limit = "invalid"
        elif image_order_limit == 3.5:
            with pytest.raises(TypeError, match="Image order limit must be an positive integer if not None."):
                self.HotRegion.image_order_limit = 3.5

    @pytest.mark.parametrize(
        "sqrt_num_cells, min_sqrt_num_cells, max_sqrt_num_cells, fast_sqrt_num_cells, fast_min_sqrt_num_cells, fast_max_sqrt_num_cells, expect_error",
        [
            # Valid cases
            (32, 10, 80, 16, 4, 16, None),           # Standard valid input
            (4, 4, 4, 4, 4, 4, None),               # Minimal valid input
            (80, 10, 80, 16, 4, 16, None),          # Edge of max valid

            # Invalid cases
            (3, 10, 80, 16, 4, 16, ValueError),     # sqrt_num_cells < 4
            ("invalid", 10, 80, 16, 4, 16, ValueError), # Non-integer sqrt_num_cells
            (32, 80, 10, 16, 4, 16, AssertionError), # min > max for regular cells
            (32, 10, 80, 16, 16, 4, AssertionError), # min > max for fast cells
        ],
    )
    def test_set_num_cells(self, sqrt_num_cells, min_sqrt_num_cells, max_sqrt_num_cells,
                           fast_sqrt_num_cells, fast_min_sqrt_num_cells, fast_max_sqrt_num_cells,
                           expect_error):
        if expect_error:
            with pytest.raises(expect_error):
                self.HotRegion.set_num_cells(
                    sqrt_num_cells, min_sqrt_num_cells, max_sqrt_num_cells,
                    fast_sqrt_num_cells, fast_min_sqrt_num_cells, fast_max_sqrt_num_cells
                )
        else:
            self.HotRegion.set_num_cells(
                sqrt_num_cells, min_sqrt_num_cells, max_sqrt_num_cells,
                fast_sqrt_num_cells, fast_min_sqrt_num_cells, fast_max_sqrt_num_cells
            )

            # Validate attributes for regular cells
            assert self.HotRegion.sqrt_num_cells == sqrt_num_cells
            assert self.HotRegion._num_cells == sqrt_num_cells**2
            assert self.HotRegion._min_sqrt_num_cells == min_sqrt_num_cells
            assert self.HotRegion._max_sqrt_num_cells == max_sqrt_num_cells

            # Validate attributes for fast cells
            assert self.HotRegion._fast_num_cells == fast_sqrt_num_cells**2
            assert self.HotRegion._fast_min_sqrt_num_cells == fast_min_sqrt_num_cells
            assert self.HotRegion._fast_max_sqrt_num_cells == fast_max_sqrt_num_cells

    def test_leaves_property(self):
        # Set a value directly to the private attribute
        self.HotRegion._leaves = "test_value"

        # Verify that the property returns the correct value
        assert self.HotRegion.leaves == "test_value", (
            f"Expected leaves to return 'test_value', but got {self.HotRegion.leaves}"
        )
        
    @pytest.mark.parametrize(
        "num_leaves, num_phases, phases, fast_num_leaves, fast_num_phases, fast_phases, do_fast, expect_error",
        [
            # Valid cases not fast
            (5, None, None, None, None, None, False, None),
            (5, None, np.linspace(0, 1, 5), None, None, None, False, None),
            (5, 6, None, None, None, None, False, None),
            (5, None, np.linspace(0, 1, 5), None, None, None, False, None),
            
            
            # Invalid cases
            (None, None, None, None, None, None, False, TypeError), # not allowed no leaves
            (5, "not_a_number", None, None, None, None, False, ValueError), # not allowed invalid num_phases
            (5, None, "not_an_array", None, None, None, False, TypeError), # not allowed invalid phases
            (5, None, np.array([0, 2, 1]), None, None, None, False, TypeError), # non monotonical
            (5, None, np.array([-0.1, 0.5, 1]), None, None, None, False, TypeError), # negative
            (5, None, np.array([0, 0.5, 1.1]), None, None, None, False, TypeError), # more than 1

            
            # Fast valid cases
            (5, None, None, 5, None, None, True, None),
            (5, None, np.linspace(0, 1, 5), 5, None, np.linspace(0, 1, 5), True, None),
            (5, 6, None, 5, 6, None, True, None),
            (5, None, np.linspace(0, 1, 5), 5, 6, np.linspace(0, 1, 5), True, None),
            
            # Fast invalid cases
            (5, None, None, None, None, None, True, TypeError), # not allowed no leaves
            (5, None, None, 5, "not_a_number", None, True, ValueError), # not allowed invalid num_phases
            (5, None, None, 5, None, "not_an_array", True, TypeError), # not allowed invalid phases
            (5, None, None, 5, None, np.array([0, 2, 1]), True, TypeError), # non monotonical
            (5, None, None, 5, None, np.array([-0.1, 0.5, 1]), True, TypeError), # negative
            (5, None, None, 5, None, np.array([0, 0.5, 1.1]), True, TypeError), # more than 1
        ]
    )
    def test_set_phases(self,num_leaves, num_phases, phases, fast_num_leaves, fast_num_phases, fast_phases, do_fast, expect_error):    
        HotRegion.do_fast = do_fast
        # Expect error for invalid cases
        if expect_error:
            with pytest.raises(expect_error):
                self.HotRegion.set_phases(num_leaves, num_phases=num_phases, phases=phases, fast_num_leaves=fast_num_leaves, fast_num_phases=fast_num_phases, fast_phases=fast_phases)
        else:
            # Test valid cases
            self.HotRegion.set_phases(num_leaves, num_phases=num_phases, phases=phases, fast_num_leaves=fast_num_leaves, fast_num_phases=fast_num_phases, fast_phases=fast_phases)
            if phases is None and num_phases is None:
                num_temp = num_leaves
                linspace_temp = np.linspace(0, 1, num_temp)
            elif phases is None:
                num_temp = num_phases
                linspace_temp = np.linspace(0, 1, num_temp)
            else:
                linspace_temp = phases
            
            assert np.allclose(self.HotRegion._phases_cycles, linspace_temp), f"Expected phases {linspace_temp}, but got {self.HotRegion._phases_cycles}"
            assert np.allclose(self.HotRegion._leaves, np.linspace(0.0, 2*np.pi, 5)), f"Expected leaves to be np.linspace(0.0, 2*np.pi, 5), but got {self.HotRegion._leaves}"
            
        #HotRegion.do_fast = False # restore
        
        
    def test_phases_in_cycles_property(self):
        # Set a value directly to the private attribute
        self.HotRegion._phases_cycles = "test_value"
        self.HotRegion._fast_phases_cycles = "test_value"

        # Verify that the property returns the correct value
        assert self.HotRegion.phases_in_cycles == "test_value", (
            f"Expected phases_in_cycles to return 'test_value', but got {self.HotRegion.phases_in_cycles}"
        )
        assert self.HotRegion.fast_phases_in_cycles == "test_value", (
            f"Expected fast_phases_in_cycles to return 'test_value', but got {self.HotRegion.fast_phases_in_cycles}"
        )
        
    def test_num_cells_property(self):
        # Set a value directly to the private attribute
        self.HotRegion._num_cells = "test_value"

        # Verify that the property returns the correct value
        assert self.HotRegion.num_cells == "test_value", (
            f"Expected num_cells to return 'test_value', but got {self.HotRegion.num_cells}"
        )
        
        
    @pytest.mark.parametrize(
        "input_value, expected_value, expect_error",
        [
            # Valid cases
            (True, True, None),  # Set True
            (False, False, None),  # Set False
            
            # Invalid cases
            (1, None, TypeError),  # Non-boolean integer
            ("True", None, TypeError),  # Non-boolean string
            (None, None, TypeError),  # None is invalid
        ]
    )
    def test_is_antiphased(self, input_value, expected_value, expect_error):
    
        # Handle invalid cases
        if expect_error:
            with pytest.raises(expect_error):
                self.HotRegion.is_antiphased = input_value
        else:
            # Set the value
            self.HotRegion.is_antiphased = input_value
            
            # Assert getter returns the correct value
            assert self.HotRegion.is_antiphased == expected_value, f"Expected {expected_value}, but got {self.HotRegion.is_antiphased}"

    @pytest.mark.parametrize(
        "input_extension, expected_value, expect_error",
        [
            # Valid cases
            ("BB", 1, None),
            ("Num4D", 2, None),
            ("Pol_BB_Burst", 3, None),
            ("Pol_Num2D", 4, None),
            ("user", 5, None),
            ("Num5D", 6, None),
    
            # Invalid cases
            ("InvalidExt", None, TypeError),  # Not a recognized string
            (123, None, TypeError),          # Non-string input
            (None, None, TypeError),         # None is invalid
        ]
    )
    def test_atm_ext(self, input_extension, expected_value, expect_error):
        # For invalid cases, assert that the proper exception is raised
        if expect_error:
            with pytest.raises(expect_error):
                self.HotRegion.atm_ext = input_extension
        else:
            # Set the value
            self.HotRegion.atm_ext = input_extension
            
            # Assert the getter returns the correct value
            assert self.HotRegion.atm_ext == expected_value, f"Expected {expected_value}, but got {self.HotRegion.atm_ext}"
    
    
    @pytest.mark.parametrize(
        "input_option, expected_value, expect_error",
        [
            # Valid cases
            (0, 0, None),
            (1, 1, None),
            (2, 2, None),
            (3, 3, None),
    
            # Invalid cases
            (4, None, TypeError),      # Not in the valid range
            (-1, None, TypeError),     # Negative value
            ("1", None, TypeError),    # Non-integer input
            (None, None, TypeError),   # None is invalid
        ]
    )
    def test_beam_opt(self, input_option, expected_value, expect_error):
        # For invalid cases, assert that the proper exception is raised
        if expect_error:
            with pytest.raises(expect_error):
                self.HotRegion.beam_opt = input_option
        else:
            # Set the value
            self.HotRegion.beam_opt = input_option
            
            # Assert the getter returns the correct value
            assert self.HotRegion.beam_opt == expected_value, f"Expected {expected_value}, but got {self.HotRegion.beam_opt}"

    def test_print_settings(self, capfd):
        # Call the method
        self.HotRegion.print_settings()
        
        # Capture the output
        captured = capfd.readouterr()
        
        # Assert that something was printed
        assert captured.out != "", "Expected some output, but nothing was printed."

    def test_construct_cellmesh_runs(self):
        # This function calls multiple underlying functions, so it is more of an integration test.
        # Example valid case, not asserting the actual outcome yet.
        self.HotRegion.fast_mode = False
        self.HotRegion.set_num_cells(32, 10, 80, 16, 4, 16)
        self.HotRegion._HotRegion__construct_cellMesh(self.Spacetime, None, 1)
        # Not yet testing invalid cases and failing cases yet. Also: not testing yet underlying functions individually.
        
    @pytest.mark.parametrize(
        "super_cell_area, cede_cell_area, expected",
        [
            (42, 84, (42, 84)),  # Case where both attributes are defined
            (42, None, (42, None)),  # Case where _cede_cellArea is missing
        ],
    )
    def test_cell_area(self, super_cell_area, cede_cell_area, expected):        
        # Mock the _super_cellArea attribute
        self.HotRegion._super_cellArea = super_cell_area
    
        # Conditionally set or remove _cede_cellArea
        if cede_cell_area is not None:
            self.HotRegion._cede_cellArea = cede_cell_area
        elif hasattr(self.HotRegion, '_cede_cellArea'):
            del self.HotRegion._cede_cellArea
    
        # Access the private property using name mangling
        cell_area = self.HotRegion._HotRegion__cellArea
    
        # Assert the expected value
        assert cell_area == expected

    def test_compute_rays_runs(self):
        # just testing here that the __compute_rays function runs. Calibrate
        # lags is called within this function at the end.
        self.HotRegion._HotRegion__compute_rays(self.Spacetime, 
                                                self.Photosphere, 1)
        #print('superlag', self.HotRegion._super_lag)
        #print('superlag', self.HotRegion._super_lag.shape)
        
        
    @pytest.mark.parametrize('cede', [True, False])
    def test_calibrate_lag_runs(self, cede):
        # just testing here that the __calibrate_lag function runs. Not relying
        # on compute_rays to do that, since that function relies on this 
        # function.

        mock_r_s_over_r = 0.5924*np.ones(42)
        self.HotRegion._r_s_over_r = mock_r_s_over_r
        mock_super_lag = 0.1*np.ones((42, 50))
        self.HotRegion._super_lag = mock_super_lag
        if cede:
            self.HotRegion._cede_r_s_over_r = mock_r_s_over_r
            self.HotRegion._cede_lag = mock_super_lag
        elif not cede:
            pass
        self.HotRegion._HotRegion__calibrate_lag(self.Spacetime, 
                                                 self.Photosphere)
        
        
    def test_cellParamVecs_runs(self):
        # just testing if compute_cellParamVecs runs. I am relying on construct
        # cellmesh and compute rays, but I tested those separately above.
        self.HotRegion._HotRegion__construct_cellMesh(self.Spacetime,
                                  None,
                                  1)
        self.HotRegion._HotRegion__compute_rays(self.Spacetime, 
                                                self.Photosphere, 1)
        self.HotRegion._HotRegion__compute_cellParamVecs()

    def test_psi(self):
        # Define test inputs (theta, phi, colatitude)
        theta = np.array([0, np.pi/4, np.pi/2])  # Example angles in radians
        phi = np.array([0, np.pi/2, np.pi])      # Example angles in radians
        colatitude = np.array([np.pi/4, np.pi/2, 3*np.pi/4])  # Example colatitudes
    
        # Expected results (precomputed for the given inputs)
        expected_results = np.arccos(
            np.cos(colatitude) * np.cos(theta)
            + np.sin(colatitude) * np.sin(theta) * np.cos(phi)
        )
    
        # Call the psi function
        results = self.HotRegion.psi(theta, phi, colatitude)
    
        # Validate the results against the expected values
        np.testing.assert_allclose(results, expected_results, rtol=1e-6, atol=1e-8)
    
    def test_embed_runs(self):
        # Just testing that the embed function runs.
        self.HotRegion.embed(self.Spacetime, self.Photosphere, None, 1)
        # Also testing the optional args that produce a correction.           
        # Mock function to be passed as *args
        def mock_correction_function(vector):
            return np.ones((len(vector), len(vector), 1))

        self.HotRegion.embed(self.Spacetime, self.Photosphere, None, 1, mock_correction_function)
        
        
    def test_do_fast_getter_setter(self):
        # Initial state of _do_fast
        initial_value = False
        self.HotRegion.do_fast = initial_value
        
        # Test the getter
        assert self.HotRegion.do_fast == initial_value, f"Expected {initial_value}, but got {self.HotRegion.do_fast}"
    
        # Test the setter
        new_value = True
        self.HotRegion.do_fast = new_value
        assert self.HotRegion.do_fast == new_value, f"Expected {new_value}, but got {self.HotRegion.do_fast}"

    def test_fast_mode_getter_setter(self):
        # Initial state of _do_fast
        initial_value = False
        self.HotRegion.fast_mode = initial_value
        
        # Test the getter
        assert self.HotRegion.fast_mode == initial_value, f"Expected {initial_value}, but got {self.HotRegion.do_fast}"
    
        # Test the setter
        new_value = True
        self.HotRegion.fast_mode = new_value
        assert self.HotRegion.fast_mode == new_value, f"Expected {new_value}, but got {self.HotRegion.do_fast}"

    def test_cede_getter_setter(self):
        # Test with initial value of _cede
        initial_value = False
        self.HotRegion._cede = initial_value
    
        # Test the getter
        assert self.HotRegion.cede == initial_value, f"Expected {initial_value}, but got {self.HotRegion.cede}"
    
        # Test the setter with a valid boolean value (True)
        self.HotRegion.cede = True
        assert self.HotRegion.cede == True, f"Expected True, but got {self.HotRegion.cede}"
    
        # Test the setter with a valid boolean value (False)
        self.HotRegion.cede = False
        assert self.HotRegion.cede == False, f"Expected False, but got {self.HotRegion.cede}"


    def test_omit_getter_setter(self):
        # Test with initial value of _omit
        initial_value = False
        self.HotRegion._omit = initial_value
    
        # Test the getter
        assert self.HotRegion.omit == initial_value, f"Expected {initial_value}, but got {self.HotRegion.omit}"
    
        # Test the setter with a valid boolean value (True)
        self.HotRegion.omit = True
        assert self.HotRegion.omit == True, f"Expected True, but got {self.HotRegion.omit}"
    
        # Test the setter with a valid boolean value (False)
        self.HotRegion.omit = False
        assert self.HotRegion.omit == False, f"Expected False, but got {self.HotRegion.omit}"
        

    def test_concentric_getter_setter(self):
        # Test with initial value of _concentric
        initial_value = False
        self.HotRegion._concentric = initial_value
    
        # Test the getter
        assert self.HotRegion.concentric == initial_value, f"Expected {initial_value}, but got {self.HotRegion.concentric}"
    
        # Test the setter with a valid boolean value (True)
        self.HotRegion.concentric = True
        assert self.HotRegion.concentric == True, f"Expected True, but got {self.HotRegion.concentric}"
    
        # Test the setter with a valid boolean value (False)
        self.HotRegion.concentric = False
        assert self.HotRegion.concentric == False, f"Expected False, but got {self.HotRegion.concentric}"
        
        
        
    def test_integrate_runs(self):
        # Just testing that the embed function runs.
        # split=False
        # hot_region = HotRegion(self.hot_bounds, self.hot_values, split=split)
        # hot_region.symmetry = True
        # print(hot_region._integrator)
        self.HotRegion.embed(self.Spacetime, self.Photosphere, None, 1)
        

        # just testing if the integrate runs, first with the default BB 
        # atmosphere given a hotregion and no elsewhere
        hot_atmosphere = ()
        elsewhere_atmosphere = ()
        atm_ext_else = None
        energies = np.linspace(1,2,10) # keV
        self.HotRegion.integrate(self.Spacetime, energies, 1, hot_atmosphere, elsewhere_atmosphere, atm_ext_else)
        
