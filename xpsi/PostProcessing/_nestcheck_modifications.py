import getdist
import getdist.plots
getdist.chains.print_load_details = False
from ._global_imports import *
from nestcheck.data_processing import process_samples_array

def getdist_kde(x, samples, weights, **kwargs):
    """
    Implement the GetDist 1D Kernel Density Estimator.
    GetDist executes boundary correction for density estimation near
    the parameter limits. Limits are *required* for proper
    GetDist KDE at parameter boundaries, and can be passed via the kwargs.
    """
    settings = kwargs.get('settings', {'fine_bins': 1024,
                                       'smooth_scale_1D': -1.0,
                                       'boundary_correction_order': 1,
                                       'mult_bias_correction_order': 1})

    ranges = kwargs.get('ranges')
    if ranges is None:
        raise ValueError('Supply parameter bounds for KDE.')

    idx = kwargs.get('idx')

    bcknd = getdist.mcsamples.MCSamples(sampler='nested',
                                        samples=samples,
                                        weights=weights,
                                        names=['x'],
                                        ranges=dict(x=ranges[idx]),
                                        settings=settings)

    normalize = kwargs.get('normalize', False)
    if normalize:
        bcknd.get1DDensity('x').normalize(by='integral',
                                               in_place=True)
    return bcknd.get1DDensity('x').Prob(x)

def process_multinest_run(dead_points, live_points, file_root, base_dir, **kwargs):
    """Modified version of nestcheck.data_processing.process_multinest_run().
    Loads data (both .txt and .npy) from a MultiNest run into the nestcheck 
    dictionary format for analysis.

    N.B. producing required output file containing information about the
    iso-likelihood contours within which points were sampled (where they were
    "born") requies MultiNest version 3.11 or later.

    Parameters
    ----------
    dead_points: ndarray
        Dead points 
    live_points: ndarray
        Live points 
    file_root: str
        Root name for output files. When running MultiNest, this is determined
        by the nest_root parameter.
    base_dir: str
        Directory containing output files. When running MultiNest, this is
        determined by the nest_root parameter.
    kwargs: dict, optional
        Passed to ns_run_utils.check_ns_run (via process_samples_array)

    Returns
    -------
    ns_run: dict
        Nested sampling run dict (see the module docstring for more details).
    """
    # Remove unnecessary final columns
    dead_points = dead_points[:, :-2]
    live_points = live_points[:, :-1]
    assert dead_points[:, -2].max() < live_points[:, -2].min(), (
        'final live points should have greater logls than any dead point!',
        dead_points, live_points)
    ns_run = process_samples_array(_np.vstack((dead_points, live_points)), **kwargs)
    assert _np.all(ns_run['thread_min_max'][:, 0] == -_np.inf), (
        'As MultiNest does not currently perform dynamic nested sampling, all '
        'threads should start by sampling the whole prior.')
    ns_run['output'] = {}
    ns_run['output']['file_root'] = file_root
    ns_run['output']['base_dir'] = base_dir
    return ns_run
