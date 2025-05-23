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

def process_multinest_run(file_root, base_dir, filetype, **kwargs):
    """Modified version of nestcheck.data_processing.process_multinest_run().
    Loads data (both .txt and .npy) from a MultiNest run into the nestcheck 
    dictionary format for analysis.

    N.B. producing required output file containing information about the
    iso-likelihood contours within which points were sampled (where they were
    "born") requies MultiNest version 3.11 or later.

    Parameters
    ----------
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
    # Load dead and live points
    if filetype == ".npy":
        dead = _np.load(_os.path.join(base_dir, file_root) + 'dead-birth.npy')
        live = _np.load(_os.path.join(base_dir, file_root)
                        + 'phys_live-birth.npy')
    elif filetype == ".txt":
        dead = _np.loadtxt(_os.path.join(base_dir, file_root) + 'dead-birth.txt')
        live = _np.loadtxt(_os.path.join(base_dir, file_root)
                        + 'phys_live-birth.txt')   
    # Remove unnecessary final columns
    dead = dead[:, :-2]
    live = live[:, :-1]
    assert dead[:, -2].max() < live[:, -2].min(), (
        'final live points should have greater logls than any dead point!',
        dead, live)
    ns_run = process_samples_array(_np.vstack((dead, live)), **kwargs)
    assert _np.all(ns_run['thread_min_max'][:, 0] == -_np.inf), (
        'As MultiNest does not currently perform dynamic nested sampling, all '
        'threads should start by sampling the whole prior.')
    ns_run['output'] = {}
    ns_run['output']['file_root'] = file_root
    ns_run['output']['base_dir'] = base_dir
    return ns_run
