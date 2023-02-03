import getdist
import getdist.plots
getdist.chains.print_load_details = False

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
