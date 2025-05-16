import h5py 
from xpsi.global_imports import *

def load_group_from_hdf5(filepath: str, ext: str = 'post_unweighted'):
    """ Load a specific group from a hdf5 based on the extension given. The 
    root + extensions are complementary to the standard output files from MultiNest.

    :param filepath: Path to HDF5 file. 
    :type filerootpath: str

    :param ext: Extension to the root name, defaults to 'post_unweighted' 
        (which is similar to the MultiNest [root].txt output file). 
    :type ext: str, optional

    :return: N-D array with the data.  
    :rtype: np.array
    """
    with h5py.File(filepath, 'r') as f:
        # get specific data group matching to a original multinest output file 
        data_group = f[ext]

        # read each column in order into a list
        col_names = list(data_group.keys())
        data_cols = [data_group[name][()] for name in col_names]

        # stack columns into the ND array 
        data = _np.column_stack(data_cols)

    return data 

def add_transform_group_to_hdf5(filerootpath: str, data: _np.array,
                                paramnames: list, ext: str = 'post_unweighted'):
    """Used for transforming a specific group/multinest file in the hdf5 file.

    :param filerootpath: Path to HDF5 file. 
    :type filerootpath: str

    :param data: Transformed data
    :type data: np.array

    :param paramnames: list with parameter names including the derived ones! (and weights, etc?)
    :type paramnames: list

    :param ext: Extension to the root name, defaults to 'post_unweighted' 
        (which is similar to the MultiNest [root].txt output file). 
    :type ext: str, optional
    """
    # Set metadata depending on extension 
    if ext == 'post_unweighted':
        metadata = ['weights', '-2logl'] + paramnames
    if ext == 'dead-birth':
        metadata = paramnames + ['logl', 'logl_birth', 
                                'log_prior_mass', 'node_no']
    if ext == 'phys_live-birth':
        metadata = paramnames + ['logl', 'logl_birth', 
                                'node_no']
    
    # Check if metadata matches amount of columns in the transformed data 
    if len(metadata) != (data.shape[0] if data.ndim == 1 else data.shape[1]):
        raise ValueError(
            f"Length mismatch: metadata has {len(metadata)} items, "
            f"but data has {(data.shape[0] if data.ndim == 1 else data.shape[1])} columns. "
            "Check if xpsi.Runs.load_runs.names is correct."
            )
    
    with h5py.File(filerootpath, 'a') as f:
        # Overwrite transform dataset if it already exists by removing it and 
        # creating a new one
        if f'transformed-{ext}' in f:
            print("ext", f[ext])
            del f[f'transformed-{ext}']

        print("keys", f.keys())
              
        data_group = f.create_group(f"transformed-{ext}", track_order=True) 

        # Add datasets, one for each column
        for idx, name in enumerate(metadata):
            data_group.create_dataset(name, data=data[:,idx])


def check_group_exists_in_hdf5(filepath: str, ext: str = 'post_unweighted'):
    """Check whether a certain group exists inside the HDF5. 

    :param filerootpath: Path to HDF5 file. 
    :type filerootpath: str

    :param ext: Extension to the root name, defaults to 'post_unweighted' 
        (which is similar to the MultiNest [root].txt output file). 
    :type ext: str, optional

    :return: Returns True if the group exits, otherwise returns False. 
    :rtype: bool
    """
    with h5py.File(filepath, 'r') as f:
        if ext in f:
            return True
        else:
            return False

def get_param_data_from_hdf5(filepath: str, ext: str = 'post_unweighted'):
    """ Load a specific group from a hdf5 based on the extension given. The 
    root + extensions are complementary to the standard output files from MultiNest.

    :param filepath: Path to HDF5 file. 
    :type filerootpath: str

    :param ext: Extension to the root name, defaults to 'post_unweighted' 
        (which is similar to the MultiNest [root].txt output file). 
    :type ext: str, optional

    :return: N-D array with the data.  
    :rtype: np.array
    """
    with h5py.File(filepath, 'r') as f:
        # get specific data group matching to a original multinest output file 
        data_group = f['transformed-post_unweighted']

        # read each column in order into a list
        col_names = list(data_group.keys())
        col_names_param = col_names.copy()  

        # remove all non-parameter columns         
        for drop in ('weights', '-2logl'):
            if drop in col_names_param:
                col_names_param.remove(drop) 

        print("col names 2", col_names)

        data_cols = [data_group[name][()] for name in col_names_param]

        # stack columns into the ND array 
        data = _np.column_stack(data_cols)

    return data 