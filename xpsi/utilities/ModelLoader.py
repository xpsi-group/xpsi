import os, sys, importlib
from xpsi.global_imports import xpsiError

class ModelLoaderError(xpsiError):
    """ Raised if there is a problem with the module loader. """

def load_model( model_path , 
               config_path = None , 
               execution_path = None ):
    """ Function to load a model from a given python script. 
    The model can be loaded with or without a configuration file, from an execution directory or locally.

    Args:
        model_path (str): Path to the main file which loads the model.
        config_path (str|None, optional): Path to the configuration file, if needed, for loading the main. No configuration file is used if None. Defaults to None.
        execution_path (str|None, optional): Path to the execution directory, if needed, for loading the main. Execution happens in the current directory if None. Defaults to None.

    Returns:
        Namespace: Imported main with a X-PSI defined Neutron Star model. The different attributes are defined in the main.
    """

    # Check if the files exist
    assert os.path.isfile(model_path), f'File {model_path} does not exist.'
    abs_model_path = os.path.abspath(model_path)
    if config_path is not None:
        assert os.path.isfile(config_path), f'File {config_path} does not exist.'
        abs_config_path = os.path.abspath(config_path)

    # Save original values
    pwd = os.getcwd()
    original_sys_path = sys.path.copy()
    original_sys_modules = sys.modules.copy()

    # Change directory and make the symbolic link
    if execution_path is not None:
        os.chdir(execution_path)
    os.system(f'ln -s { abs_model_path } local_main.py')
    
    # Try to catch the error and make sure the cleanup happens
    try:

        spec = importlib.util.spec_from_file_location(f'local_main', 'local_main.py')
        model = importlib.util.module_from_spec(spec)
        sys.modules[f'local_main'] = model

        # Case with config_path provided
        if config_path is not None:
            os.system(f'ln -s { abs_config_path } local_config.ini')
            sys.argv = ['local_main.py', f'@local_config.ini']

        # Do the loading
        spec.loader.exec_module(model)

        # Successful
        print( f'Model located at {model_path}' + (f' with configuration {config_path}' if config_path is not None else '') + ' has been sucessfully loaded' )

        # Return the model
        return model
    
    except NameError:
        raise ModelLoaderError("The model could not be loaded because of an ill-defined import in the model at model_path. Look the traceback for more information.")

    except:
        raise ModelLoaderError("The model could not be loaded for unkown reasons. Look the traceback for more information.") 
    
    # Cleanup whatever happens
    finally:
        os.system("rm local_main.py")
        if config_path is not None:
            os.system("rm local_config.ini")
        os.chdir( pwd )
        sys.path = original_sys_path
        sys.modules = original_sys_modules