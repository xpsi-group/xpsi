import os, sys, site, sysconfig, importlib, importlib.util, importlib.machinery
from xpsi.global_imports import xpsiError

class ModelLoaderError(xpsiError):
    """ Raised if there is a problem with the module loader. """

# Name under which the main of the model is executed
_MODEL_MODULE_NAME = 'local_main'


###############################################################
#### HELPERS TO ISOLATE THE USER CODE OF EACH LOADED MODEL ####
###############################################################
def _is_within(path, root):
    """ True if `path` is `root` or lies somewhere below it. """
    try:
        return os.path.commonpath([path, root]) == root
    except ValueError:  # e.g. different drives, or mix of absolute/relative
        return False

def _library_roots():
    """ Directories holding the standard library and installed third-party packages.
    Anything living there is shared between models and is never purged. """
    roots = set()
    for key in ('stdlib', 'platstdlib', 'purelib', 'platlib'):
        path = sysconfig.get_paths().get(key)
        if path:
            roots.add(os.path.realpath(path))
    try:
        for path in site.getsitepackages():
            roots.add(os.path.realpath(path))
    except AttributeError:  # some old virtualenvs do not provide it
        pass
    try:
        roots.add(os.path.realpath(site.getusersitepackages()))
    except AttributeError:
        pass
    return roots

def _module_locations(module):
    """ Real paths of the file (or of the package directories) a module comes from. """
    file = getattr(module, '__file__', None)
    if file:
        locations = [file]
    else:  # namespace packages have no __file__ but have a __path__
        try:
            locations = list(getattr(module, '__path__', []))
        except Exception:
            locations = []
    return [os.path.realpath(location) for location in locations]

def _is_user_module(module, library_roots):
    """ A user module is a pure python module that does not come from the standard
    library or from an installed package. Compiled extensions are never considered
    as user modules since they cannot be reloaded safely within a process. """
    file = getattr(module, '__file__', None)
    if file and file.endswith(tuple(importlib.machinery.EXTENSION_SUFFIXES)):
        return False
    locations = _module_locations(module)
    if not locations:  # built-in or frozen module
        return False
    return not any(_is_within(location, root) for location in locations for root in library_roots)

def _purge_new_user_modules(modules_before):
    """ Remove from sys.modules the user modules that were imported since the snapshot
    `modules_before` (a set of module names), so that a later load re-imports them
    from its own folder. Only new top-level names are considered, which protects the
    submodules that libraries (e.g. xpsi) import lazily. Modifies sys.modules in place. """
    tops_before = {name.partition('.')[0] for name in modules_before}
    roots = _library_roots()
    for name in list(sys.modules):
        if name in modules_before or name.partition('.')[0] in tops_before:
            continue
        if _is_user_module(sys.modules[name], roots):
            del sys.modules[name]


#####################
#### MODEL LOADER ####
#####################
def load_model( model_path , 
               config_path = None , 
               execution_path = None ):
    """ Function to load a model from a given python script. 
    The model can be loaded with or without a configuration file, from an execution directory or locally.

    Several models can be loaded one after the other (and be used simultaneously) even if their
    custom modules share the same names: the modules imported by a model from its own folders are
    removed from the module cache after the load, so the next model imports its own versions.
    The models already loaded keep working with their own classes.

    Args:
        model_path (str): Path to the main file which loads the model.
        config_path (str|None, optional): Path to the configuration file, if needed, for loading the main. No configuration file is used if None. Defaults to None.
        execution_path (str|None, optional): Path to the execution directory, if needed, for loading the main. Execution happens in the current directory if None. Defaults to None.

    Returns:
        Namespace: Imported main with a X-PSI defined Neutron Star model. The different attributes are defined in the main.
    """

    # Check that the files exist, and resolve everything before any change of directory
    if not os.path.isfile(model_path):
        raise ModelLoaderError(f'File {model_path} does not exist.')
    abs_model_path = os.path.abspath(model_path)

    abs_config_path = None
    if config_path is not None:
        if not os.path.isfile(config_path):
            raise ModelLoaderError(f'File {config_path} does not exist.')
        abs_config_path = os.path.abspath(config_path)

    if execution_path is not None and not os.path.isdir(execution_path):
        raise ModelLoaderError(f'Directory {execution_path} does not exist.')

    # Save the original state. Lists and dicts modified by the main are restored in place.
    pwd = os.getcwd()
    original_sys_path = list(sys.path)
    original_sys_argv = sys.argv
    modules_before = set(sys.modules)

    # Try to catch the error and make sure the cleanup happens
    try:

        # Execution happens in the execution directory (relative paths of the main are relative to it)
        if execution_path is not None:
            os.chdir(execution_path)

        # Make sure the finders do not serve stale directory listings from a previous model
        importlib.invalidate_caches()

        # The main reads its configuration file from the command line
        sys.argv = [abs_model_path] + ([f'@{abs_config_path}'] if abs_config_path is not None else [])

        # Load the main directly from its real path (no symbolic links needed)
        loader = importlib.machinery.SourceFileLoader(_MODEL_MODULE_NAME, abs_model_path)
        spec = importlib.util.spec_from_file_location(_MODEL_MODULE_NAME, abs_model_path, loader=loader)
        model = importlib.util.module_from_spec(spec)
        sys.modules[_MODEL_MODULE_NAME] = model

        # Do the loading
        spec.loader.exec_module(model)

    except NameError as e:
        raise ModelLoaderError("The model could not be loaded because of an ill-defined import in the model at model_path. Look the traceback for more information.") from e

    except (Exception, SystemExit) as e:
        # SystemExit is included since argparse exits on a bad configuration file
        raise ModelLoaderError(f"The model could not be loaded ({type(e).__name__}: {e}). Look the traceback for more information.") from e

    # Cleanup whatever happens
    finally:
        # Purge before leaving the execution directory, in case some __file__ are relative to it
        sys.modules.pop(_MODEL_MODULE_NAME, None)
        _purge_new_user_modules(modules_before)
        importlib.invalidate_caches()
        os.chdir( pwd )
        sys.path[:] = original_sys_path
        sys.argv = original_sys_argv

    # Successful
    print( f'Model located at {model_path}' + (f' with configuration {config_path}' if config_path is not None else '') + ' has been sucessfully loaded' )
    # Return the model
    return model
