import importlib
import sys
import inspect
import functools

def cytest(func):
    """
    Wraps `func` in a plain Python function.
    """

    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        bound = inspect.signature(func).bind(*args, **kwargs)
        return func(*bound.args, **bound.kwargs)

    return wrapped

# list of Cython modules containing tests
cython_test_modules = ["mytest_hot"]

for mod in cython_test_modules:
    try:
        # For each callable in `mod` with name `test_*`,
        # set the result as an attribute of this module.
        mod = importlib.import_module(mod)
        for name in dir(mod):
            item = getattr(mod, name)
            if callable(item) and name.startswith("test_"):
                setattr(sys.modules[__name__], name, item)
    except ImportError:
        pass