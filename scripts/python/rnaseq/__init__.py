import importlib
from . import readwrite
importlib.reload(readwrite)
from . import setup
importlib.reload(setup)
from . import pp
importlib.reload(pp)
from . import dea
importlib.reload(dea)
from . import utils
importlib.reload(utils)
from . import plotting
importlib.reload(plotting)


# from os.path import dirname, basename, isfile, join
# import glob
# modules = glob.glob(join(dirname(__file__), "*.py"))
# __all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]



# __all__ = ["readwrite", "rnadata", "annot", "pp", "dea", "utils"]