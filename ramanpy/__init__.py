from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

from .specific_fit_classes import *
from ramanpy import runners
from .process_results_raman import *