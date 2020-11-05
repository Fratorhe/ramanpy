from ._version import get_versions
from .specific_fit_classes import *

__version__ = get_versions()['version']
del get_versions

from .runners import *
