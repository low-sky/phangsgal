# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
# ----------------------------------------------------------------------------

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys
from distutils.version import LooseVersion

__minimum_python_version__ = "3.6"

__all__ = []

class UnsupportedPythonError(Exception):
    pass


if LooseVersion(sys.version) < LooseVersion(__minimum_python_version__):
    raise UnsupportedPythonError("phangs-alma does not support Python < {}"
                                 .format(__minimum_python_version__))

if not _ASTROPY_SETUP_:   # noqa
    # For egg_info test builds to pass, put package imports here.
    from .phangs_galaxies import *
    from .utilities import *
    from .XCO import *
    # Then you can be explicit to control what ends up in the namespace,
    # __all__ += ['do_primes']   # noqa
    # or you can keep everything from the subpackage with the following instead
    __all__ += phangs_galaxies.__all__
    __all__ += utilities.__all__
    __all__ += XCO.__all__
    
    
