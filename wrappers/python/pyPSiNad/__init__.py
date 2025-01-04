from __future__ import absolute_import
__author__ = "Xin He"

import os, os.path
import sys
from . import version

# if sys.platform == 'win32':
#     _path = os.environ['PATH']
#     os.environ['PATH'] = '%(lib)s;%(lib)s\plugins;%(path)s' % {
#         'lib': version.pyPSiNad_library_path, 'path': _path}
#     try:
#         with os.add_dll_directory(version.pyPSiNad_library_path):
#             from . import _pyPSiNad
#     except:
#         pass

# os.path.append(version.KIDS_library_path)

# import pyPSiNad.libpyPSiNad_v1 as lib
import pyPSiNad.libpyPSiNad_v2 as lib
# import pyPSiNad._psinad as _psinad
import pyPSiNad.ext as ext
import pyPSiNad.example as example

# __version__ = Platform.KIDSVersion()

class KIDSException(Exception): # for swig wrapper only
    """This is the class used for all exceptions thrown by the C++ library."""
    pass
