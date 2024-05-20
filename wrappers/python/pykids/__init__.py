from __future__ import absolute_import
__author__ = "Xin He"

import os, os.path
import sys
from . import version

# if sys.platform == 'win32':
#     _path = os.environ['PATH']
#     os.environ['PATH'] = '%(lib)s;%(lib)s\plugins;%(path)s' % {
#         'lib': version.pykids_library_path, 'path': _path}
#     try:
#         with os.add_dll_directory(version.pykids_library_path):
#             from . import _pykids
#     except:
#         pass

os.path.append(version.KIDS_library_path)
# import pykids.libpykids_v1 as v1
import pykids.libpykids_v2 as v2
import pykids._kids as _kids
import pykids.ext as ext

if sys.platform == 'win32':
    os.environ['PATH'] = _path
    del _path

# __version__ = Platform.KIDSVersion()

class KIDSException(Exception): # for swig wrapper only
    """This is the class used for all exceptions thrown by the C++ library."""
    pass
