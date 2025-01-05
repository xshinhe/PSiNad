from __future__ import absolute_import
__author__ = "Xin He"

import os, os.path
import sys
import pkgutil
import importlib
from . import version

# import pyPSiNad.libpyPSiNad_v1 as lib
import pyPSiNad.libpyPSiNad_v2 as lib

# 合并命名空间
for attr in dir(lib):
    # print(attr)
    if not attr.startswith('_'):
        setattr(sys.modules[__name__], attr, getattr(lib, attr))

# 删除对 lib 和 libpyPSiNad_v2 的引用
modules_to_delete = ['pyPSiNad.lib', 'pyPSiNad.libpyPSiNad_v2']
for module_name in modules_to_delete:
    if module_name in sys.modules:
        del sys.modules[module_name]

# 可选：删除当前作用域中的引用
del lib

__all__ = []

# 动态导入所有子模块
for loader, module_name, is_pkg in pkgutil.walk_packages(__path__):
    if not module_name.startswith('_') and module_name not in dir():
        __all__.append(module_name)
        # print(module_name)
        importlib.import_module(f".{module_name}", package=__name__)

def onlySampling(M: Model) -> Solver:
    return defaultSolverFactory('Sampling', M)

# __version__ = Platform.KIDSVersion()

class KIDSException(Exception): # for swig wrapper only
    """This is the class used for all exceptions thrown by the C++ library."""
    pass
