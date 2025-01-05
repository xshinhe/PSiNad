import os
import shutil
import sys
import functools
import warnings
import platform
import numpy
from glob import glob
from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension
from Cython.Build import cythonize

MAJOR_VERSION_NUM = '@PSINAD_MAJOR_VERSION@'
MINOR_VERSION_NUM = '@PSINAD_MINOR_VERSION@'
BUILD_INFO = '@PSINAD_BUILD_VERSION@'
GIT_VERSION = '@PSINAD_GIT_VERSION@'
cmake_source_dir = "@CMAKE_SOURCE_DIR@"
IS_RELEASED = False

def not_recommended(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        warnings.warn(f"Function {func.__name__} is not recommended for use.", category=UserWarning)
        return func(*args, **kwargs)
    return wrapper

@not_recommended
def remove_package(mod, verbose):
    try:
        install_path = mod.__path__[0]
        if os.path.exists(install_path) and verbose:
            print(f'REMOVING "{install_path}"')
            shutil.rmtree(install_path)
    except (AttributeError, IndexError, OSError):
        pass

@not_recommended
def uninstall(verbose=True):
    try:
        import pyPSiNad
        remove_package(pyPSiNad, verbose)
    except ImportError:
        pass

def report_error(message):
    sys.stderr.write("ERROR: {}\n".format(message))
    sys.exit(1)

def write_version_py():
    GIT_REVISION = GIT_VERSION
    VERSION = FULL_VERSION = f"{MAJOR_VERSION_NUM}.{MINOR_VERSION_NUM}.{BUILD_INFO}"
    if not IS_RELEASED:
        FULL_VERSION += f'.dev-{GIT_REVISION[:7]}'

    with open('pyPSiNad/version.py', 'w') as f:
        f.write(f"""
# THIS FILE IS GENERATED FROM SETUP.PY
short_version = '{VERSION}'
version = '{VERSION}'
full_version = '{FULL_VERSION}'
git_revision = '{GIT_REVISION}'
release = {str(IS_RELEASED)}
PSINAD_library_path = r'{os.getenv('PSINAD_LIB_PATH')}'

if not release:
    version = full_version
""")

def build_setup_kwargs():
    if sys.version_info < (2, 7):
        report_error("PSINAD requires Python 2.7 or better.")
    if platform.system() in ['Windows']:
        report_error(f"Not supported on platform {platform.system()}")

    pyPSiNad_include_path = os.getenv('PSINAD_INCLUDE_PATH')
    if not pyPSiNad_include_path:
        report_error("Set PSINAD_INCLUDE_PATH to point to the include directory for PSINAD")
    pyPSiNad_lib_path = os.getenv('PSINAD_LIB_PATH')
    if not pyPSiNad_lib_path:
        report_error("Set PSINAD_LIB_PATH to point to the lib directory for PSINAD")
    library_dirs=[pyPSiNad_lib_path]

    include_dirs=pyPSiNad_include_path.split(';') + [numpy.get_include()]

    # 指定自带的 pybind11 路径
    print(cmake_source_dir)
    pybind11_path = os.path.join(cmake_source_dir, "thirdpart", "pybind11", "include")
    print(pybind11_path)
    if not os.path.exists(pybind11_path):
        raise RuntimeError(f"pybind11 not found at {pybind11_path}")
    # 指定自带的 Eigen 路径
    eigen_path = os.path.join(cmake_source_dir, "thirdpart", "Eigen")
    if not os.path.exists(eigen_path):
        raise RuntimeError(f"Eigen not found at {eigen_path}")

    include_dirs.extend([pybind11_path, eigen_path])

    define_macros = [('MAJOR_VERSION', MAJOR_VERSION_NUM),
                     ('MINOR_VERSION', MINOR_VERSION_NUM)]

    libraries = ['PSiNad_shared'] #, 'PSINADPlugin']

    extension_args = {
        "name": "pyPSiNad._psinad",
        "sources": ["swig/PSINADSwig.cxx"],
        "include_dirs": include_dirs,
        "define_macros": define_macros,
        "library_dirs": library_dirs,
        "libraries": libraries,
        "extra_compile_args": ['-std=c++17', '-Wno-c++11-extensions'],
        "extra_link_args": ['-std=c++17'],
        "runtime_library_dirs":library_dirs,
    }

    setup_kwargs = {
        "name": "pyPSINAD",
        "version": f"{MAJOR_VERSION_NUM}.{MINOR_VERSION_NUM}.{BUILD_INFO}",
        "author": "Liu-group  (Xin He et al.)",
        "author_email": "xshinhe@pku.edu.cn",
        "license": "Python Software Foundation License (BSD-like)",
        "url": "https://github.com/liugroup/PSiNad",
        "packages": [
            "pyPSiNad",
        ],
        "package_data": {
            "pyPSiNad": ["libpyPSiNad_v1.so"],
        },
        "platforms": ["Linux"],
        "description": "Python wrapper for PSiNad (Phase Space Integrated Nonadiabatic Dynamics)",
        "long_description": """
        PSINAD (Phase Space Integrated Nonadiabatic Dynamics) offers an open-source framework tailored for simulating 
        nonadiabatic dynamics based quantum phase space and advanced trajectory-based approximations.
        """,
        "ext_modules" : [
            Extension(**extension_args),
            Pybind11Extension(
                "pyPSiNad.libpyPSiNad_v2",
                sorted(glob("pybind11/libpyPSiNad_v2.cpp")),
                include_dirs=include_dirs,
                library_dirs=library_dirs,
                libraries=libraries,
                extra_compile_args=['-std=c++17'],
                extra_link_args=['-lPSiNad_shared', '-Wl,-rpath,/usr/local/psinad/lib'],
            ),
            Extension(
                "pyPSiNad.ext.examples",
                sources=[
                    "ext/examples/test.cpp",
                    # "ext/examples/test.pyx",  # Cython
                ],
                include_dirs=include_dirs + [
                    "ext/examples/include",
                    numpy.get_include(),
                ],
                language="c++",
            ),
        ],
        "install_requires": [
            "cython>=0.29.0",  # 需要 Cython 0.29.0 或更高版本
            "numpy>=1.20.0",   # 需要 NumPy 1.20.0 或更高版本
        ],
    }

    return setup_kwargs    

if __name__ == '__main__':
    setupKeywords=build_setup_kwargs()
    write_version_py()
    setup(**setupKeywords)
