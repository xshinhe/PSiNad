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

MAJOR_VERSION_NUM = '@KIDS_MAJOR_VERSION@'
MINOR_VERSION_NUM = '@KIDS_MINOR_VERSION@'
BUILD_INFO = '@KIDS_BUILD_VERSION@'
GIT_VERSION = '@KIDS_GIT_VERSION@'
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
        import pykids
        remove_package(pykids, verbose)
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

    with open('pykids/version.py', 'w') as f:
        f.write(f"""
# THIS FILE IS GENERATED FROM SETUP.PY
short_version = '{VERSION}'
version = '{VERSION}'
full_version = '{FULL_VERSION}'
git_revision = '{GIT_REVISION}'
release = {str(IS_RELEASED)}
KIDS_library_path = r'{os.getenv('KIDS_LIB_PATH')}'

if not release:
    version = full_version
""")

def build_setup_kwargs():
    if sys.version_info < (2, 7):
        report_error("KIDS requires Python 2.7 or better.")
    if platform.system() in ['Darwin', 'Windows']:
        report_error(f"Not supported on platform {platform.system()}")

    pykids_include_path = os.getenv('KIDS_INCLUDE_PATH')
    if not pykids_include_path:
        report_error("Set KIDS_INCLUDE_PATH to point to the include directory for KIDS")
    pykids_lib_path = os.getenv('KIDS_LIB_PATH')
    if not pykids_lib_path:
        report_error("Set KIDS_LIB_PATH to point to the lib directory for KIDS")
    library_dirs=[pykids_lib_path]
    include_dirs=pykids_include_path.split(';') + [numpy.get_include()]

    define_macros = [('MAJOR_VERSION', MAJOR_VERSION_NUM),
                     ('MINOR_VERSION', MINOR_VERSION_NUM)]

    libraries = ['KIDS'] #, 'KIDSPlugin']

    extension_args = {
        "name": "pykids._kids",
        "sources": ["swig/KIDSSwig.cxx"],
        "include_dirs": include_dirs,
        "define_macros": define_macros,
        "library_dirs": library_dirs,
        "libraries": libraries,
        "extra_compile_args": ['-std=c++11'],
        "extra_link_args": [],
        "runtime_library_dirs":library_dirs,
    }

    print(f"library_dirs={library_dirs}")
    print(f"include_dirs={include_dirs}")
    print(f"libraries={libraries}")

    setup_kwargs = {
        "name": "KIDS",
        "version": f"{MAJOR_VERSION_NUM}.{MINOR_VERSION_NUM}.{BUILD_INFO}",
        "author": "Xin He | Liu-group",
        "license": "Python Software Foundation License (BSD-like)",
        "url": "https://github.com/xshinhe/KIDS",
        "packages": [
            "pykids",
        ],
        "package_data": {
            "pykids": ["libpykids_v1.so"],
        },
        "platforms": ["Linux"],
        "description": "Python wrapper for KIDS",
        "long_description": """
        KIDS (Kernel Integrated Dynamics Simulator) offers an open-source framework tailored for simulating 
        chemical and physical dynamics, with a primary focus on atomic and molecular scales in condensed matter. 
        It is designed for (classical / qauntum) dynamics simulation of small system, few-body system, reduced 
        systems, and even large many-particle (i.e. molecules &amp; condensed matter) systems. It provides a versatile 
        platform for the development of advanced algorithms, offering ease of use and accessibility at a minimal 
        cost.
        """,
        "ext_modules" : [
            Extension(**extension_args),
            Pybind11Extension(
                "pykids.libpykids_v2",
                sorted(glob("pybind11/libpykids_v2.cpp")),
                include_dirs=include_dirs,
                library_dirs=library_dirs,
                libraries=libraries,
                extra_compile_args=['-lKIDS', '-Wl,-rpath,/usr/local/kids/lib'],
                extra_link_args=['-Wl,-R/usr/local/kids/lib'],
            ),
            Extension(
                "pykids.ext.examples",
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
    }

    return setup_kwargs    

if __name__ == '__main__':
    setupKeywords=build_setup_kwargs()
    write_version_py()
    setup(**setupKeywords)
