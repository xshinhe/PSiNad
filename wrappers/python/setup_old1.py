import ast
import re
import os
import sys
import platform
import numpy
from glob import glob
from setuptools import setup
from Cython.Build import cythonize

MAJOR_VERSION_NUM='@KIDS_MAJOR_VERSION@'
MINOR_VERSION_NUM='@KIDS_MINOR_VERSION@'
BUILD_INFO='@KIDS_BUILD_VERSION@'
GIT_VERSION='@KIDS_GIT_VERSION@'
IS_RELEASED=False

__author__ = "Xin He"
__version__ = "%s.%s" % (MAJOR_VERSION_NUM, MINOR_VERSION_NUM)

def reportError(message):
    sys.stdout.write("ERROR: ")
    sys.stdout.write(message)
    sys.stdout.write("\nExiting\n")
    sys.exit(1)

def removeRecursive(dir):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if os.path.isdir(path):
            removeRecursive(path)
        else:
            os.remove(path)
    os.rmdir(dir)

def removePackage(mod, verbose):
        try:
            pathList = mod.__path__
        except AttributeError:
            return
        if len(pathList) > 1:
           raise Exception("more than one item in KIDS.__path__")
        installPath = pathList[0]
        if os.path.exists(installPath):
            if verbose:
                sys.stdout.write('REMOVING "%s"\n' % installPath)
            removeRecursive(installPath)

def uninstall(verbose=True):
    save_path=sys.path[:]
    sys.path=[]
    for item in save_path:
        if item!='.' and item!=os.getcwd():
            sys.path.append(item)
    try:
        import pykids
        removePackage(pykids, verbose)
    except ImportError:
        pass
    sys.path=save_path

def writeVersionPy(
    filename="pykids/version.py", # @deprecated from cmake config other than from file
    major_version_num=MAJOR_VERSION_NUM,
    minor_version_num=MINOR_VERSION_NUM, 
    build_info=BUILD_INFO):
    cnt = """
# THIS FILE IS GENERATED FROM KIDS SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s
KIDS_library_path = r'%(path)s'

if not release:
    version = full_version
"""
    ''' @deprecated
    if os.path.exists(filename):
        with open(filename) as f:
            text = f.read()
            match = re.search(r"git_revision\s+=\s+(.*)", text, re.MULTILINE)
        try:
            git_revision = ast.literal_eval(match.group(1))
        except:
            git_revision = 'Unknown'
    else:
        git_revision = 'Unknown'
    '''
    git_revision = GIT_VERSION
    version = full_version = '%s.%s.%s' % (major_version_num, minor_version_num, build_info)
    if not IS_RELEASED:
        full_version += '.dev-' + git_revision[:7]
    with open(filename, 'w') as a:
        a.write(cnt % {'version': version,
                       'full_version' : full_version,
                       'git_revision' : git_revision,
                       'isrelease': str(IS_RELEASED),
                       'path': os.getenv('KIDS_LIB_PATH')})

def buildKeywordDictionary(major_version_num=MAJOR_VERSION_NUM,
                           minor_version_num=MINOR_VERSION_NUM,
                           build_info=BUILD_INFO):
    from setuptools import Extension
    try:
        from pybind11.setup_helpers import Pybind11Extension
        print('%'*100)
    except ImportError:
        from setuptools import Extension as Pybind11Extension

    if platform.system() == "Windows" or platform.system() == 'Darwin':
        reportError("Not support now for platform %s"%platform.system())

    setupKeywords = {}
    setupKeywords["name"]              = "KIDS"
    setupKeywords["version"]           = "%s.%s.%s" % (major_version_num,
                                                       minor_version_num,
                                                       build_info)
    setupKeywords["author"]            = "Xin He | Liu-group"
    setupKeywords["license"]           = "Python Software Foundation License (BSD-like)"
    setupKeywords["url"]               = "https://github.com/xshinhe/KIDS"
    setupKeywords["packages"]          = [
                                          "pykids",
                                          "pykids.ext",
                                          # "pykids.chem",
                                          # "pykids.core",
                                          # "pykids.models",
                                          # "pykids.solvers",
                                          # "pykids.app"
                                          ]
    setupKeywords["data_files"]        = []
    setupKeywords["package_data"]      = {"pykids" : [
                                            "libpykids_v1.so",
                                            ],
                                          # "pykids.app" : ['data/*.json', 'data/*.dat', 'data/*.ds'],
                                          }
    setupKeywords["platforms"]         = ["Linux"] #, "Mac OS X", "Windows"]
    setupKeywords["description"]       = "Python wrapper for KIDS"
    setupKeywords["long_description"]  = \
    """
    KIDS (Kernel Integrated Dynamics Simulator) offers an open-source framework tailored for simulating 
    chemical and physical dynamics, with a primary focus on atomic and molecular scales in condensed matter. 
    It is designed for (classical / qauntum) dynamics simulation of small system, few-body system, reduced 
    systems, and even large many-particle (i.e. molecules & condensed matter) systems. It provides a versatile 
    platform for the development of advanced algorithms, offering ease of use and accessibility at a minimal 
    cost.
    """

    define_macros = [('MAJOR_VERSION', major_version_num),
                     ('MINOR_VERSION', minor_version_num)]

    libraries=['KIDS'] #, 'KIDSPlugins']

    pykids_include_path = os.getenv('KIDS_INCLUDE_PATH')
    if not pykids_include_path:
        reportError("Set KIDS_INCLUDE_PATH to point to the include directory for KIDS")
    pykids_lib_path = os.getenv('KIDS_LIB_PATH')
    if not pykids_lib_path:
        reportError("Set KIDS_LIB_PATH to point to the lib directory for KIDS")


    ## ADD SWIG EXTENSION
    extra_compile_args=['-std=c++11']
    extra_link_args=[] # fix for crossing platform @todo
    library_dirs=[pykids_lib_path]
    include_dirs=pykids_include_path.split(';')
    include_dirs.append(numpy.get_include())
    extensionArgs = {"name": "pykids._kids", # swig library
                    "sources": ["swig/KIDSSwig.cxx"],
                    "include_dirs": include_dirs,
                    "define_macros": define_macros,
                    "library_dirs": library_dirs,
                    "libraries": libraries,
                    "extra_compile_args": extra_compile_args,
                    "extra_link_args": extra_link_args}
    extensionArgs["runtime_library_dirs"] = library_dirs
    setupKeywords["ext_modules"] = [Extension(**extensionArgs)]

    ## ADD PYBIND11 EXTENSION # pybind11 library
    setupKeywords["ext_modules"] += [Pybind11Extension(
                                            "pykids.libpykids_v2",
                                            sorted(glob("pybind11/libpykids_v2.cpp")),
                                            # include_dirs=include_dirs,
                                            # library_dirs=library_dirs,
                                        ),
                                    ]

    ## OTHER EXTENSION                                    
    # setupKeywords["ext_modules"] += cythonize('ext/*.pyx')
    setupKeywords["ext_modules"] += cythonize(Extension(
        "pykids.ext.examples",
        sources=[
            "ext/examples/src/test.cpp",
            "ext/examples/test.pyx",
        ],
        include_dirs=include_dirs +[
            "ext/examples/include",
            "ext/examples/",
            numpy.get_include(),
        ],
        language="c++",
    ))

    outputString = ''
    firstTab     = 40
    secondTab    = 60
    for key in sorted(iter(setupKeywords)):
         value         = setupKeywords[key]
         outputString += key.rjust(firstTab) + str( value ).rjust(secondTab) + "\n"
    # sys.stdout.write("%s" % outputString)
    return setupKeywords

def main():
    if sys.version_info < (2, 7):
        reportError("KIDS requires Python 2.7 or better.")
    if platform.system() == 'Darwin' or platform.system() == 'Windows':
        reportError("Not support now for platform %s"%platform.system())
    try:
        uninstall()
    except:
        pass
    setupKeywords=buildKeywordDictionary()
    writeVersionPy()
    print(setupKeywords)
    setup(**setupKeywords)

if __name__ == '__main__':
    main()
