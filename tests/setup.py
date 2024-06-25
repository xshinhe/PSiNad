from distutils.core import setup
from Cython.Build import cythonize
setup(
      ext_modules = cythonize("Kernel_Customized.py",language_level = "3")
)