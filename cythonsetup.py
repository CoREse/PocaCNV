from setuptools import Extension, setup
from Cython.Build import cythonize

# setup.py
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "rdprocessing",
        ["rdprocessing.pyx","getRDScores.cpp"],
        #libraries=["CGetRDScores"],
        extra_compile_args=['-fopenmp', '-O3', "-std=c++11"],
        extra_link_args=['-fopenmp'],
        language="c++"
    )
]

setup(
    name='rdprocessing',
    ext_modules=cythonize(ext_modules,language_level="3"),
)
