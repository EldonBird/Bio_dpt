from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "pcr_lib",                          # name of the generated Python module
        ["src/pcr_lib.cpp"],               # your C++ source file
        include_dirs=[pybind11.get_include()],
        language="c++"
    )
]

setup(
    name="pcr_lib",
    version="0.1",
    ext_modules=ext_modules,
)
