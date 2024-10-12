"""Setup SANSMIC"""

# Always prefer setuptools over distutils
from os import path
from glob import glob

from setuptools import find_packages, setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

__version__ = "1.0.1"

here = path.abspath(path.dirname(__file__))

ext_modules = [
    Pybind11Extension(
        "sansmic.libsansmic",
        sources=sorted(glob("src/ext_modules/libsansmic/*.cpp")),
        cxx_std=17,
        # Example: passing in the version to the compiled code
        define_macros=[
            ("VERSION_INFO", '"{}"'.format(__version__)),
        ],
    ),
]

setup(
    name="sansmic",  # Required
    version=__version__,  # Required
    description="The SANSMIC solution mining code",  # Optional
    # long_description=long_description,  # Optional
    # long_description_content_type='text/markdown',  # Optional (see note above)
    packages=[
        "sansmic",
    ],  # Required
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
