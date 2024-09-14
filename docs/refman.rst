API Reference
=============

.. toctree::
    :hidden:
    :maxdepth: 4
    :caption: Python package

    apidocs/sansmic

.. toctree::
    :hidden:
    :maxdepth: 4
    :caption: C++ library

    apidocs/index

There are two parts of the sansmic package, the python package and the
C++ library.

The python modules are documented in typical python fashion using sphinx
processing of the docstrings. The pybind11 module creates the python
wrapper functions for the :class:`sansmic.libsansmic` binary module;
these are also documented with the native python modules.

The C++ functions themselves are documented separately from the python code.
The documentation is in doxygen style, passed through breathe and exhale
prior to being processed by sphinx, which produces this documentation.
This leads to a different type of structure than sphinx normally produces,
but one which will be familiar to doxygen users.


.. HACK to create autosummary
    .. autosummary::
        :toctree: apidocs
        :recursive:
        :template: autosummary/module.rst

        sansmic
