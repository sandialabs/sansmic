Getting started
===============
Please note, sansmic does not have a graphical user interface (:term:`GUI`). This means
that to run sansmic, you will need a command prompt or terminal of some sort,
or you can run sansmic from within python or a Jupyter notebook. (Setting up
Jupyter is beyond the scope of this guide, please see https://jupyter.org/).


Installing sansmic
------------------
To install sansmic, you will need to have Python installed. Installing and setting up
Python is outside the scope of this User Guide, but https://www.python.org/ is the
official site for the software and it can also be installed for free from the Windows
Store. For Mac users, please find some good instructions online. If you use Linux or
other Unix-like systems, the authors are going to assume you already have and use
python on your system.

Once python is installed, open up a Command Prompt, Console, or xterm window, and install
sansmic using the ``pip`` command. This will download and install the latest version of sansmic
and the two (or three) required packages for sansmic to your python installation. [#]_

.. code:: bash

    python3 -m pip install sansmic


That's it! You can test if sansmic installed property by trying the command

.. code:: bash

   sansmic --help


If the installation was successful, it should print out a help screen, and you can
start running sansmic. By default, sansmic will output results to Comma Separated Values
(:term:`CSV`) formatted files, that are readable by any data processing or spreadsheet
software, or :term:`JSON` files, which can be read in by numpy and pandas, among other
packages. (If this is all you want to do, you can skip forward to the next chapter!)


Advanced installation options
-----------------------------
Because python is a modular programming language, not every installation has, or wants
to have, every python package installed that sansmic *could* make use of. With input
and output formats, specifically, some users may want some functionality and not others.


Input and output formats
~~~~~~~~~~~~~~~~~~~~~~~~
The following packages make additional file formats available to sansmic. You can
install all of them by installing using pip extras notation.

.. code:: bash

   python3 -m pip install sansmic[formats]


* either ``pyyaml`` or ``ruamel.yaml`` - use :term:`YAML` (.yaml) formatted scenario configuration files.
* ``h5py`` - use :term:`HDF5` (.h5) data files for saving results.
* ``openpyxl`` - save results in Microsoft Excel (.xlsx) formatted files.
* ``tabulate`` - print results to the screen in Markdown (.md) format, which is prettier.
* ``lasio`` - use geometry data from a :term:`LAS` (.las) formatted text file.


Building sansmic
~~~~~~~~~~~~~~~~
If you are building sansmic from source (a likelihood for \*nix and some Mac users),
you will need the following packages, but pip *should* download them for you. You
also need a C++ compiler and python header libraries.

* ``setuptools`` - provides the build backend for the package.
* ``pybind11`` - required header libraries that make the C++ libsansmic library compilable.


For developers
~~~~~~~~~~~~~~~
Rather than repeat the gory details of setting up sansmic for development,
please see the ``CONTRIBUTING.md`` file in the git repository
(https://github.com/sandialabs/sansmic).


Building documentation
~~~~~~~~~~~~~~~~~~~~~~
You have already found the documentation since you are reading this! But, if you want to
make it available offline, or in another format, there are some additional packages you
will need, and you will need the **Doxygen** program (https://www.doxygen.nl).
After you get Doxygen, you can get the necessary python packages by using the "docs"
pip extras marker.

.. code:: bash

   python3 -m pip install sansmic[docs]


The libraries this will install (if you don't already have them) are:

* ``Sphinx`` - the main documentation engine sansmic.
* ``sphinx_design`` - layout options package, required by pydata-sphinx-theme.
* ``pydata-sphinx-theme`` - the theme used for the sansmic documentation.
* ``sphinxcontrib-bibtex`` - use better bibliographic citations.
* ``sphinx-argparse`` - generate program help automatically.
* ``breathe`` - process doxygen output files.
* ``exhale`` - automate the doxygen+breathe process.


.. only:: html

   .. rubric:: Notes

.. [#] There are two packages that are always required to run sansmic, numpy and pandas.
   If your Python version is less than 3.11, sansmic will also require - and pip will
   automatically install - the tomli package for toml support, which is already included
   in Python v3.11+.
