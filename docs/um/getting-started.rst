Getting started
===============
Please note that sansmic does not have a graphical user interface (:term:`GUI`). This means
that to run sansmic, you will need a command prompt or terminal of some sort,
or you can run sansmic from within python or a Jupyter notebook. (Setting up
Jupyter is beyond the scope of this guide, please see https://jupyter.org/).

.. note::

   Based on feedback from potential users and collaborators, there is a now standalone
   executable that is built that will run on modern Windows[#]_ operating systems without
   needing to install Python at the user- or system-level. See `Standalone installation for Windows-users`_,
   below.


Standalone installation for Windows-users
-----------------------------------------
For Windows users (only), a standalone executable has been produced. This executable
is available on GitHub releases (https://github.com/sandialabs/sansmic/releases/latest)
and will be at the bottom of the page with a name like
``sansmic-v1.0.7-standalone-win_amd64.zip``.
There is also a ``[...].sigstore.json`` file that can be used to verify the zip-file.

Once the zip-file is downloaded, unzip the file (right-click, select "Extract All ...")
to a new folder. There will be a folder called "sansmic" that will contain the
sansmic executable. You will need to put this "sansmic" folder somewhere and add it
to your "PATH" environment variable. To do this:

* Click on the magnifying glass "Search" icon (it is possible that clicking the Windows button
  will also bring up the search function).
* Type "edit environment variables for your account" and click on the result that comes up
* Double-click on "Path" in the *top* section of the window that appears.
* Click "New" and type the location of the "sansmic" folder (e.g., "C:\\Users\\foo\\sansmic"
  if my username is "foo" and I dragged the "sansmic" folder to my "foo" icon in the Explorer
  sidebar).
* Click "Okay" to close each of the windows that are open. Sansmic should now be available
  for you to run from the Terminal app!

To run the examples, open the Terminal app (you can use the search function for this, too)
and, assuming a user named "foo" once again, type:

.. code:: bash

   sansmic --version


This should give you the version number (and it should match the version you downloaded!).
To run an example, continue with:

.. code:: bash

   cd C:\Users\foo\sansmic\examples
   sansmic -v baseline.toml


This should show a progress bar and create files in the examples directory. You are
ready to start using sansmic on your own problems!



Python "pip install" (preferred)
--------------------------------
Sansmic is built and tested with CPython 3.9--3.12 on 64-bit
Installing and setting up Python is outside the scope of this User Guide, but
https://www.python.org/ is the official site for the software, and it can also
be installed for free from the Windows Store. For Mac and Linux users, the authors
recommend following the appropriate instructions on the Python main site if you
do not already have Python installed.

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
~~~~~~~~~~~~~~
Rather than repeat the gory details of setting up sansmic for development,
please see the ``CONTRIBUTING.md`` file in the git repository
(https://github.com/sandialabs/sansmic).



.. only:: html

   .. rubric:: Notes

.. [#] Windows is a registered trademark of Microsoft Corp.
.. [#] There are two packages that are always required to run sansmic, numpy and pandas.
   If your Python version is less than 3.11, sansmic will also require - and pip will
   automatically install - the tomli package for toml support, which is already included
   in Python v3.11+.
