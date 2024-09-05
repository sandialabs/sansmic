Getting started
===============
Please note, sansmic does not have a graphical user interface (GUI). This means
that to run sansmic, you will need a command prompt or terminal of some sort,
or you can run sansmic from within python or jupyter/jupyterlab.



Installing sansmic
------------------
To install sansmic, you will need to have Python installed. Installing and setting up
Python is outside the scope of this User Guide, but https://www.python.org/ is the
official site for the software and it can also be installed for free from the Windows
Store. For Mac users, please find some good instructions online. If you use Linux or
other Unix-like systems, the authors are going to assume you already have and use
python on your system.

Once python is installed, open up a Command Prompt, Console, or xterm window, and install
sansmic using the ``pip`` command.

.. code:: bash

    python -m pip install sansmic


This will download and install the latest version of sansmic and the required packages
for sansmic to your python installation.


Running sansmic
---------------
There are two commands that can be run from the command line, sansmic_ and sansmic-convert_.


``sansmic``
-----------

.. argparse::
   :module: sansmic.app
   :func: _main_parser
   :prog: sansmic


``sansmic-convert``
-------------------

.. argparse::
   :module: sansmic.app
   :func: _convert_parser
   :prog: sansmic-convert


For developers
--------------
The sansmic package is a mixture of C++ and python code, so you will need to have a
C++ compiler installed if you plan on modifying any of the C++ code. If not, then
you can get away without having one and using the distributed binary library with
the wheel.
