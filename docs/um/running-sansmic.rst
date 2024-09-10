Running sansmic
===============
There are two commands that can be run from the command line, sansmic_ and sansmic-convert_.
The following include all available output options, but note that not all options may be
available if you have not installed the appropriate packages (see Advanced Installation).


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
