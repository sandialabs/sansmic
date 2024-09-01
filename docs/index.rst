=======
sansmic
=======

The sansmic leaching simulator
==============================

.. toctree::
   :hidden:
   :maxdepth: 1

   userman
   refman


**sansmic** is a package for the simulation of leaching in underground salt
caverns. The core is written in C++ for speed, but the primary API and the
command-line program are Python-based. Lower case *sansmic* is used to
differentiate it from the previous, FORTRAN-based versions of the same name.

Solution mining is the process of disolving underground mineral deposits by
pumping fresh, or at least unsaturated, water down a well, allowing it to
leach the minerals into solution, and then pumping the solution back to the
surface. Solution mining is also used to create underground caverns --
typically in salt deposits -- that are then used for energy storage; typically
this is through storage of natural gas or crude oil, but hydrogen storage is
beginning to gain momentum.

SANSMIC was developed in the early 1980s to fill a specific modeling need for
the U.S. Strategic Petroleum Reserve (:term:`SPR`) that could not be met
using the tools available at the time. That specific need was to
model the injection of both raw water, for cavern leaching, and crude
oil simultaneously in a process referred to as "leach/fill". SANSMIC has
primarily been used at the :term:`SPR` for modeling impacts due to emergency
exchanges, sales, remediations, and the original construction.



Citing sansmic
--------------



Indices and tables
------------------

* :ref:`references`
* :ref:`nomenclature`
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
