Introduction
============

The U.S. Strategic Petroleum Reserve (:term:`SPR`) is the largest government-owned
crude oil storage facility in the world. It is comprised of 62 underground solution-mined
salt caverns located across four sites in Texas and Louisiana. The SPR is operated by the
U.S. Department of Energy (:term:`DOE`). Sandia National Laboratories (:term:`Sandia` or SNL)
has served as geotechnical advisor to DOE for the SPR project since 1979.
The daily operations and maintenance of the sites is handled for DOE by a
Management and Operations (:term:`M&O`) contractor.

When the SPR was created in the late 1970s, the U.S. Government purchased
over twenty caverns to jumpstart storage.
Of those purchased, only 14 caverns were determined to be suitable for crude oil storage,
11 of which are still in service.
The remaining 49 active caverns were designed and constructed by DOE from the
early 1980s to the early 1990s specifically for use as crude oil storage caverns.

Solution mining is not like traditional mining. A well, or wells, are drilled
into an underground salt formation; fresh water is injected into the well,
which disolves the salt and creates a void space filled with brine.
Once the salt has disolved, more water is injected, forcing the brine out, and
enlarging the void space. This process is continued until either, (a) the
desired mineral has been mined out (and separated from the brine on the surface),
or, (b) the brine filled cavern is of the desired size (as was the case with
the SPR). Once the SPR caverns were created, oil was pumped in, displacing the brine,
and the caverns officially became "storage caverns".

Because of the SPR's size -- each of the SPR caverns could easily hold the Empire State Building
inside it -- there is no way to store enough brine on the surface to use
brine to get the oil back out. This means that withdrawing the oil from the reserve
requires pumping raw water back into the cavern to displace the oil.
Raw water is any fresh, brackish, or saline water that is "undersaturated" when compared to
brine; however, using raw water means that every time oil is removed from the SPR caverns the
caverns get larger. The sansmic software is used to model how operations will
affect the size and shape of the SPR caverns.


The origins of SANSMIC
----------------------
The original authors of SANSMIC are not available for the current authors
to ask questions of, and so much of the history of SANSMIC comes from
reports or notes within other documents or, frankly, rumor and second- or
third-hand knowledge. The following is the best information the current authors
have at this time.

The SANSMIC program was written between 1980 and 1983 by A.J. Russo
at :term:`Sandia` :cite:`russo1983`.
The program was a modification of SALT [1]_, by A. Saberian and A.L. Podio,
:cite:`saberian1974,saberian1976,saberian1977` which was developed from 1974 to 1977
with support from the Solution Mining Research Institute (:term:`SMRI`).
Russo's modifications were related to the addition of a moving oil-brine
interface (:term:`OBI`); while the SALT code modeled the ordinary leaching
that is used for standard solution mining operations, the U.S. Government
decided that the SPR caverns at the Bryan Mound site should be filled
while they were still being mined. This leach-fill mode, with a moving
oil-brine interface during leaching, and a raw-water-based withdrawal mode
were the primary drivers behind Russo developing SANSMIC as a separate program.

The SANSMIC program was written in FORTRAN, and it was modified at least
twice since that time, once in 1991 and again in the early 2010s. When
Weber et al. attempted to validate the code in 2014,
they found that there were multiple, previously undocumented modifications to the
code :cite:`weber2014,weber2015`.
They documented their findings in a design document and they re-validated
the code using some of the original validation tests done by T. Eyermann in the
1980s :cite:`eyermann1984`.
(A technical review of the validation document :cite:`weber2014`
was graciously provided by Thomas Eyermann during the work in 2014.)

While the FORTRAN code could likely have been maintained without modification for
some additional time, it was felt by the current Sandia team that supports the SPR
that there were outstanding needs that could not be reasonably developed
in FORTRAN at this point and that a change was needed. The SANSMIC program
was rewritten in the C++ and Python programming languages, and is being released
as the current software package, `sansmic`.
The lower-case name `sansmic` will be used throughout
this document to indicate the new version of the software, and the all-caps
name `SANSMIC` will refer to the FORTRAN versions of the software.

Also of note, the SANSMIC software used version 8.1 of the Sandia Mathematical
Program Library (Sandia National Laboratories, 1980) to solve the ODEs.
The functions used (ODE, INTRP, STEP1) were written by L.F. Shampine
and M. K. Gordon :cite:`shampine1975`. The elegance and speed of these functions
led to their continued use as the ODE solver in the current version of sansmic.


Changes in sansmic
------------------
The most significant change between SANSMIC and sansmic is obviously the
programming language used. The FORTRAN code was functional code with
numerous goto statements, hard-coded parameters and memory limits,
and required file-based scratch storage.
The new code is object-oriented, uses dynamically assigned parameters
and memory management, and allows for interactive and on-demand
queries of the simulation state.

The second largest change is in the input file format. sansmic
will read the old-style, formatted text files that were used by SANSMIC.
It will also convert them to a more human-friendly format; currently,
the TOML, YAML, and JSON formats are all supported.

Numerous functions that are no longer relevant to the use of sansmic with
the SPR caverns have also been deprecated, such as resetting the geometry
of the cavern mid-simulation.




.. [1] The original SALT code by Saberian and Podio became the :term:`SMRI` software
    SALGAS. The current version of SALGAS (v4.2)
    is available from SMRI. https://www.solutionmining.org/software
