.. _nomenclature:

Nomenclature
============

Notation
--------

The following abbreviations, acronyms, initialisms, and mathematical and unit
symbols are used in this documentation.

.. glossary::
    :sorted:

    DOE
        `U.S. Department of Energy <https://www.energy.gov/>`_

    SPR
        `U.S. Strategic Petroleum Reserve <https://www.spr.doe.gov/>`_

    Sandia
        `Sandia National Laboratories <https://www.sandia.gov/>`_

    SNL
        `Sandia National Laboratories <https://www.sandia.gov/>`_

    SMRI
        `Solution Mining Research Institute <https://www.solutionmining.org>`_

    NIST
        `National Institute of Standards and Technology <https://www.nist.gov/>`_

    sg
        specific gravity

    M&O
        maintenance and operations

    EOT
        end of tubing

    OBI
        oil-brine interface

    TD
        total depth

    ID
        inner diameter, internal diameter

    OD
        outer diameter, outside diameter

    MD
        :term:`measured depth`

    USC
        :term:`U.S. Customary<customary units>` (see also :cite:`NIST-HB-44`)

    SI
        International System of Units (see also :cite:`NIST-SP-330`)

    BPD
        barrels per day

    QA
        quality assurance

    I/O
        input/output

    wt. pct.
        weight percent

    I/F
        interface

    coeff.
        coefficient

    API (software)
        application programming interface

    API (organization)
        American Petroleum Institute

    ODE
        ordinary differential equation

    WL
        withdrawal leach mode

    OL
        ordinary leach mode

    LF
        leach/fill mode

    SF
        product storage fill mode

    ZDP
        :term:`zero-depth point`

    bbl
        (unit symbol) :term:`barrel`, oil barrel

    ft
        (unit symbol) :term:`foot`, international foot

    in
        (unit symbol) :term:`inch`, standard inch

    Mbbl
        (unit symbol) one **thousand** :term:`barrels<barrel>` (= 10³ bbl)

        **NOTE:** this is a :term:`U.S. customary<customary units>` prefixed symbol;
        this is *not* a "mega-barrel" and is not using the SI "mega-" prefix.
        **Do not** use this notation with SI units.

    MMbb
        (unit symbol) one **million** :term:`barrels<barrel>` (= 10⁶ bbl)

        **NOTE:** this is a :term:`U.S. customary<customary units>` prefixed symbol - **do not** use
        this notation with SI units.

    °F
        (unit symbol) :term:`degree Fahrenheit`

    lb
        (unit symbol) :term:`pound`, avoirdupois pound

    lbf
        (unit symbol) :term:`pound-force`

    psi
    lbf/in²
        (unit symbol) :term:`pound-force per square inch`

    GUI
        graphical user interface

    TOML
        (file format) `[Tom's Obvious Minimal Language] <https://toml.io/>`_

    YAML
        (file format) `YAML Ain't Markup Language <https://yaml.org/>`_

    JSON
        (file format) `JavaScript Object Notation <https://www.json.org/>`_

    CSV
        (file format) Comma Separated Values

    HDF5
        (file format) `Heirarchical Data Format v5 <https://docs.hdfgroup.org/hdf5/v1_14/index.html>`_

    LAS
        (file format) `Log ASCII Standard <https://www.usgs.gov/programs/national-geological-and-geophysical-data-preservation-program/las-format>`_

    MarkDown
        (file format) `MarkDown <https://github.com/skills/communicate-using-markdown>`_

    reStructuredText
        (file format) `reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_


.. _glossary:

Glossary
--------

.. glossary::
    :sorted:

    SANSMIC : Software code
        Sandia Solution Mining Code

    SALT : Software code
        :term:`SMRI` solution mining code

    customary units
        Sansmic uses :term:`U.S. Customary units` internally.
        Customary units are non-SI units still used in the U.S. and some other countries. Customary units
        still in use are nearly all based on the English :term:`foot` in some way. However,
        :term:`U.S. Customary units` *are not the same as Imperial units,* and referring to them
        that way is a recipe for confusion :cite:`NIST-HB-44`.

        Note that U.S. Customary units are also *not the same* as now-deprecated :term:`U.S. Survey units`.

        *See:* `NIST Handbook 44 appendices B and C <https://www.nist.gov/pml/owm/nist-handbook-44-current-edition>`_.

    U.S. Customary units
        The International :term:`customary units` plus additional units that are U.S. specific, such as the
        U.S. liquid gallon. The :term:`oil barrel` is not techincally a U.S. Customary unit, although it is
        defined in terms of U.S. Customary units (1 barrel = 42 U.S. gallons).
        Note that U.S. Customary units are also *not the same* as now-deprecated :term:`U.S. Survey units`.

        *See:* `NIST Handbook 44 appendices B and C <https://www.nist.gov/pml/owm/nist-handbook-44-current-edition>`_.

    U.S. Survey units
        **DEPRECATED** The old (pre-2022) U.S. Survey units are based on the foot of the late 1800s. This foot was
        renamed the "U.S. Survey Foot" in 1959, when the :term:`international foot` became the standard
        "foot" in the United States  :cite:`NIST-HB-44`. In 2022, the U.S. Government
        deprecated the U.S. Survey Foot and its derived units  :cite:`NIST-HB-44`. All survey
        units (acre, chain, etc.) are now defined in terms of the international foot :cite:`NIST-HB-44`.

        When referring to historical measurements, and
        using the deprecated units, the unit should be prefixed with the letter "s" (e.g., sft and sac instead
        of ft and ac)  :cite:`NIST-HB-44`. They should no longer be used when creating new documents or new measurements.
        The "survey foot" and "survey mile" are now the preferred terms for the old U.S. Survey
        Foot and U.S. Survey Mile with unit symbols "sft" and "smi", not "usft" or "usmi" :cite:`NIST-HB-44`.

        *See:* `NIST Handbook 44 appendices B and C <https://www.nist.gov/pml/owm/nist-handbook-44-current-edition>`_.

    barrel
    barrels
    oil barrel
        U.S. and petroleum industry unit of volume. The oil barrel is equal to 42 U.S. liquid gallons, which are
        in turn defined as exactly 231 cubic inches. Thus, the oil barrel is equal to 9702 cubic inches or
        539/96 cubic feet (9702/1728 == 539/96) :cite:`NIST-HB-44`.
        The unit symbol is ":term:`bbl`" or "bbl oil".

    foot
    feet
    international foot
        :term:`Customary unit<customary units>` of length, defined as exactly 0.3048 meters :cite:`NIST-HB-44`.
        The unit symbols ":term:`ft`" and "ift" are both acceptable and interchangable :cite:`NIST-HB-44`.

    inch
    inches
        :term:`Customary unit<customary units>` of length, defined as exactly 0.0254 meters (12 in = 1 :term:`ft`) :cite:`NIST-HB-44`.
        Unit symbol is ":term:`in`" :cite:`NIST-HB-44`.

    degree Fahrenheit
        :term:`Customary unit<customary units>` of temperature. A change of one degree Fahrenheit is exactly
        equal to a change of 5/9 degree Celsius. The unit symbol is ":term:`°F`" or "degF" if limited
        to alphanumeric characters :cite:`NIST-HB-44`. Conversion is done using the following
        formula.

        .. math::

            T~/^\circ\mathrm{F} = 32 + \frac{5}{9} T~/^\circ\mathrm{C}

    pound
    pounds
    pound avoirdupois
        :term:`Customary unit<customary units>` of mass. One pound avoirdupois is defined as exactly 0.45359237
        kilograms  :cite:`NIST-HB-44`. The unit symbol is ":term:`lb`" or "lb avdp" (to differentiate it from some other pound) :cite:`NIST-HB-44`.

    pound-force
    pounds-force
        :term:`Customary unit<customary units>` of force, equal to one :term:`pound` accelerated at
        standard gravity (gₙ ≡ 9.80665 m/s²). One pound-force is exactly equal to 4.4482216152605 newtons
        The unit symbol is ":term:`lbf`" :cite:`NIST-HB-44`.

    pound-force per square inch
    pounds-force per square inch
        A :term:`customary unit<customary units>` of pressure. The technical unit symbol is "lbf/in²",
        but the symbol ":term:`psi`" is more commonly used; "PSI" is an abbreviation appropriate for
        text, but should not be used as a unit symbol :cite:`NIST-HB-44`. The PSI does not have an exact decimal representation,
        though it can be written exactly as a fraction. 1 psi is approximately equal to 6894.757 pascal.

    zero-depth point
        The zero-depth point (:term:`ZDP`) is the point on a wellhead, rig, or other permanent datum, that is used to set

        .. math::

            z_\mathrm{zdp} := 0.0

        for :term:`measured depth` (:term:`MD`) values. This should not change over time once a well
        has been completed.

    measured depth
        The measured depth (:term:`MD`) is the positive-valued distance from a specific :term:`zero-depth point` (:term:`ZDP`) downhole within a
        well, wellbore, or cavern. It is not corrected or adjusted to another datum (e.g., sea level), and
        it is not corrected for well deviation or horizontal movements.



U.S. Customary Units
--------------------
It is assumed that any user is familiar with :term:`SI` units; however,
for those that want the most up-to-date information, NIST-SP-330 provides the U.S.
English-laguage version of the official SI definition. :cite:`NIST-SP-330`

The following :term:`U.S. Customary units` (USC) are non-SI units that are
used in SANSMIC along with the equivalent SI unit and conversion
factor. See NIST-HB-44 Appendix B and C for more information on customary units
and the legal and official conversion factors specified by the U.S.
Government. :cite:`NIST-HB-44`


.. only:: html

    .. dropdown:: USC to SI conversion factors

        +----------------+----------------------------------+
        | USC unit       | Definition (in SI units)         |
        +================+==================================+
        | :term:`ft`     | := 0.3048 m                      |
        +----------------+----------------------------------+
        | :term:`in`     | := 2.54 cm                       |
        +----------------+----------------------------------+
        | in²            | ≡ 6.4516 cm²                     |
        +----------------+----------------------------------+
        | ft³            | ≡ 0.028316846592 m³              |
        +----------------+----------------------------------+
        | :term:`Mbbl`/d | ≡ 158.987294938 m³/d             |
        +----------------+----------------------------------+
        | :term:`lb`     | := 0.45359237 kg                 |
        +----------------+----------------------------------+
        | :term:`lbf`    | ≡ 4.4482216152605 N              |
        +----------------+----------------------------------+
        | :term:`psi`    | ≡ 1 lbf/in²                      |
        |                |                                  |
        |                | ≡ (44482216152605 / 6451600) kPa |
        |                |                                  |
        |                | ≈ 6.894757 kPa                   |
        +----------------+----------------------------------+
        | :term:`sg`     | = 1000 kg/m³                     |
        +----------------+----------------------------------+


.. only:: latex

    +----------------+-------------------------------------------------+
    | USC unit       | Definition (in SI units)                        |
    +================+=================================================+
    | :term:`ft`     | :math:`:=\qty{0.3048}{m}`                       |
    +----------------+-------------------------------------------------+
    | :term:`in`     | :math:`:=\qty{2.54}{cm}`                        |
    +----------------+-------------------------------------------------+
    | in²            | :math:`\equiv\qty{6.4516}{cm^2}`                |
    +----------------+-------------------------------------------------+
    | ft³            | :math:`\equiv\qty{0.028316846592}{m^3}`         |
    +----------------+-------------------------------------------------+
    | :term:`Mbbl`/d | :math:`\equiv\qty{158.987294938}{m^3/d}`        |
    +----------------+-------------------------------------------------+
    | :term:`lb`     | :math:`:=\qty{0.45359237}{kg}`                  |
    +----------------+-------------------------------------------------+
    | :term:`lbf`    | :math:`\equiv\qty{4.4482216152605}{N}`          |
    +----------------+-------------------------------------------------+
    | :term:`psi`    | :math:`\equiv\qty{1}{lbf/in^2}`                 |
    |                |                                                 |
    |                | :math:`\equiv44482216152605/6451600\unit{kPa}`  |
    |                |                                                 |
    |                | :math:`\approx\qty{6.894757}{kPa}`              |
    +----------------+-------------------------------------------------+
    | :term:`sg`     | :math:`\sim\qty{1000}{kg/m^3}`                  |
    +----------------+-------------------------------------------------+
