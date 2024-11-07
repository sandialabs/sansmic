Scenarios
=========

sansmic runs a simulation based on a user-defined scenario.
In SANSMIC, this was defined in a .DAT file. Modern scenarios are
defined programmatically, using :class:`~sansmic.model.Scenario`,
:class:`~sansmic.model.AdvancedOptions`,
and :class:`~sansmic.model.StageDefinition` objects, or through
a dictionary representation of these classes.
The pereferred file format for a scenario file is :term:`TOML`,
although :term:`JSON`, and :term:`YAML` will be parsed if the proper
extension is used. All configuration examples will be provided in
TOML format.

The scenario file starts with options that are not contained in
a TOML table. The following is the basic layout of a scenario file

.. code-block:: toml

    title = 'A title for my scenario'   # optional but highly recommended
    num-cells = 200                     # required
    # ... additional per-scenario options will follow

    [advanced]
    # ... advanced and unusual settings, such as ...
    relative-error = 1.0e-4             # the ODE solver tolerance

    [defaults]
    # ... values which can be overridden in a specific stage, but which
    #     you may want a default value so you don't *have* to specify
    #     more than once
    solver-timestep = 0.1         # hours
    inner-tbg-inside-diam = 9.85  # inches

    [[stages]]
    # ... per-stage options go here

    # [[stages]]
    # ... define an additional stage with another [[stages]] entry


Optional but highly recommended scenario keys are

* :confval:`title`
* :confval:`comments`


Required scenario keys are

* :confval:`num-cells`
* :confval:`geometry-data`
* :confval:`geometry-format`
* :confval:`cavern-height`
* :confval:`floor-depth`
* :confval:`ullage-standoff`
* :confval:`insolubles-ratio`
* :confval:`units` (in progress)


Valid :confval:`[defaults]` keys are

* :confval:`solver-timestep`
* :confval:`save-frequency`
* :confval:`inner-tbg-inside-diam`
* :confval:`inner-tbg-outside-diam`
* :confval:`outer-csg-inside-diam`
* :confval:`outer-csg-outside-diam`


The `[advanced]`_ section keys are

* :confval:`absolute-error`
* :confval:`relative-error`
* :confval:`coallescing-wells`
* :confval:`well-separation`
* :confval:`jet-model-version`
* :confval:`dissolution-factor`
* :confval:`max-brine-sg`
* :confval:`solid-density`
* :confval:`entrainment-coeff`
* :confval:`molecular-diffusion`
* :confval:`eddy-coefficient`
* :confval:`diffusion-beta`

..
    * :confval:`plume-model-version`
    * :confval:`temperature-model-version`


The `[[stages]]`_ table require that **all** of the following
be defined for **each** stage

* :confval:`title`
* :confval:`simulation-mode`
* :confval:`injection-duration`
* :confval:`rest-duration`
* :confval:`brine-injection-sg`
* :confval:`brine-injection-depth` xor :confval:`brine-injection-height`
* :confval:`brine-production-depth` xor :confval:`brine-production-height`
* :confval:`brine-injection-rate`


The keys which are required for at least the **first** stage are

* :confval:`set-cavern-sg`
* :confval:`brine-interface-depth` xor :confval:`brine-interface-height`


The keys which can be set using the ``[defaults]`` section
can also be set in the ``[stages]``, but they must be defined
in every stage if they are not given defaults.

Finally, the following options are optional or are only required
in a certain :confval:`simulation-mode`

* :confval:`set-initial-conditions`
* :confval:`product-injection-rate`
* :confval:`stop-condition`
* :confval:`stop-value`

..
    hack - don't use these
    * :confval:`product-injection-depth`
    * :confval:`product-production-depth`


Scenario options
----------------

.. confval:: title

    A descriptive title for the simulation or for a stage.
    The simulation title is optional, but every stage *must* have
    a title defined.


.. confval:: comments

    A longer block for detailed comments, citations, etc., can be
    defined at the scenario level. Use triple-quotes to allow use
    of multiple lines, like this

    .. code-block:: toml

        comments = """This is an example
        of a comment that spans multiple lines."""

    If an old-style DAT file is used, this section will initially
    be filled in with notes about the conversion process to TOML.


.. confval:: num-cells

    The number of cells to use in the vertical domain. Each cell
    is between two nodes. The spacing of the nodes is defined by
    the :confval:`cavern-height` divided by :confval:`num-cells`.


.. confval:: geometry-format

    There are ways to define the geometry, and it is defined by this
    field. The valid options are given below. For all options except
    *radius-list*, the python module will linearly interpolate between
    depths to get the radii for each node, which is not necessarily a
    volume-conserving operation.

    "radius-list"
        The radii for each *node* are listed in order from the *bottom*
        of the cavern to the top. There are :confval:`num-cells` + 1 nodes
        that are equally spaced based on the :confval:`cavern-height` divided
        by the number of cells. The radii are provided in units of :term:`foot`.
        The radii can be listed in the scenario configuration file or in
        a separate file - see :confval:`geometry-data`.

    "radius-table"
        The radii and depths are both specified, and the radii for each
        evenly-spaced node are interpolated between depths if needed. This
        is not a volume-conserving operation [currently].

    "volume-table"
        The volume of the cavern is specified for a list of depths. The
        volumes are then converted to a cumulative volume "strapping" curve,
        and the cumulative volume is interpolated for each evenly-spaced
        node. After this interpolation is performed, the volumes are de-accumulated
        and the radii calculated. This is a volume-conserving operation.

    "volume-list"
        This is actually a list of volumes *and* depths, and is identical
        to a volume table unless the data is in a separate file, in which
        case the number of entries must be specified on the first line,
        followed by all [double check] and then all [double check].
        If depths and radii are specified in the configuration file, then
        this is identical to the *volume-table* format.

    ..
        The following new formats are planned, but not yet implemented:

        "infer"
            The format will be inferred from the fields provided in the dictionary
            which must be one of the standard fields defined in :confval:`geometry-data`.
            If a file is specified, it must be in TOML, JSON, or YAML format.

        <dictionary>
            Will allow importing a file using :mod:`pandas` or :mod:`lasio`. If a piece
            of data is provided in a field not named one of those defined in
            :confval:`geometry-data`, then providing that name as a field will allow
            the user to map that field name to the appropriate column in the file.
            The "units" key will allow for files that are not specified with the same
            units as the rest of the configuration file.


.. confval:: geometry-data
    :type: string or dict

    If the :confval:`geometry-data` field is a string, it is assumed to
    be a path to a file that contains *only* data in the format defined
    by the :confval:`geometry-format`. I.e., if the geometry format is
    *radius-list*, then it must be a list of numbers, one per line, and
    nothing else in the file.

    If the :confval:`geometry-data` is a dictionary, it can have the
    following fields; the fields needed and how they are interpreted and/or
    interpolated are determined by the :confval:`geometry-format`.

    .. confval:: geometry-data.radii

        A list of radii. If :confval:`geometry-format` = *radius-list*,
        then this should be the only field provided, and the radii must
        be appropriately specified from bottom to top.

        If :confval:`geometry-format` = *radius-table*, then this field
        should have the radii corresponding to the depths in
        :confval:`geometry-data.depths`, and the order is irrelevant provided
        the two fields are in the same order.


    .. confval:: geometry-data.volumes


    .. confval:: geometry-data.depths



.. confval:: cavern-height

    The total height of the model domain in the vertical axis.


.. confval:: floor-depth

    The measured depth of the floor of the cavern/well. This
    is a positive value below a surface datum.


.. confval:: ullage-standoff

    The ullage, or remaining usable volume, is calculated by
    removing the volume between the floor and this distance
    from the brine volume.


.. confval:: insolubles-ratio

    The volumetric fraction of the rock salt that is composed
    of insoluble material that will accumulate on the cavern floor.


.. confval:: units

    .. warning::
        Units is an open issue and is the first priority
        for updates. Currently, only `ft-in-bbl` is enabled.

    Specify the units that all options will be provided in. The
    following values are defined

    "ft-in-bbl"
        Depths, heights, and cavern radii are specified in feet;
        casing and tubing radii are specified in inches; volumes
        and flow rates are defined in barrels and barrels/day,
        respectively.

    "ft-in-ft3"
        Depths, heights, and cavern radii are specified in feet;
        casing and tubing radii are specified in inches; volumes
        and flow rates are defined in cubic feet and barrels/day,
        respectively.

    "m-cm-m3"
        Depths, heights, and cavern radii are specified in meters;
        casing and tubing radii are specified in centimeters; volumes
        and flow rates are defined in cubic meters and cubic meters/day,
        respectively.


.. confval:: [defaults]

    The following values can be set in the :confval:`[defaults]` section,
    in which case their value will be overridden by a value set in
    a stage definition. However, if they are not specified here,
    they must be specified in every stage definition.

    * :confval:`solver-timestep`
    * :confval:`save-frequency`
    * :confval:`inner-tbg-inside-diam`
    * :confval:`inner-tbg-outside-diam`
    * :confval:`outer-csg-inside-diam`
    * :confval:`outer-csg-outside-diam`


``[advanced]``
----------------
The advanced options table is not required, and all options have
default values that should not be changed unless necessary. Some
options provided here will be deprecated in the future, but because
v1.0.0 of sansmic is designed to replicate the FORTRAN code, they
are still available.

.. confval:: absolute-error

    A new option, this allows manual specification of the ODE solver
    absolute error tolerance. By default this value is 1.0e-2, and this
    appears to be generally sufficient.


.. confval:: relative-error

    A new option, this allows manual specification of the ODE solver
    relative error tolerance. By default this value is 1.0e-4, and this
    appears to be generally sufficient.


.. confval:: coallescing-wells

    This option will be deprecated. Define the number of wells that
    are being used in ordinary leaching, by default 1. This value
    should always be 1 for withdrawal or leach-fill modes and for
    ordinary leaching where the cavern has already coallesced.


.. confval:: well-separation

    This option will be deprecated. Define the distance between
    coallescing wells in ordinary leaching mode. If the number of
    wells is 1, then this value has no impact.


.. confval:: jet-model-version

    The jet model is 1 by default. To turn it off and have the
    plume start at the :term:`EOT`, use a value of 0. If and when
    the jet model is improved, it will be assigned a different version
    number, or possibly a name, that can be used to select it.

..
    future development

    `plume-model-version`
    ~~~~~~~~~~~~~~~~~~~~~~~

    `temperature-model-version`
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. confval:: dissolution-factor

    This value should nearly always be left alone with its default
    value of 1.0. This adjusts the rate of dissolution compared to
    standard dissolution of halite.


.. confval:: max-brine-sg

    Override the maximum density (in specific gravity terms) of the
    brine in the cavern. This is only needed if the salt rock is not
    halite. The default is 1.2019 sg.


.. confval:: solid-density

    Override the density of the surrounding rock salt (in g/cm3).
    This density sould be the bulk density of the salt, and should
    not include insoluble material (which will be accounted for
    volumetrically). Default value is 2.16 g/cm3.


.. confval:: entrainment-coeff

    The entrainment coefficient 'alpha' in the dissolution equations.


.. confval:: molecular-diffusion

    The molecular diffusion coefficient, D_mol.


.. confval:: eddy-coefficient

    The eddy coefficient, D_0, for diffusion calculations.


.. confval:: diffusion-beta

    The 'beta' coefficient in diffusion calculations.



``[[stages]]``
--------------
A new stage is defined by the ``[[stages]]`` table header.
At least one stage is required, and there are no limits on the
number of stages. A stage **requires** a title key.

.. confval:: simulation-mode

    Which fluid flow regime is being modeled. Valid options are

    "ordinary"
        Ordinary leaching with "brine" injection and production.

    "withdrawal"
        Withdrawal leach with "brine" injection and "product"
        production.

    "leach-fill"
        Leaching while filling, with "brine" injection and production and "prodcut" injection.


.. confval:: solver-timestep

    Can also be set in the ``[defaults]`` section; if so, it will
    be overridden by a setting in the ``[[stages]]`` definition.
    If not provided a default value, it must be specified in each stage.

    The solver timestep is specified in decimal hours.


.. confval:: save-frequency

    Can also be set in the ``[defaults]`` section; if so, it will
    be overridden by a setting in the ``[[stages]]`` definition.
    If not provided a default value, it must be specified in each stage.

    The save frequency is specified in integer number of *timesteps*.
    It can also be set with one of the following three strings

    "hourly"
        Save results every hour, where an "hour" is based on the
        rounded, integer value of the inverse of the solver timestep.
        If you want this to be exact, make sure that the timestep
        is exactly equal to an hour divided by an integer.
        Saving results every hour is likely overkill unless you are
        debugging.

    "daily"
        Save results every day, where "day" is calculated in the same
        way as an hour is for 'hourly'. *This is the recommended option
        for most simulations.*

    "bystage"
        Only save results at the end of each injection and each rest
        period. This is the most efficient method in terms of storage
        and processing time, but provides little data for debugging or
        careful analysis.


.. confval:: injection-duration

    The duration of active injection activities, defined in fractional
    hours. This value cannot be 0. This value is required for all
    stages and cannot be given a default value.


.. confval:: rest-duration

    The duration of the post-injection inactive or static period,
    defined in fractional hours. This value can be 0.
    This value is required for all stages and cannot be given a
    default value.


.. confval:: inner-tbg-inside-diam

    Can also be set in the ``[defaults]`` section; if so, it will
    be overridden by a setting in the ``[[stages]]`` definition.
    If not provided a default value, it must be specified in every stage.

    The :term:`ID` of the innermost tubing for concentric completions,
    or for the brine string for single-string completions.
    Specified in decimal :term:`inch`.


.. confval:: inner-tbg-outside-diam

    Can also be set in the ``[defaults]`` section; if so, it will
    be overridden by a setting in the ``[[stages]]`` definition.
    If not provided a default value, it must be specified in every stage.

    The :term:`OD` of the innermost tubing for concentric completions,
    or for the brine string for single-string completions.
    Specified in decimal :term:`inch`.


.. confval:: outer-csg-inside-diam

    Can also be set in the ``[defaults]`` section; if so, it will
    be overridden by a setting in the ``[[stages]]`` definition.
    If not provided a default value, it must be specified in every stage.

    The :term:`ID` of the outermost tubing for concentric completions,
    or for the production casing for single-string/slick well
    completions. Specified in decimal :term:`inch`.


.. confval:: outer-csg-outside-diam

    Can also be set in the ``[defaults]`` section; if so, it will
    be overridden by a setting in the ``[[stages]]`` definition.
    If not provided a default value, it must be specified in every stage.

    The :term:`OD` of the outermost tubing for concentric completions,
    or for the production casing for single-string/slick well
    completions. Specified in decimal :term:`inch`.


.. confval:: brine-injection-sg

    The specific gravity of the water being injected, regardless of if
    is raw water or brine. This value cannot be less than 1.0 and is
    also limited by the :confval:`max-brine-sg`.


.. confval:: brine-injection-rate

    The rate of brine injection in *barrels per day*.
    Variable injection rates are still in the process of being
    implemented. The brine/raw water is injected at a constant
    rate over the entire injection duration.


.. confval:: brine-injection-depth

    The *depth* (preferred) or *height* within the cavern
    where the EOT for the injection string is located. If above the
    production depth, it will be assigned to the outer casing,
    otherwise it will apply to the inner casing.

    The :confval:`brine-injection-depth` is specified in :term:`foot`
    below the :term:`ZDP` at the surface.

    The :confval:`brine-injection-height` is specified in :term:`foot`
    above the original :confval:`floor-depth`

    Do not use both -depth and -height options in the same stage.


.. confval:: brine-injection-height

    See :confval:`brine-injection-depth`


.. confval:: brine-production-depth

    The *depth* (preferred) or *height* within the cavern
    where the EOT for the injection string is located. If above the
    injection depth, it will be assigned to the outer casing,
    otherwise it will apply to the inner casing.

    The :confval:`brine-production-depth` is specified in :term:`foot`
    below the :term:`ZDP` at the surface.

    The :confval:`brine-production-height` is specified in :term:`foot`
    above the original :confval:`floor-depth`

    Do not use both -depth and -height options in the same stage.


.. confval:: brine-interface-height

    See :confval:`brine-interface-depth`


.. confval:: brine-interface-depth

    The initial brine interface *depth* (preferred) or *height*
    within the cavern. This is only required for the first stage, and
    should not be set in subsequent stages unless you *want* to
    manually reset the interface position. If so, you should also
    set :confval:`set-initial-conditions` to true. If left blank or set
    to 0, it will behave as if unset and use the value from the
    previous stage.

    The :confval:`brine-interface-depth` is specified in :term:`foot`
    below the :term:`ZDP` at the surface.

    The :confval:`brine-interface-height` is specified in :term:`foot`
    above the original :confval:`floor-depth`

    Do not use both -depth and -height options in the same stage.


.. confval:: brine-production-height

    See :confval:`brine-production-depth`


.. confval:: set-cavern-sg

    The initial specific gravity for the cavern brine. This is
    optional. If set for the first stage, then the cavern will
    be initialized to the specified brine sg value. If omitted
    for the first stage, it will be set to the maximum brine
    specific gravity.

    *Changed: v1.0.6*

    The SANSMIC user guide from 2015 had an error, and stated that
    subsequent stages did not use this value. In fact, this value
    was still being set unless it was set to 1.0 or less.

    sansmic will behave *as was described in the user manual*. To
    replicate the old behavior, set :confval:`set-initial-conditions`
    to true to force initialization of the cavern brine to the value
    you specify in this key. (You can also set this to 0.0 in an
    old DAT file to force this behavior in the FORTRAN version)


.. confval:: set-initial-conditions

    Force cavern specific gravity and/or the initial brine interface
    level to be set for subsequent (not the first) stages. By default
    false.


.. confval:: product-injection-rate

    This option is defined in the same way as :confval:`brine-injection-rate`,
    but applies only to the `leach-fill` simulation mode. It cannot
    be used to move the brine interface, currently, as it will cause
    significant (and non-physical) leaching at the brine injection
    location even when the brine injection rate is set to 0.0. This
    is in the list of issues to either fix or create a new simulation
    mode to handle.


.. confval:: stop-condition

    It is possible to add a stop condition to an injection period
    based on either the brine interface location or the total cavern
    size specified in :confval:`stop-value`.
    By default, the only stop condition is reaching the end of the
    injection period.

    Valid values are "duration" (or blank, the default), "depth" or
    "volume".


.. confval:: stop-value

    The value for the stop condition, which is either a *depth* below
    surface, or a volume in *barrels*.

    *Changed: v1.0.7*

    If a depth is specified as the stop value and the :confval:`simulation-mode`
    mode is "withdrawal", then the injection period will stop when
    the duration is exceeded or the interface depth becomes smaller
    (higher) than the value given. If the  :confval:`simulation-mode` is "leach-fill",
    then the simulation will stop when the interface depth becomes
    larger (deeper) than the value given. Note that this is *not* how
    the old FORTRAN version of SANSMIC behaved.
