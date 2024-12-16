Outputs
=======




Files created
-------------------






File formats
------------


HDF5-formatted files
~~~~~~~~~~~~~~~~~~~~

.. rubric:: Notation and naming

``\-- item``
    an item is either a group (folder) or a data object (array)
``@-- attribute``
    an attribute is metadata associated with an item
``\-- data: <type>[dim1][dim2]``
    data is a an N-D array of type ``<type>`` (usually string, float, or integer)
    with ``dim1`` rows and, optionally, ``dim2`` columns.
``(identifier)``
    an identifier surrounded by ``( )`` indicates that a group or data is optional
``<identifier>``
    an identifier surrounded by ``< >`` indicates that the group is indexed, e.g., 1, 2, ...


.. rubric:: General organization


.. rubric:: Spatially referenced data


.. rubric:: Time dependent data



.. rubric:: File layout


.. plantuml::

    root level
    \-- sansmic
    \-- cavern
    \-- scenario
        \-- <id>
            \-- config
            \-- simulation

.. verbatim::

    sansmic
     @-- version: Integer[2]
     \-- author
     |   @-- name: String[]
     |   @-- (email: String[])
     \-- creator
         @-- name: String[]
         @-- version: String[]

    cavern
     @-- name: String[]
     @-- wells: String[variable]
     \-- (depthzero)
     |    @-- crs: String[]
     |    \-- elev: Float[]
     \-- (centerline)
     |    @-- crs: String[]
     |    @-- axis: String[2]
     |    \-- coords: Float[2]
     \-- geometry: Float[variable][by-format]
     |    @-- format: String[]
     |    @-- num-cells: Integer[]
     |    @-- depth-range: Float[2]

    config
     @-- title: String[]
     @-- comments: String[]
     @-- units: String[]
     \-- casings
     |    \-- <id:Integer>: Float[2]
     \-- solver
     |    @-- timestep: Float[]
     |    @-- save-frequency: Integer[]
     |    @-- absolue-error: Float[]
     |    @-- relative-error: Float[]
     |    ...
     \-- stages
         \-- <id:Integer>
             @-- title: String[]

    simulation
     \-- node_depth:       Float[nn]
     \-- node_height:      Float[nn]
     \-- node_radius0:     Float[nn]
     \-- time_hours:       Float[ns]
     \-- time_days:        Float[ns]
     \-- state
     |    \-- step_id:         Integer[ns]
     |    \-- stage_id:        Integer[ns]
     |    \-- run_mode:        Integer[ns]
     |    \-- inj_phase:       Integer[ns]
     |    \-- inj_cell:        Integer[ns]
     |    \-- prod_cell:       Integer[ns]
     |    \-- plume_cell:      Integer[ns]
     |    \-- obi_cell:        Integer[ns]
     |    \-- convergence:     Float[ns]
     |    \-- stepsize:        Float[ns]
     \-- instantaneous
     |    \-- inj_depth:       Float[ns]
     |    \-- prod_depth:      Float[ns]
     |    \-- plume_depth:     Float[ns]
     |    \-- obi_depth:       Float[ns]
     |    \-- flow_out:        Float[ns]
     |    \-- flow_in:         Float[ns]
     |    \-- jet_length:      Float[ns]
     |    \-- jet_velocity:    Float[ns]
     |    \-- jet_radius:      Float[ns]
     |    \-- inj_sg:          Float[ns]
     |    \-- prod_sg:         Float[ns]
     |    \-- cavern_sg:       Float[ns]
     |    \-- cavern_volume:   Float[ns]
     \-- cumulative
     |    \-- insol_top_depth: Float[ns]
     |    \-- insol_height:    Float[ns]
     |    \-- inj_volume:      Float[ns]
     |    \-- prod_volume:     Float[ns]
     |    \-- fill_volume:     Float[ns]
     |    \-- insol_volume:    Float[ns]
     |    \-- vented_volume:   Float[ns]
     \-- nodes
          \-- radius:          Float[ns][nn]
          \-- dr0:             Float[ns][nn]
          \-- dr_dt:           Float[ns][nn]
          \-- theta:           Float[ns][nn]
          \-- xincl:           Float[ns][nn]
          \-- volume:          Float[ns][nn]
          \-- density:         Float[ns][nn]
          \-- dC_dt:           Float[ns][nn]
          \-- dC_dz:           Float[ns][nn]
          \-- dis_flag:        Integer[ns][nn]
          \-- dis_factor:      Float[ns][nn]
          \-- diff_coeff:      Float[ns][nn]
          \-- plume_radius:    Float[ns][nn]
          \-- plume_velocity:  Float[ns][nn]
          \-- plume_density:   Float[ns][nn]
          \-- inflow:          Float[ns][nn]
