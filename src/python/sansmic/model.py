# coding: utf-8
#
# Copyright (c) 2024 National Technology and Engineering Solutions of
# Sandia, LLC. Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.
#
# SPDX-License-Identifier: BSD-3-Clause.

"""The core python classes for SANSMIC."""

import warnings
from dataclasses import InitVar, asdict, dataclass, field
from enum import IntEnum
from fractions import Fraction
import logging
from math import nan
from types import NoneType
from typing import Any, Dict, List, Union

import numpy as np
import pandas as pd

logger = logging.getLogger("sansmic")

try:
    from . import libsansmic as _ext
except ImportError:
    logger.critical('The C++ library is not installed. Conversions will work, but the main program will not run.')

class SansmicConfigError(TypeError):
    """Specific error for missing configuration data."""

    pass


class Units(IntEnum):
    """The units that are used to define the scenario.

    Depths, heights, and cavern radii are in the first (larger) length unit.
    Tubing radii are in the second (smaller) length unit.
    Volumes and volumetric flow rates are in the third unit.

    Durations are in hours. Constant injection rates are in <volume-unit> per day.
    File-based injection rates are in <volume-unit> per hour.
    """

    FT_IN_BBL = 1
    """foot/inch/barrel"""
    FT_IN_FT3 = 2
    """foot/inch/cubic foot"""
    M_CM_M3 = 3
    """meter/centimeter/cubic meter"""

    @property
    def inch(self):
        """1 in ≝ 0.0254 m"""
        return Fraction(254, 10000)

    @property
    def foot(self):
        """1 ft ≝ 0.3048 m"""
        return Fraction(3048, 10000)

    @property
    def cubic_foot(self):
        """1 ft³ == 0.028316846592 m³"""
        return Fraction(3048**3, 10000**3)

    @property
    def barrel(self):
        """1 bbl == 0.158987294928 m³"""
        return Fraction(42 * 231 * 254, 10000**3)

    @property
    def centimeter(self):
        """1 cm ≝ 0.01 m"""
        return Fraction(1, 100)


class StopCondition(IntEnum):
    """Which stop condition should be used to end the simulation."""

    DEPTH = -1
    """Stop when the interface reaches a specified height above the cavern floor"""
    DURATION = 0
    """Stop only based on duration"""
    VOLUME = 1
    """Stop when the total cavern volume reaches a specified value."""


class GeometryFormat(IntEnum):
    """The format for the initial cavern geometry."""

    RADIUS_LIST = 0
    """The radius is provided from the bottom of the cavern to the top;
    :attr:`~Scenario.num_cells` + 1 values, equally spaced"""
    VOLUME_LIST = 1
    VOLUME_TABLE = -1
    RADIUS_TABLE = 2
    LAYER_CAKE = 5
    """The geometry is provided in a 'layer-cake' style LAS file"""


class SimulationMode(IntEnum):
    """The simulation mode determines which options are active for injection."""

    ORDINARY = 0
    """Ordinary leaching, with raw water or undersaturated brine injected through
    the inner tubing or outer casing and brine produced from the other."""
    WITHDRAWAL = 1
    """Withdrawal leach, with brine injected through the suspended tubing (hanging 
    string) and product produced from top of cavern."""
    LEACH_FILL = 2
    """Simultaneous leaching (water/brine injection and brine production) and 
    product injection from the top."""
    STORAGE_FILL = -1
    """Storage fill, with product injected from the top of the cavern and brine
    produced from the suspended tubing/hanging string."""


class RateScheduleType(IntEnum):
    """The way that the rate of injection is specified."""

    CONSTANT_RATE = 1
    """Injection occurs at a constant rate of volume :class:`Units` per day."""
    DATAFILE = 2
    """Injection rate is specified in :class:`Units` per hour in a file."""


@dataclass
class StageDefinition:
    """Defines a SANSMIC simulation stage which is an injection period followed by a rest period."""

    title: str = None
    """Title for the stage."""
    simulation_mode: SimulationMode = None
    """The simulation mode used in this stage."""
    solver_timestep: float = None
    """The solver timestep in hours."""
    save_frequency: int = None
    """The save frequency in number of timesteps"""
    injection_duration: float = None
    """The duration of the injection phase of the stage."""
    rest_duration: float = None
    """The duration of the post-injection rest phase of the stage."""

    inner_tbg_inside_diam: float = None
    """The inner tubing inner diameter (ID)."""
    inner_tbg_outside_diam: float = None
    """The inner tubing outer diameter (OD)."""
    outer_csg_inside_diam: float = None
    """The outer tubing/casing inner diameter (ID)."""
    outer_csg_outside_diam: float = None
    """The outer tubing/casing outer diameter (OD)."""

    brine_injection_sg: float = None
    """The water/brine injection specific gravity."""
    brine_injection_depth: float = None
    """The depth below surface/:term:`ZDP` where the end of the injection string/casing/tubing is positioned."""
    brine_production_depth: float = None
    """The depth below surface/:term:`ZDP` where the end of the production string/casing/tubing is positioned."""
    brine_injection_rate: Union[float, str] = 0
    """The meaning of this field is based on the type of data that is provided (a value of 0 means off).

    float : constant rate
        Injection occurs at a constant rate for the duration of the
        :attr:`injection_duration`.
    
    str : file path
        Injection rates should be read from a file. The file can be any type that can
        be read into a pd DataFrame. An 'index' column is optional, but if not provided, then all data is assumed to be pre-sorted. In addition to the 
        'index' column, files may contain exactly one set of the following columns:
        
        ['rate', 'duration']
            This style gives a list of rates and how long they should continue for.
            Durations are provided in decimal hours. If the sum of durations is greater than :attr:`injection_duration`, an error will be raised. If the 
            sum of durations is less than :attr:`injection_duration`, the remaining
            time will be assumed to be no-flow.
        
        ['time', 'rate']
            This style gives a list of new rate values that start at the given time 
            after the start of the stage. The 'time' values must be increasing and
            never duplicate. Times are in decimal hours since the start of the stage.
            The last rate will be used up until the end of the 
            :attr:`injection_duration`. If a time is given that is after the 
            :attr:`injection_duration`, then an error will be raised.

        ['hourly']
            This style specifies a new rate for each hour of the stage. There must be
            a value specified for each hour of the :attr`injection_duration`, or an 
            error will be raised.

        The optional column ['sg'] can be provided to change the injection water/brine
        specific gravity; however, if provided, it must contain a value for every row.

        
        If 'duration' and 'time' occur in the same file, an error will be raised.
        If 'time', 'rate', or 'duration' occur in the same file as 'hourly', then
        an error will be raised. 
        
        Columns with names other than those described above will be ignored.

        .. attention::

            In all cases, a no-flow, leaching-only period will occur for
            :attr:`rest_duration` hours **after** all injection data
            is processed -- this means that if an injection file finishes with 0s
            that there will be an additional :attr:`rest_duration` hours added to the 
            end of the stage.

        .. warning::

            The same file cnnot be used for both brine injection and product
            injection data.
    """

    set_initial_conditions: bool = None
    """Unlink initial cavern brine gravity and interface level from previous stage. 
    Automatically set to True for the first stage added to a model."""
    set_cavern_sg: float = None
    """Set the initial specific gravity for all brine-filled cells of the cavern to 
    this value. If set_initial_conditions is False and this is not None (or 0) an error
    will be raised upon scenario validation."""
    set_interface_level: float = None
    """Set the initial oil-brine interface or blanket level. By default None, which
    will link to the previous stage (as will a value of 0)."""
    # set_cavern_temp: float = None
    # """Set the initial temperature for all brine-filled cells of the cavern; by
    # default None, which will link to the previous stage."""

    product_injection_rate: Union[float, str] = 0.0
    """Either a constant rate of product injection or a file with an injection schedule, by default 0."""
    product_water_content: float = None
    """The volume-percent water within the product; a value of None turns this option off."""
    product_injection_depth: float = None
    """The depth below surface/:term:`ZDP` where the product injection string -- or casing or chimney bottom -- is positioned; a value 0 or None sets this to the roof of the cavern."""
    product_production_depth: float = None
    """The depth below surface/:term:`ZDP` where the product production string -- or casing or chimney bottom -- is positioned; a value 0 or None sets this to the roof of the cavern."""

    stop_condition: StopCondition = StopCondition.DURATION
    """Set an additional stop condition for the stage, see :class:`StopCondition`."""
    stop_value: float = 0
    """If the :attr:`stop_condition` is not :attr:`~StopCondition.DURATION`, then the depth or volume value."""

    defaults: InitVar[Dict[str, Any]] = None

    valid_default_keys = [
        "solver_timestep",
        "save_frequency",
        "inner_tbg_inside_diam",
        "inner_tbg_outside_diam",
        "outer_csg_inside_diam",
        "outer_csg_outside_diam",
    ]

    def __post_init__(self, defaults=None):
        if defaults is None:
            return
        if not isinstance(defaults, dict):
            raise TypeError("defaults must be a dictionary")
        for k, v in defaults.items():
            if k not in self.valid_default_keys:
                logger.warning(
                    "Ignoring non-defaultable or unknown setting {} = {}".format(k, repr(v))
                )
            elif getattr(self, k) is None:
                setattr(self, k, v)

    def __setattr__(self, name, value):
        if isinstance(value, str) and value.strip() == "":
            value = None
        if name == "simulation_mode" and not isinstance(value, SimulationMode):
            if isinstance(value, int):
                value = SimulationMode(value)
            elif isinstance(value, str):
                value = SimulationMode[value.upper().replace(" ", "_").replace("-", "_")]
            elif isinstance(value, _ext.CRunMode):
                value = SimulationMode(int(value))
            else:
                TypeError("simulation_mode cannot be of type {}".format(type(value)))
        elif name == "stop_condition" and not isinstance(value, StopCondition):
            if isinstance(value, int):
                value = StopCondition(value)
            elif isinstance(value, str):
                value = StopCondition[value.upper().replace(" ", "_").replace("-", "_")]
            elif isinstance(value, _ext.CGeometryFormat):
                value = StopCondition(int(value))
            else:
                SansmicConfigError("stop_condition cannot be of type {}".format(type(value)))
        elif name == "set_cavern_sg" and value is not None and value < 1.0:
            value = None
        elif name == "set_initial_conditions" and value is False:
            self.set_cavern_sg = None
            self.set_interface_level = None
        super().__setattr__(name, value)

    def validate(self):
        """Validate that all required options have been selected.

        Raises
        ------
        SansmicConfigError
            if there are missing or misconfigured options
        ValueError
            if timing options are invalid
        """
        if self.solver_timestep <= 0:
            raise ValueError("Timestep must be greater than 0 hours")
        if self.injection_duration <= 0:
            raise ValueError("Injection duration must be greater than 0 hours")
        if self.rest_duration < 0:
            raise ValueError("Rest duration must be positive")

        # Validate appropriate initial conditions settings
        if self.set_initial_conditions and (
            self.set_cavern_sg is None or self.set_interface_level is None
        ):
            raise SansmicConfigError(
                "An initial stage must have both set_cavern_sg and set_interface_level."
            )
        elif not self.set_initial_conditions and self.set_cavern_sg:
            # raise SansmicConfigError(
            warnings.warn(
                "Setting the starting cavern sg ought to use 'set_initial_conditions' to be set to True"
            )
        elif not self.set_initial_conditions and self.set_interface_level:
            warnings.warn(
                "Make sure you meant to reset the interface level; use 0.0 or None to continue from the last stage."
            )

        # Validate appropriate brine injection and production settings
        if (
            self.simulation_mode
            in [SimulationMode.ORDINARY, SimulationMode.LEACH_FILL, SimulationMode.STORAGE_FILL]
            and self.brine_production_depth is None
        ):
            raise SansmicConfigError(
                "Missing required 'brine_production_depth' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode
            in [SimulationMode.ORDINARY, SimulationMode.LEACH_FILL, SimulationMode.WITHDRAWAL]
            and self.brine_injection_depth is None
        ):
            raise SansmicConfigError(
                "Missing required 'brine_injection_depth' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode
            in [SimulationMode.ORDINARY, SimulationMode.LEACH_FILL, SimulationMode.WITHDRAWAL]
            and self.brine_injection_rate is None
        ):
            raise SansmicConfigError(
                "Missing required 'brine_injection_rate' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode
            in [SimulationMode.ORDINARY, SimulationMode.LEACH_FILL, SimulationMode.WITHDRAWAL]
            and self.brine_injection_sg is None
        ):
            raise SansmicConfigError(
                "Missing required 'brine_injection_sg' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )

        # Validate casing size information
        if (
            self.simulation_mode
            in [
                SimulationMode.ORDINARY,
                SimulationMode.LEACH_FILL,
                SimulationMode.WITHDRAWAL,
                SimulationMode.STORAGE_FILL,
            ]
            and self.inner_tbg_inside_diam is None
        ):
            raise SansmicConfigError(
                "Missing required 'inner_tbg_inside_diam' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode
            in [
                SimulationMode.ORDINARY,
                SimulationMode.LEACH_FILL,
                SimulationMode.WITHDRAWAL,
                SimulationMode.STORAGE_FILL,
            ]
            and self.inner_tbg_outside_diam is None
        ):
            raise SansmicConfigError(
                "Missing required 'inner_tbg_outside_diam' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode in [SimulationMode.ORDINARY, SimulationMode.LEACH_FILL]
            and self.outer_csg_inside_diam is None
        ):
            raise SansmicConfigError(
                "Missing required 'outer_csg_inside_diam' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode in [SimulationMode.ORDINARY, SimulationMode.LEACH_FILL]
            and self.outer_csg_outside_diam is None
        ):
            raise SansmicConfigError(
                "Missing required 'outer_csg_outside_diam' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )

        # Product injection depth an rate validation
        if (
            self.simulation_mode in [SimulationMode.LEACH_FILL, SimulationMode.STORAGE_FILL]
            and self.product_injection_depth is None
        ):
            self.product_injection_depth = 0.0
            warnings.warn(
                "Missing 'product_injection_depth' options for {} simulation mode. Setting to 0.0".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode in [SimulationMode.LEACH_FILL, SimulationMode.STORAGE_FILL]
            and self.product_injection_rate is None
        ):
            raise SansmicConfigError(
                "Missing required 'product_injection_rate' options for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )

    def to_dict(self) -> Dict[str, Any]:
        """Convert the stage definition to a dictionary"""
        ret = asdict(self)
        return ret

    def _to_cstage(self) -> _ext.CStage:
        """Create a CStage object for the C++ interface."""
        self.validate()
        stage = _ext.CStage()
        stage.title = self.title
        stage.mode = _ext.CRunMode(int(self.simulation_mode))
        stage.print_interval = self.save_frequency
        stage.subsequent = 0 if self.set_initial_conditions else 1
        stage.iResetGeo = 0
        stage.rest_duration = self.rest_duration
        stage.stop_value = (
            int(self.stop_condition) * abs(self.stop_value) if self.stop_condition else 0
        )
        stage.injection_depth = (
            self.brine_injection_depth if self.brine_injection_depth is not None else 0
        )
        stage.production_depth = (
            self.brine_production_depth if self.brine_production_depth is not None else 0
        )
        stage.interface_depth = (
            self.set_interface_level if self.set_interface_level is not None else 0
        )
        stage.injection_rate = self.brine_injection_rate
        stage.inn_tbg_inside_radius = self.inner_tbg_inside_diam / 2.0
        stage.inn_tbg_outside_radius = self.inner_tbg_outside_diam / 2.0
        stage.out_csg_inside_radius = (
            self.outer_csg_inside_diam / 2.0
            if self.outer_csg_inside_diam
            else self.inner_tbg_inside_diam / 2.0
        )
        stage.out_csg_outside_radius = (
            self.outer_csg_outside_diam / 2.0
            if self.outer_csg_outside_diam
            else self.inner_tbg_outside_diam / 2.0
        )
        stage.injection_fluid_sg = self.brine_injection_sg
        stage.cavern_sg = (
            0.0
            if not self.set_initial_conditions or self.set_cavern_sg is None
            else self.set_cavern_sg
        )
        stage.timestep = self.solver_timestep
        stage.injection_duration = self.injection_duration
        stage.fill_rate = (
            self.product_injection_rate if self.product_injection_rate is not None else 0.0
        )
        stage.tDelay = 0.0
        return stage


@dataclass
class Scenario:
    """A SANSMIC scenario definition used to run a simulation."""

    title: str = None
    """General title for the scenario."""
    num_cells: int = None
    """Number of cells to use in the model."""
    geometry_data: Union[List[float], List[List[float]], str] = None
    """The geometry data in a list or list of lists, or a filename."""
    geometry_format: GeometryFormat = GeometryFormat.RADIUS_LIST
    """The format for the geometry data."""
    cavern_height: float = None
    """The total cavern height."""
    floor_depth: float = None
    """The depth from surface datum of the floor."""
    ullage_standoff: float = 40
    """The standoff distance of the lowest brine interface depth."""
    insolubles_ratio: float = 0.04
    """The volume ratio of insoluble material within the salt."""
    dissolution_factor: float = 1.0
    """**Danger** modify the dissolution factor by a certain amount; leave this as 1.0."""
    coallescing_wells: int = 1
    """Number of coallescing wells for cavern development modeling; leave as 1 for completed caverns."""
    well_separation: float = 0.0
    """The separation distance between wells for cavern development modeling."""
    units: Units = Units.FT_IN_BBL
    """The units used in describing the scenario."""
    stages: List[StageDefinition] = field(default_factory=list)
    """The activity stages to simulate."""

    defaults: Dict[str, Union[int, float]] = field(default_factory=dict)
    """Default values for a subset of stage attributes, see :attr:`StageDefinition.valid_default_keys`."""

    def __post_init__(self):
        converted_stages = list()
        for stage in self.stages:
            if isinstance(stage, dict):
                stage = StageDefinition(**stage, defaults=self.defaults)
            converted_stages.append(stage)
        self.stages = converted_stages

    def __setattr__(self, name, value):
        if isinstance(value, str) and value.strip() == "":
            value = None
        if name == "geometry_format" and not isinstance(value, GeometryFormat):
            if isinstance(value, int):
                value = GeometryFormat(value)
            elif isinstance(value, str):
                value = GeometryFormat[value.upper().replace(" ", "_").replace("-", "_")]
            elif isinstance(value, _ext.CGeometryFormat):
                value = GeometryFormat(int(value))
            else:
                SansmicConfigError("geometry_format cannot be of type {}".format(type(value)))
        elif name == "units" and not isinstance(value, Units):
            if isinstance(value, int):
                value = Units(value)
            elif isinstance(value, str):
                value = Units[value.upper().replace(" ", "_").replace("-", "_")]
            else:
                SansmicConfigError("Units cannot be of type {}".format(type(value)))
        super().__setattr__(name, value)

    def _to_cmodel(self, prefix="temp"):
        """Create a C++ model object; in general, this should only be called internally."""
        cmodel = _ext.CModel(prefix)
        for stage_num, stage in enumerate(self.stages):
            if stage_num == 0:
                stage.set_initial_conditions = True
            cstage = stage._to_cstage()
            cstage.num_cells = self.num_cells
            cstage.geometry_format = _ext.CGeometryFormat(int(self.geometry_format))
            cstage.num_coallescing = self.coallescing_wells
            cstage.cavern_height = self.cavern_height
            cstage.ullage_standoff = self.ullage_standoff
            cstage.coallescing_well_separation = self.well_separation
            cstage.dissolution_factor = self.dissolution_factor
            cstage.insoluble_fraction = self.insolubles_ratio
            cstage.refDep = self.floor_depth
            cstage.depth = self.floor_depth
            if stage_num == 0:
                # Add the geometry data
                if self.geometry_format is GeometryFormat.RADIUS_LIST and isinstance(
                    self.geometry_data, str
                ):
                    radii = list()
                    with open(self.geometry_data, "r") as fin:
                        for line in fin.readlines():
                            if len(line.strip()) > 0:
                                radii.append(float(line.strip()))
                    cstage.radius_vector = radii
                elif isinstance(self.geometry_data, dict):
                    if "radii" in self.geometry_data:
                        tmp = self.geometry_data["radii"].copy()
                        tmp.insert(0, 0)
                        cstage.radius_vector = tmp
                        cstage.nData = len(self.geometry_data["radii"])
                    if "depths" in self.geometry_data:
                        tmp = self.geometry_data["depths"]
                        tmp.insert(0, 0)
                        cstage.depth_vector = tmp
                        cstage.nData = len(self.geometry_data["depths"])
                    if "volumes" in self.geometry_data:
                        tmp = self.geometry_data["volumes"]
                        tmp.insert(0, 0)
                        cstage.volume_vector = tmp
                        cstage.nData = len(self.geometry_data["volumes"])
            # add the stage to the c++ model
            cmodel.add_stage(cstage)
        return cmodel

    def new_stage(self, pos: int = None, **kwargs) -> StageDefinition:
        """Add a new stage to the scenario in the optionally specified `pos` position
        based on keyword arguments. Passes :attr:`Scenario.defaults` unless a separate
        `defaults` dictionary is passed as one of the keyword arguments.

        Parameters
        ----------
        pos : int or None, keyword only
            The position in the stages list to insert the stage,
            by default None which will append to the end.
        kwargs : keyword arguments
            Any valid keyword argument for a StageDefinition object can be passed.

        Returns
        -------
        StageDefinition
            new stage created from the keyword arguments which has already been
            inserted into the proper location
        """
        defaults = kwargs.pop("defaults", self.defaults)
        stage = StageDefinition(defaults=defaults, **kwargs)
        if pos is None:
            self.stages.append(stage)
        else:
            self.stages.insert(pos, stage)
        return stage

    def to_dict(self) -> Dict[str, Any]:
        """Convert the scenario definition to a dictionary"""
        ret = asdict(self)
        return ret

    def new_simulation(self, prefix="temp") -> "Simulator":
        """Create a new :class:`Simulator` object.

        Parameters
        ----------
        prefix : str
            The prefix to use when creating output files.
        """
        return Simulator(self, prefix)


class StageIterator:
    """Iterator over each stage within an open simulation.

    Examples
    --------
    >>> stages: StageIterator = sim.stages()
    >>> for stage in stages:
    ...     pass
    >>> res = sim.get_results()
    """

    def __init__(self, sim: "Simulator"):
        self._sim = sim
        self.__model_stages = self._sim.num_stages
        self.__current_stage = 0

    def __iter__(self):
        self.__model_stages = self._sim.num_stages
        self.__current_stage = 0
        return self

    def __next__(self):
        self._sim._has_run = True
        if self.__current_stage < self.__model_stages:
            self._sim._run_stage(self.__current_stage)
            self.__current_stage = self.__current_stage + 1
            return self.__current_stage
        else:
            raise StopIteration


class StepIterator:
    """Iterator over each step within an open simulation.

    Examples
    --------
    >>> steps: StepIterator = sim.steps()
    >>> for step in steps:
    ...     if step % 100 == 0:
    ...         step_res = sim.get_current_state()
    >>> res = sim.get_results()
    """

    def __init__(self, sim: "Simulator"):
        self._sim = sim
        self.__model_stages = self._sim.num_stages
        self.__current_stage = 0
        self.__current_step = 0

    def __iter__(self):
        self.__model_stages = self._sim.num_stages
        self.__current_stage = 0
        self.__current_step = 0
        return self

    def __next__(self):
        self._sim._has_run = True
        if self.__current_stage < self.__model_stages:
            cs = self._sim._run_step(self.__current_stage)
            self.__current_stage = cs
            self.__current_step = self.__current_step + 1
            return self.__current_step
        else:
            raise StopIteration

    @property
    def stage(self) -> int:
        """The current stage number."""
        return self.__current_stage


class Simulator:
    """The SANSMIC Python simulator.

    Parameters
    ----------
    scenario : Scenario
        The scenario to run
    prefix : str
        The output file prefix to use


    Examples
    --------
    Like Python file objects, the Simulator class can be used with a context
    manager -- the "with-as" syntax -- or by explicitly calling the
    :meth:`~Simulator.open` and :meth:`~Simulator.close` methods.

    """

    def __init__(self, scenario: Union[Scenario, _ext.CModel], prefix="temp"):
        self._scenario = None
        self._prefix = prefix
        self._cmodel = None
        if isinstance(scenario, Scenario):
            self._scenario = scenario
        else:
            raise TypeError("Invalid scenario type")
        self._has_run = False
        self._is_open = False
        self._is_initialized = False
        self._is_finalized = False

    def __enter__(self):
        self.open(self._prefix)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        pass

    def __del__(self):
        if self._is_open or self._cmodel is not None:
            try:
                self._cmodel.close_outfiles()
            except Exception as e:
                warnings.warn(e)
        del self._cmodel
        self._cmodel = None

    @property
    def is_open(self) -> bool:
        """Is this simulator attached to a C++ SANSMIC model?"""
        return self._is_open

    @property
    def current_stage(self) -> int:
        """The current stage of the simulation"""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        return self._cmodel._get_current_stage()

    @property
    def current_time(self) -> float:
        """The current time (within the stage) of the simulation"""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        return self._cmodel._get_current_time()

    @property
    def current_volume(self) -> float:
        """The current simulated cavern volume"""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        return self._cmodel._get_current_volume()

    @property
    def is_running(self):
        """Is the simulation mid-stage?"""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        return self._cmodel.is_running()

    @property
    def num_stages(self):
        """The total number of stages that have been defined"""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        return self._cmodel.num_stages()

    def _initialize(self):
        """Initialize the simulation (for iterators)"""
        if self._is_initialized:
            raise RuntimeError("The simulation has already been initialized")
        self._cmodel.initialize()
        self._cmodel.open_outfiles(False)
        self._is_initialized = True
        self._has_run = False
        self._is_finalized = False

    def _run_stage(self, stage_num: int):
        """Run the specified stage, all steps."""
        stage_num = self._cmodel.run_stage(stage_num)
        return stage_num

    def _run_step(self, stage_num: int):
        """Run the next timestep in the specified stage."""
        step_num = self._cmodel.run_step(stage_num)
        return step_num

    def _finalize(self):
        """Finalize the simulation (for iterators)"""
        if self._is_finalized:
            raise RuntimeError("The simulation has already been finalized")
        self._cmodel.close_outfiles()
        self._is_finalized = True
        self._is_initialized = True

    def open(self, prefix=None):
        """Open a simulation.

        Parameters
        ----------
        prefix : str
            The output filename prefix to use, by default None which
            uses the value specified when the Simulator was created.
        """
        if self._is_open or self._cmodel is not None:
            raise RuntimeError("The simulation is already open")
        if prefix is None:
            prefix = self._prefix
        else:
            self._prefix = prefix
        if self._cmodel is None and self._scenario is not None:
            self._cmodel = self._scenario._to_cmodel(prefix)
        self._is_open = True

    def close(self):
        """Close and garbage collect the C++ model object."""
        if not self._is_finalized and self._is_initialized:
            self._finalize()
        self._cmodel.close_outfiles()
        del self._cmodel
        self._cmodel = None
        self._is_open = False
        self._is_finalized = False
        self._is_initialized = False

    def run(self):
        """Run the complete simulation; requires the Simulator to have been opened first."""
        if self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        self._is_finalized = False
        self._is_initialized = True
        self._cmodel.run()
        self._has_run = True
        self._is_initialized = True
        self._is_finalized = True

    def stages(self) -> StageIterator:
        """Get a generator that will run each stage in turn."""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        self._initialize()
        return StageIterator(self)

    def steps(self) -> StepIterator:
        """Get a generator that will run each step of each stage in turn."""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        self._initialize()
        return StepIterator(self)

    def get_current_state(self) -> "Results":
        """Get the current state of the model as a single-timestep results object."""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        data = self._cmodel._get_current_state()
        results = Results(data)
        return results

    def get_results(self) -> "Results":
        """Get the results for an entire simulation."""
        if self._cmodel is None:
            raise RuntimeError(
                "The simulation is already closed; results must be retrieved prior to closing the simulation"
            )
        if not self._has_run:
            raise RuntimeError("No results exist - the simulation has not been run/stages/steps.")
        data = self._cmodel._get_results()
        results = Results(data)
        return results


class Results:
    """
    Python container for the results of a sansmic simulation.

    Parameters
    ----------
    data : :class:`sansmic.library.libsansmic.CResults`
        results object from the C++ library

    """

    def __init__(self, data: _ext.CResults) -> None:
        self._data = data
        self._t_index = np.arange(1, len(self._data.t) + 1)
        self._n_index = np.arange(1, len(self._data.r[0]) + 1)

    def __repr__(self):
        return "<Results: from {:.3g} to {:.3g} d ({:d} stages)>".format(
            self._data.t[0] / 24.0,
            self._data.t[-1] / 24.0,
            max(self._data.stage) - min(self._data.stage) + 1,
        )

    def get_timeseries_data(self) -> pd.DataFrame:
        """Get the timeseries of summary results as a dataframe"""
        return pd.DataFrame.from_dict(
            {
                "time": self.t(),
                "days": self.t() / 24.0,
                "convergence": self.err(),
                "cavern volume": self.V_cavTot(),
                "insolubles volume": self.V_insolTot(),
                "insolubles vented": self.V_insolVent(),
                "water injected": self.V_injTot(),
                "oil injected": self.V_fillTot(),
                "output sg": self.sg_out(),
                "cavern sg": self.sg_cavAve(),
                "insolubles depth": self.z_insol(),
                "injection depth": self.z_inj(),
                "production depth": self.z_prod(),
                "plume depth": self.z_plm(),
                "obi depth": self.z_obi(),
                "production rate": self.Q_out(),
                "stage": self.stage(),
                "phase": self.phase(),
                "time step": self.dt(),
                "injection cell": self.injCell(),
                "production cell": self.prodCell(),
                "plume cell": self.plmCell(),
                "obi cell": self.obiCell(),
                "jet length": self.l_jet(),
                "jet velocity": self.u_jet(),
                "jet radius": self.r_jet(),
            }
        )

    def z_0(self) -> pd.Series:
        """node depth, in ft"""
        return pd.Series(self._data.z_0, index=self._n_index, name="depth")

    def h_0(self) -> pd.Series:
        """node height above original floor, in ft"""
        return pd.Series(self._data.h_0, index=self._n_index, name="height")

    def r_0(self) -> pd.Series:
        """initial cavern radius, in ft"""
        return pd.Series(self._data.r_0, index=self._n_index, name="radius")

    def t(self) -> pd.Series:
        """time of each row of time-based results"""
        return pd.Series(self._data.t, index=self._t_index, name="time")

    def stage(self) -> pd.Series:
        """stage number"""
        return pd.Series(self._data.stage, index=self._t_index, dtype=int, name="stage")

    def phase(self) -> pd.Series:
        """phase (0=static, 1=dynamic)"""
        return pd.Series(self._data.phase, index=self._t_index, dtype=int, name="phase")

    def injCell(self) -> pd.Series:
        """injection cell index"""
        return pd.Series(self._data.injCell, index=self._t_index, dtype=int, name="injection cell")

    def prodCell(self) -> pd.Series:
        """production cell index"""
        return pd.Series(
            self._data.prodCell, index=self._t_index, dtype=int, name="production cell"
        )

    def obiCell(self) -> pd.Series:
        """cell containing oil blanket interface"""
        return pd.Series(self._data.obiCell, index=self._t_index, dtype=int, name="obi cell")

    def plmCell(self) -> pd.Series:
        """cell containing plume stagnation level"""
        return pd.Series(self._data.plumCell, index=self._t_index, dtype=int)

    def err(self) -> pd.Series:
        """mass ballance error calculation"""
        return pd.Series(self._data.err, index=self._t_index, name="convergence")

    def z_obi(self) -> pd.Series:
        """depth of blanket interface"""
        return pd.Series(self._data.z_obi, index=self._t_index, name="obi depth")

    def z_inj(self) -> pd.Series:
        """depth of injection cell"""
        return pd.Series(self._data.z_inj, index=self._t_index, name="injection depth")

    def z_prod(self) -> pd.Series:
        """depth of production cell"""
        return pd.Series(self._data.z_prod, index=self._t_index, name="production depth")

    def z_plm(self) -> pd.Series:
        """depth of plume stagnation cell"""
        return pd.Series(self._data.z_plm, index=self._t_index, name="plume depth")

    def z_insol(self) -> pd.Series:
        """depth of top of insolubles"""
        return pd.Series(self._data.z_insol, index=self._t_index, name="insolubles depth")

    def h_insol(self) -> pd.Series:
        """height of insolubles above original floor"""
        return pd.Series(self._data.h_insol, index=self._t_index, name="insolubles height")

    def l_jet(self) -> pd.Series:
        """length of injection jet"""
        return pd.Series(self._data.l_jet, index=self._t_index, name="jet length")

    def r_jet(self) -> pd.Series:
        """velocity of injection jet"""
        return pd.Series(self._data.r_jet, index=self._t_index, name="jet radius")

    def u_jet(self) -> pd.Series:
        """radius of injection jet"""
        return pd.Series(self._data.u_jet, index=self._t_index, name="jet velocity")

    def V_injTot(self) -> pd.Series:
        """total injected water volume (bbl)"""
        return pd.Series(self._data.V_injTot, index=self._t_index, name="water injected")

    def V_fillTot(self) -> pd.Series:
        """total injected (+) or withdrawn (-) oil volume (bbl)"""
        return pd.Series(self._data.V_fillTot, index=self._t_index, name="oil injected")

    def V_cavTot(self) -> pd.Series:
        """total cavern volume (bbl)"""
        return pd.Series(self._data.V_cavTot, index=self._t_index, name="cavern volume")

    def V_insolTot(self) -> pd.Series:
        """total volume of insolubles (bbl)"""
        return pd.Series(self._data.V_insolTot, index=self._t_index, name="insolubles volume")

    def V_insolVent(self) -> pd.Series:
        """volume of insolubles vented (bbl)"""
        return pd.Series(self._data.V_insolVent, index=self._t_index, name="insolubles vented")

    def Q_out(self) -> pd.Series:
        """current brine exiting cavern (bbl/d)"""
        return pd.Series(self._data.Q_out, index=self._t_index, name="production rate")

    def sg_out(self) -> pd.Series:
        """current brine concentration exiting cavern (sg)"""
        return pd.Series(self._data.sg_out, index=self._t_index, name="output sg")

    def sg_cavAve(self) -> pd.Series:
        """current average concentration in brine (sg)"""
        return pd.Series(self._data.sg_cavAve, index=self._t_index, name="cavern sg")

    def dt(self) -> pd.Series:
        """the calculation timestep used"""
        return pd.Series(self._data.dt, index=self._t_index, name="time step")

    def r(self) -> pd.DataFrame:
        """current radius (ft)"""
        return pd.DataFrame(self._data.r, index=self._t_index, columns=self._n_index)

    def dr(self) -> pd.DataFrame:
        """change in radius since start of simulation (ft)"""
        return pd.DataFrame(self._data.dr_0, index=self._t_index, columns=self._n_index)

    def sg(self) -> pd.DataFrame:
        """current concentration (sg)"""
        return pd.DataFrame(self._data.sg, index=self._t_index, columns=self._n_index)

    def theta(self) -> pd.DataFrame:
        """current wall angle (deg)"""
        return pd.DataFrame(self._data.theta, index=self._t_index, columns=self._n_index)

    def Q_inj(self) -> pd.DataFrame:
        """current inflow (bbl/d)"""
        return pd.DataFrame(self._data.Q_inj, index=self._t_index, columns=self._n_index)

    def V(self) -> pd.DataFrame:
        """current volume of cell defined by node (bbl)"""
        return pd.DataFrame(self._data.V, index=self._t_index, columns=self._n_index)

    def f_dis(self) -> pd.DataFrame:
        """current dissolution factor, including adjustments"""
        return pd.DataFrame(self._data.f_dis, index=self._t_index, columns=self._n_index)

    def f_flag(self) -> pd.DataFrame:
        """current dissolution process adjustments flag"""
        return pd.DataFrame(
            self._data.f_flag, index=self._t_index, columns=self._n_index, dtype=int
        )

    def xincl(self) -> pd.DataFrame:
        """current wall angle adjustment factor"""
        return pd.DataFrame(self._data.xincl, index=self._t_index, columns=self._n_index)

    def amd(self) -> pd.DataFrame:
        """debug variable"""
        return pd.DataFrame(self._data.amd, index=self._t_index, columns=self._n_index)

    def D_coeff(self) -> pd.DataFrame:
        """current diffusion coefficient in use at node depth"""
        return pd.DataFrame(self._data.D_coeff, index=self._t_index, columns=self._n_index)

    def dC_dz(self) -> pd.DataFrame:
        """vertical concentration gradient (sg)"""
        return pd.DataFrame(self._data.dC_dz, index=self._t_index, columns=self._n_index)

    def C_old(self) -> pd.DataFrame:
        """previous step concentration (sg)"""
        return pd.DataFrame(self._data.C_old, index=self._t_index, columns=self._n_index)

    def C_new(self) -> pd.DataFrame:
        """current concentration (sg)"""
        return pd.DataFrame(self._data.C_new, index=self._t_index, columns=self._n_index)

    def dC_dt(self) -> pd.DataFrame:
        """change in concentration over time"""
        return pd.DataFrame(self._data.dC, index=self._t_index, columns=self._n_index)

    def dr_dt(self) -> pd.DataFrame:
        """current wall recession rate (ft)"""
        return pd.DataFrame(self._data.dr, index=self._t_index, columns=self._n_index)

    def C_plm(self) -> pd.DataFrame:
        """concentration within plume at node depth"""
        return pd.DataFrame(self._data.C_plm, index=self._t_index, columns=self._n_index)

    def u_plm(self) -> pd.DataFrame:
        """plume velocity at node depth"""
        return pd.DataFrame(self._data.u_plm, index=self._t_index, columns=self._n_index)

    def r_plm(self) -> pd.DataFrame:
        """plume radius at node depth"""
        return pd.DataFrame(self._data.r_plm, index=self._t_index, columns=self._n_index)


class _OutDataBlock:
    def __init__(self):
        self.t = None
        self.dt = None
        self.t_prime = None
        self.phase = None
        self.stage = None
        self.injCell = None
        self.prodCell = None
        self.obiCell = None
        self.plumCell = None
        self.z_inj = None
        self.z_prod = None
        self.z_obi = None
        self.z_plm = None
        self.l_jet = None
        self.r_jet = None
        self.u_jet = None
        self.idx = list()
        self.h = list()
        self.r = list()
        self.dr_0 = list()
        self.sg = list()
        self.theta = list()
        self.Q_inj = list()
        self.V = list()
        self.f_dis = list()
        self.f_flag = list()
        self.xincl = list()
        self.amd = list()
        self.D_coeff = list()
        self.dC_dz = list()
        self.C_old = list()
        self.C_new = list()
        self.dC = list()
        self.dr = list()
        self.C_plm = list()
        self.u_plm = list()
        self.r_plm = list()
        self.V_cavTot = None
        self.Q_out = None
        self.sg_out = None
        self.V_insolTot = None
        self.z_insol = None
        self.z_obi = None
        self.V_insolVent = None
        self.V_brine = None
        self.V_ullage = None
        self.V_usable = None
        self.err = None


class _OutputData(_OutDataBlock):
    def __init__(self):
        super().__init__()
        self.timeTotal = None
        self.r_0 = None
        self.z_0 = None
        self.h_0 = None
        self.V_injTot = None
        self.V_fillTot = None
        self.sg_cavAve = None
