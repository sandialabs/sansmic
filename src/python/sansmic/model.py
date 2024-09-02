# coding: utf-8
#
# Copyright (c) 2024 National Technology and Engineering Solutions of
# Sandia, LLC. Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.
#
# SPDX-License-Identifier: BSD-3-Clause.

"""The core python classes for SANSMIC."""

import logging
import warnings
from copy import deepcopy
from dataclasses import InitVar, asdict, dataclass, field
from enum import IntEnum
from fractions import Fraction
from math import nan
from types import NoneType
from typing import Any, Dict, List, Union

import numpy as np
import pandas as pd

logger = logging.getLogger("sansmic")

has_ext = False
try:
    from . import libsansmic as _ext

    has_ext = True


except ImportError:
    logger.critical(
        "The C++ library is not installed. Conversions will work, but the main program will not run."
    )
    has_ext = False

    class _ext:
        class _CModel:
            pass

        class _CStage:
            pass

        class _CResults:
            pass

        class _CGeometryFormat(IntEnum):
            blah = 1

        class _CRunMode(IntEnum):
            blah = 1

        CStage = _CStage
        CModel = _CModel
        CGeometryFormat = _CGeometryFormat
        CRunMode = _CRunMode
        CResults = _CResults


def _rename_with_underscore(orig: dict, new: dict):
    """Rename all keys in a dictionary with underscores instead of spaces or dashes"""
    for k, v in orig.items():
        k = k.lower().replace("-", "_")
        new[k] = v


def _rename_with_dash(orig: dict, new: dict):
    """Rename all keys in a dictionary with dashes instead of underscores"""
    for k, v in orig.items():
        k = k.lower().replace("_", "-")
        new[k] = v


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
        """1 in ≔ 0.0254 m"""
        return Fraction(254, 10000)

    @property
    def foot(self):
        """1 ft ≔ 0.3048 m"""
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
        """1 cm ≔ 0.01 m"""
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
class AdvancedOptions:
    """Advanced configuration options. Most of these can and should be left as None which
    will use the default C++ library values."""

    absolute_error: float = None  #  1.0e-2
    """ODE solver absolute error tolerance; CModel default is 0.01"""
    relative_error: float = None  #  1.0e-4
    """ODE solver absolute error tolerance; CModel default is 0.0001"""
    jet_model_version: int = None  #  1
    """Jet model version; CModel default is 1"""
    plume_model_version: int = None  #  1
    """Plume model version; CModel default is 1"""
    temperature_model_version: int = None  #  0
    """Temperature model version; CModel default is 0"""
    dissolution_factor: float = None  # 1.0
    """Dissolution factor; CModel default is 1.0 - this should not be changed unless you are sure you know the effects"""
    max_brine_sg: float = None  #  1.2019
    """Maximum brine specific gravity; CSalt default is 1.2019"""
    solid_density: float = None  #  2.16
    """Rock density in solid form; CSalt default is 2.16 g/cc"""
    entrainment_coeff: float = None  #  0.09
    """Dissolution entrainment coefficient; CModel default is 0.09"""
    molecular_diffusion: float = None  #  5.03e-5
    """Molecular diffusion coefficient; CModel default is 5.03e-5"""
    eddy_coefficient: float = None  #  1.142e5
    """Eddy coefficient; CModel default is 1.142e5"""
    diffusion_beta: float = None  #  0.147
    """Diffusion beta coefficient; default is 0.147"""

    @classmethod
    def from_dict(cls, opts: dict) -> "AdvancedOptions":
        """Create a new object from a dictionary of options.

        This method differs from the __init__ constructor by automatically
        converting non-underscore characters - e.g., ``-`` or ``.`` or `` ``
        to underscores and changing keys to lower-case prior to creating
        the object.

        See also: :meth:`~StageDefinition.to_dict`

        Parameters
        ----------
        opts : dict
            The initialization values

        Returns
        -------
        StageDefinition
            the new stage
        """
        new_opts = dict()
        _rename_with_underscore(opts, new_opts)
        return cls(**new_opts)

    def to_dict(self, keep_empty: bool = False):
        """Convert the object's data to a dictionary of options.

        This method differs from the :meth:`~dataclasses.asdict`
        method by automatically converting underscore characters to
        hyphens for a more readable dictionary. Specifically used
        when creating TOML, JSON and YAML files.

        See also: :meth:`~StageDefinition.from_dict`

        Parameters
        ----------
        keep_empty : bool
            Keep values set to None and set them to ``""``

        Returns
        -------
        dict
            the options dictionary
        """
        ret = dict()
        _rename_with_dash(asdict(self), ret)
        keys = list(ret.keys())
        for k in keys:
            if ret[k] is None:
                if keep_empty:
                    ret[k] = ""
                else:
                    del ret[k]
        return ret


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
    """The meaning of this field is based on the type of data that is provided (a value of 0 means off)."""
    # """
    # float : constant rate
    #     Injection occurs at a constant rate for the duration of the
    #     :attr:`injection_duration`.

    # str : file path
    #     Injection rates should be read from a file. The file can be any type that can
    #     be read into a pd DataFrame. An 'index' column is optional, but if not provided,
    #     then all data is assumed to be pre-sorted. In addition to the
    #     'index' column, files may contain exactly one set of the following columns:

    #     ['rate', 'duration']
    #         This style gives a list of rates and how long they should continue for.
    #         Durations are provided in decimal hours. If the sum of durations is greater
    #         than :attr:`injection_duration`, an error will be raised. If the
    #         sum of durations is less than :attr:`injection_duration`, the remaining
    #         time will be assumed to be no-flow.

    #     ['time', 'rate']
    #         This style gives a list of new rate values that start at the given time
    #         after the start of the stage. The 'time' values must be increasing and
    #         never duplicate. Times are in decimal hours since the start of the stage.
    #         The last rate will be used up until the end of the
    #         :attr:`injection_duration`. If a time is given that is after the
    #         :attr:`injection_duration`, then an error will be raised.

    #     ['hourly']
    #         This style specifies a new rate for each hour of the stage. There must be
    #         a value specified for each hour of the :attr`injection_duration`, or an
    #         error will be raised.

    #     The optional column ['sg'] can be provided to change the injection water/brine
    #     specific gravity; however, if provided, it must contain a value for every row.

    #     If 'duration' and 'time' occur in the same file, an error will be raised.
    #     If 'time', 'rate', or 'duration' occur in the same file as 'hourly', then
    #     an error will be raised.

    #     Columns with names other than those described above will be ignored.

    #     .. attention::

    #         In all cases, a no-flow, leaching-only period will occur for
    #         :attr:`rest_duration` hours **after** all injection data
    #         is processed -- this means that if an injection file finishes with 0s
    #         that there will be an additional :attr:`rest_duration` hours added to the
    #         end of the stage.

    #     .. warning::

    #         The same file cnnot be used for both brine injection and product
    #         injection data.
    # """

    set_initial_conditions: bool = None
    """Unlink initial cavern brine gravity and interface level from previous stage.
    Automatically set to True for the first stage added to a model."""
    set_cavern_sg: float = None
    """Set the initial specific gravity for all brine-filled cells of the cavern to
    this value. If set_initial_conditions is False and this is not None (or 0) an error
    will be raised upon scenario validation."""
    brine_interface_depth: float = None
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
                    "Ignoring non-defaultable or unknown setting {} = {}".format(
                        k, repr(v)
                    )
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
                value = SimulationMode[
                    value.upper().replace(" ", "_").replace("-", "_")
                ]
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
                TypeError("stop_condition cannot be of type {}".format(type(value)))
        elif name == "set_cavern_sg" and value is not None and value < 1.0:
            value = None
        elif name == "set_initial_conditions" and value is False:
            self.set_cavern_sg = None
            self.brine_interface_depth = None
        super().__setattr__(name, value)

    @classmethod
    def from_dict(cls, opts: dict) -> "StageDefinition":
        """Create a new object from a dictionary of options.

        This method differs from the __init__ constructor by automatically
        converting non-underscore characters - e.g., ``-`` or ``.`` or `` ``
        to underscores and changing keys to lower-case prior to creating
        the object.

        See also: :meth:`~StageDefinition.to_dict`

        Parameters
        ----------
        opts : dict
            The initialization values

        Returns
        -------
        StageDefinition
            the new stage
        """
        new_opts = dict()
        _rename_with_underscore(opts, new_opts)
        return cls(**new_opts)

    def to_dict(self, keep_empty: bool = False):
        """Convert the object's data to a dictionary of options.

        This method differs from the :meth:`~dataclasses.asdict`
        method by automatically converting underscore characters to
        hyphens for a more readable dictionary. Specifically used
        when creating TOML, JSON and YAML files.

        See also: :meth:`~StageDefinition.from_dict`

        Parameters
        ----------
        keep_empty : bool
            Keep values set to None and set them to ``""``

        Returns
        -------
        dict
            the options dictionary
        """
        ret = dict()
        _rename_with_dash(asdict(self), ret)
        keys = list(ret.keys())
        for k in keys:
            if ret[k] is None:
                if keep_empty:
                    ret[k] = ""
                else:
                    del ret[k]
        return ret

    def validate(self):
        """Validate that all required options have been selected.

        Raises
        ------
        TypeError
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
            self.set_cavern_sg is None or self.brine_interface_depth is None
        ):
            raise TypeError(
                "An initial stage must have both set_cavern_sg and set_interface_level."
            )
        elif not self.set_initial_conditions and self.set_cavern_sg:
            # raise TypeError(
            warnings.warn(
                "Setting the starting cavern sg ought to use 'set_initial_conditions' to be set to True"
            )
        elif not self.set_initial_conditions and self.brine_interface_depth:
            warnings.warn(
                "Make sure you meant to reset the interface level; use 0.0 or None to continue from the last stage."
            )

        # Validate appropriate brine injection and production settings
        if (
            self.simulation_mode
            in [
                SimulationMode.ORDINARY,
                SimulationMode.LEACH_FILL,
                SimulationMode.STORAGE_FILL,
            ]
            and self.brine_production_depth is None
        ):
            raise TypeError(
                "Missing required 'brine_production_depth' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode
            in [
                SimulationMode.ORDINARY,
                SimulationMode.LEACH_FILL,
                SimulationMode.WITHDRAWAL,
            ]
            and self.brine_injection_depth is None
        ):
            raise TypeError(
                "Missing required 'brine_injection_depth' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode
            in [
                SimulationMode.ORDINARY,
                SimulationMode.LEACH_FILL,
                SimulationMode.WITHDRAWAL,
            ]
            and self.brine_injection_rate is None
        ):
            raise TypeError(
                "Missing required 'brine_injection_rate' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode
            in [
                SimulationMode.ORDINARY,
                SimulationMode.LEACH_FILL,
                SimulationMode.WITHDRAWAL,
            ]
            and self.brine_injection_sg is None
        ):
            raise TypeError(
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
            raise TypeError(
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
            raise TypeError(
                "Missing required 'inner_tbg_outside_diam' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode in [SimulationMode.ORDINARY, SimulationMode.LEACH_FILL]
            and self.outer_csg_inside_diam is None
        ):
            raise TypeError(
                "Missing required 'outer_csg_inside_diam' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode in [SimulationMode.ORDINARY, SimulationMode.LEACH_FILL]
            and self.outer_csg_outside_diam is None
        ):
            raise TypeError(
                "Missing required 'outer_csg_outside_diam' for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )

        # Product injection depth an rate validation
        if (
            self.simulation_mode
            in [SimulationMode.LEACH_FILL, SimulationMode.STORAGE_FILL]
            and self.product_injection_depth is None
        ):
            self.product_injection_depth = 0.0
            warnings.warn(
                "Missing 'product_injection_depth' options for {} simulation mode. Setting to 0.0".format(
                    self.simulation_mode.name
                )
            )
        if (
            self.simulation_mode
            in [SimulationMode.LEACH_FILL, SimulationMode.STORAGE_FILL]
            and self.product_injection_rate is None
        ):
            raise TypeError(
                "Missing required 'product_injection_rate' options for {} simulation mode.".format(
                    self.simulation_mode.name
                )
            )

    def _to_cstage(self) -> _ext.CStage:
        """Create a CStage object for the C++ interface.
        This method is protected because, in general, it is better not to
        try to operate directly on the C++ stage object.
        """
        self.validate()
        stage = _ext.CStage()
        stage.title = self.title
        stage.mode = _ext.CRunMode(int(self.simulation_mode))
        stage.print_interval = self.save_frequency
        stage.subsequent = 0 if self.set_initial_conditions else 1
        stage.rest_duration = self.rest_duration
        stage.stop_value = (
            int(self.stop_condition) * abs(self.stop_value)
            if self.stop_condition
            else 0
        )
        stage.injection_depth = (
            self.brine_injection_depth if self.brine_injection_depth is not None else 0
        )
        stage.production_depth = (
            self.brine_production_depth
            if self.brine_production_depth is not None
            else 0
        )
        stage.interface_depth = (
            self.brine_interface_depth if self.brine_interface_depth is not None else 0
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
            self.product_injection_rate
            if self.product_injection_rate is not None
            else 0.0
        )
        return stage


@dataclass
class Scenario:
    """A SANSMIC scenario definition used to run a simulation."""

    title: str = None
    """General title for the scenario."""
    comments: str = """"""
    """Comments about the scenario"""
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
    """**warning** modify the dissolution factor by a certain amount; leave this as 1.0."""
    coallescing_wells: int = 1
    """**warning** number of coallescing wells for cavern development modeling; leave as 1 for completed caverns."""
    well_separation: float = 0.0
    """**warning** the separation distance between wells for cavern development modeling."""
    units: Units = Units.FT_IN_BBL
    """The units used in describing the scenario."""
    defaults: Dict[str, Union[int, float]] = field(default_factory=dict)
    """Default values for a subset of stage attributes, see :attr:`StageDefinition.valid_default_keys`."""
    advanced: AdvancedOptions = field(default_factory=AdvancedOptions)
    """Advanced and/or uncommonly used options"""
    stages: List[StageDefinition] = field(default=None)
    """The activity stages to simulate."""

    def __setattr__(self, name, value):
        """The setattr method is overloaded. IntEnum parameters are
        automatically converted from strings or integers and stages and
        advanced options are converted to the appropriate classes if a
        dictionary is passed."""
        if isinstance(value, str) and value.strip() == "":
            value = None
        if name == "geometry_format" and not isinstance(value, GeometryFormat):
            if isinstance(value, int):
                value = GeometryFormat(value)
            elif isinstance(value, str):
                value = GeometryFormat[
                    value.upper().replace(" ", "_").replace("-", "_")
                ]
            elif isinstance(value, _ext.CGeometryFormat):
                value = GeometryFormat(int(value))
            else:
                TypeError("geometry_format cannot be of type {}".format(type(value)))
        elif name == "units" and not isinstance(value, Units):
            if isinstance(value, int):
                value = Units(value)
            elif isinstance(value, str):
                value = Units[value.upper().replace(" ", "_").replace("-", "_")]
            else:
                TypeError("Units cannot be of type {}".format(type(value)))
        elif (
            name == "advanced"
            and value is not None
            and not isinstance(value, AdvancedOptions)
        ):
            if isinstance(value, dict):
                value = AdvancedOptions.from_dict(value)
            else:
                raise TypeError("advanced cannot be a {}".format(type(value)))
        # elif name == "stages" and (not isinstance(value, list) or value is None):
        #     raise TypeError("stages cannot be a {}".format(type(value)))
        elif name == "stages":
            if not hasattr(self, "stages") or value is None:
                value = list()
            elif not hasattr(self, name) or self.stages is None:
                converted_stages = list()
                for stage in value:
                    if isinstance(stage, dict):
                        stage = StageDefinition(**stage, defaults=self.defaults)
                    converted_stages.append(stage)
                value = converted_stages
            else:
                TypeError("stages cannot set directly - please use add_stages")
        super().__setattr__(name, value)

    @classmethod
    def from_dict(cls, opts: dict) -> "Scenario":
        """Create a new object from a dictionary of options.

        This method differs from the __init__ constructor by automatically
        converting non-underscore characters - e.g., ``-`` or ``.`` or `` ``
        to underscores and changing keys to lower-case prior to creating
        the object.

        See also: :meth:`~StageDefinition.to_dict`

        Parameters
        ----------
        opts : dict
            The initialization values

        Returns
        -------
        StageDefinition
            the new stage
        """
        new_opts = dict()
        _rename_with_underscore(opts, new_opts)
        return cls(**new_opts)

    def to_dict(self, keep_empty: bool = False):
        """Convert the object's data to a dictionary of options.

        This method differs from the :meth:`~dataclasses.asdict`
        method by automatically converting underscore characters to
        hyphens for a more readable dictionary. Specifically used
        when creating TOML, JSON and YAML files.

        See also: :meth:`~StageDefinition.from_dict`

        Parameters
        ----------
        keep_empty : bool
            Keep values set to None and set them to ``""``

        Returns
        -------
        dict
            the options dictionary
        """
        ret = dict()
        _rename_with_dash(asdict(self), ret)
        ret["advanced"] = self.advanced.to_dict(keep_empty)
        if len(ret["advanced"]) == 0:
            del ret["advanced"]
        del ret["stages"]
        ret["stages"] = list()
        for stage in self.stages:
            ret["stages"].append(stage.to_dict(keep_empty))
        keys = list(ret.keys())
        for k in keys:
            if ret[k] is None:
                if keep_empty:
                    ret[k] = ""
                else:
                    del ret[k]
        return ret

    def new_stage(self, pos: int = None, **kwargs) -> StageDefinition:
        """Add a new stage in the optionally-specified `pos` position, and create it
        based on keyword arguments. Passes existing :attr:`~Scenario.defaults`
        unless a separate `defaults` dictionary is passed as one of the keyword arguments.

        Parameters
        ----------
        pos : int or None, keyword only
            The position in the stages list to insert the stage,
            by default None which will append to the end.
        kwargs : keyword arguments
            Any valid keyword argument for a :class:`StageDefinition` object

        Returns
        -------
        StageDefinition
            new stage created from the keyword arguments. It will have been
            added into the stage list in the proper position.
        """
        defaults = kwargs.pop("defaults", self.defaults)
        stage = StageDefinition(defaults=defaults, **kwargs)
        if pos is None:
            self.stages.append(stage)
        else:
            self.stages.insert(pos, stage)
        return stage

    def new_simulation(self, prefix="temp", verbosity=0) -> "Simulator":
        """Create a new :class:`Simulator` object.

        Parameters
        ----------
        prefix : str
            The prefix to use when creating output files, by default temp.
        verbosity : int
            A verbosity level to pass to the C++ model, by default 0.
        """
        return Simulator(self, prefix, verbosity)

    def _to_cscenario(self):
        """Create a C++ model object; in general, this should only be called internally."""
        cscenario = _ext.CScenario()
        cscenario.fraction_insolubles = self.insolubles_ratio
        cscenario.floor_depth = self.floor_depth
        cscenario.geometry_format = _ext.CGeometryFormat(int(self.geometry_format))
        cscenario.num_cells = self.num_cells
        cscenario.ullage_standoff = self.ullage_standoff
        cscenario.cavern_height = self.cavern_height
        cscenario.coallescing_wells = self.coallescing_wells
        cscenario.well_separation = self.well_separation
        if self.advanced.absolute_error is not None:
            cscenario.absolute_error = self.advanced.absolute_error

        if self.advanced.diffusion_beta is not None:
            cscenario.diffusion_beta = self.advanced.diffusion_beta

        if self.advanced.dissolution_factor is not None:
            cscenario.dissolution_factor = self.advanced.dissolution_factor

        if self.advanced.eddy_coefficient is not None:
            cscenario.eddy_coefficient = self.advanced.eddy_coefficient

        if self.advanced.entrainment_coeff is not None:
            cscenario.entrainment_coeff = self.advanced.entrainment_coeff

        if self.advanced.jet_model_version is not None:
            cscenario.jet_model_version = self.advanced.jet_model_version

        if self.advanced.max_brine_sg is not None:
            cscenario.max_brine_sg = self.advanced.max_brine_sg

        if self.advanced.molecular_diffusion is not None:
            cscenario.molecular_diffusion = self.advanced.molecular_diffusion

        if self.advanced.plume_model_version is not None:
            cscenario.plume_model_version = self.advanced.plume_model_version

        if self.advanced.relative_error is not None:
            cscenario.relative_error = self.advanced.relative_error

        if self.advanced.solid_density is not None:
            cscenario.solid_density = self.advanced.solid_density

        if self.advanced.temperature_model_version is not None:
            cscenario.temperature_model_version = (
                self.advanced.temperature_model_version
            )

        if self.geometry_format is GeometryFormat.RADIUS_LIST and isinstance(
            self.geometry_data, str
        ):
            radii = list()
            with open(self.geometry_data, "r") as fin:
                for line in fin.readlines():
                    if len(line.strip()) > 0:
                        radii.append(float(line.strip()))
            cscenario.radius_vector = radii
        elif isinstance(self.geometry_data, dict):
            if "radii" in self.geometry_data:
                tmp = self.geometry_data["radii"].copy()
                tmp.insert(0, 0)
                cscenario.radius_vector = tmp
            if "depths" in self.geometry_data:
                tmp = self.geometry_data["depths"]
                tmp.insert(0, 0)
                cscenario.depth_vector = tmp
            if "volumes" in self.geometry_data:
                tmp = self.geometry_data["volumes"]
                tmp.insert(0, 0)
                cscenario.volume_vector = tmp
        for stage_num, stage in enumerate(self.stages):
            if stage_num == 0:
                stage.set_initial_conditions = True
            cstage = stage._to_cstage()
            cscenario.add_stage(cstage)
        return cscenario


class StepwiseIterator:
    """Iterator over each step within an open simulation.

    Examples
    --------
    >>> for step in sim.steps:
    ...     if step % 100 == 0:
    ...         step_res = sim.get_current_state()
    >>> res = sim.results

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
            stage_done = self._sim._run_step()
            self.__current_stage = self.__current_stage + stage_done
            self.__current_step = self.__current_step + 1
            return self.__current_stage, self.__current_step
        else:
            raise StopIteration


class Simulator:
    """The SANSMIC Python simulator.

    Parameters
    ----------
    scenario : Scenario
        The scenario to run
    prefix : str, optional
        The output file prefix to use, by default "temp"


    Examples
    --------
    The simulator can be created using context manager style "with-as"
    statements to automatically close the C++ model when the run is completed.
    This can be accomplished by directly creating the simulator object
    from the class or by using the :meth:`~Scenario.new_simulation` method.

    *Running in batch mode without context management*

    .. code:: python

        sim = sansmic.Simulator(scenario, prefix='run13')
        sim.open()
        sim.run_sim()
        sim.close()
        res13 = sim.results


    *Run in batch mode with context management*

    .. code:: python

        with scenario.new_simulation('run14') as sim:
            sim.run_sim()
        res14 = sim.results


    *Run in stepwise mode with context management*

    .. code:: python

        with scenario.new_simulation('run15') as sim:
            for stage, step in sim.steps:
                pass
        res15 = sim.results

    """

    def __init__(
        self, scenario: Union[Scenario, _ext.CModel], prefix="temp", verbosity=0
    ):
        self._scenario = None
        self._prefix = prefix
        self._cmodel = None
        self._verbosity = verbosity
        if isinstance(scenario, Scenario):
            self._scenario = scenario
        else:
            raise TypeError("Invalid scenario")
        self._has_run = False
        self._is_open = False
        self._is_initialized = False
        self._is_finalized = False
        self.__stepwise = None
        self.__results = None

    def __enter__(self):
        self.open(self._prefix)
        self._cmodel.verbosity = self._verbosity
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @property
    def steps(self) -> StepwiseIterator:
        """Provides an iterator that will run each step of each stage in turn.

        Examples
        --------
        .. code:: python

            # Generate and step through the model
            with scenario.new_simulation('run14') as sim:
                for stage, step in sim.steps:
                    pass
            results = sim.results

        """
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        if self.__stepwise is not None:
            raise RuntimeError("The simulation is already running")
        self._initialize()
        self.__stepwise = StepwiseIterator(self)
        return self.__stepwise

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
    def is_running(self) -> bool:
        """Is the simulation mid-stage?"""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        return self._cmodel.is_running()

    @property
    def num_stages(self) -> int:
        """The total number of stages that have been defined"""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        return self._cmodel.num_stages()

    def _initialize(self):
        """Initialize the simulation (for iterators)"""
        if self._is_initialized:
            raise RuntimeError("The simulation has already been initialized")
        self._cmodel.open_outfiles(False)
        self._is_initialized = True
        self._has_run = False
        self._is_finalized = False

    def _run_step(self):
        """Run the next timestep in the specified stage."""
        step_num = self._cmodel.run_next_step()
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
            cscenario = self._scenario._to_cscenario()
            self._cmodel = _ext.CModel(prefix)
            self._cmodel.configure(cscenario)
        self.__results = None
        self._is_open = True

    def close(self):
        """Close and garbage collect the C++ model object."""
        if not self._is_finalized and self._is_initialized:
            self._finalize()
        data = self._cmodel._get_results()
        self.__results = Results(data)
        self.__stepwise = None
        del self._cmodel
        self._cmodel = None
        self._is_open = False
        self._is_finalized = False
        self._is_initialized = False

    def run_sim(self):
        """Run the complete simulation; requires the Simulator to have been opened first."""
        if self._cmodel is None:
            raise RuntimeError("The simulation is not open")
        if self.__stepwise is not None:
            raise RuntimeError("The simulation is already running in stepwise mode")
        self._is_finalized = False
        self._is_initialized = True
        self._cmodel.run_sim()
        self._has_run = True
        self._is_initialized = True
        self._is_finalized = True

    def get_current_state(self) -> "Results":
        """Get the current state of the model as a single-timestep results object."""
        if not self._is_open or self._cmodel is None:
            raise RuntimeError("The simulation has not been started")
        data = self._cmodel._get_current_state()
        results = Results(data)
        return results

    @property
    def results(self) -> "Results":
        """The results for an entire simulation."""
        if not self._has_run and self.__results is None:
            warnings.warn(
                "Simulation incomplete - use get_current_state to get partial results"
            )
            return None
        return self.__results


class Results:
    """
    Python container for the results of a sansmic simulation.

    Parameters
    ----------
    data : :class:`~sansmic.libsansmic.CResults`
        results object from the C++ library
    """

    def __init__(self, data: _ext.CResults) -> None:
        self._data = data
        self.z = pd.Series(data.z_0, name="depth")
        """The depth (:term:`MD`) of the bottom of each cell"""
        self.h = pd.Series(data.h_0, name="height")
        """The height above the original :term:`TD` of the bottom of each cell"""
        self.r_0 = pd.Series(data.r_0, name="radius")
        """The initial radius of each node"""
        self.t = pd.Series(data.t, name="time")
        """The times when data were saved"""
        self.step = pd.Series(data.step, name="step")
        """The steps when the data were saved"""
        # self.summary.set_index(self.t, inplace=True, drop=True)
        self.r = pd.DataFrame(data.r, columns=self.z).T.sort_index(ascending=True)
        """The radius of each node with respect to time"""
        self.dr = pd.DataFrame(data.dr_0, columns=self.z).T.sort_index(ascending=True)
        """The change in radius :term:`wrt` time"""
        self.sg = pd.DataFrame(data.sg, columns=self.z).T.sort_index(ascending=True)
        """The specific gravity of brine in the cell :term:`wrt` time"""
        self.theta = pd.DataFrame(data.theta, columns=self.z).T.sort_index(
            ascending=True
        )
        """The angle of cell wall :term:`wrt` time"""
        self.Q_inj = pd.DataFrame(data.Q_inj, columns=self.z).T.sort_index(
            ascending=True
        )
        """The instantaneous injection rate at each time"""
        self.V = pd.DataFrame(data.V, columns=self.z).T.sort_index(ascending=True)
        """The volume of each cell through time"""
        self.f_dis = pd.DataFrame(data.f_dis, columns=self.z).T.sort_index(
            ascending=True
        )
        """The modified dissolution factor :term:`wrt` time"""
        self.f_flag = pd.DataFrame(data.f_flag, columns=self.z, dtype=int).T.sort_index(
            ascending=True
        )
        """The dissolution status variable indicating if dissolution is occuring."""
        self.xincl = pd.DataFrame(data.xincl, columns=self.z).T.sort_index(
            ascending=True
        )
        """The inclination of the cell wall"""
        self.amd = pd.DataFrame(data.amd, columns=self.z).T.sort_index(ascending=True)
        """Temporary variable for debugging"""
        self.D_coeff = pd.DataFrame(data.D_coeff, columns=self.z).T.sort_index(
            ascending=True
        )
        """The calculated effective diffusion coefficient"""
        self.dC_dz = pd.DataFrame(data.dC_dz, columns=self.z).T.sort_index(
            ascending=True
        )
        """The change in concentration with respect to depth in each cell at each time"""
        self.C_old = pd.DataFrame(data.C_old, columns=self.z).T.sort_index(
            ascending=True
        )
        """The previous concentration in the cell"""
        self.C_new = pd.DataFrame(data.C_new, columns=self.z).T.sort_index(
            ascending=True
        )
        """The current concentration in the cell"""
        self.dC_dt = pd.DataFrame(data.dC, columns=self.z).T.sort_index(ascending=True)
        """The change in concentration with respect to time, related to the rate of dissolution"""
        self.dr_dt = pd.DataFrame(data.dr, columns=self.z).T.sort_index(ascending=True)
        """The change in radius with respect to time"""
        self.C_plm = pd.DataFrame(data.C_plm, columns=self.z).T.sort_index(
            ascending=True
        )
        """The concentration within the injection plume"""
        self.u_plm = pd.DataFrame(data.u_plm, columns=self.z).T.sort_index(
            ascending=True
        )
        """The velocity of the fluid in the injection plume"""
        self.r_plm = pd.DataFrame(data.r_plm, columns=self.z).T.sort_index(
            ascending=True
        )
        """The radius of the injection plume"""
        self.summary = pd.DataFrame.from_dict(
            {
                "t_h": self.t,
                "t_d": self.t / 24.0,
                "step": data.step,
                "stage": data.stage,
                "phase": data.phase,
                "i_inj": data.injCell,
                "i_prod": data.prodCell,
                "i_plume": data.plumCell,
                "i_obi": data.obiCell,
                "err_ode": data.err,
                "z_inj": data.z_inj,
                "z_prod": data.z_prod,
                "z_plume": data.z_plm,
                "z_obi": data.z_obi,
                "z_insol": data.z_insol,
                "h_insol": data.h_insol,
                "l_jet": data.l_jet,
                "u_jet": data.u_jet,
                "r_jet": data.r_jet,
                "V_inj": data.V_injTot,
                "V_fill": data.V_fillTot,
                "V_cav": data.V_cavTot,
                "V_insol": data.V_insolTot,
                "V_vented": data.V_insolVent,
                "Q_out": data.Q_out,
                "sg_out": data.sg_out,
                "sg_ave": data.sg_cavAve,
                "dt_h": data.dt,
            },
        )
        """
        The cavern-wide data summary from the simulation. Contains the following elements:

        .. rubric:: summary DataFrame columns

        t_h
            the time in hours
        t_d
            the time in days
        step
            the step number
        stage
            the stage number
        phase
            the leaching phase, where 1 is active and 0 is passive
        i_inj
            the cell where injection is occurring
        i_prod
            the cell where production is occurring
        i_plume
            the plume stagnation cell index
        i_obi
            the OBI cell index
        err_ode
            the ODE convergence value (closer to 1.0 is better)
        z_inj
            the injection depth
        z_prod
            the production depth
        z_plume
            the plume stagnation depth
        z_obi
            the OBI depth
        z_insol
            the depth to the top of the insoluble pile
        h_insol
            the height of the insolubles above the original TD
        l_jet
            the jet length
        u_jet
            the jet velocity
        r_jet
            the jet radius
        V_inj
            the cumulative injection volume
        V_fill
            the cumulative fill volume
        V_cav
            the total cavern volume
        V_insol
            the volume of insolubles that have fallen
        V_vented
            the volume of insolubles that were sucked out with the brine (ordinary and leach/fill modes only)
        Q_out
            the instantaneous outflow in bbl/d
        sg_out
            the instantaneous outflow specific gravity
        sg_ave
            the average cavern brine concentration
        dt_h
            the timestep size

        """

    def __repr__(self) -> str:
        return "<Results: from {:.3g} to {:.3g} d ({:d} stages)>".format(
            self.t.iloc[0] / 24.0,
            self.t.iloc[-1] / 24.0,
            max(self.summary.stage) - min(self.summary.stage) + 1,
        )


class _OutDataBlock:
    """Useful for processing old SANSMIC .out files"""

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
        self.h_insol = None
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
        self.step = None


class _OutputData(_OutDataBlock):
    """Useful for processing old SANSMIC .out files"""

    def __init__(self):
        super().__init__()
        self.timeTotal = None
        self.r_0 = None
        self.z_0 = None
        self.h_0 = None
        self.V_injTot = None
        self.V_fillTot = None
        self.sg_cavAve = None


if not has_ext:
    _ext.CGeometryFormat = GeometryFormat
    _ext.CRunMode = SimulationMode
    _ext.CModel = Scenario
    _ext.CStage = StageDefinition
