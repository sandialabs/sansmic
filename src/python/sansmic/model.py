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
from dataclasses import InitVar, asdict, dataclass, field
from enum import IntEnum
from fractions import Fraction
from typing import Any, Dict, List, Literal, Union

import numpy as np
import pandas as pd

try:
    import h5py
except ImportError as e:
    h5py = e

from . import libsansmic as _ext

logger = logging.getLogger("sansmic")


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

    @classmethod
    def inch(self):
        """1 in ≔ 0.0254 m"""
        return Fraction(254, 10000)

    @classmethod
    def foot(self):
        """1 ft ≔ 0.3048 m"""
        return Fraction(3048, 10000)

    @classmethod
    def cubic_foot(self):
        """1 ft³ == 0.028316846592 m³"""
        return Fraction(3048**3, 10000**3)

    @classmethod
    def barrel(self):
        """1 bbl == 0.158987294928 m³"""
        return Fraction(42 * 231 * 254, 10000**3)

    @classmethod
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
    """ODE solver absolute error tolerance; CScenario default is 0.01"""
    relative_error: float = None  #  1.0e-4
    """ODE solver absolute error tolerance; CScenario default is 0.0001"""
    coallescing_wells: int = None  #  1
    """The number of wells for multi-well cavern construction; CScenario default is 1"""
    well_separation: float = None  # 0.0
    """The separation distance for multi-well cavern construction; CScenario default is 0.0"""
    jet_model_version: int = None  #  1
    """Jet model version; CScenario default is 1"""
    plume_model_version: int = None  #  1
    """Plume model version; CScenario default is 1"""
    temperature_model_version: int = None  #  0
    """Temperature model version; CScenario default is 0"""
    dissolution_factor: float = None  # 1.0
    """Dissolution factor; CScenario default is 1.0 - this should not be changed unless you are sure you know the effects"""
    max_brine_sg: float = None  #  1.2019
    """Maximum brine specific gravity; CSalt default is 1.2019"""
    solid_density: float = None  #  2.16
    """Rock density in solid form; CSalt default is 2.16 g/cc"""
    entrainment_coeff: float = None  #  0.09
    """Dissolution entrainment coefficient; CScenario default is 0.09"""
    molecular_diffusion: float = None  #  5.03e-5
    """Molecular diffusion coefficient; CScenario default is 5.03e-5"""
    eddy_coefficient: float = None  #  1.142e5
    """Eddy coefficient; CScenario default is 1.142e5"""
    diffusion_beta: float = None  #  0.147
    """Diffusion beta coefficient; default is 0.147"""

    def __setattr__(self, name, value):
        if isinstance(value, str) and value.strip() == "":
            value = None
        super().__setattr__(name, value)

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
    save_frequency: Union[int, Literal["hourly", "daily", "bystage"]] = "daily"
    """The save frequency in number of timesteps or one of {"hourly", "daily", "bystage"}, by default "daily"."""
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
    this value."""
    brine_interface_depth: float = None
    """Set the initial oil-brine interface or blanket level. By default None, which
    will link to the previous stage (as will a value of 0)."""
    # set_cavern_temp: float = None
    # """Set the initial temperature for all brine-filled cells of the cavern; by
    # default None, which will link to the previous stage."""

    product_injection_rate: Union[float, str] = 0.0
    """Either a constant rate of product injection or a file with an injection schedule, by default 0."""
    # product_water_content: float = None
    # """The volume-percent water within the product; a value of None turns this option off."""
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
            k2 = k.strip().replace("-", "_").replace(" ", "_").replace(".", "_")
            if k2 not in self.valid_default_keys:
                logger.warning(  # pragma: no cover
                    "Ignoring non-defaultable or unknown setting {} = {}".format(
                        k, repr(v)
                    )
                )
            elif getattr(self, k2) is None:
                setattr(self, k2, v)

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
        # elif name == "set_initial_conditions" and value is False:
        #     self.set_cavern_sg = None
        #     self.brine_interface_depth = None
        super().__setattr__(name, value)

    @classmethod
    def from_dict(cls, opts: dict, defaults=None) -> "StageDefinition":
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
        return cls(**new_opts, defaults=defaults)

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

        if self.solver_timestep and self.solver_timestep <= 0:
            raise ValueError("Timestep must be greater than 0 hours")
        if self.injection_duration <= 0:
            raise ValueError("Injection duration must be greater than 0 hours")
        if self.rest_duration < 0:
            raise ValueError("Rest duration must be positive")

        # Validate appropriate initial conditions settings
        if self.set_initial_conditions and self.set_cavern_sg is None:
            logger.warning(  # pragma: no cover
                "Setting the initial conditions without setting cavern sg -- cavern sg will be set to fully saturated brine"
            )
            self.set_cavern_sg = 10.0
        # elif (
        #     not self.set_initial_conditions
        #     and self.set_cavern_sg is not None
        #     and self.set_cavern_sg >= 1.0
        # ):
        #     # raise TypeError(
        #     logger.warning( # pragma: no cover
        #         "Setting the starting cavern sg requires set_initial_conditions to be set to True -- setting to 0.0"
        #     )
        #     self.set_cavern_sg = 0.0
        elif not self.set_initial_conditions and self.brine_interface_depth:
            logger.warning(  # pragma: no cover
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

        # Product injection depth and rate validation
        # if (
        #     self.simulation_mode in [SimulationMode.LEACH_FILL, SimulationMode.STORAGE_FILL]
        #     and self.product_injection_depth is None
        # ):
        #     self.product_injection_depth = 0.0
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

    def _to_cstage(self, defaults=None) -> _ext.CStage:
        """Create a CStage object for the C++ interface.
        This method is protected because, in general, it is better not to
        try to operate directly on the C++ stage object.
        """

        self.validate()
        stage = _ext.CStage()
        if defaults is None:
            defaults = dict()
        stage.timestep = (
            self.solver_timestep
            if self.solver_timestep
            else defaults.get("solver_timestep", 0.1)
        )
        stage.title = self.title
        stage.mode = _ext.CRunMode(int(self.simulation_mode))
        if isinstance(self.save_frequency, str):
            if self.save_frequency == "hourly":
                stage.print_interval = int(np.round(1 / stage.timestep))
            elif self.save_frequency == "daily":
                stage.print_interval = int(np.round(24.0 / stage.timestep))
            elif self.save_frequency == "bystage":
                stage.print_interval = int(
                    np.round(
                        (self.injection_duration + self.rest_duration) / stage.timestep
                    )
                )
            else:
                stage.print_interval = int(self.save_frequency)
        elif isinstance(self.save_frequency, (int, float)):
            stage.print_interval = int(self.save_frequency)
        else:
            stage.print_interval = defaults.get(
                "save_frequency", int(np.round(24.0 / stage.timestep))
            )
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
        stage.cavern_sg = self.set_cavern_sg if self.set_cavern_sg is not None else 0.0
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
    units: Units = Units.FT_IN_BBL
    """The units used in describing the scenario."""
    defaults: Dict[str, Union[int, float, str]] = field(default_factory=dict)
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
        elif name == "stages":
            if not hasattr(self, "stages") or value is None:
                value = list()
            elif not hasattr(self, name) or self.stages is None:
                converted_stages = list()
                for ct, stage in enumerate(value):
                    if isinstance(stage, dict):
                        if (
                            "set-initial-conditions" in stage
                            or "set_initial_conditions" in stage
                        ):
                            pass
                        else:
                            stage["set-initial-conditions"] = ct == 0
                        stage = StageDefinition.from_dict(stage, defaults=self.defaults)
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
            stg_dict = stage.to_dict(keep_empty)
            for k in ret["defaults"].keys():
                if (
                    k in stg_dict
                    and k in ret["defaults"]
                    and stg_dict[k] == ret["defaults"][k]
                ):
                    del stg_dict[k]
            ret["stages"].append(stg_dict)
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

    def new_simulation(
        self,
        prefix="temp",
        verbosity=0,
        generate_tst_file=True,
        generate_out_file=False,
    ) -> "Simulator":
        """Create a new :class:`Simulator` object.

        Parameters
        ----------
        prefix : str
            The prefix to use when creating output files, by default temp.
        verbosity : int
            A verbosity level to pass to the C++ model, by default 0.
        generate_tst_file : bool, optional
            Generate an old-style .TST file, by default True
        generate_out_file : bool, optional
            Generate an old-style .OUT file, by default False
        """

        return Simulator(self, prefix, verbosity, generate_tst_file, generate_out_file)

    def _to_cscenario(self):
        """Create a C++ model object; in general, this should only be called internally."""

        cscenario = _ext.CScenario()
        cscenario.fraction_insolubles = self.insolubles_ratio
        cscenario.floor_depth = self.floor_depth
        cscenario.geometry_format = _ext.CGeometryFormat(int(self.geometry_format))
        cscenario.num_cells = self.num_cells
        cscenario.ullage_standoff = self.ullage_standoff
        cscenario.cavern_height = self.cavern_height

        # Use defaults for all advanced settings unless they have been changed
        if self.advanced.coallescing_wells is not None:
            cscenario.coallescing_wells = self.advanced.coallescing_wells

        if self.advanced.well_separation is not None:
            cscenario.well_separation = self.advanced.well_separation

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

        # TODO: This needs to be far more robust
        if self.geometry_format is GeometryFormat.RADIUS_LIST and isinstance(
            self.geometry_data, str
        ):
            radii = list([0])  # Have to pad the data by 1
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
            cstage = stage._to_cstage(defaults=self.defaults)
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
    vebosity : int, optional
        The verbosity level for console output, by default 0
    generate_tst_file : bool, optional
        Generate an old-style .TST file, by default True
    generate_out_file : bool, optional
        Generate an old-style .OUT file, by default False
    """

    def __init__(
        self,
        scenario: Union[Scenario, _ext.CModel],
        prefix="temp",
        verbosity=0,
        generate_tst_file=True,
        generate_out_file=False,
    ):
        self._scenario = None
        self._prefix = prefix
        self._cmodel = None
        self._verbosity = verbosity
        self._b_use_tstfile = generate_tst_file
        self._b_use_outfile = generate_out_file
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
        self._cmodel.set_verbosity_level(self._verbosity)
        self._cmodel.generate_tst_file(self._b_use_tstfile)
        self._cmodel.generate_out_file(self._b_use_outfile)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @property
    def steps(self) -> StepwiseIterator:
        """Provides an iterator that will run each step of each stage in turn."""

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
        self.df_t_1D: pd.DataFrame = None
        """Results that are 1D indexed in the D(time) domain."""
        self.df_z_1D: pd.DataFrame = None
        """Results that are 1D indexed in the D(space) domain."""
        self.df_t_z_2D: pd.DataFrame = None
        """Results that are 2D indexed in the D(time, space) domain."""

        if data is None:
            return
        t = pd.Series(data.t, name="time")
        r = pd.DataFrame(data.r)
        dr = pd.DataFrame(data.dr_0)
        sg = pd.DataFrame(data.sg)
        theta = pd.DataFrame(data.theta)
        Q_inj = pd.DataFrame(data.Q_inj)
        V = pd.DataFrame(data.V)
        f_dis = pd.DataFrame(data.f_dis)
        f_flag = pd.DataFrame(data.f_flag, dtype=int)
        xincl = pd.DataFrame(data.xincl)
        amd = pd.DataFrame(data.amd)
        D_coeff = pd.DataFrame(data.D_coeff)
        dC_dz = pd.DataFrame(data.dC_dz)
        C_old = pd.DataFrame(data.C_old)
        C_new = pd.DataFrame(data.C_new)
        dC_dt = pd.DataFrame(data.dC)
        dr_dt = pd.DataFrame(data.dr)
        C_plm = pd.DataFrame(data.C_plm)
        u_plm = pd.DataFrame(data.u_plm)
        r_plm = pd.DataFrame(data.r_plm)

        self.df_z_1D = pd.DataFrame.from_dict(
            dict(
                z=data.z_0,
                h=data.h_0,
                r=data.r_0,
            )
        )
        self.df_z_1D.index.set_names(["node_idx"], inplace=True)

        self.df_t_1D = pd.DataFrame.from_dict(
            dict(
                t_h=t,
                t_d=t / 24.0,
                step=data.step,
                stage=data.stage,
                phase=data.phase,
                i_inj=data.injCell,
                i_prod=data.prodCell,
                i_plume=data.plumCell,
                i_obi=data.obiCell,
                err_ode=data.err,
                z_inj=data.z_inj,
                z_prod=data.z_prod,
                z_plume=data.z_plm,
                z_obi=data.z_obi,
                z_insol=data.z_insol,
                h_insol=data.h_insol,
                l_jet=data.l_jet,
                u_jet=data.u_jet,
                r_jet=data.r_jet,
                V_inj=data.V_injTot,
                V_fill=data.V_fillTot,
                V_cav=data.V_cavTot,
                V_insol=data.V_insolTot,
                V_vented=data.V_insolVent,
                Q_out=data.Q_out,
                sg_out=data.sg_out,
                sg_ave=data.sg_cavAve,
                dt_h=data.dt,
            ),
        )
        self.df_t_1D.index.set_names(["save_idx"], inplace=True)

        self.df_t_z_2D = pd.DataFrame.from_dict(
            dict(
                r=r.stack(),
                dr=dr.stack(),
                dr_dt=dr_dt.stack(),
                theta=theta.stack(),
                xincl=xincl.stack(),
                V=V.stack(),
                sg=sg.stack(),
                dC_dt=dC_dt.stack(),
                dC_dz=dC_dz.stack(),
                f_flag=f_flag.stack(),
                f_dis=f_dis.stack(),
                D_coeff=D_coeff.stack(),
                r_plm=r_plm.stack(),
                u_plm=u_plm.stack(),
                C_plm=C_plm.stack(),
                Q_inj=Q_inj.stack(),
            )
        )
        self.df_t_z_2D.index.set_names(["save_idx", "node_idx"], inplace=True)

    def __repr__(self) -> str:
        return "<Results: from {:.3g} to {:.3g} d ({:d} stages)>".format(
            self.df_t_1D.t_d.iloc[0],
            self.df_t_1D.t_d.iloc[-1],
            max(self.df_t_1D.stage) - min(self.df_t_1D.stage) + 1,
        )

    @property
    def time(self) -> pd.Series:
        """Time domain: time in days, by timestep"""
        return self.df_t_1D.t_d

    @property
    def injection_cell(self) -> pd.Series:
        """The injection cell, by timestep"""
        return self.df_t_1D.i_inj

    @property
    def production_cell(self) -> pd.Series:
        """The production cell, by timestep"""
        return self.df_t_1D.i_prod

    @property
    def plume_stagnation_cell(self) -> pd.Series:
        """The plume stagnation cell, by timestep"""
        return self.df_t_1D.i_plume

    @property
    def interface_cell(self) -> pd.Series:
        """The brine interface cell, by timestep"""
        return self.df_t_1D.i_obi

    @property
    def ode_convergence(self) -> pd.Series:
        """The ODE convergence factor, by timestep"""
        return self.df_t_1D.err_ode

    @property
    def injection_depth(self) -> pd.Series:
        """The water/brine injection depth, by timestep"""
        return self.df_t_1D.z_inj

    @property
    def production_depth(self) -> pd.Series:
        """The depth of brine production, by timestep"""
        return self.df_t_1D.z_prod

    @property
    def plume_stagnation_depth(self) -> pd.Series:
        """The plume stagnation depth, by timestep"""
        return self.df_t_1D.z_plm

    @property
    def interface_depth(self) -> pd.Series:
        """The brine interface depth, by timestep"""
        return self.df_t_1D.z_obi

    @property
    def cavern_total_depth(self) -> pd.Series:
        """The depth to the new cavern TD, by timestep"""
        return self.df_t_1D.z_insol

    @property
    def insoluble_material_height(self) -> pd.Series:
        """The height of the insoluble materials deposited, by timestep"""
        return self.df_t_1D.h_insol

    @property
    def jet_length(self) -> pd.Series:
        """The length of the injection fluid jet's penetration below EOT, by timestep"""
        return self.df_t_1D.l_jet

    @property
    def jet_velocity(self) -> pd.Series:
        """The velocity of the injection fluid, by timestep"""
        return self.df_t_1D.u_jet

    @property
    def jet_radius(self) -> pd.Series:
        """The radius of the injection fluid jet and starting plume radius, by timestep"""
        return self.df_t_1D.r_jet

    @property
    def cumulative_injection_volume(self) -> pd.Series:
        """Total injected water/brine volume, by timestep"""
        return self.df_t_1D.V_inj

    @property
    def cumulative_fill_volume(self) -> pd.Series:
        """Total product fill volume, by timestep"""
        return self.df_t_1D.V_fill

    @property
    def cavern_volume(self) -> pd.Series:
        """Total cavern volume, by timestep"""
        return self.df_t_1D.V_cav

    @property
    def cumulative_insoluble_material_volume(self) -> pd.Series:
        """Total volume of insoluble materials deposited, by timestep"""
        return self.df_t_1D.V_insol

    @property
    def cumulative_insoluble_material_vented(self) -> pd.Series:
        """Total volume of insoluble materials vented, by timestep"""
        return self.df_t_1D.V_vented

    @property
    def brine_production_rate(self) -> pd.Series:
        """Instantaneous brine production rate, by timestep"""
        return self.df_t_1D.Q_out

    @property
    def brine_production_sg(self) -> pd.Series:
        """Instantaneous brine production cell sg, by timestep"""
        return self.df_t_1D.sg_out

    @property
    def cavern_average_brine_sg(self) -> pd.Series:
        """Instantaneous average cavern brine sg, by timestep"""
        return self.df_t_1D.sg_ave

    @property
    def depths(self) -> pd.Series:
        """Vertical domain: depths, by node"""
        return self.df_z_1D.z

    @property
    def node_heights(self) -> pd.Series:
        """Vertical domain: heights above initial TD, by node"""
        return self.df_z_1D.h

    @property
    def hours(self) -> pd.Series:
        """Time domain: time in hours, by timestep"""
        return self.df_t_1D.t_h

    @property
    def step_size(self) -> pd.Series:
        """Simulation step size, by timestep"""
        return self.df_t_1D.dt

    @property
    def step_number(self) -> pd.Series:
        """Simulation step number, by timestep"""
        return self.df_t_1D.step

    @property
    def stage_number(self) -> pd.Series:
        """Simulation stage number, by timestep"""
        return self.df_t_1D.stage

    @property
    def injection_phase(self) -> pd.Series:
        """Simulation injection(1) or static(0) phase, by timestep"""
        return self.df_t_1D.phase

    @property
    def radius(self) -> pd.DataFrame:
        """Cavern radius, by node and timestep"""
        return self.df_t_z_2D.r.unstack().T

    @property
    def change_in_radius(self) -> pd.DataFrame:
        """Change in cavern radius since previous time, by node and timestep"""
        return self.df_t_z_2D.dr.unstack().T

    @property
    def rate_of_change_in_radius(self) -> pd.DataFrame:
        """Rate of change in cavern radius, by node and timestep"""
        return self.df_t_z_2D.dr_dt.unstack().T

    @property
    def wall_angle(self) -> pd.DataFrame:
        """Dissolution wall angle, by node and timestep"""
        return self.df_t_z_2D.theta.unstack().T

    @property
    def wall_factor(self) -> pd.DataFrame:
        """Dissolution wall factor, by node and timestep"""
        return self.df_t_z_2D.xincl.unstack().T

    @property
    def cell_volume(self) -> pd.DataFrame:
        """Cell volume above node, by node and timestep"""
        return self.df_t_z_2D.V.unstack().T

    @property
    def plume_radius(self) -> pd.DataFrame:
        """Plume radius at node depth, by node and timestep"""
        return self.df_t_z_2D.r_plm.unstack().T

    @property
    def plume_velocity(self) -> pd.DataFrame:
        """Plume velocity at node depth, by node and timestep"""
        return self.df_t_z_2D.u_plm.unstack().T

    @property
    def plume_sg(self) -> pd.DataFrame:
        """Plume concentration(sg) at node depth, by node and timestep"""
        return self.df_t_z_2D.C_plm.unstack().T

    @property
    def cell_injection_rate(self) -> pd.DataFrame:
        """Injection rate within cell above node, by node and timestep"""
        return self.df_t_z_2D.Q_inj.unstack().T

    @property
    def cell_sg(self) -> pd.DataFrame:
        """Concentration(sg) in cell above node, by node and timestep"""
        return self.df_t_z_2D.sg.unstack().T

    @property
    def rate_of_change_in_sg(self) -> pd.DataFrame:
        """Rate of change in concentration in cell above node, by node and timestep"""
        return self.df_t_z_2D.dC_dt.unstack().T

    @property
    def vertical_diffusion_rate(self) -> pd.DataFrame:
        """Vertical change in concentration across node, by node and timestep"""
        return self.df_t_z_2D.dC_dz.unstack().T

    @property
    def state_indicator(self) -> pd.DataFrame:
        """Dissolution state flag, by node and timestep"""
        return self.df_t_z_2D.f_flag.unstack().T

    @property
    def effective_dissolution_factor(self) -> pd.DataFrame:
        """Effective dissolution factor at node, by node and timestep"""
        return self.df_t_z_2D.f_dis.unstack().T

    @property
    def effective_diffusion_coefficient(self) -> pd.DataFrame:
        """Effective diffusion coefficient at node, by node and timestep"""
        return self.df_t_z_2D.D_coeff.unstack().T

    def to_dict(self) -> dict:
        """Convert the results object to a dictionary format."""
        ret = dict()
        for k in ["df_t_1D", "df_z_1D", "df_t_z_2D"]:
            ret[k] = getattr(self, k).to_dict("tight")
        return ret

    @classmethod
    def from_dict(cls, d) -> "Results":
        """Create a new object from a tight dictionary representation."""
        new = cls(None)
        for k in ["df_t_1D", "df_z_1D", "df_t_z_2D"]:
            setattr(new, k, pd.DataFrame.from_dict(d[k], orient="tight"))
        return new

    def to_hdf(
        self,
        filename: str,
        *,
        compression: str = "gzip",
        compression_opts: int = 9,
        shuffle: bool = True,
        fletcher32: bool = True,
        **kwargs,
    ):
        """Write results to an HDF5 file.

        Parameters
        ----------
        filename : str
            A filename with extension ".h5", or filename prefix without extension, to write to

        Keyword Arguments
        -----------------
        compression : str
            The HDF5 compression filter library to use, by default 'gzip'
        compression_opts : int
            The compression level to use, by default 9 (highest)
        shuffle : bool
            Whether to use the shuffle filter, by default True
        fletcher32 : bool
            Whether to use checksum filter, by default True
        kwargs : additional keyword arguments
            Any additional arguments to be passed to the :mod:`h5py` module ``create_dataset`` call.

        """

        if isinstance(h5py, ImportError):
            raise RuntimeError("Optional dependency not installed: h5py") from h5py
        if not filename.lower().endswith(".h5"):
            filename = filename + ".h5"
        with h5py.File(filename, "w") as f:
            f.create_dataset(
                "df_t_1D",
                data=self.df_t_1D.to_records(),
                dtype=[(k, v) for (k, v) in self.df_t_1D.dtypes.items()],
                compression=compression,
                compression_opts=compression_opts,
                shuffle=shuffle,
                fletcher32=fletcher32,
                **kwargs,
            )
            f.create_dataset(
                "df_z_1D",
                data=self.df_z_1D.to_records(),
                dtype=[(k, v) for (k, v) in self.df_z_1D.dtypes.items()],
                compression=compression,
                compression_opts=compression_opts,
                shuffle=shuffle,
                fletcher32=fletcher32,
                **kwargs,
            )
            f.create_dataset(
                "df_t_z_2D",
                data=self.df_t_z_2D.to_records(),
                dtype=[(k, v) for (k, v) in self.df_t_z_2D.dtypes.items()],
                compression=compression,
                compression_opts=compression_opts,
                shuffle=shuffle,
                fletcher32=fletcher32,
                **kwargs,
            )

    @classmethod
    def from_hdf(cls, filename: str) -> "Results":
        """Read results from an HDF5 file.

        Parameters
        ----------
        filename : str
            A filename with extension ".h5", or filename prefix without extension, to read from.

        """

        results = cls(None)
        if isinstance(h5py, ImportError):
            raise RuntimeError("Optional dependency not installed: h5py") from h5py
        if not filename.lower().endswith(".h5"):
            filename = filename + ".h5"
        with h5py.File(filename, "r") as f:
            results.df_t_1D = pd.DataFrame(
                np.array(f["df_t_1D"], dtype=f["df_t_1D"].dtype).view(np.recarray)
            )
            results.df_z_1D = pd.DataFrame(
                np.array(f["df_z_1D"], dtype=f["df_z_1D"].dtype).view(np.recarray)
            )
            results.df_t_z_2D = pd.DataFrame(
                np.array(f["df_t_z_2D"], dtype=f["df_t_z_2D"].dtype).view(np.recarray)
            )
        return results


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
