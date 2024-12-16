# coding: utf-8
#
# Copyright (c) 2024 National Technology and Engineering Solutions of
# Sandia, LLC. Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.
#
# SPDX-License-Identifier: BSD-3-Clause.

"""Enums for sansmic configuration options.

Also includes a metaclass (:class:`CaseNormalizedEnumType`) and a
subclass of :class:`enum.IntEnum` that are used to normalize case
of members of the enum to `UPPER_CASE_UNDERSCORED` when used in code
and to `lower-case-dashed` when output as a string value.
"""

from enum import EnumType, IntEnum, Enum
from fractions import Fraction


class CaseNormalizedEnumType(EnumType):
    """Change to a case-normalized name handling."""

    def __getitem__(cls, name):
        """
        Return the member matching `name`.
        """
        return cls._member_map_[name.replace("-", "_").replace(" ", "_").upper()]

    def __setattr__(cls, name, value):
        """
        Block attempts to reassign Enum members.

        A simple assignment to the class namespace only changes one of the
        several possible ways to get an Enum member from the Enum class,
        resulting in an inconsistent Enumeration.
        """
        member_map = cls.__dict__.get("_member_map_", {})
        if name.upper().replace("-", "_").replace(" ", "_") in member_map:
            raise AttributeError("cannot reassign member %r" % (name,))
        super().__setattr__(name, value)


class CaseNormalizedIntEnum(IntEnum, metaclass=CaseNormalizedEnumType):
    """A case-normalized name handling subclass of IntEnum."""

    def __str__(self):
        return "%s" % (self._name_.replace("_", "-").lower())


class CaseNormalizedEnum(Enum, metaclass=CaseNormalizedEnumType):
    """A case-normalized name handling subclass of IntEnum."""

    def __str__(self):
        return "%s" % (self._name_.replace("_", "-").lower())


class Units(CaseNormalizedIntEnum):
    """The units that are used to define the scenario.

    Depths, heights, and cavern radii are in the first (larger) length unit.
    Tubing radii are in the second (smaller) length unit.
    Volumes and volumetric flow rates are in the third unit.

    Durations are in hours. Constant injection rates are in <volume-unit> per day.
    File-based injection rates are in <volume-unit> per hour.

    WARNING: This class has been modified to be case and separator (-, _, and space)
    insensitive. I.e., ``Units['ft-IN bBl']`` is
    normalized and will return ``Units['FT_IN_BBL']``.
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
        return Fraction(3048, 10000) ** 3

    @classmethod
    def barrel(self):
        """1 bbl == 0.158987294928 m³"""
        return 42 * 231 * Fraction(254, 10000) ** 3

    @classmethod
    def centimeter(self):
        """1 cm ≔ 0.01 m"""
        return Fraction(1, 100)


class StopCondition(CaseNormalizedIntEnum):
    """Which stop condition should be used to end the simulation.

    WARNING: This class has been modified to be case and separator (-, _, and space)
    insensitive. I.e., ``StopCondition['dePTH']`` is
    normalized and will return ``StopCondition['DEPTH']``.
    """

    DEPTH = -1
    """Stop when the interface reaches a specified height above the cavern floor"""
    DURATION = 0
    """Stop only based on duration"""
    VOLUME = 1
    """Stop when the total cavern volume reaches a specified value."""

    def __getitem__(cls, name):
        """
        Return the member matching `name`.
        """
        return cls._member_map_[name.replace("-", "_").replace(" ", "_").upper()]


class GeometryFormat(CaseNormalizedIntEnum):
    """The format for the initial cavern geometry.

    WARNING: This class has been modified to be case and separator (-, _, and space)
    insensitive. I.e., ``GeometryFormat['radius LISt']`` is
    normalized and will return ``GeometryFormat['RADIUS_LIST']``.
    """

    RADIUS_LIST = 0
    """The radius is provided from the bottom of the cavern to the top;
    :attr:`~Scenario.num_cells` + 1 values, equally spaced"""
    VOLUME_LIST = 1
    VOLUME_TABLE = -1
    RADIUS_TABLE = 2
    LAYER_CAKE = 5
    """The geometry is provided in a 'layer-cake' style LAS file"""

    def __getitem__(cls, name):
        """
        Return the member matching `name`.
        """
        return cls._member_map_[name.replace("-", "_").replace(" ", "_").upper()]


class SimulationMode(CaseNormalizedIntEnum):
    """The simulation mode determines which options are active for injection.

    WARNING: This class has been modified to be case and separator (-, _, and space)
    insensitive. I.e., ``SimulationMode['leACH-FiLL']`` is
    normalized and will return ``SimulationMode['LEACH_FILL']``.
    """

    ORDINARY = 0
    """Ordinary leaching, with raw water or undersaturated brine injected through
    the inner tubing or outer casing and brine produced from the other."""
    WITHDRAWAL = 1
    """Withdrawal leach, with brine injected through the suspended tubing (hanging
    string) and product produced from top of cavern."""
    LEACH_FILL = 2
    """Simultaneous leaching (water/brine injection and brine production) and
    product injection from the top."""

    def __getitem__(cls, name):
        """
        Return the member matching `name`.
        """
        return cls._member_map_[name.replace("-", "_").replace(" ", "_").upper()]


class RateScheduleType(CaseNormalizedIntEnum):
    """The way that the rate of injection is specified.

    WARNING: This class has been modified to be case and separator (-, _, and space)
    insensitive. I.e., ``RateScheduleType['Constant Rate']`` is
    normalized and will return ``RateScheduleType['CONSTANT_RATE']``.
    """

    CONSTANT_RATE = 1
    """Injection occurs at a constant rate of volume :class:`Units` per day."""
    DATAFILE = 2
    """Injection rate is specified in :class:`Units` per hour in a file."""

    def __getitem__(cls, name):
        """
        Return the member matching `name`.
        """
        return cls._member_map_[name.replace("-", "_").replace(" ", "_").upper()]


class SaveFrequency(CaseNormalizedEnum):
    """The frequency at which data should be saved by the underlying C++ module."""

    HOURLY = 1
    """Save results as close to every hour as machine precision will allow."""
    DAILY = 2
    """Save results as close to every day as machine precision will allow."""
    STAGE = 3
    """Save results at the end of the injection phase and at the end of the stage (if rest-duration > 0)."""
    BYSTAGE = STAGE
    """To-be-deprecated alias"""
