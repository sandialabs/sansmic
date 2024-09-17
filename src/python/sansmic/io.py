# coding: utf-8
#
# Copyright (c) 2024 National Technology and Engineering Solutions of
# Sandia, LLC. Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.
#
# SPDX-License-Identifier: BSD-3-Clause.

"""Provide input and output functions for sansmic."""

import json
import logging
import os
import sys
import warnings
from dataclasses import fields
from enum import IntEnum

if sys.version_info[1] < 11:
    import tomli as tomllib
else:
    import tomllib

import numpy as np
import pandas as pd

try:
    import yaml
except ImportError as e:
    try:
        import ruamel.yaml as _yaml

        yaml = _yaml.YAML
    except ImportError:
        yaml = e

try:
    import h5py
except ImportError as e:
    h5py = e

from .model import (
    GeometryFormat,
    Results,
    Scenario,
    SimulationMode,
    StageDefinition,
    StopCondition,
    _OutDataBlock,
    _OutputData,
)

logger = logging.getLogger("sansmic")


def read_scenario(
    config_file: str, warn=True, strict=False, *, format=None
) -> Scenario:
    """Load a sansmic scenario file.

    Parameters
    ----------
    config_file : str
        filename or full path to a new-style configuration file
    warn : bool, optional
        emit warnings if there are extra, unrecognized options in the file, by default True
    strict : bool, optional
        raise errors if there are extra, unrecognized options in the file, by default False

    Returns
    -------
    Scenario
        the sansmic scenario

    Raises
    ------
    RuntimeError
        if there is an unrecognized option and ``strict`` is True, or if there are no
        stages defined.
    """
    logger.debug(
        "Attempting to read scenario {} in {} mode".format(
            config_file, "strict" if strict else "permissive"
        )
    )
    try:
        if config_file.lower().endswith(".json") or format == "json":
            with open(config_file, "r") as fin:
                data = json.load(fin)
        elif config_file.lower().endswith(".yaml") or format == "yaml":
            with open(config_file, "r") as fin:
                data = yaml.safe_load(fin)
        elif config_file.lower().endswith(".toml") or format == "toml":
            with open(config_file, "rb") as fin:
                data = tomllib.load(fin)
        elif config_file.lower().endswith(".dat") or format == "dat":
            with open(config_file, "r") as fin:
                data = read_dat(fin)
            return data
        else:
            with open(config_file, "rb") as fin:
                data = tomllib.load(fin)
    except:
        logger.error("Error reading scenario file {}".format(config_file))
        raise

    if "stages" not in data or len(data["stages"]) < 1:
        logger.error("No stages provided. Failed to load valid scenario.")
        raise RuntimeError("No stages provided.")
    return Scenario.from_dict(data)


def read_dat(str_or_buffer, *, ignore_errors=False) -> Scenario:
    """Read an old-style SANSMIC input DAT file.

    Parameters
    ----------
    str_or_buffer : str or stream
        The file or stream buffer to process.

    Returns
    -------
    Scenario
        the scenario derived from the DAT file

    """
    scenario = Scenario()
    with open(str_or_buffer, "r") if isinstance(
        str_or_buffer, str
    ) else str_or_buffer as fin:
        prefix = fin.name
        scenario.title = "Converted from {}".format(prefix)
        first = True
        logger.debug("Reading old DAT format file: {}".format(prefix))
        while True:
            stage = StageDefinition()
            # Stage - block 1
            title = fin.readline().strip()
            if title == "END":
                break
            logger.debug(" Record 1: Stage title")
            logger.debug("  title         = {}".format(title))
            # Stage - block 2
            (
                num_cells,
                mode,
                print_interval,
                subsequent,
                iResetGeo,
                rest_duration,
                num_coallescing,
                geometry_format,
                stop_value,
            ) = (
                fin.readline().strip().split()
            )
            logger.debug(" Record 2: General")
            logger.debug("  ndiv          = {}".format(num_cells))
            logger.debug("  leachtype     = {}".format(mode))
            logger.debug("  iprint        = {}".format(print_interval))
            logger.debug("  repeat        = {}".format(subsequent))
            logger.debug("  resetgeo[depr]= {}".format(iResetGeo))
            logger.debug("  iWait         = {}".format(rest_duration))
            logger.debug("  nco           = {}".format(num_coallescing))
            logger.debug("  idata         = {}".format(geometry_format))
            logger.debug("  ivol          = {}".format(stop_value))
            # Stage - block 3
            (
                cavern_height,
                injection_height,
                production_height,
                interface_height,
                ullage_standoff,
            ) = (
                fin.readline().strip().split()
            )
            logger.debug(" Record 3: Heights")
            logger.debug("  zmax          = {}".format(cavern_height))
            logger.debug("  zi            = {}".format(injection_height))
            logger.debug("  zp            = {}".format(production_height))
            logger.debug("  zb            = {}".format(interface_height))
            logger.debug("  zu            = {}".format(ullage_standoff))
            # Stage - block 4
            injection_rate = fin.readline().strip()
            logger.debug(" Record 4: Injection flow rates")
            logger.debug("  QI            = {}".format(injection_rate))
            # TODO: FIXME: handle injection tables? Or do we ignore for old
            # style inputs?
            #
            # Stage - block 5
            (
                inn_tbg_inside_radius,
                inn_tbg_outside_radius,
                out_csg_inside_radius,
                out_csg_outside_radius,
            ) = (
                fin.readline().strip().split()
            )
            logger.debug(" Record 5: Casing and tubing")
            logger.debug("  rpi           = {}".format(inn_tbg_inside_radius))
            logger.debug("  rpo           = {}".format(inn_tbg_outside_radius))
            logger.debug("  rcasi         = {}".format(out_csg_inside_radius))
            logger.debug("  rcaso         = {}".format(out_csg_outside_radius))
            # Stage - block 6
            injection_fluid_sg, cavern_sg = fin.readline().strip().split()
            logger.debug(" Record 6: Water and brine")
            logger.debug("  sginj         = {}".format(injection_fluid_sg))
            logger.debug("  sgcav         = {}".format(cavern_sg))
            # Stage - block 7
            timestep, injection_duration = fin.readline().strip().split()
            logger.debug(" Record 7: Timing")
            logger.debug("  dt            = {}".format(timestep))
            logger.debug("  tend          = {}".format(injection_duration))
            # Stage - block 8
            fill_rate, tDelay, coallescing_well_separation = (
                fin.readline().strip().split()
            )
            logger.debug(" Record 8: Oil fill rates")
            logger.debug("  qfil          = {}".format(fill_rate))
            logger.debug("  tdlay  [depr.]= {}".format(tDelay))
            logger.debug("  sep           = {}".format(coallescing_well_separation))
            if float(tDelay) != 0:
                logger.warning("The TDLAY option should not be used. Ignoring.")
            # First stage - block 9
            if first:
                logger.debug("First stage initialization")
                logger.debug(" Record 9: Geometry")
                geometry_data = list()
                first = False
                geometry_format = GeometryFormat(int(geometry_format))
                if geometry_format == GeometryFormat.RADIUS_LIST:
                    radii = list()
                    for _ in range(int(num_cells) + 1):
                        r = float(fin.readline().strip())
                        radii.append(r)
                    geometry_data = dict(radii=radii)
                elif geometry_format == GeometryFormat.VOLUME_LIST:
                    ndata = int(fin.readline().strip())
                    depths = list()
                    volumes = list()
                    for _ in range(ndata):
                        z = float(fin.readline().strip())
                        depths.append(z)
                    for _ in range(ndata):
                        v = float(fin.readline().strip())
                        volumes.append(v)
                    geometry_data = dict(depths=depths, volumes=volumes)
                elif geometry_format == GeometryFormat.VOLUME_TABLE:
                    ndata = int(fin.readline().strip())
                    depths = list()
                    volumes = list()
                    for _ in range(ndata):
                        z, v = fin.readline().strip().split()
                        depths.append(float(z))
                        volumes.append(float(v))
                    geometry_data = dict(depths=depths, volumes=volumes)
                elif geometry_format == GeometryFormat.VOLUME_TABLE:
                    ndata = int(fin.readline().strip())
                    depths = list()
                    radii = list()
                    for _ in range(ndata):
                        z, v = fin.readline().strip().split()
                        depths.append(float(z))
                        radii.append(float(v))
                    geometry_data = dict(depths=depths, radii=radii)

                dissolution_factor, insoluble_fraction, refdep, depth = (
                    fin.readline().strip().split()
                )
                logger.debug(" Record 10: Miscellaneous")
                logger.debug("  zdis          = {}".format(dissolution_factor))
                logger.debug("  zfin          = {}".format(insoluble_fraction))
                logger.debug("  refdep [depr.]= {}".format(refdep))
                logger.debug("  depth         = {}".format(depth))
                if refdep != depth:
                    logger.warning(
                        "The REFDEP is no longer used, only DEPTH. Ignoring REFDEP."
                    )
                if float(dissolution_factor) != 1.0:
                    logger.warning(
                        "The ZDIS should always be 1.0. This is a dangerous choice."
                    )
                scenario.geometry_format = geometry_format
                scenario.geometry_data = geometry_data
                scenario.num_cells = int(num_cells)
                scenario.coallescing_wells = int(num_coallescing)
                scenario.well_separation = float(coallescing_well_separation)
                df = float(dissolution_factor)
                if df != 1.0:
                    scenario.advanced.dissolution_factor = df
                scenario.insolubles_ratio = float(insoluble_fraction)
                scenario.floor_depth = float(depth)
                scenario.cavern_height = float(cavern_height)
                scenario.ullage_standoff = float(ullage_standoff)
            else:
                logger.debug("Subsequent stage initialization")
            if not first and scenario.num_cells != int(num_cells) and not ignore_errors:
                logger.critical("Old-style DAT has changes in NDIV throughout the file")
                raise TypeError("Invalid data in DAT file: NDIV not constant.")
            elif not first and int(iResetGeo) != 0 and not ignore_errors:
                logger.critical(
                    "The RESETGEO flag (item 5 on line 2 of stage {}) is set to 1. This is not valid in this version of SANSMIC. Please correct this and rerun.".format(
                        len(scenario.stages) + 1
                    )
                )
                raise TypeError("Invalid data in DAT file: RESETGEO not 0")
            if not first and isinstance(cavern_sg, (float, int)) and cavern_sg > 1.0:
                logger.warning(
                    "The REPEAT option was supposed to turn off the cavern SG; it did not do so. sansmic is currently mimicing SANSMIC behavior and resetting the cavern brine to {} sg in stage {}. \n\nIf this is not what is intended, please manually remove the 'set-cavern-sg' entry from stages after stage 1. This behavior will change in future releases.".format(cavern_sg, len(scenario.stages)+1)
                )
            stage.title = title
            stage.simulation_mode = SimulationMode(int(mode))
            stage.save_frequency = int(print_interval)
            stage.set_initial_conditions = bool(1 - int(subsequent))
            stage.rest_duration = float(rest_duration)
            stage.stop_value = float(stop_value)
            if stage.stop_value < 0:
                stage.stop_condition = StopCondition.DEPTH
            if stage.stop_value > 0:
                stage.stop_condition = StopCondition.VOLUME
            stage.brine_injection_depth = float(injection_height)
            stage.brine_production_depth = float(production_height)
            stage.brine_interface_depth = float(interface_height)
            stage.brine_injection_rate = float(injection_rate)
            stage.inner_tbg_inside_diam = float(inn_tbg_inside_radius) * 2
            stage.inner_tbg_outside_diam = float(inn_tbg_outside_radius) * 2
            stage.outer_csg_inside_diam = float(out_csg_inside_radius) * 2
            stage.outer_csg_outside_diam = float(out_csg_outside_radius) * 2
            stage.brine_injection_sg = float(injection_fluid_sg)
            stage.set_cavern_sg = float(cavern_sg)
            stage.solver_timestep = float(timestep)
            stage.injection_duration = float(injection_duration)
            stage.product_injection_rate = float(fill_rate)
            scenario.stages.append(stage)
            logger.debug("Finished reading stage")
    return scenario


def write_scenario(scenario: Scenario, filename: str, *, redundant=False, format=None):
    """Write a new-style SANSMIC scenario file (preferred extension is .toml)"""
    sdict = scenario.to_dict(redundant)
    keys = [k for k in sdict.keys()]
    for k in keys:
        if not redundant and k == "coallescing-wells" and sdict[k] == 1:
            sdict[k] = None
            sdict["well-separation"] = None
        if sdict[k] is None:
            if redundant:
                sdict[k] = ""
            else:
                del sdict[k]
        elif isinstance(sdict[k], IntEnum):
            sdict[k] = sdict[k].name.lower().replace("_", "-")
    for ct, s in enumerate(sdict["stages"]):
        if (
            not redundant
            and scenario.stages[ct].stop_condition == StopCondition.DURATION
        ):
            del s["stop-condition"]
            del s["stop-value"]
        keys = [k for k in s.keys()]
        for k in keys:
            if (
                not redundant
                and k in scenario.defaults
                and scenario.defaults.get(k, None) == s[k]
            ):
                del s[k]
                continue
            if s[k] is None:
                if redundant:
                    s[k] = ""
                else:
                    del s[k]
            elif isinstance(s[k], IntEnum):
                s[k] = s[k].name.lower().replace("_", "-")
        if not redundant and ct == 0 and "set-initial-conditions" in s:
            del s["set-initial-conditions"]
    name, ext = os.path.splitext(filename)
    with open(filename, "w") as fout:
        if ext.lower() == ".toml" or format == "toml":
            for k, v in sdict.items():
                if k in ["stages", "defaults", "advanced"]:
                    continue
                elif k == "geometry-data" and isinstance(v, dict):
                    for k2, v2 in v.items():
                        fout.write("geometry-data.{} = {}\n".format(k2, v2))
                    continue
                if isinstance(v, bool):
                    v = str(v).lower()
                elif isinstance(v, str):
                    v = repr(v)
                elif isinstance(v, IntEnum):
                    v = repr(v.name.lower().replace("_", "-"))
                fout.write("{} = {}\n".format(k, v))

            if "advanced" in sdict:
                fout.write("\n[advanced]\n")
                for k, v in sdict["advanced"].items():
                    if isinstance(v, bool):
                        v = str(v).lower()
                    elif isinstance(v, str):
                        v = repr(v)
                    elif isinstance(v, IntEnum):
                        v = repr(v.name.lower().replace("_", "-"))
                    fout.write("{} = {}\n".format(k, v))

            if len(sdict["defaults"]) > 0:
                fout.write("\n[defaults]\n")
                for k, v in sdict["defaults"].items():
                    if isinstance(v, bool):
                        v = str(v).lower()
                    elif isinstance(v, str):
                        v = repr(v)
                    elif isinstance(v, IntEnum):
                        v = repr(v.name.lower().replace("_", "-"))
                    fout.write("{} = {}\n".format(k, v))

            for s in sdict["stages"]:
                fout.write("\n[[stages]]\n")
                for k, v in s.items():
                    if isinstance(v, bool):
                        v = str(v).lower()
                    elif isinstance(v, str):
                        v = repr(v)
                    elif isinstance(v, IntEnum):
                        v = repr(v.name.lower().replace("_", "-"))
                    fout.write("{} = {}\n".format(k, v))
        elif ext.lower() == ".json" or format == "json":
            json.dump(sdict, fout)
        elif ext.lower() == ".yaml" or format == "yaml":
            yaml.dump(sdict, fout)
        else:
            logger.critical(
                'Unknown file extension ({}) and invalid format ({}). Valid extensions/formats are "toml", "yaml" and "json".'.format(
                    ext, format
                )
            )
            raise RuntimeError("Unknown file format for scenario output")


def read_classic_out_ddl(file_prefix):
    """
    Read in the data from ".out" and ".ddl" files from a classic SANSMIC run.

    Parameters
    ----------
    file_prefix : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    started = False
    dataline = False
    skip1 = False
    steps = list()
    step = None
    phaseline = False

    r_0 = list()
    z_0 = list()
    h_0 = list()

    with open(file_prefix + ".ddl", "r") as fddl:
        for ct, line in enumerate(fddl.readlines()):
            if ct == 0:
                continue
            words = line.split()
            idx = int(words[0])
            z = float(words[1])
            r = float(words[2])
            r_0.insert(0, r)
            z_0.insert(0, z)
            if idx == 1:
                break
    for z in z_0:
        h = z - z_0[0]
        h_0.append(h)

    with open(file_prefix + ".out", "r") as fout:
        for line in fout.readlines():
            words = line.split()
            if skip1:
                skip1 = False
                continue
            if not skip1 and len(words) == 0:
                continue
            if not started and words[0] == "TIME=":
                started = True
                step = _OutDataBlock()
            elif not started:
                continue
            if words[0] == "TIME=":
                step.t = float(words[1]) * 24.0
                step.dt = float(words[7]) * 24.0
                step.t_prime = float(words[11]) * 24.0
                phaseline = True
            elif phaseline and words[0] == "INJECTION":
                step.phase = 1
                step.stage = int(words[2])
                phaseline = False
            elif phaseline and words[0] == "WORKOVER":
                step.phase = 0
                step.stage = int(words[2])
                phaseline = False
            elif words[0] == "I(INJ),":
                step.injCell = int(words[4])
                step.prodCell = int(words[5])
                step.obiCell = int(words[6])
                step.plumCell = int(words[7])
            elif words[0] == "H(INJ),":
                step.z_inj = z_0[0] - float(words[4])
                step.z_prod = z_0[0] - float(words[5])
                # step.z_obi = z_0[0] - float(words[6])
                step.z_plm = z_0[0] - float(words[7])
            elif words[0] == "Ljet(ft),":
                step.l_jet = float(words[4])
                step.r_jet = float(words[5])
                step.u_jet = float(words[6])
                skip1 = True
            elif words[0] == "HEIGHT":
                skip1 = True
                dataline = True
            elif dataline:
                cellid = int(line[0:5])
                if cellid == 1:
                    dataline = False
                    skip1 = True
                step.idx.insert(0, cellid)
                step.h.insert(0, float(line[5:16]))
                step.r.insert(0, float(line[16:27]))
                step.dr_0.insert(0, float(line[27:38]))
                step.sg.insert(0, float(line[38:49]))
                step.theta.insert(0, float(line[49:60]))
                step.Q_inj.insert(0, float(line[60:71]))
                step.V.insert(0, float(line[71:82]))
                step.f_dis.insert(0, float(line[82:93]))
                step.f_flag.insert(0, int(line[93:98]))
                step.xincl.insert(0, float(line[98:109]))
                step.amd.insert(0, float(line[109:120]))
                step.D_coeff.insert(0, float(line[120:131]))
                step.dC_dz.insert(0, float(line[131:142]))
                step.C_old.insert(0, float(line[142:153]))
                step.C_new.insert(0, float(line[153:164]))
                step.dC.insert(0, float(line[164:175]))
                step.dr.insert(0, float(line[175:186]))
                step.C_plm.insert(0, float(line[186:197]))
                step.u_plm.insert(0, float(line[197:208]))
                step.r_plm.insert(0, float(line[208:]))
            elif words[0] == "TOTAL":
                step.V_cavTot = float(line[25:36])
            elif words[0] == "BRINE" and words[1] == "OUT":
                step.Q_out = float(line[25:36])
            elif words[0] == "OUTLET":
                step.sg_out = float(line[25:36])
            elif words[0] == "VOLUME":
                step.V_insolTot = float(line[25:36])
            elif words[0] == "INSOL":
                step.z_insol = z_0[0] - float(line[25:36])
            elif words[0] == "BLANKET":
                step.z_obi = z_0[0] - float(line[25:36])
            elif words[0] == "VOL":
                step.V_insolVent = float(line[25:36])
            elif words[0] == "BRINE" and words[1] == "VOLUME":
                step.V_brine = float(line[25:36])
            elif words[0] == "ULLAGE":
                step.V_ullage = float(line[25:36])
            elif words[0] == "USEABLE":
                step.V_usable = float(line[25:36])
            elif words[0] == "CFAC":
                step.err = float(line[25:36])
            elif line.startswith("----------"):
                # TODO: add datablock to list
                steps.append(step)
                started = False
            else:
                print("unrecognized line: ", line)
    combo = _OutputData()
    combo.z_0 = z_0
    combo.h_0 = h_0
    combo.r_0 = r_0
    combo.V_injTot = list()
    combo.V_fillTot = list()
    combo.sg_cavAve = list()

    for step in steps:
        combo.V_injTot.append(np.nan)
        combo.V_fillTot.append(np.nan)
        combo.sg_cavAve.append(np.nan)
        for key in step.__dict__.keys():
            if not hasattr(combo, key):
                continue
            else:
                tmp = getattr(combo, key)
                if tmp is None:
                    setattr(combo, key, list())
                if key == "t":
                    val = step.t + step.t_prime
                else:
                    val = getattr(step, key)
                getattr(combo, key).append(val)

    res = Results(combo)
    return res


def read_tst_file(file_prefix: str):
    """Read a .tst output file into a DataFrame.

    Parameters
    ----------
    file_prefix : str
        Filename without extension.
    """
    Q_fill = list()
    t = list()
    V_cavTot = list()
    err = list()
    sg_out = list()
    sg_cavAve = list()
    V_insolTot = list()
    z_insol = list()
    z_obi = list()
    V_insolVent = list()
    V_ullage = list()
    V_usable = list()
    Q_inj = list()
    Q_fill = list()
    V_injTot = list()
    V_fillTot = list()
    with open(file_prefix + ".tst", "r") as fout:
        for line in fout.readlines():
            words = line.split()
            if len(words) == 0:
                continue
            if words[0] in ["#", "t", "(d)", "File=", "TIME", "Days", "---", "==="]:
                continue
            t.append(float(words[0]))
            V_cavTot.append(float(words[1]))
            err.append(float(words[2]))
            sg_out.append(float(words[3]))
            sg_cavAve.append(float(words[4]))
            V_insolTot.append(float(words[5]))
            z_insol.append(float(words[6]))
            z_obi.append(float(words[7]))
            V_insolVent.append(float(words[8]))
            V_ullage.append(float(words[9]))
            V_usable.append(float(words[10]))
            Q_inj.append(float(words[11]))
            V_injTot.append(float(words[12]))
            Q_fill.append(float(words[13]))
            V_fillTot.append(float(words[14]))
    return pd.DataFrame.from_dict(
        dict(
            t_d=t,
            V_cav=V_cavTot,
            err_ode=err,
            sg_out=sg_out,
            sg_ave=sg_cavAve,
            V_insol=V_insolTot,
            z_insol=z_insol,
            z_obi=z_obi,
            V_vented=V_insolVent,
            Q_inj=Q_inj,
            Q_fill=Q_fill,
            V_inj=V_injTot,
            V_fill=V_fillTot,
        )
    )


def write_hdf_results(results: Results, filename: str, **kwargs):
    """Write results to an HDF5 file.

    Parameters
    ----------
    results : Results
        The results object to write
    filename : str
        The filename to write the results to.

    Keyword Arguments
    -----------------
    kwargs : additional keyword arguments
        See :meth:`~sansmic.model.Results.to_hdf` for keyword arguments.

    """
    results.to_hdf(filename)


def read_hdf_results(filename: str) -> Results:
    """Read results from an HDF5 file.

    Parameters
    ----------
    filename : str
        Filename to read the results from.

    """
    results = Results.from_hdf(filename)
    return results


def write_json_results(results: Results, filename: str, **kwargs):
    """Write results to a JSON file.

    Parameters
    ----------
    results : Results
        The results object to write
    filename : str
        The filename to write the results to.

    Keyword Arguments
    -----------------
    kwargs : additional keyword arguments
        See :meth:`json.dump` for valid keyword arguments

    """
    with open(filename, "w") as f:
        json.dump(results.to_dict(), f, **kwargs)


def read_json_results(filename: str):
    """Read results from a JSON file.

    Parameters
    ----------
    filename : str
        Filename to read the results from.
    """
    with open(filename, "r") as f:
        d = json.load(f)
    return Results.from_dict(d)


def write_csv_results(results: Results, prefix: str):
    """Write results to CSV files.

    Parameters
    ----------
    results : Results
        The results to write out.
    prefix : str
        Write results files in CSV formats.
    """
    with open(prefix + "-summary.csv", "w") as f:
        results.df_t_1D.to_csv(f, lineterminator="\n", index=False)
    with open(prefix + "-radius.csv", "w") as f:
        df = results.radius
        df.index = results.depths
        df = df.T
        df.index = results.time
        df = df.T
        df.to_csv(f, lineterminator="\n")
    with open(prefix + "-density.csv", "w") as f:
        df = results.cell_sg
        df.index = results.depths
        df = df.T
        df.index = results.time
        df = df.T
        df.to_csv(f, lineterminator="\n")
    with open(prefix + "-wall-angle.csv", "w") as f:
        df = results.wall_angle
        df.index = results.depths
        df = df.T
        df.index = results.time
        df = df.T
        df.to_csv(f, lineterminator="\n")
    with open(prefix + "-dr_dt.csv", "w") as f:
        df = results.rate_of_change_in_radius
        df.index = results.depths
        df = df.T
        df.index = results.time
        df = df.T
        df.to_csv(f, lineterminator="\n")
    with open(prefix + "-dC_dt.csv", "w") as f:
        df = results.rate_of_change_in_sg
        df.index = results.depths
        df = df.T
        df.index = results.time
        df = df.T
        df.to_csv(f, lineterminator="\n")
    with open(prefix + "-dC_dz.csv", "w") as f:
        df = results.vertical_diffusion_rate
        df.index = results.depths
        df = df.T
        df.index = results.time
        df = df.T
        df.to_csv(f, lineterminator="\n")
