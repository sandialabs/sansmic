# coding: utf-8
#
# Copyright (c) 2024 National Technology and Engineering Solutions of
# Sandia, LLC. Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.
#
# SPDX-License-Identifier: BSD-3-Clause.

"""
The SANSMIC main program.

This submodule provides the hook for the command line programs ``sansmic``
and ``sansmic-convert``.

"""

import logging
from argparse import ArgumentParser
from os.path import splitext

from numpy import round
from pip._vendor.rich.progress import Progress
from sansmic import __version__
from sansmic import io as sio

try:
    import h5py
except ImportError as e:
    h5py = e

logging.basicConfig()
logger = logging.getLogger("sansmic")


def _get_verbosity(args):
    """Set the logger level for sansmic."""
    level = logging.WARNING
    if args.verbose:
        level = logging.INFO
    elif args.quiet:
        level = logging.ERROR
    if args.debug:
        level = logging.DEBUG
    logger.setLevel(level)
    if args.debug:
        args.verbose = 3
        return 3
    return args.verbose


def _main_parser(default_hdf5=None):
    if default_hdf5 is None:
        kwargs = dict(default=None)
    else:
        kwargs = dict()
    parser = ArgumentParser(
        prog="sansmic",
        description=f"A solution mining code. v{__version__}",
    )
    parser.add_argument(
        "datafile",
        metavar="INPUTFILE",
        help="the SANSMIC input file to use; if extension is '.dat', --toml is implied",
    )

    outp = parser.add_argument_group(
        "Output options",
        """Default options:  ``{}  {}  --tst  --no-old-out``""".format(
            "--hdf" if default_hdf5 or default_hdf5 is None else "--no-hdf",
            "--json" if not default_hdf5 and default_hdf5 is not None else "--no-json",
        ),
    )
    outp.add_argument(
        "-o",
        "--output",
        dest="prefix",
        help="change output file prefix; if used, --toml is implied",
        default=None,
    )

    wwo = outp.add_mutually_exclusive_group()
    wwo.add_argument(
        "--toml", default=None, action="store_true", help="Create a TOML scenario file."
    )
    wwo.add_argument(
        "--no-toml",
        dest="toml",
        action="store_false",
        help="Don't create a TOML scenario file.",
        **kwargs,
    )

    wwo = outp.add_mutually_exclusive_group()
    wwo.add_argument(
        "--hdf",
        default=default_hdf5,
        action="store_true",
        help="Create an HDF5 formatted results file.",
    )
    wwo.add_argument(
        "--no-hdf",
        dest="hdf",
        action="store_false",
        help="Don't create an HDF5 results file.",
        **kwargs,
    )

    wwo = outp.add_mutually_exclusive_group()
    wwo.add_argument(
        "--json",
        default=not default_hdf5 if default_hdf5 is not None else None,
        action="store_true",
        help="Create a JSON formatted results file.",
    )
    wwo.add_argument(
        "--no-json",
        dest="json",
        action="store_false",
        help="Don't create a JSON results file.",
        **kwargs,
    )

    wwo = outp.add_mutually_exclusive_group()
    wwo.add_argument(
        "--tst",
        default=True,
        action="store_true",
        help="Create a daily summary TST text file.",
    )
    wwo.add_argument(
        "--no-tst",
        dest="tst",
        action="store_false",
        help="Don't create a TST summary file.",
    )

    wwo = outp.add_mutually_exclusive_group()
    wwo.add_argument(
        "--old-out",
        dest="old_out",
        action="store_true",
        default=False,
        help="Create an old-style SANSMIC OUT file.",
    )
    wwo.add_argument(
        "--no-old-out",
        dest="old_out",
        action="store_false",
        help="Don't create an old-style OUT file.",
    )

    outp = parser.add_argument_group(
        "Console/stdout/stderr reporting",
        "The following options are mutually exclusive",
    )
    verb = outp.add_mutually_exclusive_group()
    verb.add_argument(
        "-v",
        "--verbose",
        help="increase reporting details",
        action="count",
        default=0,
    )
    verb.add_argument(
        "-q",
        "--quiet",
        help="turn off runtime screen output",
        action="store_true",
        default=False,
    )
    verb.add_argument(
        "--debug",
        help="turn on debug messages and set max verbosity",
        action="store_true",
        default=False,
    )
    return parser


def main(args=None, ret=False):
    """Command line function to run sansmic.

    If this function is run from somewhere other than the command line, then
    arguments can be passed in to control execution as a list of strings.
    The ``ret`` argument can be set to ``True`` to return the results when
    using this as a function, also.

    Parameters
    ----------
    args : argparse.Namespace or list[str]
        Arguments as a list of strings that would be used on the command-line.
    ret : bool
        Should the function return a results object, by default False.

    """
    if ret or args is not None:
        default_hdf5 = False
    else:
        default_hdf5 = not isinstance(h5py, Exception)
    parser = _main_parser(default_hdf5)
    args = parser.parse_args(args=args)
    datafile = args.datafile
    prefix = splitext(datafile)[0] if args.prefix is None else args.prefix
    verbosity = _get_verbosity(args)
    logger.debug("Running sansmic with {}".format(args))
    model = sio.read_scenario(datafile, warn=not args.quiet)
    logger.info("Successfully created scenario from {}".format(datafile))
    if (args.toml is None and args.prefix is not None) or datafile.lower().endswith(
        ".dat"
    ):
        args.toml = True
    elif args.toml is None and args.prefix is None:
        args.toml = False
    if args.toml:
        sio.write_scenario(model, prefix + ".toml")
    logger.debug("Running simulation")
    with model.new_simulation(prefix, verbosity, args.tst, args.old_out) as sim:
        logger.debug("Created new simulation")
        if args.quiet or not args.verbose:
            logger.debug("Running in batch mode")
            if not args.quiet:
                print("Running sansmic in batch mode from ".format(repr(datafile)))
            sim.run_sim()
        else:
            logger.debug("Running in stepwise mode")
            last_stage = 0
            last_step = 0
            n_stages = len(model.stages)
            dt_stage = [s.solver_timestep for s in model.stages]
            day_size = [int(round(24.0 / dt)) for dt in dt_stage]
            t_stage = [s.injection_duration + s.rest_duration for s in model.stages]
            t_inject = [s.injection_duration for s in model.stages]
            stage_sizes = [
                int(round(t_stage[ct] / dt_stage[ct])) for ct in range(n_stages)
            ]
            n_steps = sum(stage_sizes)
            p_freq = day_size[0]
            print(
                "Running sansmic scenario: {}".format(
                    datafile if not model.title else model.title
                )
            )
            stage = 0
            with Progress() as progress:
                if args.verbose >= 1:
                    task = progress.add_task("Progress...", total=n_steps)
                if args.verbose >= 2:
                    task_S = progress.add_task(
                        "[red]Stage {}...".format(stage + 1), total=stage_sizes[stage]
                    )
                for stage, step in sim.steps:
                    if last_stage != stage:
                        if stage >= len(model.stages):
                            if args.verbose >= 1:
                                progress.update(
                                    task,
                                    completed=n_steps,
                                )
                            if args.verbose >= 2:
                                progress.update(
                                    task_S,
                                    completed=n_steps,
                                )
                                print("All stages complete.")
                        else:
                            last_stage = stage
                            last_step = step
                            p_freq = day_size[stage]
                            if args.verbose >= 2:
                                progress.update(
                                    task_S,
                                    completed=n_steps,
                                )
                                task_S = progress.add_task(
                                    "[red]Stage {}...".format(stage + 1),
                                    total=stage_sizes[stage],
                                )
                    else:
                        if args.verbose >= 1 and (step - last_step) % p_freq == 0:
                            progress.update(
                                task,
                                advance=p_freq,
                            )
                        if args.verbose >= 2 and (step - last_step) % p_freq == 0:
                            progress.update(
                                task_S,
                                advance=p_freq,
                            )
        logger.debug("Simulation complete")
    res = sim.results
    if not args.quiet:
        print("Initial and final results:")
        print(
            (res.df_t_1D.iloc[[0, -1], [1, 3, 13, 15, 19, 20, 21, 26]]).to_string(
                index=False
            )
        )
    if args.hdf:
        sio.write_hdf(res, prefix + ".h5")
    if args.json:
        sio.write_json(res, prefix + ".json")
    logger.debug("Sansmic complete")
    if ret:
        return res


def _convert_parser():
    parser = ArgumentParser(
        prog="sansmic-convert",
        description=f"Convert old SANSMIC DAT files to new sansmic scenario files. v{__version__}",
    )
    parser.add_argument(
        "infile", metavar="OLD_FILE", help="the SANSMIC input file to convert"
    )
    parser.add_argument(
        "outfile",
        metavar="NEW_FILE",
        help="the new scenario file to create [extension choices: .toml, .json, .yaml]",
        default=None,
    )
    parser.add_argument(
        "-f",
        "--full",
        action="store_true",
        default=False,
        help="Keep all entries, even if they are blank",
    )
    return parser


def convert():
    """Command line function to convert a scenario/dat to a new format."""

    parser = _convert_parser()
    args = parser.parse_args()
    infile = args.infile
    logger.debug("Running sansmic-convert")
    model = sio.read_scenario(infile)
    logger.debug("Successfully created scenario from {}".format(infile))
    sio.write_scenario(model, args.outfile, redundant=args.full)
    logger.debug("Successfully wrote scenario to {}".format(args.outfile))
