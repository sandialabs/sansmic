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
from argparse import ArgumentParser, Action
from os.path import splitext
import sys

from numpy import round
from pip._vendor.rich.progress import Progress
import sansmic.io

try:
    import h5py
except ImportError as e:
    h5py = e

logging.basicConfig()
logger = logging.getLogger("sansmic")


class LicenseAction(Action):
    def __call__(default, parser, namespace, values, option_string=None):
        if values != "bar":
            from sansmic import __license__

            print(__license__)
            parser.exit(0)


class CopyrightAction(Action):
    def __call__(default, parser, namespace, values, option_string=None):
        if values != "bar":
            from sansmic import __copyright__

            print(__copyright__)
            parser.exit(0)


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


def _main_parser(defaults=False):
    from sansmic import __version__

    if not defaults:
        kwargs = dict(default=None)
    else:
        kwargs = dict()
    parser = ArgumentParser(
        prog=f"sansmic v.{__version__}",
        description="Simulate leaching in an underground salt cavern due to ordinary solution mining, product withdrawal, or product fill.",
        epilog=f"sansmic (c) 2024 NTESS. Use the --copyright or --license flags for full details.",
    )
    parser.add_argument(
        "--copyright",
        default=None,
        nargs=0,
        action=CopyrightAction,
        help="Display full copyright details and exit.",
    )
    parser.add_argument(
        "--license",
        default=None,
        nargs=0,
        action=LicenseAction,
        help="Display full license details and exit.",
    )
    parser.add_argument(
        "datafile",
        metavar="INPUTFILE",
        help="the SANSMIC input file to use; if extensansmic.ion is '.dat', --toml is implied",
    )

    outp = parser.add_argument_group(
        "Output options",
        """Default options:  ``--csv  --tst  --no-hdf  --no-json  --no-old-out``""",
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
        "--toml",
        default=None,
        action="store_true",
        help="Create a TOML scenario file.",
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
        "--csv",
        action="store_true",
        default=True if defaults else None,
        help="Create CSV results files.",
    )
    wwo.add_argument(
        "--no-csv",
        dest="csv",
        action="store_false",
        help="Don't create CSV results files.",
        **kwargs,
    )

    wwo = outp.add_mutually_exclusive_group()
    wwo.add_argument(
        "--hdf",
        default=False if defaults else None,
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
        default=False if defaults else None,
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
        default=True if defaults else None,
        action="store_true",
        help="Create a daily summary TST text file.",
    )
    wwo.add_argument(
        "--no-tst",
        dest="tst",
        action="store_false",
        help="Don't create a TST summary file.",
        **kwargs,
    )

    wwo = outp.add_mutually_exclusive_group()
    wwo.add_argument(
        "--old-out",
        dest="old_out",
        action="store_true",
        default=False if defaults else None,
        help="Create an old-style SANSMIC OUT file.",
    )
    wwo.add_argument(
        "--no-old-out",
        dest="old_out",
        action="store_false",
        help="Don't create an old-style OUT file.",
        **kwargs,
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
    extra_args = args is not None
    if ret or args is not None:
        defaults = False
    else:
        defaults = True
    parser = _main_parser(defaults)
    args = parser.parse_args(args=args)

    # Wrap the different sections in try/except blocks
    try:
        datafile = args.datafile
        prefix = splitext(datafile)[0] if args.prefix is None else args.prefix
        verbosity = _get_verbosity(args)
        logger.debug("Running sansmic with {}".format(args))
        model = sansmic.io.read_scenario(datafile, warn=not args.quiet)
        logger.info("Successfully created scenario from {}".format(datafile))
        if not args.toml:
            pass
        elif (
            args.toml is None and args.prefix is not None
        ) or datafile.lower().endswith(".dat"):
            args.toml = True
        elif args.toml is None and args.prefix is None:
            args.toml = False
        if args.toml:
            sansmic.io.write_scenario(model, prefix + ".toml")
    except Exception as e:
        if extra_args:
            raise e
        parser.error(str(e))

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

    # Wrap the outputs in a try/except block
    try:
        if args.csv:
            sansmic.io.write_csv_results(res, prefix)
        if args.json:
            sansmic.io.write_json_results(res, prefix + ".json")
        if args.hdf:
            sansmic.io.write_hdf_results(res, prefix + ".h5")
        logger.debug("Sansmic complete")
    except Exception as e:
        if extra_args:
            raise e
        logger.critical("Error while writing results - some results may be missing")
        parser.error(str(e))

    if ret:
        return res


def _convert_parser():
    from sansmic import __version__

    parser = ArgumentParser(
        prog=f"sansmic-convert v.{__version__}",
        description="Convert from an old-style DAT file to the new TOML format.",
        epilog=f"sansmic (c) 2024 NTESS. Use the --copyright or --license flags for full details.",
    )
    parser.add_argument(
        "--copyright",
        default=None,
        nargs=0,
        action=CopyrightAction,
        help="Display full copyright details and exit.",
    )
    parser.add_argument(
        "--license",
        default=None,
        nargs=0,
        action=LicenseAction,
        help="Display full license details and exit.",
    )
    parser.add_argument(
        "infile", metavar="OLD_FILE", help="the SANSMIC input file to convert"
    )
    parser.add_argument(
        "outfile",
        metavar="NEW_FILE",
        help="the new scenario file to create [extensansmic.ion choices: .toml, .json, .yaml]",
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


def convert(args=None):
    """Command line function to convert a scenario/dat to a new format."""
    extra_args = args is not None
    parser = _convert_parser()
    args = parser.parse_args(args=args)
    if args.license:
        from sansmic import __license__

        print(__license__)
        parser.exit(0)
    elif args.copyright:
        from sansmic import __copyright__

        print(__copyright__)
        parser.exit(0)
    infile = args.infile
    logger.debug("Running sansmic-convert")
    try:
        model = sansmic.io.read_scenario(infile)
        logger.debug("Successfully created scenario from {}".format(infile))
        sansmic.io.write_scenario(model, args.outfile, redundant=args.full)
    except Exception as e:
        if extra_args:
            raise e
        parser.error(str(e))
    logger.debug("Successfully wrote scenario to {}".format(args.outfile))
