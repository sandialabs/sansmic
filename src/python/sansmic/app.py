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

from sansmic import __version__
from sansmic import io as sio
from pip._vendor.rich.progress import Progress

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


def main(args=None):
    """Command line function to run sansmic."""
    parser = ArgumentParser(
        prog="sansmic",
        description=f"A solution mining code. v{__version__}",
    )
    parser.add_argument(
        "datafile", metavar="INPUTFILE", help="the SANSMIC input file to use"
    )
    outp = parser.add_argument_group(
        "Output options",
        """By default, only warnings will be output to the console, using
--quiet will turn these off. Using --verbose will change the behavior up to
three times. 1: turn on logging INFO level output. 2: output a summary of
results at the end of each stage, and use a progress bar. 3: output summary
results every simulation day. Using --debug will set the logging level to
DEBUG instead of INFO and automatically set verbosity to 3. These options
are mutually exclusive.
""",
    )
    outp.add_argument(
        "-o",
        "--output",
        dest="prefix",
        help="change result files' prefix [default: from datafile]",
        default=None,
    )
    verb = outp.add_mutually_exclusive_group()
    verb.add_argument(
        "-v",
        "--verbose",
        help="turn on info-level messages, increase the output verbosity up to 3 times",
        action="count",
        default=0,
    )
    verb.add_argument(
        "-q",
        "--quiet",
        help="set log level to errors-only",
        action="store_true",
        default=False,
    )
    verb.add_argument(
        "--debug",
        help="turn on debug-level messages and set verbosity to 3",
        action="store_true",
        default=False,
    )
    args = parser.parse_args(args=args)
    datafile = args.datafile
    prefix = splitext(datafile)[0] if args.prefix is None else args.prefix
    verbosity = _get_verbosity(args)
    logger.debug("Running sansmic with {}".format(args))
    model = sio.read_scenario(datafile, warn=not args.quiet)
    logger.info("Successfully created scenario from {}".format(datafile))
    logger.debug("Running simulation")
    with model.new_simulation(prefix, verbosity) as sim:
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
            p_freq = model.stages[last_stage].save_frequency
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
                            p_freq = model.stages[last_stage].save_frequency
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
        print("Final results:")
        print(
            (res.summary.iloc[[0, -1], [1, 3, 13, 15, 19, 20, 21, 26]]).to_markdown(
                index=False
            )
        )
        summary = res.summary.to_markdown(index=False)
        print(res.summary)
    logger.debug("Sansmic complete")


def convert():
    """Command line function to convert a scenario/dat to a new format."""
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
    verb = parser.add_argument_group(
        "Output options",
        """By default, only warnings will be output to the console
        The following options are mutually exclusive.""",
    )
    verb = verb.add_mutually_exclusive_group()
    verb.add_argument(
        "-v",
        "--verbose",
        help="turn on info-level messages, increase the output verbosity up to 3 times",
        action="count",
        default=0,
    )
    verb.add_argument(
        "-q",
        "--quiet",
        help="set log level to errors-only",
        action="store_true",
        default=False,
    )
    verb.add_argument(
        "--debug",
        help="turn on debug-level messages and set verbosity to 3",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()
    if args.debug and args.quiet:
        parser.error("argument --debug: not allowed with argument -q/--quiet")
    infile = args.infile
    verbosity = _get_verbosity(args)
    logger.debug("Running sansmic-convert")
    model = sio.read_scenario(infile, warn=not args.quiet)
    logger.debug("Successfully created scenario from {}".format(infile))
    sio.write_scenario(model, args.outfile, redundant=args.full)
    logger.debug("Successfully wrote scenario to {}".format(args.outfile))
