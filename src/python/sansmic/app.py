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

from sansmic import __version__, io

logging.basicConfig()
logger = logging.getLogger("sansmic")


def main():
    """Command line function to run sansmic."""
    parser = ArgumentParser(
        prog="sansmic",
        description=f"A solution mining code. v{__version__}",
    )
    parser.add_argument("datafile", help="the SANSMIC input file to use")
    parser.add_argument(
        "-o",
        "--output",
        dest="prefix",
        help="change result files' prefix [default: from datafile]",
        default=None,
    )
    verb = parser.add_mutually_exclusive_group()
    verb.add_argument(
        "-v",
        "--verbose",
        help="Increase the output verbosity.",
        action="count",
        default=0,
    )
    verb.add_argument(
        "-q",
        "--quiet",
        help="Only report errors, not warnings.",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()
    datafile = args.datafile
    prefix = splitext(datafile)[0] if args.prefix is None else args.prefix
    level = logging.WARNING
    if args.verbose >= 2:
        level = logging.DEBUG
    elif args.verbose:
        level = logging.INFO
    elif args.quiet:
        level = logging.ERROR
    logger.setLevel(level)
    logger.debug("Running sansmic with {}".format(args))
    model = io.read_scenario(datafile, warn=not args.quiet)
    logger.debug("Successfully created scenario from {}".format(datafile))
    logger.info("Running simulation")
    with model.new_simulation(prefix) as sim:
        sim.run()


def convert():
    """Command line function to convert a scenario/dat to a new format."""
    parser = ArgumentParser(
        prog="sansmic-convert",
        description=f"A solution mining code. v{__version__}",
    )
    parser.add_argument("infile", help="the original input file to convert")
    parser.add_argument(
        dest="outfile",
        help="the new scenario file to create",
        default=None,
    )
    verb = parser.add_mutually_exclusive_group()
    verb.add_argument(
        "-v",
        "--verbose",
        help="increase the output verbosity: once -> info, twice -> debug.",
        action="count",
        default=0,
    )
    verb.add_argument(
        "-q",
        "--quiet",
        help="only report errors, not warnings.",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()
    infile = args.infile
    level = logging.WARNING
    if args.verbose >= 2:
        level = logging.DEBUG
    elif args.verbose:
        level = logging.INFO
    elif args.quiet:
        level = logging.ERROR
    logger.setLevel(level)
    logger.debug("Running sansmic-convert")
    model = io.read_scenario(infile, warn=not args.quiet)
    logger.debug("Successfully created scenario from {}".format(infile))
    io.write_scenario(model, args.outfile)
    logger.debug("Successfully wrote scenario to {}".format(args.outfile))
