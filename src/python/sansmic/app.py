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

import click
from numpy import round
from pip._vendor.rich.progress import Progress
import sansmic.io

try:  # pragma: no cover
    import h5py
except ImportError as e:
    h5py = e

logging.basicConfig()
logger = logging.getLogger("sansmic")


def print_license(ctx: click.Context, param, value):
    """
    Print the license for sansmic.

    This is a click option callback function.
    """
    if not value or ctx.resilient_parsing:
        return
    from sansmic import __license__

    click.echo_via_pager(__license__)
    ctx.exit()


def print_copyright(ctx: click.Context, param, value):
    """
    Print the copyright for sansmic.

    This is a click option callback function.
    """
    if not value or ctx.resilient_parsing:
        return
    from sansmic import __copyright__

    click.echo_via_pager(__copyright__)
    ctx.exit()


@click.group()
def cli():
    pass


@click.command()
@click.argument("scenario_file", type=click.Path(exists=True))
@click.option(
    "-o",
    "--prefix",
    default=None,
    help="Prefix for output files.  [default: from SCENARIO_FILE]",
)
@click.option(
    "--toml/--no-toml",
    default=None,
    help="Create a TOML scenario file.  [default: iff SCENARIO_FILE is .DAT format or PREFIX is set]",
)
@click.option(
    "--csv/--no-csv", default=True, help="Create CSV output files.", show_default=True
)
@click.option(
    "--hdf/--no-hdf",
    default=False,
    help=(
        "Create an HDF5 output file.  [default: no-hdf]"
        if not isinstance(h5py, Exception)
        else "[DISABLED: h5py not found]"
    ),
)
@click.option(
    "--tst/--no-tst",
    default=True,
    help="Create an old-style TST file.",
    show_default=True,
)
@click.option(
    "--json/--no-json",
    default=False,
    help="Create a JSON output data file.",
    show_default=True,
)
@click.option(
    "--oldout/--no-oldout",
    default=False,
    help="Create an old-style OUT file.",
    show_default=True,
)
@click.option(
    "-v", "--verbose", count=True, help="Increase log details.", show_default=True
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    default=False,
    help="Suppress all screen output.",
    show_default=True,
)
@click.option(
    "--debug",
    is_flag=True,
    default=False,
    help="Output debug information.",
    show_default=True,
)
@click.version_option(package_name="sansmic", message="%(version)s")
@click.option(
    "--license",
    is_flag=True,
    callback=print_license,
    expose_value=False,
    is_eager=True,
    help="Show license and exit.",
)
@click.option(
    "--copyright",
    is_flag=True,
    callback=print_copyright,
    expose_value=False,
    is_eager=True,
    help="Print copyright and exit.",
)
def run(
    scenario_file,
    *,
    prefix=None,
    toml=None,
    csv=True,
    hdf=False,
    tst=True,
    oldout=False,
    verbose=0,
    quiet=False,
    debug=False,
    json=False,
):
    from sansmic import __version__

    if debug:
        logger.setLevel(logging.DEBUG)
        verbose = 3
    elif verbose == 2:
        logger.setLevel(logging.INFO)
    elif quiet:
        logger.setLevel(logging.ERROR)
    try:
        if prefix is None:
            prefix = splitext(scenario_file)[0]
        verbosity = verbose
        model = sansmic.io.read_scenario(scenario_file, warn=not quiet)
        logger.info("Successfully created scenario from {}".format(scenario_file))
        if toml is False:
            pass
        elif (toml is None and prefix is not None) or scenario_file.lower().endswith(
            ".dat"
        ):
            toml = True
        elif toml is None and prefix is None:
            toml = False
        if toml:
            sansmic.io.write_scenario(model, prefix + ".toml")
            logger.info("Wrote scenario to {}".format(prefix + ".toml"))
    except Exception as e:
        logger.critical(str(e))
        raise e

    logger.debug("Running simulation")
    with model.new_simulation(prefix, verbosity, tst, oldout) as sim:
        logger.debug("Created new simulation")
        if quiet or not verbose:
            logger.debug("Running in batch mode")
            if not quiet:
                click.echo(
                    "Running sansmic in batch mode from {}".format(repr(scenario_file))
                )
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
            click.echo(
                "Running sansmic scenario: {}".format(
                    scenario_file if not model.title else model.title
                )
            )
            stage = 0
            with Progress() as progress:
                if verbose >= 1:
                    task = progress.add_task("Progress...", total=n_steps)
                if verbose >= 2:
                    task_S = progress.add_task(
                        "[red]Stage {}...".format(stage + 1), total=stage_sizes[stage]
                    )
                for stage, step in sim.steps:
                    if last_stage != stage:
                        if stage >= len(model.stages):
                            if verbose >= 1:
                                progress.update(
                                    task,
                                    completed=n_steps,
                                )
                            if verbose >= 2:
                                progress.update(
                                    task_S,
                                    completed=n_steps,
                                )
                                print("All stages complete.")
                        else:
                            last_stage = stage
                            last_step = step
                            p_freq = day_size[stage]
                            if verbose >= 2:
                                progress.update(
                                    task_S,
                                    completed=n_steps,
                                )
                                task_S = progress.add_task(
                                    "[red]Stage {}...".format(stage + 1),
                                    total=stage_sizes[stage],
                                )
                    else:
                        if verbose >= 1 and (step - last_step) % p_freq == 0:
                            progress.update(
                                task,
                                advance=p_freq,
                            )
                        if verbose >= 2 and (step - last_step) % p_freq == 0:
                            progress.update(
                                task_S,
                                advance=p_freq,
                            )
        logger.debug("Simulation complete")
    res = sim.results
    if not quiet:
        click.echo("Initial and final results:")
        click.echo(
            (res.df_t_1D.iloc[[0, -1], [1, 3, 13, 15, 19, 20, 21, 26]]).to_string(
                index=False
            )
        )

    # Wrap the outputs in a try/except block
    try:
        if csv:
            sansmic.io.write_csv_results(res, prefix)
        if json:
            sansmic.io.write_json_results(res, prefix + ".json")
        if hdf:
            sansmic.io.write_hdf_results(res, prefix + ".h5")
        logger.debug("Sansmic complete")
    except Exception as e:  # pragma: no cover
        logger.critical("Error while writing results - some results may be missing")
        raise e

    return res


@click.command()
@click.argument("dat_file", type=click.Path(exists=True))
@click.argument("out_file", type=click.Path(), required=False)
@click.option(
    "--full",
    is_flag=True,
    default=False,
    help="Keep all entries even blanks.",
    show_default=True,
)
@click.version_option(package_name="sansmic", message="%(version)s")
@click.option(
    "--license",
    is_flag=True,
    callback=print_license,
    expose_value=False,
    is_eager=True,
    help="Show license and exit.",
)
@click.option(
    "--copyright",
    is_flag=True,
    callback=print_copyright,
    expose_value=False,
    is_eager=True,
    help="Print copyright and exit.",
)
def convert(dat_file, out_file=None, full=False):
    """Command line function to convert a scenario/dat to a new format."""
    logger.debug("Running sansmic-convert")
    if out_file is None:
        out_file = click.prompt("Please enter the name of the file to create")
    try:
        model = sansmic.io.read_scenario(dat_file)
        logger.debug("Successfully created scenario from {}".format(dat_file))
        sansmic.io.write_scenario(model, out_file, redundant=full)
    except Exception as e:
        logger.critical(str(e))
        raise e
    logger.debug("Successfully wrote scenario to {}".format(out_file))
