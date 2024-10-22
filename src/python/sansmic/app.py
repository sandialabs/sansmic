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
from datetime import datetime as dt
from pathlib import Path

import click
from numpy import round
from pip._vendor.rich.progress import Progress

import sansmic.io
from ._version import __version__, copyright, license

try:  # pragma: no cover
    import h5py
except ImportError as e:
    h5py = e

startup = dt.now().isoformat()
logger = logging.getLogger("sansmic")
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.DEBUG)


def print_license(ctx: click.Context, param, value):
    """
    Print the license for sansmic.

    This is a click option callback function.
    """
    if not value or ctx.resilient_parsing:
        return
    click.echo_via_pager(license)
    ctx.exit()


def print_copyright(ctx: click.Context, param, value):
    """
    Print the copyright for sansmic.

    This is a click option callback function.
    """
    if not value or ctx.resilient_parsing:
        return
    click.echo_via_pager(copyright)
    ctx.exit()


class _Echo:
    """Filter messages output using ``click.secho``.

    Output is determined by comparing the verbosity level and quiet-mode flag
    to the rules for for the method used.

    .. rubric:: Rules

    ``__call__(...)``
        Output message if not in quiet mode.
    ``verbose(...)``
        Output message only if verbose is non-zero.
    ``force(...)``
        Always output the message.

    The parameters for each of the calls is the same as :meth:`click.echo`.


    Parameters
    ----------
    verbosity : int
        The verbosity level, by default 0.
    quiet : bool
        Quiet-mode flag, by default False.
    """

    def __init__(self, verbosity: int = 0, quiet: bool = False):
        self.verbosity = verbosity
        self.quiet = quiet

    def __call__(
        self,
        message=None,
        file=None,
        nl: bool = True,
        err: bool = False,
        color: bool = None,
        **styles,
    ):
        """Output message unless :attr:`quiet` was set to True."""
        if self.quiet:
            return
        click.secho(message=message, file=file, nl=nl, err=err, color=color, **styles)

    def verbose(
        self,
        message=None,
        file=None,
        nl: bool = True,
        err: bool = False,
        color: bool = None,
        **styles,
    ):
        """Output message iff :attr:`verbose` >= 1."""
        if not self.verbosity:
            return
        click.secho(message=message, file=file, nl=nl, err=err, color=color, **styles)

    def force(
        self,
        message=None,
        file=None,
        nl: bool = True,
        err: bool = False,
        color: bool = None,
        **styles,
    ):
        """Output message even if :attr:`quiet` was set to True."""
        click.secho(message=message, file=file, nl=nl, err=err, color=color, **styles)


@click.group()
def cli():
    pass


@click.command("run")
@click.argument(
    "scenario_file",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
)
@click.option(
    "-o",
    "--output-name",
    "prefix",
    default=None,
    help="Filename stem for output files. Should NOT include a directory path (use -O to specify the output directory).  [default: stem of SCENARIO_FILE]",
)
@click.option(
    "-O",
    "--output-path",
    "opath",
    default=".",
    type=click.Path(
        exists=True,
        writable=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    help="Directory for output files. Directory must already exist.  [default: current directory]",
)
@click.option(
    "--csv/--no-csv",
    default=True,
    help="Create CSV output files.",
    show_default=True,
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
    "-v",
    "--verbose",
    count=True,
    help="Increase log details.",
    show_default=True,
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
    scenario_file: Path,
    *,
    prefix: str = None,
    opath: Path = Path("."),
    toml: bool = True,
    csv: bool = True,
    hdf: bool = False,
    tst: bool = True,
    oldout: bool = False,
    verbose: int = 0,
    quiet: bool = False,
    debug: bool = False,
    json: bool = False,
):

    echo = _Echo(verbosity=verbose, quiet=quiet)
    echo(f"sansmic v{__version__}")

    if prefix is None:
        prefix = scenario_file.stem
    pprefix = opath.joinpath(prefix)

    log_file = opath.joinpath(prefix).with_suffix(".log")
    with open(log_file, "w") as flog:
        flog.write(f"program:  sansmic v{__version__}\n")
        flog.write(f"startup:  {startup}\n")
        flog.write(f"input:    {scenario_file}\n")
        flog.write(f"output:   {pprefix}\n")
        flog.write("run-log:\n")

    file_log = logging.FileHandler(log_file)
    file_log.setLevel(logging.INFO)
    if debug:
        file_log.setFormatter(
            logging.Formatter(
                "- time: %(asctime)s\n  level: %(levelname)s\n  file: %(filename)s:%(lineno)d\n  funcName: %(funcName)s\n  message: %(message)s",
                datefmt="%Y-%m-%dT%H:%M:%S",
            )
        )
    elif verbose < 2:
        file_log.setFormatter(logging.Formatter("- %(message)s"))
    else:
        file_log.setFormatter(logging.Formatter("- message: %(message)s"))

    if debug:
        file_log.setLevel(logging.DEBUG)
        verbose = 4
        echo.verbosity = 3
    elif verbose > 2:
        file_log.setLevel(logging.DEBUG)
        echo.verbosity = 2
    elif verbose == 2:
        file_log.setLevel(logging.INFO)
        echo.verbosity = 2
    elif verbose == 1:
        file_log.setLevel(logging.INFO)
    elif quiet:
        file_log.setLevel(logging.WARNING)

    logger.addHandler(file_log)

    toml_file = opath.joinpath(prefix).with_suffix(".toml")
    if toml_file.exists() and toml_file == scenario_file:
        toml = False

    try:
        model = sansmic.io.read_scenario(scenario_file, warn=not quiet)
        echo.verbose(f'Loaded scenario "{scenario_file}"')
        if toml:
            sansmic.io.write_scenario(model, toml_file)
            echo.verbose(f'Wrote scenario "{toml_file}"')
        file_log.flush()
    except Exception as e:
        logger.critical(str(e))
        file_log.flush()
        file_log.stream.close()
        raise e

    with model.new_simulation(str(pprefix.absolute()), verbose, tst, oldout) as sim:
        logger.debug("Created new simulation object")
        if verbose >= 2 and not debug:
            file_log.setFormatter(
                logging.Formatter(
                    "- message: %(message)s\n  elapsed: %(relativeCreated).0f ms"
                )
            )
        logger.info("Running simulation")
        if not verbose:
            echo("Running sansmic simulation...")
            sim.run_sim()
        else:
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
            echo(
                "Running {}".format(
                    "sansmic simulation..." if not model.title else model.title
                )
            )
            stage = 0
            with Progress() as progress:
                file_log.flush()
                if verbose >= 1:
                    task = progress.add_task("Progress...", total=n_steps)
                if verbose >= 2:
                    logger.info(f"Starting stage {stage+1}")
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
                                logger.info(f"Stage {stage} complete")
                                progress.update(
                                    task_S,
                                    completed=n_steps,
                                )
                            logger.info("All stages complete")
                            file_log.flush()
                        else:
                            last_stage = stage
                            last_step = step
                            p_freq = day_size[stage]
                            if verbose >= 2:
                                progress.update(
                                    task_S,
                                    completed=n_steps,
                                )
                                logger.info(f"Stage {stage} complete")
                                logger.info(f"Starting stage {stage+1}")
                                file_log.flush()
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
        echo("Simulation complete.")
    logger.info("Simulation complete - writing results")
    file_log.flush()

    res = sim.results
    echo("Starting and ending cavern state: ")
    echo(
        (res.df_t_1D.iloc[[0, -1], [1, 3, 13, 15, 19, 20, 21, 26]]).to_string(
            index=False
        )
    )

    # Wrap the outputs in a try/except block
    try:
        if csv:
            sansmic.io.write_csv_results(res, opath.joinpath(prefix))
            logger.info(f'Wrote results to "{prefix}.*.csv"')
            echo("CSV result files written to {}.*.csv".format(opath.joinpath(prefix)))
        if json:
            sansmic.io.write_json_results(
                res, opath.joinpath(prefix).with_suffix(".json")
            )
            logger.info(f'Wrote results to "{prefix}.json"')
            echo(
                "JSON results file written to {}".format(
                    opath.joinpath(prefix).with_suffix(".json")
                )
            )
        if hdf:
            sansmic.io.write_hdf_results(res, opath.joinpath(prefix).with_suffix(".h5"))
            logger.info(f'Wrote results to "{prefix}.h5"')
            echo(
                "HDF5 results file written to {}".format(
                    opath.joinpath(prefix).with_suffix(".h5")
                )
            )
    except Exception as e:  # pragma: no cover
        logger.critical("Error while writing results - some results may be missing")
        file_log.stream.close()
        raise e
    logger.info("All processes complete")
    file_log.flush()
    file_log.stream.write(f"shutdown: {dt.now().isoformat()}\n")
    file_log.stream.close()
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


if __name__ == "__main__":
    run()
