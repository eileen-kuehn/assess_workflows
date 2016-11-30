import os
import click
import logging

from assess_workflows.generic.structure import Structure
from assess_workflows.utils.utils import do_multicore
from gnmutils.sources.datasource import DataSource
from gnmutils.utils import relevant_directories
from utility.exceptions import ExceptionFrame
from utility.report import LVL


@click.group()
@click.option("--basepath", "basepath", multiple=False, required=True)
@click.option("--workflow-name", "workflow_name", multiple=False, required=True)
@click.option("--step", "step", default=1, multiple=False)
@click.pass_context
def cli(ctx, basepath, workflow_name, step):
    ctx.obj["structure"] = Structure(basepath=basepath, name=workflow_name, step=step)


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.option("--output_path", "output_path", multiple=False, required=True)
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def prepare_raw_data(ctx, paths, output_path, pcount):
    data = []
    for path in paths:
        # prepare data
        for folder, workernode_subdir, run_subdir, _ in relevant_directories(path):
            data.append({
                "path": os.path.join(os.path.join(folder, workernode_subdir), run_subdir),
                "output_path": output_path
            })
    if pcount > 1:
        do_multicore(
            count=pcount,
            target=_prepare_raw_data,
            data=data
        )
    else:
        _prepare_raw_data(data)


def _prepare_raw_data(kwargs):
    """
    Method processes raw data from a given path and saves results to a given output_path.

    :param path: the path to process
    :param output_path: the path to save the output to
    """
    with ExceptionFrame():
        path = kwargs.get("path", None)
        output_path = kwargs.get("output_path", None)
        data_source = DataSource.best_available_data_source()
        for job in data_source.jobs(
                source="raw", path=path, data_path=output_path, stateful=True):
            data_source.write_job(data=job, path=output_path)
        for traffic in data_source.traffics(
                source="raw", path=path, data_path=output_path, stateful=True):
            data_source.write_traffic(data=traffic, path=output_path)


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.option("--output-path", "output_path", multiple=False, required=True)
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def create_payloads(ctx, paths, output_path, pcount):
    pass

cli.add_command(prepare_raw_data)
cli.add_command(create_payloads)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix="DISS")
