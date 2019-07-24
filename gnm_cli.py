import os
import glob
import click
import logging

from assess_workflows.generic.structure import Structure
from assess_workflows.utils.utils import do_multicore
from gnmutils.pilot import Pilot
from gnmutils.sources.datasource import DataSource
from gnmutils.utils import relevant_directories


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
        # TODO: is this called for every filename?!
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
        _prepare_raw_data(data[0])


@click.command()
@click.option("--paths", "paths", multiple=True, required=True)
@click.option("--output_path", "output_path", multiple=False, required=True)
@click.option("--pcount", "pcount", type=int, default=1)
@click.pass_context
def create_payloads(ctx, paths, output_path, pcount):
    data = []
    for path in paths:
        # prepare data
        for folder, workernode_subdir, run_subdir, _ in relevant_directories(path):
            # get all relevant files
            current_path = os.path.join(os.path.join(folder, workernode_subdir), run_subdir)
            data.extend([{
                "path": filename,
                "output_path": output_path
            } for filename in glob.glob("%s/*-process.csv" % current_path)])
    if pcount > 1:
        do_multicore(
            count=pcount,
            target=_create_payloads,
            data=data
        )
    else:
        for element in data:
            _create_payloads(element)


def _prepare_raw_data(kwargs):
    """
    Method processes raw data from a given path and saves results to a given output_path.

    :param path: the path to process
    :param output_path: the path to save the output to
    """
    path = kwargs.get("path", None)
    output_path = kwargs.get("output_path", None)
    data_source = DataSource.best_available_data_source()
    for job in data_source.jobs(
            source="raw", path=path, data_path=output_path, stateful=False):
        data_source.write_job(data=job, path=output_path)
    for traffic in data_source.traffics(
            source="raw", path=path, data_path=output_path, stateful=False):
        data_source.write_traffic(data=traffic, path=output_path)


def _create_payloads(kwargs):
    path = kwargs.get("path", None)
    output_path = kwargs.get("output_path", None)
    data_source = DataSource.best_available_data_source()
    for pilot in data_source.jobs(path=path):
        if pilot is not None:
            pilot.__class__ = Pilot
            if pilot.is_cms_pilot():
                pilot.prepare_traffic()
                for payload, _ in pilot.payloads():
                    # write file per payload
                    data_source.write_payload(path=output_path,
                                              data=payload)
                else:
                    logging.info("current pilot is not a CMS pilot")


cli.add_command(prepare_raw_data)
cli.add_command(create_payloads)

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.WARNING)
    logging.getLogger("EXCEPTION").setLevel(logging.INFO)
    cli(obj={}, auto_envvar_prefix="DISS")
