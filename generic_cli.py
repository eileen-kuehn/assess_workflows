"""
This file shall try to standardise my current workflows. It should provide a fixed structure for
contents that are created and used within the workflow.

To initialise the structure, the command init_environment can be used. It will create a structure
looking like this:

* base_folder
    * ./data
        * ./raw
        * ./processed
        * ./conferences <--- those should not be touched because they might already be published
    * ./figures
        * ./conferences <--- those should not be touched because they might already be published
        * ./exploratory
    * ./models
    * ./scripts
        * ./helper
            * ./conferences
        * ./md
        * ./python
        * ./R
    * ./text
    * ./workflows       <--- workflows are attached here...

There are specific workflow for tasks that are done. Each of those workflow gets its own folder
and also substructure. It should be self-contained by also has to be able to rely on the outer
folder structure, so that no raw or processed data needs to be copied. Only a specification of
which files to be used is needed here.

* ./workflows
    * ./<unique_workflow_id>
        * configuration.py <--- contains the configuration to run
        * run_workflow.sh  <--- contains the actual steps to be executed to get the results
        * README.md        <--- steps, thoughts, etc. are described in here
        * ./input          <--- either links or accumulated paths to data are stored in here
        * ./exploratory    <--- exploratory results are stored here
        * ./intermediate   <--- intermediate results are stored here
        * ./final          <--- final results are stored here

But how to integrate the workflows? How to integrate the actual research log? How to integrate
the stuff you are doing exploratory? Still no idea... but a first step taken here
"""
import os
import stat
import json
import click
import shutil
import logging
import datetime

from assess_workflows.generic.structure import Structure

from utility.exceptions import ExceptionFrame
from utility.report import LVL


@click.group()
@click.option("--basepath", "basepath", required=True, multiple=False)
@click.option("--workflow-name", "workflow_name", multiple=False, default=None)
@click.option("--step", "step", default=1, multiple=False)
@click.pass_context
def cli(ctx, basepath, workflow_name, step):
    ctx.obj["structure"] = Structure(basepath=basepath, name=workflow_name)


@click.command()
@click.option("--name", "name", default="another")
@click.pass_context
def create_workflow(ctx, name):
    workflow_name = "%s_%s_%s" % (datetime.date.today(), name, "workflow")
    structure = ctx.obj.get("structure", None)
    structure.name = workflow_name
    if not os.path.exists(structure.workflow_path()):
        try:
            structure.input_path()
            structure.exploratory_path()
            structure.intermediate_path()
            structure.final_path()
            with open(structure.workflow_config_file_path(), "w") as config_file:
                # create config object and write in json format
                config_object = {
                    "DISS_WORKFLOW_NAME": workflow_name
                }
                config_file.write(json.dumps(config_object))
            with open(structure.configuration_file_path(), "w") as configuration_file:
                configurations = [{"algorithms": [None], "signatures": [None]}]
                configuration_file.write("configurations = %s\n" % configurations)
            with open(structure.execution_file_path(), "w") as run_file:
                file_content = """#!/usr/local/bin/python

import os
import json
import click
import logging

from utility.exceptions import ExceptionFrame
from utility.report import LVL

from assess_workflows.generic.workflow import Workflow
from assess_workflows.generic.structure import Structure


@click.command()
@click.option("--start", "start", required=False, default=0)
@click.option("--end", "end", required=False, type=int)
def start(start, end):
    env_vars = json.load(open(os.path.join(
        os.path.join(os.path.dirname(__file__), os.pardir), "base_configuration.json"), "r"))
    env_vars.update(json.load(open("workflow_configuration.json", "r")))
    workflow = Workflow()
    # TODO: add your tasks to the workflow
    workflow.execute(env_vars, start=start, end=end)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger('EXCEPTION').setLevel(LVL.INFO)
    with ExceptionFrame():
        start()
"""
                run_file.write(file_content)
            # make file executable
            st = os.stat(structure.execution_file_path())
            os.chmod(structure.execution_file_path(), st.st_mode | stat.S_IEXEC)
        except OSError:
            click.echo("The workflow has already been created before, skipping here...")
        else:
            click.echo("created workflow structure at %s" % structure.workflow_path(workflow_name))


@click.command()
@click.option("--from_step", "from_steps", required=True, multiple=True, type=int)
@click.option("--to_step", "to_steps", required=True, multiple=True, type=int)
@click.pass_context
def intermediate_as_input(ctx, from_steps, to_steps):
    structure = ctx.obj.get("structure", None)
    for from_step in from_steps:
        file_path = structure.intermediate_file_path(step=from_step)
        for to_step in to_steps:
            try:
                os.symlink(file_path, structure.input_file_path(step=to_step))
            except OSError:
                pass


@click.command()
@click.option("--from_step", "from_steps", required=True, multiple=True, type=int)
@click.option("--file_type", "file_type", type=str, default="json")
@click.pass_context
def finalise(ctx, from_steps, file_type):
    structure = ctx.obj.get("structure", None)
    for from_step in from_steps:
        file_path = structure.intermediate_file_path(step=from_step, file_type=file_type)
        shutil.copyfile(file_path, structure.final_file_path(step=from_step, file_type=file_type))
        # TODO: maybe write protect the file?

cli.add_command(create_workflow)
cli.add_command(intermediate_as_input)
cli.add_command(finalise)

if __name__ == '__main__':
    logging.getLogger().setLevel(LVL.WARNING)
    logging.getLogger("EXCEPTION").setLevel(LVL.INFO)
    with ExceptionFrame():
        cli(obj={}, auto_envvar_prefix="DISS")
